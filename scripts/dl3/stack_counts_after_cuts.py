import logging
from argparse import ArgumentParser

import numpy as np
from astropy import units as u
from astropy.table import Table

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(input_paths, output_path, norm):
    cuts_after_trigger = []
    cuts_after_gh = []
    cuts_after_gh_theta = []
    t_effective = []
    t_elapsed = []

    for path in input_paths:
        cuts = Table.read(path, path="cuts")
        cuts_after_trigger.append(cuts["after_trigger"])
        cuts_after_gh.append(cuts["after_gh"])
        cuts_after_gh_theta.append(cuts["after_gh_theta"])
        t_effective.append(cuts.meta["t_effective"])
        t_elapsed.append(cuts.meta["t_elapsed"])

    t_eff = u.Quantity(t_effective).reshape(-1, 1)
    if norm == "none":
        norm = u.Quantity(1)
    elif norm == "eff":
        norm = np.sum(t_eff)
    else:
        raise NotImplementedError(f"Unsupported norm {norm}")

    cuts_after_trigger = np.sum(
        np.array(cuts_after_trigger) / norm,
        axis=0,
    )
    cuts_after_gh = np.sum(
        np.array(cuts_after_gh) / norm,
        axis=0,
    )

    table = Table(
        {
            "after_trigger": cuts_after_trigger,
            "after_gh": cuts_after_gh,
            "center": cuts["center"],
            "high": cuts["high"],
            "low": cuts["low"],
        },
    )
    table.write(output_path, path="cuts", overwrite=True, serialize_meta=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-paths", required=True, nargs="+")
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("--norm", default="none")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.input_paths, args.output_path, args.norm)
