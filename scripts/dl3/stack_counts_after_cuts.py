import logging
from argparse import ArgumentParser

import numpy as np
from astropy import units as u
from astropy.table import Table

log = logging.getLogger(__name__)


def main(input_paths, output_path, norm):
    cuts_after_trigger = []
    cuts_after_gh = []
    intensities_after_trigger = []
    intensities_after_gh = []
    cog_after_trigger = []
    cog_after_gh = []
    t_effective = []
    t_elapsed = []

    for path in input_paths:
        cuts = Table.read(path, path="cuts")
        cuts_after_trigger.append(cuts["after_trigger"])
        cuts_after_gh.append(cuts["after_gh"])
        t_effective.append(cuts.meta["t_effective"])
        t_elapsed.append(cuts.meta["t_elapsed"])

        intensity = Table.read(path, path="intensity")
        intensities_after_trigger.append(intensity["after_trigger"])
        intensities_after_gh.append(intensity["after_gh"])

        cog = Table.read(path, path="cog")
        cog_after_trigger.append(cog["after_trigger"])
        cog_after_gh.append(cog["after_gh"])

    t_eff = u.Quantity(t_effective).reshape(-1, 1)
    if norm == "none":
        norm = u.Quantity(1)
    elif norm == "eff":
        norm = np.sum(t_eff)
    else:
        raise NotImplementedError(f"Unsupported norm {norm}")

    # e reco
    cuts_after_trigger = np.sum(
        np.array(cuts_after_trigger),
        axis=0,
    )
    cuts_after_gh = np.sum(
        np.array(cuts_after_gh),
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
        meta={"t_elapsed": sum(t_elapsed), "t_effective": sum(t_effective)},
    )
    table.write(output_path, path="cuts", overwrite=True, serialize_meta=True)

    # intensity
    intensities_after_trigger = np.sum(
        np.array(intensities_after_trigger),
        axis=0,
    )
    intensities_after_gh = np.sum(
        np.array(intensities_after_gh),
        axis=0,
    )
    table = Table(
        {
            "after_trigger": intensities_after_trigger,
            "after_gh": intensities_after_gh,
            "center": intensity["center"],
            "high": intensity["high"],
            "low": intensity["low"],
        },
        meta={"t_elapsed": sum(t_elapsed), "t_effective": sum(t_effective)},
    )
    table.write(output_path, path="intensity", append=True, serialize_meta=True)

    # cog
    cog_after_trigger = np.sum(
        np.array(cog_after_trigger),
        axis=0,
    )
    cog_after_gh = np.sum(
        np.array(cog_after_gh),
        axis=0,
    )
    table = Table(
        {
            "after_trigger": cog_after_trigger,
            "after_gh": cog_after_gh,
            "bin_centers": cog["bin_centers"],
        },
        meta={"t_elapsed": sum(t_elapsed), "t_effective": sum(t_effective)},
    )
    table.write(output_path, path="cog", append=True, serialize_meta=True)





if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-paths", required=True, nargs="+")
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("--norm", default="none")
    args = parser.parse_args()
    main(args.input_paths, args.output_path, args.norm)
