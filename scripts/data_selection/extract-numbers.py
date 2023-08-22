import logging
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from astropy.table import Table

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def to_num(n: int) -> str:
    return r"\num{" + str(n) + "}"


def main():
    parser = ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("output_directory")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(logfile=args.log_file, verbose=args.verbose)
    outdir = Path(args.output_directory)

    tbl = Table.read("build/dl1-datachecks-masked.h5")

    mask_source = tbl["mask_run_id"]
    mask_time_pedestal_pointing = tbl["mask_run_selection"]
    mask_ped_charge = tbl["mask_pedestal_charge"]
    mask_cosmics = mask_ped_charge & tbl["mask_cosmics"] & tbl["mask_cosmics_above"]

    with open(outdir / "runselection-01-observing-source.tex", "w") as f:
        f.write(to_num(np.count_nonzero(mask_source)))

    with open(outdir / "runselection-02-ok-during-timeframe.tex", "w") as f:
        f.write(to_num(np.count_nonzero(mask_time_pedestal_pointing)))

    with open(outdir / "runselection-03-pedestal-charge.tex", "w") as f:
        f.write(to_num(np.count_nonzero(mask_ped_charge)))

    with open(outdir / "runselection-04-cosmics.tex", "w") as f:
        f.write(to_num(np.count_nonzero(mask_cosmics)))


if __name__ == "__main__":
    main()
