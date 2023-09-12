import logging
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from astropy.table import Table

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


# TODO Logging (but its in the lstchain env, so just basic logs)
def main(hdu_index_path, bkg_dir, bkg_pattern):
    t = Table.read(hdu_index_path)
    ids = sorted(np.unique(t["OBS_ID"]))
    log.info(len(ids), "ids in file:", ids)

    bkg_files = list(Path(bkg_dir).glob(bkg_pattern))

    if len(bkg_files) == 1:
        bkg_files = bkg_files * len(ids)
    elif len(bkg_files) == len(ids):
        log.info(
            """
            Number of bkg files matching.
            Sorting them assuming that leads
            to the same order as the obs_ids""",
        )
        bkg_files = sorted(bkg_files)
    else:
        raise Exception(
            f"""
            Number of bkg files does not match.
            We have {len(ids)} obs,
            but {len(bkg_files)} bkg files""",
        )
    for i, bkg in zip(ids, bkg_files):
        log.info(f"Linking {bkg} to obs id {i}")
        t.add_row([i, "bkg", "bkg_3d", bkg_dir, bkg, "BACKGROUND", -1])
    log.info(t)
    t.write(hdu_index_path, overwrite=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--hdu-index-path", required=True)
    parser.add_argument("--bkg-dir", required=True)  # needs to be relative(?)
    parser.add_argument("--bkg-pattern", default="bkg*.fits.gz")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.hdu_index_path, args.bkg_dir, args.bkg_pattern)
