import logging
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.io import fits
from os.path import relpath

log = logging.getLogger(__name__)


def main(hdu_index_path, bkg_files):
    t = Table.read(hdu_index_path)
    ids = sorted(np.unique(t["OBS_ID"]))

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
    # remove old links in case of reexecution of bkg rules
    t.remove_rows(np.nonzero(t["HDU_TYPE"] == "bkg"))
    for i, b in zip(ids, bkg_files):
        bkg = Path(b)
        log.info(f"Linking {bkg} to obs id {i}")
        with fits.open(bkg) as f:
            t.add_row(
                [
                    i,
                    "bkg",
                    f[1].header['HDUCLAS4'].lower(),
                    relpath(bkg.parent, hdu_index_path.parent),
                    bkg.name,
                    "BACKGROUND",
                    -1,
                ],
            )
        log.info(t)
    t.write(hdu_index_path, overwrite=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--hdu-index-path", required=True)
    parser.add_argument("--bkg-files", required=True, nargs="+")
    parser.add_argument("--log-file")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
    )  # TODO setup logging without scriptutils...
    args = parser.parse_args()
    main(Path(args.hdu_index_path), args.bkg_files)
