from argparse import ArgumentParser

import numpy as np
from astropy.table import Table

parser = ArgumentParser()
parser.add_argument("--hdu-index-path", required=True)
parser.add_argument("--bkg-dir", required=True)  # needs to be relative(?)
parser.add_argument("--bkg-file", nargs="+")
args = parser.parse_args()


# TODO Logging (but its in the lstchain env, so just basic logs)
def main(hdu_index_path, bkg_dir, bkg_file):
    t = Table.read(hdu_index_path)
    ids = sorted(np.unique(t["OBS_ID"]))
    print(len(ids), "ids in file:", ids)
    if len(bkg_file) == 1:
        bkg_file = bkg_file * len(ids)
    elif len(bkg_file) == len(ids):
        print(
            """
            Number of bkg files matching.
            Sorting them assuming that leads
            to the same order as the obs_ids""",
        )
        bkg_file = sorted(bkg_file)
    else:
        raise Exception(
            f"""
            Number of bkg files does not match.
            We have {len(ids)} obs,
            but {len(bkg_file)} bkg files""",
        )
    for i, bkg in zip(ids, bkg_file):
        print(f"Linking {bkg} to obs id {i}")
        t.add_row([i, "bkg", "bkg_3d", bkg_dir, bkg, "BACKGROUND", -1])
    print(t)
    t.write(hdu_index_path, overwrite=True)


if __name__ == "__main__":
    main(**vars(args))
