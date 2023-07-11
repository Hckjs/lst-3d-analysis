from gammapy.irf import Background3D
import astropy.units as u
import numpy as np
from gammapy.data import DataStore
from gammapy.analysis import Analysis, AnalysisConfig
from astropy.io import fits
from argparse import ArgumentParser
from astropy.table import Table

parser = ArgumentParser()
parser.add_argument("--hdu-index-path", required=True)
parser.add_argument("--bkg-dir", required=True) # needs to be relative(?)
parser.add_argument("--bkg-file", required=True) # Note that this does only allow a single stacked bkg right now like this
args = parser.parse_args()


def main(hdu_index_path, bkg_dir, bkg_file):
    t = Table.read(hdu_index_path)
    for i in np.unique(t["OBS_ID"]):
        t.add_row([i, "bkg", "bkg_3d", bkg_dir, bkg_file, "BACKGROUND", -1])
    print(t)
    t.write(hdu_index_path, overwrite=True)



if __name__ == "__main__":
    main(**vars(args))
