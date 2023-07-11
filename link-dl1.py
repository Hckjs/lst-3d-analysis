from pathlib import Path
from argparse import ArgumentParser
from lstchain.io.io import dl1_params_lstcam_key
import pandas as pd
import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, angular_separation
from astropy.time import Time
from link_utils import build_altaz, get_pointings_of_irfs, sin_delta, cos_zenith, euclidean_distance


def get_irf_pointings(mc_dir):
    # Why is the parent needed?
    irf_paths = [p for p in mc_dir.glob("*_test.dl1.h5")]
    # The test suffix screws up the split at underscores, also there is an extra corsika 
    filelist = [p.name.split("node_")[1].split("_test.dl1.h5")[0] for p in irf_paths]
    irf_pointings: AltAz = get_pointings_of_irfs(filelist)
    irf_sindelta = sin_delta(irf_pointings)
    irf_coszenith = cos_zenith(irf_pointings)
    return irf_sindelta, irf_coszenith, irf_paths


def get_pointing_from_file(run):
    df = pd.read_hdf(run, dl1_params_lstcam_key)
    mean_alt = df.alt_tel.mean()
    mean_az = df.az_tel.mean()
    pointing = build_altaz(alt = mean_alt*u.rad, az=mean_az*u.rad)
    sindelta = sin_delta(pointing)
    coszenith = cos_zenith(pointing)
    return sindelta, coszenith

def main():
    parser = ArgumentParser()
    parser.add_argument("--dl1-dir", required=True)
    parser.add_argument("--test-link-dir", required=True)
    parser.add_argument("--mc-nodes-dir", required=True)
    args = parser.parse_args()

    dl1_files = Path(args.dl1_dir).glob("dl1*.h5")
    irf_sd, irf_cz, irf_paths = get_irf_pointings(Path(args.mc_nodes_dir))
    link_dir = Path(args.test_link_dir)

    print("dl1 files:", dl1_files)
    for f in dl1_files:
        sd, cz = get_pointing_from_file(f)
        nearest = euclidean_distance(
                    x1=sd,
                    y1=cz,
                    x2=irf_sd,
                    y2=irf_cz,
            ).argmin()
        link = link_dir / f.name
        print("Will link from", link, "to", irf_paths[nearest])
        if link.is_symlink():
            link.unlink()
        link.symlink_to(irf_paths[nearest].resolve())




if __name__ == "__main__":
    main()
