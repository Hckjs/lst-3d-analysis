# Code for bkgmodel adapted from Simone Mender
# https://github.com/cta-observatory/pybkgmodel/tree/add_exlcusion_region_method
import argparse
import logging
from pathlib import Path
import yaml
from rich import progress
import astropy.units as u
import numpy as np
import pandas as pd
import pickle
from astropy.coordinates import AltAz, Angle, EarthLocation, SkyCoord
from astropy.coordinates.erfa_astrom import ErfaAstromInterpolator, erfa_astrom
from gammapy.data import DataStore
from gammapy.maps import MapAxis, WcsGeom
from scriptutils.log import setup_logging
from scriptutils.bkg import ExclusionMapBackgroundMaker
import pickle

log = logging.getLogger(__name__)


def main():
    """
    Function running the entire background reconstruction procedure.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--cached-count-maps", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix", default="bkg")
    parser.add_argument("--dummy-output", required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    out = Path(args.output_dir)

    erfa_astrom.set(ErfaAstromInterpolator(300 * u.s))

    with open(args.config) as f:
        config = yaml.safe_load(f)
    e_binning = config["binning"]["energy"]
    fov_binning = config["binning"]["offset"]
    exclusion = config["exclusion"]
    matching = config["run_matching"]

    # TODO Define that properly somewhere
    location = EarthLocation.of_site("Roque de los Muchachos")
    e_reco = MapAxis.from_energy_bounds(
        u.Quantity(e_binning["min"]),
        u.Quantity(e_binning["max"]),
        e_binning["n_bins"],
        name="energy",
    )
    ds = DataStore.from_dir(args.input_dir)
    with open(args.cached_counts_maps, "rb") as f:
        cached_maps = pickle.load(f)

    # Select similar runs. This is only zenith right now
    # Ra/dec is assumed to be correclty incorporated earlier
    # Az is neglected for now, maybe thats fine for one LST
    # Time is neglected as well. Very different runs should be excluded before this
    # and there is no notion of MC-periods in LST so far
    criteria = pd.DataFrame(
        {
            "obs_id": ds.obs_table["OBS_ID"],
            "cos_zenith": np.cos(ds.obs_table["ZEN_PNT"].quantity),
        },
    )
    for obs_id in ds.obs_ids:
        cos_zenith_diff = np.abs(
            (criteria["cos_zenith"] - criteria["cos_zenith"][0]).values,
        )
        mask = cos_zenith_diff < matching["max_cos_zenith_diff"]
        selected_ids = criteria["obs_id"][mask].values
        # select fitting runs
        bkg_maker = ExclusionMapBackgroundMaker(
            e_reco,
            location,
            nbins=fov_binning["n_bins"],
            offset_max=u.Quantity(fov_binning["max"]),
            exclusion_radius=u.Quantity(exclusion["radius"]),
            exclusion_sources=[SkyCoord(**s) for s in exclusion["sources"]],
        )
        bkg_maker.run(ds, selected_ids, cached_maps=cached_maps)
        if config["hdu_type"] == "3D":
            bkg = bkg_maker.get_bg_3d()
        elif config["hdu_type"] == "2D":
            bkg = bkg_maker.get_bg_2d()
        else:
            raise NotImplementedError()
        bkg.write(
            out / f"{config['prefix']}_{obs_id}.fits.gz",
            overwrite=args.overwrite,
        )
    Path(args.dummy_output).touch()


if __name__ == "__main__":
    main()
