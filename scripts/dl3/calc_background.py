# Code for bkgmodel adapted from Simone Mender
# https://github.com/cta-observatory/pybkgmodel/tree/add_exlcusion_region_method
import argparse
import logging
import pickle
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import yaml
from astropy.coordinates import EarthLocation
from astropy.coordinates.erfa_astrom import ErfaAstromInterpolator, erfa_astrom
from gammapy.data import DataStore
from gammapy.maps import MapAxis
from regions import Regions

from scriptutils.bkg import ExclusionMapBackgroundMaker
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main():  # noqa
    """
    Function running the entire background reconstruction procedure.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-runs", required=True, nargs="+")
    parser.add_argument("--cached-maps", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--exclusion", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix", default="bkg")
    parser.add_argument("--dummy-output", required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    out = Path(args.output_dir)

    exclusion_regions = Regions.read(args.exclusion, format="ds9")
    erfa_astrom.set(ErfaAstromInterpolator(300 * u.s))

    with open(args.config) as f:
        config = yaml.safe_load(f)
    e_binning = config["binning"]["energy"]
    fov_binning = config["binning"]["offset"]
    matching = config["run_matching"]
    assert len(matching) == 1, "Got more than one matching key in config"
    match_on = list(matching.keys())[0]
    log.info(f"Selecting based on: {match_on}")

    # TODO Define that properly somewhere
    location = EarthLocation.of_site("Roque de los Muchachos")
    e_reco = MapAxis.from_energy_bounds(
        u.Quantity(e_binning["min"]),
        u.Quantity(e_binning["max"]),
        e_binning["n_bins"],
        name="energy",
    )
    ds = DataStore.from_events_files(args.input_runs)
    with open(args.cached_maps, "rb") as f:
        cached_maps = pickle.load(f)

    # Select similar runs. This is only zenith right now
    # Ra/dec is assumed to be correclty incorporated earlier
    # Az is neglected for now, maybe thats fine for one LST
    # Time is neglected as well. Very different runs should be excluded before this
    # and there is no notion of MC-periods in LST so far

    # looks like obs_table isbehaving inconsistenly
    zens = []
    decs = []
    # TODO Define somewhere and import
    location = EarthLocation.of_site("Roque de los Muchachos")
    for obs_id in ds.obs_ids:
        obs = ds.obs(obs_id)
        decs.append(obs.pointing.get_icrs().dec.deg)
        zens.append(
            90 - obs.pointing.get_altaz(location=location, obstime=obs.tmid).alt.deg,
        )
    criteria = pd.DataFrame(
        {
            "obs_id": ds.obs_ids,
            "zenith": zens,
            "cos_zenith": np.cos(np.deg2rad(zens)),
        },
    )

    log.info(f"Selection criteria: {criteria}")

    for i, obs_id in enumerate(ds.obs_ids):
        if match_on == "max_cos_zenith_diff":
            cos_zenith_diff = np.abs(
                (criteria["cos_zenith"] - criteria["cos_zenith"][i]).values,
            )
            mask = cos_zenith_diff < matching["max_cos_zenith_diff"]
            selected_ids = criteria["obs_id"][mask].values
        elif match_on == "zenith_bin_edges":
            # this is a bit of a waste creating the bins and indices for every obs...
            edges = matching["zenith_bin_edges"]
            # Select all runs ending up in the same bin
            idx = np.digitize(criteria["zenith"], edges)
            mask = idx == idx[i]
            selected_ids = criteria["obs_id"][mask].values
        elif match_on == "cos_zenith_bin_edges":
            # this is a bit of a waste creating the bins and indices for every obs...
            edges = matching["cos_zenith_bin_edges"]
            # Select all runs ending up in the same bin
            idx = np.digitize(criteria["cos_zenith"], edges)
            mask = idx == idx[i]
            selected_ids = criteria["obs_id"][mask].values
        # TODO n_zenith_bins? then I have bins with the same amount of runs...
        else:
            raise NotImplementedError()

        log.info(f"Selected to match {obs_id}: {len(selected_ids)} ({selected_ids})")
        log.info(f"Selected to match {criteria[mask]}")

        # select fitting runs
        bkg_maker = ExclusionMapBackgroundMaker(
            e_reco,
            location,
            exclusion_regions=exclusion_regions,
            nbins=fov_binning["n_bins"],
            n_offset_bins=fov_binning.get("n_offset_bins", 8),
            offset_max=u.Quantity(fov_binning["max"]),
        )
        bkg_maker.run(ds, selected_ids, cached_maps=cached_maps)
        if config["hdu_type"] == "3D":
            bkg = bkg_maker.get_bg_3d()
        elif config["hdu_type"] == "2D":
            bkg = bkg_maker.get_bg_2d()
        else:
            raise NotImplementedError()
        bkg.write(
            out / f"{config['prefix']}_{obs_id:05d}.fits.gz",
            overwrite=args.overwrite,
        )
    Path(args.dummy_output).touch()


if __name__ == "__main__":
    main()
