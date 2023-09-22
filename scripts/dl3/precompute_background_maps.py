# Code for bkgmodel adapted from Simone Mender
# https://github.com/cta-observatory/pybkgmodel/tree/add_exlcusion_region_method
import argparse
import logging
import pickle

import astropy.units as u
import yaml
from astropy.coordinates import EarthLocation
from astropy.coordinates.erfa_astrom import ErfaAstromInterpolator, erfa_astrom
from gammapy.data import DataStore
from gammapy.maps import MapAxis
from regions import Regions

from scriptutils.bkg import ExclusionMapBackgroundMaker
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main():
    """
    Function running the entire background reconstruction procedure.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-runs", required=True, nargs="+")
    parser.add_argument("--config", required=True)
    parser.add_argument("--exclusion", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    erfa_astrom.set(ErfaAstromInterpolator(300 * u.s))

    with open(args.config) as f:
        config = yaml.safe_load(f)
    e_binning = config["binning"]["energy"]
    fov_binning = config["binning"]["offset"]
    exclusion_regions = Regions.parse(args.exclusion, format="ds9")

    # TODO Define that properly somewhere
    location = EarthLocation.of_site("Roque de los Muchachos")
    e_reco = MapAxis.from_energy_bounds(
        u.Quantity(e_binning["min"]),
        u.Quantity(e_binning["max"]),
        e_binning["n_bins"],
        name="energy",
    )

    ds = DataStore.from_events_files(args.input_runs)
    bkg_maker = ExclusionMapBackgroundMaker(
        e_reco,
        location,
        exclusion_regions=exclusion_regions,
        nbins=fov_binning["n_bins"],
        offset_max=u.Quantity(fov_binning["max"]),
    )
    cached_maps = bkg_maker.fill_all_maps(ds, None)
    with open(args.output, "wb") as f:
        pickle.dump(cached_maps, f)


if __name__ == "__main__":
    main()
