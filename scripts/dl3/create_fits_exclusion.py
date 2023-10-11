import logging
from argparse import ArgumentParser

from astropy.coordinates import SkyCoord
from gammapy.analysis import AnalysisConfig
from gammapy.maps import WcsGeom
from regions import Regions

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    # Create geometry to convert to pixel region
    config = AnalysisConfig.read(args.config)
    wcs_geom_settings = config.datasets.geom.wcs
    geom_params = {}
    skydir_settings = wcs_geom_settings.skydir
    skydir = SkyCoord(
        skydir_settings.lon,
        skydir_settings.lat,
        frame=skydir_settings.frame,
    )
    geom_params["skydir"] = skydir
    geom_params["frame"] = skydir_settings.frame
    geom_params["binsz"] = wcs_geom_settings.binsize
    width = wcs_geom_settings.width.width.to("deg").value
    height = wcs_geom_settings.width.height.to("deg").value
    geom_params["width"] = (width, height)
    geom = WcsGeom.create(**geom_params)

    region = Regions.read(args.input)
    region = region.to_pixel(geom.wcs)
    region.write(args.output, format="fits", overwrite=True)
