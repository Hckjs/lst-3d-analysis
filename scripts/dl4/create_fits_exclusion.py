import logging
from argparse import ArgumentParser

from gammapy.analysis import AnalysisConfig, Analysis
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

    region = Regions.read(args.input)

    config = AnalysisConfig.read(args.config)
    analysis = Analysis(config)
    geom = analysis._create_geometry()

    exclusion_mask = geom.region_mask(region, inside=False)
    exclusion_mask.write(args.output, overwrite=True)

