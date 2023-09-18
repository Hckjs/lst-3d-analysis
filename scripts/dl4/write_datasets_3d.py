import logging
from argparse import ArgumentParser

from gammapy.analysis import Analysis, AnalysisConfig

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(config, output):
    # Standard high-level interface stuff
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)
    analysis.get_observations()
    analysis.get_datasets()

    analysis.datasets.write(output, overwrite=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.config, args.output)
