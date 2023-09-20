import logging
from argparse import ArgumentParser

import matplotlib.pyplot as plt
from gammapy.estimators import FluxPoints

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(flux_points, output):
    flux_points = FluxPoints.read(flux_points, format="profile")
    fig = plt.figure()
    flux_points.plot()
    fig.savefig(output)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--flux-points", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.flux_points, args.output)
