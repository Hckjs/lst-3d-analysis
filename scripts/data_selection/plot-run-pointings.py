import logging
from argparse import ArgumentParser

from astropy.table import Table
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main():
    pointings = Table.read(args.input_path)

    # TODO 2 plots (also azimuth)
    fig, ax = plt.subplots()
    ax.plot(pointings["run_id"], pointings["zen"], ls="")
    fig.savefig(args.output_path)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("-o", "--output_path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main()
