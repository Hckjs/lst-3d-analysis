import logging
from argparse import ArgumentParser

from astropy.table import Table
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(input_path, output_path):
    cuts = Table.read(input_path, path="cuts")

    x = cuts["center"]
    xerr = (cuts["high"] - x, x - cuts["low"])

    fig, ax = plt.subplots()

    ax.errorbar(
        x,
        cuts["after_trigger"],
        xerr=xerr,
        ls="",
        label="Trigger",
    )
    ax.errorbar(
        x,
        cuts["after_gh"],
        xerr=xerr,
        ls="",
        label="GH Cut",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(rf"E_{{\text{{reco}}}} \:/\: {x.unit}")
    ax.set_ylabel(cuts["after_trigger"].unit)

    ax.legend(title="Counts after")

    fig.savefig(output_path)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.input_path, args.output_path)
