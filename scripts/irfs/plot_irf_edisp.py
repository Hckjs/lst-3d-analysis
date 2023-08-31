import logging
from argparse import ArgumentParser

from gammapy.irf import load_irf_dict_from_file
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(input_path, output, log_file, verbose):
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    edisp = load_irf_dict_from_file(input_path)["edisp"]

    fig, ax = plt.subplots()
    # TODO Multiple offsets?
    edisp.plot_bias(ax=ax, offset=None, add_cbar=True)
    fig.savefig(output)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(**vars(args))
