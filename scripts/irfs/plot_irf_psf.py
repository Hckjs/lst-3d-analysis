import logging
from argparse import ArgumentParser

import numpy as np
from gammapy.irf import load_irf_dict_from_file
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(input_path, output, log_file, verbose):
    setup_logging(logfile=log_file, verbose=verbose)

    psf = load_irf_dict_from_file(input_path)["psf"]

    psf.axes["energy_true"]
    offsets = psf.axes["offset"]

    np.atleast_1d(offsets.center)
    np.atleast_1d(offsets.bin_width)
    # TODO Create my own plots, this is more of a placeholder
    # and it even breaks right now due to tight layout and colorbar ...
    # psf.peek()

    fig = plt.gcf()
    fig.savefig(output)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(**vars(args))
