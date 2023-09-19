import logging
from argparse import ArgumentParser

import astropy.units as u
import matplotlib
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


log = logging.getLogger(__name__)


def main(input_path, output_path):
    pointings = Table.read(input_path)
    figures = []

    fig, ax = plt.subplots()
    ax.plot(pointings["run_id"], pointings["zen"].quantity.to_value(u.rad), ls="")
    ax.set_xlabel("Run Id")
    ax.set_ylabel("Zenith / deg")
    figures.append(fig)

    ax.plot(
        pointings["run_id"],
        np.cos(pointings["zen"].quantity.to_value(u.rad)),
        ls="",
    )
    ax.set_xlabel("Run Id")
    ax.set_ylabel("Cos (zenith)")
    figures.append(fig)

    fig, ax = plt.subplots()
    ax.plot(pointings["run_id"], pointings["az"].quantity.to_value(u.rad), ls="")
    ax.set_xlabel("Run Id")
    ax.set_ylabel("Azimuth / deg")
    figures.append(fig)

    if output_path is None:
        plt.show()
    else:
        with PdfPages(output_path) as pdf:
            for fig in figures:
                pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_path", required=True)
    parser.add_argument("-o", "--output_path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.input_path, args.output_path)
