import logging
from argparse import ArgumentParser

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
    # TODO: This could be a h5 table with units attached
    pointings = Table.read(input_path)
    figures = []

    fig, ax = plt.subplots()
    ax.scatter(pointings["run_id"], pointings["zen"])
    ax.set_xlabel("Run Id")
    ax.set_ylabel("Zenith / deg")
    figures.append(fig)

    fig, ax = plt.subplots()
    ax.scatter(
        pointings["run_id"],
        np.cos(np.deg2rad(pointings["zen"])),
    )
    ax.set_xlabel("Run Id")
    ax.set_ylabel("Cos (zenith)")
    figures.append(fig)

    fig, ax = plt.subplots()
    ax.scatter(pointings["run_id"], pointings["az"])
    ax.set_xlabel("Run Id")
    ax.set_ylabel("Azimuth / deg")
    figures.append(fig)

    for c in pointings.colnames:
        if c.startswith("cosmic"):
            fig, ax = plt.subplots()
            ax.scatter(pointings["run_id"], pointings[c])
            ax.set_xlabel("Run Id")
            ax.set_ylabel(c)
            figures.append(fig)

            fig, ax = plt.subplots()
            ax.scatter(pointings["zen"], pointings[c])
            ax.set_xlabel("Zenith / deg")
            ax.set_ylabel(c)
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
