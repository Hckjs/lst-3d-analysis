import logging
from argparse import ArgumentParser

import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import QTable
from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)


def main(input, output):
    data = QTable.read(input)
    log.info(data.colnames)
    log.info(data)
    figures = []

    fig, ax = plt.subplots()
    ax.plot([1,2,3])
    figures.append(fig)

    with PdfPages(output) as pdf:
        for fig in figures:
            pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.input, args.output)
