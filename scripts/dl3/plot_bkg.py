import logging
from argparse import ArgumentParser

import astropy.units as u
import matplotlib
from gammapy.irf import Background3D
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


log = logging.getLogger(__name__)


def main(input_path, output):
    bkg_3d = Background3D.read(input_path, hdu="BACKGROUND")
    # TODO: Could use the bins in the file, but unsure about edges vs centers etc...
    energies = u.Quantity([20, 50, 100, 200, 500, 1000, 2000, 5000], u.GeV)
    figures = []

    # there is no ax keyword or similar, figure always gets created in the function :/
    # All in one figure
    bkg_3d.plot_at_energy(energies)
    fig = plt.gcf()
    figures.append(fig)

    # separate figures
    for e in energies:
        bkg_3d.plot_at_energy([e])
        fig = plt.gcf()
        figures.append(fig)

    if output is None:
        plt.show()
    else:
        with PdfPages(output) as pdf:
            for fig in figures:
                pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(**vars(args))
