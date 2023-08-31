import logging
from argparse import ArgumentParser

import astropy.units as u
import matplotlib
from gammapy.irf import load_irf_dict_from_file
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


log = logging.getLogger(__name__)


def main(input_path, output, log_file, verbose):
    setup_logging(logfile=log_file, verbose=verbose)

    psf = load_irf_dict_from_file(input_path)["psf"]
    figures = []

    # TODO Could create my own plots as well,
    # but I think these are fine for the time being
    fig, ax = plt.subplots()
    figures.append(fig)
    psf.plot_containment_radius(
        ax=ax,
        fraction=0.68,
    )

    fig, ax = plt.subplots()
    figures.append(fig)
    psf.plot_containment_radius_vs_energy(
        fraction=(0.68,),
        offset=[0, 0.4, 0.8, 1.2] * u.deg,
        ax=ax,
    )

    fig, ax = plt.subplots()
    figures.append(fig)
    psf.plot_containment_radius_vs_energy(
        fraction=(0.95,),
        offset=[0, 0.4, 0.8, 1.2] * u.deg,
        ax=ax,
    )

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
    main(**vars(args))
