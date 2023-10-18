import logging
from argparse import ArgumentParser

import matplotlib
import matplotlib.pyplot as plt
from gammapy.estimators import FluxMaps

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)


def main(flux_maps, output):
    ts_maps = FluxMaps.read(flux_maps)
    figures = []

    fig = plt.figure()
    ax = ts_maps.flux.plot(add_cbar=True, fig=fig)
    ax.set_title("Flux")
    figures.append(fig)

    fig = plt.figure()
    ax = ts_maps.flux_err_maps.plot(add_cbar=True, fig=fig)
    ax.set_title("Flux Error")
    figures.append(fig)

    fig = plt.figure()
    ax = ts_maps.flux_ul.plot(add_cbar=True, fig=fig)
    ax.set_title("Flux UL")
    figures.append(fig)

    fig = plt.figure()
    ax = ts_maps.sqrt_ts.plot(add_cbar=True, fig=fig)
    ax.set_title("Sqrt(TS)")
    figures.append(fig)

    fig = plt.figure()
    ax = ts_maps.npred_excess.plot(add_cbar=True, fig=fig)
    ax.set_title("N Pred Excess")
    figures.append(fig)

    fig = plt.figure()
    ax = ts_maps.npred_excess_ul.plot(add_cbar=True, fig=fig)
    ax.set_title("N Pred Excess UL")
    figures.append(fig)

    if output is None:
        plt.show()
    else:
        with PdfPages(output) as pdf:
            for fig in figures:
                pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--flux-maps", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.flux_maps, args.output)
