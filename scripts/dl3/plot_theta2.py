import logging
import warnings
from argparse import ArgumentParser

import matplotlib
import numpy as np
from astropy.io import fits
from astropy.table import Table
from gammapy.maps import MapAxis
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)
# this warning happens for empty energy bins
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message="Attempting to set identical low and high ylims makes "
    "transformation singular; automatically expanding.",
)


def plot_theta_squared_table(table, *, preliminary=False, ylim=None):
    theta2_axis = MapAxis.from_bounds(
        min(table["theta2_min"]),
        max(table["theta2_max"]),
        nbin=len(table["theta2_max"]),
        unit=table["theta2_max"].unit,
    )

    x = theta2_axis.center
    xerr = theta2_axis.bin_width

    on = table["counts"]
    off = table["counts_off"]

    fig, ax = plt.subplots()

    ax.errorbar(
        x,
        on,
        xerr=xerr / 2,
        yerr=np.sqrt(on),
        ls="",
        label="On",
        zorder=99,
    )
    ax.bar(
        x,
        off,
        width=xerr,
        color="gray",
        alpha=0.2,
    )
    ax.errorbar(
        x,
        off,
        xerr=xerr / 2,
        yerr=np.sqrt(off),
        ls="",
        label="Off",
    )

    cut = table.meta.get("CUT")
    if cut:
        ax.axvline(
            cut,
            linestyle="dashed",
            color="k",
            label="Rad Max Cut",
        )

    if ylim is None:
        ymin = min(min(on), min(off))
        ymax = max(max(on), max(off))
        ylim = (0.9 * ymin, 1.1 * ymax)
    ax.set_ylim(ylim)
    ax.set_xlim(min(table["theta2_min"]), max(table["theta2_max"]))

    ax.legend(loc="upper right")
    ax.set_xlabel(r"$\theta^2 \:\:/\:\: \mathrm{deg}^2$")
    ax.set_ylabel("Counts")

    if preliminary:
        ax.annotate(
            "PRELIMINARY",
            (1, 0.68),
            xycoords="axes fraction",
            horizontalalignment="right",
            fontsize="larger",
            color="gray",
        )

    return fig, ax


def plot_excess(table, ax):
    theta2_axis = MapAxis.from_bounds(
        min(table["theta2_min"]),
        max(table["theta2_max"]),
        nbin=len(table["theta2_max"]),
        unit=table["theta2_max"].unit,
    )
    on = table["counts"]
    off = table["counts_off"]
    x = theta2_axis.center
    xerr = theta2_axis.bin_width
    excess = on - off
    ax.bar(x = x, width=xerr, height=excess)
    return ax


def main(input_path, output, preliminary):
    figures = []
    with fits.open(input_path) as f:
        plt.rcParams["figure.max_open_warning"] = len(f)
        # Skip primary
        for hdu in f[1:]:
            table = Table.read(hdu)
            fig, ax = plot_theta_squared_table(table, preliminary=preliminary)
            low = table.meta["ELOW"]
            high = table.meta["EHI"]
            ax.set_title(f"{low} - {high}")
            figures.append(fig)

            fig, ax = plt.subplots()
            plot_excess(table, ax)
            ax.set_title(f"{low} - {high}")
            figures.append(fig)

    with PdfPages(output) as pdf:
        for fig in figures:
            pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--preliminary", action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.input_path, args.output, args.preliminary)
