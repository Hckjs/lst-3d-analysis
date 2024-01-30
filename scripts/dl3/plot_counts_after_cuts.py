import logging
from argparse import ArgumentParser

import matplotlib
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


log = logging.getLogger(__name__)


def plot_cuts(input_path):
    cuts = Table.read(input_path, path="cuts")
    time = cuts.meta["effective_time"]

    x = cuts["center"]
    xerr = (cuts["high"] - x, x - cuts["low"])

    fig, ax = plt.subplots()
    ax.errorbar(
        x,
        cuts["after_trigger"] / time,
        xerr=xerr,
        y_err=np.sqrt(cuts["after_trigger"]) / time,
        ls="",
        label="Trigger",
    )
    ax.errorbar(
        x,
        cuts["after_gh"] / time,
        y_err=np.sqrt(cuts["after_gh"]) / time,
        xerr=xerr,
        ls="",
        label="GH Cut",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(rf"E_{{\text{{reco}}}} \:/\: {x.unit}")
    ax.set_ylabel(cuts["after_trigger"].unit)
    ax.legend(title="Counts after")
    return fig


def plot_cog(input_path):
    cog = Table.read(input_path, path="cog")
    time = cog.meta["effective_time"]

    fig_trigger, ax_trigger = plt.subplots()
    im = ax_trigger.pcolormesh(cog["after_trigger"] / time)
    ax_trigger.set_title("After trigger")
    fig_trigger.colorbar(im)

    fig_gh, ax_gh = plt.subplots()
    im_gh = ax_gh.pcolormesh(cog["after_gh"] / time)
    ax_gh.set_title("After gh")
    fig_gh.colorbar(im_gh)

    for ax in (ax_trigger, ax_gh):
        ax.set_xlabel("x")
        ax.set_xlabel("y")

    return fig_trigger, fig_gh


def plot_intensity(input_path):
    intensity = Table.read(input_path, path="intensity")
    time = intensity.meta["effective_time"]

    x = intensity["center"]
    xerr = (intensity["high"] - x, x - intensity["low"])

    fig, ax = plt.subplots()
    ax.errorbar(
        x,
        intensity["after_trigger"] / time,
        xerr=xerr,
        y_err=np.sqrt(intensity["after_trigger"]) / time,
        ls="",
        label="Trigger",
    )
    ax.errorbar(
        x,
        intensity["after_gh"] / time,
        y_err=np.sqrt(intensity["after_gh"]) / time,
        xerr=xerr,
        ls="",
        label="GH Cut",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(rf"E_{{\text{{reco}}}} \:/\: {x.unit}")
    ax.set_ylabel(intensity["after_trigger"].unit)
    ax.legend(title="Counts after")
    return fig


def main(input_path, output_path):
    fig_cuts = plot_cuts(input_path)
    fig_intensity = plot_intensity(input_path)
    fig_cog_trigger, fig_cog_gh = plot_cog(input_path)
    figures = [fig_cuts, fig_intensity, fig_cog_trigger, fig_cog_gh]

    if output_path is None:
        plt.show()
    else:
        with PdfPages(output_path) as pdf:
            for fig in figures:
                pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.input_path, args.output_path)
