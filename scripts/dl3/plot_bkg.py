import logging
from argparse import ArgumentParser

import astropy.units as u
import numpy as np
import matplotlib
from gammapy.irf import Background3D
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


log = logging.getLogger(__name__)


def plot_at_energy(
    bkg3d, energy=None, add_cbar=True, ncols=3, figsize=None, **kwargs
):
    """Plot the background rate in Field of view coordinates at a given energy.

    Parameters
    ----------
    energy : `~astropy.units.Quantity`
        list of Energy
    add_cbar : bool
        Add color bar?
    ncols : int
        Number of columns to plot
    figsize : tuple
        Figure size
    **kwargs : dict
        Keyword arguments passed to `~matplotlib.pyplot.pcolormesh`.
    """
    n = len(energy)
    cols = min(ncols, n)
    rows = 1 + (n - 1) // cols
    width = 12
    cfraction = 0.0
    if add_cbar:
        cfraction = 0.15
    if figsize is None:
        figsize = (width, rows * width // (cols * (1 + cfraction)))

    fig, axes = plt.subplots(
        ncols=cols,
        nrows=rows,
        figsize=figsize,
        gridspec_kw={"hspace": 0.2, "wspace": 0.3},
    )

    x = bkg3d.axes["fov_lat"].edges.value
    y = bkg3d.axes["fov_lon"].edges.value
    X, Y = np.meshgrid(x, y)

    for i, ee in enumerate(energy):
        if len(energy) == 1:
            ax = axes
        else:
            ax = axes.flat[i]
        bkg = bkg3d.evaluate(energy=ee)
        bkg_unit = bkg.unit
        bkg = bkg.value
       # with quantity_support():
            #caxes = ax.pcolormesh(X, Y, bkg.squeeze(), **kwargs)
        caxes = ax.pcolormesh(X, Y, bkg.squeeze(), **kwargs)

        print(bkg3d.axes)
        print(bkg3d.axes["fov_lat"])
#        bkg3d.axes["fov_lat"].format_plot_xaxis(ax)
#        bkg3d.axes["fov_lon"].format_plot_yaxis(ax)
        ax.set_title(str(ee))
        if add_cbar:
            label = f"Background [{bkg_unit.to_string(UNIT_STRING_FORMAT)}]"
            cbar = ax.figure.colorbar(caxes, ax=ax, label=label, fraction=cfraction)
            cbar.formatter.set_powerlimits((0, 0))

        row, col = np.unravel_index(i, shape=(rows, cols))
        if col > 0:
            ax.set_ylabel("")
        if row < rows - 1:
            ax.set_xlabel("")
        ax.set_aspect("equal", "box")
    return fig


def main(input_path, output):
    bkg_3d = Background3D.read(input_path, hdu="BACKGROUND")
    # TODO: Could use the bins in the file, but unsure about edges vs centers etc...
    energies = u.Quantity([20, 50, 100, 200, 500, 1000, 2000, 5000], u.GeV)
    figures = []

    # there is no ax keyword or similar, figure always gets created in the function :/
    # All in one figure
    #bkg_3d.plot_at_energy(energies, add_cbar=False)
    #fig = plt.gcf()
    fig = plot_at_energy(bkg_3d, energies, add_cbar=False)
    figures.append(fig)

    # separate figures
    for e in energies:
        fig = plot_at_energy(bkg_3d, [e], add_cbar=False)
        #bkg_3d.plot_at_energy([e], add_cbar=False)
        #fig = plt.gcf()
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

    main(args.input_path, args.output)
