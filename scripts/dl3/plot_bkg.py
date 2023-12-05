import logging
from argparse import ArgumentParser

import astropy.units as u
import matplotlib
import numpy as np
from gammapy.irf import Background2D, Background3D
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


log = logging.getLogger(__name__)


def plot3d_at_energy(
    bkg3d,
    energy=None,
    add_cbar=True,
    ncols=3,
    figsize=None,
    **kwargs,
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
    lat, lon = np.meshgrid(x, y)

    for i, ee in enumerate(energy):
        if len(energy) == 1:
            ax = axes
        else:
            ax = axes.flat[i]
        bkg = bkg3d.evaluate(energy=ee)
        bkg_unit = bkg.unit
        bkg = bkg.value
        caxes = ax.pcolormesh(lat, lon, bkg.squeeze(), **kwargs)
        ax.set_title(str(ee))
        if add_cbar:
            label = f"Background [{bkg_unit}]"
            cbar = ax.figure.colorbar(caxes, ax=ax, label=label, fraction=cfraction)
            cbar.formatter.set_powerlimits((0, 0))

        row, col = np.unravel_index(i, shape=(rows, cols))
        if col > 0:
            ax.set_ylabel("")
        if row < rows - 1:
            ax.set_xlabel("")
        ax.set_aspect("equal", "box")
    return fig


def plot_at_energy(bkg, energy=None, add_cbar=True, ncols=3, figsize=None, **kwargs):
    if isinstance(bkg, Background3D):
        return plot3d_at_energy(bkg, energy, add_cbar, ncols, figsize, **kwargs)
    else:
        # Avoid https://github.com/gammapy/gammapy/issues/4950
        # Its only some plots, no need to kill the whole workflow
        try:
            fig = plot3d_at_energy(
                bkg.to_3d(),
                energy,
                add_cbar,
                ncols,
                figsize,
                **kwargs,
            )
        except ValueError as e:
            log.error(e)
            fig = plt.figure()
        return fig


def main(input_path, output):
    try:
        bkg = Background3D.read(input_path, hdu="BACKGROUND")
    except Exception:
        bkg = Background2D.read(input_path, hdu="BACKGROUND")
    # TODO: Could use the bins in the file, but unsure about edges vs centers etc...
    energies = u.Quantity([20, 50, 100, 200, 500, 1000, 2000, 5000], u.GeV)
    figures = []

    # there is no ax keyword or similar, figure always gets created in the function :/
    # All in one figure
    fig = plot_at_energy(bkg, energies, add_cbar=False)
    figures.append(fig)

    if not isinstance(bkg, Background2D):
        bkg = bkg.to_2d()

    fig, ax = plt.subplots()
    bkg.plot_energy_dependence(ax=ax)
    figures.append(fig)

    fig, ax = plt.subplots()
    bkg.plot_offset_dependence(ax=ax)
    figures.append(fig)

    fig, ax = plt.subplots()
    bkg.plot_spectrum(ax=ax)
    figures.append(fig)

    # separate figures
    for e in energies:
        fig = plot_at_energy(bkg, [e], add_cbar=False)
        # bkg_3d.plot_at_energy([e], add_cbar=False)
        # fig = plt.gcf()
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
