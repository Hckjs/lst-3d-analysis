from argparse import ArgumentParser

import matplotlib
import numpy as np
from gammapy.estimators import FluxMaps
from gammapy.maps import WcsNDMap
from matplotlib import pyplot as plt
from scipy.stats import norm

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

parser = ArgumentParser()
parser.add_argument("--lima-maps-input", required=True)
parser.add_argument("--exclusion-map-input", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()


def main(lima_maps_input, exclusion_map_input, output):
    figures = []
    lima_maps = FluxMaps.read(lima_maps_input)
    exclusion_mask = WcsNDMap.read(exclusion_map_input)

    significance_map = lima_maps["sqrt_ts"]
    excess_map = lima_maps["npred_excess"]

    fig, (ax1, ax2) = plt.subplots(
        figsize=(11, 5),
        subplot_kw={"projection": lima_maps.geom.wcs},
        ncols=2,
    )
    ax1.set_title("Significance map")
    significance_map.plot(ax=ax1, add_cbar=True)
    ax2.set_title("Excess map")
    excess_map.plot(ax=ax2, add_cbar=True)
    figures.append(fig)

    significance_map_off = significance_map * exclusion_mask
    significance_all = significance_map.data[np.isfinite(significance_map.data)]
    significance_off = significance_map_off.data[np.isfinite(significance_map_off.data)]

    fig, ax = plt.subplots()
    ax.hist(
        significance_all,
        density=True,
        alpha=0.5,
        color="red",
        label="all bins",
        bins=21,
    )

    ax.hist(
        significance_off,
        density=True,
        alpha=0.5,
        color="blue",
        label="off bins",
        bins=21,
    )

    # Now, fit the off distribution with a Gaussian
    mu, std = norm.fit(significance_off)
    x = np.linspace(-8, 8, 50)
    p = norm.pdf(x, mu, std)
    ax.plot(x, p, lw=2, color="black", label=f"mu={mu:.2f}\nstd={std:.2f}")
    ax.legend()
    ax.set_xlabel("Significance")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    xmin, xmax = np.min(significance_all), np.max(significance_all)
    ax.set_xlim(xmin, xmax)
    figures.append(fig)

    if output is None:
        plt.show()
    else:
        with PdfPages(output) as pdf:
            for fig in figures:
                pdf.savefig(fig)


if __name__ == "__main__":
    main(**vars(args))
