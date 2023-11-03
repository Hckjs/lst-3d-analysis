import logging
from argparse import ArgumentParser

import matplotlib
import numpy as np
from gammapy.estimators import FluxMaps
from gammapy.maps import WcsNDMap
from matplotlib import pyplot as plt
from scipy.stats import norm

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


# TODO Lots of hardcoded values here
def main(lima_maps_input, exclusion_mask, output):
    figures = []
    lima_maps = FluxMaps.read(lima_maps_input)
    # this might have multiple energy bins and then the shapes will not match
    exclusion_mask = WcsNDMap.read(exclusion_mask).sum_over_axes()

    significance_map = lima_maps["sqrt_ts"]
    excess_map = lima_maps["npred_excess"]

    fig, (ax1, ax2) = plt.subplots(
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
    log.info(
        f"All significances: {np.mean(significance_all)} +- {np.std(significance_all)}",
    )
    log.info(
        f"Off significances: {np.mean(significance_off)} +- {np.std(significance_off)}",
    )

    x = np.linspace(-5, 5, 50)
    fig, ax = plt.subplots()
    ax.hist(
        significance_all,
        density=True,
        alpha=0.5,
        color="red",
        label="all bins",
        bins=x,
    )

    ax.hist(
        significance_off,
        density=True,
        alpha=0.5,
        color="blue",
        label="off bins",
        bins=x,
    )

    # Now, fit the off distribution with a Gaussian
    mu, std = norm.fit(significance_off)
    p = norm.pdf(x, mu, std)
    ax.plot(x, p, lw=2, color="black", label=f"mu={mu:.2f}\nstd={std:.2f}")
    ax.legend()
    ax.set_xlabel("Significance")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    # xmin, xmax = np.min(significance_all), np.max(significance_all)
    #    ax.set_xlim(xmin, xmax)
    figures.append(fig)

    if output is None:
        plt.show()
    else:
        with PdfPages(output) as pdf:
            for fig in figures:
                pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input-maps", required=True)
    parser.add_argument("--exclusion-mask", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.input_maps, args.exclusion_mask, args.output)
