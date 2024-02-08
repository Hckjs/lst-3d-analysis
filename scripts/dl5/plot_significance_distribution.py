import logging
import pickle
from argparse import ArgumentParser

import matplotlib
import numpy as np
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
    with open(lima_maps_input, "rb") as f:
        maps = pickle.load(f)

    # this might have multiple energy bins and then the shapes will not match
    # There is no reason why the map would be energy dependent.
    # We could also take a slice instead
    exclusion_mask = WcsNDMap.read(exclusion_mask).sum_over_axes()
    excess_counts = []
    npred_counts = []

    for name, lima_maps in maps.items():
        significance_map = lima_maps["sqrt_ts"]
        excess_map = lima_maps["npred_excess"]
        npred_map = lima_maps["npred"]

        fig, (ax1, ax2) = plt.subplots(
            subplot_kw={"projection": lima_maps.geom.wcs},
            ncols=2,
        )
        ax1.set_title(f"Significance map ({name})")
        significance_map.plot(ax=ax1, add_cbar=True)
        exclusion_mask.interp_to_geom(geom=significance_map.geom).reduce_over_axes(
            func=np.logical_or,
        ).plot_mask(ax=ax1, hatches=["\\"], colors="C8")

        ax2.set_title("Excess map")
        excess_map.plot(ax=ax2, add_cbar=True)
        exclusion_mask.interp_to_geom(geom=excess_map.geom).reduce_over_axes(
            func=np.logical_or,
        ).plot_mask(ax=ax2, hatches=["\\"], colors="C8")
        figures.append(fig)

        significance_all = significance_map.data[np.isfinite(significance_map.data)]
        log.info(f"{len(significance_all)} bins in total")
        log.info(
            f"All: {np.mean(significance_all)} +- {np.std(significance_all)}",
        )

        significance_off = significance_map.data[
            np.isfinite(significance_map.data) & exclusion_mask.data.astype(bool)
        ]
        log.info(f"{len(significance_off)} bins in off")
        log.info(
            f"Off: {np.mean(significance_off)} +- {np.std(significance_off)}",
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
        ax.set_title(f"Significance ({name})")
        figures.append(fig)

        if (not "GeV" in name) and (not "TeV" in name):
            excess_counts.append(excess_map.data[
                np.isfinite(significance_map.data) & exclusion_mask.data.astype(bool)
            ].sum())
            npred_counts.append(npred_map.data[
                np.isfinite(significance_map.data) & exclusion_mask.data.astype(bool)
            ].sum())
        

    fig, ax = plt.subplots()
    ax.plot(np.cumsum(excess_counts))
    ax.set_title("Cumulative excess counts in on region")
    ax.set_ylabel("Counts")
    ax.set_xlabel("Observations")
    figures.append(fig)

    fig, ax = plt.subplots()
    ax.plot(np.cumsum(npred_counts))
    ax.set_title("Cumulative npred counts in on region")
    ax.set_ylabel("Counts")
    ax.set_xlabel("Observations")
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
