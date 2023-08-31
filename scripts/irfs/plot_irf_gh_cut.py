import logging
from argparse import ArgumentParser

from astropy import units as u
from astropy.table import Table
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(input_path, output, log_file, verbose):
    setup_logging(logfile=log_file, verbose=verbose)

    energy_unit = u.TeV

    gh_cuts = Table.read(input_path, hdu="GH_CUTS")

    fig, ax = plt.subplots()

    ax.errorbar(
        gh_cuts["center"].quantity.to_value(energy_unit),
        gh_cuts["cut"],
        xerr=(gh_cuts["center"] - gh_cuts["low"], gh_cuts["high"] - gh_cuts["center"]),
        ls="",
        label="G/H Cut",
        # lw=2,
        color="black",
    )
    ax.set_xlabel(f"$E_{{\\mathrm{{reco}}}} / {energy_unit}$")
    ax.set_ylabel("Gammaness")

    ax.bar(
        gh_cuts["center"].quantity.to_value(energy_unit),
        gh_cuts["cut"],
        bottom=0,
        width=gh_cuts["high"] - gh_cuts["low"],
        color="gray",
        alpha=0.1,
        label="Discarded Region",
    )

    ax.set_ylim(0, 1.05)

    ax.set_xscale("log")
    ax.legend()

    fig.savefig(output)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(**vars(args))
