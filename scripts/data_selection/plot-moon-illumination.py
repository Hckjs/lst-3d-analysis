import logging
from argparse import ArgumentParser

from astropy.table import Table
from matplotlib import colors
from matplotlib import pyplot as plt

from scriptutils.config import Config
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main():
    parser = ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("-o", "--output_path", required=True)
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(logfile=args.log_file, verbose=args.verbose)
    config = Config.parse_file(args.config)

    runsummary = Table.read(args.input_path)
    runsummary = runsummary[runsummary["mask_run_selection"]]

    norm = colors.Normalize(0, 1)

    ped_std = runsummary["ped_charge_stddev"]

    moon_alt = runsummary["moon_altitude"]
    mask_altitude = moon_alt < 0
    moon_light = runsummary["moon_illumination"]
    moon_light[mask_altitude] = 0

    fig, ax = plt.subplots()

    im = ax.scatter(
        moon_alt,
        ped_std,
        label="Runs",
        c=moon_light,
        norm=norm,
        edgecolor="k",
    )

    ax.set_ylim(0)
    ax.set_xlim(-90, 80)

    ax.fill_between(
        (0, 1),
        [config.pedestal.ll],
        [config.pedestal.ul],
        alpha=0.1,
        label="Selection",
        transform=ax.get_yaxis_transform(),
    )

    ax.set_ylabel("Pedestal Charge Std. Dev. / p.e.")
    ax.set_xlabel("Altitude / deg")

    fig.colorbar(im, ax=ax, label="Moon Illumination")

    ax.legend()

    fig.savefig(args.output_path)


if __name__ == "__main__":
    main()
