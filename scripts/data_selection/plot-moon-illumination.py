import logging
from argparse import ArgumentParser

from astroplan.moon import moon_illumination
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, get_moon
from astropy.table import Table
from astropy.time import Time
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
    time = Time(runsummary["time"], format="unix", scale="utc")

    norm = colors.Normalize(0, 1)
    location = EarthLocation.from_geodetic(
        u.Quantity(-17.89139, u.deg),
        u.Quantity(28.76139, u.deg),
        height=u.Quantity(2184, u.m),
    )

    ped_std = runsummary["ped_charge_stddev"]

    altaz = AltAz(obstime=time, location=location)

    moon = get_moon(time, location=location).transform_to(altaz)

    mask_altitude = moon.alt.to_value(u.deg) < 0
    moon_light = moon_illumination(time)
    moon_light[mask_altitude] = 0

    fig, ax = plt.subplots()

    im = ax.scatter(
        moon.alt.to_value(u.deg),
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
