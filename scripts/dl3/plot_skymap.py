import logging
from argparse import ArgumentParser

from astropy import units as u
from gammapy.maps import WcsMap
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)
angle = u.deg


def main(input_path, output_path):
    skymap = WcsMap.read(input_path, map_type="wcs")

    geom = skymap.geom
    source = geom.center_skydir
    pointing_ra = skymap.meta["pointing_ra_deg"]
    pointing_dec = skymap.meta["pointing_dec_deg"]

    edges = u.Quantity(
        [
            u.Quantity(geom.pix_to_coord([i, j]))
            for i, j in zip(range(int(geom.npix[0])), range(int(geom.npix[1])))
        ],
    ).T

    fig, ax = plt.subplots()

    mesh = ax.pcolormesh(
        *edges.to_value(angle),
        skymap.data[0, ...],
    )

    ax.scatter(
        source.ra.to_value(u.deg),
        source.dec.to_value(u.deg),
        ec="k",
        fc="w",
        label="Source",
    )

    ax.scatter(
        pointing_ra,
        pointing_dec,
        ec="k",
        fc="g",
        label="Pointing",
    )

    fig.colorbar(mesh, ax=ax)

    ax.set_xlabel(f"RA / {angle}")
    ax.set_ylabel(f"Dec / {angle}")

    ax.legend()

    ax.set_aspect(1)

    fig.savefig(output_path)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.input_path, args.output_path)
