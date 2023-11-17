import logging
from argparse import ArgumentParser

import astropy.units as u
import matplotlib
import numpy as np
from astropy.table import Table
from gammapy.data import DataStore
from matplotlib import pyplot as plt

from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


log = logging.getLogger(__name__)


def main(input_path, output):
    ds = DataStore.from_file(input_path)
    figures = []

    zen = []
    az = []
    ontime = []
    elapsed_time = []
    counts = []
    obs_ids = []
    deadc = []

    for o in ds.get_observations():
        i = o.obs_info
        log.info(o.obs_id)
        zen.append(90 - i["ALT_PNT"])
        az.append(i["AZ_PNT"])
        ontime.append(i["ONTIME"])
        elapsed_time.append(i["TELAPSE"])
        counts.append(len(o.events.table))
        deadc.append(i["DEADC"])
        obs_ids.append(o.obs_id)
    rate = np.array(counts) / u.Quantity(ontime, u.s)
    rate2 = np.array(counts) / u.Quantity(elapsed_time, u.s)
    rate3 = np.array(counts) / u.Quantity(elapsed_time, u.s) / np.array(deadc)

    Table(
        {
            "rate": rate,
            "counts": counts,
            "ontime": ontime,
            "elapsed_time": elapsed_time,
            "obs_id": obs_ids,
            "zen": zen,
            "az": az,
        },
    )

    for r in (rate, rate2, rate3):
        fig, ax = plt.subplots()
        ax.scatter(obs_ids, r)
        ax.set_xlabel("Obs id")
        ax.set_ylabel(f"Rate / {rate.unit}")
        figures.append(fig)

        fig, ax = plt.subplots()
        ax.scatter(zen, r)
        ax.set_xlabel("Zenith / deg")
        ax.set_ylabel(f"Rate / {rate.unit}")
        figures.append(fig)

        fig, ax = plt.subplots()
        ax.scatter(np.cos(zen), r)
        ax.set_xlabel("cos(zd)")
        ax.set_ylabel(f"Rate / {rate.unit}")
        figures.append(fig)

        fig, ax = plt.subplots()
        ax.scatter(az, r)
        ax.set_xlabel("Azimuth / deg")
        ax.set_ylabel(f"Rate / {rate.unit}")
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
