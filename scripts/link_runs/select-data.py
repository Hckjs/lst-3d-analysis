import logging
from argparse import ArgumentParser

import numpy as np
import pandas as pd
from astropy.time import Time

from ..config import Config
from ..log import setup_logging

log = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-c", "--config", required=True)
    args = parser.parse_args()

    config = Config.parse_file(args.config)
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    log.info(config)

    tables = pd.read_html(args.input_file)

    df = tables[0]  # only one table
    log.info(f"Checking a total of {len(df)} runs")

    # sanitize input
    df["Date directory"] = np.char.replace(
        np.array(df["Date directory"], dtype=str),
        "-",
        "",
    ).astype(int)
    df["Run start [UTC]"] = np.array(df["Run start [UTC]"], dtype=np.datetime64)
    df["Run ID"] = df["Run ID"].apply("{:05d}".format)

    # select runs
    mask_source = df["Source name"] == config.source
    log.info(
        f"Found {np.count_nonzero(mask_source)} runs observing {config.source}",
    )

    run_start = Time(df["Run start [UTC]"])
    mask_time = (config.time_start < run_start) & (run_start < config.time_stop)
    log.info(f"Found {np.count_nonzero(mask_time)} runs in the selected timeframe")

    mask = mask_source & mask_time

    df[mask].to_csv(args.output_file, index=False)
