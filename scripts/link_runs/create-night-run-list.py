import json
import logging
from argparse import ArgumentParser

import pandas as pd

from ..config import Config
from ..log import setup_logging

log = logging.getLogger(__name__)


def main():
    parser = ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(logfile=args.log_file, verbose=args.verbose)
    Config.parse_file(args.config)

    df = pd.read_csv(args.input_path, dtype={"Run ID": str})

    runs = {}

    for night, group in df.groupby("Date directory"):
        runs[night] = sorted(list(group["Run ID"]))
        log.info(f"Runs of night {night}: {runs}")

    with open(args.output_path, "w") as f:
        json.dump(runs, f)


if __name__ == "__main__":
    main()
