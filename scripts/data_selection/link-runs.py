import json
import logging
from argparse import ArgumentParser
from itertools import chain
from pathlib import Path

from tqdm import tqdm

from scriptutils.link_utils import link
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main() -> None:  # noqa: PLR-915
    parser = ArgumentParser()
    parser.add_argument("--runs", required=True)
    parser.add_argument("--dl1-link-dir", required=True)
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(logfile=args.log_file, verbose=args.verbose)

    with open(args.runs, "r") as f:
        runs = json.load(f)
    n_runs = len(set(chain(*runs.values())))
    log.info(f"Checking {n_runs} runs")

    outdir_dl1 = Path(args.dl1_link_dir)

    filename_dl1 = "dl1_LST-1.Run{run_id}.h5"
    template_target_dl1 = (
        Path("/fefs/aswg/data/real/DL1/{night}/v0.9/tailcut84") / filename_dl1
    ).as_posix()
    template_linkname_dl1 = (outdir_dl1 / filename_dl1).as_posix()

    progress = tqdm(total=n_runs)
    for night, run_ids in runs.items():
        for run_id in run_ids:
            log.info(f"Linking for run {run_id}")
            target_dl1 = Path(template_target_dl1.format(night=night, run_id=run_id))
            linkname_dl1 = Path(
                template_linkname_dl1.format(
                    night=night,
                    run_id=run_id,
                ),
            )

            link(target_dl1, linkname_dl1)
            progress.update()

    Path(args.output_path).touch()


if __name__ == "__main__":
    main()
