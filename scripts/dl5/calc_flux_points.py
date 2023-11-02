import logging
from argparse import ArgumentParser

from astropy.time import Time
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import Datasets

from scriptutils.io import load_datasets_with_models
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def select_timeframe(datasets, t_start, t_stop):
    t_start = Time(t_start, format="mjd") if t_start else datasets.gti.time_start[0]
    t_stop = Time(t_stop, format="mjd") if t_stop else datasets.gti.time_stop[-1]
    return datasets.select_time(t_start, t_stop)


def main(  # noqa: PLR0913
    config,
    datasets_path,
    best_model_path,
    bkg_models_path,
    output,
    t_start,
    t_stop,
    **kwargs,
):
    config = AnalysisConfig.read(config)

    analysis = Analysis(config)

    datasets = load_datasets_with_models(datasets_path, bkg_models_path)
    # Since we read the datasets here, they need not to be stacked
    # Especially our helper script to handle energy-dependent theta-cuts
    # does not stack at all.
    # This is helpful to still have per-run information, but means
    # we need to explicitly stack here, which should be cheap

    # Also if they are not stacked, we can select the ones we want here:
    datasets = select_timeframe(datasets, t_start, t_stop)

    if config.datasets.stack:
        analysis.datasets = Datasets([datasets.stack_reduce()])
    else:
        analysis.datasets = datasets

    log.info(f"Models before reading best fit: {analysis.datasets.models.names}")
    # this should expand
    analysis.read_models(best_model_path)
    log.info(f"Models after reading best fit: {analysis.datasets.models.names}")

    analysis.get_flux_points()

    analysis.flux_points.write(output, overwrite=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--datasets-path", required=True)
    parser.add_argument("--best-model-path", required=True)
    parser.add_argument("--bkg-models-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--t-start", required=False)
    parser.add_argument("--t-stop", required=False)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(**vars(args))
