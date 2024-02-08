import logging
import pickle
from argparse import ArgumentParser

from gammapy.estimators import TSMapEstimator
from gammapy.modeling.models import (
    GaussianSpatialModel,
    PowerLawSpectralModel,
    SkyModel,
)

from scriptutils.io import load_datasets_with_models
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(datasets_path, models_path, output):
    datasets = load_datasets_with_models(datasets_path, models_path)
    maps = {}
    log.info(f"Running on {len(datasets)} datasets")

    for d in datasets:
        log.info(
            f"Estimating on dataset {d.name} with models {d.models.names}",
        )
        try:
            e_reco = d.background.geom.axes["energy"].edges
            estimator = TSMapEstimator()
            name = d.name
            d = d.to_masked()
            ts_maps = estimator.run(d)
            maps[name] = ts_maps
            log.info(f"Calculating ts maps in bins of energy: {e_reco}")
            for e_min, e_max in zip(e_reco[:-1], e_reco[1:]):
                estimator = TSMapEstimator(energy_edges = [e_min, e_max])
                ts_maps = estimator.run(d)
                maps[f"{name}_{e_min:.1f}_{e_max:.1f}"] = ts_maps

        except ValueError as e:
            log.error(f"Can not compute significance for dataset {d}.")
            log.error(e)
        except AttributeError as e:
            log.error(e)

    # Stacked
    try:
        estimator = TSMapEstimator()
        stacked = datasets.stack_reduce()
        ts_maps = estimator.run(stacked)
        maps["stacked"] = ts_maps
        for e_min, e_max in zip(e_reco[:-1], e_reco[1:]):
            estimator = TSMapEstimator(energy_edges = [e_min, e_max])
            ts_maps = estimator.run(stacked)
            maps[f"stacked_{e_min:.1f}_{e_max:.1f}"] = ts_maps

    except ValueError as e:
        log.error("Can not compute significance for stacked dataset.")
        log.error(e)
    except AttributeError as e:
        log.error(e)

    with open(output, "wb") as f:
        pickle.dump(maps, f)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--datasets-path", required=True)
    parser.add_argument("--models-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    log.info(f"Estimating significances on {args.datasets_path}")
    main(args.datasets_path, args.models_path, args.output)
