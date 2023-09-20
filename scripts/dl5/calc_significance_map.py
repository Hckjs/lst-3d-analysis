import logging
from argparse import ArgumentParser

from gammapy.datasets import Datasets
from gammapy.estimators import TSMapEstimator
from gammapy.modeling.models import (
    GaussianSpatialModel,
    PowerLawSpectralModel,
    SkyModel,
)

from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(dataset_path, output):
    datasets = Datasets.read(dataset_path)
    # TS map estimator only works on a single dataset (?)
    stacked = datasets.stack_reduce()

    # TODO make the smoothing configurable
    # Understand how exactly this works
    # Do I need the best fit models attached? fit bkg? No, right?
    # Kernel width vs spatial model size?
    model = SkyModel(
        spectral_model=PowerLawSpectralModel(),
        spatial_model=GaussianSpatialModel(sigma="0.05 deg"),
    )
    estimator = TSMapEstimator(
        model=model,
        kernel_width="0.2 deg",
        selection_optional="all",
    )
    ts_maps = estimator.run(stacked)
    ts_maps.write(output, overwrite=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--dataset-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.config, args.dataset_path, args.output)
