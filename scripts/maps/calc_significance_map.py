from argparse import ArgumentParser

import astropy.units as u
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import Datasets
from gammapy.estimators import ExcessMapEstimator

parser = ArgumentParser()
parser.add_argument("-c", "--config", required=True)
parser.add_argument("--dataset-path", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()


def main(config, dataset_path, output):
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)
    analysis.get_observations()

    datasets = Datasets.read(dataset_path)

    # TODO make the smoothing configurable
    estimator = ExcessMapEstimator(0.02 * u.deg, selection_optional=[])
    lima_maps = estimator.run(datasets.stack_reduce())

    lima_maps.write(output, overwrite=True)


if __name__ == "__main__":
    main(**vars(args))
