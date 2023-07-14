

from argparse import ArgumentParser
from gammapy.analysis import Analysis, AnalysisConfig


parser = ArgumentParser()
parser.add_argument("-c", "--config", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()


def main(config, output):

    # Standard high-level interface stuff
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)
    analysis.get_observations()
    analysis.get_datasets()

    analysis.datasets.write(output, overwrite=True)


if __name__ == "__main__":
    main(**vars(args))
