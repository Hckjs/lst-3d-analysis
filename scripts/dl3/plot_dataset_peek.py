from argparse import ArgumentParser

import matplotlib
import matplotlib.pyplot as plt
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import Datasets

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

parser = ArgumentParser()
parser.add_argument("-o", "--output", required=True)
parser.add_argument("-c", "--config", required=True)
parser.add_argument("--dataset-path", required=True)
args = parser.parse_args()


def peek(data):
    data.peek()
    fig = plt.gcf()
    return fig


def counts(data):
    fig, ax = plt.subplots()
    data.background.get_spectrum().plot(label="bkg", ax=ax)
    data.counts.get_spectrum().plot(label="counts", ax=ax)
    ax.legend()
    ax.set_title(data.name)
    return fig


def main(config, dataset_path, output):
    config = AnalysisConfig.read(config)

    analysis = Analysis(config)
    # No explicit stacking here to get the plots per dataset
    # Also note, that I do NOT load the models here
    # I am not sure whether I want this with or without models (and thus background fit)
    # right now I want to see how the excess looks like WITHOUT any fitting
    analysis.datasets = Datasets.read(dataset_path)

    figures = []
    for d in analysis.datasets:
        figures.append(peek(d))
        figures.append(counts(d))

    # TODO Mark this as stacked, the name cant be set unfortunately...
    if len(analysis.datasets) > 1:
        stacked = analysis.datasets.stack_reduce()
        figures.append(peek(stacked))
        figures.append(counts(stacked))

    with PdfPages(output) as pdf:
        for fig in figures:
            pdf.savefig(fig)


if __name__ == "__main__":
    main(**vars(args))
