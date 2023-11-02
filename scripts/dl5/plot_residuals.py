import logging
from argparse import ArgumentParser

import matplotlib
import matplotlib.pyplot as plt
from gammapy.analysis import Analysis, AnalysisConfig

from scriptutils.io import load_datasets_with_models
from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)


def main(  # noqa: PLR0913
    config,
    datasets_path,
    bkg_models_path,
    best_model_path,
    output,
    t_start,
    t_stop,
    **kwargs,
):
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)

    analysis.datasets = load_datasets_with_models(datasets_path, bkg_models_path)
    # should expand
    analysis.read_models(best_model_path)

    figures = []
    fig, ax = plt.subplots()
    analysis.datasets.plot_residuals_spatial(ax=ax)

    fig, ax = plt.subplots()
    analysis.datasets.plot_residuals_spectral(ax=ax)

    with PdfPages(output) as pdf:
        for fig in figures:
            pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--datasets-path", required=True)
    parser.add_argument("--bkg-models-path", required=True)
    parser.add_argument("--best-model-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--t-start", required=False)
    parser.add_argument("--t-stop", required=False)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(**vars(args))
