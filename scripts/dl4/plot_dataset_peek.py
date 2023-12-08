import logging
from argparse import ArgumentParser

import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.maps import RegionNDMap, WcsNDMap
from gammapy.modeling.models import (
    PiecewiseNormSpectralModel,
    PowerLawNormSpectralModel,
)

from scriptutils.io import load_datasets_with_models
from scriptutils.log import setup_logging

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)


def peek(data):
    data.peek()
    fig = plt.gcf()
    return fig


def bkg_exclusion(data, exclusion_mask):
    fig, ax = plt.subplots()
    data.counts.sum_over_axes().plot(ax=ax)
    data.mask_safe_image.plot_mask(ax=ax, hatches=["///"], colors="C7")
    exclusion_mask.interp_to_geom(geom=data.counts.geom).reduce_over_axes(
        func=np.logical_or,
    ).plot_mask(ax=ax, hatches=["\\"], colors="C8")
    return fig


def counts(data):
    fig, ax = plt.subplots()

    bkg_spectrum = data.background.get_spectrum()
    bkg_spectrum.plot(label="bkg", ax=ax)

    data.counts.get_spectrum().plot(label="counts", ax=ax)

    # stacked datasets have no bkg model anymore
    if data.models:
        geom = bkg_spectrum.geom.copy()
        energy_axis = geom.axes["energy"].center
        bkg_norm_model = data.models[f"{data.name}-bkg"].spectral_model
        log.info(f"Model is an {bkg_norm_model.tag}")
        # the api for evaluate is slightly different...
        if isinstance(bkg_norm_model, PiecewiseNormSpectralModel):
            log.info("Plotting piecewise norm")
            norms = bkg_norm_model.evaluate(energy_axis)
            norms_map = RegionNDMap(geom, norms.to_value(u.one))
            (bkg_spectrum * norms_map).plot(label="bkg fit", ax=ax)
            ax_norms = ax.twinx()
            e = bkg_norm_model.energy
            n = bkg_norm_model.norms
            ax_norms.plot(e, n, color="k", ls="--", marker=".")

        elif isinstance(bkg_norm_model, PowerLawNormSpectralModel):
            log.info("Plotting global norm")
            norm = bkg_norm_model.parameters["norm"].value
            (bkg_spectrum * norm).plot(label="bkg fit", ax=ax)
    ax.legend()
    ax.set_title(data.name)
    return fig


# maybe move to utils at some point
def list_models(dataset):
    if dataset.models:
        log.info(
            f"Dataset {dataset.name} has the following models: {dataset.models.names}",
        )
    else:
        log.info(f"Dataset {dataset.name} has no models")


def main(config, datasets_path, models_path, output):
    config = AnalysisConfig.read(config)

    Analysis(config)
    datasets = load_datasets_with_models(datasets_path, models_path)

    exclusion_mask = WcsNDMap.read(config.datasets.background.exclusion)
    exclusion_mask.data = exclusion_mask.data.astype(bool)

    figures = []
    for d in datasets:
        list_models(d)
        figures.append(peek(d))
        figures.append(bkg_exclusion(d, exclusion_mask))
        figures.append(counts(d))

    # TODO Mark this as stacked, the name cant be set unfortunately...
    if len(datasets) > 1:
        stacked = datasets.stack_reduce()
        list_models(stacked)
        figures.append(peek(stacked))
        figures.append(counts(stacked))
        figures.append(bkg_exclusion(stacked, exclusion_mask))

    with PdfPages(output) as pdf:
        for fig in figures:
            pdf.savefig(fig)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--datasets-path", required=True)
    parser.add_argument("--models-path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.config, args.datasets_path, args.models_path, args.output)
