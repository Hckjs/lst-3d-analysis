import logging
from argparse import ArgumentParser

import numpy as np
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.makers import (
    DatasetsMaker,
)
from gammapy.modeling.models import PiecewiseNormSpectralModel

from scriptutils.io import save_datasets_with_models

log = logging.getLogger("__name__")


def main(config, output_datasets, output_models):
    # Standard high-level interface stuff
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)
    analysis.get_observations()

    datasets_settings = analysis.config.datasets
    offset_max = datasets_settings.geom.selection.offset_max

    log.info("Creating reference dataset and makers.")
    stacked = analysis._create_reference_dataset(name="stacked")

    maker = analysis._create_dataset_maker()
    maker_safe_mask = analysis._create_safe_mask_maker()
    bkg_maker = analysis._create_background_maker()
    energy_axis = analysis._make_energy_axis(analysis.config.datasets.geom.axes.energy)
    bkg_spectral = PiecewiseNormSpectralModel(
        energy_axis.center,
        norms=np.ones(len(energy_axis.center)),
    )
    bkg_maker.default_spectral_model = bkg_spectral

    makers = [maker, maker_safe_mask, bkg_maker]
    makers = [maker for maker in makers if maker is not None]

    log.info("Start the data reduction loop.")

    datasets_maker = DatasetsMaker(
        makers,
        stack_datasets=datasets_settings.stack,
        n_jobs=analysis.config.general.n_jobs,
        cutout_mode="trim",
        cutout_width=2 * offset_max,
    )
    analysis.datasets = datasets_maker.run(stacked, analysis.observations)

    save_datasets_with_models(
        analysis.datasets,
        output_datasets,
        output_models,
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-o", "--output-datasets", required=True)
    parser.add_argument("-m", "--output-models", required=True)
    parser.add_argument("-j", "--n-jobs", default=1, type=int)
    args = parser.parse_args()
    main(**vars(args))
