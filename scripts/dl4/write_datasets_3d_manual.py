import logging
from argparse import ArgumentParser

import numpy as np
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.makers import (
    DatasetsMaker,
)
from gammapy.modeling.models import PiecewiseNormSpectralModel

from scriptutils.io import save_datasets_with_models
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)

class DatasetsMakerWithNames(DatasetsMaker):
    tag = "DatasetsMakerWithNames"

    def make_dataset(self, dataset, observation):
        """Make single dataset.

        Parameters
        ----------
        dataset : `~gammapy.datasets.MapDataset`
            Reference dataset
        observation : `Observation`
            Observation
        """
        name=f"Run_{observation.obs_id}"
        if self._apply_cutout:
            cutouts_kwargs = {
                "position": observation.get_pointing_icrs(observation.tmid).galactic,
                "width": self.cutout_width,
                "mode": self.cutout_mode,
                "name": name,
            }
            dataset_obs = dataset.cutout(
                **cutouts_kwargs,
            )
        else:
            dataset_obs = dataset.copy(name=name)

        if dataset.models is not None:
            models = dataset.models.copy()
            models.reassign(dataset.name, dataset_obs.name)
            dataset_obs.models = models

        log.info(f"Computing dataset {dataset_obs.name} for observation {observation.obs_id}")


        for maker in self.makers:
            log.info(f"Running {maker.tag} on dataset {dataset_obs.name}")
            dataset_obs = maker.run(dataset=dataset_obs, observation=observation)

        return dataset_obs


def main(config, output_datasets, output_models, n_jobs):
    # Standard high-level interface stuff
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)
    analysis.get_observations()

    datasets_settings = analysis.config.datasets
    offset_max = datasets_settings.geom.selection.offset_max

    for obs in analysis.observations:
        pointing = obs.pointing
        log.info(type(pointing))
        log.info(f"Location: {pointing.location}")
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
    # make sure norms are all positive
    for p in bkg_spectral.parameters:
        if p.name.startswith("norm_"):
            p.min = 1e-3
            p.max = 1e3

    bkg_maker.default_spectral_model = bkg_spectral

    makers = [maker, maker_safe_mask, bkg_maker]
    makers = [maker for maker in makers if maker is not None]

    log.info("Start the data reduction loop.")

    datasets_maker = DatasetsMakerWithNames(
        makers,
        stack_datasets=datasets_settings.stack,
        n_jobs=n_jobs,
        cutout_mode="trim",
        cutout_width=2 * offset_max,
    )
    analysis.datasets = datasets_maker.run(stacked, analysis.observations)
    log.info("Datasets created")

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
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--n-jobs", default=1, type=int)
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.config, args.output_datasets, args.output_models, args.n_jobs)
