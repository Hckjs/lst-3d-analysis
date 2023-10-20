import logging
from argparse import ArgumentParser

import numpy as np
from astropy import units as u
from gammapy.estimators import FluxProfileEstimator
from gammapy.modeling.models import PowerLawSpectralModel
from gammapy.utils.regions import (
    make_concentric_annulus_sky_regions,
)

from scriptutils.io import load_datasets_with_models
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(datasets_path, models_path, output):
    datasets = load_datasets_with_models(datasets_path, models_path)
    # Should be the same geom for all of them by construction
    center = datasets[0].counts.geom.center_skydir

    # TODO Hardcoded values.
    regions = make_concentric_annulus_sky_regions(
        center=center,
        radius_max="1.5 deg",
        nbin=15,
    )
    flux_profile_estimator = FluxProfileEstimator(
        regions=regions,
        spectrum=PowerLawSpectralModel(index=2.3),
        energy_edges=[100, 5000] * u.GeV,
        selection_optional=["ul", "scan"],
        norm_values=np.linspace(-1, 5, 11),
    )
    profile = flux_profile_estimator.run(datasets=datasets)
    # TODO this seems to lose the stat profile
    profile.write(
        output,
        format="profile",
        overwrite=True,
        sed_type="dnde",
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--datasets-path", required=True)
    parser.add_argument("--models-path", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(args.dataset_path, args.output)
