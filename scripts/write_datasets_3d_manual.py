from argparse import ArgumentParser
from os import cpu_count
import numpy as np
import yaml
from astropy.coordinates import SkyCoord
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import SpectrumDataset
from gammapy.makers import (
    DatasetsMaker,
    FoVBackgroundMaker,
    MapDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom
from regions import PointSkyRegion, Regions
from gammapy.modeling.models import PiecewiseNormSpectralModel


def main(config, bkg_config, output, n_off_regions, n_jobs):
    """
    Create dl4 datasets from a dl3 datastore for a pointlike analysis
    This can currently (Gammapy 0.20.1) not be done with the high-level interface
    as energy-dependent theta-cuts have just been added.
    """
    if n_jobs == 0:
        n_jobs = cpu_count()
    # Standard high-level interface stuff
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)
    analysis.get_observations()

    # Define things for the dataset maker step
    # point sky region > circle sky region for energy dependent cuts
    on_region = analysis.config.datasets.on_region
    target_position = SkyCoord(on_region.lon, on_region.lat, frame=on_region.frame)
    on_region = PointSkyRegion(target_position)

    with open(bkg_config) as f:
        b = yaml.safe_load(f)["exclusion_regions"]
        exclusion_regions = [Regions.parse(reg,format='ds9') for reg in b["exclusion_regions"]]

    energy_axis_config = config.datasets.geom.axes.energy
    energy_axis = MapAxis.from_bounds(
        name="energy",
        lo_bnd=energy_axis_config.min.value,
        hi_bnd=energy_axis_config.max.to_value(energy_axis_config.min.unit),
        nbin=energy_axis_config.nbins,
        unit=energy_axis_config.min.unit,
        interp="log",
        node_type="edges",
    )
    energy_axis_true_config = config.datasets.geom.axes.energy_true
    energy_axis_true = MapAxis.from_bounds(
        name="energy_true",
        lo_bnd=energy_axis_true_config.min.value,
        hi_bnd=energy_axis_true_config.max.to_value(energy_axis_true_config.min.unit),
        nbin=energy_axis_true_config.nbins,
        unit=energy_axis_true_config.min.unit,
        interp="log",
        node_type="edges",
    )
    geom = RegionGeom.create(region=on_region, axes=[energy_axis])
    empty = MapDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name="empty"
    )

    dataset_maker = MapDatasetMaker()
    safe_mask_maker = analysis._create_safe_mask_maker()
    bkg_spektral = PiecewiseNormSpectralModel(energy_axis.center, norms=np.ones(len(energy_axis.center)))
    bkg_maker = FoVBackgroundMaker(
            method=config.background.parameters["method"],# Needs to be set
            exclusion_mask=~geom.region_mask(regions=[exclusion_regions]),
            spectral_model=bkg_spektral,
            )

    makers = [
        dataset_maker,
        safe_mask_maker,
        bkg_maker,
    ]

    datasets_maker = DatasetsMaker(makers, stack_datasets=False, n_jobs=n_jobs)
    datasets = datasets_maker.run(
        empty,
        analysis.observations,
    )

    datasets.write(output, overwrite=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-b", "--bkg-config", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--n-off-regions", default=1, type=int)
    parser.add_argument("-j", "--n-jobs", default=1, type=int)
    args = parser.parse_args()
    main(**vars(args))
