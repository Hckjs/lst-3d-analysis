from argparse import ArgumentParser
from astropy.coordinates import SkyCoord
from gammapy.maps import MapAxis, RegionGeom
from gammapy.analysis import Analysis, AnalysisConfig
from regions import CircleSkyRegion

parser = ArgumentParser()
parser.add_argument("-c", "--config", required=True)
parser.add_argument("-o", "--output", required=True)
#parser.add_argument("--exclusion-config", required=True)
args = parser.parse_args()


def create_region(lon, lat, frame, radius):
    on_center = SkyCoord(lon, lat, frame=frame)
    on_region = CircleSkyRegion(on_center, radius)
    return on_region


def main(config, output):
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)
    analysis.get_observations()

    g = analysis._create_geometry()
    e = g.axes["energy"]
    geom_image = g.to_image().to_cube([e.squash()])

    on_region = create_region(
            lon=config.datasets.on_region.lon,
            lat=config.datasets.on_region.lat,
            frame=config.datasets.on_region.frame,
            radius=config.datasets.on_region.radius
            )

    # TODO add more exclusions.  
    # Make the exclusion mask(s)
    # This is from a self written config
#    with open(exclusion_config) as f:
#        c = yaml.safe_load(f)
#    for s, p in c.items():
    # construct like above and add to list

    all_regions = [on_region]
    exclusion_mask = ~geom_image.region_mask(all_regions)
    exclusion_mask.write(output)


if __name__ == "__main__":
    main(**vars(args))
