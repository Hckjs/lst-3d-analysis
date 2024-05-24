import logging
from argparse import ArgumentParser

from astropy.table import Table
import arviz as az
import numpy as np
import astropy.constants as c
import astropy.units as u
import yaml
from gammapy.astro.darkmatter import JFactory
from gammapy.astro.darkmatter.profiles import BurkertProfile, EinastoProfile, NFWProfile
from gammapy.maps import WcsNDMap
from titrate.upperlimits import ULFactory

from scriptutils.io import load_datasets_with_models
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)

profiles = {
    "NFW": NFWProfile,
    "Burkert": BurkertProfile,
    "Einasto": EinastoProfile,
    #    "Zhao": ZhaoProfile,
}
parameters = {
    "NFW": ["rho_s", "r_s"],
    "Burkert": ["rho_s", "r_s"],
    "Einasto": ["rho_s", "r_s", "alpha"],
    #    "Zhao": ZhaoProfile,
}

def create_jmap(profile, distance, geom):
    log.info(profile)
    log.info(distance)
    jfactory = JFactory(
        geom=geom,
        profile=profile,
        distance=distance,
    )
    jfactor = jfactory.compute_differential_jfactor()
    jfact_map = WcsNDMap(geom=geom, data=jfactor.value, unit=jfactor.unit)
    return jfact_map


def add_units(values):
    for k, v in values.items():
        if isinstance(k, u.Quantity):
            continue
        if k == "r_s":
            values[k] = v * u.kpc
        elif k == "rho_s":
            values[k] = (v * u.M_sun / u.kpc**3 * c.c**2).to(u.GeV / u.cm**3)
        else:
            continue
    return values


def main(
    datasets_path,
    models_path,
    config_path,
    output,
    channel,
    **kwargs,
):
    with open(config_path) as f:
        config = yaml.safe_load(f)
    p = config["profile"]
    distance = u.Quantity(p["distance"])
    profile = profiles[p["name"]]
    # TODO Load arviz inference data here
    if v := p["params"].get("values"):
        values = v
    elif path := p["params"].get("table_path"):
        table = Table.read(path)
        # should match already, but just to make sure
        table["rho_s"] = table["rho_s"].quantity.to_value(u.M_sun/u.kpc**3)
        table["r_s"] = table["r_s"].quantity.to_value(u.kpc)
        if "beta0" in table:
            table.remove_column("beta0")
        table = table[parameters[p["name"]]]
        values = {c: np.median(table[c]) for c in table.colnames}
    elif path := p["params"].get("idata_path"):
        raise NotImplementedError()
    else:
        raise KeyError()
    values = add_units(values)

    log.info(profile)
    log.info(values)
    # TODO Right now we lose all model information in titrate even if we load them here
    datasets = load_datasets_with_models(datasets_path, models_path)
    # for now stack the datasets
    stacked = datasets.stack_reduce()
    geom = stacked.geoms["geom"].to_image()
    j_map = create_jmap(profile(**values), distance, geom)

    masses = config["mass_range"]

    # TODO This is only for testing purposes. of course we cannot fake data here...
    # The calculation will fail due to the bad bkg description though and I want
    # to test the pipeline
    stacked.fake()

    ulfactory = ULFactory(
        stacked,
        [channel],
        u.Quantity(masses["min"]),
        u.Quantity(masses["max"]),
        masses["n_values"],
        j_map,
        max_workers=8,
    )
    ulfactory.compute()
    ulfactory.save_results(output, overwrite=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--datasets-path", required=True)
    parser.add_argument("--models-path", required=True)
    parser.add_argument("--config-path", required=True)
    parser.add_argument("--channel", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    main(**vars(args))
