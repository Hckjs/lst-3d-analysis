import logging
from argparse import ArgumentParser

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_moon
from astropy.table import Table
from astropy.time import Time

from scriptutils.config import Config
from scriptutils.log import setup_logging
from scriptutils.stats import bounds_std as bounds

log = logging.getLogger(__name__)


def get_mask(x, le=np.inf, ge=-np.inf):
    return np.logical_and(
        np.greater_equal(x, ge),
        np.less_equal(x, le),
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("datacheck_path")
    parser.add_argument("--output-runlist", required=True)
    parser.add_argument("--output-datachecks", required=True)
    parser.add_argument("--output-config", required=True)
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    config = Config.parse_file(args.config)
    output_config = config.copy()

    runs = pd.read_csv(args.input_path, dtype={"Run ID": str})
    setup_logging(logfile=args.log_file, verbose=args.verbose)

    run_ids = np.array(runs["Run ID"], dtype=int)

    runsummary = Table.read(args.datacheck_path)
    log.info(f"Checking a total of {len(runsummary)} runs")

    time = Time(runsummary["time"], format="unix", scale="utc")

    mask_time = get_mask(
        time,
        ge=config.time_start,
        le=config.time_stop,
    )

    mask_run_id = np.in1d(np.array(runsummary["runnumber"]), run_ids)
    log.info(
        "Runs observing the source: "
        f"{np.count_nonzero(mask_run_id)} / {len(runsummary)} ",
    )

    mask_pedestals_ok = np.isfinite(runsummary["num_pedestals"])

    # Exclude runs that are too far from source

    tel_pointing = SkyCoord(
        ra=runsummary["mean_ra"],
        dec=runsummary["mean_dec"],
    )

    source_coordinates = SkyCoord(
        ra=config.source_ra_deg,
        dec=config.source_dec_deg,
        unit=u.deg,
    )

    separation = tel_pointing.separation(source_coordinates)

    mask_separation_low = np.isclose(
        u.Quantity(0.4, u.deg),
        separation,
        rtol=0,
        atol=0.15,
    )

    zenith = 90 - np.rad2deg(runsummary["mean_altitude"])
    mask_zenith = zenith < config.max_zenith_deg

    mask = (
        mask_pedestals_ok & mask_run_id & mask_time & mask_separation_low & mask_zenith
    )

    runsummary["mask_run_id"] = mask_run_id
    runsummary["mask_time"] = mask_time
    runsummary["mask_separation_low"] = mask_separation_low
    runsummary["mask_pedestals_ok"] = mask_pedestals_ok
    runsummary["mask_zenith"] = mask_zenith

    runsummary["mask_run_selection"] = mask

    after_pedestals_run_id_time_separation = np.count_nonzero(mask)
    s = (
        "After selecting the time of the dataset, "
        "after removing runs with problems in pedestals and after removing "
        "mispointed runs and those with high zenith,"
        f" {after_pedestals_run_id_time_separation} runs are kept."
    )
    log.info(s)

    mask = mask & (~np.isin(runsummary["runnumber"], config.never_include))
    log.info(
        f"Selected runs after removing from blacklist: {runsummary['runnumber'][mask]}",
    )

    location = EarthLocation.from_geodetic(
        u.Quantity(-17.89139, u.deg),
        u.Quantity(28.76139, u.deg),
        height=u.Quantity(2184, u.m),
    )

    altaz = AltAz(obstime=time[mask], location=location)

    ped_std = runsummary["ped_charge_stddev"]

    ped_ll = config.pedestal.ll
    ped_ul = config.pedestal.ul

    if config.pedestal.sigma is not None:
        sigma = config.pedestal.sigma
        log.info(
            "Calculating pedestal cuts based on configured sigma interval "
            "for runs selected before and with moon below horizon.",
        )

        altaz = AltAz(obstime=time, location=location)
        moon = get_moon(time, location=location).transform_to(altaz)

        ped_ll, ped_ul = bounds(ped_std[mask & (moon.alt.to_value(u.deg) < 0)], sigma)

        log.info("Calculated %f sigma interval is (%f, %f)", sigma, ped_ll, ped_ul)

    mask_pedestal_charge = get_mask(ped_std, ge=ped_ll, le=ped_ul)

    runsummary["mask_pedestal_charge"] = mask_pedestal_charge & mask
    mask = runsummary["mask_pedestal_charge"] & mask

    after_pedestal_charge = np.count_nonzero(mask)
    s = (
        "After checking for the deviation of the pedestal charges, "
        "possibly dependent on the moon elevation and illumination, "
        f"{after_pedestal_charge} runs are kept."
    )
    log.info(s)

    # Check cosmics rates

    cosmics_rate = runsummary["cosmics_rate"]

    cos_ll = config.cosmics.ll
    cos_ul = config.cosmics.ul

    if config.cosmics.sigma is not None:
        sigma = config.cosmics.sigma
        log.info(
            "Calculating cosmics cuts based on "
            "configured sigma interval "
            "for runs selected before.",
        )

        cos_ll, cos_ul = bounds(cosmics_rate[mask], sigma)

        log.info("Calculated %f sigma interval is (%f, %f)", sigma, cos_ll, cos_ul)

    mask_cosmics = get_mask(cosmics_rate, ge=cos_ll, le=cos_ul)

    runsummary["mask_cosmics"] = mask_cosmics

    mask = runsummary["mask_cosmics"] & mask

    after_cosmics = np.count_nonzero(mask)
    s = (
        "After checking the rate of cosmics and selecting "
        "based on AGN Zoo Paper cuts, "
        f"{after_cosmics} runs are kept."
    )
    log.info(s)

    cosmics_rate_above10 = runsummary["cosmics_rate_above10"]
    cosmics_rate_above30 = runsummary["cosmics_rate_above30"]

    cos_10_ll = config.cosmics_10.ll
    cos_10_ul = config.cosmics_10.ul

    if config.cosmics_10.sigma is not None:
        sigma = config.cosmics_10.sigma
        log.info(
            "Calculating cosmics above 10 cuts based "
            "on configured sigma interval "
            "for runs selected before.",
        )

        cos_10_ll, cos_10_ul = bounds(cosmics_rate_above10[mask], sigma)

        log.info(
            "Calculated %f sigma interval is (%f, %f)",
            sigma,
            cos_10_ll,
            cos_10_ul,
        )

    cos_30_ll = config.cosmics_30.ll
    cos_30_ul = config.cosmics_30.ul

    if config.cosmics_30.sigma is not None:
        sigma = config.cosmics_30.sigma
        log.info(
            "Calculating cosmics above 30 cuts based "
            "on configured sigma interval for "
            "runs selected before.",
        )

        cos_30_ll, cos_30_ul = bounds(cosmics_rate_above30[mask], sigma)

        log.info(
            "Calculated %f sigma interval is (%f, %f)",
            sigma,
            cos_30_ll,
            cos_30_ul,
        )

    mask_above10 = get_mask(
        cosmics_rate_above10,
        le=cos_10_ul,
        ge=cos_10_ll,
    )
    mask_above30 = get_mask(
        cosmics_rate_above30,
        le=cos_30_ul,
        ge=cos_30_ll,
    )

    mask_above = mask_above10 & mask_above30
    runsummary["mask_cosmics_above"] = mask_above
    mask = runsummary["mask_cosmics_above"] & mask

    after_cosmics_above_n = np.count_nonzero(mask)
    s = (
        "After cutting on the rate of cosmics with more than 10 p.e. "
        f"(resp. 30 p.e.) {after_cosmics_above_n} runs are kept."
    )
    log.info(s)

    duration = np.sum(runsummary["elapsed_time"][mask].quantity).to(u.h)
    s = (
        f"Selected a total of {np.count_nonzero(mask)} runs "
        f"with observation time of {duration:.2f}"
    )
    log.info(s)

    log.info(
        f"Selected runs after all selection steps: {runsummary['runnumber'][mask]}",
    )

    mask = mask | np.isin(runsummary["runnumber"], config.always_include)
    log.info(
        f"Selected runs after adding from whitelist: {runsummary['runnumber'][mask]}",
    )

    runs[mask].to_csv(args.output_runlist, index=False)

    runsummary.write(
        args.output_datachecks,
        serialize_meta=True,
        overwrite=True,
        path="data",
        compression=True,
    )

    output_config.pedestal = {"ul": ped_ul, "ll": ped_ll, "sigma": None}
    output_config.cosmics = {"ul": cos_ul, "ll": cos_ll, "sigma": None}
    output_config.cosmics_10 = {"ul": cos_10_ul, "ll": cos_10_ll, "sigma": None}
    output_config.cosmics_30 = {"ul": cos_30_ul, "ll": cos_30_ll, "sigma": None}
    with open(args.output_config, "w") as f:
        f.write(output_config.json())
