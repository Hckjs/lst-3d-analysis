import logging
from argparse import ArgumentParser

import numpy as np
from astropy.io import fits
from astropy.table import Table
from gammapy.analysis import Analysis, AnalysisConfig

from scriptutils.io import load_datasets_with_models
from scriptutils.log import setup_logging

log = logging.getLogger(__name__)


def main(config, datasets_path, models_path, output):
    config = AnalysisConfig.read(config)
    analysis = Analysis(config)
    analysis.get_observations()

    datasets = load_datasets_with_models(datasets_path, models_path)
    hdulist = [fits.PrimaryHDU()]

    # vs livetime
    info_table = datasets.info_table(cumulative=True)
    table = Table(
        {
            "excess": info_table["excess"],
            "sqrt_ts": info_table["sqrt_ts"],
            "livetime": info_table["livetime"],
        },
    )
    table.meta["CUMUL"] = True
    print(table)
    hdulist.append(fits.table_to_hdu(table))

    # vs run id and vs zenith
    # this requires the datasets to not be stacked
    if len(analysis.observations) == len(datasets):
        obs_ids = []
        excess = []
        sqrt_ts = []
        for obs, ds in zip(analysis.observations, datasets):
            obs_ids.append(obs.obs_id)
            # zenith.append(obs.get_pointing_altaz(obs.tmid).zen.to_value(u.deg))
            excess.append(ds.info_dict()["excess"])
            sqrt_ts.append(ds.info_dict()["sqrt_ts"])
        table = Table(
            {
                "excess": np.array(excess),
                "sqrt_ts": np.array(sqrt_ts),
                #   "zenith": np.array(zenith),
                "obs_id": np.array(obs_ids),
            },
        )
        table.meta["CUMUL"] = False
        hdulist.append(fits.table_to_hdu(table))

    fits.HDUList(hdulist).writeto(output, overwrite=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--datasets-path", required=True)
    parser.add_argument("--models-path", required=True, help="Bkg fit")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    setup_logging(logfile=args.log_file, verbose=args.verbose)
    main(args.config, args.dataset_path, args.output)
