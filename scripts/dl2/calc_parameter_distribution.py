import logging
from argparse import ArgumentParser

import numpy as np
from astropy.table import Table
from lstchain.io.io import read_dl2_params

log = logging.getLogger(__name__)

intensity_bins = [80, 200, 800, 3200]


# TODO add data? evaluate irf cuts? cant be the stacked mc then...
def main(input_gamma, input_proton, parameter, output):
    hists = {"gamma": {}, "proton": {}}
    for k, p in zip(("gamma", "proton"), (input_gamma, input_proton)):
        events = read_dl2_params(p)
        log.info(events.columns)
        columns = ["gh_score", "reco_energy", "intensity"]
        events = events[columns]
        for lower, upper in zip(intensity_bins[:-1], intensity_bins[1:]):
            mask = (events["intensity"] >= lower) & (events["intensity"] < upper)
            selected = events[mask]
            # TODO Log or not. save unit
            counts, edges = np.histogram(selected[parameter], bins=100)
            hists[k][f"{lower}_{upper}"] = {
                "lower": lower,
                "upper": upper,
                "bin_edges": edges,
                "counts": counts,
                "parameter": parameter,
            }

    Table(hists).write(output, overwrite=True, serialize_meta=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input-gamma", required=True)
    parser.add_argument("--input-proton", required=True)
    parser.add_argument("--parameter", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    main(
        args.input_gamma,
        args.input_proton,
        args.parameter,
        args.output,
    )
