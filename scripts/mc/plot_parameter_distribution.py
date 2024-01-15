import logging
from argparse import ArgumentParser

import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

log = logging.getLogger(__name__)


def main(input_path, output_path):
    hists = Table.read(input_path)
    for i, (k, v) in enumerate(hists.items()):
        if i == 0:
            fig, axes = plt.subplots(2, np.ceil(len(v) / 2))
        for j, (results, ax) in enumerate(zip(v.values(), axes)):
            e = results["edges"]
            x = (e[:-1] + e[1:]) / 2
            xerr = (e[1:] - x, x - e[:-1])
            ax.errorbar(x, results["counts"], xerr=xerr, ls="")
        ax.set_yscale("log")
        ax.set_xlabel(results["parameter"])
        ax.set_ylabel("Counts")
        ax.set_title(f"{results['lower']}-{results['upper']}")
        ax.legend()

    fig.savefig(output_path)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-path", required=True)
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.input_path, args.output_path)
