from argparse import ArgumentParser

import numpy as np
from gammapy.irf import load_irf_dict_from_file
from matplotlib import pyplot as plt

parser = ArgumentParser()
parser.add_argument("-i", "--input-path", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()


def main(input_path, output):
    aeff = load_irf_dict_from_file(input_path)["aeff"]

    e_true = aeff.axes["energy_true"]

    offsets = aeff.axes["offset"]
    unit = offsets.unit

    off_centers = np.atleast_1d(offsets.center)
    off_widths = np.atleast_1d(offsets.bin_width)

    fig, ax = plt.subplots()
    for i, (off_center, off_width) in enumerate(zip(off_centers, off_widths)):
        alpha = i / len(off_centers)
        label = r" \pm ".join(
            [
                f"${off_center.to_value(unit)}",
                f"{off_width.to_value(unit) / 2:.1f}$ {unit}",
            ],
        )
        ax.errorbar(
            e_true.center,
            aeff.data[:, i],
            xerr=e_true.bin_width / 2,
            ls="",
            label=f"Offset: {label}",
            color="black",
            alpha=alpha,
        )
    ax.set_xlabel(rf"$E_{{\mathrm{{true}}}}$ / {e_true.center.unit}")
    ax.set_ylabel(f"Effective Area / {aeff.unit}")

    ax.set_xscale("log")
    ax.set_yscale("log")

    fig.savefig(output)


if __name__ == "__main__":
    main(**vars(args))
