from argparse import ArgumentParser

from astropy.table import Table
from matplotlib import pyplot as plt


def main():
    pointings = Table.read(args.input_path)

    # TODO 2 plots (also azimuth)
    fig, ax = plt.subplots()
    ax.plot(pointings["run_id"], pointings["zen"], ls="")
    fig.savefig(args.output_path)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("-o", "--output_path", required=True)
    args = parser.parse_args()

    main()
