import logging
from argparse import ArgumentParser
from pathlib import Path

import joblib
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from lstchain.io import (
    read_configuration_file,
    replace_config,
    standard_config,
)
from lstchain.io.io import (
    read_dl2_params,
)
from lstchain.visualization import plot_dl2

if matplotlib.get_backend() == "pgf":
    from matplotlib.backends.backend_pgf import PdfPages
else:
    from matplotlib.backends.backend_pdf import PdfPages


def main():  # noqa
    parser = ArgumentParser()
    parser.add_argument("--gammas", required=True)
    parser.add_argument("--protons", required=True)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("-o", "--output-path", required=True)
    parser.add_argument("--log-file")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    # TODO basic logging replacememtn
    log = logging.getLogger(__name__)
    log.level = logging.DEBUG if args.verbose else logging.INFO
    if args.log_file is not None:
        file_handler = logging.FileHandler(args.log_file)
        file_formatter = logging.Formatter(
            fmt="%(asctime)s|%(levelname)s|%(name)s|%(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S",
        )
        file_handler.setFormatter(file_formatter)
        log.addHandler(file_handler)

    # TODO load data and models
    custom_config = read_configuration_file(args.config)
    config = replace_config(standard_config, custom_config)
    #
    models_dict = {"reg_energy": joblib.load(Path(args.model_dir, "reg_energy.sav"))}
    models_dict["cls_gh"] = joblib.load(Path(args.model_dir, "cls_gh.sav"))
    if config["disp_method"] == "disp_vector":
        models_dict["disp_vector"] = joblib.load(
            Path(args.model_dir, "reg_disp_vector.sav"),
        )
    elif config["disp_method"] == "disp_norm_sign":
        models_dict["disp_norm"] = joblib.load(
            Path(args.model_dir, "reg_disp_norm.sav"),
        )
        models_dict["disp_sign"] = joblib.load(
            Path(args.model_dir, "cls_disp_sign.sav"),
        )
    # data. does this work or do i need to load just the feature and fdo it by hand?
    gammas = read_dl2_params(args.gammas)
    selected_gammas = gammas.query("reco_type==0 & mc_type==0")
    protons = read_dl2_params(args.protons)
    protons.query("reco_type==0 & mc_type==0")
    dl2 = pd.concat([gammas, protons], ignore_index=True)

    # lstchain does not foresee us passing axes or figure here :(
    figures = []
    plot_dl2.plot_features(dl2)
    fig = plt.gcf()
    fig.tight_layout()
    figures.append(fig)

    plot_dl2.energy_results(selected_gammas)
    figures.append(plt.gcf())

    plot_dl2.direction_results(selected_gammas)
    figures.append(plt.gcf())

    plot_dl2.plot_disp_vector(selected_gammas)
    fig = plt.gcf()
    fig.tight_layout()
    figures.append(fig)

    plot_dl2.plot_pos(dl2)
    figures.append(plt.gcf())

    fig, ax = plt.subplots()
    plot_dl2.plot_roc_gamma(dl2, ax=ax)
    figures.append(fig)

    plot_dl2.plot_models_features_importances(args.model_dir, args.config)
    figures.append(plt.gcf())

    fig, ax = plt.subplots()
    ax.hist(
        dl2[dl2["mc_type"] == 101]["gammaness"],  # noqa
        bins=100,
        histtype="step",
        density=True,
    )
    ax.hist(
        dl2[dl2["mc_type"] == 0]["gammaness"],
        bins=100,
        histtype="step",
        density=True,
    )
    figures.append(fig)

    with PdfPages(args.output_path) as pdf:
        for fig in figures:
            pdf.savefig(fig)


if __name__ == "__main__":
    main()
