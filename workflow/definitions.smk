import json
import os
from pathlib import Path

# Unfortunately I need these as parameters for the rules, so I define them here
main_config_path = (config_dir / "lst_agn.json").resolve()
with open(main_config_path, "r") as f:
    config_agn = json.load(f)
PRODUCTION = config_agn["production"]
DECLINATION = config_agn["declination"]
# Set all configs
# .resolve(), because I am paranoid. Thats a TODO kinda
# It seems to be, that snakemake changes the current working directory and then local paths fail
# Per Analysis configs are not in here. They have fixed names in the analysis-* dirs
config_dir = Path(config.get("config_dir", "../lst-analysis-config"))
MATPLOTLIBRC = os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc")
CONFIGS = {
    "agn": (config_dir / "lst_agn.json").resolve(),
    "data_selection": (config_dir / "data-selection.json").resolve(),
    "irf_tool": (config_dir / "irf_tool_config.json").resolve(),
    "bkg_model": (config_dir / "bkgmodel.yml").resolve(),
}

build = Path("build") / config_dir.name
OUTDIRS = {
    "data_selection": (build / "data_selection").resolve(),
    "mc_nodes": (
        build / "mc_nodes"
    ).resolve(),  # This is not really configurable, its hard coded in the link script
    "mc": (build / "mc").resolve(),
    "models": (build / "models").resolve(),
    "dl1": (build / "dl1").resolve(),
    "dl2": (build / "dl2").resolve(),
    "dl3": (build / "dl3").resolve(),
    "irfs": (build / "irfs").resolve(),
    "dl4": (build / "dl4").resolve(),
    "dl5": (build / "dl5").resolve(),
}

# Set all enviroments
env_dir = Path("workflow/envs")
ENVS = {
    "lstchain": config_agn.get("lstchain_enviroment", "lstchain-v0.9.13"),
    "gammapy": (env_dir / "agn_analysis.yml").resolve(),
    "run_selection": (env_dir / "run_selection.yml").resolve(),
    "background": (env_dir / "background.yml").resolve(),
    "dark_matter": (env_dir / "dark_matter.yml").resolve(),
}

scripts_dir = Path("scripts")
SCRIPTS = {
    "data_selection": (build / "data_selection").resolve(),
    "mc": (scripts_dir / "mc").resolve(),
    "dl1": (scripts_dir / "dl1").resolve(),
    "dl2": (scripts_dir / "dl2").resolve(),
    "irfs": (scripts_dir / "irfs").resolve(),
    "dl3": (scripts_dir / "dl3").resolve(),
    "dl4": (scripts_dir / "dl4").resolve(),
    "fit": (scripts_dir / "fit").resolve(),
    "dm": (scripts_dir / "dm").resolve(),
}


# TODO This is the most critical part as the further evaluation depends on this checkpoint
# Have to make sure this works as expected
def RUN_IDS(wildcards):
    with open(checkpoints.run_ids.get(**wildcards).output, "r") as f:
        runs = json.load(f)
    return sorted(set(chain(*runs.values())))


models_to_train = [
    Path(OUTDIRS["models"]) / "reg_energy.sav",
    Path(OUTDIRS["models"]) / "cls_gh.sav",
    Path(OUTDIRS["models"]) / "reg_disp_norm.sav",
    Path(OUTDIRS["models"]) / "cls_disp_sign.sav",
]

# TODO: cuts are not really IRFs, should separate that.
# Add radmax here if 1D
irfs_to_produce = ["aeff", "gh_cut", "edisp", "psf"]  # TODO script missing
