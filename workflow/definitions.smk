import json
import os
from pathlib import Path

# "Main" paths. Everuthing else is relative to these
scripts_dir = Path("scripts")
env_dir = Path("workflow/envs")
config_dir = Path(config.get("config_dir", "../lst-analysis-config"))
main_config_path = (config_dir / "lst_agn.json").absolute()
build = Path("build") / config_dir.name
MATPLOTLIBRC = os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc")

with open(main_config_path, "r") as f:
    config_agn = json.load(f)
PRODUCTION = config_agn["production"]
DECLINATION = config_agn["declination"]

# Set all configs
# .absolute(), because I am paranoid. Thats a TODO kinda
# It seems to be, that snakemake changes the current working directory and then local paths fail
# Per Analysis configs are not in here. They have fixed names in the analysis-* dirs
CONFIGS = {
    "agn": (config_dir / "lst_agn.json").absolute(),
    "data_selection": (config_dir / "data_selection.json").absolute(),
    "irf_tool": (config_dir / "irf_tool_config.json").absolute(),
    "bkg_model": (config_dir / "bkgmodel.yml").absolute(),
}

OUTDIRS = {
    "data_selection": (build / "data_selection").absolute(),
    "mc_nodes": (
        build / "mc_nodes"
    ).absolute(),  # This is not really configurable, its hard coded in the link script
    "mc": (build / "mc").absolute(),
    "models": (build / "models").absolute(),
    "models_link": (build / "models/mcpipe").absolute(),
    "dl1": (build / "dl1").absolute(),
    "dl2": (build / "dl2").absolute(),
    "dl3": (build / "dl3").absolute(),
    "irfs": (build / "irfs").absolute(),
    "dl4": (build / "dl4").absolute(),
    "dl5": (build / "dl5").absolute(),
}

# Set all enviroments
ENVS = {
    "lstchain": config_agn.get("lstchain_enviroment", "lstchain-v0.9.13"),
    "gammapy": (env_dir / "agn_analysis.yml").absolute(),
    "data_selection": (env_dir / "data-selection.yml").absolute(),
    "background": (env_dir / "background.yml").absolute(),
    "dark_matter": (env_dir / "dark_matter.yml").absolute(),
}

SCRIPTS = {
    "data_selection": (scripts_dir / "data_selection").absolute(),
    "mc": (scripts_dir / "mc").absolute(),
    "dl1": (scripts_dir / "dl1").absolute(),
    "dl2": (scripts_dir / "dl2").absolute(),
    "irfs": (scripts_dir / "irfs").absolute(),
    "dl3": (scripts_dir / "dl3").absolute(),
    "dl4": (scripts_dir / "dl4").absolute(),
    "fit": (scripts_dir / "fit").absolute(),
    "dm": (scripts_dir / "dm").absolute(),
}


# TODO This is the most critical part as the further evaluation depends on this checkpoint
# Have to make sure this works as expected
def RUN_IDS(wildcards):
    with open(checkpoints.run_ids.get(**wildcards).output, "r") as f:
        runs = json.load(f)
    return sorted(set(chain(*runs.values())))


# SHould be the same for protons
def MC_NODES(wildcards):
    exists = Path(checkpoints.link_mc.get(**wildcards).output).exists()
    nodes = [x.name for x in (mc_nodes / "GammaDiffuse").glob("*") if x.is_dir()]
    return expand(
        mc / "{{wildcards.particle}}/{node}_{{wildcards.train_or_test}}.dl1.h5",
        node=nodes,
    )


models_to_train = [
    Path(OUTDIRS["models"]) / "reg_energy.sav",
    Path(OUTDIRS["models"]) / "cls_gh.sav",
    Path(OUTDIRS["models"]) / "reg_disp_norm.sav",
    Path(OUTDIRS["models"]) / "cls_disp_sign.sav",
]

# TODO: cuts are not really IRFs, should separate that.
# Add radmax here if 1D
irfs_to_produce = ["aeff", "gh_cut", "edisp", "psf"]  # TODO script missing
