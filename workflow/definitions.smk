import json
import os
from pathlib import Path

# "Main" paths. Everuthing else is relative to these
scripts_dir = Path("scripts")
env_dir = Path("workflow/envs")
config_dir = Path(config.get("config_dir", "../lst-analysis-config"))
main_config_path = (config_dir / "lst_agn.json").absolute()
build = Path("build") / config_dir.name
analyses = [
    x.name for x in config_dir.iterdir() if x.name.startswith("analysis") and x.is_dir()
]

# Using the one used for the dl1 MCs as I dont reprocess anything
# Best might be to reproduce dl1 mc and data?
# For now this is fine as mcpipe uses the default config from lstchain
# (at least in my test production) and lstchain versions dont differ much in the output
# (hopefully). Thats a fat TODO though...
lstchain_config = build / "lstchain_config.json"
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
    "dl1": (build / "dl1").absolute(),
    "dl2": (build / "dl2").absolute(),
    "dl3": (build / "dl3").absolute(),
    "irfs": (build / "irfs").absolute(),
    "dl4": (build / "dl4").absolute(),
    "dl5": (build / "dl5").absolute(),
}

# Set all enviroments
ENVS = {
    "lstchain": config_agn.get("lstchain_enviroment", "lstchain-v0.10.3"),
    "gammapy": (env_dir / "agn-analysis.yml").absolute(),
    "data_selection": (env_dir / "data-selection.yml").absolute(),
    "background": (env_dir / "background.yml").absolute(),
    "dark_matter": (env_dir / "dark_matter.yml").absolute(),
}

SCRIPTS = {
    "data_selection": (scripts_dir / "data_selection").absolute(),
    "mc": (scripts_dir / "mc").absolute(),
    "irfs": (scripts_dir / "irfs").absolute(),
    "dl1": (scripts_dir / "dl1").absolute(),
    "dl2": (scripts_dir / "dl2").absolute(),
    "dl3": (scripts_dir / "dl3").absolute(),
    "dl4": (scripts_dir / "dl4").absolute(),
    "dl5": (scripts_dir / "dl5").absolute(),
    "dm": (scripts_dir / "dm").absolute(),
}


# TODO This is the most critical part as the further evaluation depends on this checkpoint
# Have to make sure this works as expected
def RUN_IDS(wildcards):
    with open(checkpoints.run_ids.get(**wildcards).output.runlist, "r") as f:
        runs = json.load(f)
    return sorted(set(chain(*runs.values())))


def MC_NODES(wildcards):
    exists = Path(checkpoints.link_mc.get(**wildcards).output.dummy).exists()
    mc_nodes = Path(OUTDIRS["mc_nodes"]) / f"{wildcards.particle}"
    nodes = [x.name for x in mc_nodes.glob("*") if x.is_dir()]
    return nodes


def MC_NODES_DL1(wildcards):
    out = Path(OUTDIRS["mc"]) / f"{wildcards.particle}/dl1"
    nodes = MC_NODES(wildcards)
    return [out / f"{node}_{wildcards.train_or_test}.dl1.h5" for node in nodes]


def MC_NODES_DL2(wildcards):
    out = Path(OUTDIRS["mc"]) / f"{wildcards.particle}/dl2"
    nodes = MC_NODES(wildcards)
    return [out / f"{node}_{wildcards.train_or_test}.dl2.h5" for node in nodes]


# TODO: cuts are not really IRFs, should separate that.
# Add radmax here if 1D
irfs_to_produce = ["aeff", "gh_cut", "edisp", "psf"]  # TODO script missing


def MC_NODES_IRFs(wildcards):
    exists = Path(checkpoints.link_mc.get(**wildcards).output.dummy).exists()
    out = Path(OUTDIRS["irfs"]) / "plots"
    mc_nodes = Path(OUTDIRS["mc_nodes"]) / "GammaDiffuse"
    nodes = [x.name for x in mc_nodes.glob("*") if x.is_dir()]
    return [out / f"{irf}_{node}.pdf" for node in nodes for irf in irfs_to_produce]


models_to_train = [
    Path(OUTDIRS["models"]) / "reg_energy.sav",
    Path(OUTDIRS["models"]) / "cls_gh.sav",
    Path(OUTDIRS["models"]) / "reg_disp_norm.sav",
    Path(OUTDIRS["models"]) / "cls_disp_sign.sav",
]


def DL2_FILES(wildcards):
    ids = RUN_IDS(wildcards)
    out = Path(OUTDIRS["dl2"])
    return [out / f"LST-1.Run{run_id}.dl2.h5" for run_id in ids]


def DL3_FILES(wildcards):
    ids = RUN_IDS(wildcards)
    out = Path(OUTDIRS["dl3"])
    return [out / f"LST-1.Run{run_id}.dl3.fits.gz" for run_id in ids]


def IRF_FILES(wildcards):
    ids = RUN_IDS(wildcards)
    out = Path(OUTDIRS["dl3"])
    return [out / f"irfs_{run_id}.fits.gz" for run_id in ids]


def BKG_FILES(wildcards):
    ids = RUN_IDS(wildcards)
    out = Path(OUTDIRS["dl3"])
    return [out / f"bkg_{run_id}.fits.gz" for run_id in ids]
