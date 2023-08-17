import json
from pathlib import Path

# Unfortunately I need these as parameters for the rules, so I define them here
with open(main_config_path, "r") as f:
    config_agn = json.load(f)
PRODUCTION = config_agn["production"]
DECLINATION = config_agn["declination"]


# Set all configs
# .resolve(), because I am paranoid
# It seems to be, that snakemake changes the current working directory and then local paths fail
# Per Analysis configs are not in here. They have fixed names in the analysis-* dirs
config_dir = Path(config.get("config_dir", "../lst-analysis-config"))
CONFIGS = {
    "agn": (config_dir / "lst_agn.json").resolve(),
    "data_selection": (config_dir / "data-selection.json").resolve(),
    "irf_tool": (config_dir / "irf_tool_config.json").resolve(),
    "bkg_model": (config_dir / "bkgmodel.yml").resolve()
}

build = Path("build") / config_dir.name

# Set all enviroments
env_dir = Path("workflow/envs")
ENVS = {
    "lstchain": config_agn.get("lstchain_enviroment", "lstchain-v0.9.13"), 
    "gammapy": (env_dir / "agn_analysis.yml").resolve(), 
    "run_selection": (env_dir / "drun_selection.yml").resolve(),
    "background_creation": (env_dir / "background.yml").resolve(),
    "dark_matter": (env_dir / "dark_matter.yml").resolve(),
}

scripts_dir = Path("scripts")
SCRIPTS = {
    "run_selection": (scripts_dir / "run_selection").resolve(),
    "link_runs": (scripts_dir / "link_runs").resolve(),
    "mc": (scripts_dir / "mc").resolve(),
    "dl1": (scripts_dir / "dl1").resolve(),
    "dl2": (scripts_dir / "dl2").resolve(),
    "irfs": (scripts_dir / "irfs").resolve(),
    "dl3": (scripts_dir / "dl3").resolve(),
    "dl4": (scripts_dir / "dl4").resolve(),
    "fit": (scripts_dir / "fit").resolve(),
    "dm": (scripts_dir / "dm").resolve(),
}


