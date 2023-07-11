import json
from itertools import chain
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import AltAz
from pathlib import Path
from argparse import ArgumentParser
from tqdm import tqdm
from link_utils import build_altaz, get_pointings_of_irfs, sin_delta, cos_zenith, euclidean_distance


def main() -> None:
    prod = args.prod
    dec = args.dec

    runsummary = Table.read(args.runsummary)

    path = Path(template_irf.format(prod=prod, dec=dec, node=""))
    filelist = [p.name for p in path.parent.iterdir()]
    irf_pointings: AltAz = get_pointings_of_irfs(filelist)

    irf_sindelta = sin_delta(irf_pointings)
    irf_coszenith = cos_zenith(irf_pointings)

    progress = tqdm(total=n_runs)

    for night, run_ids in runs.items():
        for run_id in run_ids:
            target_dl1 = Path(template_target_dl1.format(night=night, run_id=run_id))

            run = runsummary[runsummary["runnumber"] == int(run_id)]

            pointing = build_altaz(
                alt=run["mean_altitude"] * u.rad,
                az=run["mean_azimuth"] * u.rad,
            )
            sindelta = sin_delta(pointing)
            coszenith = cos_zenith(pointing)

            nearest_irf = euclidean_distance(
                x1=sindelta,
                y1=coszenith,
                x2=irf_sindelta,
                y2=irf_coszenith,
            ).argmin()
            node = filelist[nearest_irf]

            linkname_dl1 = Path(
                template_linkname_dl1.format(
                    night=night,
                    run_id=run_id,
                ),
            )


            linkname_dl1.parent.mkdir(exist_ok=True, parents=True)
            if linkname_dl1.exists() and linkname_dl1.is_symlink():
                linkname_dl1.unlink()
            linkname_dl1.symlink_to(target_dl1)

            # Need to merge/split diffuse data myself and train a model!
            progress.update()

    target_nodes = Path(template_target_nodes.format(prod=prod, dec=dec, node=node))
    linkname_nodes.parent.mkdir(exist_ok=True, parents=True)
    if linkname_nodes.exists() and linkname_nodes.is_symlink():
        linkname_nodes.unlink()
    linkname_nodes.symlink_to(target_nodes)

    target_protons_dl1 = Path(template_target_protons_dl1.format(prod=prod, dec=dec))
    linkname_protons_dl1.parent.mkdir(exist_ok=True, parents=True)
    if linkname_protons_dl1.exists() and linkname_protons_dl1.is_symlink():
        linkname_protons_dl1.unlink()
    linkname_protons_dl1.symlink_to(target_protons_dl1)

    # to get the config
    target_model = Path(template_target_model.format(prod=prod, dec=dec))
    linkname_model = template_linkname_model
    linkname_model.parent.mkdir(exist_ok=True, parents=True)
    print(linkname_model, linkname_model.exists(), linkname_model.is_symlink())
    if linkname_model.exists() and linkname_model.is_symlink():
        linkname_model.unlink()
    linkname_model.symlink_to(target_model)

    Path(args.output_path).touch()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--runs", required=True)
    parser.add_argument("--prod", required=True)
    parser.add_argument("--dec", required=True)
    parser.add_argument("--runsummary", required=True)
    parser.add_argument("-o", "--output-path", required=True)
    args = parser.parse_args()

    with open(args.runs, "r") as f:
        runs = json.load(f)
    n_runs = len(set(chain(*runs.values())))

    build_dir = Path(args.output_path).parent
    outdir_dl1 = build_dir / "dl1/"
    filename_dl1 = "dl1_LST-1.Run{run_id}.h5"
    template_target_dl1 = (
        (Path("/fefs/aswg/data/real/DL1/{night}/v0.9/tailcut84") / filename_dl1)
        .resolve()
        .as_posix()
    )
    template_linkname_dl1 = (outdir_dl1 / filename_dl1).resolve().as_posix()

    # This is what changes compared to 1D
    # We link unmerged dl1 nodes here!
    outdir_nodes = build_dir / "mc_nodes"
    template_target_nodes = "/fefs/aswg/data/mc/DL1/AllSky/{prod}/TrainingDataset/{dec}"  # noqa
    linkname_nodes = outdir_nodes

    #filename_irf = "irf_Run{run_id}.fits.gz"
    template_irf = "/fefs/aswg/data/mc/IRF/AllSky/{prod}/TestingDataset/{dec}/{node}/irf_{prod}_{node}.fits.gz"  # noqa

    template_target_protons_dl1 = "/fefs/aswg/data/mc/DL1/AllSky/{prod}/TrainingDataset/{dec}/Protons/dl1_{prod}_{dec}_Protons_merged.h5"
    linkname_protons_dl1 = outdir_dl1 / "train/proton_diffuse_merged.dl1.h5"

    # Just link to get the proper configs for training
    outdir_model = build_dir / "models/mcpipe/"
    template_target_model = "/fefs/aswg/data/models/AllSky/{prod}/{dec}/"
    template_linkname_model = outdir_model

    main()
