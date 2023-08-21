env = ENVS["lstchain"]
bkg_env = ENVS["background"]

# Having these as paths makes them easier to use
dl2 = Path(OUTDIRS["dl2"])
dl3 = Path(OUTDIRS["dl3"])
irfs = Path(OUTDIRS["irfs"])
bkg = Path(OUTDIRS["bkg"])
models = Path(OUTDIRS["models"])

irf_config = CONFIGS["irf_tool"]
bkg_config = CONFIGS["bkg_model"]

scripts = Path(SCRIPTS["dl3"])

plots = Path(PLOTSDIRS["dl3"])


rule dl3:
    output:
        dl3 / "dl3_LST-1.Run{run_id}.fits.gz",
    input:
        data=dl2 / "dl2_LST-1.Run{run_id}.h5",
        irf=irfs / "irf_Run{run_id}.fits.gz",
        config=irf_config,
    conda:
        env
    resources:
        mem_mb=12000,
        time=30,
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        """
        lstchain_create_dl3_file  \
            --input-dl2 {input.data}  \
            --output-dl3-path $(dirname $(realpath {output}))  \
            --input-irf {input.irf}  \
            --config {input.config} \
            --gzip \
            --overwrite \
        """


# Using my fork here currently


# TODO Write my own script.
# no clear way to swap between runwise and stacked in my workflow :/
# maybe define a function, that returns the corresponding bkg name to a run based on
# a variable. I would need to parse that from the bkgmodel config and
# also get it into the link bkg script ...
# build_dir / "background/stacked_bkg_map.fits"
# dont ask... result of my hacks, should be solved later upstream
rule calc_background:
    output:
        expand(
            dl3 / "dl3_LST-1.Run{run_id}.fits.fits",
            run_id=RUN_IDS,
        ),
    input:
        runs=expand(
            dl3 / "dl3_LST-1.Run{run_id}.fits.gz",
            run_id=RUN_IDS,
        ),
        config=bkg_config,
    conda:
        bkg_env
    log:
        out=dl3 / "calc_bkg.log",
        err=dl3 / "calc_bkg.err",
    shell:
        """
        bkgmodel --config {input.config}
        """


# Use lstchain env here to ensure we can load it
# run id is stacked only right now, but this way it can be expanded
# data=build_dir / "background/{run_id}_bkg_map.fits", # stacked
rule plot_background:
    output:
        plots / "background/{run_id}.pdf",
    input:
        data=background / "dl3_LST-1.Run{run_id}.fits.fits",  #runwise
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
        script=scripts / "plot_bkg.py",
    conda:
        env
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


# bkg = build_dir / "background/stacked_bkg_map.fits"
# dont ask... result of my hacks, should be solved later upstream
rule dl3_hdu_index:
    output:
        dl3 / "hdu-index.fits.gz",
    input:
        runs=expand(
            dl3 / "dl3_LST-1.Run{run_id}.fits.gz",
            run_id=RUN_IDS,
        ),
        bkg=expand(
            dl3 / "dl3_LST-1.Run{run_id}.fits.fits",
            run_id=RUN_IDS,
        ),
    params:
        bkg_script=scripts / "link_bkg.py",
        bkg_dir=lambda w, input: os.path.relpath(
            Path(input.bkg[0]).parent, Path(input.runs[0]).parent
        ),
        bkg_files=lambda w, input: [Path(x).name for x in input.bkg],
    conda:
        env
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    resources:
        time=15,
    shell:
        """
        lstchain_create_dl3_index_files  \
            --input-dl3-dir {build_dir}/dl3  \
            --output-index-path {build_dir}/dl3  \
            --file-pattern 'dl3_*.fits.gz'  \
            --overwrite 

        python {params.bkg_script} \
        --hdu-index-path {output} \
        --bkg-dir {params.bkg_dir} \
        --bkg-file {params.bkg_files} 
        """


# Plots using dl3 files
rule observation_plots:
    input:
        build_dir / "dl3/hdu-index.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script="scripts/events.py",
    output:
        build_dir / "plots/{analysis}/observation_plots.pdf",
    resources:
        mem_mb=64000,
    conda:
        gammapy_env
    shell:
        """
        python {input.script} \
            -c {input.config} \
            -o {output} \
        """


rule calc_theta2_per_obs:
    output:
        build_dir / "dl3/theta2/{run_id}.fits.gz",
    input:
        data=build_dir / "dl3/dl3_LST-1.Run{run_id}.fits.gz",
        script="scripts/calc_theta2_per_obs.py",
        config=data_selection_config_path,
        index=build_dir / "dl3/hdu-index.fits.gz",
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    resources:
        mem_mb=16000,
    conda:
        gammapy_env
    log:
        build_dir / "logs/dl3/theta2/{run_id}.log",
    shell:
        "python {input.script} -i {build_dir}/dl3 -o {output} --obs-id {wildcards.run_id} --config {input.config} --log-file {log}"


rule stack_theta2:
    output:
        build_dir / "dl3/theta2/stacked.fits.gz",
    input:
        runs=expand(
            build_dir / "dl3/theta2/{run_id}.fits.gz",
            run_id=RUN_IDS,
        ),
        script="scripts/stack_theta2.py",
    conda:
        gammapy_env
    log:
        build_dir / "logs/dl3/theta2_stacked.log",
    shell:
        "python {input.script} -o {output} --input-files {input.runs} --log-file {log}"


rule plot_theta:
    output:
        build_dir / "plots/theta2/{runid}.pdf",
    input:
        data=build_dir / "dl3/theta2/{runid}.fits.gz",
        script="scripts/plot_theta2.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule bkg_exclusion:
    input:
        config=config_dir / "{analysis}/analysis.yaml",
        script="scripts/make_exclusion_region.py",
    output:
        build_dir / "{analysis}/exclusion.fits.gz",
    conda:
        gammapy_env
    shell:
        "python {input.script} -c {input.config} -o {output}"


rule calc_skymap_per_obs:
    output:
        build_dir / "dl3/skymap_dl3/{run_id}.fits",
    input:
        data=build_dir / "dl3/dl3_LST-1.Run{run_id}.fits.gz",
        script="scripts/calc_skymap_gammas.py",
        config=irf_config_path,
        index=build_dir / "dl3/hdu-index.fits.gz",
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    resources:
        # mem_mb=16000,
        time=5,
    conda:
        gammapy_env
    shell:
        "python {input.script} -i {build_dir}/dl3 -o {output} --obs-id {wildcards.run_id} --config {input.config}"


rule plot_skymap_dl3:
    output:
        build_dir / "plots/skymap_dl3/{runid}.pdf",
    input:
        data=build_dir / "dl3/skymap_dl3/{runid}.fits",
        script="scripts/plot_skymap_dl3.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    resources:
        time=5,
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule calc_skymap:
    resources:
        mem_mb="64G",
        time=10,
    conda:
        lstchain_env
    output:
        build_dir / "dl3/skymap/{run_id}.fits",
    input:
        data=build_dir / "dl2/dl2_LST-1.Run{run_id}.h5",  #?????????
        config=irf_config_path,
        script="scripts/calc_skymap.py",
    shell:
        "python {input.script} -i {input.data} -o {output} -c {input.config}"


rule plot_skymap:
    conda:
        lstchain_env
    output:
        build_dir / "plots/skymap/{run_id}.pdf",
    input:
        data=build_dir / "dl3/skymap/{run_id}.fits",
        script="scripts/plot_skymap.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule stack_skymaps:
    conda:
        lstchain_env
    output:
        build_dir / "dl3/skymap/stacked.fits",
    input:
        data=expand(
            build_dir / "dl3/skymap/{run_id}.fits",
            run_id=RUN_IDS,
        ),
        script="scripts/stack_skymap.py",
    shell:
        "python {input.script} -i {input.data} -o {output}"


rule stack_skymaps_dl3:
    conda:
        lstchain_env
    output:
        build_dir / "dl3/skymap_dl3/stacked.fits",
    input:
        data=expand(
            build_dir / "dl3/skymap_dl3/{run_id}.fits",
            run_id=RUN_IDS,
        ),
        script="scripts/stack_skymap.py",
    shell:
        "python {input.script} -i {input.data} -o {output}"
