env = ENVS["lstchain"]
bkg_env = ENVS["background"]

dl2 = Path(OUTDIRS["dl2"])
dl3 = Path(OUTDIRS["dl3"])
plots = dl3 / "plots"
irfs = Path(OUTDIRS["irfs"])

scripts = Path(SCRIPTS["dl3"])

irf_config = CONFIGS["irf_tool"]
bkg_config = CONFIGS["bkg_model"]


rule dl3:
    input:
        index=dl3 / "hdu-index.fits.gz",
        bkg=dl3 / "bkg-exists",
        plots=DL3_PLOTS,


rule dl2_to_dl3:
    output:
        run=dl3 / "LST-1.Run{run_id}.dl3.fits.gz",
    input:
        data=dl2 / "LST-1.Run{run_id}.dl2.h5",
        irfs=MC_NODES_IRFs,
        config=irf_config,
    params:
        irf_pattern="irfs_*.fits.gz",
        out=dl3,
        in_irfs=irfs,
    conda:
        env
    resources:
        mem_mb=12000,
        time=30,
    log:
        dl3 / "create_dl3_{run_id}.log",
    shell:
        """
        lstchain_create_dl3_file  \
            --input-dl2 {input.data}  \
            --output-dl3-path {params.out}  \
            --input-irf-path {params.in_irfs}  \
            --irf-file-pattern {params.irf_pattern} \
            --config {input.config} \
            --gzip \
            --use-nearest-irf-node \
            --overwrite \
            --log-file {log}
        """


rule dl3_hdu_index:
    output:
        dl3 / "hdu-index.fits.gz",
    input:
        runs=DL3_FILES,
    params:
        outdir=dl3,
    conda:
        env
    log:
        dl3 / "hdu_index.log",
    resources:
        time=15,
    shell:
        """
        lstchain_create_dl3_index_files  \
            --input-dl3-dir {params.outdir}  \
            --output-index-path {params.outdir}  \
            --file-pattern '*.dl3.fits.gz'  \
            --overwrite \
            --log-file {log}
        """


rule calc_count_maps:
    output:
        dl3 / "bkg_cached_maps.pkl",
    input:
        runs=DL3_FILES,
        index=dl3 / "hdu-index.fits.gz",
        config=bkg_config,
        script=scripts / "precompute_background_maps.py",
    params:
        obs_dir=dl3,
        bkg_dir=dl3,
    conda:
        bkg_env
    resources:
        partition="long",
        time=360,
    log:
        dl3 / "calc_count_maps.log",
    shell:
        """python {input.script} \
        --input-dir {params.obs_dir} \
        --output {output} \
        --config {input.config} \
        --log-file {log} \
        --overwrite
        """


rule calc_background:
    output:
        #        bkg=expand(dl3 / "bkg_{run_id}.fits.gz", run_id=RUN_IDS), # cant use function in output
        dummy=dl3 / "bkg-exists",
    input:
        runs=DL3_FILES,
        index=dl3 / "hdu-index.fits.gz",
        config=bkg_config,
        calc_script=scripts / "calc_background.py",
        link_script=scripts / "link_bkg.py",
        cached_maps=dl3 / "bkg_cached_maps.pkl",
    params:
        obs_dir=dl3,
        bkg_dir=dl3,
    conda:
        bkg_env
    resources:
        partition="long",
        time=360,
    log:
        dl3 / "calc_bkg.log",
    shell:
        """python {input.calc_script} \
        --input-dir {params.obs_dir} \
        --output-dir {params.bkg_dir} \
        --dummy-output {output.dummy} \
        --cached-maps {input.cached_maps} \
        --config {input.config} \
        --log-file {log} \
        --overwrite

        python {input.link_script} \
        --hdu-index-path {input.index} \
        --bkg-dir {params.bkg_dir} \
        """


# Use lstchain env here to ensure we can load it
rule plot_background:
    output:
        dl3 / "plots/bkg_{run_id}.pdf",
    input:
        data=dl3 / "bkg_{run_id}.fits.gz",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
        script=scripts / "plot_bkg.py",
    conda:
        env
    log:
        dl3 / "plots/bkg_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule observation_plots:
    input:
        dl3 / "hdu-index.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script=scripts / "obs_plots.py",
    output:
        plots / "{analysis}/obs_plots.pdf",
    resources:
        mem_mb=16000,
    conda:
        gammapy_env
    log:
        plots / "{analysis}/obs_plots.log",
    shell:
        """
        python {input.script} \
            -c {input.config} \
            -o {output} \
            --log-file {log}
        """


rule calc_theta2_per_obs:
    output:
        plots / "theta2/{run_id}.fits.gz",
    input:
        data=dl3 / "LST-1.Run{run_id}.dl3.fits.gz",
        script=scripts / "calc_theta2_per_obs.py",
        config=data_selection_config_path,  # this seems unnecessary
        index=dl3 / "hdu-index.fits.gz",
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    resources:
        mem_mb=16000,
    conda:
        gammapy_env
    log:
        dl3 / "theta2/calc_{run_id}.log",
    shell:
        "python {input.script} -i {dl3} -o {output} --obs-id {wildcards.run_id} --config {input.config} --log-file {log}"


rule stack_theta2:
    output:
        dl3 / "theta2/stacked.fits.gz",
    input:
        runs=expand(
            dl3 / "theta2/{run_id}.fits.gz",
            run_id=RUN_IDS,
        ),
        script=scripts / "stack_theta2.py",
    conda:
        gammapy_env
    log:
        dl3 / "theta2/theta2_stacked.log",
    shell:
        "python {input.script} -o {output} --input-files {input.runs} --log-file {log}"


rule plot_theta:
    output:
        plots / "theta2/{runid}.pdf",
    input:
        data=dl3 / "theta2/{runid}.fits.gz",
        script="scripts/plot_theta2.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    log:
        dl3 / "theta2/plot_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"


rule calc_skymap:
    output:
        dl3 / "skymaps/{run_id}.fits",
    input:
        data=dl3 / "LST-1.Run{run_id}.dl3.fits.gz",
        script=scripts / "calc_skymap_gammas.py",
        config=irf_config_path,
        index=dl3 / "hdu-index.fits.gz",
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    resources:
        # mem_mb=16000,
        time=5,
    conda:
        gammapy_env
    log:
        dl3 / "skymaps/calc_{run_id}.log",
    shell:
        "python {input.script} -i {dl3} -o {output} --obs-id {wildcards.run_id} --config {input.config} --log-file {log}"


rule plot_skymap:
    output:
        build_dir / "plots/skymap_dl3/{runid}.pdf",
    input:
        data=dl3 / "skymap_dl3/{runid}.fits",
        script=scripts / "plot_skymap_dl3.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    resources:
        time=5,
    log:
        dl3 / "skymaps/plot_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"


rule stack_skymaps:
    conda:
        lstchain_env
    output:
        dl3 / "skymap_dl3/stacked.fits",
    input:
        data=expand(
            dl3 / "skymap_dl3/{run_id}.fits",
            run_id=RUN_IDS,
        ),
        script="scripts/stack_skymap.py",
    log:
        dl3 / "skymaps/stack.log",
    shell:
        "python {input.script} -i {input.data} -o {output} --log-file {log}"


rule cuts_dl2_dl3:
    output:
        dl3 / "counts/after_gh_theta_cut_{run_id}.h5",
    input:
        dl2=dl2 / "dl2_LST-1.Run{run_id}.h5",
        irf=irfs / "irf_Run{run_id}.fits.gz",
        config=irf_config_path,
        script=scripts / "calc_counts_after_cuts.py",
    resources:
        mem_mb="64G",
        time=10,
    conda:
        lstchain_env
    log:
        dl3 / "counts/calc_{run_id}.log",
    shell:
        "python {input.script} --input-dl2 {input.dl2} --input-irf {input.irf} -c {input.config} -o {output} --log-file {log}"


rule stack_cuts_dl2_dl3:
    output:
        dl3 / "counts/after_gh_theta_cut_stacked.h5",
    input:
        data=expand(
            dl3 / "counts/after_gh_theta_cut_{run_id}.h5",
            run_id=RUN_IDS,
        ),
        script=scripts / "stack_counts_after_cuts.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        lstchain_env
    log:
        dl3 / "counts/stack.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"


rule plot_cuts_dl2_dl3:
    output:
        plots / "counts_after_cuts.pdf",
    input:
        data=dl3 / "counts/after_gh_theta_cut.h5",
        script=scripts / "plot_counts_after_cuts.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        lstchain_env
    log:
        dl3 / "counts/plot_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"
