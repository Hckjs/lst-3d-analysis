env = ENVS["lstchain"]
bkg_env = ENVS["background"]
gammapy_env = ENVS["gammapy"]

dl2 = Path(OUTDIRS["dl2"])
dl3 = Path(OUTDIRS["dl3"])
plots = dl3 / "plots"
irfs = Path(OUTDIRS["irfs"])

scripts = Path(SCRIPTS["dl3"])

irf_config = CONFIGS["irf_tool"]
bkg_config = CONFIGS["bkg_model"]
data_selection_config = CONFIGS["data_selection"]

dl3_plot_types = ["theta2", "skymap", "counts_after_cuts", "bkg"]
# bkg plots as well!


def DL3_PLOTS(wildcards):
    ids = RUN_IDS(wildcards)
    per_run = [dl3 / f"plots/{p}/{p}_{run}.pdf" for p in dl3_plot_types for run in ids]
    return per_run + [dl3 / f"plots/{p}/{p}_stacked.pdf" for p in dl3_plot_types[:-1]]


rule dl3:
    input:
        index=dl3 / "hdu-index.fits.gz",
        bkg=dl3 / "bkg-exists",
        runwise_plots=DL3_PLOTS,


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


rule calc_count_maps:
    output:
        dl3 / "bkg_cached_maps.pkl",
    input:
        runs=DL3_FILES,
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
        --input-runs {input.runs} \
        --output {output} \
        --config {input.config} \
        --log-file {log} \
        --overwrite
        """


rule calc_background:
    output:
        dummy=dl3 / "bkg-exists",
    input:
        runs=DL3_FILES,
        config=bkg_config,
        script=scripts / "calc_background.py",
        cached_maps=dl3 / "bkg_cached_maps.pkl",
    params:
        obs_dir=dl3,
        bkg_dir=dl3,
    conda:
        bkg_env
    resources:
        partition="short",
    log:
        dl3 / "calc_bkg.log",
    shell:
        """python {input.script} \
        --input-runs {input.runs} \
        --output-dir {params.bkg_dir} \
        --dummy-output {output.dummy} \
        --cached-maps {input.cached_maps} \
        --config {input.config} \
        --log-file {log} \
        --overwrite
        """


rule dl3_hdu_index:
    output:
        dl3 / "hdu-index.fits.gz",
    input:
        runs=DL3_FILES,
        link_script=scripts / "link_bkg.py",
        dummy=dl3 / "bkg-exists",
    params:
        outdir=dl3,
        bkg=BKG_FILES,
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

        python {input.link_script} \
        --hdu-index-path {output} \
        --bkg-files {params.bkg} \
        """


rule calc_theta2_per_obs:
    output:
        dl3 / "theta2/{run_id}.fits.gz",
    input:
        data=dl3 / "LST-1.Run{run_id}.dl3.fits.gz",
        script=scripts / "calc_theta2_per_obs.py",
        config=data_selection_config,  # this seems unnecessary
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


def dl3_all_theta_tables(wildcards):
    ids = RUN_IDS(wildcards)
    return [dl3 / f"theta2/{run}.fits.gz" for run in ids]


rule stack_theta2:
    output:
        dl3 / "theta2/stacked.fits.gz",
    input:
        runs=dl3_all_theta_tables,
        script=scripts / "stack_theta2.py",
    conda:
        gammapy_env
    log:
        dl3 / "theta2/theta2_stacked.log",
    shell:
        "python {input.script} -o {output} --input-files {input.runs} --log-file {log}"


rule plot_theta:
    output:
        plots / "theta2/theta2_{run_id}.pdf",
    input:
        data=dl3 / "theta2/{run_id}.fits.gz",
        script=scripts / "plot_theta2.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        plots / "theta2/plot_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"


rule plot_background:
    output:
        plots / "bkg/bkg_{run_id}.pdf",
    input:
        data=dl3 / "bkg-exists",
        script=scripts / "plot_bkg.py",
        rc=MATPLOTLIBRC,
    params:
        data=lambda wildcards: dl3 / f"bkg_{wildcards.run_id}.fits.gz",
    conda:
        gammapy_env
    log:
        dl3 / "plots/bkg/bkg_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {params.data} -o {output}"


rule calc_skymap:
    output:
        dl3 / "skymap/{run_id}.fits.gz",
    input:
        data=dl3 / "LST-1.Run{run_id}.dl3.fits.gz",
        script=scripts / "calc_skymap_gammas.py",
        config=irf_config,
        index=dl3 / "hdu-index.fits.gz",
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    resources:
        # mem_mb=16000,
        time=5,
    conda:
        gammapy_env
    log:
        dl3 / "skymap/calc_{run_id}.log",
    shell:
        "python {input.script} -i {dl3} -o {output} --obs-id {wildcards.run_id} --config {input.config} --log-file {log}"


def dl3_all_skymaps(wildcards):
    ids = RUN_IDS(wildcards)
    return [dl3 / f"skymap/{run}.fits.gz" for run in ids]


rule stack_skymaps:
    output:
        dl3 / "skymap/stacked.fits.gz",
    input:
        data=dl3_all_skymaps,
        script=scripts / "stack_skymap.py",
    conda:
        gammapy_env
    log:
        dl3 / "skymap/stack.log",
    shell:
        "python {input.script} -i {input.data} -o {output} --log-file {log}"


rule plot_skymap:
    output:
        plots / "skymap/skymap_{run_id}.pdf",
    input:
        data=dl3 / "skymap/{run_id}.fits.gz",
        script=scripts / "plot_skymap.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    resources:
        time=5,
    log:
        plots / "skymap/plot_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"


rule cuts_dl2_dl3:
    output:
        dl3 / "counts_after_cuts/{run_id}.h5",
    input:
        dl2=dl2 / "LST-1.Run{run_id}.dl2.h5",
        irf=dl3 / "LST-1.Run{run_id}.dl3.fits.gz",
        config=irf_config,
        script=scripts / "calc_counts_after_cuts.py",
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    resources:
        mem_mb="64G",
        time=10,
    conda:
        env
    log:
        dl3 / "counts_after_cuts/calc_{run_id}.log",
    shell:
        "python {input.script} --input-dl2 {input.dl2} --input-irf {input.irf} -c {input.config} -o {output} --log-file {log}"


def dl3_all_counts(wildcards):
    ids = RUN_IDS(wildcards)
    return [dl3 / f"counts_after_cuts/{run}.h5" for run in ids]


rule stack_cuts_dl2_dl3:
    output:
        dl3 / "counts_after_cuts/stacked.h5",
    input:
        data=dl3_all_counts,
        script=scripts / "stack_counts_after_cuts.py",
        rc=MATPLOTLIBRC,
    conda:
        env
    log:
        dl3 / "counts_after_cuts/stack.log",  # TODO use this
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule plot_cuts_dl2_dl3:
    output:
        plots / "counts_after_cuts/counts_after_cuts_{run_id}.pdf",
    input:
        data=dl3 / "counts_after_cuts/{run_id}.h5",
        script=scripts / "plot_counts_after_cuts.py",
        rc=MATPLOTLIBRC,
    conda:
        env
    log:
        plots / "counts_after_cuts/plot_counts_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"
