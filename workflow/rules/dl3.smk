lstchain_env = ENVS["lstchain"]
bkg_env = ENVS["background"]
gammapy_env = ENVS["gammapy"]

dl2 = Path(OUTDIRS["dl2"])
dl3 = Path(OUTDIRS["dl3"])
irfs = Path(OUTDIRS["irfs"])

scripts = Path(SCRIPTS["dl3"])
data_selection_config = CONFIGS["data_selection"]


rule dl3:
    input:
        index=[dl3 / f"{analysis}/hdu-index.fits.gz" for analysis in analyses],
        bkg=[dl3 / f"{analysis}/bkg-exists" for analysis in analyses],
        plots=DL3_PLOTS,


rule dl2_to_dl3:
    output:
        run=dl3 / "{analysis}/LST-1.Run{run_id}.dl3.fits.gz",
    input:
        data=dl2 / "LST-1.Run{run_id}.dl2.h5",
        irfs=MC_NODES_IRFs,  # changes baes on analysis
        config=config_dir / "{analysis}/irf_tool_config.json",
    params:
        irf_pattern="irfs_*.fits.gz",
        out=lambda wc: dl3 / wc.get("analysis"),
        in_irfs=lambda wc: irfs / wc.get("analysis"),
    conda:
        lstchain_env
    resources:
        mem_mb=12000,
        time=30,
    log:
        dl3 / "{analysis}/create_dl3_{run_id}.log",
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
        dl3 / "{analysis}/bkg_cached_maps.pkl",
    input:
        runs=DL3_FILES,
        config=config_dir / "{analysis}/bkgmodel.yml",
        script=scripts / "precompute_background_maps.py",
        bkg_exclusion_regions=config_dir / "{analysis}/bkg_exclusion",
    conda:
        bkg_env
    resources:
        partition="long",
        time=360,
    log:
        dl3 / "{analysis}/calc_count_maps.log",
    shell:
        """python {input.script} \
        --input-runs {input.runs} \
        --exclusion {input.bkg_exclusion_regions} \
        --output {output} \
        --config {input.config} \
        --log-file {log} \
        --overwrite
        """


rule calc_background:
    output:
        dummy=dl3 / "{analysis}/bkg-exists",
    input:
        runs=DL3_FILES,
        config=config_dir / "{analysis}/bkgmodel.yml",
        script=scripts / "calc_background.py",
        cached_maps=dl3 / "{analysis}/bkg_cached_maps.pkl",
        bkg_exclusion_regions=config_dir / "{analysis}/bkg_exclusion",
    params:
        bkg_dir=lambda wc: dl3 / wc.get("analysis"),
    conda:
        bkg_env
    resources:
        partition="short",
    log:
        dl3 / "{analysis}/calc_bkg.log",
    shell:
        """python {input.script} \
        --input-runs {input.runs} \
        --output-dir {params.bkg_dir} \
        --exclusion {input.bkg_exclusion_regions} \
        --dummy-output {output.dummy} \
        --cached-maps {input.cached_maps} \
        --config {input.config} \
        --log-file {log} \
        --verbose \
        --overwrite
        """


def DL3_INDEX_FILELIST(wildcards):
    files = DL3_FILES(wildcards)
    return ("--file-list " + " --file-list ".join([f.name for f in files]),)


rule dl3_hdu_index:
    output:
        dl3 / "{analysis}/hdu-index.fits.gz",
    input:
        runs=DL3_FILES,
        index_script=scripts / "create_hdu_index.py",
        link_script=scripts / "link_bkg.py",
        dummy=dl3 / "{analysis}/bkg-exists",
    params:
        outdir=lambda wc: dl3 / wc.get("analysis"),
        bkg=BKG_FILES,
        filelist=DL3_INDEX_FILELIST,
    conda:
        lstchain_env
    log:
        dl3 / "{analysis}/hdu_index.log",
    shell:
        """
        python {input.index_script}  \
            --input-dl3-dir {params.outdir}  \
            --output-index-path {params.outdir}  \
            {params.filelist} \
            --overwrite \
            --log-file {log}

        python {input.link_script} \
        --hdu-index-path {output} \
        --bkg-files {params.bkg} \
        """


rule plot_dl3_rates:
    output:
        dl3 / "{analysis}/plots/rates.pdf",
    input:
        index=dl3 / "{analysis}/hdu-index.fits.gz",
        script=scripts / "plot_rates.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        plots / "{analysis}/rates.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.index} -o {output} --log-file {log}"


rule plot_dl3_irf_comparison:
    output:
        dl3 / "{analysis}/plots/irfs.pdf",
    input:
        index=dl3 / "{analysis}/hdu-index.fits.gz",
        script=scripts / "plot_irf_comparison.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        dl3 / "{analysis}/plots/irfs.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.index} -o {output} --log-file {log}"


rule calc_theta2_per_obs:
    output:
        dl3 / "{analysis}/theta2/{run_id}.fits.gz",
    input:
        data=dl3 / "{analysis}/LST-1.Run{run_id}.dl3.fits.gz",
        script=scripts / "calc_theta2_per_obs.py",
        config=data_selection_config,  # this seems unnecessary
        index=dl3 / "{analysis}/hdu-index.fits.gz",
    params:
        outdir=lambda wc: dl3 / wc.get("analysis"),
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    resources:
        mem_mb=16000,
    conda:
        gammapy_env
    log:
        dl3 / "{analysis}/theta2/calc_{run_id}.log",
    shell:
        "python {input.script} -i {params.outdir} -o {output} --obs-id {wildcards.run_id} --config {input.config} --log-file {log}"


def dl3_all_theta_tables(wildcards):
    ids = RUN_IDS(wildcards)
    return [dl3 / f"{wildcards.analysis}/theta2/{run}.fits.gz" for run in ids]


rule stack_theta2:
    output:
        dl3 / "{analysis}/theta2/stacked.fits.gz",
    input:
        runs=dl3_all_theta_tables,
        script=scripts / "stack_theta2.py",
    conda:
        gammapy_env
    log:
        dl3 / "{analysis}/theta2/theta2_stacked.log",
    shell:
        "python {input.script} -o {output} --input-files {input.runs} --log-file {log}"


rule plot_theta:
    output:
        dl3 / "{analysis}/plots/theta2/theta2_{run_id}.pdf",
    input:
        data=dl3 / "{analysis}/theta2/{run_id}.fits.gz",
        script=scripts / "plot_theta2.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        dl3 / "{analysis}/plots/theta2/plot_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"


rule plot_background:
    output:
        dl3 / "{analysis}/plots/bkg/bkg_{run_id}.pdf",
    input:
        data=dl3 / "{analysis}/bkg-exists",
        script=scripts / "plot_bkg.py",
        rc=MATPLOTLIBRC,
    params:
        data=lambda wildcards: dl3
        / f"{wildcards.analysis}/bkg_{wildcards.run_id}.fits.gz",
    conda:
        gammapy_env
    log:
        dl3 / "{analysis}/plots/bkg/bkg_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {params.data} -o {output}"


rule calc_skymap:
    output:
        dl3 / "{analysis}/skymap/{run_id}.fits.gz",
    input:
        data=dl3 / "{analysis}/LST-1.Run{run_id}.dl3.fits.gz",
        script=scripts / "calc_skymap_gammas.py",
        config=config_dir / "{analysis}/irf_tool_config.json",
        index=dl3 / "{analysis}/hdu-index.fits.gz",
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    params:
        outdir=lambda wc: dl3 / wc.get("analysis"),
        n_bins=50,
    conda:
        gammapy_env
    log:
        dl3 / "{analysis}/skymap/calc_{run_id}.log",
    shell:
        "python {input.script} -i {params.outdir} -o {output} --obs-id {wildcards.run_id} --config {input.config} --log-file {log} --n-bins {params.n_bins}"


def dl3_all_skymaps(wildcards):
    ids = RUN_IDS(wildcards)
    return [dl3 / f"{wildcards.analysis}/skymap/{run}.fits.gz" for run in ids]


rule stack_skymaps:
    output:
        dl3 / "{analysis}/skymap/stacked.fits.gz",
    input:
        data=dl3_all_skymaps,
        script=scripts / "stack_skymap.py",
    conda:
        gammapy_env
    log:
        dl3 / "{analysis}/skymap/stack.log",
    shell:
        "python {input.script} -i {input.data} -o {output} --log-file {log}"


rule plot_skymap:
    output:
        dl3 / "{analysis}/plots/skymap/skymap_{run_id}.pdf",
    input:
        data=dl3 / "{analysis}/skymap/{run_id}.fits.gz",
        script=scripts / "plot_skymap.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    resources:
        time=5,
    log:
        dl3 / "{analysis}/plots/skymap/plot_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"


rule cuts_dl2_dl3:
    output:
        dl3 / "{analysis}/counts_after_cuts/{run_id}.h5",
    input:
        dl2=dl2 / "LST-1.Run{run_id}.dl2.h5",
        irf=dl3 / "{analysis}/LST-1.Run{run_id}.dl3.fits.gz",
        config=config_dir / "{analysis}/irf_tool_config.json",
        script=scripts / "calc_counts_after_cuts.py",
    wildcard_constraints:
        run_id="\d+",  # dont match on "stacked".
    resources:
        mem_mb="64G",
    conda:
        lstchain_env
    log:
        dl3 / "{analysis}/counts_after_cuts/calc_{run_id}.log",
    shell:
        "python {input.script} --input-dl2 {input.dl2} --input-irf {input.irf} -c {input.config} -o {output} --log-file {log}"


def dl3_all_counts(wildcards):
    ids = RUN_IDS(wildcards)
    return [dl3 / f"{wildcards.analysis}/counts_after_cuts/{run}.h5" for run in ids]


rule stack_cuts_dl2_dl3:
    output:
        dl3 / "{analysis}/counts_after_cuts/stacked.h5",
    input:
        data=dl3_all_counts,
        script=scripts / "stack_counts_after_cuts.py",
        rc=MATPLOTLIBRC,
    conda:
        lstchain_env
    log:
        dl3 / "{analysis}/counts_after_cuts/stack.log",  # TODO use this
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule plot_cuts_dl2_dl3:
    output:
        dl3 / "{analysis}/plots/counts_after_cuts/counts_after_cuts_{run_id}.pdf",
    input:
        data=dl3 / "{analysis}/counts_after_cuts/{run_id}.h5",
        script=scripts / "plot_counts_after_cuts.py",
        rc=MATPLOTLIBRC,
    conda:
        lstchain_env
    log:
        dl3 / "{analysis}/dl3/counts_after_cuts/plot_counts_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"


rule plot_run_irf:
    output:
        "{somepath}/{analysis}/plots/{irf}/{irf}_{run_id}.pdf",
    input:
        data=dl3 / "{analysis}/LST-1.Run{run_id}.dl3.fits.gz",
        script=irf_scripts / "plot_irf_{irf}.py",
        rc=MATPLOTLIBRC,
    conda:
        plot_env
    resources:
        mem_mb=1000,
        time=20,
    wildcard_constraints:
        irf="|".join(irfs_to_produce),
    log:
        "{somepath}/{analysis}/plots/{irf}/{irf}_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} \
        python {input.script} \
        -i {input.data} \
        -o {output} \
        --log-file {log}"
