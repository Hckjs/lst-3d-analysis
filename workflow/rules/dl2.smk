env = ENVS["lstchain"]
# Having these as paths makes them easier to use
scripts = Path(SCRIPTS["dl2"])
dl2 = Path(OUTDIRS["dl2"])
models = Path(OUTDIRS["models"])
plots = Path(PLOTSDIRS["dl2"])
# TODO i dont like this path too much, but for now keep it
# At least its clear where the config comes from
config = models / "mcpipe/lstchain_config.json"


rule dl2_stage:
    input:
        runs=expand(
            dl2 / "dl2_LST-1.Run{run_id}.h5",
            run_id=RUN_IDS,
        ),


rule dl2:
    resources:
        mem_mb=64000,
        cpus=4,
    output:
        dl2 / "{potentially_test}dl2_LST-1.Run{run_id}.h5",
    input:
        data=dl1 / "{potentially_test}dl1_LST-1.Run{run_id}.h5",
        config=config,
        models=models_to_train,
    conda:
        env
    # allow wildcard to be empty to also match observed runs
    wildcard_constraints:
        potentially_test=".*",
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        """
        lstchain_dl1_to_dl2  \
            --input-file {input.data}  \
            --output-dir $(dirname {output}) \
            --path-models {input.model_dir}  \
            --config {input.config} 
        """


rule cuts_dl2_dl3:
    resources:
        mem_mb="64G",
        time=10,
    conda:
        lstchain_env
    output:
        build_dir / "dl3/counts/after_gh_theta_cut_{run_id}.h5",
    input:
        dl2=build_dir / "dl2/dl2_LST-1.Run{run_id}.h5",
        irf=build_dir / "irf/irf_Run{run_id}.fits.gz",
        config=irf_config_path,
        script="scripts/calc_counts_after_cuts.py",
    shell:
        "python {input.script} --input-dl2 {input.dl2} --input-irf {input.irf} -c {input.config} -o {output}"


rule stack_cuts_dl2_dl3:
    conda:
        lstchain_env
    output:
        build_dir / "dl3/counts/after_gh_theta_cut_{norm}_stacked.h5",
    input:
        data=expand(
            build_dir / "dl3/counts/after_gh_theta_cut_{run_id}.h5",
            run_id=RUN_IDS,
        ),
        script="scripts/stack_counts_after_cuts.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --norm {wildcards.norm}"


rule plot_cuts_dl2_dl3:
    conda:
        lstchain_env
    output:
        build_dir / "plots/counts_after_gh_theta_cut_{norm}.pdf",
    input:
        data=build_dir / "dl3/counts/after_gh_theta_cut_{norm}.h5",
        script="scripts/plot_counts_after_cuts.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


# Really using the lstchain env should be fine, but I want to not
# use the lstchain env for plots at all. TODO: One Plot env?
irf_env = ENVS["lstchain"]
plot_env = ENVS["gammapy"]
# Having these as paths makes them easier to use
scripts = Path(SCRIPTS["irfs"])
out = Path(OUTDIRS["irfs"])
dl2_test_files = Path(OUTDIRS["dl2"]) / "test"
plots = Path(PLOTSDIRS["dl2"])
config = CONFIGS["irf_tool"]


rule irf_stage:
    input:
        expand(
            plots / "{irf}/{irf}_Run{run_id}.pdf", run_id=RUN_IDS, irf=irfs_to_produce
        ),


rule irf:
    output:
        out / "irfs_Run{run_id}.fits.gz",
    input:
        gammas=dl2_test_files / "dl2_LST-1.Run{run_id}.h5",
        config=config,
    conda:
        irf_env
    resources:
        mem_mb=8000,
        time=10,
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        """
        lstchain_create_irf_files \
            -o {output} \
            -g {input.gammas} \
            --config {input.config} \
        """


rule plot_irf:
    output:
        plots / "{irf}/{irf}_Run{run_id}.pdf",
    input:
        data=irf / "irfs_Run{run_id}.fits.gz",
        script=scripts / "plot_{irf}.py",
        rc=MATPLOTLIBRC,
    conda:
        plot_env
    resources:
        mem_mb=1000,
        time=20,
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"
