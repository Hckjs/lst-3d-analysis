env = ENVS["lstchain"]
irfs = Path(OUTDIRS["irfs"])
irf_config = CONFIGS["irf_tool"]
irf_scripts = Path(SCRIPTS["irfs"])
dl2 = Path(OUTDIRS["dl2"])
dl2_scripts = Path(SCRIPTS["dl2"])
models = Path(OUTDIRS["models"])
config = lstchain_config


rule dl2_stage:
    input:
        runs=expand(
            dl2 / "LST-1.Run{run_id}.dl2.h5",
            run_id=RUN_IDS,
        ),
        irfs=MC_NODES_IRFs,


rule dl1_to_dl2:
    resources:
        mem_mb=64000,
        cpus=4,
    output:
        "{somepath}/dl2" / "{base}.dl2.h5",
    input:
        data="{somepath}/dl1" / "{base}.dl1.h5",
        config=config,
        models=models_to_train,
    conda:
        env
    log:
        "{somepath}/dl2/dl1_to_dl2_{base}.log",
    shell:
        """
        lstchain_dl1_to_dl2  \
            --input-file {input.data}  \
            --output-dir $(dirname {output}) \
            --path-models {models}  \
            --config {input.config}
        """


# TODO Logging
rule irf:
    output:
        irfs / "irf_{node}.fits.gz",
    input:
        gammas=mc / "GammaDiffuse/dl1/{node}_test.dl1.h5",
        config=irf_config,
    conda:
        env
    resources:
        mem_mb=8000,
        time=10,
    log:
        irfs / "irf_{node}.log",
    shell:
        """
        lstchain_create_irf_files \
            -o {output} \
            -g {input.gammas} \
            --config {input.config}
        """


rule plot_irf:
    output:
        "{somedir}/plots/{irf}_{base}.pdf",
    input:
        data="{somedir}/irfs_{base}.fits.gz",
        script=irf_scripts / "plot_{irf}.py",
        rc=MATPLOTLIBRC,
    conda:
        plot_env
    resources:
        mem_mb=1000,
        time=20,
    log:
        "{somedir}/plots/{irf}.log",
    shell:
        "MATPLOTLIBRC={input.rc} \
        python {input.script} \
        -i {input.data} \
        -o {output} \
        --log-file {log}"
