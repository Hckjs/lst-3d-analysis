env = ENVS["lstchain"]
plot_env = ENVS["gammapy"]
irfs = Path(OUTDIRS["irfs"])
irf_config = CONFIGS["irf_tool"]
irf_scripts = Path(SCRIPTS["irfs"])
dl2 = Path(OUTDIRS["dl2"])
dl2_scripts = Path(SCRIPTS["dl2"])
models = Path(OUTDIRS["models"])
config = lstchain_config


rule dl2:
    input:
        irfs=MC_NODES_IRFs,
        runs=DL2_FILES,


rule dl1_to_dl2:
    output:
        Path("{somepath}/dl2") / "{base}.dl2.h5",
    input:
        data=Path("{somepath}/dl1") / "{base}.dl1.h5",
        config=config,
        models=models_to_train,
    conda:
        env
    resources:
        mem_mb=64000,
        cpus=4,
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


rule irf:
    output:
        irfs / "irfs_{node}.fits.gz",
    input:
        gammas=mc / "GammaDiffuse/dl2/{node}_test.dl2.h5",
        config=irf_config,
    conda:
        env
    resources:
        mem_mb=8000,
        time=10,
    log:
        irfs / "irfs_{node}.log",
    shell:
        """
        lstchain_create_irf_files \
            -o {output} \
            -g {input.gammas} \
            --config {input.config}
        """


rule plot_irf:
    output:
        "{somepath}/plots/{irf}_{base}.pdf",
    input:
        data="{somepath}/irfs_{base}.fits.gz",
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
        "{somepath}/plots/{irf}_{base}.log",
    shell:
        "MATPLOTLIBRC={input.rc} \
        python {input.script} \
        -i {input.data} \
        -o {output} \
        --log-file {log}"
