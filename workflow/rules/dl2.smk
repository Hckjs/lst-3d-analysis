lstchain_env = ENVS["lstchain"]
plot_env = ENVS["gammapy"]
irfs = Path(OUTDIRS["irfs"])
irf_scripts = Path(SCRIPTS["irfs"])
dl2 = Path(OUTDIRS["dl2"])
dl2_scripts = Path(SCRIPTS["dl2"])
models = Path(OUTDIRS["models"])
# config = lstchain_config
config = CONFIGS["lstchain"]


rule dl2:
    input:
        runs=DL2_FILES,


rule dl1_to_dl2:
    output:
        Path("{somepath}/dl2") / "{base}.dl2.h5",
    input:
        data=Path("{somepath}/dl1") / "{base}.dl1.h5",
        config=config,
        models=models_to_train,
    conda:
        lstchain_env
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
        irfs / "{analysis}/irfs_{node}.fits.gz",
    input:
        gammas=mc / "GammaDiffuse/dl2/{node}_test.dl2.h5",
        config=config_dir / "{analysis}/irf_tool_config.json",
        edisp_script=irf_scripts / "fix_edisp.py",
    conda:
        lstchain_env
    resources:
        mem_mb=8000,
        time=10,
    log:
        irfs / "{analysis}/irfs_{node}.log",
    shell:
        """
        lstchain_create_irf_files \
            -o {output} \
            -g {input.gammas} \
            --config {input.config}

        python {input.edisp_script} {output}
        """


rule plot_irf:
    output:
        "{somepath}/{analysis}/plots/{irf}/{irf}_{base}.pdf",
    input:
        data="{somepath}/{analysis}/irfs_{base}.fits.gz",
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
        "{somepath}/{analysis}/plots/{irf}/{irf}_{base}.log",
    shell:
        "MATPLOTLIBRC={input.rc} \
        python {input.script} \
        -i {input.data} \
        -o {output} \
        --log-file {log}"
