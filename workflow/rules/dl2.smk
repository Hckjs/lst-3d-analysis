lstchain_env = ENVS["lstchain"]
plot_env = ENVS["gammapy"]
scripts = Path(SCRIPTS["dl2"])
irfs = Path(OUTDIRS["irfs"])
irf_scripts = Path(SCRIPTS["irfs"])
mc = Path(OUTDIRS["mc"])
dl1 = Path(OUTDIRS["dl1"])
dl2 = Path(OUTDIRS["dl2"])
dl2_scripts = Path(SCRIPTS["dl2"])
models = Path(OUTDIRS["models"])
config = CONFIGS["lstchain"]


rule dl2:
    input:
        runs=DL2_FILES,


rule dl1_to_dl2_mc:
    output:
        mc / "{particle}/dl2/{base}.dl2.h5",
    input:
        config=config,
        models=models_to_train,
        mc_exists=mc / "mc-linked.txt",
        data=mc / "{particle}/dl1/{base}.dl1.h5",
    conda:
        lstchain_env
    resources:
        mem_mb=64000,
        cpus=4,
    log:
        mc / "{particle}/dl2/dl1_to_dl2_{base}.log",
    shell:
        """
        lstchain_dl1_to_dl2  \
            --input-file {params.data}  \
            --output-dir $(dirname {output}) \
            --path-models {models}  \
            --config {input.config}
        """


rule dl1_to_dl2_data:
    output:
        dl2 / "{base}.dl2.h5",
    input:
        config=config,
        models=models_to_train,
        dl1_exists=dl1 / "runs-linked.txt",
    params:
        data=lambda wc: dl1 / f"{wc.get('base')}.dl1.h5",
    conda:
        lstchain_env
    resources:
        mem_mb=64000,
        cpus=4,
    log:
        dl2 / "dl1_to_dl2_{base}.log",
    shell:
        """
        lstchain_dl1_to_dl2  \
            --input-file {params.data}  \
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


rule calc_parameter_histograms:
    output:
        mc / "{parameter}_distribution.h5",
    input:
        gamma=mc / "GammaDiffuse/dl2/GammaDiffuse_test.dl2.h5",
        proton=mc / "Protons/dl2/Protons_test.dl2.h5",
        script=scripts / "calc_parameter_distribution.py",
    conda:
        lstchain_env
    log:
        models / "calc_{parameter}_distribution.log",
    shell:
        """
        python {input.script} \
        --input-gamma {input.gamma} \
        --input-proton {input.proton} \
        --parameter {wildcards.parameter} \
        --output {output} \
        --log-file {log}
        """


rule plot_parameter_histograms:
    output:
        plots / "{parameter}.pdf",
    input:
        data=mc / "{parameter}_distribution.h5",
        script=scripts / "plot_parameter_distribution.py",
    conda:
        lstchain_env
    log:
        models / "plot_{parameter}_distribution.log",
    shell:
        """
        python {input.script} \
        --input-path {input.data} \
        --output-path {output} \
        --log-file {log}
        """
