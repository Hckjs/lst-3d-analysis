lstchain_env = ENVS["lstchain"]
link_env = ENVS["data_selection"]
plot_env = ENVS["gammapy"]
scripts = Path(SCRIPTS["mc"])
mc = Path(OUTDIRS["mc"])
models = mc / "models"

# Need some extra dirs
mc_nodes = Path(OUTDIRS["mc_nodes"])
dl1 = Path(OUTDIRS["dl1"])
models = Path(OUTDIRS["models"])
config = CONFIGS["lstchain"]

plots = mc / "plots"


rule mc:
    input:
        link=mc / "mc-linked.txt",
        models=models_to_train,
        plots=plots / "gh_score.pdf",


localrules:
    link_mc,


checkpoint link_mc:
    output:
        dummy=mc / "mc-linked.txt",
        config=mc / "lstchain_mcpipe.json",
    input:
        script=scripts / "link-mc.py",
    params:
        production=PRODUCTION,
        declination=DECLINATION,
        mc_nodes=mc_nodes,
    conda:
        link_env
    log:
        mc / "link_mc.log",
    shell:
        "python \
        {input.script} \
        --prod {params.production} \
        --dec {params.declination} \
        --mc-nodes-link-dir {params.mc_nodes} \
        --model-config-link-path {output.config} \
        --log-file {log} \
        --verbose \
        --output-path {output.dummy}"


rule merge_gamma_mc_per_node:
    output:
        train=mc / "GammaDiffuse/dl1/{node}_train.dl1.h5",
        test=mc / "GammaDiffuse/dl1/{node}_test.dl1.h5",
    input:
        dummy=mc / "mc-linked.txt",
        script=scripts / "merge_mc_nodes.py",
    params:
        train_size=train_size,
        directory=lambda wildcards: mc_nodes / f"GammaDiffuse/{wildcards.node}",
    conda:
        lstchain_env
    log:
        mc / "GammaDiffuse/dl1/merge_gamma_mc_{node}.log",
    shell:
        "python {input.script} \
        --input-dir {params.directory} \
        --train-size {params.train_size} \
        --output-train {output.train} \
        --output-test {output.test} \
        --pattern 'dl1_*.h5' \
        --log-file {log}"


rule merge_proton_mc_per_node:
    output:
        train=mc / "Protons/dl1/{node}_train.dl1.h5",
        train=mc / "Protons/dl1/{node}_test.dl1.h5",
    input:
        dummy=mc / "mc-linked.txt",
        script=scripts / "merge_mc_nodes.py",
    params:
        train_size=0.9,
        directory=lambda wildcards: mc_nodes / f"Protons/{wildcards.node}",
    conda:
        lstchain_env
    log:
        mc / "Protons/dl1/merge_proton_mc_{node}.log",
    shell:
        "python {input.script} \
        --input-dir {params.directory} \
        --train-size {params.train_size} \
        --output-train {output.train} \
        --output-test {output.test} \
        --pattern 'dl1_*.h5' \
        --log-file {log}"


rule merge_train_or_test_of_all_nodes:
    output:
        mc / "{particle}/{particle}_{train_or_test}.dl1.h5",
    input:
        nodes=MC_NODES_DL1,
        script=scripts / "merge_mc_nodes.py",
    params:
        directory=lambda wildcards: mc / f"{wildcards.particle}/dl1",
        pattern=lambda wildcards: f"*_{wildcards.train_or_test}.dl1.h5",
        out_type=lambda wildcards: f"output-{wildcards.train_or_test}",
    conda:
        lstchain_env
    log:
        mc / "{particle}/merge_all_{particle}_{train_or_test}.log",
    shell:
        """
        python {input.script} \
        --input-dir {params.directory} \
        --pattern {params.pattern} \
        --{params.out_type} {output} \
        --log-file {log}
        """


rule train_models:
    output:
        models_to_train,
    input:
        gamma=mc / "GammaDiffuse/GammaDiffuse_train.dl1.h5",
        proton=mc / "Protons/Protons_train.dl1.h5",
        config=config,
    resources:
        mem_mb=64000,
        cpus=8,
        partition="long",
        time=1200,
    conda:
        lstchain_env
    log:
        models / "train_models.log",
    shell:
        """
        lstchain_mc_trainpipe \
        --fg {input.gamma} \
        --fp {input.proton} \
        --config {input.config} \
        --log-file {log} \
        --output-dir {models}
        """


rule calc_parameter_histograms:
    output:
        mc / "{parameter}_distribution.h5",
    input:
        gamma=mc / "GammaDiffuse/GammaDiffuse_test.dl1.h5",
        proton=mc / "Protons/Protons_test.dl1.h5",
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
