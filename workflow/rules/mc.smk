env = ENVS["lstchain"]
link_env = ENVS["data_selection"]
scripts = Path(SCRIPTS["mc"])
mc = Path(OUTDIRS["mc"])

# Need some extra dirs
mc_nodes = Path(OUTDIRS["mc_nodes"])
dl1 = Path(OUTDIRS["dl1"])
models = Path(OUTDIRS["models"])
model_config = models / "config"

plots = mc / "plots"

# TODO Configurable
train_size = 0.4


rule mc:
    input:
        link=mc / "mc-linked.txt",
        models=models_to_train,


localrules:
    link_mc,


checkpoint link_mc:
    output:
        dummy=mc / "mc-linked.txt",
        config=model_config / "lstchain_config.json",
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
        train=mc / "GammaDiffuse/{node}_train.dl1.h5",
        test=mc / "GammaDiffuse/{node}_test.dl1.h5",
    input:
        mc / "mc-linked.txt",
    params:
        train_size=train_size,
        directory=lambda wildcards: mc_nodes / f"GammaDiffuse/{wildcards.node}",
    conda:
        env
    log:
        mc / "merge_gamma_mc_{node}.log",
    shell:
        "python scripts/merge_mc_nodes.py \
        --input-dir {params.directory} \
        --train-size {params.train_size} \
        --output-train {output.train} \
        --output-test {output.test} \
        --log-file {log}"


rule merge_proton_mc_per_node:
    output:
        train=mc / "Proton/{node}_train.dl1.h5",
    input:
        mc / "mc-linked.txt",
    params:
        train_size=1.0,
        directory=lambda wildcards: mc_nodes / f"Proton/{wildcards.node}",
    conda:
        env
    log:
        mc / "merge_proton_mc_{node}.log",
    shell:
        "python scripts/merge_mc_nodes.py \
        --input-dir {params.directory} \
        --train-size {params.train_size} \
        --output-train {output.train} \
        --log-file {log}"


rule merge_train_or_test_of_all_nodes:
    output:
        dl1 / "{train_or_test}/{particle}_{train_or_test}.dl1.h5",
    input:
        MC_NODES,
    params:
        directory=lambda wildcards: mc / f"{wildcards.particle}",
        pattern=lambda wildcards: f"*_{wildcards.train_or_test}.dl1.h5",
        out_type=lambda wildcards: f"output-{wildcards.train_or_test}",
    conda:
        env
    log:
        mc / "merge_all_{particle}_{train_or_test}.log",
    shell:
        """
        python scripts/merge_mc_nodes.py \
        --input-dir {params.directory} \
        --pattern {params.pattern} \
        --{params.out_type} {output} \
        --log-file {log}
        """


# TODO Any chance to get logging in here?
rule train_models:
    output:
        models_to_train,
    input:
        gamma=dl1 / "train/GammaDiffuse_train.dl1.h5",
        proton=dl1 / "train/Proton_train.dl1.h5",
        config=model_config / "lstchain_config.json",
    resources:
        mem_mb=64000,
        cpus=8,
        partition="long",
        time=1200,
    conda:
        env
    log:
        models / "train_models.log",
    shell:
        """
        lstchain_mc_trainpipe \
        --fg {input.gamma} \
        --fp {input.proton} \
        --config {input.config} \
        --output-dir {models}
        """
