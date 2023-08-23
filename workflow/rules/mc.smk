env = ENVS["lstchain"]
scripts = Path(SCRIPTS["mc"])
mc = Path(OUTDIRS["mc"])

# Need some extra dirs
mc_nodes = Path(OUTDIRS["mc_nodes"])
dl1 = Path(OUTDIRS["dl1"])
models = Path(OUTDIRS["models"])

plots = mc / "plots"
train_size = (0.4,)  # TODO put in config?


rule mc:
    input:
        m=models_to_train,


checkpoint link_mc:
    output:
        out / "mc-linked.txt",
    input:
        script=scripts / "link-mc.py",
    params:
        production=PRODUCTION,
        declination=DECLINATION,
        mc_nodes=mc_link_location,
        models=models_link_location,
    conda:
        env
    log:
        out / "link_mc.log",
    shell:
        "python \
        {input.script} \
        --prod {params.production} \
        --dec {params.declination} \
        --mc-nodes-link-dir {params.mc_nodes} \
        --models-link-dir {params.models} \
        --log-file {log} \
        --output-path {output}"


rule merge_gamma_mc_per_node:
    output:
        train=mc / "GammaDiffuse/{node}_train.dl1.h5",
        test=mc / "GammaDiffuse/{node}_test.dl1.h5",
    input:
        # put the all linked here and in the merge train test to amek sure
        # this is valid!
        x.name for x in (mc_nodes / "GammaDiffuse").glob("*") if x.is_dir(),
    params:
        train_size=train_size,
        directory=lambda node: mc_nodes / f"GammaDiffuse/{node}",
    conda:
        env
    log:
        mc / "merge_gamma_mc_{node}.log",
    shell:
        """
        python scripts/merge_mc_nodes.py \
        --input-dir {params.directory} \
        --train-size {params.train_size} \
        --output-train {output.train} \
        --output-test {output.test} \
        --log-file {log}
        """


# Cant do this in link script because the merge into train and test did not happen yet at that stage
rule link_test_nodes_to_dl1_runs:
    input:
        runs=expand(
            dl1 / "dl1_LST-1.Run{run_id}.h5",
            run_id=RUN_IDS,
        ),
        test_files=expand(mc / "GammaDiffuse/{node}_test.dl1.h5", node=all_gamma_nodes),
    output:
        runs=expand(
            dl1 / "test/dl1_LST-1.Run{run_id}.h5",
            run_id=RUN_IDS,
        ),
    # Just to avoid complex things in the command below
    params:
        dl1_test_dir=dl1 / "test",
        mc_nodes_dir=mc / "GammaDiffuse",  # these are the merged train/test ones
    conda:
        env
    log:
        models / "link_test_nodes_to_dl1_runs.log",
    shell:
        """
        python link-dl1.py \
        --dl1-dir {dl1} \
        --test-link-dir {params.dl1_test_dir} \
        --mc-nodes-dir {params.mc_nodes_dir} \
        --log-file {log}
        """


rule merge_train_or_test_of_all_nodes:
    output:
        dl1 / "{train_or_test}/{particle}_{train_or_test}.dl1.h5",
    input:
        files=expand(
            mc / "{{particle}}/{node}_{{train_or_test}}.dl1.h5",
            node=all_gamma_nodes,
        ),
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
        --{params.out_type} {output}
        """


rule train_models:
    output:
        models_to_train,
    input:
        gamma=dl1 / "train/GammaDiffuse_train.dl1.h5",
        proton=dl1 / "train/proton_diffuse_merged.dl1.h5",
        config=models / "mcpipe/lstchain_config.json",
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
        --output-dir {models} \
        > {log} 2>&1
        """
