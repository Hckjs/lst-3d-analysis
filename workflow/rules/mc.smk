env = ENVS["lstchain"]
# Having these as paths makes them easier to use
scripts = Path(SCRIPTS["mc"])
# Only the intermediate files end up in out, models are extra and final MCs are in dl1
mc_nodes = Path(OUTDIRS["mc_nodes"])
mc_merged = Path(OUTDIRS["mc_merged"])
dl1 = Path(OUTDIRS["dl1"])
models = Path(OUTDIRS["models"])
plots = Path(PLOTSDIRS["link_runs"])
train_size = (0.4,)  # TODO put in config?


rule mc_stage_target:
    input:
        m=models_to_train,


def all_gamma_nodes(wildcards):
    linked = checkpoints.run_ids.get(**wildcards).output
    all_gamma_nodes = [
        x.name for x in (mc_nodes / "GammaDiffuse").glob("*") if x.is_dir()
    ]
    return all_gamma_nodes


rule merge_gamma_mc_per_node:
    output:
        train=mc_merged / "GammaDiffuse/{node}_train.dl1.h5",
        test=mc_merged / "GammaDiffuse/{node}_test.dl1.h5",
    input:
        all_gamma_nodes,
    params:
        train_size=train_size,
        directory=lambda node: mc_nodes / f"GammaDiffuse/{node}",
    conda:
        env
    log:
        merged / "merge_gamma_mc_{node}.log",
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
        test_files=expand(
            merged / "GammaDiffuse/{node}_test.dl1.h5", node=all_gamma_nodes
        ),
    output:
        runs=expand(
            dl1 / "test/dl1_LST-1.Run{run_id}.h5",
            run_id=RUN_IDS,
        ),
    # Just to avoid complex things in the command below
    params:
        dl1_test_dir=dl1 / "test",
        mc_nodes_dir=merged / "GammaDiffuse",  # these are the merged train/test ones
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
        dl1 / "{train_or_test}/{particle}_train.dl1.h5",
    input:
        files=expand(
            merged / "{{particle}}/{node}_{{train_or_test}}.dl1.h5",
            node=all_gamma_nodes,
        ),
    params:
        directory=lambda wildcards: merged / f"{wildcards.particle}",
        pattern=lambda wildcards: f"*_{wildcards.train_or_test}.dl1.h5",
        out_type=lambda wildcards: f"output-{wildcards.train_or_test}",
    conda:
        env
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
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
