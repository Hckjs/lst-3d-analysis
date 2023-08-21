env = ENVS["lstchain"]
# Having these as paths makes them easier to use
scripts = Path(SCRIPTS["dl2"])
dl2 = Path(OUTDIRS["dl2"])
models = Path(OUTDIRS["models"])
plots = Path(PLOTSDIRS["dl2"])
# TODO i dont like this path too much, but for now keep it
# At least its clear where the config comes from
config = models / "mcpipe/lstchain_config.json"


# rule dl2_stage:
#    input: ???


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
