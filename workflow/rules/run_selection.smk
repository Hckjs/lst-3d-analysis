# Select runs based on datacheck values

scripts = SCRIPTS["run_selection"]
env = ENVS["run_selection"]
config = CONFIGS["run_selection"]
plots = build / "plots / run_selection"
run_selection_plots = [
    PLOTDIRS["run_selection"] / "run_selection/moon-illumination.pdf",
    PLOTDIRS["run_selection"] / "run_selection/cosmics.pdf",
    PLOTDIRS["run_selection"] / "run_selection/cosmics-above.pdf",
]


rule run_selection_stage:
    input:
        run_selection_plots,


rule gather_test_nodes:
    output:
        build / "allsky-mc/test-nodes.csv",
    input:
        script=scripts / "gather-test-nodes.py",
        run_ids=run_ids,
    params:
        production=PRODUCTION,
        declination=DECLINATION,
    conda:
        env
    log:
        out=lambda wildcards, output: Path(output).with_suffix(".log"),
        err=lambda wildcards, output: Path(output).with_suffix(".err"),
    shell:
        "python {input.script} \
        --prod {params.production} \
        --dec {params.declination} \
        -o {output}"


rule gather_run_pointings:
    output:
        build / "allsky-mc/run-pointings.csv",
    input:
        runs=build / "runs.json",
        datacheck=build / "dl1-datachecks-masked.h5",
        script=scripts / "gather-run-pointings.py",
    conda:
        env
    log:
        out=lambda wildcards, output: Path(output).with_suffix(".log"),
        err=lambda wildcards, output: Path(output).with_suffix(".err"),
    shell:
        "python {input.script} \
        --runs {input.runs} \
        --runsummary {input.datacheck} \
        -o {output}"


rule plot_data_selection:
    output:
        plots / "run_selection/{name}.pdf",
    input:
        data=build / "dl1-datachecks-masked.h5",
        config=config / "dl1-selection-cuts-config.json",
        script=scripts / "plot-{name}.py",
    conda:
        env
    log:
        out=lambda wildcards, output: Path(output).with_suffix(".log"),
        err=lambda wildcards, output: Path(output).with_suffix(".err"),
    shell:
        "python {input.script} {input.data} -c {input.config} -o {output}"


rule numbers:
    output:
        plots / "runselection-01-observing-source.tex",
        plots / "runselection-02-ok-during-timeframe.tex",
        plots / "runselection-03-pedestal-charge.tex",
        plots / "runselection-04-cosmics.tex",
    input:
        data=build / "dl1-datachecks-masked.h5",
        script=scripts / "extract-numbers.py",
    conda:
        env
    log:
        out=lambda wildcards, output: Path(output).with_suffix(".log"),
        err=lambda wildcards, output: Path(output).with_suffix(".err"),
    shell:
        "python {input.script} {input.data} {plots}"
