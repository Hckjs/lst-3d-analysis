env = ENVS["data_selection"]
config = CONFIGS["data_selection"]
scripts = Path(SCRIPTS["data_selection"])
out = Path(OUTDIRS["data_selection"])
plots = Path(PLOTSDIRS["data_selection"])


rule link_runs_stage:
    input:
        out / "all-linked.txt",


localrules:
    runlist,
    select_datasets,
    merge_datachecks,
    run_ids,
    data_check,


rule runlist:
    output:
        out / "runlist.html",
    shell:
        """
        echo 'Provide the file {output}. The command is:'
        echo 'curl --user <username>:<password> https://lst1.iac.es/datacheck/lstosa/LST_source_catalog.html -o runlist.html'
        """


rule select_datasets:
    output:
        out / "runlist.csv",
    input:
        data=out / "runlist.html",
        config=config,
        script=scripts / "select-data.py",
    conda:
        env
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        "python {input.script} {input.data} {output} -c {input.config}"


rule merge_datachecks:
    output:
        output=out / "dl1-datachecks-merged.h5",
    input:
        data=out / "runlist.csv",
        script=scripts / "merge-datachecks.py",
    conda:
        env
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        "python {input.script} {input.data} {output.output}"


rule data_check:
    output:
        runlist=out / "runlist-checked.csv",
        datachecks=out / "dl1-datachecks-masked.h5",
        config=out / "dl1-selection-cuts-config.json",
    input:
        runlist=out / "runlist.csv",
        datachecks=out / "dl1-datachecks-merged.h5",
        config=config,
        script=scripts / "data-check.py",
    conda:
        env
    log:
        out=out / "datackeck.log",
        err=out / "datackeck.err",
    shell:
        "\
        python \
            {input.script} \
            {input.runlist} \
            {input.datachecks} \
            --config {input.config} \
            --output-runlist {output.runlist} \
            --output-datachecks {output.datachecks} \
            --output-config {output.config} \
        "


rule run_ids:
    output:
        out / "runs.json",
    input:
        data=out / "runlist-checked.csv",
        config=config,
        script=scripts / "create-night-run-list.py",
    conda:
        env
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        "python {input.script} {input.data} {output} -c {input.config}"


checkpoint link_paths:
    output:
        out / "all-linked.txt",
    input:
        runs=out / "runs.json",
        datacheck=out / "dl1-datachecks-masked.h5",
        script=scripts / "link-paths.py",
    params:
        production=PRODUCTION,
        declination=DECLINATION,
    conda:
        env
    log:
        out=lambda wildcards, output: output.with_suffix(".log"),
        err=lambda wildcards, output: output.with_suffix(".err"),
    shell:
        "python {input.script} \
        --runs {input.runs} \
        --prod {params.production} \
        --dec {params.declination} \
        --runsummary {input.datacheck} \
        -o {output}"


# Select runs based on datacheck values

run_selection_plots = [
    out / "plots/run_selection/moon-illumination.pdf",
    out / "plots/run_selection/cosmics.pdf",
    out / "plots/run_selection/cosmics-above.pdf",
]


rule run_selection_stage:
    input:
        run_selection_plots,


rule gather_run_pointings:
    output:
        out / "run-pointings.csv",
    input:
        runs=out / "runs.json",
        datacheck=out / "dl1-datachecks-masked.h5",
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
        out / "plots/{name}.pdf",
    input:
        data=out / "dl1-datachecks-masked.h5",
        config=config / "dl1-selection-cuts-config.json",
        script=scripts / "plot-{name}.py",
    conda:
        env
    log:
        out=lambda wildcards, output: Path(output).with_suffix(".log"),
        err=lambda wildcards, output: Path(output).with_suffix(".err"),
    shell:
        "python {input.script} {input.data} -c {input.config} -o {output}"
