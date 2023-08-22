env = ENVS["data_selection"]
config = CONFIGS["data_selection"]
scripts = Path(SCRIPTS["data_selection"])
out = Path(OUTDIRS["data_selection"])
plots = out / "plots"


run_selection_plots = [
    plots / f"{name}.pdf"
    for name in ["moon-illumination", "cosmics", "cosmics-above", "run-pointings"]
]


rule dl1:
    input:
        out / "all-linked.txt",
        run_selection_plots,


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
        echo 'curl --user <username>:<password> https://lst1.iac.es/datacheck/lstosa/LST_source_catalog.html -o {output}'
        echo 'You might need to create the output directory first.'
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
        out=out / "select_datasets.log",
        err=out / "select_datasets.err",
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
        out=out / "merge_datacheck.log",
        err=out / "merge_datacheck.err",
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
        out=out / "datacheck.log",
        err=out / "datacheck.err",
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
        out=out / "check_runlist.log",
        err=out / "check_runlist.err",
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
        out=out / "link_paths.log",
        err=out / "link_paths.err",
    shell:
        "python {input.script} \
        --runs {input.runs} \
        --prod {params.production} \
        --dec {params.declination} \
        --runsummary {input.datacheck} \
        -o {output}"


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
        out=out / "run_pointings.log",
        err=out / "run_pointings.err",
    shell:
        "python {input.script} \
        --runs {input.runs} \
        --runsummary {input.datacheck} \
        -o {output}"


rule plot_data_selection:
    output:
        plots / "{name}.pdf",
    input:
        data=out / "dl1-datachecks-masked.h5",
        config=out / "dl1-selection-cuts-config.json",
        script=scripts / "plot-{name}.py",
    conda:
        env
    log:
        out=plots / "{name}.log",
        err=plots / "{name}.err",
    shell:
        "python {input.script} {input.data} -c {input.config} -o {output}"
