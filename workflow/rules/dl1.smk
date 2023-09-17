env = ENVS["data_selection"]
config = CONFIGS["data_selection"]
scripts = Path(SCRIPTS["data_selection"])
out = Path(OUTDIRS["data_selection"])
dl1_link_location = Path(OUTDIRS["dl1"])
plots = out / "plots"


run_selection_plots = [
    plots / f"{name}.pdf"
    for name in ["moon-illumination", "cosmics", "cosmics-above", "run-pointings"]
]


rule dl1:
    input:
        out / "runs-linked.txt",
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
        out / "merge_datacheck.log",
    shell:
        "python {input.script} {input.data} {output.output} --log-file {log}"


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
        out / "datacheck.log",
    shell:
        "python \
        {input.script} \
        {input.runlist} \
        {input.datachecks} \
        --config {input.config} \
        --output-runlist {output.runlist} \
        --output-datachecks {output.datachecks} \
        --output-config {output.config} \
        --log-file {log}"


checkpoint run_ids:
    output:
        runlist=out / "runs.json",
    input:
        data=out / "runlist-checked.csv",
        config=config,
        script=scripts / "create-night-run-list.py",
    conda:
        env
    log:
        out / "check_runlist.log",
    shell:
        "python \
        {input.script} \
        {input.data} \
        {output} \
        -c {input.config} \
        --log-file {log}"


checkpoint link_runs:
    output:
        out / "runs-linked.txt",
    input:
        runs=out / "runs.json",
        datacheck=out / "dl1-datachecks-masked.h5",
        script=scripts / "link-runs.py",
    params:
        dl1=dl1_link_location,
    conda:
        env
    log:
        out / "link_runs.log",
    shell:
        "python \
        {input.script} \
        --runs {input.runs} \
        --dl1-link-dir {params.dl1} \
        --log-file {log} \
        --output-path {output}"


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
        plots / "{name}.log",
    shell:
        "python \
        {input.script} \
        {input.data} \
        -c {input.config} \
        -o {output} \
        --log-file {log} "


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
        out / "run_pointings.log",
    shell:
        "python {input.script} \
        --runs {input.runs} \
        --runsummary {input.datacheck} \
        --output {output} \
        --log-file {log} "


rule plot_run_pointings:
    output:
        plots / "run-pointings.pdf",
    input:
        pointings=out / "run-pointings.csv",
        script=scripts / "plot-run-pointings.py",
    conda:
        env
    log:
        plots / "run_pointings.log",
    shell:
        "python {input.script} \
        --input {input.pointings} \
        --output {output} \
        --log-file {log} "
