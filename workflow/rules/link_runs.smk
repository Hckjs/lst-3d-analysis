# "Preselect" runs and link (no data quality checks)
# Choose runs observing the source and find the observations of these runs
# TODO: Merge with run_selection?

env = ENVS["link_runs"]
config = CONFIGS["run_selection"]
# Having these as paths makes them easier to use
scripts = Path(SCRIPTS["link_runs"])
out = Path(OUTDIRS["link_runs"])
plots = Path(PLOTSDIRS["link_runs"])


rule link_runs_stage:
    input:
        OUTDIRS["link_runs"] / "all-linked.txt",


localrules:
    runlist,
    select_datasets,
    merge_datachecks,
    run_ids,
    data_check,


# Can not automatize this easily as it requires a password
# Could read that from somewhere, but I think this is fine for now
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
    log:  # there is not really one good name here based on output
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
