# Create DL4 datasets, plot sensitivity, significance, ...
rule dataset_3d:
    input:
        data=build_dir / "dl3/hdu-index.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script="scripts/write_datasets_manual.py",
    output:
        build_dir / "dl4/{analysis}/datasets.fits.gz",
    conda:
        gammapy_env
    shell:
        "python {input.script} -c {input.config} -o {output}"


rule calc_dl4_diagnostics:
    output:
        build_dir / "dl4/{analysis}/dl4_diagnostics.fits.gz",
    input:
        data=build_dir / "dl4/{analysis}/datasets.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script="scripts/calc_dl4_diagnostics.py",
    resources:
        mem_mb=16000,
    conda:
        gammapy_env
    shell:
        "python {input.script} -c {input.config} -o {output} --dataset-path {input.data}"


rule peek_datasets:
    output:
        build_dir / "plots/{analysis}/dataset_peek.pdf",
    input:
        data=build_dir / "dl4/{analysis}/datasets.fits.gz",
        script="scripts/plot_dataset_peek.py",
        config=config_dir / "{analysis}/analysis.yaml",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -c {input.config} -o {output} --dataset-path {input.data}"


rule plot_dl4:
    output:
        build_dir / "plots/{analysis}/{name}.pdf",
    input:
        data=build_dir / "dl4/{analysis}/{name}.fits.gz",
        script="scripts/plot_{name}.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"
