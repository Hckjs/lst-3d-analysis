env = ENVS["gammapy"]

dl3 = Path(OUTDIRS["dl3"])
dl4 = Path(OUTDIRS["dl4"])
scripts = Path(SCRIPTS["dl4"])

plot_types = ["obs_plots", "dataset_peek"]  # , "dl4_diagnostics"]


rule dl4:
    input:
        [dl4 / "{analysis}/{plot}.pdf" for analysis in analyses for plot in plot_types],


rule observation_plots:
    output:
        plots / "{analysis}/obs_plots.pdf",
    input:
        dl3 / "hdu-index.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script=scripts / "obs_plots.py",
        rc=MATPLOTLIBRC,
    resources:
        mem_mb=16000,
    conda:
        gammapy_env
    log:
        plots / "{analysis}/obs_plots.log",
    shell:
        """
        MATPLOTLIBRC={input.rc} \
        python {input.script} \
            -c {input.config} \
            -o {output} \
            --log-file {log}
        """


rule dataset_3d:
    output:
        build_dir / "dl4/{analysis}/datasets.fits.gz",
    input:
        data=build_dir / "dl3/hdu-index.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script="scripts/write_datasets_3d.py",
    conda:
        env
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
