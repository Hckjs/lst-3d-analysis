gammapy_env = ENVS["gammapy"]
dl3 = Path(OUTDIRS["dl3"])
dl4 = Path(OUTDIRS["dl4"])
plots = dl4 / "plots"
scripts = Path(SCRIPTS["dl4"])

dl4_plot_types = ["dataset_peek"]  # , "dl4_diagnostics"]


rule dl4:
    input:
        [
            plots / f"{analysis}/{plot}.pdf"
            for analysis in analyses
            for plot in dl4_plot_types
        ],


rule create_fov_bkg_exclusion:
    output:
        dl4 / "{analysis}/bkg_exclusion.fits.gz",
    input:
        region=config_dir / "bkg_exclusion",
        script=scripts / "create_fits_exclusion.py",
        config=config_dir / "{analysis}/analysis.yaml",
    conda:
        gammapy_env
    log:
        dl4 / "{analysis}/create_exclusion.log",
    shell:
        "python {input.script}  -i {input.region} -o {output} --log-file {log} -c {input.config}"


rule create_dataset:
    output:
        dl4 / "{analysis}/datasets.fits.gz",
    input:
        data=dl3 / "hdu-index.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script=scripts / "write_datasets_3d.py",
        bkg_exclusion_regions=dl4 / "{analysis}/bkg_exclusion.fits.gz",
    conda:
        gammapy_env
    log:
        dl4 / "{analysis}/datasets.log",
    shell:
        "python {input.script} -c {input.config}  -o {output} --log-file {log}"


rule calc_dl4_diagnostics:
    output:
        dl4 / "{analysis}/dl4_diagnostics.fits.gz",
    input:
        data=dl4 / "{analysis}/datasets.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script=scripts / "calc_dl4_diagnostics.py",
    resources:
        mem_mb=16000,
    conda:
        gammapy_env
    log:
        dl4 / "{analysis}/dl4_diagnostics.log",
    shell:
        "python {input.script} -c {input.config} -o {output} --dataset-path {input.data} --log-file {log}"


rule peek_datasets:
    output:
        plots / "{analysis}/dataset_peek.pdf",
    input:
        data=dl4 / "{analysis}/datasets.fits.gz",
        script=scripts / "plot_dataset_peek.py",
        config=config_dir / "{analysis}/analysis.yaml",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        plots / "{analysis}/dataset_peek.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -c {input.config} -o {output} --dataset-path {input.data} --log-file {log}"


rule plot_dl4_dianotics:
    output:
        plots / "{analysis}/dl4_diagnostics.pdf",
    input:
        data=dl4 / "{analysis}/dl4_diagnostics.fits.gz",
        script=scripts / "plot_dl4_diagnostics.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    log:
        plots / "{analysis}/dl4_diagnostics.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"
