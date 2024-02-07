gammapy_env = ENVS["gammapy"]
dl3 = Path(OUTDIRS["dl3"])
dl4 = Path(OUTDIRS["dl4"])
scripts = Path(SCRIPTS["dl4"])

dl4_plot_types = ["dataset_peek"]  # , "dl4_diagnostics"]
dl4_plots = [
    dl4 / f"{analysis}/plots/{plot}.pdf"
    for analysis in analyses
    for plot in dl4_plot_types
]


rule dl4:
    input:
        dl4_plots,


rule create_fov_bkg_exclusion:
    output:
        dl4 / "{analysis}/bkg_exclusion.fits.gz",
    input:
        region=config_dir / "{analysis}/bkg_exclusion",
        script=scripts / "create_fits_exclusion.py",
        config=config_dir / "{analysis}/analysis.yaml",
    conda:
        gammapy_env
    log:
        dl4 / "{analysis}/create_exclusion.log",
    shell:
        "python {input.script}  -i {input.region} -o {output} --log-file {log} -c {input.config} --verbose"


rule create_dataset:
    output:
        datasets=dl4 / "{analysis}/datasets.fits.gz",
        bkg_fit=dl4 / "{analysis}/bkg_fit.yaml",
    input:
        data=dl3 / "{analysis}/hdu-index.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        script=scripts / "write_datasets_3d_manual.py",
        bkg_exclusion_regions=dl4 / "{analysis}/bkg_exclusion.fits.gz",
    conda:
        gammapy_env
    resources:
        cpus=1,
        mem_mb=32000,
    log:
        dl4 / "{analysis}/datasets.log",
    shell:
        "python {input.script} -c {input.config}  -o {output.datasets} -m {output.bkg_fit} --log-file {log} --n-jobs {resources.cpus} --verbose"


rule calc_dl4_diagnostics:
    output:
        dl4 / "{analysis}/dl4_diagnostics.fits.gz",
    input:
        data=dl4 / "{analysis}/datasets.fits.gz",
        bkg_fit=dl4 / "{analysis}/bkg_fit.yaml",
        config=config_dir / "{analysis}/analysis.yaml",
        script=scripts / "calc_dl4_diagnostics.py",
    resources:
        mem_mb=16000,
    conda:
        gammapy_env
    log:
        dl4 / "{analysis}/dl4_diagnostics.log",
    shell:
        "python {input.script} -c {input.config} -o {output} --datasets-path {input.data} --models-path {input.bkg_fit} --log-file {log}"


rule peek_datasets:
    output:
        dl4 / "{analysis}/plots/dataset_peek.pdf",
    input:
        data=dl4 / "{analysis}/datasets.fits.gz",
        bkg_fit=dl4 / "{analysis}/bkg_fit.yaml",
        script=scripts / "plot_dataset_peek.py",
        config=config_dir / "{analysis}/analysis.yaml",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        dl4 / "{analysis}/plots/dataset_peek.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -c {input.config} -o {output} --datasets-path {input.data} --models-path {input.bkg_fit} --log-file {log}"


rule plot_dl4_dianotics:
    output:
        dl4 / "{analysis}/plots/dl4_diagnostics.pdf",
    input:
        data=dl4 / "{analysis}/dl4_diagnostics.fits.gz",
        script=scripts / "plot_dl4_diagnostics.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    log:
        dl4 / "{analysis}/plots/dl4_diagnostics.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --log-file {log}"
