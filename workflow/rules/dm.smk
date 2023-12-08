dm_env = ENVS["dark_matter"]
dm = Path(OUTDIRS["dm"])
dl4 = Path(OUTDIRS["dl4"])
scripts = Path(SCRIPTS["dm"])
channels = ["b", "tau", "W", "mu"]  # configurable?


rule dm:
    input:
        ul_plots=expand(
            dm / f"{analysis}/plots/uls_{channel}.pdf",
            analysis=analyses,
            channel=channels,
        ),


rule gen_uls:
    output:
        dm / "{analysis}/uls_{channel}.fits.gz",
    input:
        dataset=dl4 / "{analysis}/datasets.fits.gz",
        dm_config=configs / "dm_config.yaml",
        script=scripts / "fit_uls.py",
        models_fit=dl4 / "{analysis}/bkg_fit.yaml",
    conda:
        dm_env
    log:
        dm / "{analysis}/{channel}.log",
    shell:
        """
        python {input.script} \
        --dataset {input.dataset} \
        --config-path {input.dm_config} \
        --models-path {input.models_fit} \
        --output {output} \
        --log-file {log} \
        --channel {channel} \
        --verbose
        """


rule plot_uls:
    output:
        dm / "{analysis}/plots/uls_{channel}.pdf",
    input:
        uls=dm / "{analysis}/uls_channel.fits.gz",
        script=scripts / "fit_uls.py",
    conda:
        dm_env
    log:
        dm / "{analysis}/plots/{channel}.log",
    shell:
        """
        python {input.script} \
        --dataset {input.dataset} \
        --output {output} \
        --log-file {log} \
        --verbose
        """
