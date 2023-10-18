gammapy_env = ENVS["gammapy"]
dl4 = Path(OUTDIRS["dl4"])
dl5 = Path(OUTDIRS["dl5"])
scripts = Path(SCRIPTS["dl5"])

dl5_plot_types = ["significance_map", "2d_flux_profile"]  # , flux_points, light_curve]


rule dl5:
    input:
        [
            dl5 / f"{analysis}/plots/{plot}.pdf"
            for analysis in analyses
            for plot in dl5_plot_types
        ],


rule calc_significance_map:
    output:
        dl5 / "{analysis}/significance_map.fits.gz",
    input:
        data=dl4 / "{analysis}/datasets.fits.gz",
        script=scripts / "calc_significance_map.py",
    conda:
        gammapy_env
    log:
        dl5 / "{analysis}/calc_significance_map.log",
    shell:
        """
        python {input.script} \
        --dataset-path {input.data} \
        --output {output} \
        --log-file {log}
        """


rule plot_significance_map:
    output:
        dl5 / "{analysis}/plots/significance_map.pdf",
    input:
        lima_map=dl5 / "{analysis}/significance_map.fits.gz",
        script=scripts / "plot_significance_map.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        dl5 / "{analysis}/plots/plot_significance_map.log",
    shell:
        """
        MATPLOTLIBRC={input.rc} \
        python {input.script} \
        --flux-maps {input.lima_map} \
        --output {output} \
        --log-file {log}
        """


rule plot_significance_distribution:
    output:
        dl5 / "{analysis}/plots/significance_distribution.pdf",
    input:
        lima_map=dl5 / "{analysis}/significance_map.fits.gz",
        script=scripts / "plot_significance_distribution.py",
        rc=MATPLOTLIBRC,
        exclusion_mask=dl4 / "{analysis}/exclusion.fits.gz",
    conda:
        gammapy_env
    log:
        dl5 / "{analysis}/plots/plot_significance_distribution.log",
    shell:
        """
        MATPLOTLIBRC={input.rc} \
        python {input.script} \
        --input-maps {input.lima_map} \
        --exclusion-mask {input.exclusion_mask} \
        --output {output} \
        --log-file {log}
        """


rule calc_2d_flux_profile:
    output:
        dl5 / "{analysis}/2d_flux_profile.fits.gz",
    input:
        data=dl4 / "{analysis}/datasets.fits.gz",
        script=scripts / "calc_2d_flux_profile.py",
    conda:
        gammapy_env
    log:
        dl5 / "{analysis}/calc_2d_flux_profile.log",
    shell:
        """
        python {input.script} \
        --dataset-path {input.data} \
        --output {output} \
        --log-file {log}
        """


rule plot_2d_flux_profile:
    output:
        dl5 / "{analysis}/plots/2d_flux_profile.pdf",
    input:
        flux_points=dl5 / "{analysis}/2d_flux_profile.fits.gz",
        script=scripts / "plot_2d_flux_profile.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        dl5 / "{analysis}/plots/plot_2d_flux_profile.log",
    shell:
        """
        MATPLOTLIBRC={input.rc} \
        python {input.script} \
        --flux-points {input.flux_points} \
        --output {output} \
        --log-file {log}
        """
