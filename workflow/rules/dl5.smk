gammapy_env = ENVS["gammapy"]
dl4 = Path(OUTDIRS["dl4"])
dl5 = Path(OUTDIRS["dl5"])
plots = dl5 / "plots"
scripts = Path(SCRIPTS["dl5"])

dl5_plot_types = [
    "significance_map",
]  # 2d_flux_profile, flux_points, light_curve]


rule dl5:
    input:
        [
            plots / f"{analysis}/{plot}.pdf"
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
        plots / "{analysis}/significance_map.pdf",
    input:
        lima_map=dl5 / "{analysis}/significance_map.fits.gz",
        script=scripts / "plot_significance_map.py",
        rc=MATPLOTLIBRC,
    conda:
        gammapy_env
    log:
        dl5 / "{analysis}/plot_significance_map.log",
    shell:
        """
        MATPLOTLIBRC={input.rc} \
        python {input.script} \
        --input {input.lima_map} \
        --output {output} \
        --log-file {log}
        """


rule plot_significance_distribution:
    output:
        plots / "{analysis}/significance_distribution.pdf",
    input:
        lima_map=dl5 / "{analysis}/significance_map.fits.gz",
        script=scripts / "plot_significance_distribution.py",
        rc=MATPLOTLIBRC,
        exclusion_mask=dl4 / "{analysis}/exclusion.fits.gz",
    conda:
        gammapy_env
    log:
        dl5 / "{analysis}/plot_significance_distribution.log",
    shell:
        """
        MATPLOTLIBRC={input.rc} \
        python {input.script} \
        --input-maps {input.lima_map} \
        --exclusion-mask {input.exclusion_mask} \
        --output {output} \
        --log-file {log}
        """


rule model_best_fit:
    input:
        config=config_dir / "{analysis}/analysis.yaml",
        dataset=build_dir / "dl4/{analysis}/datasets.fits.gz",
        model=config_dir / "{analysis}/models.yaml",
        script="scripts/fit-model.py",
    output:
        build_dir / "dl4/{analysis}/model-best-fit.yaml",
    conda:
        gammapy_env
    shell:
        """
        python {input.script} \
            -c {input.config} \
            --dataset-path {input.dataset} \
            --model-config {input.model} \
            -o {output} \
        """


# Fit flux etc.
rule calc_flux_points:
    input:
        data=build_dir / "dl4/{analysis}/datasets.fits.gz",
        config=config_dir / "{analysis}/analysis.yaml",
        model=build_dir / "dl4/{analysis}/model-best-fit.yaml",
        script="scripts/calc_flux_points.py",
    output:
        build_dir / "dl4/{analysis}/flux_points.fits.gz",
    conda:
        gammapy_env
    shell:
        """
        python {input.script} \
            -c {input.config} \
            --dataset-path {input.data} \
            --best-model-path {input.model} \
            -o {output}
        """


rule plot_flux_points:
    input:
        data=build_dir / "dl4/{analysis}/flux_points.fits.gz",
        model=build_dir / "dl4/{analysis}/model-best-fit.yaml",
        script="scripts/plot_flux_points.py",
    output:
        build_dir / "plots/{analysis}/flux_points.pdf",
    conda:
        gammapy_env
    shell:
        """
        python {input.script} \
            -i {input.data} \
            --best-model-path {input.model} \
            -o {output}
        """


rule calc_light_curve:
    input:
        model=build_dir / "dl4/{analysis}/model-best-fit.yaml",
        config=config_dir / "{analysis}/analysis.yaml",
        dataset=build_dir / "dl4/{analysis}/datasets.fits.gz",
        script="scripts/calc_light_curve.py",
    output:
        build_dir / "dl4/{analysis}/light_curve.fits.gz",
    conda:
        gammapy_env
    shell:
        """
        python {input.script} \
            -c {input.config} \
            --dataset-path {input.dataset} \
            --best-model-path {input.model} \
            -o {output} \
        """
