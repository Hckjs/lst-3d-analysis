rule calc_significance_map:
    output:
        build_dir / "dl4/{analysis}/significance_map.fits.gz",
    input:
        data=build_dir / "dl4/{analysis}/datasets.fits.gz",
        script="scripts/calc_significance_map.py",
        config=config_dir / "{analysis}/analysis.yaml",
    conda:
        gammapy_env
    shell:
        "python {input.script} -c {input.config} -o {output} --dataset-path {input.data}"


rule plot_significance_map:
    output:
        build_dir / "plots/{analysis}/significance_map.pdf",
    input:
        lima_map=build_dir / "dl4/{analysis}/significance_map.fits.gz",
        exclusion_mask=build_dir / "{analysis}/exclusion.fits.gz",
        script="scripts/plot_significance_map.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    conda:
        gammapy_env
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} --lima-maps-input {input.lima_map} --exclusion-map-input {input.exclusion_mask} -o {output}"


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
