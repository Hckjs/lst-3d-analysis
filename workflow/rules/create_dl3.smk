# Using my fork here currently
# no clear way to swap between runwise and stacked in my workflow :/
# maybe define a function, that returns the corresponding bkg name to a run based on
# a variable. I would need to parse that from the bkgmodel config and
# also get it into the link bkg script ...
# build_dir / "background/stacked_bkg_map.fits"
# dont ask... result of my hacks, should be solved later upstream
rule calc_background:
    conda:
        background_env
    output:
        expand(
            build_dir / "background/dl3_LST-1.Run{run_id}.fits.fits",
            run_id=RUN_IDS,
        ),
    input:
        runs=expand(
            build_dir / "dl3/dl3_LST-1.Run{run_id}.fits.gz",
            run_id=RUN_IDS,
        ),
        config=bkg_config_path,
    shell:
        """
    bkgmodel --config {input.config}
        """


# Use lstchain env here to ensure we can load it
# run id is stacked only right now, but this way it can be expanded
# data=build_dir / "background/{run_id}_bkg_map.fits", # stacked
rule plot_background:
    output:
        build_dir / "plots/background/{run_id}.pdf",
    conda:
        lstchain_env
    input:
        data=build_dir / "background/dl3_LST-1.Run{run_id}.fits.fits",  #runwise
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
        script="scripts/plot_bkg.py",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule dl3:
    output:
        build_dir / "dl3/dl3_LST-1.Run{run_id}.fits.gz",
    input:
        data=build_dir / "dl2/dl2_LST-1.Run{run_id}.h5",
        irf=build_dir / "irf/irf_Run{run_id}.fits.gz",
        config=irf_config_path,
    resources:
        mem_mb=12000,
        time=30,
    conda:
        lstchain_env
    log:
        build_dir / "logs/dl3/{run_id}.log",
    shell:
        """
        lstchain_create_dl3_file  \
            --input-dl2 {input.data}  \
            --output-dl3-path $(dirname $(realpath {output}))  \
            --input-irf {input.irf}  \
            --config {input.config} \
            --gzip \
            --overwrite \
        """


# bkg = build_dir / "background/stacked_bkg_map.fits"
# dont ask... result of my hacks, should be solved later upstream
rule dl3_hdu_index:
    conda:
        lstchain_env
    output:
        build_dir / "dl3/hdu-index.fits.gz",
    input:
        runs=expand(
            build_dir / "dl3/dl3_LST-1.Run{run_id}.fits.gz",
            run_id=RUN_IDS,
        ),
        bkg=expand(
            build_dir / "background/dl3_LST-1.Run{run_id}.fits.fits",
            run_id=RUN_IDS,
        ),
    params:
        bkg_script="scripts/link_bkg.py",
        bkg_dir=lambda w, input: os.path.relpath(
            Path(input.bkg[0]).parent, Path(input.runs[0]).parent
        ),
        bkg_files=lambda w, input: [Path(x).name for x in input.bkg],
    resources:
        time=15,
    shell:
        """
        lstchain_create_dl3_index_files  \
            --input-dl3-dir {build_dir}/dl3  \
            --output-index-path {build_dir}/dl3  \
            --file-pattern 'dl3_*.fits.gz'  \
            --overwrite 

    python {params.bkg_script} \
    --hdu-index-path {output} \
    --bkg-dir {params.bkg_dir} \
    --bkg-file {params.bkg_files} 
        """


rule cuts_dl2_dl3:
    resources:
        mem_mb="64G",
        time=10,
    conda:
        lstchain_env
    output:
        build_dir / "dl3/counts/after_gh_theta_cut_{run_id}.h5",
    input:
        dl2=build_dir / "dl2/dl2_LST-1.Run{run_id}.h5",
        irf=build_dir / "irf/irf_Run{run_id}.fits.gz",
        config=irf_config_path,
        script="scripts/calc_counts_after_cuts.py",
    shell:
        "python {input.script} --input-dl2 {input.dl2} --input-irf {input.irf} -c {input.config} -o {output}"


rule stack_cuts_dl2_dl3:
    conda:
        lstchain_env
    output:
        build_dir / "dl3/counts/after_gh_theta_cut_{norm}_stacked.h5",
    input:
        data=expand(
            build_dir / "dl3/counts/after_gh_theta_cut_{run_id}.h5",
            run_id=RUN_IDS,
        ),
        script="scripts/stack_counts_after_cuts.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output} --norm {wildcards.norm}"


rule plot_cuts_dl2_dl3:
    conda:
        lstchain_env
    output:
        build_dir / "plots/counts_after_gh_theta_cut_{norm}.pdf",
    input:
        data=build_dir / "dl3/counts/after_gh_theta_cut_{norm}.h5",
        script="scripts/plot_counts_after_cuts.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule calc_skymap:
    resources:
        mem_mb="64G",
        time=10,
    conda:
        lstchain_env
    output:
        build_dir / "dl3/skymap/{run_id}.fits",
    input:
        data=build_dir / "dl2/dl2_LST-1.Run{run_id}.h5",
        config=irf_config_path,
        script="scripts/calc_skymap.py",
    shell:
        "python {input.script} -i {input.data} -o {output} -c {input.config}"


rule plot_skymap:
    conda:
        lstchain_env
    output:
        build_dir / "plots/skymap/{run_id}.pdf",
    input:
        data=build_dir / "dl3/skymap/{run_id}.fits",
        script="scripts/plot_skymap.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


rule stack_skymaps:
    conda:
        lstchain_env
    output:
        build_dir / "dl3/skymap/stacked.fits",
    input:
        data=expand(
            build_dir / "dl3/skymap/{run_id}.fits",
            run_id=RUN_IDS,
        ),
        script="scripts/stack_skymap.py",
    shell:
        "python {input.script} -i {input.data} -o {output}"


rule stack_skymaps_dl3:
    conda:
        lstchain_env
    output:
        build_dir / "dl3/skymap_dl3/stacked.fits",
    input:
        data=expand(
            build_dir / "dl3/skymap_dl3/{run_id}.fits",
            run_id=RUN_IDS,
        ),
        script="scripts/stack_skymap.py",
    shell:
        "python {input.script} -i {input.data} -o {output}"
