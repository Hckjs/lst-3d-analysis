env = ENVS["lstchain"]
bkg_env = ENVS["background"]

dl2 = Path(OUTDIRS["dl2"])
dl3 = Path(OUTDIRS["dl3"])
irfs = Path(OUTDIRS["irfs"])

scripts = Path(SCRIPTS["dl3"])

irf_config = CONFIGS["irf_tool"]
bkg_config = CONFIGS["bkg_model"]


rule dl3:
    input:
        dl3 / "hdu-index.fits.gz",


rule dl2_to_dl3:
    output:
        run=dl3 / "LST-1.Run{run_id}.dl3.fits.gz",
        irf=dl3 / "irfs_{run_id}.fits.gz",
    input:
        data=dl2 / "LST-1.Run{run_id}.dl2.h5",
        irfs=MC_NODES_IRFs,
        config=irf_config,
    params:
        interp="linear",  # nope
        irf_pattern="irfs_*.fits.gz",
        out=dl3,
        in_irfs=irfs,
    conda:
        env
    resources:
        mem_mb=12000,
        time=30,
    log:
        dl3 / "create_dl3_{run_id}.log",
    shell:
        """
        lstchain_create_dl3_file  \
            --input-dl2 {input.data}  \
            --output-dl3-path {params.out}  \
            --input-irf-path {params.in_irfs}  \
            --irf-file-pattern {params.irf_pattern} \
            --final-irf-file {output.irf} \
            --config {input.config} \
            --use-nearest-irf-node \
            --gzip \
            --overwrite \
        """


rule calc_background:
    output:
        BKG_FILES,
    input:
        runs=DL3_FILES,
        config=bkg_config,  # thats not the normal pybkgmodel config, but more keys
        script=scripts / "calc_background.py",
    conda:
        bkg_env
    log:
        dl3 / "calc_bkg.log",
    shell:
        """python {input.script} \
        --input-runs {input.runs} \
        --config {input.config} \
        --output {output} \
        --log-file {log}
        """


# Use lstchain env here to ensure we can load it
rule plot_background:
    output:
        dl3 / "plots/bkg_{run_id}.pdf",
    input:
        data=dl3 / "bkg_{run_id}.fits.gz",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
        script=scripts / "plot_bkg.py",
    conda:
        env
    log:
        dl3 / "plots/bkg_{run_id}.log",
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


# bkg = build_dir / "background/stacked_bkg_map.fits"
rule dl3_hdu_index:
    output:
        dl3 / "hdu-index.fits.gz",
    input:
        runs=DL3_FILES,
        bkg=BKG_FILES,
    params:
        bkg_script=scripts / "link_bkg.py",
        bkg_dir=dl3,
        outdir=dl3,
    conda:
        env
    log:
        dl3 / "hdu_index.log",
    resources:
        time=15,
    shell:
        """
        lstchain_create_dl3_index_files  \
            --input-dl3-dir {params.outdir}  \
            --output-index-path {params.outdir}  \
            --file-pattern '*.dl3.fits.gz'  \
            --overwrite

        python {params.bkg_script} \
        --hdu-index-path {output} \
        --bkg-dir {params.bkg_dir} \
        --bkg-file {input.bkg}
        """
