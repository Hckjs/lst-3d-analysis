env = ENVS["lstchain"]
# Having these as paths makes them easier to use
scripts = Path(SCRIPTS["irfs"])
out = Path(OUTDIRS["irfs"])
dl2 = Path(OUTDIRS["dl2"])
plots = Path(PLOTSDIRS["dl2"])
config = CONFIGS["irf_tool"]


rule irf:
    output:
        out / "irf_Run{run_id}.fits.gz",
    input:
        gammas=dl2 / "test/dl2_LST-1.Run{run_id}.h5",
        config=irf_config_path,
    resources:
        mem_mb=8000,
        time=10,
    conda:
        lstchain_env
    shell:
        """
        lstchain_create_irf_files \
            -o {output} \
            -g {input.gammas} \
            --config {input.config} \
        """


rule plot_irf:
    output:
        plots / "irf/{irf}_Run{run_id}.pdf",
    input:
        data=irf / "irf_Run{run_id}.fits.gz",
        script=scripts / "plot_irf_{irf}.py",
        rc=MATPLOTLIBRC,
    resources:
        mem_mb=1000,
        time=20,  # minutes
    conda:
        gammapy_env
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"
