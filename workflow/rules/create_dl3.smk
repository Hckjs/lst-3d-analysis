rule merge_gamma_mc_per_node:
    output:
        train=build_dir / "merged/GammaDiffuse/{node}_train.dl1.h5",
        test=build_dir / "merged/GammaDiffuse/{node}_test.dl1.h5",
    params:
        train_size=0.5,
        directory=lambda node: build_dir / f"mc_nodes/GammaDiffuse/{node}"
    conda:
        lstchain_env
    shell:
        """
        python scripts/merge_mc_nodes.py \
	--input-dir {params.directory} \
	--train-size {params.train_size} \
	--output-train {output.train} \
	--output-test {output.test} 
	"""

# TODO This should be a function getting the wildcard particle, but right now its only one particle, so its fine
from pathlib import Path
all_gamma_nodes = [ x.name for x in (build_dir / "mc_nodes/GammaDiffuse").glob("*") if x.is_dir()]

rule merge_train_or_test_of_all_nodes:
    output:
        build_dir / "dl1/{train_or_test}/{particle}_train.dl1.h5"
    input:
        files=expand(build_dir / "merged/{{particle}}/{node}_{{train_or_test}}.dl1.h5", node=all_gamma_nodes),
    params:
        directory=lambda wildcards: build_dir / f"merged/{wildcards.particle}",
        pattern=lambda wildcards: f"*_{wildcards.train_or_test}.dl1.h5",
	out_type=lambda wildcards: f"output-{wildcards.train_or_test}",
    conda:
        lstchain_env
    shell:
        """
        python scripts/merge_mc_nodes.py \
	--input-dir {params.directory} \
        --pattern {params.pattern} \
	--{params.out_type} {output}
	"""

rule train_models:
    output:
        build_dir/"models/reg_energy.sav",
        build_dir/"models/cls_gh.sav",
        build_dir/"models/reg_disp_norm.sav",
        build_dir/"models/cls_disp_sign.sav",
    input:
        gamma=build_dir / "dl1/train/GammaDiffuse_train.dl1.h5",
        proton=build_dir / "dl1/train/proton_diffuse_merged.dl1.h5",
        config=build_dir / "models/mcpipe/lstchain_config.json",
    resources:
        mem_mb=64000,
        cpus=8,
	partition="long",
	time=1200,
    params:
        outdir=build_dir / "models"
    conda:
        lstchain_env
    shell:
        """
	lstchain_mc_trainpipe \
	--fg {input.gamma} \
	--fp {input.proton} \
	--config {input.config} \
	--output-dir {params.outdir}
	"""

# Cant do this in link script because the merge into train and test did not happen yet at that stage
rule link_test_nodes:
    input:
        runs=expand(
            build_dir / "dl1/dl1_LST-1.Run{run_id}.h5",
            run_id=RUN_IDS,
        ),
        dir = build_dir / "dl1",
        test_files=expand(build_dir / "merged/GammaDiffuse/{node}_test.dl1.h5", node=all_gamma_nodes),
    output:
        runs=expand(
            build_dir / "dl1/test/dl1_LST-1.Run{run_id}.h5",
            run_id=RUN_IDS,
        ),
    params:
        test_dir = build_dir / "dl1/test",
        mc_nodes_dir = build_dir / "merged/GammaDiffuse",
    conda:
        lstchain_env
    shell:
        """
	python link-dl1.py \
	--dl1-dir {input.dir} \
	--test-link-dir {params.test_dir} \
	--mc-nodes-dir {params.mc_nodes_dir}
	"""
        


# Adapt this to match for the diff gamma test set! needs a wildcard for the dir probably
rule dl2:
    resources:
        mem_mb=64000,
        cpus=4
    output:
        build_dir / "dl2/{potentially_test}dl2_LST-1.Run{run_id}.h5",
    input:
        data=build_dir / "dl1/{potentially_test}dl1_LST-1.Run{run_id}.h5",
        config=build_dir / "models/mcpipe/lstchain_config.json",
        e_model=build_dir/"models/reg_energy.sav",
        gh_model=build_dir/"models/cls_gh.sav",
        disp_model=build_dir/"models/reg_disp_norm.sav",
        sign_model=build_dir/"models/cls_disp_sign.sav",
        model_dir=build_dir / "models"
    conda:
        lstchain_env
    # allow wildcard to be empty
    wildcard_constraints:
        potentially_test = ".*"
    shell:
        """
        lstchain_dl1_to_dl2  \
            --input-file {input.data}  \
            --output-dir $(dirname {output}) \
            --path-models {input.model_dir}  \
            --config {input.config} 
        """


rule irf:
    resources:
        mem_mb=8000,
        time=10,
    output:
        build_dir / "irf/irf_Run{run_id}.fits.gz",
    input:
        gammas=build_dir / "dl2/test/dl2_LST-1.Run{run_id}.h5",
        config=irf_config_path,
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
        build_dir / "plots/irf/{irf}_Run{run_id}.pdf",
    input:
        data=build_dir / "irf/irf_Run{run_id}.fits.gz",
        script="scripts/plot_irf_{irf}.py",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
    resources:
        mem_mb=1000,
        time=5,  # minutes
    conda:
        gammapy_env
    shell:
        "MATPLOTLIBRC={input.rc} python {input.script} -i {input.data} -o {output}"


# Using my fork here currently
# One stacked map for now. Maybe per run later?
# Need to check again how this is done in pybkgmodel
rule calc_background:
    conda:
        background_env
    output:
        build_dir / "background/stacked_bkg_map.fits"
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
rule plot_background:
    output:
        build_dir / "plots/background/{run_id}.pdf",
    conda:
        lstchain_env
    input:
        data=build_dir / "background/{run_id}_bkg_map.fits",
        rc=os.environ.get("MATPLOTLIBRC", config_dir / "matplotlibrc"),
	script="scripts/plot_bkg.py"
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
	bkg_file=build_dir / "background/stacked_bkg_map.fits"
    params:
        bkg_dir = build_dir / "background",
        bkg_script = "scripts/link_bkg.py",
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
	--bkg-file {input.bkg_file}
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

