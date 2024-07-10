# ===== target ======
rule all:
    input:
        # "input/candidate_alleles.multiallelic.for_py.csv"
        expand(
            "{sample_name}/{sample_name}.mpileup.DP.merged.csv",
            sample_name=sample_names
        ),
        # merged_mpileup_vcf = expand(
        #     "{sample_name}/combined_vcf/{sample_name}-genotyped_combined.vcf.gz",
        #     sample_name = sample_names,
        # ),
        # merged_mpileup_vcf_AD_py = expand(
        #     "{sample_name}/{sample_name}-genotyped_combined_AD_for_py.txt",
        #     sample_name = sample_names
        # ),


rule 0_get_tumor_pon:
    input: 
		ss_h5s = expand(
			"/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/ss_h5/{sample_name}.h5",
			sample_name = SAMPLE_NAMES
		),
        master_sample_sheet = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/Tapestri_batch2_samples_MASTER.xlsx"
    output: 
        tumor_pon = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1_general_genetics/snv_blacklists/all_tumor_samples_N={N}_composite_snv_sheet.csv"
    params:
    log:
    conda: 
    # hardcoded
        '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/envs/mosaic-custom.yaml'
    resources:
        mem_mb = 8000,
        n_threads = 4,
        time_min = 49,
    shell:
        """
        python /data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1b_general_genetics/supp1b.0_get_tumor_pon.py \
            > {log} 2>&1
        """
rule 1_filter_snvs_get_scMAFs:
    input: 
		ss_h5s = expand(
			"/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/ss_h5/{sample_name}.h5",
			sample_name = SAMPLE_NAMES
		),
        master_sample_sheet = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/Tapestri_batch2_samples_MASTER.xlsx"
    output: 
        tumor_pon = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1_general_genetics/snv_blacklists/all_tumor_samples_N={N}_composite_snv_sheet.csv"
    params:
    log:
    conda: 
    # hardcoded
        '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/envs/mosaic-custom.yaml'
    resources:
        mem_mb = 8000,
        n_threads = 4,
        time_min = 49,
    shell:
        """
        python /lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1b_general_genetics/supp1b.1_filter_snvs.py \
            > {log} 2>&1
        """