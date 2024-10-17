# ===== target ======
rule all:
    input:

rule 0_get_tumor_pon:
    input: 
		ss_h5s = expand(
			"/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/ss_h5/{patient_name}.h5",
			patient_name = patient_names
		),
        master_sample_sheet = "Tapestri_batch2_samples_MASTER_INTERNAL.xlsx"
    output: 
        tumor_pon = "supp1_general_genetics/snv_blacklists/all_tumor_samples_composite_snv_sheet.csv"
    log:
        "analysis_log/0_get_tumor_pon.log"
    conda: "mosaic"
    resources:
        mem_mb = 8000,
        n_threads = 4,
        time_min = 49,
    shell:
        """
        python supp1b_general_genetics/supp1b.0_get_tumor_pon.py \
            > {log} 2>&1
        """

rule 1_filter_snvs_get_scMAFs:
    input: 
		ss_h5s = expand(
			"data_compiled/ss_h5/{patient_name}.h5",
			patient_name = patient_names
		),
        master_sample_sheet = "Tapestri_batch2_samples_MASTER.xlsx"
    output: 
        tumor_pon = "supp1_general_genetics/snv_blacklists/all_tumor_samples_N={N}_composite_snv_sheet.csv",
        expand("supp1b_general_genetics/all_vars_mut_prev=0.005/vois/{patient_name}_all_vars_mut_prev=0.005_voi.txt", patient_name=patient_names)
    params:
    log:
    conda: "mosaic"
    resources:
        mem_mb = 8000,
        n_threads = 4,
        time_min = 49,
    shell:
        """
        python supp1b_general_genetics/supp1b.1_filter_snvs.py \
            > {log} 2>&1
        """
