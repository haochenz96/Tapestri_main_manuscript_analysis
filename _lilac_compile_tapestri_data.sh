COMPILED_DATA_DIR=/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled
mkdir -p ${COMPILED_DATA_DIR}/fillout_h5
mkdir -p ${COMPILED_DATA_DIR}/falcon_solutions
mkdir -p ${COMPILED_DATA_DIR}/manual_annotated_snv_lists
mkdir -p ${COMPILED_DATA_DIR}/cravat

# # ===== Caitlin_Aki cases =====
# CAITLIN_AKI_CASES=(M04 M07 M11 M12 M13 RA15_06 RA15_16 RA16_08 RA16_29 RA17_13 RA17_22 RA19_02 RA19_10 RA19_21 RA20_05 RA21_17) # N=15
# CAITLIN_AKI_h5_dir=/lila/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_compiled/fillout_h5
# CAITLIN_AKI_falcon_dir=/lila/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_compiled/falcon_results
# CAITLIN_AKI_annotated_snv_lists_dir=/lila/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_compiled/manual_annotated_snv_lists
# for case_i in ${CAITLIN_AKI_CASES[@]}
# do
#     for i in $(find -L ${CAITLIN_AKI_h5_dir} -name ${case_i}*.patient_wide.genotyped.h5)
#     do
#         rsync -avzu -L $i ${COMPILED_DATA_DIR}/fillout_h5/
#     done
#     for i in $(find -L ${CAITLIN_AKI_h5_dir} -name ${case_i}*CRAVAT_output_cleaned.txt)
#     do
#         rsync -avzu -L $i ${COMPILED_DATA_DIR}/cravat
#     done

#     # for i in $(find -L ${CAITLIN_AKI_falcon_dir}/ \
#     #     -name "${case_i}*assignment.updated.csv" -o -name *unique*clone_profile*.csv -o -name *png)
#     # do
#     #     rsync -avzu -L $i ${COMPILED_DATA_DIR}/falcon_solutions/
#     # done

#     # for snv_f in $(find -L ${CAITLIN_AKI_annotated_snv_lists_dir}/ -name ${case_i}*.txt)
#     # do
#     #     rsync -avzu -L $snv_f ${COMPILED_DATA_DIR}/manual_annotated_snv_lists/
#     # done
# done


# ===== BRCA cases =====
BRCA_CASES=(BPA-1 BPA-2 BPA-3 BPA-5)
BRCA_snv_pipeline_output_dir=/lila/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/fillout_h5/
BRCA_falcon_output_dir=/lila/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/falcon_solutions/
BRCA_annotated_snv_lists_dir=/lila/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/manual_annotated_snv_lists/

for case_i in ${BRCA_CASES[@]}
do
	# rsync -avzu -L ${BRCA_snv_pipeline_output_dir} ${COMPILED_DATA_DIR}/fillout_h5/
	# rsync -avzu -L ${BRCA_falcon_output_dir} ${COMPILED_DATA_DIR}/falcon_solutions/
	rsync -avzu -L ${BRCA_annotated_snv_lists_dir} ${COMPILED_DATA_DIR}/manual_annotated_snv_lists/
done

# # rename BPA-2 and BPA-5's names
# cd /lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled
# rename --verbose BPA-2-IR BPA-2 */*
# rename --verbose BPA-5-RSX BPA-5 */*

# mv fillout_h5/BPA-2.mpileup.renamed.h5 fillout_h5/BPA-2.patient_wide.genotyped.h5
# mv fillout_h5/BPA-5.mpileup.renamed.h5 fillout_h5/BPA-5.patient_wide.genotyped.h5

# # ===== TOPCOAT cases =====
# TOPCOAT_CASES=(TP11 TP12 TP49)
# TOPCOAT_snv_pipeline_output_dir=/lila/data/iacobuzc/haochen/Tapestri_batch3/batch3_data_compiled/fillout_h5/
# TOPCOAT_falcon_output_dir=/lila/data/iacobuzc/haochen/Tapestri_batch3/batch3_data_compiled/falcon_solutions/
# TOPCOAT_annotated_snv_lists_dir=/lila/data/iacobuzc/haochen/Tapestri_batch3/batch3_data_compiled/manual_annotated_snv_lists/

# for case_i in ${TOPCOAT_CASES[@]}
# do
# 	# rsync -avzu -L ${TOPCOAT_snv_pipeline_output_dir} ${COMPILED_DATA_DIR}/fillout_h5/
# 	rsync -avzu -L ${TOPCOAT_falcon_output_dir} ${COMPILED_DATA_DIR}/falcon_solutions/
# 	# rsync -avzu -L ${TOPCOAT_annotated_snv_lists_dir} ${COMPILED_DATA_DIR}/manual_annotated_snv_lists/
# done
