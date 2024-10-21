# ===== 2 fast-condor =====
conda activate condor

FAST_CONDOR_SCRIPT=/data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/fast-ConDoR/src/fast-condor.py

cd /data/iacobuzc/haochen/condor_pipeline_trial

# manual test run
PATIENT_NAME=PC03

python /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/fast-ConDoR/src/fast-condor.py \
    -i condor_inputs/${PATIENT_NAME}/character_bin_matrix.csv\
    -r condor_inputs/${PATIENT_NAME}/total_readcounts.csv \
    -v condor_inputs/${PATIENT_NAME}/alt_readcounts.csv  \
    -s condor_inputs/${PATIENT_NAME}/germline_mutations.txt  \
    -s2 condor_inputs/${PATIENT_NAME}/somatic_mutations.txt -o condor_outputs/${PATIENT_NAME}/out \
    -m /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/references/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.gene_pos_ado_added.csv \
    -d ${PATIENT_NAME} \
    --scr \
    --cnp /data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/falcon_solutions/${PATIENT_NAME}.unique_cn_clone_profiles.csv 
    # --subclonal_mutations /data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/manual_subclonal_snvs/${PATIENT_NAME}.subclonal_mutations.yaml 

# # ===== solution heatmap =====
# conda activate condor
# cd /data/iacobuzc/haochen/condor_pipeline_trial
# python /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/scripts/3a_generate_solution_heatmaps.py \
#     -d ${PATIENT_NAME} \
#     -c condor_inputs/${PATIENT_NAME}/character_bin_matrix.csv \
#     -v condor_inputs/${PATIENT_NAME}/character_vaf_matrix.csv \
#     -s condor_outputs/${PATIENT_NAME}/out_B.csv \
#     -g condor_inputs/${PATIENT_NAME}/germline_mutations.txt \
#     -m condor_inputs/${PATIENT_NAME}/somatic_mutations.txt \
#     -t condor_inputs/${PATIENT_NAME}/total_readcounts.csv \
#     -a condor_inputs/${PATIENT_NAME}/alt_readcounts.csv \
#     -o condor_outputs/${PATIENT_NAME}/heatmaps/condor_solution_heatmap.png \
#     -p condor_outputs/${PATIENT_NAME}/heatmaps/vaf_heatmap.png \
#     -i /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/references/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.gene_pos_ado_added.csv \
#     -snvs /lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/manual_annotated_snv_lists/${PATIENT_NAME}-patient-all_vars-voi.hz_curated.txt

# # ===== post-condor sc-heatmap =====
# conda activate mosaic
# cd /data/iacobuzc/haochen/condor_pipeline_trial
# python /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/scripts/1b_generate_sc_heatmaps.py \
#     --patient_name ${PATIENT_NAME} \
#     --patient_h5 /data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/fillout_h5/${PATIENT_NAME}.patient_wide.genotyped.h5 \
#     --snv_f /lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/manual_annotated_snv_lists/${PATIENT_NAME}-patient-all_vars-voi.hz_curated.txt \
#     --clone_assignment condor_downstream/${PATIENT_NAME}/${PATIENT_NAME}_final_sc_clone_assignment.csv \
#     --output_dir post_condor_sc_heatmaps \
#     --post_condor 