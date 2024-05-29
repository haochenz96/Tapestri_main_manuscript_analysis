mamba activate samtools

case_names=(RA17_13)

for case_i in ${case_names[@]}; do
	# Subset the WGS SNP BAF data to the ROI
	awk '{print $1, $2, $2, $3, $4, $5}' OFS='\t' "/juno/work/iacobuzc/haochen/bulk_analysis/hatchet-from_KM/${case_i}/baf/tumor.1bed" | \
		bedtools intersect \
			-a "/juno/work/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp_methods/1_validating_cn_calling/${case_i}/${case_i}_roi.bed" \
			-b stdin \
			-wao \
		> "/juno/work/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp_methods/1_validating_cn_calling/${case_i}/${case_i}_roi_wgs_snp_baf.bed"
done