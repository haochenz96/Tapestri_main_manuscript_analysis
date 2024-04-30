# ===== raw data =====
# # Tapestri batch2
# rsync -avzu \
# 	--exclude "raw_read_counts/" \
# 	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_compiled \
# 	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled

# # Tapestri batch3
# rsync -avzu \
# 	-L \
# 	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch3/batch3_data_compiled/ \
# 	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled

# # Tapestri BRCA
# rsync -avzu \
# 	-L \
# 	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/ \
# 	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled \

# ===== condor_pipeline =====
# # Tapestri batch2
# rsync -avzu \
# 	-L \
# 	--exclude "TP12/"
# 	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch2/analysis/condor-pipeline/condor_inputs \
# 	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/condor_pipeline
# rsync -avzu \
# 	-L \
# 	--exclude "TP12/"
# 	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch2/analysis/condor-pipeline/condor_outputs \
# 	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/condor_pipeline
rsync -avzu \
	-L \
	--exclude "**TP12**" \
	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch2/analysis/condor-pipeline/condor_downstream \
	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/condor_pipeline 

# # Tapestri batch3
# rsync -avzu \
# 	-L \
# 	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch3/analysis/condor-pipeline/condor_inputs \
# 	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/condor_pipeline
# rsync -avzu \
# 	-L \
# 	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch3/analysis/condor-pipeline/condor_outputs \
# 	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/condor_pipeline
# rsync -avzu \
# 	-L \
# 	zhangh5@lilac-xfer01:/data/iacobuzc/haochen/Tapestri_batch3/analysis/condor-pipeline/condor_downstream \
# 	/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/condor_pipeline
