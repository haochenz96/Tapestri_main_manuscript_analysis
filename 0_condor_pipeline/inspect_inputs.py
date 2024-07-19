# %% 
import pandas as pd

og = "/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/0_condor_pipeline/condor_inputs/RA17_13/character_vaf_matrix.csv"
og = pd.read_csv(og, index_col=0)
nw = "/data/iacobuzc/haochen/Tapestri_batch2/analysis/cn_clones/RA17_13/RA17_13-cohort-cn_calling-falcon_AJ_wo_MYC/RA17_13-cn_call_with_homdel/selected_solution/RA17_13_homdel_nclones=9.sample_sc_clone_assignment.updated.csv"
nw = pd.read_csv(nw, index_col=0)

# %%