# %%
# @HZ 2024-05-04
# this should be a concatenation of unfiltered mutation lists for each SAMPLE

import mosaic.io as mio
import pandas as pd
from pathlib import Path

# %% 
MASTER_SAMPLE_SHEET = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/Tapestri_batch2_samples_MASTER.xlsx"
H5_DIR = Path("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/ss_h5")
OUTDIR=Path("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1_general_genetics/snv_blacklists")
# ----- patient-sample map -----
master_sample_df = pd.read_excel(
	MASTER_SAMPLE_SHEET, 
	sheet_name = "all_sample_clinical_info"
	)
master_sample_df = master_sample_df[master_sample_df["censor"] == 0]
sample_map = {}
sample_h5s = {}
for patient_i in master_sample_df['case'].unique(): 
	case_samples = master_sample_df[master_sample_df['case'] == patient_i]['sample']
	sample_map[patient_i] = case_samples.tolist()
	for sample_i in case_samples:
		h5 = list(Path(H5_DIR).glob(f"*{sample_i}.mpileup.h5"))
		if not len(h5) == 1:
			raise ValueError(f"!= 1 H5 found for {sample_i}")
		sample_h5s[sample_i] = h5[0]

sample_names = master_sample_df["sample"].unique()
sample_objs = {}
for sample_i in sample_names:
    h5 = str(sample_h5s[sample_i])
    sample_objs[sample_i] = mio.load(h5)
    sample_objs[sample_i].dna.genotype_variants(
        min_dp = 8,
        min_alt_read = 3,
        assign_low_conf_genotype = True,
    )

# 1. create composite variant sheet
snv_prev_sheets = {}
for sample_i in sample_names:
    snv_prev_sheets[sample_i] = sample_objs[sample_i].dna.get_attribute('mut_filtered', constraint='row').sum(axis=0)
    snv_prev_sheets[sample_i] = snv_prev_sheets[sample_i][
        snv_prev_sheets[sample_i] >= 3
    ] # take only SNVs present in at least 3 single cells
    snv_prev_sheets[sample_i] /= sample_objs[sample_i].dna.shape[0] # get mut_proportion
N=len(snv_prev_sheets)
concat_snv_sheet = pd.concat(snv_prev_sheets, axis=1)
concat_snv_sheet['all_tumor_samples_occurence'] = (concat_snv_sheet > 0).sum(axis=1)
concat_snv_sheet['all_tumor_samples_median_sc_prev'] = concat_snv_sheet.median(axis=1)
concat_snv_sheet.index.name = 'condensed_format'
# write output
concat_snv_sheet.to_csv(
    OUTDIR / f'all_tumor_samples_N={N}_composite_snv_sheet.csv', index=True, header=True
)
# # embed()
# # write filtered output
# concat_snv_sheet.loc[
#     (concat_snv_sheet['all_tumor_samples_occurence'] > 5) & 
#     (concat_snv_sheet['all_tumor_samples_median_sc_prev'] < 0.3)
# ].to_csv(
#     OUTDIR / 'all_tumor_samples_composite_snv_sheet.filtered.csv', index=True, header=True
# )





# %%
