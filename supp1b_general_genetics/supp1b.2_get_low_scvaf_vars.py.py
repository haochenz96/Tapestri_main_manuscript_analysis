import mosaic.io as mio
import pandas as pd
from pathlib import Path

from IPython import embed

wd = Path('/data/iacobuzc/haochen/Tapestri_batch2/analysis/general_genetics')
sc_maf = wd / "1_scMAF" / "pan_cohort_scMAF.csv"
sc_maf_df = pd.read_csv(sc_maf, index_col=0)
sc_maf_df = sc_maf_df[
    (sc_maf_df["CCF"] > 0.01) & 
    (sc_maf_df["gene_name"] != "KRAS")
    ] 
output_dir = wd / '2_tumors_pon'
output_dir.mkdir(exist_ok=True, parents=True)

# 1. create composite variant sheet

# %%
multilayer_snv_sheet = {}

multilayer_snv_sheet['mut_filtered_num_cells'] = sc_maf_df.pivot(columns='patient_name', values='mut_filtered_num_cells')
multilayer_snv_sheet['mean_sc_AF_in_mut_filtered'] = sc_maf_df.pivot(columns='patient_name', values='mean_sc_AF_in_mut_filtered')

# %%
recurrent_snvs = multilayer_snv_sheet['mut_filtered_num_cells'][
    ((~multilayer_snv_sheet['mut_filtered_num_cells'].isna()).sum(axis=1) >= 3)
].index

low_mean_af_snvs = multilayer_snv_sheet['mean_sc_AF_in_mut_filtered'][
    multilayer_snv_sheet['mean_sc_AF_in_mut_filtered'].mean(axis=1) < 30
].index

# %%
blacklist = list(set(recurrent_snvs) | set(low_mean_af_snvs))

# # %%
# v = 'chr10:131506192:C/T' 
# sc_maf_df.loc[v, ['bulk_tumor_fillout_VAF','bulk_normal_VAF']].any().any()

# %%
from tea.cravat import NONFUNC_SO

sc_maf_df_unique = sc_maf_df[~sc_maf_df.index.duplicated(keep='first')]
tumor_pon_df = pd.DataFrame(index = blacklist)

tumor_pon_df['num_patients'] = (~multilayer_snv_sheet['mut_filtered_num_cells'].loc[blacklist].isna()).sum(axis=1)
tumor_pon_df['mean_mut_filtered_num_cells'] = multilayer_snv_sheet['mut_filtered_num_cells'].loc[blacklist].mean(axis=1, skipna=True)
tumor_pon_df['mean_sc_AF_in_mut_filtered'] = multilayer_snv_sheet['mean_sc_AF_in_mut_filtered'].loc[blacklist].mean(axis=1, skipna=True)
tumor_pon_df['functional'] = ~sc_maf_df_unique['snv_class'].isin(NONFUNC_SO)

tumor_pon_df[['gene_name','snv_class','HGVSp']] = sc_maf_df_unique[['gene_name','snv_class','HGVSp']].loc[blacklist]

tumor_pon_df['detected_by_bulk'] = None
for var_i, _ in tumor_pon_df.iterrows():
    tumor_pon_df.loc[var_i, 'detected_by_bulk'] = sc_maf_df.loc[var_i, ['bulk_tumor_fillout_VAF','bulk_normal_VAF']].any().any()


# %%
tumor_pon_df.to_csv(
    output_dir / 'all_tumor_cases.pon_blacklist.csv', index=True, header=True
)



