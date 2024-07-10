# %% 
import pandas as pd
import os
from glob import glob
# %%
cravat_dir = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/cravat"
cravat_fs = glob(os.path.join(cravat_dir, "*.txt"))

assert len(cravat_fs) == 24

cravat_dfs = [pd.read_csv(f, sep='\t', header=[0,1], index_col=0) for f in cravat_fs]
# for df in cravat_dfs:
# 	df = df.reset_index()
cols_to_exclude = ["bulk_comparison", "Original Input", "Tapestri_result"]
master_df = pd.concat(
    [df for df in cravat_dfs], 
    axis=0)

# %% drop sample specific columns
cols_to_exclude = ["bulk_comparison", "Original Input", "Tapestri_result"]
master_df = master_df.drop(columns=[c for c in master_df.columns if c[0] in cols_to_exclude])

# drop duplicates by index name
master_df = master_df[~master_df.index.duplicated(keep='first')]
# %%
master_df.to_csv("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1b_general_genetics/master_cravat.csv", index=True)
# %%
