# %%
import pandas as pd
import numpy as np
import mosaic.io as mio
import plotly.express as px
from pathlib import Path
import os

os.chdir(Path(__file__).parent)
# patient_names = [
#     "PC01", "PC02", "PC03", "PC04", "PC05", "PC06", "PC07", "PC08", "PC09", "PC10", "PC11", "PC12",
#     "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21-IR", "PC22", "PC23", "PC24-RSX",
#     ]

patient_names = ["PC11"]
# kras_variants = {}
for patient_name in patient_names:
    pt_h5 = f"../data_compiled/fillout_h5/{patient_name}.patient_wide.genotyped.h5"
    pt = mio.load(pt_h5)
    cravat = pd.read_csv(f"../data_compiled/fillout_h5/{patient_name}_CRAVAT_output_cleaned.txt", sep = "\t", index_col = 0, header=[0,1])
    cravat = cravat[
        (~cravat[("Variant Annotation", "Sequence Ontology")].isin(NONFUNC_SO)) & 
        (~(cravat[("PoN_comparison", "PoN-superset-normals-occurence")] > 1))
        ]

# %%
pt.dna.genotype_variants()
# %% read in high-quality SNVs
snv_f = f"../data_compiled/manual_annotated_snv_lists/{patient_name}-patient-all_vars-voi.hz_curated.txt"
snv_df = pd.read_csv(snv_f, sep = '\t', index_col = 0, comment='#')
snv_df["annotation"].fillna("", inplace=True)
# filter out artifact
snv_df = snv_df[~snv_df['annotation'].str.contains("artifact")]
# filter out germline HOM
snv_df = snv_df[~snv_df['annotation'].str.contains("germline_HOM")]
snv_df.sort_values(by=["mut_prev", "HGVSp"], ascending=False)
snv_ann_map = snv_df['HGVSp'].to_dict()

# %% ----- get ME matrix -----
VOI = snv_df[snv_df["annotation"].str.contains("somatic")].index.tolist()
VOI += ["chr12:25398284:C/T"]
VOI = [v for v in VOI if v != "chr3:178936082:G/A"]
snv_ann_map["chr12:25398284:C/T"] = "KRAS p.Gly12Asp"


from tea.snv.me import get_exclusivity
E = 1 - get_exclusivity(
    pt.dna.get_attribute('DP', constraint='row', features=VOI),
    pt.dna.get_attribute('alt_read_count', constraint='row', features=VOI),
    # rm_irrelevant_cells = False
)

# %%
import matplotlib.pyplot as plt
from tea.plots import *
mpl_config_params(font_size=8)
import seaborn as sns

E.columns = [snv_ann_map[i] for i in E.columns]
E.index = [snv_ann_map[i] for i in E.index]
# %%
fig, ax = plt.subplots(figsize=(5,5))         # Sample figsize in inches
sns.heatmap(
    E, annot=True, 
    cmap='Blues', fmt='.2f', # vmin=0, vmax=1,
    ax=ax
    )

# %%
fig.savefig(f"{patient_name}_voi_me.pdf", bbox_inches='tight')

# %% Venn diagram
from matplotlib_venn import venn2

# %%
kras1_cells = set(pt.dna.barcodes()[pt.dna.get_attribute('mut_filtered', features=["chr12:25398284:C/A"])[0].astype(bool)])
kras2_cells = set(pt.dna.barcodes()[pt.dna.get_attribute('mut_filtered', features=["chr12:25398284:C/T"])[0].astype(bool)])
fig, ax = plt.subplots(figsize=(5,5))
venn2(
    [set(kras1_cells), set(kras2_cells)], 
    ('KRAS p.Gly12Asp', 'KRAS p.Gly12Asp'),
    ax=ax
    )

fig.savefig(f"{patient_name}_kras_venn.pdf", bbox_inches='tight')

# %%
kras2 = pt.dna[list(kras2_cells), pt.dna.ids()]
# what are the mutations positive in kras2?
voi = kras2.ids()[kras2.get_attribute('mut_filtered', constraint='row').sum(axis=0) > 2]
# %%
cravat.loc[[v for v in voi if v in cravat.index]].iloc[:, :20]
# %%
