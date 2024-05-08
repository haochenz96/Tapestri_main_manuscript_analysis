# %% import things
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyparsing import line
import seaborn as sns
import os
from tea.cravat import NONFUNC_SO
import pickle

# cd to the directory of the script
WD = Path(__file__).parent
os.chdir(WD)
pan_cohort_scMAF = pd.read_csv("pan_cohort_somatic_scMAF_mut_prev=0.005.csv")

def config_params(font_size=7):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'

MUT_PREV=0.005

# %%
GOI=["KRAS", "TP53","CDKN2A", "SMAD4", "SMAD2", "SMAD3", "TGFBR1", "TGFBR2", "ACVR1B", "BMPR1A", "ARID1A", "ARID2", "BRCA2", "ATM", "BAP1", "PIK3CA", "FGFR1","RNF43", "POLD1", "IRF6", "GATA6", "MYC", "MTOR"]

pan_cohort_scMAF["functional?"] = ~pan_cohort_scMAF["snv_class"].isin(NONFUNC_SO)
pan_cohort_scMAF["driver?"] = pan_cohort_scMAF["gene_name"].isin(GOI)
pan_cohort_scMAF["functional_driver?"] = (~pan_cohort_scMAF["snv_class"].isin(NONFUNC_SO)) & (pan_cohort_scMAF["gene_name"].isin(GOI))
import seaborn as sns
import matplotlib.pyplot as plt
# use Seaborn stacked barplot - make a histogram of CCF, colored by functional status
# log scale

# pan_cohort_scMAF_filtered = pan_cohort_scMAF[pan_cohort_scMAF["CCF"]>0]

config_params(font_size=12)

fig, ax= plt.subplots(1, 3, figsize=(18, 5))
# 1. functional mutations to driver genes
sns.histplot(
	pan_cohort_scMAF[pan_cohort_scMAF["functional_driver?"]],
	x="mut_filtered_sc_freq",
	multiple="stack",
	bins=100,
	color="red", linewidth=0,
	ax=ax[0],
)
ax[0].set_title("Functional SNVs to driver genes")
ax[0].set_yscale("log", base=2)
# set y axis range
ax[0].set_ylim(1, 2**11)
# 2. all mutations to driver genes
sns.histplot(
	pan_cohort_scMAF[pan_cohort_scMAF["driver?"]],
	x="mut_filtered_sc_freq",
	multiple="stack",
	bins=100,
	color="blue", linewidth=0,
	ax=ax[1],
)
ax[1].set_title("All SNVs to driver genes")
ax[1].set_yscale("log", base=2)
ax[1].set_ylim(1, 2**11)
# 3. all mutations to all genes
sns.histplot(
	pan_cohort_scMAF,
	x="mut_filtered_sc_freq",
	multiple="stack",
	bins=100,
	color="black", linewidth=0,
	ax=ax[2],
)
ax[2].set_title("All SNVs")
ax[2].set_yscale("log", base=2)
ax[2].set_ylim(1, 2**11)

# %%
fig.savefig(WD / f"pan_cohort_somatic_scMAF_mut_prev={MUT_PREV}_CCF_hist.pdf")

# %% Calculate dN/dS
# load L matrix
with open("dnds/panel3359_gene_L_mat_map.pkl", "rb") as f:
	L = pickle.load(f)
GOI_mutation_type_counts = pd.concat([L[gene_i].sum(axis=0) for gene_i in GOI],axis=1).sum(axis=1)
panel_wide_mutation_type_counts = pd.concat([L[gene_i].sum(axis=0) for gene_i in L.keys()],axis=1).sum(axis=1)
NONCODING_SO = ["intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "splice_site_variant", "2kb_upstream_variant"]

# calculate GOI dN/dS
GOI_dN = pan_cohort_scMAF[
	(pan_cohort_scMAF["gene_name"].isin(GOI)) &
	(~pan_cohort_scMAF["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
].shape[0] / GOI_mutation_type_counts[["Missense", "Nonsense"]].sum()
GOI_dS = pan_cohort_scMAF[
	(pan_cohort_scMAF["gene_name"].isin(GOI)) &
	(pan_cohort_scMAF["snv_class"].isin(["synonymous_variant"]))
].shape[0] / GOI_mutation_type_counts[["Synonymous"]].sum()

print(f"GOI dN/dS: {GOI_dN/GOI_dS}")

# calculate panel-wide dN/dS
panel_wide_dN = pan_cohort_scMAF[
	(~pan_cohort_scMAF["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
].shape[0] / panel_wide_mutation_type_counts[["Missense", "Nonsense"]].sum()
panel_wide_dS = pan_cohort_scMAF[
	(pan_cohort_scMAF["snv_class"].isin(["synonymous_variant"]))
].shape[0] / panel_wide_mutation_type_counts[["Synonymous"]].sum()

print(f"Panel-wide dN/dS: {panel_wide_dN/panel_wide_dS}")

# %%