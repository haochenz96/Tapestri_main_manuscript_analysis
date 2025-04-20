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
pan_cohort_scMAF = pan_cohort_scMAF[pan_cohort_scMAF["patient_name"] != "RA16_08"]
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
DRIVERS=["KRAS", "TP53","CDKN2A", "SMAD4", "SMAD2", "SMAD3", "TGFBR1", "TGFBR2", "ACVR1B", "BMPR1A", "ARID1A", "ARID2", "BRCA2", "ATM", "BAP1", "PIK3CA", "FGFR1","RNF43", "POLD1", "IRF6", "GATA6", "MYC", "MTOR"]

pan_cohort_scMAF["functional?"] = ~pan_cohort_scMAF["snv_class"].isin(NONFUNC_SO)
pan_cohort_scMAF["driver?"] = pan_cohort_scMAF["gene_name"].isin(DRIVERS)
pan_cohort_scMAF["functional_driver?"] = (~pan_cohort_scMAF["snv_class"].isin(NONFUNC_SO)) & (pan_cohort_scMAF["gene_name"].isin(DRIVERS))
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


dnds_result = {}
DRIVERS_mutation_type_counts = pd.concat([L[gene_i].sum(axis=0) for gene_i in DRIVERS],axis=1).sum(axis=1)
panel_wide_mutation_type_counts = pd.concat([L[gene_i].sum(axis=0) for gene_i in L.keys()],axis=1).sum(axis=1)
NONCODING_SO = ["intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "splice_site_variant", "2kb_upstream_variant"]


# %% calculate DRIVERS dN/dS
N = pan_cohort_scMAF[
	(pan_cohort_scMAF["gene_name"].isin(DRIVERS)) &
	(~pan_cohort_scMAF["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
].shape[0]
DRIVER_dN = N / DRIVERS_mutation_type_counts[["Missense", "Nonsense"]].sum()
S = pan_cohort_scMAF[
	(pan_cohort_scMAF["gene_name"].isin(DRIVERS)) &
	(pan_cohort_scMAF["snv_class"].isin(["synonymous_variant"]))
].shape[0]
DRIVER_dS = S / DRIVERS_mutation_type_counts[["Synonymous"]].sum()

dnds_result["DRIVERS"] = {
	"N": N,
	"S": S,
	"dN/dS": DRIVER_dN / DRIVER_dS,
}


# calculate panel-wide dN/dS
N = pan_cohort_scMAF[
	(~pan_cohort_scMAF["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
].shape[0]
panel_wide_dN = N / panel_wide_mutation_type_counts[["Missense", "Nonsense"]].sum()
S = pan_cohort_scMAF[
	(pan_cohort_scMAF["snv_class"].isin(["synonymous_variant"]))
].shape[0]
panel_wide_dS = S / panel_wide_mutation_type_counts[["Synonymous"]].sum()

dnds_result["PANEL_WIDE"] = {
	"N": N,
	"S": S,
	"dN/dS": panel_wide_dN / panel_wide_dS,
}

# %% stratify by CCF
ccf_cutoff = 0.005
for sub_df, ccf_cutoff in zip([pan_cohort_scMAF[pan_cohort_scMAF["CCF"]<ccf_cutoff], pan_cohort_scMAF[pan_cohort_scMAF["CCF"]>=ccf_cutoff]], ["<0.005", ">=0.005"]):	
	N = sub_df[
		(~sub_df["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
	].shape[0]
	dN = N / panel_wide_mutation_type_counts[["Missense", "Nonsense"]].sum()
	S = sub_df[
		(sub_df["snv_class"].isin(["synonymous_variant"]))
	].shape[0]
	dS = S / panel_wide_mutation_type_counts[["Synonymous"]].sum()
	dnds_result[f"CCF{ccf_cutoff}"] = {
		"N": N,
		"S": S,
		"dN/dS": dN / dS,
	}


# %% stratify by sample collection time
patient_sheet = pd.read_excel("../Tapestri_batch2_samples_MASTER_INTERNAL.xlsx", sheet_name = "all_case_genetics", skiprows=1, header=0, index_col=0)
patient_sheet["collection_time"] = patient_sheet["collection_method"].copy().map({
      "resection": "local_disease",
      "autopsy": "metastatic_disease",
      "biopsy": "metastatic_disease",
})

pan_cohort_scMAF["sample_collection_time"] = pan_cohort_scMAF["patient_name"].map(patient_sheet["collection_time"])
for group_name, sub_df in pan_cohort_scMAF.groupby("sample_collection_time"):	
	N = sub_df[
		(~sub_df["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
	].shape[0]
	dN = N / panel_wide_mutation_type_counts[["Missense", "Nonsense"]].sum()
	S = sub_df[
		(sub_df["snv_class"].isin(["synonymous_variant"]))
	].shape[0]
	dS = S / panel_wide_mutation_type_counts[["Synonymous"]].sum()
	if dS == 0:
		dS = 1e-17
		print(f"[WARNING] No synonymous mutations for {group_name}")
	dnds_result[group_name] = {
		"N": N,
		"S": S,
		"dN/dS": dN / dS,
	}

# %% focus on tgfb
TGFB_genes = ["SMAD4", "SMAD2", "SMAD3", "TGFBR1", "TGFBR2", "ACVR1B", "BMPR1A", "ARID1A", "ARID2"]
TGFB_mutation_type_counts = pd.concat([L[gene_i].sum(axis=0) for gene_i in TGFB_genes],axis=1).sum(axis=1)
# calculate dN/dS for TGFB, stratified by sample collection time
for group_name, sub_df in pan_cohort_scMAF.groupby("sample_collection_time"):
	sub_df = sub_df[sub_df["gene_name"].isin(TGFB_genes)]
	N = sub_df[
		(~sub_df["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
	].shape[0]
	dN = N / TGFB_mutation_type_counts[["Missense", "Nonsense"]].sum()
	S = sub_df[
		(sub_df["snv_class"].isin(["synonymous_variant"]))
	].shape[0]
	dS = S / TGFB_mutation_type_counts[["Synonymous"]].sum()
	if dS == 0:
		dS = 1e-17
		print(f"[WARNING] No synonymous mutations for {group_name}")
	dnds_result[f"TGFB_{group_name}"] = {
		"N": N,
		"S": S,
		"dN/dS": dN / dS,
	}

# %%
MAPK_genes = ["AKT1", "NRAS", "PIK3CA", "FGFR3", "BRAF", "FGFR1", "KRAS", "ERBB3", "MAP2K1", "MAP2K4", "NF1", "ERBB2", "GNAS"]
MAPK_mutation_type_counts = pd.concat([L[gene_i].sum(axis=0) for gene_i in MAPK_genes],axis=1).sum(axis=1)
# calculate dN/dS for MAPK, stratified by sample collection time
for group_name, sub_df in pan_cohort_scMAF.groupby("sample_collection_time"):
	sub_df = sub_df[sub_df["gene_name"].isin(MAPK_genes)]
	N = sub_df[
		(~sub_df["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
	].shape[0]
	dN = N / MAPK_mutation_type_counts[["Missense", "Nonsense"]].sum()
	S = sub_df[
		(sub_df["snv_class"].isin(["synonymous_variant"]))
	].shape[0]
	dS = S / MAPK_mutation_type_counts[["Synonymous"]].sum()
	if dS == 0:
		dS = 1e-17
		print(f"[WARNING] No synonymous mutations for {group_name}")
	dnds_result[f"MAPK_{group_name}"] = {
		"N": N,
		"S": S,
		"dN/dS": dN / dS,
	}

# %% look at the MAPK mutations in more detail
pan_cohort_scMAF[
	(pan_cohort_scMAF["sample_collection_time"] == "local_disease") & 
	(pan_cohort_scMAF["gene_name"].isin(MAPK_genes)) &
	(~pan_cohort_scMAF["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
]

# %% What's going on with PC03
import mosaic.io as mio
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
pt.dna.genotype_variants()
# %% Venn diagram
from matplotlib_venn import venn2
VOI = ["chr12:25398284:C/A", "chr3:178936082:G/A", "chr3:178936040:A/G", "chr8:38314936:C/T"]
ANCHOR = "chr3:178936082:G/A"
OTHERS = [v for v in VOI if v != ANCHOR]
figs = {}
for v in OTHERS:
	kras_cells = set(pt.dna.barcodes()[pt.dna.get_attribute('mut_filtered', features=[ANCHOR])[0].astype(bool)])
	map2k4_cells = set(pt.dna.barcodes()[pt.dna.get_attribute('mut_filtered', features=[v])[0].astype(bool)])
	fig, ax = plt.subplots(figsize=(5,5))
	venn2(
		[set(kras_cells), set(map2k4_cells)], 
		('KRAS p.G12V', v),
		ax=ax
		)
	figs[v] = fig

# %% where was chr3:178936082:G/A located?
pd.Series(pt.dna.row_attrs["sample_name"][
	(pt.dna.get_attribute('mut_filtered', features=["chr3:178936040:A/G"])==1)[0]
	]).value_counts()
# %% what's going on with PC11?
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

# %% ME
ANCHOR = 
from tea.snv.me import get_exclusivity
E = 1 - get_exclusivity(
    pt.dna.get_attribute('DP', constraint='row', features=VOI),
    pt.dna.get_attribute('alt_read_count', constraint='row', features=VOI),
    rm_irrelevant_cells = True
)
import matplotlib.pyplot as plt
from tea.plots import *
mpl_config_params(font_size=8)
import seaborn as sns

# E.columns = [snv_ann_map[i] for i in E.columns]
# E.index = [snv_ann_map[i] for i in E.index]
fig, ax = plt.subplots(figsize=(5,5))         # Sample figsize in inches
sns.heatmap(
    E, annot=True, 
    cmap='Blues', fmt='.2f', # vmin=0, vmax=1,
    ax=ax
    )

# %% Venn diagram
from matplotlib_venn import venn2

kras_cells = set(pt.dna.barcodes()[pt.dna.get_attribute('mut_filtered', features=["chr12:25398284:C/A"])[0].astype(bool)])
map2k4_cells = set(pt.dna.barcodes()[pt.dna.get_attribute('mut_filtered', features=["chr17:12016559:C/T"])[0].astype(bool)])
fig, ax = plt.subplots(figsize=(5,5))
venn2(
    [set(kras_cells), set(map2k4_cells)], 
    ('KRAS p.G12V', 'MAP2K4 p.P232L'),
    ax=ax
    )


# %%
dna_damage_repair_genes = ["ERCC3", "BARD1", "FANCD2", "POLQ", "TOPBP1", "ATR", "RNF168", "POLN", "FANCE", "MGMT", "MUS81", "ATM", "BRCA2", "HERC2", "TP53BP1", "FANCI", "MPG", "PALB2", "FANCA", "BRCA1", "BRD4", "PARP4", "TP63", "RECQL5", "TP53"]
dna_damage_repair_mutation_type_counts = pd.concat([L[gene_i].sum(axis=0) for gene_i in dna_damage_repair_genes],axis=1).sum(axis=1)
# calculate dN/dS for dna_damage_repair, stratified by sample collection time
for group_name, sub_df in pan_cohort_scMAF.groupby("sample_collection_time"):
	sub_df = sub_df[sub_df["gene_name"].isin(dna_damage_repair_genes)]
	N = sub_df[
		(~sub_df["snv_class"].isin(set(NONFUNC_SO + NONCODING_SO)))
	].shape[0]
	dN = N / dna_damage_repair_mutation_type_counts[["Missense", "Nonsense"]].sum()
	S = sub_df[
		(sub_df["snv_class"].isin(["synonymous_variant"]))
	].shape[0]
	dS = S / dna_damage_repair_mutation_type_counts[["Synonymous"]].sum()
	if dS == 0:
		dS = 1e-17
		print(f"[WARNING] No synonymous mutations for {group_name}")
	dnds_result[f"DNA_DAMAGE_REPAIR_{group_name}"] = {
		"N": N,
		"S": S,
		"dN/dS": dN / dS,
	}

# %%
dnds_result_df = pd.DataFrame(dnds_result)
dnds_result_df.to_csv("dnds_result.csv", index=True)
# %%
