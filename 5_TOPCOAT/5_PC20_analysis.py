# %% Import 
from tkinter import font
from networkx import constraint
import pandas as pd
import numpy as np
import mosaic.io as mio
import plotly.express as px
from tea.plots import plot_snv_clone
from tea.utils import sort_for_var
from pathlib import Path
import os
from tea.utils import rgb_string_to_hex

os.chdir(Path(__file__).parent)
patient_name = "PC20"

pt_h5 = list(Path("../data_compiled/fillout_h5").glob(f"{patient_name}*.h5"))[0]
pt = mio.load(pt_h5)

snv_f = list(Path("../data_compiled/manual_annotated_snv_lists").glob(f"{patient_name}*voi*.txt"))[0]
output_dir = Path(".")
condor_downstream_dir = Path("../0_condor_pipeline/condor_downstream/ete_trees_refined_subclonal_snvs")

amplicon_params = Path('/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_project/tap_cn_calling/train-normals/train-combined_8_normals/NB_train-combined_8_normals-results.gene_added.csv')
amp_params_df = pd.read_csv(
    amplicon_params,
    index_col= 0,
)

COLOR_SEQUENCE= [rgb_string_to_hex(x) for x in px.colors.qualitative.Pastel]

# %% ===== Add raw FALCON results to H5 =====
cn_assignment_f = list(condor_downstream_dir.glob(f"{patient_name}/{patient_name}*assignment*csv"))[0]
cn_assignment_df = pd.read_csv(cn_assignment_f, index_col = 0)
print(f'[INFO] Loaded CN clone assignment file {cn_assignment_f}.')
# add cn_clone info
cn_assignment_df['cell_barcode_formatted'] = cn_assignment_df['cell_barcode'] + "-" + cn_assignment_df.index
cn_assignment_df.set_index('cell_barcode_formatted', inplace=True)

# assign cells from pt that are not in cn_assignment_df to clone 0
cn_assignment_df = cn_assignment_df.reindex(pt.dna.barcodes(), fill_value=0)


cn_clone_palette = dict(zip(
    np.sort(cn_assignment_df['final_clone_id'].unique()), 
    np.array(COLOR_SEQUENCE)[np.sort(cn_assignment_df['final_clone_id'].unique())]
    ))
# rename the keys
cn_clone_palette = {f"CN_clone-{k}": v for k, v in cn_clone_palette.items()}

pt.dna.row_attrs['label'] = np.array(list(
    map(
        lambda x: f"CN_clone-{int(cn_assignment_df.loc[x, 'final_clone_id'])}", 
        pt.dna.barcodes())
    ))
pt.dna.set_palette(cn_clone_palette)
pt.cnv.row_attrs['label'] = pt.dna.row_attrs['label']
pt.cnv.set_palette(cn_clone_palette)
pt.cnv.get_gene_names(amplicon_params, gene_name_col = 'gene')
pt.cnv.var = pd.DataFrame.from_dict(pt.cnv.col_attrs).set_index("id")
raw_rc = pt.cnv.get_attribute('read_counts', constraint='row')

normalized_rc = raw_rc / amp_params_df.loc[raw_rc.columns]['amplicon_factor'] / raw_rc.sum(axis=1).values[:, None]
# add to sample
pt.cnv.add_layer('normalized_read_counts', normalized_rc.values)

normalized_rc_binary = (normalized_rc == 0).astype(int)
pt.cnv.add_layer('normalized_read_counts_binary[zero/nonzero]', normalized_rc_binary.values)


# %% read in high-quality SNVs
snv_df = pd.read_csv(snv_f, sep = '\t', index_col = 0, comment='#')
snv_df["annotation"].fillna("", inplace=True)
# filter out artifact
snv_df = snv_df[~snv_df['annotation'].str.contains("artifact")]
# filter out germline HOM
snv_df = snv_df[~snv_df['annotation'].str.contains("germline_HOM")]
snv_df.sort_values(by=["mut_prev", "HGVSp"], ascending=False)
snv_ann_map = snv_df['HGVSp'].to_dict()

# # BPA-2
# GOI=["INO80", "CDKN2A", "RREB1", "TP53", "PBRM1", "OGG1", "BRCA2", "ASCC3", "PER1", "KMT2D", "RXRA", "CHEK2", "DNMT3B", "RTEL1", "TGFBR2", "NOTCH2", "CHD3", "POLE", "KRAS"]

# voi = snv_df[snv_df['HGVSp'].str.split().str[0].isin(GOI)].index.tolist()
# for RA17_13, we take everything
voi = snv_df.index.tolist()

# ===== highlight vars with bulk annotation ====='
# highlight vars
germline_hom_var_col = '#00cc66' # green
germline_het_var_col = '#2693ff' # blue
somatic_var_col = '#ff0000' # red
# highlight_col = '#e100ff' # purple

for var_i in snv_ann_map:
    # germline
    if snv_df.loc[var_i, 'annotation'] == 'germline_HOM':
        snv_ann_map[var_i] = f'<span style="color:{germline_hom_var_col};">' + snv_ann_map[var_i] + '</span>'
    elif snv_df.loc[var_i, 'annotation'] == 'germline_HET':
        snv_ann_map[var_i] = f'<span style="color:{germline_het_var_col};">' + snv_ann_map[var_i] + '</span>'
    elif "somatic" in snv_df.loc[var_i, 'annotation']:
        snv_ann_map[var_i] = f'<span style="color:{somatic_var_col};">' + snv_ann_map[var_i] + '</span>'
    elif "artifact" in snv_df.loc[var_i, 'annotation']:
        snv_ann_map[var_i] = f'<span style="color:{highlight_col};">' + snv_ann_map[var_i] + '</span>'
    else:
        pass

# %% ===== too many normal cells, sample to 1/3 of other cells =====
normal_cells = pt.dna.barcodes()[(pt.dna.get_attribute("label") == "CN_clone-0")[0]]
other_cells = pt.dna.barcodes()[(pt.dna.get_attribute("label") != "CN_clone-0")[0]]
np.random.seed(0)
sampled_normal_cells = np.random.choice(normal_cells, size = len(other_cells) // 3, replace = False)
sampled_cells = np.concatenate([sampled_normal_cells, other_cells])
pt = pt[sampled_cells]
# %% ===== plot DNA heatmap =====
fig = plot_snv_clone(
    pt,
    sample_name=patient_name,
    story_topic = f'{patient_name}-high_conf_snvs',
    voi = voi,
    attribute = "AF_MISSING",
    barcode_sort_method = "stringsort",
    ann_map = snv_ann_map
)
# # map the x-axis ticks according to snv_ann_map
# fig.update_xaxes(ticktext = list(map(lambda x: snv_ann_map[x], snv_df.index.tolist())))

# save the figure
# write to PDF, with width proportion to the number of SNVs
fig.write_image(
    str(output_dir / f"{patient_name}_sc_heatmap_DNA.pdf"), 
    width = 500 + 10 * len(snv_df.index.tolist())
)
print(f'[INFO] Saved heatmap for {patient_name} to {output_dir / f"{patient_name}-DNA-heatmap.pdf"}.')
# %% get sorted barcodes to unify barcode order between DNA and CNV
attribute="AF_MISSING"
# also make sure to retrieve the sorted_bars 
# del pt.dna.__dict__["_Assay__heatmap"]
sorted_bars = sort_for_var(
    dna = pt.dna,
    vars = voi,
    attribute = attribute,
    method = "stringsort"
    )


# %% ===== CNV heatmap =====
sc_cnv_heatmap = pt.cnv.heatmap(
    'normalized_read_counts_binary[zero/nonzero]',
    features = ["TGFBR2", "CDKN2A", "KRAS", "TP53", "SMAD4"],
    bars_order = sorted_bars,
)

# update colorscale title to be: "homdel"
# update colorscale to be two discrete values: 0, 1
# decrease the height of the colorscale
sc_cnv_heatmap.layout.coloraxis.colorbar.title = 'zero-reads'
sc_cnv_heatmap.layout.coloraxis.colorbar.len = 0.2
sc_cnv_heatmap.layout.coloraxis.colorbar.tickvals = [0, 1]


sc_cnv_heatmap.layout.coloraxis.colorscale = [
    [0, 'rgb(0, 0, 0)'], 
    [1, 'rgb(255, 225, 0)']
]
sc_cnv_heatmap.update_layout(font_family="Arial")
sc_cnv_heatmap.show()
# %% Save sc_CNV heatmap
# sc_dna_heatmap.write_image(
#     output_dir / f'{patient_name}-snv_heatmap-AF_MISSING.png',
#     height=800, width=300, scale=3,
#     )

sc_cnv_heatmap.write_image(
    output_dir / f"{patient_name}_sc_heatmap_cnv_binarized_select_genes.pdf",
    width=500
    )

# %% ===== for each single sample, plot the DNA heatmap =====
sample_names = pd.Series(pt.dna.barcodes()).str.split("-", n=1).str[1].unique()
tumor_cells = pt.dna.barcodes()[(pt.dna.get_attribute("label") != "CN_clone-0")[0]]

ss_tumors = {}
for sample_i in sample_names:
    sample_cells = pt.dna.barcodes()[pd.Series(pt.dna.barcodes()).str.contains(sample_i)]
    sample_tumor_cells = list(set(tumor_cells) & set(sample_cells))
    ss_tumors[sample_i] = pt[sample_tumor_cells]

    # plot DNA heatmap
    fig = plot_snv_clone(
        ss_tumors[sample_i],
        sample_name=sample_i,
        story_topic = f'{patient_name}-high_conf_snvs',
        voi = voi,
        attribute = "AF_MISSING",
        barcode_sort_method = "stringsort",
        ann_map = snv_ann_map
    )
    # save the figure
    fig.write_image(
        str(output_dir / f"{sample_i}_sc_heatmap_DNA.pdf"), 
        width = 500 + 10 * len(snv_df.index.tolist())
    )
    print(f'[INFO] Saved heatmap for {sample_i} to {output_dir / f"{sample_i}_DNA_heatmap_tumor_cells.pdf"}.')
# %% ===== difference across timepoints seem to be: CHEK2 and MYC loci =====
chek2 = "chr22:29130300:C/T"
myc = "chr8:128750707:T/TC"
kras = "chr12:25380277:G/T"
# distribution of myc in those having positive read in T1
myc_af_t1 = ss_tumors["PC20-IR"].dna.get_attribute("AF_MISSING", features = [myc], constraint="row")
myc_af_t1 = myc_af_t1[myc_af_t1[myc] > 0]

import seaborn as sns
from tea.plots import mpl_config_params
import matplotlib.pyplot as plt

mpl_config_params(font_size = 12)
fig, ax = plt.subplots()

sns.histplot(
    myc_af_t1[myc], bins = 50,
    ax=ax, color="#ff0000"
    )
ax.set_xlabel("MYC sc-VAF at T1")
ax.set_ylabel("Frequency")
fig.savefig(output_dir / "PC20-IR_MYC_AF_distribution.pdf")

# %% ----- distribution of CHEK2 locus AF in all tumor cells
chek2_af_t0 = ss_tumors["PC20-RSX"].dna.get_attribute("AF_MISSING", features = [chek2], constraint="row")
chek2_af_t1 = ss_tumors["PC20-IR"].dna.get_attribute("AF_MISSING", features = [chek2], constraint="row")

mpl_config_params(font_size = 12)
fig, ax = plt.subplots(1,2, figsize = (8,5))
sns.histplot(
    chek2_af_t0[chek2], bins = 50, stat = "percent",
    ax=ax[0], multiple="dodge", color="#2693ff"
    )
sns.histplot(
    chek2_af_t1[chek2], bins = 50, stat = "percent",
    ax=ax[1], multiple="dodge", color="#2693ff"
    )
# add legend
# ax.legend(["T0", "T1"])
ax[0].set_ylabel("Frequency (%)")
ax[1].set_ylabel("")
for ax_i in ax:
    ax_i.set_xlabel("")
    ax_i.set_xlim(-10, 100)
    ax_i.set_ylim(0, 100)
# set figure title, below the x-axis
fig.text(0.5, 0.01, "CHEK2 sc-VAF at T0 (left), T1 (right)", ha="center")
# fig.suptitle("CHEK2 AF distribution", y=0)

fig.savefig(output_dir / "PC20_two_timepts_CHEK2_AF_distribution.pdf")
plt.close()
# %% ===== what's myc's sc-prev in T0? =====
ss_tumors["PC20-RSX"].dna.genotype_variants()
myc_ALT = ss_tumors["PC20-RSX"].dna.get_attribute("alt_read_count", features=[myc], constraint="row")

ss_tumors["PC20-IR"].dna.genotype_variants()
kras_myc_mut_genotype = ss_tumors["PC20-IR"].dna.get_attribute("mut_filtered", features=[kras, myc], constraint="row")
kras_mutant_cells = kras_myc_mut_genotype[kras].sum()
kras_myc_mutant_cells = (kras_myc_mut_genotype.sum(axis=1) == 2).sum()
t1_myc_ccf = (kras_myc_mut_genotype.sum(axis=1) == 2).sum() / kras_myc_mut_genotype[kras].sum()
# %%
