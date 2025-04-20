# %%
import sys
import argparse
import yaml
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
patient_name = "PC11"

pt_h5 = f"../data_compiled/fillout_h5/{patient_name}.patient_wide.genotyped.h5"
pt = mio.load(pt_h5)

snv_f = f"../data_compiled/manual_annotated_snv_lists/{patient_name}-patient-all_vars-voi.hz_curated.txt"
output_dir = Path(".")
condor_downstream_dir = Path("../0_condor_pipeline/condor_downstream")

amplicon_params = Path('/Users/haochenzhang/Iacobuzio_lab/Tapestri_batch2/tap_cn_calling/train-normals/train-combined_8_normals/NB_train-combined_8_normals-results.gene_added.csv')
amp_params_df = pd.read_csv(
    amplicon_params,
    index_col= 0,
)

COLOR_SEQUENCE= [rgb_string_to_hex(x) for x in px.colors.qualitative.Pastel]

# %% ===== Add final CN clone assignment to H5 =====
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

# RA15_06
GOI=["PBRM1", "TP53", "ARID1B", "WRN", "TGFBR2", "RTEL1", "ARID1A", "SMAD4", "KRAS", "PARP4"]
# RA17_22
GOI=["CHEK2", "PTPRS", "PARP4", "TP53", "AXIN1", "ACVR2A", "PBRM1", "TGFBR2", "RBM10", "TTK", "KRAS", "SMAD4", "PIK3CA"]
voi = snv_df[snv_df['HGVSp'].str.split().str[0].isin(GOI)].index.tolist()

# ===== highlight vars with bulk annotation ====='
# highlight vars
germline_hom_var_col = '#00cc66' # green
germline_het_var_col = '#2693ff' # blue
somatic_var_col = '#ff0000' # red
# highlight_col = '#e100ff' # yellow

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

# plot heatmap
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
	width = 10 * len(snv_df.index.tolist())
)
print(f'[INFO] Saved heatmap for {patient_name} to {output_dir / f"{patient_name}-DNA-heatmap.pdf"}.')
# %% sort barcodes
attribute="AF_MISSING"
# also make sure to retrieve the sorted_bars 
# del pt.dna.__dict__["_Assay__heatmap"]
sorted_bars = sort_for_var(
    dna = pt.dna,
    vars = voi,
    attribute = attribute,
    method = "stringsort"
    )


# %% CNV heatmap
sc_cnv_heatmap = pt.cnv.heatmap(
    'normalized_read_counts_binary[zero/nonzero]',
    features = ["KRAS", "SMAD4"],
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
# %% Save figures
# sc_dna_heatmap.write_image(
#     output_dir / f'{patient_name}-snv_heatmap-AF_MISSING.png',
#     height=800, width=300, scale=3,
#     )


sc_cnv_heatmap.write_image(
    output_dir / f"{patient_name}_sc_heatmap_cnv_binarized_select_genes.pdf",
    width=500
    )

# %% KRAS allelic imblanace in RA17_22
# plot distribution of AF at KRAS locus (chr12:25398284:C/A), split by label
kras_af = pt.dna.get_attribute("AF", constraint = "row", features = ["chr12:25398284:C/A"])
cell_labels = pt.dna.get_attribute("label", constraint = "row")
kras_af["label"] = cell_labels.loc[kras_af.index]
# filter out clone 0
kras_af = kras_af[kras_af["label"] != "CN_clone-0"]
kras_af.columns = ["KRAS_AF", "CN_clone"]
kras_af.sort_values(by="CN_clone", inplace=True)
# seaborn violin plot
import seaborn as sns
import matplotlib.pyplot as plt
from tea.plots import mpl_config_params

# sns.violinplot(
#     data=kras_af, x="CN_clone", y="KRAS_AF",
#     legend=False, density_norm="count", inner=None,
#     # palette=RESPONSE_PALETTE, order=RESPONSE_PALETTE.keys(),
#     # ax=ax[(count%2), count//2]
# )
fig, ax = plt.subplots()
mpl_config_params()
sns.boxplot(
    data=kras_af, x="CN_clone", y="KRAS_AF",
    hue = "CN_clone",
    palette = cn_clone_palette,
)
ax.set_xlabel("")
ax.set_aspect(0.05)
fig.savefig(
    output_dir / f"{patient_name}_KRAS_AF_boxplot.pdf",
    bbox_inches="tight"
)
# %% ===== amplifications analysis =====

KRAS_amps = ["AMPL38584", "AMPL115373", "AMPL38586", "AMPL38587", "AMPL115376"]

# %% Plot the raw read counts for KRAS amplification
import matplotlib.pyplot as plt
import seaborn as sns
from tea.plots import mpl_config_params
# use seaborn, for each gene (myc, mtor, gata6), for each clone, plot the normalized read counts (row-wise normalized)

raw_rc_df = pt.cnv.get_attribute(
    "normalized_read_counts", constraint="row",
    features = KRAS_amps
).T
raw_rc_df["gene"] = ["KRAS"] * len(KRAS_amps)
# take the mean of the normalized read counts for each gene
raw_rc_df = raw_rc_df.groupby("gene").mean().T
raw_rc_df["clone_id"] = pt.cnv.row_attrs["label"]

raw_rc_df_for_plotting = raw_rc_df.melt(id_vars=["clone_id"], var_name="gene", value_name="normalized_read_counts")
# raw_rc_df_for_plotting["clone_id"] = raw_rc_df_for_plotting["clone_id"].str.split("-").str[-1].astype(int)
# sort
raw_rc_df_for_plotting = raw_rc_df_for_plotting.sort_values(by=["gene", "clone_id"])
# cn_clone_palette_for_plotting = dict(zip(
#     np.sort(raw_rc_df_for_plotting["clone_id"].unique()), 
#     np.array(COLOR_SEQUENCE)[np.sort(raw_rc_df_for_plotting["clone_id"].unique())]
#     ))

# sns stripplot
mpl_config_params()
fig, ax = plt.subplots(figsize=(4, 5))
sns.stripplot(
    data=raw_rc_df_for_plotting[raw_rc_df_for_plotting["gene"] == "KRAS"],
    x="clone_id", y="normalized_read_counts", hue="clone_id", palette=cn_clone_palette,
    s=1, linewidth=0, alpha=0.8,
    ax=ax
)
ax.set_ylim([-1, 8])
ax.tick_params(axis='x', labelsize=8)
ax.tick_params(axis='y', labelsize=12)
ax.set_title("KRAS", fontsize = 15)
ax.set_xlabel(
    'Clone ID', fontsize=15, # move down a little bit
    labelpad=10
)
ax.set_ylabel('normalized read counts', fontsize=15)
# ax.get_legend().remove()
fig.savefig(
    f"{patient_name}_clone_raw_rc_stripplot_MYC_low_rc_range.pdf",
    bbox_inches='tight',
)

# %%

