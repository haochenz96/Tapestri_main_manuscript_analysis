# %% Import 
from tkinter import font
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
patient_name = "PC10"

pt_h5 = list(Path("../data_compiled/fillout_h5").glob(f"{patient_name}*.h5"))[0]
pt = mio.load(pt_h5)

snv_f = list(Path("../data_compiled/manual_annotated_snv_lists").glob(f"{patient_name}*voi*.txt"))[0]
output_dir = Path(".")
condor_downstream_dir = Path("../0_condor_pipeline/condor_downstream")

amplicon_params = Path('/Users/haochenzhang/Iacobuzio_lab/Tapestri_batch2/tap_cn_calling/train-normals/train-combined_8_normals/NB_train-combined_8_normals-results.gene_added.csv')
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

# %% ===== plot heatmap =====
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

# %% ==== Simple 2D histogram of AF distribution for the two TGFBR1 mutations =====
tgfbr1_snvs = ["chr9:101891277:C/T", "chr9:101900209:C/T"]
tgfbr1_snvs_af = pt.dna.get_attribute("AF_MISSING", features = tgfbr1_snvs, constraint="row")
tgfbr1_snvs_af.index.name = "cell_barcode_formatted"

tgfbr1_snvs_af_positive = tgfbr1_snvs_af[(tgfbr1_snvs_af > 20).any(axis=1)]

tgfbr1_snvs_af_positive_long = tgfbr1_snvs_af_positive.reset_index().melt(id_vars = "cell_barcode_formatted", var_name = "SNV", value_name = "AF_MISSING")

# ===== scatter
# scale down the dot size
scatter_2d = px.scatter(tgfbr1_snvs_af_positive, x = tgfbr1_snvs[0], y = tgfbr1_snvs[1])
scatter_2d.update_traces(marker=dict(size=2))
# set x and y ranges to [0,100], fix the ratio to 1:1
scatter_2d.update_xaxes(range=[0, 100])
scatter_2d.update_yaxes(range=[0, 100])
scatter_2d.update_layout(width=800, height=800, autosize=False, margin=dict(l=0, r=0, b=0, t=0))

# ===== density
contour_2d = px.density_contour(
    tgfbr1_snvs_af_positive, 
    x = tgfbr1_snvs[0], y = tgfbr1_snvs[1], 
    nbinsx=10, nbinsy=10
    )
contour_2d.update_layout(
    font_family="Arial",
    # title_font_color="red",
    # legend_title_font_color="green"
)
contour_2d.update_traces(
    contours_showlabels = True,
    contours_coloring="fill",
    colorbar=dict(
        len=0.3,
        tickvals=[],
        title='Number of single nuclei', 
        titlefont=dict(
        size=10,
        ))
    )

contour_2d.update_xaxes(dtick=10, range=[0, 100], title = snv_ann_map[tgfbr1_snvs[0]])
contour_2d.update_yaxes(dtick=10, range=[0, 100], title = snv_ann_map[tgfbr1_snvs[1]])
contour_2d.update_layout(
    title = f"{patient_name} - {snv_ann_map['chr9:101891277:C/T']} vs {snv_ann_map['chr9:101900209:C/T']}",
    #yaxis=dict(range=[-1.75, 1.5]),
    height = 800,
    width =800,
    bargap = 0,
    #hovermode = 'closest',
    title_text='DENSITY PLOT',
    showlegend = True
    )
contour_2d.show()

# save the contour plot
contour_2d.write_image(f"{patient_name}-TGFBR1-density.pdf")

# %% Complex 2D histogram of above
import plotly.graph_objects as go

contour_2d_complex = go.Figure()
contour_2d_complex.add_trace(go.Histogram2dContour(
    x = tgfbr1_snvs_af_positive[tgfbr1_snvs[0]], 
    y = tgfbr1_snvs_af_positive[tgfbr1_snvs[1]],
    colorscale = 'Hot',
    reversescale = True,
    xaxis = 'x',
    yaxis = 'y',
    colorbar=dict(
    len=0.3,
    tickvals=[],
    title='Number of single nuclei', 
    titlefont=dict(
    size=10,
    ))
))
contour_2d_complex.add_trace(go.Scatter(
    x = tgfbr1_snvs_af_positive[tgfbr1_snvs[0]], 
    y = tgfbr1_snvs_af_positive[tgfbr1_snvs[1]],
    xaxis = 'x',
    yaxis = 'y',
    mode = 'markers',
    marker = dict(
        color = 'rgba(0,0,0,0.2)',
        size = 3
    )
))
contour_2d_complex.add_trace(go.Histogram(
    y = tgfbr1_snvs_af_positive[tgfbr1_snvs[0]],
    xaxis = 'x2',
    marker = dict(
        color = 'rgba(0,0,0,1)'
    )
))
contour_2d_complex.add_trace(go.Histogram(
    x = tgfbr1_snvs_af_positive[tgfbr1_snvs[1]],
    yaxis = 'y2',
    marker = dict(
        color = 'rgba(0,0,0,1)'
    )
))

contour_2d_complex.update_layout(
    autosize = False,
    xaxis = dict(
        zeroline = True,
        domain = [0,0.85],
        showgrid = False,
        range=[-1, 100],
        tick0 = 0, dtick=10, 
    ),
    yaxis = dict(
        zeroline = True,
        domain = [0,0.85],
        showgrid = False,
        range=[-1, 100],
        tick0 = 0, dtick=10, 
    ),
    xaxis2 = dict(
        zeroline = False,
        domain = [0.85,1],
        showgrid = False,
        showticklabels = False,
    ),
    yaxis2 = dict(
        zeroline = False,
        domain = [0.85,1],
        showgrid = False,
        showticklabels = False,
    ),
    height = 600,
    width = 600,
    bargap = 0,
    hovermode = 'closest',
    showlegend = False
)
contour_2d_complex.update_layout(
    # title="Plot Title",
    xaxis_title=snv_df.loc[tgfbr1_snvs[0], "HGVSp"],
    yaxis_title=snv_df.loc[tgfbr1_snvs[1], "HGVSp"],
    legend_title="Legend Title",
    font=dict(
        family="Arial",
        size=12,
    ),
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
)
contour_2d_complex.show()

# %% Save complex 2D histogram
contour_2d_complex.write_image(f"{patient_name}_2xTGFBR1_sc_contour_complex.pdf")
# %% ===== amplifications analysis =====

myc_amps = ["AMPL126740","AMPL143573","AMPL257985","AMPL257986","AMPL257987","AMPL257988"]
mtor_amps = ["AMPL257566","AMPL253201","AMPL257567","AMPL250219","AMPL248729"]
gata6_amps = ["AMPL257559","AMPL257560","AMPL257561","AMPL167407","AMPL167408","AMPL167408","AMPL167409"]

# %% Plot the raw read counts for MYC, MTOR, GATA6
import matplotlib.pyplot as plt
import seaborn as sns
from tea.plots import mpl_config_params
# use seaborn, for each gene (myc, mtor, gata6), for each clone, plot the normalized read counts (row-wise normalized)

raw_rc_df = pt.cnv.get_attribute(
    "normalized_read_counts", constraint="row",
    features = myc_amps + mtor_amps + gata6_amps
).T
raw_rc_df["gene"] = ["MYC"] * len(myc_amps) + ["MTOR"] * len(mtor_amps) + ["GATA6"] * len(gata6_amps)
# take the mean of the normalized read counts for each gene
raw_rc_df = raw_rc_df.groupby("gene").mean().T
raw_rc_df["clone_id"] = pt.cnv.row_attrs["label"]

raw_rc_df_for_plotting = raw_rc_df.melt(id_vars=["clone_id"], var_name="gene", value_name="normalized_read_counts")
raw_rc_df_for_plotting["clone_id"] = raw_rc_df_for_plotting["clone_id"].str.split("-").str[-1].astype(int)
# sort
raw_rc_df_for_plotting = raw_rc_df_for_plotting.sort_values(by=["gene", "clone_id"])
cn_clone_palette_for_plotting = dict(zip(
    np.sort(raw_rc_df_for_plotting["clone_id"].unique()), 
    np.array(COLOR_SEQUENCE)[np.sort(raw_rc_df_for_plotting["clone_id"].unique())]
    ))

# sns stripplot
mpl_config_params()
# plot stripplots (1x3)
fig, ax = plt.subplots(1, 3, figsize=(12, 5))
sns.set_theme(style="whitegrid")
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', weight='bold')
sns.stripplot(
    data=raw_rc_df_for_plotting[raw_rc_df_for_plotting["gene"] == "MYC"],
    x="clone_id", y="normalized_read_counts", hue="clone_id", palette=cn_clone_palette_for_plotting,
    s=1, linewidth=0, alpha=0.8,
    ax=ax[0]
)
ax[0].set_title("MYC")
sns.stripplot(
    data=raw_rc_df_for_plotting[raw_rc_df_for_plotting["gene"] == "MTOR"],
    x="clone_id", y="normalized_read_counts", hue="clone_id", palette=cn_clone_palette_for_plotting,
    s=1, linewidth=0, alpha=0.8,
    ax=ax[1], 
)
ax[1].set_title("MTOR")
sns.stripplot(
    data=raw_rc_df_for_plotting[raw_rc_df_for_plotting["gene"] == "GATA6"],
    x="clone_id", y="normalized_read_counts", hue="clone_id", palette=cn_clone_palette_for_plotting,
    s=1, linewidth=0, alpha=0.8,
    ax=ax[2], 
)
ax[2].set_title("GATA6")

# disable all x ticks and labels
for ax_i in ax:
    # ax_i.set_xticks([])
    # ax_i.set_xticklabels([])
    ax_i.set_xlabel('')
    ax_i.set_ylabel('')
    # increase y tick font size
    ax_i.tick_params(axis='x', labelsize=12)
    ax_i.tick_params(axis='y', labelsize=12)
ax[1].set_xlabel(
    'Clone ID', fontsize=15, # move down a little bit
    labelpad=10
)
ax[0].set_ylabel('normalized read counts', fontsize=15)

# disable color legend for all subplots
for ax_i in ax:
    ax_i.get_legend().remove()
# add extra padding horizontally
fig.tight_layout(pad=2.0)

fig.savefig(f"{patient_name}_clone_raw_rc_stripplot_MYC_MTOR_GATA6.pdf")
# fig.savefig(f"{patient_name}-clone_raw_rc-stripplot-MYC_MTOR_GATA6.png", bbox_inches='tight', dpi=300)

# %% focus on MYC, show the lower read count region
mpl_config_params()
fig, ax = plt.subplots(figsize=(4, 5))
sns.stripplot(
    data=raw_rc_df_for_plotting[raw_rc_df_for_plotting["gene"] == "MYC"],
    x="clone_id", y="normalized_read_counts", hue="clone_id", palette=cn_clone_palette_for_plotting,
    s=1, linewidth=0, alpha=0.8,
    ax=ax
)
ax.set_ylim([-1, 20])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.set_title("MYC", fontsize = 15)
ax.set_xlabel(
    'Clone ID', fontsize=15, # move down a little bit
    labelpad=10
)
ax.set_ylabel('normalized read counts', fontsize=15)
ax.get_legend().remove()
fig.savefig(
    f"{patient_name}_clone_raw_rc_stripplot_MYC_low_rc_range.pdf",
    bbox_inches='tight',
)



# %%
