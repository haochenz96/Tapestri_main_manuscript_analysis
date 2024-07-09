# %%
import os
import logging
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

os.chdir(os.path.dirname(__file__))

case_name = "RA16_29"
unique_cn_clone_profiles_df = pd.read_csv(
    f"LWP_hatchet_results/{case_name}_hatchet_profiles.csv", 
    index_col=0
    )
# ----- adjust for WGD -----
unique_cn_clone_profiles_df = unique_cn_clone_profiles_df / 2

# ----- config -----
color_sequence = px.colors.qualitative.Pastel
amp_gene_map_f = "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_project/tap_cn_calling/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt"
amp_gene_map_df = pd.read_csv(amp_gene_map_f, sep="\t", index_col="amplicon_number")

# %%
# ----- plotting params -----
unique_amplicon_order = amp_gene_map_df.index
gene_names = amp_gene_map_df.loc[unique_amplicon_order, "gene_name"].values
gene_names_to_plot = pd.Series(gene_names).value_counts()[pd.Series(gene_names).value_counts() >= 3].index # plot gene names covered by at least 3 amplicons

# Draw cluster labels
# cluster_labels = np.arange(unique_cn_clone_profiles_df.shape[0])
cluster_labels = unique_cn_clone_profiles_df.index.values
cluster_labels = [i.split("_")[1] for i in cluster_labels]
cluster_labels = [int(i) for i in cluster_labels]
logging.info(f"identified {len(cluster_labels)} unique clones")

cn_clone_palette = dict(zip(cluster_labels, np.array(color_sequence)[cluster_labels]))
cluster_colors = [cn_clone_palette[i] for i in cluster_labels]
if len(cluster_labels) > len(color_sequence):
	logging.warning(f"more clusters than colors, some clusters will be colored the same")
elif len(cluster_labels) == 1: # @HZ 2024-04-26: if only one clster, go.heatmap will return error
	cluster_colors = None

###########  ----- CN cluster profiles -----  ################
# draw subplots
fig = make_subplots(
	rows=2, cols=2,
	shared_yaxes=True, shared_xaxes=True,
	horizontal_spacing=0.01,
	vertical_spacing=0.01,
	column_widths=[1 / 25, 24 / 25],
	row_heights=[1 / 25, 24 / 25],
	)
# embed()
# get labels
labs = go.Heatmap(
	z=cluster_labels,
	y=np.arange(unique_cn_clone_profiles_df.shape[0]),
	x=[0] * unique_cn_clone_profiles_df.shape[0],
	customdata=cluster_labels,
	colorscale=cluster_colors,
	hovertemplate='label: %{customdata}<extra></extra>',
	showlegend=False,
	showscale=False
	)
fig.add_trace(labs, row=2, col=1)

# @HZ 12/20/2022: both ticktext and tickvals are needed to manually draw the ticklabels, otherwise it won't show up
fig.layout.yaxis3.ticktext = cluster_labels
fig.layout.yaxis3.tickvals = np.arange(unique_cn_clone_profiles_df.shape[0])

# Draw gene names
# embed()
un_genes, idx_genes, inv_genes, cnts_genes = np.unique(gene_names, return_index=True, return_inverse=True, return_counts=True)
gene_col_binary = [0]
for i in np.arange(1,len(inv_genes)):
	# iterate through inv_genes, get connectivity of amplicons
	if inv_genes[i] == inv_genes[i-1]:
		gene_col_binary.append(gene_col_binary[i-1])
	else:
		gene_col_binary.append(abs(1-gene_col_binary[i-1]))

ticks = (idx_genes + cnts_genes / 2).astype(int)
gene_names_subplot = go.Heatmap(
	z = gene_col_binary,
	x = unique_cn_clone_profiles_df.columns,
	y = [0] * unique_cn_clone_profiles_df.shape[1],
	colorscale = [[0, 'rgb(0,0,0)'], [1, 'rgb(144,144,144)']],
	showlegend=False,
	showscale=False,
)
fig.add_trace(gene_names_subplot, row=1, col=2)

# Draw main heatmap
# labels = np.tile(labels[:, None], (1, unique_cn_clone_profiles_df.shape[1]))
vals = go.Heatmap(
	z=unique_cn_clone_profiles_df,
	y=np.arange(unique_cn_clone_profiles_df.shape[0]),
	x=unique_cn_clone_profiles_df.columns,
	# customdata=labels,
	coloraxis='coloraxis',
	hovertemplate='%{z:.2f}<br>%{x}<extra>%{customdata}</extra>',
	showlegend=False,
	showscale=False
)
fig.add_trace(vals, row=2, col=2)

# draw gene names
genes_ticks = (idx_genes + cnts_genes / 2).astype(int)
	
# embed()

fig.layout.xaxis2.ticktext = [ i if i in gene_names_to_plot else "" for i in gene_names[genes_ticks] ]
fig.layout.xaxis2.tickvals = genes_ticks
fig.layout.xaxis2.tickfont = {
	'size': 8,
}
fig.update_layout({'xaxis2': {'ticklen': 4, 'side': 'top', 'tickangle': -60, 'showticklabels': True}})

# draw chromosome numbers
chromosome_ordered = amp_gene_map_df.loc[unique_amplicon_order, "chr"].values
un, ind, cnts = np.unique(chromosome_ordered, return_index=True, return_counts=True)
ticks = (ind + cnts / 2).astype(int)

fig.layout.xaxis4.ticktext = chromosome_ordered[ticks]
fig.layout.xaxis4.tickvals = ticks
fig.update_layout({'xaxis4': {'tickangle': -45, 'showticklabels': True}})

for i in ind:
	fig.add_vline(i - 0.5, line_color='lightcyan', line_width=1, row=2, col=2)
	
######################################################
# update color schemes
num_vals = 7
colorscale = [
	(0, 'rgb(163, 163, 163)'), (1/num_vals, 'rgb(163, 163, 163)'), # NA
	(1/num_vals, 'rgb(0,0,0)'), (2/num_vals, 'rgb(0,0,0)'), # 0
	(2/num_vals, 'rgb(7, 95, 237)'), (3/num_vals, 'rgb(7, 95, 237)'), # 1
	(3/num_vals, 'rgb(146, 170, 209)'), (4/num_vals, 'rgb(146, 170, 209)'), # 2 
	(4/num_vals, 'rgb(237, 156, 57)'), (5/num_vals, 'rgb(237, 156, 57)'), # 3
	(5/num_vals, 'rgb(242, 29, 22)'), (6/num_vals, 'rgb(242, 29, 22)'), # 4
	(6/num_vals, 'rgb(202, 82, 250)'), (1, 'rgb(202, 82, 250)') # 5+
	]

colorbar_ticktext=[str(i) for i in list(range(num_vals-1))]
colorbar_ticktext.insert(0, 'NA')
colorbar_ticktext[-1] += '+'
colorbar_tickvals = [(num_vals-1)/(num_vals*2) * (2*i + 1) - 1 for i in range(num_vals)] 
fig.update_layout(
	coloraxis=dict(
		colorscale=colorscale,
		colorbar_tickvals = colorbar_tickvals,
		colorbar_ticktext = colorbar_ticktext,
		colorbar_title=dict(
			font_size = 10,
			text = 'total_copy_number',
		),
		cmax=num_vals-2,
		cmin=-1
	),
	font_family = 'Arial',
	width=1600,
    height=400,
	)

fig.write_image(
	f"LWP_hatchet_results/{case_name}_hatchet_profiles.pdf",
)

# %%
