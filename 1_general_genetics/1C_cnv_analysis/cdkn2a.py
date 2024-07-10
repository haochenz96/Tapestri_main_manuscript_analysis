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

os.chdir(Path(os.path.realpath(__file__)).parent)

# need to use the params for Tapestri V3
amplicon_params = Path('/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_project/tap_cn_calling/train-normals/train-combined_8_normals/NB_train-combined_8_normals-results.gene_added.csv')
amp_params_df = pd.read_csv(
    amplicon_params,
    index_col= 0,
)

COLOR_SEQUENCE= [rgb_string_to_hex(x) for x in px.colors.qualitative.Pastel]

sample_sheet = "../Tapestri_all_patient_sample_map.yaml"
with open(sample_sheet, "r") as f:
    patient_sample_map = yaml.load(f, Loader=yaml.FullLoader)
patients = patient_sample_map.keys()

for patient_name in patients:
    pt_h5 = list(Path("data_compiled/fillout_h5").glob(f"{patient_name}*.h5"))
    if len(pt_h5) != 1:
        raise ValueError(f"[INFO] Found {len(pt_h5)} h5 files for {patient_name}. Expected 1.")
    pt = mio.load(pt_h5)

    snv_f = list(Path("data_compiled/manual_annotated_snv_lists").glob(f"{patient_name}*voi*.txt"))[0]
    output_dir = Path(".")
    condor_downstream_dir = Path("0_condor_pipeline/condor_downstream/ete_trees_refined_subclonal_snvs")
    # %% Add raw FALCON results to H5
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

# RA21_17
GOI=["KRAS", "TP53", "ARID2", "CTNNB1", "HERC2"]

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
	width = 100 + 10 * len(snv_df.index.tolist())
)
print(f'[INFO] Saved heatmap for {patient_name} to {output_dir / f"{patient_name}-DNA-heatmap.pdf"}.')
# %%
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
    features = ["TGFBR2", "CDKN2A", "KRAS", "BRCA2", "TP53", "SMAD4"],
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


# %% Plot ordered AF distribution for the two arid2 mutations

# arid2_snvs = ['chr12:46245592:C/G', 'chr12:46243885:C/G']
arid2_snvs = snv_df[snv_df["HGVSp"].str.contains("ARID2")].index.tolist()
arid2_snvs_af = pt.dna.get_attribute("AF_MISSING", features = arid2_snvs, constraint="row")
arid2_snvs_af.index.name = "cell_barcode_formatted"
# hierarchically sort the dataframe
arid2_snvs_af = arid2_snvs_af.sort_values(by=arid2_snvs_af.columns.tolist(), ascending=False)

# wide to long, preserve the original index
arid2_snvs_af_long = arid2_snvs_af.reset_index().melt(id_vars = "cell_barcode_formatted", var_name = "SNV", value_name = "AF_MISSING")

# make the scatterplot
fig = px.scatter(
    arid2_snvs_af_long, x="cell_barcode_formatted", y="AF_MISSING", color='SNV',
    )

# scale down the dot size
scatter_2d = px.scatter(arid2_snvs_af, x = arid2_snvs[0], y = arid2_snvs[1])
scatter_2d.update_traces(marker=dict(size=2))
# set x and y ranges to [0,100], fix the ratio to 1:1
scatter_2d.update_xaxes(range=[0, 100])
scatter_2d.update_yaxes(range=[0, 100])
scatter_2d.update_layout(width=800, height=800, autosize=False, margin=dict(l=0, r=0, b=0, t=0))

# # lay over a contour plot to indicate density of the dots in the scatterplot above
# scatter_2d.add_contour(
#     x = arid2_snvs_af[arid2_snvs[0]],
#     y = arid2_snvs_af[[1]],
#     contours_coloring='fill',
#     line_width=0,
#     opacity=0.5,
#     contours=dict(
#         start=0,
#         end=arid2_snvs_af[arid2_snvs[0]].max(),
#         size=arid2_snvs_af[arid2_snvs[0]].max()/10
#     )
# )

arid2_snvs_af_positive = arid2_snvs_af[(arid2_snvs_af > 0).any(axis=1)]

contour_2d = px.density_contour(arid2_snvs_af_positive, x = arid2_snvs[0], y = arid2_snvs[1], nbinsx=10, nbinsy=10)
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

contour_2d.update_xaxes(dtick=10, range=[0, 100], title = snv_ann_map[arid2_snvs[0]])
contour_2d.update_yaxes(dtick=10, range=[0, 100], title = snv_ann_map[arid2_snvs[1]])
contour_2d.update_layout(
    title = f"{patient_name} - {snv_ann_map['chr12:46245592:C/G']} vs {snv_ann_map['chr12:46243885:C/G']}",
    #yaxis=dict(range=[-1.75, 1.5]),
    height = 800,
    width =800,
    bargap = 0,
    #hovermode = 'closest',
    title_text='DENSITY PLOT',
    showlegend = True
    )
contour_2d.show()

# save
contour_2d.write_image(f"{patient_name}-ARID2-density-contour.pdf")


# %% Density heatmap
density_heatmap = px.density_heatmap(
    arid2_snvs_af_positive, 
    x = arid2_snvs[0], y = arid2_snvs[1],
    nbinsx=20, nbinsy=20 
    )
density_heatmap.update_xaxes(dtick=10, range=[-10, 105], title = snv_ann_map[arid2_snvs[0]])
density_heatmap.update_yaxes(dtick=10, range=[-10, 105], title = snv_ann_map[arid2_snvs[1]])
density_heatmap.update_layout(width=400, height=400, autosize=False)
density_heatmap.write_image(f"{patient_name}-ARID2-density-heatmap.pdf")

# %%
