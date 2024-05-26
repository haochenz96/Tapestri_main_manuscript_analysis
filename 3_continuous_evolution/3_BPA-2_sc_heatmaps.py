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
patient_name = "BPA-2"

pt_h5 = list(Path("../data_compiled/fillout_h5").glob(f"{patient_name}*.h5"))[0]
pt = mio.load(pt_h5)

snv_f = list(Path("../data_compiled/manual_annotated_snv_lists").glob(f"{patient_name}*voi*.txt"))[0]
output_dir = Path(".")
condor_downstream_dir = Path("../0_condor_pipeline/condor_downstream/ete_trees_refined_subclonal_snvs")

# need to use the params for Tapestri V3
amplicon_params = Path('/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_project/tap_cn_calling/train-normals/v3_chem_RA18_18-11_1/NB_train-v3_chem_RA18_18-11_1-results.csv')
amp_params_df = pd.read_csv(
    amplicon_params,
    index_col= 0,
)

COLOR_SEQUENCE= [rgb_string_to_hex(x) for x in px.colors.qualitative.Pastel]

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

# BPA-2
GOI=["INO80", "CDKN2A", "RREB1", "TP53", "PBRM1", "OGG1", "BRCA2", "ASCC3", "PER1", "KMT2D", "RXRA", "CHEK2", "DNMT3B", "RTEL1", "TGFBR2", "NOTCH2", "CHD3", "POLE", "KRAS"]

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

# %% Fix barcode names
pt.cnv.row_attrs["barcode"] = pt.cnv.row_attrs["barcode"] + "-" + pt.cnv.row_attrs["sample_name"]

# %% CNV heatmap
sc_cnv_heatmap = pt.cnv.heatmap(
    'normalized_read_counts_binary[zero/nonzero]',
    features = ["CDKN2A", "BRCA2"],
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

# %% ===== mutual colocalization test for BRCA2 homdel and other voi =====
from tea.snv.me import *

# calculate presence prob of TP53 HET and BRCA2 HOMDEL
pp = pt.dna.get_attribute("AF_MISSING", constraint = "row", features = ["chr17:7573996:A/T", "chr13:32914437:GT/G"])

# for TP53 HET, get binary indicator for AF within [20, 80]
p53_pp = ((pp["chr17:7573996:A/T"] >= 20) & (pp["chr17:7573996:A/T"] <= 80)).astype(int)
# for BRCA2, get binary indicator for AF = -50 (DP=0)
brca2_pp = (pp["chr13:32914437:GT/G"] == -50).astype(int)

me_inv_pval = __compute_mutual_correlation_pval(
    brca2_pp, p53_pp, 
    # null_prob_mut_1, null_prob_mut_2, 
    gametes = [(1,0), (0,1)],
    rm_irrelevant_cells = True
    )

print("BRCA2 HOMDEL and TP53 HET mutual colocalization p-value: ", me_inv_pval)

# %%
