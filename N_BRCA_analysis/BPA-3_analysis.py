# %%
import pandas as pd
import numpy as np
from pathlib import Path
import mosaic.io as mio
import seaborn as sns
import plotly.express as px
from tea.plots import plot_snv_clone
import os

# change to script directory
# get script dir
os.chdir(Path(os.path.realpath(__file__)).parent)

data_dir = Path('/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled')
amplicon_params = Path('../../Tapestri_project/tap_cn_calling/train-normals/train-combined_8_normals/NB_train-combined_8_normals-results.gene_added.csv')

amp_params_df = pd.read_csv(
    amplicon_params,
    index_col= 0,
)


# %%
# sample_objs = {}

# for sample_i in sample_names:
#     h5_f = data_dir / f"{patient_name}-patient" / 'cn_clone_added_h5s' / f"{sample_i}_cn_clone_added.h5"
#     sample_objs[sample_i] = mio.load(h5_f)

patient_name = "BPA-3"
pt_h5 = data_dir / "fillout_h5" / f"{patient_name}.patient_wide.genotyped.h5"
# pt_h5 = "/lila/data/iacobuzc/haochen/Tapestri_batch2/pipeline_results/BPA-2/BPA-2-IR/OUTPUTS_from_mpileup/BPA-2-IR.mpileup.h5"
pt = mio.load(pt_h5)

obs_df = pd.DataFrame(index=pt.dna.barcodes())
obs_df["sample"] = obs_df.index.str.split("-",n=1).str[1]
obs_df["cell"] = obs_df.index.str.split("-",n=1).str[0]
pt.obs = obs_df
samples = pt.obs["sample"].unique()
sample_objs = {}
for sample_i in samples:
    sample_objs[sample_i] = pt[pt.obs["sample"] == sample_i]

# %% Add raw FALCON results to H5
# HZ refined
cn_assignment_f = "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/0_condor_pipeline/condor_downstream/ete_trees_refined_subclonal_snvs/BPA-3/BPA-3_final_sc_clone_assignment.csv"
cn_assignment_df = pd.read_csv(cn_assignment_f, index_col = 0)
print(f'[INFO] Loaded CN clone assignment file {cn_assignment_f}.')
# add cn_clone info
cn_assignment_df['cell_barcode_formatted'] = cn_assignment_df['cell_barcode'] + "-" + cn_assignment_df.index
cn_assignment_df.set_index('cell_barcode_formatted', inplace=True)

# assign cells from pt that are not in cn_assignment_df to clone 0
cn_assignment_df = cn_assignment_df.reindex(pt.dna.barcodes(), fill_value=0)


cn_clone_palette = dict(zip(
    np.sort(cn_assignment_df['final_clone_id'].unique()), 
    np.array(px.colors.qualitative.Set3)[np.sort(cn_assignment_df['final_clone_id'].unique())]
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
snv_f = data_dir / "manual_annotated_snv_lists" / f"{patient_name}-all_vars-voi.hz_curated.txt"
snv_df = pd.read_csv(snv_f, sep = '\t', index_col = 0, comment='#')
snv_df = snv_df[snv_df["annotation"].notna()]
snv_ann_map = snv_df['HGVSp'].to_dict()

# ===== highlight vars with bulk annotation ====='
# highlight vars
germline_hom_var_col = '#00cc66' # green
germline_het_var_col = '#2693ff' # blue
somatic_var_col = '#ff0000' # red
likely_artifact_col = '#ffcc00' # yellow
highlight_color = '#f700ff'

for var_i in snv_ann_map:
    # germline
    if var_i == "chr13:32914288:TAACC/T":
        snv_ann_map[var_i] = f'<span style="color:{highlight_color};">' + snv_ann_map[var_i] + '</span>'
        continue
    if snv_df.loc[var_i, 'annotation'] == 'germline_HOM':
        snv_ann_map[var_i] = f'<span style="color:{germline_hom_var_col};">' + snv_ann_map[var_i] + '</span>'
    elif snv_df.loc[var_i, 'annotation'] == 'germline_HET':
        snv_ann_map[var_i] = f'<span style="color:{germline_het_var_col};">' + snv_ann_map[var_i] + '</span>'
    elif "somatic" in snv_df.loc[var_i, 'annotation'].lower():
        snv_ann_map[var_i] = f'<span style="color:{somatic_var_col};">' + snv_ann_map[var_i] + '</span>'
    elif snv_df.loc[var_i, 'annotation'] == 'likely_artifact':
        snv_ann_map[var_i] = f'<span style="color:{likely_artifact_col};">' + snv_ann_map[var_i] + '</span>'
    else:
        pass

# %% plot heatmap
# del pt.dna.__dict__["_Assay__heatmap"]
fig = plot_snv_clone(
    pt,
    sample_name=patient_name,
    story_topic = f'{patient_name}-high_conf_snvs',
    voi = snv_df.index.tolist(),
    attribute = "AF_MISSING",
    ann_map = snv_ann_map
)
# # map the x-axis ticks according to snv_ann_map
# fig.update_xaxes(ticktext = list(map(lambda x: snv_ann_map[x], snv_df.index.tolist())))

# %% save the figure
# write to PDF, with width proportion to the number of SNVs
fig.write_image(
    f"{patient_name}-DNA-heatmap.pdf", 
    width = 500 + 10 * len(snv_df.index.tolist()),
    engine="kaleido"
    )
# %%
from tea.utils import sort_for_var
attribute="AF_MISSING"
# also make sure to retrieve the sorted_bars 
del pt.dna.__dict__["_Assay__heatmap"]
sorted_bars = sort_for_var(
    dna = pt.dna,
    vars = snv_df.index.tolist(),
    attribute = attribute,
    method = "hier"
    )


# %% CNV heatmap
sc_cnv_heatmap = pt.cnv.heatmap(
    'normalized_read_counts_binary[zero/nonzero]',
    features = ['CDKN2A', 'KRAS', 'BRCA2', 'TP53', 'SMAD4'],
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
output_dir = Path(patient_name)
output_dir.mkdir(exist_ok=True, parents=True)

# sc_dna_heatmap.write_image(
#     output_dir / f'{patient_name}-snv_heatmap-AF_MISSING.png',
#     height=800, width=300, scale=3,
#     )


sc_cnv_heatmap.write_image(
    output_dir / f'{patient_name}-cnv_heatmap-binarized_select_genes.png',
    height=800, width=800, scale=3,
    )

# %% What's going on with normal pancreas' BRCA2?
n_f = "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/batch2_data_compiled/fillout_h5/normal_pancreas/RA18_10-30_1.mpileup.h5"
# ptn = mio.load(n_f)
ptn.cnv.get_gene_names(amplicon_params, gene_name_col = 'gene')
ptn.cnv.var = pd.DataFrame.from_dict(ptn.cnv.col_attrs).set_index("id")
raw_rc = ptn.cnv.get_attribute('read_counts', constraint='row')

normalized_rc = raw_rc / amp_params_df.loc[raw_rc.columns]['amplicon_factor'] / raw_rc.sum(axis=1).values[:, None]
# add to sample
ptn.cnv.add_layer('normalized_read_counts', normalized_rc.values)

normalized_rc_binary = (normalized_rc == 0).astype(int)
ptn.cnv.add_layer('normalized_read_counts_binary[zero/nonzero]', normalized_rc_binary.values)

normal_cnv_heatmap = ptn.cnv.heatmap(
    'normalized_read_counts_binary[zero/nonzero]',
    features = ['CDKN2A', 'KRAS', 'BRCA2', 'TP53', 'SMAD4'],
)
normal_cnv_heatmap.layout.coloraxis.colorbar.title = 'zero-reads'
normal_cnv_heatmap.layout.coloraxis.colorbar.len = 0.2
normal_cnv_heatmap.layout.coloraxis.colorbar.tickvals = [0, 1]


normal_cnv_heatmap.layout.coloraxis.colorscale = [
    [0, 'rgb(0, 0, 0)'], 
    [1, 'rgb(255, 225, 0)']
]
normal_cnv_heatmap.update_layout(font_family="Arial")
normal_cnv_heatmap.show()


# %% [focus on BRCA2 "chr13:32914437:GT/G"]
voi = ["chr13:32914437:GT/G"]
# 1. In BPA-2-RSX (supposingly low tumor purity)
# (1) how many cells have homdel?
af_df = sample_objs["BPA-2-RSX"].dna.get_attribute("AF_MISSING", features=voi, constraint="row")
af_df[af_df[voi] == -50].shape[0]
# (2) in cells w/o homdel, what's the distribution of AF?

# %% bulk data of BPA-3
# what's going on with chr17:7579801:G/C
bulk_germline_maf_f = "/lila/data/iacobuzc/haochen/bulk_analysis/tempo_pipeline/WGS_batch2/output/germline/RSX640_normal/combined_mutations/RSX640_tumor__RSX640_normal.germline.unfiltered.maf"
gmaf = pd.read_csv(bulk_germline_maf_f, sep = '\t', comment='#')
coi = ["t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count","n_alt_count"]
coi = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"] + coi
gmaf = gmaf.loc[gmaf["Hugo_Symbol"] == "TP53", coi]
# %%
