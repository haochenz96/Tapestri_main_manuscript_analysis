# %%
import pandas as pd
from pathlib import Path
import mosaic.io as mio
import numpy as np
from glob import glob
import plotly.express as px
import seaborn as sns
from tea.utils import rgb_string_to_hex
from tea.plots import mpl_config_params
COLOR_SEQUENCE= [rgb_string_to_hex(x) for x in px.colors.qualitative.Pastel]
mpl_config_params(font_size=12)

patients = ["M04", "BPA-1", "BPA-2", "BPA-3", "BPA-4", "BPA-5"]
# patients = ["BPA-5"]
brca_vars = [
    "chr13:32907415:T/A", # M04
    "chr13:32936733:A/T", # BPA-1
    "chr13:32914437:GT/G", # BPA-2
    "chr13:32914288:TAACC/T", # BPA-3
    "chr13:32911463:A/G", # BPA-4
    "chr13:32914437:GT/G", # BPA-5
]
output_dir = Path("/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/N_BRCA_analysis")

data_dir = Path('/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled')
condor_downstream_dir = Path('/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/0_condor_pipeline/condor_downstream/ete_trees_refined_subclonal_snvs')
# %%
for patient_name in patients:
    pt_h5 = glob(str(data_dir / "fillout_h5") + f"/*{patient_name}*.h5")
    if len(pt_h5) > 1:
        raise ValueError(f"Found multiple H5 files for {patient_name}")
    else:
        pt_h5 = pt_h5[0]
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

    # Add raw FALCON results to H5
    # HZ refined
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
    # pt.cnv.get_gene_names(amplicon_params, gene_name_col = 'gene')
    # pt.cnv.var = pd.DataFrame.from_dict(pt.cnv.col_attrs).set_index("id")
    # raw_rc = pt.cnv.get_attribute('read_counts', constraint='row')

    # normalized_rc = raw_rc / amp_params_df.loc[raw_rc.columns]['amplicon_factor'] / raw_rc.sum(axis=1).values[:, None]
    # # add to sample
    # pt.cnv.add_layer('normalized_read_counts', normalized_rc.values)

    # normalized_rc_binary = (normalized_rc == 0).astype(int)
    # pt.cnv.add_layer('normalized_read_counts_binary[zero/nonzero]', normalized_rc_binary.values)

    # read in high-quality SNVs
    snv_f = list((data_dir / "manual_annotated_snv_lists").glob(f"{patient_name}*.txt"))[0]
    snv_df = pd.read_csv(snv_f, sep = '\t', index_col = 0, comment='#')
    snv_df = snv_df[snv_df["annotation"].notna()]

    # ----- filter -----
    snv_df = snv_df[~snv_df['annotation'].str.contains("artifact")]
    # filter out germline HOM
    snv_df = snv_df[~snv_df['annotation'].str.contains("germline_HOM")]

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
        if var_i in(brca_vars):
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

    #
    # plot distribution of BRCA2 locus AF, split by tumor vs non-tumor barcode
    brca_var = snv_df.index[snv_df.index.isin(brca_vars)].tolist()
    if len(brca_var) > 1:
        raise ValueError(f"Found multiple gBRCA2 loci for {patient_name}")
    brca_df = pd.DataFrame(
        {"AF": pt.dna.get_attribute("AF", features=brca_var, constraint="row")[brca_var[0]].values,
        "clond_id": pt.dna.get_attribute("label", features=brca_var, constraint="row")[brca_var[0]].values},
        index=pt.dna.barcodes()
        )
    brca_df["cell_id"] = brca_df["clond_id"].apply(lambda x: "tumor" if x != "CN_clone-0" else "non-tumor")
    # map to tumor / non-tumor
    # statistical test
    from scipy.stats import ranksums

    # Perform Wilcoxon rank-sum test
    tumor_af = brca_df[brca_df["cell_id"]=="tumor"]["AF"]
    non_tumor_af = brca_df[brca_df["cell_id"]=="non-tumor"]["AF"]

    stat, p_value = ranksums(tumor_af, non_tumor_af)

    # Print the results
    print(f"Wilcoxon rank-sum test statistic: {stat}")
    print(f"p-value: {p_value}")

    # Save the results to a file
    with open(output_dir / f"{patient_name}_wilcoxon_test_results.txt", "w") as f:
        f.write(f"Wilcoxon rank-sum test statistic: {stat}\n")
        f.write(f"p-value: {p_value}\n")

    # plot 
    # sort brca_df by non-tumor first, tumor second
    brca_df = brca_df.sort_values(by="cell_id", ascending=False)
    g = sns.FacetGrid(brca_df, col="cell_id")
    g.map(sns.histplot, "AF", stat="percent", bins=50)
    # add p-value to the plot
    # g.axes[0][1].text(1, 0.5, f"p-value: {p_value:.2e}", fontsize=12, ha='right', va='center', transform=g.axes[0][1].transAxes)
    # save figure
    g.savefig(output_dir / f"{patient_name}_BRCA2_AF_distribution.png", dpi= 400)

# %%
