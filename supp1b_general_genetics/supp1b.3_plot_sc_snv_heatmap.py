# %%
import sys
import argparse
import yaml
import pandas as pd
import mosaic.io as mio
import plotly.express as px
from tea.plots import plot_snv_clone
from pathlib import Path
import os 
os.chdir(Path(__file__).parent)
sample_name = "PC18-T0-pancreas_biopsy"

pt_h5 = "../data_compiled/fillout_h5/TP12.patient_wide.genotyped.h5"
pt = mio.load(pt_h5)
ss = pt[[i for i in pt.dna.barcodes() if sample_name in i]]

snv_f = "../data_compiled/manual_annotated_snv_lists/TP12-all_vars-voi.hz_curated.txt"
output_dir = Path(".")
# def main(args):
# %%

# # Add raw FALCON results to H5
# cn_assignment_f = list(opt_nclones_dir.glob(f"{sample_name}*assignment*csv"))[0]
# cn_assignment_df = pd.read_csv(cn_assignment_f, index_col = 0)
# print(f'[INFO] Loaded CN clone assignment file {cn_assignment_f}.')

# cn_clone_palette = dict(zip(
# 	np.sort(cn_assignment_df['clone_id'].unique()), 
# 	np.array(px.colors.qualitative.Set3)[np.sort(cn_assignment_df['clone_id'].unique())]
# 	))
# # rename the keys
# cn_clone_palette = {f"CN_clone-{k}": v for k, v in cn_clone_palette.items()}

# # add cn_clone info
# cn_assignment_df['cell_barcode_formatted'] = cn_assignment_df['cell_barcode'] + "-" + cn_assignment_df.index
# cn_assignment_df.set_index('cell_barcode_formatted', inplace=True)
# pt.dna.row_attrs['label'] = np.array(list(
# 	map(
# 		lambda x: f"CN_clone-{int(cn_assignment_df.loc[x, 'clone_id'])}", 
# 		pt.dna.barcodes())
# 	))
# pt.dna.set_palette(cn_clone_palette)
# pt.cnv.row_attrs['label'] = pt.dna.row_attrs['label']
# pt.cnv.set_palette(cn_clone_palette)
# # pt.cnv.get_gene_names()

num_cells = ss.dna.shape[0]

# read in high-quality SNVs
snv_df = pd.read_csv(snv_f, sep = '\t', index_col = 0, comment='#')
snv_df["annotation"].fillna("", inplace=True)
snv_ann_map = snv_df['HGVSp'].to_dict()

# ===== highlight vars with bulk annotation ====='
# highlight vars
germline_hom_var_col = '#00cc66' # green
germline_het_var_col = '#2693ff' # blue
somatic_var_col = '#ff0000' # red
likely_artifact_col = '#ffcc00' # yellow

for var_i in snv_ann_map:
	# germline
	if snv_df.loc[var_i, 'annotation'] == 'germline_HOM':
		snv_ann_map[var_i] = f'<span style="color:{germline_hom_var_col};">' + snv_ann_map[var_i] + '</span>'
	elif snv_df.loc[var_i, 'annotation'] == 'germline_HET':
		snv_ann_map[var_i] = f'<span style="color:{germline_het_var_col};">' + snv_ann_map[var_i] + '</span>'
	elif "somatic" in snv_df.loc[var_i, 'annotation']:
		snv_ann_map[var_i] = f'<span style="color:{somatic_var_col};">' + snv_ann_map[var_i] + '</span>'
	elif "artifact" in snv_df.loc[var_i, 'annotation']:
		snv_ann_map[var_i] = f'<span style="color:{likely_artifact_col};">' + snv_ann_map[var_i] + '</span>'
	else:
		pass

# plot heatmap
fig = plot_snv_clone(
	ss,
	sample_name=sample_name,
	story_topic = f'{sample_name}-high_conf_snvs',
	voi = snv_df.index.tolist(),
	attribute = "AF_MISSING",
	ann_map = snv_ann_map
)
# # map the x-axis ticks according to snv_ann_map
# fig.update_xaxes(ticktext = list(map(lambda x: snv_ann_map[x], snv_df.index.tolist())))

# save the figure
# write to PDF, with width proportion to the number of SNVs
fig.write_image(
	str(output_dir / f"{sample_name}-DNA-heatmap.pdf"), 
	width = 500 + 10 * len(snv_df.index.tolist())
)
print(f'[INFO] Saved heatmap for {sample_name} to {output_dir / f"{sample_name}-DNA-heatmap.pdf"}.')

# if __name__ == "__main__":
# 	main(sys.argv[1:])
# 	argparse.ArgumentParser(
# 		description="Plot SNV clone heatmap."
# 	)
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument(
# 		'--input_h5', 
# 		type=str, extension=['.h5'],
# 	)
# 	parser.add_argument(
# 		'--annotated_snv_list',
# 		type=str,
# 	)
# 	parser.add_argument(
# 		'--output_dir',
# 		type=str,
# 	)
# 	parser.add_argument(
# 		'--output_prefix',
# 		type=str,
# 	)
# 	parser.add_argument(
# 		"--cn_assignment",
# 		type=str, extension=['.csv'],
# 	)
# %%
