# %%
from pathlib import Path
import os, sys
import yaml

import mosaic.io as mio
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from tea.cravat import NONFUNC_SO
from tea.format import isNaN
from tea.plots import plot_snv_clone
from collections import Counter

from tea.cravat import get_technical_artifact_mask
from tea.format import CONDENSED_SNV_FORMAT, check_matrix_format
from tea.utils import get_simple_timestamp
timestamp = get_simple_timestamp()
from IPython import embed

	
# ===== IO =====
WD = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1_general_genetics"
os.chdir(WD)
WD = Path(".")
ANALYSIS_CONFIG = "pan_cohort_snv_heatmaps.yaml"
with open(ANALYSIS_CONFIG, 'r') as f:
	analysis_config = yaml.safe_load(f)

# ===== PARAMS =====
FILTER_SNV_ONLY = False
PON_OCCURENCE_FREQUENCY=0.5
# ==================

master_sample_sheet = analysis_config["master_sample_sheet"]
h5_dir = analysis_config["h5_dir"]
cravat_dir = analysis_config['cravat_dir']
falcon_dir = analysis_config['falcon_dir']

# ----- patient-sample map -----
master_sample_df = pd.read_excel(
	master_sample_sheet, 
	sheet_name = "all_sample_clinical_info"
	)
master_sample_df = master_sample_df[master_sample_df["censor"] == 0]
# ----- master CRAVAT -----
cravat_fs = list(Path(cravat_dir).glob("*.txt"))
assert len(cravat_fs) == master_sample_df['case'].nunique(), "Number of CRAVAT files != number of patients"
patient_cravats = {}
patient_sc_assignments = {}
# construct case to sample map
sample_map = {}
sample_h5s = {}
for patient_i in master_sample_df['case'].unique():
	patient_cravat = [f for f in cravat_fs if patient_i in f.stem]
	if not len(patient_cravat) == 1:
		raise ValueError(f"!= 1 CRAVAT found for {patient_i}")
	# patient_cravat = patient_cravat[0]
	# cravat_dfs.append(pd.read_csv(patient_cravat, sep='\t', index_col=0, header=[0.1]))
	patient_cravats[patient_i] = patient_cravat[0]
	
	patient_sc_assignment = list(Path(falcon_dir).glob(f"{patient_i}*assignment*.csv"))
	if not len(patient_sc_assignment) == 1:
		raise ValueError(f"!= 1 SC assignment found for {patient_i}")
	patient_sc_assignments[patient_i] = patient_sc_assignment[0]
 
	case_samples = master_sample_df[master_sample_df['case'] == patient_i]['sample']
	sample_map[patient_i] = case_samples.tolist()
	for sample_i in case_samples:
		h5 = list(Path(h5_dir).glob(f"*{sample_i}.mpileup.h5"))
		if not len(h5) == 1:
			raise ValueError(f"!= 1 H5 found for {sample_i}")
		sample_h5s[sample_i] = h5[0]

# # flatte patient_cravats to get master cravat
# cravat_dfs = []
# for patient_i in patient_cravats.values():
# 	for i in patient_i:
# 		cravat_dfs.append(i)

# cols_to_exclude = ["bulk_comparison", "Original Input", "Tapestri_result"]
# master_cravat = pd.concat([df.drop(columns=cols_to_exclude) for df in cravat_dfs], axis=1)
# master_cravat.to_csv("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1_general_genetics/master_cravat.csv", index=True)

# ----- blacklists -----
tumor_pon = pd.read_csv(analysis_config["snv_selection_params"]["tumor_pon"], index_col=0)
PON_OCCURENCE_FREQUENCY=0.5
threshold_sample_count = (len(tumor_pon.columns)-2) * PON_OCCURENCE_FREQUENCY
tumor_pon_blacklist = tumor_pon.index[tumor_pon["all_tumor_samples_occurence"] > threshold_sample_count]
# plot distribution of tumor_pon["all_tumor_samples_occurence"]
# px.histogram(tumor_pon["all_tumor_samples_occurence"]).show()

manual_blacklist = pd.read_csv(analysis_config["snv_selection_params"]["manual_blacklist"], index_col=0)
manual_snp_blacklist = manual_blacklist.index[manual_blacklist["type"] == "snp"]
manual_snp_loci_blacklist = [x.rsplit(":", 1)[0] for x in manual_snp_blacklist]
manual_snv_blacklist = manual_blacklist.index[manual_blacklist["type"] == "snv"]

# %% ====== technical filter SNVs, make falcon-clone clustered-sc-heatmaps ======
n=0
# per_patient_somatic_vars = {}
for patient_i, sample_names in sample_map.items():
	# if n>0:
	# 	break
	# n+=1
	if not patient_i == "BPA-3":
		continue
	sample_objs = {}
	patient_cravat = pd.read_csv(patient_cravats[patient_i], sep="\t", index_col=0, header=[0,1])

	for sample_i in sample_names:
		sample_objs[sample_i] = mio.load(sample_h5s[sample_i])
		# genotype
		if not "NGT" in sample_objs[sample_i].dna.layers:
			sample_objs[sample_i].dna.genotype_variants(
				min_dp = 8,
				min_alt_read = 3,
				assign_low_conf_genotype = True,
				)

	# ====== load cn_assignment_df ======
	cn_assignment_df = pd.read_csv(patient_sc_assignments[patient_i], index_col = 0)
	print(f'[INFO] Loaded CN clone assignment file {patient_sc_assignments[patient_i]}.')
	if not 'clone_id' in cn_assignment_df.columns:
		try: # try to add clone_id column
			cn_assignment_df['clone_id'] = cn_assignment_df['cluster_id']
		except:
			raise ValueError(f'[ERROR] `clone_id`/`cluster_id` column not found in CN clone assignment file!')
	unique_cluster_ids_sorted = np.sort(np.unique(cn_assignment_df['clone_id'].astype(int)))
	unique_cluster_ids_sorted_named = [ f"CN_clone-{clone_i}" for clone_i in unique_cluster_ids_sorted ]
	# embed()
	cn_clone_palette = dict(zip(unique_cluster_ids_sorted_named, np.array(px.colors.qualitative.Set3)[unique_cluster_ids_sorted]))

	for sample_i in sample_names:
		# add cn_clone info
		if not sample_i in cn_assignment_df.index:
			raise ValueError(f"{sample_i} not in cn_assignment_df.index")
		cn_assignment_dict = cn_assignment_df.loc[sample_i,:].set_index('cell_barcode').to_dict()['clone_id']
		try:
			sample_objs[sample_i].dna.row_attrs['label'] = np.array(list(map(lambda x: f"CN_clone-{int(cn_assignment_dict[x])}", sample_objs[sample_i].dna.barcodes())))
			sample_objs[sample_i].dna.set_palette(cn_clone_palette)
		except KeyError:
			print (f'[ERROR] {sample_i} not found in cn_assignment_df.index!')
			# embed()
			cn_err = 1
		else:
			sample_objs[sample_i].cnv.row_attrs['label'] = sample_objs[sample_i].dna.row_attrs['label']
			sample_objs[sample_i].cnv.set_palette(cn_clone_palette)
			# sample_objs[sample_i].cnv.get_gene_names()

		num_cells = sample_objs[sample_i].dna.shape[0]
		print(f'[INFO] {sample_i} has {num_cells} cells.')
		# if args.write_cn_clone_added_h5s and cn_err == 0:
		# 	out_dir = WD / 'cn_clone_added_h5s'
		# 	out_dir.mkdir(exist_ok=True, parents=True)
		# 	try:
		# 		mio.save(sample_objs[sample_i], str(out_dir / f"{sample_i}_cn_clone_added.h5"))
		# 		print(f'[INFO] {out_dir / f"{sample_i}_cn_clone_added.h5"} saved.')
		# 	except FileExistsError:
		# 		print(f'[WARNING] {out_dir / f"{sample_i}_cn_clone_added.h5"} already exists! Overwriting...')
		# 		try: 
		# 			# delete and try again
		# 			os.remove(str(out_dir / f"{sample_i}_cn_clone_added.h5"))
		# 			mio.save(sample_objs[sample_i], str(out_dir / f"{sample_i}_cn_clone_added.h5"))
		# 		except:
		# 			print(f'[ERROR] {out_dir / f"{sample_i}_cn_clone_added.h5"} already exists and cannot be overwritten!')

	# ====== set up for plotting ======
	####################################
	# ### ----- SNV selection ----- ####
	####################################

	snv_selection_params = analysis_config['snv_selection_params']

	# 1. mutational prevalence
	MUT_PREV = snv_selection_params['mut_prev_threshold']

	# 2. technical artifact filters
	bq_prev_threshold = snv_selection_params['bq_prev_threshold']
	if bq_prev_threshold is not None and type(bq_prev_threshold) is not float: # single-value
		raise ValueError(f"bq_prev_threshold must be float, not {type(bq_prev_threshold)}")

	normals_occurences = snv_selection_params['normals_occurences'] if 'normals_occurences' in snv_selection_params else 3 # <--- default to 3

	if 'ado_threshold' in snv_selection_params and snv_selection_params['ado_threshold'] is not None:
		print(f"[INFO] filtering out SNVs with ADO > {snv_selection_params['ado_threshold']} in ANY sample...")
		ado_threshold = snv_selection_params['ado_threshold']
	else:
		ado_threshold = None

	# 3. functional SNVs
	topic = snv_selection_params['topic']
	if not type(topic) is str: # single-value
		raise ValueError(f"topic must be str, not {type(topic)}")
	try: 
		func_only = snv_selection_params['func_only']
		func_only = bool(func_only)
	except KeyError:
		func_only = False
	except TypeError:
		func_only = False

	# 4. germline SNPs
	germline_attrs = {}
	for germline_attr_i in ['select_hom_germline_snps_af', 'rescue_1000genome_af']:
		if not germline_attr_i in snv_selection_params:
			germline_attrs[germline_attr_i] = False
		else:
			germline_attrs[germline_attr_i] = snv_selection_params[germline_attr_i]

	# 5. plotting params
	attribute = snv_selection_params['attribute']
	if not type(attribute) is list:
		attribute = [attribute]

	# whitelist snvs
	if 'whitelist_snvs' in snv_selection_params and snv_selection_params['whitelist_snvs']:
		whitelist_snvs = set(snv_selection_params['whitelist_snvs'])
	else:
		whitelist_snvs = []

	# blacklist snvs
	if 'blacklist_snvs' in snv_selection_params and snv_selection_params['blacklist_snvs']:
		if type(snv_selection_params['blacklist_snvs']) is list:
			blacklist_snvs = set(snv_selection_params['blacklist_snvs'])
		else:
			try:
				blacklist_snv_df = pd.read_csv(snv_selection_params['blacklist_snvs'], index_col=0)
				if not check_matrix_format(blacklist_snv_df, CONDENSED_SNV_FORMAT):
					raise ValueError(f"[ERROR] blacklist_snvs file not in the correct format (index column needs to be in condensed SNV format).")
				blacklist_snvs = set(blacklist_snv_df.index)
			except:
				raise ValueError(f"[ERROR] blacklist_snvs should either be a list of SNVs or a path to a CSV file whose index is the list of SNVs.")
	else:
		blacklist_snvs = []
		
	# # highlight snvs
	# if 'highlight_snvs' in snv_selection_params and snv_selection_params['highlight_snvs']:
	# 	highlight_snvs = set(snv_selection_params['highlight_snvs'])
	# else:
	# 	highlight_snvs = []

	print(f"""
			===== 
			[INFO] MUT_PREV = {MUT_PREV} 
			=====
			""")
	
	(WD / f"{topic}_mut_prev={MUT_PREV}" / "sc_heatmaps").mkdir(exist_ok=True, parents=True) 
	(WD / f"{topic}_mut_prev={MUT_PREV}" / "vois").mkdir(exist_ok=True, parents=True)

	# ----- @HZ 07/17/2023 filter T>C artifacts that are rare -----
	if 'filter_TtoC_artifact' in snv_selection_params:
		try:
			filter_TtoC_artifact = snv_selection_params['filter_TtoC_artifact']['filter'] 
			filter_TtoC_artifact_lower_thres = snv_selection_params['filter_TtoC_artifact']['lower_thres']
			filter_TtoC_artifact_upper_thres = snv_selection_params['filter_TtoC_artifact']['upper_thres']
		except KeyError:
			raise ValueError(f"[ERROR] filter_TtoC_artifact must have keys ['filter', 'lower_thres', 'upper_thres']")
		else:
			if filter_TtoC_artifact_lower_thres >= filter_TtoC_artifact_upper_thres:
				raise ValueError(f"[ERROR] filter_TtoC_artifact_lower_thres must be strictly smaller than filter_TtoC_artifact_upper_thres")
	else:
		if MUT_PREV < 0.01:
			print(f"[INFO] MUT_PREV is lower than default upper_thres (0.01) for T>C filter. The filter will be applied.")
			filter_TtoC_artifact = True
			filter_TtoC_artifact_lower_thres = MUT_PREV
			filter_TtoC_artifact_upper_thres = 0.01 
		else:
			print(f"[WARNING] MUT_PREV is higher than default upper_thres (0.01) for T>C filter. The filter will not be applied.")
			filter_TtoC_artifact = False
	# -----------------------------------------------------------------

	voi_union = set()
	voi_count_union = Counter()
	# ann_map_union = {}
	bulk_germline_vars = set()
	bulk_somatic_vars = set()
	TtoC_artifact_blacklist = set()
	ado_blacklist = set()

	# ====== for each sample_i, filter SNVs and get a union set for plotting  ======
	for sample_i in sample_names:
		num_cells = sample_objs[sample_i].dna.shape[0]

		snvs = patient_cravat.index
		mut_prev_series = sample_objs[sample_i].dna.get_attribute("mut_filtered", constraint="row").sum(axis=0)[snvs]
	
		technical_mask = get_technical_artifact_mask(
			snvs,
			patient_cravat, 
			num_cells = num_cells, 
			mut_prev_series = mut_prev_series,
			bq_prev_threshold = bq_prev_threshold, 
			normals_pon_occurence = normals_occurences, 
			rescue_1000genome_af = germline_attrs['rescue_1000genome_af'], 
			filter_broad_wes_pon = False,
			ado_threshold = None,
			ngt_df = None,
			filter_TtoC_artifact = None,
			blacklist = None,
			)			

		# filters on functional SNVs
		if func_only:
			mask = technical_mask & ~patient_cravat[('Variant Annotation', 'Sequence Ontology')].isin(NONFUNC_SO)
		else:
			mask = technical_mask

		voi = snvs[mask].tolist()

		# BLACKLIST T>C artifacts that are rare in EACH SAMPLE
		if filter_TtoC_artifact:
			tc_mask = (snvs.str.endswith('T/C')) & (mut_prev_series >= filter_TtoC_artifact_lower_thres * num_cells) & (mut_prev_series <= filter_TtoC_artifact_upper_thres * num_cells)
			print(f"[INFO] {sample_i} has {tc_mask.sum()} T>C artifacts that are rare. They will be added to the blacklist.")
			TtoC_artifact_blacklist = TtoC_artifact_blacklist.union(
				set(snvs[tc_mask].tolist())
			)

		# BLACKLIST SNVs that have too high ADO in any sample
		if ado_threshold is not None:
			ado_high_vars = sample_objs[sample_i].dna.ids()[
				(sample_objs[sample_i].dna.get_attribute('NGT',constraint='row') == 3).sum(axis=0) > (ado_threshold*num_cells)
			]
			print(f"[INFO] {sample_i} has {len(ado_high_vars)} SNVs with ADO > {ado_threshold}. They will be added to the blacklist.")
			ado_blacklist = ado_blacklist.union(set(ado_high_vars))

		# filters on mut_prev_threshold
		prev_filtered_vars = sample_objs[sample_i].dna.ids()[
			sample_objs[sample_i].dna.get_attribute("mut_filtered", constraint="row").sum(axis=0) >= (MUT_PREV * num_cells)
		]
		# take intersection
		voi = [ v for v in voi if v in prev_filtered_vars ]
		print(f'{sample_i}: {len(voi)} / {len(patient_cravat)} SNVs are kept (min prev: {MUT_PREV})')

		# voi = list(set(voi).union(whitelist_snvs))
		voi_mut_prev = Counter(mut_prev_series.to_dict())

		# ====== save bulk annotated vars ======

		# select SNVs detected in matched bulk normal
		# embed()
		# sys.exit()
		try:
			bulk_normal_vars = patient_cravat.index[(patient_cravat[('bulk_comparison', 'bulk-matched_bulk_normal-AF')] > 0)]
		except KeyError:
			print(f'bulk normal annotation not found in CRAVAT DF for {sample_i}')
		else:
			if len(bulk_normal_vars) == 0:
				print(f'[WARNING] No bulk normal SNVs detected in {sample_i}')
			bulk_germline_vars.update(bulk_normal_vars)

		# select SNVs detected in matched bulk cohort
		try:
			bulk_cohort_vars = patient_cravat.index[(patient_cravat[('bulk_comparison', 'bulk-matched_bulk_cohort-AF')] > 0)]
		except KeyError:
			print(f'bulk tumor annotation not found in CRAVAT DF for {sample_i}')
		else:
			if len(bulk_cohort_vars) == 0:
				print(f'[WARNING] No bulk cohort SNVs detected in {sample_i}')
			bulk_somatic_vars.update(bulk_cohort_vars)

		# ===== get a union of all samples SNVs =====
		voi_union = voi_union.union(set(voi))
		voi_count_union += voi_mut_prev
		# ann_map_union.update(ann_map)

	# select homozygous germline SNVs, if specified
	if germline_attrs['select_hom_germline_snps_af']:
		alt_combined = pd.concat(
			[ sample_objs[sample_i].dna.get_attribute('alt_read_count', constraint='row') for sample_i in sample_names ], axis=0
		)
		dp_combined = pd.concat(
			[ sample_objs[sample_i].dna.get_attribute('DP', constraint='row') for sample_i in sample_names ], axis=0 
		)
		overall_af = alt_combined.sum(axis=0) / dp_combined.sum(axis=0)
		germline_hom_snps_from_tapestri = [x for x in voi_union if overall_af[x] > germline_attrs['select_hom_germline_snps_af']]
		print(f"[INFO] selected {len(germline_hom_snps_from_tapestri)} homozygous germline SNVs (AF > {germline_attrs['select_hom_germline_snps_af']})")
	else:
		germline_hom_snps_from_tapestri = []
	# embed()
	# remove SNVs that are blacklisted
	print(f"[DEBUG] {len(voi_union)} SNVs before blacklist filtering")
	voi_union = voi_union.difference(TtoC_artifact_blacklist)
	print(f"[DEBUG] {len(voi_union)} SNVs after TtoC blacklist filtering")
	voi_union = voi_union.difference(ado_blacklist)
	print(f"[DEBUG] {len(voi_union)} SNVs after ADO blacklist filtering")
	voi_union = voi_union.difference(blacklist_snvs)
	print(f"[DEBUG] {len(voi_union)} SNVs after manual blacklist filtering")

	# add whitelist
	for v in whitelist_snvs:
		if not v in patient_cravat.index:
			raise ValueError(f"[ERROR] {v} from whitelist is not found in CRAVAT DF index!")
	voi_union = list(voi_union.union(whitelist_snvs))

	# remove SNV/SNPs that are manually blacklisted
	manual_artifact_vars = [
	 	x for x in bulk_germline_vars if x.rsplit(":", 1)[0] in manual_snp_loci_blacklist
		] + [
		x for x in voi_union if x in manual_snv_blacklist
	 	]

	bulk_germline_vars = [x for x in bulk_germline_vars if x not in manual_artifact_vars]

	# ===== annotation =====
	ann = patient_cravat.loc[voi_union, :].index.map(
		lambda x: 
		patient_cravat.loc[x, ('Variant Annotation', 'Gene')] + ' ' + patient_cravat.loc[x, ('Variant Annotation', 'Protein Change')] if not isNaN(patient_cravat.loc[x, ('Variant Annotation','Protein Change')])
		else patient_cravat.loc[x, ('Variant Annotation','Gene')] + ' ' + patient_cravat.loc[x, ('Variant Annotation','Sequence Ontology')]
	)
	ann_map_union = dict(zip(voi_union, ann))
	
	# filter somatic vars
	somatic_vars = set(voi_union) - set(germline_hom_snps_from_tapestri) - set(bulk_germline_vars) - set(manual_artifact_vars)
	# filter out PON blacklist
	somatic_vars_in_pon = [x for x in somatic_vars if x in tumor_pon_blacklist]
	all_artifacts = set(manual_artifact_vars).union(set(somatic_vars_in_pon))

	voi_sorted = sorted(voi_union, key=voi_count_union.get, reverse=True)
	# label germline and bulk somatic vars:
	def __annotate_snvs(var, germline_hom_vars, other_bulk_germline_vars, bulk_somatic_vars, artifact_vars):
		if var in artifact_vars:
			return "likely_artifact"		
		elif var in germline_hom_vars:
			return "germline_HOM"
		elif var in other_bulk_germline_vars:
			return "germline_HET"
		elif var in bulk_somatic_vars:
			return "bulk_somatic"
		else:
			return "NA"
 
	with open(WD / f"{topic}_mut_prev={MUT_PREV}" / "vois" / f'{patient_i}_{topic}_mut_prev={MUT_PREV}_voi.txt', 'w') as f:
		voi_df = pd.DataFrame.from_dict(voi_count_union, orient='index', columns=['mut_prev']).loc[voi_sorted]
		voi_df['HGVSp'] = voi_df.index.map(lambda x: ann_map_union[x])
		voi_df['annotation'] = voi_df.index.map(
			lambda x: __annotate_snvs(x, germline_hom_snps_from_tapestri, bulk_germline_vars, bulk_somatic_vars, all_artifacts)
		)
		voi_df.index.name = 'condensed_format'
		f.write(f"# {timestamp}\n")
		f.write(f"# patient_name: {patient_i}")
		f.write(f"# sample_names: {sample_names}")
		f.write('# ======  [snv_selection_params] ======\n')
		f.write(f"# attribute: {attribute}\n")
		f.write(f"# mut_prev_threshold: {MUT_PREV}\n")
		f.write(f"# bq_prev_threshold: {bq_prev_threshold}\n")
		f.write(f"# topic: {topic}\n")
		f.write(f"# func_only: {func_only}\n")
		f.write(f"# normals_occurences: {normals_occurences}\n")
		f.write(f"# ado_threshold: {ado_threshold}\n")
		f.write(f"# germline_attrs: {germline_attrs}\n")
		f.write(f"# filter_TtoC_artifact: {snv_selection_params['filter_TtoC_artifact']}\n")
		f.write(f"# whitelist_snvs: {whitelist_snvs}\n")
		if 'blacklist_snvs' in snv_selection_params:
			f.write(f"# blacklist_snvs: {snv_selection_params['blacklist_snvs']}\n")
		else:
			f.write(f"# blacklist_snvs: None\n")
		f.write('# =====================================\n')

		voi_df.to_csv(f, sep='\t', index=True, header=True)
  
		per_patient_somatic_vars[patient_i] = voi_df.index[
	  		(~voi_df["annotation"].str.contains("germline")) &
			(~voi_df["annotation"].str.contains("artifact"))
		].tolist()
		
	if FILTER_SNV_ONLY:
		# stop after filtering
		print(f"[INFO] SNV filtering done. Skipping heatmap making...")
		continue

	# ====== plot ======
	# ===== highlight vars with bulk annotation ====='
	# highlight vars
	germline_hom_var_col = '#00cc66' # green
	germline_het_var_col = '#2693ff' # blue
	somatic_var_col = '#ff0000' # red
	artifact_var_col = "#ffcc00"

	for var_i in ann_map_union:
		if var_i in manual_artifact_vars:
			ann_map_union[var_i] = f'<span style="color:{artifact_var_col};">' + ann_map_union[var_i] + '</span>'
		elif var_i in somatic_vars_in_pon:
			ann_map_union[var_i] = f'<span style="color:{artifact_var_col};">' + ann_map_union[var_i] + '</span>'
		# germline
		elif var_i in germline_hom_snps_from_tapestri:
			ann_map_union[var_i] = f'<span style="color:{germline_hom_var_col};">' + ann_map_union[var_i] + '</span>'
		elif var_i in bulk_germline_vars:
			ann_map_union[var_i] = f'<span style="color:{germline_het_var_col};">' + ann_map_union[var_i] + '</span>'
		elif var_i in bulk_somatic_vars:
			ann_map_union[var_i] = f'<span style="color:{somatic_var_col};">' + ann_map_union[var_i] + '</span>'
		else:
			pass

	for sample_i in sample_names:
		for attribute_i in attribute:
			fig = plot_snv_clone(
				sample_objs[sample_i],
				sample_name=sample_i,
				story_topic = f'high_conf_mutations-{topic}',
				voi = voi_sorted,
				attribute = attribute_i,
				ann_map = ann_map_union
			)
			fig.update_layout(
				width = max(len(voi_sorted) * 20, 600),
			)
			fig.write_image(
				WD / f"{topic}_mut_prev={MUT_PREV}" / "sc_heatmaps" / f'{sample_i}-{topic}-{attribute_i}.pdf',
			)

			# fig.show()

# %% ===== construct per-patient scMAF =====
sample_mafs = {}
DATA_DIR = Path("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled")
PATIENT_SAMPLE_MAP_F = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/1_general_genetics/Tapestri_all_patient_sample_map.yaml"
with open(PATIENT_SAMPLE_MAP_F, "r") as f:
	patient_sample_map = yaml.safe_load(f)
MASTER_CRAVAT_F="/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1_general_genetics/master_cravat.csv"
master_cravat_df = pd.read_csv(MASTER_CRAVAT_F, index_col=0, header=[0,1])

# manual tuning:
per_patient_somatic_vars["M04"] = set(per_patient_somatic_vars["M04"]) | {"chr13:32907415:T/A"}
per_patient_somatic_vars["M07"] = set(per_patient_somatic_vars["M07"]) | {"chr17:7577556:C/T"}
per_patient_somatic_vars["RA15_06"] = set(per_patient_somatic_vars["RA15_06"]) |{"chr18:48593507:C/T"}
per_patient_somatic_vars["RA16_08"] = set(per_patient_somatic_vars["RA16_08"]) |{"chr2:25457242:C/T"}
per_patient_somatic_vars["RA17_13"] = set(per_patient_somatic_vars["RA17_13"]) |{"chr18:48581173:G/T", "chr9:21974695:G/GT"}
per_patient_somatic_vars["RA19_02"] = set(per_patient_somatic_vars["RA19_02"]) |{"chr12:25398284:C/A"}
per_patient_somatic_vars["RA19_10"] = set(per_patient_somatic_vars["RA19_10"]) |{"chr12:25398284:C/A"}
per_patient_somatic_vars["RA19_21"] = set(per_patient_somatic_vars["RA19_21"]) |{"chr1:209969821:C/T", "chr3:178936091:G/A"}

per_patient_somatic_vars["RA20_05"] = set(per_patient_somatic_vars["RA20_05"]) |{"chr12:25398284:C/A", "chr17:7578496:A/T"}


# # %% RA19_10 debug
# voi_df = pd.read_csv(
# 	"/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1_general_genetics/all_vars_mut_prev=0.005/vois/RA19_10_all_vars_mut_prev=0.005_voi.txt",
# 	sep="\t", index_col=0, comment="#"
# )
# per_patient_somatic_vars["RA19_10"] = voi_df.index[
# 	(~voi_df["annotation"].str.contains("germline", na=False)) &
# 	(~voi_df["annotation"].str.contains("artifact", na=False))
# ].tolist()
# %%
for patient_i in master_sample_df['case'].unique():
	if patient_i in sample_mafs:
		continue
	print(f"[INFO] processing {patient_i}")
	voi = list(per_patient_somatic_vars[patient_i])
	mpileup_h5 = list((DATA_DIR / "fillout_h5").glob(f"{patient_i}*.h5"))
	if not len(mpileup_h5) == 1:
		raise ValueError(f"!= 1 mpileup h5 found for {patient_i}")
	sample_obj = mio.load(str(mpileup_h5[0]))
	sample_obj.dna = sample_obj.dna[sample_obj.dna.barcodes(), voi]
	sample_obj.dna.genotype_variants(
		min_dp=8,
		min_alt_read=3,
		assign_low_conf_genotype=True,
	)
	num_cells = sample_obj.dna.shape[0]


	# fetch raw layers
	mut_filtered_layer = sample_obj.dna.get_attribute(
		"mut_filtered", features=voi, constraint="row"
	)
	mut_unfiltered_layer = sample_obj.dna.get_attribute(
		"mut_unfiltered", features=voi, constraint="row"
	)
	# get AF in cells that are mutated for each variant
	af_mut_filtered_layer = sample_obj.dna.get_attribute(
		"AF", features=voi, constraint="row"
	)[mut_filtered_layer.astype(bool)]
	af_mut_unfiltered_layer = sample_obj.dna.get_attribute(
		"AF", features=voi, constraint="row"
	)[mut_unfiltered_layer.astype(bool)]

	# create sample MAF
	sample_maf = pd.DataFrame(
		columns=[
			"condensed_format",
			"patient_name",
			"sample_name",
			"gene_name",
			"snv_class",
			"HGVSp",
			"total_num_cells",
			"mut_filtered_num_cells",
			"mut_filtered_sc_freq",
			"mut_unfiltered_num_cells",
			"mut_unfiltered_sc_freq",
			"mean_sc_AF_in_mut_filtered",
			"median_sc_AF_in_mut_filtered",
			"driver_snv",
			"num_cancer_sc",
			"num_cancer_sc_w_var",
			"CCF",
			"bulk_tumor_fillout_VAF",
			"bulk_normal_VAF",
		]
	)
	sample_maf["condensed_format"] = voi
	sample_maf["patient_name"] = patient_i
	sample_maf["sample_name"] = patient_i

	sample_maf[["gene_name", "snv_class", "HGVSp"]] = master_cravat_df.loc[
		voi,
		[
			("Variant Annotation", "Gene"),
			("Variant Annotation", "Sequence Ontology"),
			("Variant Annotation", "Protein Change"),
		],
	].values

	sample_maf["total_num_cells"] = sample_obj.dna.shape[0]
	sample_maf[
		[
			"mut_filtered_num_cells",
			"mut_unfiltered_num_cells",
			"mean_sc_AF_in_mut_filtered",
			"median_sc_AF_in_mut_filtered",
		]
	] = pd.concat(
		[
			mut_filtered_layer.sum(axis=0),
			mut_unfiltered_layer.sum(axis=0),
			af_mut_filtered_layer.mean(axis=0, skipna=True),
			af_mut_filtered_layer.median(axis=0, skipna=True),
		],
		axis=1,
	).values

	sample_maf[["mut_filtered_sc_freq", "mut_unfiltered_sc_freq"]] = sample_maf[
		["mut_filtered_num_cells", "mut_unfiltered_num_cells"]
	].divide(sample_maf["total_num_cells"], axis=0)

	# # %%
	# case BPA-3 needs to be considered specifically, as we need to define cancer cells by CN calling result (CN clone 3 from the final ETE tree)
	if patient_i == "BPA-3":
		# need to use FALCON clone to define cancer cells
		cn_assignment_f = list(Path("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/1_general_genetics/1A_sc_oncoprint").glob(f"{patient_i}*assignment*csv"))
		if len(cn_assignment_f) != 1:
			raise ValueError(f"[ERROR] for patient {patient_i}, CN assignment file not found")
		cn_assignment_df = pd.read_csv(cn_assignment_f[0]).sort_values(["sample_name", "cell_barcode"])
		# cn_assignment_df["sample_name"] = cn_assignment_df["sample_name"].map({"BPA-3-RSX": "BPA-3-T", "BPA-3-N": "RSX640_N"})
		cn_assignment_df["barcode-sample"] = cn_assignment_df["cell_barcode"] + "-" + cn_assignment_df["sample_name"]
		tumor_cells = cn_assignment_df.loc[cn_assignment_df["final_clone_id"] == 3, "barcode-sample"].values
		tumor_cells_bool = mut_filtered_layer.index.isin(tumor_cells)
		sample_maf["driver_snv"] = "refined_tree_clone3"
		# embed()
	elif "driver_snv" in patient_sample_map[patient_i]:
		driver_snv = patient_sample_map[patient_i]["driver_snv"]
		sample_maf["driver_snv"] = driver_snv
		if driver_snv in mut_filtered_layer.columns:
			tumor_cells_bool = mut_filtered_layer[driver_snv].astype(bool)
		else:
			raise KeyError(
					f"[ERROR] for patient {patient_i}, driver SNV not found in mut_filtered_layer"
				)
	else:
	   raise KeyError("[ERROR] for patient {patient_i}, driver SNV not available")
	mut_filtered_layer_cancer_cells = mut_filtered_layer[tumor_cells_bool]
	num_cancer_sc = mut_filtered_layer_cancer_cells.shape[0]
	num_cancer_sc_w_var = mut_filtered_layer_cancer_cells.sum(axis=0)
	sample_maf["num_cancer_sc"] = num_cancer_sc
	sample_maf["num_cancer_sc_w_var"] = num_cancer_sc_w_var.values
	ccf_vals = (num_cancer_sc_w_var / num_cancer_sc).values
	sample_maf["CCF"] = ccf_vals
	# bulk coverage
	try:
		sample_maf["bulk_tumor_fillout_VAF"] = master_cravat_df.loc[
			voi, [("bulk_comparison", "bulk-matched_bulk_cohort-AF")]
		].values
	except:
		print(f"[ERROR] for sample {patient_i}, bulk tumor fillout VAF not available")
		e = 1
	try:
		sample_maf["bulk_normal_VAF"] = master_cravat_df.loc[
			voi, [("bulk_comparison", "bulk-matched_bulk_normal-AF")]
		].values
	except:
		print(f"[ERROR] for sample {patient_i}, bulk normal VAF not available")
		e = 1
	sample_mafs[patient_i] = sample_maf

# %%
pan_cohort_scMAF = pd.concat(sample_mafs, ignore_index=True)
pan_cohort_scMAF.to_csv(WD / f"pan_cohort_somatic_scMAF_mut_prev={MUT_PREV}.csv", index=False, header=True)

# %%
GOI=["KRAS", "TP53","CDKN2A", "SMAD4", "TGFBR1", "TGFBR2", "BRCA2", "ATM", "FGFR1",]

pan_cohort_scMAF["functional?"] = ~pan_cohort_scMAF["snv_class"].isin(NONFUNC_SO)
pan_cohort_scMAF["driver?"] = pan_cohort_scMAF["gene_name"].isin(GOI)
pan_cohort_scMAF["functional_driver?"] = (~pan_cohort_scMAF["snv_class"].isin(NONFUNC_SO)) & (pan_cohort_scMAF["gene_name"].isin(GOI))
import seaborn as sns
import matplotlib.pyplot as plt
# use Seaborn stacked barplot - make a histogram of CCF, colored by functional status
# log scale

# set default font to be Arial, font size 8
plt.rcParams['font.family'] = "serif"
plt.rcParams['font.size'] = 8

pan_cohort_scMAF_filtered = pan_cohort_scMAF[pan_cohort_scMAF["CCF"]>0]

fig, ax= plt.subplots(1, 2, figsize=(12, 5))
sns.histplot(
	pan_cohort_scMAF_filtered[pan_cohort_scMAF_filtered["functional_driver?"]],
	x="CCF",
	multiple="stack",
	bins=100,
	color="red",
	ax=ax[0],
)
ax[0].set_title("Functional driver SNVs")
ax[0].set_yscale("log", base=2)
# set y axis range
ax[0].set_ylim(1, 2**10)
sns.histplot(
	pan_cohort_scMAF_filtered[pan_cohort_scMAF_filtered["driver?"]],
	x="CCF",
	multiple="stack",
	bins=100,
	color="blue",
	ax=ax[1],
)
ax[1].set_title("All SNVs")
ax[1].set_yscale("log", base=2)
ax[1].set_ylim(1, 2**10)

# %%
fig.savefig(WD / f"pan_cohort_somatic_scMAF_mut_prev={MUT_PREV}_CCF_hist.pdf")

# %%
