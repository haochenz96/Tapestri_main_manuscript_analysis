rule all:
	input:

rule analyze_ete_trees:
	input:
		tree_f = expand(
			"ete_trees_w_trunk/{patient_name}/{patient_name}_HZ_ETE_tree.refined.subclonal_snvs.trunk.pkl",
			patient_name = PATIENT_NAMES
		),
	output:
		tree_log_df = "1B.1_tree_analysis_stats.csv"

rule ete_tree_to_nhx:
	input:
		tree_f = expand(
			"ete_trees_w_trunk/{patient_name}/{patient_name}_HZ_ETE_tree.refined.subclonal_snvs.trunk.pkl",
			patient_name = PATIENT_NAMES
		),
	output:
		tree_nhx_f = expand(
			"ete_trees_w_trunk/{patient_name}/{patient_name}_HZ_ETE_tree.refined.subclonal_snvs.trunk.nhx",
			patient_name = PATIENT_NAMES
		),

rule nhx_to_ggtrees:
