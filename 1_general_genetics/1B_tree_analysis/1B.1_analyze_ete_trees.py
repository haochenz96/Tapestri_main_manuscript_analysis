# %% 
import ete3 
import pandas as pd
import pickle
from pathlib import Path
from collections import Counter

pickle_dir = Path("../../condor_pipeline/condor_downstream/ete_trees_refined_subclonal_snvs")
goi = ["KRAS", "TP53","CDKN2A", "SMAD4", "SMAD2", "SMAD3", "TGFBR1", "TGFBR2", "ACVR1B", "BMPR1A", "ARID1A", "ARID2", "BRCA2", "ATM", "BAP1", "FGFR1","FBXW7", "RNF43", "POLD1", "IRF6", "GATA6", "MYC", "MTOR"]
sample_sheet = pd.read_excel("../../Tapestri_batch2_samples_MASTER.xlsx", sheet_name="all_case_genetics", index_col=0, skiprows=1)
sample_sheet = sample_sheet[sample_sheet["censor"] == 0]
poi = sample_sheet.index
# rename:
# BPA-2 -> BPA-2-IR
# BPA-5 -> BPA-5-RSX
poi = [x if x not in ["BPA-2", "BPA-5"] else f"{x}-IR" if x == "BPA-2" else f"{x}-RSX" for x in poi]

# for f in pickle_dir.glob("*solT_cell"):
#     print(f"Processing {f.stem}")

# %%
master_maf_dict = {}
tree_log_df = pd.DataFrame(columns=["patient_name","truncal_snvs","total_snvs","truncal_cnvs","total_cnvs"])

for patient_name in poi:
	# initiate a df in MAF format to log LOH events of genes of interest
	maf = pd.DataFrame(columns=["patient_name","gene_name","alteration_class","CCF"])
	# initiate a df that logs tree stats:
	# - number of truncal SNVs (immediately following root) against total SNVs
	# - number of truncal SNPs (immediately following root) against total SNPs

	# patient_name = "RA19_10"
	f = pickle_dir / f"{patient_name}/{patient_name}_HZ_ETE_tree.refined.subclonal_snvs.pkl"
	if patient_name in ["BPA-2-IR", "BPA-5-RSX"]:
		patient_name = patient_name.rsplit("-",1)[0]
		
	# first, we need to find the node that represents the trunk. The complexities are:
	# 1. in some cases, the branch following the root is not the trunk. Need to move down the tree until a somatic SNV event is found
	# 2. in some cases (e.g. RA19_21), the root has two children; then a selection needs to be made
	with open(f, "rb") as f:
		tree = pickle.load(f)
		# 1. find trunk
		trunk = tree.get_children()[0]
		if len(tree.get_children()) > 1:
			# if the root has two children, we need to select the one that has the most events
			event_count = Counter()
			events = {}
			for child in tree.get_children():
				event_count[child] = len(child.somatic_snv_events) + len(child.germline_snp_events)
				events[child] = (child.somatic_snv_events | child.germline_snp_events).values()
			trunk = event_count.most_common(1)[0][0]
			print(f"for patient {patient_name}, selected trunk with events: {events[trunk]}")
		if len(trunk.somatic_snv_events) == 0:
			# if the trunk has no somatic SNV events, we need to move down the tree until we find one; unless this case is BPA-3 which has no driver SNVs
			if patient_name == "BPA-3":
				print(f"for patient {patient_name}, no driver SNV present")
			else:
				for node in trunk.traverse("preorder"):
					if len(node.somatic_snv_events) > 0:
						print(f"for patient {patient_name}, moved down to {node.name} with somatic snvs: {node.somatic_snv_events}")
						trunk = node
						break
		# total tumor size is sum of all clones' clone_sizes after the trunk
		total_tumor_size = sum([d.clone_size for d in trunk.traverse("preorder") if "clone_size" in d.__dict__])
		# need to isolate somatic SNV GAIN events
		truncal_snv_events = trunk.somatic_snv_events.values()
		n_truncal_snv_GAINS = len([i for i in truncal_snv_events if i[1] == "GAIN"])
		n_truncal_cnv_events = len(trunk.germline_snp_events) + len(truncal_snv_events) - n_truncal_snv_GAINS
		tree_log_df.loc[patient_name, ["truncal_snvs", "truncal_cnvs"]] = [n_truncal_snv_GAINS, n_truncal_cnv_events]
		c = Counter()
		# 2. count events
		for node in trunk.traverse("preorder"):
			# print(node.germline_snp_events)
			# print(node.somatic_snv_events)
			descendants_size = sum([d.clone_size for d in node.traverse("preorder") if "clone_size" in d.__dict__])
			ccf = descendants_size / total_tumor_size
			# A. log the SNP events
			for _, snp_i in node.germline_snp_events.items():
				gene_name = snp_i[0].split(" ")[0]
				if gene_name in goi and gene_name not in maf.gene_name.values:
					print(f"for patient {patient_name}, found LOH event in {gene_name} with CCF {ccf}")
					maf = pd.concat(
						[maf, pd.DataFrame({"patient_name":patient_name,"gene_name":gene_name,"alteration_class":"LOH","CCF":[ccf]})],
						ignore_index=True
					)
			# B. log the events
			# need to isolate somatic SNV GAIN events
			somatic_snv_events = node.somatic_snv_events.values()
			n_somatic_snv_GAINS = len([i for i in somatic_snv_events if i[1] == "GAIN"])
			n_somatic_cnv_events = len(node.germline_snp_events) + len(somatic_snv_events) - n_somatic_snv_GAINS
			c["total_snvs"] += n_somatic_snv_GAINS
			c["total_cnvs"] += n_somatic_cnv_events
		tree_log_df.loc[patient_name, ["total_snvs", "total_cnvs"]] = [c["total_snvs"], c["total_cnvs"]]

	master_maf_dict[patient_name] = maf
# %%
master_maf = pd.concat(master_maf_dict.values(), ignore_index=True)
master_maf.to_csv("1B_pan_cohort_scMAF.LOH_events.csv", index=False)
# %%
tree_log_df["truncal_snv_density"] = tree_log_df["truncal_snvs"] / tree_log_df["total_snvs"]
tree_log_df["truncal_snp_density"] = tree_log_df["truncal_cnvs"] / tree_log_df["total_cnvs"]
mean_truncal_snv_density = tree_log_df["truncal_snv_density"].mean()
print(f"mean truncal SNV density: {mean_truncal_snv_density}")
mean_truncal_snp_density = tree_log_df["truncal_snp_density"].mean()
print(f"mean truncal SNP density: {mean_truncal_snp_density}")

# %% t-test, if the two means are significantly different
from scipy.stats import ttest_ind
ttest_ind(tree_log_df["truncal_snv_density"].values, tree_log_df["truncal_snp_density"].values)
# %%