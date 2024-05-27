# %% 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from tea.plots import mpl_config_params
mpl_config_params(font_size=12)

case = "RA17_13"
os.makedirs(case, exist_ok=True)

# %%
# read BAF
baf_f = f"{case}/{case}_roi_wgs_snp_baf.bed"
baf = pd.read_csv(baf_f, sep="\t", index_col=False, names=["roi_chr", "roi_start", "roi_end", "ann", "chr", "start", "end", "sample", "AAC", "BAC"])

# calculate BAF
baf["BAF"] = baf["BAC"] / (baf["AAC"] + baf["BAC"])

# plot BAF
goi = baf["ann"].unique()

# # RA15_16
# soi = ["MPAM02PT3"]

# RA17_13
soi = ["MPAM05PT3", "MPAM05PT6", "MPAM05PT8", "MPAM05PT9", "MPAM05PT11"]

# # RA17_22
# soi = ["MPAM06PT8"]


for g in goi:
	for s in soi:
		fig, ax = plt.subplots(figsize=(8, 3))
		sns.scatterplot(
			data=baf[
				(baf["ann"] == g) & 
				(baf["sample"] == s)
				],
			x="start", 
			y="BAF", 
			ax=ax,
			color="blue", s=15,
			# outline color
			linewidth=0.5, edgecolor="black"
			)
		ax.set(
			title=f"SNP BAF in {g} region",
			ylim=[-0.1, 1.1],
			xlabel="Genomic coordinate",
			ylabel="B allele frequency"
			)
		fig.savefig(f"{case}/{case}_BAF_{s}_{g}.pdf", bbox_inches="tight")
# %%
