# %% import things
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyparsing import line
import seaborn as sns
import os
from tea.cravat import NONFUNC_SO
import pickle
import glob
import mosaic.io as mio
# cd to the directory of the script
WD = Path(__file__).parent
os.chdir(WD)
ss_h5_dir = "../data_compiled/ss_h5"

pan_cohort_scMAF = pd.read_csv("pan_cohort_somatic_scMAF_mut_prev=0.005.csv")
case_sheet = pd.read_excel("../Tapestri_batch2_samples_MASTER_INTERNAL.xlsx", skiprows=1, index_col=0)
case_rename = case_sheet["HZ_specific_case_ID"].to_dict()
pan_cohort_scMAF["patient_name"] = pan_cohort_scMAF["patient_name"].map(case_rename)
sample_sheet = pd.read_excel("../Tapestri_batch2_samples_MASTER_INTERNAL.xlsx", index_col=1, sheet_name = "all_sample_clinical_info")
sample_rename = sample_sheet["HZ_specific_sample_ID"].to_dict()
# %% Focus on one case: PC20
h5s = glob.glob(os.path.join(ss_h5_dir, "PC20*.h5"))
samples = {}
for h5 in h5s:
	smp = mio.load(h5)
	smp.dna.genotype_variants()
	sample_name = smp.dna.metadata['sample_name'][0][0]
	sample_name = sample_rename[sample_name]
	# filter to SNVs that are in more than 1 cell
	voi = smp.dna.ids()[
		smp.dna.get_attribute("mut_filtered").sum(axis=0) > 1
    ]
	smp.dna = smp.dna[list(smp.dna.barcodes()), voi]
	samples[sample_name] = smp

# %%

# %%


# %%
