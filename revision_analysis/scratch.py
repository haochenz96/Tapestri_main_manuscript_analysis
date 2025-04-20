# %%
import pickle

f = "/Users/haochenzhang/Iacobuzio_lab/Tapestri_main_manuscript_analysis/0_condor_pipeline/condor_outputs/pickle_files/PC15_self.solT_cell"

# inspect condor SNV orders
with open(f, "rb") as f:
    nx_tree = pickle.load(f)
# %%
