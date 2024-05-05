# %% 
import mosaic.io as mio
import pandas as pd
from pathlib import Path
from tea.snv.me import get_exclusivity
import matplotlib.pyplot as plt
import seaborn as sns

patient="RA17_22"
h5= "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/fillout_h5/RA17_22.patient_wide.genotyped.h5"

# %% 
s = mio.load(h5)
s.dna.genotype_variants()
# %%
voi = "chr12:25398284:C/T"
# which barcodes are positive?
boi = s.dna.barcodes()[s.dna.get_attribute("mut_filtered",features=[voi])[0].astype(bool)]
len(boi)

# total number of barcodes
total_sample_barcode_count = pd.Series([i.rsplit("-")[2] for i in s.dna.barcodes()]).value_counts()
# barcode carrying second KRAS
boi_count = pd.Series([i.rsplit("-")[2] for i in boi]).value_counts()

# looks like the only sample with >0.5% sc-prev of KRAS is 33_2, with 15 / 1810 barcodes positive
# %%
voi = ["chr12:25398284:C/A", "chr18:48573531:G/A", "chr18:48573570:G/A", "chr12:25398284:C/T"]
ss_barcodes = [i for i in s.dna.barcodes() if i.endswith("32_1") ]
soi = s[ss_barcodes]

g00 = 0
g01 = 1 
g10 = 1
g11 = 0

E2 = get_exclusivity(
    soi.dna.get_attribute('DP', constraint='row', features=voi),
    soi.dna.get_attribute('alt_read_count', constraint='row', features=voi),
    # ADO_PRECISION, FP,
    gametes = [(g00,g01), (g10,g11)],
    # rm_irrelevant_cells = False
)

E2 = 1 - E2
fig, ax = plt.subplots(figsize=(5,5))         # Sample figsize in inches
sns.heatmap(
    E2, annot=True, 
    cmap='Blues', fmt='.2f', # vmin=0, vmax=1,
    ax=ax
    )

# %%
