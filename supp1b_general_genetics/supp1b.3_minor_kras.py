# %% 
import mosaic.io as mio
import pandas as pd
from pathlib import Path
from tea.snv.me import get_exclusivity
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

def config_params(font_size=7):

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'

patient="RA17_22"
h5= "../data_compiled/fillout_h5/RA17_22.patient_wide.genotyped.h5"

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
ann_map = {
    "chr12:25398284:C/A": "KRASp.G12V", 
    "chr12:25398284:C/T": "KRASp.G12D", 
    "chr18:48573531:G/A": "SMAD4p.A39T", 
    "chr18:48573570:G/A": "SMAD4p.D52N",
}
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
<<<<<<<< HEAD:supp1b_general_genetics/supp1b.2_minor_kras.py

# %%
config_params(font_size=8)

E2.columns = [ann_map[i] for i in E2.columns]
E2.index = [ann_map[i] for i in E2.index]
fig, ax = plt.subplots(figsize=(3,3))         # Sample figsize in inches
========
# %%
config_params(font_size=8)
fig, ax = plt.subplots(figsize=(5,5))         # Sample figsize in inches
>>>>>>>> a375bcd71c01471628d58fa77ad5e29ce33a8bd0:supp1b_general_genetics/supp1b.3_minor_kras.py
sns.heatmap(
    E2, annot=True, 
    cmap='Blues', fmt='.2f', # vmin=0, vmax=1,
    ax=ax
    )

# %%
<<<<<<<< HEAD:supp1b_general_genetics/supp1b.2_minor_kras.py
fig.savefig("supp1b.2_minor_kras.pdf", bbox_inches='tight')

# %%
========
fig.write()
>>>>>>>> a375bcd71c01471628d58fa77ad5e29ce33a8bd0:supp1b_general_genetics/supp1b.3_minor_kras.py
