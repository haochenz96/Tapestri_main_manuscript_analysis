# %% 

import mosaic.io as mio

h5 = "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled/fillout_h5/BPA-5-RSX.mpileup.renamed.h5"
pt = mio.load(h5)

# %% 
pt.dna.genotype_variants()
# %%
# get the barcodes that are mutated for KRAS
soi = ["chr18:48575152:C/T", "chr12:25398284:C/T"]
boi = pt.dna.barcodes()[
	pt.dna.get_attribute("mut_filtered", features=soi, constraint="row").any(axis=1)
]
len(boi)
sorted([int(b.split("-BPA")[0].split("cell_")[1]) for b in boi ])

# %%
