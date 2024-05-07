# %% 
import mosaic.io as mio

ORG_BARCODE="org_barcode"
BARCODE="barcode"
SAMPLE="sample_name"
# %% 
# per https://github.com/haochenz96/mosaic/blob/6217b39074d68ded75a60496a4033f6b8f271ce0/src/mosaic/assay.py#L220
h5 = "/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/SNV_FILLOUT_results/BPA-5-RSX/BPA-5-RSX.mpileup.h5"
s = mio.load(h5)
# %% 
s.dna.add_row_attr(ORG_BARCODE, s.dna.barcodes())
s.dna.add_row_attr(BARCODE, s.dna.barcodes() + '-' + s.dna.row_attrs[SAMPLE])
s.cnv.add_row_attr(ORG_BARCODE, s.dna.barcodes())
# %%
mio.save(s,"/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/SNV_FILLOUT_results/BPA-5-RSX/BPA-5-RSX.mpileup.renamed.h5")
# %%
