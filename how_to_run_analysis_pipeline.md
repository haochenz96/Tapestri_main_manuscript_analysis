<!-- 1. run supp1a to get sample-level library technical details -->
1. run supp1b to get technical aspects of general genetics, which will include the pan-cohort panel of normal and blacklist.
2. run 0_condor_pipeline
3. run 1_general_genetics
   1. `1A.1_get_scMAF.py` -- automatic. Run it to get `1A_pan_cohort_scMAF.SNVs.csv`, which includes short mutations required for making the oncoplot.
   2. `1B_pan_cohort_phylogeny_scMAF` -- automatic. Run it to get `1B_pan_cohort_phylogeny_scMAF.csv`, which includes all the LOH events needed for making the oncoplot.
   3. `1A.2_snv_cnv_combined_analysis.py` -- manual. Need to update patient/sample name to run and get sc-heatmaps as needed.
   4. `Figure1-oncoprint.r` -- automatic. Makes `1A_composite_maf.FINAL.csv`, `1A_composite_maf.TGFB.csv`, `1A_CCF_stats.csv`, `1A_ARID1A_CCF_stats.csv`. Most importantly, makes `Figure_1B_sc_oncoprint.pdf/png`
