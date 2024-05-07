# %% IMPORT MODULES
from pathlib import Path
import pandas as pd
import numpy as np

analysis_dir = Path("/data/iacobuzc/haochen/Tapestri_batch2/analysis")
snv_lists_dir = Path("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/manual_annotated_snv_lists")
# %% FILTER BLACKLIST
N=69
PON_OCCURENCE_FREQUENCY=0.5
threshold_sample_count = N * PON_OCCURENCE_FREQUENCY
tumors_pon_f = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1b_general_genetics/snv_black_whitelists/all_tumor_samples_N=69_composite_snv_sheet.csv"
tumors_pon_df = pd.read_csv(tumors_pon_f, sep=',', header=0, index_col=0)

snv_blacklist = tumors_pon_df[
    (tumors_pon_df["all_tumor_samples_occurence"] > threshold_sample_count)
].index.tolist()
# do not include KRAS
snv_blacklist = [x for x in snv_blacklist if not x.startswith('chr12:2539828')]

snp_blacklist = [
    'chr18:48584855:AT/A', # SMAD4 intron
    'chr18:48584855:A/ATT', # SMAD4 intron
    'chr11:92577070:T/C', # FAT3 intron
    'chr14:75508280:GA/G', # MLH3 intron,
    'chr10:88679200:C/T',# BMPR1A synonymous
    'chr4:96035796:GAC/G', # BMPR1B intron M12
    'chrX:53565685:A/ACACTCTCT', # HUWE1 intron RA19_10
    'chr4:96073753:AAAT/A', # BMPR1B intron RA19_10, RA19_21
    'chr7:154760608:TCTG/T', # PAXP1 Q434del RA19_21
    "chr9:101911436:A/AT", # TGFBR1 intron RA21_17
    "chr11:64572557:A/G", # MEN1 silent (noisy amplicon/LOH pattern)
    "chr2:203384776:CT/C", # BMPR2 from TP12 (noisy amplicon/LOH pattern)
    "chr11:64572557:A/G", # MEN1 His433 (noisy amplicon/LOH pattern)
    "chr14:105242966:T/C", # AKT1 intron (noisy amplicon/LOH pattern)
]
snp_loci_blacklist = [x.rsplit(":", 1)[0] for x in snp_blacklist]

manual_snv_blacklist = [
    'chr5:38953665:GA/G',
    'chr5:38953665:G/GA',
    'chr4:187538942:T/A',
    'chr4:187538942:T/G',
    'chr7:151926977:A/G',
    'chr11:64572557:A/G',
    'chr4:96035796:GAC/G', # BMPR1B intron
    'chr18:48584855:A/AT',
    'chr18:48584855:A/ATT', # SMAD4 intron
    'chr11:92577070:T/C', # FAT3 intron
    'chr18:45394653:G/A', # SMAD2 intron # for M11 only: somatic --> germline
    'chr5:68667298:G/T', # RAD17 5_prime_UTR
    'chr17:37881703:G/T', # ERBB2 intron variant RA20_05
    'chr2:241660378:G/A', # KIF1A p.Ser1615= RA20_05
    'chr3:38523684:C/T', # ACVR2B intron_variant
    "chr6:52149559:C/A", # MCM3 2kb_upstream
    "chr7:151927025:A/G", # KMT2C Y987H
    "chr7:151926945:G/A", # KMT2C intron
    "chr8:30933848:TA/T", # WRN intron
    "chr2:203379681:A/C", # BMPR2 L200=
    "chr5:68667280:G/T", # RAD17 from TP6
    "chr11:65631240:A/T", # MUS81 from TP6
]
snv_blacklist += manual_snv_blacklist
snv_loci_blacklist = [x.rsplit(":", 1)[0] for x in snv_blacklist]

# # %% Blacklist the SNPs and SNVs
# # patients_of_interest = ["BPA-5-RSX"]
# output_dir = Path("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/supp1b_general_genetics/all_vars_mut_prev=0.01/sc_heatmaps")
# output_dir.mkdir(exist_ok=True, parents=True)
# for p in snv_lists_dir.glob("**/*-voi.hz_curated.txt"):
#     # if not any([x in str(p) for x in patients_of_interest]):
#     #     continue
#     manual_snv = pd.read_csv(p, sep='\t', index_col=0, comment = '#')
#     manual_snv.fillna({'annotation': ''}, inplace=True)
    
#     # flag all loci in blacklist as "likely_artifact"
#     count = 0
#     blacklisted = []
#     for locus, row in manual_snv.iterrows():
#         if "germline" in row["annotation"]:
#             if locus.rsplit(":", 1)[0] in snp_loci_blacklist:
#                 if not manual_snv.loc[locus, 'annotation'] == "likely_artifact":
#                     manual_snv.loc[locus, 'annotation'] = 'likely_artifact'
#                     count += 1
#                     # blacklisted.append(manual_snv.loc[locus, 'HGVSp'])
#                     blacklisted.append(locus)
#         else:
#             if locus.rsplit(":", 1)[0] in snv_loci_blacklist:
#                 if not manual_snv.loc[locus, 'annotation'] == "likely_artifact":
#                     manual_snv.loc[locus, 'annotation'] = 'likely_artifact'
#                     count += 1
#                     # blacklisted.append(manual_snv.loc[locus, 'HGVSp'])
#                     blacklisted.append(locus)
#     print(f"""[INFO] Blacklisted {count} SNVs/SNPs in {p.stem}: 
#           {blacklisted}""")
#     # write output 
#     output_f = output_dir / (p.stem)
#     manual_snv.to_csv(output_f, sep='\t', index=True)


# %% Load data and plot heatmaps
import yaml
import mosaic.io as mio
import plotly.express as px
from tea.plots import plot_snv_clone

snv_analysis_dir = analysis_dir / "snv-cnv-combined"

h5_dir = Path('/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/fillout_h5')
opt_nclones_dir = Path("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/data_compiled/falcon_solutions")
output_dir = Path("/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/0_condor_pipeline/pre_condor_sc_heatmaps")
output_dir.mkdir(exist_ok=True, parents=True)
patient_info_f = "/lila/data/iacobuzc/haochen/Tapestri_main_manuscript_analysis/1_general_genetics/Tapestri_all_patient_sample_map.yaml"
with open(patient_info_f, 'r') as f:
    patient_info = yaml.safe_load(f)
patient_names = patient_info.keys()
patient_names = ["M13", "BPA-3"]
for patient_name in patient_names:
    # if (output_dir / f"{patient_name}-DNA-heatmap.pdf").is_file():
    #     print(f'[INFO] Skipping {patient_name}.')
    #     continue
    print(f'[INFO] Processing {patient_name}.')

    try:
        pt_h5 = list(h5_dir.glob(f"{patient_name}*h5"))[0]
    except:
        print(f'[ERROR] No H5 file found for {patient_name}.')
        continue
    pt = mio.load(pt_h5)

    # Add raw FALCON results to H5
    cn_assignment_f = list(opt_nclones_dir.glob(f"{patient_name}*assignment*csv"))[0]
    cn_assignment_df = pd.read_csv(cn_assignment_f, index_col = 0)
    print(f'[INFO] Loaded CN clone assignment file {cn_assignment_f}.')

    cn_clone_palette = dict(zip(
        np.sort(cn_assignment_df['clone_id'].unique()), 
        np.array(px.colors.qualitative.Set3)[np.sort(cn_assignment_df['clone_id'].unique())]
        ))
    # rename the keys
    cn_clone_palette = {f"CN_clone-{k}": v for k, v in cn_clone_palette.items()}

    # add cn_clone info
    cn_assignment_df['cell_barcode_formatted'] = cn_assignment_df['cell_barcode'] + "-" + cn_assignment_df.index
    cn_assignment_df.set_index('cell_barcode_formatted', inplace=True)
    pt.dna.row_attrs['label'] = np.array(list(
        map(
            lambda x: f"CN_clone-{int(cn_assignment_df.loc[x, 'clone_id'])}", 
            pt.dna.barcodes())
        ))
    pt.dna.set_palette(cn_clone_palette)
    pt.cnv.row_attrs['label'] = pt.dna.row_attrs['label']
    pt.cnv.set_palette(cn_clone_palette)
    # pt.cnv.get_gene_names()

    num_cells = pt.dna.shape[0]


    # read in high-quality SNVs
    snv_f = list(snv_lists_dir.glob(f"{patient_name}*voi*.txt"))[0]
    snv_df = pd.read_csv(snv_f, sep = '\t', index_col = 0, comment='#')
    snv_df["annotation"].fillna("", inplace=True)
    snv_ann_map = snv_df['HGVSp'].to_dict()

    # ===== highlight vars with bulk annotation ====='
    # highlight vars
    germline_hom_var_col = '#00cc66' # green
    germline_het_var_col = '#2693ff' # blue
    somatic_var_col = '#ff0000' # red
    likely_artifact_col = '#ffcc00' # yellow

    for var_i in snv_ann_map:
        # germline
        if snv_df.loc[var_i, 'annotation'] == 'germline_HOM':
            snv_ann_map[var_i] = f'<span style="color:{germline_hom_var_col};">' + snv_ann_map[var_i] + '</span>'
        elif snv_df.loc[var_i, 'annotation'] == 'germline_HET':
            snv_ann_map[var_i] = f'<span style="color:{germline_het_var_col};">' + snv_ann_map[var_i] + '</span>'
        elif "somatic" in snv_df.loc[var_i, 'annotation']:
            snv_ann_map[var_i] = f'<span style="color:{somatic_var_col};">' + snv_ann_map[var_i] + '</span>'
        elif "artifact" in snv_df.loc[var_i, 'annotation']:
            snv_ann_map[var_i] = f'<span style="color:{likely_artifact_col};">' + snv_ann_map[var_i] + '</span>'
        else:
            pass

    # plot heatmap
    fig = plot_snv_clone(
        pt,
        sample_name=patient_name,
        story_topic = f'{patient_name}-high_conf_snvs',
        voi = snv_df.index.tolist(),
        attribute = "AF_MISSING",
        ann_map = snv_ann_map
    )
    # # map the x-axis ticks according to snv_ann_map
    # fig.update_xaxes(ticktext = list(map(lambda x: snv_ann_map[x], snv_df.index.tolist())))

    # save the figure
    # write to PDF, with width proportion to the number of SNVs
    fig.write_image(
        str(output_dir / f"{patient_name}-DNA-heatmap.pdf"), 
        width = 500 + 10 * len(snv_df.index.tolist())
    )
    print(f'[INFO] Saved heatmap for {patient_name} to {output_dir / f"{patient_name}-DNA-heatmap.pdf"}.')



# %%
