# %%
# We want two things:
# (1) per-sample: MAF showing clonality and mean single-cell
#     columns : ['condensed_format', 'patient_name', 'sample_name', 'gene_name', 'snv_class', 'HGVSp', 'total_num_cells', 'mut_filtered_num_cells', 'mut_filtered_sc_freq', 'mut_unfiltered_num_cells', 'mut_unfiltered_sc_freq', 'mean_sc_AF_in_mut_filtered', 'median_sc_AF_in_mut_filtered', 'driver_snv', 'num_cancer_sc', 'num_cancer_sc_w_var', 'CCF', 'bulk_tumor_fillout_VAF', 'bulk_normal_VAF']
# (2) per sample: log sheet recording whether this procedure is completed or not

# io
import yaml
from pathlib import Path
import pandas as pd
import numpy as np
import mosaic.io as mio
import warnings

from IPython import embed

warnings.filterwarnings("ignore")

# SNV filtering
from tea.cravat import NONFUNC_SO, get_technical_artifact_mask
from tea.format import isNaN

# plotting
import matplotlib.pyplot as plt

data_dir = Path("../../data_compiled")
wd = Path(".")
wd.mkdir(exist_ok=True)

############### PARAMS ###############
mut_prev_threshold = 0.005
normals_pon_occurence = 3
func_vars = True
blacklist = None
######################################
master_patient_df = pd.read_excel("../../Tapestri_batch2_samples_MASTER.xlsx", sheet_name="all_case_genetics", skiprows=1, index_col=0)
cases_to_exclude = master_patient_df.index[master_patient_df["censor"] == 1].tolist()
all_patient_sample_map_f = "../Tapestri_all_patient_sample_map.yaml"
with open(all_patient_sample_map_f, "r") as f:
    patient_sample_map = yaml.safe_load(f)
# remove cases to exclude
for case in cases_to_exclude:
    if case in patient_sample_map:
        patient_sample_map.pop(case)
sample_mafs = []
processed_cases = []
# %%
for patient_i in patient_sample_map:
    # patient_i = "BPA-2"
    if patient_i in processed_cases:
        print(f"[INFO] patient {patient_i} already processed! Skipping...")
        continue
    print(f"[INFO] processing patient {patient_i}")
    # if count >2:
    #     break
    # get H5 and SNV list
    mpileup_h5 = list((data_dir / "fillout_h5").glob(f"{patient_i}*.h5"))
    if len(mpileup_h5) == 0:
        raise ValueError(
            f"[ERROR] patient {patient_i}'s H5 file does not exist! Skipping..."
        )
    if len(mpileup_h5) > 1:
        raise ValueError(f"[ERROR] patient {patient_i} has multiple H5 files! Skipping...")
    snv_list_f = list(
        (data_dir / "manual_annotated_snv_lists").glob(f"{patient_i}*voi.hz_curated.txt")
    )
    if len(snv_list_f) > 1:
        raise ValueError(f"[ERROR] patient {patient_i} has multiple SNV lists! Skipping...")
    snv_list = pd.read_csv(snv_list_f[0], sep="\t", index_col=0)
    snv_list = snv_list[
        ~(
            snv_list["annotation"].str.contains("germline")
            | snv_list["annotation"].str.contains("artifact")
        )
    ]
    snv_ann_map = snv_list["HGVSp"].to_dict()
    voi = snv_list.index.tolist()
    cravat_f = list(
        (data_dir / "fillout_h5").glob(f"{patient_i}*CRAVAT_output_cleaned.txt")
    )
    if len(snv_list_f) > 1:
        raise ValueError(f"[ERROR] patient {patient_i} has multiple SNV lists! Skipping...")
    cravat_df = pd.read_csv(cravat_f[0], sep="\t", index_col=0, header=[0, 1]).loc[voi]
    # load and filter H5
    sample_obj = mio.load(str(mpileup_h5[0]))
    sample_obj.dna = sample_obj.dna[sample_obj.dna.barcodes(), voi]
    sample_obj.dna.genotype_variants(
        min_dp=8,
        min_alt_read=3,
        assign_low_conf_genotype=True,
    )
    num_cells = sample_obj.dna.shape[0]

    # ===== variant filtering =====
    # functional
    # NONFUNC_SO.append('=')
    if func_vars:
        # filter out snvs with non-functional SO in the value of snv_ann_map
        voi = cravat_df.index[
            cravat_df[("Variant Annotation", "Sequence Ontology")].apply(
                lambda x: x not in NONFUNC_SO
            )
        ].tolist()

    # fetch raw layers
    mut_filtered_layer = sample_obj.dna.get_attribute(
        "mut_filtered", features=voi, constraint="row"
    )
    mut_unfiltered_layer = sample_obj.dna.get_attribute(
        "mut_unfiltered", features=voi, constraint="row"
    )
    # get AF in cells that are mutated for each variant
    af_mut_filtered_layer = sample_obj.dna.get_attribute(
        "AF", features=voi, constraint="row"
    )[mut_filtered_layer.astype(bool)]
    af_mut_unfiltered_layer = sample_obj.dna.get_attribute(
        "AF", features=voi, constraint="row"
    )[mut_unfiltered_layer.astype(bool)]

    # create sample MAF
    sample_maf = pd.DataFrame(
        columns=[
            "condensed_format",
            "patient_name",
            "sample_name",
            "gene_name",
            "snv_class",
            "HGVSp",
            "total_num_cells",
            "mut_filtered_num_cells",
            "mut_filtered_sc_freq",
            "mut_unfiltered_num_cells",
            "mut_unfiltered_sc_freq",
            "mean_sc_AF_in_mut_filtered",
            "median_sc_AF_in_mut_filtered",
            "driver_snv",
            "num_cancer_sc",
            "num_cancer_sc_w_var",
            "CCF",
            "bulk_tumor_fillout_VAF",
            "bulk_normal_VAF",
        ]
    )
    sample_maf["condensed_format"] = voi
    sample_maf["patient_name"] = patient_i
    sample_maf["sample_name"] = patient_i

    sample_maf[["gene_name", "snv_class", "HGVSp"]] = cravat_df.loc[
        voi,
        [
            ("Variant Annotation", "Gene"),
            ("Variant Annotation", "Sequence Ontology"),
            ("Variant Annotation", "Protein Change"),
        ],
    ].values

    sample_maf["total_num_cells"] = sample_obj.dna.shape[0]
    sample_maf[
        [
            "mut_filtered_num_cells",
            "mut_unfiltered_num_cells",
            "mean_sc_AF_in_mut_filtered",
            "median_sc_AF_in_mut_filtered",
        ]
    ] = pd.concat(
        [
            mut_filtered_layer.sum(axis=0),
            mut_unfiltered_layer.sum(axis=0),
            af_mut_filtered_layer.mean(axis=0, skipna=True),
            af_mut_filtered_layer.median(axis=0, skipna=True),
        ],
        axis=1,
    ).values

    sample_maf[["mut_filtered_sc_freq", "mut_unfiltered_sc_freq"]] = sample_maf[
        ["mut_filtered_num_cells", "mut_unfiltered_num_cells"]
    ].divide(sample_maf["total_num_cells"], axis=0)

    # # %%
    # case BPA-3 needs to be considered specifically, as we need to define cancer cells by CN calling result (CN clone NOT 3 from the final ETE tree)
    if patient_i == "BPA-3":
        # need to use FALCON clone to define cancer cells
        cn_assignment_f = list(Path(".").glob(f"{patient_i}*assignment*csv"))
        if len(cn_assignment_f) != 1:
            raise ValueError(f"[ERROR] for patient {patient_i}, CN assignment file not found")
        cn_assignment_df = pd.read_csv(cn_assignment_f[0]).sort_values(["sample_name", "cell_barcode"])
        # cn_assignment_df["sample_name"] = cn_assignment_df["sample_name"].map({"BPA-3-RSX": "BPA-3-T", "BPA-3-N": "RSX640_N"})
        cn_assignment_df["barcode-sample"] = cn_assignment_df["cell_barcode"] + "-" + cn_assignment_df["sample_name"]
        tumor_cells = cn_assignment_df.loc[cn_assignment_df["final_clone_id"] != 3, "barcode-sample"].values
        tumor_cells_bool = mut_filtered_layer.index.isin(tumor_cells)
        sample_maf["driver_snv"] = "refined_tree_clone3"
        # embed()
    # case BPA-4 needs to be considered specifically, as we need to define cancer cells by CN calling result (CN clone not 3 from the final ETE tree)
    elif "driver_snv" in patient_sample_map[patient_i]:
        driver_snv = patient_sample_map[patient_i]["driver_snv"]
        sample_maf["driver_snv"] = driver_snv
        if driver_snv in mut_filtered_layer.columns:
            tumor_cells_bool = mut_filtered_layer[driver_snv].astype(bool)
        else:
            raise KeyError(
                    f"[ERROR] for patient {patient_i}, driver SNV not found in mut_filtered_layer"
                )
    else:
       raise KeyError("[ERROR] for patient {patient_i}, driver SNV not available")
    mut_filtered_layer_cancer_cells = mut_filtered_layer[tumor_cells_bool]
    num_cancer_sc = mut_filtered_layer_cancer_cells.shape[0]
    num_cancer_sc_w_var = mut_filtered_layer_cancer_cells.sum(axis=0)
    sample_maf["num_cancer_sc"] = num_cancer_sc
    sample_maf["num_cancer_sc_w_var"] = num_cancer_sc_w_var.values
    ccf_vals = (num_cancer_sc_w_var / num_cancer_sc).values
    sample_maf["CCF"] = ccf_vals
    # bulk coverage
    try:
        sample_maf["bulk_tumor_fillout_VAF"] = cravat_df.loc[
            voi, [("bulk_comparison", "bulk-matched_bulk_cohort-AF")]
        ].values
    except:
        print(f"[ERROR] for sample {patient_i}, bulk tumor fillout VAF not available")
        e = 1
    try:
        sample_maf["bulk_normal_VAF"] = cravat_df.loc[
            voi, [("bulk_comparison", "bulk-matched_bulk_normal-AF")]
        ].values
    except:
        print(f"[ERROR] for sample {patient_i}, bulk normal VAF not available")
        e = 1

    sample_mafs.append(sample_maf)
    processed_cases.append(patient_i)

# %%
pan_cohort_scMAF = pd.concat(sample_mafs, ignore_index=True)

pan_cohort_scMAF.to_csv(wd / f"pan_cohort.scMAF.csv", index=False, header=True)



# %%
