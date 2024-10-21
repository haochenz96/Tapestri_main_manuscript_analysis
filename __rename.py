# %% 
import os
import pandas as pd
import numpy as np
from pathlib import Path
import mosaic.io as mio
os.chdir('/Users/haochenzhang/Iacobuzio_lab/Tapestri_main_manuscript_analysis/')
def rename_in_directory(directory: str, original_name: str, mapped_name: str):
    # Walk through the directory
    for dirpath, dirnames, filenames in os.walk(directory, topdown=False):
        # Rename files
        for filename in filenames:
            if original_name in filename:
                old_file_path = os.path.join(dirpath, filename)
                new_file_name = filename.replace(original_name, mapped_name)
                new_file_path = os.path.join(dirpath, new_file_name)
                os.rename(old_file_path, new_file_path)
        
        # Rename directories
        for dirname in dirnames:
            if original_name in dirname:
                old_dir_path = os.path.join(dirpath, dirname)
                new_dir_name = dirname.replace(original_name, mapped_name)
                new_dir_path = os.path.join(dirpath, new_dir_name)
                os.rename(old_dir_path, new_dir_path)

    print(f"Renaming completed. All instances of '{original_name}' have been replaced with '{mapped_name}' in '{directory}'")

# Usage example
# rename_in_directory('/path/to/directory', 'old_name', 'new_name')
# %% load renaming map
patient_rename = pd.read_excel('Tapestri_batch2_samples_MASTER_INTERNAL.xlsx', sheet_name = "all_case_genetics", skiprows=1)
patient_rename_map = dict(zip(patient_rename['case_ID'], patient_rename['HZ_specific_case_ID']))
sample_rename = pd.read_excel('Tapestri_batch2_samples_MASTER_INTERNAL.xlsx', sheet_name='all_sample_clinical_info')
sample_rename_map = dict(zip(sample_rename['sample'], sample_rename['HZ_specific_sample_ID']))

# %%
os.chdir('/Users/haochenzhang/Iacobuzio_lab/Tapestri_main_manuscript_analysis/')
for original_name, mapped_name in sample_rename_map.items():
    rename_in_directory(
        'data_compiled/ss_h5', 
        original_name,
        mapped_name
    )
# %% ==== rename all h5 metadata =====
# @HZ 2024-10-20: NOT WORKING!!!!
data_compiled_dir = Path('data_compiled')
(data_compiled_dir / 'fillout_h5_renamed').mkdir(exist_ok=True)

renanmed_data_dir = Path()
count = 0
for patient_name in patient_rename_map.values():
    h5_file_path = data_compiled_dir / 'fillout_h5' / f'{patient_name}.patient_wide.genotyped.h5'
    renamed_h5_file_path = data_compiled_dir / 'fillout_h5_renamed' / f'{patient_name}.patient_wide.genotyped.renamed.h5'
    if os.path.exists(h5_file_path):
        print(f"found {h5_file_path}")
    else:
        raise ValueError(f"h5 file not found for {patient_name}")
    
    s = mio.load(h5_file_path)
    os.remove(h5_file_path) # remove original h5 file
    og_sample_names = s.dna.metadata["sample_name"]
    s.dna.add_metadata(
        "sample_name", 
        np.array(
            [[sample_rename_map[og_sample_names[i][0]]] for i in range(len(og_sample_names))]
        ).astype(object)
    )
    s.cnv.add_metadata(
        "sample_name", 
        s.dna.metadata["sample_name"]
    )

    # save
    mio.save(s, renamed_h5_file_path)

    count += 1
    if count > 0:
        break
