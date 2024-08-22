from collections import defaultdict
import csv
import h5py
import json
import nibabel as nib
import numpy as np
import os
from os.path import join as ospj
import pandas as pd
import re
import shutil
import sys
from tqdm import tqdm

# This script saves out tract profiles for all subjects to h5 file  
# script adapted from Chenying's concifti code :) 

########################################
# Set directories
########################################
config_file = "/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_HBN.json"

# Read config from the specified file
with open(config_file, "rb") as f:
    config = json.load(f)

dataset = config['dataset']
data_root = config['tract_profiles_data_root']
outputs_root = config['tract_profiles_outputs_root']
h5_dir = ospj(f"{outputs_root}", "h5_files")

if not os.path.exists(f"{h5_dir}"):
        os.makedirs(f"{h5_dir}")
        print(f"Directory 'h5_dir' created.")
else:
        print(f"Directory 'h5_dir' already exists.")

 
########################################
# Define functions
########################################

# write hdf5 files for csv's
def write_hdf5(scalar):
    """
    Load all gifti data.
    Parameters
    -----------
    scalar: str
        Name of scalar (e.g. dti_fa, or dti_md)
    """

    # define cohort filename and load cohort file
    cohort_file = f"{tract_profile_scalar}_cohortfile.csv"
    cohort_df = pd.read_csv(ospj(outputs_root, "cohortfiles", tract_profile_scalar, cohort_file))
            
    # upload each subject's data
    scalars = defaultdict(list)
    sources_lists = defaultdict(list)
    
    print("Loading tract profiles for each subject")
    for ix, row in tqdm(cohort_df.iterrows(), total=cohort_df.shape[0]):   # ix: index of row (start from 0); row: one row of data
        scalar_file = row['source_file']
        data = pd.read_csv(scalar_file, header=None)
        scalars[row['scalar_name']].append(data.T.to_numpy())   # append to specific scalar_name
        sources_lists[row['scalar_name']].append(row['source_file'])  # append source csv filename to specific scalar_name

    # Write the output
    output_h5 = ospj(f"{tract_profile_scalar}.h5")
    h5_finaldir = ospj(f"{h5_dir}")

    if not os.path.exists(f"{h5_finaldir}"):
            os.makedirs(f"{h5_finaldir}")
            print(f"Directory {h5_finaldir} created.")
    else:
            print(f"Directory {h5_finaldir} already exists.")

    output_file = ospj(h5_finaldir,output_h5)
    f = h5py.File(output_file, "w")

    for scalar_name in scalars.keys():  # in the cohort.csv, two or more scalars in one sheet is allowed, and they can be separated to different scalar group.
        one_scalar_h5 = f.create_dataset('scalars/{}/values'.format(scalar_name),
                         data=np.row_stack(scalars[scalar_name]))
         
        one_scalar_h5.attrs['column_names'] = list(sources_lists[scalar_name])  # column names: list of source csv filenames
    f.close()
    print(f"{tract_profile_scalar} h5 saved")
    return int(not os.path.exists(output_file))



##################################################
# Convert csv's to h5 files 
##################################################
   
# set tract profile scalars
tract_profile_scalars = ["dti_fa", "dti_md"]
for tract_profile_scalar in tract_profile_scalars:
    write_hdf5(tract_profile_scalar)