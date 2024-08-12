import csv
import glob
import json
import nibabel as nib
import numpy as np
import os
from os.path import join as ospj
import pandas as pd
import re
import shutil
import sys
 

# This script collates tract profiles across subjects into 1 large csv

########################################
# Set directories
########################################
config_file = "/cbica/projects/luo_wm_dev/code/superficialWM_analyses/config_HBN.json"

# Read config from the specified file
with open(config_file, "rb") as f:
    config = json.load(f)

dataset = config['dataset']
data_root = config['tract_profiles_data_root']
outputs_root = config['tract_profiles_outputs_root']

 
if not os.path.exists(ospj(data_root, "all_subjects")):
    os.makedirs(ospj(data_root, "all_subjects"))
    print(f"all_subjects created")
else:
    print(f"all_subjects already exists.")


if not os.path.exists(outputs_root):
    os.makedirs(outputs_root)
    print(f"{outputs_root} created")
else:
    print(f"{outputs_root}  already exists.")

  
 
# Initialize a DataFrame to store the aggregated FA and MD values for each node and tract
all_tract_profiles = pd.DataFrame(columns=['subjectID', 'tractID', 'nodeID', 'dti_fa', 'dti_md'])
subs_missing_tracts = []
# Iterate through each subject's directory
for sub in os.listdir(data_root):
    sub_dir = os.path.join(data_root, sub)
    if os.path.isdir(sub_dir):
        # Iterate through each subject's CSV files
        for filename in os.listdir(sub_dir):
            if filename.endswith("profiles_dwi.csv"):
                filepath = os.path.join(sub_dir, filename)
                # Read the CSV file
                df = pd.read_csv(filepath)
                if len(df) < 2200: # if missing tracts 
                    print(sub) 
                    subs_missing_tracts.append(sub)   
                else:
                    # Add a column for the subject ID
                    df['subjectID'] = sub
                    # Append the data to the all_tract_profiles DataFrame
                    all_tract_profiles = pd.concat([all_tract_profiles, df], ignore_index=True)

 
if dataset == "HCPD":
    print('HCPD')
    sub_0197045 = pd.read_csv("/cbica/projects/luo_wm_dev/input/HCPD/HCPD_tractprofiles/sub-0197045/sub-0197045_ses-V1_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-profiles_dwi.csv")
    sub_0197045['tractID'].unique()
    sub_0001305 = pd.read_csv("/cbica/projects/luo_wm_dev/input/HCPD/HCPD_tractprofiles/sub-0001305/sub-0001305_ses-V1_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-profiles_dwi.csv")
    sub_0001305['tractID'].unique()
    np.setdiff1d(sub_0001305['tractID'].unique(), sub_0197045['tractID'].unique()) # sub_0197045 missing "FP"

    # identify where FP should be
    fp_indices = sub_0001305[sub_0001305['tractID'] == 'FP'].index

    # create a DataFrame with NaN values
    nan_df = pd.DataFrame(np.nan, index=range(len(fp_indices)), columns=sub_0197045.columns)
    nan_df['tractID'] = "FP"
    nan_df['nodeID'] = range(0,100)

    # insert the NaN DataFrame into df1 at the specified indices
    sub_0197045 = pd.concat([sub_0197045.iloc[:fp_indices[0], :], nan_df, sub_0197045.iloc[fp_indices[0]:, :]]).reset_index(drop=True)
    sub_0197045['subjectID'] = "sub-0197045"
    all_tract_profiles = pd.concat([all_tract_profiles, sub_0197045], ignore_index=True)
    all_tract_profiles = all_tract_profiles.iloc[:, :-1]
    all_tract_profiles.to_csv(f'{data_root}/all_subjects/collated_tract_profiles.tsv', index=False)


 

if dataset == "HBN":
    sub_0001305 = pd.read_csv("/cbica/projects/luo_wm_dev/input/HCPD/HCPD_tractprofiles/sub-0001305/sub-0001305_ses-V1_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-profiles_dwi.csv")
    all_tracts = sub_0001305['tractID'].unique()
    rm_subs = []
    for idx, sub in enumerate(subs_missing_tracts): 
        sub_csv = pd.read_csv(glob.glob(f"/cbica/projects/luo_wm_dev/input/{dataset}/{dataset}_tractprofiles/{sub}/{sub}_ses-HBNsite*_acq-64dir_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-profiles_dwi.csv")[0])
        
        # Extract subject ID
        subjectID = sub
        
        # Extract unique tract IDs
        unique_tractIDs = sub_csv['tractID'].unique()
        
        print(subjectID)
        print(unique_tractIDs)
        print(len(unique_tractIDs))
       
        # if missing more than 1 tract (i.e. has less than 21 tracts), then remove
        if len(unique_tractIDs) < 21: # if missing >1 tract 
            rm_subs.append(sub) 

 

    # Which tracts are the subjects missing? 
    sub_missing_tract_list = []
    for idx, sub in enumerate(subs_missing_tracts): 
        print(sub)
        sub_csv = pd.read_csv(glob.glob(f"/cbica/projects/luo_wm_dev/input/{dataset}/{dataset}_tractprofiles/{sub}/{sub}_ses-HBNsite*_acq-64dir_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-profiles_dwi.csv")[0])
        
        # Extract subject ID
        subjectID = sub
        
        # Extract unique tract IDs
        unique_tractIDs = sub_csv['tractID'].unique()
        
        missing_tract = set(all_tracts.tolist()) - set(unique_tractIDs)

         # Append data to the list
        sub_missing_tract_list.append({'subjectID': subjectID, 'missing_tract': missing_tract})

        
    sub_missing_tract_df = pd.DataFrame(sub_missing_tract_list)
    sub_missing_tract_df.to_csv(f'{data_root}/all_subjects/subjects_with_missing_tracts.csv', index=False)
    
    # I'm just going to exclude all subjects that are missing any tracts from HBN, and do the same for HCP-D
    all_tract_profiles = all_tract_profiles.iloc[:, :-1]
    all_tract_profiles.to_csv(f'{data_root}/all_subjects/collated_tract_profiles.tsv', index=False)
     