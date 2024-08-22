import csv
import os
from os.path import join as ospj
from numpy import savetxt
import pandas as pd
import re
import sys
import json
import shutil

"""
This script collates the pyAFQ tract profiles output across subjects in preparation for visualization and fitting GAMs.
Each pyAFQ tract profiles csv should be a 1D array of length 2800 (100 nodes * 28 tracts). It also identifies which subjects 
have missing tracts and saves out sub-ids in order to exclude them in final sample construction. 

**i might be able to delete all the subject folders and the individual level tract profiles**

Output:
- collated_tract_profiles.tsv for each dataset (rows = 100 nodes * 28 tracts * X subjects; columns = 4 metrics + tractID + nodeID)
- subs_missing_tracts.txt: 1D vector of subject IDs for people who were missing entire tracts
- subs_with_nans.txt: 1D vector of subject IDs for people who had NaN data (seems like there are no NaNs in anyone's data)
"""
 
for dataset in ["HCPD", "HBN", "PNC"]:  
    print(f"Processing {dataset}")
    ########################################
    # Set directories
    ########################################
    config_file = f"/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_{dataset}.json"
    with open(config_file, "rb") as f:
        config = json.load(f)
    data_root = config['tract_profiles_root']
    if not os.path.exists(ospj(data_root, "all_subjects")):
        os.makedirs(ospj(data_root, "all_subjects"))
        print(f"all_subjects created")
    else:
        print(f"all_subjects already exists.")
    sample_selection_dir = ospj(config['data_root'], "sample_selection_files")

    ########################################
    # Reformat tract profiles
    ########################################
    # load each subject's csv. make new dataframes for dki/dti fa and md separately
    # that dataframe should have column names as tractID_nodeID and the number should correspond to fa or md
    if dataset == "PNC":
        all_tract_profiles = pd.DataFrame(columns=['sub', 'tractID', 'nodeID', 'dti_fa', 'dti_md'])
    else: 
        all_tract_profiles = pd.DataFrame(columns=['sub', 'tractID', 'nodeID', 'dki_fa', 'dki_md','dti_fa', 'dti_md'])

    # also, keep track of subjects who are missing data
    subs_missing_tracts = []
    subs_with_nans = {}

    for sub in os.listdir(data_root):
        sub_dir = os.path.join(data_root, sub)
        if os.path.isdir(sub_dir) and sub != "all_subjects":
            for filename in os.listdir(sub_dir):
                if filename.endswith("profiles_dwi.csv"):
                    filepath = os.path.join(sub_dir, filename)
                    df = pd.read_csv(filepath)
                    # check for rows with NaNs
                    nan_rows = df[df.isna().any(axis=1)]
                    if not nan_rows.empty:
                        print(f"NaNs {sub}")
                        subs_with_nans[sub] = nan_rows.index.tolist()
                    if len(df) < 2800:
                        print(f"Missing tracts across all measures {sub}")
                        subs_missing_tracts.append(sub) # subjects with missing data will NOT be included in the collated tsv
                    else:
                        # add a column for subject ID
                        df['sub'] = sub
                        # append the data to the all_tract_profiles df
                        all_tract_profiles = pd.concat([all_tract_profiles, df], ignore_index=True)
                    os.remove(filepath)
            # delete the entire subject directory to save storage
            shutil.rmtree(sub_dir)
             
    # save out collated tract profiles (subjects missing data are NOT included in this)
    all_tract_profiles = all_tract_profiles.iloc[:, :-1]
    all_tract_profiles.to_csv(f'{data_root}/all_subjects/collated_tract_profiles_nocovbat.tsv', index=False)
    
    # save the subject IDs with missing tracts to a text file
    with open(os.path.join(sample_selection_dir, 'subs_missing_tracts.txt'), 'w') as f:
        for subject in subs_missing_tracts:
            f.write(f"{subject}\n")

    # save the rows with NaNs for each subject to a text file
    with open(os.path.join(sample_selection_dir, 'subs_with_nans.txt'), 'w') as f:
        for subject, nan_rows in subs_with_nans.items():
            f.write(f"{subject}: {nan_rows}\n")

