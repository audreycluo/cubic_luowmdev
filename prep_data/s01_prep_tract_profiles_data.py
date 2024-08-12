import csv
import os
from os.path import join as ospj
from numpy import savetxt
import pandas as pd
import re
import sys
import json


# This script reformats the pyAFQ tract profiles output for each subject to have separate csv's for dki_fa, dki_md, dti_fa, dti_md in preparation for fitting GAMs.
# Tract profiles are modified to have 1 entry per unique tractID/nodeID combination (i.e. ARC_L_1) 
# Each csv is a 1D array of length 2800 (100 nodes * 28 tracts)
 

for dataset in ["HCPD", "HBN"]: 
    print(f"Processing {dataset}")
    ########################################
    # Set directories
    ########################################
    config_file = f"/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_{dataset}.json"

    # Read config from the specified file
    with open(config_file, "rb") as f:
        config = json.load(f)

    #dataset = config['dataset']
    data_root = config['tract_profiles_root']
    sample_selection_dir = ospj(config['data_root'], "sample_selection_files")
    
    demographics = config['demographics_temp']
    qc = config['qc_temp']

    
    ########################################
    # Reformat tract profiles
    ########################################

    # load each subject's csv. make a new dataframe for fa and md separately
    # that dataframe should have column names as tractID_nodeID and the number should correspond to fa or md
    # Initialize an empty list to store subject IDs with missing tracts
    subs_missing_tracts = []

    # Initialize a dictionary to store the rows with NaNs for each subject
    subs_with_nans = {}

    for sub in os.listdir(data_root):
        sub_dir = os.path.join(data_root, sub)
        if os.path.isdir(sub_dir):
            # Iterate through each subject's tract profile
            for filename in os.listdir(sub_dir):
                if filename.endswith("profiles_dwi.csv"):
                    filepath = os.path.join(sub_dir, filename)
                    # Read the csv file
                    df = pd.read_csv(filepath)
                    
                    # Check for rows with NaNs
                    nan_rows = df[df.isna().any(axis=1)]
                    if not nan_rows.empty:
                        print(f"NaNs {sub}")
                        subs_with_nans[sub] = nan_rows.index.tolist()
                    
                    if len(df) < 2800:
                        print(f"Missing tracts across all measures {sub}")
                        subs_missing_tracts.append(sub)
                    else:
                        # Create a new column tract_node by combining tractID and nodeID
                        df['tract_node'] = df['tractID'].astype(str) + '_' + df['nodeID'].astype(str)

                        # make DKI A dataframe
                        df_dki_fa = df[['tract_node', 'dki_fa']]
                        df_dki_fa = df_dki_fa.set_index('tract_node')
                        df_dki_fa = df_dki_fa.to_numpy()
                        savetxt(f'{sub_dir}/{sub}_tract_profiles_dki_fa_nocovbat.csv', df_dki_fa, delimiter=',')
                    
                        # make DKI MD dataframe
                        df_dki_md = df[['tract_node', 'dti_md']]
                        df_dki_md = df_dki_md.set_index('tract_node')
                        df_dki_md = df_dki_md.to_numpy()
                        savetxt(f'{sub_dir}/{sub}_tract_profiles_dki_md_nocovbat.csv', df_dki_md, delimiter=',')

                        # make DTI FA dataframe
                        df_fa = df[['tract_node', 'dti_fa']]
                        df_fa = df_fa.set_index('tract_node')
                        df_fa = df_fa.to_numpy()
                        savetxt(f'{sub_dir}/{sub}_tract_profiles_dti_fa_nocovbat.csv', df_fa, delimiter=',')
                    
                        # make DTI MD dataframe
                        df_md = df[['tract_node', 'dti_md']]
                        df_md = df_md.set_index('tract_node')
                        df_md = df_md.to_numpy()
                        savetxt(f'{sub_dir}/{sub}_tract_profiles_dti_md_nocovbat.csv', df_md, delimiter=',')

    # Save the subject IDs with missing tracts to a text file
    with open(os.path.join(sample_selection_dir, 'subs_missing_tracts.txt'), 'w') as f:
        for subject in subs_missing_tracts:
            f.write(f"{subject}\n")

    # Save the rows with NaNs for each subject to a text file
    with open(os.path.join(sample_selection_dir, 'subs_with_nans.txt'), 'w') as f:
        for subject, nan_rows in subs_with_nans.items():
            f.write(f"{subject}: {nan_rows}\n")

 