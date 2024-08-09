import csv
import os
from os.path import join as ospj
from numpy import savetxt
import pandas as pd
import re
import sys
import json


# This script reformats the pyAFQ tract profiles output for each subject to have separate csv's for dti_fa and dti_md in preparation for ModelArray.
# Tract profiles are modified to have 1 entry per unique tractID/nodeID combination (i.e. ARC_L_1) 
# Each csv is a 1D array of length 2200 (100 nodes * 22 tracts)

# fyi: sub-0197045 is missing tract "FP"  
 
########################################
# Set directories
########################################
config_file = "/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_HCPD.json"

# Read config from the specified file
with open(config_file, "rb") as f:
    config = json.load(f)

dataset = config['dataset']
data_root = config['tract_profiles_data_root']
demographics = config['demographics_file']
qc = config['qc_file']

 
########################################
# Reformat tract profiles
########################################

# load each subject's csv. make a new dataframe for fa and md separately
# that dataframe should have column names as tractID_nodeID and the number should correspond to fa or md
 

# Iterate through each subject's directory
if dataset == "HCPD" or dataset == "HBN":
    for sub in os.listdir(data_root):
        sub_dir = os.path.join(data_root, sub)
        if os.path.isdir(sub_dir):
            # Iterate through each subject's tract profile
            for filename in os.listdir(sub_dir):
                if filename.endswith("profiles_dwi.csv"):
                    filepath = os.path.join(sub_dir, filename)
                    # Read the csv file
                    df = pd.read_csv(filepath)
                    if len(df) < 2200:
                        print(sub)
                    else:
                        # Create a new column tract_node by combining tractID and nodeID
                        df['tract_node'] = df['tractID'].astype(str) + '_' + df['nodeID'].astype(str)
                        # make FA dataframe
                        df_fa = df[['tract_node', 'dti_fa']]
                        df_fa = df_fa.set_index('tract_node')
                        df_fa = df_fa.to_numpy()
                        savetxt(f'{sub_dir}/{sub}_tract_profiles_dti_fa_nocovbat.csv', df_fa, delimiter=',')
                    
                        # make MD dataframe
                        df_md = df[['tract_node', 'dti_md']]
                        df_md = df_md.set_index('tract_node')
                        df_md = df_md.to_numpy()
                        savetxt(f'{sub_dir}/{sub}_tract_profiles_dti_md_nocovbat.csv', df_md, delimiter=',')
else:
    for sub in os.listdir(data_root):
        sub_dir = os.path.join(data_root, sub)
        if os.path.isdir(sub_dir):
            # Iterate through each subject's tract profile
            for filename in os.listdir(sub_dir):
                if filename.endswith("profiles_dwi.csv"):
                    filepath = os.path.join(sub_dir, filename)
                    # Read the csv file
                    df = pd.read_csv(filepath)
                    if len(df) < 2200:
                        print(sub)
                    else:
                        # Create a new column tract_node by combining tractID and nodeID
                        df['tract_node'] = df['tractID'].astype(str) + '_' + df['nodeID'].astype(str)
                        # make FA dataframe
                        df_fa = df[['tract_node', 'dti_fa']]
                        df_fa = df_fa.set_index('tract_node')
                        df_fa = df_fa.to_numpy()
                        savetxt(f'{sub_dir}/{sub}_tract_profiles_dti_fa.csv', df_fa, delimiter=',')
                    
                        # make MD dataframe
                        df_md = df[['tract_node', 'dti_md']]
                        df_md = df_md.set_index('tract_node')
                        df_md = df_md.to_numpy()
                        savetxt(f'{sub_dir}/{sub}_tract_profiles_dti_md.csv', df_md, delimiter=',')
                    

## no longer want this code - removing subjects if they're missing tracts (as of 6/2024)
if dataset == "HCPD":                  
    # make files separately for sub_0197045: add NaNs for FP
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
    
    # Create a new column tract_node by combining tractID and nodeID
    sub_0197045['tract_node'] = sub_0197045['tractID'].astype(str) + '_' + sub_0197045['nodeID'].astype(str)
    # make FA dataframe
    sub_0197045_fa = sub_0197045[['tract_node', 'dti_fa']]
    sub_0197045_fa = sub_0197045_fa.set_index('tract_node')
    sub_0197045_fa = sub_0197045_fa.to_numpy()
    savetxt(f'{sub_dir}/{sub}_tract_profiles_dti_fa.csv', sub_0197045_fa, delimiter=',')
    
    # make MD dataframe
    sub_0197045_md = sub_0197045[['tract_node', 'dti_md']]
    sub_0197045_md = sub_0197045_md.set_index('tract_node')
    sub_0197045_md = sub_0197045_md.to_numpy()
    savetxt(f'{sub_dir}/{sub}_tract_profiles_dti_md.csv', sub_0197045_md, delimiter=',')
    