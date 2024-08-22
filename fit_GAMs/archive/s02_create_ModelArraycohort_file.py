import csv
import os
from os.path import join as ospj
import pandas as pd
import re
import sys
import json
 

# This script creates cohort files for ModelArray with tract profiles data

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
demographics = config['demographics_file']
qc = config['qc_file']

# make output cohortfiles directories
if not os.path.exists(ospj(f"{outputs_root}", "cohortfiles", "dti_fa")):
        os.makedirs(ospj(f"{outputs_root}", "cohortfiles", "dti_fa"))
        print(f"Directory 'cohortfiles/dti_fa' created.")
else:
        print(f"Directory 'cohortfiles/dti_fa' already exists.")

if not os.path.exists(ospj(f"{outputs_root}", "cohortfiles", "dti_md")):
        os.makedirs(ospj(f"{outputs_root}", "cohortfiles", "dti_md"))
        print(f"Directory 'cohortfiles/dti_md' created.")
else:
        print(f"Directory 'cohortfiles/dti_md' already exists.")
 

########################################
# Make Cohort Files
########################################

# load demographics and qc files
dem_df = pd.read_csv(demographics)
qc_df = pd.read_csv(qc)
merged_df = pd.merge(dem_df, qc_df, on="sub")
merged_df = merged_df.drop(columns=['t1_neighbor_corr', 'race', 'site']) 

 
# define function for making cohort files  
 
def make_cohort_df(scalars):
    
    # Iterate over scalars
    for scalar in scalars:
        scalar_dataframes = {}
        print(scalar)
        
        # Add a new column 'scalar_name' with scalar name
        scalar_df = merged_df.copy()  # Create a copy of the original dataframe
        scalar_df['scalar_name'] = scalar
        
        # Create source_file column    
        scalar_df['source_file'] = scalar_df.apply(lambda row: f"{data_root}/{row['sub']}/{row['sub']}_tract_profiles_{scalar}.csv", axis=1)

        scalar_df.to_csv(f'{outputs_root}/cohortfiles/{scalar}/{scalar}_cohortfile.csv', index=False)
     

# make cohort files
scalars = ["dti_fa", "dti_md"]
make_cohort_df(scalars)
