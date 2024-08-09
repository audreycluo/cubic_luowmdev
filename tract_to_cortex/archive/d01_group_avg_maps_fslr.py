import json
import nibabel as nib
import numpy as np
import os


import glob

import nibabel as nib
from nilearn import surface, datasets, plotting
from nilearn.plotting import show
import numpy as np
import os
from os.path import join as ospj
import re
import sys

"""
This script takes the binary, subject-level tract-to-cortex maps (in fsLR 32k) for each tract, and averages them across subjects. 
This produces a group-level tract-to-cortex map of the proportion of subjects that have a tract connecting to each cortical region.

ingredients for each group-level map:
- subject-level tract-to-cortex maps for tract of interest

"""

 
config_file = sys.argv[1]

# Read config from the specified file
with open(config_file, "rb") as f:
    config = json.load(f)

data_root = config['data_root']
dataset = config['dataset']
derivs_dir = ospj(data_root, f"derivatives/{dataset}_vol_to_surf")
out_dir = ospj(derivs_dir, subject, "native_acpc")

# Create directory for vol_to_surf outputs
if not os.path.exists(derivs_dir):
    os.makedirs(derivs_dir)

if not os.path.exists(ospj(derivs_dir, subject, "native_acpc")):
        os.makedirs(ospj(derivs_dir, subject, "native_acpc"))
        print(f"Directory derivatives/{dataset}_vol_to_surf/{subject}/native_acpc created.")
else:
        print(f"Directory derivatives/{dataset}_vol_to_surf/{subject}/native_acpc already exists.")


# Define the directories
subjects_dir = "/path/to/subjects"  # Directory containing subject folders
tracts = 22  # Number of tracts


# Initialize a dictionary to store the sum of each tract's data
tract_sums = {tract: None for tract in range(tracts)}

# Loop over each subject
for subject_id in subject_ids:
    subject_path = os.path.join(subjects_dir, subject_id)
    for tract in range(tracts):
        tract_file = os.path.join(subject_path, f"tract_{tract+1}.shape.gii")
        
        # Load the tract-to-cortex map for the current tract and subject
        gii_data = nib.load(tract_file)
        data_array = gii_data.darrays[0].data  # Assuming data is in the first darray
        
        if tract_sums[tract] is None:
            # Initialize the sum array if it's the first subject for this tract
            tract_sums[tract] = np.zeros_like(data_array, dtype=np.float64)
        
        # Accumulate the binary tract-to-cortex data
        tract_sums[tract] += data_array

# Compute the average for each tract to get the proportion of subjects
tract_averages = {tract: tract_sums[tract] / len(subject_ids) for tract in range(tracts)}

# Save the averaged data
output_dir = "/path/to/output"
os.makedirs(output_dir, exist_ok=True)

for tract in range(tracts):
    # Create a new Gifti image with the averaged tract-to-cortex data
    gii_img = nib.gifti.GiftiImage()
    gii_data_array = nib.gifti.GiftiDataArray(data=tract_averages[tract])
    gii_img.add_gifti_data_array(gii_data_array)
    
    # Save the Gifti image to the output directory
    output_file = os.path.join(output_dir, f"average_tract_{tract+1}.shape.gii")
    nib.save(gii_img, output_file)

print("Averaging completed and files saved.")