#!/bin/bash

# define variables passed from the submission script
subject=$1
config_file=$2
 
source /cbica/projects/luo_wm_dev/miniconda3/etc/profile.d/conda.sh
conda activate luo_wm_dev

# run python script
python b02_align_surfaces_with_tract_volumes.py ${subject} ${config_file}  
