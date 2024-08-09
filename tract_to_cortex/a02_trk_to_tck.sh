#!/bin/bash

# define variables passed from the submission script
subject=$1
config_file=$2
pyafq_dir=$3
 
source /cbica/projects/luo_wm_dev/miniconda3/etc/profile.d/conda.sh
conda activate luo_wm_dev

# run python script
python a02_trk_to_tck.py ${subject} ${config_file} ${pyafq_dir}
