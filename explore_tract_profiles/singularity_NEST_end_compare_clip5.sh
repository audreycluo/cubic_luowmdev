#!/bin/bash

dataset=$1
tract_name=$2
 
# Run the R script for the specific dataset and tract
singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.6.sif Rscript --save /cbica/projects/luo_wm_dev/code/tract_profiles/explore_tract_profiles/NEST_wrapper_end_compare_clip5.R ${dataset} ${tract_name}
