#!/bin/bash

dataset=$1
tract_name=$2

 

# Run the R script for the specific dataset and tract
singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.7.sif Rscript --save /cbica/projects/luo_wm_dev/code/tract_profiles/significance_testing/NEST/deep_to_peripheral/NEST_wrapper_clipEnds_clip5.R ${dataset} ${tract_name}
