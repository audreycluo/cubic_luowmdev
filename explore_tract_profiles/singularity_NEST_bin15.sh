#!/bin/bash

dataset=$1
tract_list=$2

inputarray=()
while IFS= read -r line; do
    inputarray+=("$line")
done < "${tract_list}"

# Get the tract for this SLURM array task
tract_name="${inputarray[$SLURM_ARRAY_TASK_ID]}"


# Run the R script for the specific dataset and tract
singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.6.sif Rscript --save /cbica/projects/luo_wm_dev/code/tract_profiles/explore_tract_profiles/NEST_wrapper_clipEnds_bin15.R ${dataset} ${tract_name}
