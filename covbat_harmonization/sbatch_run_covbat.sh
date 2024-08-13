#!/bin/bash

#SBATCH --job-name=covbat_harmonization     # Job name
#SBATCH --ntasks=2                          # Number of tasks (parallel processes)
#SBATCH --cpus-per-task=10                  # Number of CPU cores per task
#SBATCH --mem-per-cpu=10gb                  # Memory per CPU core
#SBATCH --array=1-2                         # Array job, one task per dataset
#SBATCH --output="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/covbat/covbat_%A_%a.out"
#SBATCH --error="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/covbat/covbat_%A_%a.err"

set -e

datasets=("HCPD" "HBN")

# get the dataset corresponding to this task
dataset=${datasets[$SLURM_ARRAY_TASK_ID-1]}
 
 
singularity run --cleanenv \
  /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.5.sif \
   Rscript --save /cbica/projects/luo_wm_dev/code/tract_profiles/covbat_harmonization/covbat_multisite_tract_profiles.R ${dataset}