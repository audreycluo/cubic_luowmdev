#!/bin/bash

# set variables
datasets=("HBN" "PNC")
scalars=("dti_md")

r_script="/cbica/projects/luo_wm_dev/code/tract_profiles/fit_GAMs/fit_GAMs_envSES.R"

# loop through each multishell dataset and scalar
for dataset in "${datasets[@]}"; do
  config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_${dataset}.json"
    
    # where to save output and error logs
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/fit_GAMs/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi
  
  for scalar in "${scalars[@]}"; do
      job_name="fit_GAMs_envSES_${dataset}_${scalar}"
      # submit the job to SLURM using my Singularity container
      sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G --time=6:00:00 --output=${logs_dir}/${job_name}_%j.out --error=${logs_dir}/${job_name}_%j.err --wrap="singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.5.sif Rscript --save ${r_script} ${dataset} ${scalar}"
      echo "Submitted ${dataset} ${scalar}"
  done
done
 