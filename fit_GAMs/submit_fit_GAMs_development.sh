#!/bin/bash

# set variables
datasets=("HCPD" "HBN")
scalars=("dki_fa" "dki_md" "dti_fa" "dti_md")
scalars_PNC=("dti_fa" "dti_md")
r_script="/cbica/projects/luo_wm_dev/code/tract_profiles/fit_GAMs/fit_GAMs_development.R"

# loop through each multishell dataset and scalar
for dataset in "${datasets[@]}"; do
  config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_${dataset}.json"
    
    # where to save output and error logs
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/fit_GAMs/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi
  
  for scalar in "${scalars[@]}"; do
      job_name="fit_GAMs_dev_${dataset}_${scalar}"
      # submit the job to SLURM using my Singularity container
      sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=25G --output=${logs_dir}/${job_name}.out --error=${logs_dir}/${job_name}.err --wrap="singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.5.sif Rscript --save ${r_script} ${dataset} ${scalar}"
      echo "Submitted ${dataset} ${scalar}"
  done
done

# submit singleshell data
dataset="PNC"
config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_${dataset}.json"
  
# where to save output and error logs
logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/fit_GAMs/${dataset}"
if [ ! -d "${logs_dir}" ]; then
  mkdir -p ${logs_dir}
fi

for scalar in "${scalars_PNC[@]}"; do
  job_name="fit_GAMs_dev_${dataset}_${scalar}"
  # submit the job to SLURM using my Singularity container
  sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=25G --output=${logs_dir}/${job_name}.out --error=${logs_dir}/${job_name}.err --wrap="singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.5.sif Rscript --save ${r_script} ${dataset} ${scalar}"
  echo "Submitted ${dataset} ${scalar}"
done
