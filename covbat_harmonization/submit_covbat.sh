#!/bin/bash

# set variables
datasets=("HCPD" "HBN")
r_script="/cbica/projects/luo_wm_dev/code/tract_profiles/covbat_harmonization/covbat_multisite_tract_profiles.R"

# loop through each multishell dataset  
for dataset in "${datasets[@]}"; do
  # where to save output and error logs
  logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/covbat/${dataset}"
  if [ ! -d "${logs_dir}" ]; then
      mkdir -p "${logs_dir}"
  fi

  job_name="covbat_${dataset}"
  # submit the job to SLURM using my Singularity container
  sbatch --job-name="${job_name}" --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=25G --output="${logs_dir}/${job_name}.out" --error="${logs_dir}/${job_name}.err" --wrap="singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.5.sif Rscript --save ${r_script} ${dataset}"
  echo "Submitted ${dataset}"
done
