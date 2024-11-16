#!/bin/bash

# set variables

r_script="/cbica/projects/luo_wm_dev/code/tract_profiles/covbat_harmonization/covbat_all_datasets.R"

# where to save output and error logs
logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/covbat/all_datasets"
if [ ! -d "${logs_dir}" ]; then
    mkdir -p "${logs_dir}"
fi
job_name="covbat_all_datasets"
sbatch --job-name="${job_name}" --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=10G --time=01:00:00 --output="${logs_dir}/${job_name}_%j.out" --error="${logs_dir}/${job_name}_%j.err" --wrap="singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.5.sif Rscript --save ${r_script}"

echo "Submitted covbat across all datasets"
 
