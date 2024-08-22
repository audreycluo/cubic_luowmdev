#!/bin/bash

# set variables
datasets=("HCPD" "HBN" "PNC")
r_script="/cbica/projects/luo_wm_dev/code/tract_profiles/explore_tract_profiles/NEST_wrapper.R"

# loop through each dataset 
for dataset in "${datasets[@]}"; do
    # where to save output and error logs
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/NEST/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi
  
  
    job_name="NEST_${dataset}"
    # submit the job to SLURM using my Singularity container
    sbatch --job-name=${job_name} --nodes=1 --ntasks=1 --cpus-per-task=8 --mem-per-cpu=3G --output=${logs_dir}/${job_name}.out --error=${logs_dir}/${job_name}.err --wrap="singularity run --cleanenv /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.6.sif Rscript --save ${r_script} ${dataset}"
    echo "Submitted NEST for ${dataset}"
   
done