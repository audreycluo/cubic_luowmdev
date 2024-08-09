#!/bin/bash

datasets=("HCPD")

for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_${dataset}.json"
    
    # Set where to save out error/output logs
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/tract_to_cortex/c_vol_to_surf/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi

    # Set dirs
    data_root=$(jq -r '.data_root' ${config_file})
     
    # Set subjects file
    subjects_file="${data_root}/subject_list/${dataset}_tempsubject_list.txt"
    #subjects_file="${data_root}/subject_list/${dataset}_subject_list_test.txt"
    #subjects_file="${data_root}/subject_list/${dataset}_subject_list_rerun.txt"
    # Loop through subjects in the file
    while IFS= read -r subject; do
        echo "Processing ${subject}"

        job_name="vol_to_surf"
        output_file="${logs_dir}/${job_name}_${subject}.out"
        error_file="${logs_dir}/${job_name}_${subject}.err"

        # Submit the first job
        jobid1=$(sbatch <<EOT | awk '{print $4}'
#!/bin/bash
#SBATCH --job-name=${job_name}_${subject}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=${output_file}
#SBATCH --error=${error_file}

source /cbica/projects/luo_wm_dev/miniconda3/etc/profile.d/conda.sh
conda activate luo_wm_dev

python c01_vol_to_surf_individual.py ${subject} ${config_file} 2.50
EOT
)
        echo "Submitted vol_to_surf (Job ID: $jobid1)"
        
        # submit second job
        job_name="native_to_fsLR"
        output_file="${logs_dir}/${job_name}_${subject}.out"
        error_file="${logs_dir}/${job_name}_${subject}.err"

        jobid2=$(sbatch <<EOT | awk '{print $4}'
#!/bin/bash
#SBATCH --job-name=${job_name}_${subject}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=${output_file}
#SBATCH --error=${error_file}
#SBATCH --dependency=afterok:${jobid1}

./c02_native_to_fsLR.sh ${subject} ${config_file} 
EOT
)
        echo "Submitted native_to_fsLR.sh with dependency on vol_to_surf.sh (Job ID: $jobid2)"

    done < ${subjects_file}
done
