#!/bin/bash

# submit this with bash b00_wrapper_transforms.sh

datasets=("HCPD")  

for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_${dataset}.json"
    
    # where to save output and error logs
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/tract_to_cortex/b_transforms/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi

    # set dir
    data_root=$(jq -r '.data_root' ${config_file})
    
    # subjects file
    #subjects_file="${data_root}/subject_list/${dataset}_tempsubject_list.txt"
    #subjects_file="${data_root}/subject_list/${dataset}_subject_list_test.txt"
    subjects_file="${data_root}/subject_list/${dataset}_subject_list_rerun.txt"
    
    # Loop through subjects in the file
    while IFS= read -r subject; do
        echo "Processing ${subject}"

        # define variables for the SLURM directive
        job_name="determine_freesurfer-to-native_acpc_volume_xfm"
        output_file="${logs_dir}/${job_name}_${subject}.out"
        error_file="${logs_dir}/${job_name}_${subject}.err"
        
        # create a unique temporary file for the first job
        temp_file1=$(mktemp temp_sbatch_script_${subject}_XXXX.sh)
        
        # create and submit the first job
        cat <<EOT > ${temp_file1}
#!/bin/bash
#SBATCH --job-name=${job_name}_${subject}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=${output_file}
#SBATCH --error=${error_file}

./b01_determine_freesurfer-to-native_acpc_volume_xfm.sh ${subject} ${config_file} 
EOT
        jobid1=$(sbatch ${temp_file1} | awk '{print $4}')
        echo "Submitted determine_freesurfer-to-native_acpc_volume_xfm.sh (Job ID: $jobid1)"

        # create a unique temporary file for the second job
        temp_file2=$(mktemp temp_sbatch_script_${subject}_XXXX.sh)
        
        # submit the second job
        job_name="align_surfaces_with_tract_volumes"
        output_file="${logs_dir}/${job_name}_${subject}.out"
        error_file="${logs_dir}/${job_name}_${subject}.err"

        cat <<EOT > ${temp_file2}
#!/bin/bash
#SBATCH --job-name=${job_name}_${subject}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=${output_file}
#SBATCH --error=${error_file}
#SBATCH --dependency=afterok:$jobid1

./b02_align_surfaces_with_tract_volumes.sh ${subject} ${config_file} 
EOT
        jobid2=$(sbatch ${temp_file2} | awk '{print $4}')
        echo "Submitted align_surfaces_with_tract_volumes.sh with dependency on determine_freesurfer-to-native_acpc_volume_xfm.sh (Job ID: $jobid2)"
        
        # remove the temporary files after submission
        rm -f ${temp_file1} ${temp_file2}
        
    done < ${subjects_file}
done
