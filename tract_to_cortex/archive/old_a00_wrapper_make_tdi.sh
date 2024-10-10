#!/bin/bash

# submit this with bash a00_wrapper_make_tdi.sh

datasets=("HCPD" "HBN" "PNC")  

for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_${dataset}.json"
    
    # where to save output and error logs
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/tract_to_cortex/a_make_tdi/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi

    # set dir
    data_root=$(jq -r '.data_root' ${config_file})
    qsiprep_dir="${data_root}/raw/datalad_qsiprep"
    pyafq_dir="${data_root}/babs_qsirecon_pyafq_act/merge_ds"

    # subjects file
    subjects_file="${data_root}/subject_list/${dataset}_tempsubject_list.txt"
    #subjects_file="${data_root}/subject_list/${dataset}_subject_list_test.txt"
    
    
    # Loop through subjects in the file
    while IFS= read -r subject; do
        echo "Processing ${subject}"

        # i want my error and out files to have custom names (job_name + subjectID) and be saved into my logs folder
        # which is why i need to create these temporary sbatch scripts
        # its just way easier to track my progress and troubleshoot later on
        # being extra now = being happier later muahah

        # define variables for the SLURM directive
        job_name="datalad_get_trks"
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

./a01_datalad_get_trks.sh ${subject} ${config_file} ${pyafq_dir} ${qsiprep_dir}
EOT
        jobid1=$(sbatch ${temp_file1} | awk '{print $4}')
        echo "Submitted datalad_get_trks.sh (Job ID: $jobid1)"

        # create a unique temporary file for the second job
        temp_file2=$(mktemp temp_sbatch_script_${subject}_XXXX.sh)

        # submit the second job
        job_name="trk_to_tck"
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

./a02_trk_to_tck.sh ${subject} ${config_file} ${pyafq_dir} 
EOT
        jobid2=$(sbatch ${temp_file2} | awk '{print $4}')
        echo "Submitted trk_to_tck.sh with dependency on datalad_get_trks.sh (Job ID: $jobid2)"

        # create a unique temporary file for the third job
        temp_file3=$(mktemp temp_sbatch_script_${subject}_XXXX.sh)

        # submit the third job
        job_name="delete_trks"
        output_file="${logs_dir}/${job_name}_${subject}.out"
        error_file="${logs_dir}/${job_name}_${subject}.err"

        cat <<EOT > ${temp_file3}
#!/bin/bash
#SBATCH --job-name=${job_name}_${subject}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=${output_file}
#SBATCH --error=${error_file}
#SBATCH --dependency=afterok:$jobid2

./a03_delete_trks.sh ${subject} ${pyafq_dir}
EOT
        jobid3=$(sbatch ${temp_file3} | awk '{print $4}')
        echo "Submitted delete_trks.sh with dependency on trk_to_tck.sh (Job ID: $jobid3)"

        # create a unique temporary file for the fourth job
        temp_file4=$(mktemp temp_sbatch_script_${subject}_XXXX.sh)

        # submit the fourth job
        job_name="pyafq_bundles_to_tdi"
        output_file="${logs_dir}/${job_name}_${subject}.out"
        error_file="${logs_dir}/${job_name}_${subject}.err"

        cat <<EOT > ${temp_file4}
#!/bin/bash
#SBATCH --job-name=${job_name}_${subject}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=${output_file}
#SBATCH --error=${error_file}
#SBATCH --dependency=afterok:$jobid3

./a04_pyafq-bundles_to_tdi.sh ${subject} ${config_file}
EOT
        jobid4=$(sbatch ${temp_file4} | awk '{print $4}')
        echo "Submitted pyafq-bundles_to_tdi.sh with dependency on delete_trks.sh (Job ID: $jobid4)"

        # create a unique temporary file for the fifth job
        temp_file5=$(mktemp temp_sbatch_script_${subject}_XXXX.sh)

        # submit the last job
        job_name="cleanup"
        output_file="${logs_dir}/${job_name}_${subject}.out"
        error_file="${logs_dir}/${job_name}_${subject}.err"

        cat <<EOT > ${temp_file5}
#!/bin/bash
#SBATCH --job-name=${job_name}_${subject}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=${output_file}
#SBATCH --error=${error_file}
#SBATCH --dependency=afterok:$jobid4

./a05_cleanup.sh ${subject} ${config_file}
EOT
        jobid5=$(sbatch ${temp_file5} | awk '{print $4}')
        echo "Submitted cleanup.sh with dependency on pyafq_bundles_to_tdi.sh (Job ID: $jobid5)"

        # remove the temporary files after submission
        rm -f ${temp_file1} ${temp_file2} ${temp_file3} ${temp_file4} ${temp_file5}

    done < ${subjects_file}
done
