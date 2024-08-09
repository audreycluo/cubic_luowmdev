#!/bin/bash

# this can just be run by doing ./wrapper_tract_to_cortex_jobs.sh

datasets=("HCPD")

for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_${dataset}.json"
    
    # set where to save out error/output logs
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/tract_to_cortex/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi

    # set dirs
    data_root=$(jq -r '.data_root' ${config_file})
    qsiprep_dir="${data_root}/datalad_qsiprep"
    pyafq_dir="${data_root}/babs_qsirecon_pyafq/merge_ds"

    # set subjects file
    #subjects_file="${data_root}/subject_list/${dataset}_tempsubject_list.txt"
    subjects_file="${data_root}/subject_list/${dataset}_subject_list_test.txt"
    
    # loop through subjects in the file
    while IFS= read -r subject; do
        echo "Processing ${subject}"

        # submit first job
        jobid1=$(sbatch datalad_get_trks.sh ${subject} ${config_file} ${logs_dir} ${pyafq_dir} ${qsiprep_dir} | awk '{print $4}')
        echo "Submitted datalad_get_trks.sh (Job ID: $jobid1)"

        # submit my second job with dependency on the first job - i.e. it'll be submitted only after first job is done!
        jobid2=$(sbatch --dependency=afterok:$jobid1 trk_to_tck.sh ${subject} ${config_file} ${logs_dir} ${pyafq_dir} | awk '{print $4}')
        echo "Submitted trk_to_tck.sh with dependency on datalad_get_trks.sh (Job ID: $jobid2)"

        # submit my third job with dependency on the second job 
        jobid3=$(sbatch --dependency=afterok:$jobid2 delete_trks.sh ${subject} ${config_file} ${logs_dir} ${pyafq_dir} | awk '{print $4}')
        echo "Submitted delete_trks.sh with dependency on trk_to_tck.sh (Job ID: $jobid3)"

        # submit my fourth job with dependency on the third job 
        jobid4=$(sbatch --dependency=afterok:$jobid3 pyafq-bundles_to_tdi.sh ${subject} ${config_file} ${logs_dir} ${pyafq_dir} | awk '{print $4}')
        echo "Submitted pyafq-bundles_to_tdi.sh with dependency on delete_trks.sh (Job ID: $jobid4)"

        # submit my last job with dependency on the fourth job 
        jobid5=$(sbatch --dependency=afterok:$jobid4 cleanup.sh ${subject} ${config_file} ${logs_dir} ${pyafq_dir} | awk '{print $4}')
        echo "Submitted cleanup.sh with dependency on pyafq-bundles_to_tdi.sh (Job ID: $jobid5)"

    done < ${subjects_file}

done
# yayyy
       