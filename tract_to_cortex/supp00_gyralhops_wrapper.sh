#!/bin/bash

datasets=( "PNC" "HCPD" "HBN") 

for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_${dataset}.json"
    
    # where to save output and error logs
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/tract_to_cortex/supp_gyralhops/${dataset}"
    if [ ! -d "${logs_dir}" ]; then
        mkdir -p ${logs_dir}
    fi

    # set dir
    data_root=$(jq -r '.data_root' ${config_file})
    qsiprep_dir="${data_root}/raw/datalad_qsiprep"
    # make probseg dir
    if [ ! -d "${data_root}/derivatives/gyral_hops" ]; then
        mkdir -p ${data_root}/derivatives/gyral_hops
    fi
    
    # subjects file
    subjects_file="${data_root}/subject_list/final_sample/${dataset}_WMDev_FinalSample.txt"
    mapfile -t subjects_array < <(tail -n +2 ${subjects_file}) # skip header
    for i in "${!subjects_array[@]}"; do
        subjects_array[$i]=$(echo "${subjects_array[$i]}" | tr -d '"') # remove quotes
    done
    subject_count=${#subjects_array[@]} 
    
    jid_supp01=$(sbatch --parsable \
        --array=0-$((subject_count-1))%2 \
        --job-name="supp01_gyralhops_datalad_${dataset}" \
        ./supp01_gyralhops_datalad.sh ${subjects_file} ${qsiprep_dir} ${logs_dir} "supp01_gyralhops_datalad_${dataset}")

   jid_supp02=$(sbatch --parsable \
        --dependency=afterok:${jid_supp01} \
        --array=0-$((subject_count-1)) \
        --job-name="supp02_gyralhops_probseg_to_acpc_${dataset}" \
        ./supp02_gyralhops_probseg_to_acpc.sh ${subjects_file} ${config_file} ${logs_dir} "supp02_gyralhops_probseg_to_acpc_${dataset}")

    jid_supp03=$(sbatch --parsable \
        --dependency=afterok:${jid_supp02} \
        --array=0-$((subject_count-1)) \
        --job-name="supp03_gyralhops_vol_to_surf_${dataset}" \
        ./supp03_gyralhops_vol_to_surf.sh ${subjects_file} ${config_file} ${logs_dir} "supp03_gyralhops_vol_to_surf_${dataset}")
    
    sbatch --parsable \
            --dependency=afterok:${jid_supp03} \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task=1 \
            --time=3:00:00 \
            --output=${logs_dir}/supp04_gyralhops_best_depths_${dataset}_%j.out \
            --error=${logs_dir}/supp04_gyralhops_best_depths_${dataset}_%j.err \
            --job-name="supp04_gyralhops_best_depths_${dataset}" \
            --wrap="python ./supp04_gyralhops_best_depths.py ${config_file}"

done
