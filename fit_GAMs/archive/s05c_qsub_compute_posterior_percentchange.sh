#!/bin/bash

datasets=("HCPD" "HBN")
scalars=("dti_fa" "dti_md")

 
functions_dir="/cbica/projects/luo_wm_dev/code/tract_profiles"

####################################
# Compute Posterior Percent Change #
####################################
# Loop through each dataset
for dataset in "${datasets[@]}"; do
    config=${functions_dir}/config/config_${dataset}.json
	outputs_root=$(jq -r '.tract_profiles_outputs_root' ${config})
    outputs_dir_logs="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/fit_GAMs/posterior_percentchange/${dataset}"

        if [ ! -d ${outputs_dir_logs} ]; then
            mkdir -p ${outputs_dir_logs} 
        fi
        if [ ! -d ${outputs_dir_logs} ]; then
            mkdir -p ${outputs_dir_logs} 
        fi

    for scalar in "${scalars[@]}"; do
         
        qsub -l h_vmem=128G,s_vmem=128G \
                    -N posterior_percentchange_${scalar}_${dataset} \
                    -b y \
                    -V \
                    -j n \
                    -o ${outputs_dir_logs} \
                    -e ${outputs_dir_logs} \
                    ${functions_dir}/fit_GAMs/s05b_singularity_compute_percentchange.sh ${config} ${scalar}
    done
done
