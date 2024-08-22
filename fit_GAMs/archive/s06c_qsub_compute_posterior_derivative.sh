#!/bin/bash

datasets=("HCPD" "HBN")
scalars=("dti_md" "dti_fa")

functions_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/fit_GAMs"

#################################
# Compute Posterior Derivatives #
#################################
# Loop through each dataset
for dataset in "${datasets[@]}"; do
    for scalar in "${scalars[@]}"; do
        config_file="/cbica/projects/luo_wm_dev/code/superficialWM_analyses/config_${dataset}.json"
        outputs_dir_logs="/cbica/projects/luo_wm_dev/output/${dataset}/tract_profiles/GAM/logs"
        qsub -l h_vmem=128G,s_vmem=128G \
                    -N posterior_derivative_${scalar}_${dataset} \
                    -b y \
                    -V \
                    -j n \
                    -o ${outputs_dir_logs} \
                    -e ${outputs_dir_logs} \
                    ${functions_dir}/s06b_singularity_compute_posterior_derivative.sh ${config_file} ${scalar}
    done
done
