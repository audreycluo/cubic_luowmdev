#!/bin/bash

# set variables
#datasets=("HCPD" "HBN" "PNC")
datasets=("HBN")
r_script="/cbica/projects/luo_wm_dev/code/tract_profiles/explore_tract_profiles/NEST_wrapper_end_compare_clip5.R"
#tract_list=("Arcuate" "Inferior_Fronto-occipital" "Callosum_Motor")
 
#inputarray=("${tract_list[@]}")
#tract_count=${#inputarray[@]}

# submit job array for each dataset with elements in array being tracts
for dataset in "${datasets[@]}"; do
    logs_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/NEST/${dataset}"
    mkdir -p ${logs_dir}

    sbatch --job-name=NEST_${dataset}_end_compare_low_p \
           --nodes=1 \
           --ntasks=1 \
           --cpus-per-task=8 \
           --mem-per-cpu=4G \
           --time=48:00:00  \
           --output=${logs_dir}/NEST_${dataset}_end_compare_%j_ifo_covbat_all.out \
           --error=${logs_dir}/NEST_${dataset}_end_compare_%j_ifo_covbat_all.err \
           singularity_NEST_end_compare_clip5.sh ${dataset} "Inferior_Fronto-occipital"
done
