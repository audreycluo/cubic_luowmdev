#!/bin/bash
#SBATCH --job-name=covbat_harmonization
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/covbat/covbat_HCPD_HBN.out"
#SBATCH --error="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/covbat/covbat_HCPD_HBN.err"

datasets=("HCPD" "HBN")

functions_dir="/cbica/projects/luo_wm_dev/code/tract_profiles/covbat_harmonization"
outputs_dir_logs="/cbica/projects/luo_wm_dev/code/tract_profiles/logs"

########################### 
# Do Covbat Harmonization #
########################### 
# Loop through each dataset
for dataset in "${datasets[@]}"; do
    qsub -l h_vmem=64G,s_vmem=64G \
                -N covbat_tractprofiles_${dataset} \
                -b y \
                -V \
                -j n \
                -o ${outputs_dir_logs} \
                -e ${outputs_dir_logs} \
                ${functions_dir}/singularity_covbat_tractprofiles.sh ${dataset} 
done
