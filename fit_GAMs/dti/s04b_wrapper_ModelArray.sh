#!/bin/bash


datasets=("HCPD" "HBN")
functions_dir="/cbica/projects/luo_wm_dev/code/tract_profiles"

for dataset in "${datasets[@]}"; do
	config=${functions_dir}/config/config_${dataset}.json
	data_root=$(jq -r '.tract_profiles_data_root' ${config})
	outputs_root=$(jq -r '.tract_profiles_outputs_root' ${config})
	
	config_csv=${functions_dir}/config/config_${dataset}.csv
	########################################
	# Check for required directories
	########################################

	# Check for GAM directory
	if [ ! -d ${outputs_root}/GAM/ ]; then
		echo "Skipping - no GAM directory"
		continue
	fi

	########################################
	# Create log directory if not already made
	########################################
	outputs_dir_logs="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/fit_GAMs/ModelArray/${dataset}"

	if [ ! -d ${outputs_dir_logs} ]; then
		mkdir -p ${outputs_dir_logs} 
	fi

	########################################
	# Submit job
	########################################

	qsub -l h_vmem=64G,s_vmem=64G \
		-N qsub_ModelArray_tract_profiles_${dataset} \
		-b y \
		-V \
		-j n \
		-o ${outputs_dir_logs} \
		-e ${outputs_dir_logs} \
		./s04a_singularity_ModelArray.sh ${config_csv}

done 
 