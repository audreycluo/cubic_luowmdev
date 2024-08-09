#!/bin/bash

config_file="../config/config_HCPD.json"
dataset=$(jq -r '.dataset' ${config_file})

########################################
# Set log directory
########################################

if [ ! -d /cbica/projects/luo_wm_dev/code/tract_profiles/logs ]; then
	mkdir -p /cbica/projects/luo_wm_dev/code/tract_profiles/logs
fi
outputs_dir_logs="/cbica/projects/luo_wm_dev/code/tract_profiles/logs"

########################################
# Submit job
########################################
qsub -l h_vmem=60G,s_vmem=60G \
    -N datalad_get_zips_${dataset} \
    -b y \
    -V \
	-j n \
	-o ${outputs_dir_logs} \
	-e ${outputs_dir_logs} \
	./datalad_qsiprep_zips.sh -c ${config_file}
	