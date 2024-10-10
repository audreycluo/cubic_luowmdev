#!/bin/bash

# Define variables passed from the submission script
subject=$1
config_file=$2
data_root=$(jq -r '.data_root' ${config_file})
dataset=$(jq -r '.dataset' ${config_file})

rm -rf ${data_root}/derivatives/${dataset}_tck_temp/${subject} # delete the subject's tck folder since we don't need it anymore
rm -rf ${data_root}/derivatives/${dataset}_tdi_maps/${subject}/mgz # delete temporary mgz folder
rm -rf ${data_root}/derivatives/${dataset}_tdi_maps/${subject}/nifti # delete temporary nifti folder
 