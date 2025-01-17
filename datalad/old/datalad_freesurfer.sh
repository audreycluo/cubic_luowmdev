#!/bin/bash

# submit this script with something like: 
# qsub -l h_vmem=60G,s_vmem=60G -o ~/jobs -e ~/jobs datalad_freesurfer.sh -c ../config_<dataset>.json 

##############################################
# Set config file (based on dataset to datalad get)
##############################################
# Parse command-line arguments
while getopts ":c:" opt; do
  case $opt in
    c) config_file="$OPTARG"
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
    ;;
  esac
done

# Check if the config_file variable is set
if [ -z "$config_file" ]; then
  echo "Error: Please specify a config file using the -c option." >&2
  exit 1
fi
########################################
# Set directories
########################################
dataset=$(jq -r '.dataset' "$config_file")
data_root=$(jq -r '.data_root' "$config_file")
freesurfer_datalad_dir="${data_root}/datalad_freesurfer"
freesurfer_dir="${data_root}/${dataset}_freesurfer"
if [ ! -d ${freesurfer_dir} ]; then
    mkdir -p ${freesurfer_dir} 
fi

########################################
# Set ssh key 
########################################
# to avoid putting in password for each datalad get
eval $(ssh-agent)
ssh-add ~/.ssh/id_rsa
  
########################################
# freesurfer
########################################
# clone fstabulate 
datalad clone ria+ssh://audluo@bblsub.pmacs.upenn.edu:/static/LINC_${dataset}#~fstabulate datalad_freesurfer

cd ${freesurfer_datalad_dir}/inputs/data
datalad get -n . # to get input zips 
 