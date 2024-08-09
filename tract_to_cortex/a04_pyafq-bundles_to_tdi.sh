#!/bin/bash

module load gcc/9.4.0
module load mrtrix/3.0.4 


# Define variables passed from the submission script
subject=$1
config_file=$2
data_root=$(jq -r '.data_root' ${config_file})
dataset=$(jq -r '.dataset' ${config_file})

########################################
# Check for required files
########################################
# Check if converted .tck files exist
if ! find "${data_root}/derivatives/${dataset}_tck_temp/${subject}" -name "*.tck" -print -quit | grep -q .; then
    echo "No tract (.tck) files for ${subject}"
fi

########################################
# Create output directories
########################################
# Create output directories for tract density maps
if [ ! -d ${data_root}/derivatives/${dataset}_tdi_maps/${subject} ]; then
    mkdir -p ${data_root}/derivatives/${dataset}_tdi_maps/${subject}/mgz
    mkdir -p ${data_root}/derivatives/${dataset}_tdi_maps/${subject}/nifti
    mkdir -p ${data_root}/derivatives/${dataset}_tdi_maps/${subject}/tdi_binarized
fi
outputs_dir_mgz=${data_root}/derivatives/${dataset}_tdi_maps/${subject}/mgz
outputs_dir_nifti=${data_root}/derivatives/${dataset}_tdi_maps/${subject}/nifti
outputs_dir_tdi=${data_root}/derivatives/${dataset}_tdi_maps/${subject}/tdi_binarized

########################################
# Convert tck to TDI maps
########################################
# Iterate over tract (.tck) files
for tract in ${data_root}/derivatives/${dataset}_tck_temp/${subject}/*.tck; do
    
    # Extract file name (without .tck extension)
    tract_fname=$(basename ${tract} | sed 's/.tck//g')
    tract_label=$(echo ${tract_fname} | cut -d'-' -f6 | cut -d'_' -f1)

    template=${data_root}/raw/datalad_qsiprep/qsiprep/${subject}/*/dwi/${subject}_*-T1w_dwiref.nii.gz

    # Convert tck files to TDI maps
    if [ ! -f ${outputs_dir_mgz}/${tract_fname}.mgz ]; then
        tckmap ${tract} -template ${template} -contrast tdi ${outputs_dir_mgz}/${subject}_${tract_label}.mgz
    fi

    # Convert TDI maps to nifti
    mri_convert --in_type mgz \
        --out_type nii \
        ${outputs_dir_mgz}/${subject}_${tract_label}.mgz \
        ${outputs_dir_nifti}/${subject}_${tract_label}.LPS.nii.gz

    # Change orientation of NIFTIs to LAS+  
    mri_convert --in_type nii \
                --out_type nii \
                --out_orientation LAS+ \
                ${outputs_dir_nifti}/${subject}_${tract_label}.LPS.nii.gz \
                ${outputs_dir_nifti}/${subject}_${tract_label}.LAS.nii.gz

    mrcalc ${outputs_dir_nifti}/${subject}_${tract_label}.LAS.nii.gz 0 -gt 1 0 -if ${outputs_dir_tdi}/${subject}_${tract_label}_binarized.nii.gz
    # all voxels with nonzero TDI are set to 1. All voxels with TDI = 0 are set to 0. 
done
 
