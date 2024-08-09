#!/bin/bash


# Adapted from Marc Jaskir's excellent code: https://github.com/marcjaskir/tracts/blob/main/04-map_tract_volumes_to_surfaces.sh

module load connectome_workbench/1.4.2

########################################
# Define variables passed from the submission script
########################################
subject=$1
config_file=$2

########################################
# Set directories
########################################
proj_root="/cbica/projects/luo_wm_dev/"
dataset=$(jq -r '.dataset' "$config_file")
data_root=$(jq -r '.data_root' "$config_file")
derivs_root="${data_root}/derivatives"
xfms_root="${derivs_root}/${dataset}_fs_qsiprep_xfm"
cortical_maps_root="${derivs_root}/${dataset}_vol_to_surf"
 

# Create output directory for fslr_32 outputs
if [ ! -d ${cortical_maps_root}/${subject}/fslr_32k ]; then
    mkdir -p ${cortical_maps_root}/${subject}/fslr_32k
fi
outputs_dir_native_acpc=${cortical_maps_root}/${subject}/native_acpc # already exists from c01
outputs_dir_fslr_32k=${cortical_maps_root}/${subject}/fslr_32k
 

########################################
# Create midthickness surfaces from Freesurfer data (used for resampling to fsLR)
########################################
# fsLR 32k

# Left hemisphere
wb_shortcuts -freesurfer-resample-prep \
    ${xfms_root}/${subject}/surfaces/native_acpc/${subject}.lh.white.native_acpc.surf.gii \
    ${xfms_root}/${subject}/surfaces/native_acpc/${subject}.lh.pial.native_acpc.surf.gii \
    ${xfms_root}/${subject}/surfaces/freesurfer/${subject}.lh.sphere.freesurfer.surf.gii \
    ${proj_root}/templates/fslr_32k/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii \
    ${outputs_dir_fslr_32k}/lh.midthickness.surf.gii \
    ${outputs_dir_fslr_32k}/lh.midthickness.fslr_32k.surf.gii \
    ${outputs_dir_fslr_32k}/lh.sphere.reg.surf.gii

# Right hemisphere
wb_shortcuts -freesurfer-resample-prep \
    ${xfms_root}/${subject}/surfaces/native_acpc/${subject}.rh.white.native_acpc.surf.gii \
    ${xfms_root}/${subject}/surfaces/native_acpc/${subject}.rh.pial.native_acpc.surf.gii \
    ${xfms_root}/${subject}/surfaces/freesurfer/${subject}.rh.sphere.freesurfer.surf.gii \
    ${proj_root}/templates/fslr_32k/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii \
    ${outputs_dir_fslr_32k}/rh.midthickness.surf.gii \
    ${outputs_dir_fslr_32k}/rh.midthickness.fslr_32k.surf.gii \
    ${outputs_dir_fslr_32k}/rh.sphere.reg.surf.gii


########################################
# Warp tract surfaces to fsLR
########################################

for tract_file in ${outputs_dir_native_acpc}/*.shape.gii; do

    # Extract tract label
    tract_fname=$(basename ${tract_file})
    tract_label=$(echo ${tract_fname} | cut -d '_' -f2)
    depth=$(echo "${tract_fname}" | sed -n 's/.*_\([0-9]*\.[0-9]*\)\.shape\.gii/\1/p')

    echo "Warping ${tract_label} ${depth} to fsLR"

    # Extract hemisphere
    hemi=$(echo ${tract_label} | rev | cut -c1-1 | rev)
   
    ###############
    # fsLR 32k
    ###############
        # Left hemisphere
    if [ ${hemi} == 'L' ]; then
        wb_command -metric-resample \
            ${tract_file} \
            ${xfms_root}/${subject}/surfaces/freesurfer/${subject}.lh.sphere.freesurfer.surf.gii \
            ${proj_root}/templates/fslr_32k/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii \
            ADAP_BARY_AREA \
            ${outputs_dir_fslr_32k}/${subject}.${tract_label}_${depth}.shape.gii \
            -area-surfs ${outputs_dir_fslr_32k}/lh.midthickness.surf.gii ${outputs_dir_fslr_32k}/lh.midthickness.fslr_32k.surf.gii

    # Right hemisphere
    elif [ ${hemi} == 'R' ]; then
        wb_command -metric-resample \
            ${tract_file} \
            ${xfms_root}/${subject}/surfaces/freesurfer/${subject}.rh.sphere.freesurfer.surf.gii \
            ${proj_root}/templates/fslr_32k/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii \
            ADAP_BARY_AREA \
            ${outputs_dir_fslr_32k}/${subject}.${tract_label}_${depth}.shape.gii \
            -area-surfs ${outputs_dir_fslr_32k}/rh.midthickness.surf.gii ${outputs_dir_fslr_32k}/rh.midthickness.fslr_32k.surf.gii
    fi       
done

# Clean up intermediary files
rm ${outputs_dir_fslr_32k}/lh.midthickness.surf.gii
rm ${outputs_dir_fslr_32k}/rh.midthickness.surf.gii
rm ${outputs_dir_fslr_32k}/lh.sphere.reg.surf.gii
rm ${outputs_dir_fslr_32k}/rh.sphere.reg.surf.gii
rm ${outputs_dir_fslr_32k}/lh.midthickness.fslr_32k.surf.gii
rm ${outputs_dir_fslr_32k}/rh.midthickness.fslr_32k.surf.gii
 