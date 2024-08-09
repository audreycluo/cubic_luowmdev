#!/bin/bash


##############################################
# Set config file (based on dataset to analyze)
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

freesurfer_dir="${data_root}/datalad_freesurfer/inputs/data"
qsiprep_dir="${data_root}/datalad_qsiprep"
SUBJECTS_DIR="${freesurfer_dir}/freesurfer"


########################################
# Read in subject ID
########################################
subject=${@:OPTIND}


########################################
# Get required files
########################################
# Freesurfer
cd ${freesurfer_dir}
datalad get ${subject}*.zip
mkdir -p "freesurfer/${subject}/mri"
echo "Freesurfer ${subject} zip successfully gotten"
unzip -j ${subject}*.zip "freesurfer/${subject}/mri/brain.mgz" -d "freesurfer/${subject}/mri" # unzip and copy over just this file
unzip -j ${subject}*.zip "freesurfer/${subject}/mri/ribbon.mgz" -d "freesurfer/${subject}/mri" 
find . -name "${subject}*.zip" -exec datalad drop --nocheck {} \; # drop zip to save space
echo "Freesurfer ${subject} zip successfully dropped"

# QSIPrep
cd ${qsiprep_dir}
datalad get ${subject}*.zip
echo "Qsiprep ${subject} zip successfully gotten"
mkdir -p "qsiprep/${subject}/anat/"
unzip -j ${subject}*.zip "qsiprep/${subject}/anat/${subject}_desc-preproc_T1w.nii.gz" -d "qsiprep/${subject}/anat"
unzip -j ${subject}*.zip "qsiprep/${subject}/anat/${subject}_desc-brain_mask.nii.gz" -d "qsiprep/${subject}/anat"
find . -name "${subject}*.zip" -exec datalad drop --nocheck {} \; # drop zip to save space
echo "Qsiprep ${subject} zip successfully dropped" 


########################################
# Create output directories
########################################   
rm -rf ${data_root}/freesurfer_qsiprep_xfm/${subject} ## to remove
# Create output directory for transforms (including reference volumes) and transformed volumes
if [ ! -d ${data_root}/freesurfer_qsiprep_xfm/${subject} ]; then
    mkdir -p ${data_root}/freesurfer_qsiprep_xfm/${subject}
fi


########################################
# Transformation time!
########################################  
# set variables
run_qsiprep_cmd="singularity exec /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-0.21.4.sif"
workdir="${data_root}/freesurfer_qsiprep_xfm"
 
# Convert FS brain from .mgz to .nii
${run_qsiprep_cmd} mrconvert -strides -1,-2,3 \
    ${SUBJECTS_DIR}/${subject}/mri/brain.mgz ${workdir}/${subject}/fs_brain.nii


# Register FreeSurfer brain to QSIPrep T1w
${run_qsiprep_cmd} antsRegistration --collapse-output-transforms 1 \
    --dimensionality 3 --float 0 \
    --initial-moving-transform [ ${qsiprep_dir}/qsiprep/${subject}/anat/${subject}_desc-preproc_T1w.nii.gz, ${workdir}/${subject}/fs_brain.nii, 1 ] \
    --initialize-transforms-per-stage 0 --interpolation BSpline \
    --output [ ${workdir}/${subject}/transform, ${workdir}/${subject}/transform_Warped.nii.gz ] \
    --transform Rigid[ 0.1 ] \
    --metric Mattes[ ${qsiprep_dir}/qsiprep/${subject}/anat/${subject}_desc-preproc_T1w.nii.gz, ${workdir}/${subject}/fs_brain.nii, 1, 32, Random, 0.25 ] \
    --convergence [ 1000x500x250x100, 1e-06, 10 ] \
    --smoothing-sigmas 3.0x2.0x1.0x0.0mm --shrink-factors 8x4x2x1 \
    --use-histogram-matching 0 \
    --masks [ ${qsiprep_dir}/qsiprep/${subject}/anat/${subject}_desc-brain_mask.nii.gz, NULL ] \
    --winsorize-image-intensities [ 0.002, 0.998 ] \
    --write-composite-transform 0
echo "${subject} registered FreeSurfer brain to QSIPrep T1w"

# Convert ANTs .mat transform to .txt, and rename it
${run_qsiprep_cmd} ConvertTransformFile 3 ${workdir}/${subject}/transform0GenericAffine.mat ${workdir}/${subject}/${subject}_from-FS_to-T1wACPC_mode-image_xfm.txt

# Convert ribbon.mgz to .nii.gz
${run_qsiprep_cmd} mrconvert -strides -1,-2,3 \
    ${SUBJECTS_DIR}/${subject}/mri/ribbon.mgz ${workdir}/${subject}/ribbon.nii.gz

# can delete brain.mgz and subject folder to save space
#rm -rf ${SUBJECTS_DIR}/${subject}
 
${run_qsiprep_cmd} antsApplyTransforms \
    -d 3 \
    -i ${workdir}/${subject}/ribbon.nii.gz \
    -r ${qsiprep_dir}/qsiprep/${subject}/anat/${subject}_desc-preproc_T1w.nii.gz \
    -o ${workdir}/${subject}/ribbon_transformed.nii.gz \
    -n NearestNeighbor \
    -t ${workdir}/${subject}/transform0GenericAffine.mat
echo "${subject} ribbon.mgz transformed to to QSIPrep T1w"


# can delete qsiprep subject folder to save space
#rm -rf ${qsiprep_dir}/qsiprep/${subject}