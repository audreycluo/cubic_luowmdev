#!/bin/bash

# define variables passed from the submission script
subject=$1
config_file=$2
pyafq_dir=$3
qsiprep_dir=$4



######################
# get required files #
######################
cd ${qsiprep_dir}
if find "qsiprep/${subject}" -maxdepth 1 -type d -name "ses-*" | read; then
    echo "qsiprep dir for subject exists"
else
    # get the dwiref image to use as template for DIPY's load_tractogram
    datalad get ${subject}*.zip
    echo "Qsiprep ${subject} zip successfully gotten"
    unzip -o ${subject}*.zip "qsiprep/${subject}/ses*/dwi/${subject}_ses-V1_space-T1w_dwiref.nii.gz" # -o flag to automatically overwrite
    find . -name "${subject}*.zip" -exec datalad drop --nocheck {} \;
    echo "Qsiprep ${subject} zip successfully dropped"
fi

 

cd ${pyafq_dir}
if [ ! -d "qsirecon/${subject}" ]; then
    datalad get ${subject}*.zip
    echo "PyAFQ ${subject} zip successfully gotten"
    unzip -o ${subject}*.zip # -o flag automatically overwrites existing files (like CITATION etc that we don't care about)
    find . -name "${subject}*.zip" -exec datalad drop --nocheck {} \; # drop zip to save space
    echo "PyAFQ ${subject} zip successfully dropped"
fi
 


