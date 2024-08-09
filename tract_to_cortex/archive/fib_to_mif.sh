#!/bin/bash

export FOD_TEMPLATE_FIB="/cbica/projects/luo_wm_dev/input/HCP1065/fib_MNI/HCP1065.1.25mm.odf.fib"
export MNI_TEMPLATE_T1="/cbica/projects/luo_wm_dev/input/HCP1065/T1_MNI/mni_icbm152_t1_tal_nlin_asym_09a.nii"
export OUT_TEMPLATE_MIF="/cbica/projects/luo_wm_dev/input/HCP1065/mif_MNI/HCP1065.1.25mm.fod.mif"
export QSIPREP_SIF="/cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-0.21.4.sif"
 
 
singularity exec --writable-tmpfs \
    -B $FOD_TEMPLATE_FIB \
    -B $MNI_TEMPLATE_T1 \
    $QSIPREP_SIF \
    fib2mif --fib $FOD_TEMPLATE_FIB --mif $OUT_TEMPLATE_MIF --ref_image $MNI_TEMPLATE_T1

    