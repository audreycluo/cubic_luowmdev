import nibabel as nib
import numpy as np


# HCP1650 t1
t1_img = nib.load('/cbica/projects/luo_wm_dev/input/HCP1065/T1_MNI/mni_icbm152_t1_tal_nlin_asym_09a.nii')

# load ROI from the default bundle dict
roi_img = nib.load('/Users/audluo/AFQ_data/templates/ATR_roi1_L.nii.gz')

# Print header information
print("T1 Image Header:")
print(t1_img.header)
print(t1_img.affine)

print("\nROI Image Header:")
print(roi_img.header)
print(roi_img.affine)

# nooo they're not aligned lol