#!/bin/bash

singularity run --cleanenv \
  /cbica/projects/luo_wm_dev/software/ModelArray/modelarray_confixel_latest.sif \
   Rscript --save /cbica/projects/luo_wm_dev/code/tract_profiles/fit_GAMs/s04_ModelArray_fitGAMs_tract_profiles.R $1 
