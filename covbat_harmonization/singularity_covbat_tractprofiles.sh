#!/bin/bash

singularity run --cleanenv \
  /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.5.sif \
   Rscript --save /cbica/projects/luo_wm_dev/code/tract_profiles/covbat_harmonization/covbat_tractprofiles.R $1