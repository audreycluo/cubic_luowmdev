#!/bin/bash

singularity run --cleanenv \
  /cbica/projects/luo_wm_dev/software/docker/r-packages-for-cubic_0.0.4.sif \
   Rscript --save /cbica/projects/luo_wm_dev/code/tract_profiles/fit_GAMs/s05a_compute_posterior_percentchange.R $1 $2