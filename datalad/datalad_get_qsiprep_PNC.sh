#!/bin/bash

# qsub -l h_vmem=64G,s_vmem=64G datalad_get_qsiprep_PNC.sh  
cd /cbica/projects/luo_wm_dev/input/PNC/datalad_qsiprep

missing=$(git annex find --not --in here)
datalad get $missing 

 