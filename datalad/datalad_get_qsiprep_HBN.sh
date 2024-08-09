#!/bin/bash

# qsub -l h_vmem=60G,s_vmem=60G datalad_get_qsiprep_HBN.sh  
cd /cbica/projects/luo_wm_dev/input/HBN/datalad_qsiprep
missing=$(git annex find --not --in here)
datalad get $missing 