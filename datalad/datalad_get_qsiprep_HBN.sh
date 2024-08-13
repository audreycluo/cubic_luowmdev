#!/bin/bash
#SBATCH --job-name=datalad_get_HBN
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/HBN/qsiprep.out
#SBATCH --error=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/HBN/qsiprep.err

cd /cbica/projects/luo_wm_dev/input/HBN/raw/datalad_qsiprep
missing=$(git annex find --not --in here)
datalad get $missing 

 