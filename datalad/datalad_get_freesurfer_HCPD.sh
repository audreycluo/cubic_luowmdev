#!/bin/bash
#SBATCH --job-name=datalad_get_freesurfer_HCPD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/HCPD/freesurfer.out
#SBATCH --error=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/HCPD/freesurfer.err

cd /cbica/projects/luo_wm_dev/input/HCPD/raw/datalad_freesurfer/inputs/data/
missing=$(git annex find --not --in here)
datalad get $missing 

 