#!/bin/bash
#SBATCH --job-name=datalad_get_HCPD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/HCPD/qsiprep.out
#SBATCH --error=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/HCPD/qsiprep.err

cd /cbica/projects/luo_wm_dev/input/HCPD/raw/datalad_qsiprep
missing=$(git annex find --not --in here)
datalad get $missing 

 