#!/bin/bash
#SBATCH --job-name=datalad_get_freesurfer_PNC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/PNC/freesurfer_%j.out
#SBATCH --error=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/PNC/freesurfer_%j.err

cd /cbica/projects/luo_wm_dev/input/PNC/raw/datalad_freesurfer/inputs/data/
# datalad get -n .
missing_files=$(git annex find --not --in here)

# some participants in my list don't have freesurfer outputs
no_zip_file_list="/cbica/projects/luo_wm_dev/input/PNC/sample_selection_files/missing_freesurfer.txt"
> $no_zip_file_list

while IFS= read -r participant; do
    zip_file_pattern="${participant}_freesurfer-20.2.3.zip"
    existing_zip_file=$(find . -name "$zip_file_pattern" | head -n 1)

    if [ -z "$existing_zip_file" ]; then 
        echo "${participant} has no zip file in the folder." >> $no_zip_file_list; 
    elif echo "$missing_files" | grep -q "$zip_file_pattern"; then
        echo "${participant}'s needs to be gotten. Attempting to download..."
        datalad get "$zip_file_pattern"
    fi
done < /cbica/projects/luo_wm_dev/input/PNC/subject_list/PNC_subject_list_babs.txt

echo "Participants with no zip file have been logged in: $no_zip_file_list"
 
