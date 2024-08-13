#!/bin/bash
#SBATCH --job-name=babs_mergeds_tractprofiles
#SBATCH --ntasks=3                          # Number of tasks (parallel processes)
#SBATCH --cpus-per-task=10                  # Number of CPU cores per task
#SBATCH --array=1-3                         # Array job, one task per dataset
#SBATCH --output="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/babs_mergeds_tractprofiles_%A_%a.out"
#SBATCH --error="/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/babs_mergeds_tractprofiles_%A_%a.err"

datasets=("HCPD" "HBN" "PNC")

# get the dataset corresponding to this task
dataset=${datasets[$SLURM_ARRAY_TASK_ID-1]}
 
dest_dir="/cbica/projects/luo_wm_dev/input/${dataset}/derivatives/tract_profiles"
src_dir="/cbica/projects/luo_wm_dev/input/${dataset}/derivatives/babs_qsirecon_pyafq_dki/merge_ds"

echo "Processing dataset: ${dataset}"
echo "Source directory: ${src_dir}"
echo "Destination directory: ${dest_dir}"

# create directory if it doesn't exist
mkdir -p "${dest_dir}" || { echo "Failed to create destination directory"; exit 1; }

cd "${src_dir}" || { echo "Failed to change directory to ${src_dir}"; exit 1; }

export dest_dir

# parallel processing 
find . -name 'sub*zip' | parallel -j 10 '
    file={}
    sub=$(basename ${file%_*})

    echo "Processing file: $file"
    echo "Destination directory: ${dest_dir}"
    echo "Subject: ${sub}"

    if ! [ -d "${dest_dir}/${sub}" ]; then
        echo "Getting ${file}"
        datalad get "${file}"
        mkdir -p "${dest_dir}/${sub}"
        
        # unzip only the tract profiles csv
        unzip -j "${file}" "qsirecon-PYAFQ/${sub}/*/dwi/*/*profiles_dwi.csv" -d "${dest_dir}/${sub}"
        
        datalad drop "${file}"  # drop the zip
    fi
' dest_dir="${dest_dir}"
 
 