#!/bin/bash
#SBATCH --job-name=datalad_get_babs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --output=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/babs_mergeds_tractprofiles.out
#SBATCH --error=/cbica/projects/luo_wm_dev/code/tract_profiles/logs/datalad/babs_mergeds_tractprofiles.err

set -x

datasets=("HCPD" "HBN")

for dataset in "${datasets[@]}"; do 
    dest_dir="/cbica/projects/luo_wm_dev/input/${dataset}/derivatives/tract_profiles_dki"
    src_dir="/cbica/projects/luo_wm_dev/input/${dataset}/derivatives/babs_qsirecon_pyafq_dki/merge_ds"

    echo "Processing dataset: ${dataset}"
    echo "Source directory: ${src_dir}"
    echo "Destination directory: ${dest_dir}"

    # create directory if it doesn't exist
    mkdir -p "${dest_dir}" || { echo "Failed to create destination directory"; exit 1; }

    cd "${src_dir}" || { echo "Failed to change directory to ${src_dir}"; exit 1; }

    export dest_dir

    # parallel processing 
    find . -name 'sub*zip' | parallel -j 8 '
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
done
 