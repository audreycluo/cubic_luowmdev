#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1                           
#SBATCH --cpus-per-task=1                   
#SBATCH --time=00:30:00                    
#SBATCH --output=/dev/null    # suppress default output file
#SBATCH --error=/dev/null     # suppress default error file
 
# note: for whatever reason, datalad get takes way longer to run when i run a job array on all participants. 
# If i just do a small array of 3 people, the time limit can be way shorter 
# (despite 1 job array element = 1 participant in both cases). *shrug* weird datalad/cubic quirks i guess.

######################
# setup for job array 
######################
# define variables passed from the submission script
subjects_file=$1  
pyafq_dir=$2
qsiprep_dir=$3
logs_dir=$4
job_name=$5

mapfile -t subjects_array < <(tail -n +2 ${subjects_file}) # skip header
for i in "${!subjects_array[@]}"; do
    subjects_array[$i]=$(echo "${subjects_array[$i]}" | tr -d '"') # remove quotes
done

# get the subject corresponding to my job array element
subject=${subjects_array[$SLURM_ARRAY_TASK_ID]}

# define the log filenames dynamically based on subject
output_file="${logs_dir}/${job_name}_${subject}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
error_file="${logs_dir}/${job_name}_${subject}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# Redirect stdout and stderr to the desired log files
exec > "${output_file}"
exec 2> "${error_file}"
 
echo "Datalad getting qsiprep and pyafq for ${subject}"

######################
# get required files #
######################
cd ${qsiprep_dir}
if find "qsiprep/${subject}" -maxdepth 1 -type d -name "ses-*" | grep -q .; then
    echo "qsiprep dir for subject exists"
else
    # get the dwiref image to use as template for DIPY's load_tractogram
    if datalad get ${subject}*.zip; then
        echo "Qsiprep ${subject} zip successfully gotten"
        unzip -o ${subject}*.zip "qsiprep/${subject}/ses*/dwi/${subject}*T1w_dwiref.nii.gz"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to unzip ${subject} files"
            exit 1   
        fi
        find . -name "${subject}*.zip" -exec datalad drop --nocheck {} \;
        echo "Qsiprep ${subject} zip successfully dropped"
    else
        echo "Error: Failed to get qsiprep ${subject} zip"
        exit 1  
    fi
fi

cd ${pyafq_dir}
if [ ! -d "qsirecon-PYAFQ/${subject}" ]; then
    if datalad get ${subject}*.zip; then
        echo "PyAFQ ${subject} zip successfully gotten"
        unzip -o ${subject}*.zip
        if [ $? -ne 0 ]; then
            echo "Error: Failed to unzip ${subject} files"
            exit 1  
        fi
        find . -name "${subject}*.zip" -exec datalad drop --nocheck {} \;
        echo "PyAFQ ${subject} zip successfully dropped"
    else
        echo "Error: Failed to get PyAFQ ${subject} zip"
        exit 1   
    fi
fi

