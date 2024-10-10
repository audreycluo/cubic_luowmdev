#!/bin/bash

# define variables passed from the submission script
subject=$1
pyafq_dir=$2


rm -rf ${pyafq_dir}/qsirecon-PYAFQ/${subject}
rm ${pyafq_dir}/qsirecon-PYAFQ/${subject}.html