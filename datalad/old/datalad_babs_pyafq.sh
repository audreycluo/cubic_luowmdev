#!/bin/bash

dataset="HBN"
# qsub -l h_vmem=60G,s_vmem=60G datalad_babs_pyafq.sh  
# clone from babs project
#cd /cbica/projects/luo_wm_dev/input/${dataset}

#datalad clone \
#   ria+file://${PWD}/babs_qsirecon_pyafq/output_ria#~data \
#   qsirecon_pyafq
# i dont think i need this step. i can just use the merge-ds folder
    
#mkdir /cbica/projects/luo_wm_dev/input/${dataset}/${dataset}_tractprofiles

cd /cbica/projects/luo_wm_dev/input/${dataset}/babs_qsirecon_pyafq/merge_ds

# for now, unzip and copy just the tract profiles csv's
for file in sub*zip ; do 
    sub=${file%_*}
    if ! [ -d /cbica/projects/luo_wm_dev/input/${dataset}/${dataset}_tractprofiles/$sub ]; then

        datalad get $file  
        mkdir /cbica/projects/luo_wm_dev/input/${dataset}/${dataset}_tractprofiles/$sub  
        sub=${file%_*}   

        
        datalad get $file  

        # this unzips just the one file I want
        unzip -j "$file" "qsirecon/$sub/*/dwi/*/*profiles_dwi.csv" -d /cbica/projects/luo_wm_dev/input/${dataset}/${dataset}_tractprofiles/$sub 
        datalad drop $file   # drop the zip
    fi
done
 