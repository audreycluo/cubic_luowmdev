
###################
#### Setup BABS ###
###################
 

# singularity pull docker://pennbbl/qsiprep:0.22.0
 
# Prepare containerized BIDS App as a DataLad dataset 
cd ~/software/qsiprep

conda activate babs
datalad create -D "qsiprep container 0.22.0" qsiprep-container-0-22-0
cd qsiprep-container-0-22-0
datalad containers-add \
    --url /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep_0.22.0.sif \
    qsiprep-0-22-0


# need to datalad get all the qsiprep zips 
cd /cbica/projects/luo_wm_dev/input/HCPD/datalad_qsiprep
datalad get *zip # i ended up dropping these files after babs-init bc they take up so much space. NOPE that was a bad idea, need the files datlad getted

cd ..
conda activate babs


###################
# Initialize BABS #
###################

# HCPD
babs-init --where_project /cbica/projects/luo_wm_dev/input/HCPD \
    --project_name babs_qsirecon_pyafq --input qsiprep /cbica/projects/luo_wm_dev/input/HCPD/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/HCPD/subject_list/HCPD_subject_list_babs.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container \
    --container_name qsiprep-0-19-1 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-19-1_qsirecon_mrtrix_pyafq.yaml \
    --type_session single-ses --type_system sge 

# HCPD - redo with dki
babs-init --where_project /cbica/projects/luo_wm_dev/input/HCPD \
    --project_name babs_qsirecon_pyafq_dki --input qsiprep /cbica/projects/luo_wm_dev/input/HCPD/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/HCPD/subject_list/HCPD_subject_list_babs.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container-0-21-4 \
    --container_name qsiprep-0-21-4 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-21-4_qsirecon_mrtrix_pyafq_dki.yaml \
    --type_session single-ses --type_system sge 

# HCPD - regular test subject
babs-init --where_project /cbica/projects/luo_wm_dev/input/HCPD \
    --project_name babs_qsirecon_pyafq_test --input qsiprep /cbica/projects/luo_wm_dev/input/HCPD/raw/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/HCPD/subject_list/HCPD_subject_list_babs_test.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container \
    --container_name qsiprep-0-19-1 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-19-1_qsirecon_mrtrix_pyafq.yaml \
    --type_session single-ses --type_system sge 


# HCPD - dki with test subjects, slurm
babs-init --where_project /cbica/projects/luo_wm_dev/input/HCPD/derivatives \
    --project_name babs_qsirecon_pyafq_dki_test --input qsiprep /cbica/projects/luo_wm_dev/input/HCPD/raw/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/HCPD/subject_list/HCPD_subject_list_babs_test.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container-0-22-0 \
    --container_name qsiprep-0-22-0 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-22-0_qsirecon_mrtrix_pyafq_dki.yaml \
    --type_session single-ses --type_system slurm 

# HCPD - dki with test subjects, sge TT_TT
babs-init --where_project /cbica/projects/luo_wm_dev/input/HCPD/derivatives \
    --project_name babs_qsirecon_pyafq_dki --input qsiprep /cbica/projects/luo_wm_dev/input/HCPD/raw/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/HCPD/subject_list/HCPD_subject_list_babs.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container-0-22-0 \
    --container_name qsiprep-0-22-0 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-22-0_qsirecon_mrtrix_pyafq_dki.yaml \
    --type_session single-ses --type_system sge 
 

 
##########

# PNC - doing this with qsiprep 0.21.4
babs-init --where_project /cbica/projects/luo_wm_dev/input/PNC \
    --project_name babs_qsirecon_pyafq --input qsiprep /cbica/projects/luo_wm_dev/input/PNC/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/PNC/sample_selection_files/PNC_WMDev_FinalSample_N1153_age8to22.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container-0-21-4 \
    --container_name qsiprep-0-21-4 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-21-4_qsirecon_mrtrix_pyafq_singleshell.yaml \
    --type_session single-ses --type_system sge # note that yaml code will have to be edited after babs init to include the correct json
 
# PNC 0.22.0
babs-init --where_project /cbica/projects/luo_wm_dev/input/PNC/derivatives \
    --project_name babs_qsirecon_pyafq_dki --input qsiprep /cbica/projects/luo_wm_dev/input/PNC/raw/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/PNC/subject_list/PNC_subject_list_babs.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container-0-22-0 \
    --container_name qsiprep-0-22-0 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-22-0_qsirecon_mrtrix_pyafq_singleshell_slurm.yaml \
    --type_session single-ses --type_system slurm 

# HBN
babs-init --where_project /cbica/projects/luo_wm_dev/input/HBN \
    --project_name babs_qsirecon_pyafq --input qsiprep /cbica/projects/luo_wm_dev/input/HBN/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/HBN/sample_selection_files/HBN_WMDev_TempSample_N1276_age5to22.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container \
    --container_name qsiprep-0-19-1 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-19-1_qsirecon_mrtrix_pyafq.yaml \
    --type_session single-ses --type_system sge 
 

babs-init --where_project /cbica/projects/luo_wm_dev/input/HBN/derivatives \
    --project_name babs_qsirecon_pyafq_dki --input qsiprep /cbica/projects/luo_wm_dev/input/HBN/raw/datalad_qsiprep \
    --list_sub_file /cbica/projects/luo_wm_dev/input/HBN/sample_selection_files/HBN_WMDev_TempSample_N1276_age5to22.txt \
    --container_ds /cbica/projects/luo_wm_dev/software/qsiprep/qsiprep-container-0-22-0 \
    --container_name qsiprep-0-22-0 --container_config_yaml_file /cbica/projects/luo_wm_dev/code/tract_profiles/run_babs_qsirecon/babs_qsiprep-0-22-0_qsirecon_mrtrix_pyafq_dki_slurm.yaml \
    --type_session single-ses --type_system slurm 
 

 

# note:
# "Registering the input dataset(s)...
# Cloning input dataset #1: 'qsiprep'
# [INFO   ] RIA store unavailable.""
# this error is okay on cubic for cloning pmacs static datasets

cd babs_qsirecon_pyafq_test/
babs-status --project-root $PWD # check status - must be in root of babs project
head analysis/code/participant_job.sh
babs-check-setup --project-root ${PWD} --job-test # success!

babs-submit --project-root $PWD   # submit 1 participant as test
babs-status --project-root $PWD  # to check status of that job - this takes several hours

# after test job finishes successfully: 
babs-submit --project-root $PWD --all 

# submit failed
babs-status \
    --project-root $PWD \
    --resubmit failed  

babs-merge --project-root $PWD

# resubmit a test subject - only works if the subject failed though
babs-status \
    --project-root $PWD   \
    --resubmit-job sub-NDARPT777FDA

babs-status \
    --project-root $PWD   \
    --resubmit-job sub-NDARPF682GDC

########################
# Consume results!!!!! #
########################
cd /cbica/projects/luo_wm_dev/input/HCPD
datalad clone \
    ria+file://${PWD}/babs_qsirecon_pyafq/output_ria#~data \
    qsirecon_pyafq


########################
# Failed BABS projects #
########################
cd /cbica/projects/luo_wm_dev/input/HBN/

 
chmod -R +w babs_qsirecon_pyafq 
rm -rf babs_qsirecon_pyafq