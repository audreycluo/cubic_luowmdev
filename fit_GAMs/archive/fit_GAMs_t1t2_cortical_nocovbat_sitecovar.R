library(data.table)
library(dplyr)
library(mgcv)
library(parallel)
library(rjson)
library(stringr)
library(tidyr)
source("/cbica/projects/luo_wm_dev/code/tract_profiles/fit_GAMs/gam_functions/GAM_functions_tractprofiles.R")

# This script fits developmental regional GAMs on T1/T2 data from Graham's paper using functions from GAM_functions_tractprofiles.R.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9302463/ 
# fyi: GAMs are fit for 1 specified scalar (i.e. dti_md)
# Specifically this script does the following:

# 1) run gam.fit.smooth to compute developmental measures + model checks:
# - Age effect: partial Rsq and delta adjusted Rsq
# - Age of maximal developmental change
# - Age of decrease offset: age of maturation if DWI measure decreases with age (e.g. MD)
# - Age of increase offset: age of maturation if DWI measure increases with age (e.g. FA)
# - Age of last change
# - Age of first developmental slowing
# - Model AIC, k-index for checking the model

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
scalar = args[2]
print(paste("Fitting t1/t2 cortical GAMs for", dataset))
print(paste("Scalar:", scalar))

################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_%1$s.json", dataset))
demographics <- read.csv(config_data$demo_qc)
data_root <- config_data$tract_profiles_root
outputs_root <- config_data$outputs_root
GAM_outputs_dir <- paste0(outputs_root, "/GAM/", scalar)

if (dir.exists(GAM_outputs_dir)) {
  print(paste(GAM_outputs_dir, "already exists"))
} else {
  dir.create(GAM_outputs_dir, recursive = TRUE)
  print(paste(GAM_outputs_dir, "created"))
}

################### 
# Define functions 
###################
# 1) run gam.fit.smooth to compute developmental measures + model checks
run_gam.fit.smooth <- function(gam_df, smooth.var, covs, k, set.fx) {
  GAM_dev_measures <- map_dfr(glasser_labels$region, 
                                    function(x){as.data.frame(gam.fit.smooth.cortex(gam.data = gam_df, 
                                                                             region = as.character(x), 
                                                                             smooth_var = smooth.var, 
                                                                             covariates = covs, 
                                                                             knots = k, 
                                                                             set_fx = set.fx))}) 
  #GAM_dev_measures <- do.call(rbind, GAM_dev_measures)
  write.csv(GAM_dev_measures, sprintf("%1$s/%2$s_GAM_t1t2_measures_age_mat.csv", GAM_outputs_dir, dataset), quote = F, row.names =F)
}
 
################### 
# Load files  
###################
glasser_labels <- read.csv("/cbica/projects/luo_wm_dev/atlases/glasser/HCP-MMP1_UniqueRegionList.csv")
glasser_labels$regionID <- c(1:360)
glasser_labels$region <- gsub("7Pl", "7PL", glasser_labels$region)   
glasser_labels$region <- gsub("-", ".", glasser_labels$region) # remove hyphens for formatting reasons
glasser_labels$region[1:180] <- paste0("L_", glasser_labels$region[1:180]) # label left hemi
glasser_labels$region[181:360] <- paste0("R_", glasser_labels$region[181:360]) # label right hemi


t1t2 <- read.csv("/cbica/projects/luo_wm_dev/input/cortical_data/t1t2_hcpd/n628_hcpd_newCorr_myelin_Aug2021.csv")
t1t2$subject_id <- gsub("HCD", "sub-", t1t2$subject_id)
t1t2 <- t1t2 %>% rename(sub=subject_id)
t1t2_glasser <- t1t2[,c(1:361)]
names(t1t2_glasser)[2:181] <- glasser_labels$region[1:180]
names(t1t2_glasser)[182:361] <- glasser_labels$region[181:360]
t1t2_covariates <- t1t2 %>% select(sub, Sex, Scanner, Reference.Voltage, Mean.Pseudo.Transmit.Map, T2..Dropout.Threshold, Smoothing.FWHM.mm, Pseudotransmit.Reference.Value, Correction.Equation.Slope, Corrected.CSF.Regressor)

######################################################### 
# Merge demographics and qc with tract profiles data
#########################################################
# merge
demographics <- demographics %>% select(sub, age, mean_fd)
gam_df <- left_join(demographics, t1t2_glasser, by = "sub")
gam_df <- left_join(gam_df, t1t2_covariates, by = "sub")
gam_df <- gam_df %>% relocate(age, mean_fd, .after = last_col())
gam_df <- gam_df[-c(which(rowSums(is.na(gam_df)) > 1)),] # exclude 6 subjects for missing T1w/T2w data
gam_df$Sex <- as.factor(gam_df$Sex)

################### 
# Fit GAMs
###################
# set variables for GAM
smooth.var = "age"
# for covariates, see https://pmc.ncbi.nlm.nih.gov/articles/PMC9302463/#sec22 and https://github.com/andlab-harvard/hcpd_myelin/blob/main/07_brms_regional_posterior_deriv_analysis.R#L74C31-L74C222
covs = "Sex + Scanner + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"

k = 3 
set.fx = TRUE

# 1) run gam.fit.smooth to compute temporal measures + model checks
print("Running gam.fit.smooth")
run_gam.fit.smooth(gam_df, smooth.var, covs, k, set.fx)

print("Script finished!")
