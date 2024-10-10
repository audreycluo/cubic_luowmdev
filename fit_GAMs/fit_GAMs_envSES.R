library(data.table)
library(dplyr)
library(mgcv)
library(parallel)
library(rjson)
library(stringr)
library(tidyr)
source("/cbica/projects/luo_wm_dev/code/tract_profiles/fit_GAMs/gam_functions/GAM_functions_tractprofiles.R")

# This script fits developmental nodewise GAMs on tract profiles data using functions from GAM_functions_tractprofiles.R.
# fyi: GAMs are fit for 1 specified scalar (i.e. dti_md)
# Specifically this script does the following:

# 1) run gam.smooth.predict.covariateinteraction to look at developmental trajectory for 10th and 90th envSES percentiles   
# 2) run gam.fit.covariate to look at the main effect of envSES  

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
scalar = args[2]
print(paste("Fitting age by envSES GAMs for", dataset))
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
# 1) run gam.smooth.predict.covariateinteraction to look at developmental trajectory for 10th and 90th p-factor percentiles
run_gam.smooth.predict.covariateinteraction <- function(gam_df, smooth.var, int.var, covs, k, set.fx, filename, age1, age2){
  smooth.fittedvals.lowerpercentile <- mclapply(tract_node_labels, 
                                                function(x){as.data.frame(gam.smooth.predict.covariateinteraction(gam.data = gam_df,  
                                                                                                                  tract_node = as.character(x), 
                                                                                                                  smooth_var = smooth.var, 
                                                                                                                  int_var = int.var, 
                                                                                                                  int_var.predict = quantile(gam_df[[int.var]], c(0.1), na.rm=T),
                                                                                                                  covariates = covs, 
                                                                                                                  knots = k, 
                                                                                                                  set_fx = set.fx,
                                                                                                                  increments = 200,
                                                                                                                  age1 = age1,
                                                                                                                  age2 = age2)) %>% mutate(tract = x)},  mc.cores = 4) 
  smooth.fittedvals.lowerpercentile <- do.call(rbind, smooth.fittedvals.lowerpercentile)
  smooth.fittedvals.lowerpercentile$envSES <- c("low")
  
  smooth.fittedvals.upperpercentile <- mclapply(tract_node_labels, 
                                                function(x){as.data.frame(gam.smooth.predict.covariateinteraction(gam.data = gam_df,  
                                                                                                                  tract_node = as.character(x), 
                                                                                                                  smooth_var = smooth.var, 
                                                                                                                  int_var = int.var, 
                                                                                                                  int_var.predict = quantile(gam_df[[int.var]], c(0.9), na.rm=T),
                                                                                                                  covariates = covs, 
                                                                                                                  knots = k, 
                                                                                                                  set_fx = set.fx,
                                                                                                                  increments = 200,
                                                                                                                  age1 = age1,
                                                                                                                  age2 = age2)) %>% mutate(tract = x)}, mc.cores = 4) 
  smooth.fittedvals.upperpercentile <- do.call(rbind, smooth.fittedvals.upperpercentile)
  smooth.fittedvals.upperpercentile$envSES <- c("high")
  
  smooth.fittedvals.byenvSES <- rbind(smooth.fittedvals.lowerpercentile, smooth.fittedvals.upperpercentile)
  write.csv(smooth.fittedvals.byenvSES, sprintf("%1$s/%2$s_GAM_ageby%3$s_interaction.csv", GAM_outputs_dir, dataset, filename), quote = F, row.names =F)
}

# 2) run gam.fit.covariate to look at the main effect of envSES (sensitivity analysis)
run_gam.fit.covariate <- function(gam_df, smooth.var, covariate.interest, covariates.noninterest, k, set.fx, filename){
  envSES.maineffects <- mclapply(tract_node_labels, 
                                  function(x){as.data.frame(gam.fit.covariate(gam.data = gam_df,  
                                                                              tract_node = as.character(x), 
                                                                              smooth_var = smooth.var, 
                                                                              covariate.interest = covariate.interest, 
                                                                              covariates.noninterest = covariates.noninterest, 
                                                                              knots = k, 
                                                                              set_fx = set.fx)) %>% mutate(tract = x)}, mc.cores = 4) 
  envSES.maineffects <- do.call(rbind, envSES.maineffects)
  write.csv(envSES.maineffects, sprintf("%1$s/%2$s_GAM_%3$s_maineffects.csv", GAM_outputs_dir, dataset, filename), quote = F, row.names =F)
}

################### 
# Load files  
###################
if (dataset == "PNC" & !file.exists(sprintf("%1$s/all_subjects/collated_tract_profiles_final.RData", data_root))) {
  print("Formatting PNC collated profiles tsv")
  all_subjects <- fread(sprintf("%1$s/all_subjects/collated_tract_profiles_nocovbat.tsv", data_root))
  all_subjects$tractID <- gsub("Fronto-occipital", "Fronto.occipital", all_subjects$tractID)
  all_subjects <- all_subjects %>% mutate(hemi = ifelse(grepl("Left", tractID), "Left", "Right")) %>% 
    mutate(tract_node = gsub(" ", "_", paste0(tractID, "_", nodeID)))
  all_subjects$sub <- as.factor(all_subjects$sub)
  all_subjects <- all_subjects %>% 
    mutate(nodeID = str_extract(tract_node, "[0-9]+")) %>%
    mutate(tractID = gsub("_[0-9]+", "", tract_node)) %>%
    mutate(hemi = str_extract(tractID, "Left|Right"))
  all_subjects$nodeID <- as.numeric(all_subjects$nodeID)
  all_subjects$sub <- as.factor(all_subjects$sub)
  all_subjects <- all_subjects %>% select(sub, tract_node, nodeID, tractID, hemi, dti_fa, dti_md) %>% arrange(sub, tractID, nodeID, hemi)
  saveRDS(all_subjects, sprintf("%1$s/all_subjects/collated_tract_profiles_final.RData", data_root))
  
  tract_profiles_long <- all_subjects
  
} else {
  print("Loading RData file")
  tract_profiles_long <- readRDS(sprintf("%1$s/all_subjects/collated_tract_profiles_final.RData", data_root))
}

################################################################# 
# Merge demographics and qc with tract profiles data and envSES
################################################################# 

# load envSES csv
envSES <- fread(sprintf("/cbica/projects/luo_wm_dev/input/%s/sample_selection_files/%s_envSES.csv", dataset, dataset))
if(dataset=="PNC") {
  bblid_scanid <- read.csv("/cbica/projects/luo_wm_dev/input/PNC/sample_selection_files/bblid_scanid_sub.csv")
  envSES <- merge(envSES, bblid_scanid, by ="scanid")
  envSES$sub <- paste0("sub-", envSES$rbcid)
  demographics <- left_join(demographics, envSES[, c("sub", "envSES")])
} else if (dataset=="HBN") {
  envSES$sub <- paste0("sub-", envSES$GUID)
  demographics <- left_join(demographics, envSES[,c("sub", "General_SES")], by = "sub")
  demographics <- demographics %>% rename(envSES = General_SES)
  
}

# first format and make tract profiles data wide
# such that columns = tract_node, rows = subjects
tract_profiles_wide <- tract_profiles_long %>% select(sub, tract_node, all_of(scalar)) %>% pivot_wider(names_from = "tract_node", values_from = !!sym(scalar))
tract_node_labels = names(tract_profiles_wide)[-1]

# merge
gam_df <- merge(tract_profiles_wide, demographics, by = "sub")
gam_df <- na.omit(gam_df) # exclude subjects if missing p-factor data (0 in PNC, 8 in HBN)

# check envSES and mean_fd
cor.test(gam_df$envSES, gam_df$mean_fd, method = c("pearson")) #PNC: rho = 0.009300684, p-value = 0.7579; HBN: rho = -0.01118563, p-value = 0.7368 

smooth.var = "age"
k = 3
set.fx = T
covs = "mean_fd + sex"

if (dataset == "HBN") {
  age1 = 5
  age2 = 22
} else if (dataset == "HCPD") {
  age1 = 8
  age2 = 22
} else if (dataset == "PNC") {
  age1 = 8
  age2 = 23
}

################### 
# Fit da GAMz 
###################
# 1) run gam.smooth.predict.covariateinteraction to look at developmental trajectory for 10th and 90th p-factor percentiles
# envSES
int.var = "envSES"
filename = "envSES"
run_gam.smooth.predict.covariateinteraction(gam_df, smooth.var, int.var, covs, k, set.fx, filename, age1, age2)

# 2) run gam.fit.covariate to look at the main effect of envSES (sensitivity analysis)
# envSES
covariate.interest = "envSES"
covariates.noninterest = "mean_fd + sex"
filename = "envSES"
run_gam.fit.covariate(gam_df, smooth.var, covariate.interest, covariates.noninterest, k, set.fx, filename)
 
print("Script finished!")


