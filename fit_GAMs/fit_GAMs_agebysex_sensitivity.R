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

# 1) run gam.factorsmooth.interaction to look at the sex-by-age interaction for each node (sensitivity analysis)
# 2) run gam.fit.covariate to look at the main effect of sex (sensitivity analysis)
# 3) run gam.smooth.predict.covariateinteraction.factor to look at developmental trajectory for female and males

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
scalar = args[2]
print(paste("Fitting age by sex GAMs for", dataset))
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
# 1) run gam.factorsmooth.interaction to look at the sex-by-age interaction for each node (sensitivity analysis)
run_gam.factorsmooth.interaction <- function(gam_df, smooth.var, int.var, covs, k, set.fx){
  agebysex.interactioneffects <- mclapply(tract_node_labels, 
                                         function(x){as.data.frame(gam.factorsmooth.interaction(gam.data = gam_df,  
                                                                                                tract_node = as.character(x), 
                                                                                                smooth_var = smooth.var, 
                                                                                                int_var = int.var, 
                                                                                                covariates = covs, 
                                                                                                knots = k, 
                                                                                                set_fx = set.fx)) %>% mutate(tract = x)}, mc.cores = 4) 
  agebysex.interactioneffects <- do.call(rbind, agebysex.interactioneffects)
  write.csv(agebysex.interactioneffects, sprintf("%1$s/%2$s_GAM_agebysex_interaction.csv", GAM_outputs_dir, dataset), quote = F, row.names =F)
}

# 2) run gam.fit.covariate to look at the main effect of sex (sensitivity analysis)
run_gam.fit.covariate <- function(gam_df, smooth.var, covariate.interest, covariates.noninterest, k, set.fx){
  sex.maineffects <- mclapply(tract_node_labels, 
                                          function(x){as.data.frame(gam.fit.covariate(gam.data = gam_df,  
                                                                                      tract_node = as.character(x), 
                                                                                      smooth_var = smooth.var, 
                                                                                      covariate.interest = covariate.interest, 
                                                                                      covariates.noninterest = covariates.noninterest, 
                                                                                      knots = k, 
                                                                                      set_fx = set.fx)) %>% mutate(tract = x)}, mc.cores = 4) 
  sex.maineffects <- do.call(rbind, sex.maineffects)
  write.csv(sex.maineffects, sprintf("%1$s/%2$s_GAM_sex_maineffects.csv", GAM_outputs_dir, dataset), quote = F, row.names =F)
}


# 3) run gam.smooth.predict.covariateinteraction.factor to look at developmental trajectory for female and males
run_gam.smooth.predict.sex_interaction <- function(gam_df, smooth_var, int_var, covs, k, set_fx, age1, age2){
  sex_levels <- c("F", "M")   
  
  results_list <- list()
  for (sex in sex_levels) {
    results <- lapply(tract_node_labels, function(x) {
      as.data.frame(gam.smooth.predict.covariateinteraction.factor(
        gam.data = gam_df, 
        tract_node = as.character(x), 
        smooth_var = smooth_var, 
        int_var = int_var, 
        int_var.predict = sex,
        covariates = covs, 
        knots = k, 
        set_fx = set_fx,
        increments = 200,
        age1 = age1,
        age2 = age2)) %>% mutate(tract = x, sex = sex)
    })
    
    results <- do.call(rbind, results)
    results_list[[sex]] <- results
  }
  
  all_results <- do.call(rbind, results_list)
  write.csv(all_results, sprintf("%1$s/%2$s_GAM_agebysex_interaction_fits.csv", GAM_outputs_dir, dataset), quote = F, row.names = F)
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

######################################################### 
# Merge demographics and qc with tract profiles data
#########################################################
# first format and make tract profiles data wide
# such that columns = tract_node, rows = subjects
tract_profiles_wide <- tract_profiles_long %>% select(sub, tract_node, all_of(scalar)) %>% pivot_wider(names_from = "tract_node", values_from = !!sym(scalar))
tract_node_labels = names(tract_profiles_wide)[-1]
# merge
gam_df <- merge(tract_profiles_wide, demographics, by = "sub")
gam_df$sex.ordered <- factor(gam_df$sex, levels = c("M", "F"), ordered = TRUE)
smooth.var = "age"
k = 3
set.fx = T
int.var = "sex.ordered" # might need to edit this
covs = "mean_fd"


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
# 1) run gam.factorsmooth.interaction to look at the sex-by-age interaction for each node (sensitivity analysis)
run_gam.factorsmooth.interaction(gam_df, smooth.var, int.var, covs, k, set.fx)

# 2) run gam.fit.covariate to look at the main effect of sex (sensitivity analysis)
covariate.interest = "sex"
covariates.noninterest = "mean_fd"
run_gam.fit.covariate(gam_df, smooth.var, covariate.interest, covariates.noninterest, k, set.fx)

# 3) run gam.smooth.predict.covariateinteraction.factor to look at developmental trajectories for female and males
covs = "sex.ordered + mean_fd"
run_gam.smooth.predict.sex_interaction(gam_df, smooth.var, int.var, covs, k, set.fx, age1, age2)


print("Script finished!")