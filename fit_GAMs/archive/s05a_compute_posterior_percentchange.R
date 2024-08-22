library(data.table)
library(dplyr)
library(magrittr)
library(rjson)
library(stringr)
source("/cbica/projects/luo_wm_dev/code/functions/GAM_functions.R")


################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly=TRUE)
config_file <- args[1]
scalar <- args[2]
  
################## 
# Set Directories 
################## 
config <- fromJSON(file = config_file)
dataset <- config$dataset
 
outputs_root <- config$tract_profiles_outputs_root
outdir_posterior <- paste0(config$tract_profiles_outputs_root, "/GAM/", scalar, "/posterior_percentchange")
 
if (!dir.exists(outdir_posterior)) {
  # If directory doesn't exist, create it
  dir.create(outdir_posterior, recursive = TRUE)
  print(paste("Directory", outdir_posterior, "created."))
} else {
  print(paste("Directory", outdir_posterior, "already exists."))
}

################## 
# Define Functions
################## 

## Function to calculate percent change of GAM Fitted Value Predictions at min and max age (both raw and normalized)
#function to estimate npd posterior draw fitted values for each tract 
#and compute the percent change of diffusion scalar between min and max age for each tract, for each draw
#executed via rerun below

np <- 200 #number of ages to get the fitted values at
npd <- 10000 #number of posterior draws  

compute_percent_change <- function(gam_df){ 
   
  all_tract_names <- unique(gam_df$tract_hemi)
   
  print("Get posterior predicted fitted values")
  gam.fits.tracts <- map_dfr(all_tract_names, function(x){
    gam.smooth.predict_posterior_tractprofiles(dataframe = gam_df, tract_name = as.character(x), smooth_var = "age", covariates = "sex + mean_fd", set_fx = TRUE, 
                                 draws = npd, increments = np, return_posterior_fits = TRUE, all_ages = FALSE)}) #run gam.fits to get simulated fits
  
  # Compute percent change in DTI measure between min and max age for each node at each draw
  print("Computing percent change")
  percent_change_values <- gam.fits.tracts %>%
    group_by(draw,nodeID,tract_hemi) %>%  
    summarise(percent_change = (posterior.fits[which.max(age)]-posterior.fits[which.min(age)])/(posterior.fits[which.min(age)]))  # percent change for each draw for each tract 
    # percent change is computed  by computing the difference in DTI metric at min and max of a given draw . 
    # So for draw 1 of node 0, compute percent change in DTI metric at min and max age. Then for draw 2 of node 0, compute percent change again, until draw 1000 of node 0.
    # then repeat computations for draw 1 of node 1, draw 2 of node 2, etc. 
   
  percent_change_values.wide <- percent_change_values %>% pivot_wider(names_from = "draw", values_from = "percent_change", names_sort = FALSE)
  print("Finished computing percent change")
  
  print("Computing normalized percent change")
  percent_change_values_norm <- gam.fits.tracts %>%
    group_by(draw,nodeID,tract_hemi) %>%  
    summarise(percent_change = (posterior.fits[which.max(age)]-posterior.fits[which.min(age)])/(nodewise_mean_all_age[which.min(age)]))  # percent change for each draw for each tract 
  # percent change is computed relative to tractwise-nodewise mean (also fyi nodewise_mean_all_age[which.min(age)] and nodewise_mean_all_age[which.max(age)] are equal)
  # For draw 1 of node 0, compute percent change in DTI metric at min and max age. Then for draw 2 of node 0, compute percent change again, until draw 1000 of node 0.
  # then repeat computations for draw 1 of node 1, draw 2 of node 2, etc. 
  
  percent_change_values_norm.wide <- percent_change_values_norm %>% pivot_wider(names_from = "draw", values_from = "percent_change", names_sort = FALSE)
  print("Finished computing normalized percent change")
  
  rm(gam.fits.tracts)
  gc()
   
  saveRDS(percent_change_values.wide, sprintf("%1$s/%2$s_posterior_percentchange.RData", outdir_posterior, dataset))
  saveRDS(percent_change_values_norm.wide, sprintf("%1$s/%2$s_posterior_percentchange_normalized.RData", outdir_posterior, dataset))
   
}

 
 


################## 
# Load files
################## 

cohortfile <- read.csv(sprintf("/cbica/projects/luo_wm_dev/output/%1$s/tract_profiles/cohortfiles/dti_fa/dti_fa_cohortfile.csv", dataset)) # just for demographics, scalar doesn't matter
cohortfile <- cohortfile %>% rename(subjectID = sub)
all_subjects <- fread(sprintf("/cbica/projects/luo_wm_dev/input/%1$s/%1$s_tractprofiles/all_subjects/collated_tract_profiles_%1$s_reoriented.tsv", dataset))

gam_df <- merge(all_subjects, cohortfile, by="subjectID")
gam_df <- gam_df %>% mutate(tract_node = paste0(tract_hemi, "_", nodeID))
gam_df$sex <- as.factor(gam_df$sex)

print("Files loaded, computing gams next")


#########################################################################
# Calculate Posterior Percent Change in DTI Metric (at min and max age) #
#########################################################################
compute_percent_change(gam_df)

