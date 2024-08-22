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
outdir_posterior <- paste0(config$tract_profiles_outputs_root, "/GAM/", scalar, "/posterior_derivs/")
 
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

##################################################################################
#################### Compute Posterior Derivatives on Each Node ################## 
################################################################################## 

np <- 200 #number of ages to get the derivatives at
npd <- 1000 #number of posterior draws   

compute_nodewise_derivs <- function(df, scalar) {
  tract_nodes <- df$tract_node
  gam.derivs <- lapply(tract_nodes, function(x){
    gam.smooth.deriv_posterior(dataframe=df, tract_node = as.character(x), smooth_var = "age", covariates = "sex + mean_fd", knots = 3, set_fx = TRUE, 
                          draws = npd, increments = np, return_posterior_derivatives = TRUE)}) #run gam.derivs to get simulated derivs
  names(gam.derivs) <- tract_nodes
  gam.derivs.wide <- lapply(1:length(gam.derivs), function(x) pivot_wider(gam.derivs[[x]], names_from = "draw", values_from = "posterior.derivative", names_sort = FALSE)) # pivot wider for each tract_node
  names(gam.derivs.wide) <- tract_nodes
  gam.derivs.wide <- lapply(gam.derivs.wide, function(x) x[,-c(1:2)]) # remove age and tract_node column
  saveRDS(gam.derivs.wide, sprintf("%1$sderivs_200ages_%2$s.RData", outdir_posterior, scalar), row.names=F)
  
  rm(gam.derivs.wide)
  gc()
  print(paste("Posterior derivatives for", tract_node, "finished"))
}
 


################## 
# Load files
################## 

cohortfile <- read.csv(sprintf("/cbica/projects/luo_wm_dev/output/%1$s/tract_profiles/cohortfiles/dti_fa/dti_fa_cohortfile.csv", dataset)) # just for demographics, scalar doesn't matter
cohortfile <- cohortfile %>% rename(subjectID = sub)
all_subjects <- fread(sprintf("/cbica/projects/luo_wm_dev/input/%1$s/%1$s_tractprofiles/all_subjects/collated_tract_profiles.tsv", dataset))

gam_df <- merge(all_subjects, cohortfile, by="subjectID")
gam_df <- gam_df %>% mutate(tract_node = paste0(tract_hemi, "_", nodeID))
gam_df$sex <- as.factor(gam_df$sex)

print("Files loaded, computing gams next")

#################################################
# Calculate Posterior Derivatives for each Node #
#################################################
compute_nodewise_derivs(gam_df, scalar)

