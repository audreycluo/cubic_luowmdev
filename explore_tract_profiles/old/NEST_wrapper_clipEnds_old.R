library(dplyr)
library(gratia) 
library(mgcv)
library(parallel)
library(rjson)
library(stringr)
library(tidyr)
library(NEST)
source("/cbica/projects/luo_wm_dev/code/tract_profiles/explore_tract_profiles/NEST_custom_functions.R")

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
print(paste("Running NEST for", dataset))

################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_%1$s.json", dataset))
demographics <- read.csv(config_data$demo_qc)
outputs_root <- "/cbica/projects/luo_wm_dev/output/"
tract_prof_output_root <- config_data$outputs_root
NEST_outputs_dir <- paste0(tract_prof_output_root, "/NEST/")

if (dir.exists(NEST_outputs_dir)) {
  print(paste(NEST_outputs_dir, "already exists"))
} else {
  dir.create(NEST_outputs_dir, recursive = TRUE)
  print(paste(NEST_outputs_dir, "created"))
}

################### 
# Define functions 
###################

# a wrapper function for running NEST on each tract
NEST_wrapper <- function(tract, dataset, bin_size, clipEnds) {
  #clipEnds <- clipEnds - 1 # if want to remove the first 3 nodes, would want to remove nodes 0, 1, and 2. Thus add "-1" for easier coding down the line
  df <- get(paste0(dataset))
  df_tract <- df %>% filter(tract_label == tract)
  df_wide <- df_tract %>% select(sub, tract_node, all_of("dti_md")) %>% pivot_wider(names_from = "tract_node", values_from = "dti_md") 
  tract_filename <- gsub(" ", "_", tract)
  if(!identical(demographics$sub, as.character(df_wide$sub))) {
    stop("Error: The 'sub' columns in 'demographics' and 'df_wide' do not match.")
  }
  
  df_wide <- df_wide %>% select(-sub)
  # set phenotype of interest (age) and covariates (sex + mean_fd)
  demographics$sex <- as.numeric(as.factor(demographics$sex))
  # 1 = F, 2 = M
  covs <- as.matrix(cbind( demographics$age, demographics$sex, as.numeric(demographics$mean_fd))) # matrix dimensions = N x 3
  colnames(covs) <- c("age", "sex", "mean_fd")
  
  
  # idx = identify which locations should be ignored. E.g. nodes that are neither deep nor peripheral
  # i'm keeping the deep wm nodes as the central 10 nodes 
  # but varying which nodes we consider 'peripheral'
  if(tract == "Corticospinal") { # for CST, peripheral = tract by motor cortex, deep = the other end of the tract
    if(bin_size == 5) {
      idx <- c(1:(clipEnds), 11:90, (100-clipEnds+1):100, # only keep the tract ends for CST
               101:(100+clipEnds), 111:190, (200-clipEnds+1):200) 
      
      # indicate where peripheral WM is (0 = deep, 1 = peripheral) ## start here
      peripheral_idx <-  as.numeric(grepl("_([1-9]|1[0-4])$", names(df_wide[, -idx]))) # nodes 5-9 = peripheral
      
    } else if (bin_size == 10) {
      idx <- c(1:(clipEnds), (11+clipEnds):(90-clipEnds), (100-clipEnds+1):100, # only keep the tract ends for CST
               101:(100+clipEnds), (111+clipEnds):(190-clipEnds), (200-clipEnds+1):200) 
      
      # indicate where peripheral WM is (0 = deep, 1 = peripheral) ## start here
      peripheral_idx <-  as.numeric(grepl("_([1-9]|1[0-4])$", names(df_wide[, -idx]))) # nodes 5-14 = peripheral
      
    } else if (bin_size == 15) { 
      idx <- c(1:(clipEnds), (16+clipEnds):(85-clipEnds), (100-clipEnds+1):100, # only keep the tract ends for CST
               101:(100+clipEnds), (116+clipEnds):(185-clipEnds), (200-clipEnds+1):200) 
      
      # indicate where peripheral WM is (0 = deep, 1 = peripheral) ## start here
      peripheral_idx <-  as.numeric(grepl("_([1-9]|1[0-9])$", names(df_wide[, -idx]))) # nodes 5-29 = peripheral
      
    }
  } else if (str_detect(tract, "Callosum")) {
    if(bin_size == 5) {
      idx <- c(1:(clipEnds), 11:40, 61:90, (100-clipEnds+1):100) 
      # 'network of interest': indicate where peripheral WM is (0 = deep, 1 = peripheral)
      peripheral_idx <-  as.numeric(grepl("_([5-9]|9[0-4])$", names(df_wide[, -idx])))  # nodes 5-9 and 90-94 = peripheral
      
    } else if (bin_size == 10) {
      idx <- c(1:(clipEnds), 16:40, 61:85, (100-clipEnds+1):100) 
      # 'network of interest': indicate where peripheral WM is (0 = deep, 1 = peripheral)
      peripheral_idx <-  as.numeric(grepl("_([5-9]|1[0-4]|8[5-9]|9[0-4])$", names(df_wide[, -idx])))  # nodes 5-14 and 85-94 = peripheral
      
    } else if (bin_size == 15) {
      idx <- c(1:(clipEnds), 21:40, 61:80, (100-clipEnds+1):100) 
      # 'network of interest': indicate where peripheral WM is (0 = deep, 1 = peripheral)
      peripheral_idx <-  as.numeric(grepl("_([5-9]|1[0-9]|8[0-9]|9[0-4])$", names(df_wide[, -idx])))  # nodes 5-14 and 80-94 = peripheral
      
    }
  } else {
    if(bin_size == 5) {
      idx <- c(1:(clipEnds), 11:40, 61:90, (100-clipEnds+1):100, 
               101:(100+clipEnds), 111:140, 161:190, (200-clipEnds+1):200) 
      # 'network of interest': indicate where peripheral WM is (0 = deep, 1 = peripheral)
      peripheral_idx <-  as.numeric(grepl("_([5-9]|9[0-4])$", names(df_wide[, -idx])))  # nodes 5-9 and 90-94 = peripheral
      
    } else if (bin_size == 10) {
      idx <- c(1:(clipEnds), 16:40, 61:85, (100-clipEnds+1):100, 
               101:(100+clipEnds), 116:140, 161:185, (200-clipEnds+1):200) 
      # 'network of interest': indicate where peripheral WM is (0 = deep, 1 = peripheral)
      peripheral_idx <-  as.numeric(grepl("_([5-9]|1[0-4]|8[5-9]|9[0-4])$", names(df_wide[, -idx])))  # nodes 5-14 and 85-94 = peripheral
      
    } else if (bin_size == 15) {
      idx <- c(1:(clipEnds), 21:40, 61:80, (100-clipEnds+1):100, 
               101:(100+clipEnds), 121:140, 161:180, (200-clipEnds+1):200) 
      # 'network of interest': indicate where peripheral WM is (0 = deep, 1 = peripheral)
      peripheral_idx <-  as.numeric(grepl("_([5-9]|1[0-9]|8[0-9]|9[0-4])$", names(df_wide[, -idx])))  # nodes 5-14 and 80-94 = peripheral
      
    }
  }
  
  # final variables for NEST
  X <- as.matrix(df_wide[, -idx])
  colnames(X) <- NULL
  dat <- covs
  net <- list(peripheral_nodes <- peripheral_idx)
  
  print(paste("Running NEST for", tract))
  print(paste("Number of nodes included as peripheral nodes:", bin_size))
  
  # one-sided: is delta adj rsq more extreme peripheral vs. deep nodes?
  print("Running one-sided NEST")
  result_onesided <- NEST_audrey(
    statFun = "gam.deltaRsq",
    args = list(X = X, 
                dat = dat, 
                gam.formula = as.formula(X.v ~ s(age, k = 3, fx = TRUE) + sex + mean_fd), 
                lm.formula = as.formula(X.v ~ age + sex + mean_fd),  
                y.in.gam = "s(age)", 
                y.in.lm = "age", 
                y.permute = "age",  
                n.perm = 999),
    net.maps = net,
    one.sided = TRUE,
    n.cores = 12, 
    seed = 123, 
    what.to.return = "everything"
  )
  saveRDS(result_onesided, paste0(NEST_outputs_dir, tract_filename, "_numNodes", bin_size, "_onesided_clipEnds", clipEnds, ".RData"))
  
  # two-sided: is delta adj rsq different peripheral vs. deep nodes?
  print("Running two-sided NEST")
  result_twosided <- NEST_audrey(
    statFun = "gam.deltaRsq",
    args = list(X = X, 
                dat = dat, 
                gam.formula = as.formula(X.v ~ s(age, k = 3, fx = TRUE) + sex + mean_fd), 
                lm.formula = as.formula(X.v ~ age + sex + mean_fd),  
                y.in.gam = "s(age)", 
                y.in.lm = "age", 
                y.permute = "age",  
                n.perm = 999), 
    net.maps = net,
    one.sided = FALSE,
    n.cores = 4, 
    seed = 123, 
    what.to.return = "everything"
  )
  
  saveRDS(result_onesided, paste0(NEST_outputs_dir, tract_filename, "_numNodes", bin_size, "_twosided_clipEnds", clipEnds, ".RData"))
}

df <- readRDS(sprintf("%1$s/%2$s/tract_profiles/all_subjects/tract_profiles_for_viz.RData", outputs_root, dataset))
assign(dataset, df)
print(paste0(dataset, " df loaded"))

tracts <- unique(df$tract_label)  

# Run NEST
lapply(tracts, NEST_wrapper, dataset = dataset, bin_size = 5, clipEnds = 5)
lapply(tracts, NEST_wrapper, dataset = dataset, bin_size = 10, clipEnds = 5)


# # may want to do 9999 for n.perm
# will need to request like 72 hours to run this stuff
