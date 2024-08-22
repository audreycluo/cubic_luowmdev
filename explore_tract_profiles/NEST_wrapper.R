library(dplyr)
library(gratia) 
library(mgcv)
library(parallel)
library(rjson)
library(stringr)
library(tidyr)
library(NEST)


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
# this is a custom statFun NEST function for computing delta adjusted rsq
#' statFun.gam.deltaRsq() will be called by NEST_audrey() if NEST argument statFun=="gam.deltaRsq"
#'@param X n x p matrix of image measurements (n = number of subjects, p = number of image locations)
#'@param dat matrix including phenotype and any covariates. make sure all columns in dat are names
#'@param gam.formula user-specified formula for gam fit at each vertex (in each gam, outcome will be a different column of X and the right side of the formula should include variables included in dat matrix)
#'@param lm.formula user-specified formula for lm fit at each vertex (purpose is to get the sign of the delta adjusted R-squared)
#'@param y.in.gam character specifying how the phenotype variable can be identified in gam. if phenotype is included as a linear term in the gam formula, this could just be the variable name. if there's a smooth term for the phenotype, it may be something like s(y)
#'@param y.in.lm character specifying how the phenotype variable can be identified in lm. again, might just be the name of the variable, but could be different (e.g., if testing an interaction term)
#'@param y.permute specify names of columns in dat matrix that should be permuted when getNull = TRUE
#'@param n.cores for parallelization, number of cores to use (default is 1 / no parallelization)
#'@param seed optional to set seed
#'@param n.perm number of permutations to conduct for inference. default is 999 (i.e. minimum p-value 1/1000)
#'@param getNull whether to obtain null distribution vs. just get observed map of statistics. default will be TRUE inside NEST function, but the statFun function will then be recursively called to get null distribution and getNull will then switch to FALSE
#'@export
statFun.gam.deltaRsq = function(X, dat, gam.formula, lm.formula, y.in.gam, y.in.lm, y.permute = NULL, n.cores = 1, seed = NULL, n.perm = 999, getNull = TRUE){
  
  xyz = cbind(X.v = X[,1], # just use the first image location for building model.matrix (will end up refitting at every location though)
              dat)
  
  # initial gam (to get model matrix)
  gam.init = gam(formula = gam.formula, data = data.frame(xyz))
  gam.model.matrix = model.matrix(gam.init)
  
  # weird bug: model matrix is not full rank; remove duplicated column
  if (sum(duplicated(t(gam.model.matrix))) > 0){
    gam.model.matrix = gam.model.matrix[,-which(duplicated(t(gam.model.matrix)))]
  }
  
  columns.to.test = colnames(gam.model.matrix)[grepl(y.in.gam, colnames(gam.model.matrix), fixed = TRUE)]
  
  # pre-specify the permutations (i.e., do the same permutations of people across vertices)
  if (getNull == TRUE){
    if (is.null(y.permute)){
      message("need to specify which columns of `dat` should be permuted")
      return(NULL)
    }
    if (!is.null(seed)){set.seed(seed)}
    perm.ind = lapply(1:n.perm, FUN = function(k){
      sample(1:nrow(X), replace = FALSE)
    })
  }
  
  # Calculate delta adjusted R-squared (ΔR^2adj) for observed data:
  stat.obs = unlist(mclapply(1:ncol(X), FUN = function(v){
    xyz[,"X.v"] = X[,v]
    
    # Full model (with smooth term)
    gam.fullmodel = gam(formula = gam.formula, data = data.frame(xyz))
    gam.fullmodel.results = summary(gam.fullmodel)
    
    # Null model (without smooth term)
    nullmodel.formula = update(gam.formula, . ~ . - s(age, k = 3, fx = TRUE))  # Remove the smooth term
    gam.nullmodel = gam(formula = nullmodel.formula, data = data.frame(xyz))
    gam.nullmodel.results = summary(gam.nullmodel)
    
    # Calculate delta adjusted R-squared (ΔR^2adj)
    adjRsq = gam.fullmodel.results$r.sq - gam.nullmodel.results$r.sq
    
    # Determine the sign using the linear model
    lm.model = lm(lm.formula, data = data.frame(xyz))
    lm.model.t = summary(lm.model)$coefficients[2, 3]  # t-value for the first covariate (assumed to be the smooth term)
    
    if(lm.model.t < 0){
      adjRsq = adjRsq * -1
    }
    
    return(adjRsq)
  }, mc.cores = n.cores))
  
  # If getNull is TRUE, calculate the null distribution of ΔR^2adj
  if (getNull == TRUE){
    stat.null = lapply(1:n.perm, FUN = function(k){
      
      dat.k = dat
      dat.k[,y.permute] = dat.k[perm.ind[[k]], y.permute]
      
      statFun.gam.deltaRsq(X = X,
                           dat = dat.k,
                           gam.formula = gam.formula,
                           lm.formula = lm.formula,
                           y.in.gam = y.in.gam,
                           y.in.lm = y.in.lm, n.cores = n.cores, seed = seed,
                           n.perm = NULL,  # no permutations to do with getNull = FALSE
                           getNull = FALSE)
      
    })
    return(list(T.obs = stat.obs,
                T.null = stat.null))
  } else{
    return(stat.obs)
  }
}


#' Updated NEST function that accepts statFun.gam.deltaRsq
#'@param statFun Function name for the statistical function to be used (e.g., "lm", "gam.mvwald", "gam.deltaRsq").
#'@param args List of arguments required by the selected statFun function.
#'@param net.maps List of vectors (one for each network), where each vector contains 0s and 1s indicating network membership.
#'@param one.sided Logical. If TRUE, performs a one-sided test.
#'@param n.cores Number of cores to use for parallel processing.
#'@param seed Optional seed for reproducibility.
#'@param what.to.return Character vector indicating what outputs to return (e.g., "pval", "ES", "T.obs", "T.null").
#'@export
NEST_audrey = function(statFun, args, net.maps, one.sided = TRUE, n.cores = 1, seed = NULL, what.to.return = c("pval")){
  
  # Check input requirements
  if (!is.list(net.maps)){
    message("net.maps argument should be a list, of vectors (one for each network)")
    return(NULL)
  }
  
  if (ncol(args$X) != length(net.maps[[1]])){
    message("each network should be specified as a vector with length equal to number of columns of X.")
    return(NULL)
  }
  
  if (!identical(as.numeric(sort(unique(unlist(net.maps)))), as.numeric(c(0,1)))){ 
    message("net.maps should include only 0's and 1's (1 = in network/ROI; 0 = outside network/ROI)")
    return(NULL)
  }
  
  # Map brain-phenotype associations
  if (statFun == "lm"){
    required.args = c("X", "y")
    optional.args = c("Z", "type", "FL", "getNull","n.perm")
    
    args = checkArgs(args = args, required.args = required.args, optional.args = optional.args)
    if (!isFALSE(args)){
      statFun.out = statFun.lm(X = args$X,
                               y = args$y,
                               Z = args$Z,
                               type = args$type,
                               n.cores = n.cores,
                               seed = seed,
                               FL = args$FL,
                               n.perm = args$n.perm,
                               getNull = args$getNull)
    } else{
      message("fix args!")
      return(NULL)
    }
  }
  
  if (statFun == "gam.mvwald" || statFun == "gam.deltaRsq"){
    
    # Check args:
    required.args = c("X","dat","gam.formula","lm.formula","y.in.gam","y.in.lm")
    optional.args = c("n.perm")
    
    args = checkArgs(args = args, required.args = required.args, optional.args = optional.args)
    
    if (!isFALSE(args)){
      if (statFun == "gam.mvwald") {
        statFun.out = statFun.gam.mvwald(X = args$X,
                                         dat = args$dat,
                                         gam.formula = args$gam.formula,
                                         lm.formula = args$lm.formula,
                                         y.in.gam = args$y.in.gam,
                                         y.in.lm = args$y.in.lm,
                                         y.permute = args$y.permute,
                                         n.cores = n.cores, seed = seed,
                                         n.perm = args$n.perm,
                                         getNull = TRUE)
      } else if (statFun == "gam.deltaRsq") {
        statFun.out = statFun.gam.deltaRsq(X = args$X,
                                           dat = args$dat,
                                           gam.formula = args$gam.formula,
                                           lm.formula = args$lm.formula,
                                           y.in.gam = args$y.in.gam,
                                           y.in.lm = args$y.in.lm,
                                           y.permute = args$y.permute,
                                           n.cores = n.cores, seed = seed,
                                           n.perm = args$n.perm,
                                           getNull = TRUE)
      }
    } else{
      message("fix args!")
      return(NULL)
    }
  }
  
  ES.list = mclapply(net.maps, FUN = function(net.map){
    ES.obs = enrichScore(stat.map = statFun.out$T.obs,
                         net.map = net.map,
                         one.sided = one.sided,
                         save.detail = F)
    
    ES.null = lapply(1:args$n.perm, FUN = function(k){
      enrichScore(stat.map = statFun.out$T.null[[k]],
                  net.map = net.map,
                  one.sided = one.sided,
                  save.detail = F)
    })
    
    return(list(ES.obs = ES.obs,
                ES.null.dist = ES.null))
  }, mc.cores = n.cores)
  
  pval.list = mclapply(1:length(net.maps), FUN = function(net){
    pvalFun(obs = ES.list[[net]]$ES.obs, null.dist = ES.list[[net]]$ES.null.dist)
  })
  
  out = list()
  
  if ("everything" %in% what.to.return){
    out$pval = pval.list
    out$ES = ES.list
    out$T.obs = statFun.out$T.obs
    out$T.null = statFun.out$T.null
    return(out)
  }else{
    if ("pval" %in% what.to.return){ 
      out$pval = pval.list
    }
    
    if ("ES" %in% what.to.return){
      out$ES = ES.list
    }
    
    if ("ES.obs" %in% what.to.return){
      out$ES.obs = sapply(ES.list, FUN = function(x){x$ES.obs})
    }
    
    if ("T.obs" %in% what.to.return){
      out$T.obs = statFun.out$T.obs
    }
    
    if ("T.null" %in% what.to.return){
      out$T.null = statFun.out$T.null
    }
    return(out)
  }
}

# a wrapper function for running NEST on each tract
NEST_wrapper <- function(tract, dataset, bin_num_nodes) {
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
    if(bin_num_nodes == 10) {
      idx <- c(11:90, 111:190) 
      # 'network of interest': indicate where peripheral WM is (0 = deep, 1 = peripheral)
      peripheral_idx <-  as.numeric(grepl("_(0|[1-9])$", names(df_wide[, -idx]))) # nodes 0-9 = peripheral
    } else if (bin_num_nodes == 15) {
      idx <- c(16:85, 116:185) 
      peripheral_idx <- as.numeric(grepl("_(0?[0-9]|1[0-5])$", names(df_wide[, -idx]))) # nodes 0-15 = peripheral
    } else if (bin_num_nodes == 20) { 
      idx <- c(21:80, 121:180) 
      peripheral_idx <- as.numeric(grepl("_(0?[0-9]|1[0-9]|20)$", names(df_wide[, -idx]))) # nodes 0-20 = peripheral
    }
  } else {
    if(bin_num_nodes == 10) {
      idx <- c(11:40, 61:90, 111:140, 161:190) 
      # 'network of interest': indicate where peripheral WM is (0 = deep, 1 = peripheral)
      peripheral_idx <-  as.numeric(grepl("_(0|[1-9]|9[0-9])$", names(df_wide[, -idx])))  # nodes 0-9 and 90-99 = peripheral
    } else if (bin_num_nodes == 15) {
      idx <- c(16:40, 61:85, 116:140, 161:185) 
      peripheral_idx <- as.numeric(grepl("_(0?[0-9]|1[0-4]|8[5-9]|9[0-9])$", names(df_wide[, -idx]))) # nodes 0-15 and 85-99 = peripheral
    } else if (bin_num_nodes == 20) {
      idx <- c(21:40, 61:80, 121:140, 161:180) 
      peripheral_idx <- as.numeric(grepl("_(0?[0-9]|1[0-9]|8[0-9]|9[0-9])$", names(df_wide[, -idx])))  # nodes 0-20 and 80-99 = peripheral
    }
  }
  
  # final variables for NEST
  X <- as.matrix(df_wide[, -idx])
  colnames(X) <- NULL
  dat <- covs
  net <- list(peripheral_nodes <- peripheral_idx)
  
  print(paste("Running NEST for", tract))
  print(paste("Number of nodes included as peripheral nodes:", bin_num_nodes))
  
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
    n.cores = 8, 
    seed = 123, 
    what.to.return = "everything"
  )
  saveRDS(result_onesided, paste0(NEST_outputs_dir, tract_filename, "_numNodes", bin_num_nodes, "_onesided.RData"))
  
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
    n.cores = 8, 
    seed = 123, 
    what.to.return = "everything"
  )
  
  saveRDS(result_onesided, paste0(NEST_outputs_dir, tract_filename, "_numNodes", bin_num_nodes, "_twosided.RData"))
}

df <- readRDS(sprintf("%1$s/%2$s/tract_profiles/all_subjects/tract_profiles_for_viz.RData", outputs_root, dataset))
assign(dataset, df)
print(paste0(dataset, " df loaded"))

tracts <- unique(df$tract_label)[!str_detect(unique(df$tract_label), "Callosum")] # excluding callosal tracts for now

# Run NEST
lapply(tracts, NEST_wrapper, dataset = dataset, bin_num_nodes = 10)
lapply(tracts, NEST_wrapper, dataset = dataset, bin_num_nodes = 15)
lapply(tracts, NEST_wrapper, dataset = dataset, bin_num_nodes = 20)
 
