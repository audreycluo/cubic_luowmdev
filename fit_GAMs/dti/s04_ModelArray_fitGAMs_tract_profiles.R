
library(ModelArray)

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly=TRUE)
config_file <- args[1]



################## 
# Set Directories 
################## 
config <- read.csv(config_file)
outputs_root <- config$tract_profiles_outputs_root
outdir_fa <-  paste0(outputs_root, "/GAM/dti_fa")
outdir_md <-  paste0(outputs_root, "/GAM/dti_md")

if (!dir.exists(outdir_fa)) {
  # If directory doesn't exist, create it
  dir.create(outdir_fa, recursive = TRUE)
  print(paste("Directory", outdir_fa, "created."))
} else {
  print(paste("Directory", outdir_fa, "already exists."))
}

if (!dir.exists(outdir_md)) {
  # If directory doesn't exist, create it
  dir.create(outdir_md, recursive = TRUE)
  print(paste("Directory", outdir_md, "created."))
} else {
  print(paste("Directory", outdir_md, "already exists."))
}

################## 
# Define Functions
################## 
fitGAMs_ModelArray <- function(scalar, outdir) {
    
  # Create a ModelArray-class object
  # filename of dsistudio scalar data (.h5 file):
  h5_path <- sprintf("%1$s/h5_files/%2$s.h5", outputs_root, scalar)
   
  
  # create a ModelArray-class object:
  modelarray <- ModelArray(h5_path, scalar_types = c(scalar))
  # let's check what's in it:
  #modelarray
  #scalars(modelarray)[["dti_fa"]]
  
  phenotypes <- read.csv(paste0(outputs_root, "/cohortfiles/", scalar, "/", scalar, "_cohortfile.csv"))
  
  # formula:
  smooth_var <- "age"
  covariates <- "sex + mean_fd"
  formula.gam <- as.formula(sprintf("%s ~ s(%s, k=3, fx=T) + %s", scalar, smooth_var, covariates))
 
  
  print(paste0("fitting GAMs for ", scalar))
  
  mygam.try <- ModelArray.gam(
    formula.gam,
    modelarray,
    phenotypes,
    scalar,
    changed.rsq.term.index = c(1)
  )
   
  saveRDS(mygam.try, sprintf("%1$s/GAMresults.%2$s.age.RData", outdir, scalar))
  write.csv(mygam.try, sprintf("%1$s/GAMresults.%2$s.age.csv", outdir, scalar))
}

 
fitLinearModel_ModelArray <- function(scalar, outdir) {
  
  # Create a ModelArray-class object
  # filename of dsistudio scalar data (.h5 file):
  h5_path <- sprintf("%1$s/h5_files/%2$s.h5", outputs_root, scalar)
  
  # create a ModelArray-class object:
  modelarray <- ModelArray(h5_path, scalar_types = c(scalar))
  phenotypes <- read.csv(paste0(outputs_root, "/cohortfiles/", scalar, "/", scalar, "_cohortfile.csv"))
  
  # formula:
  covariates <- "age + sex + mean_fd"
  formula.lm <- as.formula(sprintf("%s ~  %s", scalar, covariates))
  
  print(paste0("fitting linear model for ", scalar))
  
  mylm.try <- ModelArray.lm(
    formula.lm,
    modelarray,
    phenotypes,
    scalar 
  )
   
  write.csv(mylm.try, sprintf("%1$s/LinearModelresults.%2$s.age.csv", outdir, scalar))
}


  
################## 
# Fit GAMs!
################## 
fitGAMs_ModelArray("dti_fa", outdir_fa)
fitGAMs_ModelArray("dti_md", outdir_md)

################################################################################# 
# Fit linear model to determine sign of age effect, as computed by delta adj R2 #
#################################################################################
fitLinearModel_ModelArray("dti_fa", outdir_fa)
fitLinearModel_ModelArray("dti_md", outdir_md)
 