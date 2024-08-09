library(ComBatFamily)
library(data.table)
library(dplyr)
library(mgcv)
library(rjson)
library(stringr)
library(tidyr)


################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]


################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_%1$s.json", dataset))
demographics <- read.csv(config_data$demographics_file)
qc_file <- read.csv(config_data$qc_file)
demographics <- merge(demographics, qc_file, by = "sub")

data_root <- config_data$tract_profiles_data_root
  

################## 
# Load files 
################## 
all_subjects <- fread(sprintf("%1$s/all_subjects/collated_tract_profiles_nocovbat.tsv", data_root))

# df needs to have subjectID as rows and tract_node as columns. 
all_subjects <- all_subjects %>%
  mutate(
    hemi = ifelse(grepl("_L", tractID), "Left", "Right"),
    tractID = gsub("_[LR]", "", tractID),
    tractID = case_when(
      tractID == "ATR" ~ "Anterior Thalamic Radiation",
      tractID == "CGC" ~ "Cingulum Cingulate",
      tractID == "CST" ~ "Corticospinal Tract",
      tractID == "IFO" ~ "Inferior Fronto-occipital Fasciculus",
      tractID == "ILF" ~ "Inferior Longitudinal Fasciculus",
      tractID == "SLF" ~ "Superior Longitudinal Fasciculus",
      tractID == "ARC" ~ "Arcuate Fasciculus",
      tractID == "UNC" ~ "Uncinate Fasciculus",
      tractID == "FA" ~ "Forceps Minor",
      tractID == "FP" ~ "Forceps Major",
      tractID == "pARC" ~ "Posterior Arcuate",
      tractID == "VOF" ~ "Vertical Occipital Fasciculus",
      TRUE ~ tractID
    ),
    tract_hemi = paste(hemi, tractID),
    tract_hemi = gsub("Right Forceps", "Forceps", tract_hemi),
    tract_hemi = str_replace_all(tract_hemi, " ", "_")
  )

all_subjects$hemi[all_subjects$tract_hemi == "Forceps_Minor" | all_subjects$tract_hemi == "Forceps_Major"] <- NA
all_subjects$subjectID <- as.factor(all_subjects$subjectID)
all_subjects <- all_subjects %>% mutate(tract_node = paste0(tract_hemi, "_", nodeID))

md <- all_subjects %>% select(subjectID, dti_md, tract_node)
md_wide <- md %>% pivot_wider(names_from = "tract_node", values_from = "dti_md")
md_wide <- data.frame(md_wide)

fa <- all_subjects %>% select(subjectID, dti_fa, tract_node)
fa_wide <- fa %>% pivot_wider(names_from = "tract_node", values_from = "dti_fa")
fa_wide <- data.frame(fa_wide)

# set rownames = subjectID; remove subjectID column 
row.names(md_wide) <- md_wide$subjectID
md_to_harmonize <- md_wide %>% select(-subjectID)

row.names(fa_wide) <- fa_wide$subjectID
fa_to_harmonize <- fa_wide %>% select(-subjectID)

# if dataset = HCPD, will need to temporarily remove sub-0197045's FP data since they were missing FP (NA's), and covbat does not like NAs!!
if(dataset == "HCPD") {
  NA_indices <- which(is.na(md_to_harmonize), arr.ind= TRUE)
  md_to_harmonize <- na.omit(md_to_harmonize)
  fa_to_harmonize <- na.omit(fa_to_harmonize)
  md_wide <- md_wide[-c(which(md_wide$subjectID == "sub-0197045")),]
 
}

# reorder demographics to match md_wide and fa_wide
demographics <- demographics %>% rename(subjectID=sub)
demographics <- left_join(md_wide[,c(1,2)], demographics, by="subjectID")
demographics <- demographics %>% select(subjectID, age, sex, race, site, mean_fd)

# set vectors
age_vec <- demographics$age 
sex_vec <- as.factor(demographics$sex) 

mean_fd_vec <- demographics$mean_fd 

covar_df <- bind_cols(demographics$subjectID, as.numeric(age_vec), as.factor(sex_vec), as.numeric(mean_fd_vec))
covar_df <- dplyr::rename(covar_df, subjectID=...1,
                          age = ...2,
                          sex = ...3,
                          mean_fd = ...4)

################## 
# Harmonize data 
################## 

data.harmonized_fa <- covfam(fa_to_harmonize, bat = as.factor(demographics$site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd))
print("FA harmonized")
data.harmonized_md <- covfam(md_to_harmonize, bat = as.factor(demographics$site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd))
print("MD harmonized") 

fa_covbat <- data.frame(data.harmonized_fa$dat.covbat)
fa_covbat$subjectID <- rownames(fa_covbat)
fa_covbat <- fa_covbat %>% relocate(subjectID)
rownames(fa_covbat) <- NULL
 
md_covbat <- data.frame(data.harmonized_md$dat.covbat)
md_covbat$subjectID <- rownames(md_covbat)
md_covbat <- md_covbat %>% relocate(subjectID)
rownames(md_covbat) <- NULL

saveRDS(fa_covbat, sprintf("%1$s/all_subjects/collated_tract_profiles_fa_covbat.RData", data_root))
saveRDS(md_covbat, sprintf("%1$s/all_subjects/collated_tract_profiles_md_covbat.RData", data_root))


##########
md <- all_subjects %>% 
  group_by(subjectID, tract) %>% 
  summarize(dti_md = mean(dti_md, na.rm = TRUE))

fa <- all_subjects %>% 
  group_by(subjectID, tract) %>% 
  summarize(dti_fa = mean(dti_fa, na.rm = TRUE))

md_wide <- md %>% pivot_wider(names_from = "tract", values_from = "dti_md")
fa_wide <- fa %>% pivot_wider(names_from = "tract", values_from = "dti_fa")

row.names(md_wide) <- md_wide$subjectID
md_to_harmonize <- md_wide %>% select(-subjectID)

row.names(fa_wide) <- fa_wide$subjectID
fa_to_harmonize <- fa_wide %>% select(-subjectID)

if(dataset == "HCPD") {
  md_to_harmonize <- na.omit(md_to_harmonize)
  fa_to_harmonize <- na.omit(fa_to_harmonize)
  md_wide <- md_wide[-c(which(md_wide$subjectID == "sub-0197045")),]
}

demographics <- left_join(md_wide[,c(1)], demographics, by="subjectID")
demographics <- demographics %>% select(subjectID, age, sex, race, site, mean_fd)

age_vec <- demographics$age 
sex_vec <- as.factor(demographics$sex) 
mean_fd_vec <- demographics$mean_fd 

covar_df <- bind_cols(demographics$subjectID, as.numeric(age_vec), as.factor(sex_vec), as.numeric(mean_fd_vec))
covar_df <- dplyr::rename(covar_df, subjectID=...1, age = ...2, sex = ...3, mean_fd = ...4)




harmonize_tractwise <- function(data, bat, covar, model, formula, score.model = NULL, score.args = NULL, eb = TRUE, robust.LS = FALSE, ref.batch = NULL, percent.var = 0.95, n.pc = NULL, std.var = TRUE, ...) {
  results <- list()
  
  for (tract in colnames(data)) {
    cat("Harmonizing tract:", tract, "\n")
    
    tract_data <- data[[tract]]
    harmonized <- covfam(tract_data, bat = as.factor(bat), covar = covar, model = model, formula = formula, score.model = score.model, score.args = score.args, eb = eb, robust.LS = robust.LS, ref.batch = ref.batch, percent.var = percent.var, n.pc = n.pc, std.var = std.var, ...)
    
    results[[tract]] <- harmonized$dat.covbat
  }
  
  harmonized_data <- do.call(cbind, results)
  return(harmonized_data)
}

data.harmonized_fa <- harmonize_tractwise(fa_to_harmonize, bat = demographics$site, covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd) + as.factor(site))
print("FA harmonized")
data.harmonized_md <- harmonize_tractwise(md_to_harmonize, bat = demographics$site, covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd) + as.factor(site))
print("MD harmonized")
