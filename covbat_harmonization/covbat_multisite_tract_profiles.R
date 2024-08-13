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
print(paste("Processing", dataset))

################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_%1$s.json", dataset))
demographics <- read.csv(config_data$demo_qc)
data_root <- config_data$tract_profiles_root
  
################## 
# Define function 
##################
# @param df, dataframe of covbat harmonized data
# @param scalar, string of name of scalar, e.g. "dki_fa"
# this function makes the covbat df long and the output should have columns "sub", "tract_node", "nodeID", "tractID", "hemi", and the scalar
format_covbat <- function(df, scalar) {
  df_long <- df %>% pivot_longer(cols = -sub, names_to = "tract_node")
  df_long <- df_long %>% 
    mutate(nodeID = str_extract(tract_node, "[0-9]+")) %>%
    mutate(tractID = gsub("_[0-9]+", "", tract_node)) %>%
    mutate(hemi = str_extract(tractID, "Left|Right"))
  df_long$nodeID <- as.numeric(df_long$nodeID)
  df_long$sub <- as.factor(df_long$sub)
  df_long$tractID <- gsub("Fronto.occipital", "Fronto-occipital", df_long$tractID)
  names(df_long)[which(names(df_long) == "value")] <- paste0(scalar)
  return(df_long)
}

################## 
# Load files 
################## 
all_subjects <- fread(sprintf("%1$s/all_subjects/collated_tract_profiles_nocovbat.tsv", data_root))
all_subjects <- all_subjects %>% mutate(hemi = ifelse(grepl("Left", tractID), "Left", "Right")) %>% 
  mutate(tract_node = gsub(" ", "_", paste0(tractID, "_", nodeID)))
all_subjects$sub <- as.factor(all_subjects$sub)

# df needs to have sub as rows and tract_node as columns. 
dki_fa <- all_subjects %>% select(sub, dki_fa, tract_node)
dki_fa_wide <- dki_fa %>% pivot_wider(names_from = "tract_node", values_from = "dki_fa")
dki_fa_wide <- data.frame(dki_fa_wide)

dti_fa <- all_subjects %>% select(sub, dti_fa, tract_node)
dti_fa_wide <- dti_fa %>% pivot_wider(names_from = "tract_node", values_from = "dti_fa")
dti_fa_wide <- data.frame(dti_fa_wide)

dki_md <- all_subjects %>% select(sub, dki_md, tract_node)
dki_md_wide <- dki_md %>% pivot_wider(names_from = "tract_node", values_from = "dki_md")
dki_md_wide <- data.frame(dki_md_wide)

dti_md <- all_subjects %>% select(sub, dti_md, tract_node)
dti_md_wide <- dti_md %>% pivot_wider(names_from = "tract_node", values_from = "dti_md")
dti_md_wide <- data.frame(dti_md_wide)

# set rownames = sub; remove sub column for proper covbat formatting
row.names(dki_fa_wide) <- dki_fa_wide$sub
dki_fa_to_harmonize <- dki_fa_wide %>% select(-sub)

row.names(dti_fa_wide) <- dti_fa_wide$sub
dti_fa_to_harmonize <- dti_fa_wide %>% select(-sub)

row.names(dki_md_wide) <- dki_md_wide$sub
dki_md_to_harmonize <- dki_md_wide %>% select(-sub)

row.names(dti_md_wide) <- dti_md_wide$sub
dti_md_to_harmonize <- dti_md_wide %>% select(-sub)

# reorder demographics row to match the wide diffusion metric df's 
demographics <- left_join(dti_md_wide[,c(1,2)], demographics, by="sub")
demographics <- demographics %>% select(sub, age, sex, race, site, mean_fd)

# set covariate vectors
age_vec <- demographics$age 
sex_vec <- as.factor(demographics$sex) 
mean_fd_vec <- demographics$mean_fd 
covar_df <- bind_cols(demographics$sub, as.numeric(age_vec), as.factor(sex_vec), as.numeric(mean_fd_vec))
covar_df <- dplyr::rename(covar_df, sub=...1,
                          age = ...2,
                          sex = ...3,
                          mean_fd = ...4)

################## 
# Harmonize data 
################## 
data.harmonized_dki_fa <- covfam(dki_fa_to_harmonize, bat = as.factor(demographics$site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd))
print("DKI FA harmonized")

data.harmonized_dti_fa <- covfam(dti_fa_to_harmonize, bat = as.factor(demographics$site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd))
print("DTI FA harmonized")
 
data.harmonized_dki_md <- covfam(dki_md_to_harmonize, bat = as.factor(demographics$site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd))
print("DKI MD harmonized")

data.harmonized_dti_md <- covfam(dti_md_to_harmonize, bat = as.factor(demographics$site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd))
print("DTI MD harmonized")


# clean up covbat output for saving to RData
dki_fa_covbat <- data.frame(data.harmonized_dki_fa$dat.covbat)
dki_fa_covbat$sub <- rownames(dki_fa_covbat)
dki_fa_covbat <- dki_fa_covbat %>% relocate(sub)
rownames(dki_fa_covbat) <- NULL

dti_fa_covbat <- data.frame(data.harmonized_dti_fa$dat.covbat)
dti_fa_covbat$sub <- rownames(dti_fa_covbat)
dti_fa_covbat <- dti_fa_covbat %>% relocate(sub)
rownames(dti_fa_covbat) <- NULL
 
dki_md_covbat <- data.frame(data.harmonized_dki_md$dat.covbat)
dki_md_covbat$sub <- rownames(dki_md_covbat)
dki_md_covbat <- dki_md_covbat %>% relocate(sub)
rownames(dki_md_covbat) <- NULL

dti_md_covbat <- data.frame(data.harmonized_dti_md$dat.covbat)
dti_md_covbat$sub <- rownames(dti_md_covbat)
dti_md_covbat <- dti_md_covbat %>% relocate(sub)
rownames(dti_md_covbat) <- NULL

# final formatting...
final_dki_fa_covbat <- format_covbat(dki_fa_covbat, "dki_fa")
final_dti_fa_covbat <- format_covbat(dti_fa_covbat, "dti_fa")
final_dki_md_covbat <- format_covbat(dki_md_covbat, "dki_md")
final_dti_md_covbat <- format_covbat(dti_md_covbat, "dti_md")
 
# merge all the covbat harmonized data together
merged_fa <- merge(final_dki_fa_covbat, final_dti_fa_covbat)
merged_md <- merge(final_dki_md_covbat, final_dti_md_covbat)
merged_covbat_all <- merge(merged_fa, merged_md)
merged_covbat_all <- merged_covbat_all %>% arrange(sub, tractID, nodeID, hemi) 

# save out!
saveRDS(merged_covbat_all, sprintf("%1$s/all_subjects/collated_tract_profiles_final.RData", data_root))
 


