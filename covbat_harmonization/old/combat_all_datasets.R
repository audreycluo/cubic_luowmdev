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
#args <- commandArgs(trailingOnly = TRUE) 
#dataset = args[1]
#print(paste("Processing", dataset))

################## 
# Set Directories 
################## 
HCPD_config_data <- fromJSON(file=sprintf("/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_%1$s.json", "HCPD"))
HCPD_demographics <- read.csv(HCPD_config_data$demo_qc)
HCPD_demographics$site <- paste0("HCPD", HCPD_demographics$site)
HCPD_data_root <- HCPD_config_data$tract_profiles_root
  
HBN_config_data <- fromJSON(file=sprintf("/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_%1$s.json", "HBN"))
HBN_demographics <- read.csv(HBN_config_data$demo_qc)
HBN_data_root <- HBN_config_data$tract_profiles_root

PNC_config_data <- fromJSON(file=sprintf("/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_%1$s.json", "PNC"))
PNC_demographics <- read.csv(PNC_config_data$demo_qc)
PNC_demographics <- PNC_demographics %>% mutate(site = "PNC")
PNC_data_root <- PNC_config_data$tract_profiles_root

################## 
# Define function 
##################
# @param df, dataframe of combat harmonized data
# @param scalar, string of name of scalar, e.g. "dti_fa"
# this function makes the combat df long and the output should have columns "sub", "tract_node", "nodeID", "tractID", "hemi", and the scalar
format_combat <- function(df, scalar) {
  df_long <- df %>% pivot_longer(cols = -sub, names_to = "tract_node")
  df_long <- df_long %>% 
    mutate(nodeID = str_extract(tract_node, "[0-9]+")) %>%
    mutate(tractID = gsub("_[0-9]+", "", tract_node)) %>%
    mutate(hemi = str_extract(tractID, "Left|Right"))
  df_long$nodeID <- as.numeric(df_long$nodeID)
  df_long$sub <- as.factor(df_long$sub)
  names(df_long)[which(names(df_long) == "value")] <- paste0(scalar)
  return(df_long)
}

################## 
# Load files 
################## 
HCPD_all_subjects <- fread(sprintf("%1$s/all_subjects/collated_tract_profiles_nocovbat.tsv", HCPD_data_root))
HCPD_all_subjects$tractID <- gsub("Fronto-occipital", "Fronto.occipital", HCPD_all_subjects$tractID)
HCPD_all_subjects <- HCPD_all_subjects %>% mutate(hemi = ifelse(grepl("Left", tractID), "Left", "Right")) %>% 
  mutate(tract_node = gsub(" ", "_", paste0(tractID, "_", nodeID))) 
HCPD_all_subjects$sub <- as.factor(HCPD_all_subjects$sub)

HBN_all_subjects <- fread(sprintf("%1$s/all_subjects/collated_tract_profiles_nocovbat.tsv", HBN_data_root))
HBN_all_subjects$tractID <- gsub("Fronto-occipital", "Fronto.occipital", HBN_all_subjects$tractID)
HBN_all_subjects <- HBN_all_subjects %>% mutate(hemi = ifelse(grepl("Left", tractID), "Left", "Right")) %>% 
  mutate(tract_node = gsub(" ", "_", paste0(tractID, "_", nodeID))) 
HBN_all_subjects$sub <- as.factor(HBN_all_subjects$sub)

PNC_all_subjects <- fread(sprintf("%1$s/all_subjects/collated_tract_profiles_nocovbat.tsv", PNC_data_root))
PNC_all_subjects$tractID <- gsub("Fronto-occipital", "Fronto.occipital", PNC_all_subjects$tractID)
PNC_all_subjects <- PNC_all_subjects %>% mutate(hemi = ifelse(grepl("Left", tractID), "Left", "Right")) %>% 
  mutate(tract_node = gsub(" ", "_", paste0(tractID, "_", nodeID))) 
PNC_all_subjects$sub <- as.factor(PNC_all_subjects$sub)

all_subjects <- rbind(HCPD_all_subjects, HBN_all_subjects, PNC_all_subjects)


# df needs to have sub as rows and tract_node as columns. 
dti_md <- all_subjects %>% select(sub, dti_md, tract_node)
dti_md_wide <- dti_md %>% pivot_wider(names_from = "tract_node", values_from = "dti_md")
dti_md_wide <- data.frame(dti_md_wide)

# set rownames = sub; remove sub column for proper combat formatting
row.names(dti_md_wide) <- dti_md_wide$sub
dti_md_to_harmonize <- dti_md_wide %>% select(-sub)

# reorder demographics row to match the wide diffusion metric df's 
demographics <- rbind(HCPD_demographics, HBN_demographics, PNC_demographics)
demographics <- left_join(dti_md_wide[,c(1,2)], demographics, by="sub")
demographics <- demographics %>% select(sub, age, sex, race, mean_fd, site)
# sanity check:
# identical(row.names(dti_md_wide), demographics$sub)


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

data.harmonized_dti_md <- comfam(dti_md_to_harmonize, bat = as.factor(demographics$site), covar = covar_df, model = gam, formula = y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(mean_fd), ref.batch = "PNC")
print("DTI MD harmonized")

# clean up combat output for saving to RData
dti_md_combat <- data.frame(data.harmonized_dti_md$dat.combat)
dti_md_combat$sub <- rownames(dti_md_combat)
dti_md_combat <- dti_md_combat %>% relocate(sub)
rownames(dti_md_combat) <- NULL

# final formatting...
final_dti_md_combat <- format_combat(dti_md_combat, "dti_md")
 
# merge all the combat harmonized data together
merged_combat_all <- final_dti_md_combat %>% arrange(sub, tractID, nodeID, hemi) 

# save out!
saveRDS(merged_combat_all, "/cbica/projects/luo_wm_dev/input/covbat_all_datasets/collated_tract_profiles_combat_all_datasets_site.RData")

