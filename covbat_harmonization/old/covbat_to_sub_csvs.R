 
library(data.table)
library(dplyr)
library(rjson)
library(stringr)
library(tidyr)

dataset="HCPD"
config <- fromJSON(file=sprintf("/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_%1$s.json", dataset))
outputs_root <- config$tract_profiles_outputs_root
 
format_covbat <- function(scalar) {
  df <- readRDS(sprintf("/cbica/projects/luo_wm_dev/input/%1$s/%1$s_tractprofiles/all_subjects/collated_tract_profiles_%2$s_covbat.RData",dataset, scalar))
  df <- df %>% pivot_longer(cols = -subjectID, names_to = "tract_node")
  df <- df %>% 
    mutate(tract_hemi = gsub("_[0-9]+", "", tract_node)) %>% 
    mutate(nodeID = str_extract(tract_node, "[0-9]+")) %>%
    mutate(tractID = gsub("_", " ", gsub("Left_|Right_", "", tract_hemi))) %>%
    mutate(hemi = str_extract(tract_hemi, "Left|Right"))
  df$nodeID <- as.numeric(df$nodeID)
  df$subjectID <- as.factor(df$subjectID)
  df <- df %>% select(-tract_node)
  df$tractID <- gsub("Fronto.occipital", "Fronto-occipital", df$tractID)
  names(df)[which(names(df) == "value")] <- paste0("dti_", scalar)
  
  return(df)
}
all_subjects_covbat_fa <- format_covbat("fa")
all_subjects_covbat_md <- format_covbat("md")

# save out individual covbat harmonized csv's for each subject
for(sub in unique(all_subjects_covbat_fa$subjectID)) {
  
  to_save_md <- all_subjects_covbat_md %>% filter(subjectID == sub)
  write.table(to_save_md$dti_md, sprintf("/cbica/projects/luo_wm_dev/input/%1$s/%1$s_tractprofiles/%2$s/%2$s_tract_profiles_dti_md.csv", dataset, sub), col.names=FALSE, row.names=F, sep=",")
  
  to_save_fa <- all_subjects_covbat_fa %>% filter(subjectID == sub)
  write.table(to_save_fa$dti_fa, sprintf("/cbica/projects/luo_wm_dev/input/%1$s/%1$s_tractprofiles/%2$s/%2$s_tract_profiles_dti_fa.csv", dataset, sub), col.names=FALSE, row.names=F, sep=",")
  print(paste(sub, "saved"))
}

all_subjects_covbat <- merge(all_subjects_covbat_fa, all_subjects_covbat_md)
all_subjects_covbat <- all_subjects_covbat %>% arrange(subjectID, tractID, tract_hemi, nodeID) 
write.table(all_subjects_covbat, sprintf("/cbica/projects/luo_wm_dev/input/%1$s/%1$s_tractprofiles/all_subjects/collated_tract_profiles.tsv",dataset), row.names=F, sep=',') 
