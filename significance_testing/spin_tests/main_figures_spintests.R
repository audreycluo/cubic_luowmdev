library(dplyr)
library(parallel)
library(tidyr)
source("/cbica/projects/luo_wm_dev/code/tract_profiles/results/main_figures_functions.R")

# Spin tests for Figures 6 and 7
## saves out the spun p-values, empirical test statistic, and dof where appropriate in a csv for each dataset

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
print(paste("Running spin tests for", dataset))

################# 
# set directories
################# 
proj_root <- "/cbica/projects/luo_wm_dev/"
input_root <- paste0(proj_root, "input")
output_root <- paste0(proj_root, "output")

output_dir <- paste0(output_root, "/", dataset, "/tract_profiles/fig6_fig7_spintests")
 
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message("Folder created: ", output_dir)
} else {
  message("Folder already exists: ", output_dir)
}

scalar = "dti_md"
ageeffects <- read.csv(sprintf("%1$s/%2$s/tract_profiles/GAM/%3$s/%2$s_GAM_dev_measures.csv", output_root, dataset, scalar))
invisible(assign(paste0(dataset, "_ageeffects"), ageeffects, envir = .GlobalEnv))

# format age effect df's (the functions were originally formatted for all 3 datasets, hence the lapply's)
ageeffect_df_names <- c()
ageeffect_df_names <- append(ageeffect_df_names, paste0(dataset, "_ageeffects"))
ageeffect_dfs <- lapply(ageeffect_df_names, format_ageeffect)
ageeffect.fdr_dfs <- lapply(ageeffect_dfs, sig_nodes)
names(ageeffect.fdr_dfs) <- ageeffect_df_names
tracts = unique(ageeffect.fdr_dfs[[paste0(dataset, "_ageeffects")]]$tract_label)

# load glasser labels
glasser_labels <- read.csv("/cbica/projects/luo_wm_dev/atlases/glasser/HCP-MMP1_UniqueRegionList.csv")
glasser_labels$regionID <- c(1:360)
glasser_labels$region <- gsub("7Pl", "7PL", glasser_labels$region)   

# load maps for depth = 1.5mm
lapply("1.5", load_maps, paste0(dataset)) # makes lh_maps_[dataset], rh_maps_[dataset]

# load maps for depth = 1.0mm
#lapply("1.0", load_maps, paste0(dataset)) # makes lh_maps_[dataset], rh_maps_[dataset]

# load S-A axis
glasser_SAaxis <- read.csv("/cbica/projects/luo_wm_dev/SAaxis/glasser_SAaxis.csv")
glasser_SAaxis <- glasser_SAaxis %>% select(SA.axis_rank, label)
glasser_SAaxis$regionName <- gsub("_ROI", "", glasser_SAaxis$label)
glasser_SAaxis$regionName <- gsub("^(.)_(.*)$", "\\2_\\1", glasser_SAaxis$region)
 
# make cortical endpoint maps for age of maturation
threshold=0.3
deveffects_5_agemat <- ageeffect.fdr_dfs[[paste0(dataset, "_ageeffects")]] %>% filter((nodeID < 10 & nodeID > 4) | (nodeID < 95 & nodeID > 89)) %>% mutate(node_position = case_when((nodeID < 10 & nodeID > 4) ~ "end1", (nodeID < 95 & nodeID > 89) ~ "end2")) %>%
  select(tract_label, tractID, nodeID, node_position, hemi, smooth.peak.change, smooth.decrease.offset, smooth.last.change, smooth.slowing.onset) %>%
  group_by(tractID, node_position, hemi) %>% summarise(mean_peak_change = mean(smooth.peak.change, na.rm = T),
                                                       mean_ageeffect = mean(smooth.decrease.offset, na.rm = T),
                                                       mean_last_change = mean(smooth.last.change, na.rm = T),
                                                       mean_dev_slowing = mean(smooth.slowing.onset, na.rm = T))
invisible(assign(paste0(dataset, "_deveffects_5_agemat"), deveffects_5_agemat, envir = .GlobalEnv))
make_maps(paste0(dataset), "5_agemat") # makes [bundle_name]_deveffect_[dataset] for each dataset and bundle. e.g. IFOL_deveffect_PNC
 
# aggregate age maps
region_all <- aggregate_age_maps(paste0(dataset))  
lh_by_region <- region_all[[1]]
rh_by_region <- region_all[[2]]
invisible(assign(paste0("region_all_", dataset), region_all, envir = .GlobalEnv))
invisible(assign(paste0("lh_by_region_", dataset), lh_by_region, envir = .GlobalEnv))
invisible(assign(paste0("rh_by_region_", dataset), rh_by_region, envir = .GlobalEnv))

# compute SA rank for endpoints
all_endpoints <- compute_mean_SA(paste0(dataset))  
lh_all_endpoints <- all_endpoints[[1]]
rh_all_endpoints <- all_endpoints[[2]]
all_endpoints <- all_endpoints[[3]]
invisible(assign(paste0("all_endpoints_", dataset), all_endpoints, envir = .GlobalEnv))
invisible(assign(paste0("lh_all_endpoints_", dataset), lh_all_endpoints, envir = .GlobalEnv))
invisible(assign(paste0("rh_all_endpoints_", dataset), rh_all_endpoints, envir = .GlobalEnv))

###############################################################
# Fig. 6: 
# For individual datasets: compute mean difference in S-A axis 
# rank between endpoints vs. Difference in age of maturation 
# of endpoints (delta-delta plot)
###############################################################
lh_names <- c("ARCL", "ILFL", "IFOL", "SLFL", "pARCL", "UNCL", "VOFL", "COrbL", "CAntFrL" ,"CSupFrL", "CMotL", "CSupParL", "CPostParL", "CTempL", "COccL")
rh_names <- c("ARCR", "ILFR", "IFOR", "SLFR", "pARCR", "UNCR", "VOFR", "COrbR", "CAntFrR" ,"CSupFrR", "CMotR", "CSupParR", "CPostParR", "CTempR", "COccR")

diffs <- compute_diffs_wrapper(lh_all_endpoints, rh_all_endpoints)  
diffs <- diffs %>% mutate(group = ifelse(str_detect(bundle_name, "IFO") | str_detect(bundle_name, "ILF"), 
                                                 "Large Difference\n in S-A Rank", 
                                                 "Small Difference\n in S-A Rank"))
diffs$group <- factor(diffs$group, levels = c("Small Difference\n in S-A Rank", "Large Difference\n in S-A Rank"))
invisible(assign(paste0("diffs_", dataset), diffs, envir = .GlobalEnv))

print(paste("Delta-delta spin test running for", dataset))
delta_p <- perm.sphere.SAaxis_delta(SAaxis = glasser_SAaxis$SA.axis_rank, perm.id = perm.id.full, dataset) 
delta_p_df <- as.data.frame(delta_p)
write.csv(delta_p_df, paste0(output_dir, "/fig6_pspin_depth1.0.csv"), row.names=F)

########################################################################  
# Fig 7: Spearman's correlation for fully matured endpoints vs. S-A rank
########################################################################  
# binarize age of maturation df's to get matured and not-yet-matured regions
aggregated_axis <- merge_SA_parcel(paste0(dataset)) # 360 glasser parcels with mean age maturation and SA rank

# Remove parcels where age of maturation is the max age (meaning that node/parcel never reached maturation in our age window)
aggregated_axis_binary <- aggregated_axis
max_value <- max(aggregated_axis_binary$regional_mean_ageeffect, na.rm = TRUE)
aggregated_axis_binary$regional_mean_ageeffect[aggregated_axis_binary$regional_mean_ageeffect == max_value] <- NA
invisible(assign(paste0("aggregated_axis_", dataset, "_binary"), aggregated_axis_binary, envir = .GlobalEnv))

# spin for Spearmans
print(paste("Age of maturation vs. SA rank spin test running for", dataset))
parcelSA_p <- perm.sphere.p(x = aggregated_axis_binary$SA.axis_rank, y = aggregated_axis_binary$regional_mean_ageeffect, perm.id = perm.id.full, corr.type = "spearman") 
rho.emp <- cor.test(aggregated_axis_binary$SA.axis_rank, aggregated_axis_binary$regional_mean_ageeffect, method = "spearman", use = "complete.obs")$estimate 
ageMat_SA_p_df <- data.frame(list(p.perm = parcelSA_p, rho.emp = rho.emp), row.names = NULL)
ageMat_SA_p_df <- as.data.frame(ageMat_SA_p_df)
write.csv(ageMat_SA_p_df, paste0(output_dir, "/fig7_spearmans_pspin_depth1.5.csv"), row.names=F)

# spin t-test
aggregated_axis_binary <- aggregated_axis_binary %>% mutate(maturation_status = ifelse(is.na(regional_mean_ageeffect), 0, 1))
parcelSA_p_ttest <- perm.sphere.p.ttest(SAaxis = glasser_SAaxis$SA.axis_rank, perm.id = perm.id.full, dataset = paste0(dataset), alternative = "greater", var.equal = FALSE)  
parcelSA_p_ttest <- as.data.frame(parcelSA_p_ttest)
write.csv(parcelSA_p_ttest, paste0(output_dir, "/fig7_ttest_pspin_depth1.5.csv"), row.names=F)

print("Script finished!")
