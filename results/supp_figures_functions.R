library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(gratia) 
library(knitr)
#library(kableExtra)
library(mgcv)
library(RColorBrewer)
library(scales)
library(stringr)
library(rjson)
library(tidyr)

########################
# Supplementary Figs
########################
# spin test for tract-level comparison to the S-A axis. This function spins the S-A axis then recomputes average S-A rank for each end. 
# Option to do pearson's (all endpoints) or spearman's (only matured included in correlation, plus a spun t-test)
# @param SAaxis, vector of S-A axis ranks
perm.sphere.SAaxis <- function(SAaxis, perm.id, dataset, spun_ttest = FALSE, alternative = "greater", var.equal = FALSE) {
  
  nroi = dim(perm.id)[1]  # number of regions
  nperm = dim(perm.id)[2] # number of permutations
  
  # spin SA axis: permutation of measures
  SAaxis.perm = array(NA,dim=c(nroi,nperm))
  for (r in 1:nperm) {
    for (i in 1:nroi) {
      SAaxis.perm[i,r] = SAaxis[perm.id[i,r]]
    }
  }
  # if just want to compute spun p-value for age of maturation vs mean S-A rank for all cortical endpoints: 
  if (spun_ttest == FALSE) { 
    # empirical correlation
    all_endpoints <- get(paste0("all_endpoints_", dataset))
    rho.emp = cor(all_endpoints$mean_SA, all_endpoints$age_effect, method = "pearson", use="complete.obs")  
    
    # correlation between permuted S-A axis and age of maturation
    rho.null.SAaxis = vector(length=nperm)
    
    # when averaging across datasets, use PNC's bundle-to-cortex probability map
    if(dataset=="avg_datasets") {
      dataset <- "PNC"
    }
    
    for (r in 1:nperm) {
      print(r)
      # merging of permuted S-A axis to age of maturation df
      perm_mean_SAaxis <- compute_mean_SA(dataset, SAaxis.perm[,r]) # merge permuted S-A rank with age of maturation AND compute mean permuted S-A ranks
      perm_mean_SAaxis <- perm_mean_SAaxis[[3]] # get all_endpoints df
      rho.null.SAaxis[r] = cor(perm_mean_SAaxis$mean_SA, perm_mean_SAaxis$age_effect, method="pearson", use="complete.obs")
    }
    # p-value definition depends on the sign of the empirical correlation
    if (rho.emp>0) {
      p.perm = sum(rho.null.SAaxis>rho.emp)/nperm
    } else { 
      p.perm = sum(rho.null.SAaxis<rho.emp)/nperm
    } 
    return(list(p.perm = p.perm, rho.emp = rho.emp))
    # return p-value for the correlation between age of maturation vs. mean S-A rank for cortical endpoint
    
  } else if(spun_ttest == TRUE) {
    
    all_endpoints <- get(paste0("all_endpoints_", dataset))
    all_endpoints_binary <- all_endpoints
    max_value <- max(all_endpoints_binary$age_effect, na.rm = TRUE)
    all_endpoints_binary$age_effect[all_endpoints_binary$age_effect == max_value] <- NA
    all_endpoints_binary <- all_endpoints_binary %>% mutate(maturation_status = ifelse(is.na(age_effect), 0, 1)) # 0 = not matured, 1 = matured
    
    # empirical correlation
    rho.emp = cor(all_endpoints_binary$mean_SA, all_endpoints_binary$age_effect, method = "spearman", use="complete.obs")  
    
    # empirical t value
    t.emp <- t.test(all_endpoints_binary$mean_SA[which(all_endpoints_binary$maturation_status == 0)], # compare differences in sa rank
                    all_endpoints_binary$mean_SA[which(all_endpoints_binary$maturation_status == 1)], 
                    alternative = alternative, var.equal = var.equal)$statistic # hypothesis is non-matured tract ends have greater S-A rank
    
    df.emp <- t.test(all_endpoints_binary$mean_SA[which(all_endpoints_binary$maturation_status == 0)], # compare differences in sa rank
                    all_endpoints_binary$mean_SA[which(all_endpoints_binary$maturation_status == 1)], 
                    alternative = alternative, var.equal = var.equal)$parameter 
    
    # set vector for correlations between permuted S-A axis and age of maturation
    rho.null.SAaxis = vector(length=nperm)
    # set vector for t-tests between permuted S-A axis differences and age of maturation differences between endpoints
    t.null.SAaxis = vector(length=nperm)
    
    # when averaging across datasets, use PNC's bundle-to-cortex probability map
    if(dataset=="avg_datasets") {
      dataset <- "PNC"
    }
    
    for (r in 1:nperm) {
      print(r)
      # merging of permuted S-A axis to age of maturation df
      perm_mean_SAaxis <- compute_mean_SA(dataset, SAaxis.perm[,r]) # merge permuted S-A rank with age of maturation AND compute mean permuted S-A ranks
      perm_mean_SAaxis <- perm_mean_SAaxis[[3]] # get all_endpoints df
      max_value_perm <- max(perm_mean_SAaxis$age_effect, na.rm = TRUE)
      perm_mean_SAaxis$age_effect[perm_mean_SAaxis$age_effect == max_value_perm] <- NA
      rho.null.SAaxis[r] = cor(perm_mean_SAaxis$mean_SA, perm_mean_SAaxis$age_effect, method="spearman", use="complete.obs")
      
      perm_mean_SAaxis <- perm_mean_SAaxis %>% mutate(maturation_status = ifelse(is.na(age_effect), 0, 1)) # 0 = not matured, 1 = matured
      
      # calculate t-test and permuted t-test: do tract ends that have NOT matured between endpoints actually have greater s-a ranks than tract ends that have matured?? 
      t.null.SAaxis[r] <- t.test(perm_mean_SAaxis$mean_SA[which(perm_mean_SAaxis$maturation_status == 0)], 
                                 perm_mean_SAaxis$mean_SA[which(perm_mean_SAaxis$maturation_status == 1)], 
                                 alternative = alternative, var.equal = var.equal)$statistic
    }
    # p-value definition depends on the sign of the empirical correlation
    if (rho.emp>0) {
      p.perm.cor = sum(rho.null.SAaxis>rho.emp)/nperm
    } else { np
      p.perm.cor = sum(rho.null.SAaxis<rho.emp)/nperm
    } 
    
    # p-value definition for t-test
    if (t.emp>0) {
      p.perm.t = sum(t.null.SAaxis>t.emp)/nperm
    } else { 
      p.perm.t= sum(t.null.SAaxis<t.emp)/nperm
    } 
    
    return(list(rho.emp = rho.emp, p.perm.cor = p.perm.cor, t.emp = t.emp, p.perm.t = p.perm.t, df.emp = df.emp))
  }
}

# make summary df's: this creates summary_[dataset]_[scalar] dataframes that include mean, sd, se for each scalar/dataset
make_summary_dfs <- function(scalar, dataset_name, df) {
  # set variables
  mean_scalar <- paste0("mean_", scalar)
  # Process the given dataset df
  scalar_data <- df %>% select(sub, tractID, tract_label, nodeID, hemi, all_of(scalar))
  sqrt_n <- sqrt(length(unique(scalar_data$sub)))
  setDT(scalar_data)
  summary_data <- scalar_data[, .(
    mean_scalar = mean(.SD[[1]], na.rm = TRUE),
    sd = sd(.SD[[1]], na.rm = TRUE),
    se = sd(.SD[[1]], na.rm = TRUE) / sqrt_n,
    ymin_sd = mean(.SD[[1]], na.rm = TRUE) - sd(.SD[[1]], na.rm = TRUE),
    ymax_sd = mean(.SD[[1]], na.rm = TRUE) + sd(.SD[[1]], na.rm = TRUE),
    ymin_se = mean(.SD[[1]], na.rm = TRUE) - (sd(.SD[[1]], na.rm = TRUE) / sqrt_n),
    ymax_se = mean(.SD[[1]], na.rm = TRUE) + (sd(.SD[[1]], na.rm = TRUE) / sqrt_n),
    cv = (sd(.SD[[1]], na.rm = TRUE) / (mean(.SD[[1]], na.rm = TRUE))) * 100
  ), by = .(tract_label, tractID, nodeID, hemi), .SDcols = scalar]
  
  setnames(summary_data, "mean_scalar", paste0("mean_", scalar))
  summary_data[, Dataset := dataset_name]
  summary_dataset_scalar <- paste0("summary_", dataset_name, "_", scalar)
  assign(summary_dataset_scalar, summary_data, envir = .GlobalEnv)
}

# hex plot function for correlations between datasets
hex_plot <- function(df, x, y, text, ylim1, ylim2, xlim1, xlim2, x_text, y_text) {
  plot <- ggplot(df, aes_string(x = x, y = y)) +
    geom_hex(bins = 7) +
    paletteer::scale_fill_paletteer_c("ggthemes::Blue-Green Sequential", direction = -1) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +   
    annotate("text", x = x_text, y = y_text, 
             label = text, 
             hjust = 0, vjust = 1, size = 7, color = "black") + 
    theme(legend.position = "bottom", 
          legend.key.width = unit(2, 'cm'), legend.key.height = unit(1, 'cm'),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.margin = unit(c(0.2, 0.5, 0.2, 1), "cm")) + 
    labs(x = x, y = y, fill = "Count") + ylim(ylim1, ylim2) + xlim(xlim1, xlim2)
  return(plot)
}


# function for permutation-based correlation test
perm_cor_test <- function(vector1, vector2, n_perm = 10000, seed = 123) {
  # empirical correlation
  rho.emp <- cor(vector1, vector2, use = "complete.obs")
  
  # null correlation vector
  rho.perm <- numeric(n_perm)
  
  # set seed
  set.seed(seed)
  
  # permute
  for (i in 1:n_perm) {
    permuted_vector2 <- sample(vector2)  # shuffle one vector
    rho.perm[i] <- cor(vector1, permuted_vector2, use = "complete.obs")
  }
  
  # compute p-value (two-tailed test)
  exceed_count <- sum(abs(rho.perm) >= abs(rho.emp))
  p_value <- (1 + exceed_count) / (1 + n_perm)
  
  return(list(rho.emp = rho.emp, p_value = p_value))
}



plot_meanSA_by_age_mat <- function(dataset, annot_text, binary=NULL) {
  title = gsub("PD", "P-D", dataset) # for HCPD
  
  if(is.null(binary)) {
    all_endpoints <- get(paste0("all_endpoints_", dataset))
    
    SA_plot <- ggplot(all_endpoints, aes(x = mean_SA, y = age_effect, label = bundle_name)) +
      geom_point(aes(fill = mean_SA, color = mean_SA, shape = end), size = 4) + 
      paletteer::scale_fill_paletteer_c("grDevices::ag_Sunset", direction = -1, 
                                        limits = c(min(all_endpoints$mean_SA), max(all_endpoints$mean_SA)), oob = squish) + 
      paletteer::scale_color_paletteer_c("grDevices::ag_Sunset", direction = -1,
                                         limits = c(min(all_endpoints$mean_SA), max(all_endpoints$mean_SA)), oob = squish) +
      scale_shape_manual(values = c(19,1)) + 
      geom_smooth(data = all_endpoints, method='lm', se=TRUE, fill=alpha(c("gray70"),.9), col="black") + 
      ggrepel::geom_text_repel(
        size = 5, 
        position = position_jitter(seed = 3), 
        box.padding = 0.5,  # Increase space around labels
        point.padding = 0.5,  # Increase space around points
        segment.color = "gray",  # Line color connecting point to label
        segment.size = 0.5,  # Line thickness
        max.overlaps = Inf 
      )  + 
      annotate(geom="text", x=180, y=27, label=annot_text, color="black", size=8) + 
      labs(title = title) + theme_classic() +
      theme(legend.position = "none",
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24),
            plot.title = element_text(hjust = 0.5, size = 24),
            axis.title = element_blank(),
            axis.text.x = element_text(size = 24),
            axis.text.y = element_text(size = 24)) + xlim(15, 340) + ylim(13, 27) + guides(shape = guide_legend("Endpoint"), color = FALSE, fill = FALSE) 
  } else {
    all_endpoints <- get(paste0("all_endpoints_", dataset))
    all_endpoints_binary <- all_endpoints
    max_value <- max(all_endpoints_binary$age_effect, na.rm = TRUE)
    all_endpoints_binary$age_effect[all_endpoints_binary$age_effect == max_value] <- NA
    all_endpoints_binary <- all_endpoints_binary %>% mutate(maturation_status = ifelse(is.na(age_effect), 0, 1)) # 0 = not matured, 1 = matured
    
    SA_plot <- ggplot(all_endpoints_binary, aes(x = mean_SA, y = age_effect, label = bundle_name)) +
      geom_point(aes(fill = mean_SA, color = mean_SA, shape = end), size = 4) + 
      paletteer::scale_fill_paletteer_c("grDevices::ag_Sunset", direction = -1, 
                                        limits = c(min(all_endpoints$mean_SA), max(all_endpoints$mean_SA)), oob = squish) + 
      paletteer::scale_color_paletteer_c("grDevices::ag_Sunset", direction = -1,
                                         limits = c(min(all_endpoints$mean_SA), max(all_endpoints$mean_SA)), oob = squish) +
      scale_shape_manual(values = c(19,1)) + 
      geom_smooth(data = all_endpoints_binary, method='lm', se=TRUE, fill=alpha(c("gray70"),.9), col="black") + 
      ggrepel::geom_text_repel(
        size = 5, 
        position = position_jitter(seed = 3), 
        box.padding = 0.5,  # space around labels
        point.padding = 0.5,  # space around points
        segment.color = "gray",  # line color connecting point to label
        segment.size = 0.5,  # line thickness
        max.overlaps = Inf 
      )  + 
      annotate(geom="text", x=180, y=27, label=annot_text, color="black", size=8) + 
      labs(title = title) + theme_classic() +
      theme(legend.position = "none",
            legend.text = element_text(size = 24),
            legend.title = element_text(size = 24),
            plot.title = element_text(hjust = 0.5, size = 24),
            axis.title = element_blank(),
            axis.text.x = element_text(size = 24),
            axis.text.y = element_text(size = 24)) + xlim(15, 340) + ylim(13, 27) + guides(shape = guide_legend("Endpoint"), color = FALSE, fill = FALSE) 
    
  }
  
  return(SA_plot)
}



plot_ageeffects <- function(bundle_name, dataset, ylim1, ylim2) {
  df <- get(paste0(bundle_name, "_deveffect_", dataset))
  if(grepl("L$", bundle_name)) {
    hemi = "left"
    cortical_pos1 <- "left lateral" 
    cortical_pos2 <- "left medial"
  } else{
    hemi = "right"
    cortical_pos1 <- "right lateral" 
    cortical_pos2 <- "right medial"
  }
  plot_lateral <- ggplot() + 
    geom_brain(data = df, atlas= glasser, 
               mapping=aes(fill=mean_age_effect), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos1)) +
    paletteer::scale_fill_paletteer_c("grDevices::RdYlBu", direction = -1 , limits = c(ylim1, ylim2), oob = squish, na.value = "white") +
    #scale_fill_gradientn(colors = aquamarine, na.value = "white", limits = c(ylim1, ylim2), oob=squish) +
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -1, 0.5, -1), "cm"),
          plot.title = element_blank()) 
  
  plot_medial <- ggplot() + 
    geom_brain(data = df, atlas= glasser, 
               mapping=aes(fill=mean_age_effect), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos2)) +
    paletteer::scale_fill_paletteer_c("grDevices::RdYlBu", direction = -1 , limits = c(ylim1, ylim2), oob = squish, na.value = "white") +
    #scale_fill_gradientn(colors = aquamarine, na.value = "white", limits = c(ylim1, ylim2), oob=squish) +
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -1, 0.5, -1), "cm"),
          plot.title = element_blank())
  
  return( plot_grid(plot_lateral, plot_medial))
  
}

arrange_ageeffect_plots <- function(lh_plots, rh_plots) {
  # arrange maps together  
  COrb_plots <- ggarrange(lh_plots$COrbL, rh_plots$COrbR, ncol=1, nrow = 2)
  COrb_final <- annotate_figure(COrb_plots, top = text_grob("Callosum Orbital", 
                                                            color = "black", size = 20))
  
  CAntFr_plots <- ggarrange(lh_plots$CAntFrL, rh_plots$CAntFrR, ncol=1, nrow = 2)
  CAntFr_final <- annotate_figure(CAntFr_plots, top = text_grob("Callosum Anterior Frontal", 
                                                                color = "black", size = 20))
  
  CSupFr_plots <- ggarrange(lh_plots$CSupFrL, rh_plots$CSupFrR, ncol=1, nrow = 2)
  CSupFr_final <- annotate_figure(CSupFr_plots, top = text_grob("Callosum Superior Frontal", 
                                                                color = "black", size = 20))
  
  CMot_plots <- ggarrange(lh_plots$CMotL, rh_plots$CMotR, ncol=1, nrow = 2)
  CMot_final <- annotate_figure(CMot_plots, top = text_grob("Callosum Motor", 
                                                            color = "black", size = 20))
  
  CSupPar_plots <- ggarrange(lh_plots$CSupParL, rh_plots$CSupParR, ncol=1, nrow = 2)
  CSupPar_final <- annotate_figure(CSupPar_plots, top = text_grob("Callosum Superior Parietal", 
                                                                  color = "black", size = 20))
  
  CPostPar_plots <- ggarrange(lh_plots$CPostParL, rh_plots$CPostParR, ncol=1, nrow = 2)
  CPostPar_final <- annotate_figure(CPostPar_plots, top = text_grob("Callosum Posterior Parietal", 
                                                                    color = "black", size = 20))
  
  CTemp_plots <- ggarrange(lh_plots$CTempL, rh_plots$CTempR, ncol=1, nrow = 2)
  CTemp_final <- annotate_figure(CTemp_plots, top = text_grob("Callosum Temporal", 
                                                              color = "black", size = 20))
  
  
  COcc_plots <- ggarrange(lh_plots$COccL, rh_plots$COccR, ncol=1, nrow = 2)
  COcc_final <- annotate_figure(COcc_plots, top = text_grob("Callosum Occipital", 
                                                            color = "black", size = 20))
  
  ARC_plots <- ggarrange(lh_plots$ARCL, rh_plots$ARCR, ncol=1, nrow = 2)
  ARC_final <- annotate_figure(ARC_plots, top = text_grob("Arcuate Fasciculus", 
                                                          color = "black", size = 20))
  
  
  CST_plots <- ggarrange(lh_plots$CSTL, rh_plots$CSTR, ncol=1, nrow = 2)
  CST_final <- annotate_figure(CST_plots, top = text_grob("Corticospinal Tract", 
                                                          color = "black", size = 20))
  
  IFO_plots <- ggarrange(lh_plots$IFOL, rh_plots$IFOR, ncol=1, nrow = 2)
  IFO_final <- annotate_figure(IFO_plots, top = text_grob("Inferior Fronto-occipital Fasciculus", 
                                                          color = "black", size = 20))
  
  ILF_plots <- ggarrange(lh_plots$ILFL, rh_plots$ILFR, ncol=1, nrow = 2)
  ILF_final <- annotate_figure(ILF_plots, top = text_grob("Inferior Longitudinal Fasciculus", 
                                                          color = "black", size = 20))
  
  pARC_plots <- ggarrange(lh_plots$pARCL, rh_plots$pARCR, ncol=1, nrow = 2)
  pARC_final <- annotate_figure(pARC_plots, top = text_grob("Posterior Arcuate Fasciculus", 
                                                            color = "black", size = 20))
  
  SLF_plots <- ggarrange(lh_plots$SLFL, rh_plots$SLFR, ncol=1, nrow = 2)
  SLF_final <- annotate_figure(SLF_plots, top = text_grob("Superior Longitudinal Fasciculus", 
                                                          color = "black", size = 20))
  
  UNC_plots <- ggarrange(lh_plots$UNCL, rh_plots$UNCR, ncol=1, nrow = 2)
  UNC_final <- annotate_figure(UNC_plots, top = text_grob("Uncinate Fasciculus", 
                                                          color = "black", size = 20))
  
  VOF_plots <- ggarrange(lh_plots$VOFL, rh_plots$VOFR, ncol=1, nrow = 2)
  VOF_final <- annotate_figure(VOF_plots, top = text_grob("Vertical Occipital Fasciculus", 
                                                          color = "black", size = 20))
  
  all_plots <- ggarrange(COrb_final, CAntFr_final, CSupFr_final, CMot_final, CSupPar_final, CPostPar_final, CTemp_final, COcc_final, ARC_final, CST_final, IFO_final, ILF_final, pARC_final, SLF_final, UNC_final, VOF_final, ncol = 4, nrow = 4)
  
  # stitch it all together
  legend <- get_legend(lh_plots$VOFL + theme(legend.position = "bottom",
                                             legend.key.height = unit(1.5, 'cm'),
                                             legend.key.width = unit(3.5, 'cm'),
                                             legend.margin=margin(-1,0,-1,0),
                                             legend.text = element_text(size=28),
                                             legend.title = element_blank()))
  
  
  legend <- ggplot() + 
    geom_brain(data = VOFL_deveffect_HCPD, atlas= glasser, 
               mapping=aes(fill=mean_age_effect), 
               show.legend=TRUE, 
               hemi = 'right',
               position = position_brain("right lateral")) +
    paletteer::scale_fill_paletteer_c("grDevices::RdYlBu", direction = -1 , limits = c(15, 23), oob = squish, na.value = "white") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.key.height = unit(0.5, 'cm'),
      legend.key.width = unit(1.5, 'cm'),
      legend.margin=margin(0,0,0,0),
      legend.text = element_text(size = 28),
      legend.title = element_blank()
    ) 
  legend <- get_legend(legend)
  legend_title <- textGrob(expression(paste("Age of Maturation")), 
                           gp=gpar(col="black", fontsize=28), hjust = 0.5)
  
  legend_with_title <- plot_grid(legend, legend_title, ncol = 1, rel_heights = c(0.8, 0.2))
  
  final_plot <- plot_grid(all_plots, legend_with_title, ncol = 1, rel_heights = c(1, 0.12))
  return(final_plot)
}

# plot S-A rank of tract endpoints on cortical surface
plot_SA_surface <- function(bundle_name, df) {
  if(grepl("L$", bundle_name)) {
    hemi = "left"
    cortical_pos1 <- "left lateral" 
    cortical_pos2 <- "left medial" 
  } else{
    hemi = "right"
    cortical_pos1 <- "right lateral" 
    cortical_pos2 <- "right medial" 
  }
  plot_lateral <- ggplot() + 
    geom_brain(data = df, atlas= glasser, 
               mapping=aes(fill=SA.axis_rank), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos1)) +
    scale_fill_viridis_c(option = "magma", na.value = "white", limits = c(1, 360), direction = -1 ) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin =  margin(0, -1, -5, 0),
          plot.title = element_blank()) 
  
  plot_medial <- ggplot() + 
    geom_brain(data = df, atlas= glasser, 
               mapping=aes(fill=SA.axis_rank), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos2)) +
    scale_fill_viridis_c(option = "magma", na.value = "white", limits = c(1, 360), direction = -1 ) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin =  margin(0, -1, -5, 0),
          plot.title = element_blank()) 
  return(plot_grid(plot_lateral, plot_medial, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr"))
}

wrapper_plot_SA <- function(bundle_name, dataset) {
  L_deveffect <- get(paste0(bundle_name, "L_deveffect_", dataset))
  R_deveffect <- get(paste0(bundle_name, "R_deveffect_", dataset))
  
  merged_L <- merge(L_deveffect, glasser_SAaxis,  by = "regionName") %>% select(-label)
  merged_R <- merge(R_deveffect, glasser_SAaxis,  by = "regionName") %>% select(-label)
  
  merged_L_mean <- merged_L %>% group_by(end) %>% 
    mutate(mean_SA = mean(SA.axis_rank, na.rm = TRUE)) %>% 
    ungroup() %>%
    mutate(SA.axis_rank = ifelse(is.na(mean_age_effect), NA, mean_SA)) %>%  # Replace SA.axis_rank where mean_age_effect is NA
    select(-mean_SA)
  
  merged_R_mean <- merged_R %>% group_by(end) %>% 
    mutate(mean_SA = mean(SA.axis_rank, na.rm = TRUE)) %>% 
    ungroup() %>%
    mutate(SA.axis_rank = ifelse(is.na(mean_age_effect), NA, mean_SA)) %>%  # Replace SA.axis_rank where mean_age_effect is NA
    select(-mean_SA)
  
  cortical_pos1 <- "left lateral" 
  plot_legend <- ggplot() + 
    geom_brain(data = merged_L_mean, atlas= glasser, 
               mapping=aes(fill=SA.axis_rank), 
               show.legend=TRUE, 
               hemi = "left",
               position = position_brain(cortical_pos1)) +
    scale_fill_viridis_c(option = "magma", na.value = "white", limits = c(1, 360), direction = -1) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.key.width = unit(2, 'cm'), legend.key.height = unit(0.5, 'cm'),
          plot.margin = margin(0, -1, -5, 0),
          
          plot.title = element_blank()) 
  
  legend <- get_legend(plot_legend)
  
  map_L_mean <- plot_SA_surface("L", merged_L_mean)
  map_R_mean <- plot_SA_surface("R", merged_R_mean)
  
  map_merged <- ggarrange(map_L_mean, map_R_mean, ncol = 2)
  return(list(map_merged, legend))
}


arrange_SA_plots <- function(plots) {
  # arrange maps together  
  COrb_final <- annotate_figure(plots$COrb[[1]], top = text_grob("Callosum Orbital", 
                                                                 color = "black", size = 20))
  
  CAntFr_final <- annotate_figure(plots$CAntFr[[1]], top = text_grob("Callosum Anterior Frontal", 
                                                                     color = "black", size = 20))
  
  CSupFr_final <- annotate_figure(plots$CSupFr[[1]], top = text_grob("Callosum Superior Frontal", 
                                                                     color = "black", size = 20))
  
  CMot_final <- annotate_figure(plots$CMot[[1]], top = text_grob("Callosum Motor", 
                                                                 color = "black", size = 20))
  
  CSupPar_final <- annotate_figure(plots$CSupPar[[1]], top = text_grob("Callosum Superior Parietal", 
                                                                       color = "black", size = 20))
  
  CPostPar_final <- annotate_figure(plots$CPostPar[[1]], top = text_grob("Callosum Posterior Parietal", 
                                                                         color = "black", size = 20))
  
  CTemp_final <- annotate_figure(plots$CTemp[[1]], top = text_grob("Callosum Temporal", 
                                                                   color = "black", size = 20))
  
  COcc_final <- annotate_figure(plots$COcc[[1]], top = text_grob("Callosum Occipital", 
                                                                 color = "black", size = 20))
  
  ARC_final <- annotate_figure(plots$ARC[[1]], top = text_grob("Arcuate Fasciculus", 
                                                               color = "black", size = 20))
  
  CST_final <- annotate_figure(plots$CST[[1]], top = text_grob("Corticospinal Tract", 
                                                               color = "black", size = 20))
  
  IFO_final <- annotate_figure(plots$IFO[[1]], top = text_grob("Inferior Fronto-occipital Fasciculus", 
                                                               color = "black", size = 20))
  
  ILF_final <- annotate_figure(plots$ILF[[1]], top = text_grob("Inferior Longitudinal Fasciculus", 
                                                               color = "black", size = 20))
  
  pARC_final <- annotate_figure(plots$pARC[[1]], top = text_grob("Posterior Arcuate Fasciculus", 
                                                                 color = "black", size = 20))
  
  SLF_final <- annotate_figure(plots$SLF[[1]], top = text_grob("Superior Longitudinal Fasciculus", 
                                                               color = "black", size = 20))
  
  UNC_final <- annotate_figure(plots$UNC[[1]], top = text_grob("Uncinate Fasciculus", 
                                                               color = "black", size = 20))
  
  VOF_final <- annotate_figure(plots$VOF[[1]], top = text_grob("Vertical Occipital Fasciculus", 
                                                               color = "black", size = 20))
  
  all_plots <- ggarrange(COrb_final, CAntFr_final, CSupFr_final, CMot_final, CSupPar_final, CPostPar_final, CTemp_final, COcc_final, ARC_final, CST_final, IFO_final, ILF_final, pARC_final, SLF_final, UNC_final, VOF_final, ncol = 4, nrow = 4)
  
  # stitch it all together
  legend <- plots$VOF[[2]]
  legend_title <- textGrob(expression(paste("Mean Sensorimotor-Association Axis Rank")), 
                           gp=gpar(col="black", fontsize=20), hjust = 0.5)
  legend_with_title <- plot_grid(legend, legend_title, ncol = 1, rel_heights = c(0.8, 0.2))
  final_plot <- plot_grid(all_plots, legend_with_title, ncol = 1, rel_heights = c(1, 0.12))
  return(final_plot)
}
