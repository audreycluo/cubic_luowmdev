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

# function for plotting tractprofiles mean and sd or se for a given scalar
plot_mean_var <- function(tract, scalar, variance_measure, color1, color2, color3, ylim1 = 0, ylim2 = 1, ylim3 = NULL, ylim4 = NULL, clipEnds) {
  mean_scalar <- paste0("mean_", scalar)
  summary_HCPD <- get(paste0("summary_HCPD_", scalar)) %>% filter(nodeID > (clipEnds-1) & nodeID < (99-clipEnds+1))
  summary_HBN <- get(paste0("summary_HBN_", scalar)) %>% filter(nodeID > (clipEnds-1) & nodeID < (99-clipEnds+1))
  summary_PNC <- get(paste0("summary_PNC_", scalar)) %>% filter(nodeID > (clipEnds-1) & nodeID < (99-clipEnds+1))
  
  df_profiles <- rbind(summary_HCPD, summary_HBN, summary_PNC)
  df <- df_profiles %>% filter(tract_label == tract)
  df$Dataset <- factor(df$Dataset, levels = c("PNC", "HCPD", "HBN"))
  # plot
  if(str_detect(tract, "Callosum")) {
    plot <- ggplot(data = df, aes(x = nodeID, y = get(mean_scalar))) +
      geom_ribbon(data = df, aes_string(x = "nodeID", ymin = paste0("ymin_", variance_measure), ymax = paste0("ymax_", variance_measure), fill = "Dataset"), alpha = .3) + 
      geom_line(data = df, aes(x = nodeID, y = get(mean_scalar), color = Dataset),size = 1) +
      scale_fill_manual(values = c("HCPD" = color1, "HBN" = color2, "PNC" = color3)) +
      scale_color_manual(values = c("HCPD" = color1, "HBN" = color2, "PNC" = color3)) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=20, hjust=0.5),
            plot.margin = unit(c(1, 1, 0.2, 0.2), "cm")) + labs(title = tract) + ylim(ylim3, ylim4) + 
      scale_x_continuous(
        breaks = c(10, 50, 90),                 
        labels = c("Left", "Deep", "Right")  
      )
    
  } else {
    if(tract %in% c("Posterior Arcuate", "Vertical Occipital")) {
      x_scale <- scale_x_continuous(breaks = c(15, 52, 86), labels = c("Superior", "Deep", "Inferior"))
    } else if(tract %in% c("Corticospinal")) {
      x_scale <- scale_x_continuous(breaks = c(15, 75), labels = c("Superior", "Deep"))
    } else {
      x_scale <- scale_x_continuous(breaks = c(15, 49, 86), labels = c("Anterior", "Deep", "Posterior")) 
    }
    plot <- ggplot(data = df, aes(x = nodeID, y = get(mean_scalar), linetype = hemi)) +
      geom_ribbon(data = df, aes_string(x = "nodeID", ymin = paste0("ymin_", variance_measure), ymax = paste0("ymax_", variance_measure), fill = "Dataset"), alpha = .3) + 
      geom_line(data = df, aes(x = nodeID, y = get(mean_scalar), color = Dataset),size = 1) +
      guides(linetype = guide_legend("Hemi", override.aes = list(fill = "white"))) + 
      scale_fill_manual(values = c("HCPD" = color1, "HBN" = color2, "PNC" = color3)) +
      scale_color_manual(values = c("HCPD" = color1, "HBN" = color2, "PNC" = color3)) +
      x_scale + 
      theme(legend.position = "bottom",
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            legend.box = "vertical",
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=20, hjust=0.5),
            plot.margin = unit(c(1, 1, 0.2, 0.2), "cm")) + labs(title = tract) + ylim(if (str_detect(tract, "Corticospinal") | str_detect(tract, "Inferior Fronto")) 
              ylim3 else ylim1, if (str_detect(tract, "Corticospinal") | str_detect(tract, "Inferior Fronto")) ylim4 else ylim2)
  }
  return(plot)
}

# plot coefficient of variance along tracts
plot_cv <- function(tract, scalar, variance_measure, color1, color2, color3, ylim1 = 0, ylim2 = 1, ylim3 = NULL, ylim4 = NULL, clipEnds) {
  mean_scalar <- paste0("mean_", scalar)
  summary_HCPD <- get(paste0("summary_HCPD_", scalar)) %>% filter(nodeID > (clipEnds-1) & nodeID < (99-clipEnds+1))
  summary_HBN <- get(paste0("summary_HBN_", scalar)) %>% filter(nodeID > (clipEnds-1) & nodeID < (99-clipEnds+1))
  summary_PNC <- get(paste0("summary_PNC_", scalar)) %>% filter(nodeID > (clipEnds-1) & nodeID < (99-clipEnds+1))
  
  df_profiles <- rbind(summary_HCPD, summary_HBN, summary_PNC)
  df <- df_profiles %>% filter(tract_label == tract)
  df$Dataset <- factor(df$Dataset, levels = c("PNC", "HCPD", "HBN"))
  
  # plot
  if(str_detect(tract, "Callosum")) {
    plot <- ggplot(data = df, aes(x = nodeID, y = get(variance_measure))) +
      geom_line(data = df, aes(x = nodeID, y = get(variance_measure), color = Dataset),size = 1) +
      scale_fill_manual(values = c("HCPD" = color1, "HBN" = color2, "PNC" = color3)) +
      scale_color_manual(values = c("HCPD" = color1, "HBN" = color2, "PNC" = color3)) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=20, hjust=0.5),
            plot.margin = unit(c(1, 1, 0.2, 0.2), "cm")) + labs(title = tract) + ylim(ylim3, ylim4) +
      scale_x_continuous(
        breaks = c(10, 50, 90),                 
        labels = c("Left", "Deep", "Right")  
      )
    
  } else {
    if(tract %in% c("Posterior Arcuate", "Vertical Occipital")) {
      x_scale <- scale_x_continuous(breaks = c(15, 52, 86), labels = c("Superior", "Deep", "Inferior"))
    } else if(tract %in% c("Corticospinal")) {
      x_scale <- scale_x_continuous(breaks = c(15, 75), labels = c("Superior", "Deep"))
    } else {
      x_scale <- scale_x_continuous(breaks = c(15, 49, 86), labels = c("Anterior", "Deep", "Posterior")) 
    }
    
    plot <- ggplot(data = df, aes(x = nodeID, y = get(variance_measure), linetype = hemi)) +
      geom_line(data = df, aes(x = nodeID, y = get(variance_measure), color = Dataset),size = 1) +
      guides(linetype = guide_legend("Hemi", override.aes = list(fill = "white"))) + 
      scale_fill_manual(values = c("HCPD" = color1, "HBN" = color2, "PNC" = color3)) +
      scale_color_manual(values = c("HCPD" = color1, "HBN" = color2, "PNC" = color3)) +
      x_scale + 
      theme(legend.position = "bottom",
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            legend.box = "vertical",
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=20, hjust=0.5),
            plot.margin = unit(c(1, 1, 0.2, 0.2), "cm")) + labs(title = tract) + ylim(if (str_detect(tract, "Corticospinal") | str_detect(tract, "Inferior Fronto")) 
              ylim3 else ylim1, if (str_detect(tract, "Corticospinal") | str_detect(tract, "Inferior Fronto")) ylim4 else ylim2)
  }
  return(plot)
}


# arrange callosum 
arrange_callosum_var_plots <- function(list_figures) {
  callosum_plot_final <- ggarrange(list_figures$`Callosum Orbital`, list_figures$`Callosum Anterior Frontal`, list_figures$`Callosum Superior Frontal`, list_figures$`Callosum Motor`, 
                                   list_figures$`Callosum Superior Parietal`, list_figures$`Callosum Posterior Parietal`, list_figures$`Callosum Temporal`, list_figures$`Callosum Occipital`, 
                                   nrow = 2, ncol = 4, common.legend = T, legend = "bottom") +
    bgcolor("white") + theme_void()
  return(callosum_plot_final)
}


# arrange association bundles + CST 
arrange_all_association_var_plots <- function(list_figures) {
  association_plots_final <- ggarrange(list_figures$`Arcuate`, list_figures$`Inferior Fronto-occipital`, list_figures$`Inferior Longitudinal`, list_figures$`Superior Longitudinal`, 
                                       list_figures$`Posterior Arcuate`, list_figures$`Uncinate`, list_figures$`Vertical Occipital`, list_figures$`Corticospinal`, 
                                       nrow = 2, ncol = 4, common.legend = T, legend = "bottom") +
    bgcolor("white") + theme_void()
  return(association_plots_final) 
}

