library(ggpubr)
library(ggplot2)
library(dplyr)
library(patchwork)

#utils
source("src/utils.R")


# area maps
matched_map_severity<-readRDS("./plots/maps/study_area_severity_map.rds")
matched_map_recovery<-readRDS("./plots/maps/study_area_recovery_map.rds")

# spatial residual plots
spat.resid.sev_plot<-readRDS("./plots/maps/residuals_map_severity.rds")
spat.resid.rec_plot<-readRDS("./plots/maps/residuals_map_recovery.rds")

# treatment effects plots
treat.plots_sev <-readRDS("./plots/treat_effects/fig_sev_treat_effect_all.rds")
treat.plots_rec <-readRDS("./plots/treat_effects/fig_rec_treat_effect_all.rds")

#Panels of all main plots from the general psm analysis

# recovery
publication_panel_all_fires(matched_map_recovery, treat.plots_rec, spat.resid.rec_plot, "recovery") 

#severity
publication_panel_all_fires(matched_map_severity, treat.plots_sev, spat.resid.sev_plot, "severity")


# set out the pplots
# for the tags - make sure the annotation is done for each plot first
# because theme elements are being set at the end in the tight layout - 
# the tags need to be set individually first
# First set tags individually on each plot

