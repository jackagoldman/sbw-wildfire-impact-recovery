# Propensity score analysis 


# load required packages
library(tidyverse)
library(ggplot2)
library(sf)
library(viridis)
library(cowplot)
library(MatchIt)
library(sensemakr)
library(MuMIn)
source("src/best_model_functions.R")
source("src/utils.R")


# Load your data 
source("/home/goldma34/sbw-wildfire-impact-recovery/src/load_data.R")  # Replace with actual script that loads hist_gt90_1

# Part 1: SEVERITY ==============================

## Propensity score matching for full severity Data ---------------------------

# Full matching on a probit PS
m.out2 <- matchit(history ~  host_pct +isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri ,
                  data = history_gt90, 
                  method = "nearest",
                  distance = "glm",
                  link = "probit",
                  m.order = "random",
                  ### set a maximum distance for the matches
                  caliper = 0.25)



m.data <- match_data(m.out2)


m.data %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",
                             history ==0 ~"Non-Defoliated"))


#save matched data
write.csv(m.data, "/home/goldma34/sbw-wildfire-impact-recovery/data/on_sev_match_data.csv")


# Convert the data frame to an sf object
m.data.sf <- st_as_sf(m.data, coords = c("x", "y"), crs = 4326)


## Assess balance for severity ======================
v <- data.frame(old = c("host_pct", "isi_90", "dc_90", 
"dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri"),
new = c("Host Species Percentage", "Intial Spread Index", "Drought Code", 
        "Duff Moisture Code", "Fine Fuel Moisture Code", "Build Up Index", "Fire Weather Index", "Mean Terrain Ruggedness Index"))

cobalt::love.plot(m.out2, stats = c("mean.diffs"), 
                  threshold = c(m = .25), 
                  binary = "std",
                  abs = TRUE,
                  var.order = "unadjusted",
                  var.names = v,
                  limits = c(0, 1),
                  grid = FALSE,
                  wrap = 20,
                  sample.names = c("Unmatched", "Matched"),
                  position = "top",
                  shapes = c("circle", "triangle"),
                  colors = c("#FF8C00A0", "#8B0000A0"))


host_plot <- cobalt::bal.plot(m.out2, 
                              var.name = "host_pct", 
                              which = "both",
                              colors = c("#FF8C00A0", "#8B0000A0"),
                              sample.names = c("Unmatched", "Matched")) +scale_fill_manual(
                                values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                                labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Host Species Percentage") + ggtitle(NULL) +theme_bw() + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

isi_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "isi_90", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Initial Spread Index") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

dc_plot <- cobalt::bal.plot(m.out2, 
                            var.name = "dc_90", 
                            which = "both", 
                            colors = c("#FF8C00A0", "#8B0000A0"),
                            sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                              values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                              labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Drought Code") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

dmc_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "dmc_90", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Duff Moisture Code") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

ffmc_plot <- cobalt::bal.plot(m.out2, 
                              var.name = "ffmc_90", 
                              which = "both", 
                              colors = c("#FF8C00A0", "#8B0000A0"),
                              sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                                values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                                labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Fine Fuel Moisture Code") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

bui_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "bui_90", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Build Up Index") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

fwi_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "fwi_90", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Fire Weather Index") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

tri_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "mean_tri", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Fire Weather Index") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# Arrange the plots in a 2-column layout
sev_combined <- plot_grid(host_plot, isi_plot, dc_plot, dmc_plot, ffmc_plot, bui_plot,fwi_plot, ncol = 1, align = 'v', rel_heights = c(1, 1, 1, 1, 1, 1))

# save sev_combined
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/balance/fig_sev_all_balance.png",
  plot = sev_combined,
  width = 14, height = 5 * ceiling(length(sev_combined)/2), # adjust height based on number of plots
  dpi = 300
)


## Severity average treatment effect estimate =====================

# fit model
fit_sev <- lm(rbr_w_offset  ~ history + host_pct+ isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,  data = m.data,weights = weights)

# avg_comparisons for treatment effect
fit_att_sev <- as.data.frame(marginaleffects::avg_comparisons(fit_sev,
                                                          variables = "history",
                                                          vcov = ~subclass,
                                                          newdata = subset(history == 1))) %>%  mutate(Fires= "all")

# results
fit_att_sev

#rsquared
r.squaredGLMM(fit_sev)

# Create a severity treatment effect data frame with better formatted column names
sev_treatment_effects <- fit_att_sev %>%
  rename(
    Treatment_Effect = estimate,
    SE = std.error,
    CI_Lower = conf.low,
    CI_Upper = conf.high,
    P_Value = p.value
  ) %>%
  mutate(
    Model = "Severity",
    .before = 1
  ) %>%
  select(Model, Treatment_Effect, SE, CI_Lower, CI_Upper, P_Value, Fires)

# save rds
saveRDS(sev_treatment_effects, "/home/goldma34/sbw-wildfire-impact-recovery/results/all_fires/treatment_effects_sev.rds")


# Save this for later when we'll combine with recovery effects

# effects plot
p.sev_all <- marginaleffects::plot_predictions(fit_sev, condition = c("history"), draw = FALSE) 

p.sev_all<- p.sev_all %>% 
  mutate(Fires = "All")  %>% 
  mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))

# Get max y value for positioning
sev_y_max <- max(p.sev_all$conf.high) * 1.1

# Create base plot
plot.p.sev_all <- ggplot(p.sev_all, aes(x = history, y = estimate, color = history)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Burn Severity", x = "Defoliation History", color = "History") +
  scale_x_discrete(labels = c("Defoliated", "Non-Defoliated")) +
  scale_color_manual(values = c("Defoliated" = "#8B0000A0", "Non-Defoliated" = "#FF8C00A0")) +
  theme_bw()+
  theme(legend.position = "none")

# Add significance annotation
plot.p.sev_all <- add_significance_annotation(plot.p.sev_all, fit_att_sev, y_position = sev_y_max)

# save rds
saveRDS(plot.p.sev_all, "/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/fig_sev_treat_effect_all.rds")

# Save plot
ggsave(plot = plot.p.sev_all, filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/fig_severity_treat_effect_all.png", width = 6, height = 4, dpi = 300)

## Sensitivity analysis for severity #####

# benchmark covariates
covariates <- c("host_pct", "isi_90", "dc_90", "dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri")

# run sensitivity analysis
sev.sensitivity <- sensemakr(model = fit_sev, 
                                treatment = "history",
                                benchmark_covariates = covariates,
                                alpha = 0.05, 
                                reduce = TRUE)

#outcome
summary(sev.sensitivity)


#The results suggest that the 'history' variable has a moderate coefficient estimate with a relatively 
#low standard error, leading to a high t-value and indicating that the effect is statistically significant. 
#The sensitivity statistics imply that the observed effect is reasonably robust and would 
#require relatively strong unobserved confounders to nullify or make the effect statistically insignificant.

#plot bias contour of point estimate
plot(sev.sensitivity)

# Convert sensitivity plot to ggplot and save
ggplot_sev_sensitivity <- plot(sev.sensitivity, type = "contour") + theme_bw()
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/sensitivity/fig_sev_sensitivity.png",
  plot = ggplot_sev_sensitivity,
  width = 8, height = 6, dpi = 300
)

# plot extreme scenario
plot(sev.sensitivity, type = "extreme")

#####################################
# Part 2: RECOVERY ==============================
######################################

## Propensity score matching for full recovery Data --------------------------

# Full matching on a probit PS
m.out.rec <- matchit(history ~  host_pct + rbr_w_offset + mean_temperature + sum_precipitation_mm + mean_tri,
                     data = history_gt90,
                     method = "nearest",
                     distance = "glm",
                     link = "probit",
                     m.order = "random",
                     ### set a maximum distance for the matches
                     caliper = 0.25)



m.data.rec <- match_data(m.out.rec)


m.data.rec %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",
                             history ==0 ~"Non-Defoliated"))


write.csv(m.data.rec, "/home/goldma34/sbw-wildfire-impact-recovery/data/on_rec_match_data.csv")

# Convert the data frame to an sf object
m.data.rec.sf <- st_as_sf(m.data.rec, coords = c("x", "y"), crs = 4326)


## Assess balance for recovery ==============
v.r <- data.frame(old = c("host_pct", "rbr_w_offset", "mean_temperature", "sum_precipitation_mm", "mean_tri"),
                new = c("Host Species Percentage", "burn severity", "post-fire mean temperature", 
                        "post-fire total precipitation (mm)", "mean terrain ruggedness index (m)"))

cobalt::love.plot(m.out.rec, stats = c("mean.diffs"), 
                  threshold = c(m = .25), 
                  binary = "std",
                  abs = TRUE,
                  var.order = "unadjusted",
                  var.names = v.r,
                  limits = c(0, 1),
                  grid = FALSE,
                  wrap = 20,
                  sample.names = c("Unmatched", "Matched"),
                  position = "top",
                  shapes = c("circle", "triangle"),
                  colors = c("#FF8C00A0", "#8B0000A0"))


host_plot <- cobalt::bal.plot(m.out.rec, 
                              var.name = "host_pct", 
                              which = "both",
                              colors = c("#FF8C00A0", "#8B0000A0"),
                              sample.names = c("Unmatched", "Matched")) +scale_fill_manual(
                                values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                                labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Host Species Percentage") + ggtitle(NULL) +theme_bw() + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

rbr_plot <- cobalt::bal.plot(m.out.rec, 
                             var.name = "rbr_w_offset", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Burn Severity (RBR)") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

temp_plot <- cobalt::bal.plot(m.out.rec, 
                            var.name = "mean_temperature", 
                            which = "both", 
                            colors = c("#FF8C00A0", "#8B0000A0"),
                            sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                              values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                              labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Mean Temperature (°C)") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

ppt_plot <- cobalt::bal.plot(m.out.rec, 
                             var.name = "sum_precipitation_mm", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Total Precipitation (mm)") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

tri_plot <- cobalt::bal.plot(m.out.rec, 
                             var.name = "mean_tri", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Mean Terrain Ruggedness Index (m)") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))



# Arrange the plots in a 2-column layout
rec_combined <- plot_grid(host_plot, rbr_plot, temp_plot, ppt_plot, tri_plot, ncol = 1, align = 'v', rel_heights = c(1, 1, 1, 1, 1, 1))

# save sev_combined
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/balance/fig_rec_all_balance.png",
  plot = rec_combined,
  width = 14, height = 5 * ceiling(length(rec_combined)/2), # adjust height based on number of plots
  dpi = 300
)


## Recovery average treatment estimate ========
m.data.rec_test  <- m.data.rec %>% 
  mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))

fit_rec_baseline <- lm(recovery  ~ history + host_pct+rbr_w_offset + mean_temperature + sum_precipitation_mm + mean_tri,  data = m.data.rec,weights = weights)

fit_att_rec <- as.data.frame(marginaleffects::avg_comparisons(fit_rec_baseline,
                                                              variables = "history",
                                                              vcov = ~subclass,
                                                              newdata = subset(history == 1))) %>% 
  mutate(Fires= "all")

# summary
fit_att_rec

#r squared
r.squaredGLMM(fit_rec_baseline)

# Create a recovery treatment effect data frame with better formatted column names
rec_treatment_effects <- fit_att_rec %>%
  rename(
    Treatment_Effect = estimate,
    SE = std.error,
    CI_Lower = conf.low,
    CI_Upper = conf.high,
    P_Value = p.value
  ) %>%
  mutate(
    Model = "Recovery",
    .before = 1
  ) %>%
  select(Model, Treatment_Effect, SE, CI_Lower, CI_Upper, P_Value, Fires)

#save rds
saveRDS(rec_treatment_effects, "/home/goldma34/sbw-wildfire-impact-recovery/results/all_fires/treatment_effects_rec.rds")

# Combine both treatment effects
combined_effects <- bind_rows(sev_treatment_effects, rec_treatment_effects)

# Export combined treatment effects to CSV
write.csv(
  combined_effects,
  "/home/goldma34/sbw-wildfire-impact-recovery/results/all_fires/treatment_effects_all_fires.csv",
  row.names = FALSE
)

# Print a message confirming export
cat("Treatment effects exported to results/treatment_effects.csv\n")


# effects plot
p.rec_all <- marginaleffects::plot_predictions(fit_rec_baseline, condition = c("history"), draw = FALSE) 

p.rec_all <- p.rec_all %>% 
  mutate(Fires = "All")  %>% 
  mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))

# Get max y value for positioning
rec_y_max <- max(p.rec_all$conf.high) * 1.1

# Create base plot
plot.p.rec_all <- ggplot(p.rec_all, aes(x = history, y = estimate, color = history)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Recovery Magnitude (%)", x = "Defoliation History", color = "History") +
  scale_x_discrete(labels = c("Defoliated", "Non-Defoliated")) +
  scale_color_manual(values = c("Defoliated" = "#8B0000A0", "Non-Defoliated" = "#FF8C00A0")) +
  theme_bw()+
  theme(legend.position = "none")

# Add significance annotation
plot.p.rec_all <- add_significance_annotation(plot.p.rec_all, fit_att_rec, y_position = rec_y_max)

#save rds
saveRDS(plot.p.rec_all, "/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/fig_rec_treat_effect_all.rds")


# Save plot
ggsave(plot = plot.p.rec_all, filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/fig_rec_treat_effect_all.png", width = 6, height = 4, dpi = 300)


##  sensitivity analysis for recovery  -----------

# benchmark covariates
covariates_rec <- c("host_pct", "rbr_w_offset", "mean_temperature", "sum_precipitation_mm", "mean_tri")

# run sensitivity analysis
rec.sensitivity <- sensemakr(model = fit_rec_baseline, 
                             treatment = "history",
                             benchmark_covariates = covariates_rec,
                             alpha = 0.05, 
                             reduce = TRUE)



#outcome
summary(rec.sensitivity)

#The results suggest that the 'history' variable has a moderate coefficient estimate with a relatively 
#low standard error, leading to a high t-value and indicating that the effect is statistically significant. 
#The sensitivity statistics imply that the observed effect is reasonably robust and would require relatively 
#strong unobserved confounders to nullify or make the effect statistically insignificant.
#Overall, these results indicate that 'history' has a significant and robust effect on the outcome.

#plot bias contour of point estimate
plot(rec.sensitivity)

# Convert sensitivity plot to ggplot and save
ggplot_rec_sensitivity <- plot(rec.sensitivity, type = "contour") + theme_bw()
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/sensitivity/fig_rec_sensitivity.png",
  plot = ggplot_rec_sensitivity,
  width = 8, height = 6, dpi = 300
)
# plot extreme scenario
plot(rec.sensitivity, type = "extreme")


# Export R-square values to csv
# Create dataframe with R-squared values for both models
r_squared_table <- data.frame(
  Model = c("Severity", "Recovery"),
  R2_marginal = c(r.squaredGLMM(fit_sev)[1, "R2m"], r.squaredGLMM(fit_rec_baseline)[1, "R2m"]),
  R2_conditional = c(r.squaredGLMM(fit_sev)[1, "R2c"], r.squaredGLMM(fit_rec_baseline)[1, "R2c"])
)

# Export R-squared table to CSV
write.csv(
  r_squared_table,
  "/home/goldma34/sbw-wildfire-impact-recovery/results/all_fires/r_squared_values_all_fires.csv",
  row.names = FALSE
)

# Print a message confirming export
cat("R-squared values exported to results/r_squared_values.csv\n")