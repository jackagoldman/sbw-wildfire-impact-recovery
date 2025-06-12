### Subgroup spatial autocorrelation ###

# load required packages
library(lme4)
library(MuMIn)
library(lmerTest)
library(nlme)
library(sp)
library(gstat)
library(ggplot2)
library(car)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(cowplot)
library(spdep)
library(adespatial)


#check and set working directory
# Function to set appropriate path based on working directory
set_appropriate_path <- function() {
  current_wd <- getwd()
  cat("Current working directory:", current_wd, "\n")
  
  # Check if working directory contains '/goldma34/'
  if (grepl("/goldma34/", current_wd)) {
    base_path <- "/home/goldma34/sbw-wildfire-impact-recovery/"
    cat("Using server path:", base_path, "\n")
  } else {
    # Use current working directory as base
    base_path <- file.path(getwd())
    cat("Using local path:", base_path, "\n")
  }
  
  return(base_path)
}

# Set the base path
base_path <- set_appropriate_path()


#utils
source(file.path(base_path, "/src/utils.R"))


# Load your data 
source(file.path(base_path, "/src/load_data.R"))

#load
source(file.path(base_path, "/src/best_model_functions.R"))

# Load Ontario map data
ontario <- ne_states(country = "Canada", returnclass = "sf") %>%
  filter(name == "Ontario")


######################################
# SPATIAL DISTRIBUTION MODELS
######################################

# defol dataset - models are for defolaited only fires, make two dataframes for attach residuals and plotting
defol_only_sev <- subset(history_gt90, history == 1)
defol_only_rec <- subset(history_gt90, history == 1)


# window of opportunity subset 1-4
defol_only_sev_1 <- subset(defol_only_sev, window_opp == 1)
defol_only_sev_2 <- subset(defol_only_sev, window_opp == 2)
defol_only_sev_3 <- subset(defol_only_sev, window_opp == 3)
defol_only_sev_4 <- subset(defol_only_sev, window_opp == 4)


########################################
# Part 1: SEVERITY ======================
########################################

# sample sizes for different windows of opportunity
cat("Sample size for Window of Opportunity 1:", nrow(defol_only_sev_1), "\n")
cat("Sample size for Window of Opportunity 2:", nrow(defol_only_sev_2), "\n")
cat("Sample size for Window of Opportunity 3:", nrow(defol_only_sev_3), "\n")
cat("Sample size for Window of Opportunity 4:", nrow(defol_only_sev_4), "\n")


### Window of Opportunity 1 =========================
# Define the model formula
formula.sev <- rbr_w_offset ~host_pct + Cumulative_Years_Defol + isi_90 + 
  dc_90 + dmc_90 + ffmc_90 + bui_90+ fwi_90 + mean_tri+  x + y 


# Fit the model using nlme with spatial correlation, where correlation decreases with distance
model_gls_sev_wo1 <- gls(formula.sev,
                 correlation = corGaus(form = ~ x + y), 
                 data = defol_only_sev_1,
                 control = list(singular.ok = TRUE))

# Identify variables with high VIF
high_vif_results <- identify_high_vif(model_gls_sev_wo1, threshold = 10)

# Get the list of variables with high VIF
high_vif_vars <- high_vif_results$high_vif_vars

# Print the variables with high VIF
print(high_vif_vars)

# remove high VIF variables from the model
formula.wo1 <- update(formula.sev, . ~ . - isi_90 - dc_90 - dmc_90  - bui_90 - fwi_90)
# Fit the model again without high VIF variables
model_gls_sev_wo1 <- gls(formula.wo1,
                         correlation = corGaus(form = ~ x + y), 
                         data = defol_only_sev_1,
                         control = list(singular.ok = TRUE))    

summary(model_gls_sev_wo1)

# get list of variables in model
formula_vars_wo1 <- all.vars(formula.wo1)

# identify high VIF variables in the new model
high_vif_results_wo1 <- identify_high_vif(model_gls_sev_wo1, threshold = 10)
# Get the list of variables with high VIF in the new model
high_vif_vars_wo1 <- high_vif_results_wo1$high_vif_vars
# Print the variables with high VIF in the new model
print(high_vif_vars_wo1)

# Residuals vs. Fitted Values Plot

# With these lines:
resid_data <- data.frame(
  fitted = fitted(model_gls_sev_wo1),
  residuals = resid(model_gls_sev_wo1)
)

resid_plot_sev_wo1 <- ggplot(resid_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid") +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw()


#save residuals plot 
ggsave(
  filename = file.path(base_path, "/plots/spatial_analysis/residuals_fitted_severity_wo1.png"),
  plot = resid_plot_severity,
  width = 12,
  height = 8,
  dpi = 300
)

# Display the model summary
mod.sum_sev_wo1 <- summary(model_gls_sev_wo1)
print(mod.sum_sev_wo1)
saveRDS(mod.sum_sev_wo1, file.path(base_path, "/results/spatial_analysis/model_summary_severity_wo1.rds"))

# r2
r2.sev_wo1 <- rsquared.gls(model_gls_sev_wo1, formula.wo1)
print(r2.sev_wo1)
saveRDS(r2.sev_wo1, file.path(base_path, "/results/spatial_analysis/r2_severity.rds"))

# extract residuals
defol_only_sev_1$residuals_gls <- residuals(model_gls_sev_wo1)

colnames(defol_only_sev_1) <- make.names(colnames(defol_only_sev_1))

# get coords
coordinates(defol_only_sev_1) <- ~ x + y

# Compute the variogram
variogram_gls_wo1 <- gstat::variogram(residuals_gls ~ 1, data = defol_only_sev_1)

# Plot the variogram
var_plot <- plot(variogram_gls_wo1, main = "Variogram of Residuals (Gaussian)")

# confirm defol_only as dataframe
df.sev <- as.data.frame(defol_only_sev_1)

# Create a spatial plot of residuals
spat.resid.sev_plot_wo1 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev, aes(x = x, y = y, color = residuals_gls)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

#save rds
saveRDS(spat.resid.sev_plot_wo1, file.path(base_path, "/plots/maps/residuals_map_severity_wo1.rds"))


# save spatial plot of residuals
ggsave(
  filename = file.path(base_path, "/plots/maps/residuals_map_severity_wo1.png"),
  plot = spat.resid.sev_plot_wo1,
  width = 12,
  height = 8,
  dpi = 300
)

# window opportunity 2 =========================

# Fit the model using nlme with spatial correlation, where correlation decreases with distance
model_gls_sev_wo2 <- gls(formula.sev,
                 correlation = corGaus(form = ~ x + y), 
                 data = defol_only_sev_2,
                 control = list(singular.ok = TRUE))

# Identify variables with high VIF
high_vif_results_2 <- identify_high_vif(model_gls_sev_wo2, threshold = 10)

# Get the list of variables with high VIF
high_vif_vars_2 <- high_vif_results_2$high_vif_vars_2

# Print the variables with high VIF
print(high_vif_vars_2)

# remove high VIF variables from the model
# keep fwi_90 and cumulative_Years_Defol
formula.wo2 <- update(formula.sev, . ~ . - isi_90 - dc_90 - dmc_90  - bui_90 - ffmc_90)
# Fit the model again without high VIF variables
model_gls_sev_wo2 <- gls(formula.wo2,
                         correlation = corGaus(form = ~ x + y), 
                         data = defol_only_sev_2,
                         control = list(singular.ok = TRUE))    

summary(model_gls_sev_wo1)

# get list of variables in model
formula_vars_wo2 <- all.vars(formula.wo2)

# identify high VIF variables in the new model
high_vif_results_wo2 <- identify_high_vif(model_gls_sev_wo2, threshold = 10)
# Get the list of variables with high VIF in the new model
high_vif_vars_wo2 <- high_vif_results_wo2$high_vif_vars
# Print the variables with high VIF in the new model
print(high_vif_vars_wo1)

# Residuals vs. Fitted Values Plot

# With these lines:
resid_data_wo2 <- data.frame(
  fitted = fitted(model_gls_sev_wo2),
  residuals = resid(model_gls_sev_wo2)
)

resid_plot_sev_wo2 <- ggplot(resid_data_wo2, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid") +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw()


# save residuals plot 
ggsave(
  filename = file.path(base_path, "/plots/spatial_analysis/residuals_fitted_severity_wo2.png"),
  plot = resid_plot_sev_wo2,
  width = 12,
  height = 8,
  dpi = 300
)

# Display the model summary
mod.sum_sev_wo2 <- summary(model_gls_sev_wo2)
print(mod.sum_sev_wo2)
saveRDS(mod.sum_sev_wo1, file.path(base_path, "/results/spatial_analysis/model_summary_severity_wo2.rds"))

# r2
r2.sev_wo2 <- rsquared.gls(model_gls_sev_wo2, formula.wo2)
print(r2.sev_wo2)
saveRDS(r2.sev_wo2, file.path(base_path, "/results/spatial_analysis/r2_severity.rds"))

# extract residuals
defol_only_sev_2$residuals_gls <- residuals(model_gls_sev_wo2)

colnames(defol_only_sev_2) <- make.names(colnames(defol_only_sev_2))

# get coords
coordinates(defol_only_sev_2) <- ~ x + y

# Compute the variogram
variogram_gls_wo2 <- gstat::variogram(residuals_gls ~ 1, data = defol_only_sev_2)

# Plot the variogram
var_plot_wo2 <- plot(variogram_gls_wo2, main = "Variogram of Residuals (Gaussian)")

# confirm defol_only as dataframe
df.sev_wo2 <- as.data.frame(defol_only_sev_2)

# Create a spatial plot of residuals
spat.resid.sev_plot_wo2 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_wo2, aes(x = x, y = y, color = residuals_gls)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

#save rds
saveRDS(spat.resid.sev_plot_wo1, file.path(base_path, "/plots/maps/residuals_map_severity_wo1.rds"))


# save spatial plot of residuals
ggsave(
  filename = file.path(base_path, "/plots/maps/residuals_map_severity_wo1.png"),
  plot = spat.resid.sev_plot_wo1,
  width = 12,
  height = 8,
  dpi = 300
)

# test MSR

# Create spatial weights
coords <- cbind(defol_only_sev_2$x, defol_only_sev_2$y)

#
nb <- try(dnearneigh(defol_only_sev_2, 0, 0.5), silent = TRUE)
if (inherits(nb, "try-error") || any(card(nb) == 0)) {
  # If some points have no neighbors, try increasing the distance
  cat("Some points have no neighbors with distance 0.5, trying 1.0...\n")
  nb <- dnearneigh(defol_only_sev_2, 0, 1.0)
  
  # Check if we still have isolated points
  if (any(card(nb) == 0)) {
    cat("Some points still have no neighbors with distance 1.0, trying 2.0...\n")
    nb <- dnearneigh(defol_only_sev_2, 0, 2.0)
    
    # If still problematic, use k-nearest neighbors instead
    if (any(card(nb) == 0)) {
      cat("Using k-nearest neighbors (k=5) instead of distance-based...\n")
      nb <- knn2nb(knearneigh(coordinates(defol_only_sev_2), k=5))
    }
  }
}


listw <- nb2listw(nb, style = "W")

# Compute Moran's Eigenvector Maps (MEMs)
mem <- scores.listw(listw, MEM.autocor = "positive")


# Assume 'resid' is a vector of residuals from a model
resid <- residuals(model_gls_sev_wo2, type = "normalized")

# Moran Spectral Randomization test
msr_test <- moran.randtest(resid, listw)

# Plot and interpret
plot(msr_test)

#A significant result indicates non-random spatial structure in the residuals.


# window opportunity 3 =========================

# Fit the model using nlme with spatial correlation, where correlation decreases with distance
model_gls_sev_wo3 <- gls(formula.sev,
                 correlation = corGaus(form = ~ x + y), 
                 data = defol_only_sev_3,
                 control = list(singular.ok = TRUE))

# Identify variables with high VIF
high_vif_results_3 <- identify_high_vif(model_gls_sev_wo3, threshold = 10)

# Get the list of variables with high VIF
high_vif_vars_3 <- high_vif_results_3$high_vif_vars

# Print the variables with high VIF
print(high_vif_vars_3)

# remove high VIF variables from the model
# keep fwi_90 and cumulative_Years_Defol
formula.wo3 <- update(formula.sev, . ~ . - isi_90 - dc_90 - dmc_90  - bui_90 - fwi_90)
# Fit the model again without high VIF variables
model_gls_sev_wo3 <- gls(formula.wo3,
                         correlation = corGaus(form = ~ x + y), 
                         data = defol_only_sev_3,
                         control = list(singular.ok = TRUE))    

summary(model_gls_sev_wo3)

# get list of variables in model
formula_vars_wo3 <- all.vars(formula.wo3)

# identify high VIF variables in the new model
high_vif_results_wo3 <- identify_high_vif(model_gls_sev_wo3, threshold = 10)
# Get the list of variables with high VIF in the new model
high_vif_vars_wo3 <- high_vif_results_wo3$high_vif_vars
# Print the variables with high VIF in the new model
print(high_vif_vars_wo1)

# Residuals vs. Fitted Values Plot

# With these lines:
resid_data_wo3 <- data.frame(
  fitted = fitted(model_gls_sev_wo3),
  residuals = resid(model_gls_sev_wo3)
)

resid_plot_sev_wo3 <- ggplot(resid_data_wo3, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid") +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw()


#save residuals plot 
ggsave(
  filename = file.path(base_path, "/plots/spatial_analysis/residuals_fitted_severity_wo3.png"),
  plot = resid_plot_sev_wo3,
  width = 12,
  height = 8,
  dpi = 300
)

# Display the model summary
mod.sum_sev_wo3 <- summary(model_gls_sev_wo3)
print(mod.sum_sev_wo3)
saveRDS(mod.sum_sev_wo3, file.path(base_path, "/results/spatial_analysis/model_summary_severity_wo3.rds"))

# r2
r2.sev_wo3 <- rsquared.gls(model_gls_sev_wo3, formula.wo3)
print(r2.sev_wo3)
saveRDS(r2.sev_wo3, file.path(base_path, "/results/spatial_analysis/r2_severity.rds"))

# extract residuals
defol_only_sev_3$residuals_gls <- residuals(model_gls_sev_wo3)

colnames(defol_only_sev_3) <- make.names(colnames(defol_only_sev_3))

# get coords
coordinates(defol_only_sev_3) <- ~ x + y

# Compute the variogram
variogram_gls_wo3 <- gstat::variogram(residuals_gls ~ 1, data = defol_only_sev_3)

# Plot the variogram
var_plot_wo3 <- plot(variogram_gls_wo3, main = "Variogram of Residuals (Gaussian)")

# confirm defol_only as dataframe
df.sev_wo3 <- as.data.frame(defol_only_sev_3)

# Create a spatial plot of residuals
spat.resid.sev_plot_wo3 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_wo3, aes(x = x, y = y, color = residuals_gls)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

#save rds
saveRDS(spat.resid.sev_plot_wo1, file.path(base_path, "/plots/maps/residuals_map_severity_wo1.rds"))


# save spatial plot of residuals
ggsave(
  filename = file.path(base_path, "/plots/maps/residuals_map_severity_wo1.png"),
  plot = spat.resid.sev_plot_wo1,
  width = 12,
  height = 8,
  dpi = 300
)

# test MSR

# Create spatial weights
coords <- cbind(defol_only_sev_3$x, defol_only_sev_3$y)

#
nb <- try(dnearneigh(defol_only_sev_3, 0, 0.5), silent = TRUE)
if (inherits(nb, "try-error") || any(card(nb) == 0)) {
  # If some points have no neighbors, try increasing the distance
  cat("Some points have no neighbors with distance 0.5, trying 1.0...\n")
  nb <- dnearneigh(defol_only_sev_3, 0, 1.0)
  
  # Check if we still have isolated points
  if (any(card(nb) == 0)) {
    cat("Some points still have no neighbors with distance 1.0, trying 2.0...\n")
    nb <- dnearneigh(defol_only_sev_3, 0, 2.0)
    
    # If still problematic, use k-nearest neighbors instead
    if (any(card(nb) == 0)) {
      cat("Using k-nearest neighbors (k=5) instead of distance-based...\n")
      nb <- knn2nb(knearneigh(coordinates(defol_only_sev_3), k=5))
    }
  }
}


listw <- nb2listw(nb, style = "W")

# Compute Moran's Eigenvector Maps (MEMs)
mem <- scores.listw(listw, MEM.autocor = "positive")


# Assume 'resid' is a vector of residuals from a model
resid <- residuals(model_gls_sev_wo3, type = "normalized")

# Moran Spectral Randomization test
msr_test <- moran.randtest(resid, listw)

# Plot and interpret
plot(msr_test)

#A significant result indicates non-random spatial structure in the residuals.








# Part 2: RECOVERY ==================

# Define the model formula
formula.rec <- recovery ~host_pct +  Cumulative_Years_Defol:window_opp + rbr_w_offset +
  mean_temperature + sum_precipitation_mm + mean_tri+ x + y 

# Fit the model using nlme with spatial correlation, where correlation decreases with distance
# gls

model_gls_rec <- gls(formula.rec,
                     correlation = corGaus(form = ~ x + y), 
                     data = defol_only_rec)

# Residuals vs. Fitted Values Plot
# With these lines:
resid.rec_data <- data.frame(
  fitted = fitted(model_gls_rec),
  residuals = resid(model_gls_rec)
)

resid_plot_recovery <- ggplot(resid.rec_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid") +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw()



#save residuals plot 
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/spatial_analysis/residuals_fitted_recovery.png",
  plot = resid_plot_recovery,
  width = 12,
  height = 8,
  dpi = 300
)

# Display the model summary
mod.sum_rec <-summary(model_gls_rec)
print(mod.sum_rec)
saveRDS(mod.sum_rec, "/home/goldma34/sbw-wildfire-impact-recovery/results/spatial_analysis/model_summary_recovery.rds")
# r2
r2.rec <-rsquared.gls(model_gls_rec, formula.rec)
print(r2.rec)
saveRDS(r2.rec, "/home/goldma34/sbw-wildfire-impact-recovery/results/spatial_analysis/r2_recovery.rds")

# extract residuals
defol_only_rec$residuals_gls_rec <- residuals(model_gls_rec)

colnames(defol_only_rec) <- make.names(colnames(defol_only_rec))

# get coords
coordinates(defol_only_rec) <- ~ x + y

# Compute the variogram

variogram_gls_rec <- variogram(residuals_gls_rec ~ 1, data = defol_only_rec)

# Plot the variogram
var_plot_rec <- plot(variogram_gls_rec, main = "Variogram of Residuals (Gaussian)")

# confirm defol_only as dataframe
df_rec <- as.data.frame(defol_only_rec)

# Create a spatial plot of residuals
spat.resid.rec_plot <-ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df_rec, aes(x = x, y = y, color = residuals_gls_rec)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

#save rds
saveRDS(spat.resid.rec_plot, "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/residuals_map_recovery.rds")

# save spatial plot of residuals
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/residuals_map_recovery.png",
  plot = spat.resid.rec_plot,
  width = 12,
  height = 8,
  dpi = 300
)


# summary tables
export_gls_coefficient_tables(model_gls, model_gls_rec)




