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

options(scipen = 999)  # Disable scientific notation for better readability


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


# Window of Opportunity 1 =========================

# Define the model formula
formula.sev <- rbr_w_offset ~host_pct + Cumulative_Years_Defol + isi_90 + 
  dc_90 + dmc_90 + ffmc_90 + bui_90+ fwi_90 + mean_tri

 # Identify mems
 window_opp_1_mems <- get_mems(defol_only_sev_1)

#get data
window_opp_1_mems_data <- window_opp_1_mems$data_w_first_mem
# get weights
window_opp_1_weights <- window_opp_1_mems$weights
# view the first MEM
print(first_mem <- window_opp_1_mems$first_mem)

# LM model
window_opp_1_mod.1 <- lm_model(window_opp_1_mems_data, window_opp_1_weights, formula.sev, null_model = TRUE)

# check moran's I - is it significantly different from 0?
morans_i_wo1 <-window_opp_1_mod.1$morans_i
morans_i_wo1

# check model summary
mod.sum_sev_1 <- summary(window_opp_1_mod.1$model)

# check adj r-squared
print(r2_wo1 <- window_opp_1_mod.1$adj_r2_model)
print(r2_w1_null <-window_opp_1_mod.1$adj_r2_null)

# whats the formula 
eval(mod.sum_sev_1$call[[2]])

# check variogram
variogram_wo1 <- window_opp_1_mod.1$vario
plot(variogram_wo1, main = "Variogram of Residuals (Gaussian)")

#identify when gamma begins to go down
max_gamma_row_wo1 <- variogram_wo1[which.max(variogram_wo1$gamma), ]
(max_gamma_row_wo1$dist / 1000) # convert to km

# Create a dataframe with residuals and coordinates
df.sev_1 <- data.frame(
  x = defol_only_sev_1$x,
  y = defol_only_sev_1$y,
  residuals_lm_mem = residuals(window_opp_1_mod.1$model)
)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo1 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_1, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

spat.resid.sev_plot_wo1

#
# Window of Opptunity 2 ==========================


 # Identify mems
 window_opp_2_mems <- get_mems(defol_only_sev_2)

#get data
window_opp_2_mems_data <- window_opp_2_mems$data_w_first_mem
# get weights
window_opp_2_weights <- window_opp_2_mems$weights
# view the first MEM
print(first_mem.2 <- window_opp_2_mems$first_mem)

# LM model
window_opp_2_mod.1 <- lm_model(window_opp_2_mems_data, window_opp_2_weights, formula.sev, null_model = TRUE)

# check moran's I - is it significantly different from 0?
morans_i_wo2 <-window_opp_2_mod.1$morans_i
morans_i_wo2

# check model summary
mod.sum_sev_2 <- summary(window_opp_2_mod.1$model)

# check adj r-squared
print(r2_wo2 <- window_opp_2_mod.1$adj_r2_model)
print(r2_w2_null <-window_opp_2_mod.1$adj_r2_null)

# whats the formula 
eval(mod.sum_sev_2$call[[2]])

# check variogram
variogram_wo2 <- window_opp_2_mod.1$vario
plot(variogram_wo2, main = "Variogram of Residuals (Gaussian)")

#identify when gamma begins to go down
max_gamma_row <- variogram_wo2[which.max(variogram_wo2$gamma), ]
(max_gamma_row$dist / 1000) # convert to km

# Create a dataframe with residuals and coordinates
df.sev_2 <- data.frame(
  x = defol_only_sev_2$x,
  y = defol_only_sev_2$y,
  residuals_lm_mem = residuals(window_opp_2_mod.1$model)
)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo2 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_2, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

spat.resid.sev_plot_wo2

# Window of Opptunity 3 ==========================
# MEM as spatial predictor

 #dentify mems
 window_opp_3_mems <- get_mems(defol_only_sev_3)

#get data
window_opp_3_mems_data <- window_opp_3_mems$data_w_first_mem
# get weights
window_opp_3_weights <- window_opp_3_mems$weights
# view the first MEM
print(first_mem.3 <- window_opp_3_mems$first_mem)

# LM model
window_opp_3_mod.1 <- lm_model(window_opp_3_mems_data, window_opp_3_weights, formula.sev, null_model = TRUE)

# check moran's I - is it significantly different from 0?
morans_i_wo3 <-window_opp_3_mod.1$morans_i
morans_i_wo3

# check model summary
mod.sum_sev_3 <- summary(window_opp_3_mod.1$model)

# check adj r-squared
print(r2_wo3 <- window_opp_3_mod.1$adj_r2_model)
print(r2_w3_null <-window_opp_3_mod.1$adj_r2_null)

# whats the formula 
eval(mod.sum_sev_3$call[[2]])

# check variogram
variogram_wo3 <- window_opp_3_mod.1$vario
plot(variogram_wo3, main = "Variogram of Residuals (Gaussian)")

#identify when gamma begins to go down
max_gamma_row3<- variogram_wo3[which.max(variogram_wo3$gamma), ]
(max_gamma_row3$dist / 1000) # convert to km

# Create a dataframe with residuals and coordinates
df.sev_3 <- data.frame(
  x = defol_only_sev_3$x,
  y = defol_only_sev_3$y,
  residuals_lm_mem = residuals(window_opp_3_mod.1$model)
)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo3 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_3, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

spat.resid.sev_plot_wo3

# Window of Opptunity 4 ==========================
#dentify mems
window_opp_4_mems <- get_mems(defol_only_sev_4)

#get data
window_opp_4_mems_data <- window_opp_4_mems$data_w_first_mem
# get weights
window_opp_4_weights <- window_opp_4_mems$weights
# view the first MEM
print(first_mem.4 <- window_opp_4_mems$first_mem)

# LM model
window_opp_4_mod.1 <- lm_model(window_opp_4_mems_data, window_opp_4_weights, formula.sev, null_model = TRUE)

# check moran's I - is it significantly different from 0?
morans_i_wo4 <-window_opp_4_mod.1$morans_i
morans_i_wo4

# check model summary
mod.sum_sev_4 <- summary(window_opp_4_mod.1$model)

# check adj r-squared
print(r2_wo4 <- window_opp_4_mod.1$adj_r2_model)
print(r2_w4_null <-window_opp_4_mod.1$adj_r2_null)

# whats the formula 
eval(mod.sum_sev_4$call[[2]])

# check variogram
variogram_wo4 <- window_opp_4_mod.1$vario
plot(variogram_wo4, main = "Variogram of Residuals (Gaussian)")

#identify when gamma begins to go down
max_gamma_row4<- variogram_wo4[which.max(variogram_wo4$gamma), ]
(max_gamma_row4$dist / 1000) # convert to km

# fit model with second mem
# get second row order 
second_mem <- window_opp_4_mems$mems_selected[2,"order"]
second_mem <- window_opp_4_mems$mems[, second_mem]
#add second mem to data
# as data.frame
second_mem <- as.data.frame(second_mem)
# rename the column to first_mem
colnames(second_mem) <- "second_mem"
window_opp_4_mems_data <- cbind(window_opp_4_mems$data_w_first_mem, second_mem =second_mem)

# Fit the model with the second MEM
formula <- update(formula, . ~ . + first_mem + second_mem)
  
mem_model <- lm(formula, data = window_opp_4_mems_data)
# vif
mem_model_vif <- identify_high_vif(mem_model, threshold = 10)
  
# check vifs
mem_list_vif <- mem_model_vif$high_vif_vars
  
# Get all terms from the right side of the formula
all_terms_mem_model <- attr(terms(formula), "term.labels")
  
# Filter out the variables to remove
terms_to_keep_mem_model <- setdiff(all_terms_mem_model, mem_list_vif)
  
# Create a new formula
formula_vif_check <- reformulate(termlabels = terms_to_keep_mem_model, response = response_var)
# Fit the large scale model without high VIF variables
model_vif_removed <- lm(formula_vif_check, data = window_opp_4_mems_data)
# get adjusted r2
adj_r2_model <- summary(model_vif_removed)$adj.r.squared

# get moran's I
moran_i <- moranNP.randtest(residuals(model_vif_removed), listw, nrepet = 999, alter = "two-sided") 

# get variogram

vario <- variogram(residuals(model_vif_removed) ~ 1, data = window_opp_4_mems_data)


# Create a dataframe with residuals and coordinates
df.sev_4 <- data.frame(
  x = defol_only_sev_4$x,
  y = defol_only_sev_4$y,
  residuals_lm_mem = residuals(model_vif_removed)
)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo4 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_4, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

spat.resid.sev_plot_wo4
































