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
library(adesgraphics)


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
  dc_90 + dmc_90 + ffmc_90 + bui_90+ fwi_90 + mean_tri+  x + y 


# MEM as spatial predictor

# Create spatial weights
# coordinates as matrix
coords_1 <- as.matrix(coordinates(defol_only_sev_1))

# k near neighbor to find max distance for d2
knn1 <- knearneigh(coords_1, k = 2)
nbknn1 <- knn2nb(knn1, sym = TRUE)
# Calculate distances
knn_dists <- nbdists(nbknn1, coords_1)

# Upper threshold (d2)
d2 <- max(unlist(knn_dists))

# Use dnearneigh with the calculated d2
d <- dnearneigh(coords_1, d1 = 0, d2 = d2)

# transform coords  to data frame for plotting
coords_df <- as.data.frame(coords_1)

# Create a plot to visualize the neighbors
try({
  plot(d, coords_df, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d2))
}, silent = TRUE)

# how many disjointed subgraphs
n.comp.nb(d)

# Create neighbors using different methods
nbtri <- tri2nb(coords_1)
nbgab <- graph2nb(gabrielneigh(coords_1), sym = TRUE)
nbrel <- graph2nb(relativeneigh(coords_1), sym = TRUE)

# Create a plot to visualize the neighbors
try({
  plot(nbtri, coords_df, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d2))
},silent = TRUE)

try({
  plot(nbgab, coords_df, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d2))
},silent = TRUE)

try({
  plot(nbrel, coords_df, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d2))
},silent = TRUE)

# how many disjointed subgraphs
n.comp.nb(nbtri)
n.comp.nb(nbgab)
n.comp.nb(nbrel)

# define spatial weight matrix
#weight edges between neighbours as a function of spatial distances
distgab <- nbdists(nbgab, coords_1)
fdist <- lapply(distgab, function(x) 1 - x/max(dist(coords_1)))
# Create a listw object
listw_gab <- nb2listw(nbgab, glist = fdist, style = "W", zero.policy = TRUE)

#Create MEM
mem.gab <- mem(listw_gab)
mem.gab

# select best mem
selected <- forward.sel(Y = as.data.frame(defol_only_sev_1$rbr_w_offset), X = as.data.frame(mem.gab))
selected_mem <- mem.gab[, selected$order]
# Print selected MEMs
print(selected_mem)
# selecte first MEM
first_mem <- selected_mem[, 1]
# Print selected MEM
print(first_mem)

# cbind selected MEMs to the original data
defol_only_sev_1 <- cbind(defol_only_sev_1, first_mem)

# LM
# modify formula to include selected MEMs
formula.sev_mem <- update(formula.sev, . ~ . + first_mem)
# remove x and y from the formula
formula.sev_mem <- update(formula.sev_mem, . ~ . - x - y)
# fit the model using LM
model_lm_sev_mem <- lm(formula.sev_mem, data = defol_only_sev_1)
# Display the model summary
mod.sum_sev_mem <- summary(model_lm_sev_mem)
# print
print(mod.sum_sev_mem)

#check multicollinearity
vif_results <- vif(model_lm_sev_mem)
# Print the VIF results
print(vif_results)
# Identify variables with high VIF
identify_high_vif <- function(model, threshold = 10) {
  vif_values <- vif(model)
  high_vif_vars <- names(vif_values[vif_values > threshold])
  
  if (length(high_vif_vars) > 0) {
    cat("Variables with VIF greater than", threshold, ":\n")
    print(high_vif_vars)
  } else {
    cat("No variables with VIF greater than", threshold, "\n")
  }
  
  return(list(high_vif_vars = high_vif_vars))
}
# Identify variables with high VIF
high_vif_results <- identify_high_vif(model_lm_sev_mem, threshold = 10)

# Get the list of variables with high VIF
high_vif_vars <- high_vif_results$high_vif_vars
# Print the variables with high VIF
print(high_vif_vars)

# remove high VIF variables from the model
formula.wo1 <- update(formula.sev_mem, . ~ . - isi_90 - dc_90 - dmc_90  - bui_90 - fwi_90)

#refit the model again without high VIF variables
model_lm_sev_mem_wo1 <- lm(formula.wo1, data = defol_only_sev_1)
# Display the model summary
mod.sum_sev_mem_wo1 <- summary(model_lm_sev_mem_wo1)
# print
print(mod.sum_sev_mem_wo1)
#check multicollinearity
vif_results_wo1 <- vif(model_lm_sev_mem_wo1)
# Print the VIF results
print(vif_results_wo1)

#get residuals
defol_only_sev_1$residuals_lm_mem <- residuals(model_lm_sev_mem)

## Compute the variogram
variogram_lm_wo1 <- gstat::variogram(residuals_lm_mem ~ 1, data = defol_only_sev_1)

# Plot the variogram
lm_var_plot <- plot(variogram_lm_wo1, main = "Variogram of Residuals (Gaussian)")

print(lm_var_plot)

# confirm defol_only as dataframe
df.sev_1 <- as.data.frame(defol_only_sev_1)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo1 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_1, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

# moran's I
NP.rbr_w1 <- moranNP.randtest(defol_only_sev_1$rbr_w_offset, listw_gab, nrepet = 999, alter = "two-sided") 
NP.rbr_w1

# visualize moran's I
plot(NP.rbr_w1, main = "Moran's I for Residuals (Gaussian)")

# moran's I
NP.host_w1 <- moranNP.randtest(defol_only_sev_1$host_pct ,listw_gab, nrepet = 999, alter = "two-sided") 
NP.host_w1

# visualize moran's I
plot(NP.host_w1, main = "Moran's I for Residuals (Gaussian)")

# Window of Opptunity 2 ==========================
# MEM as spatial predictor

# Create spatial weights
# coordinates as matrix
coords_2 <- as.matrix(coordinates(defol_only_sev_2))

# k near neighbor to find max distance for d2
knn2 <- knearneigh(coords_2, k = 2)
nbknn2 <- knn2nb(knn2, sym = TRUE)
# Calculate distances
knn_dists2 <- nbdists(nbknn2, coords_2)

# Upper threshold (d2)
d2.2 <- max(unlist(knn_dists2))

# Use dnearneigh with the calculated d2
d.2 <- dnearneigh(coords_2, d1 = 0, d2 = d2.2)

# transform coords  to data frame for plotting
coords_df2 <- as.data.frame(coords_2)

# Create a plot to visualize the neighbors
try({
  plot(d.2, coords_df2, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d2.2))
}, silent = TRUE)

# how many disjointed subgraphs
n.comp.nb(d.2)

# Create neighbors using different methods
nbtri.2 <- tri2nb(coords_2)
nbgab.2 <- graph2nb(gabrielneigh(coords_2), sym = TRUE)
nbrel.2 <- graph2nb(relativeneigh(coords_2), sym = TRUE)

# Create a plot to visualize the neighbors
try({
  plot(nbtri.2, coords_df2, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d2))
},silent = TRUE)

try({
  plot(nbgab.2, coords_df2, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d2))
},silent = TRUE)

try({
  plot(nbrel.2, coords_df2, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d2))
},silent = TRUE)

# how many disjointed subgraphs
n.comp.nb(nbtri.2)
n.comp.nb(nbgab.2)
n.comp.nb(nbrel.2)

# define spatial weight matrix
#weight edges between neighbours as a function of spatial distances
distgab.2 <- nbdists(nbgab.2, coords_2)
fdist.2 <- lapply(distgab.2, function(x) 1 - x/max(dist(coords_2)))
# Create a listw object
listw_gab.2 <- nb2listw(nbgab.2, glist = fdist.2, style = "W", zero.policy = TRUE)

#Create MEM
mem.gab.2 <- mem(listw_gab.2)
mem.gab.2
# select best mem
selected.2 <- forward.sel(Y = as.data.frame(defol_only_sev_2$rbr_w_offset), X = as.data.frame(mem.gab.2))
selected_mem.2 <- mem.gab.2[, selected.2$order]
# Print selected MEMs
print(selected_mem)
# selecte first MEM
first_mem <- selected_mem[, 1]
# Print selected MEM
print(first_mem)

# define spatial weight matrix
#weight edges between neighbours as a function of spatial distances
distd.2 <- nbdists(d.2, coords_2)
fdistd.2 <- lapply(distd.2, function(x) 1 - x/max(dist(coords_2)))
# Create a listw object
listw_d.2 <- nb2listw(d.2, glist = fdistd.2, style = "W", zero.policy = TRUE)

#Create MEM
mem.d.2 <- mem(listw_d.2)
mem.d.2
# select best mem
selected.d2 <- forward.sel(Y = as.data.frame(defol_only_sev_2$rbr_w_offset), X = as.data.frame(mem.d.2))
selected_mem.d2 <- mem.d.2[, selected.d2$order]
# Print selected MEMs
print(selected_mem.d2)
# selecte first MEM
first_mem.d2<- selected_mem.d2[, 1]
# Print selected MEM
print(first_mem.d2)

# cbind selected MEMs to the original data
defol_only_sev_2 <- cbind(defol_only_sev_2, first_mem.d2)

# LM
# modify formula to include selected MEMs
formula.sev_mem.d2 <- update(formula.sev, . ~ . + first_mem.d2)
# remove x and y from the formula
formula.sev_mem.d2 <- update(formula.sev_mem.d2, . ~ . - x - y)
# fit the model using LM
model_lm_sev_mem.d2 <- lm(formula.sev_mem.d2, data = defol_only_sev_2)
# Display the model summary
mod.sum_sev_mem.d2 <- summary(model_lm_sev_mem.d2)
# print
print(mod.sum_sev_mem.d2)

#check multicollinearity
vif_results.d2 <- vif(model_lm_sev_mem.d2)
# Print the VIF results
print(vif_results.d2)


# Identify variables with high VIF
high_vif_results.d2 <- identify_high_vif(model_lm_sev_mem.d2, threshold = 10)

# Get the list of variables with high VIF
high_vif_vars.d2 <- high_vif_results.d2$high_vif_vars
# Print the variables with high VIF
print(high_vif_vars.d2)

# remove high VIF variables from the model
formula.wo2 <- update(formula.sev_mem.d2, . ~ . - dc_90 - dmc_90  - bui_90)

#refit the model again without high VIF variables
model_lm_sev_mem_wo2 <- lm(formula.wo2, data = defol_only_sev_2)
# Display the model summary
mod.sum_sev_mem_wo2 <- summary(model_lm_sev_mem_wo2)
# print
print(mod.sum_sev_mem_wo2)
#check multicollinearity
vif_results_wo2 <- vif(model_lm_sev_mem_wo2)
# Print the VIF results
print(vif_results_wo2)

#get residuals
defol_only_sev_2$residuals_lm_mem <- residuals(model_lm_sev_mem_wo2)

## Compute the variogram
variogram_lm_wo2 <- gstat::variogram(residuals_lm_mem ~ 1, data = defol_only_sev_2)

# Plot the variogram
lm_var_plot_wo2 <- plot(variogram_lm_wo2, main = "Variogram of Residuals (Gaussian)")

print(lm_var_plot_wo2)

# confirm defol_only as dataframe
df.sev_2 <- as.data.frame(defol_only_sev_2)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo2 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_2, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

# moran's I
NP.rbr_w2 <- moranNP.randtest(defol_only_sev_2$rbr_w_offset, listw_d.2, nrepet = 999, alter = "two-sided") 
NP.rbr_w2

# visualize moran's I
plot(NP.rbr_w1, main = "Moran's I for Residuals (Gaussian)")

# Window of Opptunity 3 ==========================
# MEM as spatial predictor

# Create spatial weights
# coordinates as matrix
coords_3 <- as.matrix(coordinates(defol_only_sev_3))

# k near neighbor to find max distance for d3
knn3 <- knearneigh(coords_3, k = 2)
nbknn3 <- knn2nb(knn3, sym = TRUE)
# Calculate distances
knn_dists3 <- nbdists(nbknn3, coords_3)

# Upper threshold (d3)
d3.3 <- max(unlist(knn_dists3))

# Use dnearneigh with the calculated d3
d.3 <- dnearneigh(coords_3, d1 = 0, d2 = d3.3)

# transform coords  to data frame for plotting
coords_df3 <- as.data.frame(coords_3)

# Create a plot to visualize the neighbors
try({
  plot(d.3, coords_df3, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d3.3))
}, silent = TRUE)

# how many disjointed subgraphs
n.comp.nb(d.3)

# Create neighbors using different methods
nbtri.3 <- tri3nb(coords_3)
nbgab.3 <- graph3nb(gabrielneigh(coords_3), sym = TRUE)
nbrel.3 <- graph3nb(relativeneigh(coords_3), sym = TRUE)

# Create a plot to visualize the neighbors
try({
  plot(nbtri.3, coords_df3, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d3))
},silent = TRUE)

try({
  plot(nbgab.3, coords_df3, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d3))
},silent = TRUE)

try({
  plot(nbrel.3, coords_df3, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d3))
},silent = TRUE)

# how many disjointed subgraphs
n.comp.nb(nbtri.3)
n.comp.nb(nbgab.3)
n.comp.nb(nbrel.3)

# define spatial weight matrix
#weight edges between neighbours as a function of spatial distances
distgab.3 <- nbdists(nbgab.3, coords_3)
fdist.3 <- lapply(distgab.3, function(x) 1 - x/max(dist(coords_3)))
# Create a listw object
listw_gab.3 <- nb3listw(nbgab.3, glist = fdist.3, style = "W", zero.policy = TRUE)

#Create MEM
mem.gab.3 <- mem(listw_gab.3)
mem.gab.3
# select best mem
selected.3 <- forward.sel(Y = as.data.frame(defol_only_sev_3$rbr_w_offset), X = as.data.frame(mem.gab.3))
selected_mem.3 <- mem.gab.3[, selected.3$order]
# Print selected MEMs
print(selected_mem)
# selecte first MEM
first_mem <- selected_mem[, 1]
# Print selected MEM
print(first_mem)

# define spatial weight matrix
#weight edges between neighbours as a function of spatial distances
distd.3 <- nbdists(d.3, coords_3)
fdistd.3 <- lapply(distd.3, function(x) 1 - x/max(dist(coords_3)))
# Create a listw object
listw_d.3 <- nb2listw(d.3, glist = fdistd.3, style = "W", zero.policy = TRUE)

#Create MEM
mem.d.3 <- mem(listw_d.3)
mem.d.3
# select best mem
selected.d3 <- forward.sel(Y = as.data.frame(defol_only_sev_3$rbr_w_offset), X = as.data.frame(mem.d.3))
selected_mem.d3 <- mem.d.3[, selected.d3$order]
# Print selected MEMs
print(selected_mem.d3)
# selecte first MEM
first_mem.d3<- selected_mem.d3[, 1]
# Print selected MEM
print(first_mem.d3)

# cbind selected MEMs to the original data
defol_only_sev_3 <- cbind(defol_only_sev_3, first_mem.d3)

# LM
# modify formula to include selected MEMs
formula.sev_mem.d3 <- update(formula.sev, . ~ . + first_mem.d3)
# remove x and y from the formula
formula.sev_mem.d3 <- update(formula.sev_mem.d3, . ~ . - x - y)
# fit the model using LM
model_lm_sev_mem.d3 <- lm(formula.sev_mem.d3, data = defol_only_sev_3)
# Display the model summary
mod.sum_sev_mem.d3 <- summary(model_lm_sev_mem.d3)
# print
print(mod.sum_sev_mem.d3)

#check multicollinearity
vif_results.d3 <- vif(model_lm_sev_mem.d3)
# Print the VIF results
print(vif_results.d3)


# Identify variables with high VIF
high_vif_results.d3 <- identify_high_vif(model_lm_sev_mem.d3, threshold = 10)

# Get the list of variables with high VIF
high_vif_vars.d3 <- high_vif_results.d3$high_vif_vars
# Print the variables with high VIF
print(high_vif_vars.d3)

# remove high VIF variables from the model
formula.wo3 <- update(formula.sev_mem.d3, . ~ . - isi_90 - fwi_90 - dc_90 - dmc_90  - bui_90)

#refit the model again without high VIF variables
model_lm_sev_mem_wo3 <- lm(formula.wo3, data = defol_only_sev_3)
# Display the model summary
mod.sum_sev_mem_wo3 <- summary(model_lm_sev_mem_wo3)
# print
print(mod.sum_sev_mem_wo3)
#check multicollinearity
vif_results_wo3 <- vif(model_lm_sev_mem_wo3)
# Print the VIF results
print(vif_results_wo3)

#get residuals
defol_only_sev_3$residuals_lm_mem <- residuals(model_lm_sev_mem_wo3)

## Compute the variogram
variogram_lm_wo3 <- gstat::variogram(residuals_lm_mem ~ 1, data = defol_only_sev_3)

# Plot the variogram
lm_var_plot_wo3 <- plot(variogram_lm_wo3, main = "Variogram of Residuals (Gaussian)")

print(lm_var_plot_wo3)

# confirm defol_only as dataframe
df.sev_3 <- as.data.frame(defol_only_sev_3)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo3 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_3, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

# moran's I
NP.rbr_w3 <- moranNP.randtest(defol_only_sev_3$rbr_w_offset, listw_d.3, nrepet = 999, alter = "two-sided") 
NP.rbr_w3

# visualize moran's I
plot(NP.rbr_w3, main = "Moran's I for Residuals (Gaussian)")

# Window of Opptunity 4 ==========================
# MEM as spatial predictor

# Create spatial weights

coordinates(defol_only_sev_4)<-~ x + y


# coordinates as matrix
coords_4 <- as.matrix(coordinates(defol_only_sev_4))

# k near neighbor to find max distance for d4
knn4 <- knearneigh(coords_4, k = 2)
nbknn4 <- knn2nb(knn4, sym = TRUE)
# Calculate distances
knn_dists4 <- nbdists(nbknn4, coords_4)

# Upper threshold (d4)
d4.4 <- max(unlist(knn_dists4))

# Use dnearneigh with the calculated d4
d.4 <- dnearneigh(coords_4, d1 = 0, d2 = d4.4)

# transform coords  to data frame for plotting
coords_df4 <- as.data.frame(coords_4)

# Create a plot to visualize the neighbors
try({
  plot(d.4, coords_df4, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d4.4))
}, silent = TRUE)

# how many disjointed subgraphs
n.comp.nb(d.4)

# Create neighbors using different methods
nbtri.4 <- tri2nb(coords_4)
nbgab.4 <- graph2nb(gabrielneigh(coords_4), sym = TRUE)
nbrel.4 <- graph2nb(relativeneigh(coords_4), sym = TRUE)

# Create a plot to visualize the neighbors
try({
  plot(nbtri.4, coords_df4, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d4))
},silent = TRUE)

try({
  plot(nbgab.4, coords_df4, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d4))
},silent = TRUE)

try({
  plot(nbrel.4, coords_df4, col = "red", main = "Nearest neighbors")
  title(main = paste("Nearest neighbors d=", d4))
},silent = TRUE)

# how many disjointed subgraphs
n.comp.nb(nbtri.4)
n.comp.nb(nbgab.4)
n.comp.nb(nbrel.4)

# define spatial weight matrix
#weight edges between neighbours as a function of spatial distances
distgab.4 <- nbdists(nbgab.4, coords_4)
fdist.4 <- lapply(distgab.4, function(x) 1 - x/max(dist(coords_4)))
# Create a listw object
listw_gab.4 <- nb2listw(nbgab.4, glist = fdist.4, style = "W", zero.policy = TRUE)

#Create MEM
mem.gab.4 <- mem(listw_gab.4)
mem.gab.4
# select best mem
selected.4 <- forward.sel(Y = as.data.frame(defol_only_sev_4$rbr_w_offset), X = as.data.frame(mem.gab.4))
selected_mem.4 <- mem.gab.4[, selected.4$order]
# Print selected MEMs
print(selected_mem.4)
# selecte first MEM
first_mem.d4 <- selected_mem.4[, 1]
# Print selected MEM
print(first_mem.d4)

# define spatial weight matrix
#weight edges between neighbours as a function of spatial distances
distd.4 <- nbdists(d.4, coords_4)
fdistd.4 <- lapply(distd.4, function(x) 1 - x/max(dist(coords_4)))
# Create a listw object
listw_d.4 <- nb4listw(d.4, glist = fdistd.4, style = "W", zero.policy = TRUE)

#Create MEM
mem.d.4 <- mem(listw_d.4)
mem.d.4
# select best mem
selected.d4 <- forward.sel(Y = as.data.frame(defol_only_sev_4$rbr_w_offset), X = as.data.frame(mem.d.4))
selected_mem.d4 <- mem.d.4[, selected.d4$order]
# Print selected MEMs
print(selected_mem.d4)
# selecte first MEM
first_mem.d4<- selected_mem.d4[, 1]
# Print selected MEM
print(first_mem.d4)

# cbind selected MEMs to the original data
defol_only_sev_4 <- cbind(defol_only_sev_4, first_mem.d4)

# LM
# modify formula to include selected MEMs
formula.sev_mem.d4 <- update(formula.sev, . ~ . + first_mem.d4)
# remove x and y from the formula
formula.sev_mem.d4 <- update(formula.sev_mem.d4, . ~ . - x - y)
# fit the model using LM
model_lm_sev_mem.d4 <- lm(formula.sev_mem.d4, data = defol_only_sev_4)
# Display the model summary
mod.sum_sev_mem.d4 <- summary(model_lm_sev_mem.d4)
# print
print(mod.sum_sev_mem.d4)

#check multicollinearity
vif_results.d4 <- vif(model_lm_sev_mem.d4)
# Print the VIF results
print(vif_results.d4)


# Identify variables with high VIF
high_vif_results.d4 <- identify_high_vif(model_lm_sev_mem.d4, threshold = 10)

# Get the list of variables with high VIF
high_vif_vars.d4 <- high_vif_results.d4$high_vif_vars
# Print the variables with high VIF
print(high_vif_vars.d4)

# remove high VIF variables from the model
formula.wo4 <- update(formula.sev_mem.d4, . ~ . - dc_90 - dmc_90  - fwi_90 - bui_90)

#refit the model again without high VIF variables
model_lm_sev_mem_wo4 <- lm(formula.wo4, data = defol_only_sev_4)
# Display the model summary
mod.sum_sev_mem_wo4 <- summary(model_lm_sev_mem_wo4)
# print
print(mod.sum_sev_mem_wo4)
#check multicollinearity
vif_results_wo4 <- vif(model_lm_sev_mem_wo4)
# Print the VIF results
print(vif_results_wo4)

#get residuals
defol_only_sev_4$residuals_lm_mem <- residuals(model_lm_sev_mem_wo4)

## Compute the variogram
variogram_lm_wo4 <- gstat::variogram(residuals_lm_mem ~ 1, data = defol_only_sev_4)

# Plot the variogram
lm_var_plot_wo4 <- plot(variogram_lm_wo4, main = "Variogram of Residuals (Gaussian)")

print(lm_var_plot_wo4)

# confirm defol_only as dataframe
df.sev_4 <- as.data.frame(defol_only_sev_4)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo4 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_4, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

# moran's I
NP.rbr_w4 <- moranNP.randtest(defol_only_sev_4$rbr_w_offset, listw_gab.4, nrepet = 999, alter = "two-sided") 
NP.rbr_w4

# visualize moran's I
plot(NP.rbr_w4, main = "Moran's I for Residuals (Gaussian)")




























# GLS ==============================

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


# save residuals plot 
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




# window opportunity 2 ==============================

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




