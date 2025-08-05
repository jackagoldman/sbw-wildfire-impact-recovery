# Spatial Analysis

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

#utils
source("/home/goldma34/sbw-wildfire-impact-recovery/src/utils.R")


# Load your data 
source("/home/goldma34/sbw-wildfire-impact-recovery/src/load_data.R") 

#load
source("/home/goldma34/sbw-wildfire-impact-recovery/src/best_model_functions.R") 


# Load Ontario map data
ontario <- ne_states(country = "Canada", returnclass = "sf") %>%
  filter(name == "Ontario")

# Create the inset pie chart for m.data.sf
inset_data1 <- h90.sf %>%
  group_by(history) %>%
  summarise(count = n())

inset_pie1 <- ggplot(inset_data1, aes(x = "", y = count, fill = factor(history))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 3)  # Adjust the size here

# Create the inset pie chart for m.data.sf2
inset_data2 <- m.data.sev_sf %>%
  group_by(history) %>%
  summarise(count = n())

inset_pie2 <- ggplot(inset_data2, aes(x = "", y = count, fill = factor(history))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 3)  # Adjust the size here

# Create the first map with inset
map1 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = h90.sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                     name = NULL, 
                     labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                     guide = guide_legend(title.position = "top", title.hjust = 0.4)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = "right",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) +
  labs(title = "A)") +
  annotation_custom(ggplotGrob(inset_pie1), xmin = -78, xmax = -73, ymin = 42, ymax = 52)

  #save rds
saveRDS(map1, "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/study_area_unmatched_map.rds")

# Create the second map with inset using m.data.sf2
map2 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = m.data.sev_sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                     name = NULL, 
                     labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                     guide = guide_legend(title.position = "top", title.hjust = 0.4)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) +
  labs(title = "B)") +
  annotation_custom(ggplotGrob(inset_pie2), xmin = -78, xmax = -73, ymin = 42, ymax = 52)

map2
#save rds
saveRDS(map2, "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/study_area_severity_map.rds")

# Create the inset pie chart for m.data.rec_sf
inset_data3 <- m.data.rec_sf %>%
  group_by(history) %>%
  summarise(count = n())

inset_pie3 <- ggplot(inset_data3, aes(x = "", y = count, fill = factor(history))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 3)  # Adjust the size here

# recovery map using m.data.rec_sf
map3 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = m.data.rec_sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                     name = NULL, 
                     labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                     guide = guide_legend(title.position = "top", title.hjust = 0.4)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) +
  labs(title = "C)") +
  annotation_custom(ggplotGrob(inset_pie3), xmin = -78, xmax = -73, ymin = 42, ymax = 52)

# save rds
saveRDS(map3, "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/study_area_recovery_map.rds")

# Extract the legend from one of the maps
legend <- cowplot::get_legend(map1)

# Combine the legend and the two maps
combined_plot <- plot_grid(
  legend,
  plot_grid(
    map1 + theme(legend.position = "none"),
    map2 + theme(legend.position = "none"),
    ncol = 2,
    rel_widths = c(1, 1)
  ),
  ncol = 1,
  rel_heights = c(0.1, 1)
)


# save each plot inidividually
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/study_area_unmatched.png",
  plot = map1,
  width = 12,
  height = 8,
  dpi = 300
)
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/study_area_matched_severity.png",
  plot = map2,
  width = 12,
  height = 8,
  dpi = 300
)

# Save the combined plot    
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/study_area_and_sev_matched.png",
  plot = combined_plot,
  width = 12,
  height = 8,
  dpi = 300
)

######################################
# SPATIAL DISTRIBUTION MODELS
######################################

# defol dataset - models are for defolaited only fires, make two dataframes for attach residuals and plotting
defol_only_sev <- subset(history_gt90, history == 1)
defol_only_rec <- subset(history_gt90, history == 1)



########################################
# Part 1: SEVERITY ======================
########################################
# MEM ================


# Define the model formula
formula.sev <- rbr_w_offset ~ history + host_pct + Cumulative_Years_Defol:window_opp + isi_90 + 
  dc_90 + dmc_90 + ffmc_90 + bui_90+ fwi_90 + mean_tri


# Identify mems
full_sev_mems <- get_mems(m.data)

#get data
full_sev_mems_data <- full_sev_mems$data_w_ten_mems
# get weights
full_sev_weights <- full_sev_mems$weights

# standardize the data



cols_to_scale <- c("host_pct",  "isi_90", "dc_90", 
                   "dmc_90", "ffmc_90", "bui_90", "fwi_90")

# Scale the selected columns in the @data slot
full_sev_mems_data@data[cols_to_scale] <- scale(full_sev_mems_data@data[cols_to_scale])


# LM model
full_sev.1 <- lm_model(full_sev_mems_data, full_sev_weights, formula.sev, null_model = TRUE, mems= c("first_mem", "second_mem", "third_mem", "fourth_mem"))

# check moran's I - is it significantly different from 0?
morans_i_wo1 <-full_sev.1$morans_i
morans_i_wo1

# check model summary
full.sum_sev_1 <- summary(full_sev.1$model)

# avg_comparisons for treatment effect
fit_att_sev_test <- as.data.frame(marginaleffects::avg_comparisons(full_sev.1$model,
                                                          variables = "history",
                                                          vcov = ~subclass,
                                                          newdata = subset(history == 1))) %>%  mutate(Fires= "all")

# check adj r-squared
print(r2_wo1 <- full_sev.1$adj_r2_model)
print(r2_w1_null <-full_sev.1$adj_r2_null)

# whats the formula 
eval(mod.sum_sev_1$call[[2]])

# check variogram
variogram_wo1 <- full_sev.1$vario
plot(variogram_wo1, main = "Variogram of Residuals (Gaussian)")

#identify when gamma begins to go down
max_gamma_row_wo1 <- variogram_wo1[which.max(variogram_wo1$gamma), ]
(max_gamma_row_wo1$dist / 1000) # convert to km

# Create a dataframe with residuals and coordinates
df.sev_1 <- data.frame(
  x = m.data$x,
  y = m.data$y,
  residuals_lm_mem = residuals(full_sev.1$model)
)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_wo1 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_1, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

spat.resid.sev_plot_wo1





# Define the model formula =====================
formula.sev <- rbr_w_offset ~ history + host_pct + Cumulative_Years_Defol:window_opp + isi_90 + 
  dc_90 + dmc_90 + ffmc_90 + bui_90+ fwi_90 + mean_tri+  x + y 

# Fit the model using nlme with spatial correlation, where correlation decreases with distance
model_gls <- gls(formula.sev,
                 correlation = corGaus(form = ~ x + y), 
                 data = m.data,
                 control = list(singular.ok = TRUE))

# Check for multicollinearity - remove anything > 10
print(vif(model_gls))

# Residuals vs. Fitted Values Plot

# With these lines:
resid_data <- data.frame(
  fitted = fitted(model_gls),
  residuals = resid(model_gls)
)

resid_plot_severity <- ggplot(resid_data, aes(x = fitted, y = residuals)) +
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
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/spatial_analysis/residuals_fitted_severity.png",
  plot = resid_plot_severity,
  width = 12,
  height = 8,
  dpi = 300
)

# Display the model summary
mod.sum_sev <- summary(model_gls)
print(mod.sum_sev)
saveRDS(mod.sum_sev, "/home/goldma34/sbw-wildfire-impact-recovery/results/spatial_analysis/model_summary_severity.rds")
# r2
r2.sev <- rsquared.gls(model_gls, formula.sev)
print(r2.sev)
saveRDS(r2.sev, "/home/goldma34/sbw-wildfire-impact-recovery/results/spatial_analysis/r2_severity.rds")

# extract residuals
defol_only_sev$residuals_gls <- residuals(model_gls)

colnames(defol_only_sev) <- make.names(colnames(defol_only_sev))

# get coords
coordinates(defol_only_sev) <- ~ x + y

# Compute the variogram
variogram_gls <- variogram(residuals_gls ~ 1, data = defol_only_sev)

# Plot the variogram
var_plot <- plot(variogram_gls, main = "Variogram of Residuals (Gaussian)")

# confirm defol_only as dataframe
df.sev <- as.data.frame(defol_only_sev)

# Create a spatial plot of residuals
spat.resid.sev_plot <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev, aes(x = x, y = y, color = residuals_gls)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

#save rds
saveRDS(spat.resid.sev_plot, "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/residuals_map_severity.rds")


# save spatial plot of residuals
ggsave(
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/residuals_map_severity.png",
  plot = spat.resid.sev_plot ,
  width = 12,
  height = 8,
  dpi = 300
)
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



# Part 3 - NON defoliated fires =====================
#non-defol dataset - models are for defolaited only fires, make two dataframes for attach residuals and plotting
non_defol_only_sev <- subset(history_gt90, history == 0)
non_defol_only_rec <- subset(history_gt90, history == 0)


# remove CYD from the formula for non-defoliated fires
formula.non.sev <- rbr_w_offset ~host_pct +  isi_90 + 
  dc_90 + dmc_90 + ffmc_90 + bui_90+ fwi_90 + mean_tri




#dentify mems
non_defol_mems <- get_mems(non_defol_only_sev)

#get data
non_defol_mems_data <- defol_mems$data_w_first_mem
# get weights
non_defol_weights <- defol_mems$weights
# view the first MEM
print(first_mem.ndefol <- non_defol_mems$first_mem)


# LM model
ndefol_mod.1 <- lm_model(non_defol_mems_data, non_defol_weights, formula.non.sev, null_model = TRUE)



# check moran's I - is it significantly different from 0?
morans_i_nd <-ndefol_mod.1$morans_i
morans_i_nd

# check model summary
mod.sum_sev_nd <- summary(ndefol_mod.1$model)

# check adj r-squared
print(r2_nd <- ndefol_mod.1$adj_r2_model)
print(r2_nd_null <-ndefol_mod.1$adj_r2_null)

# whats the formula 
eval(mod.sum_sev_nd$call[[2]])

# check variogram
variogram_nd <- ndefol_mod.1$vario
plot(variogram_nd, main = "Variogram of Residuals (Gaussian)")

##add second mem to data
# fit model with second mem
# get second row order 
second_mem <- non_defol_mems$mems_selected[2,"order"]
second_mem <- non_defol_mems$mems[, second_mem]
#add second mem to data
# as data.frame
second_mem <- as.data.frame(second_mem)
# rename the column to first_mem
colnames(second_mem) <- "second_mem"
nd_mems_data <- cbind(non_defol_mems$data_w_first_mem, second_mem =second_mem)

# Fit the model with the second MEM
formula <- update(formula.non.sev, . ~ . + first_mem + second_mem)
  
mem_model <- lm(formula, data = nd_mems_data)
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
model_vif_removed <- lm(formula_vif_check, data = nd_mems_data)
# get adjusted r2
adj_r2_model <- summary(model_vif_removed)$adj.r.squared

# get moran's I
moran_i <- moranNP.randtest(residuals(model_vif_removed), listw, nrepet = 999, alter = "two-sided") 

# get variogram

vario <- variogram(residuals(model_vif_removed) ~ 1, data = nd_mems_data)

#identify when gamma begins to go down
max_gamma_row_nd<- vario[which.max(vario$gamma), ]
(max_gamma_row_nd$dist / 1000) # convert to km

# Create a dataframe with residuals and coordinates
df.sev_nd <- data.frame(
  x = non_defol_sev$x,
  y = non_defol_sev$y,
  residuals_lm_mem = residuals(mem_model)
)

# Create a spatial plot of residuals with lm
spat.resid.sev_plot_nd <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.sev_nd, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

spat.resid.sev_plot_nd

#save rds
saveRDS(spat.resid.sev_plot_nd, "./plots/maps/residuals_map_nd_sev.rds")

# save spatial plot of residuals
ggsave(
  filename = "./plots/maps/residuals_map_nd_sev.png",
  plot = spat.resid.sev_plot_nd,
  width = 12,
  height = 8,
  dpi = 300
)


#Part 4 nd Recovery Models =========================

# remove CYD from the formula for non-defoliated fires
formula.non.rec <- recovery ~host_pct  + rbr_w_offset + mean_temperature + sum_precipitation_mm + mean_tri 




#dentify mems
non_defol_mems.rec <- get_mems(non_defol_only_rec)

#get data
non_defol_mems_data.rec <- non_defol_mems.rec$data_w_first_mem
# get weights
non_defol_weights.rec <- non_defol_mems.rec$weights
# view the first MEM
print(first_mem.ndefol.rec <- non_defol_mems.rec$first_mem)


# LM model
ndefol_mod.1.rec <- lm_model(non_defol_mems_data.rec, non_defol_weights.rec, formula.non.rec, null_model = TRUE)



# check moran's I - is it significantly different from 0?
morans_i_nd.rec <-ndefol_mod.1.rec$morans_i
morans_i_nd.rec

# check model summary
mod.sum_rec_nd.rec <- summary(ndefol_mod.1.rec$model)

# check adj r-squared
print(r2_nd <- ndefol_mod.1$adj_r2_model)
print(r2_nd_null <-ndefol_mod.1$adj_r2_null)

# whats the formula 
eval(mod.sum_rec_nd$call[[2]])

# check variogram
variogram_nd <- ndefol_mod.1$vario
plot(variogram_nd, main = "Variogram of Residuals (Gaussian)")

##add second mem to data
# fit model with second mem
# get second row order 
second_mem.rec <- non_defol_mems.rec$mems_selected[2,"order"]
second_mem.rec <- non_defol_mems.rec$mems[, second_mem]
#add second mem to data
# as data.frame
second_mem.rec <- as.data.frame(second_mem.rec)
# rename the column to first_mem
colnames(second_mem.rec) <- "second_mem"
nd_mems_data.rec <- cbind(non_defol_mems_data.rec, second_mem =second_mem.rec)

# Fit the model with the second MEM
formula.rec <- update(formula.non.rec, . ~ . + first_mem + second_mem)
  
mem_model.rec <- lm(formula.rec, data = nd_mems_data.rec)
# vif
mem_model_vif.rec <- identify_high_vif(mem_model.rec, threshold = 10)
  
# check vifs
mem_list_vif.rec <- mem_model_vif.rec$high_vif_vars
  
# Get all terms from the right side of the formula
all_terms_mem_model <- attr(terms(formula), "term.labels")
  
# Filter out the variables to remove
terms_to_keep_mem_model <- setdiff(all_terms_mem_model, mem_list_vif)
  
# Create a new formula
formula_vif_check <- reformulate(termlabels = terms_to_keep_mem_model, response = response_var)
# Fit the large scale model without high VIF variables
model_vif_removed <- lm(formula_vif_check, data = nd_mems_data)
# get adjusted r2
adj_r2_model <- summary(model_vif_removed)$adj.r.squared

# get moran's I
moran_i <- moranNP.randtest(residuals(model_vif_removed), listw, nrepet = 999, alter = "two-sided") 

# get variogram

vario <- variogram(residuals(model_vif_removed) ~ 1, data = nd_mems_data)

#identify when gamma begins to go down
max_gamma_row_nd<- vario[which.max(vario$gamma), ]
(max_gamma_row_nd$dist / 1000) # convert to km

# Create a dataframe with residuals and coordinates
df.rec_nd <- data.frame(
  x = non_defol_rec$x,
  y = non_defol_rec$y,
  residuals_lm_mem = residuals(mem_model)
)

# Create a spatial plot of residuals with lm
spat.resid.rec_plot_nd <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df.rec_nd, aes(x = x, y = y, color = residuals_lm_mem)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

spat.resid.rec_plot_nd

#save rds
saveRDS(spat.resid.rec_plot_nd, "./plots/maps/residuals_map_nd_rec.rds")

# save spatial plot of residuals
ggsave(
  filename = "./plots/maps/residuals_map_nd_rec.png",
  plot = spat.resid.rec_plot_nd,
  width = 12,
  height = 8,
  dpi = 300
)



























