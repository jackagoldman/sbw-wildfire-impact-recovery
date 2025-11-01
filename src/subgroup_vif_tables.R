# Load required libraries
library(MatchIt)
library(dplyr)
library(marginaleffects)
library(MuMIn)
library(pwr)

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

# Source the file with function definitions
source(file.path(base_path, "/src/best_model_functions.R"))

#source the file with covariate balance plots
source(file.path(base_path, "/src/covariate_balance_plots.R"))

# Load your data 
source(file.path(base_path, "/src/load_data.R"))


# read in best models for all subgroups (1-4) and intermediate based on following path  saveRDS(best_model_sev_1, file.path(base_path,"/results/subgroup/best_model_subgroup1_severity.RDS")) # nolint: line_length_linter.
best_model_sev_1 <- readRDS(file.path(base_path,"/results/subgroup/best_model_subgroup1_severity.RDS"))
best_model_rec_1 <- readRDS(file.path(base_path,"/results/subgroup/best_model_subgroup1_recovery.RDS"))
best_model_sev_2 <- readRDS(file.path(base_path,"/results/subgroup/best_model_subgroup2_severity.RDS"))
best_model_rec_2 <- readRDS(file.path(base_path,"/results/subgroup/best_model_subgroup2_recovery.RDS"))
best_model_sev_3 <- readRDS(file.path(base_path,"/results/subgroup/best_model_subgroup3_severity.RDS"))
best_model_rec_3 <- readRDS(file.path(base_path,"/results/subgroup/best_model_subgroup3_recovery.RDS"))
best_model_sev_4 <- readRDS(file.path(base_path,"/results/subgroup/best_model_subgroup4_severity.RDS"))
best_model_rec_4 <- readRDS(file.path(base_path,"/results/subgroup/best_model_subgroup4_recovery.RDS"))
best_model_sev_intermediate <- readRDS(file.path(base_path,"/results/subgroup/best_model_intermediate_severity.RDS"))
best_model_rec_intermediate <- readRDS(file.path(base_path,"/results/subgroup/best_model_intermediate_recovery.RDS"))


# For each best model, match the data
m.data.sev_1 <- match_data(best_model_sev_1)
m.data.rec_1 <- match_data(best_model_rec_1)
m.data.sev_2 <- match_data(best_model_sev_2)
m.data.rec_2 <- match_data(best_model_rec_2)
m.data.sev_3 <- match_data(best_model_sev_3)
m.data.rec_3 <- match_data(best_model_rec_3)
m.data.sev_4 <- match_data(best_model_sev_4)
m.data.rec_4 <- match_data(best_model_rec_4)
m.data.sev_intermediate <- match_data(best_model_sev_intermediate)
m.data.rec_intermediate <- match_data(best_model_rec_intermediate)

# fit severity models with appropriate matched data
 fit.sev_1 <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_1, weights = weights)
 fit.sev_2 <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_2, weights = weights)
 fit.sev_3 <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_3, weights = weights)
 fit.sev_4 <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_4, weights = weights)

#calculate multicollinearity using  vif function from car package
library(car)
vif_sev_1 <- vif(fit.sev_1)
vif_sev_2 <- vif(fit.sev_2)
vif_sev_3 <- vif(fit.sev_3)
vif_sev_4 <- vif(fit.sev_4)

# make tables
vif_table_1 <- data.frame(
  feature = names(vif_sev_1), # remove intercept 
  VIF = as.numeric(vif_sev_1)
)

vif_table_2 <- data.frame(
  feature = names(vif_sev_2), # remove intercept 
  VIF = as.numeric(vif_sev_2)
)

vif_table_3 <- data.frame(
  feature = names(vif_sev_3), # remove intercept 
  VIF = as.numeric(vif_sev_3)
)

vif_table_4 <- data.frame(
  feature = names(vif_sev_4), # remove intercept 
  VIF = as.numeric(vif_sev_4)
)

# write to csv
write.csv(vif_table_1, file.path(base_path,"/results/subgroup/vif_severity_subgroup1.csv"), row.names = FALSE)
write.csv(vif_table_2, file.path(base_path,"/results/subgroup/vif_severity_subgroup2.csv"), row.names = FALSE)
write.csv(vif_table_3, file.path(base_path,"/results/subgroup/vif_severity_subgroup3.csv"), row.names = FALSE)
write.csv(vif_table_4, file.path(base_path,"/results/subgroup/vif_severity_subgroup4.csv"), row.names = FALSE)

# Vif workflow for recovery
fit.rec_1 <- lm(recovery ~ history + host_pct + rbr_w_offset+ mean_temperature + sum_precipitation_mm + mean_tri,
                data = m.data.rec_1, weights = weights)
fit.rec_2 <- lm(recovery ~ history + host_pct + rbr_w_offset+ mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec_2, weights = weights)
fit.rec_3 <- lm(recovery ~ history + host_pct + rbr_w_offset+ mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec_3, weights = weights)
fit.rec_4 <- lm(recovery ~ history + host_pct + rbr_w_offset+ mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec_4, weights = weights)

#calculate multicollinearity using  vif function from car package
vif_rec_1 <- vif(fit.rec_1)
vif_rec_2 <- vif(fit.rec_2)
vif_rec_3 <- vif(fit.rec_3)
vif_rec_4 <- vif(fit.rec_4) 

# make tables
vif_table_rec_1 <- data.frame(
  feature = names(vif_rec_1), # remove intercept 
  VIF = as.numeric(vif_rec_1)
)           
vif_table_rec_2 <- data.frame(
  feature = names(vif_rec_2), # remove intercept 
  VIF = as.numeric(vif_rec_2)
)
vif_table_rec_3 <- data.frame(
  feature = names(vif_rec_3), # remove intercept 
  VIF = as.numeric(vif_rec_3)
)
vif_table_rec_4 <- data.frame(
  feature = names(vif_rec_4), # remove intercept 
  VIF = as.numeric(vif_rec_4)
)
# write to csv
write.csv(vif_table_rec_1, file.path(base_path,"/results/subgroup/vif_recovery_subgroup1.csv"), row.names = FALSE)
write.csv(vif_table_rec_2, file.path(base_path,"/results/subgroup/vif_recovery_subgroup2.csv"), row.names = FALSE)
write.csv(vif_table_rec_3, file.path(base_path,"/results/subgroup/vif_recovery_subgroup3.csv"), row.names = FALSE)
write.csv(vif_table_rec_4, file.path(base_path,"/results/subgroup/vif_recovery_subgroup4.csv"), row.names = FALSE)

# do the same for intermediate subgroup
fit.sev_intermediate <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev_intermediate, weights = weights)
vif_sev_intermediate <- vif(fit.sev_intermediate)
vif_table_intermediate <- data.frame(
  feature = names(vif_sev_intermediate), # remove intercept 
  VIF = as.numeric(vif_sev_intermediate)
)
write.csv(vif_table_intermediate, file.path(base_path,"/results/subgroup/vif_severity_intermediate.csv"), row.names = FALSE)    

fit.rec_intermediate <- lm(recovery ~ history + host_pct + rbr_w_offset+ mean_temperature + sum_precipitation_mm + mean_tri,
                data = m.data.rec_intermediate, weights = weights)
vif_rec_intermediate <- vif(fit.rec_intermediate)
vif_table_rec_intermediate <- data.frame(
  feature = names(vif_rec_intermediate), # remove intercept 
  VIF = as.numeric(vif_rec_intermediate)
)
write.csv(vif_table_rec_intermediate, file.path(base_path,"/results/subgroup/vif_recovery_intermediate.csv"), row.names = FALSE)
