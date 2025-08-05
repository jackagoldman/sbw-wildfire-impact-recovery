# Load required libraries
library(sensemakr)
library(dplyr)
library(ggplot2)

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

# Define output directory
result_dir <- file.path(base_path, "results/subgroup/")

# Function to safely read RDS files
safe_read_rds <- function(file_path) {
  if (file.exists(file_path)) {
    tryCatch({
      readRDS(file_path)
    }, error = function(e) {
      cat("Error reading", file_path, ":", e$message, "\n")
      NULL
    })
  } else {
    cat("File not found:", file_path, "\n")
    NULL
  }
}

# Load all saved model fits
cat("Loading models...\n")
fit_sev_1 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup1_severity.RDS"))
fit_sev_2 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup2_severity.RDS"))
fit_sev_3 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup3_severity.RDS"))
fit_sev_4 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup4_severity.RDS"))

fit_rec_1 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup1_recovery.RDS"))
fit_rec_2 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup2_recovery.RDS"))
fit_rec_3 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup3_recovery.RDS"))
fit_rec_4 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup4_recovery.RDS"))

# set treatment variable
treatment_var <- "history"

# get benchmark covariates based on covariates used in each model
get_model_covariates <- function(model, treatment_var = "history", 
                                exclude_treatment = TRUE, include_details = FALSE) {
  # Check if model is NULL
  if (is.null(model)) {
    cat("Error: Model is NULL\n")
    return(NULL)
  }
  
  # Extract coefficient names (exclude intercept)
  coef_names <- names(coef(model))
  coef_names <- setdiff(coef_names, "(Intercept)")
  
  # Exclude treatment variable if requested
  if (exclude_treatment && treatment_var %in% coef_names) {
    coef_names <- setdiff(coef_names, treatment_var)
  }
  
  # If no details requested, just return the covariate names
  if (!include_details) {
    return(coef_names)
  }
  
  # Get model summary for detailed information
  model_summary <- summary(model)
  coef_table <- coef(model_summary)
  
  # Extract detailed information for each covariate
  covariate_details <- list(
    names = coef_names,
    coefficients = coef(model)[coef_names],
    p_values = coef_table[coef_names, "Pr(>|t|)"],
    t_values = coef_table[coef_names, "t value"],
    std_errors = coef_table[coef_names, "Std. Error"]
  )
  
  # Try to calculate VIF values if possible
  tryCatch({
    if (requireNamespace("car", quietly = TRUE)) {
      vif_values <- car::vif(model)
      covariate_details$vif <- vif_values[names(vif_values) %in% coef_names]
    }
  }, error = function(e) {
    cat("Note: Could not calculate VIF values:", e$message, "\n")
  })
  
  # Sort covariates by significance (absolute t-value)
  if (length(coef_names) > 0) {
    importance_order <- order(abs(coef_table[coef_names, "t value"]), decreasing = TRUE)
    covariate_details$by_importance <- coef_names[importance_order]
  }
  
  return(covariate_details)
}

# Create a list of all models
all_models <- list(
  "Subgroup 1 (0-2 years) - Severity" = fit_sev_1,
  "Subgroup 2 (3-5 years) - Severity" = fit_sev_2,
  "Subgroup 3 (6-9 years) - Severity" = fit_sev_3, 
  "Subgroup 4 (10+ years) - Severity" = fit_sev_4,
  "Subgroup 1 (0-2 years) - Recovery" = fit_rec_1,
  "Subgroup 2 (3-5 years) - Recovery" = fit_rec_2,
  "Subgroup 3 (6-9 years) - Recovery" = fit_rec_3,
  "Subgroup 4 (10+ years) - Recovery" = fit_rec_4
)

# Get covariates for each model
all_covariates <- lapply(names(all_models), function(model_name) {
  model <- all_models[[model_name]]
  if (!is.null(model)) {
    covs <- get_model_covariates(model)
    cat(model_name, "covariates:", paste(covs, collapse=", "), "\n")
    return(covs)
  } else {
    cat(model_name, "is NULL\n")
    return(NULL)
  }
})
names(all_covariates) <- names(all_models)


# Use the covariates for sensitivity analysis
all_sens_results <- list()  # Initialize list to store all sensitivity results

for (model_name in names(all_models)) {
  model <- all_models[[model_name]]
  covs <- all_covariates[[model_name]]
  
  if (!is.null(model) && !is.null(covs) && length(covs) > 0) {
    cat("\nRunning sensitivity analysis for", model_name, "\n")
    
    # Run the sensitivity analysis with model-specific covariates
    sens_result <- tryCatch({
      sensemakr(
        model = model,
        treatment = treatment_var,
        benchmark_covariates = covs,
        q = 1,
        alpha = 0.05
      )
    }, error = function(e) {
      cat("Error in sensitivity analysis:", e$message, "\n")
      NULL
    })
    
    # Store the result in the list (even if NULL)
    all_sens_results[[model_name]] <- sens_result
    
    # Process or save the individual results as needed 
    if (!is.null(sens_result)) {
      # Print sensitivity analysis results
      cat("Sensitivity analysis result for", model_name, ":\n")
      print(sens_result)
      
      # Optionally save individual result to a file
      output_file <- paste0(result_dir, "sensitivity_analysis_", gsub(" ", "_", model_name), ".RDS")
      saveRDS(sens_result, output_file)
      cat("Saved individual sensitivity analysis result to", output_file, "\n")
    } else {
      cat("No valid sensitivity analysis result for", model_name, "\n")
    }  
  }
}

all_sens_results

#  access all models stats
sev1_stats <- as.data.frame(all_sens_results[["Subgroup 1 (0-2 years) - Severity"]]$sensitivity_stats)
sev2_stats <- as.data.frame(all_sens_results[["Subgroup 2 (3-5 years) - Severity"]]$sensitivity_stats)
sev3_stats <- as.data.frame(all_sens_results[["Subgroup 3 (6-9 years) - Severity"]]$sensitivity_stats)
sev4_stats <- as.data.frame(all_sens_results[["Subgroup 4 (10+ years) - Severity"]]$sensitivity_stats)
rec1_stats <- as.data.frame(all_sens_results[["Subgroup 1 (0-2 years) - Recovery"]]$sensitivity_stats)
rec2_stats <- as.data.frame(all_sens_results[["Subgroup 2 (3-5 years) - Recovery"]]$sensitivity_stats)
rec3_stats <- as.data.frame(all_sens_results[["Subgroup 3 (6-9 years) - Recovery"]]$sensitivity_stats)
rec4_stats <- as.data.frame(all_sens_results[["Subgroup 4 (10+ years) - Recovery"]]$sensitivity_stats)
# Combine all stats into a single data frame
all_stats <- rbind(
  sev1_stats %>% mutate(Subgroup = "0-2 years", Response = "Severity"),
  sev2_stats %>% mutate(Subgroup = "3-5 years", Response = "Severity"),
  sev3_stats %>% mutate(Subgroup = "6-9 years", Response = "Severity"),
  sev4_stats %>% mutate(Subgroup = "10+ years", Response = "Severity"),
  rec1_stats %>% mutate(Subgroup = "0-2 years", Response = "Recovery"),
  rec2_stats %>% mutate(Subgroup = "3-5 years", Response = "Recovery"),
  rec3_stats %>% mutate(Subgroup = "6-9 years", Response = "Recovery"),
  rec4_stats %>% mutate(Subgroup = "10+ years", Response = "Recovery")
  )  %>% relocate(
    Subgroup,
    Response,
    .before = "treatment")

# Save all stats to a CSV file
output_stats_file <- file.path(result_dir, "sensitivity_analysis_stats.csv")
write.csv(all_stats, output_stats_file, row.names = FALSE)
