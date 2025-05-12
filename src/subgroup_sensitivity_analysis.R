# Load required libraries
library(sensemakr)
library(dplyr)
library(ggplot2)

# Define output directory
result_dir <- "/home/goldma34/sbw-wildfire-impact-recovery/results/subgroup/"

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

# Define covariates for each model type
covariates_sev <- c("host_pct", "isi_90", "dc_90", "dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri")
covariates_rec <- c("host_pct", "mean_temperature", "sum_precipitation_mm", "mean_tri")

# Define benchmark covariates (variables that are likely to be important confounders)
benchmark_covariates_sev <- c("host_pct", "fwi_90")  # Select key variables as benchmarks
benchmark_covariates_rec <- c("host_pct", "mean_temperature")

# Function to analyze one model and return results using sensemakr
analyze_one_model_sensemakr <- function(model, model_name, treatment_var, covariates, benchmark_covs) {
  if (is.null(model)) {
    cat("Model is NULL for", model_name, "\n")
    return(NULL)
  }
  
  cat("Analyzing", model_name, "...\n")
  
  # Try to run sensemakr analysis
  tryCatch({
    # Run sensemakr analysis
    sens_result <- sensemakr(
      model = model,
      treatment = treatment_var,
      benchmark_covariates = benchmark_covs,
      kd = c(1, 2, 3),    # Multiples of the benchmark for confounder strength
      ky = c(1, 2, 3),    # Multiples of the benchmark for confounder strength
      q = 1,             # Proportion of treatment effect 
      alpha = 0.05       # Significance level
    )
    
    # Extract summary information
    summary_stats <- summary(sens_result)
    
    # Extract treatment effect information
    treatment_effect <- sens_result$estimate
    std_error <- sens_result$se
    t_value <- sens_result$t_statistic
    p_value <- 2 * pt(-abs(t_value), sens_result$dof)
    
    # Extract robustness values
    rv_qa <- sens_result$sensitivity_stats$rv_qa[1]  # RV for bringing effect to zero
    rv_qa_alpha <- sens_result$sensitivity_stats$rv_qa_alpha[1]  # RV for bringing to insignificance
    
    # Extract partial R² values
    partial_r2_t_x <- sens_result$sensitivity_stats$r2yd.x[1]
    partial_r2_t_y <- sens_result$sensitivity_stats$r2yd.x[1]
    
    # Create result row
    result <- data.frame(
      Model = model_name,
      Treatment_Effect = treatment_effect,
      Standard_Error = std_error,
      T_Value = t_value,
      P_Value = p_value,
      Partial_R2_Treatment = partial_r2_t_x,
      RV_q1 = rv_qa,             # Robustness value for bringing effect to zero
      RV_q1_alpha = rv_qa_alpha, # Robustness value for bringing to insignificance
      Effect_Direction = ifelse(treatment_effect > 0, "Positive", "Negative"),
      Significance = ifelse(p_value < 0.05, "Significant", "Not Significant"),
      stringsAsFactors = FALSE
    )
    
    # Add robustness ratio and assessment
    typical_r2 <- 0.15  # A moderate benchmark partial R2
    result$Robustness_Ratio <- result$RV_q1 / typical_r2
    
    result$Robustness_Assessment <- case_when(
      result$Robustness_Ratio > 2 ~ "Very Strong",
      result$Robustness_Ratio > 1 ~ "Strong",
      result$Robustness_Ratio > 0.5 ~ "Moderate",
      result$Robustness_Ratio > 0.25 ~ "Weak",
      TRUE ~ "Very Weak"
    )
    
    # Return both the result data and the sensemakr object
    list(result = result, sens_object = sens_result)
    
  }, error = function(e) {
    cat("Error analyzing", model_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Analyze all models
cat("Starting sensitivity analysis using sensemakr...\n")

# Create list of all models, names, etc.
models <- list(fit_sev_1, fit_sev_2, fit_sev_3, fit_sev_4,
              fit_rec_1, fit_rec_2, fit_rec_3, fit_rec_4)

model_names <- c("Subgroup 1 (0-2 years) - Severity",
                "Subgroup 2 (3-5 years) - Severity",
                "Subgroup 3 (6-9 years) - Severity",
                "Subgroup 4 (10+ years) - Severity",
                "Subgroup 1 (0-2 years) - Recovery",
                "Subgroup 2 (3-5 years) - Recovery",
                "Subgroup 3 (6-9 years) - Recovery",
                "Subgroup 4 (10+ years) - Recovery")

treatment_var <- "history"

cov_list <- list(
  covariates_sev, covariates_sev, covariates_sev, covariates_sev,
  covariates_rec, covariates_rec, covariates_rec, covariates_rec
)

benchmark_list <- list(
  benchmark_covariates_sev, benchmark_covariates_sev, benchmark_covariates_sev, benchmark_covariates_sev,
  benchmark_covariates_rec, benchmark_covariates_rec, benchmark_covariates_rec, benchmark_covariates_rec
)

# Run analysis for each model and collect results
all_results <- list()
all_sens_objects <- list()

for (i in 1:length(models)) {
  result <- analyze_one_model_sensemakr(
    models[[i]], 
    model_names[i], 
    treatment_var, 
    cov_list[[i]], 
    benchmark_list[[i]]
  )
  
  if (!is.null(result)) {
    all_results[[i]] <- result$result
    all_sens_objects[[i]] <- result$sens_object
    names(all_sens_objects)[i] <- model_names[i]
  } else {
    all_results[[i]] <- NULL
    all_sens_objects[[i]] <- NULL
  }
}

# Combine all non-null results
valid_results <- all_results[!sapply(all_results, is.null)]
valid_sens_objects <- all_sens_objects[!sapply(all_sens_objects, is.null)]

if (length(valid_results) > 0) {
  combined_results <- do.call(rbind, valid_results)
  
  # Save full results
  write.csv(combined_results, paste0(result_dir, "sensitivity_analysis_sensemakr.csv"), row.names = FALSE)
  
  # Create interpretable summary
  interpretable_summary <- combined_results %>%
    select(Model, Treatment_Effect, P_Value, Partial_R2_Treatment, RV_q1, Robustness_Assessment) %>%
    mutate(
      Effect_Size = round(Treatment_Effect, 2),
      P_Value = round(P_Value, 4),
      Partial_R2 = round(Partial_R2_Treatment, 3),
      Robustness_Value = round(RV_q1, 3)
    ) %>%
    select(Model, Effect_Size, P_Value, Partial_R2, Robustness_Value, Robustness_Assessment)
  
  # Save interpretable summary
  write.csv(interpretable_summary, paste0(result_dir, "sensitivity_analysis_interpretable_sensemakr.csv"), row.names = FALSE)
  
  # Display results
  cat("\nSensitivity Analysis Results:\n")
  print(interpretable_summary)
  
  cat("\nSensitivity analysis complete. Results saved to:",
      paste0(result_dir, "sensitivity_analysis_sensemakr.csv"), "and",
      paste0(result_dir, "sensitivity_analysis_interpretable_sensemakr.csv"), "\n")
  
  # Create and save contour plots for each model
  cat("\nGenerating contour plots for each model...\n")
  
  for (i in seq_along(valid_sens_objects)) {
    model_name <- names(valid_sens_objects)[i]
    sens_obj <- valid_sens_objects[[i]]
    
    # Clean model name for file naming
    clean_name <- gsub(" |\\(|\\)|\\+|-", "_", model_name)
    
    # Create plot
    plot_file <- paste0(result_dir, "contour_plot_", clean_name, ".png")
    
    # Generate and save contour plot
    png(plot_file, width = 8, height = 6, units = "in", res = 300)
    plot(sens_obj, 
         type = "contour", 
         main = paste("Sensitivity Analysis for", model_name),
         xlab = "Partial R² of confounder(s) with treatment",
         ylab = "Partial R² of confounder(s) with outcome")
    dev.off()
    
    cat("Saved contour plot for", model_name, "to", plot_file, "\n")
    
    # Generate and save extreme scenarios plot
    plot_file2 <- paste0(result_dir, "extreme_scenarios_", clean_name, ".png")
    
    png(plot_file2, width = 8, height = 6, units = "in", res = 300)
    plot(sens_obj, 
         type = "extreme", 
         main = paste("Extreme Scenarios for", model_name))
    dev.off()
    
    cat("Saved extreme scenarios plot for", model_name, "to", plot_file2, "\n")
  }
  
  # Generate combined plots for severity and recovery models
  # Group objects by severity and recovery
  severity_objects <- valid_sens_objects[grep("Severity", names(valid_sens_objects))]
  recovery_objects <- valid_sens_objects[grep("Recovery", names(valid_sens_objects))]
  
  # Create comparison tables
  if (length(severity_objects) > 0) {
    severity_table <- ovb_bounds(severity_objects)
    write.csv(severity_table, paste0(result_dir, "sensitivity_bounds_severity.csv"), row.names = TRUE)
  }
  
  if (length(recovery_objects) > 0) {
    recovery_table <- ovb_bounds(recovery_objects)
    write.csv(recovery_table, paste0(result_dir, "sensitivity_bounds_recovery.csv"), row.names = TRUE)
  }
  
} else {
  cat("\nNo valid sensitivity analysis results were generated.\n")
}

# Save the sensitivity objects for later use
if (length(valid_sens_objects) > 0) {
  saveRDS(valid_sens_objects, paste0(result_dir, "sensitivity_objects.RDS"))
  cat("\nSensitivity objects saved to:", paste0(result_dir, "sensitivity_objects.RDS"), "\n")
}

cat("\nSensitivity analysis complete.\n")