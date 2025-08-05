#=============================
# FIRE INSECT CO-OCCURRENCE 3-9 YEARS ANALYSIS
#=============================


# Load required libraries
library(MatchIt)
library(dplyr)
library(marginaleffects)
library(MuMIn)
library(sensemakr)

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

# Source the file with covariate balance plots
source(file.path(base_path, "/src/covariate_balance_plots.R"))

# Load your data
source(file.path(base_path, "/src/load_data.R"))

# Check if data is loaded
if (!exists("hist_gt90_2_2")) {
  stop("Data not loaded. Make sure hist_gt90_2_2 exists.")
}

# Define parameters
methods <- c("nearest", "optimal")
distances <- c("glm", "mahalanobis")
links <- c("logit", "probit")
m_orders <- c("random")
calipers <- c(0.1, 0.2, 0.3)
replace_list <- c(FALSE)
use_mahvars_list <- c(FALSE, TRUE)

#rename df
intermediate <- hist_gt90_2_2

#============================
# SEVERITY ANALYSIS FOR HIST_GT90_2_2
#============================
cat("\n\n======= SEVERITY ANALYSIS FOR intermediate =======\n")
cat("Running model selection...\n")
best_model_info_sev <- find_best_model(intermediate, response = "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_sev <- best_model_info_sev$best_model
all_results_sev <- best_model_info_sev$results

# Check if we found a best model
if (!is.null(best_model_sev)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.sev <- match_data(best_model_sev)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.sev %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data
  cat("\nFitting linear model...\n")
  fit.sev <- lm(rbr_w_offset ~ history + host_pct + isi_90 + dc_90 + dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,
               data = m.data.sev, weights = weights)

  # check vifs 
  vif_values <- car::vif(fit.sev)

  # if any values in vif_values list is greater than 10 remove and update model
  if (any(vif_values > 10)) { 
    cat("High VIF detected. Removing variables with VIF > 10...\n")
    high_vif_vars <- names(vif_values[vif_values > 10])
    cat("Removing variables:", paste(high_vif_vars, collapse = ", "), "\n")
   # Create a formula to remove those variables
   remove_formula <- as.formula(
   paste(". ~ . -", paste(high_vif_vars, collapse = " - ")))
   # Update the model
   fit.sev <- update(fit.sev, remove_formula)
  } else {
    cat("No high VIF detected. Proceeding with original model.\n")
  }             
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.sev))

  # Calculate R-squared
  cat("\nR-squared:\n")
  fit.sev_r2 <- r.squaredGLMM(fit.sev)
  print(fit.sev_r2)
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_sev, file.path(base_path, "/results/subgroup/best_model_intermediate_severity.RDS"))
  saveRDS(fit.sev, file.path(base_path,"/results/subgroup/fit_model_intermediate_severity.RDS"))
  saveRDS(fit.sev_r2,file.path(base_path, "/results/subgroup/fit_model_intermediate_severity_r2.RDS"))
  
  # Calculate treatment effects
  cat("\nCalculating treatment effects...\n")
  avg_effects <- avg_comparisons(
    fit.sev, 
    variables = "history",
    newdata = subset(fit.sev$model, history == 1)
  )
  print(avg_effects)
  
  # Save treatment effects
  saveRDS(avg_effects,file.path(base_path, "/results/subgroup/treatment_effects_intermediate_severity.RDS"))
  
  # Generate balance plots
  cat("\nGenerating balance plots...\n")
  plot_dir <- file.path(base_path,"/plots/balance/intermediate/")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  model_list <- list("3-9 years - Severity" = best_model_sev)
  balance_plots <- generate_all_balance_plots(model_list, plot_dir)
  
  # Power analysis
  cat("\nPerforming power analysis...\n")
  library(pwr)
  
  # Extract observed effect
  model_summary <- summary(fit.sev)
  treatment_coef <- coef(model_summary)["history", ]
  observed_effect <- treatment_coef["Estimate"]
  observed_se <- treatment_coef["Std. Error"]
  
  # Calculate Cohen's d effect size
  d <- observed_effect / sd(m.data.sev$rbr_w_offset)
  cat("\nObserved effect size (Cohen's d):", round(d, 3), "\n")
  
  # Calculate power for different sample sizes
  sample_sizes <- seq(20, 300, by = 10)
  powers <- sapply(sample_sizes, function(n) {
    pwr.t.test(d = d, n = n/2, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
  })
  
  # Create power curve plot
  library(ggplot2)
  power_data <- data.frame(sample_size = sample_sizes, power = powers)
  power_plot <- ggplot(power_data, aes(x = sample_size, y = power)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Total Sample Size", 
      y = "Statistical Power",
      title = paste("Power Analysis for intermediate (d =", round(d, 2), ")")
    ) +
    theme_minimal()
  print(power_plot)
  
  # Save power plot
  power_dir <- file.path(base_path,"/plots/power_analysis/intermediate/")
  dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(power_dir, "power_curve_severity.png"), power_plot, width = 8, height = 6, dpi = 300)
  
  # Find required sample size for power = 0.8
  min_sample_size <- sample_sizes[min(which(powers >= 0.8))]
  cat("Minimum sample size required for 80% power:", min_sample_size, "\n")
  
  # Calculate minimum detectable effect size with current sample
  current_sample_size <- nrow(m.data.sev)
  min_d <- pwr.t2n.test(
    n1 = sum(m.data.sev$history == 1), 
    n2 = sum(m.data.sev$history == 0), 
    power = 0.8,
    sig.level = 0.05,
    alternative = "two.sided"
  )$d
  
  min_raw_effect <- min_d * sd(m.data.sev$rbr_w_offset)
  cat("\nMinimum detectable effect size (Cohen's d) with current sample size:", round(min_d, 3), "\n")
  cat("Minimum detectable effect size (raw units) with current sample size:", round(min_raw_effect, 3), "\n")
  
} else {
  cat("No suitable severity model found. Try different parameters.\n")
}

#============================
# RECOVERY ANALYSIS FOR intermediate
#============================
cat("\n\n======= RECOVERY ANALYSIS FOR intermediate =======\n")
cat("Running model selection...\n")
best_model_info_rec <- find_best_model(intermediate, response = "recovery", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_rec <- best_model_info_rec$best_model
all_results_rec <- best_model_info_rec$results

# Check if we found a best model
if (!is.null(best_model_rec)) {
  # Continue with analysis using best model
  cat("Getting matched data...\n")
  m.data.rec <- match_data(best_model_rec)
  
  # Print summary statistics
  cat("\nSummary of matched data:\n")
  summary_table <- m.data.rec %>% 
    group_by(history) %>% 
    summarise("Number of Fires" = n()) %>% 
    mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))
  print(summary_table)
  
  # Fit model with matched data - using recovery covariates
  cat("\nFitting linear model...\n")
  fit.rec <- lm(recovery ~ history + host_pct + rbr_w_offset + mean_temperature + sum_precipitation_mm + mean_tri,
               data = m.data.rec, weights = weights)
  
  # Print model summary
  cat("\nModel summary:\n")
  print(summary(fit.rec))
  
  # Calculate R-squared
  cat("\nR-squared:\n")
  fit.rec_r2 <- r.squaredGLMM(fit.rec)
  print(fit.rec_r2)
  
  # Save results
  cat("\nSaving results...\n")
  saveRDS(best_model_rec, file.path(base_path,"/results/subgroup/best_model_intermediate_recovery.RDS"))
  saveRDS(fit.rec,file.path(base_path, "/results/subgroup/fit_model_intermediate_recovery.RDS"))
  saveRDS(fit.rec_r2, file.path(base_path,"/results/subgroup/fit_model_intermediate_recovery_r2.RDS"))
  
  # Calculate treatment effects
  cat("\nCalculating treatment effects...\n")
  avg_effects_rec <- avg_comparisons(
    fit.rec, 
    variables = "history",
    newdata = subset(fit.rec$model, history == 1)
  )
  print(avg_effects_rec)
  
  # Save treatment effects
  saveRDS(avg_effects_rec, file.path(base_path,"/results/subgroup/treatment_effects_intermediate_recovery.RDS"))
  
  # Generate balance plots
  cat("\nGenerating balance plots...\n")
  plot_dir <- file.path(base_path,"/plots/balance/intermediate/")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  model_list <- list("3-9 years - Recovery" = best_model_rec)
  balance_plots <- generate_all_balance_plots(model_list, plot_dir)
  
  # Power analysis
  cat("\nPerforming power analysis...\n")
  
  # Extract observed effect
  model_summary <- summary(fit.rec)
  treatment_coef <- coef(model_summary)["history", ]
  observed_effect <- treatment_coef["Estimate"]
  observed_se <- treatment_coef["Std. Error"]
  
  # Calculate Cohen's d effect size
  d <- observed_effect / sd(m.data.rec$recovery)
  cat("\nObserved effect size (Cohen's d):", round(d, 3), "\n")
  
  # Calculate power for different sample sizes
  sample_sizes <- seq(20, 300, by = 10)
  powers <- sapply(sample_sizes, function(n) {
    pwr.t.test(d = d, n = n/2, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
  })
  
  # Create power curve plot
  power_data <- data.frame(sample_size = sample_sizes, power = powers)
  power_plot_rec <- ggplot(power_data, aes(x = sample_size, y = power)) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Total Sample Size", 
      y = "Statistical Power",
      title = paste("Power Analysis for intermediate Recovery (d =", round(d, 2), ")")
    ) +
    theme_minimal()
  print(power_plot_rec)
  
  # Save power plot
  power_dir <- file.path(base_path,"/plots/power_analysis/intermediate/")
  dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(power_dir, "power_curve_recovery.png"), power_plot_rec, width = 8, height = 6, dpi = 300)
  
  # Find required sample size for power = 0.8
  min_sample_size <- sample_sizes[min(which(powers >= 0.8))]
  cat("Minimum sample size required for 80% power:", min_sample_size, "\n")
  
  # Calculate minimum detectable effect size with current sample
  current_sample_size <- nrow(m.data.rec)
  min_d <- pwr.t2n.test(
    n1 = sum(m.data.rec$history == 1), 
    n2 = sum(m.data.rec$history == 0), 
    power = 0.8,
    sig.level = 0.05,
    alternative = "two.sided"
  )$d
  
  min_raw_effect <- min_d * sd(m.data.rec$recovery)
  cat("\nMinimum detectable effect size (Cohen's d) with current sample size:", round(min_d, 3), "\n")
  cat("Minimum detectable effect size (raw units) with current sample size:", round(min_raw_effect, 3), "\n")
  
} else {
  cat("No suitable recovery model found. Try different parameters.\n")
}

#============================
# SUMMARY
#============================
cat("\n\n======= SUMMARY OF ANALYSIS FOR intermediate =======\n")

# Create a summary dataframe
summary_df <- data.frame(
  Analysis = c("Severity", "Recovery"),
  Model_Found = c(!is.null(best_model_sev), !is.null(best_model_rec)),
  stringsAsFactors = FALSE
)

# Add sample sizes if models were found
if (!is.null(best_model_sev)) {
  summary_df$Total_Sample_Size[summary_df$Analysis == "Severity"] <- nrow(m.data.sev)
  summary_df$Treated_Sample[summary_df$Analysis == "Severity"] <- sum(m.data.sev$history == 1)
  summary_df$Control_Sample[summary_df$Analysis == "Severity"] <- sum(m.data.sev$history == 0)
}

if (!is.null(best_model_rec)) {
  summary_df$Total_Sample_Size[summary_df$Analysis == "Recovery"] <- nrow(m.data.rec)
  summary_df$Treated_Sample[summary_df$Analysis == "Recovery"] <- sum(m.data.rec$history == 1)
  summary_df$Control_Sample[summary_df$Analysis == "Recovery"] <- sum(m.data.rec$history == 0)
}

# Print summary
print(summary_df)

# Save summary
write.csv(summary_df, file.path(base_path,"results/subgroup/intermediate_analysis_summary.csv"), row.names = FALSE)

cat("\nAnalysis complete. Results saved to results directory.\n")

#============================
# GENERATE TREATMENT EFFECT PLOTS
#============================
cat("\n\n======= GENERATING TREATMENT EFFECT PLOTS =======\n")

# Source dedicated intermediate treatment effects plotting script
source(file.path(base_path, "/src/treatment_effect_plots.R"))
generate_all_intermediate_treatment_effect_plots()
generate_all_combined_intermediate_treatment_effect_plots()


#============================
# GENERATE TREATMENT EFFECT TABLES
#============================
cat("\n\n======= GENERATING TREATMENT EFFECT TABLES =======\n")

# Source dedicated intermediate treatment effects export script
source(file.path(base_path,"/src/export_treatment_effects.R"))
export_intermediate_treatment_effects()

cat("\nAnalysis complete. Results saved to results directory.\n")



#============================
# GENERATE SENSITIVITY ANALYSIS

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

# intermediate severity and recovery models
fit_int_sev <- safe_read_rds(paste0(result_dir, "fit_model_intermediate_severity.RDS"))
fit_int_rec <- safe_read_rds(paste0(result_dir, "fit_model_intermediate_recovery.RDS"))

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

# set treatment variable
treatment_var <- "history"

# Create a list of all models
all_models <- list(
  "intermediate severity" = fit_int_sev,
  "intermediate recovery" = fit_int_rec
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
sev_stats <- as.data.frame(all_sens_results[["intermediate severity"]]$sensitivity_stats)
rec_stats <- as.data.frame(all_sens_results[["intermediate recovery"]]$sensitivity_stats)


all_stats <- rbind(
  sev_stats %>% mutate(Subgroup = "Intermediate", Response = "Severity"),
  rec_stats %>% mutate(Subgroup = "Intermediate", Response = "Recovery")
  )  %>% relocate(
    Subgroup,
    Response,
    .before = "treatment")


# Save all stats to a CSV file
output_stats_file <- file.path(result_dir, "intermediate_sensitivity_analysis_stats.csv")
write.csv(all_stats, output_stats_file, row.names = FALSE)
