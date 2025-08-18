library(ggplot2)
library(tidyverse)


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

# read in all subgroup fits
fit_sev_1 <- readRDS(paste0(base_path, "/results/subgroup/fit_model_subgroup1_severity.RDS"))
fit_sev_2 <- readRDS(paste0(base_path, "/results/subgroup/fit_model_subgroup2_severity.RDS"))
fit_sev_3 <- readRDS(paste0(base_path, "/results/subgroup/fit_model_subgroup3_severity.RDS"))
fit_sev_4 <- readRDS(paste0(base_path, "/results/subgroup/fit_model_subgroup4_severity.RDS"))
fit_rec_1 <- readRDS(paste0(base_path, "/results/subgroup/fit_model_subgroup1_recovery.RDS"))
fit_rec_2 <- readRDS(paste0(base_path, "/results/subgroup/fit_model_subgroup2_recovery.RDS"))
fit_rec_3 <- readRDS(paste0(base_path, "/results/subgroup/fit_model_subgroup3_recovery.RDS"))
fit_rec_4 <- readRDS(paste0(base_path, "/results/subgroup/fit_model_subgroup4_recovery.RDS"))    

# Convert summary statistics to a data frame with consistent scientific notation
sev_model_stats <- function(model_summary) {
  # Extract coefficient table
  coef_table <- as.data.frame(model_summary$coefficients)
  
  # Add term names as a column
  coef_table$term <- rownames(coef_table)
  
  # Rename columns for clarity
  names(coef_table) <- c("estimate", "std_error", "t_value", "p_value", "term")
  
  # Reorder columns
  coef_table <- coef_table[, c("term", "estimate", "std_error", "t_value", "p_value")]
  
  # Sort: intercept first, then by significance (p-value)
  intercept_row <- coef_table[coef_table$term == "(Intercept)", ]
  other_rows <- coef_table[coef_table$term != "(Intercept)", ]
  other_rows <- other_rows[order(other_rows$p_value), ]
  
  # Combine rows back together
  coef_table <- rbind(intercept_row, other_rows)
  
  # Add significance stars
  coef_table$significance <- ""
  coef_table$significance[coef_table$p_value < 0.05] <- "*"
  coef_table$significance[coef_table$p_value < 0.01] <- "**"
  coef_table$significance[coef_table$p_value < 0.001] <- "***"
  
  # Format numeric columns with scientific notation
  coef_table$estimate <- formatC(coef_table$estimate, format = "e", digits = 4)
  coef_table$std_error <- formatC(coef_table$std_error, format = "e", digits = 4)
  coef_table$t_value <- formatC(coef_table$t_value, format = "e", digits = 4)
  coef_table$p_value <- formatC(coef_table$p_value, format = "e", digits = 4)
  
  # Return the data frame
  return(coef_table)
}

# Create the data frame from severity model summary
sev_stats_1 <- sev_model_stats(summary(fit_sev_1))
sev_stats_2 <- sev_model_stats(summary(fit_sev_2))
sev_stats_3 <- sev_model_stats(summary(fit_sev_3))
sev_stats_4 <- sev_model_stats(summary(fit_sev_4))      

#create the data frame from recovery model summary
rec_stats_1 <- sev_model_stats(summary(fit_rec_1))
rec_stats_2 <- sev_model_stats(summary(fit_rec_2))
rec_stats_3 <- sev_model_stats(summary(fit_rec_3))
rec_stats_4 <- sev_model_stats(summary(fit_rec_4))

# Print the results
print(sev_stats_1)
print(sev_stats_2)
print(sev_stats_3)
print(sev_stats_4)
print(rec_stats_1)
print(rec_stats_2)
print(rec_stats_3)
print(rec_stats_4)

# Save to CSV files
write.csv(sev_stats_1, file.path(base_path, "results/subgroup/severity_model_statistics_subgroup1.csv"), 
          row.names = FALSE)
write.csv(sev_stats_2, file.path(base_path, "results/subgroup/severity_model_statistics_subgroup2.csv"), 
          row.names = FALSE)                    
write.csv(sev_stats_3, file.path(base_path, "results/subgroup/severity_model_statistics_subgroup3.csv"), 
          row.names = FALSE)
write.csv(sev_stats_4, file.path(base_path, "results/subgroup/severity_model_statistics_subgroup4.csv"), 
          row.names = FALSE)
write.csv(rec_stats_1, file.path(base_path, "results/subgroup/recovery_model_statistics_subgroup1.csv"), 
          row.names = FALSE)
write.csv(rec_stats_2, file.path(base_path, "results/subgroup/recovery_model_statistics_subgroup2.csv"), 
          row.names = FALSE)
write.csv(rec_stats_3, file.path(base_path, "results/subgroup/recovery_model_statistics_subgroup3.csv"), 
          row.names = FALSE)
write.csv(rec_stats_4, file.path(base_path, "results/subgroup/recovery_model_statistics_subgroup4.csv"), 
          row.names = FALSE)


# get sample sizes for each subgroup
get_sample_size <- function(model) {
    if (is.null(model)) {
        return(0)
    }
    return(length(model$residuals))
}   


# Sample sizes for each subgroup
sample_size_sev_1 <- get_sample_size(fit_sev_1)
sample_size_sev_2 <- get_sample_size(fit_sev_2)
sample_size_sev_3 <- get_sample_size(fit_sev_3)
sample_size_sev_4 <- get_sample_size(fit_sev_4)
sample_size_rec_1 <- get_sample_size(fit_rec_1)
sample_size_rec_2 <- get_sample_size(fit_rec_2)
sample_size_rec_3 <- get_sample_size(fit_rec_3)
sample_size_rec_4 <- get_sample_size(fit_rec_4)

# Create a data frame with sample sizes
sample_sizes <- data.frame(
    Subgroup = c("Subgroup 1", "Subgroup 2", "Sub
group 3", "Subgroup 4"),
    Sample_Size_Severity = c(sample_size_sev_1, sample_size_sev_2, sample_size_sev_3, sample_size_sev_4),
    Sample_Size_Recovery = c(sample_size_rec_1, sample_size_rec_2, sample_size_rec_3, sample_size_rec_4)
)   
# Print the sample sizes
print(sample_sizes)

# Save sample sizes to CSV
write.csv(sample_sizes, file.path(base_path, "results/subgroup/sample_sizes_subgroups.csv"), 
          row.names = FALSE)

# get data for intermediate
fit_sev_int <- readRDS(paste0(base_path, "/results/subgroup/fit_model_intermediate_severity.RDS"))
fit_rec_int <- readRDS(paste0(base_path, "/results/subgroup/fit_model_intermediate_recovery.RDS"))

#get sample sizes for intermediate
sample_size_sev_int <- get_sample_size(fit_sev_int)
sample_size_rec_int <- get_sample_size(fit_rec_int)

# add to sample sizes data frame with intermediate in subgroup column
sample_sizes <- rbind(sample_sizes,
                        data.frame(Subgroup = "Intermediate",
                                   Sample_Size_Severity = sample_size_sev_int,
                                   Sample_Size_Recovery = sample_size_rec_int))     
# Print the updated sample sizes
print(sample_sizes)


# get sample sizes for each level of history
get_sample_size_by_history <- function(model, history_value) {
  if (is.null(model)) {
    return(0)
  }
  
  # Extract history column from model data
  history_col <- model$model$history  # Using model$model instead of model$data
  
  # Convert both the column and the input value to character for consistent comparison
  history_col <- as.character(history_col)
  history_value <- as.character(history_value)
  
  # Count matches
  return(sum(history_col == history_value))
}
# Sample sizes for each history level in severity models
sample_size_sev_1_hist_0 <- get_sample_size_by_history(fit_sev_1, as.factor(0))
sample_size_sev_1_hist_1 <- get_sample_size_by_history(fit_sev_1, 1)
sample_size_sev_2_hist_0 <- get_sample_size_by_history(fit_sev_2, 0)
sample_size_sev_2_hist_1 <- get_sample_size_by_history(fit_sev_2, 1)
sample_size_sev_3_hist_0 <- get_sample_size_by_history(fit_sev_3, 0)
sample_size_sev_3_hist_1 <- get_sample_size_by_history(fit_sev_3, 1)
sample_size_sev_4_hist_0 <- get_sample_size_by_history(fit_sev_4, 0)
sample_size_sev_4_hist_1 <- get_sample_size_by_history(fit_sev_4, 1)
# Sample sizes for each history level in recovery models
sample_size_rec_1_hist_0 <- get_sample_size_by_history(fit_rec_1, 0)
sample_size_rec_1_hist_1 <- get_sample_size_by_history(fit_rec_1, 1)
sample_size_rec_2_hist_0 <- get_sample_size_by_history(fit_rec_2, 0)
sample_size_rec_2_hist_1 <- get_sample_size_by_history(fit_rec_2, 1)
sample_size_rec_3_hist_0 <- get_sample_size_by_history(fit_rec_3, 0)
sample_size_rec_3_hist_1 <- get_sample_size_by_history(fit_rec_3, 1)
sample_size_rec_4_hist_0 <- get_sample_size_by_history(fit_rec_4, 0)
sample_size_rec_4_hist_1 <- get_sample_size_by_history(fit_rec_4, 1)    

# get for intermediate
sample_size_sev_int_hist_0 <- get_sample_size_by_history(fit_sev_int, 0)
sample_size_sev_int_hist_1 <- get_sample_size_by_history(fit_sev_int, 1)
sample_size_rec_int_hist_0 <- get_sample_size_by_history(fit_rec_int, 0)
sample_size_rec_int_hist_1 <- get_sample_size_by_history(fit_rec_int, 1)    

# create a dataframe with samplesizes by history for eeach subgroup and intermediate and add as columns to samples_sizes
sample_sizes <- sample_sizes %>% 
  mutate(
    Sample_Size_Severity_History_0 = c(sample_size_sev_1_hist_0, sample_size_sev_2_hist_0, 
                                        sample_size_sev_3_hist_0, sample_size_sev_4_hist_0, 
                                        sample_size_sev_int_hist_0),
    Sample_Size_Severity_History_1 = c(sample_size_sev_1_hist_1, sample_size_sev_2_hist_1, 
                                        sample_size_sev_3_hist_1, sample_size_sev_4_hist_1, 
                                        sample_size_sev_int_hist_1),
    Sample_Size_Recovery_History_0 = c(sample_size_rec_1_hist_0, sample_size_rec_2_hist_0, 
                                        sample_size_rec_3_hist_0, sample_size_rec_4_hist_0, 
                                        sample_size_rec_int_hist_0),
    Sample_Size_Recovery_History_1 = c(sample_size_rec_1_hist_1, sample_size_rec_2_hist_1, 
                                        sample_size_rec_3_hist_1, sample_size_rec_4_hist_1, 
                                        sample_size_rec_int_hist_1)
  ) 
# Print the updated sample sizes with history
print(sample_sizes)

# save the updated sample sizes to CSV
write.csv(sample_sizes, file.path(base_path, "results/subgroup/sample_sizes_subgroups_with_history.csv"), 
          row.names = FALSE)
