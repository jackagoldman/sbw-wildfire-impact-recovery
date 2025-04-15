#' Calculate Minimum Detectable Effect Size for Propensity Score Matched Data
#'
#' This function calculates the minimum effect size that can be detected with
#' 80% power given your current sample size or other sample sizes of interest.
#' 
#' @param model_fit A fitted model object (from lm)
#' @param matched_data The matched data used in the model
#' @param power_target The power level to target (default 0.8)
#' @param sig_level The significance level (default 0.05)
#' @param sample_sizes Optional vector of sample sizes to evaluate
#' @param plot_filename Optional filename to save the plot
#' @return A list with minimum detectable effect sizes and a plot
min_detectable_effect_size <- function(model_fit, matched_data, 
                                       power_target = 0.8, 
                                       sig_level = 0.05,
                                       sample_sizes = NULL,
                                       plot_filename = NULL) {
  
  require(pwr)
  require(ggplot2)
  require(dplyr)
  
  # Get current sample size from matched data
  current_sample_size <- nrow(matched_data)
  cat("Current total sample size:", current_sample_size, "\n")
  
  # Get treated and control sample sizes
  treated_size <- sum(matched_data$history == 1)
  control_size <- sum(matched_data$history == 0)
  cat("Treatment group:", treated_size, "fires\n")
  cat("Control group:", control_size, "fires\n")
  
  # Calculate Cohen's d for the model
  model_summary <- summary(model_fit)
  treatment_coef <- coef(model_summary)["history", ]
  observed_effect <- treatment_coef["Estimate"]
  observed_se <- treatment_coef["Std. Error"]
  
  outcome_var <- all.vars(formula(model_fit))[1]
  outcome_sd <- sd(matched_data[[outcome_var]])
  cohens_d <- observed_effect / outcome_sd
  
  cat("Observed effect size (raw):", round(observed_effect, 3), "\n")
  cat("Outcome standard deviation:", round(outcome_sd, 3), "\n")
  cat("Observed Cohen's d:", round(cohens_d, 3), "\n")
  
  # Calculate minimum detectable effect size for current sample size
  min_d <- pwr.t2n.test(n1 = treated_size, 
                        n2 = control_size, 
                        power = power_target,
                        sig.level = sig_level,
                        alternative = "two.sided")$d
  
  min_raw <- min_d * outcome_sd
  
  cat("Minimum detectable effect size (Cohen's d):", round(min_d, 3), "\n")
  cat("Minimum detectable effect size (raw units):", round(min_raw, 3), "\n")
  
  # If sample sizes are provided, create a plot of effect size vs. sample size
  if (!is.null(sample_sizes)) {
    # Calculate minimum detectable effect sizes for various sample sizes
    effect_sizes <- sapply(sample_sizes, function(n) {
      # Assume the same treated/control ratio as in the current data
      n1 <- round(n * (treated_size / current_sample_size))
      n2 <- n - n1
      pwr.t2n.test(n1 = n1, n2 = n2, power = power_target, 
                   sig.level = sig_level, alternative = "two.sided")$d
    })
    
    # Convert to raw units
    raw_effects <- effect_sizes * outcome_sd
    
    # Create data frame for plotting
    plot_data <- data.frame(
      sample_size = sample_sizes,
      cohens_d = effect_sizes,
      raw_effect = raw_effects
    )
    
    # Mark current sample size
    current_line <- geom_vline(xintercept = current_sample_size, 
                               linetype = "dashed", color = "red")
    
    # Create plot
    p1 <- ggplot(plot_data, aes(x = sample_size, y = cohens_d)) +
      geom_line(color = "blue", size = 1) +
      current_line +
      labs(
        x = "Total Sample Size", 
        y = "Minimum Detectable Effect Size (Cohen's d)",
        title = paste0("Minimum Effect Size Detectable with ", 
                       power_target*100, "% Power"),
        subtitle = paste0("Current sample size = ", current_sample_size,
                          ", Minimum d = ", round(min_d, 3))
      ) +
      theme_minimal() +
      geom_text(aes(x = current_sample_size + 10, 
                    y = max(cohens_d) * 0.9,
                    label = paste("Current n =", current_sample_size)),
                hjust = 0, vjust = 1)
    
    # Create plot with raw effect units
    p2 <- ggplot(plot_data, aes(x = sample_size, y = raw_effect)) +
      geom_line(color = "darkgreen", size = 1) +
      current_line +
      labs(
        x = "Total Sample Size", 
        y = paste("Minimum Detectable Effect Size (", outcome_var, ")"),
        title = paste0("Minimum Effect Size Detectable with ", 
                       power_target*100, "% Power"),
        subtitle = paste0("Current sample size = ", current_sample_size,
                          ", Min. effect = ", round(min_raw, 3))
      ) +
      theme_minimal() +
      geom_text(aes(x = current_sample_size + 10, 
                    y = max(raw_effect) * 0.9,
                    label = paste("Current n =", current_sample_size)),
                hjust = 0, vjust = 1)
    
    # Add observed effect as a horizontal line if available
    if (!is.na(observed_effect)) {
      p1 <- p1 + geom_hline(yintercept = cohens_d, 
                             linetype = "dotted", color = "darkred") +
        geom_text(aes(x = min(sample_sizes) + 10, 
                      y = cohens_d + 0.02,
                      label = paste("Observed effect d =", round(cohens_d, 3))),
                  hjust = 0, color = "darkred")
      
      p2 <- p2 + geom_hline(yintercept = observed_effect, 
                             linetype = "dotted", color = "darkred") +
        geom_text(aes(x = min(sample_sizes) + 10, 
                      y = observed_effect + 0.05 * max(raw_effect),
                      label = paste("Observed effect =", round(observed_effect, 3))),
                  hjust = 0, color = "darkred")
    }
    
    # Save plots if filename provided
    if (!is.null(plot_filename)) {
      dir.create(dirname(plot_filename), recursive = TRUE, showWarnings = FALSE)
      
      # Save Cohen's d plot
      d_filename <- gsub("\\.png$", "_cohens_d.png", plot_filename)
      ggsave(d_filename, p1, width = 8, height = 6, dpi = 300)
      
      # Save raw effect plot
      raw_filename <- gsub("\\.png$", "_raw.png", plot_filename)
      ggsave(raw_filename, p2, width = 8, height = 6, dpi = 300)
      
      cat("Plots saved to", d_filename, "and", raw_filename, "\n")
    }
    
    return(list(
      current_sample_size = current_sample_size,
      treated_size = treated_size,
      control_size = control_size,
      observed_effect = observed_effect,
      observed_cohens_d = cohens_d,
      min_detectable_d = min_d,
      min_detectable_raw = min_raw,
      plot_cohens_d = p1,
      plot_raw = p2,
      data = plot_data
    ))
  } else {
    return(list(
      current_sample_size = current_sample_size,
      treated_size = treated_size,
      control_size = control_size,
      observed_effect = observed_effect,
      observed_cohens_d = cohens_d,
      min_detectable_d = min_d,
      min_detectable_raw = min_raw
    ))
  }
}

