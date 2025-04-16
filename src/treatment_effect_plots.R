#' Generate Treatment Effect Plots for Subgroup Models
#'
#' This function creates standardized treatment effect plots showing
#' the difference between defoliated and non-defoliated areas for
#' all your subgroup models.
#'
#' @param model_list List of fitted model objects
#' @param model_names Character vector of model names (must match model_list)
#' @param response_type Character either "Severity" or "Recovery"
#' @param output_dir Directory to save plots
#' @param combined_plot Logical; whether to also create a combined plot of all models
#' @return List of ggplot objects
#' @examples
#' # For severity models
#' sev_models <- list(fit.sev_1, fit.sev_2, fit.sev_3, fit.sev_4)
#' names <- c("0-2 years", "3-5 years", "6-9 years", "10+ years")
#' plot_treatment_effects(sev_models, names, "Severity", "plots/treat_effects/")
generate_treatment_effect_plots <- function(model_list, model_names, response_type,
                                           output_dir = "/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/",
                                           combined_plot = TRUE,
                                           all_categories_plot = TRUE) {
  
  # Load required libraries
  require(ggplot2)
  require(dplyr)
  require(marginaleffects)
  require(patchwork)
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Validate inputs
  if (length(model_list) != length(model_names)) {
    stop("Length of model_list and model_names must be the same")
  }
  
  if (!response_type %in% c("Severity", "Recovery")) {
    stop("response_type must be either 'Severity' or 'Recovery'")
  }
  
  # Set y-axis label based on response type
  y_label <- if(response_type == "Severity") {
    "Burn Severity (RBR)"
  } else {
    "Recovery Magnitude (%)"
  }
  
  # Set plot colors 
  colors <- c("Defoliated" = "#8B0000A0", "Non-Defoliated" = "#FF8C00A0")  
  
  
  # List to store plot objects
  plot_list <- list()
  all_predictions <- data.frame()
  
  # Generate plots for each model
  for (i in 1:length(model_list)) {
    # Skip if model is NULL
    if (is.null(model_list[[i]])) {
      cat("Skipping", model_names[i], "- model is NULL\n")
      next
    }
    
    tryCatch({
      # Create predictions
      preds <- plot_predictions(model_list[[i]], condition = c("history"), draw = FALSE)
      
      # Add model info
      preds <- preds %>% 
        mutate(Model = model_names[i]) %>%
        mutate(history = case_when(history == 1 ~ "Defoliated",
                                   history == 0 ~ "Non-Defoliated"))
      
      # Store for combined plot
      all_predictions <- bind_rows(all_predictions, preds)
      
      # Create individual plot
      p <- ggplot(preds, aes(x = history, y = estimate, color = history)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
        labs(
          title = paste(response_type, "-", model_names[i]),
          y = y_label, 
          x = "Defoliation History", 
          color = "History"
        ) +
        scale_color_manual(values = colors) +
        theme_bw() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 10)
        )
      
      # Store plot object
      plot_list[[model_names[i]]] <- p
      
      # Save individual plot
      plot_filename <- paste0("fig_", 
                             tolower(response_type), 
                             "_treat_effect_", 
                             gsub(" ", "_", tolower(model_names[i])), 
                             ".png")
      
      ggsave(
        plot = p,
        filename = file.path(output_dir, plot_filename),
        width = 6, 
        height = 4, 
        dpi = 300
      )
      
      cat("Saved plot for", model_names[i], "to", file.path(output_dir, plot_filename), "\n")
      
    }, error = function(e) {
      cat("Error creating plot for", model_names[i], ":", e$message, "\n")
    })
  }
  
  # Create combined plot if requested and if we have data
  if (combined_plot && nrow(all_predictions) > 0) {
    # Create faceted plot
    p_combined <- ggplot(all_predictions, aes(x = history, y = estimate, color = history)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
      facet_wrap(~ Model, scales = "free_y", ncol = 2) +
      labs(
        title = paste(response_type, "by Time Since Defoliation"),
        y = y_label, 
        x = "Defoliation History", 
        color = "History"
      ) +
      scale_color_manual(values = colors) +
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 11),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    # Save combined plot
    combined_filename <- paste0("fig_", tolower(response_type), "_treat_effects_combined.png")
    ggsave(
      plot = p_combined,
      filename = file.path(output_dir, combined_filename),
      width = 10, 
      height = 8, 
      dpi = 300
    )
    
    cat("Saved combined plot to", file.path(output_dir, combined_filename), "\n")
    
    # Add to return list
    plot_list[["combined"]] <- p_combined
  }
  if (all_categories_plot && nrow(all_predictions) > 0) {
    # Create plot with all categories on x axis
    p_all_categories <- ggplot(all_predictions, aes(x = Model, 
                        y = estimate, 
                        color = history)) +
      geom_point(size = 3,
                position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = conf.low, 
                        ymax = conf.high),
                        position = position_dodge(width = 0.5),
                        width = 0.2) +
      labs(
        y = y_label, 
        x = "Defoliation History", 
        color = "History"
      ) +
      scale_color_manual(values = colors) +
      scale_x_discrete(limits = c("0-2 years","3-5 years","6-9 years", "10+ years")) + 
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 11),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14)
      )

      # save all categories plot
      all_categories_filename <- paste0("fig_", tolower(response_type), "_treat_effects_all_categories.png")
      ggsave(
        plot = p_all_categories,
        filename = file.path(output_dir, all_categories_filename),
        width = 10, 
        height = 8, 
        dpi = 300
      )
      cat("Saved all categories plot to", file.path(output_dir, all_categories_filename), "\n")
      # Add to return list    
      plot_list[["all_categories"]] <- p_all_categories
  }
  
  return(plot_list)
}

#' Run Treatment Effect Plots for All Subgroup Models
#'
#' This function loads all the subgroup models and generates treatment effect plots.
generate_all_treatment_effect_plots <- function() {
  # Result directory
  result_dir <- "/home/goldma34/sbw-wildfire-impact-recovery/results/subgroup/"
  output_dir <- "/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/"
  
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
  
  # Load all model fits
  cat("Loading models...\n")
  fit_sev_1 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup1_severity.RDS"))
  fit_sev_2 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup2_severity.RDS"))
  fit_sev_3 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup3_severity.RDS"))
  fit_sev_4 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup4_severity.RDS"))
  
  fit_rec_1 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup1_recovery.RDS"))
  fit_rec_2 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup2_recovery.RDS"))
  fit_rec_3 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup3_recovery.RDS"))
  fit_rec_4 <- safe_read_rds(paste0(result_dir, "fit_model_subgroup4_recovery.RDS"))
  
  # Define model groups
  sev_models <- list(fit_sev_1, fit_sev_2, fit_sev_3, fit_sev_4)
  rec_models <- list(fit_rec_1, fit_rec_2, fit_rec_3, fit_rec_4)
  
   model_names <- c("0-2 years", "3-5 years", "6-9 years", "10+ years")
  
  # Generate severity plots
  cat("Generating severity treatment effect plots...\n")
  sev_plots <- generate_treatment_effect_plots(
    sev_models, 
    model_names, 
    "Severity", 
    output_dir
  )
  
  # Generate recovery plots
  cat("Generating recovery treatment effect plots...\n")
  rec_plots <- generate_treatment_effect_plots(
    rec_models, 
    model_names, 
    "Recovery", 
    output_dir
  )
  
  # Generate comprehensive plot with both responses if all models loaded
  if (all(c(!is.null(fit_sev_1), !is.null(fit_rec_1)))) {
    # Compare severity and recovery for each time period
    for (i in 1:length(model_names)) {
      if (!is.null(sev_models[[i]]) && !is.null(rec_models[[i]])) {
        comparison_plot <- sev_plots[[model_names[i]]] + rec_plots[[model_names[i]]] +
          plot_layout(ncol = 2, guides = "collect") & 
          theme(legend.position = "bottom")
        
        comparison_filename <- paste0("fig_comparison_", 
                                    gsub(" ", "_", tolower(model_names[i])), 
                                    ".png")
        
        ggsave(
          plot = comparison_plot,
          filename = file.path(output_dir, comparison_filename),
          width = 12, 
          height = 5, 
          dpi = 300
        )
      }
    }
  }
  cat("All treatment effect plots generated!\n")
}

# Run the function to generate all plots when this script is sourced
if (interactive()) {
  cat("To generate all treatment effect plots, run:\n")
  cat("generate_all_treatment_effect_plots()\n")
} else {
  # Auto-run when sourced non-interactively
  generate_all_treatment_effect_plots()
}


#' Generate Treatment Effect Plots for Intermediate Lag Analysis
#'
#' This function generates plots showing the estimated treatment effects
#' from propensity score matched models of defoliation on fire outcomes.
#' 
#' @param response_type Character, either "Severity" or "Recovery"
#' @param output_dir Directory where plots should be saved
#' @param combined_plot Logical, whether to create a combined plot
#' @param colors Vector of colors for the plots
#'
#' @return List of generated plot objects
generate_intermediate_treatment_effect_plots <- function(
  response_type = c("Severity", "Recovery"),
  output_dir = "/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/intermediate/",
  combined_plot = TRUE,
  colors = c("Defoliated" = "#8B0000A0", "Non-Defoliated" = "#FF8C00A0")) {
  
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  library(marginaleffects)
  
  # Validate response type
  response_type <- match.arg(response_type)
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Determine y-axis label based on response type
  y_label <- if(response_type == "Severity") {
    "Burn Severity (RBR)"
  } else {
    "Recovery (%)"
  }
  
  # Set up paths to model files based on response type
  if (response_type == "Severity") {
    model_file <- "/home/goldma34/sbw-wildfire-impact-recovery/results/subgroup/fit_model_intermediate_severity.RDS"
    model_names <- "3-9 years"
  } else {
    model_file <- "/home/goldma34/sbw-wildfire-impact-recovery/results/subgroup/fit_model_intermediate_recovery.RDS"
    model_names <- "3-9 years"
  }
  
  # Load model
  model <- readRDS(model_file)
  
  # Initialize output list
  plot_list <- list()
  all_predictions <- data.frame()
  
  # Create plot for the model
  tryCatch({
    # Generate treatment effect estimates
    treatment_effects <- avg_comparisons(
      model, 
      variables = "history",
      newdata = subset(model$model, history == 1)
    )
    
    # Convert to data frame for plotting
    pred_df <- as.data.frame(treatment_effects) %>%
      mutate(
        history = ifelse(contrast == "1 - 0", "Defoliated", "Non-Defoliated"),
        Model = model_names
      )
    
    # Create plot
    p <- ggplot(pred_df, aes(x = history, y = estimate, color = history)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
      labs(
        title = paste("Treatment Effect on", response_type, "for Intermediate Lag (3-9 years)"),
        y = y_label, 
        x = "Defoliation History", 
        color = "History"
      ) +
      scale_color_manual(values = colors) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    # Add to return list
    plot_list[[model_names]] <- p
    
    # Combine with other predictions
    all_predictions <- bind_rows(all_predictions, pred_df)
    
    # Save individual plot
    plot_filename <- paste0("fig_intermediate_", tolower(response_type), "_treatment_effect.png")
    ggsave(
      plot = p,
      filename = file.path(output_dir, plot_filename),
      width = 6, 
      height = 4, 
      dpi = 300
    )
    
    cat("Saved plot to", file.path(output_dir, plot_filename), "\n")
    
  }, error = function(e) {
    cat("Error creating plot for intermediate model:", e$message, "\n")
  })
  
  # Return the plot list
  return(plot_list)
}

#' Generate plots for both severity and recovery
#'
#' This function calls generate_intermediate_treatment_effect_plots 
#' for both severity and recovery.
#'
#' @return NULL
generate_all_intermediate_treatment_effect_plots <- function() {
  sev_plots <- tryCatch({
    generate_intermediate_treatment_effect_plots("Severity")
  }, error = function(e) {
    cat("Error generating severity plots:", e$message, "\n")
    list()
  })
  
  rec_plots <- tryCatch({
    generate_intermediate_treatment_effect_plots("Recovery")
  }, error = function(e) {
    cat("Error generating recovery plots:", e$message, "\n")
    list()
  })
  
  cat("Treatment effect plot generation complete.\n")
  
  # Return both sets of plots
  return(list(severity = sev_plots, recovery = rec_plots))
}

# Run the function when this script is sourced (if not being sourced for just the functions)
if (!interactive()) {
  generate_all_intermediate_treatment_effect_plots()
}

#'Generate comparison intermediate treatment effects comparison plots
#' 
#' This function generates a comparison plot of treatment effects
#' for intermediate lag (3-9 years) and the 3-5 and 6-9 years subgroup models.
#' @param response_type Character, either "Severity" or "Recovery"
#' @param output_dir Directory where plots should be saved
#' @param combined_plot Logical, whether to create a combined plot
#' @param colors Vector of colors for the plots
#'
#' @return List of generated plot objects
generate_combined_intermediate_treatment_effect_plots <- function(
  response_type = c("Severity", "Recovery"),
  output_dir = "/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/intermediate/",
  combined_plot = TRUE,
  colors = c("Defoliated" = "#8B0000A0", "Non-Defoliated" = "#FF8C00A0")) {
  
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  library(marginaleffects)

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
  
  # Function to get significance stars
  get_significance_stars <- function(p_value) {
    if(is.na(p_value)) return("n.s.")
    if(p_value > 0.05) return("n.s.") 
    else if(p_value <= 0.05 & p_value > 0.01) return("*")
    else if(p_value <= 0.01 & p_value > 0.001) return("**")
    else return("***")
  }
  
  result_dir <- "/home/goldma34/sbw-wildfire-impact-recovery/results/subgroup/"

  # Validate response type
  response_type <- match.arg(response_type)
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Determine y-axis label based on response type
  y_label <- if(response_type == "Severity") {
    "Burn Severity (RBR)"
  } else {
    "Recovery (%)"
  }
  
  # Set up paths to model files based on response type
  if (response_type == "Severity") {
    int_model_file <- paste0(result_dir, "fit_model_intermediate_severity.RDS")
    early_model_file <- paste0(result_dir, "fit_model_subgroup2_severity.RDS")
    peak_model_file <- paste0(result_dir, "fit_model_subgroup3_severity.RDS")
  } else {
    int_model_file <- paste0(result_dir, "fit_model_intermediate_recovery.RDS")
    early_model_file <- paste0(result_dir, "fit_model_subgroup2_recovery.RDS")
    peak_model_file <- paste0(result_dir, "fit_model_subgroup3_recovery.RDS")
  }
  
  # Load models
  int_model <- safe_read_rds(int_model_file)
  early_model <- safe_read_rds(early_model_file)
  peak_model <- safe_read_rds(peak_model_file)
  
  # Initialize plot list
  plot_list <- list()
  all_predictions <- data.frame()
  significance_labels <- data.frame()
  
  # Define model groups and model names
  models <- list(int_model, early_model, peak_model)
  model_names <- c("Intermediate (3-9 years)", "Early (3-5 years)", "Peak (6-9 years)")
  
  # Process each model
  for (i in 1:length(models)) {
    # Skip if model is NULL
    if (is.null(models[[i]])) {
      cat("Skipping", model_names[i], "- model is NULL\n")
      next
    }
    
    tryCatch({
      # Get sample size
      n_size <- nrow(models[[i]]$model)
      
      # Generate predictions for this model
      preds <- plot_predictions(models[[i]], condition = c("history"), draw = FALSE)
      
      # Get treatment effect and p-value
      treatment_effects <- avg_comparisons(
        models[[i]], 
        variables = "history",
        newdata = subset(models[[i]]$model, history == 1)
      )
      
      # Extract p-value for significance stars
      p_value <- as.data.frame(treatment_effects)$p.value[1]
      sig_stars <- get_significance_stars(p_value)
      
      # Create significance label data 
      sig_label <- data.frame(
        Model = model_names[i],
        label = paste0(sig_stars, "\n", "n=", n_size),
        y_pos = max(preds$conf.high) + 0.05 * diff(range(preds$estimate))
      )
      
      significance_labels <- bind_rows(significance_labels, sig_label)
      
      # Add model info
      preds <- preds %>% 
        mutate(Model = model_names[i]) %>%
        mutate(history = case_when(history == 1 ~ "Defoliated",
                                  history == 0 ~ "Non-Defoliated"))
      
      # Store for combined plot
      all_predictions <- bind_rows(all_predictions, preds)
      
      # Create individual plot
      p <- ggplot(preds, aes(x = history, y = estimate, color = history)) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
        labs(
          title = paste0(response_type, " - ", model_names[i], " (", sig_stars, ", n=", n_size, ")"),
          y = y_label, 
          x = "Defoliation History", 
          color = "History"
        ) +
        scale_color_manual(values = colors) +
        theme_bw() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12)
        )
      
      # Store in plot list
      plot_list[[model_names[i]]] <- p
      
      # Save individual plot
      plot_filename <- paste0("fig_", 
                             tolower(response_type), 
                             "_treat_effect_", 
                             gsub(" ", "_", tolower(gsub("[()]", "", model_names[i]))), 
                             ".png")
      
      ggsave(
        plot = p,
        filename = file.path(output_dir, plot_filename),
        width = 6, 
        height = 4, 
        dpi = 300
      )
      
      cat("Saved plot for", model_names[i], "to", file.path(output_dir, plot_filename), "\n")
      
    }, error = function(e) {
      cat("Error creating plot for", model_names[i], ":", e$message, "\n")
    })
  }
  
  # Create combined plot if we have data
  if (combined_plot && nrow(all_predictions) > 0) {
    # Calculate overall y-range for proper placement of significance labels
    y_range <- range(c(all_predictions$conf.low, all_predictions$conf.high), na.rm = TRUE)
    y_buffer <- 0.08 * diff(y_range)
    
    # Update y positions for significance labels based on actual data
    significance_labels <- significance_labels %>%
      group_by(Model) %>%
      mutate(
        y_pos = max(all_predictions$conf.high[all_predictions$Model == Model]) + y_buffer
      )
    
    p_combined <- ggplot(all_predictions, 
      aes(x = Model, y = estimate, color = history)) +
      geom_point(size = 3, position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                   width = 0.2,
                   position = position_dodge(width = 0.5)) +
      # Add significance and sample size labels
      geom_text(
        data = significance_labels,
        aes(x = Model, y = y_pos, label = label),
        color = "black",
        inherit.aes = FALSE,
        size = 3.5
      ) +
      labs(
           y = y_label, 
           x = "Time Window", 
           color = "History") +
      scale_color_manual(values = colors) +
      theme_bw() +
      # Expand y-axis to make room for significance labels
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      theme(
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 11),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14),
        plot.caption = element_text(hjust = 0, size = 9, face = "italic")
      )
    
    # Save combined plot
    combined_filename <- paste0("fig_", tolower(response_type), "_treat_effects_intermediate_comparison.png")
    ggsave(
      plot = p_combined,
      filename = file.path(output_dir, combined_filename),
      width = 8, 
      height = 6, 
      dpi = 300
    )
    
    cat("Saved combined plot to", file.path(output_dir, combined_filename), "\n")
    
    # Add to return list
    plot_list[["combined"]] <- p_combined
  }
  
  # Return the plot list
  return(plot_list)
}

#' Generate plots for comparison of intermediate treatment effects
#'
#' This function calls generate_combined_intermediate_treatment_effect_plots
#'
#' @return NULL
generate_all_combined_intermediate_treatment_effect_plots <- function() {
  sev_plots <- tryCatch({
    generate_combined_intermediate_treatment_effect_plots("Severity")
  }, error = function(e) {
    cat("Error generating severity plots:", e$message, "\n")
    list()
  })
  
  rec_plots <- tryCatch({
    generate_combined_intermediate_treatment_effect_plots("Recovery")
  }, error = function(e) {
    cat("Error generating recovery plots:", e$message, "\n")
    list()
  })
  
  cat("Treatment effect plot generation complete.\n")
  
  # Return both sets of plots
  return(list(severity = sev_plots, recovery = rec_plots))
}

# Run the function when this script is sourced (if not being sourced for just the functions)
if (!interactive()) { 
  cat("To generate all combined intermediate treatment effect plots, run:\n")
  cat("generate_all_combined_intermediate_treatment_effect_plots()\n") 
  } else {
  generate_all_combined_intermediate_treatment_effect_plots()
}