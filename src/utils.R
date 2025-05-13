#' Remove title from a ggplot object
#' 
#' This function takes a ggplot object and removes its title while preserving
#' all other plot elements.
#'
#' @param plot A ggplot object
#' @return A ggplot object with the title removed
#' @examples
#' # p <- ggplot(mtcars, aes(x = mpg, y = disp)) + 
#' #      geom_point() +
#' #      labs(title = "My Title")
#' # p_without_title <- remove_title(p)
remove_title <- function(plot) {
  # Create a new ggplot object with the title removed
  plot + 
    labs(title = NULL) +
    theme(plot.title = element_blank())
}




#' Calculate R-squared for GLS Model
#'
#' This function calculates the R-squared value for a Generalized Least Squares (GLS) model.
#'
#' @param model A GLS model object of class `gls`.
#' @param formula A formula specifying the model.
#' @return The R-squared value for the GLS model.
#' @examples
#' library(nlme)
#' gls_model <- gls(y ~ x, data = your_data, correlation = corExp(form = ~ x))
#' rsquared.gls(gls_model, y ~ x)
#' @export
rsquared.gls <- function(model, formula) {
  # Extract the data used to fit the model
  data <- getData(model)
  
  # Extract the response variable
  response <- model$terms[[2]]
  y <- data[[as.character(response)]]
  
  # Create the model matrix
  X <- model.matrix(formula, data = data)
  
  # Calculate the variance of the fitted values
  sigmaF <- var(as.vector(model$coefficients %*% t(X)))
  
  # Calculate the variance of the residuals
  sigmaE <- var(resid(model))
  
  # Calculate R-squared
  R_squared <- sigmaF / (sigmaF + sigmaE)
  
  print(R_squared)
}

#' Create coefficient summary tables for GLS models
#'
#' This function extracts coefficient estimates, standard errors, t-values, 
#' and p-values from GLS models and exports them to a CSV file.
#'
#' @param severity_model The GLS model for severity analysis
#' @param recovery_model The GLS model for recovery analysis
#' @param output_path Path where the CSV file should be saved
#' @param rename_vars Logical, whether to rename variables for better readability
#'
#' @return A data frame containing the coefficient tables (invisibly)
export_gls_coefficient_tables <- function(
  severity_model,
  recovery_model,
  output_path = "/home/goldma34/sbw-wildfire-impact-recovery/results/spatial_analysis/gls_model_coefficients.csv",
  rename_vars = TRUE
) {
  library(dplyr)
  
  # Extract coefficients from severity model
  sev_coef <- as.data.frame(summary(severity_model)$tTable)
  sev_coef$Variable <- rownames(sev_coef)
  sev_coef$Model <- "Severity"
  # Keep original column names from tTable
  
  # Extract coefficients from recovery model
  rec_coef <- as.data.frame(summary(recovery_model)$tTable)
  rec_coef$Variable <- rownames(rec_coef)
  rec_coef$Model <- "Recovery"
  # Keep original column names from tTable
  
  # Combine both models
  all_coef <- bind_rows(sev_coef, rec_coef)
  
  # Optionally rename variables for better readability
  if (rename_vars) {
    all_coef <- all_coef %>%
      mutate(Variable = case_when(
        Variable == "(Intercept)" ~ "Intercept",
        Variable == "host_pct" ~ "Host Percentage",
        Variable == "Cumulative_Years_Defol:window_opp1" ~ "Years Defoliated x Immediate Post-Outbreak",
        Variable == "Cumulative_Years_Defol:window_opp2" ~ "Years Defoliated x Early Window",
        Variable == "Cumulative_Years_Defol:window_opp3" ~ "Years Defoliated x Peak Window",
        Variable == "Cumulative_Years_Defol:window_opp4" ~ "Years Defoliated x Late Window",
        Variable == "isi_90" ~ "Initial Spread Index (ISI)",
        Variable == "dc_90" ~ "Drought Code (DC)",
        Variable == "dmc_90" ~ "Duff Moisture Code (DMC)",
        Variable == "ffmc_90" ~ "Fine Fuel Moisture Code (FFMC)",
        Variable == "bui_90" ~ "Buildup Index (BUI)",
        Variable == "fwi_90" ~ "Fire Weather Index (FWI)",
        Variable == "mean_tri" ~ "Mean Terrain Ruggedness Index",
        Variable == "rbr_w_offset" ~ "Burn Severity (RBR)",
        Variable == "mean_temperature" ~ "Mean Temperature",
        Variable == "sum_precipitation_mm" ~ "Total Precipitation (mm)",
        Variable == "x" ~ "Longitude",
        Variable == "y" ~ "Latitude",
        TRUE ~ Variable
      ))
  }
  # rename pvalues and t value columns
    colnames(all_coef)[colnames(all_coef) == "p-value"] <- "p.value"
    colnames(all_coef)[colnames(all_coef) == "t-value"] <- "t.value"
  
  # Format p-values for readability and add significance
  all_coef <- all_coef %>%
    mutate(
      # Create a formatted p-value column
      p_formatted = format.pval(p.value, digits = 3, eps = 0.001),
      # Add significance stars
      Significance = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      ),
      # Round numeric columns to 3 decimal places
      Value = round(Value, 3),
      Std.Error = round(Std.Error, 3),
      t.value = round(t.value, 3)
    )
  
  # First, reorganize the columns
  all_coef <- all_coef %>%
    select(Model, Variable, Value, Std.Error, t.value, p.value, p_formatted, Significance)
  
  # Rename columns for better readability
  colnames(all_coef) <- c("Model", "Variable", "Estimate", "Std.Error", "t_value", "p_value", "p_formatted", "Significance")
  
  # Write to CSV
  write.csv(all_coef, file = output_path, row.names = FALSE)
  
  cat("Coefficient tables exported to:", output_path, "\n")
  
  # Return the table invisibly
  invisible(all_coef)
}


#' Add significance annotation to treatment effect plots
#'
#' This function adds a significance annotation (stars or "n.s.") between two points
#' on a treatment effect plot based on the p-value from model results.
#'
#' @param plot The ggplot object to modify
#' @param model_results The results data frame containing p.value
#' @param y_position Vertical position for the annotation
#' @param x_positions Vector of x positions that the line should connect
#' @param bracket_height Height of the bracket lines
#' @param sig_cex Size of the text
#'
#' @return A modified ggplot with significance annotation added
add_significance_annotation <- function(plot, 
                                       model_results, 
                                       y_position = NULL,
                                       x_positions = c(1, 2),
                                       bracket_height = 0.1,
                                       sig_cex = 4) {
  
  # Extract p-value
  p_value <- model_results$p.value[1]  # Assuming the first row has the p-value
  
  # Determine significance symbol
  sig_symbol <- if (is.na(p_value)) {
    "n.s."
  } else if (p_value > 0.05) {
    "n.s."
  } else if (p_value <= 0.001) {
    "***"
  } else if (p_value <= 0.01) {
    "**"
  } else if (p_value <= 0.05) {
    "*"
  }
  
  # If y_position not provided, calculate based on plot data
  if (is.null(y_position)) {
    # Get data from the plot
    plot_data <- ggplot2::ggplot_build(plot)$data[[1]]
    # Calculate max y value plus 10% margin
    y_max <- max(plot_data$ymax, na.rm = TRUE)
    y_position <- y_max * 1.1
  }
  
  # Add significance annotation
  plot + 
    # Add line connecting points
    annotate("segment", 
             x = x_positions[1], 
             xend = x_positions[2], 
             y = y_position, 
             yend = y_position,
             color = "black", 
             size = 0.5) +
    # Add bracket end lines
    annotate("segment", 
             x = x_positions[1], 
             xend = x_positions[1], 
             y = y_position - bracket_height, 
             yend = y_position,
             color = "black", 
             size = 0.5) +
    annotate("segment", 
             x = x_positions[2], 
             xend = x_positions[2], 
             y = y_position - bracket_height, 
             yend = y_position,
             color = "black", 
             size = 0.5) +
    # Add significance text
    annotate("text", 
             x = mean(x_positions), 
             y = y_position + bracket_height, 
             label = sig_symbol,
             size = sig_cex) +
    # Expand plot limits to make room for annotation
    scale_y_continuous(limits = c(NA, y_position * 1.2), expand = expansion(mult = c(0.05, 0.2)))
}

#' Adjust Plot Size
#'
#' This function adjusts the first plot to match the size of the second plot.
#' It extracts the dimensions of the second plot and applies them to the first plot.
#'
#' @param plot_a The ggplot object to be resized
#' @param plot_b The reference ggplot object whose size will be used
#' @param preserve_aspect_ratio Logical, whether to preserve the aspect ratio of plot_a (default: TRUE)
#'
#' @return A ggplot object (plot_a) with adjusted size
#'
#' @examples
#' # p1 <- ggplot(mtcars, aes(mpg, disp)) + geom_point()
#' # p2 <- ggplot(mtcars, aes(wt, hp)) + geom_point() + theme(aspect.ratio = 0.8)
#' # p1_resized <- match_plot_size(p1, p2)
match_plot_size <- function(plot_a, plot_b, preserve_aspect_ratio = TRUE) {
  # Extract theme elements from plot_b
  b_theme <- ggplot_build(plot_b)$layout$panel_params[[1]]
  
  # Extract aspect ratio from plot_b if it exists
  if ("aspect.ratio" %in% names(plot_b$theme)) {
    target_aspect_ratio <- plot_b$theme$aspect.ratio
  } else {
    # If no aspect ratio is explicitly set, calculate it from the plot dimensions
    b_width <- diff(b_theme$x.range)
    b_height <- diff(b_theme$y.range)
    target_aspect_ratio <- b_height / b_width
  }
  
  # Apply the extracted dimensions to plot_a
  if (preserve_aspect_ratio) {
    # Only set the aspect ratio
    plot_a + theme(aspect.ratio = target_aspect_ratio)
  } else {
    # Get current dimensions of plot_a
    a_theme <- ggplot_build(plot_a)$layout$panel_params[[1]]
    a_width <- diff(a_theme$x.range)
    a_height <- diff(a_theme$y.range)
    
    # Calculate required scaling factors
    x_scale <- a_width / b_width
    y_scale <- a_height / b_height
    
    # Apply scaling to plot_a (this is more complex and may not work for all plots)
    plot_a + 
      coord_cartesian(
        xlim = c(a_theme$x.range[1], a_theme$x.range[1] + a_width/x_scale),
        ylim = c(a_theme$y.range[1], a_theme$y.range[1] + a_height/y_scale)
      ) +
      theme(aspect.ratio = target_aspect_ratio)
  }
}

# Alternative simpler function that just matches the aspect ratio
match_aspect_ratio <- function(plot_a, plot_b) {
  # Get the aspect ratio from plot_b
  if ("aspect.ratio" %in% names(plot_b$theme)) {
    # If plot_b has explicit aspect ratio, use it
    aspect_ratio <- plot_b$theme$aspect.ratio
  } else {
    # Try to compute from plot dimensions
    b_built <- ggplot_build(plot_b)
    if (length(b_built$layout$coord$aspect()) > 0) {
      aspect_ratio <- b_built$layout$coord$aspect()
    } else {
      # Default if we can't determine
      aspect_ratio <- 1
    }
  }
  
  # Apply the aspect ratio to plot_a
  plot_a + theme(aspect.ratio = aspect_ratio)
}

#' Trim Y-Axis Limits
#'
#' This function reduces the y-axis limits of a plot to the maximum height of 
#' confidence intervals or points plus a specified buffer.
#'
#' @param plot A ggplot object, typically with error bars
#' @param buffer Amount to add above the maximum value (absolute or percentage)
#' @param buffer_is_percent Logical, whether buffer is a percentage (TRUE) or absolute value (FALSE)
#' @param extract_from_data Logical, whether to extract limit from plot data (recommended)
#'
#' @return A ggplot object with adjusted y-axis limits
#'
#' @examples
#' # p <- ggplot(df, aes(x=group, y=mean)) + 
#' #      geom_point() + 
#' #      geom_errorbar(aes(ymin=lower, ymax=upper))
#' # p_trimmed <- trim_y_axis(p, buffer=50)
trim_y_axis <- function(plot, buffer=50, buffer_is_percent=FALSE, extract_from_data=TRUE) {
  
  if (extract_from_data) {
    # Extract the data from the plot
    plot_data <- ggplot_build(plot)$data
    
    # Look for error bars or confidence intervals in various formats
    # Different geoms store data differently, so we need to check multiple places
    y_values <- numeric()
    
    # Loop through all layers in the plot data
    for (i in seq_along(plot_data)) {
      layer_data <- plot_data[[i]]
      
      # Check for different column names that might contain max y values
      possible_columns <- c("ymax", "conf.high", "upper", "y")
      
      for (col in possible_columns) {
        if (col %in% colnames(layer_data)) {
          y_values <- c(y_values, layer_data[[col]])
        }
      }
    }
    
    # Find the maximum value
    if (length(y_values) > 0) {
      y_max <- max(y_values, na.rm = TRUE)
    } else {
      warning("Could not extract y values from plot. Using default limits.")
      return(plot)
    }
  } else {
    # Extract y range from plot limits
    y_range <- ggplot_build(plot)$layout$panel_params[[1]]$y.range
    y_max <- y_range[2]
  }
  
  # Calculate new upper limit
  if (buffer_is_percent) {
    y_upper <- y_max * (1 + buffer/100)
  } else {
    y_upper <- y_max + buffer
  }
  
  # Apply new limits to the plot
  plot + 
    coord_cartesian(ylim = c(NA, y_upper)) +
    # Ensure the scale doesn't add more expansion
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
}



# Function to safely write CSV files
safe_write_csv <- function(data, filepath) {
  if (!file.exists(filepath)) {
    # Create directory if it doesn't exist
    dir_path <- dirname(filepath)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Write the file
    write.csv(data, filepath)
    cat(sprintf("File created: %s\n", filepath))
  } else {
    cat(sprintf("File already exists (skipping): %s\n", filepath))
  }
}

# Function to safely write RDS files
safe_write_rds <- function(data, filepath) {
  if (!file.exists(filepath)) {
    # Create directory if it doesn't exist
    dir_path <- dirname(filepath)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    
    # Write the file
    saveRDS(data, filepath)
    cat(sprintf("File created: %s\n", filepath))
  } else {
    cat(sprintf("File already exists (skipping): %s\n", filepath))
  }
}



# publication plot main psm results
publication_panel_all_fires <- function(study_area, treatment_effects, residuals_map, outcome_variable) {
  output_dir <- "./plots/maps/panels/" 
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  if (outcome_variable == "severity") {
    output_file <- paste0(output_dir, "fig_", outcome_variable, "_panel_all_fires.png")          
  } else if (outcome_variable == "recovery") {
    output_file <- paste0(output_dir, "fig_", outcome_variable, "_panel_all_fires.png")          
  } 

  study_area_mod <- study_area + 
  theme(
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 16),
    legend.margin = margin(t = 10, r = 0, b = 0, l = 0),
    plot.margin = margin(0, 0, 0, 0),  # Zero margin on all sides
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.08, 0.98)  # Position tag further to the right
  ) + guides(colour = guide_legend(override.aes = list(size=10)))

  treatment_effects_mod <- treatment_effects +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      plot.margin = margin(0, 0, 0, 0),  # Zero margin on all sides
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0.08, 0.98)  # Position tag further to the right
    )
  if (outcome_variable == "severity") {
    treatment_effects_mod <- treatment_effects_mod + labs(y = "Median Burn Severity")
  } else if (outcome_variable == "recovery") {
    treatment_effects_mod <- treatment_effects_mod + labs(y = "Recovery Magnitude %")
  }

  if (outcome_variable == "severity"){
  residuals_mod <- residuals_map + 
      theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "right",
      legend.key.size = unit(1.5, "cm"),  # Increase overall key size
      legend.key.height = unit(1.2, "cm"), # Set specific height
      legend.key.width = unit(1.5, "cm"),  # Set specific width
      legend.margin = margin(t = 0, r = 0, b = 0, l = 5),
      legend.box.margin = margin(0, 0, 0, 10),  # Add space around legend box
      legend.text = element_text(size = 16),  # Increase text size
      plot.margin = margin(0, 0, 0, 0),  # Zero margin on all sides
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0.08, 0.98)  # Position tag further to the right
    ) + 
    guides(
      # For color legend (adjust based on your aesthetic mappings)
      color = guide_colourbar(
        barwidth = unit(1, "cm"),
        barheight = unit(5, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        label.hjust = 0.5,
      ),
      # For fill legend if present
      fill = guide_legend(
        override.aes = list(size = 8),
        keywidth = unit(2, "cm"),
        keyheight = unit(1, "cm")
      ),
      # For other aesthetics if present
      size = guide_legend(keywidth = unit(2, "cm"), keyheight = unit(1, "cm"))
    ) +
    # Apply theme to all guides
    theme(
      legend.title = element_text(size = 15, face = "italic"),
      legend.spacing.y = unit(0.5, "cm")  # Add vertical spacing between legend items
    ) } else if (outcome_variable == "recovery") {
      residuals_mod <- residuals_map + 
      scale_color_gradient2(low = "#88CCEEA0", 
        mid = "white", high = "#8B0000A0", 
        midpoint = 50) +
        theme(
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.position = "right",
          legend.key.size = unit(1.5, "cm"),  # Increase overall key size
          legend.key.height = unit(1.2, "cm"), # Set specific height
          legend.key.width = unit(1.5, "cm"),  # Set specific width
          legend.margin = margin(t = 0, r = 0, b = 0, l = 5),
          legend.box.margin = margin(0, 0, 0, 10),  # Add space around legend box
          legend.text = element_text(size = 16),  # Increase text size
          plot.margin = margin(0, 0, 0, 0),  # Zero margin on all sides
          plot.tag = element_text(size = 14, face = "bold"),
          plot.tag.position = c(0.08, 0.98)  # Position tag further to the right
        ) + 
        guides(
          # For color legend (adjust based on your aesthetic mappings)
          color = guide_colourbar(
            barwidth = unit(1, "cm"),
            barheight = unit(5, "cm"),
            title.position = "top",
            title.hjust = 0.5,
            label.position = "right",
            label.hjust = 0.5,
          ),
          # For fill legend if present
          fill = guide_legend(
            override.aes = list(size = 8),
            keywidth = unit(2, "cm"),
            keyheight = unit(1, "cm")
          ),
          # For other aesthetics if present
          size = guide_legend(keywidth = unit(2, "cm"), keyheight = unit(1, "cm"))
        ) +
        # Apply theme to all guides
        theme(
          legend.title = element_text(size = 15, face = "italic"),
          legend.spacing.y = unit(0.5, "cm")  # Add vertical spacing between legend items
        ) 
    }


    



  # Build layout with tags already positioned
  bottom_row <- treatment_effects_mod + residuals_mod + 
    plot_layout(
      guides = "keep",
      widths = c(1, 1.3)
    ) &
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.spacing = unit(0, "pt")
    )

  tight_layout <- study_area_mod / bottom_row +
    plot_layout(
      heights = c(2, 1),
      guides = "keep"
    ) &
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.spacing = unit(0, "cm"),
      panel.spacing.y = unit(-0.5, "lines")  # Negative value to remove white space
    ) &
    plot_annotation(tag_levels = 'A', tag_suffix = ")")
          
                              
    print(tight_layout)
    ggsave(output_file, plot = tight_layout, width = 12, height = 10, dpi = 300)
    cat(sprintf("Figure saved to: %s\n", output_file))
    return(tight_layout)                          
                              
}
