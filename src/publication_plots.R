library(ggpubr)
library(ggplot2)
library(dplyr)
library(patchwork)

#utils
source("/home/goldma34/sbw-wildfire-impact-recovery/src/utils.R")



# area maps
matched_map_severity<-readRDS("/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/study_area_severity_map.rds")
matched_map_recovery<-readRDS("/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/study_area_recovery_map.rds")

# spatial residual plots
spat.resid.sev_plot<-readRDS("/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/residuals_map_severity.rds")
spat.resid.rec_plot<-readRDS("/home/goldma34/sbw-wildfire-impact-recovery/plots/maps/residuals_map_recovery.rds")

# treatment effects plots
treat.plots_sev <-readRDS("/home/goldma34/sbw-wildfire-impact-recovery/plots/treat_effects/fig_sev_treat_effect_all.rds")
treat.plots_rec <-readRDS("/home/goldma34/sbw-wildfire-impact-recovery/plots/fig_rec_treat_effect_all.rds")


#check plot margin 
matched_map_severity + 
  theme(
    plot.margin = margin(0, 0, 0, 0),  # No margins
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98)  # Position in top-left corner
  )


# Use the function to create your severity layout
severity_layout <- create_balanced_layout(
  left_plot = matched_map_severity,
  top_right_plot = treat.plots_sev,
  bottom_right_plot = spat.resid.sev_plot,
  width_ratio = c(1.5, 1),
  height_ratio = c(2, 1),
  tags = c( "a", "b", "c")
)

#
severity_layout <- create_balanced_layout_simple(
  left_plot = matched_map_severity,
  top_right_plot = treat.plots_sev,
  bottom_right_plot = spat.resid.sev_plot,
  width_ratio = c(1.5, 1),
  height_ratio = c(1.5, 1),
  tags = c("a", "b", "c")
)

severity_layout

# Save the layout with cropping
save_layout_cropped(
  layout = severity_layout,
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/publication/severity_panel_2.png",
  width = 10,
  height = 7,
  dpi = 300,
  crop = TRUE,
  border = 10  # Small border in pixels
)

# Create the recovery layout with the same settings
recovery_layout <- create_balanced_layout(
  left_plot = matched_map_recovery,
  top_right_plot = spat.resid.rec_plot,
  bottom_right_plot = treat.plots_rec,
  width_ratio = c(1.5, 1),
  height_ratio = c(2, 1),
  tags = c("d", "e", "f")
)

# Save the recovery layout with cropping
save_layout_cropped(
  layout = recovery_layout,
  filename = "/home/goldma34/sbw-wildfire-impact-recovery/plots/publication/recovery_panel.png",
  width = 10,
  height = 7,
  dpi = 300,
  crop = TRUE,
  border = 10
)





### testing

layout <- c(
  area(2, 1, 3, 2),
  area(1, 3, 2, 4),
  area(3, 3, 4, 4)
)

# Show the layout to make sure it looks as it should
plot(layout)

matched_map_severity + treat.plots_sev + spat.resid.sev_plot + plot_layout(design = layout)


ggplot_build(treat.plots_sev) 


ggplot_build(spat.resid.sev_plot)





# Function to make a plot window smaller
shrink_plot_window <- function(plot, 
                              margin_reduction = margin(2, 2, 2, 2, "pt"),
                              aspect_ratio = 0.7) {
  # Reduce the plot margins
  shrunk_plot <- plot +
    theme(
      # Reduce margins around the plot
      plot.margin = margin_reduction,
      # Adjust aspect ratio (width/height) - lower values make plot wider and less tall
      aspect.ratio = aspect_ratio,
      # Reduce space between plot elements
      panel.spacing = unit(0.5, "lines"),
      # Make text smaller to fit in smaller space
      text = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      # Reduce legend size
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      # Reduce title size if it exists
      plot.title = element_text(size = 9),
      # Make sure plot takes up all available space
      plot.background = element_rect(color = NA)
    )
  
  return(shrunk_plot)
}

# Apply the function to treat.plots_sev
treat.plots_sev_small <- shrink_plot_window(
  treat.plots_sev,
  margin_reduction = margin(1, 2, 1, 1, "pt"),  # Very small margins
  aspect_ratio = 0.6  # Make it wider and less tall
)

# Check the result
treat.plots_sev_small


grouped <-matched_map_severity + treat.plots_sev_small + spat.resid.sev_plot + plot_layout(design = layout)


treat.plots_sev_compact <- treat.plots_sev_small +
  theme(
    # Reduce top margin to almost zero to bring plot down
    plot.margin = margin(0, 2, 1, 1, "pt"),
    # Further reduce aspect ratio if needed
    aspect.ratio = 0.55
  )

# Create a custom layout with specific positioning
# This uses the area() function to define precise grid positions
compact_layout <- c(
  # Map takes the left side (columns 1-2) and full height (rows 1-4)
  area(1, 1, 4, 2),
  
  # Treatment plot is positioned in top right, but moved down by reducing its row span
  # Start at row 1.5 instead of 1 to move it down
  area(1.5, 3, 3, 4),
  
  # Residual plot at bottom right
  area(3, 3, 4, 4)
)

# Visualize the layout to check positioning
# plot(compact_layout)

# Create the grouped layout with the compact treatment plot and custom layout
grouped_compact <- matched_map_severity + 
                  treat.plots_sev_compact + 
                  spat.resid.sev_plot + 
                  plot_layout(design = compact_layout)

# Display the result
grouped_compact

# Alternative approach using relative heights and grid nesting
# This might give more precise control over positioning
right_column <- (treat.plots_sev_compact / spat.resid.sev_plot) + 
  plot_layout(
    heights = c(1.5, 1),  # Adjust this ratio to control vertical positioning
    guides = "collect"
  ) &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = unit(0, "pt")  # Remove space between plots
  )

# Combine with the left plot
grouped_nested <- (matched_map_severity | right_column) +
  plot_layout(
    widths = c(1.5, 1),
    guides = "collect"
  ) &
  theme(plot.margin = margin(0, 0, 0, 0))

# Display this alternative layout
grouped_nested
