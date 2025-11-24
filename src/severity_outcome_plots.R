# Script to create violin/box plots for severity outcome distributions
# between treated and control groups
#Compare outcome distributions between treated and matched control groups

library(tidyverse)
library(ggplot2)

# Read the matched data
matched_data <- read.csv("data/on_sev_match_data.csv")

# Check the data
head(matched_data)

# Create violin plot with boxplot overlay for rbr_w_offset by history
plot <- ggplot(matched_data, aes(x = factor(history), y = rbr_w_offset, fill = factor(history))) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  labs(
    title = "All Wildfires",
    x = NULL,
    y = "Relative Burn Ratio (RBR)",
    fill = "History"
  ) +
  scale_x_discrete(labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  coord_cartesian(ylim = c(-100, 700)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "none"
  ) +  geom_text(label = "**", x = 1.5, y = 500, size = 6, color = "black") +
      geom_segment(x = 1.2, xend = 1.8, y = 490, yend = 490, color = "black", size = 1)
print(plot)

# Save the plot
ggsave(
  filename = "plots/severity_outcome_violin.png",
  plot = plot,
  width = 14,
  height = 20,
  dpi = 300
)

# Read the subgroup models and create plots
subgroup_files <- c(
  "results/subgroup/fit_model_subgroup1_severity.RDS",
  "results/subgroup/fit_model_subgroup2_severity.RDS",
  "results/subgroup/fit_model_subgroup3_severity.RDS",
  "results/subgroup/fit_model_subgroup4_severity.RDS"
)

subgroup_plots <- list()

for (i in 1:length(subgroup_files)) {
  model <- readRDS(subgroup_files[i])
  subgroup_data <- model.frame(model)  # Ensure it's a data frame
  
  # Calculate and print ymin and ymax for each subgroup
  ymin <- min(subgroup_data$rbr_w_offset, na.rm = TRUE)
  ymax <- max(subgroup_data$rbr_w_offset, na.rm = TRUE)
  print(paste("Subgroup", i, "ymin:", ymin, "ymax:", ymax))
  
  # Create violin plot
  p <- ggplot(subgroup_data, aes(x = factor(history), y = rbr_w_offset, fill = factor(history))) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
    labs(
      title = paste("Subgroup", i),
      x = "Defoliation History",
      y = "Relative Burn Ratio (RBR)",
      fill = "History"
    ) +
    scale_x_discrete(labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
    scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                      labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none"  # Remove legend for panels
    )
  
  # Add annotation for subgroup 1 and 4
  if (i %in% c(1, 4)) {
    p <- p + 
      geom_text(label = "***", x = 1.5, y = 500, size = 6, color = "red") +
      geom_segment(x = 0.6, xend = 2.4, y = 490, yend = 490, color = "red", size = 1)
  }
  
  subgroup_plots[[i]] <- p
}

# remake the subset plot but with ymin = -100 and ymax = 700
for (i in 1:length(subgroup_files)) {
  model <- readRDS(subgroup_files[i])
  subgroup_data <- model.frame(model)  # Ensure it's a data frame

  if (i == 1) {
    title_name <- paste("Subgroup", "(0-2 years)")
  } else if (i == 2) {
      title_name <- paste("Subgroup", "(3-5 years)")
    } else if (i == 3) {
      title_name <- paste("Subgroup", "(6-9 years)")
    } else if (i == 4) {
      title_name <- paste("Subgroup", "(10-15 years)")
  }
  
  # Create violin plot with fixed y-axis limits
  p <- ggplot(subgroup_data, aes(x = factor(history), y = rbr_w_offset, fill = factor(history))) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
    labs(
      title = title_name,
      x = NULL,
      y = "Relative Burn Ratio (RBR)",
      fill = "History"
    ) +
    scale_x_discrete(labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
    scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                      labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
    coord_cartesian(ylim = c(-100, 700)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      legend.position = "none"  # Remove legend for panels
    )
  
  # Add annotation for subgroup 1 and 4
  if (i %in% c(1, 4)) {
    p <- p + 
      geom_text(label = "***", x = 1.5, y = 500, size = 6, color = "black") +
      geom_segment(x = 1.2, xend = 1.8, y = 490, yend = 490, color = "black", size = 1)
  }
  
  subgroup_plots[[i]] <- p
}

# Arrange in 2x2 grid
library(cowplot)
combined_plot <- plot_grid(plotlist = subgroup_plots, ncol = 2, nrow = 2, align = 'v')
print(combined_plot)




# Save the combined plot
ggsave(
  filename = "plots/severity_subgroup_violins.png",
  plot = combined_plot,
  width = 14,
  height = 20,
  dpi = 300
)

# Print the combined plot
print(combined_plot)

# Print summary statistics
summary_stats <- matched_data %>%
  group_by(history) %>%
  summarise(
    n = n(),
    mean_rbr = mean(rbr_w_offset, na.rm = TRUE),
    median_rbr = median(rbr_w_offset, na.rm = TRUE),
    sd_rbr = sd(rbr_w_offset, na.rm = TRUE)
  )

print(summary_stats)
#save summary states to results/subgroup/severity_summary_stats.csv
write.csv(summary_stats, "results/subgroup/severity_summary_stats.csv", row.names = FALSE)

# make a plot for intermediate subgroup severity outcome
int = readRDS("results/subgroup/fit_model_intermediate_severity.RDS")
int_data <- model.frame(int)

plot_int <- ggplot(int_data, aes(x = factor(history), y = rbr_w_offset, fill = factor(history))) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  labs(
    title = "Intermediate Subgroup (3-9 years)",
    x = NULL,
    y = "Relative Burn Ratio (RBR)",
    fill = "History"
  ) +
  scale_x_discrete(labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  coord_cartesian(ylim = c(-100, 700)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "none"  # Remove legend for panels
  ) 
print(plot_int)

# Save the plot
ggsave(
  filename = "plots/severity_intermediate_subgroup_violin.png",  
  plot = plot_int,
  width = 14,
  height = 20,
  dpi = 300
)

# make a plot that compbines all 6 severity plots
# put the first plot names "plot" in the top left, then subgroup_plotlist, then plot_int in the bottom right
# add a) or similar labels to each subplot
# Modify y-axis for second column plots
subgroup_plots[[2]] <- subgroup_plots[[2]] + labs(y = NULL)
subgroup_plots[[4]] <- subgroup_plots[[4]] + labs(y = NULL)
plot_int <- plot_int + labs(y = NULL)

combined_plot_all <- plot_grid(plot, plotlist = subgroup_plots, plot_int, ncol = 2, nrow = 3, align = 'v', labels = c("a)", "b)", "c)", "d)", "e)", "f)"), label_size = 20)
print(combined_plot_all)

# save the combined plot to plots/severity_all_violins.png
ggsave(
  filename = "plots/severity_all_violins.png",
  plot = combined_plot_all,
  width = 16,
  height = 22,
  dpi = 300
)

########## RECOVERY OUTCOME PLOTS ##########

matched_data_rec <- read.csv("data/on_rec_match_data.csv")

# Check the data
head(matched_data_rec)

# Create violin plot with boxplot overlay for recovery by history
max_recovery <- max(matched_data_rec$recovery, na.rm = TRUE)
plot_rec <- ggplot(matched_data_rec, aes(x = factor(history), y = recovery, fill = factor(history))) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  labs(
    title = "All Wildfires",
    x = NULL,
    y = "Recovery Magnitude (%)",
    fill = "History"
  ) +
  scale_x_discrete(labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "none"
  ) + geom_text(label = "***", x = 1.5, y = max_recovery * 1.05, size = 6, color = "black") +
      geom_segment(x = 1.2, xend = 1.8, y = max_recovery * 1.03, yend = max_recovery * 1.03, color = "black", size = 1)


print(plot_rec)
# Save the plot
ggsave(
  filename = "plots/recovery_outcome_violin.png",
  plot = plot_rec,
  width = 14,
  height = 20,
  dpi = 300
)

# Print summary statistics for recovery
summary_stats_rec <- matched_data_rec %>%
  group_by(history) %>%
  summarise(
    n = n(),
    mean_recovery = mean(recovery, na.rm = TRUE),
    median_recovery = median(recovery, na.rm = TRUE),
    sd_recovery = sd(recovery, na.rm = TRUE)
  )

print(summary_stats_rec)

# make the subgroup plot for recovery outcome
subgroup_files_rec <- c(
  "results/subgroup/fit_model_subgroup1_recovery.RDS",
  "results/subgroup/fit_model_subgroup2_recovery.RDS",
  "results/subgroup/fit_model_subgroup3_recovery.RDS",
  "results/subgroup/fit_model_subgroup4_recovery.RDS"
) 
subgroup_plots_rec <- list()

for (i in 1:length(subgroup_files_rec)) {
  model <- readRDS(subgroup_files_rec[i])
  subgroup_data <- model.frame(model)  # Ensure it's a data frame

  if (i == 1) {
    title_name <- paste("Subgroup", "(0-2 years)")
  } else if (i == 2) {
      title_name <- paste("Subgroup", "(3-5 years)")
    } else if (i == 3) {
      title_name <- paste("Subgroup", "(6-9 years)")
    } else if (i == 4) {
      title_name <- paste("Subgroup", "(10-15 years)")
  }
  
  # Create violin plot with fixed y-axis limits
  p_rec <- ggplot(subgroup_data, aes(x = factor(history), y = recovery, fill = factor(history))) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
    labs(
      title = title_name,
      x = NULL,
      y = "Recovery Magnitude (%)",
      fill = "History"
    ) +
    scale_x_discrete(labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
    scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                      labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
    coord_cartesian(ylim = c(-100, 700)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      legend.position = "none"  # Remove legend for panels
    )
  
  # Add annotation for subgroup 1 and 4
  if (i %in% c(1, 4)) {
    p_rec <- p_rec + 
      geom_text(label = "***", x = 1.5, y = 500, size = 6, color = "black") +
      geom_segment(x = 1.2, xend = 1.8, y = 490, yend = 490, color = "black", size = 1)
  } else if (i == 3) {
     p_rec <- p_rec + 
      geom_text(label = "*", x = 1.5, y = 500, size = 6, color = "black") +
      geom_segment(x = 1.2, xend = 1.8, y = 490, yend = 490, color = "black", size = 1)
  }
  
  subgroup_plots_rec[[i]] <- p_rec
}

# Arrange in 2x2 grid
combined_plot_rec <- plot_grid(plotlist = subgroup_plots_rec, ncol = 2, nrow = 2, align = 'v')
print(combined_plot_rec)

# save the combined plot to plots/recovery_subgroup_violins.png
ggsave(
  filename = "plots/recovery_subgroup_violins.png",
  plot = combined_plot_rec,
  width = 14,
  height = 20,
  dpi = 300
)

# make a plot for intermediate subgroup severity outcome
int_rec = readRDS("results/subgroup/fit_model_intermediate_recovery.RDS")
int_data_rec <- model.frame(int_rec)

plot_int_rec <- ggplot(int_data_rec, aes(x = factor(history), y = recovery, fill = factor(history))) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  labs(
    title = "Intermediate Subgroup (3-9 years)",
    x = NULL,
    y = "Magnitude of Recovery (%)",
    fill = "History"
  ) +
  scale_x_discrete(labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  coord_cartesian(ylim = c(-100, 700)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "none"  # Remove legend for panels
  ) +  geom_text(label = "***", x = 1.5, y = max_recovery * 1.05, size = 6, color = "black") +
      geom_segment(x = 1.2, xend = 1.8, y = max_recovery * 1.03, yend = max_recovery * 1.03, color = "black", size = 1)



print(plot_int_rec)

# make a plot that compbines all 6 recovery plots
# put the first plot names "plot" in the top left, then subgroup_plotlist, then plot_int in the bottom right
# add a) or similar labels to each subplot
# Modify y-axis for second column plots
subgroup_plots_rec[[2]] <- subgroup_plots_rec[[2]] + labs(y = NULL)
subgroup_plots_rec[[4]] <- subgroup_plots_rec[[4]] + labs(y = NULL)
plot_int_rec <- plot_int_rec + labs(y = NULL)

combined_plot_all <- plot_grid(plot_rec, plotlist = subgroup_plots_rec, plot_int_rec, ncol = 2, nrow = 3, align = 'v', labels = c("a)", "b)", "c)", "d)", "e)", "f)"), label_size = 20)
print(combined_plot_all)

# gg save the combined plot
ggsave(
  filename = "plots/recovery_all_violins.png",
  plot = combined_plot_all,
  width = 16, 
  height = 24,
  dpi = 300
)
