# Script to create a forest plot for ATE estimates with confidence intervals

library(tidyverse)
library(ggplot2)

# Read the subgroup treatment effects data
treatment_effects <- read.csv("results/subgroup/all_treatment_effects.csv")

# Filter for severity models (assuming there's a column indicating model type)
# If the file has a 'Model' column, filter for 'Severity'
if ("Model" %in% colnames(treatment_effects)) {
  treatment_effects <- treatment_effects %>% filter(Model == "Severity")
}

# Check the data
print(head(treatment_effects))
print(colnames(treatment_effects))

# Create forest plot
# Assuming columns: Subgroup, estimate, conf.low, conf.high, p.value
forest_plot <- ggplot(treatment_effects, aes(x = Estimate, y = Subgroup)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbarh(aes(xmin = Conf_Low, xmax = Conf_High), height = 0.2, color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Forest Plot of Severity ATE across Subgroups",
    x = "Treatment Effect (ATE)",
    y = "Subgroup"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Add p-value annotations
forest_plot <- forest_plot +
  geom_text(
    aes(label = sprintf("p = %.3f", P_Value)),
    hjust = -0.1,
    vjust = -1,
    size = 3
  )

# Save the plot
ggsave(
  filename = "plots/forest_plot_severity_subgroups.png",
  plot = forest_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Print the plot
print(forest_plot)