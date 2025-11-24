library(ggplot2)
library(dplyr)
library(paletteer)
library(MatchIt)

# Load the data
data_sev <- read.csv("data/on_sev_match_data.csv")
data_rec <- read.csv("data/on_rec_match_data.csv")

# Create a factor for history with labels for both datasets
data_sev$history_label <- factor(data_sev$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))
data_rec$history_label <- factor(data_rec$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

# Load the best models for subgroups 1-4 severity
best_model_sev_1 <- readRDS("results/subgroup/best_model_subgroup1_severity.RDS")
best_model_sev_2 <- readRDS("results/subgroup/best_model_subgroup2_severity.RDS")
best_model_sev_3 <- readRDS("results/subgroup/best_model_subgroup3_severity.RDS")
best_model_sev_4 <- readRDS("results/subgroup/best_model_subgroup4_severity.RDS")

# Extract matched data
m.data.sev_1 <- match_data(best_model_sev_1)
m.data.sev_2 <- match_data(best_model_sev_2)
m.data.sev_3 <- match_data(best_model_sev_3)
m.data.sev_4 <- match_data(best_model_sev_4)

# Add subgroup column and history_label
m.data.sev_1$subgroup <- "Subgroup 1 (0-2 years)"
m.data.sev_1$history_label <- factor(m.data.sev_1$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

m.data.sev_2$subgroup <- "Subgroup 2 (3-5 years)"
m.data.sev_2$history_label <- factor(m.data.sev_2$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

m.data.sev_3$subgroup <- "Subgroup 3 (6-9 years)"
m.data.sev_3$history_label <- factor(m.data.sev_3$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

m.data.sev_4$subgroup <- "Subgroup 4 (10+ years)"
m.data.sev_4$history_label <- factor(m.data.sev_4$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

# Create the box plot for rbr_w_offset using data_sev
p <- ggplot(data_sev, aes(x = history_label, y = rbr_w_offset, fill = history_label)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Defoliated fires" = "#E84A5FFF", "Non-defoliated fires" = "#9EBCDAFF")) +
  labs(x = "History", y = "Median Burn Severity (RBR)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = "none")

# Save the plot
ggsave("plots/boxplot_rbr_offset.png", p, width = 8, height = 6)

# Display the plot
print(p)

# Create the box plot for Average.Recovery using data_rec
p2 <- ggplot(data_rec, aes(x = history_label, y = Average.Recovery, fill = history_label)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Defoliated fires" = "#E84A5FFF", "Non-defoliated fires" = "#9EBCDAFF")) +
  labs(x = "History", y = "Median Burn Severity (RBR)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = "none")

# Save the plot
ggsave("plots/boxplot_average_recovery.png", p2, width = 8, height = 6)

# Display the plot
print(p2)

# Filter defoliated fires
data_sev_defol <- data_sev %>% filter(history == 1)

# Create the box plot for rbr_w_offset by window_opp for defoliated fires
p3 <- ggplot(data_sev_defol, aes(x = as.factor(window_opp), y = rbr_w_offset, fill = as.factor(window_opp))) +
  geom_boxplot() +
  scale_fill_paletteer_d("MetBrewer::OKeeffe1") +
  labs(x = "Time Since Defoliation", y = "Median Burn Severity (RBR)") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = "none")

# Save the plot
ggsave("plots/boxplot_rbr_offset_defol_by_window_opp.png", p3, width = 8, height = 6)

# Display the plot
print(p3)

# Combine data
combined_data <- bind_rows(m.data.sev_1, m.data.sev_2, m.data.sev_3, m.data.sev_4)

# Create a factor for history with labels
combined_data$history_label <- factor(combined_data$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

# Calculate max y per subgroup for annotations
max_y_data <- combined_data %>%
  group_by(subgroup) %>%
  summarise(max_y = max(rbr_w_offset))

# Create annotation data for subgroups 1 and 4
annotation_data <- max_y_data %>%
  filter(subgroup %in% c("Subgroup 1 (0-2 years)", "Subgroup 4 (10+ years)")) %>%
  mutate(x = 1.5, y = max_y * 0.75, label = "***")

# Create the single plot with 4 panels in 2x2 grid
p4 <- ggplot(combined_data, aes(x = history_label, y = rbr_w_offset, fill = history_label)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_text(data = annotation_data, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 6) +
  geom_segment(data = annotation_data, aes(x = 1.2, xend = 1.8, y = y - 5, yend = y - 5), inherit.aes = FALSE) +
  facet_wrap(~subgroup, ncol = 2) +
  scale_fill_manual(values = c("Defoliated fires" = "#E84A5FFF", "Non-defoliated fires" = "#9EBCDAFF")) +
  labs(x = "Defoliation History", y = "Median Burn Severity (RBR)", fill = "History") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = "top")

# Save the plot
ggsave("plots/boxplot_severity_by_history_all_subgroups.png", p4, width = 8, height = 10)

# Load the best models for subgroups 1-4 recovery
best_model_rec_1 <- readRDS("results/subgroup/best_model_subgroup1_recovery.RDS")
best_model_rec_2 <- readRDS("results/subgroup/best_model_subgroup2_recovery.RDS")
best_model_rec_3 <- readRDS("results/subgroup/best_model_subgroup3_recovery.RDS")
best_model_rec_4 <- readRDS("results/subgroup/best_model_subgroup4_recovery.RDS")

# Extract matched data
m.data.rec_1 <- match_data(best_model_rec_1)
m.data.rec_2 <- match_data(best_model_rec_2)
m.data.rec_3 <- match_data(best_model_rec_3)
m.data.rec_4 <- match_data(best_model_rec_4)

# Add subgroup column and history_label
m.data.rec_1$subgroup <- "Subgroup 1 (0-2 years)"
m.data.rec_1$history_label <- factor(m.data.rec_1$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

m.data.rec_2$subgroup <- "Subgroup 2 (3-5 years)"
m.data.rec_2$history_label <- factor(m.data.rec_2$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

m.data.rec_3$subgroup <- "Subgroup 3 (6-9 years)"
m.data.rec_3$history_label <- factor(m.data.rec_3$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

m.data.rec_4$subgroup <- "Subgroup 4 (10+ years)"
m.data.rec_4$history_label <- factor(m.data.rec_4$history, levels = c(1, 0), labels = c("Defoliated fires", "Non-defoliated fires"))

# Combine data for recovery
combined_data_rec <- bind_rows(m.data.rec_1, m.data.rec_2, m.data.rec_3, m.data.rec_4)

# Calculate max y per subgroup for annotations (for recovery)
max_y_data_rec <- combined_data_rec %>%
  group_by(subgroup) %>%
  summarise(max_y = max(recovery))

# Create annotation data for subgroups 1, 3 and 4
annotation_data_rec <- max_y_data_rec %>%
  filter(subgroup %in% c("Subgroup 1 (0-2 years)", "Subgroup 3 (6-9 years)", "Subgroup 4 (10+ years)")) %>%
  mutate(x = 1.5, y = ifelse(subgroup %in% c("Subgroup 3 (6-9 years)", "Subgroup 4 (10+ years)"), 200, max_y * 0.75), 
         label = ifelse(subgroup == "Subgroup 3 (6-9 years)", "*", "***"))

# Create the recovery plot
p5 <- ggplot(combined_data_rec, aes(x = history_label, y = recovery, fill = history_label)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_text(data = annotation_data_rec, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 6) +
  geom_segment(data = annotation_data_rec, aes(x = 1.2, xend = 1.8, y = y - 5, yend = y - 5), inherit.aes = FALSE) +
  facet_wrap(~subgroup, ncol = 2) +
  scale_fill_manual(values = c("Defoliated fires" = "#E84A5FFF", "Non-defoliated fires" = "#9EBCDAFF")) +
  labs(x = "Defoliation History", y = "Recovery", fill = "History") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = "top")

# Save the plot
ggsave("plots/boxplot_recovery_by_history_all_subgroups.png", p5, width = 8, height = 10)

# Display the plot
print(p5)