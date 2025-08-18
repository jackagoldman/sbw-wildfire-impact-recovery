# Figures script for SBW wildfire impact and recovery project

# load required packages
library(tidyverse)
library(sf)
library(sp)
library(gstat)
library(ggplot2)
library(car)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(cowplot)

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


#utils
source(file.path(base_path, "/src/utils.R"))


# Load  data 
source(file.path(base_path,"/src/load_data.R"))

#load
source(file.path(base_path, "/src/best_model_functions.R"))

# Load Canada map data
canada <- ne_countries(country = "Canada", scale = "medium", returnclass = "sf")

# Create a small inset map of Canada with Ontario highlighted
canada_inset <- ggplot() +
  # Plot all of Canada in grey
  geom_sf(data = canada, fill = "grey80", color = "grey50") +
  # Highlight Ontario in navajowhite1
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  # Remove all unnecessary elements
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
  )

# Load Ontario map data
ontario <- ne_states(country = "Canada", returnclass = "sf") %>%
  filter(name == "Ontario")

# Create the inset pie chart for m.data.sf
inset_data1 <- h90.sf %>%
  group_by(history) %>%
  summarise(count = n())

inset_pie1 <- ggplot(inset_data1, aes(x = "", y = count, fill = factor(history))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 8)  # Adjust the size here


# Create the first map with inset
map1 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = h90.sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                     name = NULL, 
                     labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                     guide = guide_legend(title.position = "top", 
                                        title.hjust = 0.4,
                                         key.size = 3,        # Increase key size (default is ~1)
                                         override.aes = list(size = 4))) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = c(0.75, 0.85),  # Move legend further left (x=0.15 instead of right)
    legend.direction = "vertical",
    legend.text = element_text(size = 16),
    legend.justification = c(0, 1),   # Anchor at top-left
  ) +
  #labs(title = "A)") +
  annotation_custom(ggplotGrob(inset_pie1), xmin = -79, xmax = -75, ymin = 48, ymax = 58) +
   # Add the Canada inset map in the bottom left corner
  annotation_custom(ggplotGrob(canada_inset), xmin = -95, xmax = -87, ymin = 42, ymax = 47)


print(map1)

# ggsave
ggsave(
  filename = file.path(base_path, "plots/maps/figure1_study_area_map_matched_pairs_all_fires.png"),
  plot = map1,
  width = 10, height = 8, dpi = 300
)
###### Severity and Recovery ==========


# Create the inset pie chart for m.data.sf2 - severity
inset_data2 <- m.data.sev_sf %>%
  group_by(history) %>%
  summarise(count = n())

inset_pie2 <- ggplot(inset_data2, aes(x = "", y = count, fill = factor(history))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 6)  # Adjust the size here


# Create the second map with inset using m.data.sf2
severity_map <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = m.data.sev_sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                     name = NULL, 
                     labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                     guide = guide_legend(title.position = "top", 
                                        title.hjust = 0.4,
                                         key.size = 3,        # Increase key size (default is ~1)
                                         override.aes = list(size = 4))) + 
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) + labs(title = "A)") +
  annotation_custom(ggplotGrob(inset_pie2), xmin = -78, xmax = -73, ymin = 42, ymax = 52)

print(severity_map)


# Create the inset pie chart for m.data.rec_sf - recovery
inset_data3 <- m.data.rec_sf %>%
  group_by(history) %>%
  summarise(count = n())

inset_pie3 <- ggplot(inset_data3, aes(x = "", y = count, fill = factor(history))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 6)  # Adjust the size here

# recovery map using m.data.rec_sf
recovery_map <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = m.data.rec_sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                     name = NULL, 
                     labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                     guide = guide_legend(title.position = "top", title.hjust = 0.4)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) +
  labs(title = "B)") +
  annotation_custom(ggplotGrob(inset_pie3), xmin = -78, xmax = -73, ymin = 42, ymax = 52)

print(recovery_map)

# combine the two maps into one figure
# A completely different approach to create a combined plot with legend
library(gridExtra)  # Make sure this is loaded

# Set up individual maps without legends
severity_map_no_legend <- severity_map + theme(legend.position = "none")
recovery_map_no_legend <- recovery_map + theme(legend.position = "none")


#
recovery_map_with_legend <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = m.data.rec_sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                    guide = guide_legend(title.position = "top", 
                                        title.hjust = 0.4,
                                         key.size = 3,        # Increase key size (default is ~1)
                                         override.aes = list(size = 4))) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = c(0.2, 0.2),  # Position legend inside the plot at bottom left
    legend.direction = "horizontal",
    legend.text = element_text(size = 12),
    legend.margin = margin(6, 6, 6, 6)
  ) +
  labs(title = "B)") +
  annotation_custom(ggplotGrob(inset_pie3), xmin = -78, xmax = -73, ymin = 42, ymax = 52)



# Arrange both maps side by side without separate legend
combined_plot_with_internal_legends <- gridExtra::grid.arrange(
  severity_map_no_legend, 
  recovery_map_with_legend,
  ncol = 2
)

print(combined_plot_with_internal_legends)

# Save this version
ggsave(
  filename = file.path(base_path, "plots/maps/figure2_study_area_severity_recovery_maps.png"),
  plot = combined_plot_with_internal_legends,
  width = 12, height = 8, dpi = 300
)
