
library(dplyr)

data <- read.csv("data/on_rec_match_data.csv")

fires <- data %>% filter(subclass >= 1 & subclass <= 20)


# get the names of the fires from Fire_ID
fires <- fires %>% filter(history == 1)
fire_names <- unique(fires$Fire_ID)

(fire_names)

subclass = data |> select(Fire_ID, subclass)
# read in file from /Users/jgoldman/Desktop/thesis_presentation/recovery_trajectories.csv
recovery_trajectories <- read.csv("/Users/jgoldman/Desktop/thesis_presentation/recovery_trajectories_all.csv")

library(tidyr)

recovery_trajectories <- recovery_trajectories %>%
    pivot_longer(
        cols = c(layer_1_mean:layer_13_mean, layer_1_sd:layer_13_sd),
        names_to = c("layer", ".value"),
        names_pattern = "layer_(\\d+)_(mean|sd)"
    ) %>%
    rename(nbr = mean)
recovery_trajectories <- recovery_trajectories %>%
    pivot_longer(cols = layer_1_mean:layer_13_mean, names_to = "layer", names_pattern = "layer_(\\d+)_mean", values_to = "nbr")

    recovery_trajectories <- recovery_trajectories %>% rename(year = layer)

        library(ggplot2)

        # Assuming recovery_trajectories has columns: Fire_ID, treatment, year, nbr
        # Plot time series for each Fire_ID, with lines colored by treatment
    ggplot(recovery_trajectories, aes(x = as.numeric(year), y = nbr, color = Treatment, group = interaction(Fire_ID, Treatment))) +
            geom_line() +
            facet_wrap(~ Fire_ID) +
            labs(x = "Time (Year)", y = "NBR", title = "Recovery Trajectories by Fire_ID") +
            theme_minimal()


    average_trajectories <- recovery_trajectories %>%
                group_by(Treatment, year) %>%
                summarise(avg_nbr = mean(nbr, na.rm = TRUE))

    average_trajectories <- average_trajectories %>%
                    mutate(Treatment = factor(Treatment, levels = c(TRUE, FALSE), labels = c("Defoliated", "Non-defoliated")))


                        # Define custom colors for treatments
                        treatment_colors <- c("Non-defoliated" =  "#FF8C00A0", "Defoliated" = "#8B0000A0")


                     
                            # Create the publication-quality plot
                            ggplot(average_trajectories, aes(x = as.numeric(year), y = avg_nbr, color = Treatment, group = Treatment)) +
                                geom_line(size = 1.5, alpha = 0.8) +  # Thicker lines with slight transparency
                                geom_point(size = 10, alpha = 0.7, aes(color = Treatment), stroke = 1.2) +  # Bigger points with alpha, bold stroke
                                scale_color_manual(values = treatment_colors) +
                                scale_x_continuous(breaks = 1:13) +
                                labs(x = "Year", y = "Mean Normalized Burn Ratio (NBR)", color = "Defoliation History") +
                                theme_bw(base_size = 16) +  # Larger base font size for publication quality
                                theme(
                                    panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2),
                                    panel.grid = element_blank(),  # Remove gridlines
                                    axis.title = element_text(size = 18, face = "bold"),
                                    axis.text = element_text(size = 14),
                                    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                                    legend.position = c(0.8, 0.2),  # Legend inside the plot frame
                                    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
                                    legend.title = element_text(size = 14, face = "bold"),
                                    legend.text = element_text(size = 12)
                                )

    ggsave("plots/average_recovery_trajectories_pub.png", width = 10, height = 6, dpi = 300)
    library(ggplot2)

    # Assuming recovery_trajectories has columns: Fire_ID, treatment, layer, nbr
    # Plot time series for each Fire_ID, with lines colored by treatment
ggplot(recovery_trajectories, aes(x = as.numeric(layer), y = nbr, color = Treatment, group = interaction(Fire_ID, Treatment))) +
        geom_line() +
        facet_wrap(~ Fire_ID) +
        labs(x = "Time (Layer)", y = "NBR", title = "Recovery Trajectories by Fire_ID") +
        theme_minimal()


average_trajectories <- recovery_trajectories %>%
            group_by(Treatment, layer) %>%
            summarise(avg_nbr = mean(nbr, na.rm = TRUE))

average_trajectories <- average_trajectories %>%
                mutate(Treatment = factor(Treatment, levels = c(TRUE, FALSE), labels = c("Defoliated", "Non-defoliated")))


                    # Define custom colors for treatments
                    treatment_colors <- c("Non-defoliated" = "#9EBCDAFF", "Defoliated" = "#E84A5FFF")

                 
                        # Create the publication-quality plot
                        ggplot(average_trajectories, aes(x = as.numeric(layer), y = avg_nbr, color = Treatment, group = Treatment)) +
                            geom_line(size = 1.5, alpha = 0.8) +  # Thicker lines with slight transparency
                            geom_point(size = 10, alpha = 0.7, aes(color = Treatment), stroke = 1.2) +  # Bigger points with alpha, bold stroke
                            scale_color_manual(values = treatment_colors) +
                            scale_x_continuous(breaks = 1:13) +
                            labs(x = "Layer", y = "Average NBR", title = "Average Recovery Trajectories by Treatment n=40") +
                            theme_minimal(base_size = 16) +  # Larger base font size for publication quality
                            theme(
                                panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2),
                                panel.grid = element_blank(),  # Remove gridlines
                                axis.title = element_text(size = 18, face = "bold"),
                                axis.text = element_text(size = 14),
                                plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                                legend.position = c(0.8, 0.2),  # Legend inside the plot frame
                                legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
                                legend.title = element_text(size = 14, face = "bold"),
                                legend.text = element_text(size = 12)
                            )

ggsave("plots/average_recovery_trajectories.png", width = 10, height = 6, dpi = 300)



library(dplyr)
library(ggplot2)
library(rlang)


library(dplyr)
library(ggplot2)
library(tidyr)

# recovery_trajectories: columns = Treatment (0/1), subclass, year, nbr

#combine subclass to recovery_trajectories
recovery_trajectories <- recovery_trajectories %>%
  left_join(subclass, by = "Fire_ID")

# 1) Build ATT-consistent weights from subclass composition

subclass_counts <- recovery_trajectories %>%
  group_by(subclass) %>%
  summarise(n_treated  = sum(Treatment == 1),
            n_control  = sum(Treatment == 0),
            .groups = "drop")

weights_df <- recovery_trajectories %>%
  left_join(subclass_counts, by = "subclass") %>%
  mutate(w_att = ifelse(Treatment == 1, 1,
                        ifelse(n_control > 0, n_treated / n_control, 0)))
# Note: If some subclasses have only treated or only control, handle sensibly:
# - If n_control==0 and Treatment==0 -> weight 0 (excluded).
# - If n_treated==0 and Treatment==1 -> (rare) weight 0 or drop subclass.

# 2) Compute weighted means per year and per treatment
att_lines <- weights_df %>%
  group_by(year, Treatment) %>%
  summarise(wmean = if (all(is.na(nbr))) NA_real_
            else weighted.mean(nbr, w_att, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Treatment = factor(Treatment, levels = c(1,0),
                            labels = c("Defoliated","Non-defoliated")))

# 3) Plot two lines (ATT-consistent)
treatment_colors <- c("Non-defoliated" =  "#FF8C00A0",
                      "Defoliated"      = "#8B0000A0")

ggplot(att_lines, aes(x = as.numeric(year),
                      y = wmean,
                      color = Treatment,
                      group = Treatment)) +
  geom_line(size = 1.5, alpha = 0.9) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = treatment_colors, name = "Defoliation History") +
  scale_x_continuous(breaks = sort(unique(as.numeric(att_lines$year)))) +
  labs(x = "Year", y = "Weighted mean NBR",
       title = "Two-line ATT-consistent plot (subclass weights)") +
  theme_bw(base_size = 16) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2),
    panel.grid = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text  = element_text(size = 14),
    legend.position = c(0.82, 0.2),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  )
