 # read in the results from /results/all_fires/treatment_effects_all_fires_final.csv
    results <- read.csv("results/all_fires/treatment_effects_all_fires_final.csv")
    
    # make a box plot with treatment effect on the y , with CI_lower and CI_upper as the error bars and Model on the x
    library(ggplot2)
    p <- ggplot(results, aes(x = Model, y = Treatment_Effect)) +
        geom_boxplot() +
        geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
        theme_minimal() +
        labs(title = "Treatment Effects by Model",
                x = "Model",
                y = "Treatment Effect") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))    
    print(p)

# format boxplot to have colour within the bars
    p <- ggplot(results, aes(x = Model, y = Treatment_Effect, fill = Model)) +
        geom_boxplot() +
        geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
        theme_minimal() +
        labs(title = "Treatment Effects by Model",
                x = "Model",
                y = "Treatment Effect") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_brewer(palette = "Set3")    
    print(p)
