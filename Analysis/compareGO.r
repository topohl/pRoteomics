#' Compare Enrichment Analysis and Visualization
#'
#' @description
#' This script performs a comparative enrichment analysis by:
#'
#' - Reading multiple CSV files from a specified directory containing enrichment data.  
#'   The CSV filenames (excluding the extension) are used as identifiers for different comparisons.
#'
#' - Combining the individual CSV data into a single data frame with an additional 
#'   "Comparison" column indicating the source file.
#'
#' - Selecting the top enriched terms per comparison based on the highest absolute 
#'   Normalized Enrichment Score (NES), and compiling a master list of these top terms.
#'
#' - Filtering the combined data to include only these top terms, across all comparisons.
#'
#' - Reordering comparisons in the final visualization based on the maximum absolute NES.
#'
#' - Generating a significance label for each term when the adjusted p-value (p.adjust) is less than 0.05.
#'
#' @details
#' The script creates two types of visualizations:
#'
#' 1. A professional heatmap where:
#' 2. A point plot that visualizes:
#' 
#'
#' @section File Inputs:
#' CSV files are read from the directory "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment"
#' based on the pattern "*.csv". Each file serves as input for one comparison.
#'
#' @section Output:
#' The script outputs two plots:
#' \enumerate{
#'   \item A heatmap visualizing the top differential enrichments.
#'   \item A point plot emphasizing the size (via -log10(p.adjust)) and color (via NES) of each term.
#' }
#'
#' @note 
#' Ensure that the file paths and package dependencies are properly set up before running the script.
#'
#' @author Tobias Pohl

# -----------------------------------------------------
# Load Libraries
# -----------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, stringr, ggpubr, ggthemes, dplyr)

# -----------------------------------------------------
# Define Paths and Project Directories
# -----------------------------------------------------

file_paths <- list.files(path = "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment", pattern = "*.csv", full.names = TRUE)
names(file_paths) <- basename(file_paths) %>% str_remove(".csv")

# -----------------------------------------------------
# Read and Combine Data
# -----------------------------------------------------

enrichment_list <- lapply(file_paths, read.csv)
names(enrichment_list) <- names(file_paths)

combined_df <- bind_rows(
  lapply(names(enrichment_list), function(name) {
    df <- enrichment_list[[name]]
    df$Comparison <- name
    return(df)
  })
)

# -----------------------------------------------------
# Select Top Terms
# -----------------------------------------------------

# For each comparison, select the top 10 terms with the highest absolute NES value
top10_df <- combined_df %>%
    group_by(Comparison) %>%
    slice_max(order_by = abs(NES), n = 15) %>% 
    ungroup()

# Create a master list of top terms from all comparisons
top_terms <- unique(top10_df$Description)

# -----------------------------------------------------
# Filter Data for Heatmap
# -----------------------------------------------------

# For every comparison, get the rows corresponding to these top terms, even if they're not in the top 10 there
lookup_df <- combined_df %>%
    filter(Description %in% top_terms)

# Compute an ordering of comparisons by the maximum absolute NES value in the lookup dataframe
comparison_order <- lookup_df %>%
    group_by(Comparison) %>%
    summarize(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_NES)) %>%
    pull(Comparison)

# Reorder the Comparison factor using the computed order
lookup_df <- lookup_df %>%
    mutate(Comparison = factor(Comparison, levels = comparison_order))

# -----------------------------------------------------
# Create Heatmap
# -----------------------------------------------------

# Create a significance label if a p-value adjustment column is available (using p.adjust below).
# Adjust the threshold as necessary.
lookup_df <- lookup_df %>%
    mutate(sig_label = ifelse(p.adjust < 0.05, "*", ""))
    # Create a professional heatmap with Description on y-axis and Comparison on x-axis.
    # Significant changes are marked with a star.
    ggplot(lookup_df, aes(x = Comparison, y = reorder(Description, NES), fill = NES)) +
        geom_tile(color = "white", size = 0.5) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        geom_text(aes(label = sig_label), color = "black", size = 5, vjust = 0.5) +
        labs(
            title = "Heatmap of Top Differential Enrichment",
            x = "Comparison",
            y = "Enriched Terms",
            fill = "Normalized Enrichment Score (NES)"
        ) +
        theme_minimal(base_family = "Arial", base_size = 12) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(face = "bold", size = 10),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
            plot.margin = margin(10, 10, 10, 10)
        )

# -----------------------------------------------------
# Create Dot Plot
# -----------------------------------------------------

ggplot(lookup_df, aes(x = Comparison, y = reorder(Description, NES), color = NES, size = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_pubclean()