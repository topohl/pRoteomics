#' compareGO.r - Comparative Gene Ontology Enrichment Analysis and Visualization
#'
#' @description
#' This script performs a comparative analysis of GO enrichment across multiple experiments.
#' It:
#'   - Reads multiple CSV files containing enrichment data from a specified directory. The filenames (without
#'     extensions) serve as labels for individual comparisons.
#'
#'   - Combines the data into a single data frame with an added "Comparison" column indicating the source experiment.
#'
#'   - Selects the top enriched terms for each comparison based on the highest absolute Normalized Enrichment Score (NES)
#'     and compiles a master list of these terms.
#'
#'   - Filters the combined dataset to include only these top terms, ensuring consistency across all comparisons.
#'
#'   - Reorders the comparisons in the final visualizations based on the maximum absolute NES observed.
#'
#'   - Generates a significance label ("✱") for each gene set if the adjusted p-value (p.adjust) is below 0.05.
#'
#' @details
#' The script creates two primary visualizations:
#'
#'   1. A heatmap that displays differential enrichment across comparisons.
#'   2. A dot plot that highlights both the NES (color) and the significance (-log10(p.adjust)) of each term.
#'
#' Additional outputs include:
#'   - Core gene lists per term across comparisons.
#'   - A binary matrix of core gene presence/absence for computing Jaccard similarity between comparisons.
#'   - Expanded core enrichment heatmaps for individual terms.
#'
#' @section File Inputs:
#'   - CSV files from the directory: "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment"
#'     (matching the pattern "*.csv").
#'
#' @section Outputs:
#'   - Heatmaps and dot plots detailing enrichment profiles.
#'   - Core gene tables saved in CSV format.
#'   - SVG and PNG files for individual and overall enrichments.
#'
#' @note 
#'   Ensure that all file paths and package dependencies are correctly installed and set before running this script.
#'
#' @author
#'   Tobias Pohl

# -----------------------------------------------------
# Load Libraries
# -----------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, stringr, ggpubr, ggthemes, dplyr, tidyr, purrr, readr, pheatmap, tibble, tidyverse)

# -----------------------------------------------------
# Define Paths and Project Directories
# -----------------------------------------------------

# Set working directory
setwd("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler")
file_paths <- list.files(path = "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment", pattern = "*.csv", full.names = TRUE)
names(file_paths) <- basename(file_paths) %>% str_remove(".csv")

output_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results"
dir.create(output_dir, showWarnings = FALSE)

core_enrichment_dir <- file.path(output_dir, "core_enrichment")
dir.create(core_enrichment_dir, showWarnings = FALSE, recursive = TRUE)

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

# Prepare significance labels for the enrichment heatmap
lookup_df <- lookup_df %>%
  mutate(sig_label = ifelse(p.adjust < 0.05, "✱", ""))

# Reshape the data to create a matrix of NES values (rows: Description, columns: Comparison)
heatmap_data <- lookup_df %>%
  select(Description, Comparison, NES) %>%
  pivot_wider(names_from = Comparison, values_from = NES) %>%
  column_to_rownames("Description") %>%
  as.matrix()

# Similarly, create a corresponding matrix for significance labels
heatmap_labels <- lookup_df %>%
  select(Description, Comparison, sig_label) %>%
  pivot_wider(names_from = Comparison, values_from = sig_label) %>%
  column_to_rownames("Description") %>%
  as.matrix()

# Dynamically adjust the plot height to accommodate the number of rows and prevent label cutoff
plot_height <- max(8, nrow(heatmap_data) * 0.3)

# Generate heatmap
heatmap_plot <- pheatmap(
  heatmap_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = heatmap_labels,
  number_color = "black",
  color = colorRampPalette(c("#3b4cc0", "#ffffff", "#c03b3b"))(100),
  main = "Differential Enrichment Heatmap",
  fontsize = 12,
  fontsize_number = 10,
  border_color = NA,
  cellwidth = 15,
  cellheight = 15,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 90,
  silent = TRUE
)

# Render the heatmap ensuring that no elements are cut off
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)

# Save the heatmap to an SVG file with the calculated height
ggsave(output_file, heatmap_plot, width = 10, height = plot_height, dpi = 300)

# -----------------------------------------------------
# Create Dot Plot
# -----------------------------------------------------

ggplot(lookup_df, aes(x = Comparison, y = reorder(Description, NES), color = NES, size = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient2(low = "#3b4cc0", mid = "white", high = "#c03b3b", midpoint = 0) +
  theme_pubclean()

# -----------------------------------------------------
# Extract Core Genes per Comparison
# -----------------------------------------------------

# Create a named list of data frames: core genes per Comparison
core_gene_sets <- lapply(names(enrichment_list), function(name) {
  df <- enrichment_list[[name]]
  # Split semicolon-separated genes per row and unnest
  df_long <- df %>%
    select(Description, core_enrichment) %>%
    mutate(core_enrichment = str_split(core_enrichment, "/")) %>%
    unnest(core_enrichment) %>%
    mutate(Comparison = name)
  return(df_long)
})
names(core_gene_sets) <- names(enrichment_list)

# Combine all into one long dataframe
core_genes_df <- bind_rows(core_gene_sets)

# Optional: write core gene list for each Comparison & Description combo
output_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_genes_per_term/"
dir.create(output_dir, showWarnings = FALSE)

# Save full core gene table across all comparisons
write.csv(core_genes_df,
          file = file.path(output_dir, "core_genes_all_comparisons.csv"),
          row.names = FALSE)

core_genes_df %>%
  group_by(Comparison, Description) %>%
  summarise(Genes = list(unique(core_enrichment)), .groups = "drop") %>%
  rowwise() %>%
  mutate(file_name = paste0(output_dir, Comparison, "_", make.names(Description), ".csv")) %>%
  pwalk(~ write.csv(data.frame(Gene = ..3), file = ..4, row.names = FALSE))

# -----------------------------------------------------
# Compare Shared Core Genes Between Comparisons
# -----------------------------------------------------

# Aggregate all genes per Comparison (ignore term specificity for now)
core_gene_sets_aggregated <- core_genes_df %>%
  group_by(Comparison) %>%
  summarise(Genes = list(unique(core_enrichment)))

# Create a binary presence/absence matrix
all_genes <- unique(unlist(core_gene_sets_aggregated$Genes))
binary_matrix <- sapply(core_gene_sets_aggregated$Genes, function(gene_list) all_genes %in% gene_list)
rownames(binary_matrix) <- all_genes
colnames(binary_matrix) <- core_gene_sets_aggregated$Comparison

# Compute a Jaccard similarity matrix between conditions
jaccard_similarity <- function(x, y) {
  intersect = sum(x & y)
  union = sum(x | y)
  if (union == 0) return(NA)
  return(intersect / union)
}
jaccard_matrix <- outer(1:ncol(binary_matrix), 1:ncol(binary_matrix), Vectorize(function(i, j) {
  jaccard_similarity(binary_matrix[, i], binary_matrix[, j])
}))
rownames(jaccard_matrix) <- colnames(binary_matrix)
colnames(jaccard_matrix) <- colnames(binary_matrix)

# Plot heatmap of shared gene similarity
pheatmap(jaccard_matrix, main = "Jaccard Similarity of Core Genes per Comparison",
         color = colorRampPalette(c("white", "blue"))(100),
         display_numbers = TRUE, cluster_rows = TRUE, cluster_cols = TRUE)

# -----------------------------------------------------
# Expand Core Enrichment Genes for Heatmap
# -----------------------------------------------------

# Expand each enrichment file to get one gene per row
core_long_df <- bind_rows(
  lapply(names(enrichment_list), function(name) {
    df <- enrichment_list[[name]]
    df %>%
      select(Description, NES, core_enrichment) %>%
      mutate(core_enrichment = str_split(core_enrichment, "/")) %>%
      unnest(core_enrichment) %>%
      mutate(Comparison = name)
  })
)

# Rename for clarity
core_long_df <- core_long_df %>%
  rename(Gene = core_enrichment)

# First, aggregate duplicate Gene-Comparison pairs using mean NES
heatmap_df <- core_long_df %>%
  group_by(Gene, Comparison) %>%
  summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Comparison, values_from = NES)

# Then create the matrix
heatmap_matrix <- heatmap_df %>%
  column_to_rownames("Gene") %>%
  as.matrix()

heatmap_matrix[is.na(heatmap_matrix)] <- 0 # Replace NAs with 0 for heatmap visualization

# Now plot the heatmap
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Core Enrichment Gene Heatmap",
         fontsize_row = 8,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# Save the overall heatmap plot to a PNG file
output_file <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment", "core_enrichment_heatmap.png")
ggsave(output_file, width = 10, height = 8, dpi = 300)

# save to results dir

ggsave(file.path(core_enrichment_dir, "core_enrichment_heatmap.png"), width = 10, height = 8, dpi = 300)

# -----------------------------------------------------
# Generate Individual Core Enrichment Heatmaps
# -----------------------------------------------------

# Read log2fc data for each comparison from the mapped directory
log2fc_files <- list.files(
  path = "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/mapped",
  pattern = "*.csv",
  full.names = TRUE
)

names(log2fc_files) <- basename(log2fc_files) %>% str_remove(".csv")
log2fc_list <- lapply(names(log2fc_files), function(comp) {
  df <- read.csv(log2fc_files[comp])
  df$Comparison <- comp
  df
})
log2fc_df <- bind_rows(log2fc_list)

# Create output directory for individual core enrichment heatmaps
output_dir <- file.path(getwd(), "Results", "core_enrichment_heatmaps")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Loop over each enriched term and generate a corresponding heatmap
core_long_df %>%
  group_by(Description) %>%
  group_split() %>%
  walk(function(df_term) {
    description <- unique(df_term$Description)

    # Merge log2fc values based on matching Gene and Comparison
    df_term_log2fc <- df_term %>%
      left_join(log2fc_df, by = c("Gene" = "gene_symbol", "Comparison"))

    # Create a matrix of log2fc values (genes x comparisons)
    matrix_df <- df_term_log2fc %>%
      select(Gene, Comparison, log2fc) %>%
      group_by(Gene, Comparison) %>%
      summarize(log2fc = mean(log2fc, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Comparison, values_from = log2fc) %>%
      column_to_rownames("Gene")

    heatmap_matrix <- as.matrix(matrix_df)
    heatmap_matrix[is.na(heatmap_matrix)] <- 0

    # Define output filename with a safe name
    filename <- file.path(output_dir, paste0(make.names(description), ".svg"))

    # Set clustering options based on matrix dimensions
    cluster_rows_option <- nrow(heatmap_matrix) > 1
    cluster_cols_option <- ncol(heatmap_matrix) > 1

    # Set a symmetric color scale with white representing 0
    max_val <- max(abs(heatmap_matrix), na.rm = TRUE)
    breaks <- seq(-max_val, max_val, length.out = 101)

    # Dynamically calculate plot height to avoid label overlap
    row_height <- 0.05
    plot_height <- max(8, nrow(heatmap_matrix) * row_height + 0.2 * nrow(heatmap_matrix))

    # Save the heatmap as an SVG file
    svg(filename, width = 7, height = plot_height)
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(c("#3b4cc0", "#ffffff", "#c03b3b"))(100),
      breaks = breaks,
      main = description,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 8,
      cluster_rows = cluster_rows_option,
      cluster_cols = cluster_cols_option,
      border_color = NA,
      angle_col = 90,
      cellwidth = 15,
      cellheight = 15,
      legend_breaks = c(-max_val, 0, max_val),
      legend_labels = c(paste0("-", round(max_val, 2)), "0", round(max_val, 2))
    )
    dev.off()
  })