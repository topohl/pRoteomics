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
pacman::p_load(ggplot2, stringr, ggpubr, ggthemes, dplyr, tidyr, purrr, readr, pheatmap, tibble, tidyverse, RColorBrewer, writexl)

# -----------------------------------------------------
# Define Paths and Project Directories
# -----------------------------------------------------

# define ensemble profiling
# eg. baseline_cell_type_profiling,
# effects_chemogenetic_inhibition,
# effects_inhibition_memory_ensemble,
# learning_signature
# "interaction_with_learning"

# set ensemble prfiling
ensemble_profiling <- "baseline_cell_type_profiling"

# either CNO, VEH, CS or US or effects_inhibition_memory_ensemble or learning_signature
condition <- "US"

# Set working directory
setwd("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler")
# List CSV files from the core_enrichment subfolder specified by ensemble_profiling
file_paths <- list.files(
  path = file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)
names(file_paths) <- basename(file_paths) %>% str_remove(".csv")

output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/compareGO", ensemble_profiling, condition)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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

significant_only <- TRUE
top10_terms <- TRUE

if (top10_terms) {
  if (significant_only) {
    # Filter significant terms and then select top 10 by absolute NES per comparison
    top_df <- combined_df %>% 
      filter(p.adjust < 0.05) %>% 
      group_by(Comparison) %>% 
      slice_max(order_by = abs(NES), n = 10) %>% 
      ungroup()
  } else {
    # Select top 10 terms by absolute NES per comparison (without filtering for significance)
    top_df <- combined_df %>% 
      group_by(Comparison) %>% 
      slice_max(order_by = abs(NES), n = 10) %>% 
      ungroup()
  }
} else {
  if (significant_only) {
    # Keep only significant terms (p.adjust < 0.05)
    top_df <- combined_df %>% filter(p.adjust < 0.05)
  } else {
    # Use all terms if no filtering is applied
    top_df <- combined_df
  }
}

# Create a master list of top terms
top_terms <- unique(top_df$Description)

# -----------------------------------------------------
# Filter Data for Heatmap
# -----------------------------------------------------

combined_df <- combined_df %>%
  mutate(Comparison = str_trim(Comparison))

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
  mutate(sig_label = ifelse(p.adjust < 0.05, "*", ""))

# Select relevant columns and reshape the data
heatmap_data <- lookup_df %>%
  dplyr::select(Description, Comparison, NES) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = NES)

# Convert to data.frame before using column_to_rownames()
heatmap_data <- as.data.frame(heatmap_data)

# Set rownames based on the 'Description' column
rownames(heatmap_data) <- heatmap_data$Description

# Drop the 'Description' column, now that it's set as rownames
heatmap_data <- heatmap_data[, -which(names(heatmap_data) == "Description")]

# Convert to matrix
heatmap_data <- as.matrix(heatmap_data)

# save the heatmap data to a CSV file
file_name <- paste0("heatmap_data_", ensemble_profiling, "_", condition, ".xlsx")
heatmap_data_export <- tibble::rownames_to_column(as.data.frame(heatmap_data), var = "RowNames")
# Replace empty fields with NA
heatmap_data_export[heatmap_data_export == ""] <- NA
writexl::write_xlsx(heatmap_data_export, path = file.path(core_enrichment_dir, file_name))

# Similarly, create a corresponding matrix for significance labels
heatmap_labels <- lookup_df %>%
  dplyr::select(Description, Comparison, sig_label) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = sig_label) %>%
  tibble::column_to_rownames("Description")

as.matrix(heatmap_labels)

# Dynamically adjust the plot height to accommodate the number of rows and prevent label cutoff
plot_height <- max(8, nrow(heatmap_data) * 0.4)

# Use a clean, diverging palette (reversed RdBu from RColorBrewer)
my_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# Generate the heatmap with the ensemble_profiling title as headline
heatmap_plot <- pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = heatmap_labels,
  number_color = "black",
  color = my_colors,
  main = paste("Ensemble Profiling:", ensemble_profiling, "\nCondition:", condition),
  fontsize = 14,
  fontsize_number = 10,
  border_color = NA,
  cellwidth = 20,
  cellheight = 20,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 45,
  silent = TRUE
)

# Render the heatmap in VSCode/RStudio viewer
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)

# Define output directory and filename
output_file <- file.path(output_dir, paste0("enrichment_heatmap_", ensemble_profiling, "_", condition, ".svg"))

# Save the heatmap to an SVG file with clean font and correct height
svg(output_file, width = 10, height = plot_height, family = "Arial")
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# -----------------------------------------------------
# Create Dot Plot
# -----------------------------------------------------

dotplot <- ggplot(lookup_df, aes(
  x = Comparison,
  y = reorder(Description, NES, FUN = median),
  color = NES,
  size = -log10(p.adjust)
)) +
  geom_point(alpha = 0.85) +
  scale_color_gradientn(
    colours = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
    name = "NES",
    limits = c(min(lookup_df$NES, na.rm = TRUE), max(lookup_df$NES, na.rm = TRUE)),
    values = scales::rescale(c(min(lookup_df$NES, na.rm = TRUE), 0, max(lookup_df$NES, na.rm = TRUE)))
  ) +
  scale_size_continuous(
    name = expression(-log[10](p.adjust)),
    range = c(3, 10)
  ) +
  labs(
    title = "Comparative Gene Ontology Enrichment Dot Plot",
    subtitle = paste(ensemble_profiling, "under", condition, "condition"),
    x = paste("Comparison (", ensemble_profiling, ")", sep = ""),
    y = "Gene Set Description"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "black", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

# Render the dotplot
print(dotplot)

# Adjust the width dynamically based on the number of comparisons (with a minimum width of 10)
num_comparisons <- length(unique(lookup_df$Comparison))
dynamic_width <- max(8, num_comparisons * 1.2)

# Save the dotplot to an SVG file using dynamic width and the same dynamic plot_height
output_dotplot <- file.path(output_dir, paste0("enrichment_dotplot_", ensemble_profiling, "_", condition, ".svg"))
ggsave(output_dotplot, plot = dotplot, width = dynamic_width, height = plot_height, dpi = 300)

# -----------------------------------------------------
# Define genes to look for enrichment
# -----------------------------------------------------
# Define a list of genes to look for enrichment

# gene_list <- c("Tac2", "Tac3", "Tacr3")

gene_list <- c("P55099", "Q6NXX1")

# Unnest genes in core_enrichment
long_df <- combined_df %>%
  mutate(core_gene = strsplit(as.character(core_enrichment), "/|;|,|\\s+")) %>%
  unnest(core_gene)

# Filter for genes of interest
filtered_df <- long_df %>%
  filter(core_gene %in% gene_list) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

windows()

# Plot NES by Comparison, now colored and shaped by gene identity
ggplot(filtered_df, aes(
  x = Comparison,
  y = Description,
  color = NES,
  size = -log10(p.adjust),
  shape = core_gene
)) +
  geom_point(alpha = 0.9) +
  scale_color_gradientn(
    colours = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
    name = "NES"
  ) +
  scale_size_continuous(
    name = expression(-log[10](p.adjust)),
    range = c(3, 10)
  ) +
  labs(
    title = paste("Enrichment for Selected Genes:", paste(gene_list, collapse = ", ")),
    subtitle = paste("Ensemble profiling:", ensemble_profiling, "| Condition:", condition),
    x = "Comparison",
    y = "Enriched Gene Set",
    shape = "Gene"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# -----------------------------------------------------
# Extract Core Genes per Comparison
# -----------------------------------------------------

# Create a named list of data frames: core genes per Comparison
core_gene_sets <- lapply(names(enrichment_list), function(name) {
  df <- enrichment_list[[name]]
  # Split semicolon-separated genes per row and unnest
  df_long <- df %>%
    dplyr::select(Description, core_enrichment) %>%
    mutate(core_enrichment = str_split(core_enrichment, "/")) %>%
    unnest(core_enrichment) %>%
    mutate(Comparison = name)
  return(df_long)
})
names(core_gene_sets) <- names(enrichment_list)

# Combine all into one long dataframe
core_genes_df <- bind_rows(core_gene_sets)

# Optional: write core gene list for each Comparison & Description combo
output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_genes_per_term", ensemble_profiling, condition)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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

windows()

# Plot heatmap of shared gene similarity
jaccard_heatmap <- pheatmap(
  jaccard_matrix,
  main = "Jaccard Similarity of Core Genes per Comparison",
  color = colorRampPalette(c("white", "blue"))(100),
  display_numbers = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

# Render the plot in the VS Code R plot viewer
grid::grid.newpage()
grid::grid.draw(jaccard_heatmap$gtable)

# -----------------------------------------------------
# Expand Core Enrichment Genes for Heatmap
# -----------------------------------------------------

# Expand each enrichment file to get one gene per row
core_long_df <- bind_rows(
  lapply(names(enrichment_list), function(name) {
    df <- enrichment_list[[name]]
    df %>%
      dplyr::select(Description, NES, core_enrichment) %>%
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

# Read log2fc data for each comparison from the mapped/ensemble_profiling directory
log2fc_files <- list.files(
  path = file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

# Set names for the log2fc files with trimmed whitespace (directly check and trim here)
names(log2fc_files) <- basename(log2fc_files) %>% 
  str_remove(".csv") %>% 
  str_trim()  # Ensure no spaces before using them

# Print the names after trimming to check if it resolved the space issue
#print(names(log2fc_files))

# Proceed with the rest of your code and check the output again
log2fc_list <- lapply(names(log2fc_files), function(comp) {
  df <- read.csv(log2fc_files[comp])
  df$Comparison <- comp
  df
})

log2fc_df <- bind_rows(log2fc_list)

# print head of mcherry2_mcherry4 case in the Comparison column
# head(log2fc_df[log2fc_df$Comparison == "mcherry2_mcherry4", ])

# Create output directory for individual core enrichment heatmaps
output_dir <- file.path(getwd(), "Results", "core_enrichment_heatmaps", ensemble_profiling, condition)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Loop over each enriched term and generate a corresponding heatmap
core_long_df %>%
  group_by(Description) %>%
  group_split() %>%
  walk(function(df_term) {
    description <- unique(df_term$Description)

    # print description
    # print(description)

    # Merge log2fc values based on matching Gene and Comparison
    df_term_log2fc <- df_term %>%
      left_join(log2fc_df, by = c("Gene" = "gene_symbol", "Comparison"))
    
    write.csv(df_term_log2fc, file = "debug_neuron1_neuron2.csv", row.names = FALSE)

    #print(summary(core_long_df$Gene))
    #print(summary(log2fc_df$gene_symbol))
    #print(summary(core_long_df$Comparison))
    #print(summary(log2fc_df$Comparison))
    # check df_term_log2fc
    #print(head(df_term_log2fc))

    # Create a matrix of log2fc values (genes x comparisons)
    matrix_df <- df_term_log2fc %>%
      dplyr::select(Gene, Comparison, log2fc) %>%
      group_by(Gene, Comparison) %>%
      summarize(log2fc = mean(log2fc, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Comparison, values_from = log2fc) %>%
      column_to_rownames("Gene")

    # Clean up column names by trimming any leading/trailing spaces
    colnames(matrix_df) <- str_trim(colnames(matrix_df))

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

    # check matrix_df
    #print(head(matrix_df))

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