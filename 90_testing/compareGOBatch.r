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
pacman::p_load(ggplot2, stringr, ggpubr, ggthemes, dplyr, tidyr, purrr,
               readr, pheatmap, tibble, tidyverse, RColorBrewer, writexl)

uniprot_mapping_file_path <- file.path(
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler",
  "Datasets",
  "MOUSE_10090_idmapping.dat"
)

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
ensemble_profiling <- "interaction_with_learning"

# either CNO, VEH, CS or US or effects_inhibition_memory_ensemble or learning_signature
condition <- "effects_inhibition_memory_ensemble"

ont <- "BP"

# Set working directory
setwd("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler")
# List CSV files from the core_enrichment subfolder specified by ensemble_profiling
file_paths <- list.files(
  path = file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment", ont, ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)
names(file_paths) <- basename(file_paths) %>% str_remove(".csv")

output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/compareGO", ont, ensemble_profiling, condition)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

core_enrichment_dir <- file.path(output_dir, "core_enrichment")
dir.create(core_enrichment_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------
# Read and Combine Data
# -----------------------------------------------------
#' Read and Combine Enrichment Data
#'
#' @description
#' This script reads CSV files containing enrichment data from specified file paths, 
#' and combines them into a single data frame. An additional column, "Comparison", is added 
#' to indicate the source for each data frame.
#'
#' @param file_paths A named character vector, where the names represent comparison labels and the values are the respective CSV file paths.
#'
#' @return A data frame (combined_df) that aggregates all the individual data frames, 
#' with an extra column "Comparison" that holds the name corresponding to each CSV file.

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
#' Top Terms Selection for GO Batch Comparison
#' @description Filters and selects top terms from gene ontology comparisons based on significance and normalized enrichment scores.
#' @details 
#' The code uses two flags, 'top10_terms' and 'significant_only', to determine:
#'   - Whether to select only the top 10 terms per comparison or use all terms.
#'   - Whether to filter terms for significance (p.adjust < 0.05) prior to selection.
#' After filtering and ranking, top terms are aggregated into a master unique list.
#' @return A data frame ('top_df') with the selected terms per comparison and a vector ('top_terms') 
#'         containing unique term descriptions.

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
#' Process and Reorder Comparison Data for Heatmap
#' @description
#' This block filters the original data for heatmap visualization by:
#'   - Cleaning up the Comparison variable.
#'   - Subsetting the data to include rows with significance based on predefined top terms.
#'   - Computing an ordering of comparisons based on the maximum absolute NES values.
#'   - Reordering the Comparison factor to reflect the computed order.
#'
#' @return A modified data frame (lookup_df) with ordered factor levels for comparisons, 
#'         ready for visualization in a heatmap.

combined_df <- combined_df %>%
  mutate(Comparison = str_trim(Comparison))

# Include all rows for selected top terms
lookup_df <- combined_df %>%
  filter(Description %in% top_terms)

# Order comparisons by max |NES|
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
#' Generate Enrichment Heatmap
#'
#' This script processes enrichment analysis data to generate a heatmap.
#'
#' The workflow includes:
#' - Preparing significance labels based on adjusted p-values.
#' - Reshaping the data to create matrices for normalized enrichment scores (NES) and significance annotations.
#' - Saving the processed heatmap data to an Excel file.
#' - Dynamically adjusting the heatmap plot height to prevent label cutoff.
#' - Generating and rendering the heatmap with a clean, diverging RdBu color palette.
#' - Saving the final heatmap as an SVG file.
#'
#' @details
#' The input data (lookup_df) is expected to include the columns 'Description', 'Comparison', 'NES', and 'p.adjust'.
#' External packages used include dplyr, tidyr, tibble, writexl, pheatmap, RColorBrewer, and grid.
#'
#' @note Ensure that the required directories (core_enrichment_dir, output_dir) and variables (ensemble_profiling,
#'   condition, ont) are defined in the environment.

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

# Save the heatmap data to a CSV file
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
  main = paste(ont, "Ensemble Profiling:", ensemble_profiling, "\nCondition:", condition),
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
output_file <- file.path(output_dir, paste0("enrichment_heatmap_", ont, ensemble_profiling, "_", condition, ".svg"))

# Save the heatmap to an SVG file with clean font and correct height
svg(output_file, width = 10, height = plot_height, family = "Arial")
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# -----------------------------------------------------
# Create Dot Plot
# -----------------------------------------------------
#' Comparative Gene Ontology Enrichment Dot Plot
#'
#' Generates and saves a dot plot visualizing gene ontology enrichment comparisons.
#'
#' @details
#' This code creates a dot plot using ggplot2 where gene set descriptions are ordered by the median normalized enrichment score (NES).
#' Dots are sized by -log10(p.adjust) and colored according to the NES, using a gradient that spans negative to positive values.
#' The plot width is dynamically adjusted based on the number of unique comparisons.
#'
#' @param lookup_df A data frame containing GO enrichment results with columns including NES, p.adjust, Comparison, and Description.
#' @param ont A string representing the ontology category (e.g., "BP", "CC", or "MF").
#' @param ensemble_profiling A string specifying the ensemble profiling method.
#' @param condition A string indicating the experimental condition.
#' @param output_dir A string containing the file path to the directory where the plot will be saved.
#' @param plot_height Numeric value specifying the height of the plot when saved (in inches).
#'
#' @return The function renders a dot plot and saves it as an SVG file.

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
    subtitle = paste(ont, ensemble_profiling, "under", condition, "condition"),
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
dynamic_width <- max(10, num_comparisons * 4.5)

# Save the dotplot to an SVG file using dynamic width and the same dynamic plot_height
output_dotplot <- file.path(output_dir, paste0("enrichment_dotplot_", ont, "_", ensemble_profiling, "_", condition, ".svg"))
ggsave(output_dotplot, plot = dotplot, width = dynamic_width, height = plot_height, dpi = 300)

# -----------------------------------------------------
# Specify the target genes for the enrichment analysis
# -----------------------------------------------------
# Define a list of genes to look for enrichment
#'
#' @description Generates a scatter plot displaying the Normalized Enrichment Scores (NES)
#' for selected genes across different comparisons. The plot visualizes NES by color, 
#' -log10(p.adjust) by point size, and gene identity by point shape.
#'
#' @details
#' This script:
#'   - Splits the 'core_enrichment' field by various delimiters to unnest individual genes.
#'   - Filters the collapsed data for target genes specified in 'gene_list'.
#'   - Constructs a ggplot showing enrichment metrics and saves the resulting plot as an SVG file.
#'
#' @note
#' Requires pre-defined variables: 'combined_df', 'ont', 'ensemble_profiling', 'condition', 
#' and 'core_enrichment_dir'.
#'
#' @return
#' A ggplot object visualizing the enrichment of specified genes across comparisons,

gene_list <- c("P55099", "Q6NXX1")

# Unnest genes in core_enrichment
long_df <- combined_df %>%
  mutate(core_gene = strsplit(as.character(core_enrichment), "/|;|,|\\s+")) %>%
  unnest(core_gene)

# Filter for genes of interest
filtered_df <- long_df %>%
  filter(core_gene %in% gene_list) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Get min and max NES
nes_min <- min(filtered_df$NES, na.rm = TRUE)
nes_max <- max(filtered_df$NES, na.rm = TRUE)

# Use the larger absolute value for symmetric limits
max_abs_nes <- max(abs(nes_min), abs(nes_max))

windows()

# Plot NES by Comparison, now colored and shaped by gene identity
plot_selected_genes <- ggplot(filtered_df, aes(
  x = Comparison,
  y = Description,
  color = NES,
  size = -log10(p.adjust),
  shape = core_gene
)) +
  geom_point(alpha = 0.9) +
  scale_color_gradientn(
    colours = c("#6698CC", "white", "#F08C21"),
    limits = c(-max_abs_nes, max_abs_nes),
    values = c(0, 0.5, 1),  # Ensures NES = 0 is white
    name = "NES"
  ) +
  scale_size_continuous(
    name = expression(-log[10](p.adjust)),
    range = c(3, 10)
  ) +
  labs(
    title = paste("Enrichment for Selected Genes:", paste(gene_list, collapse = ", ")),
    subtitle = paste(ont, "Ensemble profiling:", ensemble_profiling, "| Condition:", condition),
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

# Save the plot to an SVG file
output_plot_file <- file.path(core_enrichment_dir, paste0("selected_genes_plot_", ont, "_", ensemble_profiling, "_", condition, ".svg"))
ggsave(output_plot_file, plot = plot_selected_genes, width = 10, height = 7, dpi = 300)

# Render the plot in the viewer
print(plot_selected_genes)

# -----------------------------------------------------
# Extract Core Genes per Comparison
# -----------------------------------------------------
#' Core Gene Extraction and CSV File Export
#' @description Extracts core genes for each comparison from an enrichment list, aggregates them,
#' and writes both a combined CSV file and individual CSV files for each unique combination
#' of Comparison and Description.
#'
#' @details The code processes a named list of enrichment data frames by:
#'   - Splitting the 'core_enrichment' column (which contains semicolon-separated gene names)
#'     into individual genes and reshaping the data to a long format.
#'   - Combining all processed data into a single data frame.
#'   - Saving the overall core gene table to a CSV file.
#'   - Grouping the data by Comparison and Description, aggregating unique core genes,
#'     and writing each group to a separate CSV file in a directory structure based on
#'     provided parameters (ont, ensemble_profiling, condition).
#'
#' @note The output directory is created if it does not exist.
#'
#' @param enrichment_list A named list of data frames; each data frame contains enrichment results
#'        with at least 'Description' and 'core_enrichment' columns.
#' @param base_dir Base directory defining where to store the output CSV files.
#' @param ont Ontology identifier used in the construction of the output directory path.
#' @param ensemble_profiling Profiling parameter used in the output directory path.
#' @param condition Condition parameter used in formulating the output directory path.
#'
#' @return No return value; CSV files are written to disk.

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

# Optional: Write core gene list for each Comparison & Description combo

# Define base directory for core genes output
base_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/core_genes"

# Construct output directory path using ont, ensemble_profiling, and condition
output_dir <- file.path(base_dir, ont, ensemble_profiling, condition)

# Create the output directory if it doesn't exist (including any necessary parent directories)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save the full core gene table across all comparisons to a CSV file
write.csv(core_genes_df,
      file = file.path(output_dir, "core_genes_all_comparisons.csv"),
      row.names = FALSE)

# For each unique combination of Comparison and Description,
# aggregate the unique core genes and write them into separate CSV files
core_genes_df %>%
  group_by(Comparison, Description) %>%                             # Group the data by Comparison and Description
  summarise(Genes = list(unique(core_enrichment)), .groups = "drop") %>% # Aggregate unique core genes into a list for each group
  rowwise() %>%                                                      # Process each row individually
  mutate(file_name = file.path(                                    # Create a file name for each group file
  output_dir,
  paste0(Comparison, "_", substr(make.names(Description), 1, 50), ".csv")  # Use Comparison and a cleaned, truncated Description
  )) %>%
  pwalk(~ write.csv(data.frame(Gene = ..3), file = ..4, row.names = FALSE)) # Write each group's gene list to its CSV file

# -----------------------------------------------------
# Compare Shared Core Genes Between Comparisons
# -----------------------------------------------------
#' Compare Shared Core Genes Between Comparisons
#'
#' Aggregates core gene sets per comparison, creates a binary presence/absence matrix,
#' computes the Jaccard similarity matrix, and visualizes the similarity as a heatmap.
#'
#' @details
#' This script performs the following steps:
#'   - Groups core genes by comparison and aggregates unique gene identifiers.
#'   - Constructs a binary matrix indicating the presence or absence of each gene across comparisons.
#'   - Calculates the Jaccard similarity index between all pairs of comparisons.
#'   - Visualizes the Jaccard similarity as a heatmap using the 'pheatmap' package.
#'
#' @note Customize the gene aggregation or similarity threshold as needed.

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
#' Generate and Save Core Enrichment Heatmap
#'
#' This section of code processes core enrichment gene data by expanding enrichment lists, 
#' aggregating duplicate gene-comparison pairs using mean NES, and pivoting the data into 
#' a wide-format matrix. A heatmap is then generated with hierarchical clustering and 
#' missing values replaced with zeros. Finally, the heatmap is saved in both SVG and PNG formats.
#'
#' @details The workflow includes:
#'   - Splitting core enrichment entries into one gene per row.
#'   - Averaging NES values for duplicate Gene-Comparison pairs.
#'   - Converting the aggregated data to a matrix for heatmap visualization.
#'   - Clustering rows and columns for pattern identification.
#'   - Saving the generated heatmap to specified file paths.
#'
#' @return A visual heatmap saved in specified file formats.

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

# Get min and max NES
nes_min <- min(heatmap_matrix, na.rm = TRUE)
nes_max <- max(heatmap_matrix, na.rm = TRUE)

# Use the larger absolute value for symmetric limits
max_abs_nes <- max(abs(nes_min), abs(nes_max))

# Create symmetric breaks for the color scale
breaks <- seq(-max_abs_nes, max_abs_nes, length.out = 101)

# Now plot the heatmap with adjusted color scale
heatmap_plot <- pheatmap(heatmap_matrix,
     color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
     breaks = breaks,
     main = "Core Enrichment Gene Heatmap",
     fontsize_row = 8,
     cluster_rows = TRUE,
     cluster_cols = TRUE)

# Save the overall heatmap plot to SVG file
output_file <- file.path(core_enrichment_dir, "core_enrichment_heatmap.svg")
svg(output_file, width = 10, height = 8)
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# Get top 25 and bottom 25 genes based on absolute NES values across comparisons
top_bottom_genes <- heatmap_df %>%
  pivot_longer(-Gene, names_to = "Comparison", values_to = "NES") %>%
  group_by(Gene) %>%
  summarize(max_abs_nes = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_nes)) %>%
  slice_head(n = 50) %>%
  bind_rows(
    arrange(., max_abs_nes) %>% slice_head(n = 50)
  ) %>%
  distinct(Gene)

# Filter the heatmap matrix to only include the top and bottom 25 genes
top_bottom_matrix <- heatmap_matrix[rownames(heatmap_matrix) %in% top_bottom_genes$Gene, , drop = FALSE]

# Order rows by max abs NES again for nice visual
gene_order <- heatmap_df %>%
  filter(Gene %in% rownames(top_bottom_matrix)) %>%
  pivot_longer(-Gene, names_to = "Comparison", values_to = "NES") %>%
  group_by(Gene) %>%
  summarize(order_metric = max(abs(NES), na.rm = TRUE)) %>%
  arrange(desc(order_metric)) %>%
  pull(Gene)

top_bottom_matrix <- top_bottom_matrix[gene_order, , drop = FALSE]

# Get min and max NES
nes_min <- min(top_bottom_matrix, na.rm = TRUE)
nes_max <- max(top_bottom_matrix, na.rm = TRUE)

# Use the larger absolute value for symmetric limits
max_abs_nes <- max(abs(nes_min), abs(nes_max))

# Create symmetric breaks for the color scale
breaks <- seq(-max_abs_nes, max_abs_nes, length.out = 101)

top_bottom_plot <- pheatmap(
  top_bottom_matrix,
  color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
  breaks = breaks,
  main = "Top & Bottom 50 Core Enrichment Genes",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 8,
  fontsize_number = 6,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  border_color = "#e7e7e7",
  border_width = 0.1,
  cellwidth = 12,
  cellheight = 8,
  angle_col = 45,
  treeheight_row = 20,
  treeheight_col = 20,
  legend = TRUE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2", "-1", "0", "1", "2"),
  annotation_legend = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE
)

# Save to SVG
svg(file.path(core_enrichment_dir, "top_bottom_core_enrichment_heatmap.svg"), width = 10, height = 10)
grid::grid.draw(top_bottom_plot$gtable)
dev.off()

# -----------------------------------------------------
# Generate Individual Core Enrichment Heatmaps
# -----------------------------------------------------
#' Generate and Save Core Enrichment Heatmaps
#'
#' @description
#' Reads log2fc CSV files for different comparisons, consolidates them, merges the data
#' with core enrichment results, and generates heatmaps for each enriched term. The heatmaps 
#' are exported as SVG files.
#'
#' @details
#' The script performs the following steps:
#' \itemize{
#'   \item Reads log2fc data from CSV files residing in a specified folder.
#'   \item Trims file names to avoid issues with whitespace.
#'   \item Consolidates individual log2fc datasets into a single dataframe.
#'   \item Merges the log2fc values with a core enrichment dataset based on gene and comparison.
#'   \item Constructs a matrix for each enriched term and generates a heatmap using the 'pheatmap' package.
#'   \item Dynamically calculates plot dimensions and saves the output as SVG.
#' }
#'
#' @param ensemble_profiling A character string indicating the ensemble profiling folder segment.
#' @param condition A character string specifying the experimental condition used in file paths.
#' @param core_long_df A dataframe containing core enrichment information (long format).
#'
#' @return No value is returned. The function creates SVG files with heatmaps as a side effect.

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
output_dir <- file.path(getwd(), "Results", "core_enrichment_heatmaps", ont, ensemble_profiling, condition)
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

    write.csv(df_term_log2fc, file = "debug_neuron1_neuron2.csv", row.names = FALSE)

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
    filename <- file.path(output_dir, paste0(substr(make.names(description), 1, 20), ".svg"))

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