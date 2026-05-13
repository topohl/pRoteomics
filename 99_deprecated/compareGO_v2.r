#' compareGO.r - Comparative Gene Ontology Enrichment Analysis and Visualization
#'
#' @description
#' This script performs a comparative analysis of GO enrichment across multiple experiments.
#' It:
#'   - Reads multiple CSV files containing enrichment data from a specified directory. 
#'     The filenames (without extensions) serve as labels for individual comparisons.
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
#'   2. A dot plot that highlights both the NES (color) and
#'      the significance (-log10(p.adjust)) of each term.
#'
#' Additional outputs include:
#'   - Core gene lists per term across comparisons.
#'   - A binary matrix of core gene presence/absence
#'     for computing Jaccard similarity between comparisons.
#'   - Expanded core enrichment heatmaps for individual terms.
#'
#' @section File Inputs:
#'   - CSV files from the directory:
#'     "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment"
#'     (matching the pattern "*.csv").
#'
#' @section Outputs:
#'   - Heatmaps and dot plots detailing enrichment profiles.
#'   - Core gene tables saved in CSV format.
#'   - SVG and PNG files for individual and overall enrichments.
#'
#' @note
#'   Ensure that all file paths and package dependencies are
#'   correctly installed and set before running this script.
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
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics",
  "Datasets",
  "MOUSE_10090_idmapping.dat"
)

# Read the UniProt mapping file and extract only Uniprot-to-GeneName mappings
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  filter(V2 == "Gene_Name") %>%
  select(UniprotID = V1, Gene_Name = V3)

# -----------------------------------------------------
# Define Paths and Project Directories
# -----------------------------------------------------

#' Configure Analysis Parameters and Directory Structure
#'
#' @description
#' This section sets up the analysis parameters and directory structure for
#' Gene Ontology (GO) enrichment analysis using clusterProfiler.
#'
#' @section Analysis Parameters:
#' - ensemble_profiling: Type of ensemble profiling analysis
#'   (e.g., "baseline_cell_type_profiling", "effects_chemogenetic_inhibition",
#'   "interaction_with_learning", etc.)
#' - condition: Experimental condition being analyzed (e.g., "CS", "US",
#'   "effects_inhibition_memory_ensemble", "learning_signature" etc.)
#' - ont: Gene Ontology domain (MF/BP/CC)
#'
#' @section Directory Structure:
#' - Input:  Datasets/core_enrichment/{ont}/{ensemble_profiling}/{condition}/
#' - Output: Results/compareGO/{ont}/{ensemble_profiling}/{condition}/
# Set analysis parameters ---------------------------------------------

# Define the ensemble profiling method used in the analysis
ensemble_profiling <- "phenotype_within_unit_simplified"  # e.g., "baseline_cell_type_profiling", "effects_chemogenetic_inhibition", "interaction_with_learning", "phenotype_within_unit_simplified"

# Specify the experimental condition (e.g., CNO, VEH, CS, US, effects_inhibition_memory_ensemble, or learning_signature)
condition <- "CA2"

# Define the Gene Ontology domain (e.g., MF, BP, or CC)
ont <- "BP"  # Biological Process

# Set up working environment ------------------------------------------

# Set the working directory to the project’s base folder
setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics")

# Import input files --------------------------------------------------

# List all CSV files from the specified core_enrichment subfolder
# The folder structure is based on the ontology, profiling method, and condition
file_paths <- list.files(
  path = file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/core_enrichment", ont, ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)
# Name each file entry using the base filename (without the .csv extension)
names(file_paths) <- basename(file_paths) %>% str_remove(".csv")

# Define output directories -------------------------------------------

# Set the output directory for comparative GO analysis results
output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/compareGO", ont, ensemble_profiling, condition)
# Create the output directory if it doesn't already exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define and create a subdirectory for core enrichment outputs
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

# Read in go_ids 
gene_go_ids <- combined_df %>%
  # Select just the columns needed
  select(GO_ID = ID, core_enrichment) %>%
  # Split the slash-separated gene accessions into rows
  mutate(core_enrichment = strsplit(core_enrichment, "/")) %>%
  unnest(core_enrichment) %>%
  rename(Gene = core_enrichment) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(GO_IDs = paste(unique(GO_ID), collapse = "; "), .groups = "drop")

# read in p-values
p_values <- combined_df %>%
  select(Comparison, Description, p.adjust) %>%
  distinct()

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

# Define filtering parameters for selecting top GO terms
significant_only <- TRUE  # If TRUE, only keep terms with p.adjust < 0.05
top10_terms <- TRUE       # If TRUE, select the top 10 terms based on absolute NES for each comparison

# Filter and select top terms from the combined data frame
if (top10_terms) {
  if (significant_only) {
    # Only include significant terms (p.adjust < 0.05), then select the top 10 with the highest absolute NES per Comparison
    top_df <- combined_df %>%
      filter(p.adjust < 0.05) %>%
      group_by(Comparison) %>%
      slice_max(order_by = abs(NES), n = 10) %>%
      ungroup()
  } else {
    # Without significance filtering, select the top 10 terms with the highest absolute NES per Comparison
    top_df <- combined_df %>%
      group_by(Comparison) %>%
      slice_max(order_by = abs(NES), n = 10) %>%
      ungroup()
  }
} else {
  if (significant_only) {
    # If not selecting top10, just retain all significant terms (p.adjust < 0.05)
    top_df <- combined_df %>%
      filter(p.adjust < 0.05)
  } else {
    # Include all terms without any filtering
    top_df <- combined_df
  }
}

# -----------------------------------------------------------------------------
# Create a master list of GO term descriptions from the filtered data
# -----------------------------------------------------------------------------

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

# Remove any extra whitespace from the Comparison names to ensure consistency
combined_df <- combined_df %>%
  mutate(Comparison = str_trim(Comparison))

# Filter the combined data to include only rows with GO term descriptions that are in the top_terms list
lookup_df <- combined_df %>%
  filter(Description %in% top_terms)

# Compute the order of comparisons by determining the maximum absolute NES for each, then sort in descending order
comparison_order <- lookup_df %>%
  group_by(Comparison) %>%
  summarize(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_NES)) %>%
  pull(Comparison)

# Recast the Comparison column as a factor using the determined order to ensure proper ordering in plots
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

# Replace NA with 0
heatmap_data[is.na(heatmap_data)] <- 0

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
  cluster_cols = FALSE,
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

# Construct dot plot for comparative GO enrichment analysis.
dotplot <- ggplot(lookup_df, aes(
  x = Comparison,
  y = reorder(Description, NES, FUN = median),  # Order gene sets by median NES
  color = NES,
  size = -log10(p.adjust)
)) +
  geom_point(alpha = 0.85) +  # Slight transparency for overlapping points
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

# Dynamically adjust the plot width based on the number of unique comparisons.
# A minimum width is set to ensure clarity when few comparisons are present.
num_comparisons <- length(unique(lookup_df$Comparison))
dynamic_width <- max(10, num_comparisons * 4.5)

# Save the dot plot as an SVG file using the dynamic width and predefined plot height.
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

# Define the target gene identifiers for enrichment analysis
gene_list <- c("P55099", "Q6NXX1")

# Split the 'core_enrichment' field into individual genes by various delimiters,
# then expand the data frame so that each gene is represented in its own row.
long_df <- combined_df %>%
  mutate(core_gene = strsplit(as.character(core_enrichment), "/|;|,|\\s+")) %>%
  unnest(core_gene)

# Retain only the rows corresponding to the genes of interest,
# and convert the 'Description' column into an ordered factor for consistent plotting.
filtered_df <- long_df %>%
  filter(core_gene %in% gene_list) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Compute the minimum and maximum Normalized Enrichment Scores (NES) across the filtered data
nes_min <- min(filtered_df$NES, na.rm = TRUE)
nes_max <- max(filtered_df$NES, na.rm = TRUE)

# Determine the symmetric limit for color scaling by taking the larger absolute NES value
max_abs_nes <- max(abs(nes_min), abs(nes_max))

# Open a new graphics window (useful for interactive sessions)
windows()

# Construct a ggplot object to visualize enrichment metrics:
# - X-axis: Comparisons
# - Y-axis: Enriched gene set descriptions
# - Point color represents the NES,
# - Point size is scaled by -log10(p.adjust)
# - Different point shapes denote distinct genes
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
    values = c(0, 0.5, 1),  # Ensures that a NES of 0 appears white for neutral contrast
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

# Define the output file path to save the generated plot in SVG format
output_plot_file <- file.path(core_enrichment_dir, paste0("selected_genes_plot_", ont, "_", ensemble_profiling, "_", condition, ".svg"))

# Export the plot to the specified SVG file with fixed dimensions and resolution
ggsave(output_plot_file, plot = plot_selected_genes, width = 10, height = 7, dpi = 300)

# Display the plot
print(plot_selected_genes)

# -----------------------------------------------
# Gene-centric Enrichment: Frequency & NES summary
# -----------------------------------------------

# Count occurrences of each gene in each Comparison and calculate mean NES
gene_enrichment_summary <- long_df %>%
  filter(core_gene %in% gene_list) %>%
  group_by(core_gene, Comparison) %>%
  summarise(
    count = n(),
    mean_NES = mean(NES, na.rm = TRUE),
    max_padj = min(p.adjust, na.rm = TRUE),  # Best (smallest) adjusted p-value across gene sets
    .groups = "drop"
  )

# Plot: Each point is a gene in a comparison
plot_gene_summary <- ggplot(gene_enrichment_summary, aes(
  x = Comparison,
  y = core_gene,
  size = count,
  fill = mean_NES
)) +
  geom_point(shape = 21, stroke = 0) +
  scale_size_continuous(
    name = "Occurrences\nin Gene Sets",
    range = c(3, 10)
  ) +
  scale_fill_gradient2(
    name = "Mean NES",
    low = "#6698CC", mid = "white", high = "#F08C21",
    midpoint = 0, limits = c(-max_abs_nes, max_abs_nes)
  ) +
  labs(
    title = "Gene-Centric Enrichment Across Comparisons",
    subtitle = paste("Ontology:", ont, "| Profiling:", ensemble_profiling, "| Condition:", condition),
    x = "Comparison",
    y = "Gene"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right"
  )

# Save plot
gene_summary_file <- file.path(core_enrichment_dir, paste0("gene_centric_enrichment_", ont, "_", ensemble_profiling, "_", condition, ".svg"))
ggsave(gene_summary_file, plot = plot_gene_summary, width = 5, height = 4, dpi = 300)

# Show plot
print(plot_gene_summary)

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
base_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/core_genes"

# Construct output directory path using ont, ensemble_profiling, and condition
output_dir <- file.path(base_dir, ont, ensemble_profiling, condition)

# Create the output directory if it doesn't exist (including any necessary parent directories)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save the full core gene table across all comparisons to a CSV file
write.csv(core_genes_df,
      file = file.path(output_dir, "core_genes_all_comparisons.csv"),
      row.names = FALSE)

# Group the core genes by Comparison and Description, and then export each group as a separate CSV file.
core_genes_df %>%
  group_by(Comparison, Description) %>%
  summarise(Genes = list(unique(core_enrichment)), .groups = "drop") %>%
  rowwise() %>%
  mutate(file_name = file.path(
    output_dir,
    paste0(Comparison, "_", substr(make.names(Description), 1, 50), ".csv")
  )) %>%
  # For each group (defined by Comparison and Description), write the unique core genes (in column 3)
  # to a CSV file. Here, pwalk iterates over each group’s elements, creating a one-column data frame
  # containing the unique gene IDs, and saves it to the file path provided in element 4.
  pwalk(~ write.csv(data.frame(Gene = ..3), file = ..4, row.names = FALSE))

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

# Save the Jaccard heatmap to an SVG file
jaccard_sim_matrix <- file.path(core_enrichment_dir, "jaccard_similarity_heatmap.svg")
svg(jaccard_sim_matrix, width = 10, height = 8)
grid::grid.draw(jaccard_heatmap$gtable)
dev.off()

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

# Reshape to long format
long_heatmap_df <- heatmap_df %>%
  pivot_longer(-Gene, names_to = "Comparison", values_to = "NES")

# Top 25 upregulated genes (highest positive NES)
top25_up <- long_heatmap_df %>%
  group_by(Gene) %>%
  summarize(max_nes = max(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_nes)) %>%
  slice_head(n = 25) %>%
  mutate(Direction = "Up")

# Top 25 downregulated genes (lowest negative NES)
top25_down <- long_heatmap_df %>%
  group_by(Gene) %>%
  summarize(min_nes = min(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_nes) %>%
  slice_head(n = 25) %>%
  mutate(Direction = "Down")

# Combine without deduplication
top_bottom_genes <- bind_rows(top25_up, top25_down)

# Create a mapping of each gene to its GO terms (Description)
gene_go <- core_long_df %>%
  select(Gene, Description) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(GO_Terms = paste(unique(Description), collapse = "; "), .groups = "drop")

# Add the GO term column to the gene list
top_bottom_genes <- top_bottom_genes %>%
  left_join(gene_go, by = "Gene")

# Save gene list with direction and GO term
write_xlsx(
  top_bottom_genes,
  file.path(core_enrichment_dir, "top_bottom_gene_list.xlsx")
)

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

# Load UniProt mapping (only Gene_Name)
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  filter(V2 == "Gene_Name") %>%
  select(UniprotID = V1, Gene_Name = V3)

# Update rownames in heatmap
mapped_names <- left_join(
  tibble(Gene = rownames(top_bottom_matrix)),
  uniprot_df,
  by = c("Gene" = "UniprotID")
)

rownames(top_bottom_matrix) <- ifelse(
  is.na(mapped_names$Gene_Name),
  mapped_names$Gene,
  mapped_names$Gene_Name
)

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
  main = "Core Enrichment Overview",
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
svg(file.path(core_enrichment_dir, "top_bottom_core_enrichment_heatmap.svg"), width = 5, height = 10)
grid::grid.draw(top_bottom_plot$gtable)
dev.off()

# -------------------------------------------------------
# Create Excel file with top/bottom genes and Uniprot IDs
# -------------------------------------------------------

#log2fc_files <- list.files(
#  path = file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/mapped", ensemble_profiling, condition),
#  pattern = "*.csv",
#  full.names = TRUE
#)

# Save top_bottom_matrix to excel with Uniprot gene names
# Read the Uniprot mapping file (adjust parameters as needed depending on file format)
#uniprot_df <- read.delim(uniprot_mapping_file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#
# Assume the first two columns are Uniprot ID and Gene Name; rename them accordingly
#colnames(uniprot_df)[1:2] <- c("UniprotID", "Gene_Name", "UniProtKB-ID")

# Convert top_bottom_matrix to a data frame with a Gene column
#top_bottom_df <- tibble::rownames_to_column(as.data.frame(top_bottom_matrix), var = "Gene")
#
# Merge gene names again for Excel (in case they weren't used as rownames above)
#top_bottom_df <- left_join(top_bottom_df, uniprot_df, by = c("Gene" = "UniprotID")) %>%
#  relocate(Gene_Name, .before = Gene) %>%
#  rename(Uniprot_ID = Gene)

#write_xlsx(
#  top_bottom_df,
#  file.path(core_enrichment_dir, "top_bottom_core_enrichment_genes.xlsx")
#)

# -------------------------------------------------------
# Create per-comparison heatmaps and Excel files
# -------------------------------------------------------

# Load UniProt mapping file (dat format)
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Extract relevant mappings
uniprot_subset <- uniprot_df %>%
  filter(V2 %in% c("Gene_Name", "UniProtKB-ID")) %>%
  pivot_wider(names_from = V2, values_from = V3, values_fn = list) %>%
  unnest(cols = everything()) %>%
  distinct(V1, .keep_all = TRUE) %>%
  rename(UniprotID = V1)

# Extract Gene Synonyms
gene_synonyms <- uniprot_df %>%
  filter(V2 == "Gene_Synonym") %>%
  distinct(V1, .keep_all = TRUE) %>%
  rename(UniprotID = V1, Gene_Synonym = V3)

# Optional: Gene Descriptions (from GSEA if available)
gene_descriptions <- core_long_df %>%
  select(Gene, Description) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(Description = paste(unique(Description), collapse = "; "), .groups = "drop")

# List log2fc input files
log2fc_files <- list.files(
  path = file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

# Read and combine all log2fc files
log2fc_df_list <- lapply(log2fc_files, function(f) {
  df <- read_csv(f, col_types = cols())
  df$Comparison <- tools::file_path_sans_ext(basename(f))
  return(df)
})

log2fc_long_df <- bind_rows(log2fc_df_list)

# Unique comparisons
unique_comparisons <- unique(log2fc_long_df$Comparison)

# Loop over each comparison
for (comp in unique_comparisons) {
  full_comp_df <- log2fc_long_df %>%
  filter(Comparison == comp, !is.na(log2fc), !is.na(padj))

# Volcano plot section (always runs, even if no sig genes)
volcano_df <- full_comp_df %>%
  filter(!str_detect(gene_symbol, "_")) %>%
  left_join(uniprot_subset %>% select(UniprotID, Gene_Name),
            by = c("gene_symbol" = "UniprotID")) %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2fc > 0 ~ "up",
    padj < 0.05 & log2fc < 0 ~ "down",
    TRUE ~ "n.s."
  ))
  
volcano_plot <- ggplot(volcano_df, aes(x = log2fc, y = -log10(padj), label = Gene_Name, color = Significance)) +
    geom_point(shape = 16, alpha = 0.6, size = 4) +
    ggrepel::geom_text_repel(
      data = subset(volcano_df, padj < 0.05),
      aes(label = Gene_Name),
      size = 4,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("up" = "#f36d07", "down" = "#455A64", "n.s." = "#c7c7c7"),
      breaks = c("up", "down", "n.s.")
    ) +
    scale_x_continuous(breaks = seq(floor(min(volcano_df$log2fc, na.rm = TRUE)),
                                    ceiling(max(volcano_df$log2fc, na.rm = TRUE)),
                                    by = 1)) +
    labs(
      title = paste(strsplit(comp, "_")[[1]], collapse = " over "),
      x = "log2 Fold Change",
      y = expression(-log[10]~(p~adj.))
    ) +
    coord_fixed(ratio = 2) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#455A64")
  
  # Save volcano plot
  # Create subfolder for volcano plots if it doesn't exist
  vol_dir <- file.path(core_enrichment_dir, "volcano_plots")
  if (!dir.exists(vol_dir)) {
    dir.create(vol_dir, recursive = TRUE)
  }
  
  ggsave(
    filename = file.path(vol_dir, paste0("log2fc_", comp, "_volcano.svg")),
    plot = volcano_plot,
    width = 5,
    height = 5
  )
  
  # Now check for significant data for heatmap + Excel
comp_df <- full_comp_df %>% filter(padj < 0.05)

if (nrow(comp_df) == 0) {
  message(paste("No significant genes found for", comp, "- skipping heatmap and Excel."))
  next
}

  # Select top and bottom 25 based on log2fc direction
  comp_df <- comp_df %>%
    mutate(Direction = case_when(
      log2fc > 0 ~ "Up",
      log2fc < 0 ~ "Down",
      TRUE ~ "Neutral"
    ))

  # Combine top 25 up and down (before deduplication)
  top_bottom_genes <- bind_rows(
    comp_df %>%
      filter(Direction == "Up") %>%
      arrange(desc(log2fc)) %>%
      slice_head(n = 25),
    comp_df %>%
      filter(Direction == "Down") %>%
      arrange(log2fc) %>%
      slice_head(n = 25)
  ) %>%
    distinct(gene_symbol, Comparison, .keep_all = TRUE)

  # Join UniProt annotations
  mapped_top_bottom <- top_bottom_genes %>%
    left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
    left_join(gene_synonyms, by = c("gene_symbol" = "UniprotID")) %>%
    left_join(gene_descriptions, by = c("gene_symbol" = "Gene")) %>%
    left_join(gene_go_ids, by = c("gene_symbol" = "Gene")) %>%
    mutate(
      Gene_Name = as.character(Gene_Name),
      `UniProtKB-ID` = as.character(`UniProtKB-ID`),
      Gene_Synonym = as.character(Gene_Synonym),
      Gene_Label = ifelse(is.na(Gene_Name), gene_symbol, Gene_Name),
      Gene_Label = make.unique(as.character(Gene_Label)),
      abs_log2fc = abs(log2fc),
      UniProt_Link = paste0("https://www.uniprot.org/uniprot/", gene_symbol)
    )

  # Add rank within direction
  mapped_top_bottom <- mapped_top_bottom %>%
    group_by(Direction) %>%
    arrange(Direction, desc(abs_log2fc)) %>%
    mutate(Rank = row_number()) %>%
    ungroup()

  # Create heatmap matrix
  heatmap_matrix <- mapped_top_bottom %>%
    select(Gene_Label, log2fc) %>%
    mutate(Gene_Label = make.unique(as.character(Gene_Label))) %>%
    tibble::column_to_rownames("Gene_Label") %>%
    as.matrix()

  if (nrow(heatmap_matrix) < 2) {
    message(paste("Skipping heatmap for", comp, "- fewer than 2 unique genes."))
  } else {
    max_abs <- max(abs(heatmap_matrix), na.rm = TRUE)
    breaks <- seq(-max_abs, max_abs, length.out = 101)

    heatmap_plot <- pheatmap(
      heatmap_matrix,
      color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
      breaks = breaks,
      main = paste("Top/Bottom 25 Genes by log2FC in", comp),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize_row = 6,
      show_colnames = FALSE,
      border_color = "#e7e7e7"
    )

    # Save heatmap
    heatmap_dir <- file.path(core_enrichment_dir, "heatmap_plots")
    if (!dir.exists(heatmap_dir)) {
      dir.create(heatmap_dir, recursive = TRUE)
    }
    svg(
      file.path(heatmap_dir, paste0("log2fc_", comp, "_heatmap.svg")),
      width = 3, height = 8
    )
    grid::grid.draw(heatmap_plot$gtable)
    dev.off()
  }

  # Save Excel with annotation
  excel_out <- mapped_top_bottom %>%
    select(
      Gene_Name,
      `UniProtKB-ID`,
      Uniprot_Accession = gene_symbol,
      Gene_Synonym,
      Comparison,
      log2fc,
      abs_log2fc,
      padj,
      Direction,
      Rank,
      GO_IDs,
      Description,
      UniProt_Link
    ) %>%
    mutate(across(everything(), as.character)) %>%
    distinct()

  proteins_sig_reg_dir <- file.path(core_enrichment_dir, "significant_proteins")
  dir.create(proteins_sig_reg_dir, showWarnings = FALSE, recursive = TRUE)
  write_xlsx(
    excel_out,
    file.path(proteins_sig_reg_dir, paste0("log2fc_", comp, "_genes.xlsx"))
  )
}


# -----------------------------------------------------
# Generate Individual Core Enrichment Heatmaps
# -----------------------------------------------------
#' Generate and Save Core Enrichment Heatmaps
#'
#' @description
#' This section reads log2fc data from CSV files and merges it with core enrichment results,
#' then generates heatmaps for each enriched term. Each heatmap is dynamically sized and saved as an SVG file.
#'
#' @details
#' The following steps are performed:
#'   - Log2fc data is read from CSV files located in the mapped/ensemble profiling directory.
#'   - File names are trimmed to remove any extraneous whitespace.
#'   - Individual log2fc datasets are consolidated into a single dataframe.
#'   - The consolidated log2fc data is merged with the core enrichment dataframe based on matching gene identifiers and comparisons.
#'   - For each enriched term, a matrix is created from the merged data.
#'   - A heatmap is generated using the 'pheatmap' package with dynamically determined dimensions to prevent label overlap.
#'   - Each heatmap is exported and saved as an SVG file.
#'
#' @param ensemble_profiling A character string indicating the ensemble profiling method used in forming file paths.
#' @param condition A character string specifying the experimental condition.
#' @param core_long_df A dataframe in long format containing core enrichment results.
#'
#' @return No return value. The function produces SVG files of heatmaps as a side effect.
#
# Read log2fc data for each comparison from the mapped/ensemble_profiling directory

log2fc_files <- list.files(
  path = file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/Datasets/mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

# Set names for the log2fc files with trimmed whitespace (directly check and trim here)
names(log2fc_files) <- basename(log2fc_files) %>%
  str_remove(".csv") %>%
  str_trim()  # Ensure no spaces before using them

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
      color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
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

# -----------------------------------------------------
# Gene-Centric Log2 Fold-Change (log2FC) Analysis
# -----------------------------------------------------
#'
#' @description
#' Processes multiple CSV files containing gene expression data (log2FC) to:
#'   - Add comparison information derived from file names.
#'   - Filter for a specific gene list.
#'   - Summarize data by gene and comparison.
#'   - Generate a bubble plot that visualizes the mean log2FC values.
#'
#' @details
#' The analysis workflow is as follows:
#'   1. Read CSV files from a designated directory.
#'   2. Extract a "Comparison" field from the filename (removing its extension) and append it to the data.
#'   3. Combine the individual data frames into one long-form data frame.
#'   4. Optionally filter the combined data to include only genes from a provided gene list.
#'   5. Group the data by gene and comparison to compute:
#'      - The count of occurrences.
#'      - The mean log2FC.
#'      - The minimum adjusted p-value (padj).
#'   6. Create a bubble plot with ggplot2 where:
#'      - The point size reflects the occurrence count.
#'      - The fill color represents the mean log2FC on a custom gradient scale.
#'   7. Save the generated plot as an SVG file and display it.
#'#'
#' @section Inputs:
#'   - CSV files containing gene expression data (with columns "gene_symbol", "log2fc", "padj") located 
#'     in a directory formed using the variables ensemble_profiling and condition.
#'   - gene_list: A vector of gene symbols for filtering the data; customize as necessary.
#'
#' @section Outputs:
#'   - A summarized data frame (gene_fc_summary) that includes, per gene and comparison, the count,
#'     mean log2FC, and the minimum adjusted p-value.
#'   - A ggplot2 bubble plot that visualizes gene-centric log2FC data.
#'   - An SVG file containing the bubble plot, saved to a path defined by core_enrichment_dir, ensemble_profiling, and condition.
#'
#' @usage
#' Source this script in an R environment after ensuring that:
#'   - All required variables (ensemble_profiling, condition, gene_list) are properly set.
#'   - The necessary dependencies are available in the environment.

# LTP gene list
#gene_list <- c(
#  "P35438", "P23818", "P11798", "P05132", "Q01147", "P21237",
#  "Q62108", "O88935", "P63085", "Q9JLN9", "P31750", "P23819",
#  "Q63844", "P15209", "P18760", "O88602"
#)

# gene list for plc
#gene_list <- c("A3KGF7", "P51432", "Q62077", "Q8CIH5", "Q8K2J0",
#  "Q8K3R3", "Q8K4S1", "Q8R3B1", "Q9Z1B3")

# gene list for tacr3 signalling pathway
#gene_list <- c("Q8K3R3", "Q8R071", "P11798", "P28652", "P20444", "P68404", "P05132", "P68181")
# TACR3-signalling broad
#gene_list <- sort(c(
#  "P21279", "P21278", "P51432", "P11881", "P68404",
#  "P63318", "P0DP26", "P0DP27", "P0DP28", "P11798",
#  "P28652", "Q61411", "Q99N57", "P31938", "P63085",
#  "Q63844", "Q01147"
#))

# TACR3-signalling CORE
gene_list <- sort(unique(c(
  "P21279", "P21278", "P51432", "P11881", "P63318",
  "P68404", "P0DP26", "P0DP27", "P0DP28", "P11798",
  "P28652", "P47937", "P47713"
)))

log2fc_files <- list.files(
  path = file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/Datasets/mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

log2fc_list <- lapply(log2fc_files, function(file) {
  read_csv(file) %>% 
    mutate(
      Comparison = tools::file_path_sans_ext(basename(file))
    )
})

# Combine all into one long data frame
log2fc_long <- bind_rows(log2fc_list)

# Step 2: Filter for genes of interest (optional)
log2fc_filtered <- log2fc_long %>%
  filter(gene_symbol %in% gene_list)  # replace with your list of genes

# Step 3: Summarize per gene per comparison
gene_fc_summary <- log2fc_filtered %>%
  group_by(gene_symbol, Comparison) %>%
  summarise(
    count = n(),
    mean_log2fc = mean(log2fc, na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: Plot
max_abs_fc <- max(abs(gene_fc_summary$mean_log2fc), na.rm = TRUE)

plot_gene_fc <- ggplot(gene_fc_summary, aes(
  x = Comparison,
  y = gene_symbol,
  size = count,
  fill = mean_log2fc
)) +
  geom_point(shape = 21, stroke = 0) +
  scale_size_continuous(guide = "none", range = c(3, 10)) +  # Remove Occurrences legend
  scale_fill_gradient2(
    name = "Mean log2FC",
    low = "#6698CC", mid = "white", high = "#F08C21",
    midpoint = 0, limits = c(-max_abs_fc, max_abs_fc)
  ) +
  labs(
    title = "Gene-Centric log2FC Across Comparisons",
    subtitle = paste("Profiling:", ensemble_profiling, "| Condition:", condition),
    x = "Comparison",
    y = "Gene"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right"
  )

# Save and show
fc_plot_file <- file.path(core_enrichment_dir, paste0("NK3R_gene_centric_log2fc_", ensemble_profiling, "_", condition, ".svg"))
ggsave(fc_plot_file, plot = plot_gene_fc, width = 5, height = 7, dpi = 300)
print(plot_gene_fc)

# save the summary table to a CSV file
summary_file <- file.path(core_enrichment_dir, paste0("NK3R_gene_centric_log2fc_summary_", ensemble_profiling, "_", condition, ".xlsx"))
writexl::write_xlsx(gene_fc_summary, path = summary_file)