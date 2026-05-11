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
  dplyr::select(UniprotID = V1, Gene_Name = V3)

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
condition <- "DG"

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

# Read files safely and build a combined data frame with a "Comparison" column
if (length(file_paths) == 0) {
  stop("No input files found in the specified core_enrichment directory: ",
       file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/core_enrichment", ont, ensemble_profiling, condition))
}

# Use purrr::imap_dfr to read each CSV and tag with the file base name (the Comparison)
combined_df <- purrr::imap_dfr(file_paths, function(path, name) {
  df <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE),
    error = function(e) {
      warning("Failed to read: ", path, "; error: ", e$message)
      return(tibble::tibble())
    }
  )
  if (nrow(df) == 0) return(df)
  df$Comparison <- name
  df
})

# Coerce p.adjust and NES to numeric if present, replacing non-numeric values with NA
if ("p.adjust" %in% names(combined_df)) {
  combined_df$p.adjust <- suppressWarnings(as.numeric(as.character(combined_df$p.adjust)))
  if (all(is.na(combined_df$p.adjust))) {
    warning("All 'p.adjust' values are NA after coercion; significance filtering will be skipped.")
  }
}

if ("NES" %in% names(combined_df)) {
  combined_df$NES <- suppressWarnings(as.numeric(as.character(combined_df$NES)))
  if (all(is.na(combined_df$NES))) {
    warning("All 'NES' values are NA after coercion; ranking by NES may not work as expected.")
  }
}

# Read in go_ids 
# Safely build gene_go_ids only if the required columns exist; otherwise return an empty tibble and warn.
if (all(c("ID", "core_enrichment") %in% names(combined_df))) {
  gene_go_ids <- combined_df %>%
    # Select just the columns needed
    dplyr::select(GO_ID = ID, core_enrichment) %>%
    # Split the slash-separated gene accessions into rows
    tidyr::separate_rows(core_enrichment, sep = "/") %>%
    dplyr::rename(Gene = core_enrichment) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarize(GO_IDs = paste(unique(GO_ID), collapse = "; "), .groups = "drop")
} else {
  warning("combined_df is missing columns 'ID' and/or 'core_enrichment'; 'gene_go_ids' will be an empty tibble.")
  gene_go_ids <- tibble::tibble(Gene = character(), GO_IDs = character())
}

# read in p-values
# Safely build p_values only if the required columns exist; otherwise return an empty tibble and warn.
if (all(c("Comparison", "Description", "p.adjust") %in% names(combined_df))) {
  p_values <- combined_df %>% dplyr::select(Comparison, Description, p.adjust)
} else {
  warning("combined_df is missing one of 'Comparison', 'Description' or 'p.adjust'; 'p_values' will be an empty tibble.")
  p_values <- tibble::tibble(Comparison = character(), Description = character(), p.adjust = numeric())
}

# Parameters controlling selection of top terms (scalars; change if you want different behavior)
top10_terms <- TRUE
significant_only <- FALSE

# Filter and select top terms from the combined data frame
if (top10_terms) {
    if (significant_only && "p.adjust" %in% names(combined_df)) {
        # only significant terms, then top 10 by absolute NES per Comparison
        top_df <- combined_df %>%
            filter(!is.na(p.adjust) & p.adjust < 0.05) %>%
            group_by(Comparison) %>%
            slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
            ungroup()
    } else {
        # top 10 by absolute NES per Comparison (no significance filter)
        top_df <- combined_df %>%
            group_by(Comparison) %>%
            slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
            ungroup()
    }
} else {
    if (significant_only && "p.adjust" %in% names(combined_df)) {
        # keep all significant terms
        top_df <- combined_df %>%
            filter(!is.na(p.adjust) & p.adjust < 0.05)
    } else {
        # keep all terms
        top_df <- combined_df
    }
}

# Exclude any GO term (Description) for which none of the comparisons are significant.
# Only apply this filtering when p.adjust is available and at least one non-NA value exists.
if ("p.adjust" %in% names(combined_df) && any(!is.na(combined_df$p.adjust))) {
    sig_terms <- combined_df %>%
        group_by(Description) %>%
        summarize(any_significant = any(!is.na(p.adjust) & p.adjust < 0.05), .groups = "drop") %>%
        filter(any_significant) %>%
        pull(Description)

    if (length(sig_terms) == 0) {
        warning("No descriptions have any significant comparisons; keeping top_df unchanged.")
    } else {
        top_df <- top_df %>% filter(Description %in% sig_terms)
    }
} else {
    # p.adjust missing or all NA: cannot determine significance across comparisons, skip exclusion
}

# Fallback: if top_df is empty, fall back to using all available descriptions (prevents downstream empty joins)
if (!exists("top_df") || nrow(top_df) == 0) {
  warning("No top terms were selected (top_df is empty); falling back to all available terms from combined_df.")
  top_terms <- unique(na.omit(combined_df$Description))
} else {
  top_terms <- unique(top_df$Description)
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
