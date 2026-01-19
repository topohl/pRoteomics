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

# -----------------------------------------------------
# Define Paths and Project Directories
# -----------------------------------------------------

#' Configure Analysis Parameters and Directory Structure
#'
#' @description
#' This section sets up the analysis parameters and directory structure for
#' Gene Ontology (GO) enrichment analysis using clusterProfiler.

# Set analysis parameters ---------------------------------------------

# Define the Gene Ontology domain (e.g., MF, BP, or CC)
ont <- "BP"  # Biological Process

# Define the ensemble profiling method used in the analysis
ensemble_profiling <- "learning_signature"

# Specify the experimental condition (e.g., CNO, VEH, CS, US, effects_inhibition_memory_ensemble, or learning_signature)
condition <- "effects_inhibition_memory_ensemble"

# Base Path Definition (Centralized for easy modification)
base_project_path <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"

# Set up working environment ------------------------------------------

# Set the working directory to the project’s base folder
setwd(base_project_path)

uniprot_mapping_file_path <- file.path(
  base_project_path,
  "Datasets",
  "MOUSE_10090_idmapping.dat"
)

# Load required packages
if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
library(clusterProfiler)
library(org.Mm.eg.db)

# Read the UniProt mapping file and extract only Uniprot-to-GeneName mappings
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  filter(V2 == "Gene_Name") %>%
  dplyr::select(UniprotID = V1, Gene_Name = V3)

# Import input files --------------------------------------------------

# List all CSV files from the specified core_enrichment subfolder
input_dir <- file.path(base_project_path, "Datasets", "core_enrichment", ont, ensemble_profiling, condition)

file_paths <- list.files(
  path = input_dir,
  pattern = "*.csv",
  full.names = TRUE
)
# Name each file entry using the base filename (without the .csv extension)
names(file_paths) <- basename(file_paths) %>% str_remove(".csv")

# Define output directories -------------------------------------------

# Main Output Directory
main_output_dir <- file.path(base_project_path, "Results", "compareGO", ont, ensemble_profiling, condition)
dir.create(main_output_dir, showWarnings = FALSE, recursive = TRUE)

# Subdirectories for organized storage
subdirs <- list(
  tables = file.path(main_output_dir, "01_Tables_and_Data"),
  plots_main = file.path(main_output_dir, "02_Main_Plots"),
  core_enrichment_plots = file.path(main_output_dir, "03_Core_Enrichment_Heatmaps"),
  gene_lists = file.path(main_output_dir, "04_Gene_Lists"),
  volcanoes = file.path(main_output_dir, "05_Volcano_Plots"),
  sig_proteins = file.path(main_output_dir, "06_Significant_Proteins"),
  go_enrichment = file.path(main_output_dir, "07_Regulated_Protein_GO"),
  gene_centric = file.path(main_output_dir, "08_Gene_Centric_Analysis")
)

# Create all subdirectories
lapply(subdirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------
# Read and Combine Data
# -----------------------------------------------------
#' Read and Combine Enrichment Data

enrichment_list <- lapply(file_paths, read.csv)
names(enrichment_list) <- names(file_paths)

combined_df <- bind_rows(
  lapply(names(enrichment_list), function(name) {
  df <- enrichment_list[[name]]
  df$Comparison <- name
  return(df)
  })
)

# Save combined dataframe to Excel
writexl::write_xlsx(
  combined_df, 
  path = file.path(subdirs$tables, paste0("Combined_Enrichment_Data.xlsx"))
)

# Read in go_ids 
gene_go_ids <- combined_df %>%
  # Select just the columns needed
  dplyr::select(GO_ID = ID, core_enrichment) %>%
  # Split the slash-separated gene accessions into rows
  mutate(core_enrichment = strsplit(core_enrichment, "/")) %>%
  unnest(core_enrichment) %>%
  dplyr::rename(Gene = core_enrichment) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(GO_IDs = paste(unique(GO_ID), collapse = "; "), .groups = "drop")

# read in p-values
p_values <- combined_df %>%
  dplyr::select(Comparison, Description, p.adjust) %>%
  distinct()

# -----------------------------------------------------
# Select Top Terms
# -----------------------------------------------------
#' Top Terms Selection for GO Batch Comparison

# Define filtering parameters for selecting top GO terms
significant_only <- TRUE  # If TRUE, only keep terms with p.adjust < 0.05
top10_terms <- TRUE
top10_per_comp <- TRUE    # If TRUE, select the top 10 terms based on absolute NES for each comparison

# Filter and select top terms from the combined data frame for the NON-REFINED (Standard) version
if (top10_per_comp) {
  if (significant_only) {
    # Only include significant terms (p.adjust < 0.05), then select the top 10 per comparison with the highest absolute NES
    top_df_per_comp <- combined_df %>%
      filter(p.adjust < 0.05) %>%
      group_by(Comparison) %>%
      arrange(desc(abs(NES))) %>%
      dplyr::slice_head(n = 10) %>%
      ungroup()
  } else {
    # Without significance filtering, select the top 10 terms per comparison with the highest absolute NES
    top_df_per_comp <- combined_df %>%
      group_by(Comparison) %>%
      arrange(desc(abs(NES))) %>%
      dplyr::slice_head(n = 10) %>%
      ungroup()
  }
} else {
  if (significant_only) {
    # If not selecting top10, just retain all significant terms (p.adjust < 0.05)
    top_df_per_comp <- combined_df %>%
      filter(p.adjust < 0.05)
  } else {
    # Include all terms without any filtering
    top_df_per_comp <- combined_df
  }
}

# -------------------------------------------------------------------------
# Define Candidate Pool for Overall Top Terms
# -------------------------------------------------------------------------

# Create a master candidate pool sorted by absolute NES
if (significant_only) {
  candidate_pool_df <- combined_df %>%
    filter(p.adjust < 0.05) %>%
    arrange(desc(abs(NES)))
} else {
  candidate_pool_df <- combined_df %>%
    arrange(desc(abs(NES)))
}

# -------------------------------------------------------------------------
# 1. Non-Refined Selection (Standard Top N)
# -------------------------------------------------------------------------
if (top10_terms) {
    top_df_standard <- candidate_pool_df %>%
      dplyr::slice_head(n = 10)
} else {
    top_df_standard <- candidate_pool_df
}

# Define the standard list of top terms
top_terms_standard <- unique(top_df_standard$Description)

# -------------------------------------------------------------------------
# 2. Refined Selection: Check Gene Overlap & Remove Redundancy
# -------------------------------------------------------------------------
#' Check for complete protein overlap in Top Terms and fill up from candidates.

# Define target pool size
target_n_terms <- 10 

# Select unique descriptions in order of valid candidates
unique_candidates <- unique(candidate_pool_df$Description)

refined_terms <- character(0)
term_signatures <- list() # Store sorted gene vectors for accepted terms
redundancy_log <- data.frame(
  Rejected_Term = character(),
  Overlapping_With = character(),
  Jaccard_Index = numeric(),
  stringsAsFactors = FALSE
)

# Helper: Get genes from the most significant instance of a term
get_term_genes <- function(desc, df) {
  # Take the row with max abs(NES) for this term to represent its core contents
  row_data <- df %>% 
    filter(Description == desc) %>% 
    arrange(desc(abs(NES))) %>% 
    dplyr::slice(1)
  
  genes <- unlist(strsplit(as.character(row_data$core_enrichment), "/"))
  return(sort(unique(genes)))
}

for (term in unique_candidates) {
  # Stop if we hit the target
  if (length(refined_terms) >= target_n_terms) break
  
  # Get genes
  current_genes <- get_term_genes(term, candidate_pool_df)
  
  is_redundant <- FALSE
  overlap_term <- NA
  score <- 0
  
  # Compare with accepted terms
  for (accepted in refined_terms) {
    accepted_genes <- term_signatures[[accepted]]
    
    # Calculate overlap (Jaccard)
    intersect_len <- length(intersect(current_genes, accepted_genes))
    union_len <- length(union(current_genes, accepted_genes))
    jaccard <- ifelse(union_len > 0, intersect_len / union_len, 0)
    
    # Check for heavy overlap (e.g. > 0.9 or 1.0)
    if (jaccard > 0.7) {
      is_redundant <- TRUE
      overlap_term <- accepted
      score <- jaccard
      break
    }
  }
  
  if (is_redundant) {
    redundancy_log <- rbind(redundancy_log, data.frame(
      Rejected_Term = term,
      Overlapping_With = overlap_term,
      Jaccard_Index = score
    ))
  } else {
    refined_terms <- c(refined_terms, term)
    term_signatures[[term]] <- current_genes
  }
}

# Save redundancy log
writexl::write_xlsx(
  redundancy_log, 
  path = file.path(subdirs$tables, "Redundant_Genesets_Log.xlsx")
)

# Define the refined list of top terms
top_terms_refined <- refined_terms

message(paste("Refined selection: Kept", length(top_terms_refined), "non-redundant terms."))
message(paste("Standard selection: Kept", length(top_terms_standard), "terms."))

# -------------------------------------------------------------------------
# Set Default "top_terms" for downstream plotting
# -------------------------------------------------------------------------

# You can choose here which one becomes the 'default' for the main heatmaps/dotplots below.
# Setting it to refined for cleaner plots:
top_terms <- top_terms_refined

# -------------------------------------------------------------------------
# Optional: Generate Comparison Plot (Refined vs Standard)
# -------------------------------------------------------------------------

# Filter data for Standard
plot_data_standard <- combined_df %>%
  filter(Description %in% top_terms_standard) %>%
  mutate(Selection_Type = "Standard")

# Filter data for Refined
plot_data_refined <- combined_df %>%
  filter(Description %in% top_terms_refined) %>%
  mutate(Selection_Type = "Refined")

# Combine for a check plot
comparison_plot_data <- bind_rows(plot_data_standard, plot_data_refined) %>%
  mutate(Comparison = factor(Comparison)) # Ensure factor

# Create a comparison dotplot
# Only plot unique descriptions per method to keep it readable
comp_plot <- ggplot(comparison_plot_data, aes(
    x = Comparison,
    y = reorder(Description, NES, FUN = mean),
    color = NES,
    size = -log10(p.adjust)
  )) +
  geom_point() +
  facet_grid(Selection_Type ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradientn(
    colours = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
    name = "NES"
  ) +
  labs(
    title = "Comparison of Top Term Selection Methods",
    y = "GO Term",
    x = "Comparison"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the comparison plot
ggsave(
  filename = file.path(subdirs$plots_main, "Comparison_Selection_Methods_Dotplot.svg"),
  plot = comp_plot,
  width = 12,
  height = 10
)

# Create a separate dotplot for Refined Selection only
refined_plot <- ggplot(plot_data_refined, aes(
  x = Comparison,
  y = reorder(Description, NES, FUN = mean),
  color = NES,
  size = -log10(p.adjust)
  )) +
  geom_point(alpha = 0.85) +
  scale_color_gradientn(
  colours = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
  name = "NES",
  limits = c(min(plot_data_refined$NES, na.rm = TRUE), max(plot_data_refined$NES, na.rm = TRUE)),
  values = scales::rescale(c(min(plot_data_refined$NES, na.rm = TRUE), 0, max(plot_data_refined$NES, na.rm = TRUE)))
  ) +
  scale_size_continuous(
  name = expression(-log[10](p.adjust)),
  range = c(1, 10),
  limits = c(0, NA)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.29, 0.29)), drop = FALSE) +
  scale_y_discrete(
  expand = expansion(add = c(0.6, 0.6)),
  labels = scales::label_wrap(40)
  ) +
  labs(
  title = "Refined Top Terms Selection - GO Enrichment Dot Plot",
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
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
  legend.title = element_text(face = "bold", size = 14),
  legend.text = element_text(size = 12),
  legend.position = "right",
  plot.margin = margin(5, 5, 5, 5),
  panel.spacing.x = unit(0.5, "lines")
  ) +
  coord_cartesian(clip = 'off')

# Dynamically adjust height relative to number of terms
num_refined_terms <- length(unique(plot_data_refined$Description))
dynamic_height_refined <- max(5, num_refined_terms * 0.7)
dynamic_width_refined <- max(6, length(unique(plot_data_refined$Comparison)) * 1.9)

# Save the refined plot separately
ggsave(
  filename = file.path(subdirs$plots_main, "Refined_Selection_Dotplot.svg"),
  plot = refined_plot,
  width = dynamic_width_refined,
  height = dynamic_height_refined,
  dpi = 300
)

# -----------------------------------------------------
# Filter Data for Heatmap
# -----------------------------------------------------
#' Process and Reorder Comparison Data for Heatmap

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
heatmap_data_export <- tibble::rownames_to_column(as.data.frame(heatmap_data), var = "RowNames")

# Replace empty fields with NA
heatmap_data_export[heatmap_data_export == ""] <- NA
writexl::write_xlsx(heatmap_data_export, path = file.path(subdirs$tables, "Matrix_Heatmap_NES.xlsx"))

# Similarly, create a corresponding matrix for significance labels
heatmap_labels <- lookup_df %>%
  dplyr::select(Description, Comparison, sig_label) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = sig_label) %>%
  tibble::column_to_rownames("Description")

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

# Define output directory and filename
output_file <- file.path(subdirs$plots_main, "Heatmap_Enrichment_Comparisons.svg")

# Save the heatmap to an SVG file with clean font and correct height
svg(output_file, width = 10, height = plot_height, family = "Arial")
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# -----------------------------------------------------
# Create Dot Plot
# -----------------------------------------------------
#' Comparative Gene Ontology Enrichment Dot Plot

# Construct dot plot for comparative GO enrichment analysis.
# Create full dotplot using top_df_per_comp (filtered per comparison)
top_df_per_comp_plot <- top_df_per_comp %>%
  mutate(Comparison = factor(Comparison, levels = comparison_order))

# Get all comparisons from comparison_order
all_comparisons <- comparison_order

# Filter combined_df to include all data points for ALL comparisons
full_data_for_plot <- combined_df %>%
  filter(Description %in% unique(top_df_per_comp_plot$Description)) %>%
  mutate(Comparison = factor(Comparison, levels = all_comparisons))

dotplot <- ggplot(full_data_for_plot, aes(
  x = Comparison,
  y = reorder(Description, NES, FUN = median),
  color = NES,
  size = -log10(p.adjust)
)) +
  geom_point(alpha = 0.85) +
  scale_color_gradientn(
  colours = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
  name = "NES",
  limits = c(min(full_data_for_plot$NES, na.rm = TRUE), max(full_data_for_plot$NES, na.rm = TRUE)),
  values = scales::rescale(c(min(full_data_for_plot$NES, na.rm = TRUE), 0, max(full_data_for_plot$NES, na.rm = TRUE)))
  ) +
  scale_size_continuous(
  name = expression(-log[10](p.adjust)),
  range = c(1, 10),
  limits = c(0, NA)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.29, 0.29)), drop = FALSE) +
  scale_y_discrete(
  expand = expansion(add = c(0.6, 0.6)),
  labels = scales::label_wrap(40)
  ) +
  labs(
  title = "Top Terms Per Comparison - GO Enrichment Dot Plot",
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
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
  legend.title = element_text(face = "bold", size = 14),
  legend.text = element_text(size = 12),
  legend.position = "right",
  plot.margin = margin(5, 5, 5, 5),
  panel.spacing.x = unit(0.5, "lines")
  ) +
  coord_cartesian(clip = 'off')

# Create a second dot plot showing only top terms using lookup_df (filtered)
top_terms_all_data <- combined_df %>%
  filter(Description %in% top_terms) %>%
  mutate(Comparison = factor(Comparison, levels = all_comparisons))

dotplot_top <- ggplot(top_terms_all_data, aes(
  x = Comparison,
  y = reorder(Description, NES, FUN = median),
  color = NES,
  size = -log10(p.adjust)
)) +
  geom_point(alpha = 0.85) +
  scale_color_gradientn(
  colours = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
  name = "NES",
  limits = c(min(top_terms_all_data$NES, na.rm = TRUE), max(top_terms_all_data$NES, na.rm = TRUE)),
  values = scales::rescale(c(min(top_terms_all_data$NES, na.rm = TRUE), 0, max(top_terms_all_data$NES, na.rm = TRUE)))
  ) +
  scale_size_continuous(
  name = expression(-log[10](p.adjust)),
  range = c(1, 10),
  limits = c(0, NA)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.29, 0.29)), drop = FALSE) +
  scale_y_discrete(
  expand = expansion(add = c(0.6, 0.6)),
  labels = scales::label_wrap(40)
  ) +
  labs(
  title = "Top Terms - Comparative GO Enrichment Dot Plot",
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
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
  legend.title = element_text(face = "bold", size = 14),
  legend.text = element_text(size = 12),
  legend.position = "right",
  plot.margin = margin(5, 5, 5, 5),
  panel.spacing.x = unit(0.5, "lines")
  ) +
  coord_cartesian(clip = 'off')

# Dynamically adjust the plot width based on the number of unique comparisons.
num_comparisons <- length(unique(lookup_df$Comparison))
dynamic_width <- max(5, num_comparisons * 1.9)

# Dynamically adjust plot height for each dotplot
num_gene_sets_per_comp <- length(unique(top_df_per_comp_plot$Description))
dynamic_height_per_comp <- max(3.7, num_gene_sets_per_comp * 0.5)

num_gene_sets_top <- length(unique(lookup_df$Description))
dynamic_height_top <- max(3.7, num_gene_sets_top * 0.7)

# Save the dot plots as SVG files
output_dotplot <- file.path(subdirs$plots_main, "Dotplot_Enrichment_TopGenes_PerComp.svg")
output_dotplot_top <- file.path(subdirs$plots_main, "Dotplot_Enrichment_TopGenes_Overall.svg")
ggsave(output_dotplot, plot = dotplot, width = dynamic_width, height = dynamic_height_per_comp, dpi = 300)
ggsave(output_dotplot_top, plot = dotplot_top, width = dynamic_width, height = dynamic_height_top, dpi = 300)

# -----------------------------------------------------
# Specify the target genes for the enrichment analysis
# -----------------------------------------------------

# Define the target gene identifiers for enrichment analysis
gene_list <- c("P55099", "Q6NXX1")

# Split the 'core_enrichment' field
long_df <- combined_df %>%
  mutate(core_gene = strsplit(as.character(core_enrichment), "/|;|,|\\s+")) %>%
  unnest(core_gene)

# Retain only the rows corresponding to the genes of interest
filtered_df <- long_df %>%
  filter(core_gene %in% gene_list) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# Compute the minimum and maximum NES
nes_min <- min(filtered_df$NES, na.rm = TRUE)
nes_max <- max(filtered_df$NES, na.rm = TRUE)
max_abs_nes <- max(abs(nes_min), abs(nes_max))

# Construct a ggplot object
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
  values = c(0, 0.5, 1),
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

# Define the output file path
output_plot_file <- file.path(subdirs$plots_main, "Plot_Selected_Genes.svg")
ggsave(output_plot_file, plot = plot_selected_genes, width = 10, height = 7, dpi = 300)

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
  max_padj = min(p.adjust, na.rm = TRUE),
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
gene_summary_file <- file.path(subdirs$gene_centric, "Gene_Centric_Summary.svg")
ggsave(gene_summary_file, plot = plot_gene_summary, width = 5, height = 4, dpi = 300)


# -----------------------------------------------------
# Extract Core Genes per Comparison
# -----------------------------------------------------
#' Core Gene Extraction and CSV File Export

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

# Save the full core gene table across all comparisons to a CSV file
write.csv(core_genes_df,
  file = file.path(subdirs$gene_lists, "Core_Genes_All_Comparisons.csv"),
  row.names = FALSE)

# Group the core genes by Comparison and Description, and then export each group as a separate CSV file.
core_genes_df %>%
  group_by(Comparison, Description) %>%
  summarise(Genes = list(unique(core_enrichment)), .groups = "drop") %>%
  rowwise() %>%
  mutate(file_name = file.path(
  subdirs$gene_lists,
  paste0("CoreGene_", Comparison, "_", substr(make.names(Description), 1, 50), ".csv")
  )) %>%
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
jaccard_heatmap <- pheatmap(
  jaccard_matrix,
  main = "Jaccard Similarity of Core Genes per Comparison",
  color = colorRampPalette(c("white", "blue"))(100),
  display_numbers = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

# Save the Jaccard heatmap to an SVG file
jaccard_sim_matrix <- file.path(subdirs$plots_main, "Heatmap_Jaccard_Similarity.svg")
svg(jaccard_sim_matrix, width = 10, height = 8)
grid::grid.draw(jaccard_heatmap$gtable)
dev.off()

# -----------------------------------------------------
# Expand Core Enrichment Genes for Heatmap
# -----------------------------------------------------
#' Generate and Save Core Enrichment Heatmap

# Expand each enrichment file to get one gene per row
core_long_df <- bind_rows(
  lapply(names(enrichment_list), function(name) {
  df <- enrichment_list[[name]]
  df %>%
  dplyr::select(Description, NES, core_enrichment) %>%
  mutate(core_enrichment = strsplit(core_enrichment, "/")) %>%
  tidyr::unnest(core_enrichment) %>%
  dplyr::rename(Gene = core_enrichment) %>%
  mutate(Comparison = name)
  })
)

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
max_abs_nes <- max(abs(nes_min), abs(nes_max))
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
output_file <- file.path(subdirs$plots_main, "Heatmap_Overall_Core_Enrichment.svg")
svg(output_file, width = 10, height = 8)
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# Reshape to long format
long_heatmap_df <- heatmap_df %>%
  pivot_longer(-Gene, names_to = "Comparison", values_to = "NES")

# Top 25 upregulated genes
top25_up <- long_heatmap_df %>%
  group_by(Gene) %>%
  summarize(max_nes = max(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_nes)) %>%
  dplyr::slice_head(n = 25) %>%
  mutate(Direction = "Up")

# Top 25 downregulated genes
top25_down <- long_heatmap_df %>%
  group_by(Gene) %>%
  summarize(min_nes = min(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_nes) %>%
  dplyr::slice_head(n = 25) %>%
  mutate(Direction = "Down")

# Combine without deduplication
top_bottom_genes <- bind_rows(top25_up, top25_down)

# Create a mapping of each gene to its GO terms (Description)
gene_go <- core_long_df %>%
  dplyr::select(Gene, Description) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(GO_Terms = paste(unique(Description), collapse = "; "), .groups = "drop")

# Add the GO term column to the gene list
top_bottom_genes <- top_bottom_genes %>%
  left_join(gene_go, by = "Gene")

# Save gene list with direction and GO term
write_xlsx(
  top_bottom_genes,
  file.path(subdirs$gene_lists, "TopBottom_GeneList.xlsx")
)

# Filter the heatmap matrix to only include the top and bottom 25 genes
top_bottom_matrix <- heatmap_matrix[rownames(heatmap_matrix) %in% top_bottom_genes$Gene, , drop = FALSE]

# Order rows
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
  dplyr::select(UniprotID = V1, Gene_Name = V3)

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
max_abs_nes <- max(abs(nes_min), abs(nes_max))
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
svg(file.path(subdirs$plots_main, "Heatmap_TopBottom_Enrichment.svg"), width = 5, height = 10)
grid::grid.draw(top_bottom_plot$gtable)
dev.off()

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
  dplyr::rename(UniprotID = V1)

# Extract Gene Synonyms
gene_synonyms <- uniprot_df %>%
  filter(V2 == "Gene_Synonym") %>%
  distinct(V1, .keep_all = TRUE) %>%
  dplyr::rename(UniprotID = V1, Gene_Synonym = V3)

# Optional: Gene Descriptions
gene_descriptions <- core_long_df %>%
  dplyr::select(Gene, Description) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(Description = paste(unique(Description), collapse = "; "), .groups = "drop")

# List log2fc input files
log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

# Read and combine all log2fc files
log2fc_df_list <- lapply(log2fc_files, function(f) {
  df <- read_csv(f, col_types = cols())
  df$Comparison <- tools::file_path_sans_ext(basename(f))
  df <- df %>%
  dplyr::rename(
  padj = adj.P.Val,
  log2fc = logFC
  ) %>%
  as.data.frame()
  return(df)
})

log2fc_long_df <- bind_rows(log2fc_df_list)

# Unique comparisons
unique_comparisons <- unique(log2fc_long_df$Comparison)

# Loop over each comparison
for (comp in unique_comparisons) {
  full_comp_df <- log2fc_long_df %>%
  filter(Comparison == comp, !is.na(log2fc), !is.na(padj))

# Volcano plot section
volcano_df <- full_comp_df %>%
  filter(!str_detect(gene_symbol, "_")) %>%
  left_join(uniprot_subset %>% dplyr::select(UniprotID, Gene_Name),
    by = c("gene_symbol" = "UniprotID")) %>%
  mutate(Significance = case_when(
  padj < 0.05 & log2fc > 0 ~ "up",
  padj < 0.05 & log2fc < 0 ~ "down",
  TRUE ~ "n.s."
  ))
  
# Identify top 5 up and down regulated genes
top5_up <- volcano_df %>%
  filter(Significance == "up") %>%
  arrange(desc(log2fc)) %>%
  dplyr::slice_head(n = 5)

top5_down <- volcano_df %>%
  filter(Significance == "down") %>%
  arrange(log2fc) %>%
  dplyr::slice_head(n = 5)

top_genes_to_label <- bind_rows(top5_up, top5_down)
# Calculate max vals
max_val <- max(
  abs(volcano_df$log2fc),
  abs(-log10(volcano_df$padj)),
  0,
  na.rm = TRUE
)
limit <- ceiling(max_val) + 1
max_x_val <- max(abs(volcano_df$log2fc), 0, na.rm = TRUE)
x_limit <- ceiling(max_x_val) + 1
max_y_val <- max(abs(-log10(volcano_df$padj)), 0, na.rm = TRUE)
y_limit <- ceiling(max_y_val) + 1

volcano_plot <- ggplot(volcano_df, aes(x = log2fc, y = -log10(padj), label = Gene_Name, color = Significance)) +
  geom_point(shape = 16, alpha = 0.6, size = 6) +
  ggrepel::geom_text_repel(
  data = top_genes_to_label,
  aes(label = Gene_Name),
  size = 9,
  max.overlaps = Inf,
  show.legend = FALSE
  ) +
  scale_color_manual(
  values = c("up" = "#f36d07", "down" = "#455A64", "n.s." = "#c7c7c7"),
  breaks = c("up", "down", "n.s.")
  ) +
  scale_x_continuous(limits = c(-x_limit, x_limit), breaks = seq(-x_limit, x_limit, by = 1)) +
  scale_y_continuous(limits = c(0, y_limit), breaks = scales::pretty_breaks()) +
  labs(
  title = paste(strsplit(comp, "_")[[1]], collapse = " over "),
  x = "log2 Fold Change",
  y = expression(-log[10]~(p~adj.))
  ) +
  theme_classic(base_size = 24) +
  theme(
  legend.position = c(1, 1),
  legend.justification = c(1, 1),
  legend.background = element_rect(fill = "white", color = NA),
  legend.title = element_blank(),
  legend.text = element_text(size = 24),
  plot.title = element_text(size = 28, face = "bold"),
  axis.title = element_text(face = "bold", size = 26),
  axis.text = element_text(size = 24, color = "black"),
  axis.line = element_blank(),
  axis.ticks = element_line(color = "black", size = 1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#455A64")

# Save volcano plot
volcano_file <- file.path(subdirs$volcanoes, paste0("Volcano_", comp, ".svg"))
message(paste("Saving volcano plot to:", volcano_file))
ggsave(
  filename = volcano_file,
  plot = volcano_plot,
  width = 6,
  height = 6
)

# Now check for significant data for heatmap + Excel
comp_df <- full_comp_df %>% 
  filter(padj < 0.05) %>%
  ungroup()

if (nrow(comp_df) == 0) {
  message(paste("  No significant genes found for", comp, "- skipping heatmap and Excel."))
  next
}

# Select top and bottom 25 based on log2fc direction
comp_df <- comp_df %>%
  mutate(Direction = case_when(
  log2fc > 0 ~ "Up",
  log2fc < 0 ~ "Down",
  TRUE ~ "Neutral"
  ))

# Combine top 25 up and down
top_bottom_genes <- bind_rows(
  comp_df %>%
  filter(Direction == "Up") %>%
  arrange(desc(log2fc)) %>%
  dplyr::slice_head(n = 25),
  comp_df %>%
  filter(Direction == "Down") %>%
  arrange(log2fc) %>%
  dplyr::slice_head(n = 25)
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
  Gene_Label = ifelse(is.na(Gene_Name) | Gene_Name == "", gene_symbol, Gene_Name),
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
  dplyr::select(Gene_Label, log2fc) %>%
  tibble::column_to_rownames("Gene_Label") %>%
  as.matrix()

if (nrow(heatmap_matrix) >= 2) {
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
  border_color = "#e7e7e7",
  silent = TRUE
  )

  # Save heatmap
  heatmap_file <- file.path(subdirs$plots_main, paste0("Heatmap_Log2FC_", comp, ".svg"))
  message(paste("  Saving heatmap to:", heatmap_file))
  svg(heatmap_file, width = 3, height = 8)
  grid::grid.draw(heatmap_plot$gtable)
  dev.off()
}

# Save Excel with annotation (Top/Bottom summary)
excel_out <- mapped_top_bottom %>%
  dplyr::select(
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

excel_file <- file.path(subdirs$sig_proteins, paste0("SigProteins_TopBottom_", comp, ".xlsx"))
write_xlsx(excel_out, excel_file)

# Save separate Excel files for ALL significantly upregulated and downregulated proteins
# Prepare full annotated data for all significant proteins
all_sig_proteins <- comp_df %>%
  left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
  left_join(gene_synonyms, by = c("gene_symbol" = "UniprotID")) %>%
  left_join(gene_descriptions, by = c("gene_symbol" = "Gene")) %>%
  left_join(gene_go_ids, by = c("gene_symbol" = "Gene")) %>%
  mutate(
  Gene_Name = as.character(Gene_Name),
  `UniProtKB-ID` = as.character(`UniProtKB-ID`),
  Gene_Synonym = as.character(Gene_Synonym),
  abs_log2fc = abs(log2fc),
  UniProt_Link = paste0("https://www.uniprot.org/uniprot/", gene_symbol)
  )

# Save upregulated proteins
upregulated <- all_sig_proteins %>%
  filter(Direction == "Up") %>%
  arrange(desc(log2fc)) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(
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

if (nrow(upregulated) > 0) {
  up_file <- file.path(subdirs$sig_proteins, paste0("SigProteins_Upregulated_", comp, ".xlsx"))
  write_xlsx(upregulated, up_file)
}

# Save downregulated proteins
downregulated <- all_sig_proteins %>%
  filter(Direction == "Down") %>%
  arrange(log2fc) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(
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

if (nrow(downregulated) > 0) {
  down_file <- file.path(subdirs$sig_proteins, paste0("SigProteins_Downregulated_", comp, ".xlsx"))
  write_xlsx(downregulated, down_file)
}
}

# -----------------------------------------------------
# Run GO Enrichment on Up/Downregulated Proteins
# -----------------------------------------------------
#' Run GO Enrichment Analysis on Significantly Regulated Proteins

# Initialize a data frame to log the results
enrichment_summary_log <- data.frame(
  Comparison = character(),
  Direction = character(),
  Input_Genes_Count = integer(),
  Enriched_Terms_Count = integer(),
  Status = character(),
  Saved_CSV = character(),
  Saved_Plot = character(),
  Timestamp = character(),
  stringsAsFactors = FALSE
)

# Get list of comparisons from the significant proteins directory
comparison_files <- list.files(
  path = subdirs$sig_proteins,
  pattern = "SigProteins_(Upregulated|Downregulated)_.*\\.xlsx$",
  full.names = TRUE
)

if (length(comparison_files) > 0) {
  # Extract unique comparison names
  comparisons <- unique(gsub("SigProteins_(Upregulated|Downregulated)_|\\.xlsx", "", 
     basename(comparison_files)))
  
  message("Found ", length(comparisons), " comparisons to process: ", paste(comparisons, collapse = ", "))

  # Loop through each comparison
  for (comp in comparisons) {
  message(paste("Processing GO enrichment for comparison:", comp))
  
  custom_theme <- function() {
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(color = "black", size = 16),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14)
  )
  }
  
  # --- Process Upregulated Proteins ---
  up_file <- file.path(subdirs$sig_proteins, paste0("SigProteins_Upregulated_", comp, ".xlsx"))
  
  log_genes_count <- 0
  log_terms_count <- 0
  log_status <- "File not found"
  log_csv <- NA
  log_plot <- NA
  
  if (file.exists(up_file)) {
  up_genes <- readxl::read_xlsx(up_file) %>%
  pull(Uniprot_Accession) %>%
  unique()
  
  log_genes_count <- length(up_genes)
  
  if (log_genes_count > 0) {
  ego_up_bp <- enrichGO(
    gene = up_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "UNIPROT",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if (!is.null(ego_up_bp) && nrow(as.data.frame(ego_up_bp)) > 0) {
    log_terms_count <- nrow(as.data.frame(ego_up_bp))
    log_status <- "Success"
    
    out_csv <- file.path(subdirs$go_enrichment, paste0("GO_", comp, "_Upregulated_BP.csv"))
    write.csv(as.data.frame(ego_up_bp), out_csv, row.names = FALSE)
    log_csv <- basename(out_csv)
    
    p_dot <- dotplot(ego_up_bp, showCategory = 15) +
    scale_color_gradientn(colours = c("#6698CC", "white", "#F08C21"), name = "p.adjust") +
    scale_size(range = c(5, 12)) +
    labs(title = "GO BP Enrichment - Upregulated", subtitle = paste("Comparison:", comp)) +
    custom_theme()
    
    out_plot <- file.path(subdirs$go_enrichment, paste0("GO_Dotplot_", comp, "_Upregulated_BP.svg"))
    ggsave(out_plot, plot = p_dot, width = 8.5, height = 8)
    log_plot <- basename(out_plot)
  } else {
    log_status <- "No significant enrichment"
  }
  } else {
  log_status <- "No genes in input file"
  }
  }
  
  enrichment_summary_log <- rbind(enrichment_summary_log, data.frame(
  Comparison = comp, Direction = "Upregulated", Input_Genes_Count = log_genes_count,
  Enriched_Terms_Count = log_terms_count, Status = log_status, Saved_CSV = log_csv,
  Saved_Plot = log_plot, Timestamp = as.character(Sys.time()), stringsAsFactors = FALSE
  ))
  
  # --- Process Downregulated Proteins ---
  down_file <- file.path(subdirs$sig_proteins, paste0("SigProteins_Downregulated_", comp, ".xlsx"))
  
  log_genes_count <- 0
  log_terms_count <- 0
  log_status <- "File not found"
  log_csv <- NA
  log_plot <- NA
  
  if (file.exists(down_file)) {
  down_genes <- readxl::read_xlsx(down_file) %>%
  pull(Uniprot_Accession) %>%
  unique()
  
  log_genes_count <- length(down_genes)
  
  if (log_genes_count > 0) {
  ego_down_bp <- enrichGO(
    gene = down_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "UNIPROT",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if (!is.null(ego_down_bp) && nrow(as.data.frame(ego_down_bp)) > 0) {
    log_terms_count <- nrow(as.data.frame(ego_down_bp))
    log_status <- "Success"
    
    out_csv <- file.path(subdirs$go_enrichment, paste0("GO_", comp, "_Downregulated_BP.csv"))
    write.csv(as.data.frame(ego_down_bp), out_csv, row.names = FALSE)
    log_csv <- basename(out_csv)
    
    p_dot <- dotplot(ego_down_bp, showCategory = 15) +
    scale_color_gradientn(colours = c("#6698CC", "white", "#F08C21"), name = "p.adjust") +
    scale_size(range = c(5, 12)) +
    labs(title = "GO BP Enrichment - Downregulated", subtitle = paste("Comparison:", comp)) +
    custom_theme()
    
    out_plot <- file.path(subdirs$go_enrichment, paste0("GO_Dotplot_", comp, "_Downregulated_BP.svg"))
    ggsave(out_plot, plot = p_dot, width = 8.5, height = 8)
    log_plot <- basename(out_plot)
  } else {
    log_status <- "No significant enrichment"
  }
  } else {
  log_status <- "No genes in input file"
  }
  }
  
  enrichment_summary_log <- rbind(enrichment_summary_log, data.frame(
  Comparison = comp, Direction = "Downregulated", Input_Genes_Count = log_genes_count,
  Enriched_Terms_Count = log_terms_count, Status = log_status, Saved_CSV = log_csv,
  Saved_Plot = log_plot, Timestamp = as.character(Sys.time()), stringsAsFactors = FALSE
  ))
  }
  
  # Save the summary log
  log_file_path <- file.path(subdirs$go_enrichment, "Summary_Log_GO_Enrichment.csv")
  write.csv(enrichment_summary_log, log_file_path, row.names = FALSE)
}

# -----------------------------------------------------
# Generate Individual Core Enrichment Heatmaps
# -----------------------------------------------------
#' Generate and Save Core Enrichment Heatmaps

# Read log2fc data for each comparison from the mapped/ensemble_profiling directory

log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

names(log2fc_files) <- basename(log2fc_files) %>%
  str_remove(".csv") %>%
  str_trim()

# Read and combine all log2fc files into one data frame
log2fc_list <- lapply(names(log2fc_files), function(comp) {
  df <- read.csv(log2fc_files[comp])
  df$Comparison <- comp
  df
})

log2fc_df <- bind_rows(log2fc_list)

# Standardize column names
log2fc_df <- log2fc_df %>%
  dplyr::rename(
  log2fc = logFC,
  pvalue = P.Value,
  padj = adj.P.Val
  )

# Generate heatmaps for each GO term in core_long_df
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
  dplyr::select(Gene, Comparison, log2fc) %>%
  group_by(Gene, Comparison) %>%
  summarize(log2fc = mean(log2fc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Comparison, values_from = log2fc) %>%
  column_to_rownames("Gene")

  colnames(matrix_df) <- str_trim(colnames(matrix_df))

  heatmap_matrix <- as.matrix(matrix_df)
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  # Define output filename with a safe name
  safe_desc <- substr(make.names(description), 1, 30)
  filename <- file.path(subdirs$core_enrichment_plots, paste0("Heatmap_Core_", safe_desc, ".svg"))

  cluster_rows_option <- nrow(heatmap_matrix) > 1
  cluster_cols_option <- ncol(heatmap_matrix) > 1

  max_val <- max(abs(heatmap_matrix), na.rm = TRUE)
  breaks <- seq(-max_val, max_val, length.out = 101)

  row_height <- 0.05
  plot_height <- max(8, nrow(heatmap_matrix) * row_height + 0.2 * nrow(heatmap_matrix))

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

# TACR3-signalling CORE
gene_list <- sort(unique(c(
  "P21279", "P21278", "P51432", "P11881", "P63318",
  "P68404", "P0DP26", "P0DP27", "P0DP28", "P11798",
  "P28652", "P47937", "P47713"
)))

log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
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
log2fc_long <- bind_rows(log2fc_list) %>%
  dplyr::rename(
  log2fc = logFC,
  pvalue = P.Value,
  padj = adj.P.Val
  )

# Step 2: Filter for genes of interest (optional)
log2fc_filtered <- log2fc_long %>%
  filter(gene_symbol %in% gene_list)

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
fc_plot_file <- file.path(subdirs$gene_centric, paste0("GeneCentric_Log2FC_Bubble_", condition, ".svg"))
ggsave(fc_plot_file, plot = plot_gene_fc, width = 5, height = 7, dpi = 300)

# save the summary table to a CSV file
summary_file <- file.path(subdirs$gene_centric, paste0("GeneCentric_Log2FC_Summary.xlsx"))
writexl::write_xlsx(gene_fc_summary, path = summary_file)

# -----------------------------------------------------
# Generate Supplementary Tables (Merged)
# -----------------------------------------------------
#' @description
#' Aggregates key findings from different stages of the analysis into a single
#' multi-sheet Excel file suitable for supplements.

# 1. Enrichment Results (All combined data)
supp_enrichment <- combined_df %>%
  dplyr::select(Comparison, ID, Description, setSize, 
     pvalue, p.adjust, qvalue, NES, core_enrichment, rank, leading_edge)

# 2. Significant Proteins (Master list across all comparisons)
# Re-filter the raw full data for significant hits and join gene names, synonyms, and descriptions
supp_sig_proteins <- log2fc_long %>%
  filter(padj < 0.05) %>%
  left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
  left_join(gene_synonyms, by = c("gene_symbol" = "UniprotID")) %>% 
  left_join(gene_descriptions, by = c("gene_symbol" = "Gene")) %>%
  mutate(Direction = ifelse(log2fc > 0, "Up", "Down")) %>%
  dplyr::select(Comparison, UniprotID = gene_symbol, Gene_Name, Gene_Synonym, Description,
     log2fc, pvalue, padj, Direction) %>%
  arrange(Comparison, padj)

# 3. Gene Centric View (from the last section used for bubble plot)
supp_gene_centric <- gene_fc_summary %>%
  left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
  dplyr::rename(Gene_Name = Gene_Name)

# 4. Matrix for Heatmap (NES) - Specific to the Top Terms visualized
supp_nes_matrix <- lookup_df %>%
  dplyr::select(Description, Comparison, NES) %>%
  pivot_wider(names_from = Comparison, values_from = NES)

# 5. Core Genes Binary Matrix (Presence/Absence)
supp_binary_matrix <- as.data.frame(binary_matrix) %>%
  tibble::rownames_to_column(var = "Gene")

# Combine into a list
supp_list <- list(
  "GO_Enrichment_Results" = supp_enrichment,
  "Significant_Proteins" = supp_sig_proteins,
  "Gene_Centric_Log2FC" = supp_gene_centric,
  "NES_Matrix_TopTerms" = supp_nes_matrix,
  "Core_Genes_Binary_Matrix" = supp_binary_matrix
)

# Write to a single Excel file
supp_output_path <- file.path(subdirs$tables, "Supplementary_Data.xlsx")
writexl::write_xlsx(supp_list, path = supp_output_path)

message("Supplementary tables saved to: ", supp_output_path)
