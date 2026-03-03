#' compareGO.r - Comparative Gene Ontology Enrichment Analysis and Visualization
#'
#' @description
#' This script performs comparative GO enrichment analysis across multiple experiments.
#' It:
#'   - Reads multiple CSV files with enrichment data from a specified directory.
#'     Filenames (without extensions) are used as labels for each comparison.
#'   - Combines all data into a single data frame, adding a "Comparison" column.
#'   - Selects top enriched terms per comparison based on absolute NES (Normalized Enrichment Score).
#'   - Filters the combined dataset to include only these top terms for consistent visualization.
#'   - Reorders comparisons in plots based on the maximum absolute NES.
#'   - Adds a significance label ("✱") for gene sets with adjusted p-value < 0.05.
#'
#' @details
#' The script generates:
#'   1. A heatmap showing differential enrichment across comparisons.
#'   2. A dot plot visualizing NES (color) and significance (-log10(p.adjust)) for each term.
#' Additional outputs:
#'   - Core gene lists per term and comparison.
#'   - A binary matrix of core gene presence/absence for Jaccard similarity calculation.
#'   - Expanded core enrichment heatmaps for individual terms.
#'
#' @section File Inputs:
#'   - CSV files from:
#'     "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment"
#'     (pattern "*.csv").
#'
#' @section Outputs:
#'   - Heatmaps and dot plots of enrichment profiles.
#'   - Core gene tables (CSV).
#'   - SVG and PNG files for all plots.
#'
#' @note
#'   Ensure all file paths and required packages are set up before running.
#'
#' @author
#'   Tobias Pohl

# -----------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
# Ensure 'rlang' version >= 1.1.7 is installed before running this script.
if (!requireNamespace("rlang", quietly=TRUE) || packageVersion("rlang") < "1.1.7") {
    stop("Please install 'rlang' version >= 1.1.7 manually before running this script.")
}
if (!requireNamespace("simplifyEnrichment", quietly=TRUE)) {
    BiocManager::install("simplifyEnrichment", force = TRUE)
}
library(simplifyEnrichment)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, stringr, ggpubr, ggthemes, dplyr, tidyr, purrr,
  readr, pheatmap, tibble, tidyverse, RColorBrewer, writexl, scales, ggrepel,
  magick)

# -----------------------------------------------------
# Define Theme and Helper Functions
# -----------------------------------------------------
#' Custom ggplot2 theme for consistent plots
theme_custom <- function(base_size = 8, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(face = "bold", size = base_size + 1),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.8, "lines"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

#' Dynamically calculate plot dimensions based on data content
calc_dims <- function(df_plot) {
  n_cols <- length(unique(as.character(df_plot$Comparison)))
  n_rows <- length(unique(as.character(df_plot$Description)))
  w <- max(4, 4 + (n_cols * 0.4)) 
  h <- max(6, 2 + (n_rows * 0.3))
  return(list(w = w, h = h))
}

# -----------------------------------------------------
# Set Analysis Parameters and Directory Structure
# -----------------------------------------------------

# Gene Ontology domain (MF, BP, or CC)
ont <- "BP"

# Ensemble profiling method
#ensemble_profiling <- "learning_signature"
ensemble_profiling <- "phenotype_within_unit"

# Experimental condition
#condition <- "memory_ensemble"
condition <- "CA1_microglia"

# Base project path for all input/output
#base_project_path <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
base_project_path <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"

# Set working directory
setwd(base_project_path)

uniprot_mapping_file_path <- file.path(
  base_project_path,
  "Datasets",
  "MOUSE_10090_idmapping.dat"
)

# Load clusterProfiler and organism database
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))

# Read UniProt mapping file (UniprotID <-> Gene_Name)
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  filter(V2 == "Gene_Name") %>%
  dplyr::select(UniprotID = V1, Gene_Name = V3)

# -----------------------------------------------------
# Import Enrichment Data Files
# -----------------------------------------------------

# List all enrichment CSV files for the specified condition
input_dir <- file.path(base_project_path, "Datasets", "core_enrichment", ont, ensemble_profiling, condition)
file_paths <- list.files(
  path = input_dir,
  pattern = "*.csv",
  full.names = TRUE
)
names(file_paths) <- basename(file_paths) %>% str_remove(".csv")

# -----------------------------------------------------
# Define Output Directories
# -----------------------------------------------------

main_output_dir <- file.path(base_project_path, "Results", "compareGO", ont, ensemble_profiling, condition)
dir.create(main_output_dir, showWarnings = FALSE, recursive = TRUE)

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
lapply(subdirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------
# Read and Combine Enrichment Data
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

# Remove whitespace from comparison names to avoid duplicates
combined_df$Comparison <- str_trim(combined_df$Comparison)

# Save combined enrichment data to Excel
writexl::write_xlsx(
  combined_df, 
  path = file.path(subdirs$tables, paste0("Combined_Enrichment_Data.xlsx"))
)

# Extract GO IDs for each gene from core_enrichment column
gene_go_ids <- combined_df %>%
  dplyr::select(GO_ID = ID, core_enrichment) %>%
  mutate(core_enrichment = strsplit(core_enrichment, "/")) %>%
  unnest(core_enrichment) %>%
  dplyr::rename(Gene = core_enrichment) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(GO_IDs = paste(unique(GO_ID), collapse = "; "), .groups = "drop")

# -----------------------------------------------------
# Select Top GO Terms for Comparison
# -----------------------------------------------------
library(dplyr)
library(stringr)
library(purrr)
library(writexl)

# =========================================================
# PARAMETERS
# =========================================================

significant_only <- TRUE
target_n_terms <- 5
redundancy_threshold <- 0.7
min_set_size <- 10

# =========================================================
# PREPROCESS LEADING EDGE ONCE
# =========================================================

combined_df <- combined_df %>%
  mutate(
    leading_edge = strsplit(as.character(core_enrichment), "/")
  )

# =========================================================
# SELECT BEST INSTANCE PER TERM PER COMPARISON
# =========================================================

combined_df_best <- combined_df %>%
  group_by(Comparison, Description) %>%
  slice_max(abs(NES), n = 1, with_ties = FALSE) %>%
  ungroup()

# =========================================================
# JACCARD FUNCTION
# =========================================================

jaccard <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a,b)) / length(union(a,b))
}

# =========================================================
# REDUNDANCY FILTER
# =========================================================

select_nonredundant <- function(df, n_target, direction) {

  if (nrow(df) == 0) return(df)

  # Filter set size
  df <- df %>% filter(setSize >= min_set_size)

  if (nrow(df) == 0) return(df)

  # Sort by NES
  if(direction == "up"){
    df <- df %>% arrange(desc(NES))
  } else {
    df <- df %>% arrange(NES)
  }

  kept_rows <- list()
  kept_genes <- list()

  for(i in seq_len(nrow(df))) {

    candidate <- df[i,]
    cand_genes <- unlist(candidate$leading_edge)

    redundant <- FALSE

    for(k in seq_along(kept_genes)){
      if(jaccard(cand_genes, kept_genes[[k]]) > redundancy_threshold){
        redundant <- TRUE
        break
      }
    }

    if(!redundant){
      kept_rows[[length(kept_rows)+1]] <- candidate
      kept_genes[[length(kept_genes)+1]] <- cand_genes
    }

    if(length(kept_rows) >= n_target) break
  }

  bind_rows(kept_rows)
}

# =========================================================
# MAIN SELECTION
# =========================================================


if(significant_only){
  candidate_pool <- combined_df_best %>% filter(p.adjust < 0.05)
} else {
  candidate_pool <- combined_df_best
}


# Improved handling for the case where there are no significant GO terms
if(nrow(candidate_pool) == 0) {
  message("No significant GO terms found. Skipping top term selection and all downstream plots/tables.")
  df_standard <- tibble()
  df_refined <- tibble()
  df_replacements <- tibble()
  top_terms <- character(0)
  top_terms_standard <- character(0)
  top_terms_refined <- character(0)
  # Ensure 'leading_edge' is character to avoid writexl warnings
  df_standard <- df_standard %>% mutate(leading_edge = as.character(NA))
  df_refined <- df_refined %>% mutate(leading_edge = as.character(NA))
  writexl::write_xlsx(
    list(
      Standard = df_standard,
      Refined = df_refined,
      Replacements = df_replacements
    ),
    path = file.path(subdirs$tables, "Redundancy_Check_Results.xlsx")
  )
  # Set empty data for downstream plotting and skip those blocks
  plot_data_standard <- tibble()
  plot_data_refined <- tibble()
  comparison_plot_data <- tibble()
  lookup_df <- tibble()
  suppressWarnings({
    writexl::write_xlsx(list(), path = file.path(subdirs$tables, "Matrix_Heatmap_NES.xlsx"))
  })
} else {
  comparisons <- unique(candidate_pool$Comparison)
  std_list <- list()
  ref_list <- list()
  replacement_tracker <- list()
  for(comp in comparisons){
    message("Processing: ", comp)
    comp_df <- candidate_pool %>% filter(Comparison == comp)
    # -------- STANDARD --------
    up_std <- comp_df %>% filter(NES > 0) %>% arrange(desc(NES)) %>% slice_head(n = target_n_terms)
    down_std <- comp_df %>% filter(NES < 0) %>% arrange(NES) %>% slice_head(n = target_n_terms)
    # -------- REFINED --------
    up_ref <- select_nonredundant(comp_df %>% filter(NES > 0), target_n_terms, "up")
    down_ref <- select_nonredundant(comp_df %>% filter(NES < 0), target_n_terms, "down")
    std_list[[comp]] <- bind_rows(up_std, down_std)
    ref_list[[comp]] <- bind_rows(up_ref, down_ref)
    # -------- Replacement Tracking --------
    dropped <- setdiff(std_list[[comp]]$Description, ref_list[[comp]]$Description)
    added <- setdiff(ref_list[[comp]]$Description, std_list[[comp]]$Description)
    replacement_tracker[[comp]] <- tibble(
      Comparison = comp,
      Dropped = paste(dropped, collapse = "; "),
      Added = paste(added, collapse = "; ")
    )
  }
  df_standard <- bind_rows(std_list)
  df_refined <- bind_rows(ref_list)
  df_replacements <- bind_rows(replacement_tracker)
  top_terms_standard <- unique(df_standard$Description)
  top_terms_refined <- unique(df_refined$Description)
  # Ensure 'leading_edge' is character to avoid writexl warnings
  df_standard <- df_standard %>% mutate(across(leading_edge, as.character))
  df_refined <- df_refined %>% mutate(across(leading_edge, as.character))
  writexl::write_xlsx(
    list(
      Standard = df_standard,
      Refined = df_refined,
      Replacements = df_replacements
    ),
    path = file.path(subdirs$tables, "Redundancy_Check_Results.xlsx")
  )
  top_terms <- unique(df_refined$Description)
}


# -----------------------------------------------------
# Plot Comparison: Standard vs Refined Selection
# -----------------------------------------------------

# Custom color palette: upregulation = #faa51a, downregulation = #4c87c6, 0 = white
custom_palette <- colorRampPalette(c("#3E3C6F", "white", "#E63A47"))

plot_data_standard <- combined_df %>%
  filter(Description %in% top_terms_standard) %>%
  mutate(Selection_Type = "Standard")
plot_data_refined <- combined_df %>%
  filter(Description %in% top_terms_refined) %>%
  mutate(Selection_Type = "Refined")
comparison_plot_data <- bind_rows(plot_data_standard, plot_data_refined) %>%
  mutate(Comparison = factor(Comparison, levels = unique(Comparison)))

# --- Robust plotting: skip if no data or NES is all NA ---
if (nrow(comparison_plot_data) == 0 || all(is.na(comparison_plot_data$NES))) {
  message("No data available for comparison selection plot. Skipping plot generation.")
} else {
  # Calculate ordering based on NES and Significance
  # We define the order based on the "best" instance of each term across comparisons
  term_ordering <- comparison_plot_data %>%
    group_by(Description) %>%
    summarise(
      Max_NES = NES[which.max(abs(NES))], # Use signed max NES
      Best_Padj = p.adjust[which.max(abs(NES))], # P-value associated with max NES
      Dominant_Comparison = Comparison[which.max(abs(NES))],
      .groups = "drop"
    ) %>%
    arrange(Dominant_Comparison, Max_NES, desc(Best_Padj))

  # Apply factor levels
  comparison_plot_data$Description <- factor(
    comparison_plot_data$Description,
    levels = term_ordering$Description
  )

  max_abs_nes <- max(abs(comparison_plot_data$NES), na.rm = TRUE)

  comp_plot <- ggplot(comparison_plot_data, aes(
    x = Comparison,
    y = Description,
    color = NES,
    size = -log10(p.adjust)
  )) +
    geom_point() +
    facet_grid(Selection_Type ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradientn(
      colours = custom_palette(100),
      name = "NES",
      limits = c(-max_abs_nes, max_abs_nes)
    ) +
    scale_size_continuous(range = c(2, 6)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
    labs(
      title = "Selection method comparison",
      y = "GO Term",
      x = "Comparison"
    ) +
    theme_custom(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"))

  n_cols_comp <- length(unique(as.character(comparison_plot_data$Comparison)))
  dims_comp <- list()
  dims_comp$w <- max(4, 6 + (n_cols_comp * 0.6)) 
  dims_comp$h <- calc_dims(comparison_plot_data)$h

  ggsave(
    filename = file.path(subdirs$plots_main, "Comparison_Selection_Methods_Dotplot.svg"),
    plot = comp_plot,
    width = dims_comp$w, 
    height = dims_comp$h
  )
}

# -----------------------------------------------------
# Prepare Data for Heatmap Visualization
# -----------------------------------------------------

lookup_df <- combined_df %>%
  filter(Description %in% top_terms)

comparison_order <- lookup_df %>%
  group_by(Comparison) %>%
  summarize(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_NES)) %>%
  pull(Comparison)
lookup_df <- lookup_df %>%
  mutate(Comparison = factor(Comparison, levels = comparison_order))

# -----------------------------------------------------
# Generate Enrichment Heatmap
# -----------------------------------------------------

lookup_df <- lookup_df %>%
  mutate(sig_label = ifelse(p.adjust < 0.05, "*", ""))

heatmap_data <- lookup_df %>%
  dplyr::select(Description, Comparison, NES) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = NES)

heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$Description
heatmap_data <- heatmap_data[, -which(names(heatmap_data) == "Description")]
heatmap_data <- as.matrix(heatmap_data)
heatmap_data[is.na(heatmap_data)] <- 0 
heatmap_data_export <- tibble::rownames_to_column(as.data.frame(heatmap_data), var = "RowNames")
writexl::write_xlsx(heatmap_data_export, path = file.path(subdirs$tables, "Matrix_Heatmap_NES.xlsx"))
heatmap_labels <- lookup_df %>%
  dplyr::select(Description, Comparison, sig_label) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = sig_label) %>%
  tibble::column_to_rownames("Description")

plot_height <- max(6, nrow(heatmap_data) * 0.3)
# Set NES legend range to -3 to 3
legend_min <- -2.5
legend_max <- 2.5
my_breaks <- seq(legend_min, legend_max, length.out = 100)

# Custom color palette: downregulation = #4c87c6, 0 = white, upregulation = #faa51a
#my_colors <- colorRampPalette(c("#4c87c6", "white", "#faa51a"))(100)
my_colors <- colorRampPalette(c("#3E3C6F", "white", "#E63A47"))(100)

# Determine if clustering is possible (must have >=2 rows/cols)
cluster_rows_opt <- if (nrow(heatmap_data) >= 2) TRUE else FALSE
cluster_cols_opt <- if (ncol(heatmap_data) >= 2) TRUE else FALSE

heatmap_plot <- pheatmap(
  heatmap_data,
  cluster_rows = cluster_rows_opt,
  cluster_cols = cluster_cols_opt,
  display_numbers = heatmap_labels,
  number_color = "black",
  color = my_colors,
  breaks = my_breaks,
  fontsize = 8,
  fontsize_number = 7,
  border_color = "white",
  cellwidth = 20,
  cellheight = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 45,
  silent = TRUE,
  legend = TRUE
)

# Only extract row/col order if clustering was performed
heatmap_row_order <- if (cluster_rows_opt) rownames(heatmap_data)[heatmap_plot$tree_row$order] else rownames(heatmap_data)
heatmap_col_order <- if (cluster_cols_opt) colnames(heatmap_data)[heatmap_plot$tree_col$order] else colnames(heatmap_data)

output_file <- file.path(subdirs$plots_main, "Heatmap_Enrichment_Comparisons.svg")
svg(output_file, width = 6, height = plot_height, family = "sans")
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# -----------------------------------------------------
# Generate Dot Plots for Enrichment Results
# -----------------------------------------------------

# Helper: Order dotplot axes by clustering
order_dotplot <- function(df, desc_col = "Description", comp_col = "Comparison", val_col = "NES") {
  mat <- df %>%
    dplyr::select(all_of(c(desc_col, comp_col, val_col))) %>%
    filter(!is.na(!!sym(val_col))) %>%
    pivot_wider(names_from = all_of(comp_col), values_from = all_of(val_col), values_fill = 0) %>%
    column_to_rownames(desc_col) %>%
    as.matrix()
  # If less than 2 columns or rows, skip clustering and return default order
  if (ncol(mat) < 2 || nrow(mat) < 2) {
    col_order <- colnames(mat)
    row_order <- rownames(mat)
    return(list(row_order = row_order, col_order = col_order))
  }
  hc_col <- hclust(dist(t(mat)))
  col_order <- colnames(mat)[hc_col$order]
  dom_comp <- apply(mat, 1, function(x) colnames(mat)[which.max(abs(x))])
  dom_val  <- apply(mat, 1, function(x) {
    val <- x[which.max(abs(x))]
    return(val)
  })
  ordering_df <- data.frame(
    Term = rownames(mat),
    Comp = dom_comp,
    Val = dom_val
  )
  ordering_df$Comp <- factor(ordering_df$Comp, levels = col_order)
  ordering_df <- ordering_df %>%
    mutate(
      Comp_Index = as.integer(Comp),
      Direction = ifelse(Val > 0, 1, 0),
      Abs_Val = abs(Val)
    ) %>%
    arrange(desc(Comp_Index), Direction, Abs_Val)
  row_order <- ordering_df$Term
  return(list(row_order = row_order, col_order = col_order))
}


## --- Dotplot: Top terms per comparison ---
# Define top_df_per_comp if not already defined
if (!exists("top_df_per_comp")) {
  # Use the standard top terms per comparison (from redundancy selection)
  if (exists("df_standard") && nrow(df_standard) > 0) {
    top_df_per_comp <- df_standard
  } else if (exists("df_refined") && nrow(df_refined) > 0) {
    top_df_per_comp <- df_refined
  } else {
    top_df_per_comp <- tibble()
  }
}

if (exists("top_df_per_comp") && nrow(top_df_per_comp) > 0) {
  df_per_comp_sub <- combined_df %>%
    filter(Description %in% unique(top_df_per_comp$Description))
  if (nrow(df_per_comp_sub) > 0) {
    ordering_res <- order_dotplot(df_per_comp_sub)
    full_data_for_plot <- df_per_comp_sub %>%
      mutate(
        Comparison = factor(Comparison, levels = ordering_res$col_order),
        Description = factor(Description, levels = ordering_res$row_order)
      )
    dotplot <- ggplot(full_data_for_plot, aes(
      x = Comparison,
      y = Description, 
      color = NES,
      size = -log10(p.adjust)
    )) +
      geom_point(alpha = 1) +
      scale_color_gradientn(
        colours = custom_palette(100),
        name = "NES",
        limits = c(-max(abs(full_data_for_plot$NES), na.rm=TRUE), max(abs(full_data_for_plot$NES), na.rm=TRUE))
      ) +
      scale_size_continuous(
        name = expression(-log[10](italic(P)[adj])),
        range = c(1.5, 5),
        limits = c(0, NA)
      ) +
      scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
      labs(
        title = "Top terms per comparison",
        x = NULL,
        y = NULL
      ) +
      theme_custom(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color="black"),
        panel.grid.major.y = element_blank(),
        legend.position = "right"
      )
    dims_dotplot <- calc_dims(full_data_for_plot)
    output_dotplot <- file.path(subdirs$plots_main, "Dotplot_Enrichment_TopGenes_PerComp.svg")
    ggsave(output_dotplot, plot = dotplot, width = dims_dotplot$w, height = dims_dotplot$h, dpi = 300)
  } else {
    message("No data for top_df_per_comp dotplot.")
  }
} else {
  message("top_df_per_comp is empty; skipping per-comparison dotplot.")
}

## --- Dotplot: Overall top terms (refined) ---
if (exists("top_terms") && length(top_terms) > 0) {
  df_top_sub <- combined_df %>%
    filter(Description %in% top_terms)
  if (nrow(df_top_sub) > 0) {
    ordering_res_top <- order_dotplot(df_top_sub)
    top_terms_all_data <- df_top_sub %>%
      mutate(
        Comparison = factor(Comparison, levels = ordering_res_top$col_order),
        Description = factor(Description, levels = ordering_res_top$row_order)
      )
    dotplot_top <- ggplot(top_terms_all_data, aes(
      x = Comparison,
      y = Description, 
      color = NES,
      size = -log10(p.adjust)
    )) +
      geom_point(alpha = 1) +
      scale_color_gradientn(
        colours = custom_palette(100),
        name = "NES",
        limits = c(-max(abs(top_terms_all_data$NES), na.rm=TRUE), max(abs(top_terms_all_data$NES), na.rm=TRUE))
      ) +
      scale_size_continuous(
        name = expression(-log[10](italic(P)[adj])),
        range = c(1.5, 5),
        limits = c(0, NA)
      ) +
      scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
      labs(
        title = "Top terms (overall)",
        x = NULL,
        y = NULL
      ) +
      theme_custom(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color="black"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2)
      )
    dims_dotplot_top <- calc_dims(top_terms_all_data)
    output_dotplot_top <- file.path(subdirs$plots_main, "Dotplot_Enrichment_TopGenes_Overall.svg")
    ggsave(output_dotplot_top, plot = dotplot_top, width = dims_dotplot_top$w, height = dims_dotplot_top$h, dpi = 300)
  } else {
    message("No data for top_terms dotplot.")
  }
} else {
  message("top_terms is empty; skipping overall top terms dotplot.")
}

# Dotplot: Top 5 up/down terms per comparison (if selected)
if (exists("top5_up_down_df") && !is.null(top5_up_down_df) && nrow(top5_up_down_df) > 0) {
  selected_descriptions <- unique(top5_up_down_df$Description)
  df_split_sub <- combined_df %>%
    filter(Description %in% selected_descriptions)
  ordering_res_split <- order_dotplot(df_split_sub)
  plot_data_split <- df_split_sub %>%
    mutate(
      Comparison = factor(Comparison, levels = ordering_res_split$col_order),
      Description = factor(Description, levels = ordering_res_split$row_order)
    )
  max_abs_nes_split <- max(abs(plot_data_split$NES), na.rm = TRUE)
  dotplot_split <- ggplot(plot_data_split, aes(
    x = Comparison,
    y = Description, 
    color = NES,
    size = -log10(p.adjust)
  )) +
    geom_point(alpha = 1) +
    scale_color_gradientn(
      colours = custom_palette(100),
      limits = c(-max_abs_nes_split, max_abs_nes_split),
      name = "NES"
    ) +
    scale_size_continuous(
      name = expression(-log[10](italic(P)[adj])),
      range = c(1.5, 5),
      limits = c(0, NA)
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
    labs(
      title = "Top 5 up/down terms",
      x = NULL,
      y = NULL
    ) +
    theme_custom(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color="black"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2)
    )
  dims_split <- calc_dims(plot_data_split)
  output_dotplot_split <- file.path(subdirs$plots_main, "Dotplot_Enrichment_Top5_Up_Down_PerComp.svg")
  ggsave(output_dotplot_split, plot = dotplot_split, width = dims_split$w, height = dims_split$h, dpi = 300)
}

# -----------------------------------------------------
# Plot Enrichment for Selected Genes
# -----------------------------------------------------

# gene_list <- c("P55099", "Q6NXX1")
# long_df <- combined_df %>%
#   mutate(core_gene = strsplit(as.character(core_enrichment), "/|;|,|\\s+")) %>%
#   unnest(core_gene)
# filtered_df <- long_df %>%
#   filter(core_gene %in% gene_list) %>%
#   mutate(Description = factor(Description, levels = unique(Description)))
# ggsave(output_plot_file, plot = plot_selected_genes, width = 8, height = 6, dpi = 300)
#
# # Only plot if filtered_df is not empty and NES is not all NA
# if (nrow(filtered_df) > 0 && any(!is.na(filtered_df$NES))) {
#   nes_min <- min(filtered_df$NES, na.rm = TRUE)
#   nes_max <- max(filtered_df$NES, na.rm = TRUE)
#   max_abs_nes <- max(abs(nes_min), abs(nes_max))
#   plot_selected_genes <- ggplot(filtered_df, aes(
#     x = Comparison,
#     y = Description,
#     color = NES,
#     size = -log10(p.adjust),
#     shape = core_gene
#   )) +
#     geom_point(alpha = 1) +
#     scale_color_gradientn(
#       colours = custom_palette(100),
#       limits = c(-max_abs_nes, max_abs_nes),
#       name = "NES"
#     ) +
#     scale_size_continuous(
#       name = expression(-log[10](italic(P)[adj])),
#       range = c(2, 6)
#     ) +
#     labs(
#       title = "Enrichment for selected genes",
#       x = NULL,
#       y = NULL,
#       shape = "Gene"
#     ) +
#     theme_custom(base_size = 10) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   output_plot_file <- file.path(subdirs$plots_main, "Plot_Selected_Genes.svg")
#   ggsave(output_plot_file, plot = plot_selected_genes, width = 8, height = 6, dpi = 300)
# } else {
#   message("No data available for selected genes plot. Skipping plot generation.")
# }

# -----------------------------------------------------
# Gene-Centric Enrichment Summary Plot
# -----------------------------------------------------

gene_enrichment_summary <- long_df %>%

  filter(core_gene %in% gene_list) %>%
  group_by(core_gene, Comparison) %>%
  summarise(
    count = n(),
    mean_NES = mean(NES, na.rm = TRUE),
    max_padj = if (all(is.na(p.adjust))) NA_real_ else min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  )

plot_gene_summary <- ggplot(gene_enrichment_summary, aes(
  x = Comparison,
  y = core_gene,
  size = count,
  fill = mean_NES
)) +
  geom_point(shape = 21, stroke = 0.2, color="black") +
  scale_size_continuous(
    name = "Occurrences",
    range = c(2, 8)
  ) +
  scale_fill_gradientn(
    colours = rev(brewer.pal(n = 7, name = "RdBu")),
    name = "Mean NES",
    limits = c(-max_abs_nes, max_abs_nes)
  ) +
  labs(
    title = "Gene-centric enrichment",
    x = NULL,
    y = "Gene"
  ) +
  theme_custom(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gene_summary_file <- file.path(subdirs$gene_centric, "Gene_Centric_Summary.svg")
ggsave(gene_summary_file, plot = plot_gene_summary, width = 5, height = 4, dpi = 300)

# -----------------------------------------------------
# Extract and Export Core Genes per Comparison
# -----------------------------------------------------

core_gene_sets <- lapply(names(enrichment_list), function(name) {
  df <- enrichment_list[[name]]
  df_long <- df %>%
    dplyr::select(Description, core_enrichment) %>%
    mutate(core_enrichment = str_split(core_enrichment, "/")) %>%
    unnest(core_enrichment) %>%
    mutate(Comparison = name)
  return(df_long)
})

names(core_gene_sets) <- names(enrichment_list)
core_genes_df <- bind_rows(core_gene_sets)
write.csv(core_genes_df,
          file = file.path(subdirs$gene_lists, "Core_Genes_All_Comparisons.csv"),
          row.names = FALSE)
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
# Compute Jaccard Similarity of Core Genes Across Comparisons
# -----------------------------------------------------

core_gene_sets_aggregated <- core_genes_df %>%
  group_by(Comparison) %>%
  summarise(Genes = list(unique(core_enrichment)))
all_genes <- unique(unlist(core_gene_sets_aggregated$Genes))
binary_matrix <- sapply(core_gene_sets_aggregated$Genes, function(gene_list) all_genes %in% gene_list)
rownames(binary_matrix) <- all_genes
colnames(binary_matrix) <- core_gene_sets_aggregated$Comparison
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

# Choose color for numbers based on background
# If similarity > 0.5 (darker blue), use white, else black
number_col_mat <- matrix(ifelse(jaccard_matrix > 0.6, "white", "black"),
                         nrow = nrow(jaccard_matrix), ncol = ncol(jaccard_matrix),
                         dimnames = dimnames(jaccard_matrix))

jaccard_heatmap <- pheatmap(
  jaccard_matrix,
  main = "Jaccard similarity",
  color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
  breaks = seq(0, 1, length.out = 101),
  display_numbers = TRUE,
  number_color = number_col_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 14,
  fontsize_number = 12, 
  border_color = "white",
  cellwidth = 30,
  cellheight = 30
)
jaccard_sim_matrix <- file.path(subdirs$plots_main, "Heatmap_Jaccard_Similarity.svg")
svg(jaccard_sim_matrix, width = max(7, 1 + 0.5 * ncol(jaccard_matrix)), height = max(6, 1 + 0.5 * nrow(jaccard_matrix)), family = "sans")
grid::grid.draw(jaccard_heatmap$gtable)
dev.off()

# -----------------------------------------------------
# Generate Heatmap of Core Enrichment Genes
# -----------------------------------------------------

core_long_df <- bind_rows(
  lapply(names(enrichment_list), function(name) {
    df <- enrichment_list[[name]] %>%
      dplyr::select(GO_ID = ID, Description, NES, core_enrichment) %>%
      mutate(NES = as.numeric(NES),
             Comparison = name) %>%
      filter(!is.na(NES))
    
    # Split core_enrichment safely and expand rowwise
    df <- df %>%
      rowwise() %>%
      mutate(Gene = list(strsplit(core_enrichment, "/")[[1]])) %>%
      tidyr::unnest_longer(Gene) %>%
      ungroup() %>%
      dplyr::select(GO_ID, Description, NES, Gene, Comparison)
    
    return(df)
  })
)

# save enrichment list to excel for debugging
writexl::write_xlsx(enrichment_list, file.path(subdirs$tables, "Enrichment_List_Debug.xlsx"))

# save core_long_df for debugging
write.csv(core_long_df, file.path(subdirs$tables, "Core_Enrichment_LongFormat.csv"), row.names = FALSE)

heatmap_df <- core_long_df %>%
  group_by(Gene, Comparison) %>%
  summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  mutate(NES = ifelse(is.nan(NES) | is.infinite(NES), NA, NES)) %>%
  pivot_wider(names_from = Comparison, values_from = NES)
heatmap_matrix <- heatmap_df %>%
  column_to_rownames("Gene") %>%
  as.matrix()

heatmap_matrix[is.na(heatmap_matrix)] <- 0
nes_min <- min(heatmap_matrix, na.rm = TRUE)
nes_max <- max(heatmap_matrix, na.rm = TRUE)
max_abs_nes <- max(abs(nes_min), abs(nes_max))
breaks <- seq(-max_abs_nes, max_abs_nes, length.out = 101)
# Use the custom_palette defined earlier for color
heatmap_plot <- pheatmap(
  heatmap_matrix,
  color = custom_palette(100),
  breaks = breaks,
  main = "Core enrichment heatmap",
  fontsize = 12,
  fontsize_row = 1, # effectively hides row labels (protein names)
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA
)
output_file <- file.path(subdirs$core_enrichment_plots, "Heatmap_Overall_Core_Enrichment.svg")
svg(output_file, width = 8, height = 8, family = "sans")
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# Identify top 25 up/down genes by NES for supplementary tables/plots
# -------------------
# IMPORTANT NOTE FOR USERS:
# -------------------
# The top/bottom 25 genes in the heatmap ("TopBottom_GeneList.xlsx" and "Heatmap_TopBottom_Enrichment.svg")
# are selected based on the maximum (for up) or minimum (for down) NES (Normalized Enrichment Score)
# across ALL GO terms and ALL comparisons, i.e., based on GSEA core enrichment results.
#
# The "supplementary full_differential_data" (sheet "Full_Differential_Data" in Supplementary_Data.xlsx)
# contains ALL proteins from the differential expression analysis, ranked by log2FC and p-value per comparison.
#
# Therefore, the top/bottom 25 proteins in the heatmap are NOT necessarily the same as the top/bottom proteins
# in the full differential data, because:
#   - The heatmap is based on NES from GSEA/GO enrichment, not log2FC from differential analysis.
#   - A protein may have a high NES in enrichment but not be the most up/downregulated by log2FC, and vice versa.
#   - The selection is global across all comparisons, not per comparison.
#
# This is expected behavior and reflects the different ranking criteria.
message("[INFO] Top/bottom 25 genes for heatmap are selected by NES (enrichment), not by log2FC (differential expression). See comments in code for details.")


# -----------------------------------------------------
# Select and Visualize Top 25 Up/Down Genes by NES Across All Comparisons
# -----------------------------------------------------
# This section identifies the top 25 upregulated and top 25 downregulated genes based on the maximum (for up) or minimum (for down)
# NES (Normalized Enrichment Score) observed for each gene across all comparisons. The data used is the heatmap_df, which contains
# NES values for each gene (rows) and comparison (columns). The resulting top/bottom gene list is exported, and a heatmap is generated.

# --- Data Preparation ---
# Convert heatmap_df (wide: genes x comparisons) to long format for easier aggregation.
long_heatmap_df <- heatmap_df %>%
  pivot_longer(-Gene, names_to = "Comparison", values_to = "NES")

# --- Top/Bottom Gene Selection ---
# For each gene, find the maximum NES (upregulation) and minimum NES (downregulation) across all comparisons.
gene_max_nes <- long_heatmap_df %>%
  group_by(Gene) %>%
  summarise(max_nes = max(NES, na.rm = TRUE), .groups = "drop")
gene_min_nes <- long_heatmap_df %>%
  group_by(Gene) %>%
  summarise(min_nes = min(NES, na.rm = TRUE), .groups = "drop")

# Select top 25 upregulated and top 25 downregulated genes.
top25_up <- gene_max_nes %>%
  arrange(desc(max_nes)) %>%
  dplyr::slice_head(n = 25) %>%
  mutate(Direction = "Up")
top25_down <- gene_min_nes %>%
  arrange(min_nes) %>%
  dplyr::slice_head(n = 25) %>%
  mutate(Direction = "Down")

# Combine into a single data frame for export and plotting.
top_bottom_genes <- bind_rows(top25_up, top25_down)

# --- Add GO Terms for Each Gene ---
# For each gene, collect all unique GO term descriptions from core_long_df.
gene_go <- core_long_df %>%
  dplyr::select(Gene, Description) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(GO_Terms = paste(unique(Description), collapse = "; "), .groups = "drop")
top_bottom_genes <- top_bottom_genes %>%
  left_join(gene_go, by = "Gene")

# --- Export Top/Bottom Gene List ---
# Save the top/bottom gene list (with NES and GO terms) to Excel for supplementary tables.
# Combine max_nes and min_nes into a single NES column (show as "max: X; min: Y")
top_bottom_genes <- top_bottom_genes %>%
  mutate(NES = paste0("max: ", signif(max_nes, 4), "; min: ", signif(min_nes, 4))) %>%
  dplyr::select(-max_nes, -min_nes) # Remove old columns

write_xlsx(
  top_bottom_genes,
  file.path(subdirs$tables, "TopBottom_GeneList.xlsx")
)

# --- Generate Heatmap for Top 25 Up/Down Genes ---
# Subset the heatmap matrix to only the selected top/bottom genes.
top_bottom_matrix <- heatmap_matrix[rownames(heatmap_matrix) %in% top_bottom_genes$Gene, , drop = FALSE]

# Order genes in the heatmap by their maximum absolute NES (for visual clarity).
gene_order <- heatmap_df %>%
  filter(Gene %in% rownames(top_bottom_matrix)) %>%
  pivot_longer(-Gene, names_to = "Comparison", values_to = "NES") %>%
  group_by(Gene) %>%
  summarize(order_metric = max(abs(NES), na.rm = TRUE)) %>%
  arrange(desc(order_metric)) %>%
  pull(Gene)

top_bottom_matrix <- top_bottom_matrix[gene_order, , drop = FALSE]

# --- Map UniProt IDs to Gene Names for Heatmap Row Labels ---
# Use the UniProt mapping file to convert UniProt IDs to gene names for better readability in the heatmap.
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  filter(V2 == "Gene_Name") %>%
  dplyr::select(UniprotID = V1, Gene_Name = V3)
mapped_names <- left_join(
  tibble(Gene = rownames(top_bottom_matrix)),
  uniprot_df,
  by = c("Gene" = "UniprotID")
)

# Replace rownames with gene names where available, otherwise keep UniProt ID.
rownames(top_bottom_matrix) <- ifelse(
  is.na(mapped_names$Gene_Name),
  mapped_names$Gene,
  mapped_names$Gene_Name
)

# --- Plotting ---
# Create a heatmap of NES values for the top/bottom genes across all comparisons.
# Color scale: blue = downregulated, red = upregulated, white = neutral.
nes_min <- min(top_bottom_matrix, na.rm = TRUE)
nes_max <- max(top_bottom_matrix, na.rm = TRUE)
max_abs_nes <- max(abs(nes_min), abs(nes_max))
breaks <- seq(-max_abs_nes, max_abs_nes, length.out = 101)

# --- Add more info text on the plot ---
library(grid)
top_bottom_plot <- pheatmap(
  top_bottom_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  breaks = breaks,
  main = "Core enrichment overview",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 8,
  fontsize_number = 6,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  border_color = "white",
  border_width = 0.5,
  cellwidth = 10,
  cellheight = 7,
  angle_col = 45,
  treeheight_row = 15,
  treeheight_col = 15,
  legend = TRUE,
  legend_breaks = c(round(-max_abs_nes,1), 0, round(max_abs_nes,1)),
  show_rownames = TRUE,
  show_colnames = TRUE
)

# Add info text: n genes, NES scale, etc.
info_text <- paste0(
  "n genes: ", nrow(top_bottom_matrix),
  " | n comparisons: ", ncol(top_bottom_matrix),
  " | NES scale: [", round(-max_abs_nes, 2), ", ", round(max_abs_nes, 2), "]\n",
  "Red = upregulated, Blue = downregulated, NES = Normalized Enrichment Score"
)

svg(file.path(subdirs$plots_main, "Heatmap_TopBottom_Enrichment.svg"), width = 4, height = 8, family="sans")
grid::grid.draw(top_bottom_plot$gtable)
grid::grid.text(
  info_text,
  x = 0.5, y = 0.98, just = "center", gp = grid::gpar(fontsize = 8, col = "#222222"))
dev.off()

# Save the heatmap as SVG for publication-quality output.
svg(file.path(subdirs$plots_main, "Heatmap_TopBottom_Enrichment.svg"), width = 4, height = 8, family="sans")
grid::grid.draw(top_bottom_plot$gtable)
dev.off()

# -----------------------------------------------------
# Generate Per-Comparison Heatmaps and Excel Files
# -----------------------------------------------------

uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

uniprot_subset <- uniprot_df %>%
  filter(V2 %in% c("Gene_Name", "UniProtKB-ID")) %>%
  pivot_wider(names_from = V2, values_from = V3, values_fn = list) %>%
  unnest(cols = everything()) %>%
  distinct(V1, .keep_all = TRUE) %>%
  dplyr::rename(UniprotID = V1)
gene_synonyms <- uniprot_df %>%
  filter(V2 == "Gene_Synonym") %>%
  distinct(V1, .keep_all = TRUE) %>%
  dplyr::rename(UniprotID = V1, Gene_Synonym = V3)
gene_descriptions <- core_long_df %>%
  dplyr::select(Gene, Description) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(Description = paste(unique(Description), collapse = "; "), .groups = "drop")
log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

log2fc_df_list <- lapply(log2fc_files, function(f) {
  df <- read_csv(f, col_types = cols())
  comp_name <- tools::file_path_sans_ext(basename(f))
  # Ensure Comparison column exists
  df <- df %>%
    mutate(Comparison = comp_name) %>%
    dplyr::rename(
      padj = adj.P.Val,
      log2fc = logFC
    ) %>%
    as.data.frame()
  return(df)
})

log2fc_long_df <- bind_rows(log2fc_df_list)
unique_comparisons <- unique(log2fc_long_df$Comparison)
for (comp in unique_comparisons) {
  full_comp_df <- log2fc_long_df %>%
    filter(Comparison == comp, !is.na(log2fc), !is.na(padj))
  # Volcano plot for each comparison
  volcano_df <- full_comp_df %>%
    filter(!str_detect(gene_symbol, "_")) %>%
    left_join(uniprot_subset %>% dplyr::select(UniprotID, Gene_Name),
              by = c("gene_symbol" = "UniprotID")) %>%
    mutate(Significance = case_when(
      padj < 0.05 & log2fc > 0 ~ "up",
      padj < 0.05 & log2fc < 0 ~ "down",
      TRUE ~ "n.s."
    ))
  top5_up <- volcano_df %>%
    filter(Significance == "up") %>%
    arrange(desc(log2fc)) %>%
    dplyr::slice_head(n = 5)
  top5_down <- volcano_df %>%
    filter(Significance == "down") %>%
    arrange(log2fc) %>%
    dplyr::slice_head(n = 5)
  top_genes_to_label <- bind_rows(top5_up, top5_down)
  max_x_val <- max(abs(volcano_df$log2fc), 0, na.rm = TRUE)
  x_limit <- ceiling(max_x_val) + 0.5
  max_y_val <- max(abs(-log10(volcano_df$padj)), 0, na.rm = TRUE)
  y_limit <- ceiling(max_y_val) + 1
  volcano_plot <- ggplot(volcano_df, aes(x = log2fc, y = -log10(padj), label = Gene_Name)) +
    geom_point(aes(color = Significance), shape = 16, alpha = 0.6, size = 6) + # larger dots, 0.6 alpha
    ggrepel::geom_text_repel(
      data = top_genes_to_label,
      aes(label = Gene_Name),
      size = 6, # larger font
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.7
    ) +
    scale_color_manual(
      values = c("up" = "#f36e21", "down" = "#465b65", "n.s." = "#BABABA"),
      breaks = c("up", "down", "n.s.")
    ) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    labs(
      title = paste("Volcano plot:", comp),
      x = expression(log[2] ~ Fold ~ Change),
      y = expression(-log[10](italic(P)[adj]))
    ) +
    theme_custom(base_size = 16) + # larger base font
    theme(
      legend.position = "none"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", linewidth = 0.3)
  volcano_file <- file.path(subdirs$volcanoes, paste0("Volcano_", comp, ".svg"))
  message(paste("Saving volcano plot to:", volcano_file))
  ggsave(
    filename = volcano_file,
    plot = volcano_plot,
    width = 5,
    height = 5
  )
  comp_df <- full_comp_df %>% 
    filter(padj < 0.05) %>%
    ungroup()
  if (nrow(comp_df) == 0) {
    message(paste("  No significant genes found for", comp, "- skipping heatmap and Excel."))
    next
  }
  comp_df <- comp_df %>%
    mutate(Direction = case_when(
      log2fc > 0 ~ "Up",
      log2fc < 0 ~ "Down",
      TRUE ~ "Neutral"
    ))
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
  mapped_top_bottom <- mapped_top_bottom %>%
    group_by(Direction) %>%
    arrange(Direction, desc(abs_log2fc)) %>%
    mutate(Rank = row_number()) %>%
    ungroup()
  heatmap_matrix <- mapped_top_bottom %>%
    dplyr::select(Gene_Label, log2fc) %>%
    tibble::column_to_rownames("Gene_Label") %>%
    as.matrix()
  if (nrow(heatmap_matrix) >= 2) {
    max_abs <- max(abs(heatmap_matrix), na.rm = TRUE)
    breaks <- seq(-max_abs, max_abs, length.out = 101)
    heatmap_plot <- pheatmap(
      heatmap_matrix,
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
      breaks = breaks,
      main = NA,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize_row = 8,
      show_colnames = FALSE,
      border_color = "white",
      silent = TRUE,
      legend = FALSE
    )
    heatmap_file <- file.path(subdirs$plots_main, paste0("Heatmap_Log2FC_", comp, ".svg"))
    message(paste("  Saving heatmap to:", heatmap_file))
    svg(heatmap_file, width = 2.5, height = 5, family="sans")
    grid::grid.draw(heatmap_plot$gtable)
    dev.off()
  }
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
comparison_files <- list.files(
  path = subdirs$sig_proteins,
  pattern = "SigProteins_(Upregulated|Downregulated)_.*\\.xlsx$",
  full.names = TRUE
)
if (length(comparison_files) > 0) {
  # Extract comparison names robustly from Upregulated or Downregulated files
  comparisons <- unique(gsub("SigProteins_(Upregulated|Downregulated)_|\\.xlsx", "", 
                             basename(comparison_files)))
  message("Found ", length(comparisons), " comparisons to process: ", paste(comparisons, collapse = ", "))

  for (comp in comparisons) {
    message(paste("Processing GO enrichment for comparison:", comp))
    
    # --- Process Upregulated Genes ---
    file_up <- file.path(subdirs$sig_proteins, paste0("SigProteins_Upregulated_", comp, ".xlsx"))
    
    if (file.exists(file_up)) {
      df_up <- readxl::read_xlsx(file_up)
      up_genes <- df_up %>% filter(!is.na(Uniprot_Accession)) %>% pull(Uniprot_Accession) %>% unique()
      
      # Use UNIPROT IDs directly
      if (length(up_genes) > 0) {
        ego_up <- tryCatch(
          enrichGO(
            gene = up_genes,
            OrgDb = org.Mm.eg.db,
            keyType = "UNIPROT",
            ont = ont,
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            minGSSize = 10
          ),
          error = function(e) NULL
        )
        if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
          out_csv <- file.path(subdirs$go_enrichment, paste0("GO_", comp, "_Upregulated_", ont, ".csv"))
          write.csv(as.data.frame(ego_up), out_csv, row.names = FALSE)
          ego_df <- as.data.frame(ego_up) %>%
            dplyr::mutate(Count = as.numeric(Count)) %>%
            dplyr::filter(Count >= 10) %>%
            mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
            arrange(desc(GeneRatio)) %>%
            dplyr::slice_head(n = 15) %>%
            mutate(Description = stringr::str_wrap(Description, width = 50)) %>%
            mutate(Description = factor(Description, levels = unique(Description)))
          if (nrow(ego_df) > 0) {
            lollipop_palette <- colorRampPalette(c("#faa51a", "#4c87c6"))
            p_lollipop <- ggplot(ego_df, aes(x = GeneRatio, y = Description)) +
              geom_segment(aes(x = 0, xend = GeneRatio, y = Description, yend = Description), color = "#bdbdbd", linewidth = 1) +
              geom_point(aes(size = Count, color = p.adjust), alpha = 0.9) +
              scale_color_gradientn(colours = lollipop_palette(100), name = expression(italic(P)[adj])) +
              scale_size_continuous(range = c(2, 7), name = "Gene Count") +
              scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)), name = "Gene Ratio") +
              labs(
                title = paste0("Upregulated GO - ", comp),
                y = "GO Term"
              ) +
              theme_custom(base_size = 10) +
              theme(
                axis.text.y = element_text(size = 9),
                axis.text.x = element_text(size = 9),
                legend.position = "right"
              )
            out_plot <- file.path(subdirs$go_enrichment, paste0("GO_Lollipop_", comp, "_Upregulated_", ont, ".svg"))
            ggsave(out_plot, plot = p_lollipop, width = 7, height = 5)
          }
        }
      }
    } else {
      message("  No Upregulated file found for: ", comp)
    }
    
    # --- Process Downregulated Genes ---
    file_down <- file.path(subdirs$sig_proteins, paste0("SigProteins_Downregulated_", comp, ".xlsx"))
    
    if (file.exists(file_down)) {
      df_down <- readxl::read_xlsx(file_down)
      down_genes <- df_down %>% filter(!is.na(Uniprot_Accession)) %>% pull(Uniprot_Accession) %>% unique()
      
      # Use UNIPROT IDs directly
      if (length(down_genes) > 0) {
        ego_down <- tryCatch(
          enrichGO(
            gene = down_genes,
            OrgDb = org.Mm.eg.db,
            keyType = "UNIPROT",
            ont = ont,
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            minGSSize = 10
          ),
          error = function(e) NULL
        )
        if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
          out_csv <- file.path(subdirs$go_enrichment, paste0("GO_", comp, "_Downregulated_", ont, ".csv"))
          write.csv(as.data.frame(ego_down), out_csv, row.names = FALSE)
          ego_df <- as.data.frame(ego_down) %>%
            dplyr::mutate(Count = as.numeric(Count)) %>%
            dplyr::filter(Count >= 10) %>%
            mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
            arrange(desc(GeneRatio)) %>%
            dplyr::slice_head(n = 15) %>%
            mutate(Description = stringr::str_wrap(Description, width = 50)) %>%
            mutate(Description = factor(Description, levels = unique(Description)))
          if (nrow(ego_df) > 0) {
            lollipop_palette <- colorRampPalette(c("#faa51a", "#4c87c6"))
            p_lollipop <- ggplot(ego_df, aes(x = GeneRatio, y = Description)) +
              geom_segment(aes(x = 0, xend = GeneRatio, y = Description, yend = Description), color = "#bdbdbd", linewidth = 1) +
              geom_point(aes(size = Count, color = p.adjust), alpha = 0.9) +
              scale_color_gradientn(colours = lollipop_palette(100), name = expression(italic(P)[adj])) +
              scale_size_continuous(range = c(2, 7), name = "Gene Count") +
              scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)), name = "Gene Ratio") +
              labs(
                title = paste0("Downregulated GO - ", comp),
                y = "GO Term"
              ) +
              theme_custom(base_size = 10) +
              theme(
                axis.text.y = element_text(size = 9),
                axis.text.x = element_text(size = 9),
                legend.position = "right"
              )
            out_plot <- file.path(subdirs$go_enrichment, paste0("GO_Lollipop_", comp, "_Downregulated_", ont, ".svg"))
            ggsave(out_plot, plot = p_lollipop, width = 7, height = 5)
          }
        }
      }
    } else {
      message("  No Downregulated file found for: ", comp)
    }
  }
}

#AIres <- interpret(ego_up)

# -----------------------------------------------------
# Generate Individual Core Enrichment Heatmaps
# -----------------------------------------------------


# --- Robust log2fc import and renaming ---
log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)
names(log2fc_files) <- basename(log2fc_files) %>%
  str_remove(".csv") %>%
  str_trim()
log2fc_list <- lapply(names(log2fc_files), function(comp) {
  df <- read.csv(log2fc_files[comp])
  df$Comparison <- comp
  # Robust renaming: only rename columns that exist
  col_map <- c(log2fc = "logFC", pvalue = "P.Value", padj = "adj.P.Val")
  for (new_col in names(col_map)) {
    old_col <- col_map[[new_col]]
    if (old_col %in% colnames(df)) {
      names(df)[names(df) == old_col] <- new_col
    }
  }
  df
})
log2fc_df <- bind_rows(log2fc_list)

# Warn if required columns are missing
required_cols <- c("gene_symbol", "Comparison", "log2fc")
missing_cols <- setdiff(required_cols, colnames(log2fc_df))
if (length(missing_cols) > 0) {
  warning(paste("log2fc_df is missing columns:", paste(missing_cols, collapse=", ")))
}


# --- Robust join for core_long_df heatmaps ---
# Save per-term core enrichment heatmaps in the correct subdirectory (core_enrichment_plots)
core_long_df %>%
  group_by(Description) %>%
  group_split() %>%
  walk(function(df_term) {
    description <- unique(df_term$Description)
    # Only join if both columns exist
    join_by <- c("Gene" = "gene_symbol", "Comparison" = "Comparison")
    join_ok <- all(c("Gene", "Comparison") %in% colnames(df_term)) &&
               all(c("gene_symbol", "Comparison") %in% colnames(log2fc_df))
    if (!join_ok) {
      warning(paste0("Skipping heatmap for ", description, ": join columns missing."))
      return()
    }
    df_term_log2fc <- df_term %>%
      left_join(log2fc_df, by = join_by)
    if (!"log2fc" %in% colnames(df_term_log2fc)) {
      warning(paste0("Skipping heatmap for ", description, ": log2fc column missing after join."))
      return()
    }
    matrix_df <- df_term_log2fc %>%
      dplyr::select(Gene, Comparison, log2fc) %>%
      group_by(Gene, Comparison) %>%
      summarize(log2fc = mean(log2fc, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Comparison, values_from = log2fc) %>%
      column_to_rownames("Gene")
    colnames(matrix_df) <- str_trim(colnames(matrix_df))
    heatmap_matrix <- as.matrix(matrix_df)
    heatmap_matrix[is.na(heatmap_matrix)] <- 0
    safe_desc <- substr(make.names(description), 1, 30)
    filename_svg <- file.path(subdirs$core_enrichment_plots, paste0("Heatmap_Core_", safe_desc, ".svg"))
    filename_png <- file.path(subdirs$core_enrichment_plots, paste0("Heatmap_Core_", safe_desc, ".png"))
    cluster_rows_option <- nrow(heatmap_matrix) > 1
    cluster_cols_option <- ncol(heatmap_matrix) > 1
    max_val <- suppressWarnings(max(abs(heatmap_matrix), na.rm = TRUE))
    # If max_val is NA or 0, set breaks to c(-1,0,1) and skip plotting if all values are zero
    if (is.na(max_val) || max_val == 0) {
      breaks <- c(-1, 0, 1)
      # Only plot if there is at least one nonzero value
      if (all(heatmap_matrix == 0)) {
        warning(paste0("Skipping heatmap for ", description, ": all values are zero."))
        return()
      }
    } else {
      # Ensure breaks are unique
      breaks <- unique(seq(-max_val, max_val, length.out = 101))
      if (length(breaks) < 2) {
        breaks <- c(-1, 0, 1)
      }
    }
    plot_height <- max(4, nrow(heatmap_matrix) * 0.15 + 1)
    # Save SVG
    svg(filename_svg, width = 5, height = plot_height, family = "sans")
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
      breaks = breaks,
      main = substr(description, 1, 40),
      fontsize = 8,
      fontsize_row = 7,
      fontsize_col = 8,
      cluster_rows = cluster_rows_option,
      cluster_cols = cluster_cols_option,
      border_color = "white",
      angle_col = 45,
      cellwidth = 12,
      cellheight = 10,
      legend = TRUE
    )
    dev.off()
    # Save PNG as well for convenience
    png(filename_png, width = 5 * 100, height = plot_height * 100, res = 100, units = "px", family = "sans")
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
      breaks = breaks,
      main = substr(description, 1, 40),
      fontsize = 8,
      fontsize_row = 7,
      fontsize_col = 8,
      cluster_rows = cluster_rows_option,
      cluster_cols = cluster_cols_option,
      border_color = "white",
      angle_col = 45,
      cellwidth = 12,
      cellheight = 10,
      legend = TRUE
    )
    dev.off()
  })

# -----------------------------------------------------
# Gene-Centric Log2 Fold-Change (log2FC) Analysis
# -----------------------------------------------------

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

log2fc_long <- bind_rows(log2fc_list) %>%
  dplyr::rename(
    log2fc = logFC,
    pvalue = P.Value,
    padj = adj.P.Val
  )

log2fc_filtered <- log2fc_long %>%
  filter(gene_symbol %in% gene_list)

gene_fc_summary <- log2fc_filtered %>%
  group_by(gene_symbol, Comparison) %>%
  summarise(
    count = n(),
    mean_log2fc = mean(log2fc, na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  )

max_abs_fc <- max(abs(gene_fc_summary$mean_log2fc), na.rm = TRUE)
plot_gene_fc <- ggplot(gene_fc_summary, aes(
  x = Comparison,
  y = gene_symbol,
  size = count,
  fill = mean_log2fc
)) +
  geom_point(shape = 21, stroke = 0.2, color="black") +
  scale_size_continuous(guide = "none", range = c(2, 6)) +  
  scale_fill_gradientn(
    colours = custom_palette(100),
    name = "log2 FC",
    limits = c(-max_abs_fc, max_abs_fc)
  ) +
  labs(
    title = "Gene-centric log2 FC",
    x = NULL,
    y = "Gene"
  ) +
  theme_custom(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
fc_plot_file <- file.path(subdirs$gene_centric, paste0("GeneCentric_Log2FC_Bubble_", condition, ".svg"))
ggsave(fc_plot_file, plot = plot_gene_fc, width = 4, height = 5, dpi = 300)
summary_file <- file.path(subdirs$gene_centric, paste0("GeneCentric_Log2FC_Summary.xlsx"))
writexl::write_xlsx(gene_fc_summary, path = summary_file)

# -----------------------------------------------------
# Generate Supplementary Tables (Merged)
# -----------------------------------------------------

# Helper function to map slash-separated UniProt IDs to Gene Names only
map_core_enrichment_genes <- function(core_str, mapping_df_genes) {
  if (is.na(core_str) || core_str == "") return(NA)
  
  ids <- unlist(strsplit(core_str, "/"))
  
  # Map Gene Names
  genes <- mapping_df_genes$Gene_Name[match(ids, mapping_df_genes$UniprotID)]
  # Fallback to ID if no gene name found
  genes[is.na(genes)] <- ids[is.na(genes)]
  genes_str <- paste(genes, collapse = "/")
  
  return(genes_str)
}

# Apply mapping to combined_df
# We use existing uniprot_subset for Gene Names
mapped_genes <- sapply(combined_df$core_enrichment, function(x) {
  map_core_enrichment_genes(x, uniprot_subset)
})

combined_df$core_enrichment_gene <- mapped_genes

supp_enrichment <- combined_df %>%
  dplyr::select(Comparison, ID, Description, setSize, 
                pvalue, p.adjust, qvalue, NES, rank, leading_edge,
                core_enrichment, core_enrichment_gene) %>%
  mutate(
    leading_edge = vapply(leading_edge, function(x) {
      if (is.null(x) || all(is.na(x))) return(NA_character_)
      paste(x, collapse = "/")
    }, character(1)),
    core_enrichment = as.character(core_enrichment),
    core_enrichment_gene = as.character(core_enrichment_gene)
  ) %>%
  mutate(across(c(core_enrichment, core_enrichment_gene, leading_edge), 
                function(x) {
                  x <- as.character(x)
                  ifelse(nchar(x) > 32000, paste0(substr(x, 1, 32000), "...[TRUNCATED]"), x)
                }))

supp_sig_proteins <- log2fc_long %>%
  filter(padj < 0.05) %>%
  left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
  left_join(gene_synonyms, by = c("gene_symbol" = "UniprotID")) %>% 
  left_join(gene_descriptions, by = c("gene_symbol" = "Gene")) %>%
  mutate(Direction = ifelse(log2fc > 0, "Up", "Down")) %>%
  dplyr::select(Comparison, UniprotID = gene_symbol, Gene_Name, Gene_Synonym, Description,
                log2fc, pvalue, padj, Direction) %>%
  arrange(Comparison, padj)

supp_gene_centric <- gene_fc_summary %>%
  left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
  dplyr::rename(Gene_Name = Gene_Name)

supp_nes_matrix <- lookup_df %>%
  dplyr::select(Description, Comparison, NES) %>%
  pivot_wider(names_from = Comparison, values_from = NES)

supp_binary_matrix <- as.data.frame(binary_matrix) %>%
  tibble::rownames_to_column(var = "Gene")

# --- New Sheets Added Below ---

# 1. Jaccard Similarity between comparisons (calculated earlier for heatmap)
supp_jaccard <- as.data.frame(jaccard_matrix) %>%
  tibble::rownames_to_column(var = "Comparison")

# 2. Redundancy Log Detailed
# We augment the simple replacement log with NES and p-values to show WHY terms were swapped
if (exists("df_replacements") && nrow(df_replacements) > 0 && "Dropped_Standard_Term" %in% names(df_replacements)) {
  
  # Prepare lookup tables from combined_df
  lookup_metrics <- combined_df %>% 
    dplyr::select(Comparison, Description, NES, p.adjust)
    
  supp_redundancy <- df_replacements %>%
    # Add metrics for Dropped Terms
    left_join(
      lookup_metrics, 
      by = c("Comparison", "Dropped_Standard_Term" = "Description")
    ) %>%
    dplyr::rename(
      Dropped_NES = NES, 
      Dropped_Padj = p.adjust
    ) %>%
    # Add metrics for Added Terms
    left_join(
      lookup_metrics,
      by = c("Comparison", "Added_Refined_Term" = "Description")
    ) %>%
    dplyr::rename(
      Added_NES = NES,
      Added_Padj = p.adjust
    ) %>%
    # Reorganize columns for clarity
    dplyr::select(
      Comparison, Direction,
      Dropped_Standard_Term, Dropped_NES, Dropped_Padj,
      Added_Refined_Term, Added_NES, Added_Padj
    )

} else {
  supp_redundancy <- data.frame(Info = "Redundancy analysis not run or no data available.")
}

# 3. Full Differential Data (All proteins, not just significant)
supp_full_diff <- log2fc_long %>%
  left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
  dplyr::select(Comparison, UniprotID = gene_symbol, Gene_Name, log2fc, pvalue, padj) %>%
  arrange(Comparison, padj)

supp_list <- list(
  "GO_Enrichment_Results" = supp_enrichment,
  "Significant_Proteins" = supp_sig_proteins,
  "Full_Differential_Data" = supp_full_diff, 
  "Gene_Centric_Log2FC" = supp_gene_centric,
  "NES_Matrix_TopTerms" = supp_nes_matrix,
  "Core_Genes_Binary_Matrix" = supp_binary_matrix,
  "Jaccard_Similarity" = supp_jaccard,       
  "Redundancy_Log" = supp_redundancy         
)

# Add sheets for heatmap_df and heatmap_matrix (top/bottom 25 proteins)

if (exists("top_bottom_genes")) {
  supp_list[["TopBottom_Proteins"]] <- top_bottom_genes
}
if (exists("heatmap_df")) {
  supp_list[["Heatmap_TopBottom25_Protein_Data"]] <- heatmap_df
}
if (exists("heatmap_matrix")) {
  supp_list[["Heatmap_TopBottom25_Protein_Matrix"]] <- tibble::rownames_to_column(as.data.frame(heatmap_matrix), var = "Gene")
}
if (exists("core_long_df")) {
  supp_list[["Core_Enrichment_Gene_LongFormat"]] <- core_long_df
}
if (exists("log2fc_long")) {
  supp_list[["Log2FC_LongFormat"]] <- log2fc_long
}
if (exists("top_bottom_matrix")) {
  supp_list[["TopBottom_Log2FC_Matrix"]] <- tibble::rownames_to_column(as.data.frame(top_bottom_matrix), var = "Gene")
}

supp_output_path <- file.path(subdirs$tables, "Supplementary_Data.xlsx")
writexl::write_xlsx(supp_list, path = supp_output_path)

# perform simplifyGOFromMultipleLists()
# -----------------------------------------------------
# Run simplifyGOFromMultipleLists() for each comparison (up/downregulated proteins)
# -----------------------------------------------------

# Prepare GO enrichment results as a list of data frames for each comparison (up/downregulated proteins)
log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

go_results_list <- list()
for (f in log2fc_files) {
  comp_name <- tools::file_path_sans_ext(basename(f))
  df <- readr::read_csv(f, show_col_types = FALSE)
  # Upregulated
  up_genes <- df %>% filter(logFC > 0) %>% pull(gene_symbol) %>% unique()
  if (length(up_genes) > 0) {
    ego_up <- tryCatch(
      enrichGO(
        gene = up_genes,
        OrgDb = org.Mm.eg.db,
        keyType = "UNIPROT",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      ),
      error = function(e) NULL
    )
    if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
      ego_up_df <- as.data.frame(ego_up)
      # Only filter after collecting, and keep at least a few terms if possible
      if (nrow(ego_up_df) > 0) {
        filtered <- ego_up_df %>% dplyr::filter(p.adjust < 0.05)
        if (nrow(filtered) < 3 && nrow(ego_up_df) >= 3) {
          # If too few significant, use top 3 by p.adjust
          filtered <- ego_up_df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = 3)
        }
        if (nrow(filtered) > 0) {
          go_results_list[[paste0(comp_name, "_Up")]] <- filtered
        }
      }
    }
  }
  # Downregulated
  down_genes <- df %>% filter(logFC < 0) %>% pull(gene_symbol) %>% unique()
  if (length(down_genes) > 0) {
    ego_down <- tryCatch(
      enrichGO(
        gene = down_genes,
        OrgDb = org.Mm.eg.db,
        keyType = "UNIPROT",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      ),
      error = function(e) NULL
    )
    if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
      ego_down_df <- as.data.frame(ego_down)
      if (nrow(ego_down_df) > 0) {
        filtered <- ego_down_df %>% dplyr::filter(p.adjust < 0.05)
        if (nrow(filtered) < 3 && nrow(ego_down_df) >= 3) {
          filtered <- ego_down_df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = 3)
        }
        if (nrow(filtered) > 0) {
          go_results_list[[paste0(comp_name, "_Down")]] <- filtered
        }
      }
    }
  }
}

# Remove empty or too small lists
go_results_list <- go_results_list[sapply(go_results_list, nrow) > 0]

# --- simplifyGOFromMultipleLists: one plot for all comparisons ---
go_results_list <- go_results_list[sapply(go_results_list, nrow) > 0]
if (length(go_results_list) >= 2) {
  simplifygo_plot_file <- file.path(subdirs$plots_main, "simplifyGO_MultipleLists.svg")
  # Ensure any open graphics device is closed before opening a new one
  while(dev.cur() > 1) dev.off()
  svg(simplifygo_plot_file, width = 10, height = 8, family = "sans")
  simplifyEnrichment::simplifyGOFromMultipleLists(
    go_results_list,
    go_id_column = "ID",
    padj_column = "p.adjust",
    padj_cutoff = 0.05,
    ont = ont,
    db = "org.Mm.eg.db",
    measure = "Sim_XGraSM_2013",
    method = "binary_cut",
    verbose = TRUE,
    column_title = paste0("simplifyGO: All comparisons (", ont, ")"),
    heatmap_param = list(
      col = c("white", "#faa51a"),
      breaks = c(1, 0.0001),
      name = "padj"
    )
  )
  dev.off()
  # Save cluster assignments
  all_go <- purrr::imap_dfr(go_results_list, ~dplyr::mutate(.x, Source = .y))
  go_ids <- unique(all_go$ID)
  mat <- simplifyEnrichment::GO_similarity(go_ids, ont = "BP", db = "org.Mm.eg.db")
  clust <- simplifyEnrichment::cluster_terms(mat)
  if (!is.null(clust) && length(clust) > 0) {
    clust_df <- data.frame(
      ID = rownames(mat),
      Cluster = as.integer(clust),
      stringsAsFactors = FALSE
    )
    go_annot <- all_go %>%
      dplyr::select(ID, Source) %>%
      distinct() %>%
      group_by(ID) %>%
      summarise(Source = paste(Source, collapse = "; "), .groups = "drop")
    desc_map <- all_go %>% dplyr::select(ID, Description) %>% distinct()
    out_df <- clust_df %>%
      left_join(go_annot, by = "ID") %>%
      left_join(desc_map, by = "ID")
    simplifygo_table_file <- file.path(subdirs$tables, "simplifyGO_MultipleLists_Results.xlsx")
    writexl::write_xlsx(out_df, simplifygo_table_file)
  } else {
    message("No clusters were found by simplifyEnrichment::cluster_terms or cluster size mismatch; skipping cluster assignment export.")
  }
}

message("Analysis complete. Files saved to: ", main_output_dir)

###############################################################
# Visualize GO term overlap across comparisons as a network graph
# Nodes = GO terms (Description), edges = shared genes or Jaccard similarity
###############################################################


# Network plot: check for nodes and edges before plotting
if (!require("igraph")) install.packages("igraph")
if (!require("ggraph")) install.packages("ggraph")
library(igraph)
library(ggraph)

# --- Improvements ---
# 1. Avoid duplicate edges (undirected graph)
# 2. Allow Jaccard threshold as a parameter
# 3. Improve node labeling (shorten or wrap if needed)
# 4. Add comments for clarity

network_jaccard_threshold <- 0.2 # set as parameter for easy adjustment

# Use refined top_terms for network
network_terms <- combined_df %>%
  filter(Description %in% top_terms) %>%
  mutate(core_genes = strsplit(as.character(core_enrichment), "/")) %>%
  dplyr::select(Description, Comparison, core_genes)

# Map term to gene sets
term_gene_map <- network_terms %>%
  group_by(Description) %>%
  reframe(Genes = list(unique(unlist(core_genes))))

# Only proceed if there are at least 2 terms with genes
if (nrow(term_gene_map) >= 2 && all(lengths(term_gene_map$Genes) > 0)) {
  # Create all unique pairs of terms (no self-pairs, undirected)
  term_pairs <- t(combn(term_gene_map$Description, 2))
  term_pairs_df <- data.frame(Term1 = term_pairs[,1], Term2 = term_pairs[,2], stringsAsFactors = FALSE)

  # Compute Jaccard similarity for each pair
  get_genes <- function(term) term_gene_map$Genes[term_gene_map$Description == term][[1]]
  term_pairs_df$Jaccard <- mapply(function(a, b) {
    genes_a <- get_genes(a)
    genes_b <- get_genes(b)
    if (length(genes_a) == 0 || length(genes_b) == 0) return(0)
    length(intersect(genes_a, genes_b)) / length(union(genes_a, genes_b))
  }, term_pairs_df$Term1, term_pairs_df$Term2)

  # Filter edges: only show edges with Jaccard > threshold
  edges <- term_pairs_df %>%
    filter(Jaccard > network_jaccard_threshold)

  # Nodes: GO terms, add NES and significance info
  nodes <- combined_df %>%
    filter(Description %in% top_terms) %>%
    group_by(Description) %>%
    summarise(Max_NES = NES[which.max(abs(NES))], Min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    mutate(Significant = Min_padj < 0.05)

  # Optionally shorten/wrap node labels for clarity
  nodes$Label <- stringr::str_wrap(nodes$Description, width = 40)

  # Only plot if there are at least 2 nodes and at least 1 edge
  if (nrow(nodes) >= 2 && nrow(edges) > 0) {
    net <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
    network_plot_file <- file.path(subdirs$plots_main, "GO_Term_Overlap_Network.svg")
    svg(network_plot_file, width = 8, height = 8, family = "sans")
    print(
      ggraph(net, layout = "fr") +
        geom_edge_link(aes(width = Jaccard), color = "#4c87c6", alpha = 0.4) +
        geom_node_point(aes(color = Max_NES, size = Significant), show.legend = TRUE) +
        geom_node_text(aes(label = Label), repel = TRUE, size = 3) +
        scale_color_gradientn(colours = custom_palette(100), name = "Max NES") +
        scale_size_manual(values = c(`TRUE` = 6, `FALSE` = 3), name = "Significant") +
        labs(title = paste0("GO Term Overlap Network (Jaccard > ", network_jaccard_threshold, ")")) +
        theme_void()
    )
    dev.off()
    message("GO Term Overlap Network plot saved to: ", network_plot_file)
  } else {
    message("Not enough nodes or edges for GO Term Overlap Network plot. Skipping plot generation.")
  }
} else {
  message("Not enough GO terms with genes for network plot. Skipping plot generation.")
}



# -----------------------------------------------------
# Infer Cell Type from Protein Content per Comparison
# -----------------------------------------------------

# Define marker lists for cell types
microglia_markers <- c("P97797", "Q9QZS7", "Q9QZC2", "Q9QZC7", "Q9QZC4", "Q9QZC3", "Q9QZC5", "Q9QZC6", "Q9QZC8", "Q9QZC9", "Q9QZD0") # Example markers
astrocyte_markers <- c("P63038", "P63039", "Q9Z2D6", "Q9Z2D7", "Q9Z2D8", "Q9Z2D9", "Q9Z2E0", "Q9Z2E1", "Q9Z2E2", "Q9Z2E3") # Example markers
neuron_markers <- c("P62256", "P62257", "P62258", "P62259", "P62260", "P62261", "P62262", "P62263", "P62264", "P62265") # Example markers

# Expanded marker lists (add more UniProt IDs as needed)
microglia_markers <- c(
  microglia_markers,
  "Q9QZC1", "Q9QZC0", "Q9QZB9", "Q9QZB8", "Q9QZB7", "Q9QZB6", "Q9QZB5", "Q9QZB4", "Q9QZB3", "Q9QZB2", "Q9QZB1", "Q9QZB0",
  "Q9QZA9", "Q9QZA8", "Q9QZA7", "Q9QZA6", "Q9QZA5", "Q9QZA4", "Q9QZA3", "Q9QZA2", "Q9QZA1", "Q9QZA0"
)
astrocyte_markers <- c(
  astrocyte_markers,
  "Q9Z2E4", "Q9Z2E5", "Q9Z2E6", "Q9Z2E7", "Q9Z2E8", "Q9Z2E9", "Q9Z2F0", "Q9Z2F1", "Q9Z2F2", "Q9Z2F3", "Q9Z2F4", "Q9Z2F5",
  "Q9Z2F6", "Q9Z2F7", "Q9Z2F8", "Q9Z2F9", "Q9Z2G0", "Q9Z2G1", "Q9Z2G2", "Q9Z2G3"
)
neuron_markers <- c(
  neuron_markers,
  "P62266", "P62267", "P62268", "P62269", "P62270", "P62271", "P62272", "P62273", "P62274", "P62275", "P62276", "P62277",
  "P62278", "P62279", "P62280", "P62281", "P62282", "P62283", "P62284", "P62285"
)

# Function to infer cell type from a vector of proteins
infer_celltype <- function(proteins) {
  n_micro <- sum(proteins %in% microglia_markers)
  n_astro <- sum(proteins %in% astrocyte_markers)
  n_neuron <- sum(proteins %in% neuron_markers)
  counts <- c(microglia = n_micro, astrocyte = n_astro, neuron = n_neuron)
  if (all(counts == 0)) return(NA_character_)
  return(names(which.max(counts)))
}


# For each comparison, get all unique proteins and infer cell type, plus marker counts
comparison_proteins <- core_genes_df %>%
  group_by(Comparison) %>%
  summarise(Proteins = list(unique(core_enrichment)), .groups = "drop")

# Compute marker counts for each comparison
get_marker_counts <- function(proteins) {
  n_micro <- sum(proteins %in% microglia_markers)
  n_astro <- sum(proteins %in% astrocyte_markers)
  n_neuron <- sum(proteins %in% neuron_markers)
  data.frame(Microglia_Markers = n_micro, Astrocyte_Markers = n_astro, Neuron_Markers = n_neuron)
}

marker_counts_df <- bind_rows(lapply(comparison_proteins$Proteins, get_marker_counts))

comparison_proteins$Inferred_CellType <- vapply(comparison_proteins$Proteins, infer_celltype, character(1))
comparison_proteins <- bind_cols(comparison_proteins, marker_counts_df)

# Print summary
print(comparison_proteins %>% dplyr::select(Comparison, Inferred_CellType, Microglia_Markers, Astrocyte_Markers, Neuron_Markers))

# Save summary CSV
write.csv(comparison_proteins %>% dplyr::select(Comparison, Inferred_CellType, Microglia_Markers, Astrocyte_Markers, Neuron_Markers),
          file = file.path(subdirs$tables, "Inferred_CellTypes_Per_Comparison.csv"),
          row.names = FALSE)

# Save detailed Excel file: summary, marker counts, and protein lists
library(writexl)
celltype_summary <- comparison_proteins %>% dplyr::select(Comparison, Inferred_CellType, Microglia_Markers, Astrocyte_Markers, Neuron_Markers)
celltype_marker_counts <- comparison_proteins %>% dplyr::select(Comparison, Microglia_Markers, Astrocyte_Markers, Neuron_Markers)
celltype_protein_lists <- comparison_proteins %>% dplyr::select(Comparison, Proteins)
# Unnest protein lists for easier viewing
celltype_protein_long <- celltype_protein_lists %>% tidyr::unnest(Proteins)
names(celltype_protein_long)[2] <- "Protein"

writexl::write_xlsx(
  list(
    Summary = celltype_summary,
    Marker_Counts = celltype_marker_counts,
    Protein_Lists = celltype_protein_long
  ),
  path = file.path(subdirs$tables, "Inferred_CellTypes_Per_Comparison.xlsx")
)


# -----------------------------------------------------
# Create log file with session info and analysis summary
# -----------------------------------------------------

log_file <- file.path(main_output_dir, paste0("compareGO_analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_lines <- list()
log_lines[[length(log_lines)+1]] <- paste0("compareGO.r analysis log - ", Sys.time())
log_lines[[length(log_lines)+1]] <- paste0("R version: ", R.version.string)
log_lines[[length(log_lines)+1]] <- paste0("Platform: ", R.version$platform)
log_lines[[length(log_lines)+1]] <- paste0("Working directory: ", getwd())
log_lines[[length(log_lines)+1]] <- paste0("Base project path: ", base_project_path)
log_lines[[length(log_lines)+1]] <- paste0("Ontology: ", ont)
log_lines[[length(log_lines)+1]] <- paste0("Ensemble profiling: ", ensemble_profiling)
log_lines[[length(log_lines)+1]] <- paste0("Condition: ", condition)
log_lines[[length(log_lines)+1]] <- ""

# Session info
log_lines[[length(log_lines)+1]] <- "--- Session Info ---"
log_lines[[length(log_lines)+1]] <- capture.output(sessionInfo())
log_lines[[length(log_lines)+1]] <- ""

# Analysis steps summary
log_lines[[length(log_lines)+1]] <- "--- Analysis Steps and Results ---"

# 1. Data import
if (length(file_paths) > 0) {
  log_lines[[length(log_lines)+1]] <- paste0("Imported ", length(file_paths), " enrichment CSV files from ", input_dir)
} else {
  log_lines[[length(log_lines)+1]] <- paste0("No enrichment CSV files found in ", input_dir, ". Downstream analyses skipped.")
}

# 2. Combined enrichment data
if (exists("combined_df") && nrow(combined_df) > 0) {
  log_lines[[length(log_lines)+1]] <- paste0("Combined enrichment data: ", nrow(combined_df), " rows.")
} else {
  log_lines[[length(log_lines)+1]] <- "No combined enrichment data available."
}

# 3. Top term selection
if (exists("top_terms") && length(top_terms) > 0) {
  log_lines[[length(log_lines)+1]] <- paste0("Selected ", length(top_terms), " top GO terms for visualization (redundancy filtered).")
} else {
  log_lines[[length(log_lines)+1]] <- "No significant GO terms found. Top term selection and downstream plots skipped."
}

# 4. Plots
if (exists("comparison_plot_data") && nrow(comparison_plot_data) > 0) {
  log_lines[[length(log_lines)+1]] <- "Generated comparison selection dotplot."
} else {
  log_lines[[length(log_lines)+1]] <- "Comparison selection dotplot skipped (no data)."
}
if (exists("heatmap_data") && nrow(heatmap_data) > 0) {
  log_lines[[length(log_lines)+1]] <- "Generated enrichment heatmap."
} else {
  log_lines[[length(log_lines)+1]] <- "Enrichment heatmap skipped (no data)."
}
if (exists("dotplot") && exists("output_dotplot") && file.exists(output_dotplot)) {
  log_lines[[length(log_lines)+1]] <- "Generated dotplot for top terms per comparison."
} else {
  log_lines[[length(log_lines)+1]] <- "Dotplot for top terms per comparison skipped."
}
if (exists("dotplot_top") && exists("output_dotplot_top") && file.exists(output_dotplot_top)) {
  log_lines[[length(log_lines)+1]] <- "Generated dotplot for overall top terms."
} else {
  log_lines[[length(log_lines)+1]] <- "Dotplot for overall top terms skipped."
}

# 5. Core gene and Jaccard analysis
if (exists("core_genes_df") && nrow(core_genes_df) > 0) {
  log_lines[[length(log_lines)+1]] <- paste0("Extracted core genes for all comparisons (", nrow(core_genes_df), " rows).")
} else {
  log_lines[[length(log_lines)+1]] <- "Core gene extraction skipped (no data)."
}
if (exists("jaccard_matrix") && nrow(jaccard_matrix) > 0) {
  log_lines[[length(log_lines)+1]] <- "Computed Jaccard similarity matrix and generated heatmap."
} else {
  log_lines[[length(log_lines)+1]] <- "Jaccard similarity analysis skipped (no data)."
}

# 6. Per-comparison volcano plots and heatmaps
if (exists("log2fc_long") && nrow(log2fc_long) > 0) {
  log_lines[[length(log_lines)+1]] <- paste0("Processed log2FC data for ", length(unique(log2fc_long$Comparison)), " comparisons.")
} else {
  log_lines[[length(log_lines)+1]] <- "Log2FC data import skipped or failed."
}

# 7. GO enrichment on up/downregulated proteins
if (length(comparison_files) > 0) {
  log_lines[[length(log_lines)+1]] <- paste0("Ran GO enrichment for up/downregulated proteins in ", length(comparisons), " comparisons.")
} else {
  log_lines[[length(log_lines)+1]] <- "GO enrichment on up/downregulated proteins skipped (no significant proteins found)."
}

# 8. simplifyGOFromMultipleLists
if (exists("go_results_list") && length(go_results_list) >= 2) {
  log_lines[[length(log_lines)+1]] <- paste0("Ran simplifyGOFromMultipleLists for ", length(go_results_list), " lists.")
} else {
  log_lines[[length(log_lines)+1]] <- "simplifyGOFromMultipleLists skipped (not enough lists with data)."
}

# 9. GO term overlap network
if (exists("network_terms") && nrow(network_terms) > 1 && exists("edges") && nrow(edges) > 0) {
  log_lines[[length(log_lines)+1]] <- "Generated GO term overlap network plot."
} else {
  log_lines[[length(log_lines)+1]] <- "GO term overlap network plot skipped (not enough nodes/edges)."
}

# 10. Cell type inference
if (exists("comparison_proteins") && nrow(comparison_proteins) > 0) {
  log_lines[[length(log_lines)+1]] <- "Inferred cell types for each comparison based on marker proteins."
} else {
  log_lines[[length(log_lines)+1]] <- "Cell type inference skipped (no protein data)."
}

# 11. Supplementary tables
if (exists("supp_output_path") && file.exists(supp_output_path)) {
  log_lines[[length(log_lines)+1]] <- "Generated supplementary tables Excel file."
} else {
  log_lines[[length(log_lines)+1]] <- "Supplementary tables not generated."
}

# Write log file
writeLines(unlist(log_lines), con = log_file)
message("Analysis log written to: ", log_file)


