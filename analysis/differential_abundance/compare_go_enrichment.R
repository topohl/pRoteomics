# ================================================================
# Script: analysis/differential_abundance/compare_go_enrichment.R
# Stage: enrichment
# Scope: dataset_specific
# Consumes: required results/differential_abundance/run_clusterprofiler_enrichment/<dataset>/models/clusterProfiler_manifest.csv; data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv; optional none declared in pipeline.yml
# Produces: results/differential_abundance/compare_go_enrichment/<dataset>/models/compareGO_input_manifest.csv; results/differential_abundance/compare_go_enrichment/<dataset>/tables
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Consumes clusterProfiler output.
# ================================================================

#' compareGO.r - Comparative Gene Ontology Enrichment Analysis and Visualization
#'
#' Consumes:
#'   - data/processed/04_differential_expression_enrichment/clusterProfiler/clusterProfiler_manifest.csv
#'   - mapped/log2FC contrast CSVs listed in that manifest
#' Produces:
#'   - compareGO tables, figures, source data, reports and logs under canonical
#'     data/processed and results folders for 04_differential_expression_enrichment/compareGO
#' Contract:
#'   - reads clusterProfiler_manifest rather than recursively discovering CSVs
#'   - validates ID, Description, NES, p.adjust, setSize, core_enrichment
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
#'   - CSV files listed in clusterProfiler_manifest.csv.
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
#  
#  

script_file <- local({
  frame_files <- vapply(sys.frames(), function(frame) {
    file <- if (!is.null(frame$ofile)) frame$ofile else NA_character_
    as.character(file)[1]
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) > 0) normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = FALSE) else NA_character_
})
script_dir <- if (!is.na(script_file)) dirname(script_file) else getwd()
paths_candidates <- c(
  repo_path("R", "paths.R"),
  file.path(getwd(), "..", "R", "paths.R"),
  repo_path("R", "paths.R"),
  file.path(script_dir, "..", "R", "paths.R")
)
paths_file <- normalizePath(paths_candidates[file.exists(paths_candidates)][1], winslash = "/", mustWork = FALSE)
if (is.na(paths_file) || !file.exists(paths_file)) {
  stop("Could not find R/paths.R. Tried:\n", paste(paths_candidates, collapse = "\n"), call. = FALSE)
}
if (!nzchar(Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = ""))) {
  Sys.setenv(PROTEOMICS_PROJECT_ROOT = dirname(dirname(paths_file)))
}
source(paths_file)
source(repo_path("R", "plotting_nature.R"))  # NATURE_REPEL_SEED
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "schema_validation.R"))
source(repo_path("R", "enrichment_plots.R"))
source(repo_path("R", "differential_abundance_paths.R"))

# Phase 6G.4: destinations resolve through the normalized output contract,
# addressed by this analysis's own identity rather than by the historical
# 04_differential_expression_enrichment stage directory. Outputs already
# written there stay exactly where they are and are read, never rewritten.
ANALYSIS_ID <- "compare_go_enrichment"
MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- "compareGO"

# Package installation policy. Keep FALSE for reproducible, fail-fast runs.
AUTO_INSTALL_MISSING_PACKAGES <- FALSE
DRY_RUN <- is_dry_run()
LEGACY_COMPAREGO_TAIL_ENABLED <- FALSE

first_existing_path <- function(paths) {
  paths <- unique(normalizePath(paths[nzchar(paths)], winslash = "/", mustWork = FALSE))
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}
if (!isTRUE(DRY_RUN)) {
  early_config_candidates <- c(
    file.path(getwd(), "compareGO_config.yml"),
    file.path(getwd(), "compareGO_config.local.yml"),
    file.path(getwd(), "config", "compareGO_config.yml"),
    file.path(getwd(), "config", "compareGO_config.local.yml"),
    repo_path("config", "compareGO_config.local.yml"),
    repo_path("config", "compareGO_config.yml")
  )
  early_config <- early_config_candidates[file.exists(early_config_candidates)][1]
  if (!is.na(early_config)) {
    cfg_lines <- readLines(early_config, warn = FALSE)
    DRY_RUN <- any(grepl("^\\s*dry_run\\s*:\\s*true\\s*$", cfg_lines, ignore.case = TRUE))
  }
}

require_or_stop <- function(pkgs, bioc = FALSE) {
  if (isTRUE(DRY_RUN)) return(invisible(TRUE))
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  stop("Missing required packages: ", paste(missing, collapse = ", "),
       ". Install them explicitly before running this script.", call. = FALSE)
}

# -----------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------

# Ensure 'rlang' version >= 1.1.7 is installed before running this script.
if (isTRUE(LEGACY_COMPAREGO_TAIL_ENABLED) && !isTRUE(DRY_RUN) && (!requireNamespace("rlang", quietly=TRUE) || packageVersion("rlang") < "1.1.7")) {
    stop("Please install 'rlang' version >= 1.1.7 manually before running this script.")
}
if (isTRUE(LEGACY_COMPAREGO_TAIL_ENABLED) && !isTRUE(DRY_RUN) && !requireNamespace("simplifyEnrichment", quietly=TRUE)) {
    stop("Please install 'simplifyEnrichment' manually before running this script.")
}
if (isTRUE(LEGACY_COMPAREGO_TAIL_ENABLED) && !isTRUE(DRY_RUN)) {
  library(simplifyEnrichment)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(writexl)
}

cran_required <- c(
  "ggplot2", "stringr", "ggpubr", "ggthemes", "dplyr", "tidyr", "purrr",
  "readr", "pheatmap", "tibble", "RColorBrewer", "writexl", "scales",
  "ggrepel", "magick"
)
if (isTRUE(LEGACY_COMPAREGO_TAIL_ENABLED)) require_or_stop(cran_required)
if (isTRUE(LEGACY_COMPAREGO_TAIL_ENABLED) && !isTRUE(DRY_RUN)) {
  suppressPackageStartupMessages(invisible(lapply(cran_required, library, character.only = TRUE)))
}

# -----------------------------------------------------
# Define Theme and Helper Functions
# -----------------------------------------------------
#' Publication-style ggplot2 theme for publication-quality figures
theme_publication <- function(base_size = 9, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Text elements
      text = element_text(color = "#2C2C2C", family = base_family, lineheight = 1.3),
      
      # Axes: minimalist but clear
      axis.line.x = element_line(color = "#2C2C2C", linewidth = 0.6),
      axis.line.y = element_line(color = "#2C2C2C", linewidth = 0.6),
      axis.ticks = element_line(color = "#2C2C2C", linewidth = 0.5),
      axis.ticks.length = unit(3, "pt"),
      axis.text = element_text(color = "#2C2C2C", size = rel(0.95)),
      axis.title = element_text(color = "#2C2C2C", size = rel(1.0), face = "plain"),
      
      # Legend: clean and prominent
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_rect(color = NA, fill = NA),
      legend.key.size = unit(10, "pt"),
      legend.key.height = unit(10, "pt"),
      legend.title = element_text(color = "#2C2C2C", size = rel(0.95), face = "plain"),
      legend.text = element_text(color = "#2C2C2C", size = rel(0.9)),
      legend.position = "right",
      legend.justification = "top",
      legend.margin = margin(5, 5, 5, 5),
      
      # Title and subtitle
      plot.title = element_text(color = "#2C2C2C", size = rel(1.15), face = "bold", 
                                hjust = 0, vjust = 1, margin = margin(b = 8)),
      plot.subtitle = element_text(color = "#555555", size = rel(0.95), hjust = 0, 
                                   margin = margin(b = 5)),
      
      # Facets
      strip.text = element_text(color = "#2C2C2C", size = rel(0.95), face = "plain"),
      strip.background = element_rect(color = "#EEEEEE", fill = "#EEEEEE", linewidth = 0.4),
      
      # Panel
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' Dynamically calculate plot dimensions based on data content (publication-optimized)
calc_dims <- function(df_plot) {
  n_cols <- length(unique(as.character(df_plot$Comparison)))
  n_rows <- length(unique(as.character(df_plot$Description)))
  # Publication-style: optimize for single or dual column layouts (85mm or 180mm)
  w <- max(3.35, 3.35 + (n_cols * 0.5))  # Single column width ~85mm = 3.35 inches
  h <- max(5, 2.5 + (n_rows * 0.25))
  return(list(w = w, h = h))
}

# =====================================================
# UTILITY FUNCTIONS FOR ENHANCED ANALYSIS
# =====================================================

#' Data Validation Function - Check for common data issues
validate_data <- function(df, name = "data") {
  issues <- list()
  
  if (nrow(df) == 0) {
    issues[[length(issues)+1]] <- paste0(name, ": Empty dataframe")
  }
  if (all(is.na(df))) {
    issues[[length(issues)+1]] <- paste0(name, ": All NA values")
  }
  if (any(duplicated(df))) {
    n_dup <- sum(duplicated(df))
    issues[[length(issues)+1]] <- paste0(name, ": ", n_dup, " duplicate rows detected")
  }
  
  if (length(issues) > 0) {
    message("[VALIDATION] ", paste(issues, collapse = " | "))
    return(FALSE)
  }
  return(TRUE)
}

#' Generate comprehensive summary statistics
generate_summary_stats <- function(enrichment_df, comparisons) {
  stats_list <- list()
  
  for (comp in unique(comparisons)) {
    comp_df <- enrichment_df %>% filter(Comparison == comp)
    
    stats_list[[comp]] <- tibble(
      Comparison = comp,
      Total_Terms = nrow(comp_df),
      Significant_Terms = sum(comp_df$p.adjust < 0.05, na.rm = TRUE),
      Upregulated_Terms = sum(comp_df$NES > 0, na.rm = TRUE),
      Downregulated_Terms = sum(comp_df$NES < 0, na.rm = TRUE),
      Mean_NES = mean(comp_df$NES, na.rm = TRUE),
      Median_NES = median(comp_df$NES, na.rm = TRUE),
      Min_NES = min(comp_df$NES, na.rm = TRUE),
      Max_NES = max(comp_df$NES, na.rm = TRUE),
      Mean_Padj = mean(comp_df$p.adjust, na.rm = TRUE),
      Median_Padj = median(comp_df$p.adjust, na.rm = TRUE),
      Mean_SetSize = mean(comp_df$setSize, na.rm = TRUE),
      Genes_Total = length(unique(unlist(strsplit(paste(comp_df$core_enrichment, collapse = "/"), "/"))))
    )
  }
  
  bind_rows(stats_list)
}

#' Analyze term consistency across comparisons
analyze_term_consistency <- function(enrichment_df) {
  term_consistency <- enrichment_df %>%
    group_by(Description) %>%
    summarise(
      Num_Comparisons = n_distinct(Comparison),
      Comparisons = paste(unique(Comparison), collapse = "; "),
      Median_NES = median(NES, na.rm = TRUE),
      Max_NES = max(NES, na.rm = TRUE),
      Min_NES = min(NES, na.rm = TRUE),
      Mean_Padj = mean(p.adjust, na.rm = TRUE),
      Direction_Consistency = ifelse(all(NES > 0) | all(NES < 0), "Consistent", "Mixed"),
      Sig_Count = sum(p.adjust < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Num_Comparisons), desc(abs(Median_NES)))
  
  return(term_consistency)
}

#' Rank genes by frequency across enrichments
rank_gene_importance <- function(core_long_df) {
  if (!"p.adjust" %in% colnames(core_long_df)) {
    core_long_df$p.adjust <- NA_real_
  }

  core_long_df <- core_long_df %>%
    mutate(p.adjust_num = suppressWarnings(as.numeric(p.adjust)))

  gene_importance <- core_long_df %>%
    group_by(Gene) %>%
    summarise(
      Freq_Across_Terms = n_distinct(Description),
      Freq_Across_Comparisons = n_distinct(Comparison),
      Max_NES = max(abs(NES), na.rm = TRUE),
      Mean_NES_Abs = mean(abs(NES), na.rm = TRUE),
      Mean_Padj = ifelse(all(is.na(p.adjust_num)), NA_real_, mean(p.adjust_num, na.rm = TRUE)),
      Top_Terms = paste(unique(Description)[1:min(3, length(unique(Description)))], collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(Freq_Across_Terms), desc(Max_NES))
  
  return(gene_importance)
}

#' Calculate comparison similarity based on gene overlap
calc_comparison_similarity <- function(core_genes_df) {
  comparisons <- unique(core_genes_df$Comparison)
  
  comparison_genes <- comparisons %>% map(function(comp) {
    unique(core_genes_df$core_enrichment[core_genes_df$Comparison == comp])
  }) %>% setNames(comparisons)
  
  # Jaccard + Overlap coefficient
  similarity_list <- list()
  for (i in 1:length(comparisons)) {
    for (j in i:length(comparisons)) {
      comp_i <- comparisons[i]
      comp_j <- comparisons[j]
      genes_i <- comparison_genes[[i]]
      genes_j <- comparison_genes[[j]]
      
      overlap <- length(intersect(genes_i, genes_j))
      union <- length(union(genes_i, genes_j))
      jaccard <- ifelse(union > 0, overlap / union, 0)
      overlap_coeff <- overlap / min(length(genes_i), length(genes_j))
      
      similarity_list[[length(similarity_list)+1]] <- tibble(
        Comparison_1 = comp_i,
        Comparison_2 = comp_j,
        Gene_Overlap = overlap,
        Jaccard_Index = jaccard,
        Overlap_Coefficient = overlap_coeff
      )
    }
  }
  
  bind_rows(similarity_list)
}

# =====================================================
# ENHANCED LIBRARY LOADING FOR NEW ANALYSES
# =====================================================

# Load additional packages for new analyses
if (isTRUE(LEGACY_COMPAREGO_TAIL_ENABLED)) {
  require_or_stop(c("ggridges", "UpSetR", "networkD3", "alluvial", "ggalluvial"))
}
if (isTRUE(LEGACY_COMPAREGO_TAIL_ENABLED) && !isTRUE(DRY_RUN)) {
  suppressPackageStartupMessages({
    library(ggridges)
    library(UpSetR)
    library(networkD3)
    library(alluvial)
    library(ggalluvial)
  })
}

# =====================================================
# CONSOLIDATE LOG2FC DATA LOADING (Eliminate Redundancy)
# =====================================================

# -----------------------------------------------------
# Set Analysis Parameters and Directory Structure (must be defined before first use)
# -----------------------------------------------------

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

read_comparego_config <- function(config_path) {
  defaults <- list(
    dataset = current_dataset(),
    legacy_mode = FALSE,
    ontology = "BP",
    route_category = "phenotype_within_unit",
    route_unit = "",
    result_types = c("GSEA_GO"),
    run_id = "",
    clusterProfiler_config_hash = "",
    clusterProfiler_manifest = "",
    uniprot_mapping_file = path_external("MOUSE_10090_idmapping.dat"),
    significant_only = TRUE,
    target_n_terms = 5,
    redundancy_threshold = 0.7,
    min_set_size = 10,
    resume_raw_tables = TRUE,
    force_raw_recompute = FALSE
  )
  if (file.exists(config_path) && requireNamespace("yaml", quietly = TRUE)) {
    yaml_cfg <- yaml::read_yaml(config_path)
    return(utils::modifyList(defaults, yaml_cfg))
  }
  defaults
}

config_candidates <- c(
  file.path(getwd(), "compareGO_config.yml"),
  file.path(getwd(), "compareGO_config.local.yml"),
  file.path(getwd(), "config", "compareGO_config.yml"),
  file.path(getwd(), "config", "compareGO_config.local.yml"),
  repo_path("config", "compareGO_config.local.yml"),
  repo_path("config", "compareGO_config.yml")
)
comparego_config_path <- config_candidates[file.exists(config_candidates)][1] %||% config_candidates[1]
comparego_cfg <- read_comparego_config(comparego_config_path)
DRY_RUN <- is_dry_run(comparego_cfg)

as_repo_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (grepl("^([A-Za-z]:|/|~)", path)) return(path)
  repo_path(path)
}
comparego_cfg$clusterProfiler_manifest <- as_repo_path(comparego_cfg$clusterProfiler_manifest)
comparego_cfg$uniprot_mapping_file <- as_repo_path(comparego_cfg$uniprot_mapping_file)
DATASET <- current_dataset_from_cli(default = comparego_cfg$dataset %||% "neuron_neuropil")
comparego_cfg$dataset <- DATASET
CANONICAL_PATHS <- differential_abundance_dirs(ANALYSIS_ID, scope = DATASET)
if (!nzchar(as.character(comparego_cfg$clusterProfiler_manifest))) {
  comparego_cfg$clusterProfiler_manifest <- path_processed(
    MODULE_ID, "clusterProfiler", DATASET, "clusterProfiler_manifest.csv"
  )
}
comparego_cfg$clusterProfiler_manifest <- as_repo_path(comparego_cfg$clusterProfiler_manifest)
legacy_mode <- isTRUE(comparego_cfg$legacy_mode)

# Gene Ontology domain (MF, BP, CC, KEGG, custom)
ont <- as.character(comparego_cfg$ontology)
ensemble_profiling <- as.character(comparego_cfg$route_category)
condition <- as.character(comparego_cfg$route_unit)
base_project_path <- repo_root()

manifest_path <- as.character(comparego_cfg$clusterProfiler_manifest)

# Phase 1C-B1 canonical execution path. The legacy analysis tail below is retained
# only as non-executable historical code until its plotting helpers are migrated.
if (isTRUE(legacy_mode)) {
  stop("compareGO production execution requires the canonical clusterProfiler manifest; legacy_mode is disabled.", call. = FALSE)
}
canonical_manifest_path <- path_processed(
  MODULE_ID, "clusterProfiler", DATASET, "clusterProfiler_manifest.csv"
)
if (!identical(
    normalizePath(manifest_path, winslash = "/", mustWork = FALSE),
    normalizePath(canonical_manifest_path, winslash = "/", mustWork = FALSE))) {
  stop("compareGO must read the canonical dataset-specific clusterProfiler manifest: ",
    canonical_manifest_path, call. = FALSE)
}
configured_result_types <- as.character(unlist(comparego_cfg$result_types))
unsupported_configured_types <- setdiff(configured_result_types, canonical_comparego_result_types())
if (length(unsupported_configured_types)) {
  stop("compareGO configuration requests unsupported result_type value(s): ",
    paste(unsupported_configured_types, collapse = ", "), call. = FALSE)
}
if (!file.exists(manifest_path)) {
  stop("clusterProfiler manifest not found: ", manifest_path,
    "\nRun analysis/differential_abundance/run_clusterprofiler_enrichment.R first.", call. = FALSE)
}

canonical_cluster_manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
validate_clusterprofiler_manifest_contract(
  canonical_cluster_manifest, strict = TRUE, require_files = TRUE
)
canonical_scope <- canonical_cluster_manifest[
  canonical_cluster_manifest$dataset == DATASET &
    canonical_cluster_manifest$route_category == ensemble_profiling,
  , drop = FALSE
]
if (nzchar(condition)) {
  canonical_scope <- canonical_scope[
    canonical_scope$route_unit == condition, , drop = FALSE
  ]
}
if (nzchar(as.character(comparego_cfg$run_id))) {
  canonical_scope <- canonical_scope[
    canonical_scope$run_id == as.character(comparego_cfg$run_id), , drop = FALSE
  ]
}
if (nzchar(as.character(comparego_cfg$clusterProfiler_config_hash))) {
  canonical_scope <- canonical_scope[
    canonical_scope$config_hash == as.character(comparego_cfg$clusterProfiler_config_hash), , drop = FALSE
  ]
}
canonical_selected <- canonical_scope[
  canonical_scope$ontology == ont & canonical_scope$result_type %in% configured_result_types,
  , drop = FALSE
]
if (!nrow(canonical_selected)) {
  stop("No canonical clusterProfiler manifest rows matched the compareGO dataset and routing configuration.", call. = FALSE)
}
canonical_key <- paste(
  canonical_selected$dataset, canonical_selected$comparison,
  canonical_selected$result_type, canonical_selected$ontology, sep = "|"
)
if (anyDuplicated(canonical_key)) {
  stop("Duplicate canonical compareGO manifest identities detected; pin one run/config before comparison.", call. = FALSE)
}

canonical_outputs <- collect_canonical_comparego_outputs(
  rbind(
    canonical_selected,
    canonical_scope[!canonical_scope$result_type %in% canonical_comparego_result_types(), , drop = FALSE]
  ),
  strict = TRUE, require_files = TRUE
)
if (isTRUE(DRY_RUN)) {
  message("[DRY RUN] Canonical compareGO manifest and provenance validation passed for ",
    nrow(canonical_selected), " analysis row(s).")
  quit(status = 0, save = "no")
}

comparego_processed_dir <- file.path(CANONICAL_PATHS$models, DATASET)
comparego_table_dir <- file.path(
  CANONICAL_PATHS$tables, DATASET, ont, ensemble_profiling,
  if (nzchar(condition)) condition else "all_route_units"
)
dir.create(comparego_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparego_table_dir, recursive = TRUE, showWarnings = FALSE)
term_comparison_file <- file.path(comparego_table_dir, "compareGO_term_comparison.csv")
term_gene_provenance_output_file <- file.path(comparego_table_dir, "compareGO_term_gene_provenance.csv")
analysis_status_summary_file <- file.path(comparego_table_dir, "compareGO_analysis_status_summary.csv")
comparego_input_manifest_file <- file.path(comparego_processed_dir, "compareGO_input_manifest.csv")

utils::write.csv(canonical_outputs$terms, term_comparison_file, row.names = FALSE)
utils::write.csv(canonical_outputs$provenance, term_gene_provenance_output_file, row.names = FALSE)
utils::write.csv(canonical_outputs$status, analysis_status_summary_file, row.names = FALSE)

comparego_manifest <- canonical_outputs$input_manifest
comparego_manifest$input_manifest <- manifest_path
comparego_manifest$comparego_contract_version <- canonical_comparego_manifest_contract_version()
comparego_status_key <- paste(
  canonical_outputs$status$dataset, canonical_outputs$status$comparison,
  canonical_outputs$status$result_type, canonical_outputs$status$ontology, sep = "|"
)
comparego_manifest_key <- paste(
  comparego_manifest$dataset, comparego_manifest$comparison,
  comparego_manifest$result_type, comparego_manifest$ontology, sep = "|"
)
comparego_manifest$comparego_analysis_status <- canonical_outputs$status$comparego_action[
  match(comparego_manifest_key, comparego_status_key)
]
comparego_manifest$term_comparison_file <- term_comparison_file
comparego_manifest$term_gene_provenance_output_file <- term_gene_provenance_output_file
comparego_manifest$analysis_status_summary_file <- analysis_status_summary_file
comparego_manifest$output_table <- term_comparison_file
validate_comparego_manifest_contract(comparego_manifest, require_files = TRUE)
validate_table_schema(comparego_manifest, "compareGO_manifest", strict = TRUE)
utils::write.csv(comparego_manifest, comparego_input_manifest_file, row.names = FALSE)

message("[INFO] Canonical compareGO completed: ", nrow(canonical_outputs$terms),
  " term row(s), ", nrow(canonical_outputs$provenance), " term-gene provenance row(s).")
quit(status = 0, save = "no")
