#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "enrichment_io.R"))

required_pkgs <- c("dplyr", "readr", "tibble", "stringr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

claim_columns <- c(
  "claim_id", "dataset", "cell_type", "region", "layer", "region_layer", "contrast",
  "group_direction", "biological_program", "evidence_type", "statistic_name",
  "statistic_value", "p_value", "fdr", "n", "core_genes", "module_id",
  "module_label", "source_file", "figure_panel_candidate", "interpretation_strength",
  "limitations"
)

empty_claims <- function() {
  out <- as.data.frame(setNames(rep(list(character()), length(claim_columns)), claim_columns), stringsAsFactors = FALSE)
  out
}

standardize_claims <- function(df) {
  if (is.null(df) || !nrow(df)) return(empty_claims())
  for (col in claim_columns) if (!col %in% names(df)) df[[col]] <- NA
  df <- df[, claim_columns, drop = FALSE]
  df
}

latest_csv <- function(root, pattern) latest_file(root, pattern)

collect_program_claims <- function(dataset) {
  f <- path_results("tables", "04_differential_expression_enrichment", "biological_program_summary", dataset, "program_summary.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"biological_program" %in% names(df)) return(empty_claims())
  df %>%
    dplyr::transmute(
      dataset = dataset,
      cell_type = dataset,
      contrast = .data$comparison,
      group_direction = .data$direction,
      biological_program = .data$biological_program,
      evidence_type = "GSEA_program_summary",
      statistic_name = "representative_NES",
      statistic_value = .data$representative_NES,
      fdr = .data$min_fdr,
      core_genes = .data$core_genes,
      source_file = f,
      figure_panel_candidate = "program_atlas_heatmap",
      interpretation_strength = vapply(seq_along(.data$min_fdr), function(i) interpretation_strength(fdr = .data$min_fdr[[i]], effect_size = .data$representative_NES[[i]], n = NA), character(1)),
      limitations = "Program mapping is regex-based and should be interpreted as thematic synthesis."
    ) %>%
    standardize_claims()
}

collect_wgcna_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_evidence_rank.csv")
  if (!file.exists(f)) f <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_priority_summary.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"ModuleID" %in% names(df)) return(empty_claims())
  for (col in c("Supermodule", "ModuleLabel_Final", "strongest_condition_contrast", "condition_model_fdr",
                "condition_model_p", "strongest_condition_adjusted_delta", "evidence_warning")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df %>%
    dplyr::transmute(
      dataset = dataset,
      cell_type = dataset,
      contrast = .data$strongest_condition_contrast,
      biological_program = dplyr::coalesce(.data$Supermodule, .data$ModuleLabel_Final),
      evidence_type = "WGCNA_module",
      statistic_name = "condition_model_fdr",
      statistic_value = .data$condition_model_fdr,
      p_value = .data$condition_model_p,
      fdr = .data$condition_model_fdr,
      module_id = .data$ModuleID,
      module_label = .data$ModuleLabel_Final,
      source_file = f,
      figure_panel_candidate = "WGCNA_module_priority",
      interpretation_strength = vapply(seq_along(.data$condition_model_fdr), function(i) interpretation_strength(fdr = .data$condition_model_fdr[[i]], effect_size = .data$strongest_condition_adjusted_delta[[i]], n = NA), character(1)),
      limitations = dplyr::coalesce(.data$evidence_warning, "Low-n module evidence; prioritize replicated or convergent modules.")
    ) %>%
    standardize_claims()
}

collect_overlap_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "05_wgcna_de_gsea_overlap", dataset, "WGCNA_vs_DE_GSEA_overlap.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"ModuleID" %in% names(df)) return(empty_claims())
  if ("status" %in% names(df)) return(empty_claims())
  df %>%
    dplyr::arrange(.data$fisher_fdr) %>%
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      dataset = dataset,
      cell_type = dataset,
      contrast = .data$contrast,
      evidence_type = "WGCNA_DE_GSEA_overlap",
      statistic_name = "fisher_fdr",
      statistic_value = .data$fisher_fdr,
      p_value = .data$fisher_p,
      fdr = .data$fisher_fdr,
      n = .data$n_DE_overlap,
      core_genes = .data$top_overlap_proteins,
      module_id = .data$ModuleID,
      source_file = f,
      figure_panel_candidate = "WGCNA_vs_DE_GSEA_overlap",
      interpretation_strength = vapply(seq_along(.data$fisher_fdr), function(i) interpretation_strength(fdr = .data$fisher_fdr[[i]], effect_size = .data$jaccard_DE[[i]], n = .data$n_DE_overlap[[i]]), character(1)),
      limitations = "Overlap supports convergence, not causality."
    ) %>%
    standardize_claims()
}

collect_behavior_claims <- function() {
  f <- path_results("tables", "08_behavior_physio_coupling", "network_behavior_coupling", "edge_behavior_figure_ready_table.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df)) return(empty_claims())
  df %>%
    dplyr::transmute(
      dataset = current_dataset(),
      contrast = .data$Edge,
      biological_program = .data$Outcome,
      evidence_type = "network_behavior_coupling",
      statistic_name = "pearson_r",
      statistic_value = .data$estimate,
      p_value = .data$p.value,
      fdr = .data$fdr,
      n = .data$n,
      source_file = f,
      figure_panel_candidate = "edge_behavior_correlation_forest",
      interpretation_strength = .data$interpretation_strength,
      limitations = .data$limitations
    ) %>%
    standardize_claims()
}

if (is_dry_run()) {
  dry_run_line("Script", "09_export_pride_journal/06_make_biological_claims_table.R")
  dry_run_line("Datasets", paste(valid_datasets(), collapse = ", "))
  dry_run_line("Output CSV", path_results("tables", "biological_claims_table.csv"))
  dry_run_line("Output XLSX", path_results("tables", "biological_claims_table.xlsx"))
  quit(status = 0, save = "no")
}

claims <- dplyr::bind_rows(
  lapply(valid_datasets(), collect_program_claims),
  lapply(valid_datasets(), collect_wgcna_claims),
  lapply(valid_datasets(), collect_overlap_claims),
  list(collect_behavior_claims())
) %>%
  standardize_claims() %>%
  dplyr::mutate(
    claim_id = sprintf("CLAIM_%04d", dplyr::row_number()),
    interpretation_strength = dplyr::coalesce(.data$interpretation_strength, "exploratory"),
    limitations = dplyr::coalesce(.data$limitations, "Conservative evidence table; review source files before manuscript claims.")
  )

dir_create(path_results("tables"))
csv_out <- path_results("tables", "biological_claims_table.csv")
xlsx_out <- path_results("tables", "biological_claims_table.xlsx")
readr::write_csv(claims, csv_out, na = "")
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(biological_claims = claims), xlsx_out)
}

message("Biological claims table written: ", csv_out)
