#!/usr/bin/env Rscript
# ================================================================
# Script: 09_export_pride_journal/07_make_biological_claims_table.R
# Stage: export
# Scope: global
# Consumes: required results/tables/04_differential_expression_enrichment/; results/tables/06_modules_WGCNA/; optional results/tables/08_behavior_physio_coupling/network_behavior_coupling/.
# Produces: results/tables/biological_claims_table.csv; results/tables/biological_claims_table.xlsx.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Final biological claims table.
# ================================================================

# Manuscript/journal biological claims index — not a PRIDE-required deposition artifact.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "schema_validation.R"))
source(repo_path("R", "final_evidence_bundle_utils.R"))

required_pkgs <- c("dplyr", "readr", "tibble", "stringr", "tidyr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

claim_columns <- c(
  "claim_id", "dataset", "region", "layer_cell_compartment", "contrast",
  "biological_program", "direction", "key_proteins_genes", "evidence_type",
  "effect_size_NES", "raw_p", "FDR", "robustness_stability_metric",
  "source_file", "figure_table_target", "interpretation_note",
  "claim_grade", "primary_evidence", "orthogonal_support", "major_limitation",
  "safe_interpretation", "unsafe_overinterpretation",
  "claim_allowed", "claim_gate_status", "claim_downgrade_reason",
  "primary_model_status", "animal_level_gate", "qc_gate",
  "missingness_gate", "batch_confound_gate", "marker_contamination_gate",
  "microglia_roi_gate", "neuropil_independence_gate", "robustness_gate",
  "evidence_independence_gate",
  "missingness_confounded", "plate_or_batch_confounded",
  "region_layer_imbalance_risk", "animal_pseudoreplication_risk",
  "early_pc_association", "marker_contamination_risk",
  "qc_interpretation_flag", "animal_level_status"
)

numeric_claim_columns <- c("effect_size_NES", "raw_p", "FDR")
logical_claim_columns <- c("claim_allowed")
character_claim_columns <- setdiff(claim_columns, c(numeric_claim_columns, logical_claim_columns))

empty_claims <- function() {
  out <- as.data.frame(setNames(rep(list(character()), length(claim_columns)), claim_columns), stringsAsFactors = FALSE)
  for (col in numeric_claim_columns) out[[col]] <- numeric()
  out
}

standardize_claims <- function(df) {
  if (is.null(df) || !nrow(df)) return(empty_claims())
  for (col in claim_columns) if (!col %in% names(df)) df[[col]] <- NA
  df <- df[, claim_columns, drop = FALSE]
  for (col in numeric_claim_columns) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  for (col in logical_claim_columns) df[[col]] <- as.logical(df[[col]])
  for (col in character_claim_columns) df[[col]] <- as.character(df[[col]])
  df
}

grade_claim <- function(evidence_type, fdr, interpretation_note) {
  evidence_type <- as.character(evidence_type)
  note <- tolower(as.character(interpretation_note))
  has_da <- grepl("de|differential|overlap", tolower(evidence_type))
  has_enrichment <- grepl("gsea|enrichment|program", tolower(evidence_type))
  has_module <- grepl("wgcna|module|network", tolower(evidence_type))
  fdr <- suppressWarnings(as.numeric(fdr))
  if (grepl("confound|unsafe|ambiguous", note)) return("X")
  if (isTRUE(has_da && has_enrichment && has_module)) return("A")
  if (isTRUE(has_da && has_enrichment)) return("B")
  if (isTRUE(has_enrichment || has_module)) return("C")
  if (!is.na(fdr) && fdr <= 0.05) return("C")
  "D"
}

primary_evidence_label <- function(evidence_type) {
  dplyr::case_when(
    grepl("WGCNA_DE_GSEA_overlap", evidence_type, ignore.case = TRUE) ~ "DA/GSEA overlap with WGCNA module evidence",
    grepl("GSEA|program|enrichment", evidence_type, ignore.case = TRUE) ~ "Enrichment/program summary",
    grepl("WGCNA|module", evidence_type, ignore.case = TRUE) ~ "Module/network evidence",
    grepl("behavior", evidence_type, ignore.case = TRUE) ~ "Behavior/network association",
    TRUE ~ evidence_type
  )
}

blank_to_na <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- NA_character_
  x
}

contains_any <- function(x, patterns) {
  x <- tolower(paste(stats::na.omit(as.character(x)), collapse = "\n"))
  if (!nzchar(x)) return(FALSE)
  any(vapply(patterns, grepl, logical(1), x = x, fixed = FALSE))
}

read_lines_if_exists <- function(path) {
  if (!file.exists(path)) return(character())
  readLines(path, warn = FALSE)
}

claim_qc_context <- function(dataset) {
  missingness_lines <- read_lines_if_exists(path_results("reports", "03_qc_exploration", "02_missingness_diagnostics", dataset, "missingness_summary.md"))
  qc_lines <- read_lines_if_exists(path_results("reports", "03_qc_exploration", "07_qc_biology_confounding_report", dataset, "qc_biology_confounding_summary.md"))
  pca <- read_csv_if_exists(path_results("tables", "03_qc_exploration", "05_pca_confounding_qc", dataset, "PCA_confounding_summary.csv"))
  variance <- read_csv_if_exists(path_results("tables", "03_qc_exploration", "06_variance_partitioning", dataset, "group_technical_confounding_screen.csv"))
  marker <- read_csv_if_exists(path_results("tables", "03_qc_exploration", "04_marker_rank_abundance_qc", dataset, "marker_rank_summary.csv"))

  pca_flag <- "not_available"
  if (!is.null(pca) && nrow(pca)) {
    p_cols <- intersect(c("q_value", "q", "FDR", "p_adj"), names(pca))
    eta_cols <- intersect(c("eta2", "variance_explained", "r2"), names(pca))
    p_vals <- if (length(p_cols)) suppressWarnings(as.numeric(pca[[p_cols[[1]]]])) else NA_real_
    eta_vals <- if (length(eta_cols)) suppressWarnings(as.numeric(pca[[eta_cols[[1]]]])) else NA_real_
    pca_flag <- if (any(!is.na(p_vals) & p_vals <= 0.10 & (!is.na(eta_vals) & eta_vals >= 0.15), na.rm = TRUE)) "possible" else "not_detected"
  }

  batch_flag <- "not_available"
  if (!is.null(variance) && nrow(variance)) {
    v_cols <- intersect(c("cramers_v_group_by_technical", "cramers_v", "eta2"), names(variance))
    vals <- if (length(v_cols)) suppressWarnings(as.numeric(variance[[v_cols[[1]]]])) else NA_real_
    batch_flag <- if (any(vals >= 0.5, na.rm = TRUE)) "possible" else "not_detected"
  } else if (contains_any(c(missingness_lines, qc_lines), c("batch", "plate", "technical"))) {
    batch_flag <- if (contains_any(c(missingness_lines, qc_lines), c("FAIL", "WARN", "confound"))) "possible" else "not_detected"
  }

  marker_flag <- "not_available"
  if (!is.null(marker) && nrow(marker)) {
    marker_flag <- if (contains_any(unlist(marker, use.names = FALSE), c("WARN", "FAIL", "contamination", "neuropil"))) "possible" else "not_detected"
  }

  data.frame(
    dataset = dataset,
    missingness_confounded = if (contains_any(missingness_lines, c("confounding: possible", "WARN", "FAIL"))) "possible" else if (length(missingness_lines)) "not_detected" else "not_available",
    plate_or_batch_confounded = batch_flag,
    region_layer_imbalance_risk = if (dataset == "neuron_neuropil" && contains_any(qc_lines, c("region", "layer", "imbalance", "WARN", "FAIL"))) "review_qc" else if (dataset == "neuron_neuropil") "not_detected" else "not_applicable",
    animal_pseudoreplication_risk = "review_claim_animal_level_status",
    early_pc_association = pca_flag,
    marker_contamination_risk = marker_flag,
    qc_interpretation_flag = dplyr::case_when(
      contains_any(qc_lines, c("FAIL")) ~ "FAIL",
      contains_any(c(qc_lines, missingness_lines), c("WARN", "possible")) ~ "WARN",
      length(qc_lines) || length(missingness_lines) ~ "PASS",
      TRUE ~ "not_available"
    ),
    stringsAsFactors = FALSE
  )
}

infer_animal_level_status <- function(evidence_type, robustness_stability_metric, source_file) {
  txt <- tolower(paste(evidence_type, robustness_stability_metric, source_file, sep = " "))
  dplyr::case_when(
    grepl("animal_level|n_animals|animal", txt) ~ "animal_level_or_reported",
    grepl("behavior|network_behavior", txt) ~ "animal_level_expected",
    grepl("wgcna|module", txt) ~ "review_source_table",
    TRUE ~ "sample_level_or_unclear"
  )
}

gate_from_flag <- function(x, pass_values, fail_values = c("possible", "FAIL", "WARN", "review_qc")) {
  x <- as.character(x)
  dplyr::case_when(
    is.na(x) | !nzchar(trimws(x)) | x %in% c("not_available", "review_source_table") ~ "missing_evidence",
    x %in% pass_values ~ "pass",
    x %in% fail_values ~ "fail",
    TRUE ~ "missing_evidence"
  )
}

claim_text_contains <- function(df, patterns) {
  txt <- tolower(paste(
    df$biological_program, df$evidence_type, df$interpretation_note,
    df$safe_interpretation, df$primary_evidence, df$orthogonal_support,
    sep = " "
  ))
  grepl(paste(patterns, collapse = "|"), txt, perl = TRUE)
}

evidence_domain_count <- function(x) {
  vapply(as.character(x), function(z) {
    if (is.na(z) || !nzchar(trimws(z))) return(0L)
    parts <- trimws(unlist(strsplit(z, ";", fixed = TRUE), use.names = FALSE))
    length(unique(parts[nzchar(parts)]))
  }, integer(1))
}

microglia_neuropil_independence_available <- function() {
  paths <- c(
    path_results("tables", "06_modules_WGCNA", "microglia_neuropil_independence", "microglia", "microglia_neuropil_independence_effects.csv"),
    path_results("tables", "06_modules_WGCNA", "microglia_neuropil_independence", "microglia", "microglia_module_neuropil_independence_classification.csv")
  )
  any(file.exists(paths))
}

add_claim_gates <- function(claims) {
  if (is.null(claims) || !nrow(claims)) return(standardize_claims(claims))

  claims$claim_allowed <- FALSE
  claims$primary_model_status <- dplyr::case_when(
    !file.exists(as.character(claims$source_file)) ~ "missing_evidence",
    is.na(claims$FDR) ~ "missing_evidence",
    claims$FDR <= 0.10 ~ "pass",
    TRUE ~ "fail"
  )
  claims$animal_level_gate <- dplyr::case_when(
    claims$animal_level_status %in% c("animal_level_or_reported", "animal_level_expected") ~ "pass",
    claims$animal_level_status == "review_source_table" ~ "missing_evidence",
    TRUE ~ "fail"
  )
  claims$qc_gate <- gate_from_flag(claims$qc_interpretation_flag, pass_values = "PASS")
  claims$missingness_gate <- gate_from_flag(claims$missingness_confounded, pass_values = c("not_detected", "not_applicable"))
  claims$batch_confound_gate <- gate_from_flag(claims$plate_or_batch_confounded, pass_values = c("not_detected", "not_applicable"))
  claims$marker_contamination_gate <- gate_from_flag(claims$marker_contamination_risk, pass_values = c("not_detected", "not_applicable"))

  purified_or_intrinsic <- claims$dataset == "microglia" & claim_text_contains(
    claims,
    c("purified\\s+microglia", "cell[- ]intrinsic", "intrinsic\\s+microglia", "isolated\\s+microglia")
  )
  targeted_microglia_support <- claims$dataset == "microglia" &
    grepl("microglia_signature|targeted_microglia|curated_microglia", claims$evidence_type, ignore.case = TRUE)
  marker_fidelity_pass <- claims$marker_contamination_gate == "pass"
  neuropil_available <- microglia_neuropil_independence_available()

  claims$microglia_roi_gate <- dplyr::case_when(
    claims$dataset != "microglia" ~ "not_applicable",
    !purified_or_intrinsic ~ "pass",
    targeted_microglia_support & marker_fidelity_pass & neuropil_available ~ "pass",
    !targeted_microglia_support | !neuropil_available ~ "missing_evidence",
    TRUE ~ "fail"
  )
  claims$neuropil_independence_gate <- dplyr::case_when(
    claims$dataset != "microglia" ~ "not_applicable",
    purified_or_intrinsic & neuropil_available ~ "pass",
    purified_or_intrinsic & !neuropil_available ~ "missing_evidence",
    TRUE ~ "not_applicable"
  )
  claims$robustness_gate <- dplyr::case_when(
    is.na(claims$robustness_stability_metric) | !nzchar(trimws(as.character(claims$robustness_stability_metric))) ~ "missing_evidence",
    grepl("not_available|unavailable", claims$robustness_stability_metric, ignore.case = TRUE) ~ "missing_evidence",
    TRUE ~ "pass"
  )
  domains <- evidence_domain_count(claims$robustness_stability_metric)
  claims$evidence_independence_gate <- dplyr::case_when(
    grepl("overlap", claims$evidence_type, ignore.case = TRUE) ~ "pass",
    grepl("biological_integration", claims$evidence_type, ignore.case = TRUE) & domains >= 2L ~ "pass",
    grepl("WGCNA_DE_GSEA_overlap", claims$evidence_type, ignore.case = TRUE) ~ "pass",
    TRUE ~ "fail"
  )

  gate_cols <- c(
    "primary_model_status", "animal_level_gate", "qc_gate", "missingness_gate",
    "batch_confound_gate", "marker_contamination_gate", "microglia_roi_gate",
    "neuropil_independence_gate", "robustness_gate", "evidence_independence_gate"
  )
  gate_matrix <- claims[, gate_cols, drop = FALSE]
  missing_gate <- apply(gate_matrix == "missing_evidence", 1, any, na.rm = TRUE)
  fail_gate <- apply(gate_matrix == "fail", 1, any, na.rm = TRUE)
  pass_or_na <- function(x) all(x %in% c("pass", "not_applicable"))
  claims$claim_allowed <- apply(gate_matrix, 1, pass_or_na)
  claims$claim_gate_status <- dplyr::case_when(
    claims$claim_allowed ~ "allowed",
    missing_gate ~ "missing_evidence",
    fail_gate ~ "blocked",
    TRUE ~ "blocked"
  )
  claims$claim_downgrade_reason <- vapply(seq_len(nrow(claims)), function(i) {
    blockers <- gate_cols[as.character(gate_matrix[i, ]) %in% c("missing_evidence", "fail")]
    if (!length(blockers)) return("none")
    paste(paste0(blockers, "=", as.character(gate_matrix[i, blockers])), collapse = "; ")
  }, character(1))

  standardize_claims(claims)
}

supermodule_annotation_for_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "module_annotation", dataset, "WGCNA_supermodule_biological_annotation.csv")
  ann <- read_csv_if_exists(f)
  if (is.null(ann) || !nrow(ann)) return(NULL)
  for (col in c("dataset", "SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "dominant_microenvironment_class")) {
    if (!col %in% names(ann)) ann[[col]] <- NA
  }
  ann %>%
    dplyr::mutate(
      dataset = dataset,
      Supermodule_DisplayLabel = dplyr::coalesce(
        as.character(.data$Supermodule_DisplayLabel),
        as.character(.data$Supermodule_FinalLabel),
        as.character(.data$Macroprogram_Display),
        as.character(.data$SupermoduleID)
      ),
      Supermodule_DisplayLabel_annotation = .data$Supermodule_DisplayLabel,
      Supermodule_FinalLabel_annotation = .data$Supermodule_FinalLabel,
      Macroprogram_Display_annotation = .data$Macroprogram_Display
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(c("SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display")),
      names_to = "supermodule_claim_key_source",
      values_to = "supermodule_claim_key"
    ) %>%
    dplyr::filter(!is.na(.data$supermodule_claim_key), nzchar(.data$supermodule_claim_key)) %>%
    dplyr::distinct(.data$dataset, .data$supermodule_claim_key, .keep_all = TRUE) %>%
    dplyr::transmute(
      dataset = .data$dataset,
      supermodule_claim_key = .data$supermodule_claim_key,
      Supermodule_DisplayLabel = .data$Supermodule_DisplayLabel_annotation,
      Supermodule_FinalLabel = .data$Supermodule_FinalLabel_annotation,
      Macroprogram_Display = .data$Macroprogram_Display_annotation,
      dominant_microenvironment_class = .data$dominant_microenvironment_class
    )
}

latest_csv <- function(root, pattern) latest_file(root, pattern)

collect_program_claims <- function(dataset) {
  f <- path_results("tables", "04_differential_expression_enrichment", "biological_program_summary", dataset, "program_summary.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"biological_program" %in% names(df)) return(empty_claims())
  for (col in c("route_category", "route_unit", "min_raw_p", "min_fdr", "representative_NES", "key_genes", "core_genes", "n_terms", "top_term")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df %>%
    dplyr::transmute(
      dataset = dataset,
      region = .data$route_unit,
      layer_cell_compartment = .data$route_category,
      contrast = .data$comparison,
      direction = .data$direction,
      biological_program = .data$biological_program,
      key_proteins_genes = dplyr::coalesce(.data$key_genes, .data$core_genes),
      evidence_type = "GSEA_program_summary",
      effect_size_NES = .data$representative_NES,
      raw_p = .data$min_raw_p,
      FDR = .data$min_fdr,
      robustness_stability_metric = paste0("n_terms=", .data$n_terms),
      source_file = f,
      figure_table_target = "program_atlas_heatmap; manuscript program summary table",
      interpretation_note = paste0(
        vapply(seq_along(.data$min_fdr), function(i) interpretation_strength(fdr = .data$min_fdr[[i]], effect_size = .data$representative_NES[[i]], n = NA), character(1)),
        "; top term: ", dplyr::coalesce(.data$top_term, "NA"),
        "; regex-based thematic synthesis"
      )
    ) %>%
    standardize_claims()
}

collect_wgcna_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_evidence_rank.csv")
  if (!file.exists(f)) f <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_priority_summary.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"ModuleID" %in% names(df)) return(empty_claims())
  for (col in c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "SupermoduleID", "Supermodule", "ModuleLabel_Final", "strongest_condition_contrast", "condition_model_fdr",
                "condition_model_p", "strongest_condition_adjusted_delta", "evidence_warning")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  super_ann <- supermodule_annotation_for_claims(dataset)
  df <- df %>%
    dplyr::mutate(
      dataset = dataset,
      supermodule_claim_key = dplyr::coalesce(as.character(.data$SupermoduleID), as.character(.data$Supermodule), as.character(.data$Supermodule_FinalLabel), as.character(.data$Macroprogram_Display))
    )
  if (!is.null(super_ann)) {
    df <- df %>%
      dplyr::left_join(super_ann, by = c("dataset", "supermodule_claim_key"), suffix = c("", ".annotation")) %>%
      dplyr::mutate(
        Supermodule_DisplayLabel = dplyr::coalesce(.data$Supermodule_DisplayLabel.annotation, .data$Supermodule_DisplayLabel),
        Supermodule_FinalLabel = dplyr::coalesce(.data$Supermodule_FinalLabel.annotation, .data$Supermodule_FinalLabel),
        Macroprogram_Display = dplyr::coalesce(.data$Macroprogram_Display.annotation, .data$Macroprogram_Display)
      )
  }
  df %>%
    dplyr::transmute(
      dataset = dataset,
      contrast = .data$strongest_condition_contrast,
      biological_program = dplyr::coalesce(
        blank_to_na(.data$Supermodule_DisplayLabel),
        blank_to_na(.data$Supermodule_FinalLabel),
        blank_to_na(.data$Macroprogram_Display),
        blank_to_na(.data$Supermodule),
        blank_to_na(.data$ModuleLabel_Final),
        blank_to_na(.data$ModuleID)
      ),
      direction = dplyr::case_when(
        suppressWarnings(as.numeric(.data$strongest_condition_adjusted_delta)) > 0 ~ "positive_effect",
        suppressWarnings(as.numeric(.data$strongest_condition_adjusted_delta)) < 0 ~ "negative_effect",
        TRUE ~ NA_character_
      ),
      evidence_type = "WGCNA_module",
      effect_size_NES = .data$strongest_condition_adjusted_delta,
      raw_p = .data$condition_model_p,
      FDR = .data$condition_model_fdr,
      robustness_stability_metric = if ("preservation_Zsummary_median" %in% names(df)) .data$preservation_Zsummary_median else NA,
      key_proteins_genes = paste(.data$ModuleID, .data$ModuleLabel_Final, sep = ": "),
      source_file = f,
      figure_table_target = "WGCNA_module_priority; WGCNA module evidence table",
      interpretation_note = paste0(
        vapply(seq_along(.data$condition_model_fdr), function(i) interpretation_strength(fdr = .data$condition_model_fdr[[i]], effect_size = .data$strongest_condition_adjusted_delta[[i]], n = NA), character(1)),
        "; ", dplyr::coalesce(.data$evidence_warning, "Low-n module evidence; prioritize replicated or convergent modules.")
      )
    ) %>%
    standardize_claims()
}

collect_overlap_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "04_wgcna_de_gsea_overlap", dataset, "WGCNA_vs_DE_GSEA_overlap.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"ModuleID" %in% names(df)) return(empty_claims())
  if ("status" %in% names(df)) return(empty_claims())
  for (col in c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "SupermoduleID", "Supermodule", "jaccard_DE", "n_DE_overlap", "top_overlap_proteins")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  super_ann <- supermodule_annotation_for_claims(dataset)
  df <- df %>%
    dplyr::mutate(
      dataset = dataset,
      supermodule_claim_key = dplyr::coalesce(as.character(.data$SupermoduleID), as.character(.data$Supermodule), as.character(.data$Supermodule_FinalLabel), as.character(.data$Macroprogram_Display))
    )
  if (!is.null(super_ann)) {
    df <- df %>%
      dplyr::left_join(super_ann, by = c("dataset", "supermodule_claim_key"), suffix = c("", ".annotation")) %>%
      dplyr::mutate(
        Supermodule_DisplayLabel = dplyr::coalesce(.data$Supermodule_DisplayLabel.annotation, .data$Supermodule_DisplayLabel),
        Supermodule_FinalLabel = dplyr::coalesce(.data$Supermodule_FinalLabel.annotation, .data$Supermodule_FinalLabel),
        Macroprogram_Display = dplyr::coalesce(.data$Macroprogram_Display.annotation, .data$Macroprogram_Display)
      )
  }
  df %>%
    dplyr::arrange(.data$fisher_fdr) %>%
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      dataset = dataset,
      contrast = .data$contrast,
      biological_program = dplyr::coalesce(
        blank_to_na(.data$Supermodule_DisplayLabel),
        blank_to_na(.data$Supermodule_FinalLabel),
        blank_to_na(.data$Macroprogram_Display),
        blank_to_na(.data$Supermodule),
        blank_to_na(.data$ModuleID)
      ),
      direction = NA_character_,
      evidence_type = "WGCNA_DE_GSEA_overlap",
      effect_size_NES = if ("jaccard_DE" %in% names(df)) .data$jaccard_DE else NA,
      raw_p = .data$fisher_p,
      FDR = .data$fisher_fdr,
      robustness_stability_metric = paste0("n_DE_overlap=", .data$n_DE_overlap),
      key_proteins_genes = .data$top_overlap_proteins,
      source_file = f,
      figure_table_target = "WGCNA_vs_DE_GSEA_overlap",
      interpretation_note = paste0(
        vapply(seq_along(.data$fisher_fdr), function(i) interpretation_strength(fdr = .data$fisher_fdr[[i]], effect_size = .data$jaccard_DE[[i]], n = .data$n_DE_overlap[[i]]), character(1)),
        "; overlap supports convergence, not causality."
      )
    ) %>%
    standardize_claims()
}

collect_behavior_claims <- function() {
  f <- path_results("tables", "08_behavior_physio_coupling", "network_behavior_coupling", "edge_behavior_figure_ready_table.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df)) return(empty_claims())
  if (!"dataset" %in% names(df)) df$dataset <- "neuron_neuropil"
  df$dataset[is.na(df$dataset) | !nzchar(as.character(df$dataset))] <- "neuron_neuropil"
  df %>%
    dplyr::transmute(
      dataset = vapply(.data$dataset, validate_dataset, character(1), source = "behavior claim dataset"),
      contrast = .data$Edge,
      biological_program = .data$Outcome,
      direction = dplyr::case_when(.data$estimate > 0 ~ "positive_correlation", .data$estimate < 0 ~ "negative_correlation", TRUE ~ "neutral"),
      evidence_type = "network_behavior_coupling",
      effect_size_NES = .data$estimate,
      raw_p = .data$p.value,
      FDR = .data$fdr,
      robustness_stability_metric = paste0("n=", .data$n),
      source_file = f,
      figure_table_target = "edge_behavior_correlation_forest",
      interpretation_note = paste(.data$interpretation_strength, .data$limitations, sep = "; ")
    ) %>%
    standardize_claims()
}

collect_microglia_signature_claims <- function(dataset) {
  f <- path_results(
    "tables",
    "04_differential_expression_enrichment",
    "microglia_targeted_signature_enrichment",
    dataset,
    "microglia_signature_claims_ready.csv"
  )
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"signature" %in% names(df)) return(empty_claims())
  for (col in c("comparison", "left_region", "right_region", "left_unit", "left_condition", "right_condition", "claim_type", "NES", "pvalue", "padj", "matched_genes", "contrast_class")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df %>%
    dplyr::transmute(
      dataset = dataset,
      region = ifelse(!is.na(.data$left_region) & !is.na(.data$right_region), paste(.data$left_region, .data$right_region, sep = "_vs_"), .data$left_region),
      layer_cell_compartment = .data$left_unit,
      contrast = .data$comparison,
      direction = dplyr::case_when(.data$NES > 0 ~ "positive_NES", .data$NES < 0 ~ "negative_NES", TRUE ~ "neutral"),
      biological_program = .data$signature,
      key_proteins_genes = .data$matched_genes,
      evidence_type = "microglia_signature_enrichment",
      effect_size_NES = .data$NES,
      raw_p = .data$pvalue,
      FDR = .data$padj,
      robustness_stability_metric = paste0("contrast_class=", .data$contrast_class, "; claim_type=", .data$claim_type),
      source_file = f,
      figure_table_target = "microglia_signature_claims_ready",
      interpretation_note = dplyr::case_when(
        .data$claim_type == "within_region_condition_microglia_program" ~ "Within-region condition contrast; conservative microglia-supported claim.",
        .data$claim_type == "regional_microglia_program" ~ "Cross-region same-condition contrast; conservative regional microglia program claim.",
        TRUE ~ "Microglia signature enrichment claim-ready row."
      )
    ) %>%
    standardize_claims()
}

collect_integration_claims <- function() {
  f <- path_results("tables", "10_biological_integration", "manuscript_program_summary", "global", "manuscript_program_summary.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"program_key" %in% names(df)) return(empty_claims())
  for (col in c("manuscript_claim_scope", "datasets_supported", "evidence_domains", "strongest_fdr", "strongest_evidence", "safe_manuscript_sentence", "main_limitation", "qc_flag")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df <- df %>%
    dplyr::mutate(
      dataset_list = strsplit(as.character(.data$datasets_supported), ";", fixed = TRUE),
      dataset_list = lapply(.data$dataset_list, function(x) {
        x <- trimws(x)
        x <- x[x %in% valid_datasets()]
        if (length(x)) x else valid_datasets()
      })
    ) %>%
    tidyr::unnest("dataset_list")
  df %>%
    dplyr::transmute(
      dataset = .data$dataset_list,
      region = NA_character_,
      layer_cell_compartment = .data$datasets_supported,
      contrast = .data$manuscript_claim_scope,
      biological_program = .data$program_key,
      direction = NA_character_,
      key_proteins_genes = NA_character_,
      evidence_type = "biological_integration_manuscript_summary",
      effect_size_NES = NA_real_,
      raw_p = NA_real_,
      FDR = .data$strongest_fdr,
      robustness_stability_metric = .data$evidence_domains,
      source_file = f,
      figure_table_target = "cross_compartment_program_atlas; evidence_priority_matrix",
      interpretation_note = paste(.data$strongest_evidence, .data$safe_manuscript_sentence, .data$main_limitation, sep = "; ")
    ) %>%
    standardize_claims()
}

if (is_dry_run()) {
  dry_run_line("Script", "09_export_pride_journal/07_make_biological_claims_table.R")
  dry_run_line("Datasets", paste(valid_datasets(), collapse = ", "))
  dry_run_line("Integration manuscript summary", path_results("tables", "10_biological_integration", "manuscript_program_summary", "global", "manuscript_program_summary.csv"), if (file.exists(path_results("tables", "10_biological_integration", "manuscript_program_summary", "global", "manuscript_program_summary.csv"))) "PASS" else "WARN")
  dry_run_line("Output CSV", path_results("tables", "biological_claims_table.csv"))
  dry_run_line("Output XLSX", path_results("tables", "biological_claims_table.xlsx"))
  dry_run_line("Final evidence bundle", path_results("tables", "10_biological_integration", "final_evidence_bundle", "global", "final_biological_evidence_bundle.xlsx"))
  quit(status = 0, save = "no")
}

claims <- dplyr::bind_rows(
  lapply(valid_datasets(), collect_program_claims),
  lapply(valid_datasets(), collect_wgcna_claims),
  lapply(valid_datasets(), collect_overlap_claims),
  lapply(valid_datasets(), collect_microglia_signature_claims),
  list(collect_behavior_claims(), collect_integration_claims())
) %>%
  standardize_claims() %>%
  dplyr::select(-dplyr::any_of(c(
    "missingness_confounded", "plate_or_batch_confounded",
    "region_layer_imbalance_risk", "animal_pseudoreplication_risk",
    "early_pc_association", "marker_contamination_risk",
    "qc_interpretation_flag", "animal_level_status"
  ))) %>%
  dplyr::left_join(dplyr::bind_rows(lapply(valid_datasets(), claim_qc_context)), by = "dataset") %>%
  dplyr::mutate(
    claim_id = sprintf("CLAIM_%04d", dplyr::row_number()),
    interpretation_note = dplyr::coalesce(.data$interpretation_note, "No interpretation note available; review source file before manuscript use."),
    claim_grade = vapply(seq_along(.data$evidence_type), function(i) grade_claim(.data$evidence_type[[i]], .data$FDR[[i]], .data$interpretation_note[[i]]), character(1)),
    primary_evidence = primary_evidence_label(.data$evidence_type),
    orthogonal_support = dplyr::case_when(
      grepl("overlap", .data$evidence_type, ignore.case = TRUE) ~ "WGCNA module and DA/GSEA convergence",
      grepl("microglia_signature", .data$evidence_type, ignore.case = TRUE) ~ "Microglia-targeted signature evidence",
      grepl("WGCNA_module", .data$evidence_type, ignore.case = TRUE) ~ "Module evidence; seek DA/enrichment support before strong biological claims",
      TRUE ~ "Review companion DA, enrichment, module, and QC outputs"
    ),
    major_limitation = dplyr::case_when(
      .data$dataset == "microglia" ~ "Microglia dataset is microglia-enriched ROI/local microenvironment, not purified microglia; region-only interpretation.",
      .data$claim_grade %in% c("C", "D") ~ "Single evidence stream or exploratory evidence; avoid causal or cell-intrinsic claims.",
      TRUE ~ "Observational proteomics; causal direction is not established."
    ),
    safe_interpretation = dplyr::case_when(
      .data$dataset == "microglia" ~ paste0("In microglia-enriched ROIs, evidence supports a local microenvironment-associated ", .data$biological_program, " signal."),
      TRUE ~ paste0("Evidence supports an association between ", .data$biological_program, " and the specified dataset/contrast.")
    ),
    unsafe_overinterpretation = dplyr::case_when(
      .data$dataset == "microglia" ~ "Do not claim purified microglial cell-intrinsic regulation or subtractive neuropil correction from these data alone.",
      TRUE ~ "Do not claim causality, mechanism, or cell-type specificity without independent validation."
    ),
    animal_level_status = infer_animal_level_status(.data$evidence_type, .data$robustness_stability_metric, .data$source_file),
    animal_pseudoreplication_risk = dplyr::case_when(
      .data$animal_level_status %in% c("animal_level_or_reported", "animal_level_expected") ~ "lower",
      .data$animal_level_status == "review_source_table" ~ "review_source_table",
      TRUE ~ "possible"
    )
  ) %>%
  add_claim_gates()

validate_table_schema(claims, "biological_claims_table", strict = TRUE)

dir_create(path_results("tables"))
csv_out <- path_results("tables", "biological_claims_table.csv")
xlsx_out <- path_results("tables", "biological_claims_table.xlsx")
readr::write_csv(claims, csv_out, na = "")
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(biological_claims = claims), xlsx_out)
}

message("Biological claims table written: ", csv_out)

write_run_manifest(
  path_results("logs", "09_export_pride_journal", "biological_claims_table", "run_manifest.yml"),
  inputs = list(source_files = unique(claims$source_file)),
  outputs = list(csv = csv_out, xlsx = if (file.exists(xlsx_out)) xlsx_out else NA_character_),
  parameters = list(datasets = valid_datasets(), schema = "biological_claims_table"),
  notes = "Evidence-graded manuscript claim table. Missing statistics remain NA."
)

bundle <- write_final_evidence_bundle(reason = "biological_claims_table")
message("Final biological evidence bundle refreshed: ", bundle$bundle_dir)
