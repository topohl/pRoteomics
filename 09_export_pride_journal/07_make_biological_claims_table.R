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

SCRIPT_ID <- "09_export_pride_journal/07_make_biological_claims_table.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

required_pkgs <- c("dplyr", "readr", "tibble", "stringr", "tidyr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

claim_columns <- c(
  "claim_id", "dataset", "region", "layer_cell_compartment", "contrast",
  "biological_program", "direction", "key_proteins_genes", "evidence_type", "claim_type",
  "claim_use_class", "raw_top_GO_term", "representative_GO_terms",
  "semantic_parent_label", "safe_program_label", "term_label_risk",
  "label_confidence", "label_basis", "label_downgrade_reason",
  "effect_size_NES", "raw_p", "FDR", "robustness_stability_metric",
  "source_file", "figure_table_target", "interpretation_note",
  "claim_grade", "primary_evidence", "orthogonal_support", "major_limitation",
  "safe_interpretation", "unsafe_overinterpretation",
  "claim_allowed", "claim_gate_status", "claim_downgrade_reason",
  "model_fit_status", "statistical_evidence_status", "claim_gate_model_status",
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
  marker <- read_csv_if_exists(path_results("tables", "03_qc_exploration", "04_marker_rank_abundance_qc", dataset, "marker_score_summary_by_metadata.csv"))
  fidelity <- read_csv_if_exists(path_results("tables", "03_qc_exploration", "04d_compartment_marker_fidelity", dataset, "compartment_marker_fidelity_summary.csv"))

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
    marker_flag <- if (contains_any(unlist(marker, use.names = FALSE), c("WARN", "FAIL", "contamination"))) "possible" else "not_detected"
  } else if (!is.null(fidelity) && nrow(fidelity)) {
    marker_flag <- if (contains_any(unlist(fidelity, use.names = FALSE), c("WARN", "FAIL", "contamination"))) "possible" else "not_detected"
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

infer_model_fit_status <- function(claims, req_primary_model = NULL) {
  txt <- tolower(paste(claims$interpretation_note, claims$robustness_stability_metric, sep = " "))
  if (is.null(req_primary_model)) req_primary_model <- rep(TRUE, nrow(claims))
  req_primary_model[is.na(req_primary_model)] <- FALSE
  dplyr::case_when(
    !req_primary_model ~ "not_applicable",
    grepl("diagnostic_only_model_fallback|fallback_used=true|fallback_type=(?!none)", txt, perl = TRUE) ~ "fallback",
    grepl("singular_model|singular", txt) ~ "singular",
    grepl("rank_deficient_model|rank deficient", txt) ~ "rank_deficient",
    grepl("claim_allowed_model=false|emmeans_failed|emmeans failed|failed|not estimable|unavailable", txt) ~ "failed",
    TRUE ~ "pass"
  )
}

infer_statistical_evidence_status <- function(fdr, raw_p, req_primary_model = NULL) {
  fdr <- suppressWarnings(as.numeric(fdr))
  raw_p <- suppressWarnings(as.numeric(raw_p))
  if (is.null(req_primary_model)) req_primary_model <- rep(TRUE, length(fdr))
  req_primary_model[is.na(req_primary_model)] <- FALSE
  dplyr::case_when(
    !req_primary_model ~ "not_applicable",
    is.na(fdr) & is.na(raw_p) ~ "missing_p_or_FDR",
    !is.na(fdr) & fdr <= 0.10 ~ "pass",
    !is.na(raw_p) & raw_p < 0.05 ~ "nominal_only",
    !is.na(fdr) & fdr > 0.10 ~ "FDR_fail",
    TRUE ~ "missing_p_or_FDR"
  )
}

is_incomplete_wgcna_label <- function(x) {
  z <- trimws(as.character(x))
  z_low <- tolower(z)
  is.na(z) | !nzchar(z) |
    grepl("(mostly|mixed:)\\s*$", z_low) |
    grepl("[/;:]\\s*$", z_low) |
    grepl("shared roi:\\s*(mostly|mixed:)?\\s*$", z_low) |
    grepl("^sm[0-9]+\\s*(·|\\u00b7|:|-)\\s*(shared roi:\\s*)?(mostly|mixed:)?\\s*$", z_low)
}

has_incomplete_wgcna_label <- function(x) {
  z <- trimws(as.character(x))
  present <- !is.na(z) & nzchar(z)
  present & is_incomplete_wgcna_label(z)
}

contains_incomplete_wgcna_label_text <- function(x) {
  z <- tolower(as.character(x))
  z[is.na(z)] <- ""
  sm_sep <- "(?:\\s*(?:[^[:alnum:][:space:]]|:)\\s*)?"
  grepl(paste0("\\bsm[0-9]+", sm_sep, "shared\\s+roi:\\s*mostly\\b"), z, perl = TRUE) |
    grepl(paste0("\\bsm[0-9]+", sm_sep, "shared\\s+roi:\\s*mixed:\\s*(?:[\\.;,\\)]|$)"), z, perl = TRUE) |
    grepl("\\bshared\\s+roi:\\s*mostly\\b", z, perl = TRUE) |
    grepl("\\bshared\\s+roi:\\s*mixed:\\s*(?:[\\.;,\\)]|$)", z, perl = TRUE) |
    grepl("\\bshared\\s+roi:\\s*(?:mostly|mixed:)\\s+[^\\.;,\\)]*[/;:]\\s*(?:[\\.;,\\)]|$)", z, perl = TRUE) |
    grepl("\\b(?:mostly|mixed:)\\s*(?:[\\.;,\\)]|$)", z, perl = TRUE)
}

extract_wgcna_display_id <- function(...) {
  vals <- list(...)
  n <- max(vapply(vals, length, integer(1)))
  out <- rep(NA_character_, n)
  for (v in vals) {
    v <- rep(as.character(v), length.out = n)
    hit <- regexpr("\\b(?:SM[0-9]+|WGCNA_[^:;|\\s]+)", v, ignore.case = FALSE, perl = TRUE)
    found <- ifelse(hit > 0, regmatches(v, hit), NA_character_)
    take <- (is.na(out) | !nzchar(out)) & !is.na(found) & nzchar(found)
    out[take] <- found[take]
  }
  out
}

replace_incomplete_wgcna_label_text <- function(x, replacement) {
  out <- as.character(x)
  repl <- as.character(replacement)
  if (length(repl) == 1L && length(out) > 1L) repl <- rep(repl, length(out))
  repl[is.na(repl) | !nzchar(repl)] <- "shared ROI module with incomplete biological label"
  patterns <- c(
    "\\bSM[0-9]+\\s*(?:[^[:alnum:][:space:]]|:)\\s*shared\\s+ROI:\\s*(?:mostly|mixed:)\\s+[^\\.;,\\)]*[/;:]\\s*",
    "\\bshared\\s+ROI:\\s*(?:mostly|mixed:)\\s+[^\\.;,\\)]*[/;:]\\s*",
    "\\bSM[0-9]+\\s*(?:[^[:alnum:][:space:]]|:)\\s*shared\\s+ROI:\\s*mostly\\b",
    "\\bSM[0-9]+\\s*(?:[^[:alnum:][:space:]]|:)\\s*shared\\s+ROI:\\s*mixed:\\s*",
    "\\bshared\\s+ROI:\\s*mostly\\b",
    "\\bshared\\s+ROI:\\s*mixed:\\s*"
  )
  for (i in seq_along(out)) {
    if (is.na(out[i]) || !nzchar(out[i])) next
    val <- out[i]
    for (pattern in patterns) {
      val <- gsub(pattern, repl[i], val, ignore.case = TRUE, perl = TRUE)
    }
    val <- gsub("\\s+([\\.;,\\)])", "\\1", val, perl = TRUE)
    val <- gsub("\\s{2,}", " ", val, perl = TRUE)
    out[i] <- trimws(val)
  }
  out
}

first_non_incomplete <- function(...) {
  vals <- list(...)
  n <- length(vals[[1]])
  out <- rep(NA_character_, n)
  for (v in vals) {
    v <- as.character(v)
    take <- (is.na(out) | !nzchar(out)) & !is_incomplete_wgcna_label(v)
    out[take] <- v[take]
  }
  out
}

wgcna_label_lookup <- function() {
  rows <- lapply(valid_datasets(), function(ds) {
    module_path <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_module_group_effects_interpretable.csv")
    super_path <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_supermodule_group_effects_interpretable.csv")
    module <- read_csv_if_exists(module_path)
    super <- read_csv_if_exists(super_path)
    out <- list()
    if (!is.null(module) && nrow(module)) {
      for (col in c("module_id", "ModuleID", "safe_display_label", "Module_CleanPlotLabel", "module_biological_label", "module_label_display", "module_label")) if (!col %in% names(module)) module[[col]] <- NA_character_
      out$module <- module %>%
        dplyr::transmute(
          dataset = ds,
          lookup_id = dplyr::coalesce(as.character(.data$module_id), as.character(.data$ModuleID)),
          repaired_label = first_non_incomplete(.data$safe_display_label, .data$Module_CleanPlotLabel, .data$module_biological_label, .data$module_label_display, .data$module_label)
        )
    }
    if (!is.null(super) && nrow(super)) {
      for (col in c("supermodule_id", "SupermoduleID", "safe_display_label", "Supermodule_FullAnnotationLabel", "Supermodule_CompositionDisplayLabel", "Supermodule_CompositionLabel", "Supermodule_CleanPlotLabel", "Supermodule_PlotLabel", "supermodule_label")) if (!col %in% names(super)) super[[col]] <- NA_character_
      out$supermodule <- super %>%
        dplyr::transmute(
          dataset = ds,
          lookup_id = dplyr::coalesce(as.character(.data$supermodule_id), as.character(.data$SupermoduleID)),
          repaired_label = first_non_incomplete(.data$safe_display_label, .data$Supermodule_FullAnnotationLabel, .data$Supermodule_CompositionDisplayLabel, .data$Supermodule_CompositionLabel, .data$Supermodule_CleanPlotLabel, .data$Supermodule_PlotLabel, .data$supermodule_label)
        )
    }
    dplyr::bind_rows(out)
  })
  dplyr::bind_rows(rows) %>%
    dplyr::filter(!is.na(.data$lookup_id), nzchar(.data$lookup_id), !is.na(.data$repaired_label), nzchar(.data$repaired_label)) %>%
    dplyr::distinct(.data$dataset, .data$lookup_id, .keep_all = TRUE)
}

repair_incomplete_wgcna_labels <- function(claims) {
  if (is.null(claims) || !nrow(claims)) return(claims)
  is_wgcna <- grepl("WGCNA|module", claims$evidence_type, ignore.case = TRUE) | claims$claim_type == "wgcna_group_effect"
  incomplete_program <- is_wgcna & has_incomplete_wgcna_label(claims$biological_program)
  incomplete_safe_label <- is_wgcna & has_incomplete_wgcna_label(claims$safe_program_label)
  incomplete_safe_text <- is_wgcna & contains_incomplete_wgcna_label_text(claims$safe_interpretation)
  incomplete <- incomplete_program | incomplete_safe_label | incomplete_safe_text
  if (!any(incomplete, na.rm = TRUE)) return(claims)
  lookup <- wgcna_label_lookup()
  lookup_id <- sub(":.*$", "", as.character(claims$key_proteins_genes))
  lookup_id[is.na(lookup_id) | !nzchar(lookup_id)] <- as.character(claims$biological_program[is.na(lookup_id) | !nzchar(lookup_id)])
  display_id <- extract_wgcna_display_id(claims$biological_program, claims$safe_program_label, claims$safe_interpretation)
  known_lookup_id <- lookup_id %in% lookup$lookup_id
  use_display_id <- !known_lookup_id & !is.na(display_id) & nzchar(display_id)
  lookup_id[use_display_id] <- display_id[use_display_id]
  joined <- data.frame(row_id = seq_len(nrow(claims)), dataset = claims$dataset, lookup_id = lookup_id, stringsAsFactors = FALSE) %>%
    dplyr::left_join(lookup, by = c("dataset", "lookup_id"))
  has_repair <- incomplete & !is.na(joined$repaired_label) & nzchar(joined$repaired_label)
  claims$biological_program[has_repair & incomplete_program] <- joined$repaired_label[has_repair & incomplete_program]
  claims$safe_program_label[has_repair & (incomplete_program | incomplete_safe_label)] <- joined$repaired_label[has_repair & (incomplete_program | incomplete_safe_label)]
  replace_rows <- has_repair & incomplete_safe_text
  claims$safe_interpretation[replace_rows] <- replace_incomplete_wgcna_label_text(
    claims$safe_interpretation[replace_rows],
    joined$repaired_label[replace_rows]
  )
  still_bad <- incomplete & !has_repair
  claims$label_confidence[incomplete] <- "low"
  append_reason <- function(old, reason) {
    old <- as.character(old)
    old[is.na(old) | !nzchar(old) | old == "not_applicable"] <- ""
    out <- ifelse(nzchar(old), paste(old, reason, sep = "; "), reason)
    out
  }
  claims$label_downgrade_reason[incomplete] <- append_reason(claims$label_downgrade_reason[incomplete], "incomplete_wgcna_label")
  fallback_label <- "shared ROI module with incomplete biological label"
  claims$biological_program[still_bad & incomplete_program] <- fallback_label
  claims$safe_program_label[still_bad] <- fallback_label
  claims$safe_interpretation[still_bad & incomplete_safe_text] <- replace_incomplete_wgcna_label_text(
    claims$safe_interpretation[still_bad & incomplete_safe_text],
    fallback_label
  )
  unrepaired_open <- still_bad & !(claims$claim_use_class %in% "blocked")
  claims$claim_use_class[unrepaired_open] <- "annotation_only"
  claims
}

wgcna_module_claim_status <- function(dataset) {
  effects <- read_csv_if_exists(path_results("tables", "06_modules_WGCNA", "group_effects", dataset, "module_group_effects.csv"))
  if (is.null(effects) || !nrow(effects) || !"module_id" %in% names(effects)) return(NULL)
  audit <- read_csv_if_exists(path_results("reviewer_audit", "wgcna_robustness_claim_gate.csv"))
  effects$module_or_supermodule_id <- as.character(effects$module_id)
  if (!is.null(audit) && nrow(audit)) {
    audit <- audit %>%
      dplyr::filter(.data$dataset == dataset, .data$level == "module") %>%
      dplyr::select("dataset", "module_or_supermodule_id", "contrast", dplyr::any_of(c("spatial_unit", "effect_scope")), "robustness_gate", "preservation_status", "sensitivity_status", "direction_stability", "confounding_status")
    join_cols <- intersect(c("dataset", "module_or_supermodule_id", "contrast", "spatial_unit", "effect_scope"), names(audit))
    effects <- effects %>% dplyr::left_join(audit, by = join_cols)
  }
  for (col in c("claim_allowed_model", "model_downgrade_reason", "animal_level_status", "robustness_gate", "preservation_status", "sensitivity_status", "direction_stability", "confounding_status")) {
    if (!col %in% names(effects)) effects[[col]] <- NA
  }
  effects %>%
    dplyr::mutate(
      model_ok = .data$claim_allowed_model %in% TRUE,
      robust_ok = .data$robustness_gate == "pass"
    ) %>%
    dplyr::group_by(.data$dataset, .data$module_id, .data$contrast) %>%
    dplyr::summarise(
      claim_allowed_model = any(.data$model_ok, na.rm = TRUE),
      animal_level_status = dplyr::case_when(
        any(.data$animal_level_status == "repeated_sample_mixed_model", na.rm = TRUE) ~ "repeated_sample_mixed_model",
        any(.data$animal_level_status == "animal_level", na.rm = TRUE) ~ "animal_level",
        any(.data$animal_level_status == "insufficient_animals", na.rm = TRUE) ~ "insufficient_animals",
        any(.data$animal_level_status == "missing_animal_id", na.rm = TRUE) ~ "missing_animal_id",
        TRUE ~ "sample_level_or_unclear"
      ),
      robustness_gate = ifelse(any(.data$robust_ok, na.rm = TRUE), "pass", "missing_required"),
      robustness_status = paste(unique(stats::na.omit(paste0(
        "preservation_status=", .data$preservation_status,
        ";sensitivity_status=", .data$sensitivity_status,
        ";direction_stability=", .data$direction_stability,
        ";confounding_status=", .data$confounding_status
      ))), collapse = " | "),
      model_downgrade_reason = paste(unique(stats::na.omit(.data$model_downgrade_reason)), collapse = ";"),
      .groups = "drop"
    )
}

gate_from_flag <- function(x, required = TRUE, pass_values, fail_values = c("possible", "FAIL", "WARN", "review_qc")) {
  x <- as.character(x)
  req <- rep(as.logical(required), length.out = length(x))
  req[is.na(req)] <- FALSE
  missing_status <- ifelse(req, "missing_required", "missing_optional")
  dplyr::case_when(
    is.na(x) | !nzchar(trimws(x)) | x %in% c("not_available", "review_source_table") ~ missing_status,
    x %in% pass_values ~ "pass",
    x %in% fail_values ~ "fail",
    TRUE ~ missing_status
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

fdr_bin <- function(fdr) {
  fdr <- suppressWarnings(as.numeric(fdr))
  dplyr::case_when(
    is.na(fdr) ~ "missing",
    fdr <= 0.05 ~ "<=0.05",
    fdr <= 0.10 ~ "0.05-0.10",
    TRUE ~ ">0.10"
  )
}

go_label_interpretation <- function(raw_top_GO_term, representative_GO_terms, key_genes, n_terms, biological_program) {
  top <- tolower(as.character(raw_top_GO_term))
  terms <- tolower(paste(raw_top_GO_term, representative_GO_terms, collapse = " "))
  genes <- tolower(as.character(key_genes))
  n_terms <- suppressWarnings(as.numeric(n_terms))

  surface <- grepl("keratinocyte|skin|epiderm|hair cell", top)
  smooth_muscle <- grepl("smooth muscle", top)
  maternal <- grepl("maternal|pregnan|placenta|reproductive", top)
  broad_development <- grepl("development|differentiation|morphogenesis|pattern", top)
  low_specificity <- grepl("^regulation of|process$|system development|cell differentiation$", top)

  adhesion_ecm <- grepl("adhesion|extracellular matrix|\\becm\\b|barrier|junction|integrin|collagen|laminin|basement", terms)
  cytoskeleton <- grepl("cytoskeleton|actin|contract|myosin|microtubule|motility|migration|shape|morphogenesis", terms) ||
    grepl("actb|actg|tubb|tuba|myh|myl|flna|flnb|vim|vcl|itga|itgb|col|lam", genes)
  regulatory <- grepl("rna|chromatin|transcription|translation|metabolic|nucleobase|protein modification|signaling", terms)
  neural <- grepl("neuron|synap|axon|dendrit|glial|brain|hippocamp|neural", terms)

  risk <- dplyr::case_when(
    surface ~ "non_neural_surface_label",
    smooth_muscle ~ "tissue_mismatched_top_term",
    maternal ~ "tissue_mismatched_top_term",
    is.finite(n_terms) && n_terms <= 1 ~ "single_term_driven",
    low_specificity ~ "low_specificity_GO",
    broad_development ~ "broad_developmental_GO",
    TRUE ~ "none"
  )

  parent <- dplyr::case_when(
    surface & adhesion_ecm & cytoskeleton ~ "adhesion/cytoskeletal morphogenesis",
    surface & adhesion_ecm ~ "ECM/barrier-associated morphogenesis",
    surface & cytoskeleton ~ "cytoskeletal morphogenesis",
    surface ~ "surface-associated developmental GO cluster",
    smooth_muscle & cytoskeleton ~ "actin/cytoskeletal contractility or morphogenesis",
    smooth_muscle ~ "contractility-associated developmental GO cluster",
    maternal & regulatory ~ "broad regulatory/developmental process",
    maternal ~ "broad developmental process",
    adhesion_ecm & cytoskeleton ~ "adhesion/cytoskeletal remodeling",
    adhesion_ecm ~ "ECM/adhesion remodeling",
    cytoskeleton ~ "cytoskeletal remodeling",
    neural ~ "neural/synaptic developmental process",
    regulatory ~ "regulatory/metabolic process",
    broad_development ~ "broad developmental/morphogenesis process",
    TRUE ~ as.character(biological_program)
  )
  confidence <- dplyr::case_when(
    risk %in% c("none") & (adhesion_ecm | cytoskeleton | neural | regulatory) ~ "high",
    risk %in% c("non_neural_surface_label", "tissue_mismatched_top_term", "broad_developmental_GO") &
      (adhesion_ecm | cytoskeleton | regulatory | neural) & is.finite(n_terms) & n_terms >= 3 ~ "moderate",
    risk %in% c("single_term_driven", "low_specificity_GO") ~ "low",
    risk %in% c("non_neural_surface_label", "tissue_mismatched_top_term") ~ "low",
    TRUE ~ "moderate"
  )
  mixed <- !adhesion_ecm & !cytoskeleton & !regulatory & !neural & broad_development
  risk <- ifelse(mixed & risk %in% c("broad_developmental_GO", "none"), "mixed_unresolved", risk)
  confidence <- ifelse(risk == "mixed_unresolved", "low", confidence)
  basis <- paste(c(
    if (adhesion_ecm) "term_cluster:adhesion_ecm" else NULL,
    if (cytoskeleton) "term_or_gene_signal:cytoskeleton_contractility" else NULL,
    if (regulatory) "term_cluster:regulatory_metabolic" else NULL,
    if (neural) "term_cluster:neural" else NULL,
    if (is.finite(n_terms)) paste0("n_terms=", n_terms) else NULL
  ), collapse = ";")
  basis[!nzchar(basis)] <- "top_term_only_or_unresolved"
  downgrade <- dplyr::case_when(
    risk == "none" ~ "none",
    confidence == "low" ~ paste0(risk, "; use parent label only as annotation/context"),
    TRUE ~ paste0(risk, "; relabeled to conservative parent process")
  )
  data.frame(
    semantic_parent_label = parent,
    safe_program_label = parent,
    term_label_risk = risk,
    label_confidence = confidence,
    label_basis = basis,
    label_downgrade_reason = downgrade,
    stringsAsFactors = FALSE
  )
}

add_go_label_interpretation <- function(claims) {
  if (is.null(claims) || !nrow(claims)) return(claims)
  for (col in c("raw_top_GO_term", "representative_GO_terms", "semantic_parent_label", "safe_program_label", "term_label_risk", "label_confidence", "label_basis", "label_downgrade_reason")) {
    if (!col %in% names(claims)) claims[[col]] <- NA_character_
  }
  idx <- which(claims$claim_type == "enrichment_program")
  if (!length(idx)) return(claims)
  n_terms <- suppressWarnings(as.numeric(sub("^n_terms=", "", as.character(claims$robustness_stability_metric[idx]))))
  interp <- dplyr::bind_rows(lapply(seq_along(idx), function(j) {
    i <- idx[[j]]
    go_label_interpretation(
      raw_top_GO_term = claims$raw_top_GO_term[[i]],
      representative_GO_terms = claims$representative_GO_terms[[i]],
      key_genes = claims$key_proteins_genes[[i]],
      n_terms = n_terms[[j]],
      biological_program = claims$biological_program[[i]]
    )
  }))
  claims[idx, names(interp)] <- interp
  missing_safe <- is.na(claims$safe_program_label) | !nzchar(dplyr::coalesce(claims$safe_program_label, ""))
  claims$safe_program_label[missing_safe] <- claims$biological_program[missing_safe]
  for (col in c("term_label_risk", "label_confidence", "label_basis", "label_downgrade_reason")) {
    missing_col <- is.na(claims[[col]]) | !nzchar(dplyr::coalesce(claims[[col]], ""))
    claims[[col]][missing_col] <- "not_applicable"
  }
  claims
}

add_claim_use_class <- function(claims) {
  if (is.null(claims) || !nrow(claims)) return(claims)
  n_terms <- suppressWarnings(as.numeric(sub("^n_terms=", "", as.character(claims$robustness_stability_metric))))
  fdr <- suppressWarnings(as.numeric(claims$FDR))
  independent <- claims$evidence_independence_gate == "pass"
  label_risky <- claims$term_label_risk %in% c(
    "tissue_mismatched_top_term", "broad_developmental_GO",
    "low_specificity_GO", "non_neural_surface_label",
    "single_term_driven", "mixed_unresolved"
  )
  label_low <- claims$label_confidence %in% c("low", "unresolved", "not_applicable")
  wgcna_annotation_limited <- claims$evidence_type == "WGCNA_module" &
    (claims$label_confidence %in% c("low", "unresolved") |
       grepl("annotation_stable_across_thresholds=FALSE|threshold_unstable", claims$interpretation_note, ignore.case = TRUE))
  claims$claim_use_class <- dplyr::case_when(
    claims$claim_gate_status == "disallowed" ~ "blocked",
    claims$claim_type == "wgcna_group_effect" & (claims$animal_level_gate == "missing_required" | claims$robustness_gate == "missing_required") ~ "blocked",
    claims$claim_type == "microglia_signature" & claims$claim_gate_status != "allowed" ~ "blocked",
    wgcna_annotation_limited & claims$claim_allowed ~ "annotation_only",
    claims$claim_type == "enrichment_program" & !is.na(fdr) & fdr > 0.10 ~ "annotation_only",
    claims$claim_type == "enrichment_program" & !is.na(fdr) & fdr > 0.05 & fdr <= 0.10 ~ "suggestive_context",
    claims$claim_type == "enrichment_program" & (is.na(n_terms) | n_terms <= 1) & !independent ~ "annotation_only",
    claims$claim_type == "enrichment_program" & label_risky & label_low & !independent ~ "annotation_only",
    claims$claim_type == "enrichment_program" & label_risky & !label_low & !is.na(fdr) & fdr <= 0.05 ~ "supporting_claim",
    claims$claim_type == "enrichment_program" & !label_risky & !is.na(fdr) & fdr <= 0.05 ~ "supporting_claim",
    claims$claim_type == "integration_summary" & claims$claim_gate_status != "disallowed" ~ "supporting_claim",
    claims$claim_allowed ~ "suggestive_context",
    TRUE ~ "annotation_only"
  )
  claims
}

go_label_risky <- function(x) {
  x %in% c(
    "tissue_mismatched_top_term", "broad_developmental_GO",
    "low_specificity_GO", "non_neural_surface_label",
    "single_term_driven", "mixed_unresolved"
  )
}

safe_label_in_interpretation <- function(safe_label, interpretation) {
  safe_label <- tolower(trimws(as.character(safe_label)))
  interpretation <- tolower(as.character(interpretation))
  !is.na(safe_label) & nzchar(safe_label) & !is.na(interpretation) &
    mapply(grepl, pattern = safe_label, x = interpretation, fixed = TRUE)
}

sync_go_safe_interpretation <- function(claims) {
  if (is.null(claims) || !nrow(claims)) return(claims)
  risky <- claims$claim_type == "enrichment_program" & go_label_risky(claims$term_label_risk)
  use_class_limited <- claims$claim_type == "enrichment_program" &
    claims$claim_use_class %in% c("annotation_only", "suggestive_context")
  update <- risky | use_class_limited
  if (!any(update, na.rm = TRUE)) return(claims)

  raw_term <- dplyr::coalesce(claims$raw_top_GO_term, "unavailable raw top GO term")
  safe_label <- dplyr::coalesce(claims$safe_program_label, claims$semantic_parent_label, claims$biological_program)
  downgrade <- dplyr::coalesce(claims$label_downgrade_reason, "Conservative interpretation required because the GO label is low-specificity or tissue-mismatched.")
  raw_warning <- ifelse(
    risky,
    paste0(" Raw top GO term: '", raw_term, "'. Do not interpret the raw GO term literally. ", downgrade),
    ""
  )

  claims$safe_interpretation[update] <- dplyr::case_when(
    claims$claim_use_class[update] == "annotation_only" ~ paste0(
      "Annotation/audit transparency only, not a manuscript claim: retain the raw GO evidence under the conservative label '",
      safe_label[update], "'.", raw_warning[update]
    ),
    claims$claim_use_class[update] == "suggestive_context" ~ paste0(
      "Suggestive/contextual evidence only for '", safe_label[update],
      "'; do not use as a standalone manuscript claim.", raw_warning[update]
    ),
    TRUE ~ paste0(
      "Evidence supports the conservative program label '", safe_label[update],
      "' for the specified dataset/contrast.", raw_warning[update]
    )
  )
  claims
}

apply_blocked_claim_wording <- function(claims) {
  if (is.null(claims) || !nrow(claims)) return(claims)
  blocked <- !claims$claim_allowed | claims$claim_use_class == "blocked" | claims$claim_gate_status == "disallowed"
  reasons <- dplyr::coalesce(claims$claim_downgrade_reason, "unspecified gate reason")
  reasons[!nzchar(trimws(reasons))] <- "unspecified gate reason"
  claims$safe_interpretation[blocked] <- paste0(
    "Not claim-eligible; retained for audit because ",
    reasons[blocked],
    ". This result is exploratory/blocked and should not be interpreted as manuscript evidence."
  )
  claims
}

positive_blocked_wording <- function(x) {
  x <- tolower(as.character(x))
  grepl("\\bevidence supports\\b|\\bsupports an association\\b|\\bsupports a local\\b", x)
}

write_blocked_claim_wording_audit <- function(claims) {
  audit_dir <- path_results("reviewer_audit")
  dir_create(audit_dir)
  audit <- claims %>%
    dplyr::filter(!.data$claim_allowed, positive_blocked_wording(.data$safe_interpretation)) %>%
    dplyr::select(
      "claim_id", "dataset", "claim_type", "evidence_type", "claim_use_class",
      "claim_gate_status", "claim_downgrade_reason", "safe_interpretation"
    )
  readr::write_csv(audit, file.path(audit_dir, "blocked_claim_wording_audit.csv"), na = "")
  invisible(audit)
}

write_wgcna_claim_source_audit <- function(claims) {
  audit_dir <- path_results("reviewer_audit")
  dir_create(audit_dir)
  audit <- claims %>%
    dplyr::filter(grepl("WGCNA|module", .data$evidence_type, ignore.case = TRUE) | .data$claim_type == "wgcna_group_effect") %>%
    dplyr::count(
      .data$evidence_type,
      .data$source_file,
      .data$claim_use_class,
      .data$model_fit_status,
      .data$statistical_evidence_status,
      .data$robustness_gate,
      .data$claim_allowed,
      name = "n_claims"
    ) %>%
    dplyr::arrange(.data$evidence_type, .data$claim_use_class, .data$model_fit_status, .data$statistical_evidence_status)
  readr::write_csv(audit, file.path(audit_dir, "wgcna_claim_source_audit.csv"), na = "")
  invisible(audit)
}

write_wgcna_label_completeness_audit <- function(claims) {
  audit_dir <- path_results("reviewer_audit")
  dir_create(audit_dir)
  if (is.null(claims) || !nrow(claims)) {
    audit <- tibble::tibble(
      claim_id = character(), dataset = character(), evidence_type = character(),
      claim_use_class = character(), claim_allowed = logical(), label_column = character(),
      label_text = character()
    )
    readr::write_csv(audit, file.path(audit_dir, "wgcna_label_completeness_audit.csv"), na = "")
    return(invisible(audit))
  }
  is_wgcna <- grepl("WGCNA|module", claims$evidence_type, ignore.case = TRUE) | claims$claim_type == "wgcna_group_effect"
  base <- claims[is_wgcna, , drop = FALSE]
  audit <- dplyr::bind_rows(
    base %>%
      dplyr::mutate(label_column = "biological_program", label_text = as.character(.data$biological_program)) %>%
      dplyr::filter(has_incomplete_wgcna_label(.data$label_text)),
    base %>%
      dplyr::mutate(label_column = "safe_program_label", label_text = as.character(.data$safe_program_label)) %>%
      dplyr::filter(has_incomplete_wgcna_label(.data$label_text)),
    base %>%
      dplyr::mutate(label_column = "safe_interpretation", label_text = as.character(.data$safe_interpretation)) %>%
      dplyr::filter(contains_incomplete_wgcna_label_text(.data$label_text))
  ) %>%
    dplyr::select(
      "claim_id", "dataset", "evidence_type", "claim_use_class", "claim_allowed",
      "label_column", "label_text"
    ) %>%
    dplyr::arrange(.data$dataset, .data$evidence_type, .data$claim_id, .data$label_column)
  readr::write_csv(audit, file.path(audit_dir, "wgcna_label_completeness_audit.csv"), na = "")
  invisible(audit)
}

microglia_neuropil_independence_available <- function() {
  paths <- c(
    path_results("reviewer_audit", "microglia_neuropil_independence_claim_gate.csv"),
    path_results("tables", "06_modules_WGCNA", "microglia_neuropil_independence", "microglia", "microglia_neuropil_independence_effects.csv"),
    path_results("tables", "06_modules_WGCNA", "microglia_neuropil_independence", "microglia", "microglia_module_neuropil_independence_classification.csv")
  )
  any(file.exists(paths))
}

read_microglia_neuropil_claim_gate <- function() {
  f <- path_results("reviewer_audit", "microglia_neuropil_independence_claim_gate.csv")
  gate <- read_csv_if_exists(f)
  if (is.null(gate) || !nrow(gate)) return(NULL)
  for (col in c(
    "module_or_supermodule_id", "contrast", "adjustment_mode",
    "independence_classification", "claim_gate_eligible", "downgrade_reason",
    "primary_effect_claim_relevant", "primary_effect_status"
  )) {
    if (!col %in% names(gate)) gate[[col]] <- NA
  }
  gate %>%
    dplyr::filter(.data$adjustment_mode %in% c("predeclared_primary", "predeclared_secondary")) %>%
    dplyr::group_by(.data$module_or_supermodule_id, .data$contrast) %>%
    dplyr::summarise(
      neuropil_predeclared_claim_gate_eligible = any(suppressWarnings(as.logical(.data$claim_gate_eligible)), na.rm = TRUE),
      neuropil_predeclared_sensitive = any(.data$independence_classification == "neuropil_sensitive", na.rm = TRUE),
      neuropil_predeclared_inconclusive = any(.data$independence_classification %in% c("inconclusive_low_power", "inconclusive_missing_match", "inconclusive_no_primary_effect"), na.rm = TRUE) &&
        !any(suppressWarnings(as.logical(.data$claim_gate_eligible)), na.rm = TRUE),
      neuropil_predeclared_no_primary_effect = any(.data$independence_classification %in% c("diagnostic_no_primary_effect", "inconclusive_no_primary_effect"), na.rm = TRUE) |
        !any(suppressWarnings(as.logical(.data$primary_effect_claim_relevant)), na.rm = TRUE),
      neuropil_predeclared_modes = paste(unique(stats::na.omit(.data$adjustment_mode)), collapse = ";"),
      neuropil_predeclared_classification = paste(unique(stats::na.omit(.data$independence_classification)), collapse = ";"),
      neuropil_predeclared_primary_effect_status = paste(unique(stats::na.omit(.data$primary_effect_status)), collapse = ";"),
      neuropil_predeclared_downgrade_reason = paste(unique(stats::na.omit(.data$downgrade_reason)), collapse = ";"),
      .groups = "drop"
    )
}

infer_claim_type <- function(evidence_type, dataset, biological_program = NA_character_, interpretation_note = NA_character_) {
  et <- tolower(as.character(evidence_type))
  ds <- as.character(dataset)
  txt <- tolower(paste(biological_program, interpretation_note, sep = " "))
  dplyr::case_when(
    grepl("microglia_signature", et) ~ "microglia_signature",
    grepl("biological_integration", et) ~ "integration_summary",
    grepl("network_behavior|behavior", et) ~ "behavior_association",
    grepl("WGCNA_DE_GSEA_overlap", evidence_type, ignore.case = TRUE) ~ "wgcna_group_effect",
    grepl("WGCNA_module|module", evidence_type, ignore.case = TRUE) & ds == "microglia" &
      grepl("microglia|neuropil|vascular|ecm|basement|perivascular|bbb", txt) ~ "microglia_roi_context",
    grepl("WGCNA|module", evidence_type, ignore.case = TRUE) ~ "wgcna_group_effect",
    grepl("GSEA|program|enrichment", evidence_type, ignore.case = TRUE) ~ "enrichment_program",
    TRUE ~ "export_context"
  )
}

gate_required <- function(claims) {
  raw_claim_text <- tolower(paste(
    claims$biological_program, claims$evidence_type, claims$interpretation_note,
    claims$primary_evidence, claims$orthogonal_support,
    sep = " "
  ))
  purified_or_intrinsic <- claims$dataset == "microglia" & claim_text_contains(
    claims,
    c("purified\\s+microglia", "cell[- ]intrinsic", "intrinsic\\s+microglia", "isolated\\s+microglia")
  )
  microglia_specific <- claims$claim_type %in% c("microglia_signature", "microglia_roi_context") |
    (claims$dataset == "microglia" & grepl("\\bmicroglia\\b|immune|myeloid|phago|lyso|complement", raw_claim_text, perl = TRUE))
  contamination_sensitive <- claims$dataset == "microglia" &
    (claims$claim_type %in% c("wgcna_group_effect", "microglia_roi_context") |
       grepl("neuropil|vascular|ecm|basement|perivascular|bbb|shared", raw_claim_text, perl = TRUE))

  list(
    purified_or_intrinsic = purified_or_intrinsic,
    microglia_specific = microglia_specific,
    contamination_sensitive = contamination_sensitive,
    primary_model = !claims$claim_type %in% c("integration_summary", "export_context"),
    animal_level = claims$claim_type %in% c("wgcna_group_effect", "behavior_association", "microglia_signature"),
    qc = !claims$claim_type %in% c("export_context"),
    missingness = !claims$claim_type %in% c("export_context"),
    batch = !claims$claim_type %in% c("export_context"),
    marker = microglia_specific | contamination_sensitive | purified_or_intrinsic,
    microglia_roi = claims$dataset == "microglia" & (microglia_specific | contamination_sensitive | purified_or_intrinsic),
    neuropil_independence = purified_or_intrinsic | contamination_sensitive | microglia_specific,
    robustness = claims$claim_type %in% c("wgcna_group_effect", "behavior_association"),
    evidence_independence = claims$claim_type %in% c("integration_summary")
  )
}

add_claim_gates <- function(claims) {
  if (is.null(claims) || !nrow(claims)) return(standardize_claims(claims))

  claims$claim_type <- infer_claim_type(claims$evidence_type, claims$dataset, claims$biological_program, claims$interpretation_note)
  req <- gate_required(claims)
  claims$claim_allowed <- FALSE
  claims$model_fit_status <- infer_model_fit_status(claims, req$primary_model)
  claims$statistical_evidence_status <- infer_statistical_evidence_status(claims$FDR, claims$raw_p, req$primary_model)
  claims$claim_gate_model_status <- dplyr::case_when(
    !req$primary_model ~ "not_applicable",
    !file.exists(as.character(claims$source_file)) ~ "missing_required",
    claims$model_fit_status %in% c("fallback", "failed", "singular", "rank_deficient") ~ "diagnostic_only",
    claims$statistical_evidence_status == "missing_p_or_FDR" ~ "missing_required",
    claims$statistical_evidence_status == "pass" ~ "pass",
    TRUE ~ "fail"
  )
  claims$primary_model_status <- dplyr::case_when(
    claims$claim_gate_model_status == "not_applicable" ~ "not_applicable",
    claims$claim_gate_model_status == "diagnostic_only" ~ "diagnostic_only",
    claims$claim_gate_model_status == "missing_required" ~ "missing_required",
    claims$claim_gate_model_status == "pass" ~ "pass",
    TRUE ~ "fail"
  )
  claims$animal_level_gate <- dplyr::case_when(
    !req$animal_level ~ "diagnostic_only",
    claims$animal_level_status %in% c("animal_level", "repeated_sample_mixed_model") ~ "pass",
    claims$animal_level_status %in% c("missing_animal_id", "sample_level_or_unclear") ~ "missing_required",
    claims$animal_level_status %in% c("insufficient_animals") ~ "fail",
    claims$animal_level_status %in% c("animal_level_or_reported", "animal_level_expected") ~ "pass",
    claims$animal_level_status == "review_source_table" ~ "missing_required",
    TRUE ~ "fail"
  )
  claims$qc_gate <- gate_from_flag(claims$qc_interpretation_flag, required = req$qc, pass_values = "PASS")
  claims$missingness_gate <- gate_from_flag(claims$missingness_confounded, required = req$missingness, pass_values = c("not_detected", "not_applicable"))
  claims$batch_confound_gate <- gate_from_flag(claims$plate_or_batch_confounded, required = req$batch, pass_values = c("not_detected", "not_applicable"))
  claims$marker_contamination_gate <- dplyr::case_when(
    !req$marker & claims$dataset %in% c("neuron_neuropil", "neuron_soma") & claims$marker_contamination_risk %in% c("not_available", NA_character_) ~ "missing_optional",
    !req$marker ~ "diagnostic_only",
    TRUE ~ gate_from_flag(claims$marker_contamination_risk, required = TRUE, pass_values = c("not_detected", "not_applicable"))
  )

  targeted_microglia_support <- claims$dataset == "microglia" &
    grepl("microglia_signature|targeted_microglia|curated_microglia", claims$evidence_type, ignore.case = TRUE)
  marker_fidelity_pass <- claims$marker_contamination_gate == "pass"
  neuropil_available <- microglia_neuropil_independence_available()
  neuropil_gate <- read_microglia_neuropil_claim_gate()
  if (!is.null(neuropil_gate) && nrow(neuropil_gate)) {
    claim_lookup <- data.frame(
      row_id = seq_len(nrow(claims)),
      module_or_supermodule_id = sub(":.*$", "", as.character(claims$key_proteins_genes)),
      contrast = as.character(claims$contrast),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(neuropil_gate, by = c("module_or_supermodule_id", "contrast")) %>%
      dplyr::arrange(.data$row_id)
    claims$neuropil_predeclared_claim_gate_eligible <- claim_lookup$neuropil_predeclared_claim_gate_eligible
    claims$neuropil_predeclared_sensitive <- claim_lookup$neuropil_predeclared_sensitive
    claims$neuropil_predeclared_inconclusive <- claim_lookup$neuropil_predeclared_inconclusive
    claims$neuropil_predeclared_no_primary_effect <- claim_lookup$neuropil_predeclared_no_primary_effect
    claims$neuropil_predeclared_classification <- claim_lookup$neuropil_predeclared_classification
    claims$neuropil_predeclared_primary_effect_status <- claim_lookup$neuropil_predeclared_primary_effect_status
    claims$neuropil_predeclared_downgrade_reason <- claim_lookup$neuropil_predeclared_downgrade_reason
  } else {
    claims$neuropil_predeclared_claim_gate_eligible <- NA
    claims$neuropil_predeclared_sensitive <- NA
    claims$neuropil_predeclared_inconclusive <- NA
    claims$neuropil_predeclared_no_primary_effect <- NA
    claims$neuropil_predeclared_classification <- NA_character_
    claims$neuropil_predeclared_primary_effect_status <- NA_character_
    claims$neuropil_predeclared_downgrade_reason <- NA_character_
  }
  neuropil_predeclared_eligible <- suppressWarnings(as.logical(claims$neuropil_predeclared_claim_gate_eligible))
  neuropil_predeclared_sensitive <- suppressWarnings(as.logical(claims$neuropil_predeclared_sensitive))
  neuropil_predeclared_inconclusive <- suppressWarnings(as.logical(claims$neuropil_predeclared_inconclusive))
  neuropil_predeclared_no_primary_effect <- suppressWarnings(as.logical(claims$neuropil_predeclared_no_primary_effect))

  claims$microglia_roi_gate <- dplyr::case_when(
    claims$dataset != "microglia" ~ "not_applicable",
    !req$microglia_roi ~ "diagnostic_only",
    !req$purified_or_intrinsic ~ "pass",
    targeted_microglia_support & marker_fidelity_pass & isTRUE(neuropil_available) & neuropil_predeclared_eligible %in% TRUE ~ "pass",
    !targeted_microglia_support | !neuropil_available ~ "missing_required",
    TRUE ~ "fail"
  )
  claims$neuropil_independence_gate <- dplyr::case_when(
    claims$dataset != "microglia" ~ "not_applicable",
    !req$neuropil_independence ~ "not_applicable",
    neuropil_predeclared_sensitive %in% TRUE ~ "fail",
    neuropil_predeclared_eligible %in% TRUE ~ "pass",
    neuropil_predeclared_no_primary_effect %in% TRUE ~ "diagnostic_only",
    neuropil_predeclared_inconclusive %in% TRUE ~ "diagnostic_only",
    neuropil_available ~ "missing_required",
    TRUE ~ "missing_required"
  )
  claims$robustness_gate <- dplyr::case_when(
    !req$robustness ~ "missing_optional",
    grepl("robustness_gate=pass", claims$robustness_stability_metric, ignore.case = TRUE) ~ "pass",
    grepl("robustness_gate=missing_required|robustness_gate=fail", claims$robustness_stability_metric, ignore.case = TRUE) ~ "missing_required",
    is.na(claims$robustness_stability_metric) | !nzchar(trimws(as.character(claims$robustness_stability_metric))) ~ "missing_required",
    grepl("not_available|unavailable", claims$robustness_stability_metric, ignore.case = TRUE) ~ "missing_required",
    TRUE ~ "pass"
  )
  domains <- evidence_domain_count(claims$robustness_stability_metric)
  claims$evidence_independence_gate <- dplyr::case_when(
    !req$evidence_independence ~ "missing_optional",
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
  missing_required_gate <- apply(gate_matrix == "missing_required", 1, any, na.rm = TRUE)
  missing_optional_gate <- apply(gate_matrix == "missing_optional", 1, any, na.rm = TRUE)
  diagnostic_primary_gate <- req$primary_model & claims$primary_model_status == "diagnostic_only"
  fail_gate <- apply(gate_matrix == "fail", 1, any, na.rm = TRUE) | diagnostic_primary_gate
  pass_or_nonblocking <- function(x) all(x %in% c("pass", "not_applicable", "diagnostic_only", "missing_optional"))
  claims$claim_allowed <- apply(gate_matrix, 1, pass_or_nonblocking) & !diagnostic_primary_gate
  claims$claim_gate_status <- dplyr::case_when(
    fail_gate | missing_required_gate ~ "disallowed",
    claims$claim_allowed & missing_optional_gate ~ "downgraded",
    claims$claim_allowed ~ "allowed",
    TRUE ~ "disallowed"
  )
  claims$claim_downgrade_reason <- vapply(seq_len(nrow(claims)), function(i) {
    relevant <- gate_cols[as.character(gate_matrix[i, ]) %in% c("missing_required", "missing_optional", "fail")]
    if (!length(relevant)) return("none")
    paste(paste0(relevant, "=", as.character(gate_matrix[i, relevant])), collapse = "; ")
  }, character(1))

  standardize_claims(claims)
}

write_claim_gate_audits <- function(claims) {
  gate_cols <- c(
    "primary_model_status", "animal_level_gate", "qc_gate", "missingness_gate",
    "batch_confound_gate", "marker_contamination_gate", "microglia_roi_gate",
    "neuropil_independence_gate", "robustness_gate", "evidence_independence_gate"
  )
  audit_dir <- path_results("reviewer_audit")
  dir_create(audit_dir)
  req <- gate_required(claims)
  required_df <- data.frame(
    row_id = seq_len(nrow(claims)),
    primary_model_status = req$primary_model,
    animal_level_gate = req$animal_level,
    qc_gate = req$qc,
    missingness_gate = req$missingness,
    batch_confound_gate = req$batch,
    marker_contamination_gate = req$marker,
    microglia_roi_gate = req$microglia_roi,
    neuropil_independence_gate = req$neuropil_independence,
    robustness_gate = req$robustness,
    evidence_independence_gate = req$evidence_independence,
    check.names = FALSE
  ) |>
    tidyr::pivot_longer(cols = dplyr::all_of(gate_cols), names_to = "gate_name", values_to = "required_for_claim")
  availability <- claims |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    dplyr::select("row_id", "dataset", "claim_type", dplyr::all_of(gate_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(gate_cols), names_to = "gate_name", values_to = "gate_status") |>
    dplyr::left_join(required_df, by = c("row_id", "gate_name")) |>
    dplyr::count(.data$dataset, .data$claim_type, .data$gate_name, .data$gate_status, .data$required_for_claim, name = "n_claims") |>
    dplyr::arrange(.data$dataset, .data$claim_type, .data$gate_name, .data$gate_status)
  summary <- claims |>
    dplyr::group_by(.data$dataset, .data$claim_type) |>
    dplyr::summarise(
      n_claims = dplyr::n(),
      allowed = sum(.data$claim_gate_status == "allowed", na.rm = TRUE),
      downgraded = sum(.data$claim_gate_status == "downgraded", na.rm = TRUE),
      disallowed = sum(.data$claim_gate_status == "disallowed", na.rm = TRUE),
      missing_required = sum(if_any(dplyr::all_of(gate_cols), ~ .x == "missing_required"), na.rm = TRUE),
      missing_optional = sum(if_any(dplyr::all_of(gate_cols), ~ .x == "missing_optional"), na.rm = TRUE),
      not_applicable = sum(if_any(dplyr::all_of(gate_cols), ~ .x == "not_applicable"), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$dataset, .data$claim_type)
  readr::write_csv(availability, file.path(audit_dir, "claim_gate_evidence_availability.csv"), na = "")
  readr::write_csv(summary, file.path(audit_dir, "claim_gate_summary.csv"), na = "")
  invisible(list(availability = availability, summary = summary))
}

write_go_label_audits <- function(claims) {
  audit_dir <- path_results("reviewer_audit")
  dir_create(audit_dir)
  enrichment <- claims |>
    dplyr::filter(.data$claim_type == "enrichment_program") |>
    dplyr::mutate(
      FDR_bin = fdr_bin(.data$FDR),
      safe_program_label_in_safe_interpretation = safe_label_in_interpretation(.data$safe_program_label, .data$safe_interpretation),
      risky_allowed_safe_label_missing = .data$claim_allowed &
        go_label_risky(.data$term_label_risk) &
        !.data$safe_program_label_in_safe_interpretation
    )
  risk <- enrichment |>
    dplyr::count(
      .data$dataset, .data$claim_type, .data$raw_top_GO_term,
      .data$safe_program_label, .data$term_label_risk,
      .data$label_confidence, .data$FDR_bin, .data$claim_use_class,
      .data$safe_program_label_in_safe_interpretation,
      name = "n_claims"
    ) |>
    dplyr::arrange(.data$dataset, .data$claim_type, .data$term_label_risk, .data$raw_top_GO_term)
  use_summary <- claims |>
    dplyr::mutate(FDR_bin = fdr_bin(.data$FDR)) |>
    dplyr::count(
      .data$dataset, .data$claim_type, .data$term_label_risk,
      .data$label_confidence, .data$FDR_bin, .data$claim_use_class,
      name = "n_claims"
    ) |>
    dplyr::arrange(.data$dataset, .data$claim_type, .data$claim_use_class)
  safe_alignment <- enrichment |>
    dplyr::filter(.data$claim_allowed, go_label_risky(.data$term_label_risk)) |>
    dplyr::group_by(.data$dataset, .data$claim_type, .data$term_label_risk, .data$claim_use_class) |>
    dplyr::summarise(
      n_allowed_go_label_risk_rows = dplyr::n(),
      n_safe_program_label_absent_from_safe_interpretation = sum(.data$risky_allowed_safe_label_missing, na.rm = TRUE),
      examples_missing = paste(utils::head(unique(.data$claim_id[.data$risky_allowed_safe_label_missing]), 8), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_safe_program_label_absent_from_safe_interpretation), .data$dataset, .data$term_label_risk)
  readr::write_csv(risk, file.path(audit_dir, "go_label_risk_audit.csv"), na = "")
  readr::write_csv(use_summary, file.path(audit_dir, "claim_use_class_summary.csv"), na = "")
  readr::write_csv(safe_alignment, file.path(audit_dir, "go_label_safe_interpretation_audit.csv"), na = "")
  invisible(list(
    go_label_risk = risk,
    claim_use_class_summary = use_summary,
    go_label_safe_interpretation = safe_alignment
  ))
}

supermodule_annotation_for_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "module_annotation", dataset, "WGCNA_supermodule_biological_annotation.csv")
  ann <- read_csv_if_exists(f)
  if (is.null(ann) || !nrow(ann)) return(NULL)
  for (col in c("dataset", "SupermoduleID", "safe_display_label", "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "dominant_microenvironment_class", "annotation_confidence", "annotation_stable_across_thresholds", "label_basis", "label_downgrade_reason")) {
    if (!col %in% names(ann)) ann[[col]] <- NA
  }
  ann %>%
    dplyr::mutate(
      dataset = dataset,
      Supermodule_DisplayLabel = dplyr::coalesce(
        as.character(.data$safe_display_label),
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
      dominant_microenvironment_class = .data$dominant_microenvironment_class,
      annotation_confidence = .data$annotation_confidence,
      annotation_stable_across_thresholds = .data$annotation_stable_across_thresholds,
      label_basis = .data$label_basis,
      label_downgrade_reason = .data$label_downgrade_reason
    )
}

latest_csv <- function(root, pattern) latest_file(root, pattern)

collect_program_claims <- function(dataset) {
  f <- resolve_input_path(
    input_name = "biological_program_summary",
    expected_path = path_results("tables", "04_differential_expression_enrichment", "biological_program_summary", dataset, "program_summary.csv"),
    required = TRUE,
    script = SCRIPT_ID,
    dataset = dataset,
    stage = "export",
    producer_script_or_artifact_id = "04_differential_expression_enrichment/03_biological_program_summary.r"
  )
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"biological_program" %in% names(df)) return(empty_claims())
  for (col in c("route_category", "route_unit", "min_raw_p", "min_fdr", "representative_NES", "key_genes", "core_genes", "n_terms", "top_term", "strongest_positive_term", "strongest_negative_term")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df %>%
    dplyr::mutate(
      representative_GO_terms = vapply(seq_len(dplyr::n()), function(i) {
        vals <- unique(stats::na.omit(c(.data$top_term[[i]], .data$strongest_positive_term[[i]], .data$strongest_negative_term[[i]])))
        vals <- vals[nzchar(as.character(vals))]
        paste(vals, collapse = "; ")
      }, character(1))
    ) %>%
    dplyr::transmute(
      dataset = dataset,
      region = .data$route_unit,
      layer_cell_compartment = .data$route_category,
      contrast = .data$comparison,
      direction = .data$direction,
      biological_program = .data$biological_program,
      key_proteins_genes = dplyr::coalesce(.data$key_genes, .data$core_genes),
      evidence_type = "GSEA_program_summary",
      claim_type = "enrichment_program",
      raw_top_GO_term = .data$top_term,
      representative_GO_terms = .data$representative_GO_terms,
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
  f <- resolve_input_path(
    input_name = "wgcna_module_evidence_rank",
    expected_path = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_evidence_rank.csv"),
    fallback_paths = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_priority_summary.csv"),
    required = TRUE,
    script = SCRIPT_ID,
    dataset = dataset,
    stage = "export",
    producer_script_or_artifact_id = "06_modules_WGCNA/01_WGCNA.r"
  )
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
  for (col in c("annotation_confidence", "annotation_stable_across_thresholds", "label_basis", "label_downgrade_reason")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  if (!is.null(super_ann)) {
    df <- df %>%
      dplyr::left_join(super_ann, by = c("dataset", "supermodule_claim_key"), suffix = c("", ".annotation")) %>%
      dplyr::mutate(
        Supermodule_DisplayLabel = dplyr::coalesce(.data$Supermodule_DisplayLabel.annotation, .data$Supermodule_DisplayLabel),
        Supermodule_FinalLabel = dplyr::coalesce(.data$Supermodule_FinalLabel.annotation, .data$Supermodule_FinalLabel),
        Macroprogram_Display = dplyr::coalesce(.data$Macroprogram_Display.annotation, .data$Macroprogram_Display),
        annotation_confidence = dplyr::coalesce(.data$annotation_confidence.annotation, .data$annotation_confidence),
        annotation_stable_across_thresholds = dplyr::coalesce(.data$annotation_stable_across_thresholds.annotation, .data$annotation_stable_across_thresholds),
        label_basis = dplyr::coalesce(.data$label_basis.annotation, .data$label_basis),
        label_downgrade_reason = dplyr::coalesce(.data$label_downgrade_reason.annotation, .data$label_downgrade_reason)
      )
  }
  status <- wgcna_module_claim_status(dataset)
  if (!is.null(status)) {
    df <- df %>%
      dplyr::left_join(status, by = c("dataset", "ModuleID" = "module_id", "strongest_condition_contrast" = "contrast"), suffix = c("", ".group_effect"))
  }
  for (col in c("claim_allowed_model", "animal_level_status", "robustness_gate", "robustness_status", "model_downgrade_reason", "annotation_confidence", "annotation_stable_across_thresholds", "label_basis", "label_downgrade_reason")) {
    if (!col %in% names(df)) df[[col]] <- NA
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
      label_confidence = dplyr::coalesce(.data$annotation_confidence, "not_applicable"),
      label_basis = dplyr::coalesce(.data$label_basis, "not_applicable"),
      label_downgrade_reason = dplyr::coalesce(.data$label_downgrade_reason, "not_applicable"),
      effect_size_NES = .data$strongest_condition_adjusted_delta,
      raw_p = .data$condition_model_p,
      FDR = .data$condition_model_fdr,
      robustness_stability_metric = dplyr::coalesce(
        paste0("robustness_gate=", dplyr::coalesce(.data$robustness_gate, "missing_required"), "; ", dplyr::coalesce(.data$robustness_status, "not_available")),
        if ("preservation_Zsummary_median" %in% names(df)) as.character(.data$preservation_Zsummary_median) else NA_character_
      ),
      key_proteins_genes = paste(.data$ModuleID, .data$ModuleLabel_Final, sep = ": "),
      source_file = f,
      figure_table_target = "WGCNA_module_priority; WGCNA module evidence table",
      interpretation_note = paste0(
        vapply(seq_along(.data$condition_model_fdr), function(i) interpretation_strength(fdr = .data$condition_model_fdr[[i]], effect_size = .data$strongest_condition_adjusted_delta[[i]], n = NA), character(1)),
        "; ", dplyr::coalesce(.data$evidence_warning, "Low-n module evidence; prioritize replicated or convergent modules."),
        "; claim_allowed_model=", dplyr::coalesce(as.character(.data$claim_allowed_model), "FALSE"),
        "; model_downgrade_reason=", dplyr::coalesce(.data$model_downgrade_reason, "not_available"),
        "; annotation_confidence=", dplyr::coalesce(.data$annotation_confidence, "not_available"),
        "; annotation_stable_across_thresholds=", dplyr::coalesce(as.character(.data$annotation_stable_across_thresholds), "not_available")
      ),
      animal_level_status = .data$animal_level_status
    ) %>%
    standardize_claims()
}

collect_overlap_claims <- function(dataset) {
  f <- resolve_input_path(
    input_name = "wgcna_de_gsea_overlap",
    expected_path = path_results("tables", "06_modules_WGCNA", "04_wgcna_de_gsea_overlap", dataset, "WGCNA_vs_DE_GSEA_overlap.csv"),
    required = TRUE,
    script = SCRIPT_ID,
    dataset = dataset,
    stage = "export",
    producer_script_or_artifact_id = "06_modules_WGCNA/04_wgcna_de_gsea_overlap.r"
  )
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
  status <- wgcna_module_claim_status(dataset)
  if (!is.null(status)) {
    df <- df %>%
      dplyr::left_join(status, by = c("dataset", "ModuleID" = "module_id", "contrast" = "contrast"), suffix = c("", ".group_effect"))
  }
  for (col in c("claim_allowed_model", "animal_level_status", "robustness_gate", "robustness_status", "model_downgrade_reason")) {
    if (!col %in% names(df)) df[[col]] <- NA
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
      robustness_stability_metric = paste0("n_DE_overlap=", .data$n_DE_overlap, "; robustness_gate=", dplyr::coalesce(.data$robustness_gate, "missing_required"), "; ", dplyr::coalesce(.data$robustness_status, "not_available")),
      key_proteins_genes = .data$top_overlap_proteins,
      source_file = f,
      figure_table_target = "WGCNA_vs_DE_GSEA_overlap",
      interpretation_note = paste0(
        vapply(seq_along(.data$fisher_fdr), function(i) interpretation_strength(fdr = .data$fisher_fdr[[i]], effect_size = .data$jaccard_DE[[i]], n = .data$n_DE_overlap[[i]]), character(1)),
        "; overlap supports convergence, not causality.",
        "; claim_allowed_model=", dplyr::coalesce(as.character(.data$claim_allowed_model), "FALSE"),
        "; model_downgrade_reason=", dplyr::coalesce(.data$model_downgrade_reason, "not_available")
      ),
      animal_level_status = .data$animal_level_status
    ) %>%
    standardize_claims()
}

collect_behavior_claims <- function() {
  f <- resolve_input_path(
    input_name = "network_behavior_coupling_claims",
    expected_path = path_results("tables", "08_behavior_physio_coupling", "network_behavior_coupling", "edge_behavior_figure_ready_table.csv"),
    required = FALSE,
    script = SCRIPT_ID,
    dataset = "global",
    stage = "export",
    producer_script_or_artifact_id = "08_behavior_physio_coupling/03_module_behavior_coupling.r"
  )
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
  f <- resolve_input_path(
    input_name = "microglia_signature_claims_ready",
    expected_path = path_results(
      "tables",
      "04_differential_expression_enrichment",
      "microglia_targeted_signature_enrichment",
      dataset,
      "microglia_signature_claims_ready.csv"
    ),
    required = identical(dataset, "microglia"),
    script = SCRIPT_ID,
    dataset = dataset,
    stage = "export",
    producer_script_or_artifact_id = "04_differential_expression_enrichment/05_microglia_targeted_signature_enrichment.r"
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

collect_wgcna_group_effect_claims <- function(dataset, level = c("module", "supermodule")) {
  level <- match.arg(level)
  filename <- paste0(level, "_group_effects.csv")
  f <- resolve_input_path(
    input_name = paste0("wgcna_", level, "_group_effects"),
    expected_path = path_results("tables", "06_modules_WGCNA", "group_effects", dataset, filename),
    required = FALSE,
    script = SCRIPT_ID,
    dataset = dataset,
    stage = "export",
    producer_script_or_artifact_id = "06_modules_WGCNA/05_module_supermodule_group_effects.r"
  )
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"contrast" %in% names(df)) return(empty_claims())
  for (col in c(
    "endpoint_id", "endpoint_label", "module_id", "supermodule_id", "module_label", "supermodule_label",
    "spatial_unit", "effect_scope", "estimate", "p_value", "FDR_global", "FDR_within_dataset_level",
    "direction", "evidence_status", "model_family", "model_formula", "primary_model_stable",
    "claim_allowed_model", "model_downgrade_reason", "fallback_used", "fallback_type",
    "animal_level_status", "pseudoreplication_guard", "biological_replicate_unit",
    "n_animals_total", "n_animals_per_group", "min_animals_per_group", "n_samples_total",
    "n_samples_per_group", "model_warning"
  )) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  id_col <- if (identical(level, "module")) "module_id" else "supermodule_id"
  df$module_or_supermodule_id <- dplyr::coalesce(as.character(df[[id_col]]), as.character(df$endpoint_id))
  audit <- read_csv_if_exists(path_results("reviewer_audit", "wgcna_robustness_claim_gate.csv"))
  if (!is.null(audit) && nrow(audit)) {
    audit <- audit %>%
      dplyr::filter(.data$dataset == dataset, .data$level == level) %>%
      dplyr::select("dataset", "level", "module_or_supermodule_id", "contrast", dplyr::any_of(c("spatial_unit", "effect_scope")), "robustness_gate", "preservation_status", "sensitivity_status", "direction_stability", "confounding_status", "robustness_downgrade_reason")
    join_cols <- intersect(c("dataset", "level", "module_or_supermodule_id", "contrast", "spatial_unit", "effect_scope"), names(audit))
    df <- df %>% dplyr::left_join(audit, by = join_cols)
  }
  for (col in c("robustness_gate", "preservation_status", "sensitivity_status", "direction_stability", "confounding_status", "robustness_downgrade_reason")) {
    if (!col %in% names(df)) df[[col]] <- NA_character_
  }
  df %>%
    dplyr::filter(!is.na(.data$contrast), nzchar(as.character(.data$contrast))) %>%
    dplyr::transmute(
      dataset = dataset,
      contrast = .data$contrast,
      biological_program = dplyr::coalesce(blank_to_na(.data$endpoint_label), blank_to_na(if (identical(level, "module")) .data$module_label else .data$supermodule_label), blank_to_na(.data$module_or_supermodule_id)),
      direction = .data$direction,
      key_proteins_genes = .data$module_or_supermodule_id,
      evidence_type = paste0("WGCNA_", level, "_group_effect"),
      effect_size_NES = .data$estimate,
      raw_p = .data$p_value,
      FDR = dplyr::coalesce(.data$FDR_global, .data$FDR_within_dataset_level),
      robustness_stability_metric = paste0(
        "robustness_gate=", dplyr::coalesce(.data$robustness_gate, "missing_required"),
        "; preservation_status=", dplyr::coalesce(.data$preservation_status, "not_available"),
        "; sensitivity_status=", dplyr::coalesce(.data$sensitivity_status, "not_available"),
        "; direction_stability=", dplyr::coalesce(.data$direction_stability, "not_available"),
        "; confounding_status=", dplyr::coalesce(.data$confounding_status, "not_available")
      ),
      source_file = f,
      figure_table_target = paste0("WGCNA_", level, "_group_effects_interpretable; ", filename),
      interpretation_note = paste0(
        "Primary WGCNA group-effect inference; level=", level,
        "; spatial_unit=", .data$spatial_unit,
        "; effect_scope=", .data$effect_scope,
        "; evidence_status=", .data$evidence_status,
        "; model_family=", .data$model_family,
        "; model_formula=", .data$model_formula,
        "; primary_model_stable=", .data$primary_model_stable,
        "; claim_allowed_model=", .data$claim_allowed_model,
        "; model_downgrade_reason=", .data$model_downgrade_reason,
        "; fallback_used=", .data$fallback_used,
        "; fallback_type=", .data$fallback_type,
        "; animal_level_status=", .data$animal_level_status,
        "; pseudoreplication_guard=", .data$pseudoreplication_guard,
        "; biological_replicate_unit=", .data$biological_replicate_unit,
        "; n_animals_total=", .data$n_animals_total,
        "; n_animals_per_group=", .data$n_animals_per_group,
        "; min_animals_per_group=", .data$min_animals_per_group,
        "; n_samples_total=", .data$n_samples_total,
        "; n_samples_per_group=", .data$n_samples_per_group,
        "; robustness_downgrade_reason=", dplyr::coalesce(.data$robustness_downgrade_reason, "not_available"),
        "; model_warning=", .data$model_warning
      ),
      animal_level_status = .data$animal_level_status,
      animal_pseudoreplication_risk = dplyr::case_when(
        .data$animal_level_status %in% c("animal_level", "repeated_sample_mixed_model") ~ "lower",
        TRUE ~ "possible"
      )
    ) %>%
    standardize_claims()
}

collect_integration_claims <- function() {
  f <- resolve_input_path(
    input_name = "manuscript_program_summary",
    expected_path = path_results("tables", "10_biological_integration", "manuscript_program_summary", "global", "manuscript_program_summary.csv"),
    required = FALSE,
    script = SCRIPT_ID,
    dataset = "global",
    stage = "export",
    producer_script_or_artifact_id = "10_biological_integration/02_manuscript_program_summary.r"
  )
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
  dry_run_line("Claim gate evidence availability audit", path_results("reviewer_audit", "claim_gate_evidence_availability.csv"))
  dry_run_line("Claim gate summary audit", path_results("reviewer_audit", "claim_gate_summary.csv"))
  dry_run_line("GO label risk audit", path_results("reviewer_audit", "go_label_risk_audit.csv"))
  dry_run_line("Claim use class summary", path_results("reviewer_audit", "claim_use_class_summary.csv"))
  dry_run_line("GO safe interpretation audit", path_results("reviewer_audit", "go_label_safe_interpretation_audit.csv"))
  dry_run_line("Blocked claim wording audit", path_results("reviewer_audit", "blocked_claim_wording_audit.csv"))
  dry_run_line("WGCNA claim source audit", path_results("reviewer_audit", "wgcna_claim_source_audit.csv"))
  dry_run_line("WGCNA label completeness audit", path_results("reviewer_audit", "wgcna_label_completeness_audit.csv"))
  dry_run_line("Final evidence bundle", path_results("tables", "10_biological_integration", "final_evidence_bundle", "global", "final_biological_evidence_bundle.xlsx"))
  quit(status = 0, save = "no")
}

claims <- dplyr::bind_rows(
  lapply(valid_datasets(), collect_program_claims),
  lapply(valid_datasets(), collect_wgcna_group_effect_claims, level = "module"),
  lapply(valid_datasets(), collect_wgcna_group_effect_claims, level = "supermodule"),
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
    "qc_interpretation_flag"
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
    animal_level_status = dplyr::coalesce(
      blank_to_na(.data$animal_level_status),
      infer_animal_level_status(.data$evidence_type, .data$robustness_stability_metric, .data$source_file)
    ),
    animal_pseudoreplication_risk = dplyr::case_when(
      .data$animal_level_status %in% c("animal_level_or_reported", "animal_level_expected") ~ "lower",
      .data$animal_level_status == "review_source_table" ~ "review_source_table",
      TRUE ~ "possible"
    )
  ) %>%
  add_claim_gates() %>%
  add_go_label_interpretation() %>%
  add_claim_use_class() %>%
  repair_incomplete_wgcna_labels() %>%
  sync_go_safe_interpretation() %>%
  apply_blocked_claim_wording() %>%
  standardize_claims()

validate_table_schema(claims, "biological_claims_table", strict = TRUE)

dir_create(path_results("tables"))
csv_out <- path_results("tables", "biological_claims_table.csv")
xlsx_out <- path_results("tables", "biological_claims_table.xlsx")
readr::write_csv(claims, csv_out, na = "")
if (requireNamespace("writexl", quietly = TRUE)) {
  tryCatch(
    writexl::write_xlsx(list(biological_claims = claims), xlsx_out),
    error = function(e) warning("Could not write biological_claims_table.xlsx: ", conditionMessage(e), call. = FALSE)
  )
}
claim_gate_audits <- write_claim_gate_audits(claims)
go_label_audits <- write_go_label_audits(claims)
blocked_wording_audit <- write_blocked_claim_wording_audit(claims)
wgcna_claim_source_audit <- write_wgcna_claim_source_audit(claims)
wgcna_label_completeness_audit <- write_wgcna_label_completeness_audit(claims)

message("Biological claims table written: ", csv_out)

write_run_manifest(
  path_results("logs", "09_export_pride_journal", "biological_claims_table", "run_manifest.yml"),
  inputs = list(source_files = unique(claims$source_file)),
  outputs = list(
    csv = csv_out,
    xlsx = if (file.exists(xlsx_out)) xlsx_out else NA_character_,
    claim_gate_evidence_availability = path_results("reviewer_audit", "claim_gate_evidence_availability.csv"),
    claim_gate_summary = path_results("reviewer_audit", "claim_gate_summary.csv"),
    go_label_risk_audit = path_results("reviewer_audit", "go_label_risk_audit.csv"),
    claim_use_class_summary = path_results("reviewer_audit", "claim_use_class_summary.csv"),
    go_label_safe_interpretation_audit = path_results("reviewer_audit", "go_label_safe_interpretation_audit.csv"),
    blocked_claim_wording_audit = path_results("reviewer_audit", "blocked_claim_wording_audit.csv"),
    wgcna_claim_source_audit = path_results("reviewer_audit", "wgcna_claim_source_audit.csv"),
    wgcna_label_completeness_audit = path_results("reviewer_audit", "wgcna_label_completeness_audit.csv")
  ),
  parameters = list(datasets = valid_datasets(), schema = "biological_claims_table"),
  notes = "Reviewer-facing manuscript claim gate. claim_grade is descriptive; claim_allowed and claim_gate_status determine eligibility. Missing statistics remain NA."
)

bundle <- write_final_evidence_bundle(reason = "biological_claims_table")
message("Final biological evidence bundle refreshed: ", bundle$bundle_dir)
