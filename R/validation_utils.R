# Lightweight validation and naming helpers shared across active scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

safe_name <- function(x, max_chars = 180) {
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>|]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- ifelse(nzchar(x), x, "unnamed")
  substr(x, 1, max_chars)
}

detect_column <- function(df, candidates, required = FALSE, context = "data frame") {
  nms <- names(df)
  nms_clean <- tolower(gsub("[^a-z0-9]", "", nms))
  cand_clean <- tolower(gsub("[^a-z0-9]", "", candidates))
  idx <- match(cand_clean, nms_clean)
  idx <- idx[!is.na(idx)]
  if (length(idx)) return(nms[idx[[1]]])
  if (isTRUE(required)) {
    stop("Could not find required column in ", context, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  NA_character_
}

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- basename(x)
  x <- sub("\\.d$", "", x, ignore.case = TRUE)
  x <- trimws(tolower(x))
  x <- gsub("[[:space:]]+", "_", x)
  x
}

sample_overlap_summary <- function(matrix_samples, metadata_samples) {
  matrix_norm <- normalize_sample_id(matrix_samples)
  metadata_norm <- normalize_sample_id(metadata_samples)
  overlap <- intersect(matrix_norm, metadata_norm)
  data.frame(
    n_matrix_samples = length(unique(matrix_norm[!is.na(matrix_norm) & nzchar(matrix_norm)])),
    n_metadata_samples = length(unique(metadata_norm[!is.na(metadata_norm) & nzchar(metadata_norm)])),
    n_overlap = length(unique(overlap)),
    overlap_fraction_matrix = if (length(unique(matrix_norm))) length(unique(overlap)) / length(unique(matrix_norm)) else NA_real_,
    overlap_fraction_metadata = if (length(unique(metadata_norm))) length(unique(overlap)) / length(unique(metadata_norm)) else NA_real_,
    stringsAsFactors = FALSE
  )
}

validate_manifest_paths <- function(manifest, path_cols = c("input_gene_file", "output_table", "output_plot"), allow_missing = TRUE) {
  if (is.null(manifest) || !nrow(manifest)) {
    return(data.frame(path_column = character(), path = character(), exists = logical(), stringsAsFactors = FALSE))
  }
  cols <- intersect(path_cols, names(manifest))
  out <- do.call(rbind, lapply(cols, function(col) {
    vals <- unique(as.character(manifest[[col]]))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    data.frame(path_column = col, path = vals, exists = file.exists(vals), stringsAsFactors = FALSE)
  }))
  if (is.null(out)) out <- data.frame(path_column = character(), path = character(), exists = logical(), stringsAsFactors = FALSE)
  if (!allow_missing && any(!out$exists)) {
    stop("Manifest contains missing paths:\n", paste(out$path[!out$exists], collapse = "\n"), call. = FALSE)
  }
  out
}

duplicate_key_summary <- function(df, keys) {
  keys <- intersect(keys, names(df))
  if (!length(keys) || is.null(df) || !nrow(df)) {
    return(data.frame(n_duplicate_keys = 0L, stringsAsFactors = FALSE))
  }
  key <- do.call(paste, c(df[keys], sep = "\r"))
  data.frame(n_duplicate_keys = sum(duplicated(key)), stringsAsFactors = FALSE)
}

interpretation_strength <- function(fdr = NA_real_, effect_size = NA_real_, n = NA_integer_, bootstrap_support = NA_real_) {
  fdr <- suppressWarnings(as.numeric(fdr))
  effect_size <- suppressWarnings(abs(as.numeric(effect_size)))
  n <- suppressWarnings(as.numeric(n))
  bootstrap_support <- suppressWarnings(as.numeric(bootstrap_support))
  if (!is.na(fdr) && fdr <= 0.05 && (is.na(effect_size) || effect_size >= 0.3) && (is.na(n) || n >= 6) && (is.na(bootstrap_support) || bootstrap_support >= 0.7)) {
    return("strong")
  }
  if (!is.na(fdr) && fdr <= 0.10 && (is.na(n) || n >= 6)) return("moderate")
  "exploratory"
}

known_pipeline_output_specs <- function() {
  group_effect_cols <- c(
    "dataset", "level", "endpoint_id", "endpoint_label", "contrast",
    "estimate", "SE", "p_value", "FDR_within_dataset_level", "FDR_global",
    "evidence_status", "n_samples", "n_animals", "n_animals_total",
    "n_animals_per_group", "min_animals_per_group", "n_samples_total",
    "n_samples_per_group", "animal_level_status", "pseudoreplication_guard",
    "model_type", "model_family", "model_formula", "primary_model_stable",
    "claim_allowed_model", "model_downgrade_reason", "fallback_used",
    "fallback_type", "formula_used", "rank_deficient_model", "singular_model",
    "emmeans_success", "animal_random_effect_used", "biological_replicate_unit",
    "model_warning"
  )
  list(
    module_group_effects.csv = list(required_columns = group_effect_cols),
    supermodule_group_effects.csv = list(required_columns = group_effect_cols),
    WGCNA_parameter_audit.csv = list(
      required_columns = c(
        "dataset", "n_samples", "n_animals", "n_features_input",
        "n_features_after_filter", "missingness_filter", "imputation_source",
        "normalization_source", "correlation_method", "network_type", "TOM_type",
        "soft_power", "scale_free_R2", "mean_connectivity", "min_module_size",
        "merge_cut_height", "deepSplit", "pamRespectsDendro", "random_seed",
        "supermodule_cut_height", "input_hashes", "provenance"
      )
    ),
    wgcna_robustness_claim_gate.csv = list(
      required_columns = c(
        "dataset", "module_or_supermodule_id", "level", "contrast",
        "robustness_gate", "preservation_status", "sensitivity_status",
        "direction_stability", "confounding_status", "claim_gate_eligible",
        "robustness_downgrade_reason"
      )
    ),
    microglia_neuropil_independence_claim_gate.csv = list(
      required_columns = c(
        "module_or_supermodule_id", "endpoint_type", "endpoint_scope", "source_level",
        "direct_independence_tested", "direct_supermodule_test", "contrast", "biological_program",
        "microenvironment_class", "adjustment_mode", "covariate_family",
        "primary_effect_status", "primary_effect_claim_relevant",
        "primary_effect_threshold",
        "independence_classification", "claim_gate_eligible", "downgrade_reason",
        "n_matched_animals", "min_animals_per_group", "effect_before",
        "effect_after", "effect_before_abs", "effect_before_near_zero",
        "percent_attenuation", "percent_attenuation_reliable",
        "direction_preserved", "audit_group_n",
        "eligible_without_primary_effect_count"
      )
    ),
    microglia_neuropil_covariate_selection_audit.csv = list(
      required_columns = c(
        "candidate_covariates", "predeclared_covariates_available",
        "selected_primary_covariate", "selected_secondary_covariate",
        "selected_exploratory_covariate", "selection_rule",
        "primary_claim_gate_eligible", "secondary_claim_gate_eligible",
        "exploratory_claim_gate_eligible"
      )
    ),
    microglia_neuropil_independence_endpoint_scope_audit.csv = list(
      required_columns = c(
        "endpoint_type", "endpoint_scope", "source_level", "n_endpoints",
        "direct_independence_tested", "consumed_by_claim_type", "notes"
      )
    ),
    claim_use_class_wording_audit.csv = list(
      required_columns = c("check_name", "severe_issue_count", "status", "example_claim_ids", "notes")
    ),
    final_claim_gate_summary.csv = list(
      required_columns = c(
        "dataset", "claim_type", "evidence_type", "claim_allowed", "claim_gate_status",
        "claim_use_class", "primary_model_status", "animal_level_gate", "microglia_roi_gate",
        "neuropil_independence_gate", "robustness_gate", "evidence_independence_gate", "n_claims"
      )
    ),
    final_reviewer_audit_manifest.csv = list(
      required_columns = c("audit_file", "exists", "n_rows", "schema_validated", "produced_by_script", "reviewer_use", "manuscript_use_allowed", "notes")
    ),
    final_evidence_bundle_validation.csv = list(
      required_columns = c("validation_check", "status", "n_violations", "details")
    ),
    WGCNA_label_candidates.csv = list(
      required_columns = c(
        "dataset", "level", "entity_id", "candidate_label", "candidate_source",
        "evidence_strength", "hub_support", "marker_context_support", "genericity_penalty",
        "ontology_mismatch_flag", "conflict_penalty", "final_label_score", "selected_label", "rejection_reason"
      )
    ),
    WGCNA_final_label_lookup.csv = list(
      required_columns = c(
        "dataset", "level", "entity_id", "parent_entity_id", "n_member_modules",
        "raw_top_GO_label", "best_data_driven_label", "parent_program", "microenvironment_context",
        "final_plot_label", "label_confidence", "label_basis", "label_score",
        "manual_review_required", "ontology_mismatch_flag", "label_rationale", "unsafe_interpretation"
      )
    ),
    wgcna_microenvironment_threshold_sensitivity.csv = list(
      required_columns = c(
        "dataset", "module_or_supermodule_id", "level",
        "class_at_0.05", "class_at_0.10", "class_at_0.20",
        "primary_class", "annotation_stable_across_thresholds",
        "classification_flip_reason", "marker_fraction_primary",
        "marker_panels_supporting", "claim_relevance"
      )
    ),
    wgcna_label_confidence_audit.csv = list(
      required_columns = c(
        "dataset", "module_or_supermodule_id", "level",
        "raw_annotation_label", "cleaned_annotation_label", "safe_display_label",
        "label_confidence", "label_basis", "label_downgrade_reason",
        "annotation_stable_across_thresholds", "unsafe_interpretation"
      )
    ),
    wgcna_annotation_source_audit.csv = list(
      required_columns = c(
        "dataset", "panel_version", "panel_id", "source_type", "source_reference",
        "allowed_use", "claim_role", "caution_note", "n_markers",
        "config_file", "config_hash"
      )
    ),
    WGCNA_module_biological_annotation.csv = list(
      required_columns = c(
        "dataset", "ModuleID", "ModuleColor", "microenvironment_class",
        "microenvironment_label", "classification_rationale", "interpretation_note"
      ),
      recommended_columns = c(
        "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms",
        "raw_top_GO_label", "raw_module_label", "raw_hub_proteins",
        "raw_marker_or_signature_label", "cleaned_biological_label",
        "cleaned_biological_label_short", "cleaned_biological_label_source",
        "cleaned_biological_label_confidence", "GO_label_relevance_flag",
        "GO_label_relevance_rationale", "microenvironment_caution_label",
        "microenvironment_caution_class", "microenvironment_caution_rationale",
        "Module_CleanPlotLabel", "module_biological_label",
        "module_biological_label_short", "module_label_display",
        "annotation_confidence", "annotation_basis", "annotation_downgrade_reason",
        "annotation_stable_across_thresholds", "unsafe_interpretation",
        "raw_annotation_label", "cleaned_annotation_label", "safe_display_label",
        "label_confidence", "label_basis", "label_downgrade_reason",
        "marker_fraction_primary", "marker_panels_supporting"
      )
    ),
    WGCNA_supermodule_biological_annotation.csv = list(
      required_columns = c(
        "dataset", "SupermoduleID", "Supermodule_DisplayLabel",
        "dominant_microenvironment_class", "dominant_module_labels",
        "Supermodule_LabelRationale", "interpretation_note"
      ),
      recommended_columns = c(
        "Supermodule_ConservativeLabel", "Supermodule_CompositionLabel",
        "Supermodule_CompositionShortLabel", "Supermodule_CompositionDisplayLabel",
        "Supermodule_CleanPlotLabel", "Supermodule_CompositionLabelSource",
        "Supermodule_CompositionConfidence", "Supermodule_CompositionRationale",
        "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms",
        "raw_top_GO_label", "raw_module_label", "raw_hub_proteins",
        "raw_marker_or_signature_label", "cleaned_biological_label",
        "cleaned_biological_label_short", "cleaned_biological_label_source",
        "cleaned_biological_label_confidence", "GO_label_relevance_flag",
        "GO_label_relevance_rationale", "microenvironment_caution_label",
        "microenvironment_caution_class", "microenvironment_caution_rationale",
        "DominantMemberTheme", "DominantMemberThemeFraction", "SecondMemberTheme",
        "SecondMemberThemeFraction", "TopMemberModuleLabels", "TopMemberGOTerms",
        "MemberThemeCounts", "MemberThemeFractions", "n_distinct_member_themes",
        "is_multi_theme_supermodule", "themes_above_display_threshold",
        "themes_omitted_from_display_label",
        "n_member_modules_with_informative_labels",
        "fraction_member_modules_with_informative_labels",
        "annotation_confidence", "annotation_basis", "annotation_downgrade_reason",
        "annotation_stable_across_thresholds", "unsafe_interpretation",
        "raw_annotation_label", "cleaned_annotation_label", "safe_display_label",
        "label_confidence", "label_basis", "label_downgrade_reason",
        "marker_fraction_primary", "marker_panels_supporting"
      )
    ),
    WGCNA_module_group_effects_interpretable.csv = list(
      required_columns = c(
        "dataset", "level", "contrast", "estimate", "p_value", "FDR_global",
        "interpretation_sentence", "ModulePlotLabel", "Supermodule_PlotLabel",
        "Supermodule_FullAnnotationLabel"
      ),
      recommended_columns = c(
        "module_label", "module_biological_label", "module_biological_label_short",
        "module_label_display", "supermodule_id", "Module_CleanPlotLabel",
        "Supermodule_CleanPlotLabel", "Supermodule_CompositionDisplayLabel",
        "Supermodule_CompositionLabel", "cleaned_biological_label",
        "microenvironment_caution_label", "GO_label_relevance_flag",
        "GO_label_relevance_rationale"
      )
    ),
    WGCNA_supermodule_group_effects_interpretable.csv = list(
      required_columns = c(
        "dataset", "level", "contrast", "estimate", "p_value", "FDR_global",
        "interpretation_sentence", "Supermodule_PlotLabel",
        "Supermodule_FullAnnotationLabel", "dominant_microenvironment_class",
        "Supermodule_LabelRationale"
      ),
      recommended_columns = c(
        "Supermodule_CleanPlotLabel", "Supermodule_CompositionDisplayLabel",
        "Supermodule_CompositionLabel", "cleaned_biological_label",
        "microenvironment_caution_label", "GO_label_relevance_flag",
        "GO_label_relevance_rationale", "TopMemberGOTerms",
        "MemberThemeCounts", "MemberThemeFractions", "n_distinct_member_themes",
        "is_multi_theme_supermodule", "themes_above_display_threshold",
        "themes_omitted_from_display_label", "supermodule_theme_label_qc_warning"
      )
    ),
    biological_claims.csv = list(
      required_columns = c(
        "claim_id", "dataset", "biological_program", "evidence_type", "claim_type",
        "claim_use_class", "raw_top_GO_term", "representative_GO_terms",
        "semantic_parent_label", "safe_program_label", "term_label_risk",
        "label_confidence", "label_basis", "label_downgrade_reason",
        "claim_grade", "primary_evidence", "orthogonal_support",
        "major_limitation", "safe_interpretation", "unsafe_overinterpretation",
        "claim_allowed", "claim_gate_status", "claim_downgrade_reason",
        "model_fit_status", "statistical_evidence_status", "claim_gate_model_status",
        "primary_model_status", "animal_level_gate", "qc_gate",
        "missingness_gate", "batch_confound_gate", "marker_contamination_gate",
        "microglia_roi_gate", "neuropil_independence_gate", "robustness_gate",
        "evidence_independence_gate",
        "missingness_confounded", "plate_or_batch_confounded",
        "region_layer_imbalance_risk", "animal_pseudoreplication_risk",
        "marker_contamination_risk", "qc_interpretation_flag"
      )
    )
  )
}

validation_displayed_supermodule_theme_count <- function(label) {
  vapply(as.character(label), function(z) {
    z <- trimws(gsub("\\s+", " ", z))
    if (is.na(z) || !nzchar(z)) return(0L)
    z <- sub("^SM[0-9]+\\s*(Â·|·|:|-)\\s*", "", z, ignore.case = TRUE)
    z <- sub("^[^:]{1,80}:\\s*(?=(mostly\\s+|mixed:|mixed\\s+multi-program|mixed /))", "", z, perl = TRUE, ignore.case = TRUE)
    if (grepl("^mixed\\s+multi-program$", z, ignore.case = TRUE)) return(0L)
    if (grepl("^mixed:\\s*", z, ignore.case = TRUE)) {
      parts <- trimws(unlist(strsplit(sub("^mixed:\\s*", "", z, ignore.case = TRUE), "\\s*;\\s*", perl = TRUE), use.names = FALSE))
      parts <- parts[nzchar(parts) & !grepl("mixed / low-specificity|mixed / unresolved|unresolved / mixed", parts, ignore.case = TRUE)]
      return(max(2L, length(unique(parts))))
    }
    z <- sub("^mostly\\s+", "", z, ignore.case = TRUE)
    z <- sub("^mixed:\\s*", "", z, ignore.case = TRUE)
    parts <- trimws(unlist(strsplit(z, "\\s*;\\s*", perl = TRUE), use.names = FALSE))
    parts <- parts[nzchar(parts) & !grepl("mixed / low-specificity|mixed / unresolved|unresolved / mixed", parts, ignore.case = TRUE)]
    length(unique(parts))
  }, integer(1))
}

validation_member_theme_fraction_values <- function(x) {
  parts <- trimws(unlist(strsplit(as.character(x), "\\s*;\\s*", perl = TRUE), use.names = FALSE))
  parts <- parts[nzchar(parts) & grepl("=", parts, fixed = TRUE)]
  vals <- suppressWarnings(as.numeric(sub("^.*=", "", parts)))
  names(vals) <- trimws(sub("=.*$", "", parts))
  vals[is.finite(vals)]
}

append_supermodule_theme_audit_messages <- function(messages, df, context = "Supermodule") {
  if (!nrow(df)) return(messages)
  if ("MemberThemeFractions" %in% names(df)) {
    empty_frac <- is.na(df$MemberThemeFractions) | !nzchar(trimws(as.character(df$MemberThemeFractions)))
    if (any(empty_frac, na.rm = TRUE)) {
      messages <- c(messages, paste0(context, " MemberThemeFractions is empty for one or more rows."))
    }
  }
  if ("themes_omitted_from_display_label" %in% names(df)) {
    omitted <- !is.na(df$themes_omitted_from_display_label) & nzchar(trimws(as.character(df$themes_omitted_from_display_label)))
    if (any(omitted, na.rm = TRUE)) {
      messages <- c(messages, paste0(context, " themes_omitted_from_display_label is non-empty for one or more rows."))
    }
  }
  if (all(c("n_distinct_member_themes", "MemberThemeFractions", "Supermodule_PlotLabel") %in% names(df))) {
    n_distinct <- suppressWarnings(as.integer(df$n_distinct_member_themes))
    shown_n <- validation_displayed_supermodule_theme_count(df$Supermodule_PlotLabel)
    hidden_three <- vapply(seq_len(nrow(df)), function(i) {
      vals <- validation_member_theme_fraction_values(df$MemberThemeFractions[[i]])
      vals <- vals[names(vals) != "mixed / low-specificity"]
      is.finite(n_distinct[[i]]) && n_distinct[[i]] >= 3L && length(vals) >= 3L && all(vals >= 0.20) && shown_n[[i]] < 3L
    }, logical(1))
    if (any(hidden_three, na.rm = TRUE)) {
      messages <- c(messages, paste0(context, " plot label shows fewer than 3 themes despite >=3 member themes all at fraction >=0.20."))
    }
  }
  if (all(c("DominantMemberThemeFraction", "Supermodule_PlotLabel") %in% names(df))) {
    dom <- suppressWarnings(as.numeric(df$DominantMemberThemeFraction))
    shown_n <- validation_displayed_supermodule_theme_count(df$Supermodule_PlotLabel)
    single_theme <- shown_n <= 1L &
      !grepl("^SM[0-9]+\\s*(Â·|·|:|-)\\s*(mixed:|mixed multi-program|mixed / low-specificity|mixed / unresolved)", as.character(df$Supermodule_PlotLabel), ignore.case = TRUE)
    if (any((!is.finite(dom) | dom < 0.60) & single_theme, na.rm = TRUE)) {
      messages <- c(messages, paste0(context, " without dominant theme >=0.60 is displayed as a single-theme program."))
    }
  }
  messages
}

validate_known_pipeline_output <- function(path, dataset = NULL) {
  specs <- known_pipeline_output_specs()
  filename <- basename(path)
  if (!filename %in% names(specs)) {
    return(data.frame(path = path, validation_status = "not_applicable", validation_message = "", stringsAsFactors = FALSE))
  }
  if (!file.exists(path)) {
    return(data.frame(path = path, validation_status = "missing", validation_message = "Expected output file does not exist.", stringsAsFactors = FALSE))
  }

  df <- tryCatch(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) e)
  if (inherits(df, "error")) {
    return(data.frame(path = path, validation_status = "error", validation_message = conditionMessage(df), stringsAsFactors = FALSE))
  }

  messages <- character()
  required <- specs[[filename]]$required_columns
  missing <- setdiff(required, names(df))
  if (length(missing)) messages <- c(messages, paste0("Missing required column(s): ", paste(missing, collapse = ", ")))
  recommended <- specs[[filename]]$recommended_columns %||% character()
  missing_recommended <- setdiff(recommended, names(df))
  if (length(missing_recommended)) messages <- c(messages, paste0("Missing recommended optional column(s): ", paste(missing_recommended, collapse = ", ")))

  if ("dataset" %in% names(df) && exists("valid_datasets", mode = "function")) {
    observed <- unique(as.character(df$dataset[!is.na(df$dataset) & nzchar(as.character(df$dataset))]))
    invalid <- setdiff(observed, valid_datasets())
    if (length(invalid)) messages <- c(messages, paste0("Invalid dataset value(s): ", paste(invalid, collapse = ", ")))
  }
  if (!is.null(dataset) && "dataset" %in% names(df)) {
    observed <- unique(as.character(df$dataset[!is.na(df$dataset) & nzchar(as.character(df$dataset))]))
    mismatch <- setdiff(observed, dataset)
    if (length(mismatch)) messages <- c(messages, paste0("Dataset values do not match selected dataset ", dataset, ": ", paste(mismatch, collapse = ", ")))
  }

  for (col in intersect(c("p_value", "FDR_within_dataset_level", "FDR_global"), names(df))) {
    values <- suppressWarnings(as.numeric(df[[col]]))
    bad <- !is.na(values) & (values < 0 | values > 1)
    if (any(bad)) messages <- c(messages, paste0(col, " has value(s) outside [0, 1]."))
  }
  for (col in intersect(c("n_samples", "n_animals"), names(df))) {
    values <- suppressWarnings(as.numeric(df[[col]]))
    bad <- !is.na(values) & values < 0
    if (any(bad)) messages <- c(messages, paste0(col, " has negative value(s)."))
  }

  if (identical(filename, "WGCNA_supermodule_biological_annotation.csv")) {
    if (all(c("dataset", "microenvironment_caution_label") %in% names(df))) {
      ds_vals <- unique(as.character(df$dataset[!is.na(df$dataset) & nzchar(as.character(df$dataset))]))
      non_micro <- ds_vals[ds_vals %in% c("neuron_soma", "neuron_neuropil")]
      if (length(non_micro) && any(grepl("shared microglia-neuropil ROI", as.character(df$microenvironment_caution_label), ignore.case = TRUE), na.rm = TRUE)) {
        messages <- c(messages, "Neuron dataset outputs must not contain shared microglia-neuropil ROI caution labels.")
      }
    }
    if (all(c("Supermodule_CleanPlotLabel", "Supermodule_CompositionLabel") %in% names(df))) {
      clean_available <- !is.na(df$Supermodule_CleanPlotLabel) & nzchar(as.character(df$Supermodule_CleanPlotLabel)) |
        !is.na(df$Supermodule_CompositionLabel) & nzchar(as.character(df$Supermodule_CompositionLabel))
      if (nrow(df) && mean(clean_available, na.rm = TRUE) < 0.80) {
        messages <- c(messages, "Supermodule_CleanPlotLabel or Supermodule_CompositionLabel is available for <=80% of supermodules.")
      }
    }
    if ("TopMemberGOTerms" %in% names(df) && all(is.na(df$TopMemberGOTerms) | !nzchar(as.character(df$TopMemberGOTerms)))) {
      messages <- c(messages, "TopMemberGOTerms is entirely empty; check WGCNA GO parsing if enrichment input exists.")
    }
    if (all(c("Supermodule_PlotLabel", "Supermodule_CleanPlotLabel", "Supermodule_CompositionLabel") %in% names(df))) {
      plot_mixed <- grepl("^SM[0-9]+\\s*(·|:|-)\\s*(mixed / unresolved|mixed / low-specificity)$", as.character(df$Supermodule_PlotLabel), ignore.case = TRUE)
      has_clean <- !is.na(df$Supermodule_CleanPlotLabel) & nzchar(as.character(df$Supermodule_CleanPlotLabel)) |
        (!is.na(df$Supermodule_CompositionLabel) & nzchar(as.character(df$Supermodule_CompositionLabel)) & !grepl("mixed / low-specificity|mixed / unresolved", as.character(df$Supermodule_CompositionLabel), ignore.case = TRUE))
      if (any(plot_mixed & has_clean, na.rm = TRUE)) {
        messages <- c(messages, "Some Supermodule_PlotLabel values are only SMxx mixed despite available cleaned/composition labels.")
      }
    }
    if ("Supermodule_CompositionLabel" %in% names(df)) {
      comp <- trimws(as.character(df$Supermodule_CompositionLabel))
      comp_present <- !is.na(comp) & nzchar(comp)
      unresolved <- comp_present & grepl("mixed|unresolved|low-specificity", comp, ignore.case = TRUE)
      if (any(comp_present) && all(unresolved[comp_present])) {
        messages <- c(messages, "All supermodule composition labels are mixed/unresolved; member-module biological labels may not be propagating.")
      }
    } else {
      label_cols <- intersect(c("Supermodule_FinalLabel", "Supermodule_DisplayLabel", "dominant_module_labels"), names(df))
      if (length(label_cols)) {
        lab <- apply(df[label_cols], 1, paste, collapse = " ")
        if (length(lab) && all(grepl("mixed|unresolved|low-specificity", lab, ignore.case = TRUE))) {
          messages <- c(messages, "All supermodules are mixed/unresolved and no composition labels are available.")
        }
      }
    }
    if ("Supermodule_CompositionLabel" %in% names(df) &&
        "Supermodule_CompositionConfidence" %in% names(df)) {
      comp <- trimws(as.character(df$Supermodule_CompositionLabel))
      conf <- trimws(as.character(df$Supermodule_CompositionConfidence))
      if (any(!is.na(comp) & nzchar(comp)) && !any(!is.na(conf) & nzchar(conf))) {
        messages <- c(messages, "Composition labels are present but all confidence values are missing.")
      }
    }
    if (all(c("dataset", "dominant_microenvironment_class", "Supermodule_CompositionLabel") %in% names(df))) {
      micro_rows <- as.character(df$dataset) == "microglia"
      caution_dom <- as.character(df$dominant_microenvironment_class) %in% c(
        "shared_microenvironment", "neuropil_sensitive",
        "vascular_basement_membrane_ecm", "vascular_bbb_mural"
      )
      overclaim <- grepl("microglia activation|microglia state", as.character(df$Supermodule_CompositionLabel), ignore.case = TRUE)
      if (any(micro_rows & caution_dom & overclaim, na.rm = TRUE)) {
        messages <- c(messages, "Microglia composition labels overclaim microglia specificity despite shared/neuropil/vascular dominance.")
      }
      conflated <- micro_rows & grepl("shared microglia-neuropil ROI|neuropil-sensitive ROI|perivascular/ECM ROI|vascular/BBB ROI", as.character(df$Supermodule_CompositionLabel), ignore.case = TRUE)
      if (any(conflated, na.rm = TRUE)) {
        messages <- c(messages, "Microglia composition labels appear to contain ROI caution wording; composition and caution should be separate.")
      }
    }
    if (all(c("GO_label_relevance_flag", "cleaned_biological_label", "raw_module_label") %in% names(df))) {
      suspicious_raw <- grepl("skin development|dorsal ventral pattern|binding sperm|zona pellucida|fertilization", as.character(df$raw_module_label), ignore.case = TRUE)
      not_cleaned <- suspicious_raw & !as.character(df$GO_label_relevance_flag) %in% c("ontology_context_mismatch", "cleaned_or_context_checked", "cleaned_go_display")
      if (any(not_cleaned, na.rm = TRUE)) {
        messages <- c(messages, "Suspicious GO labels are present but not flagged as cleaned/context-checked.")
      }
    }
    messages <- append_supermodule_theme_audit_messages(messages, df, context = "Supermodule annotation")
  }

  if (identical(filename, "WGCNA_module_group_effects_interpretable.csv")) {
    if ("module_label" %in% names(df) && all(is.na(df$module_label) | !nzchar(as.character(df$module_label)))) {
      messages <- c(messages, "Canonical module_label is entirely empty.")
    }
    if ("supermodule_id" %in% names(df) && all(is.na(df$supermodule_id) | !nzchar(as.character(df$supermodule_id)))) {
      messages <- c(messages, "Canonical module-level supermodule_id is entirely empty.")
    }
    if (all(c("cleaned_biological_label", "module_label_display", "Module_CleanPlotLabel") %in% names(df))) {
      has_clean <- !is.na(df$cleaned_biological_label) & nzchar(as.character(df$cleaned_biological_label))
      has_display <- (!is.na(df$module_label_display) & nzchar(as.character(df$module_label_display))) |
        (!is.na(df$Module_CleanPlotLabel) & nzchar(as.character(df$Module_CleanPlotLabel)))
      if (any(has_clean & !has_display, na.rm = TRUE)) {
        messages <- c(messages, "module_label_display or Module_CleanPlotLabel is empty despite cleaned_biological_label.")
      }
      mixed_user_label <- has_clean & grepl("mixed / low-specificity$", as.character(df$module_label_display), ignore.case = TRUE)
      if (any(mixed_user_label, na.rm = TRUE)) {
        messages <- c(messages, "User-facing module label is only mixed / low-specificity despite cleaned_biological_label.")
      }
    }
    if (all(c("dataset", "microenvironment_caution_label") %in% names(df))) {
      non_micro_rows <- as.character(df$dataset) %in% c("neuron_soma", "neuron_neuropil")
      if (any(non_micro_rows & grepl("shared microglia-neuropil ROI", as.character(df$microenvironment_caution_label), ignore.case = TRUE), na.rm = TRUE)) {
        messages <- c(messages, "Neuron module outputs must not contain shared microglia-neuropil ROI caution labels.")
      }
    }
    if (all(c("Supermodule_PlotLabel", "Supermodule_CleanPlotLabel", "Supermodule_CompositionLabel") %in% names(df))) {
      plot_mixed <- grepl("^SM[0-9]+\\s*(·|:|-)\\s*(mixed / unresolved|mixed / low-specificity)$", as.character(df$Supermodule_PlotLabel), ignore.case = TRUE)
      has_clean <- !is.na(df$Supermodule_CleanPlotLabel) & nzchar(as.character(df$Supermodule_CleanPlotLabel)) |
        (!is.na(df$Supermodule_CompositionLabel) & nzchar(as.character(df$Supermodule_CompositionLabel)) & !grepl("mixed / low-specificity|mixed / unresolved", as.character(df$Supermodule_CompositionLabel), ignore.case = TRUE))
      if (any(plot_mixed & has_clean, na.rm = TRUE)) {
        messages <- c(messages, "Module-level Supermodule_PlotLabel values include SMxx mixed despite available cleaned/composition labels.")
      }
    }
  }

  if (identical(filename, "WGCNA_supermodule_group_effects_interpretable.csv")) {
    if (all(c("supermodule_id", "SupermoduleID") %in% names(df))) {
      missing_super_id <- (!is.na(df$supermodule_id) & nzchar(as.character(df$supermodule_id))) &
        (is.na(df$SupermoduleID) | !nzchar(as.character(df$SupermoduleID)))
      if (any(missing_super_id, na.rm = TRUE)) {
        messages <- c(messages, "SupermoduleID is empty where supermodule_id exists.")
      }
    }
    if (all(c("dataset", "microenvironment_caution_label") %in% names(df))) {
      non_micro_rows <- as.character(df$dataset) %in% c("neuron_soma", "neuron_neuropil")
      if (any(non_micro_rows & grepl("shared microglia-neuropil ROI", as.character(df$microenvironment_caution_label), ignore.case = TRUE), na.rm = TRUE)) {
        messages <- c(messages, "Neuron supermodule outputs must not contain shared microglia-neuropil ROI caution labels.")
      }
    }
    messages <- append_supermodule_theme_audit_messages(messages, df, context = "Supermodule interpretable output")
  }

  data.frame(
    path = path,
    validation_status = if (length(messages)) "warning" else "ok",
    validation_message = paste(messages, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

validate_pipeline_outputs <- function(paths, dataset = NULL) {
  paths <- unique(paths[file.exists(paths) & !dir.exists(paths)])
  if (!length(paths)) {
    return(data.frame(path = character(), validation_status = character(), validation_message = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(paths, validate_known_pipeline_output, dataset = dataset))
}
