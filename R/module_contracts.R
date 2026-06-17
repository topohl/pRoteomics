# Shared module contracts and identifier helpers.

normalize_module_identifier <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(";.*$", "", x)
  x <- gsub("\\|.*$", "", x)
  x[is.na(x)] <- ""
  x
}

normalize_gene_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x)] <- ""
  x
}

require_module_contract_columns <- function(df, cols, artifact = "artifact") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing)) {
    stop(
      artifact, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_wgcna_module_definitions <- function(df, artifact = "WGCNA module definitions") {
  require_module_contract_columns(
    df,
    c("ModuleSet", "ModuleID", "ModuleColor", "ProteinID", "UniProt", "GeneSymbol"),
    artifact
  )
  if (!any(c("kME", "Weight") %in% colnames(df))) {
    stop(artifact, " must contain kME or Weight.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_curated_overlap_programs <- function(df, artifact = "curated overlap programs") {
  require_module_contract_columns(
    df,
    c("ModuleSet", "ModuleID", "UniProt", "GeneSymbol", "Source"),
    artifact
  )
  invisible(TRUE)
}

validate_module_score_output <- function(df, artifact = "module score output") {
  require_module_contract_columns(
    df,
    c(
      "dataset", "module_definition_source", "ModuleID",
      "ModuleScore", "ScoreType", "n_found_in_matrix", "coverage_fraction"
    ),
    artifact
  )
  if (!any(c("Sample", "AnimalID") %in% colnames(df))) {
    stop(
      artifact, " is missing required key column(s): expected Sample or AnimalID.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_wgcna_group_effects <- function(df, artifact = "WGCNA group effects") {
  require_module_contract_columns(
    df,
    c(
      "dataset", "level", "spatial_unit", "contrast", "estimate", "SE",
      "p_value", "FDR_within_dataset_level", "FDR_global", "direction",
      "effect_scope", "SpatialUnitType", "model_type", "has_repeated_animals",
      "n_animals", "n_samples", "formula_requested", "formula_used", "dropped_covariates",
      "rank_deficient_model", "model_warning"
    ),
    artifact
  )
  if (!any(c("module_id", "supermodule_id") %in% colnames(df))) {
    stop(artifact, " must contain module_id or supermodule_id.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_wgcna_module_annotation <- function(df, artifact = "WGCNA module biological annotation") {
  require_module_contract_columns(
    df,
    c(
      "dataset", "ModuleID", "ModuleColor", "n_proteins", "microenvironment_class",
      "microglia_evidence", "neuropil_evidence", "other_cellular_evidence",
      "canonical_microglia_evidence", "empirical_microglia_roi_evidence",
      "canonical_neuropil_evidence", "empirical_neuropil_evidence",
      "empirical_shared_microenvironment_evidence",
      "microglia_state_or_activation_evidence",
      "peripheral_myeloid_caution_evidence",
      "classification_threshold", "classification_rationale",
      "marker_registry_version", "empirical_marker_set_version",
      "interpretation_note"
    ),
    artifact
  )
  invisible(TRUE)
}

validate_wgcna_interpretable_summary <- function(df, artifact = "WGCNA interpretable summary") {
  require_module_contract_columns(
    df,
    c("dataset", "level", "contrast", "estimate", "p_value", "FDR_global", "interpretation_sentence"),
    artifact
  )
  invisible(TRUE)
}

validate_cross_compartment_program_atlas <- function(df, artifact = "cross-compartment program atlas") {
  require_module_contract_columns(
    df,
    c(
      "dataset", "evidence_domain", "evidence_id", "program_label",
      "entity_type", "entity_id", "source_file", "evidence_status",
      "interpretation_note", "qc_flag"
    ),
    artifact
  )
  invisible(TRUE)
}

validate_manuscript_program_summary <- function(df, artifact = "manuscript program summary") {
  require_module_contract_columns(
    df,
    c(
      "program_key", "manuscript_claim_scope", "datasets_supported",
      "evidence_domains", "strongest_evidence", "safe_manuscript_sentence",
      "main_limitation", "qc_flag"
    ),
    artifact
  )
  invisible(TRUE)
}

validate_evidence_priority_matrix <- function(df, artifact = "evidence priority matrix") {
  require_module_contract_columns(
    df,
    c(
      "priority_id", "program_key", "dataset", "priority_tier",
      "evidence_domain_count", "strongest_fdr", "robustness_flag",
      "behavior_flag", "qc_flag", "recommended_use"
    ),
    artifact
  )
  invisible(TRUE)
}

validate_final_evidence_bundle <- function(bundle_dir, artifact = "final biological evidence bundle") {
  required_sheets <- c(
    "README", "input_status", "manuscript_program_summary",
    "evidence_priority_matrix", "cross_compartment_program_atlas",
    "wgcna_key_modules", "wgcna_key_supermodules",
    "microglia_roi_signature_drivers", "qc_flags", "biological_claims"
  )
  if (!dir.exists(bundle_dir)) {
    stop(artifact, " directory does not exist: ", bundle_dir, call. = FALSE)
  }
  missing <- required_sheets[!file.exists(file.path(bundle_dir, paste0(required_sheets, ".csv")))]
  if (length(missing)) {
    stop(artifact, " is missing CSV sheet mirror(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  read_sheet <- function(sheet) {
    utils::read.csv(file.path(bundle_dir, paste0(sheet, ".csv")), check.names = FALSE, stringsAsFactors = FALSE)
  }

  readme <- read_sheet("README")
  require_module_contract_columns(readme, c("sheet", "produced_from", "meaning", "manuscript_safe_columns"), paste(artifact, "README sheet"))
  missing_readme_rows <- setdiff(required_sheets, as.character(readme$sheet))
  if (length(missing_readme_rows)) {
    stop(artifact, " README sheet does not document: ", paste(missing_readme_rows, collapse = ", "), call. = FALSE)
  }

  input_status <- read_sheet("input_status")
  require_module_contract_columns(input_status, c("input_name", "path", "status", "n_rows"), paste(artifact, "input_status sheet"))

  claims <- read_sheet("biological_claims")
  require_module_contract_columns(
    claims,
    c(
      "claim_id", "dataset", "biological_program", "evidence_type",
      "claim_grade", "primary_evidence", "orthogonal_support",
      "major_limitation", "safe_interpretation", "unsafe_overinterpretation",
      "claim_allowed", "claim_gate_status", "claim_downgrade_reason",
      "primary_model_status", "animal_level_gate", "qc_gate",
      "missingness_gate", "batch_confound_gate", "marker_contamination_gate",
      "microglia_roi_gate", "neuropil_independence_gate", "robustness_gate",
      "evidence_independence_gate",
      "missingness_confounded", "batch_or_plate_confounded",
      "region_layer_imbalance_risk", "animal_pseudoreplication_risk",
      "marker_contamination_or_roi_mixture_flag", "qc_interpretation_flag"
    ),
    paste(artifact, "biological_claims sheet")
  )

  modules <- read_sheet("wgcna_key_modules")
  require_module_contract_columns(
    modules,
    c(
      "dataset", "ModuleID", "ModuleColor", "targeted_signature_primary_driver",
      "targeted_signature_driver_class", "targeted_signature_driver_signature",
      "targeted_signature_driver_padj", "targeted_signature_driver_NES",
      "targeted_signature_driver_overlap_proteins"
    ),
    paste(artifact, "wgcna_key_modules sheet")
  )

  drivers <- read_sheet("microglia_roi_signature_drivers")
  require_module_contract_columns(
    drivers,
    c(
      "ModuleID", "ModuleColor", "microenvironment_label",
      "targeted_signature_primary_driver", "targeted_signature_driver_class",
      "targeted_signature_driver_signature", "targeted_signature_driver_padj",
      "targeted_signature_driver_NES", "targeted_signature_driver_overlap_proteins"
    ),
    paste(artifact, "microglia_roi_signature_drivers sheet")
  )

  invisible(TRUE)
}

write_contract_validation_status <- function(path, artifact, ok, message = "") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      artifact = artifact,
      ok = isTRUE(ok),
      message = as.character(message),
      stringsAsFactors = FALSE
    ),
    path,
    row.names = FALSE
  )
  invisible(path)
}
