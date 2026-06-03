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
      "classification_threshold", "classification_rationale", "interpretation_note"
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
