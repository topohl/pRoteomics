# Helpers for the manuscript-facing final biological evidence bundle.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "wgcna_claim_readiness_utils.R"))
source(repo_path("R", "wgcna_group_effect_consumer_utils.R"))
source(repo_path("R", "wgcna_stage07_semantic_utils.R"))

final_dataset_terminology <- function(dataset) {
  x <- as.character(dataset)
  out <- dplyr::case_when(
    x == "neuron_neuropil" ~ "region/layer-resolved neuron neuropil",
    x == "neuron_soma" ~ "region-resolved neuronal soma-enriched",
    x == "microglia" ~ "region-resolved microglia-enriched ROI/local microenvironment, not purified microglia",
    TRUE ~ x
  )
  out[is.na(out) | !nzchar(out)] <- NA_character_
  out
}

read_final_csv <- function(path) {
  record_input_resolution(
    script = Sys.getenv("PROTEOMICS_SCRIPT_ID", unset = NA_character_),
    dataset = "global",
    stage = "final_evidence_bundle",
    input_name = basename(path),
    expected_path = path,
    resolved_path = path,
    resolution_mode = if (file.exists(path)) "canonical" else "missing_optional",
    strict_mode = strict_inputs_enabled(),
    allowed_in_strict_mode = TRUE,
    producer_script_or_artifact_id = "final_evidence_bundle_source"
  )
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::read_csv(
      path, show_col_types = FALSE, progress = FALSE,
      guess_max = Inf
    )
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

empty_bundle_table <- function(message) {
  data.frame(status = "missing_optional", message = message, stringsAsFactors = FALSE)
}

add_dataset_terminology <- function(df) {
  if (is.null(df) || !nrow(df) || !"dataset" %in% names(df)) return(df)
  if (!"dataset_terminology" %in% names(df)) {
    df$dataset_terminology <- final_dataset_terminology(df$dataset)
  }
  df
}

bundle_input_status <- function(inputs) {
  data.frame(
    input_name = names(inputs),
    path = normalizePath(unname(unlist(inputs)), winslash = "/", mustWork = FALSE),
    status = ifelse(file.exists(unname(unlist(inputs))), "present", "missing_optional"),
    n_rows = vapply(unname(unlist(inputs)), function(path) {
      df <- read_final_csv(path)
      if (is.null(df)) 0L else nrow(df)
    }, integer(1)),
    stringsAsFactors = FALSE
  )
}

select_existing <- function(df, cols) {
  if (is.null(df) || !nrow(df)) return(df)
  df[, intersect(cols, names(df)), drop = FALSE]
}

read_all_dataset_tables <- function(path_fun) {
  pieces <- lapply(valid_datasets(), function(ds) {
    df <- read_final_csv(path_fun(ds))
    if (is.null(df) || !nrow(df)) return(NULL)
    if (!"dataset" %in% names(df)) df$dataset <- ds
    df
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(NULL)
  dplyr::bind_rows(pieces)
}

validated_neuronal_wgcna_key_row <- function(
    tier_specific_fdr, claim_gate, model_stability_status,
    model_type, support_class, model_warning
) {
  fdr <- suppressWarnings(as.numeric(tier_specific_fdr))
  claim_gate <- as.character(claim_gate)
  model_stability <- as.character(model_stability_status)
  model_text <- tolower(paste(
    as.character(model_type), as.character(support_class),
    as.character(model_warning)
  ))
  !is.na(fdr) & fdr <= 0.10 &
    claim_gate == "eligible_for_readiness_assessment" &
    !grepl("diagnostic|unstable|failed", model_stability) &
    !grepl("fallback|diagnostic|model_unstable", model_text)
}

append_no_validated_neuronal_key_status <- function(df) {
  neuronal <- c("neuron_soma", "neuron_neuropil")
  present <- unique(as.character(df$dataset[df$dataset %in% neuronal]))
  missing <- setdiff(neuronal, present)
  if (!length(missing)) return(df)
  status <- data.frame(
    dataset = missing,
    status = "no_validated_neuronal_wgcna_key_rows",
    reason = "neuronal readiness contract unavailable and all current Stage 07 inferential rows are claim-ineligible",
    stringsAsFactors = FALSE
  )
  dplyr::bind_rows(df, add_dataset_terminology(status))
}

build_wgcna_key_modules <- function() {
  modules <- read_all_dataset_tables(function(ds) {
    path_results(
      "tables", "06_modules_WGCNA", "interpretable_summary", ds,
      "WGCNA_inferential_handoff.csv"
    )
  })
  if (is.null(modules) || !nrow(modules)) {
    return(empty_bundle_table("WGCNA inferential handoffs were not available."))
  }
  modules <- wgcna_inferential_handoff_normalize_csv_types(
    modules, "combined module inferential handoffs"
  )
  wgcna_stage07_validate_inferential_handoff(
    modules, "combined module inferential handoffs"
  )
  modules <- modules |>
    dplyr::filter(
      .data$entity_level == "module",
      .data$analysis_tier == "primary_wgcna_global",
      .data$contrast == "SUS - RES",
      .data$result_scope == "primary_global_vulnerability"
    )
  if (!"module_id" %in% names(modules)) modules$module_id <- NA_character_
  if (!"ModuleID" %in% names(modules)) modules$ModuleID <- NA_character_
  micro_independence <- read_final_csv(path_results("tables", "06_modules_WGCNA", "microglia_neuropil_independence", "microglia", "microglia_module_neuropil_independence_classification.csv"))
  if (!is.null(micro_independence) && nrow(micro_independence)) {
    if (!"module_id" %in% names(micro_independence)) micro_independence$module_id <- NA_character_
    if (!"endpoint_id" %in% names(micro_independence)) micro_independence$endpoint_id <- NA_character_
    if (!"module_id" %in% names(modules)) modules$module_id <- NA_character_
    if (!"ModuleID" %in% names(modules)) modules$ModuleID <- NA_character_
    micro_independence <- micro_independence |>
      dplyr::select(dplyr::any_of(c(
        "dataset", "module_id", "endpoint_id", "neuropil_independence_classification",
        "best_contrast", "min_p_before", "min_p_after", "max_percent_attenuation",
        "matched_neuropil_covariate", "neuropil_covariate_source",
        "neuropil_independence_note"
      ))) |>
      dplyr::mutate(module_join_id = dplyr::coalesce(.data$module_id, .data$endpoint_id)) |>
      dplyr::distinct(.data$dataset, .data$module_join_id, .keep_all = TRUE)
    modules <- modules |>
      dplyr::mutate(module_join_id = dplyr::coalesce(.data$module_id, .data$ModuleID)) |>
      dplyr::left_join(
        micro_independence |> dplyr::select(-dplyr::any_of(c("module_id", "endpoint_id"))),
        by = c("dataset", "module_join_id")
      ) |>
      dplyr::select(-"module_join_id")
  }
  for (col in c("tier_specific_fdr", "targeted_signature_driver_padj", "estimate")) {
    if (!col %in% names(modules)) modules[[col]] <- NA_real_
    modules[[col]] <- suppressWarnings(as.numeric(modules[[col]]))
  }
  if (!"claim_gate" %in% names(modules)) {
    modules$claim_gate <- "not_claim_allowed_model"
  }
  for (col in c("model_type", "support_class", "model_warning", "model_stability_status")) {
    if (!col %in% names(modules)) modules[[col]] <- NA_character_
  }
  for (col in c(
    "canonical_claim_entity_id", "claim_entity_role",
    "primary_architecture_status", "spatial_dependence_class",
    "animal_stability_status", "group_effect_status",
    "allowed_claim_scope", "prohibited_claim_scope",
    "readiness_contract_version"
  )) {
    if (!col %in% names(modules)) modules[[col]] <- NA_character_
  }
  if (!"separate_manuscript_claim_allowed" %in% names(modules)) {
    modules$separate_manuscript_claim_allowed <- NA
  }
  contract <- load_microglia_wgcna_claim_readiness()
  stage13_modules <- contract$canonical_modules |>
    dplyr::transmute(
      dataset = as.character(.data$dataset),
      stage13_module_id = as.character(.data$entity_id),
      stage13_canonical_claim_entity_id = as.character(.data$canonical_claim_entity_id),
      stage13_claim_entity_role = as.character(.data$claim_entity_role),
      stage13_separate_manuscript_claim_allowed = as.logical(.data$separate_manuscript_claim_allowed),
      stage13_primary_architecture_status = as.character(.data$primary_architecture_status),
      stage13_spatial_dependence_class = as.character(.data$spatial_dependence_class),
      stage13_animal_stability_status = as.character(.data$animal_stability_status),
      stage13_group_effect_status = as.character(.data$group_effect_status),
      stage13_allowed_claim_scope = as.character(.data$allowed_claim_scope),
      stage13_prohibited_claim_scope = as.character(.data$prohibited_claim_scope),
      stage13_readiness_contract_version = as.character(.data$readiness_contract_version),
      stage13_source_file = contract$source_path
    )
  modules <- modules |>
    dplyr::mutate(stage13_module_id = dplyr::coalesce(as.character(.data$module_id), as.character(.data$ModuleID))) |>
    dplyr::left_join(stage13_modules, by = c("dataset", "stage13_module_id"), relationship = "many-to-one") |>
    dplyr::mutate(
      canonical_claim_entity_id = dplyr::coalesce(.data$stage13_canonical_claim_entity_id, as.character(.data$canonical_claim_entity_id)),
      claim_entity_role = dplyr::coalesce(.data$stage13_claim_entity_role, as.character(.data$claim_entity_role)),
      separate_manuscript_claim_allowed = dplyr::coalesce(.data$stage13_separate_manuscript_claim_allowed, as.logical(.data$separate_manuscript_claim_allowed)),
      primary_architecture_status = dplyr::coalesce(.data$stage13_primary_architecture_status, as.character(.data$primary_architecture_status)),
      spatial_dependence_class = dplyr::coalesce(.data$stage13_spatial_dependence_class, as.character(.data$spatial_dependence_class)),
      animal_stability_status = dplyr::coalesce(.data$stage13_animal_stability_status, as.character(.data$animal_stability_status)),
      group_effect_status = dplyr::coalesce(.data$stage13_group_effect_status, as.character(.data$group_effect_status)),
      allowed_claim_scope = dplyr::coalesce(.data$stage13_allowed_claim_scope, as.character(.data$allowed_claim_scope)),
      prohibited_claim_scope = dplyr::coalesce(.data$stage13_prohibited_claim_scope, as.character(.data$prohibited_claim_scope)),
      readiness_contract_version = dplyr::coalesce(.data$stage13_readiness_contract_version, as.character(.data$readiness_contract_version))
    ) |>
    dplyr::filter(
        (.data$dataset == "microglia" & .data$claim_entity_role == "canonical_module" & .data$separate_manuscript_claim_allowed %in% TRUE) |
        (.data$dataset != "microglia" & validated_neuronal_wgcna_key_row(
          .data$tier_specific_fdr, .data$claim_gate,
          .data$model_stability_status, .data$model_type,
          .data$support_class, .data$model_warning
        ))
    ) |>
    dplyr::mutate(
      wgcna_claim_semantic_status = dplyr::case_when(
        .data$dataset == "microglia" & .data$group_effect_status == "FDR_supported" ~ "stress_effect_eligible",
        .data$dataset == "microglia" & (!is.na(.data$targeted_signature_primary_driver) | !is.na(.data$best_wgcna_de_gsea_overlap)) ~ "convergent",
        .data$dataset == "microglia" ~ "architecture_only",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::arrange(.data$dataset, .data$tier_specific_fdr)
  cols <- c(
    "dataset", "dataset_terminology", "status", "reason", "module_id", "ModuleID", "ModuleColor",
    "canonical_claim_entity_id", "claim_entity_role", "separate_manuscript_claim_allowed",
    "primary_architecture_status", "spatial_dependence_class", "animal_stability_status",
    "group_effect_status", "allowed_claim_scope", "prohibited_claim_scope",
    "readiness_contract_version", "stage13_source_file", "wgcna_claim_semantic_status",
    "ModulePlotLabel", "Module_CleanPlotLabel", "module_biological_label",
    "module_biological_label_short", "module_label_display", "module_label", "module_display_label",
    "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms", "raw_top_GO_label",
    "raw_module_label", "raw_hub_proteins", "raw_marker_or_signature_label",
    "cleaned_biological_label", "cleaned_biological_label_short",
    "cleaned_biological_label_source", "cleaned_biological_label_confidence",
    "GO_label_relevance_flag", "GO_label_relevance_rationale",
    "microenvironment_caution_label", "microenvironment_caution_class",
    "microenvironment_caution_rationale", "Supermodule_CleanPlotLabel",
    "Supermodule_CompositionDisplayLabel", "Supermodule_CompositionLabel",
    "Supermodule_PlotLabel",
    "Supermodule_FullAnnotationLabel", "Supermodule_DisplayShort",
    "supermodule_id", "contrast", "estimate", "p_value",
    "tier_specific_fdr", "tier_specific_family_id",
    "tier_specific_family_size", "direction", "support_class",
    "claim_gate", "animal_level_status",
    "n_animals", "model_warning", "microenvironment_class", "microenvironment_label",
    "microenvironment_confidence", "dominant_microenvironment_class",
    "targeted_signature_primary_driver", "targeted_signature_driver_class",
    "targeted_signature_driver_signature", "targeted_signature_driver_padj",
    "targeted_signature_driver_NES", "targeted_signature_driver_overlap_proteins",
    "targeted_signature_driver_evidence_tier",
    "n_unique_targeted_signatures", "n_unique_targeted_overlap_proteins",
    "n_unique_curated_program_signatures", "n_unique_curated_program_overlap_proteins",
    "curated_program_overlap_proteins", "curated_program_overlap_warning",
    "targeted_signature_overlap_interpretation_note",
    "best_targeted_microglia_enriched_signatures",
    "best_targeted_curated_microglia_programs",
    "best_targeted_neuropil_shared_signatures",
    "neuropil_independence_classification", "best_contrast",
    "min_p_before", "min_p_after", "max_percent_attenuation",
    "matched_neuropil_covariate", "neuropil_covariate_source",
    "neuropil_independence_note",
    "best_wgcna_de_gsea_overlap", "interpretation_sentence"
  )
  append_no_validated_neuronal_key_status(select_existing(add_dataset_terminology(modules), cols))
}

build_wgcna_key_supermodules <- function() {
  supers <- read_all_dataset_tables(function(ds) {
    path_results(
      "tables", "06_modules_WGCNA", "interpretable_summary", ds,
      "WGCNA_inferential_handoff.csv"
    )
  })
  if (is.null(supers) || !nrow(supers)) {
    return(empty_bundle_table("WGCNA inferential handoffs were not available."))
  }
  supers <- wgcna_inferential_handoff_normalize_csv_types(
    supers, "combined supermodule inferential handoffs"
  )
  wgcna_stage07_validate_inferential_handoff(
    supers, "combined supermodule inferential handoffs"
  )
  supers <- supers |>
    dplyr::filter(
      .data$entity_level == "supermodule",
      .data$analysis_tier == "primary_wgcna_global",
      .data$contrast == "SUS - RES",
      .data$result_scope == "primary_global_vulnerability"
    )
  for (col in c("tier_specific_fdr", "estimate")) {
    if (!col %in% names(supers)) supers[[col]] <- NA_real_
    supers[[col]] <- suppressWarnings(as.numeric(supers[[col]]))
  }
  if (!"claim_gate" %in% names(supers)) {
    supers$claim_gate <- "not_claim_allowed_model"
  }
  for (col in c("model_type", "support_class", "model_warning", "model_stability_status")) {
    if (!col %in% names(supers)) supers[[col]] <- NA_character_
  }
  for (col in c(
    "canonical_claim_entity_id", "claim_entity_role",
    "primary_architecture_status", "spatial_dependence_class",
    "animal_stability_status", "group_effect_status",
    "allowed_claim_scope", "prohibited_claim_scope",
    "readiness_contract_version"
  )) {
    if (!col %in% names(supers)) supers[[col]] <- NA_character_
  }
  if (!"separate_manuscript_claim_allowed" %in% names(supers)) {
    supers$separate_manuscript_claim_allowed <- NA
  }
  contract <- load_microglia_wgcna_claim_readiness()
  stage13_supers <- contract$all |>
    dplyr::filter(.data$level == "supermodule") |>
    dplyr::transmute(
      dataset = as.character(.data$dataset),
      stage13_supermodule_id = as.character(.data$entity_id),
      stage13_canonical_claim_entity_id = as.character(.data$canonical_claim_entity_id),
      stage13_claim_entity_role = as.character(.data$claim_entity_role),
      stage13_separate_manuscript_claim_allowed = as.logical(.data$separate_manuscript_claim_allowed),
      stage13_primary_architecture_status = as.character(.data$primary_architecture_status),
      stage13_spatial_dependence_class = as.character(.data$spatial_dependence_class),
      stage13_animal_stability_status = as.character(.data$animal_stability_status),
      stage13_group_effect_status = as.character(.data$group_effect_status),
      stage13_allowed_claim_scope = as.character(.data$allowed_claim_scope),
      stage13_prohibited_claim_scope = as.character(.data$prohibited_claim_scope),
      stage13_readiness_contract_version = as.character(.data$readiness_contract_version),
      stage13_source_file = contract$source_path
    )
  supers <- supers |>
    dplyr::mutate(stage13_supermodule_id = as.character(.data$supermodule_id)) |>
    dplyr::left_join(stage13_supers, by = c("dataset", "stage13_supermodule_id"), relationship = "many-to-one") |>
    dplyr::mutate(
      canonical_claim_entity_id = dplyr::coalesce(.data$stage13_canonical_claim_entity_id, as.character(.data$canonical_claim_entity_id)),
      claim_entity_role = dplyr::coalesce(.data$stage13_claim_entity_role, as.character(.data$claim_entity_role)),
      separate_manuscript_claim_allowed = dplyr::coalesce(.data$stage13_separate_manuscript_claim_allowed, as.logical(.data$separate_manuscript_claim_allowed)),
      primary_architecture_status = dplyr::coalesce(.data$stage13_primary_architecture_status, as.character(.data$primary_architecture_status)),
      spatial_dependence_class = dplyr::coalesce(.data$stage13_spatial_dependence_class, as.character(.data$spatial_dependence_class)),
      animal_stability_status = dplyr::coalesce(.data$stage13_animal_stability_status, as.character(.data$animal_stability_status)),
      group_effect_status = dplyr::coalesce(.data$stage13_group_effect_status, as.character(.data$group_effect_status)),
      allowed_claim_scope = dplyr::coalesce(.data$stage13_allowed_claim_scope, as.character(.data$allowed_claim_scope)),
      prohibited_claim_scope = dplyr::coalesce(.data$stage13_prohibited_claim_scope, as.character(.data$prohibited_claim_scope)),
      readiness_contract_version = dplyr::coalesce(.data$stage13_readiness_contract_version, as.character(.data$readiness_contract_version))
    ) |>
    dplyr::filter(
        (.data$dataset == "microglia" & .data$claim_entity_role == "higher_order_block" & .data$separate_manuscript_claim_allowed %in% TRUE) |
        (.data$dataset != "microglia" & validated_neuronal_wgcna_key_row(
          .data$tier_specific_fdr, .data$claim_gate,
          .data$model_stability_status, .data$model_type,
          .data$support_class, .data$model_warning
        ))
    ) |>
    dplyr::mutate(
      wgcna_claim_semantic_status = dplyr::case_when(
        .data$dataset == "microglia" & .data$group_effect_status == "FDR_supported" ~ "stress_effect_eligible",
        .data$dataset == "microglia" ~ "architecture_only",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::arrange(.data$dataset, .data$tier_specific_fdr)
  cols <- c(
    "dataset", "dataset_terminology", "status", "reason", "supermodule_id", "Supermodule_PlotLabel",
    "canonical_claim_entity_id", "claim_entity_role", "separate_manuscript_claim_allowed",
    "primary_architecture_status", "spatial_dependence_class", "animal_stability_status",
    "group_effect_status", "allowed_claim_scope", "prohibited_claim_scope",
    "readiness_contract_version", "stage13_source_file", "wgcna_claim_semantic_status",
    "Supermodule_CleanPlotLabel", "Supermodule_CompositionDisplayLabel",
    "Supermodule_FullAnnotationLabel", "Supermodule_DisplayShort",
    "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display",
    "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms", "raw_top_GO_label",
    "raw_module_label", "raw_hub_proteins", "raw_marker_or_signature_label",
    "cleaned_biological_label", "cleaned_biological_label_short",
    "cleaned_biological_label_source", "cleaned_biological_label_confidence",
    "GO_label_relevance_flag", "GO_label_relevance_rationale",
    "microenvironment_caution_label", "microenvironment_caution_class",
    "microenvironment_caution_rationale",
    "Supermodule_ConservativeLabel", "Supermodule_CompositionLabel",
    "Supermodule_CompositionShortLabel", "Supermodule_CompositionLabelSource",
    "Supermodule_CompositionConfidence", "Supermodule_CompositionRationale",
    "DominantMemberTheme", "DominantMemberThemeFraction", "SecondMemberTheme",
    "SecondMemberThemeFraction", "TopMemberModuleLabels", "TopMemberGOTerms",
    "n_member_modules_with_informative_labels",
    "fraction_member_modules_with_informative_labels",
    "dominant_microenvironment_class", "dominant_module_labels",
    "Supermodule_LabelRationale", "contrast", "estimate", "p_value",
    "tier_specific_fdr", "tier_specific_family_id",
    "tier_specific_family_size", "direction", "support_class", "claim_gate",
    "animal_level_status", "n_animals", "model_warning", "interpretation_sentence"
  )
  append_no_validated_neuronal_key_status(select_existing(add_dataset_terminology(supers), cols))
}

build_microglia_roi_signature_drivers <- function() {
  modules <- read_final_csv(path_results("tables", "06_modules_WGCNA", "module_annotation", "microglia", "WGCNA_module_biological_annotation.csv"))
  if (is.null(modules) || !nrow(modules)) {
    return(empty_bundle_table("Microglia ROI WGCNA module annotation was not available."))
  }
  cols <- c(
    "dataset", "dataset_terminology", "ModuleID", "ModuleColor", "module_display_label",
    "microenvironment_class", "microenvironment_label", "microenvironment_confidence",
    "n_targeted_microglia_enriched_empirical_overlaps",
    "n_targeted_microglia_enriched_reference_supported_overlaps",
    "n_targeted_curated_microglia_program_overlaps",
    "n_targeted_neuropil_shared_overlaps", "n_targeted_ambiguous_overlaps",
    "n_unique_targeted_signatures", "n_unique_targeted_overlap_proteins",
    "n_unique_curated_program_signatures", "n_unique_curated_program_overlap_proteins",
    "curated_program_overlap_proteins", "curated_program_overlap_warning",
    "targeted_signature_overlap_interpretation_note",
    "targeted_signature_driver_evidence_tier",
    "targeted_signature_primary_driver", "targeted_signature_driver_class",
    "targeted_signature_driver_signature", "targeted_signature_driver_padj",
    "targeted_signature_driver_NES", "targeted_signature_driver_overlap_proteins",
    "best_targeted_microglia_enriched_signatures",
    "best_targeted_curated_microglia_programs",
    "best_targeted_neuropil_shared_signatures",
    "best_targeted_microglia_signatures",
    "classification_rationale", "interpretation_note"
  )
  out <- select_existing(add_dataset_terminology(modules), cols)
  out <- out[!is.na(out$targeted_signature_primary_driver) & nzchar(as.character(out$targeted_signature_primary_driver)), , drop = FALSE]
  if (!nrow(out)) empty_bundle_table("No targeted microglia ROI signature drivers were available.") else out
}

build_microglia_neuropil_independence <- function() {
  tbl <- read_final_csv(path_results("tables", "06_modules_WGCNA", "microglia_neuropil_independence", "microglia", "microglia_neuropil_independence_effects.csv"))
  if (is.null(tbl) || !nrow(tbl)) {
    return(empty_bundle_table("Microglia neuropil-independence sensitivity analysis was not available."))
  }
  cols <- c(
    "dataset", "dataset_terminology", "endpoint_type", "endpoint_id", "endpoint_label",
    "module_id", "contrast", "estimate_before", "estimate_after",
    "percent_attenuation", "p_before", "p_after", "FDR_before", "FDR_after",
    "neuropil_covariate_beta", "neuropil_covariate_p", "classification",
    "matched_neuropil_covariate", "matched_neuropil_label",
    "neuropil_covariate_source", "neuropil_selection_spearman",
    "n_samples", "n_matched_samples", "n_animals", "model_warning"
  )
  select_existing(add_dataset_terminology(tbl), cols)
}

build_microglia_wgcna_claim_readiness <- function() {
  contract <- load_microglia_wgcna_claim_readiness()
  out <- as.data.frame(contract$all, stringsAsFactors = FALSE)
  out$stage13_source_file <- contract$source_path
  out
}

build_qc_flags <- function(claims) {
  if (is.null(claims) || !nrow(claims) || !"dataset" %in% names(claims)) {
    return(empty_bundle_table("Biological claims QC flags were not available."))
  }
  if (!"batch_or_plate_confounded" %in% names(claims)) {
    claims$batch_or_plate_confounded <- if ("plate_or_batch_confounded" %in% names(claims)) claims$plate_or_batch_confounded else NA_character_
  }
  if (!"marker_contamination_or_roi_mixture_flag" %in% names(claims)) {
    claims$marker_contamination_or_roi_mixture_flag <- if ("marker_contamination_risk" %in% names(claims)) claims$marker_contamination_risk else NA_character_
  }
  cols <- c(
    "dataset", "dataset_terminology", "missingness_confounded",
    "plate_or_batch_confounded", "batch_or_plate_confounded",
    "region_layer_imbalance_risk", "animal_pseudoreplication_risk",
    "early_pc_association", "marker_contamination_risk",
    "marker_contamination_or_roi_mixture_flag", "qc_interpretation_flag",
    "animal_level_status", "claim_type", "claim_use_class",
    "raw_top_GO_term", "representative_GO_terms", "semantic_parent_label",
    "safe_program_label", "term_label_risk", "label_confidence",
    "label_basis", "label_downgrade_reason",
    "claim_allowed", "claim_gate_status",
    "claim_downgrade_reason", "model_fit_status", "statistical_evidence_status",
    "claim_gate_model_status", "primary_model_status", "animal_level_gate",
    "qc_gate", "missingness_gate", "batch_confound_gate",
    "marker_contamination_gate", "microglia_roi_gate",
    "neuropil_independence_gate", "robustness_gate",
    "evidence_independence_gate"
  )
  add_dataset_terminology(claims) |>
    select_existing(cols) |>
    dplyr::distinct()
}

build_bundle_readme <- function() {
  data.frame(
    sheet = c(
      "README", "input_status", "reviewer_audit_manifest", "final_validation", "manuscript_program_summary", "evidence_priority_matrix",
      "cross_compartment_program_atlas", "wgcna_key_modules", "wgcna_key_supermodules",
      "microglia_wgcna_claim_readiness",
      "microglia_roi_signature_drivers", "microglia_neuropil_independence",
      "qc_flags", "biological_claims"
    ),
    produced_from = c(
      "10_biological_integration/03_evidence_priority_matrix.r and 09_export_pride_journal/07_make_biological_claims_table.R",
      "Final bundle input scanner",
      "09_export_pride_journal/07_make_biological_claims_table.R",
      "09_export_pride_journal/07_make_biological_claims_table.R",
      "10_biological_integration/02_manuscript_program_summary.r",
      "10_biological_integration/03_evidence_priority_matrix.r",
      "10_biological_integration/01_cross_compartment_program_atlas.r",
      "06_modules_WGCNA/07_wgcna_interpretable_summary.r",
      "06_modules_WGCNA/07_wgcna_interpretable_summary.r",
      "06_modules_WGCNA/13_wgcna_claim_readiness.R",
      "06_modules_WGCNA/06_annotate_module_microenvironment.r",
      "06_modules_WGCNA/09_microglia_neuropil_independence.R",
      "09_export_pride_journal/07_make_biological_claims_table.R",
      "09_export_pride_journal/07_make_biological_claims_table.R"
    ),
    meaning = c(
      "Plain-language index for the workbook.",
      "Which optional upstream tables were present when the bundle was assembled.",
      "Index of final reviewer-facing claim and audit tables, their availability, row counts, and schema-validation status.",
      "Machine-readable PASS/FAIL checks for final claim wording, gate consistency, endpoint scope, schemas, and manifest completeness.",
      "High-level program synthesis for manuscript drafting.",
      "Program-level priority tiers based on evidence convergence and QC flags.",
      "Long evidence atlas across enrichment, WGCNA, microenvironment, robustness, behavior, and QC streams.",
      "Module-level WGCNA rows selected for manuscript inspection; includes targeted-signature driver columns when available.",
      "Supermodule-level WGCNA rows selected for manuscript inspection; includes original annotation and plotting labels.",
      "Authoritative microglia WGCNA Stage 13 handoff preserving all 22 technical identities; aliases are compatibility-only and architecture eligibility is distinct from FDR-supported stress-effect eligibility.",
      "Microglia ROI/local microenvironment targeted-signature summary; curated single-protein or neuropil-shared evidence is cautionary, not purified microglia evidence.",
      "Sensitivity model comparing microglia group effects before and after matched region-level neuron_neuropil covariate adjustment.",
      "Dataset-level QC and confounding flags propagated from the claims table.",
      "Reviewer-facing claim gate. claim_type controls gate applicability; claim_use_class and safe_program_label separate manuscript use from raw GO labels; claim_grade is descriptive."
    ),
    manuscript_safe_columns = c(
      "Use this sheet as documentation only.",
      "status, n_rows, path.",
      "audit_file, exists, n_rows, schema_validated, produced_by_script, reviewer_use, manuscript_use_allowed, notes.",
      "validation_check, status, n_violations, details.",
      "program_key, manuscript_claim_scope, strongest_evidence, safe_manuscript_sentence, main_limitation, qc_flag.",
      "priority_tier, independent_evidence_family_count, evidence_source_families, strongest_fdr, robustness_flag, behavior_flag, qc_flag, recommended_use.",
      "program_label, evidence_domain, evidence_source_family, contrast, effect_size, fdr, support_count, evidence_strength, interpretation_note, qc_flag.",
      "ModulePlotLabel, Supermodule_PlotLabel, Supermodule_FullAnnotationLabel, microenvironment_label, targeted_signature_* driver columns, FDR columns, interpretation_sentence.",
      "Supermodule_PlotLabel, conservative/composition labels, dominant_microenvironment_class, dominant_module_labels, Supermodule_LabelRationale, FDR columns, interpretation_sentence.",
      "dataset, level, entity_id, canonical_claim_entity_id, claim_entity_role, separate_manuscript_claim_allowed, primary_architecture_status, group_effect_status, allowed_claim_scope, prohibited_claim_scope, readiness_contract_version, stage13_source_file.",
      "targeted_signature_primary_driver, targeted_signature_driver_class, targeted_signature_driver_signature, targeted_signature_driver_padj, targeted_signature_driver_NES, targeted_signature_driver_overlap_proteins, targeted_signature_driver_evidence_tier, curated_program_overlap_warning.",
      "classification, estimate_before, estimate_after, percent_attenuation, FDR_before, FDR_after, neuropil_covariate_beta, neuropil_covariate_p.",
      "missingness_confounded, batch_or_plate_confounded, region_layer_imbalance_risk, animal_pseudoreplication_risk, marker_contamination_or_roi_mixture_flag, qc_interpretation_flag.",
      "claim_type, claim_use_class, raw_top_GO_term, representative_GO_terms, safe_program_label, term_label_risk, label_confidence, claim_allowed, claim_gate_status, claim_downgrade_reason, model_fit_status, statistical_evidence_status, claim_gate_model_status, gate columns, claim_grade, primary_evidence, orthogonal_support, major_limitation, safe_interpretation, unsafe_overinterpretation, QC flag columns."
    ),
    notes = c(
      "Final-facing terminology: neuron_neuropil = region/layer-resolved neuron neuropil; neuron_soma = region-resolved neuronal soma-enriched; microglia = region-resolved microglia-enriched ROI/local microenvironment, not purified microglia.",
      "Missing optional inputs are expected in some partial reruns and should limit manuscript strength.",
      "Start the reviewer audit here; manuscript_use_allowed refers to the table type, while row-level eligibility remains in biological_claims.",
      "Every validation row should report PASS before final reviewer handoff.",
      "Use for narrative synthesis, then verify source rows.",
      "Tiering is a prioritization aid, not a new statistical test.",
      "This is an evidence index, not a causal model.",
      "Microglia/immune wording should be used only where explicit immune-state evidence supports it.",
      "Both broad plot labels and full annotation labels are retained.",
      "All 22 Stage 13 identities are retained for audit. Six singleton supermodules are compatibility aliases, SM01/SM03 are descriptive blocks only, and SM09 is the only independently eligible higher-order block; none has an FDR-supported Stage 13 stress effect.",
      "Neuropil-shared and curated-program evidence should not be described as purified microglial activation.",
      "Adjustment is a sensitivity analysis, not causal subtraction of neuropil signal.",
      "QC flags are conservative and should be reviewed before strong wording.",
      "claim_grade is not an eligibility decision. GO-derived labels preserve raw terms but use safe_program_label and claim_use_class for manuscript wording; tissue-mismatched GO terms are flagged and conservatively relabeled rather than discarded."
    ),
    stringsAsFactors = FALSE
  )
}

write_bundle_csvs <- function(sheets, bundle_dir) {
  dir_create(bundle_dir)
  for (nm in names(sheets)) {
    path <- file.path(bundle_dir, paste0(nm, ".csv"))
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_csv(sheets[[nm]], path, na = "")
    } else {
      utils::write.csv(sheets[[nm]], path, row.names = FALSE, na = "")
    }
  }
}

write_final_evidence_bundle <- function(reason = "integration") {
  bundle_dir <- path_results("tables", "10_biological_integration", "final_evidence_bundle", "global")
  source_bundle_dir <- path_results("source_data", "10_biological_integration", "final_evidence_bundle", "global")
  report_dir <- path_results("reports", "10_biological_integration", "final_evidence_bundle", "global")
  log_dir <- path_results("logs", "10_biological_integration", "final_evidence_bundle", "global")
  dir_create(bundle_dir); dir_create(source_bundle_dir); dir_create(report_dir); dir_create(log_dir)

  inputs <- list(
    manuscript_program_summary = path_results("tables", "10_biological_integration", "manuscript_program_summary", "global", "manuscript_program_summary.csv"),
    evidence_priority_matrix = path_results("tables", "10_biological_integration", "evidence_priority_matrix", "global", "evidence_priority_matrix.csv"),
    cross_compartment_program_atlas = path_results("tables", "10_biological_integration", "cross_compartment_program_atlas", "global", "cross_compartment_program_atlas_long.csv"),
    biological_claims = path_results("tables", "biological_claims_table.csv"),
    microglia_module_annotation = path_results("tables", "06_modules_WGCNA", "module_annotation", "microglia", "WGCNA_module_biological_annotation.csv"),
    microglia_targeted_signature_details = path_results("tables", "06_modules_WGCNA", "module_annotation", "microglia", "WGCNA_module_targeted_signature_overlap_details.csv"),
    microglia_neuropil_independence = path_results("tables", "06_modules_WGCNA", "microglia_neuropil_independence", "microglia", "microglia_neuropil_independence_effects.csv"),
    reviewer_audit_manifest = path_results("reviewer_audit", "final_reviewer_audit_manifest.csv"),
    final_validation = path_results("reviewer_audit", "final_evidence_bundle_validation.csv"),
    microglia_wgcna_claim_readiness = microglia_wgcna_claim_readiness_path()
  )

  manuscript <- add_dataset_terminology(read_final_csv(inputs$manuscript_program_summary) %||% empty_bundle_table("Manuscript program summary was not available."))
  priority <- add_dataset_terminology(read_final_csv(inputs$evidence_priority_matrix) %||% empty_bundle_table("Evidence priority matrix was not available."))
  atlas <- add_dataset_terminology(read_final_csv(inputs$cross_compartment_program_atlas) %||% empty_bundle_table("Cross-compartment program atlas was not available."))
  claims <- add_dataset_terminology(read_final_csv(inputs$biological_claims) %||% empty_bundle_table("Biological claims table was not available."))
  reviewer_audit_manifest <- read_final_csv(inputs$reviewer_audit_manifest) %||% empty_bundle_table("Final reviewer audit manifest was not available.")
  final_validation <- read_final_csv(inputs$final_validation) %||% empty_bundle_table("Final evidence validation was not available.")
  if ("plate_or_batch_confounded" %in% names(claims) && !"batch_or_plate_confounded" %in% names(claims)) {
    claims$batch_or_plate_confounded <- claims$plate_or_batch_confounded
  }
  if ("marker_contamination_risk" %in% names(claims) && !"marker_contamination_or_roi_mixture_flag" %in% names(claims)) {
    claims$marker_contamination_or_roi_mixture_flag <- claims$marker_contamination_risk
  }

  sheets <- list(
    README = build_bundle_readme(),
    input_status = bundle_input_status(inputs),
    reviewer_audit_manifest = reviewer_audit_manifest,
    final_validation = final_validation,
    manuscript_program_summary = manuscript,
    evidence_priority_matrix = priority,
    cross_compartment_program_atlas = atlas,
    wgcna_key_modules = build_wgcna_key_modules(),
    wgcna_key_supermodules = build_wgcna_key_supermodules(),
    microglia_wgcna_claim_readiness = build_microglia_wgcna_claim_readiness(),
    microglia_roi_signature_drivers = build_microglia_roi_signature_drivers(),
    microglia_neuropil_independence = build_microglia_neuropil_independence(),
    qc_flags = build_qc_flags(claims),
    biological_claims = claims
  )

  write_bundle_csvs(sheets, bundle_dir)
  write_bundle_csvs(sheets, source_bundle_dir)
  xlsx_out <- file.path(bundle_dir, "final_biological_evidence_bundle.xlsx")
  if (requireNamespace("writexl", quietly = TRUE)) {
    tryCatch(
      writexl::write_xlsx(sheets, xlsx_out),
      error = function(e) warning("Could not write final_biological_evidence_bundle.xlsx: ", conditionMessage(e), call. = FALSE)
    )
  }
  write_run_manifest(
    file.path(log_dir, "run_manifest.yml"),
    inputs = inputs,
    outputs = list(
      bundle_dir = bundle_dir,
      source_bundle_dir = source_bundle_dir,
      workbook = if (file.exists(xlsx_out)) xlsx_out else NA_character_
    ),
    parameters = list(reason = reason),
    notes = "Final manuscript-facing biological evidence bundle assembled from existing downstream outputs only."
  )
  invisible(list(bundle_dir = bundle_dir, workbook = xlsx_out, sheets = names(sheets)))
}
