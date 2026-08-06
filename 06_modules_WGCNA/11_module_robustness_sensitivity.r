#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/11_module_robustness_sensitivity.r
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: Stage 07 inferential handoff, preservation summaries, and QC/confounding outputs.
# Produces: module robustness/sensitivity evidence.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "wgcna_group_effect_consumer_utils.R"))

SCRIPT_ID <- "06_modules_WGCNA/11_module_robustness_sensitivity.r"
run <- integration_cli(allow_all = TRUE)

claim_gate_rows <- list()

status_from_effects <- function(effects) {
  if (is.null(effects) || !nrow(effects)) return(rep("missing_group_effects", 0))
  dplyr::case_when(
    effects$claim_gate == "diagnostic_only_model" ~
      "diagnostic_only_model",
    effects$claim_gate == "not_claim_allowed_model" ~
      "model_not_claim_allowed",
    TRUE ~ "effect_available"
  )
}

make_claim_gate_audit <- function(ds, effects, level, source_file, preservation = NULL, status = NULL) {
  if (is.null(effects) || !nrow(effects)) {
    return(data.frame(
      dataset = ds,
      module_or_supermodule_id = NA_character_,
      level = level,
      contrast = NA_character_,
      spatial_unit = NA_character_,
      effect_scope = NA_character_,
      robustness_gate = "missing_required",
      preservation_status = "not_available",
      sensitivity_status = "missing_group_effects",
      direction_stability = "not_available",
      confounding_status = "not_available",
      claim_gate_eligible = FALSE,
      robustness_downgrade_reason = paste0("missing_", level, "_group_effects"),
      source_file = source_file,
      source_key = NA_character_,
      tier_specific_family_id = NA_character_,
      tier_specific_family_size = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  id_col <- "entity_id"
  contrast_col <- first_col(effects, c("contrast"))
  spatial_col <- first_col(effects, c("spatial_unit", "SpatialLabel"))
  scope_col <- first_col(effects, c("effect_scope"))
  direction_col <- first_col(effects, c("direction"))
  preservation_status <- rep("not_available", nrow(effects))
  if (identical(level, "module") && !is.null(preservation) && nrow(preservation) && !is.na(id_col)) {
    pres_key <- first_col(preservation, c("ModuleID", "module_id", "ModuleColor", "module"))
    pres_col <- first_col(preservation, c("preservation_Zsummary_median", "Zsummary", "medianRank"))
    if (!is.na(pres_key) && !is.na(pres_col)) {
      pres <- suppressWarnings(as.numeric(preservation[[pres_col]][match(effects[[id_col]], preservation[[pres_key]])]))
      preservation_status <- dplyr::case_when(
        is.na(pres) ~ "not_available",
        pres >= 10 ~ "strong_preservation",
        pres >= 2 ~ "moderate_preservation",
        TRUE ~ "weak_or_unpreserved"
      )
    }
  }
  sensitivity_status <- status_from_effects(effects)
  direction_stability <- if (!is.na(direction_col)) {
    ave(as.character(effects[[direction_col]]), effects[[id_col]], effects[[contrast_col]], FUN = function(x) {
      x <- unique(x[!is.na(x) & nzchar(x)])
      if (!length(x)) "not_available" else if (length(x) == 1L) "stable" else "mixed_direction"
    })
  } else {
    rep("not_available", nrow(effects))
  }
  confounding_status <- if (!is.null(status) && any(status$status %in% c("missing_optional", "read_error"))) "qc_incomplete" else "qc_inputs_available"
  gate <- effects$claim_gate == "eligible_for_readiness_assessment" &
    sensitivity_status == "effect_available" &
    preservation_status %in% c("not_available", "moderate_preservation", "strong_preservation") &
    direction_stability %in% c("stable", "not_available") &
    confounding_status == "qc_inputs_available"
  reason <- vapply(seq_len(nrow(effects)), function(i) {
    reasons <- c(
      if (!sensitivity_status[[i]] %in% "effect_available") sensitivity_status[[i]] else NULL,
      if (preservation_status[[i]] %in% "weak_or_unpreserved") "weak_or_unpreserved" else NULL,
      if (direction_stability[[i]] %in% "mixed_direction") "mixed_direction" else NULL,
      if (confounding_status != "qc_inputs_available") confounding_status else NULL
    )
    if (length(reasons)) paste(unique(reasons), collapse = ";") else "none"
  }, character(1))
  data.frame(
    dataset = ds,
    module_or_supermodule_id = as.character(effects[[id_col]]),
    level = level,
    contrast = if (!is.na(contrast_col)) as.character(effects[[contrast_col]]) else NA_character_,
    spatial_unit = if (!is.na(spatial_col)) as.character(effects[[spatial_col]]) else NA_character_,
    effect_scope = if (!is.na(scope_col)) as.character(effects[[scope_col]]) else NA_character_,
    robustness_gate = ifelse(gate, "pass", "missing_required"),
    preservation_status = preservation_status,
    sensitivity_status = sensitivity_status,
    direction_stability = direction_stability,
    confounding_status = confounding_status,
    claim_gate_eligible = gate,
    robustness_downgrade_reason = reason,
    source_file = as.character(effects$source_artifact),
    source_key = as.character(effects$source_key),
    tier_specific_family_id = as.character(
      effects$tier_specific_family_id
    ),
    tier_specific_family_size = as.integer(
      effects$tier_specific_family_size
    ),
    stringsAsFactors = FALSE
  )
}

make_dataset <- function(ds) {
  paths <- create_module_dirs("06_modules_WGCNA", file.path("module_robustness_sensitivity", ds))
  inputs <- list(
    inferential_handoff = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_inferential_handoff.csv"),
    preservation = path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds, "modules", "WGCNA_module_preservation_summary.csv"),
    pca_qc = path_results("tables", "03_qc_exploration", "05_pca_confounding_qc", ds, "PCA_confounding_summary.csv"),
    variance_qc = path_results("tables", "03_qc_exploration", "06_variance_partitioning", ds, "group_technical_confounding_screen.csv")
  )
  if (run$dry_run) {
    dry_run_inputs(paste(SCRIPT_ID, ds), inputs)
    handoff_ready <- file.exists(inputs$inferential_handoff)
    dry_run_line(
      "Required inferential handoff",
      inputs$inferential_handoff,
      if (handoff_ready) "PASS" else "FAIL"
    )
    dry_run_line("WGCNA robustness claim gate audit", path_results("reviewer_audit", "wgcna_robustness_claim_gate.csv"))
    return(handoff_ready)
  }
  loaded <- lapply(names(inputs), function(nm) read_csv_optional(inputs[[nm]], ds, "robustness_sensitivity", nm, required = FALSE))
  names(loaded) <- names(inputs)
  status <- do.call(rbind, lapply(loaded, `[[`, "status"))
  effect_data <- loaded$inferential_handoff$data
  if (!is.null(effect_data) && nrow(effect_data)) {
    wgcna_inferential_handoff_validate(effect_data)
  }
  module_effects <- if (is.null(effect_data)) NULL else
    effect_data[effect_data$entity_level == "module", , drop = FALSE]
  supermodule_effects <- if (is.null(effect_data)) NULL else
    effect_data[effect_data$entity_level == "supermodule", , drop = FALSE]
  effects <- module_effects
  if (is.null(effects) || !nrow(effects)) {
    out <- availability_evidence(ds, "module_robustness_sensitivity", inputs$inferential_handoff, "Module inferential handoff unavailable.")
  } else {
    mod_col <- "entity_id"
    fdr_col <- "tier_specific_fdr"
    p_col <- first_col(effects, c("p_value", "p.value"))
    est_col <- first_col(effects, c("estimate", "effect_size"))
    contrast_col <- first_col(effects, c("contrast"))
    spatial_col <- first_col(effects, c("spatial_unit", "SpatialLabel"))
    preservation <- loaded$preservation$data
    pres_key <- first_col(preservation, c("ModuleID", "module_id"))
    pres_col <- first_col(preservation, c("preservation_Zsummary_median", "Zsummary", "medianRank"))
    qc_flag <- ifelse(any(status$status %in% c("missing_optional", "read_error")), "WARN", "PASS")
    out <- data.frame(
      dataset = ds,
      evidence_domain = "module_robustness_sensitivity",
      evidence_id = as.character(effects$source_key),
      program_label = effects[[mod_col]],
      entity_type = "module",
      entity_id = effects[[mod_col]],
      contrast = if (!is.na(contrast_col)) effects[[contrast_col]] else NA_character_,
      spatial_unit = if (!is.na(spatial_col)) effects[[spatial_col]] else NA_character_,
      effect_direction = if (!is.na(est_col)) ifelse(num_or_na(effects[[est_col]]) > 0, "positive", "negative") else NA_character_,
      effect_size = if (!is.na(est_col)) num_or_na(effects[[est_col]]) else NA_real_,
      p_value = if (!is.na(p_col)) num_or_na(effects[[p_col]]) else NA_real_,
      fdr = if (!is.na(fdr_col)) num_or_na(effects[[fdr_col]]) else NA_real_,
      support_count = NA_real_,
      source_file = as.character(effects$source_artifact),
      evidence_status = NA_character_,
      interpretation_note = paste0(
        "Stage 07 inferential handoff row; family=",
        effects$tier_specific_family_id,
        "; family_size=", effects$tier_specific_family_size,
        "; claim_gate=", effects$claim_gate,
        "; source_key=", effects$source_key,
        "; preservation/QC evidence appended when available."
      ),
      qc_flag = qc_flag,
      stringsAsFactors = FALSE
    )
    if (!is.null(preservation) && nrow(preservation) && !is.na(pres_key) && !is.na(pres_col)) {
      out$support_count <- num_or_na(preservation[[pres_col]][match(out$entity_id, preservation[[pres_key]])])
      out$interpretation_note <- paste0(out$interpretation_note, " preservation_metric=", out$support_count)
    }
    out$evidence_status <- vapply(seq_len(nrow(out)), function(i) evidence_strength(out$fdr[[i]], out$support_count[[i]]), character(1))
    out <- standardize_evidence(out)
  }
  write_integration_table(out, paths, "module_robustness_sensitivity.csv")
  claim_gate_rows[[ds]] <<- dplyr::bind_rows(
    make_claim_gate_audit(ds, module_effects, "module", inputs$inferential_handoff, loaded$preservation$data, status),
    make_claim_gate_audit(ds, supermodule_effects, "supermodule", inputs$inferential_handoff, loaded$preservation$data, status)
  )
  write_csv_safe(status, file.path(paths$reports, "input_status.csv"))
  write_csv_safe(status, file.path(paths$source_data, "module_robustness_sensitivity_input_status.csv"))
  write_integration_manifest(paths, inputs, list(tables = paths$tables, source_data = paths$source_data), list(dataset = ds), "Robustness/sensitivity synthesis from the validated Stage 07 inferential handoff, preservation summaries, and QC/confounding screens. Tier-specific FDR family identity, claim gate, and exact Stage 05 source key are retained.")
  out
}

dataset_results <- lapply(integration_datasets(run$dataset), make_dataset)
if (run$dry_run) {
  required_ready <- all(unlist(dataset_results, use.names = FALSE))
  dry_run_line(
    "Status",
    if (required_ready) "ready" else "missing_required_input",
    if (required_ready) "PASS" else "FAIL"
  )
  quit(status = if (required_ready) 0 else 1, save = "no")
}
invisible(dataset_results)
audit <- dplyr::bind_rows(claim_gate_rows)
dir_create(path_results("reviewer_audit"))
write_csv_safe(audit, path_results("reviewer_audit", "wgcna_robustness_claim_gate.csv"))
message("Module robustness/sensitivity complete.")
