#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/09_module_robustness_sensitivity.r
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: WGCNA group effects, preservation summaries, and QC/confounding outputs.
# Produces: module robustness/sensitivity evidence.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))

SCRIPT_ID <- "06_modules_WGCNA/09_module_robustness_sensitivity.r"
run <- integration_cli(allow_all = TRUE)

make_dataset <- function(ds) {
  paths <- create_module_dirs("06_modules_WGCNA", file.path("module_robustness_sensitivity", ds))
  inputs <- list(
    module_effects = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"),
    supermodule_effects = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"),
    preservation = path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds, "modules", "WGCNA_module_preservation_summary.csv"),
    pca_qc = path_results("tables", "03_qc_exploration", "05_pca_confounding_qc", ds, "PCA_confounding_summary.csv"),
    variance_qc = path_results("tables", "03_qc_exploration", "06_variance_partitioning", ds, "group_technical_confounding_screen.csv")
  )
  if (run$dry_run) {
    dry_run_inputs(paste(SCRIPT_ID, ds), inputs)
    return(NULL)
  }
  loaded <- lapply(names(inputs), function(nm) read_csv_optional(inputs[[nm]], ds, "robustness_sensitivity", nm, required = FALSE))
  names(loaded) <- names(inputs)
  status <- do.call(rbind, lapply(loaded, `[[`, "status"))
  effects <- loaded$module_effects$data
  if (is.null(effects) || !nrow(effects)) {
    out <- availability_evidence(ds, "module_robustness_sensitivity", inputs$module_effects, "Module group effects unavailable.")
  } else {
    mod_col <- first_col(effects, c("module_id", "ModuleID"))
    fdr_col <- first_col(effects, c("FDR_global", "FDR_within_dataset_level", "q_value"))
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
      evidence_id = paste(ds, effects[[mod_col]], seq_len(nrow(effects)), sep = "::"),
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
      source_file = inputs$module_effects,
      evidence_status = NA_character_,
      interpretation_note = "Group-effect robustness row; preservation/QC evidence appended when available.",
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
  write_csv_safe(status, file.path(paths$reports, "input_status.csv"))
  write_csv_safe(status, file.path(paths$source_data, "module_robustness_sensitivity_input_status.csv"))
  write_integration_manifest(paths, inputs, list(tables = paths$tables, source_data = paths$source_data), list(dataset = ds), "Robustness/sensitivity synthesis from existing module effects, preservation summaries, and QC/confounding screens.")
  out
}

invisible(lapply(integration_datasets(run$dataset), make_dataset))
if (run$dry_run) quit(status = 0, save = "no")
message("Module robustness/sensitivity complete.")
