#!/usr/bin/env Rscript
# ================================================================
# Script: 08_behavior_physio_coupling/03_module_behavior_coupling.r
# Stage: coupling
# Scope: dataset_specific
# Consumes: module activity/effect tables and behavior/network coupling summaries.
# Produces: module-behavior coupling evidence.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))

SCRIPT_ID <- "08_behavior_physio_coupling/03_module_behavior_coupling.r"
run <- integration_cli(allow_all = TRUE)

make_dataset <- function(ds) {
  paths <- create_module_dirs("08_behavior_physio_coupling", file.path("module_behavior_coupling", ds))
  inputs <- list(
    module_effects = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"),
    interpretable = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_module_group_effects_interpretable.csv"),
    network_behavior = path_results("tables", "08_behavior_physio_coupling", "network_behavior_coupling", "edge_behavior_figure_ready_table.csv")
  )
  if (run$dry_run) {
    dry_run_inputs(paste(SCRIPT_ID, ds), inputs)
    return(NULL)
  }
  loaded <- lapply(names(inputs), function(nm) read_csv_optional(inputs[[nm]], ds, "module_behavior_coupling", nm, required = FALSE))
  names(loaded) <- names(inputs)
  status <- do.call(rbind, lapply(loaded, `[[`, "status"))
  modules <- loaded$interpretable$data %||% loaded$module_effects$data
  behavior <- loaded$network_behavior$data
  if (is.null(modules) || !nrow(modules)) {
    out <- availability_evidence(ds, "module_behavior_coupling", inputs$interpretable, "Module evidence unavailable.")
  } else {
    mod_col <- first_col(modules, c("module_id", "ModuleID"))
    prog_col <- first_col(modules, c("ModuleLabel_Final", "module_label", "Supermodule_PlotLabel", "Macroprogram_Display", "module_id"))
    fdr_col <- first_col(modules, c("FDR_global", "FDR_within_dataset_level"))
    est_col <- first_col(modules, c("estimate"))
    if (is.null(behavior) || !nrow(behavior)) {
      out <- availability_evidence(ds, "module_behavior_coupling", inputs$network_behavior, "Behavior/network coupling table unavailable; module rows are not behavior-coupled.")
      if (!is.na(mod_col)) {
        out <- rbind(out, standardize_evidence(data.frame(
          dataset = ds,
          evidence_domain = "module_behavior_coupling",
          evidence_id = paste(ds, modules[[mod_col]], "behavior_unavailable", sep = "::"),
          program_label = if (!is.na(prog_col)) modules[[prog_col]] else modules[[mod_col]],
          entity_type = "module",
          entity_id = modules[[mod_col]],
          effect_size = if (!is.na(est_col)) num_or_na(modules[[est_col]]) else NA_real_,
          fdr = if (!is.na(fdr_col)) num_or_na(modules[[fdr_col]]) else NA_real_,
          source_file = inputs$interpretable,
          evidence_status = "module_effect_without_behavior_input",
          interpretation_note = "Module evidence present; behavior coupling optional input missing.",
          qc_flag = "WARN",
          stringsAsFactors = FALSE
        )))
      }
    } else {
      outcome_col <- first_col(behavior, c("Outcome", "biological_program"))
      fdr_b <- first_col(behavior, c("fdr", "p.adj_BH_all_edge_phenotype_tests", "FDR"))
      est_b <- first_col(behavior, c("estimate", "effect_size_NES"))
      top_behavior <- behavior[order(num_or_na(behavior[[fdr_b]]), decreasing = FALSE, na.last = TRUE), , drop = FALSE]
      top_behavior <- utils::head(top_behavior, 10)
      module_rows <- utils::head(modules, 200)
      out <- standardize_evidence(data.frame(
        dataset = ds,
        evidence_domain = "module_behavior_coupling",
        evidence_id = paste(ds, module_rows[[mod_col]], "behavior_context", sep = "::"),
        program_label = if (!is.na(prog_col)) module_rows[[prog_col]] else module_rows[[mod_col]],
        entity_type = "module",
        entity_id = module_rows[[mod_col]],
        effect_size = if (!is.na(est_col)) num_or_na(module_rows[[est_col]]) else NA_real_,
        fdr = if (!is.na(fdr_col)) num_or_na(module_rows[[fdr_col]]) else NA_real_,
        support_count = nrow(top_behavior),
        source_file = paste(inputs$interpretable, inputs$network_behavior, sep = ";"),
        evidence_status = ifelse(nrow(top_behavior) > 0, "behavior_context_available", "behavior_context_empty"),
        interpretation_note = paste0("Top available behavior outcomes: ", paste(unique(as.character(top_behavior[[outcome_col]])), collapse = ";")),
        qc_flag = "PASS",
        stringsAsFactors = FALSE
      ))
    }
  }
  write_integration_table(out, paths, "module_behavior_coupling.csv")
  write_csv_safe(status, file.path(paths$reports, "input_status.csv"))
  write_csv_safe(status, file.path(paths$source_data, "module_behavior_coupling_input_status.csv"))
  write_integration_manifest(paths, inputs, list(tables = paths$tables, source_data = paths$source_data), list(dataset = ds), "Module-behavior coupling synthesis. Missing behavior inputs are reported as status/evidence rows.")
  out
}

invisible(lapply(integration_datasets(run$dataset), make_dataset))
if (run$dry_run) quit(status = 0, save = "no")
message("Module behavior coupling complete.")
