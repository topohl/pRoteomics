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
source(repo_path("R", "wgcna_group_effect_consumer_utils.R"))

SCRIPT_ID <- "08_behavior_physio_coupling/03_module_behavior_coupling.r"
run <- integration_cli(allow_all = TRUE)

make_dataset <- function(ds) {
  paths <- create_module_dirs("08_behavior_physio_coupling", file.path("module_behavior_coupling", ds))
  inputs <- list(
    inferential_handoff = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_inferential_handoff.csv"),
    network_behavior = path_results("tables", "08_behavior_physio_coupling", "network_behavior_coupling", "edge_behavior_figure_ready_table.csv")
  )
  if (run$dry_run) {
    dry_run_inputs(paste(SCRIPT_ID, ds), inputs)
    handoff_ready <- file.exists(inputs$inferential_handoff)
    dry_run_line(
      "Required inferential handoff",
      inputs$inferential_handoff,
      if (handoff_ready) "PASS" else "FAIL"
    )
    return(handoff_ready)
  }
  loaded <- lapply(names(inputs), function(nm) read_csv_optional(inputs[[nm]], ds, "module_behavior_coupling", nm, required = FALSE))
  names(loaded) <- names(inputs)
  status <- do.call(rbind, lapply(loaded, `[[`, "status"))
  modules <- loaded$inferential_handoff$data
  if (!is.null(modules) && nrow(modules)) {
    wgcna_inferential_handoff_validate(modules)
    modules <- modules[
      modules$entity_level == "module", , drop = FALSE
    ]
  }
  behavior <- loaded$network_behavior$data
  if (is.null(modules) || !nrow(modules)) {
    out <- availability_evidence(ds, "module_behavior_coupling", inputs$inferential_handoff, "Module inferential handoff unavailable.")
  } else {
    mod_col <- "entity_id"
    prog_col <- "display_label"
    fdr_col <- "tier_specific_fdr"
    est_col <- first_col(modules, c("estimate"))
    if (is.null(behavior) || !nrow(behavior)) {
      out <- availability_evidence(ds, "module_behavior_coupling", inputs$network_behavior, "Behavior/network coupling table unavailable; module rows are not behavior-coupled.")
      if (!is.na(mod_col)) {
        out <- rbind(out, standardize_evidence(data.frame(
          dataset = ds,
          evidence_domain = "module_behavior_coupling",
          evidence_id = paste(
            modules$source_key, "behavior_unavailable", sep = "::"
          ),
          program_label = if (!is.na(prog_col)) modules[[prog_col]] else modules[[mod_col]],
          entity_type = "module",
          entity_id = modules[[mod_col]],
          effect_size = if (!is.na(est_col)) num_or_na(modules[[est_col]]) else NA_real_,
          fdr = if (!is.na(fdr_col)) num_or_na(modules[[fdr_col]]) else NA_real_,
          source_file = modules$source_artifact,
          evidence_status = "module_effect_without_behavior_input",
          interpretation_note = paste0(
            "Module inference present; optional behavior input missing; ",
            "claim_gate=", modules$claim_gate,
            "; family=", modules$tier_specific_family_id,
            "; family_size=", modules$tier_specific_family_size,
            "; source_key=", modules$source_key
          ),
          qc_flag = ifelse(
            modules$claim_gate ==
              "eligible_for_readiness_assessment", "WARN", "DIAGNOSTIC"
          ),
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
        evidence_id = paste(
          module_rows$source_key, "behavior_context", sep = "::"
        ),
        program_label = if (!is.na(prog_col)) module_rows[[prog_col]] else module_rows[[mod_col]],
        entity_type = "module",
        entity_id = module_rows[[mod_col]],
        effect_size = if (!is.na(est_col)) num_or_na(module_rows[[est_col]]) else NA_real_,
        fdr = if (!is.na(fdr_col)) num_or_na(module_rows[[fdr_col]]) else NA_real_,
        support_count = nrow(top_behavior),
        source_file = paste(
          module_rows$source_artifact,
          inputs$network_behavior, sep = ";"
        ),
        evidence_status = ifelse(
          module_rows$claim_gate ==
            "eligible_for_readiness_assessment",
          ifelse(
            nrow(top_behavior) > 0,
            "behavior_context_available", "behavior_context_empty"
          ),
          "diagnostic_module_effect_behavior_context"
        ),
        interpretation_note = paste0(
          "Top available behavior outcomes: ",
          paste(unique(as.character(top_behavior[[outcome_col]])),
                collapse = ";"),
          "; claim_gate=", module_rows$claim_gate,
          "; family=", module_rows$tier_specific_family_id,
          "; family_size=", module_rows$tier_specific_family_size,
          "; source_key=", module_rows$source_key
        ),
        qc_flag = ifelse(
          module_rows$claim_gate ==
            "eligible_for_readiness_assessment", "PASS", "DIAGNOSTIC"
        ),
        stringsAsFactors = FALSE
      ))
    }
  }
  write_integration_table(out, paths, "module_behavior_coupling.csv")
  write_csv_safe(status, file.path(paths$reports, "input_status.csv"))
  write_csv_safe(status, file.path(paths$source_data, "module_behavior_coupling_input_status.csv"))
  write_integration_manifest(paths, inputs, list(tables = paths$tables, source_data = paths$source_data), list(dataset = ds), "Module-behavior coupling synthesis. WGCNA inference comes only from the validated Stage 07 handoff; tier-specific family identity, claim gate, and exact Stage 05 source key are retained. Missing behavior inputs are reported as status/evidence rows.")
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
message("Module behavior coupling complete.")
