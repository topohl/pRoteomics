#!/usr/bin/env Rscript
# ================================================================
# Script: 10_biological_integration/03_evidence_priority_matrix.r
# Stage: integration
# Scope: global
# Consumes: cross-compartment atlas, manuscript summary, and biological claims table when present.
# Produces: evidence priority matrix.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "final_evidence_bundle_utils.R"))

SCRIPT_ID <- "10_biological_integration/03_evidence_priority_matrix.r"
run <- integration_cli(default_dataset = "all", allow_all = TRUE)
paths <- integration_paths("evidence_priority_matrix", "global")
inputs <- list(
  atlas = path_results("tables", "10_biological_integration", "cross_compartment_program_atlas", "global", "cross_compartment_program_atlas_long.csv"),
  manuscript_summary = path_results("tables", "10_biological_integration", "manuscript_program_summary", "global", "manuscript_program_summary.csv"),
  claims = path_results("tables", "biological_claims_table.csv")
)

if (run$dry_run) {
  dry_run_inputs(SCRIPT_ID, inputs)
  quit(status = 0, save = "no")
}

loaded <- lapply(names(inputs), function(nm) read_csv_optional(inputs[[nm]], "global", "integration_priority", nm, required = FALSE))
names(loaded) <- names(inputs)
status <- do.call(rbind, lapply(loaded, `[[`, "status"))
atlas <- loaded$atlas$data
summary <- loaded$manuscript_summary$data

if (is.null(atlas) || !nrow(atlas)) {
  priority <- data.frame(
    priority_id = "PRIORITY_0001",
    program_key = "Unavailable optional evidence",
    dataset = NA_character_,
    priority_tier = "defer",
    evidence_domain_count = 0L,
    strongest_fdr = NA_real_,
    robustness_flag = "unavailable",
    behavior_flag = "unavailable",
    qc_flag = "WARN",
    recommended_use = "Do not use for biological claims until integration inputs are available.",
    stringsAsFactors = FALSE
  )
} else {
  if (!"program_key" %in% names(atlas)) atlas$program_key <- program_key(atlas$program_label)
  groups <- split(atlas, paste(atlas$dataset, atlas$program_key, sep = "||"))
  priority <- do.call(rbind, lapply(seq_along(groups), function(i) {
    d <- groups[[i]]
    domains <- unique(d$evidence_domain)
    fdr <- suppressWarnings(min(num_or_na(d$fdr), na.rm = TRUE))
    if (!is.finite(fdr)) fdr <- NA_real_
    has_robust <- any(grepl("robustness", domains))
    has_behavior <- any(grepl("behavior", domains))
    qc <- if (any(d$qc_flag %in% c("FAIL", "WARN"), na.rm = TRUE)) "WARN" else "PASS"
    tier <- if (length(domains) >= 5 && (is.na(fdr) || fdr <= 0.10) && qc == "PASS") {
      "tier_1_manuscript_anchor"
    } else if (length(domains) >= 3) {
      "tier_2_supporting"
    } else if (length(domains) >= 2) {
      "tier_3_context"
    } else {
      "defer"
    }
    data.frame(
      priority_id = sprintf("PRIORITY_%04d", i),
      program_key = d$program_key[[1]],
      dataset = d$dataset[[1]],
      priority_tier = tier,
      evidence_domain_count = length(domains),
      evidence_domains = paste(sort(domains), collapse = ";"),
      strongest_fdr = fdr,
      robustness_flag = ifelse(has_robust, "available", "missing_or_not_applicable"),
      behavior_flag = ifelse(has_behavior, "available", "missing_or_not_applicable"),
      qc_flag = qc,
      recommended_use = ifelse(tier == "defer", "Exploratory context only.", "Candidate manuscript synthesis row; verify source evidence and limitations."),
      stringsAsFactors = FALSE
    )
  }))
}

if (!is.null(summary) && nrow(summary) && "program_key" %in% names(summary)) {
  priority$manuscript_claim_scope <- summary$manuscript_claim_scope[match(priority$program_key, summary$program_key)]
}

validate_evidence_priority_matrix(priority)
invisible(write_integration_table(priority, paths, "evidence_priority_matrix.csv"))
write_csv_safe(status, file.path(paths$reports, "input_status.csv"))
write_csv_safe(status, file.path(paths$source_data, "evidence_priority_matrix_input_status.csv"))
write_integration_manifest(paths, inputs, list(tables = paths$tables, source_data = paths$source_data), list(dataset = run$dataset), "Priority matrix ranks program evidence by convergence, robustness, behavior coupling, and QC flags.")
bundle <- write_final_evidence_bundle(reason = "integration_priority_matrix")
message("Final biological evidence bundle refreshed: ", bundle$bundle_dir)
message("Evidence priority matrix complete: ", paths$tables)
