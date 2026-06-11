#!/usr/bin/env Rscript
# ================================================================
# Script: 10_biological_integration/02_manuscript_program_summary.r
# Stage: integration
# Scope: global
# Consumes: cross-compartment program atlas.
# Produces: manuscript program summary.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))

SCRIPT_ID <- "10_biological_integration/02_manuscript_program_summary.r"
run <- integration_cli(default_dataset = "all", allow_all = TRUE)
paths <- integration_paths("manuscript_program_summary", "global")
inputs <- list(
  atlas_long = path_results("tables", "10_biological_integration", "cross_compartment_program_atlas", "global", "cross_compartment_program_atlas_long.csv"),
  atlas_summary = path_results("tables", "10_biological_integration", "cross_compartment_program_atlas", "global", "cross_compartment_program_atlas.csv")
)

if (run$dry_run) {
  dry_run_inputs(SCRIPT_ID, inputs)
  quit(status = 0, save = "no")
}

loaded_long <- read_csv_optional(inputs$atlas_long, "global", "integration", "atlas_long", required = FALSE)
loaded_summary <- read_csv_optional(inputs$atlas_summary, "global", "integration", "atlas_summary", required = FALSE)
status <- rbind(loaded_long$status, loaded_summary$status)
atlas <- loaded_long$data

if (is.null(atlas) || !nrow(atlas)) {
  summary <- data.frame(
    program_key = "Unavailable optional evidence",
    manuscript_claim_scope = "no_claim",
    datasets_supported = NA_character_,
    evidence_domains = NA_character_,
    strongest_evidence = "unavailable",
    safe_manuscript_sentence = "Integration atlas was unavailable; no manuscript-level biological synthesis should be made from this table.",
    main_limitation = "Missing optional integration inputs.",
    qc_flag = "WARN",
    stringsAsFactors = FALSE
  )
} else {
  if (!"program_key" %in% names(atlas)) atlas$program_key <- program_key(atlas$program_label)
  split_prog <- split(atlas, atlas$program_key)
  summary <- do.call(rbind, lapply(names(split_prog), function(pk) {
    d <- split_prog[[pk]]
    domains <- sort(unique(d$evidence_domain[!is.na(d$evidence_domain)]))
    datasets <- sort(unique(d$dataset[d$dataset %in% valid_datasets()]))
    fdr <- suppressWarnings(min(num_or_na(d$fdr), na.rm = TRUE))
    if (!is.finite(fdr)) fdr <- NA_real_
    qc <- if (any(d$qc_flag %in% c("FAIL", "WARN"), na.rm = TRUE)) "WARN" else "PASS"
    data.frame(
      program_key = pk,
      manuscript_claim_scope = ifelse(length(domains) >= 4 && length(datasets) >= 2, "primary_systems_claim", ifelse(length(domains) >= 2, "supporting_claim", "exploratory_context")),
      datasets_supported = paste(datasets, collapse = ";"),
      evidence_domains = paste(domains, collapse = ";"),
      strongest_fdr = fdr,
      strongest_evidence = evidence_strength(fdr, length(domains)),
      safe_manuscript_sentence = paste0("Evidence supports a ", pk, " program across ", paste(datasets, collapse = ", "), " with ", length(domains), " evidence stream(s)."),
      main_limitation = ifelse(qc == "WARN", "At least one supporting input is missing or carries a QC warning; review source rows before strong wording.", "Observational proteomics; causal direction is not established."),
      qc_flag = qc,
      stringsAsFactors = FALSE
    )
  }))
}

validate_manuscript_program_summary(summary)
invisible(write_integration_table(summary, paths, "manuscript_program_summary.csv"))
write_csv_safe(status, file.path(paths$reports, "input_status.csv"))
write_csv_safe(status, file.path(paths$source_data, "manuscript_program_summary_input_status.csv"))
write_integration_manifest(paths, inputs, list(tables = paths$tables, source_data = paths$source_data), list(dataset = run$dataset), "Manuscript-facing biological program summary built only from integration atlas outputs.")
message("Manuscript program summary complete: ", paths$tables)
