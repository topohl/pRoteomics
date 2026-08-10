#!/usr/bin/env Rscript
# ================================================================
# Script: 10_biological_integration/02_manuscript_program_summary.r
# Stage: integration
# Scope: global
# Consumes: cross-compartment program atlas.
# Produces: manuscript program summary.
# Dataset behavior: global summary across the canonical integration atlas.
# Notes: Counts distinct evidence_source_family lineages among eligible evidence rows.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))

SCRIPT_ID <- "10_biological_integration/02_manuscript_program_summary.r"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
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

loaded_long <- read_csv_optional(inputs$atlas_long, "global", "integration", "atlas_long", required = TRUE)
loaded_summary <- read_csv_optional(inputs$atlas_summary, "global", "integration", "atlas_summary", required = FALSE)
status <- rbind(loaded_long$status, loaded_summary$status)
atlas <- loaded_long$data

if (is.null(atlas) || !nrow(atlas)) {
  stop("Required cross-compartment program atlas is unavailable or empty: ", inputs$atlas_long, call. = FALSE)
} else {
  if (!"program_key" %in% names(atlas)) atlas$program_key <- program_key(atlas$program_label)
  split_prog <- split(atlas, atlas$program_key)
  summary <- do.call(rbind, lapply(names(split_prog), function(pk) {
    d <- split_prog[[pk]]
    counting <- d[d$counts_toward_convergence %in% TRUE, , drop = FALSE]
    domains <- sort(unique(counting$evidence_domain[!is.na(counting$evidence_domain)]))
    families <- independent_evidence_families(counting)
    datasets <- sort(unique(counting$dataset[counting$dataset %in% valid_datasets()]))
    fdr <- suppressWarnings(min(num_or_na(counting$fdr), na.rm = TRUE))
    if (!is.finite(fdr)) fdr <- NA_real_
    qc <- if (nrow(counting) && any(counting$qc_flag %in% c("FAIL", "WARN"), na.rm = TRUE)) "WARN" else if (nrow(counting)) "PASS" else "WARN"
    semantic <- as.character(counting$evidence_semantic_class)
    has_architecture <- any(semantic == "wgcna_architecture", na.rm = TRUE)
    has_stress_effect <- any(semantic == "wgcna_stress_group_effect", na.rm = TRUE)
    has_non_wgcna <- any(semantic == "non_wgcna_evidence", na.rm = TRUE)
    semantic_scope <- dplyr::case_when(
      !nrow(counting) ~ "no_claim",
      length(families) >= 2L && !is.na(fdr) && fdr <= 0.05 ~ "convergent_stress_associated_evidence",
      has_architecture && !has_stress_effect && !has_non_wgcna ~ "wgcna_architecture_covariance_context",
      length(families) >= 2L ~ "supporting_biological_context",
      qc == "WARN" ~ "exploratory_context",
      TRUE ~ "supporting_biological_context"
    )
    scope <- dplyr::case_when(
      semantic_scope == "convergent_stress_associated_evidence" && length(datasets) >= 2L ~ "primary_systems_claim",
      semantic_scope %in% c("convergent_stress_associated_evidence", "supporting_biological_context") ~ "supporting_claim",
      semantic_scope == "wgcna_architecture_covariance_context" ~ "architecture_context",
      semantic_scope == "exploratory_context" ~ "exploratory_context",
      TRUE ~ "no_claim"
    )
    sentence <- dplyr::case_when(
      semantic_scope == "convergent_stress_associated_evidence" ~ paste0("Convergent stress-associated evidence supports a ", pk, " program across ", paste(datasets, collapse = ", "), "."),
      semantic_scope == "wgcna_architecture_covariance_context" ~ paste0("WGCNA provides ", pk, " covariance and spatial-organization context; it does not establish stress-dependent regulation or network remodelling."),
      semantic_scope == "supporting_biological_context" ~ paste0("Evidence provides supporting biological context for a ", pk, " program; verify the individual statistical endpoints before stress-effect wording."),
      semantic_scope == "exploratory_context" ~ paste0("Evidence provides exploratory context for a ", pk, " program and is not sufficient for an independent manuscript claim."),
      TRUE ~ paste0("No independently countable evidence supports a manuscript claim for the ", pk, " program.")
    )
    data.frame(
      program_key = pk,
      manuscript_claim_scope = scope,
      datasets_supported = paste(datasets, collapse = ";"),
      evidence_domains = paste(domains, collapse = ";"),
      evidence_source_families = paste(families, collapse = ";"),
      n_independent_evidence_families = length(families),
      n_evidence_rows_total = nrow(d),
      n_evidence_rows_counting_toward_convergence = nrow(counting),
      n_wgcna_architecture_rows = sum(d$evidence_semantic_class == "wgcna_architecture", na.rm = TRUE),
      n_wgcna_stress_effect_rows = sum(d$evidence_semantic_class == "wgcna_stress_group_effect", na.rm = TRUE),
      n_wgcna_alias_rows_excluded = sum(d$claim_entity_role == "compatibility_alias" & !(d$counts_toward_convergence %in% TRUE), na.rm = TRUE),
      claim_semantic_scope = semantic_scope,
      strongest_fdr = fdr,
      strongest_evidence = evidence_strength(fdr, length(families)),
      safe_manuscript_sentence = sentence,
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
