#!/usr/bin/env Rscript
# ================================================================
# Script: 10_biological_integration/01_cross_compartment_program_atlas.r
# Stage: integration
# Scope: global
# Consumes: enrichment, WGCNA, microenvironment, complex/organelle, robustness,
#           spatial architecture, behavior-coupling, and QC evidence tables.
# Produces: cross-compartment program atlas.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))

SCRIPT_ID <- "10_biological_integration/01_cross_compartment_program_atlas.r"
run <- integration_cli(default_dataset = "all", allow_all = TRUE)
paths <- integration_paths("cross_compartment_program_atlas", "global")

dataset_inputs <- function(ds) {
  list(
    enrichment_program = path_results("tables", "04_differential_expression_enrichment", "biological_program_summary", ds, "program_summary.csv"),
    external_signature = path_results("tables", "04_differential_expression_enrichment", "external_stress_disease_signature_overlap", "global", "external_stress_disease_signature_overlap.csv"),
    wgcna_interpretable = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_supermodule_group_effects_interpretable.csv"),
    microenvironment = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"),
    complex_architecture = path_results("tables", "06_modules_WGCNA", "module_complex_architecture", ds, "module_complex_architecture.csv"),
    robustness = path_results("tables", "06_modules_WGCNA", "module_robustness_sensitivity", ds, "module_robustness_sensitivity.csv"),
    module_behavior = path_results("tables", "08_behavior_physio_coupling", "module_behavior_coupling", ds, "module_behavior_coupling.csv"),
    qc_report = path_results("reports", "03_qc_exploration", "07_qc_biology_confounding_report", ds, "qc_biology_confounding_summary.md"),
    spatial_program = path_results("tables", "04_differential_expression_enrichment", "compareGO_spatial_atlas", "spatial_program_summary.csv")
  )
}

all_inputs <- unlist(lapply(integration_datasets(run$dataset), dataset_inputs), use.names = TRUE)
if (run$dry_run) {
  dry_run_inputs(SCRIPT_ID, all_inputs)
  quit(status = 0, save = "no")
}

program_evidence <- function(ds, file) {
  loaded <- read_csv_optional(file, ds, "enrichment_program", "program_summary", required = FALSE)
  df <- loaded$data
  if (is.null(df) || !nrow(df)) return(list(evidence = availability_evidence(ds, "enrichment_program", file, "Program summary unavailable."), status = loaded$status))
  prog <- first_col(df, c("biological_program", "program", "program_class", "top_term", "Description"))
  contrast <- first_col(df, c("comparison", "contrast", "phenotype_contrast"))
  spatial <- first_col(df, c("route_unit", "spatial_unit", "region"))
  est <- first_col(df, c("representative_NES", "mean_NES", "NES"))
  fdr <- first_col(df, c("min_fdr", "FDR", "p.adjust", "fdr"))
  p <- first_col(df, c("min_raw_p", "pvalue", "p.value", "p"))
  genes <- first_col(df, c("key_genes", "core_genes", "top_driver_genes"))
  ev <- standardize_evidence(data.frame(
    dataset = ds, evidence_domain = "enrichment_program",
    evidence_id = paste(ds, "enrichment", seq_len(nrow(df)), sep = "::"),
    program_label = if (!is.na(prog)) df[[prog]] else "Unlabelled enrichment program",
    entity_type = "program", entity_id = if (!is.na(prog)) df[[prog]] else NA_character_,
    contrast = if (!is.na(contrast)) df[[contrast]] else NA_character_,
    spatial_unit = if (!is.na(spatial)) df[[spatial]] else NA_character_,
    effect_size = if (!is.na(est)) num_or_na(df[[est]]) else NA_real_,
    p_value = if (!is.na(p)) num_or_na(df[[p]]) else NA_real_,
    fdr = if (!is.na(fdr)) num_or_na(df[[fdr]]) else NA_real_,
    support_count = if (!is.na(genes)) lengths(lapply(df[[genes]], split_gene_tokens)) else NA_real_,
    source_file = file,
    evidence_status = "program_evidence",
    interpretation_note = "Differential enrichment/program evidence.",
    qc_flag = "PASS",
    stringsAsFactors = FALSE
  ))
  list(evidence = ev, status = loaded$status)
}

evidence_file <- function(ds, domain, file, input_type) {
  loaded <- read_csv_optional(file, ds, domain, input_type, required = FALSE)
  df <- loaded$data
  if (is.null(df) || !nrow(df)) return(list(evidence = availability_evidence(ds, domain, file, paste(input_type, "unavailable.")), status = loaded$status))
  if (all(names(empty_evidence()) %in% names(df))) {
    return(list(evidence = standardize_evidence(df), status = loaded$status))
  }
  prog <- first_col(df, c("program_label", "Supermodule_PlotLabel", "Supermodule_DisplayLabel", "Macroprogram_Display", "biological_program"))
  id <- first_col(df, c("supermodule_id", "SupermoduleID", "module_id", "ModuleID"))
  fdr <- first_col(df, c("FDR_global", "fdr", "FDR"))
  est <- first_col(df, c("estimate", "effect_size"))
  ev <- standardize_evidence(data.frame(
    dataset = ds, evidence_domain = domain,
    evidence_id = paste(ds, domain, seq_len(nrow(df)), sep = "::"),
    program_label = if (!is.na(prog)) df[[prog]] else domain,
    entity_type = ifelse(grepl("module", domain), "module_or_supermodule", "program"),
    entity_id = if (!is.na(id)) df[[id]] else NA_character_,
    effect_size = if (!is.na(est)) num_or_na(df[[est]]) else NA_real_,
    fdr = if (!is.na(fdr)) num_or_na(df[[fdr]]) else NA_real_,
    source_file = file,
    evidence_status = domain,
    interpretation_note = paste(domain, "evidence imported for atlas."),
    qc_flag = "PASS",
    stringsAsFactors = FALSE
  ))
  list(evidence = ev, status = loaded$status)
}

all_ev <- list()
all_status <- empty_status()
for (ds in integration_datasets(run$dataset)) {
  inputs <- dataset_inputs(ds)
  pieces <- list(
    program_evidence(ds, inputs$enrichment_program),
    evidence_file(ds, "external_signature_overlap", inputs$external_signature, "external_signature_overlap"),
    evidence_file(ds, "wgcna_supermodule", inputs$wgcna_interpretable, "wgcna_interpretable"),
    evidence_file(ds, "microenvironment_marker", inputs$microenvironment, "microenvironment"),
    evidence_file(ds, "complex_organelle", inputs$complex_architecture, "complex_architecture"),
    evidence_file(ds, "robustness_sensitivity", inputs$robustness, "robustness"),
    evidence_file(ds, "behavior_coupling", inputs$module_behavior, "module_behavior")
  )
  all_ev <- c(all_ev, lapply(pieces, `[[`, "evidence"))
  all_status <- rbind(all_status, do.call(rbind, lapply(pieces, `[[`, "status")))
}

atlas_long <- standardize_evidence(do.call(rbind, all_ev))
atlas_long$program_key <- program_key(atlas_long$program_label)
atlas_long$dataset_label <- vapply(atlas_long$dataset, function(x) if (x %in% valid_datasets()) dataset_contracts()[[x]]$label else x, character(1))
atlas_long$evidence_strength <- vapply(seq_len(nrow(atlas_long)), function(i) evidence_strength(atlas_long$fdr[[i]], atlas_long$support_count[[i]]), character(1))

summary_rows <- aggregate(
  evidence_domain ~ dataset + program_key,
  atlas_long,
  function(x) paste(sort(unique(x)), collapse = ";")
)
names(summary_rows)[names(summary_rows) == "evidence_domain"] <- "evidence_domains"
summary_rows$n_evidence_domains <- lengths(strsplit(summary_rows$evidence_domains, ";", fixed = TRUE))
fdr_summary <- aggregate(fdr ~ dataset + program_key, atlas_long, function(x) suppressWarnings(min(num_or_na(x), na.rm = TRUE)))
names(fdr_summary)[names(fdr_summary) == "fdr"] <- "strongest_fdr"
summary_rows <- merge(summary_rows, fdr_summary, by = c("dataset", "program_key"), all.x = TRUE)
summary_rows$strongest_fdr[is.infinite(summary_rows$strongest_fdr)] <- NA_real_
summary_rows$atlas_status <- ifelse(summary_rows$n_evidence_domains >= 4, "multi_stream_supported", ifelse(summary_rows$n_evidence_domains >= 2, "convergent_support", "single_stream_or_missing"))

validate_cross_compartment_program_atlas(atlas_long)
invisible(write_integration_table(atlas_long, paths, "cross_compartment_program_atlas_long.csv"))
invisible(write_integration_table(summary_rows, paths, "cross_compartment_program_atlas.csv"))
write_csv_safe(all_status, file.path(paths$reports, "input_status.csv"))
write_csv_safe(all_status, file.path(paths$source_data, "cross_compartment_program_atlas_input_status.csv"))
write_integration_manifest(paths, as.list(all_inputs), list(tables = paths$tables, source_data = paths$source_data), list(dataset = run$dataset), "Manuscript-level atlas integrating enrichment, modules, microenvironment, complexes/organelles, robustness, behavior, and QC flags from existing outputs.")
message("Cross-compartment program atlas complete: ", paths$tables)
