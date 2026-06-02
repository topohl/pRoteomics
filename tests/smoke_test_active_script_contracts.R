#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

active_scripts <- c(
  "03_qc_exploration/00_dataset_qc_report.r",
  "04_differential_expression_enrichment/05_microglia_targeted_signature_enrichment.r",
  "04_differential_expression_enrichment/03_biological_program_summary.r",
  "06_modules_WGCNA/01_WGCNA.r",
  "06_modules_WGCNA/02_curated_overlap_programs.r",
  "06_modules_WGCNA/03_score_module_activity.R",
  "06_modules_WGCNA/04_wgcna_de_gsea_overlap.r",
  "09_export_pride_journal/07_make_biological_claims_table.R",
  "run_dataset_pipeline.R"
)

fail <- character()
for (script in active_scripts) {
  path <- repo_path(script)
  if (!file.exists(path)) {
    fail <- c(fail, paste("Missing active script:", script))
    next
  }
  txt <- readLines(path, warn = FALSE)
  joined <- paste(txt, collapse = "\n")
  if (!grepl("--dry-run|is_dry_run\\(|dry_run", joined)) {
    fail <- c(fail, paste("No dry-run handling detected:", script))
  }
  if (grepl("/Users/|/Volumes/|C:\\\\Users\\\\|~/Desktop|~/Documents", joined)) {
    fail <- c(fail, paste("Hard-coded local path detected:", script))
  }
  if (!grepl("create_module_dirs\\(|module_paths\\(|qc_paths\\(|path_results\\(|path_processed\\(", joined)) {
    fail <- c(fail, paste("No canonical path helper detected:", script))
  }
}

pipeline_txt <- paste(readLines(repo_path("pipeline.yml"), warn = FALSE), collapse = "\n")
for (dataset in c("neuron_neuropil", "neuron_soma", "microglia")) {
  ok <- tryCatch({
    source(repo_path("R", "dataset_config.R"), local = TRUE)
    validate_dataset(dataset) == dataset
  }, error = function(e) FALSE)
  if (!ok) fail <- c(fail, paste("Dataset validation failed:", dataset))
}
for (needle in c("03_qc_exploration/00_dataset_qc_report.r", "06_modules_WGCNA/01_WGCNA.r")) {
  if (!grepl(needle, pipeline_txt, fixed = TRUE)) {
    fail <- c(fail, paste("Pipeline does not include expected script:", needle))
  }
}

enrichment_order <- c(
  "04_differential_expression_enrichment/01_clusterProfiler.r",
  "04_differential_expression_enrichment/02_compareGO.r",
  "04_differential_expression_enrichment/04_neuropil_reference_annotation.r",
  "04_differential_expression_enrichment/05_microglia_targeted_signature_enrichment.r",
  "04_differential_expression_enrichment/03_biological_program_summary.r"
)
order_idx <- vapply(enrichment_order, function(x) {
  m <- regexpr(x, pipeline_txt, fixed = TRUE)
  as.integer(m[[1]])
}, integer(1))
if (any(order_idx <= 0) || is.unsorted(order_idx, strictly = TRUE)) {
  fail <- c(fail, "Pipeline enrichment stage order does not match required clusterProfiler -> compareGO -> neuropil -> microglia signature -> biological program summary sequence")
}

micro_sig_txt <- paste(readLines(repo_path("04_differential_expression_enrichment", "05_microglia_targeted_signature_enrichment.r"), warn = FALSE), collapse = "\n")
for (needle in c(
  "left_unit", "right_unit", "left_region", "right_region", "left_condition", "right_condition", "contrast_class",
  "microglia_signature_enrichment_with_contrast_class.csv",
  "microglia_signature_within_region_condition.csv",
  "microglia_signature_cross_region_same_condition.csv",
  "microglia_signature_cross_region_cross_condition.csv",
  "microglia_signature_leading_edge_recurrence.csv",
  "microglia_signature_claims_ready.csv"
)) {
  if (!grepl(needle, micro_sig_txt, fixed = TRUE)) {
    fail <- c(fail, paste("Microglia signature script missing expected contract marker:", needle))
  }
}

claims_export_txt <- paste(readLines(repo_path("09_export_pride_journal", "07_make_biological_claims_table.R"), warn = FALSE), collapse = "\n")
for (needle in c("microglia_signature_claims_ready.csv", "collect_microglia_signature_claims", "microglia_signature_enrichment")) {
  if (!grepl(needle, claims_export_txt, fixed = TRUE)) {
    fail <- c(fail, paste("Claims export missing expected microglia signature integration marker:", needle))
  }
}

if (length(fail)) {
  message("FAIL smoke_test_active_script_contracts")
  message(paste(fail, collapse = "\n"))
  quit(status = 1, save = "no")
}

message("PASS smoke_test_active_script_contracts")
