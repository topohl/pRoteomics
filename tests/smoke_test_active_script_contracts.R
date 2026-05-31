#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

active_scripts <- c(
  "03_qc_exploration/00_dataset_qc_report.r",
  "04_differential_expression_enrichment/03_biological_program_summary.r",
  "06_modules_WGCNA/01_WGCNA.r",
  "06_modules_WGCNA/05_wgcna_de_gsea_overlap.r",
  "09_export_pride_journal/06_make_biological_claims_table.R",
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

pipeline_txt <- paste(readLines(repo_path("run_dataset_pipeline.R"), warn = FALSE), collapse = "\n")
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

if (length(fail)) {
  message("FAIL smoke_test_active_script_contracts")
  message(paste(fail, collapse = "\n"))
  quit(status = 1, save = "no")
}

message("PASS smoke_test_active_script_contracts")
