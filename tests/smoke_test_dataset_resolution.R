#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "validation_utils.R"))

fail <- character()
expected <- c("neuron_neuropil", "neuron_soma", "microglia")
if (!identical(valid_datasets(), expected)) {
  fail <- c(fail, paste("valid_datasets changed unexpectedly:", paste(valid_datasets(), collapse = ", ")))
}

fixture_root <- tempfile("dataset_resolution_contract_")
dir.create(fixture_root, recursive = TRUE, showWarnings = FALSE)

for (dataset in expected) {
  fixture_paths <- c(
    wgcna_expression = file.path(
      fixture_root,
      paste0("20260826_pgmatrix_imputed_", dataset, "_fixture_missing70pct.xlsx")
    ),
    wgcna_metadata = file.path(fixture_root, paste0(dataset, "_wgcna_metadata.xlsx")),
    module_expression = file.path(
      fixture_root,
      paste0("20260826_pgmatrix_imputed_", dataset, "_fixture_missing70pct_with_metadata.xlsx")
    ),
    module_metadata = file.path(fixture_root, paste0(dataset, "_module_metadata.xlsx"))
  )
  if (!all(file.create(fixture_paths))) {
    stop("Could not create dataset-resolution smoke fixtures.", call. = FALSE)
  }
  do.call(Sys.setenv, as.list(c(
    PROTEOMICS_WGCNA_UPSTREAM_XLSX = fixture_paths[["wgcna_expression"]],
    PROTEOMICS_WGCNA_UPSTREAM_META_XLSX = fixture_paths[["wgcna_metadata"]],
    PROTEOMICS_MODULE_SCORE_PROTEIN_FILE = fixture_paths[["module_expression"]],
    PROTEOMICS_MODULE_SCORE_METADATA_FILE = fixture_paths[["module_metadata"]]
  )))

  w <- resolve_dataset_inputs(dataset, "wgcna", record_resolution = FALSE)
  m <- resolve_dataset_inputs(dataset, "module_score", record_resolution = FALSE)
  if (!identical(w$dataset, dataset)) fail <- c(fail, paste("WGCNA resolver returned wrong dataset for", dataset))
  if (!identical(m$dataset, dataset)) fail <- c(fail, paste("module_score resolver returned wrong dataset for", dataset))
  if (!length(w$diagnostics) || !length(m$diagnostics)) fail <- c(fail, paste("Missing diagnostics for", dataset))
  if (!identical(w$expression_file, normalizePath(fixture_paths[["wgcna_expression"]], winslash = "/"))) {
    fail <- c(fail, paste("WGCNA explicit input override not honored for", dataset))
  }
  if (!identical(m$expression_file, normalizePath(fixture_paths[["module_expression"]], winslash = "/"))) {
    fail <- c(fail, paste("module_score explicit input override not honored for", dataset))
  }
}

Sys.unsetenv(c(
  "PROTEOMICS_WGCNA_UPSTREAM_XLSX",
  "PROTEOMICS_WGCNA_UPSTREAM_META_XLSX",
  "PROTEOMICS_MODULE_SCORE_PROTEIN_FILE",
  "PROTEOMICS_MODULE_SCORE_METADATA_FILE"
))

metadata_fixture <- data.frame(
  CellTypeLayer = c("Neuron_Neuropil", "Neuron_Soma", "Microglia"),
  stringsAsFactors = FALSE
)
for (dataset in expected) {
  keep <- metadata_matches_dataset(metadata_fixture, dataset)
  if (!identical(which(keep), match(dataset, expected))) {
    fail <- c(fail, paste("Dataset metadata filter leakage for", dataset))
  }
}

if (length(fail)) {
  message("FAIL smoke_test_dataset_resolution")
  message(paste(fail, collapse = "\n"))
  quit(status = 1, save = "no")
}

message("PASS smoke_test_dataset_resolution")
