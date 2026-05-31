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

for (dataset in expected) {
  w <- resolve_dataset_inputs(dataset, "wgcna")
  m <- resolve_dataset_inputs(dataset, "module_score")
  if (!identical(w$dataset, dataset)) fail <- c(fail, paste("WGCNA resolver returned wrong dataset for", dataset))
  if (!identical(m$dataset, dataset)) fail <- c(fail, paste("module_score resolver returned wrong dataset for", dataset))
  if (!length(w$diagnostics) || !length(m$diagnostics)) fail <- c(fail, paste("Missing diagnostics for", dataset))

  if (file.exists(w$expression_file) && file.exists(w$metadata_file) && requireNamespace("readxl", quietly = TRUE)) {
    expr_head <- readxl::read_excel(w$expression_file, n_max = 3)
    md_head <- readxl::read_excel(w$metadata_file, n_max = 1000)
    sample_col <- detect_column(md_head, w$sample_id_col_candidates)
    if (is.na(sample_col)) {
      fail <- c(fail, paste("No metadata sample column detected for", dataset))
    } else {
      overlap <- sample_overlap_summary(setdiff(names(expr_head), w$protein_id_col_candidates), md_head[[sample_col]])
      if (overlap$n_overlap == 0) fail <- c(fail, paste("No sample overlap for", dataset))
      keep <- metadata_matches_dataset(md_head, dataset)
      if (any(keep) && "CellTypeLayer" %in% names(md_head)) {
        leaked <- unique(tolower(as.character(md_head$CellTypeLayer[keep])))
        if (dataset == "neuron_neuropil" && any(grepl("soma|microglia", leaked))) fail <- c(fail, "neuron_neuropil filter leakage")
        if (dataset == "neuron_soma" && any(grepl("neuropil|microglia", leaked))) fail <- c(fail, "neuron_soma filter leakage")
        if (dataset == "microglia" && any(grepl("neuron|soma|neuropil", leaked))) fail <- c(fail, "microglia filter leakage")
      }
    }
  } else {
    message("SKIP sample overlap for ", dataset, ": local matrix/metadata not both available")
  }
}

if (length(fail)) {
  message("FAIL smoke_test_dataset_resolution")
  message(paste(fail, collapse = "\n"))
  quit(status = 1, save = "no")
}

message("PASS smoke_test_dataset_resolution")
