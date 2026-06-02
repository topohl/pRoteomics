#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

message("Deprecated entrypoint: 06_modules_WGCNA/02_wgcna_de_gsea_overlap.r. Use 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r instead.")
source(repo_path("06_modules_WGCNA", "04_wgcna_de_gsea_overlap.r"))

dataset <- arg_value("--dataset", default = current_dataset())
if (nzchar(dataset)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset, source = "--dataset"))
run_wgcna_de_gsea_overlap(current_dataset(), dry_run = "--dry-run" %in% args || is_dry_run())
