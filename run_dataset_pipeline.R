#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
setwd(repo_root())

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}

has_flag <- function(flag) flag %in% args

dataset_arg <- arg_value("--dataset", default = current_dataset())
dataset <- validate_dataset(dataset_arg, source = "--dataset")
dry_run <- has_flag("--dry-run") || tolower(Sys.getenv("PROTEOMICS_DRY_RUN", unset = "")) %in% c("1", "true", "yes")

Sys.setenv(PROTEOMICS_DATASET = dataset)
if (isTRUE(dry_run)) Sys.setenv(PROTEOMICS_DRY_RUN = "true")

steps <- c(
  "01_preprocessing/03_gct_extractR.r",
  "02_id_mapping/01_MapThatProt_batch.r",
  "04_differential_expression_enrichment/01_clusterProfiler.r",
  "04_differential_expression_enrichment/02_compareGO.r"
)

cat("Dataset pipeline\n")
cat("Resolved dataset:", dataset, "\n")
cat("Dry run:", dry_run, "\n")
cat("Project root:", repo_root(), "\n\n")

for (step in steps) {
  script_path <- repo_path(step)
  if (!file.exists(script_path)) stop("Pipeline step not found: ", script_path, call. = FALSE)

  step_args <- if (isTRUE(dry_run)) c(script_path, "--dry-run") else script_path
  cat("==> Running ", step, "\n", sep = "")
  status <- system2("Rscript", step_args)
  if (!identical(status, 0L)) {
    stop("Pipeline step failed with status ", status, ": ", step, call. = FALSE)
  }
}

cat("\nDataset pipeline finished successfully.\n")
