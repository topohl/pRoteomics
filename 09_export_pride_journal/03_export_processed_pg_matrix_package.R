#!/usr/bin/env Rscript
# ================================================================
# Script: 09_export_pride_journal/03_export_processed_pg_matrix_package.R
# Stage: export
# Scope: global
# Consumes: required data/raw/pg_matrix/; data/processed/01_preprocessing/; +1 more; optional none.
# Produces: pride_submission/processed_data/; pride_submission/processed_data/data_dictionary_processed_matrices.tsv.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Exports processed pg/mapped matrix package.
# ================================================================

# Stages processed quantitative data from pg_matrix onward into pride_submission/processed_data/.
# Does not claim raw-MS/search reproducibility.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "export_helpers.R"))

ensure_pride_dirs()
config <- load_export_config()
cli <- export_cli_args()
datasets <- resolve_export_datasets(cli$dataset, config)
out_root <- pride_submission_dir("processed_data")

collect_dictionary_rows <- function() {
  data.frame(
    filename = character(),
    biological_meaning = character(),
    statistical_level = character(),
    upstream_dependency = character(),
    downstream_usage = character(),
    intended_for_PRIDE = character(),
    transient_or_canonical = character(),
    stringsAsFactors = FALSE
  )
}

dictionary <- collect_dictionary_rows()
staged <- character()
planned_count <- 0L

for (ds in datasets) {
  ds_out <- file.path(out_root, ds)
  dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

  pg_paths <- pg_matrix_input_paths(config)
  for (pg in pg_paths) {
    dest <- file.path(ds_out, "pg_matrix_reference", basename(pg))
    if (isTRUE(cli$dry_run)) {
      dry_run_line("pg_matrix reference", paste(pg, "->", dest))
    } else if (file.exists(pg)) {
      copy_export_file(pg, dest, dry_run = FALSE)
      staged <- c(staged, dest)
      dictionary <- rbind(dictionary, data.frame(
        filename = basename(pg),
        biological_meaning = "Original protein-group quantification matrix (analysis entry point).",
        statistical_level = "raw_pg_matrix",
        upstream_dependency = "external MS search (not in repo)",
        downstream_usage = "01_preprocessing imputation and dataset splits",
        intended_for_PRIDE = "yes",
        transient_or_canonical = "canonical_input",
        stringsAsFactors = FALSE
      ))
    } else {
      warning("pg_matrix input missing (expected for full export): ", pg, call. = FALSE)
    }
  }

  proc_files <- processed_files_for_dataset(ds, config, include_derived = cli$include_derived_results)
  planned_count <- planned_count + length(proc_files) + length(pg_paths)
  for (src in proc_files) {
    rel <- relative_to(src, repo_root())
    dest <- file.path(ds_out, gsub("^data/processed/", "", rel))
    dest <- gsub("^data/raw/", "pg_matrix_reference/", dest)
    if (isTRUE(cli$dry_run)) {
      dry_run_line("processed file", paste(rel, "->", dest))
    } else {
      copy_export_file(src, dest, dry_run = FALSE)
      staged <- c(staged, dest)
    }
    dictionary <- rbind(dictionary, data.frame(
      filename = basename(src),
      biological_meaning = paste0("Processed artifact for dataset ", ds, " (", classify_export_category(src, config), ")."),
      statistical_level = "processed_matrix_or_table",
      upstream_dependency = "pg_matrix onward pipeline",
      downstream_usage = "downstream mapping, DE, enrichment, modules",
      intended_for_PRIDE = "yes",
      transient_or_canonical = "canonical",
      stringsAsFactors = FALSE
    ))
  }
}

dict_path <- pride_submission_dir("processed_data", "data_dictionary_processed_matrices.tsv")
if (isTRUE(cli$dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/03_export_processed_pg_matrix_package.R")
  dry_run_line("Datasets", paste(datasets, collapse = ", "))
  dry_run_line("Staged file count (planned)", planned_count)
  dry_run_line("Dictionary target", dict_path)
  quit(status = 0, save = "no")
}

if (nrow(dictionary)) {
  write_tsv(dictionary, dict_path)
}
message("Processed pg_matrix-onward package staged under: ", out_root)
message("Staged files: ", length(staged))
message("Data dictionary: ", dict_path)
