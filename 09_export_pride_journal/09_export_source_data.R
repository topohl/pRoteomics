#!/usr/bin/env Rscript

# Collect journal source-data tables without recomputing analyses (not raw-MS/search outputs).

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "validation_utils.R"))

dry_run <- is_dry_run()
target_source <- path_results("manuscript", "source_data")
target_supp <- path_results("manuscript", "supplementary_tables")
dir_create(target_source)
dir_create(target_supp)

table_roots <- c(path_results("tables"), path_results("source_data"))

if (isTRUE(dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/09_export_source_data.R")
  dry_run_line("Candidate table/source-data roots", paste(table_roots, collapse = "; "))
  dry_run_line("Source data output", target_source)
  dry_run_line("Supplementary table output", target_supp)
  quit(status = 0, save = "no")
}

candidates <- unlist(lapply(table_roots[dir.exists(table_roots)], list.files, pattern = "\\.(csv|tsv|xlsx)$", recursive = TRUE, full.names = TRUE), use.names = FALSE)
candidates <- candidates[!grepl("results/manuscript", normalizePath(candidates, winslash = "/", mustWork = FALSE), fixed = TRUE)]

manifest <- data.frame(
  source_file = candidates,
  target_file = file.path(ifelse(grepl("source_data", candidates, fixed = TRUE), target_source, target_supp), basename(candidates)),
  stringsAsFactors = FALSE
)
if (nrow(manifest)) {
  file.copy(manifest$source_file, manifest$target_file, overwrite = TRUE)
}

manifest_path <- path_results("manuscript", "source_data_export_manifest.csv")
utils::write.csv(manifest, manifest_path, row.names = FALSE)
write_run_manifest(
  path_results("logs", "09_export_pride_journal", "source_data", "run_manifest.yml"),
  inputs = list(tables = candidates),
  outputs = list(manifest = manifest_path, source_data = target_source, supplementary_tables = target_supp),
  notes = "Collect-only manuscript source-data export; no analyses are recomputed."
)
message("Source-data export manifest written: ", manifest_path)
