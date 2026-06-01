#!/usr/bin/env Rscript

# Collect final manuscript figures from module outputs without recomputing analyses.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
dry_run <- is_dry_run()

manuscript_dirs <- path_results("manuscript", c("figure_1", "figure_2", "extended_data"))
invisible(lapply(manuscript_dirs, dir_create))

candidate_roots <- c(
  path_results("figures", "03_qc_exploration"),
  path_results("figures", "04_differential_expression_enrichment"),
  path_results("figures", "06_modules_WGCNA"),
  path_results("figures", "07_spatial_networks"),
  path_results("figures", "08_behavior_physio_coupling")
)

if (isTRUE(dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/07_export_manuscript_figures.R")
  dry_run_line("Candidate figure roots", paste(candidate_roots, collapse = "; "))
  dry_run_line("Output root", path_results("manuscript"))
  quit(status = 0, save = "no")
}

candidates <- unlist(lapply(candidate_roots[dir.exists(candidate_roots)], list.files, pattern = "\\.(svg|pdf|png)$", recursive = TRUE, full.names = TRUE), use.names = FALSE)

manifest <- data.frame(
  source_file = candidates,
  target_file = file.path(path_results("manuscript", "extended_data"), basename(candidates)),
  stringsAsFactors = FALSE
)
if (nrow(manifest)) {
  file.copy(manifest$source_file, manifest$target_file, overwrite = TRUE)
}

manifest_path <- path_results("manuscript", "figure_export_manifest.csv")
utils::write.csv(manifest, manifest_path, row.names = FALSE)
write_run_manifest(
  path_results("logs", "09_export_pride_journal", "manuscript_figures", "run_manifest.yml"),
  inputs = list(figures = candidates),
  outputs = list(manifest = manifest_path, manuscript_root = path_results("manuscript")),
  notes = "Collect-only manuscript figure export; no analyses are recomputed."
)
message("Manuscript figure export manifest written: ", manifest_path)
