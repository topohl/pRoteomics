#!/usr/bin/env Rscript

# ================================================================
# Script: 09_export_pride_journal/09_export_source_data.R
# Stage: export
# Scope: global
# Consumes: required results/tables/; optional results/source_data/.
# Produces: results/manuscript/source_data/; results/manuscript/supplementary_tables/; results/manuscript/source_data_export_manifest.csv.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Manuscript source-data export.
# ================================================================

# Collect journal source-data tables without recomputing analyses (not raw-MS/search outputs).

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "export_helpers.R"))

dry_run <- is_dry_run()
target_source <- path_results("manuscript", "source_data")
target_supp <- path_results("manuscript", "supplementary_tables")

table_roots <- c(path_results("tables"), path_results("source_data"))

# Selection is read-only, so it is computed before the dry-run guard and
# reported by it. Nothing below this point may mutate the filesystem until
# after the guard: --dry-run must be side-effect free.
candidates <- unlist(lapply(table_roots[dir.exists(table_roots)], list.files, pattern = "\\.(csv|tsv|xlsx)$", recursive = TRUE, full.names = TRUE), use.names = FALSE)
candidates <- candidates[is_exportable_result_path(candidates)]

# Journal source-data scope: drop diagnostic, superseded, proposed and
# intermediate families that are not manuscript source data. Applies to this
# export only -- analysis outputs are untouched and PRIDE packaging is
# unaffected. See the scope block in R/export_helpers.R for per-family
# justification.
scope_before <- length(candidates)
scope_reasons <- source_data_scope_exclusion_reasons(candidates)
candidates <- candidates[is.na(scope_reasons)]
scope_dropped <- table(scope_reasons[!is.na(scope_reasons)])

# Targets are budgeted and routed read-only, so the dry-run can report the same
# numbers the real run will enforce. Nothing here mutates the filesystem.
target_paths <- manuscript_table_target_paths(candidates, target_source, target_supp)
selected_info <- file.info(candidates)
n_inaccessible <- sum(!file.exists(candidates) | is.na(selected_info$size))
n_over_budget <- sum(nchar(target_paths) > manuscript_figure_path_budget())
n_duplicate_targets <- sum(duplicated(target_paths))

if (isTRUE(dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/09_export_source_data.R")
  dry_run_line("Candidate table/source-data roots", paste(table_roots, collapse = "; "))
  dry_run_line("Source data output", target_source)
  dry_run_line("Supplementary table output", target_supp)
  dry_run_line("Raw candidates before journal scope", scope_before)
  for (reason in names(scope_dropped)) {
    dry_run_line(paste0("Excluded (", reason, ")"), scope_dropped[[reason]])
  }
  dry_run_line("Scoped selected files", length(candidates))
  dry_run_line("Selected bytes", format(sum(selected_info$size, na.rm = TRUE), big.mark = ","))
  dry_run_line("Selected GB", round(sum(selected_info$size, na.rm = TRUE) / 1e9, 3))
  dry_run_line("Inaccessible selected sources", n_inaccessible)
  dry_run_line(paste0("Target paths over ", manuscript_figure_path_budget(), " chars"), n_over_budget)
  dry_run_line("Duplicate target paths", n_duplicate_targets)
  dry_run_line("Longest target path", if (length(target_paths)) max(nchar(target_paths)) else 0L)
  quit(status = 0, save = "no")
}

manifest <- data.frame(
  source_file = candidates,
  target_file = target_paths,
  stringsAsFactors = FALSE
)

# Fail closed before touching the payload: every source must be readable, every
# target must fit the budget, and no target may be claimed twice. Only then are
# the directories created and the copy verified file by file, so neither the
# manifest nor the run manifest can describe a partial export.
assert_export_sources_accessible(manifest$source_file)
assert_unique_table_export_targets(manifest$target_file, manifest$source_file)
over_budget <- manifest$target_file[nchar(manifest$target_file) > manuscript_figure_path_budget()]
if (length(over_budget)) {
  stop(
    length(over_budget), " target path(s) exceed the ",
    manuscript_figure_path_budget(), "-character budget; first: ",
    basename(over_budget[[1]]), call. = FALSE
  )
}

dir_create(target_source)
dir_create(target_supp)
copy_export_targets(manifest$source_file, manifest$target_file)

manifest_path <- path_results("manuscript", "source_data_export_manifest.csv")
utils::write.csv(manifest, manifest_path, row.names = FALSE)
write_run_manifest(
  path_results("logs", "09_export_pride_journal", "source_data", "run_manifest.yml"),
  inputs = list(tables = candidates),
  outputs = list(manifest = manifest_path, source_data = target_source, supplementary_tables = target_supp),
  notes = "Collect-only manuscript source-data export; no analyses are recomputed."
)
message("Source-data export manifest written: ", manifest_path)
