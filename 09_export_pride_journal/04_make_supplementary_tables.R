#!/usr/bin/env Rscript
# Journal/source-data supplementary tables (not PRIDE-required raw/search outputs).
# Copies only config-listed globs by default; no recursive results/ scan.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "export_helpers.R"))

ensure_pride_dirs()
config <- load_export_config()
cli <- export_cli_args()
datasets <- resolve_export_datasets(cli$dataset, config)
out_dir <- pride_submission_dir("supplementary_tables")

globs <- config$supplementary_table_globs
if (length(datasets) < length(config$datasets)) {
  globs <- globs[vapply(globs, function(g) any(vapply(datasets, function(ds) grepl(ds, g, fixed = TRUE), logical(1))) || !grepl("<dataset>", g, fixed = TRUE), logical(1))]
}
files <- supplementary_candidate_files(config, datasets, include_derived = TRUE)

if (isTRUE(cli$dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/04_make_supplementary_tables.R")
  dry_run_line("Datasets", paste(datasets, collapse = ", "))
  dry_run_line("Supplementary candidate count", length(files), if (length(files) > 0) "PASS" else "WARN")
  dry_run_line("Supplementary target", out_dir)
  quit(status = 0, save = "no")
}

if (length(files) == 0) {
  warning(
    "No supplementary candidate tables matched config globs. ",
    "Run with --include-derived-results after generating results/tables outputs.",
    call. = FALSE
  )
} else {
  for (f in files) {
    rel <- safe_filename(relative_to(f, repo_root()))
    copy_export_file(f, file.path(out_dir, rel), dry_run = FALSE)
  }
  message("Staged ", length(files), " journal supplementary table(s) in: ", out_dir)
}
