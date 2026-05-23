# Copies selected canonical tables into pride_submission/supplementary_tables.
# Consumes results/tables and results/source_data; produces journal supplementary table staging.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

ensure_pride_dirs()
DRY_RUN <- is_dry_run()
roots <- c(path_results("tables"), path_results("source_data"))
roots <- roots[dir.exists(roots)]
files <- unique(unlist(lapply(roots, list.files, pattern = "\\.(csv|tsv|xlsx)$", recursive = TRUE, full.names = TRUE)))
if (isTRUE(DRY_RUN)) {
  dry_run_line("Script", "09_export_pride_journal/03_make_supplementary_tables.R")
  dry_run_line("Candidate roots", paste(roots, collapse = "; "))
  dry_run_line("Supplementary candidate count", length(files), if (length(files) > 0) "PASS" else "WARN")
  dry_run_line("Supplementary target", pride_submission_dir("supplementary_tables"))
  quit(status = 0, save = "no")
}
if (length(files) == 0) {
  warning("No supplementary candidate tables found under results/tables or results/source_data.")
} else {
  out_dir <- pride_submission_dir("supplementary_tables")
  for (f in files) {
    rel <- safe_filename(relative_to(f, repo_root()))
    file.copy(f, file.path(out_dir, rel), overwrite = TRUE)
  }
  message("Staged ", length(files), " supplementary table candidates in: ", out_dir)
}
