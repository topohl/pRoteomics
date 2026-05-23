# Creates a lightweight methods/provenance summary for PRIDE and journal handoff.
# Consumes manifests/configs/session info; produces pride_submission/methods/methods_summary.md.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

ensure_pride_dirs()
DRY_RUN <- is_dry_run()
out_file <- pride_submission_dir("methods", "methods_summary.md")
if (isTRUE(DRY_RUN)) {
  dry_run_line("Script", "09_export_pride_journal/05_make_methods_summary.R")
  dry_run_line("Methods summary target", out_file)
  quit(status = 0, save = "no")
}
lines <- c(
  "# Proteomics Methods and Provenance Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "",
  "## Analysis Layers",
  "",
  "1. Raw/vendor MS files and search outputs are staged locally under pride_submission/.",
  "2. Canonical processed quantitative matrices live under data/processed/01_preprocessing/.",
  "3. Mapped protein/contrast tables live under data/processed/02_id_mapping/.",
  "4. Differential-expression and enrichment outputs live under data/processed and results for module 04.",
  "5. Figure source data live under results/source_data/.",
  "",
  "## Reproducibility Metadata",
  "",
  "- Input hashes are recorded where scripts use manifest-backed contracts.",
  "- Config snapshots and sessionInfo.txt are written by major refactored scripts.",
  "- Database versions and statistical thresholds should be completed in the data dictionary before submission."
)
writeLines(lines, out_file)
message("Wrote methods summary: ", out_file)
