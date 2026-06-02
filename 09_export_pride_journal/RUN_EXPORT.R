#!/usr/bin/env Rscript
# Orchestrates pg_matrix-onward PRIDE/journal export.
#
# Usage:
#   Rscript 09_export_pride_journal/RUN_EXPORT.R --dataset all --export-level pg_matrix_onward
#   Rscript 09_export_pride_journal/RUN_EXPORT.R --dataset microglia --export-level pg_matrix_onward
#   Rscript 09_export_pride_journal/RUN_EXPORT.R --dataset microglia --dry-run

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "export_helpers.R"))

setwd(repo_root())
cli <- export_cli_args()
extra <- commandArgs(trailingOnly = TRUE)
extra <- extra[!extra %in% c("--skip-validation", "--skip-supplementary", "--skip-manuscript", "--skip-claims")]

run_script <- function(rel, required = TRUE) {
  script <- repo_path(rel)
  if (!file.exists(script)) {
    msg <- paste("Export script missing:", rel)
    if (isTRUE(required)) stop(msg, call. = FALSE) else { warning(msg, call. = FALSE); return(0L) }
  }
  args <- c(script, extra)
  cat("\n==> ", rel, "\n", sep = "")
  status <- system2("Rscript", args)
  if (is.null(status)) status <- 0L
  if (!identical(status, 0L) && isTRUE(required) && !isTRUE(cli$dry_run)) {
    stop("Export step failed (exit ", status, "): ", rel, call. = FALSE)
  }
  as.integer(status)
}

cat("pRoteomics export runner\n")
cat("Scope:", export_scope_label(cli$export_level), "\n")
cat("Dataset:", cli$dataset, "\n")
cat("Dry run:", cli$dry_run, "\n\n")

steps <- c(
  "09_export_pride_journal/02_make_sample_metadata.R",
  "09_export_pride_journal/03_export_processed_pg_matrix_package.R"
)
if (!cli$skip_supplementary) steps <- c(steps, "09_export_pride_journal/04_make_supplementary_tables.R")
if (!cli$skip_manuscript) {
  steps <- c(steps, "09_export_pride_journal/09_export_source_data.R", "09_export_pride_journal/08_export_manuscript_figures.R")
}
if (!cli$skip_claims) steps <- c(steps, "09_export_pride_journal/07_make_biological_claims_table.R")
steps <- c(steps, "09_export_pride_journal/06_make_methods_summary.R", "09_export_pride_journal/01_make_pride_manifest.R")
if (!cli$skip_validation) steps <- c(steps, "09_export_pride_journal/05_validate_pride_submission.R")

statuses <- vapply(steps, run_script, integer(1))
cat("\nExport runner finished. Steps:", length(steps), " Failures:", sum(statuses != 0L), "\n")
quit(status = if (any(statuses != 0L)) 1 else 0, save = "no")
