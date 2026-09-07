#!/usr/bin/env Rscript

# Write docs/publication_freeze_manifest.yml from the current repository state.
#
# Read-only with respect to every scientific output: it hashes and reads
# existing artifacts and never recomputes, regenerates or rewrites them.
#
# Usage:
#   Rscript tools/generate_publication_freeze_manifest.R
#   Rscript tools/generate_publication_freeze_manifest.R --test-log <path>
#
# --test-log points at a persisted full-suite log to record as the freeze
# validation state. Without it the validation_state block records that no
# authoritative persisted result was available.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "publication_freeze_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1L]]
}

test_log <- value_after("--test-log")
test_state <- if (!is.na(test_log) && file.exists(test_log)) {
  log_lines <- readLines(test_log, warn = FALSE)
  # testthat's summary reporter marks a failure or error with a NUMBERED digit
  # in the per-file symbol run, not with the letters F or E. Counting letters
  # reports zero on a failing log, so the digits are what must be counted.
  file_lines <- grep("^[a-z0-9][a-z0-9-]*: ", log_lines, value = TRUE)
  symbols <- sub("^[a-z0-9][a-z0-9-]*: ", "", file_lines)
  failure_markers <- sum(vapply(gregexpr("[0-9]", symbols),
                                function(m) sum(m > 0), integer(1)))
  failing_files <- sub(":.*$", "", file_lines[grepl("[0-9]", symbols)])
  # Only record a path if the log lives inside the repository; a transient
  # out-of-tree log would make the manifest machine-specific.
  rel_log <- freeze_rel(test_log)
  in_repo <- !identical(rel_log, gsub("\\\\", "/", test_log))
  list(
    command = "Rscript tests/testthat.R",
    log_path = if (in_repo) rel_log else NA_character_,
    log_sha256 = file_hash_sha256(test_log),
    log_retained_in_repository = in_repo,
    test_files = length(file_lines),
    failures = sum(grepl("FAILURE", log_lines)),
    failure_or_error_markers = failure_markers,
    failing_files = if (length(failing_files)) failing_files else NULL,
    clean = failure_markers == 0L &&
      !any(grepl("FAILURE", log_lines)) &&
      !any(grepl("^(Failed|Error)", log_lines)),
    skipped = sum(grepl("Reason:", log_lines)),
    warnings = sum(grepl("problems\\(dat\\)", log_lines)),
    known_skip = "test-wgcna-group-effects-phase2b.R:248 optional lme4 zero-variance fixture",
    known_warnings = "two vroom parsing warnings in test-gsea-wgcna-ontology-theme-integration.R",
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    note = paste("counts parsed from a testthat summary log produced at the freeze",
                 "commit; the log itself was transient and is identified by sha256")
  )
} else {
  list(
    command = "Rscript tests/testthat.R",
    log_path = NA_character_,
    note = paste("no authoritative persisted full-suite log was supplied;",
                 "record one with --test-log to populate this block")
  )
}

generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

manifest <- build_publication_freeze_manifest(
  generated_at = generated_at,
  test_state = test_state,
  strict = TRUE
)

out <- write_publication_freeze_manifest(manifest)
cat("Publication freeze manifest written:", freeze_rel(out), "\n")
cat("  freeze commit :", manifest$freeze_identity$freeze_git_commit, "\n")
cat("  freeze tag    :", manifest$freeze_identity$freeze_git_tag, "\n")
cat("  sha256        :", file_hash_sha256(out), "\n")
