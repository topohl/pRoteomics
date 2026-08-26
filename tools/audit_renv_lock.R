#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tools/audit_renv_lock.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "R", "renv_lock_audit.R"))
audit <- audit_renv_lock(file.path(repo_root, "renv.lock"))

cat("renv.lock package records:", audit$package_count, "\n")
cat("Recorded packages:", paste(audit$packages, collapse = ", "), "\n")
if (!audit$plausibly_full_scientific_lock) {
  cat(
    "Missing scientific sentinel packages:",
    paste(audit$missing_scientific_sentinels, collapse = ", "),
    "\n"
  )
  cat("STATUS: INCOMPLETE_FOR_FULL_SCIENTIFIC_REPRODUCTION\n")
  quit(save = "no", status = 1L)
}

cat("STATUS: SCIENTIFIC_SENTINELS_PRESENT\n")
