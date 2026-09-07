#!/usr/bin/env Rscript

# Write renv.lock from the installed analysis library.
#
# This records the environment that exists. It never installs, updates or loads
# a package, and it never creates renv project infrastructure or a project
# library. Versions and source metadata come from installed DESCRIPTION files
# only.
#
# Usage:
#   Rscript tools/generate_renv_lockfile.R            # write renv.lock
#   Rscript tools/generate_renv_lockfile.R --dry-run  # report, write nothing

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tools/generate_renv_lockfile.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

source(file.path(repo_root, "R", "renv_lock_audit.R"))

cli <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% cli
lockfile <- file.path(repo_root, "renv.lock")

installed <- utils::installed.packages()
before <- if (file.exists(lockfile)) renv_lock_package_names(lockfile) else character()

lock <- build_renv_lockfile(root = repo_root, installed = installed)
direct <- attr(lock, "direct_dependencies")
recorded <- names(lock$Packages)

cat("renv lockfile generation\n")
cat("  R version            :", lock$R$Version, "\n")
cat("  Bioconductor version :", if (is.null(lock$Bioconductor)) "unavailable" else lock$Bioconductor$Version, "\n")
cat("  repositories declared:", length(lock$R$Repositories), "\n")
cat("  active files scanned :", length(renv_lock_active_source_files(repo_root)), "\n")
cat("  direct dependencies  :", length(direct), "\n")
cat("  closure (recorded)   :", length(recorded), "\n")
cat("  previous records     :", length(before), "\n")
cat("  added                :", length(setdiff(recorded, before)), "\n")
cat("  removed              :", length(setdiff(before, recorded)), "\n")
removed <- setdiff(before, recorded)
if (length(removed)) {
  cat("  removed names        :", paste(removed, collapse = ", "), "\n")
}

if (isTRUE(dry_run)) {
  cat("\n[DRY-RUN] renv.lock not written\n")
  quit(status = 0, save = "no")
}

write_renv_lockfile(lock, lockfile)
cat("\nrenv.lock written:", lockfile, "\n")

audit <- audit_renv_lock_completeness(lockfile, root = repo_root, installed = installed)
cat("  recorded_count       :", audit$recorded_count, "\n")
cat("  missing direct       :", length(audit$missing_direct), "\n")
cat("  unresolved reqs      :", length(audit$unresolved_requirements), "\n")
cat("  duplicate records    :", length(audit$duplicate_records), "\n")
cat("  complete             :", audit$complete, "\n")
if (!isTRUE(audit$complete)) {
  quit(status = 1, save = "no")
}
