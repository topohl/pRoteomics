renv_lock_package_names <- function(lockfile) {
  if (!file.exists(lockfile)) {
    stop("renv lockfile does not exist: ", lockfile, call. = FALSE)
  }

  lines <- readLines(lockfile, warn = FALSE)
  packages_start <- grep('^  "Packages": \\{$', lines)
  if (length(packages_start) != 1L) {
    stop("renv lockfile must contain exactly one Packages object", call. = FALSE)
  }

  package_lines <- lines[seq.int(packages_start + 1L, length(lines))]
  package_keys <- regexec('^    "([^"]+)": \\{$', package_lines)
  matches <- regmatches(package_lines, package_keys)
  sort(unique(vapply(matches[lengths(matches) == 2L], `[[`, character(1), 2L)))
}

renv_lock_scientific_sentinels <- function() {
  c(
    "clusterProfiler", "dplyr", "fgsea", "ggplot2", "limma", "lme4",
    "readr", "WGCNA"
  )
}

audit_renv_lock <- function(
    lockfile,
    scientific_sentinels = renv_lock_scientific_sentinels()) {
  packages <- renv_lock_package_names(lockfile)
  missing <- setdiff(scientific_sentinels, packages)

  list(
    lockfile = normalizePath(lockfile, winslash = "/", mustWork = TRUE),
    package_count = length(packages),
    packages = packages,
    scientific_sentinels = scientific_sentinels,
    missing_scientific_sentinels = missing,
    plausibly_full_scientific_lock = length(missing) == 0L
  )
}
