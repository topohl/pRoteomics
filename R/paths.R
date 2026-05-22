# Shared path helpers for pRoteomics
#
# These helpers avoid hard-coded machine-specific paths and are used by the
# PRIDE/ProteomeXchange packaging scripts. They intentionally depend only on
# base R so they can run in a minimal R installation.

repo_root <- function() {
  env_root <- Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  cur <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    markers <- c(".git", "README.md", "01_preprocessing")
    if (any(file.exists(file.path(cur, markers)))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not locate repository root. Set PROTEOMICS_PROJECT_ROOT.", call. = FALSE)
    }
    cur <- parent
  }
}

repo_path <- function(...) {
  file.path(repo_root(), ...)
}

pride_package_dir <- function() {
  file.path(repo_root(), "PRIDE_package")
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

ensure_pride_dirs <- function() {
  dirs <- file.path(
    pride_package_dir(),
    c(
      "00_metadata",
      "01_raw",
      "02_search_results",
      "02_search_results/search_parameters",
      "02_search_results/fasta_database",
      "03_processed_quantification",
      "04_downstream_analysis",
      "05_scripts"
    )
  )
  invisible(lapply(dirs, ensure_dir))
  invisible(dirs)
}
