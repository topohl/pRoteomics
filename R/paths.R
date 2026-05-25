# Shared path helpers for pRoteomics.
#
# All committed scripts should default to repository-relative paths. For local
# machines, set PROTEOMICS_PROJECT_ROOT or provide a config override outside Git.

repo_root <- function() {
  env_root <- Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  if (requireNamespace("rprojroot", quietly = TRUE)) {
    root <- tryCatch(
      rprojroot::find_root(rprojroot::has_file("README.md") | rprojroot::is_git_root),
      error = function(e) NULL
    )
    if (!is.null(root)) return(normalizePath(root, winslash = "/", mustWork = FALSE))
  }

  cur <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    markers <- c(".git", "README.md", "01_preprocessing")
    if (any(file.exists(file.path(cur, markers)))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not locate repository root. Set PROTEOMICS_PROJECT_ROOT.", call. = FALSE)
    }
    cur <- parent
  }
}

repo_path <- function(...) file.path(repo_root(), ...)

path_raw <- function(...) repo_path("data", "raw", ...)
path_metadata <- function(...) repo_path("data", "metadata", ...)
path_external <- function(...) repo_path("data", "external", ...)
path_processed <- function(...) repo_path("data", "processed", ...)
path_results <- function(...) repo_path("results", ...)

dir_create <- function(...) {
  path <- file.path(...)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ensure_dir <- dir_create

safe_filename <- function(x, max_chars = 180) {
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>|]+", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- ifelse(nzchar(x), x, "unnamed")
  substr(x, 1, max_chars)
}

file_hash <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

path_or_env <- function(env, default, must_exist = FALSE, kind = c("file", "dir", "any")) {
  kind <- match.arg(kind)
  value <- Sys.getenv(env, unset = "")
  path <- if (nzchar(value)) value else default
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)

  if (isTRUE(must_exist)) {
    exists_ok <- switch(
      kind,
      file = file.exists(path),
      dir = dir.exists(path),
      any = file.exists(path) || dir.exists(path)
    )
    if (!exists_ok) {
      stop(
        "Required ", kind, " path does not exist: ", path,
        ". Set ", env, " to override this location.",
        call. = FALSE
      )
    }
  }

  path
}

relative_to <- function(path, root = repo_root()) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}])", "\\\\\\1", root), "/?"), "", path)
}

write_session_info <- function(path = path_results("logs", "sessionInfo.txt")) {
  dir_create(dirname(path))
  capture.output(utils::sessionInfo(), file = path)
  invisible(path)
}

write_config_snapshot <- function(config, path) {
  dir_create(dirname(path))
  if (requireNamespace("yaml", quietly = TRUE)) {
    writeLines(yaml::as.yaml(config), path)
  } else {
    capture.output(str(config), file = path)
  }
  invisible(path)
}

write_run_manifest <- function(path, inputs = list(), outputs = list(), parameters = list(), notes = NULL) {
  dir_create(dirname(path))

  flatten_paths <- function(x) {
    vals <- unlist(x, recursive = TRUE, use.names = TRUE)
    vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
    as.character(vals)
  }

  input_paths <- flatten_paths(inputs)
  input_hashes <- lapply(input_paths, file_hash)
  names(input_hashes) <- names(input_paths)

  session_path <- file.path(dirname(path), "sessionInfo.txt")
  capture.output(utils::sessionInfo(), file = session_path)

  manifest <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    repo_root = repo_root(),
    inputs = inputs,
    input_hashes = input_hashes,
    outputs = outputs,
    parameters = parameters,
    notes = notes,
    session_info = relative_to(session_path)
  )

  if (requireNamespace("yaml", quietly = TRUE)) {
    writeLines(yaml::as.yaml(manifest), path)
  } else {
    capture.output(str(manifest, max.level = 4), file = path)
  }

  invisible(path)
}

is_dry_run <- function(config = NULL) {
  from_args <- "--dry-run" %in% commandArgs(trailingOnly = FALSE) ||
    "--dry-run" %in% commandArgs(trailingOnly = TRUE)
  from_env <- tolower(Sys.getenv("PROTEOMICS_DRY_RUN", unset = "")) %in% c("1", "true", "yes")
  from_config <- FALSE
  if (is.list(config)) {
    from_config <- isTRUE(config$dry_run) ||
      isTRUE(config$runtime$dry_run)
  }
  isTRUE(from_args || from_env || from_config)
}

dry_run_line <- function(label, value = "", status = NULL) {
  prefix <- if (is.null(status)) "[DRY-RUN]" else paste0("[DRY-RUN ", status, "]")
  message(prefix, " ", label, if (nzchar(as.character(value))) paste0(": ", value) else "")
}

module_paths <- function(module, substep = NULL) {
  tail <- c(module, substep)
  tail <- tail[nzchar(tail)]
  list(
    processed = do.call(path_processed, as.list(tail)),
    figures = do.call(path_results, as.list(c("figures", tail))),
    tables = do.call(path_results, as.list(c("tables", tail))),
    source_data = do.call(path_results, as.list(c("source_data", tail))),
    logs = do.call(path_results, as.list(c("logs", tail))),
    reports = do.call(path_results, as.list(c("reports", tail)))
  )
}

create_module_dirs <- function(module, substep = NULL) {
  paths <- module_paths(module, substep)
  invisible(lapply(paths, dir_create))
  paths
}

pride_submission_dir <- function(...) repo_path("pride_submission", ...)

ensure_pride_dirs <- function() {
  dirs <- pride_submission_dir(c(
    "metadata",
    "processed_data",
    "supplementary_tables",
    "methods",
    "manifests",
    "validation"
  ))
  invisible(lapply(dirs, dir_create))
  invisible(dirs)
}

# Backward-compatible alias used by earlier PRIDE helper scripts.
pride_package_dir <- function() pride_submission_dir()
