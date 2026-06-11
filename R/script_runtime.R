# Shared runtime helpers for numbered pipeline entrypoint scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("current_dataset_from_cli", mode = "function")) {
  source(repo_path("R", "dataset_config.R"))
}
if (!exists("write_run_manifest", mode = "function")) {
  source(repo_path("R", "validation_utils.R"))
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

script_has_flag <- function(flag, args = commandArgs(trailingOnly = TRUE)) {
  flag %in% args
}

script_arg_value <- function(flag, default = "", args = commandArgs(trailingOnly = TRUE)) {
  idx <- which(args == flag)
  if (!length(idx) || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

init_script_runtime <- function(script,
                                stage,
                                default_dataset = "neuron_neuropil",
                                allow_all = FALSE,
                                args = commandArgs(trailingOnly = TRUE)) {
  dataset_raw <- script_arg_value("--dataset", Sys.getenv("PROTEOMICS_DATASET", unset = default_dataset), args = args)
  if (isTRUE(allow_all) && identical(tolower(dataset_raw), "all")) {
    dataset <- "all"
    Sys.setenv(PROTEOMICS_DATASET = "all")
  } else {
    dataset <- current_dataset_from_cli(default = default_dataset, args = args)
  }
  list(
    script = script,
    stage = stage,
    dataset = dataset,
    args = args,
    dry_run = script_has_flag("--dry-run", args = args),
    started_at = Sys.time()
  )
}

input_status_row <- function(input_name,
                             path,
                             dataset = "global",
                             required = FALSE,
                             status = NULL,
                             message = NULL,
                             n_rows = NA_integer_) {
  exists <- file.exists(path) || dir.exists(path) || length(Sys.glob(path)) > 0L
  status <- status %||% if (exists) "present" else if (isTRUE(required)) "missing_required" else "missing_optional"
  message <- message %||% if (exists) "input available" else "input not available"
  data.frame(
    dataset = dataset,
    input_name = input_name,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    required = isTRUE(required),
    status = status,
    message = message,
    n_rows = n_rows,
    stringsAsFactors = FALSE
  )
}

write_input_status <- function(rows, path, dry_run = FALSE) {
  if (is.null(rows) || !nrow(rows)) return(invisible(path))
  if (isTRUE(dry_run)) {
    message("[dry-run] would write input status: ", path)
    return(invisible(path))
  }
  dir_create(dirname(path))
  utils::write.csv(rows, path, row.names = FALSE, na = "")
  invisible(path)
}

finish_script_runtime <- function(runtime,
                                  manifest_path,
                                  outputs = character(),
                                  inputs = character(),
                                  status = "completed",
                                  notes = character()) {
  if (isTRUE(runtime$dry_run)) {
    message("[dry-run] would write run manifest: ", manifest_path)
    return(invisible(manifest_path))
  }
  write_run_manifest(
    manifest_path,
    inputs = inputs,
    outputs = outputs,
    parameters = list(
      script = runtime$script,
      stage = runtime$stage,
      dataset = runtime$dataset,
      dry_run = runtime$dry_run,
      status = status,
      started_at = format(runtime$started_at, "%Y-%m-%d %H:%M:%S %Z"),
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ),
    notes = notes
  )
  invisible(manifest_path)
}
