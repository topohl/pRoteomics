# Shared runtime helpers for numbered pipeline entrypoint scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("input_addressability", mode = "function")) {
  source(repo_path("R", "paths.R"))
}
if (!exists("current_dataset_from_cli", mode = "function")) {
  source(repo_path("R", "dataset_config.R"))
}
if (!exists("write_run_manifest", mode = "function")) {
  source(repo_path("R", "validation_utils.R"))
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
  ## The ledger now records WHY an input is unusable, not just that it is.
  ## The previous vocabulary folded three different conditions into
  ## missing_required / missing_optional, which made the recorded state
  ## disagree with the run: an input behind an unmounted root, or one whose
  ## path is past the limit R can open, was written down as simply missing.
  ## `required` is already its own column, so the failure classes no longer
  ## have to carry it.
  addressability <- input_addressability(path)[[1]]
  ## Some callers pass a glob rather than a literal path. A pattern that
  ## matches something is present whatever the literal string resolves to.
  if (!identical(addressability, INPUT_STATUS_PRESENT) &&
      length(Sys.glob(path)) > 0L) {
    addressability <- INPUT_STATUS_PRESENT
  }
  present <- identical(addressability, INPUT_STATUS_PRESENT)
  data.frame(
    dataset = dataset,
    input_name = input_name,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    required = isTRUE(required),
    status = status %||% addressability,
    ## Derived, by definition, from status == present. It exists so a caller
    ## that only needs a yes/no cannot reintroduce the four-into-one collapse.
    input_present = present,
    message = message %||% input_status_message(addressability),
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
