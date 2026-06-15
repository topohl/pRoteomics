#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "pipeline_registry.R"))
setwd(repo_root())

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}

has_flag <- function(flag) flag %in% args

split_arg <- function(value) {
  value <- trimws(as.character(value))
  if (!nzchar(value) || identical(value, "none")) return(character())
  unique(trimws(unlist(strsplit(value, ","), use.names = FALSE)))
}

split_env_arg <- function(value) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) return(character())
  pairs <- unlist(strsplit(value, "\\|", fixed = FALSE), use.names = FALSE)
  pairs <- pairs[nzchar(pairs)]
  stats::setNames(
    sub("^[^=]*=", "", pairs),
    sub("=.*$", "", pairs)
  )
}

stage_arg <- arg_value("--stage", default = "all")
from_stage <- arg_value("--from-stage", default = "")
dataset_arg <- normalize_dataset(arg_value("--dataset", default = current_dataset()))
exclude_stages <- split_arg(arg_value("--exclude-stage", default = ""))
exclude_scripts <- gsub("\\\\", "/", split_arg(arg_value("--exclude-script", default = "")))
dry_run <- has_flag("--dry-run") || has_flag("--plan") ||
  tolower(Sys.getenv("PROTEOMICS_DRY_RUN", unset = "")) %in% c("1", "true", "yes")
strict_outputs <- has_flag("--strict-outputs") ||
  tolower(Sys.getenv("PROTEOMICS_STRICT_OUTPUTS", unset = "")) %in% c("1", "true", "yes")
list_stages <- has_flag("--list-stages")

if (!identical(dataset_arg, "all")) dataset_arg <- validate_dataset(dataset_arg, source = "--dataset")
if (isTRUE(dry_run)) Sys.setenv(PROTEOMICS_DRY_RUN = "true")

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
validate_pipeline_scripts_exist(registry, fail = TRUE)
validate_run_order_against_registry(registry)

stage_names <- pipeline_stage_names(registry)
valid_stages <- c(stage_names, "all")

selected_stages <- if (identical(stage_arg, "all")) {
  stage_names
} else {
  split_arg(stage_arg)
}

unknown_stages <- setdiff(c(selected_stages, exclude_stages, from_stage[nzchar(from_stage)]), valid_stages)
if (length(unknown_stages)) {
  stop("Unsupported stage(s): ", paste(unknown_stages, collapse = ", "), ". Use one of: ", paste(valid_stages, collapse = ", "), call. = FALSE)
}

if (nzchar(from_stage)) {
  if (identical(from_stage, "all")) {
    selected_stages <- stage_names
  } else {
    start <- match(from_stage, stage_names)
    selected_stages <- stage_names[start:length(stage_names)]
  }
  if (!identical(stage_arg, "all")) selected_stages <- intersect(selected_stages, split_arg(stage_arg))
}

selected_stages <- setdiff(selected_stages, exclude_stages)
if (!length(selected_stages) && !isTRUE(list_stages)) {
  stop("No stages selected after applying --stage/--from-stage/--exclude-stage.", call. = FALSE)
}

if (isTRUE(list_stages)) {
  cat("Available stages from pipeline.yml:\n")
  for (nm in stage_names) {
    cat("\n", nm, "\n", sep = "")
    print(
      pipeline_steps(registry, nm, dataset = "all", include_unsupported = TRUE)[,
        c("stage", "script", "scope", "required", "supported_datasets", "recomputes_core_state", "safe_downstream_rerun")
      ],
      row.names = FALSE
    )
  }
  quit(status = 0, save = "no")
}

steps <- pipeline_steps(registry, selected_stages, dataset = dataset_arg, include_unsupported = TRUE)
if (length(exclude_scripts)) {
  steps <- steps[!steps$script %in% exclude_scripts, , drop = FALSE]
}
if (!nrow(steps)) {
  stop("No scripts selected after applying exclusions.", call. = FALSE)
}

manifest_dir <- path_results("logs", "pipeline")
dir_create(manifest_dir)
selected_label <- paste(selected_stages, collapse = "+")
manifest_path <- file.path(
  manifest_dir,
  paste0("pipeline_manifest_", safe_name(dataset_arg), "_", safe_name(selected_label), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
)

fmt_time <- function(x) format(x, "%Y-%m-%d %H:%M:%S %Z")

path_status <- function(paths, dataset) {
  paths <- split_registry_paths(paths)
  if (!length(paths)) return(data.frame(path = character(), exists = logical()))
  data.frame(
    path = gsub("<dataset>", dataset, paths, fixed = TRUE),
    exists = vapply(paths, registry_path_exists, logical(1), dataset = dataset),
    stringsAsFactors = FALSE
  )
}

join_paths <- function(x) paste(x, collapse = "|")

compact_paths <- function(x, max_items = 4L) {
  x <- x[nzchar(x)]
  if (!length(x)) return("none")
  shown <- head(x, max_items)
  suffix <- if (length(x) > max_items) paste0(" | ... +", length(x) - max_items, " more") else ""
  paste0(paste(shown, collapse = " | "), suffix)
}

expand_existing_output_paths <- function(paths, dataset) {
  paths <- split_registry_paths(paths)
  if (!length(paths)) return(character())
  resolved <- unlist(lapply(paths, function(path) {
    p <- resolve_registry_path(path, dataset = dataset)
    if (grepl("[*?\\[]", p)) return(Sys.glob(p))
    p
  }), use.names = FALSE)
  unique(resolved[file.exists(resolved) & !dir.exists(resolved)])
}

planned_status <- function(step, missing_required, missing_optional, missing_outputs) {
  if (!isTRUE(step$supported)) return("skipped_unsupported_dataset")
  if (length(missing_required) && isTRUE(step$required)) return(if (isTRUE(dry_run)) "would_fail_missing_input" else "missing_required_input")
  if (length(missing_required)) return("skipped_missing_input")
  if (isTRUE(dry_run)) return("would_run")
  "pending"
}

make_result <- function(step, status, exit_code, started_at, finished_at, missing_required, missing_optional, missing_after, produced_before, produced_after, output_validation = NULL) {
  if (is.null(output_validation) || !nrow(output_validation)) {
    output_validation <- data.frame(path = character(), validation_status = character(), validation_message = character(), stringsAsFactors = FALSE)
  }
  context <- run_context_metadata()
  data.frame(
    dataset = step$dataset,
    selected_stage = selected_label,
    stage = step$stage,
    script = step$script,
    scope = step$scope,
    required = step$required,
    supported = step$supported,
    consumes_required = step$consumes_required,
    missing_required_inputs = join_paths(missing_required),
    consumes_optional = step$consumes_optional,
    missing_optional_inputs = join_paths(missing_optional),
    produces = step$produces,
    produced_outputs_before_run = join_paths(produced_before),
    produced_outputs_after_run = join_paths(produced_after),
    missing_expected_outputs_after_run = join_paths(missing_after),
    output_validation_status = join_paths(unique(output_validation$validation_status)),
    output_validation_messages = join_paths(output_validation$validation_message[nzchar(output_validation$validation_message)]),
    env_overrides = step$env %||% "",
    recomputes_core_state = step$recomputes_core_state,
    safe_downstream_rerun = step$safe_downstream_rerun,
    notes = step$notes,
    status = status,
    exit_code = as.integer(exit_code),
    started_at = fmt_time(started_at),
    finished_at = fmt_time(finished_at),
    git_commit = context$git_commit,
    r_version = context$r_version,
    platform = context$platform,
    env_PROTEOMICS_DATASET = Sys.getenv("PROTEOMICS_DATASET", unset = NA_character_),
    env_PROTEOMICS_DRY_RUN = Sys.getenv("PROTEOMICS_DRY_RUN", unset = NA_character_),
    env_PROTEOMICS_RECOMPUTE = Sys.getenv("PROTEOMICS_RECOMPUTE", unset = NA_character_),
    env_PROTEOMICS_WGCNA_FORCE_FULL = Sys.getenv("PROTEOMICS_WGCNA_FORCE_FULL", unset = NA_character_),
    env_PROTEOMICS_MODULE_DEFINITION_SOURCE = Sys.getenv("PROTEOMICS_MODULE_DEFINITION_SOURCE", unset = NA_character_),
    stringsAsFactors = FALSE
  )
}

cat("Dataset pipeline\n")
cat("Resolved dataset:", dataset_arg, "\n")
if (!identical(dataset_arg, "all")) cat("Dataset interpretation:", dataset_interpretation(dataset_arg), "\n")
cat("Selected stages:", paste(selected_stages, collapse = ", "), "\n")
cat("Excluded stages:", ifelse(length(exclude_stages), paste(exclude_stages, collapse = ", "), "none"), "\n")
cat("Excluded scripts:", ifelse(length(exclude_scripts), paste(exclude_scripts, collapse = ", "), "none"), "\n")
cat("Dry run / plan:", dry_run, "\n")
cat("Strict expected-output policy:", strict_outputs, "\n")
cat("Project root:", repo_root(), "\n")
cat("Registry:", repo_path("pipeline.yml"), "\n")
cat("Manifest:", manifest_path, "\n\n")

results <- data.frame()
for (i in seq_len(nrow(steps))) {
  step <- steps[i, , drop = FALSE]
  started_at <- Sys.time()
  Sys.setenv(PROTEOMICS_DATASET = if (identical(step$dataset, "global")) "all" else step$dataset)

  required_status <- path_status(step$consumes_required, step$dataset)
  optional_status <- path_status(step$consumes_optional, step$dataset)
  output_status_before <- path_status(step$produces, step$dataset)
  missing_required <- required_status$path[!is.na(required_status$exists) & !required_status$exists]
  missing_optional <- optional_status$path[!is.na(optional_status$exists) & !optional_status$exists]
  produced_before <- output_status_before$path[!is.na(output_status_before$exists) & output_status_before$exists]
  missing_outputs_before <- output_status_before$path[!is.na(output_status_before$exists) & !output_status_before$exists]

  status <- planned_status(step, missing_required, missing_optional, missing_outputs_before)
  if (isTRUE(dry_run) || !identical(status, "pending")) {
    cat(
      "[PLAN] ", step$stage, " :: ", step$script,
      " | dataset=", step$dataset,
      " | scope=", step$scope,
      " | status=", status,
      " | recomputes_core_state=", step$recomputes_core_state,
      "\n",
      sep = ""
    )
    cat("       required inputs missing: ", compact_paths(missing_required), "\n", sep = "")
    cat("       optional inputs missing: ", compact_paths(missing_optional), "\n", sep = "")
    cat("       expected outputs: ", compact_paths(path_status(step$produces, step$dataset)$path), "\n", sep = "")
    cat("       expected outputs missing now: ", compact_paths(missing_outputs_before), "\n", sep = "")
    res <- make_result(step, status, 0L, started_at, Sys.time(), missing_required, missing_optional, missing_outputs_before, produced_before, produced_before)
    results <- rbind(results, res)
    if (identical(status, "missing_required_input")) break
    next
  }

  script_path <- repo_path(step$script)
  if (!file.exists(script_path)) {
    status <- if (isTRUE(step$required)) "missing_required_script" else "missing_optional_script"
    res <- make_result(step, status, if (isTRUE(step$required)) 127L else 0L, started_at, Sys.time(), missing_required, missing_optional, missing_outputs_before, produced_before, produced_before)
    results <- rbind(results, res)
    if (isTRUE(step$required)) break
    next
  }

  if (length(missing_optional)) {
    message("[WARN] Optional input(s) missing for ", step$script, ": ", join_paths(missing_optional))
  }

  cat("==> Running ", step$stage, " :: ", step$script, " [", step$dataset, "]\n", sep = "")
  step_env <- split_env_arg(step$env)
  old_env <- stats::setNames(vapply(names(step_env), Sys.getenv, character(1), unset = NA_character_), names(step_env))
  if (length(step_env)) {
    do.call(Sys.setenv, as.list(step_env))
    cat("    env overrides: ", paste(paste(names(step_env), step_env, sep = "="), collapse = ", "), "\n", sep = "")
  }
  exit_code <- tryCatch(system2("Rscript", script_path), finally = {
    if (length(old_env)) {
      for (nm in names(old_env)) {
        if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, stats::setNames(as.list(old_env[[nm]]), nm))
      }
    }
  })
  if (is.null(exit_code)) exit_code <- 0L

  output_status_after <- path_status(step$produces, step$dataset)
  produced_after <- output_status_after$path[!is.na(output_status_after$exists) & output_status_after$exists]
  missing_after <- output_status_after$path[!is.na(output_status_after$exists) & !output_status_after$exists]
  output_validation <- validate_pipeline_outputs(expand_existing_output_paths(step$produces, step$dataset), dataset = step$dataset)
  validation_failed <- any(output_validation$validation_status %in% c("warning", "error", "missing"))
  status <- if (identical(exit_code, 0L)) "passed" else if (isTRUE(step$required)) "failed_required" else "failed_optional"
  if (length(missing_after)) {
    if (isTRUE(step$required)) {
      status <- "failed_required_missing_outputs"
      exit_code <- 1L
    } else if (isTRUE(strict_outputs)) {
      status <- "failed_optional_missing_outputs"
      if (identical(exit_code, 0L)) exit_code <- 1L
    } else if (identical(exit_code, 0L)) {
      status <- "passed_with_missing_outputs"
    } else {
      status <- "failed_optional_missing_outputs"
    }
    message("[WARN] Expected output(s) missing after ", step$script, ": ", join_paths(missing_after))
  }
  if (identical(status, "passed") && isTRUE(validation_failed)) {
    if (isTRUE(step$required)) {
      status <- "failed_required_output_validation"
      exit_code <- 1L
    } else {
      status <- "passed_with_output_warnings"
    }
    message("[WARN] Output validation for ", step$script, ": ", join_paths(output_validation$validation_message[nzchar(output_validation$validation_message)]))
  }
  res <- make_result(step, status, exit_code, started_at, Sys.time(), missing_required, missing_optional, missing_after, produced_before, produced_after, output_validation)
  results <- rbind(results, res)
  if (status %in% c("failed_required", "failed_required_missing_outputs", "failed_required_output_validation")) {
    message("[FAIL] Required step failed; stopping pipeline before downstream stages.")
    break
  }
}

utils::write.csv(results, manifest_path, row.names = FALSE)

cat("\nPipeline summary\n")
print(as.data.frame(table(results$stage, results$status)), row.names = FALSE)
cat("\nManifest written:", manifest_path, "\n")
cat("Registry audit tables:\n")
cat(" - ", file.path(manifest_dir, "unregistered_scripts.csv"), "\n", sep = "")
cat(" - ", file.path(manifest_dir, "missing_registered_scripts.csv"), "\n", sep = "")

failed_required <- results[results$status %in% c("missing_required_input", "missing_required_script", "failed_required", "failed_required_missing_outputs", "failed_required_output_validation"), , drop = FALSE]
if (nrow(failed_required)) {
  cat("\nRequired step failures:\n")
  print(failed_required[, c("dataset", "stage", "script", "status", "exit_code", "missing_required_inputs")], row.names = FALSE)
  quit(status = 1, save = "no")
}

failed_optional_strict <- results[results$status %in% c("failed_optional_missing_outputs") & isTRUE(strict_outputs), , drop = FALSE]
if (nrow(failed_optional_strict)) {
  cat("\nOptional step output failures under --strict-outputs:\n")
  print(failed_optional_strict[, c("dataset", "stage", "script", "status", "exit_code", "missing_expected_outputs_after_run")], row.names = FALSE)
  quit(status = 1, save = "no")
}

if (isTRUE(dry_run)) {
  cat("\nDry-run plan completed; no scripts were executed.\n")
} else {
  cat("\nDataset pipeline finished successfully for required steps.\n")
}
