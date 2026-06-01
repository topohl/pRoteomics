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

stage_arg <- arg_value("--stage", default = "all")
dataset_arg <- arg_value("--dataset", default = current_dataset())
dataset <- validate_dataset(dataset_arg, source = "--dataset")
dry_run <- has_flag("--dry-run") || tolower(Sys.getenv("PROTEOMICS_DRY_RUN", unset = "")) %in% c("1", "true", "yes")
list_stages <- has_flag("--list-stages")

Sys.setenv(PROTEOMICS_DATASET = dataset)
if (isTRUE(dry_run)) Sys.setenv(PROTEOMICS_DRY_RUN = "true")

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
validate_pipeline_scripts_exist(registry)
validate_run_order_against_registry(registry)

valid_stages <- c(pipeline_stage_names(registry), "all")
if (!stage_arg %in% valid_stages) {
  stop("Unsupported --stage: ", stage_arg, ". Use one of: ", paste(valid_stages, collapse = ", "), call. = FALSE)
}

if (isTRUE(list_stages)) {
  cat("Available stages from pipeline.yml:\n")
  for (nm in pipeline_stage_names(registry)) {
    cat("\n", nm, "\n", sep = "")
    print(pipeline_steps(registry, nm, dataset = dataset)[, c("script", "required", "supported_datasets")], row.names = FALSE)
  }
  quit(status = 0, save = "no")
}

selected_stages <- if (identical(stage_arg, "all")) pipeline_stage_names(registry) else stage_arg
steps <- pipeline_steps(registry, selected_stages, dataset = dataset)

manifest_dir <- path_results("logs", "pipeline")
dir_create(manifest_dir)
manifest_path <- file.path(
  manifest_dir,
  paste0("pipeline_manifest_", dataset, "_", safe_name(stage_arg), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
)

cat("Dataset pipeline\n")
cat("Resolved dataset:", dataset, "\n")
cat("Dataset interpretation:", dataset_interpretation(dataset), "\n")
cat("Stage:", stage_arg, "\n")
cat("Dry run:", dry_run, "\n")
cat("Project root:", repo_root(), "\n")
cat("Registry:", repo_path("pipeline.yml"), "\n")
cat("Manifest:", manifest_path, "\n\n")

run_one_step <- function(stage, script, required, supported) {
  started_at <- Sys.time()
  if (!isTRUE(supported)) {
    message("[SKIP] ", script, " does not support dataset ", dataset)
    return(data.frame(
      dataset = dataset,
      selected_stage = stage_arg,
      stage = stage,
      script = script,
      required = required,
      status = "skipped_unsupported_dataset",
      exit_code = 0L,
      started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
      finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    ))
  }

  script_path <- repo_path(script)
  if (!file.exists(script_path)) {
    status <- if (isTRUE(required)) "missing_required" else "missing_optional"
    message(if (isTRUE(required)) "[FAIL] " else "[WARN] ", "Pipeline script missing: ", script)
    return(data.frame(
      dataset = dataset,
      selected_stage = stage_arg,
      stage = stage,
      script = script,
      required = required,
      status = status,
      exit_code = if (isTRUE(required)) 127L else 0L,
      started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
      finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    ))
  }

  step_args <- c(script_path)
  if (isTRUE(dry_run)) step_args <- c(step_args, "--dry-run")
  cat("==> Running ", stage, " :: ", script, "\n", sep = "")
  exit_code <- system2("Rscript", step_args)
  if (is.null(exit_code)) exit_code <- 0L
  status <- if (identical(exit_code, 0L)) "passed" else if (isTRUE(required)) "failed_required" else "failed_optional"
  data.frame(
    dataset = dataset,
    selected_stage = stage_arg,
    stage = stage,
    script = script,
    required = required,
    status = status,
    exit_code = as.integer(exit_code),
    started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
    finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
}

results <- data.frame()
for (i in seq_len(nrow(steps))) {
  res <- run_one_step(steps$stage[[i]], steps$script[[i]], steps$required[[i]], steps$supported[[i]])
  results <- rbind(results, res)
  if (identical(res$status, "failed_required") || identical(res$status, "missing_required")) {
    message("[FAIL] Required step failed; stopping pipeline before downstream stages.")
    break
  }
}

utils::write.csv(results, manifest_path, row.names = FALSE)

cat("\nPipeline summary\n")
print(as.data.frame(table(results$stage, results$status)), row.names = FALSE)
cat("\nManifest written:", manifest_path, "\n")

failed_required <- results[results$status %in% c("missing_required", "failed_required"), , drop = FALSE]
if (nrow(failed_required)) {
  cat("\nRequired step failures:\n")
  print(failed_required[, c("stage", "script", "status", "exit_code")], row.names = FALSE)
  quit(status = 1, save = "no")
}

cat("\nDataset pipeline finished successfully for required steps.\n")
