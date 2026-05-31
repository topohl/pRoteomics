#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
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

stage_registry <- list(
  core = data.frame(
    script = c(
      "01_preprocessing/03_gct_extractR.r",
      "02_id_mapping/01_MapThatProt_batch.r"
    ),
    required = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  ),
  qc = data.frame(
    script = c(
      "03_qc_exploration/00_dataset_qc_report.r",
      "03_qc_exploration/02_missingness_diagnostics.r",
      "03_qc_exploration/05_pca_confounding_qc.r",
      "03_qc_exploration/07_qc_biology_confounding_report.r"
    ),
    required = c(FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  ),
  enrichment = data.frame(
    script = c(
      "04_differential_expression_enrichment/01_clusterProfiler.r",
      "04_differential_expression_enrichment/02_compareGO.r",
      "04_differential_expression_enrichment/03_biological_program_summary.r",
      "05_celltype_enrichment_EWCE/01_EWCE_E9.r"
    ),
    required = c(TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  ),
  modules = data.frame(
    script = c(
      "06_modules_WGCNA/01_WGCNA.r",
      "06_modules_WGCNA/05_wgcna_de_gsea_overlap.r",
      "06_modules_WGCNA/91_module_score.r"
    ),
    required = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  ),
  networks = data.frame(
    script = c(
      "07_spatial_networks/01_network_spatial_relations.r",
      "07_spatial_networks/02_differential_networks.r",
      "07_spatial_networks/03_bootstrap_network_stability.r",
      "07_spatial_networks/04_bootstrap_differential_network_stability.r",
      "07_spatial_networks/05_bootstrap_differential_network_figures.r",
      "07_spatial_networks/06_chord_diagram.r"
    ),
    required = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  ),
  behavior = data.frame(
    script = c("08_behavior_physio_coupling/02_network_behavior_coupling.r"),
    required = c(FALSE),
    stringsAsFactors = FALSE
  ),
  export = data.frame(
    script = c(
      "09_export_pride_journal/01_make_pride_manifest.R",
      "09_export_pride_journal/02_make_sample_metadata.R",
      "09_export_pride_journal/03_make_supplementary_tables.R",
      "09_export_pride_journal/04_validate_pride_submission.R",
      "09_export_pride_journal/05_make_methods_summary.R",
      "09_export_pride_journal/06_make_biological_claims_table.R"
    ),
    required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
)

valid_stages <- c(names(stage_registry), "all")
if (!stage_arg %in% valid_stages) {
  stop("Unsupported --stage: ", stage_arg, ". Use one of: ", paste(valid_stages, collapse = ", "), call. = FALSE)
}

if (isTRUE(list_stages)) {
  cat("Available stages:\n")
  for (nm in names(stage_registry)) {
    cat("\n", nm, "\n", sep = "")
    print(stage_registry[[nm]], row.names = FALSE)
  }
  quit(status = 0, save = "no")
}

selected_stages <- if (identical(stage_arg, "all")) names(stage_registry) else stage_arg
steps <- do.call(rbind, lapply(selected_stages, function(stage) {
  transform(stage_registry[[stage]], stage = stage)
}))
steps <- steps[, c("stage", "script", "required")]

manifest_dir <- path_results("logs", "pipeline")
dir_create(manifest_dir)
manifest_path <- file.path(
  manifest_dir,
  paste0("pipeline_manifest_", dataset, "_", safe_name(stage_arg), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
)

cat("Dataset pipeline\n")
cat("Resolved dataset:", dataset, "\n")
cat("Stage:", stage_arg, "\n")
cat("Dry run:", dry_run, "\n")
cat("Project root:", repo_root(), "\n")
cat("Manifest:", manifest_path, "\n\n")

run_one_step <- function(stage, script, required) {
  script_path <- repo_path(script)
  started_at <- Sys.time()
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

results <- do.call(rbind, lapply(seq_len(nrow(steps)), function(i) {
  res <- run_one_step(steps$stage[[i]], steps$script[[i]], steps$required[[i]])
  res
}))

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
