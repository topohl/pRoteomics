#!/usr/bin/env Rscript

# Prepare dataset-specific ProTigy inputs at the biological-replicate level.
# Existing imputed log2 values are averaged across available valid Left/Right
# hemispheres within AnimalID x canonical spatial unit. A single valid hemisphere
# is retained without hemisphere imputation. No transformation, normalization,
# protein remapping, or downstream analysis occurs.

suppressPackageStartupMessages({
  required_packages <- c("readxl", "writexl")
  missing_packages <- required_packages[!vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]
  if (length(missing_packages)) {
    stop(
      "Missing required package(s): ", paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
})

source(if (file.exists(file.path("R", "script_runtime.R"))) {
  file.path("R", "script_runtime.R")
} else {
  file.path("..", "R", "script_runtime.R")
})
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "protigy_input_utils.R"))

SCRIPT_ID <- "01_preprocessing/02a_prepare_animal_level_protigy_input.r"
MODULE_ID <- "01_preprocessing/02a_prepare_animal_level_protigy_input"
runtime <- init_script_runtime(
  script = SCRIPT_ID,
  stage = "external_protigy_input_handoff",
  default_dataset = "neuron_neuropil",
  allow_all = TRUE
)
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

datasets <- if (identical(runtime$dataset, "all")) valid_datasets() else runtime$dataset
if (
  identical(runtime$dataset, "all") &&
    nzchar(Sys.getenv("PROTEOMICS_PROTIGY_INPUT_EXPRESSION_XLSX", unset = ""))
) {
  stop(
    "The general PROTEOMICS_PROTIGY_INPUT_EXPRESSION_XLSX override cannot be used with --dataset all. ",
    "Use dataset-specific overrides ending in _NEURON_NEUROPIL, _NEURON_SOMA, and _MICROGLIA.",
    call. = FALSE
  )
}

output_paths <- function(dataset) {
  output_dir <- path_processed(
    "01_preprocessing",
    "protigy_input_animal_level",
    dataset
  )
  table_dir <- path_results("tables", MODULE_ID, dataset)
  log_dir <- path_results("logs", MODULE_ID, dataset)
  list(
    output_dir = output_dir,
    table_dir = table_dir,
    log_dir = log_dir,
    gct = file.path(output_dir, paste0(dataset, "_animal_level.gct")),
    strict_gct = file.path(
      output_dir,
      paste0(dataset, "_animal_level_complete_bilateral_sensitivity.gct")
    ),
    xlsx = file.path(output_dir, paste0(dataset, "_animal_level.xlsx")),
    audit = file.path(table_dir, "aggregation_audit.csv"),
    summary = file.path(table_dir, "aggregation_summary.csv"),
    assignment = file.path(table_dir, "source_sample_assignment.csv"),
    manifest = file.path(log_dir, "run_manifest.yml")
  )
}

resolve_one <- function(dataset, record_resolution) {
  resolve_dataset_inputs(
    dataset = dataset,
    purpose = "protigy_input",
    script = SCRIPT_ID,
    stage = runtime$stage,
    record_resolution = record_resolution
  )
}

dry_run_one <- function(dataset) {
  inputs <- resolve_one(dataset, record_resolution = FALSE)
  paths <- output_paths(dataset)
  expression_exists <- file.exists(inputs$expression_file)
  metadata_exists <- file.exists(inputs$metadata_file)
  matching_samples <- NA_integer_
  included_matching_samples <- NA_integer_
  if (expression_exists && metadata_exists) {
    expression_header <- as.data.frame(
      readxl::read_excel(inputs$expression_file, n_max = 0L),
      check.names = FALSE
    )
    expression_header <- protigy_prepare_legacy_expression(expression_header)
    sample_cols <- setdiff(names(expression_header), c("id", "Description"))
    metadata <- as.data.frame(readxl::read_excel(inputs$metadata_file))
    matching_samples <- sum(sample_cols %in% metadata$sample_id)
    matched <- metadata[match(sample_cols, metadata$sample_id), , drop = FALSE]
    excluded <- if ("exclude" %in% names(matched)) {
      protigy_excluded(matched$exclude)
    } else if ("Exclude" %in% names(matched)) {
      protigy_excluded(matched$Exclude)
    } else {
      rep(FALSE, nrow(matched))
    }
    included_matching_samples <- sum(!is.na(matched$sample_id) & !excluded)
  }
  dry_run_line("Resolved dataset", dataset)
  dry_run_line(
    "Imputed expression source",
    inputs$expression_file,
    if (expression_exists) "PASS" else "FAIL"
  )
  dry_run_line(
    "Metadata source",
    inputs$metadata_file,
    if (metadata_exists) "PASS" else "FAIL"
  )
  dry_run_line("Matching expression/metadata samples", matching_samples)
  dry_run_line("Matching samples after exclusions", included_matching_samples)
  dry_run_line("Canonical spatial unit", dataset_spatial_unit(dataset))
  dry_run_line("Expected GCT", paths$gct)
  dry_run_line("Expected strict-bilateral sensitivity GCT", paths$strict_gct)
  dry_run_line("Expected XLSX", paths$xlsx)
  dry_run_line("Expected audit directory", paths$table_dir)
  dry_run_line("Expected provenance directory", paths$log_dir)
  expression_exists && metadata_exists && !is.na(matching_samples)
}

if (isTRUE(runtime$dry_run)) {
  ok <- vapply(datasets, dry_run_one, logical(1))
  quit(status = if (all(ok)) 0L else 1L, save = "no")
}

make_summary <- function(bundle, inputs) {
  retained_assignment <- bundle$source_sample_assignment[
    bundle$source_sample_assignment$inclusion_status != "excluded_by_metadata",
    ,
    drop = FALSE
  ]
  animal_groups <- unique(retained_assignment[, c("AnimalID", "ExpGroup"), drop = FALSE])
  group_count <- function(group) sum(animal_groups$ExpGroup == group)
  audit <- bundle$aggregation_audit
  complete <- audit$hemisphere_status == "bilateral_complete"
  left_only <- audit$hemisphere_status == "left_only_observed"
  right_only <- audit$hemisphere_status == "right_only_observed"
  missing_both <- audit$hemisphere_status == "missing_both"
  data.frame(
    dataset = bundle$dataset,
    source_expression_file = normalizePath(inputs$expression_file, winslash = "/", mustWork = FALSE),
    source_expression_hash = file_hash_sha256(inputs$expression_file),
    metadata_file = normalizePath(inputs$metadata_file, winslash = "/", mustWork = FALSE),
    metadata_hash = file_hash_sha256(inputs$metadata_file),
    n_source_expression_rows = bundle$n_source_expression_rows,
    n_protigy_proteins_before_aggregation = bundle$n_proteins,
    n_proteins = bundle$n_proteins,
    n_protigy_proteins_after_aggregation = bundle$n_proteins,
    row_descriptor_name = "Description",
    description_source = bundle$description_source,
    n_description_fallback_to_id = bundle$description_fallback_count,
    n_input_samples = bundle$n_input_samples,
    n_output_animal_spatial_units = if (isTRUE(bundle$valid)) nrow(bundle$output_metadata) else 0L,
    n_unique_animals = length(unique(retained_assignment$AnimalID)),
    animals_group_1 = group_count("1"),
    animals_group_2 = group_count("2"),
    animals_group_3 = group_count("3"),
    n_complete_bilateral_units = sum(complete),
    n_incomplete_bilateral_units = sum(!complete),
    n_single_hemisphere_units = sum(left_only | right_only),
    n_left_only_units = sum(left_only),
    n_right_only_units = sum(right_only),
    n_missing_both_units = sum(missing_both),
    n_duplicate_left_units = sum(audit$n_left_source_samples > 1L),
    n_duplicate_right_units = sum(audit$n_right_source_samples > 1L),
    validation_status = bundle$validation_status,
    stringsAsFactors = FALSE
  )
}

run_one <- function(dataset) {
  Sys.setenv(PROTEOMICS_DATASET = dataset)
  inputs <- resolve_one(dataset, record_resolution = TRUE)
  paths <- output_paths(dataset)
  expression <- as.data.frame(
    readxl::read_excel(inputs$expression_file),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  metadata <- as.data.frame(
    readxl::read_excel(inputs$metadata_file),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  bundle <- protigy_prepare_animal_level(
    expression = expression,
    metadata = metadata,
    dataset = dataset,
    validate_e9_design = TRUE,
    fail_on_invalid = FALSE
  )
  summary <- make_summary(bundle, inputs)

  dir_create(paths$output_dir)
  dir_create(paths$table_dir)
  dir_create(paths$log_dir)
  utils::write.csv(bundle$aggregation_audit, paths$audit, row.names = FALSE, na = "")
  utils::write.csv(summary, paths$summary, row.names = FALSE, na = "")
  utils::write.csv(
    bundle$source_sample_assignment,
    paths$assignment,
    row.names = FALSE,
    na = ""
  )

  if (isTRUE(bundle$valid)) {
    writexl::write_xlsx(bundle$annotated_table, paths$xlsx)
    write_gct_v1.3(bundle$annotated_table, paths$gct, bundle$output_metadata)
    write_gct_v1.3(
      bundle$strict_bilateral_annotated_table,
      paths$strict_gct,
      bundle$strict_bilateral_metadata
    )
  }

  outputs <- c(
    aggregation_audit = paths$audit,
    aggregation_summary = paths$summary,
    source_sample_assignment = paths$assignment
  )
  if (isTRUE(bundle$valid)) {
    outputs <- c(
      outputs,
      animal_level_gct = paths$gct,
      complete_bilateral_sensitivity_gct = paths$strict_gct,
      animal_level_xlsx = paths$xlsx
    )
  }
  write_run_manifest(
    paths$manifest,
    inputs = list(
      imputed_expression = inputs$expression_file,
      sample_metadata = inputs$metadata_file
    ),
    outputs = as.list(outputs),
    parameters = list(
      dataset = dataset,
      canonical_spatial_unit = bundle$spatial_unit,
      biological_replicate = "AnimalID",
      protigy_row_descriptor = "Description",
      protigy_row_descriptor_count = 1L,
      protein_description_source = bundle$description_source,
      protein_description_fallback_to_id_count = bundle$description_fallback_count,
      primary_aggregation_policy = paste(
        "mean both valid hemispheres when bilateral; retain one valid hemisphere",
        "without hemisphere imputation when one-sided"
      ),
      bilateral_aggregation_method = "equal_weight_mean_LR_on_existing_imputed_log2_values",
      single_hemisphere_aggregation_method = "single_observed_hemisphere_no_imputation",
      fail_on_duplicate_same_side_samples = TRUE,
      strict_complete_bilateral_sensitivity = basename(paths$strict_gct),
      validation_status = bundle$validation_status
    ),
    notes = c(
      "Dataset-specific 01_impute.r output; joint/global QC GCT is not used.",
      "ProTigy handoff deliberately emits one Description row descriptor because the current Broad ProTigy/cmapR workflow is not compatible downstream with nrhd=0.",
      "GCT v1.3 can structurally permit nrhd=0; this is a ProTigy interoperability constraint, not a global GCT validity claim.",
      "Description uses First.Protein.Description when available; missing or absent descriptions fall back to id without remapping the id.",
      "No normalization, imputation, transformation, protein remapping, or downstream analysis is performed.",
      "Primary input retains one-sided observed units; the sensitivity GCT contains complete bilateral pairs only.",
      "ProTigy remains an external/manual boundary; this stage is not registered in pipeline.yml."
    )
  )

  message(
    "Dataset ", dataset,
    ": proteins=", bundle$n_proteins,
    ", input samples=", bundle$n_input_samples,
    ", primary animal/spatial units=", nrow(bundle$output_metadata),
    ", complete bilateral sensitivity units=", nrow(bundle$strict_bilateral_metadata),
    ", validation=", bundle$validation_status
  )
  list(
    dataset = dataset,
    success = isTRUE(bundle$valid),
    validation_status = bundle$validation_status,
    gct_path = paths$gct,
    gct_hash = if (isTRUE(bundle$valid)) file_hash_sha256(paths$gct) else NA_character_,
    strict_gct_path = paths$strict_gct,
    strict_gct_hash = if (isTRUE(bundle$valid)) file_hash_sha256(paths$strict_gct) else NA_character_,
    summary = summary
  )
}

results <- lapply(datasets, function(dataset) {
  tryCatch(
    run_one(dataset),
    error = function(e) {
      list(
        dataset = dataset,
        success = FALSE,
        validation_status = paste0("FAIL: ", conditionMessage(e)),
        gct_path = output_paths(dataset)$gct,
        gct_hash = NA_character_,
        strict_gct_path = output_paths(dataset)$strict_gct,
        strict_gct_hash = NA_character_,
        summary = NULL
      )
    }
  )
})

if (identical(runtime$dataset, "all")) {
  handoff <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      dataset = x$dataset,
      gct_path = normalizePath(x$gct_path, winslash = "/", mustWork = FALSE),
      gct_sha256 = x$gct_hash,
      complete_bilateral_sensitivity_gct_path = normalizePath(
        x$strict_gct_path,
        winslash = "/",
        mustWork = FALSE
      ),
      complete_bilateral_sensitivity_gct_sha256 = x$strict_gct_hash,
      validation_status = x$validation_status,
      stringsAsFactors = FALSE
    )
  }))
  handoff_path <- path_results("tables", MODULE_ID, "animal_level_protigy_handoff_manifest.csv")
  dir_create(dirname(handoff_path))
  utils::write.csv(handoff, handoff_path, row.names = FALSE, na = "")
  message("Wrote global handoff manifest: ", handoff_path)
}

failed <- results[!vapply(results, function(x) isTRUE(x$success), logical(1))]
if (length(failed)) {
  stop(
    "Animal-level ProTigy preparation failed validation for: ",
    paste(vapply(failed, function(x) {
      paste0(x$dataset, " [", x$validation_status, "]")
    }, character(1)), collapse = "; "),
    call. = FALSE
  )
}

message("Animal-level ProTigy inputs completed for: ", paste(datasets, collapse = ", "))
