#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "module_contracts.R"))
source(repo_path("R", "schema_validation.R"))

setwd(repo_root())

fail <- character()
check <- function(expr, label) {
  ok <- tryCatch({
    force(expr)
    TRUE
  }, error = function(e) {
    fail <<- c(fail, paste0(label, ": ", conditionMessage(e)))
    FALSE
  })
  invisible(ok)
}

read_csv_required <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

skip_missing_output <- function(path, label) {
  strict <- tolower(Sys.getenv("PROTEOMICS_STRICT_OUTPUT_CONTRACTS", unset = "false")) %in% c("1", "true", "yes")
  exists <- if (grepl("[.]csv$", path, ignore.case = TRUE)) file.exists(path) else dir.exists(path)
  if (!exists && isTRUE(strict)) stop("Missing ", label, ": ", path, call. = FALSE)
  if (!exists) {
    message("SKIP output contract not present on this runner: ", label)
    return(TRUE)
  }
  FALSE
}

registry <- NULL
check({
  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
}, "pipeline.yml parses and validates")

if (!is.null(registry)) {
  steps <- pipeline_steps(registry, pipeline_stage_names(registry), dataset = "all", include_unsupported = TRUE)
  registry_entries <- do.call(rbind, lapply(names(registry$stages), function(stage) {
    data.frame(
      stage = stage,
      script = vapply(registry$stages[[stage]]$scripts, function(x) as.character(x$script), character(1)),
      stringsAsFactors = FALSE
    )
  }))
  registry_scripts_in_order <- unique(registry_entries$script)
  duplicate_script_exemptions <- "06_modules_WGCNA/03_score_module_activity.R"
  prefix_exemptions <- c(
    "03_qc_exploration/04_marker_rank_abundance_qc.r",
    "03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r",
    "03_qc_exploration/04d_compartment_marker_fidelity.r"
  )
  required_registry_fields <- c(
    "script", "stage", "scope", "supported_datasets", "consumes_required",
    "consumes_optional", "produces", "recomputes_core_state", "safe_downstream_rerun"
  )
  missing_step_cols <- setdiff(required_registry_fields, names(steps))
  if (length(missing_step_cols)) {
    fail <- c(fail, paste("pipeline_steps() missing expected column(s):", paste(missing_step_cols, collapse = ", ")))
  }

  check(validate_pipeline_scripts_exist(registry, fail = TRUE), "all registered active scripts exist")
  check(validate_run_order_against_registry(registry), "RUN_ORDER.md references only active or legacy scripts")

  duplicate_scripts <- names(which(table(registry_entries$script) > 1L))
  unexpected_duplicate_scripts <- setdiff(duplicate_scripts, duplicate_script_exemptions)
  if (length(unexpected_duplicate_scripts)) {
    fail <- c(fail, paste("Duplicate active filenames without exemption:", paste(unexpected_duplicate_scripts, collapse = ", ")))
  }

  prefix_entries <- unique(registry_entries)
  prefix_entries$folder <- dirname(prefix_entries$script)
  prefix_entries$prefix <- suppressWarnings(as.integer(sub("^([0-9]+).*$", "\\1", basename(prefix_entries$script))))
  prefix_entries <- prefix_entries[!prefix_entries$script %in% prefix_exemptions, , drop = FALSE]
  prefix_keys <- paste(prefix_entries$stage, prefix_entries$folder, sep = "::")
  for (key in unique(prefix_keys)) {
    idx <- which(prefix_keys == key)
    prefixes <- prefix_entries$prefix[idx]
    if (anyDuplicated(prefixes)) {
      fail <- c(fail, paste("Duplicate active numeric prefixes in", key, ":", paste(prefix_entries$script[idx], collapse = ", ")))
    }
    if (length(prefixes) > 1L && any(diff(prefixes) < 0L)) {
      fail <- c(fail, paste("Active numeric prefixes contradict pipeline.yml order in", key, ":", paste(prefix_entries$script[idx], collapse = ", ")))
    }
  }

  run_order_lines <- readLines(repo_path("RUN_ORDER.md"), warn = FALSE)
  direct_matches <- regexec("^Rscript[[:space:]]+([0-9]{2}_[^[:space:]]+\\.[Rr])(?:[[:space:]]|$)", run_order_lines)
  direct_scripts <- unique(vapply(regmatches(run_order_lines, direct_matches), function(x) if (length(x) >= 2L) x[[2]] else NA_character_, character(1)))
  direct_scripts <- direct_scripts[!is.na(direct_scripts)]
  if (!identical(direct_scripts, registry_scripts_in_order)) {
    fail <- c(fail, paste0(
      "RUN_ORDER.md direct commands do not match pipeline.yml order. Expected: ",
      paste(registry_scripts_in_order, collapse = ", "), " | Found: ", paste(direct_scripts, collapse = ", ")
    ))
  }

  audit_path <- repo_path("docs", "active_script_io_audit.tsv")
  if (file.exists(audit_path)) {
    documented <- utils::read.delim(audit_path, check.names = FALSE, stringsAsFactors = FALSE)
    documented_active <- unique(documented[[1]][grepl("^active", documented[[2]])])
    hidden_documented_active <- setdiff(documented_active, registry_scripts_in_order)
    if (length(hidden_documented_active)) {
      fail <- c(fail, paste("Scripts documented as active but absent from pipeline.yml:", paste(hidden_documented_active, collapse = ", ")))
    }
  }

  registered <- unique(steps$script)
  missing_active <- registered[!file.exists(repo_path(registered))]
  if (length(missing_active)) {
    fail <- c(fail, paste("Registered active scripts missing from working tree:", paste(missing_active, collapse = ", ")))
  }

  for (field in c("consumes_required", "consumes_optional", "produces")) {
    empty_required_field <- steps$script[is.na(steps[[field]])]
    if (length(empty_required_field)) {
      fail <- c(fail, paste("Registry field unexpectedly NA for", field, ":", paste(unique(empty_required_field), collapse = ", ")))
    }
  }

  global_stages <- names(registry$stages)[vapply(registry$stages, function(stage) {
    any(vapply(stage$scripts, function(step) identical(as.character(step$scope), "global") || "global" %in% unlist(step$datasets), logical(1)))
  }, logical(1))]
  for (stage in global_stages) {
    cmd <- c("run_dataset_pipeline.R", "--dataset", "all", "--stage", stage, "--dry-run")
    status <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)
    exit <- attr(status, "status")
    if (!is.null(exit) && exit != 0L) {
      fail <- c(fail, paste0("Global/all dry-run failed for stage ", stage, ":\n", paste(status, collapse = "\n")))
    }
  }
}

for (dataset in valid_datasets()) {
  check(validate_dataset(dataset) == dataset, paste("Dataset validation", dataset))
}

module_ann_path <- path_results("tables", "06_modules_WGCNA", "module_annotation", "microglia", "WGCNA_module_biological_annotation.csv")
if (!skip_missing_output(module_ann_path, "microglia WGCNA module annotation")) {
  check({
    module_ann <- read_csv_required(module_ann_path, "microglia WGCNA module annotation")
    validate_wgcna_module_annotation(module_ann)
  }, "WGCNA module annotation contract")
}

super_ann_path <- path_results("tables", "06_modules_WGCNA", "module_annotation", "microglia", "WGCNA_supermodule_biological_annotation.csv")
if (!skip_missing_output(super_ann_path, "microglia WGCNA supermodule annotation")) {
  check({
    super_ann <- read_csv_required(super_ann_path, "microglia WGCNA supermodule annotation")
    require_module_contract_columns(
      super_ann,
      c(
        "dataset", "SupermoduleID", "dominant_microenvironment_class", "dominant_module_labels",
        "Supermodule_LabelRationale", "Supermodule_CompositionLabel",
        "Supermodule_CompositionDisplayLabel", "Supermodule_CleanPlotLabel",
        "microenvironment_caution_label", "raw_GO_BP_terms", "cleaned_biological_label"
      ),
      "WGCNA supermodule biological annotation"
    )
  }, "WGCNA supermodule annotation contract")
}

module_interp_path <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", "microglia", "WGCNA_module_group_effects_interpretable.csv")
if (!skip_missing_output(module_interp_path, "microglia WGCNA module interpretable summary")) {
  check({
    module_interp <- read_csv_required(module_interp_path, "microglia WGCNA module interpretable summary")
    validate_wgcna_interpretable_summary(module_interp)
    require_module_contract_columns(
      module_interp,
      c(
        "ModulePlotLabel", "Supermodule_PlotLabel", "Supermodule_FullAnnotationLabel",
        "module_label", "supermodule_id", "Module_CleanPlotLabel",
        "Supermodule_CleanPlotLabel", "cleaned_biological_label",
        "microenvironment_caution_label",
        "targeted_signature_primary_driver", "targeted_signature_driver_class",
        "targeted_signature_driver_signature", "targeted_signature_driver_padj",
        "targeted_signature_driver_NES", "targeted_signature_driver_overlap_proteins"
      ),
      "WGCNA module interpretable summary"
    )
  }, "WGCNA module interpretable summary contract")
}

super_interp_path <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", "microglia", "WGCNA_supermodule_group_effects_interpretable.csv")
if (!skip_missing_output(super_interp_path, "microglia WGCNA supermodule interpretable summary")) {
  check({
    super_interp <- read_csv_required(super_interp_path, "microglia WGCNA supermodule interpretable summary")
    validate_wgcna_interpretable_summary(super_interp)
    require_module_contract_columns(
      super_interp,
      c(
        "Supermodule_PlotLabel", "Supermodule_FullAnnotationLabel",
        "Supermodule_CleanPlotLabel", "Supermodule_CompositionDisplayLabel",
        "Supermodule_CompositionLabel", "cleaned_biological_label",
        "microenvironment_caution_label",
        "Supermodule_DisplayShort", "dominant_microenvironment_class",
        "dominant_module_labels", "Supermodule_LabelRationale"
      ),
      "WGCNA supermodule interpretable summary"
    )
  }, "WGCNA supermodule interpretable summary contract")
}

score_pub_validation_path <- path_results("tables", "06_modules_WGCNA", "score_publication_summary", "microglia", "WGCNA_score_publication_validation.csv")
if (!skip_missing_output(score_pub_validation_path, "microglia WGCNA score-publication validation")) {
  check({
    score_pub_validation <- read_csv_required(score_pub_validation_path, "microglia WGCNA score-publication validation")
    require_module_contract_columns(
      score_pub_validation,
      c("validation_check", "validation_status", "validation_message"),
      "WGCNA score-publication validation"
    )
    if (any(score_pub_validation$validation_status != "ok")) {
      stop("Score-publication validation has non-ok rows.", call. = FALSE)
    }
  }, "WGCNA score-publication validation contract")
}

bundle_path <- path_results("tables", "10_biological_integration", "final_evidence_bundle", "global")
if (!skip_missing_output(bundle_path, "final biological evidence bundle")) {
  check({
    validate_final_evidence_bundle(bundle_path)
  }, "final biological evidence bundle contract")
}

claims_path <- path_results("tables", "biological_claims_table.csv")
if (!skip_missing_output(claims_path, "biological claims table")) {
  check({
    claims <- read_csv_required(claims_path, "biological claims table")
    validate_table_schema(claims, "biological_claims_table", strict = TRUE)
  }, "biological claims table schema")
}

if (length(fail)) {
  message("FAIL smoke_test_active_script_contracts")
  message(paste(fail, collapse = "\n"))
  quit(status = 1, save = "no")
}

message("PASS smoke_test_active_script_contracts")
