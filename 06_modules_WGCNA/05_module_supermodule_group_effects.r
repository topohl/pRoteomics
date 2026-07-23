#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/05_module_supermodule_group_effects.r
# Stage: modules_downstream
# Scope: dataset_specific
# Purpose: canonical identity-contract-based module and supermodule
#          group-effect modelling. Phase 2 quantitative outputs only.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_identity_contract_utils.R"))
source(repo_path("R", "wgcna_group_effects_utils.R"))

run <- wgcna_cli(default_dataset = "neuron_neuropil", allow_all = FALSE)
DATASET <- wgcna_identity_validate_dataset(run$dataset)
LEVEL <- match.arg(tolower(run$level), c("module", "supermodule", "both"))
STATE_PATH <- path_processed(
  "06_modules_WGCNA", "01_WGCNA", DATASET,
  "wgcna_final_model_state.rds"
)
CONTRACT_PATHS <- wgcna_group_contract_paths(DATASET)

if (!run$dry_run && !identical(LEVEL, "both")) {
  stop(
    "Canonical Stage 05 publication requires --level both so that the ",
    "dataset-wide multiple-testing family is invariant. Single-level runs ",
    "are permitted only with --dry-run.",
    call. = FALSE
  )
}

if (!file.exists(STATE_PATH)) {
  stop("Missing frozen WGCNA state: ", STATE_PATH, call. = FALSE)
}

contract <- wgcna_group_load_identity_contract(DATASET, STATE_PATH)
state <- readRDS(STATE_PATH)
bridge <- wgcna_group_build_module_bridge(DATASET, state, contract)
module_eigengenes <- extract_module_eigengenes(state)
sample_audit <- wgcna_group_sample_audit(
  DATASET, state, module_eigengenes
)

if (run$dry_run) {
  dry_run_line(
    "Script",
    "06_modules_WGCNA/05_module_supermodule_group_effects.r"
  )
  dry_run_line("Dataset", DATASET)
  dry_run_line("Level", LEVEL)
  dry_run_line("Identity status", contract$status$status[[1]], "PASS")
  dry_run_line(
    "Identity contract SHA256",
    contract$identity_contract_sha256,
    "PASS"
  )
  dry_run_line(
    "Membership version",
    contract$membership_version,
    "PASS"
  )
  dry_run_line(
    "Frozen state SHA256",
    contract$frozen_state_sha256,
    "PASS"
  )
  dry_run_line(
    "Canonical endpoints",
    paste0(
      sum(contract$entity$level == "module"), " modules; ",
      sum(contract$entity$level == "supermodule"), " supermodules"
    ),
    "PASS"
  )
  dry_run_line(
    "Module-to-eigengene bridge",
    paste0(nrow(bridge), " complete one-to-one rows"),
    "PASS"
  )
  dry_run_line(
    "Sample/AnimalID audit",
    paste0(
      nrow(sample_audit), " samples; ",
      length(unique(sample_audit$AnimalID)), " animals; ",
      paste(unique(sample_audit$AnimalID_source), collapse = ";")
    ),
    "PASS"
  )
  dry_run_line(
    "Effect scopes",
    paste(
      c(
        "within_spatial_unit", "spatial_adjusted_global",
        "stress_by_spatial_interaction"
      ),
      collapse = ";"
    )
  )
  dry_run_line(
    "Primary packages",
    paste(wgcna_group_required_packages(), collapse = ";"),
    if (all(vapply(
      wgcna_group_required_packages(),
      requireNamespace,
      logical(1),
      quietly = TRUE
    ))) "PASS" else "FAIL"
  )
  quit(status = 0, save = "no")
}

wgcna_group_require_primary_packages()
PATHS <- wgcna_downstream_paths("group_effects", DATASET)
run_id <- paste0(
  format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
  "_", Sys.getpid()
)
failed_root <- path_results(
  "logs", "06_modules_WGCNA", "group_effects_failed",
  DATASET, run_id
)

write_failure <- function(error, validation = NULL) {
  dir_create(failed_root)
  failure <- data.frame(
    dataset = DATASET,
    status = "canonical_publication_failed",
    error_message = conditionMessage(error),
    identity_contract_sha256 = contract$identity_contract_sha256,
    membership_version = contract$membership_version,
    frozen_state_sha256 = contract$frozen_state_sha256,
    run_id = run_id,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    failure,
    file.path(failed_root, "WGCNA_group_effect_failed_run.csv"),
    row.names = FALSE,
    na = ""
  )
  if (!is.null(validation)) {
    utils::write.csv(
      validation,
      file.path(failed_root, "WGCNA_group_effect_model_validation.csv"),
      row.names = FALSE,
      na = ""
    )
  }
  writeLines(
    c(
      paste0("dataset: ", DATASET),
      paste0("run_id: ", run_id),
      paste0("error: ", conditionMessage(error)),
      "canonical_outputs_changed: false"
    ),
    file.path(failed_root, "failure_report.txt")
  )
}

result <- tryCatch({
  bundle <- wgcna_group_build_bundle(DATASET, state, contract)
  legacy_audit <- wgcna_group_legacy_staleness_audit(DATASET, PATHS)
  outputs <- list(
    module_group_effects.csv = bundle$module_group_effects,
    supermodule_group_effects.csv = bundle$supermodule_group_effects,
    supermodule_composition.csv = bundle$supermodule_composition,
    WGCNA_group_effect_endpoint_provenance.csv =
      bundle$endpoint_provenance,
    WGCNA_group_effect_model_validation.csv =
      bundle$model_validation,
    WGCNA_group_effect_sample_inclusion_audit.csv =
      bundle$sample_inclusion_audit,
    supermodule_pca_input_audit.csv =
      bundle$supermodule_pca_input_audit,
    supermodule_pca_eigenvalues.csv =
      bundle$supermodule_pca_eigenvalues,
    supermodule_pca_member_loadings.csv =
      bundle$supermodule_pca_member_loadings,
    all_supermodule_eigengene_group_values.csv =
      bundle$all_supermodule_eigengene_group_values,
    all_supermodule_eigengene_spatial_group_values.csv =
      bundle$all_supermodule_eigengene_spatial_group_values,
    WGCNA_group_effect_legacy_output_staleness_audit.csv =
      legacy_audit
  )
  wgcna_group_validate_output_bundle(outputs, contract)

  staged <- wgcna_group_prepare_stage(PATHS, DATASET, run_id)
  for (filename in names(outputs)) {
    write_csv_safe2(
      outputs[[filename]],
      file.path(staged$directories$tables, filename)
    )
    write_csv_safe2(
      outputs[[filename]],
      file.path(staged$directories$source_data, filename)
    )
  }

  manifest_inputs <- c(
    list(frozen_state = STATE_PATH),
    stats::setNames(as.list(unname(unlist(CONTRACT_PATHS))), names(CONTRACT_PATHS))
  )
  write_run_manifest(
    file.path(staged$directories$logs, "run_manifest.yml"),
    inputs = manifest_inputs,
    outputs = list(
      tables = PATHS$tables,
      source_data = PATHS$source_data,
      logs = PATHS$logs
    ),
    parameters = list(
      dataset = DATASET,
      level = LEVEL,
      group_effects_contract_version =
        wgcna_group_effects_contract_version(),
      identity_contract_version = contract$contract_version,
      identity_contract_sha256 = contract$identity_contract_sha256,
      membership_version = contract$membership_version,
      frozen_state_sha256 = contract$frozen_state_sha256,
      effect_scopes = c(
        "within_spatial_unit", "spatial_adjusted_global",
        "stress_by_spatial_interaction"
      ),
      fdr_within_level = "BH over dataset + level primary rows",
      fdr_dataset = "BH over dataset module + supermodule primary rows",
      fallback_enabled = FALSE,
      atomic_publication = TRUE
    ),
    notes = paste(
      "Phase 2 quantitative publication only.",
      "The Phase 1 identity contract is the sole membership authority.",
      "Legacy figures, labels, selected interpretations, and auxiliary outputs",
      "were preserved and recorded in the staleness audit."
    )
  )
  writeLines(
    capture.output(utils::sessionInfo()),
    file.path(staged$directories$logs, "sessionInfo.txt")
  )

  for (filename in names(outputs)) {
    table_path <- file.path(staged$directories$tables, filename)
    source_path <- file.path(staged$directories$source_data, filename)
    if (!file.exists(table_path) || !file.exists(source_path) ||
        !identical(
          file_hash_sha256(table_path),
          file_hash_sha256(source_path)
        )) {
      stop("Staged table/source mirror validation failed for ", filename,
           call. = FALSE)
    }
  }

  wgcna_group_atomic_publish(
    staged$directories,
    PATHS[c("tables", "source_data", "logs")],
    run_id = run_id
  )
  if (dir.exists(staged$root)) {
    wgcna_group_safe_remove_dir(staged$root, dirname(staged$root))
  }
  list(bundle = bundle, legacy_audit = legacy_audit)
}, error = function(e) {
  if (exists("staged", inherits = TRUE) && dir.exists(staged$root)) {
    wgcna_group_safe_remove_dir(staged$root, dirname(staged$root))
  }
  validation <- if (exists("bundle", inherits = TRUE)) {
    bundle$model_validation
  } else {
    NULL
  }
  write_failure(e, validation)
  stop(
    conditionMessage(e),
    "\nExisting canonical Stage 05 outputs were left unchanged. ",
    "Failure audit: ", failed_root,
    call. = FALSE
  )
})

message(
  "Canonical WGCNA group effects published for ", DATASET,
  ": module rows=", nrow(result$bundle$module_group_effects),
  "; supermodule rows=", nrow(result$bundle$supermodule_group_effects),
  "; model attempts=", nrow(result$bundle$model_validation),
  "; legacy outputs preserved=", nrow(result$legacy_audit), "."
)
