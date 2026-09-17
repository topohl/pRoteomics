# ================================================================
# Script: 06_modules_WGCNA/00_wgcna_identity_contract.R
# Stage: modules_downstream
# Scope: dataset_specific
# Purpose: publish a read-only canonical identity contract for a frozen WGCNA
#          module system and its provenance-selected Stage 01 supermodules.
# ================================================================

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_identity_contract_utils.R"))

run <- wgcna_cli(default_dataset = "neuron_neuropil", allow_all = FALSE)
DATASET <- wgcna_identity_validate_dataset(run$dataset)
OUTPUTS <- module_paths(
  "06_modules_WGCNA",
  file.path("identity_contract", DATASET)
)

bundle <- wgcna_identity_build_contract_bundle(DATASET)

if (run$dry_run) {
  dry_run_line("Script", "06_modules_WGCNA/00_wgcna_identity_contract.R")
  dry_run_line("Dataset", DATASET)
  dry_run_line(
    "Selected membership authority",
    wgcna_identity_relative_path(bundle$selected$membership_path),
    if (bundle$publishable) "PASS" else "FAIL"
  )
  dry_run_line("Selected cut height", bundle$selected$selected_cut_height)
  dry_run_line("Current modules", length(bundle$state_module_ids))
  dry_run_line(
    "Current supermodules",
    length(unique(bundle$membership_contract$supermodule_id))
  )
  dry_run_line(
    "Membership version",
    unique(bundle$membership_contract$membership_version)
  )
  dry_run_line(
    "Required validation failures",
    sum(bundle$validation$status == "fail" & bundle$validation$severity == "error"),
    if (bundle$publishable) "PASS" else "FAIL"
  )
  dry_run_line("Output tables", OUTPUTS$tables)
  quit(status = if (bundle$publishable) 0 else 1, save = "no")
}

invisible(lapply(
  OUTPUTS[c("tables", "source_data", "reports", "logs")],
  dir_create
))

write_identity_csv <- function(x, path) {
  dir_create(dirname(path))
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(
    x,
    temporary,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8",
    quote = TRUE
  )
  if (file.exists(path)) {
    removed <- unlink(path)
    if (!identical(removed, 0L)) {
      stop("Could not replace identity-contract output: ", path, call. = FALSE)
    }
  }
  if (!file.rename(temporary, path)) {
    stop("Could not publish identity-contract output: ", path, call. = FALSE)
  }
  invisible(path)
}

write_mirrored_csv <- function(x, filename) {
  paths <- c(
    file.path(OUTPUTS$tables, filename),
    file.path(OUTPUTS$source_data, filename)
  )
  invisible(lapply(paths, function(path) write_identity_csv(x, path)))
  paths
}

status <- data.frame(
  dataset = DATASET,
  status = if (bundle$publishable) {
    "identity_contract_publishable"
  } else {
    "identity_contract_not_publishable"
  },
  reason = if (bundle$publishable) {
    "All required identity validations passed."
  } else {
    paste(
      bundle$validation$details[
        bundle$validation$status == "fail" &
          bundle$validation$severity == "error"
      ],
      collapse = " | "
    )
  },
  selected_membership_source =
    wgcna_identity_relative_path(bundle$selected$membership_path),
  selected_cut_height = bundle$selected$selected_cut_height,
  n_modules = length(bundle$state_module_ids),
  n_supermodules = length(unique(bundle$membership_contract$supermodule_id)),
  membership_version = if (nrow(bundle$membership_contract)) {
    unique(bundle$membership_contract$membership_version)
  } else {
    NA_character_
  },
  contract_version = wgcna_identity_contract_version(),
  stringsAsFactors = FALSE
)

diagnostic_outputs <- c(
  write_mirrored_csv(
    bundle$source_hashes,
    "WGCNA_identity_source_hashes.csv"
  ),
  write_mirrored_csv(
    bundle$validation,
    "WGCNA_identity_validation.csv"
  ),
  write_mirrored_csv(
    bundle$compatibility_audit,
    "WGCNA_downstream_identity_compatibility_audit.csv"
  ),
  write_mirrored_csv(
    status,
    "WGCNA_identity_contract_status.csv"
  )
)

contract_outputs <- character()
entity_contract_path <- file.path(
  OUTPUTS$tables, "WGCNA_entity_identity_contract.csv"
)
membership_contract_path <- file.path(
  OUTPUTS$tables, "WGCNA_module_supermodule_membership_contract.csv"
)
entity_source_path <- file.path(
  OUTPUTS$source_data, "WGCNA_entity_identity_contract.csv"
)
membership_source_path <- file.path(
  OUTPUTS$source_data, "WGCNA_module_supermodule_membership_contract.csv"
)

if (bundle$publishable) {
  contract_outputs <- c(
    write_mirrored_csv(
      bundle$entity_contract,
      "WGCNA_entity_identity_contract.csv"
    ),
    write_mirrored_csv(
      bundle$membership_contract,
      "WGCNA_module_supermodule_membership_contract.csv"
    )
  )
} else {
  stale_contracts <- c(
    entity_contract_path, membership_contract_path,
    entity_source_path, membership_source_path
  )
  invisible(lapply(
    stale_contracts[file.exists(stale_contracts)],
    function(path) {
      expected_parent <- normalizePath(
        if (grepl("/source_data/|\\\\source_data\\\\", path)) {
          OUTPUTS$source_data
        } else {
          OUTPUTS$tables
        },
        winslash = "/",
        mustWork = TRUE
      )
      actual_parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
      if (!identical(tolower(actual_parent), tolower(expected_parent))) {
        stop("Refusing to invalidate an identity contract outside its dataset output directory.")
      }
      unlink(path)
    }
  ))
}

validation_failures <- bundle$validation[
  bundle$validation$status == "fail" &
    bundle$validation$severity == "error",
  ,
  drop = FALSE
]
validation_warnings <- bundle$validation[
  bundle$validation$status == "warning",
  ,
  drop = FALSE
]
compatibility_counts <- if (nrow(bundle$compatibility_audit)) {
  sort(table(bundle$compatibility_audit$compatibility_status), decreasing = TRUE)
} else {
  integer()
}

report_lines <- c(
  paste0("# WGCNA Identity Contract Audit — ", DATASET),
  "",
  paste0("- Status: `", status$status, "`"),
  paste0(
    "- Selected membership source: `",
    status$selected_membership_source,
    "`"
  ),
  paste0("- Selected cut height: `", status$selected_cut_height, "`"),
  paste0("- Current modules: `", status$n_modules, "`"),
  paste0("- Current supermodules: `", status$n_supermodules, "`"),
  paste0("- Membership version: `", status$membership_version, "`"),
  "",
  "## Validation",
  "",
  paste0("- Blocking failures: `", nrow(validation_failures), "`"),
  paste0("- Warnings: `", nrow(validation_warnings), "`"),
  "",
  "## Compatibility classifications",
  "",
  if (length(compatibility_counts)) {
    paste0(
      "- `", names(compatibility_counts), "`: ",
      as.integer(compatibility_counts)
    )
  } else {
    "- No downstream identity-bearing rows were available."
  },
  "",
  "Biological labels were audited as metadata only and were never used to establish identity."
)
report_path <- file.path(
  OUTPUTS$reports,
  "WGCNA_identity_contract_audit.md"
)
writeLines(report_lines, report_path, useBytes = TRUE)

manifest_outputs <- list(
  tables = OUTPUTS$tables,
  source_data = OUTPUTS$source_data,
  audit_report = report_path,
  published_contracts = contract_outputs,
  diagnostic_outputs = diagnostic_outputs
)
write_run_manifest(
  file.path(OUTPUTS$logs, "run_manifest.yml"),
  inputs = bundle$paths,
  outputs = manifest_outputs,
  parameters = list(
    dataset = DATASET,
    contract_version = wgcna_identity_contract_version(),
    selected_membership_source =
      wgcna_identity_relative_path(bundle$selected$membership_path),
    selected_cut_height = bundle$selected$selected_cut_height,
    membership_version = status$membership_version,
    publishable = bundle$publishable,
    recomputes_core_state = FALSE
  ),
  notes = paste(
    "Read-only WGCNA identity publication.",
    "Stage 05-13 outputs are compatibility-audit inputs only and never membership authority.",
    "Labels do not determine identity."
  )
)

if (!bundle$publishable) {
  stop(
    "WGCNA identity contract is not publishable for ", DATASET,
    ". Diagnostic hashes, validation, compatibility audit, and status were written; ",
    "entity and membership contracts were not published.",
    call. = FALSE
  )
}

message(
  "Published WGCNA identity contract for ", DATASET,
  ": ", status$n_modules, " modules, ", status$n_supermodules,
  " supermodules; membership_version=", status$membership_version
)
