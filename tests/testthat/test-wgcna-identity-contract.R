testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_identity_contract_utils.R"))

identity_test_datasets <- c("neuron_soma", "neuron_neuropil", "microglia")

identity_output_file <- function(dataset, filename) {
  path_results(
    "tables", "06_modules_WGCNA", "identity_contract", dataset, filename
  )
}

read_identity_output <- function(dataset, filename) {
  path <- identity_output_file(dataset, filename)
  testthat::expect_true(
    file.exists(path),
    info = paste(
      "Generated identity-contract integration artifact is required:",
      wgcna_identity_relative_path(path)
    )
  )
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )
}

identity_absolute_path <- function(relative_path) {
  do.call(
    repo_path,
    as.list(strsplit(relative_path, "/", fixed = TRUE)[[1]])
  )
}

testthat::test_that("all canonical datasets are accepted and all/unsupported fail", {
  testthat::expect_identical(
    wgcna_identity_datasets(),
    identity_test_datasets
  )
  for (dataset in identity_test_datasets) {
    testthat::expect_identical(
      wgcna_identity_validate_dataset(dataset),
      dataset
    )
  }
  testthat::expect_error(
    wgcna_identity_validate_dataset("all"),
    "Unsupported WGCNA identity-contract dataset"
  )
  testthat::expect_error(
    wgcna_identity_validate_dataset("neuron"),
    "Unsupported WGCNA identity-contract dataset"
  )
})

testthat::test_that("hexadecimal module syntax normalization is exact and narrow", {
  observed <- wgcna_identity_normalize_hex_module_id(c(
    "ME#abcdef", "#ABCDEF", "WGCNA_#aBcDeF",
    "MEforestgreen", "WGCNA_forestgreen", "#ABCDE", "#ABCDEFG",
    "prefix#ABCDEF", NA_character_
  ))
  testthat::expect_identical(
    observed[1:3],
    rep("WGCNA_#ABCDEF", 3)
  )
  testthat::expect_true(all(is.na(observed[4:9])))
})

testthat::test_that("live source selection is publishable for all three datasets", {
  for (dataset in identity_test_datasets) {
    bundle <- wgcna_identity_build_contract_bundle(dataset)
    testthat::expect_true(bundle$publishable, info = dataset)
    testthat::expect_false(any(
      bundle$validation$status == "fail" &
        bundle$validation$severity == "error"
    ), info = dataset)
    expected_role <- if (dataset == "microglia") {
      "current_membership_mapping"
    } else {
      "selected_cluster_assignment"
    }
    selected_row <- bundle$source_hashes[
      bundle$source_hashes$validation_status ==
        "accepted_membership_authority",
      ,
      drop = FALSE
    ]
    testthat::expect_equal(nrow(selected_row), 1L, info = dataset)
    testthat::expect_identical(selected_row$source_role, expected_role)
  }
})

testthat::test_that("generated entity and membership keys are unique", {
  for (dataset in identity_test_datasets) {
    entities <- read_identity_output(
      dataset, "WGCNA_entity_identity_contract.csv"
    )
    membership <- read_identity_output(
      dataset, "WGCNA_module_supermodule_membership_contract.csv"
    )
    testthat::expect_false(
      anyDuplicated(entities[c("dataset", "level", "entity_id")]) > 0L,
      info = dataset
    )
    testthat::expect_false(
      anyDuplicated(membership[c("dataset", "module_membership_key")]) > 0L,
      info = dataset
    )
  }
})

testthat::test_that("every current module maps once to one current supermodule", {
  for (dataset in identity_test_datasets) {
    membership <- read_identity_output(
      dataset, "WGCNA_module_supermodule_membership_contract.csv"
    )
    testthat::expect_false(anyNA(membership$module_id), info = dataset)
    testthat::expect_false(anyNA(membership$supermodule_id), info = dataset)
    counts <- table(membership$module_id)
    testthat::expect_true(all(counts == 1L), info = dataset)
  }
})

testthat::test_that("state, definitions, and selected mapping module sets agree", {
  for (dataset in identity_test_datasets) {
    bundle <- wgcna_identity_build_contract_bundle(dataset)
    definition <- wgcna_identity_read_csv(
      bundle$paths$module_definitions
    )
    definition_ids <- sort(unique(as.character(definition$ModuleID)))
    membership_ids <- sort(unique(bundle$membership_contract$module_id))
    testthat::expect_identical(
      membership_ids,
      bundle$state_module_ids,
      info = dataset
    )
    testthat::expect_identical(
      definition_ids,
      bundle$state_module_ids,
      info = dataset
    )
  }
})

testthat::test_that("membership hashes are deterministic across row order", {
  membership <- data.frame(
    dataset = "neuron_soma",
    module_id = c("WGCNA_#ABCDEF", "WGCNA_#012345"),
    supermodule_id = c("SM02", "SM01"),
    stringsAsFactors = FALSE
  )
  hash_a <- wgcna_identity_membership_version(membership)
  hash_b <- wgcna_identity_membership_version(
    membership[2:1, , drop = FALSE]
  )
  testthat::expect_identical(hash_a, hash_b)
  testthat::expect_match(hash_a, "^[0-9a-f]{64}$")
})

testthat::test_that("same SM06 string with different members is incompatible", {
  contract_key <- wgcna_identity_supermodule_membership_key(
    "neuron_neuropil", c("WGCNA_#111111", "WGCNA_#222222")
  )
  downstream_key <- wgcna_identity_supermodule_membership_key(
    "neuron_neuropil", c("WGCNA_#111111", "WGCNA_#333333")
  )
  status <- wgcna_identity_classify_compatibility(
    level = "supermodule",
    downstream_entity_id = "SM06",
    downstream_normalized_id = "SM06",
    downstream_membership_key = downstream_key,
    contract_entity_id = "SM06",
    contract_membership_key = contract_key,
    label_conflict = FALSE
  )
  testthat::expect_identical(
    status,
    "conflicting_membership_same_entity_id"
  )
})

testthat::test_that("identical membership with different labels remains identity-compatible", {
  membership_key <- wgcna_identity_supermodule_membership_key(
    "neuron_neuropil", c("WGCNA_#111111", "WGCNA_#222222")
  )
  status <- wgcna_identity_classify_compatibility(
    level = "supermodule",
    downstream_entity_id = "SM06",
    downstream_normalized_id = "SM06",
    downstream_membership_key = membership_key,
    contract_entity_id = "SM06",
    contract_membership_key = membership_key,
    label_conflict = TRUE
  )
  testthat::expect_identical(status, "label_conflict_same_membership")
})

testthat::test_that("generated source hashes prove protected frozen states are unchanged", {
  for (dataset in identity_test_datasets) {
    hashes <- read_identity_output(
      dataset, "WGCNA_identity_source_hashes.csv"
    )
    state_row <- hashes[hashes$source_role == "frozen_state", , drop = FALSE]
    testthat::expect_equal(nrow(state_row), 1L, info = dataset)
    state_path <- identity_absolute_path(state_row$source_path)
    testthat::expect_identical(
      file_hash_sha256(state_path),
      state_row$sha256,
      info = dataset
    )
  }
})

testthat::test_that("identity execution did not alter any protected Stage 01 source", {
  protected_roles <- c(
    "frozen_state", "module_definitions", "current_membership_mapping",
    "current_supermodule_summary", "selected_cluster_assignment",
    "clustering_sensitivity", "run_manifest", "wgcna_run_manifest",
    "output_manifest", "parameter_audit", "analysis_parameters"
  )
  for (dataset in identity_test_datasets) {
    hashes <- read_identity_output(
      dataset, "WGCNA_identity_source_hashes.csv"
    )
    protected <- hashes[hashes$source_role %in% protected_roles, , drop = FALSE]
    testthat::expect_equal(
      sort(protected$source_role),
      sort(protected_roles),
      info = dataset
    )
    actual <- vapply(protected$source_path, function(relative_path) {
      file_hash_sha256(identity_absolute_path(relative_path))
    }, character(1))
    testthat::expect_identical(
      unname(actual),
      unname(protected$sha256),
      info = dataset
    )
  }
})

testthat::test_that("generated contracts match observed live cardinalities", {
  for (dataset in identity_test_datasets) {
    bundle <- wgcna_identity_build_contract_bundle(dataset)
    entities <- read_identity_output(
      dataset, "WGCNA_entity_identity_contract.csv"
    )
    membership <- read_identity_output(
      dataset, "WGCNA_module_supermodule_membership_contract.csv"
    )
    testthat::expect_equal(
      sum(entities$level == "module"),
      length(bundle$state_module_ids),
      info = dataset
    )
    testthat::expect_equal(
      sum(entities$level == "supermodule"),
      length(unique(bundle$selected$membership$supermodule_id)),
      info = dataset
    )
    testthat::expect_equal(
      nrow(membership),
      length(bundle$state_module_ids),
      info = dataset
    )
  }
})

testthat::test_that("required output schemas and validations are complete", {
  entity_columns <- c(
    "dataset", "level", "entity_id", "canonical_entity_id",
    "identity_role", "n_member_modules", "membership_version",
    "state_sha256", "module_definitions_sha256",
    "supermodule_mapping_sha256", "identity_source", "identity_status",
    "selected_cut_height", "selected_cut_height_source", "contract_version"
  )
  membership_columns <- c(
    "dataset", "module_id", "supermodule_id", "module_membership_key",
    "supermodule_membership_key", "membership_version", "state_sha256",
    "mapping_sha256", "membership_status", "contract_version"
  )
  audit_columns <- c(
    "dataset", "downstream_stage", "source_file", "level",
    "downstream_entity_id", "contract_entity_id",
    "downstream_membership_key", "contract_membership_key",
    "exact_id_match", "exact_membership_match", "identity_compatible",
    "label_conflict_detected", "duplicate_downstream_key",
    "compatibility_status", "exclusion_or_warning_reason",
    "downstream_entity_id_original", "downstream_entity_id_normalized",
    "downstream_member_modules_original",
    "downstream_member_modules_normalized"
  )
  for (dataset in identity_test_datasets) {
    entities <- read_identity_output(
      dataset, "WGCNA_entity_identity_contract.csv"
    )
    membership <- read_identity_output(
      dataset, "WGCNA_module_supermodule_membership_contract.csv"
    )
    audit <- read_identity_output(
      dataset, "WGCNA_downstream_identity_compatibility_audit.csv"
    )
    validation <- read_identity_output(
      dataset, "WGCNA_identity_validation.csv"
    )
    testthat::expect_true(all(entity_columns %in% names(entities)), info = dataset)
    testthat::expect_true(all(membership_columns %in% names(membership)), info = dataset)
    testthat::expect_true(all(audit_columns %in% names(audit)), info = dataset)
    testthat::expect_false(any(
      validation$status == "fail" & validation$severity == "error"
    ), info = dataset)
  }
})

testthat::test_that("fail-closed validation and status semantics are explicit", {
  validation <- data.frame(
    status = c("pass", "fail"),
    severity = c("error", "error"),
    stringsAsFactors = FALSE
  )
  testthat::expect_false(wgcna_identity_contract_publishable(validation))
  script <- paste(
    readLines(
      testthat::test_path(
        "..", "..", "06_modules_WGCNA", "00_wgcna_identity_contract.R"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  testthat::expect_match(script, "identity_contract_not_publishable", fixed = TRUE)
  testthat::expect_match(
    script,
    "entity and membership contracts were not published",
    fixed = TRUE
  )
})

testthat::test_that("dry-run preserves generated outputs and --dataset all fails", {
  dataset <- "neuron_soma"
  output_dir <- path_results(
    "tables", "06_modules_WGCNA", "identity_contract", dataset
  )
  files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
  before <- stats::setNames(vapply(files, file_hash_sha256, character(1)), files)
  repo_root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  script <- file.path(
    "06_modules_WGCNA", "00_wgcna_identity_contract.R"
  )
  rscript <- Sys.which("Rscript")
  testthat::expect_true(nzchar(rscript))
  dry_output <- withr::with_dir(
    repo_root,
    system2(
      rscript,
      c(script, "--dataset", dataset, "--dry-run"),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  testthat::expect_null(attr(dry_output, "status"))
  after <- stats::setNames(vapply(files, file_hash_sha256, character(1)), files)
  testthat::expect_identical(after, before)

  all_output <- suppressWarnings(withr::with_dir(
    repo_root,
    system2(
      rscript,
      c(script, "--dataset", "all", "--dry-run"),
      stdout = TRUE,
      stderr = TRUE
    )
  ))
  testthat::expect_true(!is.null(attr(all_output, "status")))
  testthat::expect_match(
    paste(all_output, collapse = "\n"),
    "Expected one of|Unsupported"
  )
})

testthat::test_that("pipeline registration remains read-only and precedes Stage 05", {
  pipeline <- readLines(
    testthat::test_path("..", "..", "pipeline.yml"),
    warn = FALSE
  )
  identity_line <- grep(
    'script: "06_modules_WGCNA/00_wgcna_identity_contract.R"',
    pipeline
  )
  stage05_line <- grep(
    'script: "06_modules_WGCNA/05_module_supermodule_group_effects.r"',
    pipeline
  )
  testthat::expect_length(identity_line, 1L)
  testthat::expect_length(stage05_line, 1L)
  testthat::expect_lt(identity_line, stage05_line)
  block <- paste(
    pipeline[identity_line:min(stage05_line - 1L, length(pipeline))],
    collapse = "\n"
  )
  testthat::expect_match(block, "recomputes_core_state: false", fixed = TRUE)
  testthat::expect_match(block, "safe_downstream_rerun: true", fixed = TRUE)
})
