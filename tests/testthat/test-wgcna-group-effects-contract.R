testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_identity_contract_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_group_effects_utils.R"))

group_test_datasets <- c("neuron_soma", "neuron_neuropil", "microglia")

group_state_path <- function(dataset) {
  path_processed(
    "06_modules_WGCNA", "01_WGCNA", dataset,
    "wgcna_final_model_state.rds"
  )
}

group_output_path <- function(dataset, filename) {
  path_results(
    "tables", "06_modules_WGCNA", "group_effects", dataset, filename
  )
}

read_group_output <- function(dataset, filename) {
  path <- group_output_path(dataset, filename)
  testthat::expect_true(
    file.exists(path),
    info = paste("Generated Phase 2 artifact is required:", path)
  )
  utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )
}

load_group_fixture <- function(dataset) {
  state_path <- group_state_path(dataset)
  contract <- wgcna_group_load_identity_contract(dataset, state_path)
  state <- readRDS(state_path)
  bridge <- wgcna_group_build_module_bridge(dataset, state, contract)
  list(
    state = state,
    contract = contract,
    bridge = bridge,
    endpoints = wgcna_group_construct_endpoints(
      dataset, state, contract, bridge
    )
  )
}

testthat::test_that("all three publishable identity contracts load and validate", {
  for (dataset in group_test_datasets) {
    contract <- wgcna_group_load_identity_contract(
      dataset, group_state_path(dataset)
    )
    testthat::expect_true(
      contract$status$status[[1]] %in%
        c("publishable", "identity_contract_publishable"),
      info = dataset
    )
    testthat::expect_identical(
      contract$frozen_state_sha256,
      file_hash_sha256(group_state_path(dataset)),
      info = dataset
    )
  }
})

testthat::test_that("missing or stale membership versions fail closed", {
  testthat::expect_error(
    wgcna_group_assert_membership_version(
      c("membership_v1", "membership_v2")
    ),
    "missing or inconsistent"
  )
  testthat::expect_error(
    wgcna_group_assert_membership_version(c(NA_character_, "")),
    "missing or inconsistent"
  )
  testthat::expect_identical(
    wgcna_group_assert_membership_version(rep("membership_v1", 3)),
    "membership_v1"
  )
})

testthat::test_that("frozen-state module bridges are complete and one-to-one", {
  for (dataset in group_test_datasets) {
    fixture <- load_group_fixture(dataset)
    expected <- sort(unique(fixture$contract$membership$module_id))
    testthat::expect_identical(
      sort(fixture$bridge$module_id), expected, info = dataset
    )
    testthat::expect_false(anyDuplicated(fixture$bridge$module_id) > 0L)
    testthat::expect_false(
      anyDuplicated(fixture$bridge$state_eigengene_col_raw) > 0L
    )
  }
})

testthat::test_that("endpoint cardinalities and memberships equal the contracts", {
  for (dataset in group_test_datasets) {
    fixture <- load_group_fixture(dataset)
    provenance <- fixture$endpoints$provenance
    entity <- fixture$contract$entity
    testthat::expect_equal(
      sum(provenance$level == "module"),
      sum(entity$level == "module"),
      info = dataset
    )
    testthat::expect_equal(
      sum(provenance$level == "supermodule"),
      sum(entity$level == "supermodule"),
      info = dataset
    )
    observed <- fixture$endpoints$composition |>
      dplyr::select("dataset", "supermodule_id", "member_modules") |>
      tidyr::separate_longer_delim("member_modules", delim = ";") |>
      dplyr::rename(module_id = "member_modules") |>
      dplyr::arrange(.data$dataset, .data$module_id, .data$supermodule_id) |>
      dplyr::select("dataset", "module_id", "supermodule_id")
    expected <- fixture$contract$membership |>
      dplyr::arrange(.data$dataset, .data$module_id, .data$supermodule_id) |>
      dplyr::select("dataset", "module_id", "supermodule_id")
    testthat::expect_identical(
      as.data.frame(observed), as.data.frame(expected), info = dataset
    )
  }
})

testthat::test_that("supermodule construction and orientation are deterministic", {
  for (dataset in group_test_datasets) {
    fixture <- load_group_fixture(dataset)
    second <- wgcna_group_construct_endpoints(
      dataset, fixture$state, fixture$contract, fixture$bridge
    )
    testthat::expect_equal(
      fixture$endpoints$supermodule_eigengenes,
      second$supermodule_eigengenes,
      tolerance = 0,
      info = dataset
    )
    super_provenance <- fixture$endpoints$provenance[
      fixture$endpoints$provenance$level == "supermodule", ,
      drop = FALSE
    ]
    testthat::expect_true(all(
      super_provenance$orientation_correlation_after > 0
    ))
    testthat::expect_true(all(
      super_provenance$pc1_sign_after_orientation == 1
    ))
  }
})

testthat::test_that("singleton supermodules exactly equal member eigengenes", {
  for (dataset in group_test_datasets) {
    fixture <- load_group_fixture(dataset)
    singletons <- fixture$endpoints$composition[
      fixture$endpoints$composition$n_member_modules == 1L, ,
      drop = FALSE
    ]
    for (i in seq_len(nrow(singletons))) {
      member <- singletons$member_eigengene_columns[[i]]
      super_col <- singletons$supermodule_eigengene[[i]]
      testthat::expect_equal(
        fixture$endpoints$supermodule_eigengenes[[super_col]],
        fixture$endpoints$module_eigengenes[[member]],
        tolerance = 0,
        info = paste(dataset, singletons$supermodule_id[[i]])
      )
    }
  }
})

testthat::test_that("sample joins are exact and AnimalID provenance is explicit", {
  for (dataset in group_test_datasets) {
    state <- readRDS(group_state_path(dataset))
    eigengenes <- extract_module_eigengenes(state)
    audit <- wgcna_group_sample_audit(dataset, state, eigengenes)
    testthat::expect_equal(nrow(audit), nrow(eigengenes), info = dataset)
    testthat::expect_false(anyDuplicated(audit$Sample) > 0L)
    testthat::expect_true(all(audit$metadata_match_status == "exact_one_to_one_match"))
    testthat::expect_true(all(audit$inclusion_status == "included"))
    testthat::expect_true(all(audit$AnimalID_mapping_status == "resolved"))
    testthat::expect_true(all(!is.na(audit$AnimalID_source)))
  }
})

testthat::test_that("ambiguous AnimalID fails and replicate structure selects model type", {
  meta <- data.frame(
    AnimalID = c("A1", "A2"),
    MouseID = c("A3", "A2"),
    stringsAsFactors = FALSE
  )
  animal <- wgcna_group_animal_provenance(
    meta, c("sample_without_parse_1", "sample_without_parse_2")
  )
  testthat::expect_identical(animal$AnimalID_mapping_status[[1]], "ambiguous")
  testthat::expect_error(
    wgcna_group_model_type_for_scope(c("A1", NA_character_)),
    "missing or ambiguous"
  )
  testthat::expect_identical(
    wgcna_group_model_type_for_scope(c("A1", "A1", "A2")),
    "lmerTest_lmer"
  )
  testthat::expect_identical(
    wgcna_group_model_type_for_scope(c("A1", "A2", "A3")),
    "lm"
  )
})

testthat::test_that("primary workflow packages are installed", {
  testthat::expect_true(all(vapply(
    wgcna_group_required_packages(),
    requireNamespace,
    logical(1),
    quietly = TRUE
  )))
  testthat::expect_silent(wgcna_group_require_primary_packages())
})

testthat::test_that("named contrasts are level-order invariant and absent groups stay absent", {
  first <- wgcna_group_predeclared_contrasts(c("CON", "RES", "SUS"))
  second <- wgcna_group_predeclared_contrasts(c("SUS", "CON", "RES"))
  semantic <- function(x) {
    lapply(x, function(weights) sort(weights))
  }
  testthat::expect_identical(semantic(first), semantic(second))
  testthat::expect_identical(
    names(first[["SUS - RES"]]),
    c("CON", "RES", "SUS")
  )
  absent <- wgcna_group_predeclared_contrasts(c("CON", "SUS"))
  testthat::expect_identical(names(absent), "SUS - CON")
})

testthat::test_that("generated canonical tables contain only valid finite tests", {
  for (dataset in group_test_datasets) {
    module <- read_group_output(dataset, "module_group_effects.csv")
    supermodule <- read_group_output(dataset, "supermodule_group_effects.csv")
    primary <- dplyr::bind_rows(module, supermodule)
    testthat::expect_gt(nrow(primary), 0L)
    testthat::expect_true(all(primary$claim_allowed_model), info = dataset)
    testthat::expect_true(
      all(primary$model_valid_for_inference), info = dataset
    )
    testthat::expect_true(all(
      primary$model_stability_status %in% c(
        "stable_mixed_model", "boundary_random_intercept_zero",
        "stable_animal_level_lm"
      )
    ), info = dataset)
    testthat::expect_false(any(primary$fallback_used), info = dataset)
    named <- primary$test_type == "named_contrast"
    omnibus <- primary$test_type == "interaction_omnibus"
    testthat::expect_true(
      all(is.finite(primary$estimate[named])), info = dataset
    )
    testthat::expect_true(all(is.finite(primary$SE[named])), info = dataset)
    testthat::expect_true(all(is.finite(primary$statistic)), info = dataset)
    testthat::expect_true(all(is.finite(primary$p_value)), info = dataset)
    testthat::expect_true(all(
      primary$identical_rows_verified[omnibus] %in% TRUE
    ), info = dataset)
    testthat::expect_true(all(
      primary$min_unique_animals_compared_group >= 3L
    ), info = dataset)
    testthat::expect_true(all(
      primary$manuscript_claim_ready == "not_assessed_stage05"
    ), info = dataset)
  }
})

testthat::test_that("Phase 2B endpoint keys and FDR families are exact", {
  for (dataset in group_test_datasets) {
    primary <- dplyr::bind_rows(
      read_group_output(dataset, "module_group_effects.csv"),
      read_group_output(dataset, "supermodule_group_effects.csv")
    )
    keys <- primary[c(
      "dataset", "level", "endpoint_id", "effect_scope",
      "spatial_unit", "contrast"
    )]
    testthat::expect_false(anyDuplicated(keys) > 0L, info = dataset)
    aliases <- primary$hypothesis_level == "compatibility_alias"
    independent <- !aliases
    primary_global <- independent &
      primary$analysis_tier == "primary_wgcna_global"
    contextual <- independent &
      primary$analysis_tier == "secondary_contextual_global"
    testthat::expect_true(all(
      primary$contrast[primary_global] == "SUS - RES"
    ), info = dataset)
    testthat::expect_true(all(is.finite(
      primary$FDR_primary_global[primary_global]
    )), info = dataset)
    testthat::expect_true(all(is.na(
      primary$FDR_secondary_global[primary_global]
    )), info = dataset)
    testthat::expect_true(all(
      primary$contrast[contextual] %in% c("RES - CON", "SUS - CON")
    ), info = dataset)
    testthat::expect_true(all(is.na(
      primary$FDR_primary_global[contextual]
    )), info = dataset)
    testthat::expect_true(all(is.finite(
      primary$FDR_secondary_global[contextual]
    )), info = dataset)
    testthat::expect_true(all(vapply(
      c(
        "FDR_primary_global", "FDR_secondary_global",
        "FDR_interaction_omnibus", "FDR_local_exploratory",
        "FDR_conservative_all_tests", "FDR_global"
      ),
      function(nm) all(is.na(primary[[nm]][aliases])),
      logical(1)
    )), info = dataset)
    conservative <- independent &
      is.finite(primary$FDR_conservative_all_tests)
    testthat::expect_equal(
      primary$FDR_global[conservative],
      primary$FDR_conservative_all_tests[conservative],
      tolerance = 0, info = dataset
    )
  }
})

testthat::test_that("invalid attempts remain outside primary inference", {
  for (dataset in group_test_datasets) {
    validation <- read_group_output(
      dataset, "WGCNA_group_effect_model_validation.csv"
    )
    primary <- dplyr::bind_rows(
      read_group_output(dataset, "module_group_effects.csv"),
      read_group_output(dataset, "supermodule_group_effects.csv")
    )
    testthat::expect_true(all(
      validation$n_contrasts_emitted >= 0L
    ), info = dataset)
    failed <- validation[!grepl("^success", validation$attempt_status), ,
                         drop = FALSE]
    testthat::expect_true(all(failed$n_contrasts_emitted == 0L), info = dataset)
    testthat::expect_false(any(
      primary$statistical_support_status %in%
        c("FDR_supported", "suggestive_FDR10") &
        !primary$claim_allowed_model
    ), info = dataset)
  }
})

testthat::test_that("published endpoint and sample cardinalities derive from contracts", {
  for (dataset in group_test_datasets) {
    contract <- wgcna_group_load_identity_contract(
      dataset, group_state_path(dataset)
    )
    provenance <- read_group_output(
      dataset, "WGCNA_group_effect_endpoint_provenance.csv"
    )
    sample_audit <- read_group_output(
      dataset, "WGCNA_group_effect_sample_inclusion_audit.csv"
    )
    state <- readRDS(group_state_path(dataset))
    testthat::expect_equal(nrow(provenance), nrow(contract$entity))
    testthat::expect_equal(
      nrow(sample_audit), nrow(extract_module_eigengenes(state))
    )
    testthat::expect_true(all(
      provenance$identity_contract_sha256 ==
        contract$identity_contract_sha256
    ))
    testthat::expect_true(all(
      provenance$frozen_state_sha256 == contract$frozen_state_sha256
    ))
  }
})

testthat::test_that("legacy outputs are preserved at their audited hashes", {
  for (dataset in group_test_datasets) {
    audit <- read_group_output(
      dataset, "WGCNA_group_effect_legacy_output_staleness_audit.csv"
    )
    testthat::expect_gt(nrow(audit), 0L)
    for (i in seq_len(nrow(audit))) {
      path <- do.call(
        repo_path,
        as.list(strsplit(audit$source_file[[i]], "/", fixed = TRUE)[[1]])
      )
      testthat::expect_true(file.exists(path), info = audit$source_file[[i]])
      testthat::expect_identical(
        file_hash_sha256(path), audit$sha256[[i]],
        info = audit$source_file[[i]]
      )
    }
  }
})

testthat::test_that("atomic publication rollback restores prior directories", {
  root <- tempfile("wgcna-group-atomic-")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  roles <- c("tables", "source_data", "logs")
  targets <- stats::setNames(file.path(root, paste0("target_", roles)), roles)
  staged <- stats::setNames(file.path(root, paste0("stage_", roles)), roles)
  for (role in roles) {
    dir.create(targets[[role]], recursive = TRUE)
    dir.create(staged[[role]], recursive = TRUE)
    writeLines("old", file.path(targets[[role]], "marker.txt"))
    writeLines("new", file.path(staged[[role]], "marker.txt"))
  }
  testthat::expect_error(
    wgcna_group_atomic_publish(
      staged, targets, run_id = "injected_failure", fail_after = 1L
    ),
    "Injected atomic-publication failure"
  )
  for (role in roles) {
    testthat::expect_identical(
      readLines(file.path(targets[[role]], "marker.txt")),
      "old"
    )
  }
})

testthat::test_that("Phase 2 script is isolated from forbidden output writers", {
  script <- readLines(
    testthat::test_path(
      "..", "..", "06_modules_WGCNA",
      "05_module_supermodule_group_effects.r"
    ),
    warn = FALSE
  )
  text <- paste(script, collapse = "\n")
  testthat::expect_false(grepl("module_marker_trait_correlations.csv", text))
  testthat::expect_false(grepl("ggsave|save_plot|\\.svg", text))
  testthat::expect_false(grepl(
    "module_to_supermodule_map_with_annotations.csv", text, fixed = TRUE
  ))
  testthat::expect_true(grepl(
    "WGCNA_group_effect_legacy_output_staleness_audit.csv",
    text,
    fixed = TRUE
  ))
})
