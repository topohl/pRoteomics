testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_identity_contract_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_group_effects_utils.R"))

valid_fixture <- function(...) {
  defaults <- list(
    model_type = "lmerTest_lmer",
    fixed_effect_full_rank = TRUE,
    optimizer_success = TRUE,
    nonboundary_convergence_success = TRUE,
    finite_model_quantities = TRUE,
    contrast_estimable = TRUE,
    finite_contrast = TRUE,
    random_effect_structure = "single_animal_intercept",
    random_intercept_variance = 0.2,
    residual_variance = 1,
    is_singular_lme4 = FALSE,
    sample_contract_valid = TRUE
  )
  modifications <- list(...)
  do.call(
    wgcna_group_classify_model_diagnostics,
    utils::modifyList(defaults, modifications)
  )
}

testthat::test_that("boundary single-intercept diagnostics remain valid", {
  observed <- valid_fixture(
    random_intercept_variance = 0,
    is_singular_lme4 = TRUE
  )
  testthat::expect_true(observed$model_valid_for_inference)
  testthat::expect_identical(
    observed$model_stability_status,
    "boundary_random_intercept_zero"
  )
  testthat::expect_false(observed$primary_model_stable)
  testthat::expect_true(observed$claim_allowed_model)
  testthat::expect_match(observed$boundary_warning, "prespecified boundary")
})

testthat::test_that("variance-ratio boundary classification is scale independent", {
  tolerance <- wgcna_group_boundary_variance_ratio_tolerance()
  fixtures <- list(
    c(random = 1e-6, residual = 1e-2),
    c(random = 1, residual = 1e4),
    c(random = 2e-4, residual = 1e6),
    c(random = 0.01, residual = 1e4)
  )
  for (fixture in fixtures) {
    observed <- valid_fixture(
      random_intercept_variance = fixture[["random"]],
      residual_variance = fixture[["residual"]]
    )
    testthat::expect_lte(observed$variance_ratio, tolerance)
    testthat::expect_identical(
      observed$model_stability_status,
      "boundary_random_intercept_zero"
    )
  }
  stable <- list(
    valid_fixture(random_intercept_variance = 0.2, residual_variance = 1),
    valid_fixture(random_intercept_variance = 5e-5, residual_variance = 1e-5)
  )
  testthat::expect_true(all(vapply(
    stable,
    function(x) identical(x$model_stability_status, "stable_mixed_model"),
    logical(1)
  )))
  base <- valid_fixture(
    random_intercept_variance = 0.2, residual_variance = 1
  )
  scaled <- valid_fixture(
    random_intercept_variance = 2000, residual_variance = 10000
  )
  testthat::expect_identical(
    base$model_stability_status, scaled$model_stability_status
  )
  testthat::expect_equal(base$variance_ratio, scaled$variance_ratio)
})

testthat::test_that("invalid diagnostic fixtures fail closed", {
  rank_deficient <- valid_fixture(fixed_effect_full_rank = FALSE)
  nonconverged <- valid_fixture(
    optimizer_success = FALSE,
    nonboundary_convergence_success = FALSE
  )
  complex_singular <- valid_fixture(
    random_effect_structure = "complex_random_effects",
    is_singular_lme4 = TRUE
  )
  for (observed in list(
    rank_deficient, nonconverged, complex_singular
  )) {
    testthat::expect_false(observed$model_valid_for_inference)
    testthat::expect_identical(observed$model_stability_status, "invalid")
    testthat::expect_false(observed$claim_allowed_model)
  }
  testthat::expect_match(
    rank_deficient$failure_reason, "rank_deficient_fixed_effects"
  )
  testthat::expect_match(
    nonconverged$failure_reason,
    "optimizer_failure|nonboundary_convergence_failure"
  )
  testthat::expect_match(
    complex_singular$failure_reason, "complex_singular_random_effects"
  )
  invalid_variances <- list(
    valid_fixture(random_intercept_variance = NA_real_),
    valid_fixture(random_intercept_variance = -1),
    valid_fixture(residual_variance = 0),
    valid_fixture(residual_variance = Inf)
  )
  testthat::expect_true(all(vapply(
    invalid_variances,
    function(x) !isTRUE(x$model_valid_for_inference),
    logical(1)
  )))
})

testthat::test_that("model classification cannot depend on p or direction", {
  classifier_formals <- names(formals(wgcna_group_classify_model_diagnostics))
  testthat::expect_false(any(
    c(
      "p_value", "estimate", "effect_direction", "direction",
      "FDR", "FDR_primary_global"
    ) %in%
      classifier_formals
  ))
  first <- valid_fixture(random_intercept_variance = 0)
  second <- valid_fixture(random_intercept_variance = 0)
  testthat::expect_identical(first, second)
})

testthat::test_that("lme4 disagreement is audited without replacing ratio rule", {
  ratio_boundary_only <- valid_fixture(
    random_intercept_variance = 1e-6,
    residual_variance = 1e-2,
    is_singular_lme4 = FALSE
  )
  testthat::expect_identical(
    ratio_boundary_only$singularity_diagnostic_status,
    "variance_ratio_boundary_only"
  )
  testthat::expect_true(ratio_boundary_only$model_valid_for_inference)
  testthat::expect_identical(
    ratio_boundary_only$model_stability_status,
    "boundary_random_intercept_zero"
  )
  testthat::expect_false(ratio_boundary_only$primary_model_stable)
  testthat::expect_true(ratio_boundary_only$diagnostic_review_required)

  lme4_only <- valid_fixture(
    random_intercept_variance = 0.2,
    residual_variance = 1,
    is_singular_lme4 = TRUE
  )
  testthat::expect_identical(
    lme4_only$singularity_diagnostic_status, "lme4_singular_only"
  )
  testthat::expect_identical(
    lme4_only$model_stability_status, "stable_mixed_model"
  )
  testthat::expect_true(lme4_only$model_valid_for_inference)
  testthat::expect_false(lme4_only$primary_model_stable)
  testthat::expect_true(lme4_only$diagnostic_review_required)
})

testthat::test_that("installed lme4 boundary integration is retained when reproducible", {
  testthat::skip_if_not_installed("lmerTest")
  dat <- expand.grid(
    AnimalID = factor(sprintf("A%02d", 1:18)),
    SpatialUnit = factor(c("s1", "s2")),
    KEEP.OUT.ATTRS = FALSE
  )
  dat$StressGroup <- factor(rep(
    rep(c("CON", "RES", "SUS"), each = 6), each = 2
  ))
  dat$eigengene <- as.numeric(dat$StressGroup) +
    ifelse(dat$SpatialUnit == "s2", 0.25, -0.25) +
    rep(c(-0.1, 0.1), 18)
  fit <- suppressWarnings(lmerTest::lmer(
    eigengene ~ StressGroup + SpatialUnit + (1 | AnimalID),
    data = dat, REML = FALSE
  ))
  mm <- stats::model.matrix(
    eigengene ~ StressGroup + SpatialUnit, data = dat
  )
  diagnostics <- wgcna_group_fit_diagnostics(
    fit, "lmerTest_lmer", qr(mm)$rank, ncol(mm)
  )
  if (diagnostics$model_stability_status !=
      "boundary_random_intercept_zero") {
    testthat::skip(
      "Installed lme4 did not reproduce the optional zero-variance integration fixture."
    )
  }
  testthat::expect_true(diagnostics$model_valid_for_inference)
  testthat::expect_false(diagnostics$primary_model_stable)
  testthat::expect_true(diagnostics$claim_allowed_model)
})

make_fdr_row <- function(endpoint, level, hypothesis_level, tier, contrast,
                         scope, spatial, p_value,
                         test_type = "named_contrast") {
  row <- data.frame(
    dataset = "microglia", level = level, endpoint_id = endpoint,
    module_id = if (level == "module") endpoint else NA_character_,
    supermodule_id = if (level == "supermodule") endpoint else NA_character_,
    hypothesis_level = hypothesis_level,
    independent_hypothesis = TRUE,
    model_valid_for_inference = TRUE,
    claim_allowed_model = TRUE, fallback_used = FALSE,
    analysis_tier = tier, contrast = contrast, effect_scope = scope,
    spatial_unit = spatial, test_type = test_type, p_value = p_value,
    manuscript_claim_ready = "not_assessed_stage05",
    stringsAsFactors = FALSE
  )
  wgcna_group_bind_schema(list(row), wgcna_group_primary_schema())
}

testthat::test_that("primary and contextual global BH families are disjoint", {
  rows <- dplyr::bind_rows(
    make_fdr_row(
      "m1", "module", "module", "primary_wgcna_global",
      "SUS - RES", "spatial_adjusted_global",
      "global_spatial_adjusted", 0.01
    ),
    make_fdr_row(
      "m2", "module", "module", "primary_wgcna_global",
      "SUS - RES", "spatial_adjusted_global",
      "global_spatial_adjusted", 0.04
    ),
    make_fdr_row(
      "m1", "module", "module", "secondary_contextual_global",
      "RES - CON", "spatial_adjusted_global",
      "global_spatial_adjusted", 0.02
    ),
    make_fdr_row(
      "m2", "module", "module", "secondary_contextual_global",
      "RES - CON", "spatial_adjusted_global",
      "global_spatial_adjusted", 0.20
    ),
    make_fdr_row(
      "m1", "module", "module", "secondary_contextual_global",
      "SUS - CON", "spatial_adjusted_global",
      "global_spatial_adjusted", 0.03
    ),
    make_fdr_row(
      "m2", "module", "module", "secondary_contextual_global",
      "SUS - CON", "spatial_adjusted_global",
      "global_spatial_adjusted", 0.40
    )
  )
  adjusted <- wgcna_group_apply_fdr(
    rows, wgcna_group_primary_schema()
  )$module
  primary <- adjusted$analysis_tier == "primary_wgcna_global"
  secondary <- adjusted$analysis_tier == "secondary_contextual_global"
  testthat::expect_true(all(is.finite(
    adjusted$FDR_primary_global[primary]
  )))
  testthat::expect_true(all(is.na(
    adjusted$FDR_secondary_global[primary]
  )))
  testthat::expect_true(all(is.na(
    adjusted$FDR_primary_global[secondary]
  )))
  testthat::expect_true(all(is.finite(
    adjusted$FDR_secondary_global[secondary]
  )))
  testthat::expect_setequal(
    adjusted$contrast[secondary], c("RES - CON", "SUS - CON")
  )
  testthat::expect_identical(
    adjusted$FDR_global,
    adjusted$FDR_conservative_all_tests
  )
})

testthat::test_that("singleton aliases inherit estimates but no FDR", {
  module_rows <- dplyr::bind_rows(
    make_fdr_row(
      "m1", "module", "module", "primary_wgcna_global",
      "SUS - RES", "spatial_adjusted_global",
      "global_spatial_adjusted", 0.01
    ),
    make_fdr_row(
      "m1", "module", "module", "secondary_contextual_global",
      "RES - CON", "spatial_adjusted_global",
      "global_spatial_adjusted", 0.02
    )
  )
  module_rows$estimate <- c(0.5, 0.2)
  module_rows$SE <- c(0.1, 0.1)
  module_rows$CI_low <- module_rows$estimate - 0.196
  module_rows$CI_high <- module_rows$estimate + 0.196
  module_rows$model_formula <- "eigengene ~ StressGroup + SpatialUnit + (1 | AnimalID)"
  module_rows$formula_used <- module_rows$model_formula
  module_rows$model_stability_status <- "stable_mixed_model"
  module_rows$primary_model_stable <- TRUE
  module_rows$random_intercept_variance <- 0.2
  module_rows$residual_variance <- 1
  module_rows$variance_ratio <- 0.2
  module_rows$ICC <- 1 / 6
  module_rows$is_singular_lme4 <- FALSE
  module_rows$boundary_by_variance_ratio <- FALSE
  module_rows$boundary_variance_ratio_tolerance <- 1e-4
  module_rows$lme4_singularity_tolerance <- 1e-4
  module_rows$singularity_diagnostic_status <- "concordant_nonboundary"
  module_rows$diagnostic_review_required <- FALSE
  module_rows$singularity_class <- "stable_single_animal_intercept"
  module_rows$boundary_warning <- NA_character_
  module_rows$endpoint_construction_method <- "canonical_module_eigengene"
  module_rows$endpoint_provenance_status <- "validated"
  module_rows$membership_version <- "membership"
  module_rows$identity_contract_version <- "identity_v1"
  module_rows$identity_contract_sha256 <- "identity"
  module_rows$frozen_state_sha256 <- "state"
  module_rows <- wgcna_group_apply_fdr(
    module_rows, wgcna_group_primary_schema()
  )$module
  composition <- data.frame(
    supermodule_id = "SM01", n_member_modules = 1L,
    member_modules = "m1", stringsAsFactors = FALSE
  )
  aliases <- wgcna_group_make_singleton_alias_rows(
    module_rows, composition
  )
  testthat::expect_equal(aliases$estimate, module_rows$estimate)
  testthat::expect_equal(aliases$SE, module_rows$SE)
  testthat::expect_equal(aliases$CI_low, module_rows$CI_low)
  testthat::expect_equal(aliases$CI_high, module_rows$CI_high)
  testthat::expect_equal(aliases$p_value, module_rows$p_value)
  inherited_fields <- c(
    "model_formula", "formula_used",
    "random_intercept_variance", "residual_variance", "variance_ratio", "ICC",
    "is_singular_lme4", "boundary_by_variance_ratio",
    "boundary_variance_ratio_tolerance", "lme4_singularity_tolerance",
    "singularity_diagnostic_status", "diagnostic_review_required",
    "singularity_class", "boundary_warning",
    "model_valid_for_inference", "model_stability_status",
    "primary_model_stable", "claim_allowed_model",
    "endpoint_construction_method", "endpoint_provenance_status",
    "membership_version", "identity_contract_version",
    "identity_contract_sha256", "frozen_state_sha256"
  )
  for (nm in inherited_fields) {
    testthat::expect_identical(aliases[[nm]], module_rows[[nm]], info = nm)
  }
  testthat::expect_true(all(
    aliases$canonical_claim_entity_id == "m1"
  ))
  testthat::expect_true(all(aliases$support_source_entity_id == "m1"))
  testthat::expect_true(all(aliases$display_allowed))
  testthat::expect_true(all(!aliases$independent_hypothesis))
  testthat::expect_true(all(!aliases$enters_primary_fdr))
  testthat::expect_true(all(
    aliases$statistical_support_status ==
      "inherited_from_canonical_entity"
  ))
  testthat::expect_true(all(vapply(
    c(
      "FDR_primary_global", "FDR_secondary_global",
      "FDR_interaction_omnibus", "FDR_local_exploratory",
      "FDR_conservative_all_tests", "FDR_global",
      "FDR_within_dataset_level", "FDR_dataset_all_levels"
    ),
    function(nm) all(is.na(aliases[[nm]])),
    logical(1)
  )))
})

aggregation_fixture <- function() {
  data.frame(
    Sample = c("A1_L_1", "A1_L_2", "A1_R_1", "A2_L_1"),
    AnimalID = c("A1", "A1", "A1", "A2"),
    AnimalID_source = "synthetic_fixture",
    StressGroup = c("CON", "CON", "CON", "RES"),
    canonical_spatial_unit = "s1",
    SpatialUnitType = "region",
    Region = "s1",
    Layer = NA_character_,
    Hemisphere = c("L", "L", "R", "L"),
    endpoint = c(0, 2, 10, 4),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("hemisphere aggregation is equal-weight and content-addressed", {
  dat <- aggregation_fixture()
  observed <- wgcna_group_aggregate_endpoint(
    dat, "microglia", "module", "m1", "endpoint"
  )
  animal <- observed$animal_spatial_values
  a1 <- animal[animal$AnimalID == "A1", , drop = FALSE]
  a2 <- animal[animal$AnimalID == "A2", , drop = FALSE]
  testthat::expect_equal(a1$eigengene, 5.5)
  testthat::expect_equal(a1$n_hemispheres_observed, 2L)
  testthat::expect_identical(
    a1$hemisphere_completeness_status, "complete_bilateral_pair"
  )
  testthat::expect_equal(a2$eigengene, 4)
  testthat::expect_equal(a2$n_hemispheres_observed, 1L)
  testthat::expect_identical(
    a2$hemisphere_completeness_status,
    "one_sided_observed_no_imputation"
  )
  testthat::expect_true(all(grepl(
    "^[0-9a-f]{64}$", animal$aggregated_row_sha256
  )))

  reordered <- wgcna_group_aggregate_endpoint(
    dat[c(4, 2, 1, 3), ], "microglia", "module", "m1", "endpoint"
  )$animal_spatial_values
  testthat::expect_identical(
    animal$aggregated_row_sha256, reordered$aggregated_row_sha256
  )
  changed <- dat
  changed$endpoint[changed$Sample == "A1_L_2"] <- 3
  changed_hash <- wgcna_group_aggregate_endpoint(
    changed, "microglia", "module", "m1", "endpoint"
  )$animal_spatial_values$aggregated_row_sha256
  testthat::expect_false(identical(
    animal$aggregated_row_sha256, changed_hash
  ))
})

testthat::test_that("unequal technical counts cannot weight a hemisphere", {
  dat <- aggregation_fixture()
  extra <- dat[dat$Sample == "A1_L_1", , drop = FALSE]
  extra$Sample <- "A1_L_3"
  extra$endpoint <- 1
  dat <- dplyr::bind_rows(dat, extra)
  observed <- wgcna_group_aggregate_endpoint(
    dat, "microglia", "module", "m1", "endpoint"
  )
  a1 <- observed$animal_spatial_values[
    observed$animal_spatial_values$AnimalID == "A1", , drop = FALSE
  ]
  testthat::expect_equal(a1$eigengene, 5.5)
  testthat::expect_identical(
    a1$duplicate_source_row_status,
    "technical_duplicates_collapsed_within_hemisphere"
  )
})

testthat::test_that("ambiguous aggregation inputs fail closed", {
  dat <- aggregation_fixture()
  ambiguous <- dplyr::bind_rows(dat, dat[1, , drop = FALSE])
  ambiguous$endpoint[nrow(ambiguous)] <- 99
  testthat::expect_error(
    wgcna_group_aggregate_endpoint(
      ambiguous, "microglia", "module", "m1", "endpoint"
    ),
    "Duplicate source-sample ambiguity"
  )
  unresolved <- dat
  unresolved$Hemisphere[[1]] <- NA_character_
  testthat::expect_error(
    wgcna_group_aggregate_endpoint(
      unresolved, "microglia", "module", "m1", "endpoint"
    ),
    "unresolved hemisphere"
  )
})

testthat::test_that("interaction omnibus is a nested ML comparison on identical rows", {
  testthat::skip_if_not_installed("lmerTest")
  dat <- expand.grid(
    AnimalID = sprintf("A%02d", 1:18),
    SpatialUnit = c("s1", "s2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  dat$StressGroup <- rep(
    rep(c("CON", "RES", "SUS"), each = 6),
    times = 2
  )
  dat$eigengene <- as.numeric(factor(dat$StressGroup)) +
    ifelse(dat$SpatialUnit == "s2", 0.3, -0.3) +
    ifelse(dat$StressGroup == "SUS" & dat$SpatialUnit == "s2", 0.4, 0) +
    rep(seq(-0.2, 0.2, length.out = 18), times = 2) +
    sin(seq_len(nrow(dat))) / 50
  dat$SpatialUnitType <- "region"
  dat$AnimalID_source <- "synthetic_fixture"
  endpoint <- data.frame(
    endpoint_id = "m1",
    endpoint_construction_method = "synthetic",
    endpoint_provenance_status = "validated",
    stringsAsFactors = FALSE
  )
  contract <- list(
    membership_version = "membership_test",
    contract_version = "identity_test",
    identity_contract_sha256 = "identity",
    frozen_state_sha256 = "state"
  )
  observed <- wgcna_group_fit_interaction_omnibus(
    dat, "microglia", "module", endpoint, contract
  )
  testthat::expect_equal(nrow(observed$primary), 1L)
  testthat::expect_identical(
    observed$primary$test_type, "interaction_omnibus"
  )
  testthat::expect_true(observed$primary$identical_rows_verified)
  testthat::expect_true(is.finite(
    observed$primary$likelihood_ratio_statistic
  ))
  testthat::expect_true(is.finite(observed$primary$likelihood_ratio_df))
  testthat::expect_true(is.finite(
    observed$primary$likelihood_ratio_p_value
  ))
  required_reduced_full <- c(
    "formula", "row_hash", "fixed_effect_rank", "fixed_effect_columns",
    "optimizer_code", "optimizer_messages", "convergence_status",
    "convergence_warnings", "random_effect_structure",
    "random_intercept_variance", "residual_variance", "variance_ratio",
    "ICC", "is_singular_lme4", "boundary_by_variance_ratio",
    "singularity_class", "model_stability_status",
    "diagnostic_review_required"
  )
  testthat::expect_true(all(paste0(
    "reduced_", required_reduced_full
  ) %in% names(observed$primary)))
  testthat::expect_true(all(paste0(
    "full_", required_reduced_full
  ) %in% names(observed$primary)))
  testthat::expect_false(is.na(
    observed$primary$reduced_model_stability_status
  ))
  testthat::expect_false(is.na(
    observed$primary$full_model_stability_status
  ))
})

testthat::test_that("consumer scan is deterministic and excludes nonactive roots", {
  first <- wgcna_group_scan_downstream_consumers()
  second <- wgcna_group_scan_downstream_consumers()
  testthat::expect_identical(first, second)
  testthat::expect_true(all(first$migration_required))
  testthat::expect_identical(
    first$blocking_for_execution, first$runtime_execution_consumer
  )
  testthat::expect_identical(
    first$phase3_runtime_migration_required,
    first$runtime_execution_consumer
  )
  testthat::expect_true(all(first$consumer_kind %in% c(
    "runtime_stage", "runtime_helper", "pipeline_registry",
    "schema", "test", "other_source"
  )))
  testthat::expect_true(all(is.finite(first$line_number)))
  testthat::expect_true(all(
    first$enforcement_status == "advisory_not_runtime_enforced"
  ))
  testthat::expect_false(any(grepl(
    "^(\\.git|99_deprecated|results|data|tmp|docs)/",
    first$consumer_script
  )))
  testthat::expect_true(any(grepl(
    "^tests/testthat/", first$consumer_script
  )))
  testthat::expect_false(any(first$consumer_script %in% c(
    "R/wgcna_group_effects_utils.R",
    "06_modules_WGCNA/05_module_supermodule_group_effects.r",
    "tests/testthat/test-wgcna-group-effects-phase2b.R",
    "tests/testthat/test-wgcna-group-effects-contract.R",
    "tests/testthat/test-schema-validation.R"
  )))
  testthat::expect_setequal(
    unique(first$consumed_column),
    wgcna_group_consumer_tokens()$consumed_column
  )
  testthat::expect_true(all(nzchar(first$scan_included_path_rules)))
  testthat::expect_true(all(nzchar(first$scan_excluded_path_rules)))
})

testthat::test_that("consumer audit exhaustively matches deterministic local scan", {
  root <- normalizePath(repo_root(), winslash = "/", mustWork = TRUE)
  top_dirs <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  active_stage_dirs <- top_dirs[
    grepl("^(0[0-9]|[1-8][0-9]|9[0-8])_", top_dirs)
  ]
  roots <- c(
    file.path(root, "R"),
    file.path(root, "inst", "schemas"),
    file.path(root, "tests", "testthat"),
    file.path(root, active_stage_dirs)
  )
  files <- unlist(lapply(roots[dir.exists(roots)], function(path) {
    list.files(
      path, pattern = "\\.(R|r|ya?ml)$", recursive = TRUE,
      full.names = TRUE, include.dirs = FALSE, all.files = TRUE, no.. = TRUE
    )
  }), use.names = FALSE)
  files <- unique(c(
    files,
    list.files(
      root, pattern = "\\.(R|r|ya?ml)$", recursive = FALSE,
      full.names = TRUE, include.dirs = FALSE, all.files = TRUE, no.. = TRUE
    )
  ))
  files <- normalizePath(files, winslash = "/", mustWork = FALSE)
  relative <- gsub(
    "\\\\", "/", substring(files, nchar(root) + 2L)
  )
  explicit_exclusions <- c(
    "R/wgcna_group_effects_utils.R",
    "06_modules_WGCNA/05_module_supermodule_group_effects.r",
    "tests/testthat/test-wgcna-group-effects-phase2b.R",
    "tests/testthat/test-wgcna-group-effects-contract.R",
    "tests/testthat/test-schema-validation.R"
  )
  keep <- !relative %in% explicit_exclusions &
    !grepl("^tests/(fixtures|testthat/fixtures)/", relative)
  files <- files[keep]
  relative <- relative[keep]
  tokens <- wgcna_group_consumer_tokens()
  expected <- character()
  for (i in seq_along(files)) {
    source_lines <- readLines(files[[i]], warn = FALSE, encoding = "UTF-8")
    for (j in seq_len(nrow(tokens))) {
      if (any(grepl(tokens$scan_token[[j]], source_lines, perl = TRUE))) {
        hit_lines <- grep(
          tokens$scan_token[[j]], source_lines, perl = TRUE
        )
        expected <- c(expected, paste(
          relative[[i]], tokens$consumed_column[[j]], hit_lines,
          sep = "||"
        ))
      }
    }
  }
  observed <- wgcna_group_scan_downstream_consumers()
  observed_keys <- paste(
    observed$consumer_script, observed$consumed_column,
    observed$line_number, sep = "||"
  )
  testthat::expect_setequal(observed_keys, expected)
})

testthat::test_that("Stage 05 status separates completion from readiness", {
  status <- wgcna_group_contract_status(
    "microglia",
    c(
      repo_path("R", "wgcna_group_effects_utils.R"),
      repo_path("06_modules_WGCNA", "05_module_supermodule_group_effects.r")
    ),
    canonical_primary_outputs_complete = TRUE
  )
  testthat::expect_identical(
    status$stage05_output_status,
    "phase2b_statistical_outputs_complete"
  )
  testthat::expect_identical(
    status$publication_status, "not_assessed_stage05"
  )
  testthat::expect_false(status$downstream_compatible)
  testthat::expect_identical(
    status$downstream_contract_status, "phase3_migration_required"
  )
  testthat::expect_true(status$downstream_migration_required)
  testthat::expect_true(status$should_block_execution)
  testthat::expect_match(status$source_hashes, "wgcna_group_effects_utils")
  entries <- strsplit(status$source_hashes, ";", fixed = TRUE)[[1]]
  paths <- sub("=[0-9a-f]{64}$", "", entries)
  testthat::expect_identical(paths, sort(paths))
  testthat::expect_true(all(grepl(
    "^[^=]+=[0-9a-f]{64}$", entries
  )))
})

testthat::test_that("Stage 05 source hashes cannot omit a direct dependency", {
  required <- wgcna_group_stage05_source_dependencies(
    path_processed(
      "06_modules_WGCNA", "01_WGCNA", "microglia",
      "wgcna_final_model_state.rds"
    ),
    wgcna_group_contract_paths("microglia")
  )
  relative_required <- vapply(required, relative_to, character(1))
  expected_code_and_schemas <- c(
    "06_modules_WGCNA/05_module_supermodule_group_effects.r",
    "R/paths.R", "R/dataset_config.R", "R/dataset_inputs.R",
    "R/module_contracts.R", "R/wgcna_downstream_utils.R",
    "R/wgcna_identity_contract_utils.R",
    "R/wgcna_group_effects_utils.R",
    "inst/schemas/module_group_effects.yml",
    "inst/schemas/supermodule_group_effects.yml",
    "inst/schemas/wgcna_group_effect_model_validation.yml",
    "inst/schemas/wgcna_group_effect_animal_spatial_unit_values.yml",
    "inst/schemas/wgcna_group_effect_hemisphere_values.yml",
    paste0(
      "inst/schemas/",
      "wgcna_group_effect_downstream_consumer_migration_audit.yml"
    ),
    "inst/schemas/wgcna_group_effect_contract_status.yml"
  )
  testthat::expect_true(all(
    expected_code_and_schemas %in% relative_required
  ))
  testthat::expect_true(all(c(
    "WGCNA_entity_identity_contract.csv",
    "WGCNA_module_supermodule_membership_contract.csv",
    "WGCNA_identity_contract_status.csv",
    "wgcna_final_model_state.rds"
  ) %in% basename(required)))
  testthat::expect_silent(
    wgcna_group_assert_stage05_source_dependencies(required, required)
  )
  omitted <- required[
    !grepl("R/wgcna_downstream_utils[.]R$", gsub("\\\\", "/", required))
  ]
  testthat::expect_error(
    wgcna_group_contract_status(
      "microglia", omitted,
      canonical_primary_outputs_complete = TRUE,
      required_source_paths = required
    ),
    "omits direct dependencies"
  )
})
