testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_group_effects_utils.R"))

testthat::test_that("legacy selected-supermodule outputs are not v5 canonical outputs", {
  allowed <- wgcna_group_allowed_output_names()
  testthat::expect_false(any(grepl(
    "selected_sus|selected_supermodule|interpretation",
    allowed, ignore.case = TRUE
  )))
  testthat::expect_true(all(c(
    "module_group_effects.csv",
    "supermodule_group_effects.csv",
    "supermodule_composition.csv",
    "WGCNA_group_effect_endpoint_provenance.csv",
    "WGCNA_group_effect_model_validation.csv"
  ) %in% allowed))
  testthat::expect_match(
    wgcna_group_legacy_reason("selected_sus_res_supermodule_contents.csv"),
    "identity and claim migration deferred",
    fixed = TRUE
  )
})

testthat::test_that("v5 primary schema retains evidence and model-governance fields", {
  schema <- wgcna_group_primary_schema()
  testthat::expect_true(all(c(
    "analysis_tier", "test_type", "effect_scope", "contrast",
    "estimate", "SE", "CI_low", "CI_high", "p_value",
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "statistical_support_status", "model_valid_for_inference",
    "primary_model_stable", "claim_allowed_model", "model_downgrade_reason",
    "endpoint_construction_method", "endpoint_provenance_status"
  ) %in% names(schema)))
  testthat::expect_identical(nrow(schema), 0L)
})

testthat::test_that("legacy semantic builders cannot become active through the runner", {
  script <- paste(readLines(testthat::test_path(
    "..", "..", "06_modules_WGCNA", "05_module_supermodule_group_effects.r"
  ), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("selected_sus_res_supermodule_contents", script, fixed = TRUE))
  testthat::expect_false(grepl("add_claim_model_fields", script, fixed = TRUE))
  testthat::expect_match(script, "wgcna_group_build_bundle", fixed = TRUE)
  testthat::expect_match(script, "wgcna_group_validate_output_bundle", fixed = TRUE)
})
