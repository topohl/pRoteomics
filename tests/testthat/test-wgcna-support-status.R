source(testthat::test_path(
  "..", "..", "R", "wgcna_support_status_utils.R"
))

testthat::test_that("finite applicable FDR controls WGCNA support", {
  observed <- wgcna_group_classify_statistical_support(
    p_value = c(0.03, 0.03, 0.03, 0.20),
    applicable_fdr = c(0.40, 0.08, 0.04, 0.40)
  )
  testthat::expect_identical(
    observed,
    c(
      "not_supported", "suggestive_FDR10", "FDR_supported",
      "not_supported"
    )
  )
})

testthat::test_that("missing adjusted results fail closed", {
  observed <- wgcna_group_classify_statistical_support(
    p_value = c(0.03, 0.20, NA_real_),
    applicable_fdr = c(NA_real_, NA_real_, NA_real_)
  )
  testthat::expect_identical(
    observed, rep("not_supported", 3L)
  )
  testthat::expect_false(any(observed == "nominal_exploratory"))
})

testthat::test_that("invalid adjusted results fail closed", {
  testthat::expect_error(
    wgcna_group_classify_statistical_support(0.03, 1.1),
    "must lie in [0, 1]", fixed = TRUE
  )
})

testthat::test_that("impact audit identifies corrected-and-failed rows", {
  fixture <- data.frame(
    dataset = "microglia",
    level = "module",
    endpoint_id = c("m1", "m2"),
    canonical_claim_entity_id = c("m1", "m2"),
    contrast = "SUS - RES",
    analysis_tier = "primary_wgcna_global",
    test_type = "named_contrast",
    effect_scope = "spatial_adjusted_global",
    spatial_unit = "global_spatial_adjusted",
    p_value = c(0.03, 0.03),
    FDR_primary_global = c(0.40, 0.04),
    FDR_secondary_global = NA_real_,
    FDR_interaction_omnibus = NA_real_,
    FDR_local_exploratory = NA_real_,
    statistical_support_status = c(
      "nominal_exploratory", "FDR_supported"
    ),
    claim_entity_role = "canonical_module",
    independent_hypothesis = TRUE,
    stringsAsFactors = FALSE
  )
  observed <- wgcna_group_support_status_impact_audit(
    fixture, "fixture.csv"
  )
  testthat::expect_equal(nrow(observed), 1L)
  testthat::expect_identical(observed$endpoint_id, "m1")
  testthat::expect_equal(observed$applicable_fdr, 0.40)
  testthat::expect_identical(observed$old_status, "nominal_exploratory")
  testthat::expect_identical(observed$new_status, "not_supported")
})
