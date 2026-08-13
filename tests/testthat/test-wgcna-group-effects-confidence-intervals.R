testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_identity_contract_utils.R"))
source(testthat::test_path("..", "..", "R", "wgcna_group_effects_utils.R"))

testthat::test_that("emmeans confidence intervals use the same finite-df inference", {
  fixture <- data.frame(
    response = c(0.2, 0.4, 0.5, 1.0, 1.1, 1.4, -0.2, 0.0, 0.3),
    StressGroup = factor(rep(c("CON", "RES", "SUS"), each = 3L))
  )
  fit <- stats::lm(response ~ StressGroup, data = fixture)
  emm <- emmeans::emmeans(fit, ~ StressGroup)
  methods <- wgcna_group_predeclared_contrasts(levels(fixture$StressGroup))
  observed <- wgcna_group_emmeans_contrast_summary(emm, methods)
  direct <- as.data.frame(confint(
    emmeans::contrast(emm, method = methods, adjust = "none"),
    level = 0.95, adjust = "none"
  ))

  testthat::expect_equal(observed$CI_low, direct$lower.CL, tolerance = 1e-14)
  testthat::expect_equal(observed$CI_high, direct$upper.CL, tolerance = 1e-14)
  testthat::expect_true(all(is.finite(observed$df)))
  testthat::expect_true(all(observed$CI_method == "emmeans_contrast_summary"))
  testthat::expect_true(all(observed$CI_level == 0.95))
  testthat::expect_true(all(observed$CI_df_method == "finite_df_from_emmeans"))
})

testthat::test_that("two-sided CI exclusion reconciles with unadjusted p-values", {
  table <- data.frame(
    dataset = "fixture", level = "module", endpoint_id = paste0("m", 1:2),
    effect_scope = "spatial_adjusted_global",
    spatial_unit = "global_spatial_adjusted",
    contrast = c("RES - CON", "SUS - CON"),
    test_type = "named_contrast", estimate = c(1, 0.2), SE = c(0.2, 0.2),
    CI_low = c(0.5, -0.3), CI_high = c(1.5, 0.7),
    CI_method = "emmeans_contrast_summary", CI_level = 0.95,
    CI_df_method = "finite_df_from_emmeans", p_value = c(0.01, 0.4),
    df_den = 6
  )
  testthat::expect_equal(nrow(
    wgcna_group_ci_p_compatibility_violations(table)
  ), 0L)

  table$CI_low[[2]] <- 0.01
  testthat::expect_equal(nrow(
    wgcna_group_ci_p_compatibility_violations(table)
  ), 1L)
})

testthat::test_that("legacy normal-multiplier intervals fail provenance", {
  table <- data.frame(
    dataset = "fixture", level = "module", endpoint_id = "m1",
    effect_scope = "spatial_adjusted_global",
    spatial_unit = "global_spatial_adjusted", contrast = "RES - CON",
    test_type = "named_contrast", estimate = 1, SE = 0.2,
    CI_low = 1 - 1.96 * 0.2, CI_high = 1 + 1.96 * 0.2,
    CI_method = "normal_1.96_SE", CI_level = 0.95,
    CI_df_method = "not_used", p_value = 0.01, df_den = 4
  )
  testthat::expect_equal(nrow(
    wgcna_group_ci_p_compatibility_violations(table)
  ), 1L)
})
