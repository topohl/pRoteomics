testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tidyr")

source(testthat::test_path("..", "..", "R", "gsea_wgcna_concordance_utils.R"))
source(testthat::test_path(
  "..", "..", "R", "gsea_wgcna_concordance_diagnostic_utils.R"
))

testthat::test_that("effect percentiles retain canonical source-key grain", {
  handoff <- data.frame(
    dataset = "d",
    entity_level = "module",
    analysis_tier = "primary_wgcna_global",
    effect_scope = "spatial_adjusted_global",
    contrast = c("RES - CON", "SUS - CON", "SUS - RES"),
    estimate = c(1, -2, -3),
    CI_low = c(-1, -4, -6),
    CI_high = c(3, 0, 0),
    model_valid = TRUE,
    source_key = c("a", "b", "c")
  )
  out <- gwwd_effect_distribution_lookup(handoff)
  testthat::expect_equal(nrow(out), 3L)
  testthat::expect_equal(out$estimate_percentile, c(0, 0.5, 1))
  testthat::expect_equal(out$CI_width_percentile, c(0.25, 0.25, 1))
})

testthat::test_that("complete GSEA effects use a common three-contrast term set", {
  contrasts <- names(gww_formal_contrast_map())
  x <- tidyr::expand_grid(
    phenotype_contrast = contrasts,
    ID = c("GO:1", "GO:2")
  ) |>
    dplyr::mutate(
      dataset = "d",
      spatial_unit = "CA1",
      program_class = "Program",
      theme_id = "fixture_theme",
      manuscript_theme = "Fixture theme",
      theme_role = "primary",
      theme_claim_eligible = TRUE,
      theme_assignment_id = paste0("fixture_", dplyr::row_number()),
      Comparison = paste(.data$phenotype_contrast, "CA1", sep = "_"),
      GO_ID = .data$ID,
      GO_description = dplyr::if_else(.data$ID == "GO:1", "term one", "term two"),
      NES = c(1, 2, -1, -2, -2, -3),
      raw_p = c(0.01, 0.8, 0.02, 0.7, 0.03, 0.6),
      GSEA_FDR = c(0.04, 0.9, 0.05, 0.8, 0.06, 0.7),
      core_enrichment = "P1/P2",
      evidence_source_family = "ranked_GSEA",
      source_supplementary_file = "source.csv"
    )
  out <- gwwd_local_program_effects(x)
  testthat::expect_equal(nrow(out), 3L)
  testthat::expect_true(all(out$n_common_GO_terms == 2L))
  testthat::expect_true(all(out$anchor_GO_ID == "GO:1"))
  testthat::expect_true(any(out$gsea_min_FDR > 0.05))
  testthat::expect_match(
    unique(out$gsea_representation_method), "DESCRIPTIVE", fixed = TRUE
  )
})

testthat::test_that("directional labels are descriptive and fail closed", {
  out <- gwwd_pattern_label(
    res = c(1, -1, 1, 1),
    sus = c(-1, 1, 1, -1),
    direct = c(-1, 1, 1, -1),
    status_res = c("compatible", "compatible", "compatible", "unavailable"),
    status_sus = c("compatible", "compatible", "compatible", "compatible"),
    status_direct = c("compatible", "compatible", "incompatible", "compatible")
  )
  testthat::expect_equal(
    out,
    c(
      "RES_gt_CON_gt_SUS", "RES_lt_CON_lt_SUS", "inconsistent",
      "insufficient_data"
    )
  )
})

testthat::test_that("recurrent-cross-spatial direction counts local spatial units before aggregation", {
  local <- data.frame(
    dataset = "d",
    biological_program = "P",
    contrast = "RES - CON",
    gsea_spatial_unit = c("u1", "u2", "u3", "u4"),
    gsea_program_NES = c(1, 2, 3, -1),
    gsea_min_FDR = c(0.2, 0.3, 0.4, 0.5),
    gsea_median_FDR = c(0.4, 0.5, 0.6, 0.7),
    n_GSEA_FDR_supported_terms = 0L,
    n_common_GO_terms = 10L
  )
  out <- gwwd_recurrent_cross_spatial_program_effects(local)
  testthat::expect_equal(out$n_spatial_units, 4L)
  testthat::expect_equal(out$n_units_positive, 3L)
  testthat::expect_equal(out$n_units_negative, 1L)
  testthat::expect_equal(out$gsea_direction_sign, 1)
})

testthat::test_that("candidate mapping filters use explicit environment dataset", {
  script <- paste(
    readLines(testthat::test_path(
      "..", "..", "10_biological_integration",
      "06_gsea_wgcna_concordance_diagnostics.R"
    ), warn = FALSE),
    collapse = "\n"
  )
  testthat::expect_match(
    script, ".data$dataset == .env$dataset", fixed = TRUE
  )
  testthat::expect_false(grepl(
    ".data$dataset == dataset", script, fixed = TRUE
  ))
  testthat::expect_match(
    script, 'by = c("dataset", "ModuleID")', fixed = TRUE
  )
})

testthat::test_that("power diagnostic ranks unique endpoints without new tests", {
  x <- data.frame(
    wgcna_source_key = c("a", "a", "b"),
    dataset = "d",
    contrast = c("RES - CON", "RES - CON", "SUS - RES"),
    comparison_scope = "recurrent_cross_spatial",
    wgcna_spatial_unit = "global_spatial_adjusted",
    wgcna_entity_level = "module",
    wgcna_entity_id = c("m1", "m1", "m2"),
    wgcna_display_label = c("one", "one", "two"),
    wgcna_estimate = c(1, 1, -2),
    wgcna_SE = 1,
    wgcna_CI_low = c(-1, -1, -3),
    wgcna_CI_high = c(3, 3, -1),
    wgcna_p_value = c(0.2, 0.2, 0.01),
    wgcna_tier_specific_fdr = c(0.4, 0.4, 0.2),
    wgcna_tier_specific_family_id = "family",
    wgcna_tier_specific_family_size = 2L,
    wgcna_biological_n = 9L,
    wgcna_model_stability_status = "stable_mixed_model",
    biological_program = c("A", "B", "A"),
    gsea_evidence_id = c("g1", "g2", "g3"),
    overlap_FDR = c(0.01, 0.2, NA)
  )
  out <- gwwd_power_diagnostic(x)
  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_equal(out$wgcna_source_key[[1]], "b")
  testthat::expect_false(any(out$FDR_at_or_below_0_05))
  testthat::expect_true(out$any_overlap_FDR_supported[out$wgcna_source_key == "a"])
})

testthat::test_that("diagnostic source does not alter official classifier", {
  helper <- paste(readLines(testthat::test_path(
    "..", "..", "R", "gsea_wgcna_concordance_diagnostic_utils.R"
  )), collapse = "\n")
  script <- paste(readLines(testthat::test_path(
    "..", "..", "10_biological_integration",
    "06_gsea_wgcna_concordance_diagnostics.R"
  )), collapse = "\n")
  testthat::expect_false(grepl("gww_classify_concordance\\s*\\(", script))
  testthat::expect_false(grepl("stats::p.adjust", helper, fixed = TRUE))
  testthat::expect_match(
    script, "strict_classifications_modified = FALSE", fixed = TRUE
  )
  testthat::expect_match(
    script, "DESCRIPTIVE / HYPOTHESIS-GENERATING", fixed = TRUE
  )
})
