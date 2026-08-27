testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "wgcna_stage07_semantic_utils.R"))

plot_levels <- c("RES - CON", "SUS - CON", "SUS - RES")

# Mirrors the Stage 07 module_join shape reaching
# plot_module_effects_by_supermodule(): named contrasts plus the one
# interaction omnibus row each module carries.
module_effect_fixture <- function() {
  data.frame(
    module_id = c("M1", "M1", "M1", "M1", "M2", "M2"),
    test_type = c(
      rep("named_contrast", 3), "interaction_omnibus",
      "named_contrast", "named_contrast"
    ),
    contrast = c(
      "RES - CON", "SUS - CON", "SUS - RES",
      "StressGroup x SpatialUnit omnibus",
      "RES - CON", "SUS - CON"
    ),
    estimate = c(0.4, -0.7, 0.2, NA_real_, 1.1, -0.3),
    p_value = c(0.01, 0.02, 0.30, 0.04, 0.005, 0.44),
    statistic = c(2.1, -2.6, 1.0, 5.7, 3.3, -0.8),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("interaction_omnibus rows are excluded from the contrast dot plot input", {
  out <- wgcna_stage07_filter_named_contrast_estimable(
    module_effect_fixture(), plot_levels
  )
  testthat::expect_false("interaction_omnibus" %in% out$test_type)
  testthat::expect_false(
    "StressGroup x SpatialUnit omnibus" %in% as.character(out$contrast)
  )
  testthat::expect_equal(nrow(out), 5L)
})

testthat::test_that("named_contrast rows pass through unchanged", {
  fixture <- module_effect_fixture()
  out <- wgcna_stage07_filter_named_contrast_estimable(fixture, plot_levels)
  expected <- fixture[fixture$test_type == "named_contrast", , drop = FALSE]
  testthat::expect_equal(out, expected, ignore_attr = "row.names")
  testthat::expect_identical(out$estimate, expected$estimate)
  testthat::expect_identical(out$p_value, expected$p_value)
})

testthat::test_that("the filtered input leaves no NA aesthetic for geom_point to drop", {
  out <- wgcna_stage07_filter_named_contrast_estimable(
    module_effect_fixture(), plot_levels
  )
  # x aesthetic, fill aesthetic: the two sources of the historical
  # "Removed N rows containing missing values (geom_point())" warnings.
  testthat::expect_true(all(is.finite(out$estimate)))
  testthat::expect_false(
    any(is.na(factor(out$contrast, levels = plot_levels)))
  )
})

testthat::test_that("a non-estimable named contrast is excluded too", {
  fixture <- module_effect_fixture()
  fixture$estimate[1] <- NA_real_
  fixture$estimate[5] <- Inf
  out <- wgcna_stage07_filter_named_contrast_estimable(fixture, plot_levels)
  testthat::expect_equal(nrow(out), 3L)
  testthat::expect_true(all(is.finite(out$estimate)))
})

testthat::test_that("a contrast outside the plotted levels is excluded", {
  fixture <- module_effect_fixture()
  fixture$contrast[2] <- "RES - SUS"
  out <- wgcna_stage07_filter_named_contrast_estimable(fixture, plot_levels)
  testthat::expect_false("RES - SUS" %in% as.character(out$contrast))
  testthat::expect_equal(nrow(out), 4L)
})

testthat::test_that("filtering requires the schema columns it depends on", {
  fixture <- module_effect_fixture()
  testthat::expect_error(
    wgcna_stage07_filter_named_contrast_estimable(
      fixture[, setdiff(names(fixture), "estimate")], plot_levels
    ),
    "missing required column"
  )
})

testthat::test_that("zero-row top-supermodule input leaves no stale SVG or PDF", {
  fig_dir <- withr::local_tempdir()
  svg <- file.path(fig_dir, "top_supermodule_effects_dotplot.svg")
  pdf <- file.path(fig_dir, "top_supermodule_effects_dotplot.pdf")
  writeLines("stale historical svg", svg)
  writeLines("stale historical pdf", pdf)

  testthat::expect_message(
    removed <- wgcna_stage07_remove_stale_figure(svg),
    "Removed stale figure"
  )

  testthat::expect_false(file.exists(svg))
  testthat::expect_false(file.exists(pdf))
  testthat::expect_setequal(basename(removed), basename(c(svg, pdf)))
})

testthat::test_that("stale removal is a no-op when nothing is on disk", {
  fig_dir <- withr::local_tempdir()
  svg <- file.path(fig_dir, "top_supermodule_effects_dotplot.svg")
  testthat::expect_silent(wgcna_stage07_remove_stale_figure(svg))
  testthat::expect_false(file.exists(svg))
})

testthat::test_that("non-empty top-supermodule input still produces current outputs", {
  # The zero-row branch is the only caller of the invalidation helper, so a
  # non-empty result must keep the figure this run wrote.
  fig_dir <- withr::local_tempdir()
  svg <- file.path(fig_dir, "top_supermodule_effects_dotplot.svg")
  pdf <- file.path(fig_dir, "top_supermodule_effects_dotplot.pdf")
  writeLines("current svg", svg)
  writeLines("current pdf", pdf)

  top_super <- data.frame(
    test_type = "named_contrast",
    contrast = c("RES - CON", "SUS - CON"),
    estimate = c(0.5, -0.9),
    p_value = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  kept <- wgcna_stage07_filter_named_contrast_estimable(top_super, plot_levels)
  testthat::expect_equal(nrow(kept), 2L)

  testthat::expect_true(file.exists(svg))
  testthat::expect_true(file.exists(pdf))
})
