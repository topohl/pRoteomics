source(testthat::test_path("..", "..", "R", "plotting_nature.R"))

testthat::test_that("manuscript dimensions and millimetre conversion are stable", {
  testthat::expect_identical(
    nature_dimensions_mm(),
    c(single_column = 89, double_column = 183, maximum_height = 170)
  )
  testthat::expect_equal(mm_to_in(c(25.4, 89, 183)), c(1, 89 / 25.4, 183 / 25.4))
})

testthat::test_that("required semantic palette roles and identities are stable", {
  testthat::expect_setequal(
    names(NATURE_SEMANTIC_PALETTES),
    c("group", "dataset", "signed", "support", "jaccard")
  )
  testthat::expect_identical(
    nature_palette("group"),
    c(CON = "#3E3C6F", RES = "#C6C3BB", SUS = "#E63A48")
  )
  testthat::expect_identical(
    nature_palette("dataset"),
    c(microglia = "#A8D5CF", neuropil = "#2F6F62", soma = "#7F7F7F")
  )
})

testthat::test_that("opt-in manuscript typography stays publication-legible", {
  testthat::expect_identical(
    nature_manuscript_text_sizes_pt(),
    c(normal = 7, dense = 6.5, axis_title = 7.5, panel_letter = 9, title = 8)
  )
})

testthat::test_that("anatomical-island palette is centralized and complete", {
  palette <- nature_anatomical_island_palette()
  testthat::expect_length(palette, 18L)
  testthat::expect_identical(anyDuplicated(names(palette)), 0L)
  testthat::expect_true(all(c(
    "Microglia-enriched ROI::CA1", "Neuropil::DG-PO", "Soma::DG-SP"
  ) %in% names(palette)))
})
