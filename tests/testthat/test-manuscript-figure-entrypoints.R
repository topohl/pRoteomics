testthat::test_that("manuscript figure contract has exact unique panel identities", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  testthat::skip_if_not_installed("yaml")

  fig2 <- manuscript_figure_contract("02")
  fig3 <- manuscript_figure_contract("03")
  ids2 <- vapply(fig2$panels, function(x) as.character(x$id), character(1))
  ids3 <- vapply(fig3$panels, function(x) as.character(x$id), character(1))

  testthat::expect_identical(ids2, paste0("2", letters[1:6]))
  testthat::expect_identical(ids3, paste0("3", letters[1:5]))
  testthat::expect_false(anyDuplicated(c(ids2, ids3)) > 0L)
  testthat::expect_identical(fig2$contract_version, "manuscript_figures_v2")
  testthat::expect_identical(fig3$contract_version, "manuscript_figures_v2")
})

testthat::test_that("Figure 2a is explicitly deferred to Illustrator", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  fig2 <- manuscript_figure_contract("02")
  panel_2a <- fig2$panels[[1]]

  testthat::expect_identical(panel_2a$id, "2a")
  testthat::expect_identical(panel_2a$automation_status, "deferred_to_illustrator")
  testthat::expect_false(manuscript_figure_panel_is_automated(panel_2a))
  testthat::expect_null(panel_2a$figure_source)
  testthat::expect_null(panel_2a$primary_source)
})

testthat::test_that("contract pins canonical panels and does not use newest-file selection", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  contract_text <- paste(readLines(manuscript_figure_contract_path(), warn = FALSE), collapse = "\n")
  helper_text <- paste(readLines(repo_path("R", "manuscript_figure_utils.R"), warn = FALSE), collapse = "\n")

  testthat::expect_match(contract_text, "main_pca_by_anatomical_island_89mm.svg", fixed = TRUE)
  testthat::expect_match(contract_text, "figure2d_compact_marker_enrichment_differences_v2.svg", fixed = TRUE)
  testthat::expect_match(contract_text, "figure2f_control_anatomical_GO_regions_CA1layers_grouped.svg", fixed = TRUE)
  testthat::expect_match(contract_text, "wgcna_circular_heatmap_neuron_neuropil.svg", fixed = TRUE)
  testthat::expect_match(contract_text, "figure3b_stage07_wgcna_effects.svg", fixed = TRUE)
  testthat::expect_false(grepl("latest_input_candidate", helper_text, fixed = TRUE))
  testthat::expect_false(grepl("list.files", helper_text, fixed = TRUE))
})

testthat::test_that("hemisphere and biological-unit declarations stay explicit", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  fig2 <- manuscript_figure_contract("02")
  fig3 <- manuscript_figure_contract("03")
  by_id <- function(fig, id) fig$panels[[match(id, vapply(fig$panels, function(x) as.character(x$id), character(1)))]]

  testthat::expect_identical(by_id(fig2, "2b")$biological_unit, "measured_sample")
  testthat::expect_match(by_id(fig2, "2c")$hemisphere_handling, "not_a_stress_inference", fixed = TRUE)
  testthat::expect_identical(by_id(fig2, "2d")$biological_unit, "animal")
  testthat::expect_match(by_id(fig2, "2d")$hemisphere_handling, "equal_valid_hemispheres", fixed = TRUE)
  testthat::expect_match(by_id(fig2, "2e")$hemisphere_handling, "blocked_by_AnimalID", fixed = TRUE)
  testthat::expect_match(by_id(fig3, "3a")$hemisphere_handling, "animal_level_hemisphere_averaged", fixed = TRUE)
  testthat::expect_match(by_id(fig3, "3b")$hemisphere_handling, "equal_weight_available_LR", fixed = TRUE)
  testthat::expect_match(by_id(fig3, "3c")$hemisphere_handling, "not_applicable", fixed = TRUE)
  testthat::expect_match(by_id(fig3, "3e")$hemisphere_handling, "animal_level_hemisphere_averaged", fixed = TRUE)
})

testthat::test_that("incomplete candidates require an explicit output root", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  canonical <- manuscript_figure_args("--allow-incomplete")
  isolated <- manuscript_figure_args(c("--allow-incomplete", "--output-root", tempfile("figure-candidate-")))
  testthat::expect_true(canonical$allow_incomplete)
  testthat::expect_false(canonical$output_explicit)
  testthat::expect_true(isolated$allow_incomplete)
  testthat::expect_true(isolated$output_explicit)
})

testthat::test_that("SVG assembly is self-contained and preserves panel letters", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  testthat::skip_if_not_installed("base64enc")
  root <- tempfile("figure-assembly-")
  dir.create(root)
  a <- file.path(root, "a.svg")
  b <- file.path(root, "b.svg")
  writeLines('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"><rect width="10" height="10" fill="red"/></svg>', a)
  writeLines('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"><rect width="10" height="10" fill="blue"/></svg>', b)
  panels <- list(
    list(id = "9a", row = 1L, col = 1L),
    list(id = "9b", row = 1L, col = 2L)
  )
  target <- file.path(root, "assembled.svg")
  manuscript_figure_assemble_svg(
    c("9a" = a, "9b" = b), panels,
    list(width_mm = 183, height_mm = 80), target
  )
  text <- paste(readLines(target, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "data:image/svg+xml;base64", fixed = TRUE)
  testthat::expect_match(text, ">a</text>", fixed = TRUE)
  testthat::expect_match(text, ">b</text>", fixed = TRUE)
  testthat::expect_false(grepl(a, text, fixed = TRUE))
})

testthat::test_that("entry points are registered as downstream-only manuscript steps", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  testthat::skip_if_not_installed("yaml")
  registry <- yaml::read_yaml(repo_path("pipeline.yml"))
  entries <- unlist(lapply(registry$stages, function(stage) {
    vapply(stage$scripts %||% list(), function(x) as.character(x$script), character(1))
  }), use.names = FALSE)
  testthat::expect_true("figures/figure_02.R" %in% entries)
  testthat::expect_true("figures/figure_03.R" %in% entries)
})
