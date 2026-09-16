testthat::test_that("manuscript figure contract has exact unique panel identities", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  testthat::skip_if_not_installed("yaml")

  fig2 <- manuscript_figure_contract("02")
  fig3 <- manuscript_figure_contract("03")
  ids2 <- vapply(fig2$panels, function(x) as.character(x$id), character(1))
  ids3 <- vapply(fig3$panels, function(x) as.character(x$id), character(1))

  # Promoted in Phase 5B: Figures 2 and 3 are the final_truth_v9 generation,
  # a-h and a-i, which is the structure Results 2 and 3 were written against.
  testthat::expect_identical(ids2, paste0("2", letters[1:8]))
  testthat::expect_identical(ids3, paste0("3", letters[1:9]))
  testthat::expect_false(anyDuplicated(c(ids2, ids3)) > 0L)
  testthat::expect_identical(fig2$contract_version,
                             "manuscript_figures_v3_final_truth_v9_promoted")
  testthat::expect_identical(fig3$contract_version,
                             "manuscript_figures_v3_final_truth_v9_promoted")
})

testthat::test_that("Figure 2a is a rendered anatomy schematic, not a placeholder", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  fig2 <- manuscript_figure_contract("02")
  panel_2a <- fig2$panels[[1]]

  # Before the promotion this slot was deferred_to_illustrator with no asset at
  # all. The promoted generation draws it, so it must now behave like every
  # other automated panel.
  testthat::expect_identical(panel_2a$id, "2a")
  testthat::expect_identical(as.character(panel_2a$renderer), "f9_schematic")
  testthat::expect_true(manuscript_figure_panel_is_automated(panel_2a))
  testthat::expect_true(nzchar(as.character(panel_2a$figure_source)))
  testthat::expect_true(nzchar(as.character(panel_2a$primary_source)))
})

testthat::test_that("contract pins canonical panels and does not use newest-file selection", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  contract_text <- paste(readLines(manuscript_figure_contract_path(), warn = FALSE), collapse = "\n")
  helper_text <- paste(readLines(repo_path("R", "manuscript_figure_utils.R"), warn = FALSE), collapse = "\n")

  # Each promoted panel names one exact asset produced by the v9 entry points.
  for (asset in c("v9_schematic.svg", "v9_depth.svg", "v9_pca.svg",
                  "v9_fingerprint.svg", "v9_compartment.svg", "v9_bilateral_main.svg",
                  "v9_external_main.svg", "v9_internal_main.svg",
                  "v9_dap_track.svg", "v9_atlas.svg", "v9_bridge.svg",
                  "v9_curve_syn.svg", "v9_curve_rna.svg", "v9_curve_ox.svg",
                  "v9_prot_syn.svg", "v9_prot_rna.svg", "v9_prot_ox.svg")) {
    testthat::expect_match(contract_text, asset, fixed = TRUE)
  }
  testthat::expect_false(grepl("latest_input_candidate", helper_text, fixed = TRUE))
  testthat::expect_false(grepl("list.files", helper_text, fixed = TRUE))
})

testthat::test_that("hemisphere and biological-unit declarations stay explicit", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  # Every promoted panel must still declare both fields. The promoted generation
  # is animal-level and bilaterally aggregated throughout; the point of the test
  # is that the declaration is never silently dropped, which is how the previous
  # generation lost it on the one panel that was entirely about hemispheres.
  for (k in c("02", "03")) {
    fig <- manuscript_figure_contract(k)
    for (p in fig$panels) {
      testthat::expect_true(nzchar(as.character(p$biological_unit %||% "")),
                            info = paste("missing biological_unit:", p$id))
      testthat::expect_true(nzchar(as.character(p$hemisphere_handling %||% "")),
                            info = paste("missing hemisphere_handling:", p$id))
    }
  }
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
