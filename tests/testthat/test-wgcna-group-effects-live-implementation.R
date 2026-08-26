testthat::test_that("Phase 2B/v5 is the sole WGCNA group-effect implementation", {
  path <- testthat::test_path(
    "..", "..", "R", "wgcna_group_effects_utils.R"
  )
  expressions <- parse(file = path)
  assigned_names <- vapply(expressions, function(expression) {
    if (
      is.call(expression) &&
        identical(expression[[1]], as.name("<-")) &&
        is.symbol(expression[[2]])
    ) {
      as.character(expression[[2]])
    } else {
      NA_character_
    }
  }, character(1))
  assigned_names <- stats::na.omit(assigned_names)

  live <- c(
    "wgcna_group_fit_attempt",
    "wgcna_group_run_level",
    "wgcna_group_apply_fdr"
  )
  testthat::expect_identical(
    vapply(live, function(name) sum(assigned_names == name), integer(1)),
    stats::setNames(rep(1L, length(live)), live)
  )
  testthat::expect_false(any(grepl("phase2_legacy", assigned_names)))

  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path(
    "..", "..", "R", "wgcna_downstream_utils.R"
  ))
  source(path)
  testthat::expect_identical(
    wgcna_group_effects_contract_version(),
    "wgcna_group_effects_phase2b_v5"
  )
})
