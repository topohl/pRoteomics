source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::test_that("canonical null coalescing is scalar-safe and vector-safe", {
  testthat::expect_identical(NULL %||% "fallback", "fallback")
  testthat::expect_identical(character() %||% "fallback", "fallback")
  testthat::expect_identical(NA_character_ %||% "fallback", "fallback")
  testthat::expect_identical("value" %||% "fallback", "value")

  complete <- c("a", "b")
  partial <- c("a", NA_character_)
  all_missing <- c(NA_character_, NA_character_)
  testthat::expect_identical(complete %||% "fallback", complete)
  testthat::expect_identical(partial %||% "fallback", partial)
  testthat::expect_identical(all_missing %||% "fallback", all_missing)
})

testthat::test_that("helper source order cannot change coalescing semantics", {
  expected <- proteomics_null_coalesce
  source(repo_path("R", "validation_utils.R"))
  source(repo_path("R", "qc_exploration_utils.R"))
  testthat::expect_identical(`%||%`, expected)
  testthat::expect_identical(c("a", NA_character_) %||% "fallback",
                             c("a", NA_character_))
})

testthat::test_that("active code has one null-coalescing definition", {
  relative <- system2(
    "git", c("-C", shQuote(repo_root()), "ls-files", "--", "*.R", "*.r"),
    stdout = TRUE, stderr = TRUE
  )
  testthat::expect_null(attr(relative, "status"))
  relative <- unique(c(relative, "R/null_coalescing.R"))
  active <- !grepl(
    "(^|/)(99_deprecated|legacy)(/|$)", relative, ignore.case = TRUE
  )
  definitions <- vapply(repo_path(relative[active]), function(path) {
    any(grepl(
      "^[[:space:]]*`%\\|\\|%`[[:space:]]*<-",
      readLines(path, warn = FALSE)
    ))
  }, logical(1))
  testthat::expect_identical(
    relative[active][definitions],
    "R/null_coalescing.R"
  )
})
