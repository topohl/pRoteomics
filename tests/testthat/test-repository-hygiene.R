testthat::test_that("CSV configuration files cannot be hidden by a blanket ignore", {
  ignore_path <- testthat::test_path("..", "..", ".gitignore")
  ignore <- trimws(readLines(ignore_path, warn = FALSE))

  testthat::expect_false("*.csv" %in% ignore)
  testthat::expect_true("!config/*.csv" %in% ignore)
  testthat::expect_true("!config/**/*.csv" %in% ignore)
})

testthat::test_that("known root runtime artifacts are not tracked", {
  repo_root <- normalizePath(
    testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE
  )
  tracked <- system2(
    "git", c("-C", shQuote(repo_root), "ls-files"), stdout = TRUE
  )
  generated <- c(
    ".DS_Store",
    "control_spatial_identity_full_stderr.log",
    "control_spatial_identity_full_stdout.log",
    "control_spatial_identity_run.log",
    "permutedStats-actualModules.RData"
  )

  testthat::expect_length(intersect(tracked, generated), 0L)
})
