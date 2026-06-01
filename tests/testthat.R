#!/usr/bin/env Rscript

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Package 'testthat' is required for tests. Install it with install.packages('testthat') or renv::restore().", call. = FALSE)
}
testthat::test_dir("tests/testthat", reporter = "summary")
