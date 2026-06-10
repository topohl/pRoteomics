testthat::test_that("CI dry-run workflow does not mask failures", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  workflow <- readLines(repo_path(".github", "workflows", "r-smoke-tests.yml"), warn = FALSE)
  text <- paste(workflow, collapse = "\n")
  testthat::expect_match(text, "set -euo pipefail", fixed = TRUE)
  testthat::expect_false(grepl("set \\+e", text))
  testthat::expect_false(grepl("exit 0", text))
  for (stage in c("qc", "enrichment", "modules_wgcna", "modules_downstream")) {
    testthat::expect_match(text, stage, fixed = TRUE)
  }
})
