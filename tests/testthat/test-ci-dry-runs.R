testthat::test_that("CI workflows install dependencies fail closed", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  workflow_paths <- repo_path(
    ".github", "workflows", c("dry-run.yml", "r-smoke-tests.yml")
  )
  for (workflow_path in workflow_paths) {
    workflow <- readLines(workflow_path, warn = FALSE)
    text <- paste(workflow, collapse = "\n")
    testthat::expect_match(
      text, "uses: r-lib/actions/setup-r-dependencies@v2", fixed = TRUE,
      info = workflow_path
    )
    testthat::expect_false(
      grepl("install[.]packages[(]", text),
      info = paste(workflow_path, "must not use non-fail-fast install.packages")
    )
    dependency_step <- grep("setup-r-dependencies@v2", workflow, fixed = TRUE)
    first_validation_step <- grep(
      "Validate tests|Parse active R scripts", workflow
    )[[1]]
    testthat::expect_lt(dependency_step[[1]], first_validation_step)
  }
})

testthat::test_that("contract smoke workflow does not mask failures", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  workflow <- readLines(
    repo_path(".github", "workflows", "r-smoke-tests.yml"), warn = FALSE
  )
  text <- paste(workflow, collapse = "\n")
  testthat::expect_match(text, "set -euo pipefail", fixed = TRUE)
  testthat::expect_false(grepl("set \\+e", text))
  testthat::expect_false(grepl("exit 0", text))
  for (stage in c("qc", "enrichment", "modules_wgcna", "modules_downstream")) {
    testthat::expect_match(text, stage, fixed = TRUE)
  }
})
