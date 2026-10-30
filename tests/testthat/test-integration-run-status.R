testthat::local_edition(3)

# The run-status gate (R/statistics/integration_utils.R): a reader takes an
# integration step's tables only from a completed run, unchanged since, whose
# shared inputs have not changed. Everything is written under a temporary root.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "integration_utils.R"))

local_run_root <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root), .local_envir = env)
  root
}

testthat::test_that("a reader refuses a step without a completed run", {
  root <- local_run_root()
  paths <- integration_dirs("fixture_producer", "global", create = TRUE)
  table <- file.path(paths$tables, "out.csv")
  input <- file.path(root, "input.csv")
  writeLines("a\n1", input)
  testthat::expect_error(integration_require_completed_run("fixture_producer", consumer = "reader"),
                         "no run status")
  integration_mark_run_started(paths, "fixture_producer.R")
  writeLines("x\n1", table)
  testthat::expect_error(integration_require_completed_run("fixture_producer", table, consumer = "reader"),
                         "did not complete")
  integration_mark_run_completed(paths, "out.csv", c(shared_input = input))
  testthat::expect_silent(integration_require_completed_run("fixture_producer", table,
                                                            c(shared_input = input), consumer = "reader"))
  # a new start invalidates the earlier completion until it completes again
  integration_mark_run_started(paths, "fixture_producer.R")
  testthat::expect_error(integration_require_completed_run("fixture_producer", table, consumer = "reader"),
                         "did not complete")
})

testthat::test_that("a reader refuses a rewritten table or a changed shared input", {
  root <- local_run_root()
  paths <- integration_dirs("fixture_producer", "global", create = TRUE)
  table <- file.path(paths$tables, "out.csv")
  input <- file.path(root, "input.csv")
  writeLines("a\n1", input)
  integration_mark_run_started(paths, "fixture_producer.R")
  writeLines("x\n1", table)
  integration_mark_run_completed(paths, "out.csv", c(shared_input = input))
  writeLines("a\n2", input)
  testthat::expect_error(integration_require_completed_run("fixture_producer", table,
                                                           c(shared_input = input), consumer = "reader"),
                         "changed since its run")
  # an input the reader does not share by name is not compared
  testthat::expect_silent(integration_require_completed_run("fixture_producer", table,
                                                            c(other_input = input), consumer = "reader"))
  writeLines("x\n2", table)
  testthat::expect_error(integration_require_completed_run("fixture_producer", table, consumer = "reader"),
                         "not the table its completed run wrote")
})

testthat::test_that("a run cannot complete without its tables", {
  local_run_root()
  paths <- integration_dirs("fixture_producer", "global", create = TRUE)
  integration_mark_run_started(paths, "fixture_producer.R")
  testthat::expect_error(integration_mark_run_completed(paths, "missing.csv", character()),
                         "cannot complete without")
})
