testthat::test_that("pipeline.yml is valid and references existing active scripts", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "pipeline_registry.R"))

  testthat::skip_if_not_installed("yaml")
  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
  testthat::expect_silent(validate_pipeline_scripts_exist(registry))
  testthat::expect_silent(validate_run_order_against_registry(registry))
  steps <- pipeline_steps(registry, "enrichment", dataset = "microglia")
  testthat::expect_true("04_differential_expression_enrichment/01_clusterProfiler.r" %in% steps$script)

  network_steps <- pipeline_steps(registry, "networks", dataset = "microglia")
  testthat::expect_false(any(network_steps$supported))
  testthat::expect_false("microglia" %in% network_steps$supported_datasets)

  behavior_steps <- pipeline_steps(registry, "behavior", dataset = "microglia")
  testthat::expect_false(any(behavior_steps$supported))
})

testthat::test_that("README and RUN_ORDER do not present legacy scripts as active", {
  source(testthat::test_path("..", "..", "R", "paths.R"))

  readme <- paste(readLines(repo_path("README.md"), warn = FALSE), collapse = "\n")
  run_order <- paste(readLines(repo_path("RUN_ORDER.md"), warn = FALSE), collapse = "\n")

  active_blocks <- paste(readme, run_order, sep = "\n")
  testthat::expect_false(grepl("Backward-compatible retained names", active_blocks, fixed = TRUE))
  testthat::expect_false(grepl("04_neuropil_contamination_annotation.r", run_order, fixed = TRUE))
  testthat::expect_false(grepl("91_module_score.r", run_order, fixed = TRUE))
  testthat::expect_true(grepl("09_pride_submission/", readme, fixed = TRUE))
  testthat::expect_true(grepl("legacy", readme, ignore.case = TRUE))
})
