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
})
