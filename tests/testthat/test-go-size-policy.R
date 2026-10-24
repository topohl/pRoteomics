source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "sus_res_spatial_dap_atlas_utils.R"))

# Gene-set tests follow the main clusterProfiler GSEA/ORA term-size limits
# (config/clusterProfiler_config.yml analysis$min_gs_size / max_gs_size).

read_script <- function(name) {
  paste(readLines(repo_path("analysis/differential_abundance", name), warn = FALSE), collapse = "\n")
}

testthat::test_that("DAP-set GO ORA reads the shared configured term-size limits", {
  testthat::skip_if_not_installed("yaml")
  cfg <- yaml::read_yaml(repo_path("config", "clusterProfiler_config.yml"))$analysis
  testthat::expect_identical(
    sus_res_go_size_limits(),
    list(min_gs_size = as.integer(cfg$min_gs_size), max_gs_size = as.integer(cfg$max_gs_size))
  )
  missing_key <- tempfile(fileext = ".yml")
  yaml::write_yaml(list(analysis = list(min_gs_size = 10L)), missing_key)
  testthat::expect_error(sus_res_go_size_limits(missing_key), "max_gs_size")
  script <- read_script("build_sus_res_dap_atlas.R")
  testthat::expect_match(script, "go_size_limits <- sus_res_go_size_limits()", fixed = TRUE)
  testthat::expect_match(script, "minGSSize = go_size_limits$min_gs_size", fixed = TRUE)
  testthat::expect_match(script, "maxGSSize = go_size_limits$max_gs_size", fixed = TRUE)
  testthat::expect_false(grepl("maxGSSize = 500", script, fixed = TRUE))
})

testthat::test_that("microglia fallback methods never cap below the largest mapped signature", {
  script <- read_script("test_microglia_targeted_signatures.R")
  testthat::expect_match(script, "TARGETED_GSEA_CFG$analysis$max_gs_size", fixed = TRUE)
  testthat::expect_match(script, "max_set_size <- max(TARGETED_MAX_GS_SIZE,", fixed = TRUE)
  testthat::expect_false(grepl("max(Size|GSSize) = 500", script))
  testthat::expect_identical(
    lengths(regmatches(script, gregexpr("max(Size|GSSize) = max_set_size", script))),
    4L
  )
})
