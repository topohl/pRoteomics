testthat::test_that("dataset-scoped active scripts use shared CLI dataset resolution", {
  source(testthat::test_path("..", "..", "R", "paths.R"))

  scripts <- c(
    "01_preprocessing/03_gct_extractR.r",
    "02_id_mapping/01_MapThatProt_batch.r",
    "04_differential_expression_enrichment/01_clusterProfiler.r",
    "04_differential_expression_enrichment/02_compareGO.r",
    "04_differential_expression_enrichment/04_neuropil_contamination_annotation.r",
    "07_spatial_networks/02_differential_networks.r",
    "07_spatial_networks/03_bootstrap_network_stability.r",
    "07_spatial_networks/04_bootstrap_differential_network_stability.r",
    "07_spatial_networks/05_bootstrap_differential_network_figures.r",
    "07_spatial_networks/06_chord_diagram.r",
    "08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r",
    "08_behavior_physio_coupling/02_network_behavior_coupling.r"
  )

  for (script in scripts) {
    txt <- paste(readLines(repo_path(script), warn = FALSE), collapse = "\n")
    testthat::expect_true(
      grepl("current_dataset_from_cli\\(", txt),
      info = paste(script, "should honor --dataset before resolving paths")
    )
  }
})

testthat::test_that("dataset CLI helper prefers --dataset over environment and defaults", {
  source(testthat::test_path("..", "..", "R", "dataset_config.R"))
  old <- Sys.getenv("PROTEOMICS_DATASET", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("PROTEOMICS_DATASET") else Sys.setenv(PROTEOMICS_DATASET = old)
  }, add = TRUE)

  Sys.setenv(PROTEOMICS_DATASET = "neuron_neuropil")
  testthat::expect_identical(
    current_dataset_from_cli(args = c("--dataset", "microglia")),
    "microglia"
  )
  testthat::expect_identical(Sys.getenv("PROTEOMICS_DATASET"), "microglia")
})
