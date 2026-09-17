testthat::test_that("dataset-scoped active scripts use shared CLI dataset resolution", {
  source(testthat::test_path("..", "..", "R", "paths.R"))

  scripts <- c(
    "analysis/01_preprocessing/extract_protigy_contrasts.R",
    "analysis/01_preprocessing/map_protein_identifiers.R",
    "analysis/04_differential_abundance/run_clusterprofiler_enrichment.R",
    "analysis/04_differential_abundance/compare_go_enrichment.R",
    "analysis/04_differential_abundance/annotate_neuropil_reference.R",
    "analysis/07_spatial_networks/build_differential_networks.R",
    "analysis/07_spatial_networks/test_network_stability.R",
    "analysis/07_spatial_networks/test_differential_network_stability.R",
    "analysis/07_spatial_networks/render_differential_network_figures.R",
    "analysis/07_spatial_networks/render_network_chord_diagram.R",
    "analysis/08_integration/test_behaviour_proteomics_associations.R",
    "analysis/08_integration/test_network_behaviour_coupling.R"
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
  source(repo_path("R", "dataset_config.R"))
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
