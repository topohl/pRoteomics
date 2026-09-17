testthat::test_that("dataset-scoped active scripts use shared CLI dataset resolution", {
  source(testthat::test_path("..", "..", "R", "paths.R"))

  scripts <- c(
    "analysis/preprocessing/extract_protigy_contrasts.R",
    "analysis/preprocessing/map_protein_identifiers.R",
    "analysis/differential_abundance/run_clusterprofiler_enrichment.R",
    "analysis/differential_abundance/compare_go_enrichment.R",
    "analysis/differential_abundance/annotate_neuropil_reference.R",
    "analysis/spatial_networks/build_differential_networks.R",
    "analysis/spatial_networks/test_network_stability.R",
    "analysis/spatial_networks/test_differential_network_stability.R",
    "analysis/spatial_networks/render_differential_network_figures.R",
    "analysis/spatial_networks/render_network_chord_diagram.R",
    "analysis/integration/test_behaviour_proteomics_associations.R",
    "analysis/integration/test_network_behaviour_coupling.R"
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
