testthat::test_that("spatial network script parses --dataset before current_dataset", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  txt <- paste(readLines(repo_path("07_spatial_networks/01_network_spatial_relations.r"), warn = FALSE), collapse = "\n")

  testthat::expect_true(grepl("args <- commandArgs\\(trailingOnly = TRUE\\)", txt))
  testthat::expect_true(grepl("dataset_cli <- arg_value\\(\"--dataset\"", txt))
  testthat::expect_true(grepl("Sys.setenv\\(PROTEOMICS_DATASET = validate_dataset\\(dataset_cli, source = \"--dataset\"\\)\\)", txt))

  pos_dataset_cli <- regexpr("dataset_cli <-", txt, fixed = TRUE)[1]
  pos_current_dataset <- regexpr("SPATIAL_DATASET <- current_dataset()", txt, fixed = TRUE)[1]
  testthat::expect_gt(pos_current_dataset, pos_dataset_cli)

  testthat::expect_false(grepl('assert_dataset_capability\\(SPATIAL_DATASET, "layer"', txt))
  testthat::expect_true(grepl('assert_dataset_capability\\(SPATIAL_DATASET, "region"', txt))
  testthat::expect_true(grepl("spatial_unit <- if \\(identical\\(SPATIAL_DATASET, \"neuron_neuropil\"\\)\\) \"region_layer\" else \"region\"", txt))
  testthat::expect_true(grepl("SpatialLabel", txt, fixed = TRUE))
  testthat::expect_true(grepl('SUBSTEP_ID <- file.path\\("network_spatial_relations", SPATIAL_DATASET, spatial_unit\\)', txt))
})

testthat::test_that("dataset_config records spatial-unit contracts", {
  source(testthat::test_path("..", "..", "R", "dataset_config.R"))
  contracts <- dataset_contracts()

  testthat::expect_false(contracts$neuron_soma$layer)
  testthat::expect_false(contracts$neuron_soma$region_layer)
  testthat::expect_identical(contracts$neuron_soma$spatial_unit, "region")

  testthat::expect_false(contracts$microglia$layer)
  testthat::expect_false(contracts$microglia$region_layer)
  testthat::expect_identical(contracts$microglia$spatial_unit, "region")

  testthat::expect_true(contracts$neuron_neuropil$layer)
  testthat::expect_true(contracts$neuron_neuropil$region_layer)
  testthat::expect_identical(contracts$neuron_neuropil$spatial_unit, "region_layer")
})

testthat::test_that("spatial network dry-run honors --dataset", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  run <- function(dataset) {
    old_wd <- setwd(repo_path())
    on.exit(setwd(old_wd), add = TRUE)
    out <- suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      c("07_spatial_networks/01_network_spatial_relations.r", "--dataset", dataset, "--dry-run"),
      stdout = TRUE,
      stderr = TRUE
    ))
    paste(out, collapse = "\n")
  }

  microglia <- run("microglia")
  testthat::expect_true(grepl("Dataset: microglia", microglia, fixed = TRUE))
  testthat::expect_true(grepl("spatial_unit: region", microglia, fixed = TRUE))
  testthat::expect_true(grepl("pgmatrix_imputed_microglia_", microglia, fixed = TRUE))
  testthat::expect_true(grepl("network_spatial_relations/microglia/region", microglia, fixed = TRUE))
  testthat::expect_false(grepl("network_spatial_relations/neuron_neuropil/region_layer", microglia, fixed = TRUE))

  soma <- run("neuron_soma")
  testthat::expect_true(grepl("Dataset: neuron_soma", soma, fixed = TRUE))
  testthat::expect_true(grepl("spatial_unit: region", soma, fixed = TRUE))
  testthat::expect_true(grepl("pgmatrix_imputed_neuron_soma_", soma, fixed = TRUE))

  neuropil <- run("neuron_neuropil")
  testthat::expect_true(grepl("Dataset: neuron_neuropil", neuropil, fixed = TRUE))
  testthat::expect_true(grepl("spatial_unit: region_layer", neuropil, fixed = TRUE))
  testthat::expect_true(grepl("pgmatrix_imputed_neuron_neuropil_", neuropil, fixed = TRUE))
})

