testthat::test_that("spatial network producer is dataset and spatial-unit aware", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  txt <- paste(readLines(repo_path("07_spatial_networks/01_network_spatial_relations.r"), warn = FALSE), collapse = "\n")

  testthat::expect_false(grepl('assert_dataset_capability\\(SPATIAL_DATASET, "layer"', txt))
  testthat::expect_true(grepl('assert_dataset_capability\\(SPATIAL_DATASET, "region"', txt))
  testthat::expect_true(grepl('spatial_unit <- if \\(identical\\(SPATIAL_DATASET, "neuron_neuropil"\\)\\) "region_layer" else "region"', txt))
  testthat::expect_true(grepl('spatial_col <- if \\(identical\\(spatial_unit, "region_layer"\\)\\) "RegionLayer" else "Region"', txt))
  testthat::expect_true(grepl("SpatialLabel", txt, fixed = TRUE))
  testthat::expect_true(grepl('SUBSTEP_ID <- file.path\\("network_spatial_relations", SPATIAL_DATASET, spatial_unit\\)', txt))
  testthat::expect_true(grepl("aggregate_spatial_unit_expression", txt, fixed = TRUE))
  testthat::expect_true(grepl("dplyr::group_by\\(Protein, SpatialLabel\\)", txt))
  testthat::expect_true(grepl('SPATIAL_DATASET != "neuron_neuropil" & !is.na\\(Region\\) ~ Region', txt))
  testthat::expect_true(grepl("spatial_unit_matrix", txt, fixed = TRUE))
})

testthat::test_that("pipeline registry advertises scoped spatial network outputs", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  testthat::skip_if_not_installed("yaml")
  registry <- yaml::read_yaml(repo_path("pipeline.yml"))
  scripts <- registry$stages$networks$scripts
  producer <- scripts[[which(vapply(scripts, function(x) x$script, character(1)) == "07_spatial_networks/01_network_spatial_relations.r")]]

  testthat::expect_true("microglia" %in% producer$datasets)
  testthat::expect_true(any(grepl("network_spatial_relations/<dataset>/*/network_spatial_relations_objects.rds", producer$produces, fixed = TRUE)))
  testthat::expect_true(any(grepl("data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv", producer$consumes_required, fixed = TRUE)))
})
