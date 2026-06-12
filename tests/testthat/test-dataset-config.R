testthat::test_that("dataset normalization and capabilities are stable", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "dataset_config.R"))

  testthat::expect_identical(normalize_dataset(c("neuropil", "soma", "microglial")), c("neuron_neuropil", "neuron_soma", "microglia"))
  testthat::expect_identical(valid_datasets(), c("neuron_neuropil", "neuron_soma", "microglia"))
  testthat::expect_true(dataset_has_capability("neuron_neuropil", "layer"))
  testthat::expect_false(dataset_has_capability("microglia", "layer"))
  testthat::expect_true(dataset_has_capability("microglia", "celltype_roi"))
  testthat::expect_false(dataset_has_capability("microglia", "purified_celltype"))
  testthat::expect_false(dataset_has_capability("neuron_soma", "layer"))
  testthat::expect_silent(assert_dataset_capability("neuron_soma", "region", analysis = "region-resolved test"))
  testthat::expect_error(
    assert_dataset_capability("neuron_soma", "layer", analysis = "layer-resolved test"),
    "does not support layer-level layer-resolved test"
  )
  testthat::expect_error(
    assert_dataset_capability("microglia", "layer", analysis = "spatial network analysis"),
    "does not support layer-level spatial network analysis"
  )
  testthat::expect_error(dataset_has_capability("microglia", "cell_intrinsic"), "Unknown dataset capability")
})
