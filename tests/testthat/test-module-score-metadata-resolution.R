source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::test_that("module-score metadata merge script is dataset-aware", {
  txt <- paste(readLines(repo_path("analysis/preprocessing/build_module_score_metadata.R"), warn = FALSE), collapse = "\n")

  testthat::expect_true(grepl("--dataset", txt, fixed = TRUE))
  testthat::expect_true(grepl("PROTEOMICS_DATASET", txt, fixed = TRUE))
  testthat::expect_true(grepl("resolve_dataset_inputs(dataset_profile, purpose = \"module_score\")", txt, fixed = TRUE))
  testthat::expect_false(grepl("expected_name <- \"20260218_pgmatrix_imputed_neuron_neuropil", txt, fixed = TRUE))
  ## The property under test is that the destination is scoped to the resolved
  ## dataset rather than hard-coded. Phase 6G.7 changed how that is spelled:
  ## the scope is now the dataset argument to preprocessing_dirs() instead of a
  ## path_processed() segment, and the header names the normalized namespace.
  ## All three spellings are accepted so this locks the property, not a phase.
  testthat::expect_true(
    grepl("preprocessing_dirs(\"build_module_score_metadata\", dataset_profile)", txt, fixed = TRUE) ||
      grepl("results/preprocessing/build_module_score_metadata/<dataset>", txt, fixed = TRUE) ||
      grepl("data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>", txt, fixed = TRUE) ||
      grepl("path_processed(module_id, substep_id, dataset_profile)", txt, fixed = TRUE))
})

testthat::test_that("dataset input resolution prefers dataset-scoped module metadata", {
  txt <- paste(readLines(repo_path("R/data_contracts/dataset_inputs.R"), warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("path_processed(\n      \"01_preprocessing\",\n      \"06_merged_metadata_module_score\"", txt, fixed = TRUE) || grepl("06_merged_metadata_module_score", txt, fixed = TRUE))
  testthat::expect_true(grepl("legacy_dataset_candidates", txt, fixed = TRUE))
  testthat::expect_true(grepl("PROTEOMICS_ALLOW_GLOBAL_MODULE_SCORE_METADATA", txt, fixed = TRUE))
  testthat::expect_true(grepl("legacy global module-score metadata fallback", txt, fixed = TRUE))
})

testthat::test_that("module scoring writes overlap diagnostics and preserves spatial labels", {
  txt <- paste(readLines(repo_path("analysis/wgcna/score_module_activity.R"), warn = FALSE), collapse = "\n")

  for (needle in c(
    "module_score_sample_overlap_diagnostics.csv",
    "protein_matrix_sample_columns.csv",
    "metadata_sample_ids.csv",
    "No matching sample names between protein matrix and metadata. This usually means the module-score metadata workbook was generated for another dataset.",
    "Rscript analysis/preprocessing/build_module_score_metadata.R --dataset <dataset>",
    "SpatialUnit",
    "SpatialLabel"
  )) {
    testthat::expect_true(grepl(needle, txt, fixed = TRUE), info = needle)
  }

  testthat::expect_true(grepl("left_join(metadata, by = \"Sample\")", txt, fixed = TRUE))
  testthat::expect_true(grepl("write.xlsx(scores_df", txt, fixed = TRUE))
  testthat::expect_true(grepl("write_csv(scores_df", txt, fixed = TRUE))
})
