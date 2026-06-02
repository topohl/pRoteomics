testthat::test_that("module-score metadata merge script is dataset-aware", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  txt <- paste(readLines(repo_path("01_preprocessing/06_merged_metadata_module_score.r"), warn = FALSE), collapse = "\n")

  testthat::expect_true(grepl("--dataset", txt, fixed = TRUE))
  testthat::expect_true(grepl("PROTEOMICS_DATASET", txt, fixed = TRUE))
  testthat::expect_true(grepl("resolve_dataset_inputs(dataset_profile, purpose = \"module_score\")", txt, fixed = TRUE))
  testthat::expect_false(grepl("expected_name <- \"20260218_pgmatrix_imputed_neuron_neuropil", txt, fixed = TRUE))
  testthat::expect_true(grepl("results/module_scores/<dataset>", txt, fixed = TRUE) || grepl("path_results(\"module_scores\", dataset_profile", txt, fixed = TRUE))
})

testthat::test_that("dataset input resolution prefers dataset-scoped module metadata", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  txt <- paste(readLines(repo_path("R/dataset_inputs.R"), warn = FALSE), collapse = "\n")

  testthat::expect_true(grepl("path_results(\"module_scores\", dataset, \"sample_metadata_merged_clean_for_module_scores.xlsx\")", txt, fixed = TRUE))
  testthat::expect_true(grepl("path_processed(\"01_preprocessing\", dataset, \"sample_metadata_merged_clean_for_module_scores.xlsx\")", txt, fixed = TRUE))
  testthat::expect_true(grepl("PROTEOMICS_ALLOW_GLOBAL_MODULE_SCORE_METADATA", txt, fixed = TRUE))
  testthat::expect_true(grepl("legacy global module-score metadata fallback", txt, fixed = TRUE))
})

testthat::test_that("module scoring writes overlap diagnostics and preserves spatial labels", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  txt <- paste(readLines(repo_path("06_modules_WGCNA/03_score_module_activity.R"), warn = FALSE), collapse = "\n")

  for (needle in c(
    "module_score_sample_overlap_diagnostics.csv",
    "protein_matrix_sample_columns.csv",
    "metadata_sample_ids.csv",
    "No matching sample names between protein matrix and metadata. This usually means the module-score metadata workbook was generated for another dataset.",
    "Rscript 01_preprocessing/06_merged_metadata_module_score.r --dataset <dataset>",
    "SpatialUnit",
    "SpatialLabel"
  )) {
    testthat::expect_true(grepl(needle, txt, fixed = TRUE), info = needle)
  }

  testthat::expect_true(grepl("left_join(metadata, by = \"Sample\")", txt, fixed = TRUE))
  testthat::expect_true(grepl("write.xlsx(scores_df", txt, fixed = TRUE))
  testthat::expect_true(grepl("write_csv(scores_df", txt, fixed = TRUE))
})
