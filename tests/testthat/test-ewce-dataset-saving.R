testthat::test_that("EWCE script parses --dataset before output paths are created", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  txt <- paste(readLines(repo_path("05_celltype_enrichment_EWCE/01_EWCE_E9.r"), warn = FALSE), collapse = "\n")

  testthat::expect_true(grepl("args <- commandArgs\\(trailingOnly = TRUE\\)", txt))
  testthat::expect_true(grepl("dataset_cli <- arg_value\\(\"--dataset\"", txt))
  testthat::expect_true(grepl("Sys.setenv\\(PROTEOMICS_DATASET = validate_dataset\\(dataset_cli, source = \"--dataset\"\\)\\)", txt))
  testthat::expect_true(grepl("infer_dataset_from_path <- function\\(path\\)", txt))

  pos_dataset_cli <- regexpr("dataset_cli <-", txt, fixed = TRUE)[1]
  pos_current_dataset <- regexpr("EWCE_DATASET <- current_dataset()", txt, fixed = TRUE)[1]
  pos_create_dirs <- regexpr("CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)", txt, fixed = TRUE)[1]
  testthat::expect_gt(pos_current_dataset, pos_dataset_cli)
  testthat::expect_gt(pos_create_dirs, pos_current_dataset)
  testthat::expect_true(grepl('SUBSTEP_ID <- file.path\\("EWCE_E9", EWCE_DATASET\\)', txt))
})

testthat::test_that("EWCE dry-run honors dataset-specific output folders", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  run <- function(dataset) {
    old_wd <- setwd(repo_path())
    on.exit(setwd(old_wd), add = TRUE)
    out <- suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      c("05_celltype_enrichment_EWCE/01_EWCE_E9.r", "--dataset", dataset, "--dry-run"),
      stdout = TRUE,
      stderr = TRUE
    ))
    paste(out, collapse = "\n")
  }

  microglia <- run("microglia")
  testthat::expect_true(grepl("Dataset: microglia", microglia, fixed = TRUE))
  testthat::expect_true(grepl("pgmatrix_imputed_microglia_", microglia, fixed = TRUE))
  testthat::expect_true(grepl("EWCE_E9/microglia", microglia, fixed = TRUE))
  testthat::expect_false(grepl("EWCE_E9/neuron_neuropil", microglia, fixed = TRUE))

  soma <- run("neuron_soma")
  testthat::expect_true(grepl("Dataset: neuron_soma", soma, fixed = TRUE))
  testthat::expect_true(grepl("pgmatrix_imputed_neuron_soma_", soma, fixed = TRUE))
  testthat::expect_true(grepl("EWCE_E9/neuron_soma", soma, fixed = TRUE))
})
