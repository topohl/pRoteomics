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
  testthat::expect_true(grepl('file.path\\("EWCE_E9", EWCE_DATASET\\)', txt))
  testthat::expect_true(grepl('file.path\\("EWCE_E9_comparison", EWCE_BRANCH, EWCE_DATASET\\)', txt))
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

testthat::test_that("EWCE animal mode is isolated and reports the shared aggregation contract", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  old_unit <- Sys.getenv("PROTEOMICS_EWCE_ANALYSIS_UNIT", unset = NA_character_)
  old_branch <- Sys.getenv("PROTEOMICS_EWCE_BRANCH", unset = NA_character_)
  on.exit({
    if (is.na(old_unit)) Sys.unsetenv("PROTEOMICS_EWCE_ANALYSIS_UNIT") else Sys.setenv(PROTEOMICS_EWCE_ANALYSIS_UNIT = old_unit)
    if (is.na(old_branch)) Sys.unsetenv("PROTEOMICS_EWCE_BRANCH") else Sys.setenv(PROTEOMICS_EWCE_BRANCH = old_branch)
  }, add = TRUE)
  Sys.setenv(PROTEOMICS_EWCE_ANALYSIS_UNIT = "animal", PROTEOMICS_EWCE_BRANCH = "animal_level")
  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("05_celltype_enrichment_EWCE/01_EWCE_E9.r", "--dataset", "neuron_soma", "--dry-run"),
    stdout = TRUE,
    stderr = TRUE
  ))
  out <- paste(out, collapse = "\n")
  testthat::expect_true(grepl("Analysis unit: animal", out, fixed = TRUE))
  testthat::expect_true(grepl("Animal-level column count: 36", out, fixed = TRUE))
  testthat::expect_true(grepl("Legacy cache reuse allowed: FALSE", out, fixed = TRUE))
  testthat::expect_true(grepl("EWCE_E9_comparison/animal_level/neuron_soma", out, fixed = TRUE))
  testthat::expect_false(grepl("EWCE_E9/neuron_soma", out, fixed = TRUE))
})

testthat::test_that("EWCE animal-level safeguards and cache identity are present", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  script <- repo_path("05_celltype_enrichment_EWCE/01_EWCE_E9.r")
  txt <- paste(readLines(script, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("protigy_prepare_animal_level", txt, fixed = TRUE))
  testthat::expect_true(grepl("animal_bundle$output_metadata", txt, fixed = TRUE))
  testthat::expect_false(grepl("animal_bundle$animal_expression", txt, fixed = TRUE))
  testthat::expect_true(grepl("sample_expr_mat <- make_expr_mat(mapped_clean_df, hemisphere_sample_cols)", txt, fixed = TRUE))
  testthat::expect_true(grepl("protigy_aggregate_expression_columns", txt, fixed = TRUE))
  testthat::expect_true(grepl("identical(rownames(expr_mat), rownames(sample_expr_mat))", txt, fixed = TRUE))
  testthat::expect_true(grepl("metadata_rows_in_gene_matrix", txt, fixed = TRUE))
  testthat::expect_true(grepl("microglia_same_gene_logFC", txt, fixed = TRUE))
  testthat::expect_true(grepl("microglia_baseline_top", txt, fixed = TRUE))
  testthat::expect_true(grepl("soma_A111_CA1_single_observed_hemisphere", txt, fixed = TRUE))
  testthat::expect_true(grepl("group_by(Gene, Stratum, Cond)", txt, fixed = TRUE))
  testthat::expect_true(grepl("input_matrix_sha256", txt, fixed = TRUE))
  testthat::expect_true(grepl("hit_gene_set_sha256", txt, fixed = TRUE))
  testthat::expect_true(grepl("background_gene_set_sha256", txt, fixed = TRUE))
  testthat::expect_true(grepl("ctd_annotation_gene_set_sha256", txt, fixed = TRUE))
  testthat::expect_true(grepl('identical\\(EWCE_ANALYSIS_UNIT, "sample"\\)', txt))
  testthat::expect_true(grepl("PROTEOMICS_EWCE_BRANCH is required", txt, fixed = TRUE))
})

testthat::test_that("EWCE signature-only mode exits before parallel or bootstrap execution", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  txt <- paste(
    readLines(repo_path("05_celltype_enrichment_EWCE/01_EWCE_E9.r"), warn = FALSE),
    collapse = "\n"
  )
  testthat::expect_true(grepl("PROTEOMICS_EWCE_SIGNATURE_ONLY", txt, fixed = TRUE))
  signature_exit <- regexpr("if (EWCE_SIGNATURE_ONLY)", txt, fixed = TRUE)[1]
  future_plan <- regexpr("future::plan(future::multisession", txt, fixed = TRUE)[1]
  bootstrap_grid <- regexpr("target_grid <- target_manifest", txt, fixed = TRUE)[1]
  testthat::expect_gt(signature_exit, 0)
  testthat::expect_gt(future_plan, signature_exit)
  testthat::expect_gt(bootstrap_grid, signature_exit)
})

testthat::test_that("EWCE rejects an invalid comparison branch before analysis", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  old_branch <- Sys.getenv("PROTEOMICS_EWCE_BRANCH", unset = NA_character_)
  on.exit({
    if (is.na(old_branch)) Sys.unsetenv("PROTEOMICS_EWCE_BRANCH") else Sys.setenv(PROTEOMICS_EWCE_BRANCH = old_branch)
  }, add = TRUE)
  Sys.setenv(PROTEOMICS_EWCE_BRANCH = "bad branch")
  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("05_celltype_enrichment_EWCE/01_EWCE_E9.r", "--dataset", "neuron_soma", "--dry-run"),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  testthat::expect_true(status != 0L)
  testthat::expect_true(grepl("must match", paste(out, collapse = "\n"), fixed = TRUE))
})
