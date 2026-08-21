source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "export_helpers.R"))

testthat::test_that("canonical EWCE figure paths are exportable and comparison paths are not", {
  canonical <- file.path(
    tempdir(), "results", "figures", "05_celltype_enrichment_EWCE", "EWCE_E9",
    "microglia", "Fig1_EWCE_Summary.pdf"
  )
  comparison <- file.path(
    tempdir(), "results", "figures", "05_celltype_enrichment_EWCE", "EWCE_E9_comparison",
    "sample_level_sensitivity", "microglia", "Fig1_EWCE_Summary.pdf"
  )

  testthat::expect_true(is_exportable_result_path(canonical))
  testthat::expect_false(is_exportable_result_path(comparison))
  testthat::expect_match(canonical_ewce_figure_root(), "05_celltype_enrichment_EWCE/EWCE_E9$", fixed = FALSE)
})

testthat::test_that("source-data and table filters retain canonical EWCE only", {
  canonical <- file.path(
    tempdir(), "results", "tables", "05_celltype_enrichment_EWCE", "EWCE_E9",
    "microglia", "Supplementary_Table_EWCE.xlsx"
  )
  comparison <- file.path(
    tempdir(), "results", "source_data", "05_celltype_enrichment_EWCE", "EWCE_E9_comparison",
    "sample_level_sensitivity", "microglia", "Source_Data_EWCE_Figures.xlsx"
  )

  testthat::expect_identical(
    is_exportable_result_path(c(canonical, comparison)),
    c(TRUE, FALSE)
  )
})

testthat::test_that("supplementary EWCE glob selects canonical tables, not comparison branches", {
  root <- file.path(tempdir(), paste0("ewce_export_contract_", as.integer(stats::runif(1, 1, 1e9))))
  canonical_dir <- file.path(root, "results", "tables", "05_celltype_enrichment_EWCE", "EWCE_E9", "microglia")
  comparison_dir <- file.path(root, "results", "tables", "05_celltype_enrichment_EWCE", "EWCE_E9_comparison", "branch", "microglia")
  dir.create(canonical_dir, recursive = TRUE)
  dir.create(comparison_dir, recursive = TRUE)
  canonical_file <- file.path(canonical_dir, "Supplementary_Table_EWCE.xlsx")
  comparison_file <- file.path(comparison_dir, "Supplementary_Table_EWCE.xlsx")
  file.create(canonical_file)
  file.create(comparison_file)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  config <- list(supplementary_table_globs = file.path(
    root, "results", "tables", "05_celltype_enrichment_EWCE", "EWCE_E9", "*",
    "Supplementary_Table_EWCE.xlsx"
  ))
  selected <- supplementary_candidate_files(config, datasets = "microglia")

  testthat::expect_identical(normalize_export_path(selected), normalize_export_path(canonical_file))
  testthat::expect_false(normalize_export_path(comparison_file) %in% normalize_export_path(selected))
})

testthat::test_that("export scripts and pipeline declare the canonical-only EWCE contract", {
  figure_script <- paste(readLines(repo_path("09_export_pride_journal", "08_export_manuscript_figures.R")), collapse = "\n")
  source_data_script <- paste(readLines(repo_path("09_export_pride_journal", "09_export_source_data.R")), collapse = "\n")
  targeted_signature_script <- paste(readLines(repo_path("04_differential_expression_enrichment", "05_microglia_targeted_signature_enrichment.r")), collapse = "\n")
  pipeline <- paste(readLines(repo_path("pipeline.yml")), collapse = "\n")

  testthat::expect_match(figure_script, "canonical_ewce_figure_root\\(\\)")
  testthat::expect_match(source_data_script, "is_exportable_result_path\\(candidates\\)")
  testthat::expect_match(targeted_signature_script, "Diagnostic only")
  testthat::expect_match(pipeline, "EWCE_results_full\\.rds")
})
