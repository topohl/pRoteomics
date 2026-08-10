source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "manuscript_go_theme_utils.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "plotting_nature.R"))
source(testthat::test_path("..", "..", "R", "sus_res_biological_audit_workbook.R"))

workbook_path <- repo_path(
  "results", "reports", "04_differential_expression_enrichment",
  "sus_res_spatial_dap_atlas", "global", "sus_res_biological_audit.xlsx"
)
detail_path <- repo_path(
  "results", "source_data", "04_differential_expression_enrichment",
  "sus_res_spatial_dap_atlas", "global", "sus_res_manuscript_theme_summary.csv"
)

testthat::test_that("human-facing SUS-RES workbook has the six-sheet audit contract", {
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not(file.exists(workbook_path), "Generated SUS-RES workbook is unavailable")
  testthat::expect_identical(
    readxl::excel_sheets(workbook_path),
    c("Overview", "Spatial Map", "Theme Detail", "GO Term Audit", "DAP Summary", "Provenance")
  )
  overview <- readxl::read_excel(workbook_path, sheet = "Overview", skip = 3)
  spatial <- readxl::read_excel(workbook_path, sheet = "Spatial Map", skip = 3)
  detail <- readxl::read_excel(workbook_path, sheet = "Theme Detail", skip = 3)
  go_audit <- readxl::read_excel(workbook_path, sheet = "GO Term Audit", skip = 3)
  dap <- readxl::read_excel(workbook_path, sheet = "DAP Summary", skip = 3)
  testthat::expect_equal(nrow(overview), 8L)
  testthat::expect_equal(nrow(spatial), 8L)
  testthat::expect_equal(nrow(detail), 144L)
  testthat::expect_equal(nrow(go_audit), 336L)
  testthat::expect_equal(nrow(dap), 18L)
  testthat::expect_true(all(c(
    "Median NES (all theme terms)", "FDR support present", "Representative supported GO ID",
    "Representative supported BH FDR", "FDR-supported DAP count"
  ) %in% names(detail)))
  testthat::expect_true(all(c("Original NES", "Original BH FDR", "Assignment status") %in% names(go_audit)))
})

testthat::test_that("workbook descriptive values match the canonical CSV and contain no spreadsheet errors", {
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not(file.exists(workbook_path) && file.exists(detail_path), "Generated SUS-RES workbook inputs are unavailable")
  detail_source <- utils::read.csv(detail_path, check.names = FALSE, stringsAsFactors = FALSE)
  validation <- validate_sus_res_biological_audit_workbook(workbook_path, detail_source)
  testthat::expect_true(all(validation$validation_status == "passed"))
  testthat::expect_identical(validation$data_rows[validation$sheet == "Theme Detail"], 144L)
  testthat::expect_identical(validation$data_rows[validation$sheet == "GO Term Audit"], 336L)
})
