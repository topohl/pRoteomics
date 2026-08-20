source(testthat::test_path("..", "..", "R", "manuscript_figure3_utils.R"))

testthat::test_that("Figure 3 gene labels accept each documented symbol column", {
  for (column in figure3_gene_symbol_columns()) {
    x <- data.frame(ProteinGroupID = c("PG1", "PG2"), stringsAsFactors = FALSE)
    x[[column]] <- c("GeneA", "")
    testthat::expect_identical(
      figure3_display_gene_symbol(x), c("GeneA", "PG2"),
      info = column
    )
  }
})

testthat::test_that("Figure 3 gene labels use documented priority and ProteinGroupID fallback", {
  x <- data.frame(
    ProteinGroupID = c("PG1", "PG2", "PG3"),
    official_gene_symbol = c("Official", NA, ""),
    representative_gene_symbol = c("Representative", "Representative2", NA),
    gene_symbol = c("Legacy", "Legacy2", "Legacy3"),
    stringsAsFactors = FALSE
  )
  testthat::expect_identical(
    figure3_display_gene_symbol(x),
    c("Official", "Representative2", "Legacy3")
  )
  testthat::expect_identical(
    figure3_display_gene_symbol(data.frame(ProteinGroupID = "PG-only")),
    "PG-only"
  )
})

testthat::test_that("Figure 3 validation status is evaluated and structural failures stop", {
  pass <- figure3_validation_record(
    "unique_rows", "structural_contract", TRUE, TRUE, TRUE, "stop"
  )
  fail <- figure3_validation_record(
    "unique_rows", "structural_contract", FALSE, TRUE, FALSE, "stop"
  )
  observed <- figure3_validation_record(
    "support_count", "observed_manuscript_result", 0L,
    "observed_value_reported", NA, "report_only"
  )
  testthat::expect_identical(pass$status, "pass")
  testthat::expect_identical(fail$status, "fail")
  testthat::expect_identical(observed$status, "observed")
  testthat::expect_silent(figure3_assert_structural_checks(pass))
  testthat::expect_error(
    figure3_assert_structural_checks(fail), "unique_rows"
  )
})
