testthat::local_edition(3)

source(testthat::test_path(
  "..", "..", "R", "gsea_wgcna_concordance_utils.R"
))

base_fixture <- function() {
  data.frame(
    dataset = "microglia",
    GO_ID = c("GO:0002181", "GO:0006364", "GO:0006412"),
    NES = c(1.8, -1.2, 1.4),
    stringsAsFactors = FALSE
  )
}

CMP <- c("ca1_sus_res", "ca2_sus_res", "ca3_sus_res")

testthat::test_that("the accepted comparison fields are an exact two-name whitelist", {
  testthat::expect_identical(
    gww_comparison_source_fields(), c("source_comparison", "comparison")
  )
})

testthat::test_that("source_comparison only: accepted, preserved, normalized", {
  x <- base_fixture()
  x$source_comparison <- CMP
  out <- gww_resolve_comparison_field(x)
  testthat::expect_identical(out$comparison, CMP)
  # the raw provenance field survives untouched
  testthat::expect_identical(out$source_comparison, CMP)
})

testthat::test_that("comparison only: accepted for compatibility and preserved", {
  x <- base_fixture()
  x$comparison <- CMP
  out <- gww_resolve_comparison_field(x)
  testthat::expect_identical(out$comparison, CMP)
  testthat::expect_false("source_comparison" %in% names(out))
})

testthat::test_that("both present and identical: accepted, both preserved", {
  x <- base_fixture()
  x$source_comparison <- CMP
  x$comparison <- CMP
  out <- gww_resolve_comparison_field(x)
  testthat::expect_identical(out$comparison, CMP)
  testthat::expect_identical(out$source_comparison, CMP)
})

testthat::test_that("both present but disagreeing on a non-NA row: fails closed", {
  x <- base_fixture()
  x$source_comparison <- CMP
  x$comparison <- c(CMP[1], "ca2_res_con", CMP[3])
  testthat::expect_error(
    gww_resolve_comparison_field(x, "fixture"),
    "conflicting comparison identity"
  )
  testthat::expect_error(
    gww_resolve_comparison_field(x, "fixture"), "ca2_sus_res"
  )
})

testthat::test_that("disagreement only on NA rows is tolerated", {
  x <- base_fixture()
  x$source_comparison <- c(CMP[1], NA_character_, CMP[3])
  x$comparison <- c(CMP[1], "ca2_res_con", CMP[3])
  out <- gww_resolve_comparison_field(x)
  # source_comparison stays canonical, including its NA
  testthat::expect_identical(out$comparison, c(CMP[1], NA_character_, CMP[3]))
})

testthat::test_that("neither field present: fails closed", {
  testthat::expect_error(
    gww_resolve_comparison_field(base_fixture(), "fixture"),
    "missing a comparison identity column"
  )
})

testthat::test_that("resolution is exact-name, never loose or case-insensitive", {
  for (wrong in c("Comparison", "legacy_Comparison", "comparisons",
                  "source_Comparison", "COMPARISON")) {
    x <- base_fixture()
    x[[wrong]] <- CMP
    testthat::expect_error(
      gww_resolve_comparison_field(x, "fixture"),
      "missing a comparison identity column",
      info = wrong
    )
  }
})

testthat::test_that("source_comparison wins when both are present and agree", {
  x <- base_fixture()
  x$source_comparison <- CMP
  x$comparison <- CMP
  out <- gww_resolve_comparison_field(x)
  x2 <- x
  x2$source_comparison[2] <- NA_character_
  out2 <- gww_resolve_comparison_field(x2)
  testthat::expect_identical(out2$comparison[2], NA_character_)
  testthat::expect_identical(out$comparison[2], CMP[2])
})
