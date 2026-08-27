testthat::local_edition(3)

source(testthat::test_path(
  "..", "..", "R", "gsea_wgcna_concordance_utils.R"
))

# Minimal canonical ranked-GSEA term table as gww_validate_gsea_terms() sees it.
vocab_fixture <- function(family) {
  data.frame(
    dataset = "microglia",
    phenotype_contrast = "SUS_vs_RES",
    spatial_unit = c("CA1", "CA1", "CA2"),
    program_class = "Synapse_Vesicle_Organization",
    comparison = c("ca1_a", "ca1_b", "ca2_a"),
    ID = c("GO:0007268", "GO:0006412", "GO:0006119"),
    Description = c("first", "second", "third"),
    NES = c(1.8, -1.2, 1.4),
    pvalue = c(0.01, 0.03, 0.02),
    p.adjust = c(0.01, 0.04, 0.03),
    core_enrichment = c("P1/P2", "P3", "P2"),
    evidence_source_family = family,
    stringsAsFactors = FALSE
  )
}

CURRENT <- "canonical_compareGO_ranked_GSEA_GO_BP"
HISTORICAL <- "ranked_GSEA"

testthat::test_that("the accepted ranked-GSEA vocabulary is exactly two values", {
  testthat::expect_setequal(
    gww_ranked_gsea_source_families(), c(CURRENT, HISTORICAL)
  )
  testthat::expect_identical(gww_ranked_gsea_source_class(), "ranked_GSEA")
})

testthat::test_that("current canonical compareGO vocabulary is accepted", {
  testthat::expect_true(gww_validate_gsea_terms(vocab_fixture(CURRENT)))
})

testthat::test_that("historical ranked_GSEA vocabulary is accepted", {
  testthat::expect_true(gww_validate_gsea_terms(vocab_fixture(HISTORICAL)))
})

testthat::test_that("both vocabularies normalize to one evidence class", {
  testthat::expect_identical(
    gww_normalize_evidence_source_class(c(CURRENT, HISTORICAL)),
    c("ranked_GSEA", "ranked_GSEA")
  )
  testthat::expect_identical(
    gww_normalize_evidence_source_class(NA_character_), NA_character_
  )
})

testthat::test_that("non-ranked evidence families are rejected (fail closed)", {
  for (bad in c(
    "ORA", "over_representation", "KEGG", "canonical_compareGO_ORA_KEGG",
    "ranked_GSEA_KEGG", "", "unknown_family"
  )) {
    testthat::expect_error(
      gww_validate_gsea_terms(vocab_fixture(bad)),
      "unsupported evidence_source_family",
      info = bad
    )
  }
})

testthat::test_that("a mixed table containing one bad family is rejected", {
  mixed <- vocab_fixture(CURRENT)
  mixed$evidence_source_family[2] <- "ORA"
  testthat::expect_error(
    gww_validate_gsea_terms(mixed), "unsupported evidence_source_family"
  )
})

testthat::test_that("provenance is preserved, never rewritten to the class token", {
  x <- vocab_fixture(CURRENT)
  x$evidence_source_class <- gww_normalize_evidence_source_class(
    x$evidence_source_family
  )
  # The provenance field keeps the upstream vocabulary verbatim ...
  testthat::expect_identical(unique(x$evidence_source_family), CURRENT)
  # ... and the normalized class lives in its own explicit field.
  testthat::expect_identical(unique(x$evidence_source_class), "ranked_GSEA")
  testthat::expect_false(identical(
    unique(x$evidence_source_family), unique(x$evidence_source_class)
  ))
})

testthat::test_that("current and historical fixtures agree on everything but provenance", {
  cur <- vocab_fixture(CURRENT)
  hist <- vocab_fixture(HISTORICAL)

  # Every statistic and identity column is untouched by the vocabulary.
  stat_cols <- setdiff(names(cur), "evidence_source_family")
  testthat::expect_identical(cur[stat_cols], hist[stat_cols])

  # Both normalize to the same class, so downstream class-keyed logic is
  # identical while the provenance vocabulary stays distinct.
  testthat::expect_identical(
    gww_normalize_evidence_source_class(cur$evidence_source_family),
    gww_normalize_evidence_source_class(hist$evidence_source_family)
  )
  testthat::expect_false(identical(
    cur$evidence_source_family, hist$evidence_source_family
  ))
})
