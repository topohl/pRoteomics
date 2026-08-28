testthat::local_edition(3)

# Stage 11 consumes the canonical compareGO ranked-GSEA long atlas. That
# producer labels its output "canonical_compareGO_ranked_GSEA_GO_BP"; earlier
# canonical runs and the publication-facing Panel C source label the same
# evidence "ranked_GSEA". Both are one semantic ranked-GSEA class; anything
# else is a different inferential lineage and must fail closed.

# paths.R first: the utils header resolves its own dependencies relatively.
source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "sus_res_spatial_dap_atlas_utils.R"))

CURRENT <- "canonical_compareGO_ranked_GSEA_GO_BP"
HISTORICAL <- "ranked_GSEA"

testthat::test_that("the accepted vocabulary matches the established repository contract", {
  # The list is duplicated in sus_res_spatial_dap_atlas_utils.R to avoid
  # pulling the whole GSEA-WGCNA concordance utility into Stage 10/11. This
  # pins the duplicate to the canonical definition so the two cannot drift.
  concordance <- testthat::test_path(
    "..", "..", "R", "gsea_wgcna_concordance_utils.R"
  )
  testthat::skip_if_not(file.exists(concordance), "Concordance utility unavailable")
  env <- new.env(parent = globalenv())
  for (e in parse(concordance, keep.source = FALSE)) {
    parts <- as.list(e)
    if (length(parts) == 3L && identical(as.character(parts[[1]]), "<-") &&
        is.name(parts[[2]]) &&
        as.character(parts[[2]]) %in% c("gww_ranked_gsea_source_families",
                                        "gww_ranked_gsea_source_class")) {
      assign(as.character(parts[[2]]), eval(parts[[3]], envir = env), envir = env)
    }
  }
  testthat::expect_setequal(
    SUS_RES_RANKED_GSEA_SOURCE_FAMILIES, env$gww_ranked_gsea_source_families()
  )
  testthat::expect_identical(
    SUS_RES_RANKED_GSEA_SOURCE_CLASS, env$gww_ranked_gsea_source_class()
  )
})

testthat::test_that("1/2. both accepted vocabularies are accepted", {
  testthat::expect_silent(sus_res_normalize_ranked_gsea_source_class(CURRENT))
  testthat::expect_silent(sus_res_normalize_ranked_gsea_source_class(HISTORICAL))
})

testthat::test_that("3. both normalize to evidence_source_class ranked_GSEA", {
  testthat::expect_identical(
    sus_res_normalize_ranked_gsea_source_class(c(CURRENT, HISTORICAL)),
    c("ranked_GSEA", "ranked_GSEA")
  )
  testthat::expect_identical(
    sus_res_normalize_ranked_gsea_source_class(NA_character_), NA_character_
  )
})

testthat::test_that("4/10. raw provenance and statistics are untouched by normalization", {
  go <- data.frame(
    GO_ID = c("GO:0002181", "GO:0006397"),
    NES = c(1.84, -1.91),
    pvalue = c(1e-5, 2e-4),
    p.adjust = c(2.6e-4, 2.8e-4),
    core_enrichment = c("A/B", "C/D"),
    evidence_source_family = c(CURRENT, CURRENT),
    stringsAsFactors = FALSE
  )
  before <- go
  go$evidence_source_class <- sus_res_normalize_ranked_gsea_source_class(
    go$evidence_source_family
  )
  # raw provenance verbatim, never rewritten to the historical token
  testthat::expect_identical(go$evidence_source_family, before$evidence_source_family)
  testthat::expect_identical(unique(go$evidence_source_family), CURRENT)
  testthat::expect_identical(unique(go$evidence_source_class), "ranked_GSEA")
  # every statistical / identity column unchanged
  for (cn in c("GO_ID", "NES", "pvalue", "p.adjust", "core_enrichment")) {
    testthat::expect_identical(go[[cn]], before[[cn]], info = cn)
  }
})

testthat::test_that("5/6. unsupported families fail closed and name the offender", {
  for (bad in c("ORA", "over_representation", "KEGG",
                "canonical_compareGO_ORA_KEGG", "ranked_GSEA_KEGG",
                "", "unknown_family")) {
    testthat::expect_error(
      sus_res_normalize_ranked_gsea_source_class(bad, "fixture"),
      "unsupported evidence_source_family", info = bad
    )
  }
  testthat::expect_error(
    sus_res_normalize_ranked_gsea_source_class("ORA", "fixture"), "ORA"
  )
})

testthat::test_that("7. a mixture of the two accepted vocabularies is semantically accepted", {
  mixed <- c(CURRENT, HISTORICAL, CURRENT, HISTORICAL)
  out <- sus_res_normalize_ranked_gsea_source_class(mixed)
  testthat::expect_identical(unique(out), "ranked_GSEA")
  testthat::expect_length(out, 4L)
})

testthat::test_that("8. a mixture containing any unsupported family fails", {
  testthat::expect_error(
    sus_res_normalize_ranked_gsea_source_class(c(CURRENT, HISTORICAL, "ORA")),
    "unsupported evidence_source_family"
  )
})

testthat::test_that("9. the current canonical Stage-11 input passes the contract", {
  go_source <- testthat::test_path(
    "..", "..", "results", "source_data",
    "04_differential_expression_enrichment", "compareGO_spatial_atlas",
    "spatial_atlas_enrichment_long.csv"
  )
  testthat::skip_if_not(file.exists(go_source), "Canonical GO source unavailable")
  fams <- unique(suppressWarnings(readr::read_csv(
    go_source, col_select = "evidence_source_family",
    show_col_types = FALSE, progress = FALSE
  ))$evidence_source_family)
  testthat::expect_true(all(fams %in% SUS_RES_RANKED_GSEA_SOURCE_FAMILIES))
  testthat::expect_identical(
    unique(sus_res_normalize_ranked_gsea_source_class(fams)), "ranked_GSEA"
  )
})

testthat::test_that("the Stage-11 consumer no longer uses a literal-equality check", {
  src <- readLines(testthat::test_path(
    "..", "..", "04_differential_expression_enrichment",
    "11_stress_response_biological_audit.r"
  ), warn = FALSE)
  testthat::expect_false(any(grepl(
    'evidence_source_family) != "ranked_GSEA"', src, fixed = TRUE
  )))
  testthat::expect_false(any(grepl("GO source mixes evidence families", src, fixed = TRUE)))
  testthat::expect_true(any(grepl(
    "sus_res_normalize_ranked_gsea_source_class", src, fixed = TRUE
  )))
})
