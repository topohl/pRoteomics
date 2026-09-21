# Adjudication of the pre-existing proposed-scope figure exports.
#
# The accepted manifest predates the canonical-scope policy and still carries
# 16 proposed-scope selections. This file pins what the adjudication
# established, so the conclusion is checkable rather than remembered:
#
#   15 of 16 are byte-identical to an exported canonical counterpart.
#   The 16th, Fig_RES_SUS_divergence_publication.svg, differs from its
#   canonical namesake only in ggrepel label placement: the source-data CSV
#   behind both is byte-identical, the label text is the same 22-element set,
#   and stripping every number from the two SVGs makes them identical
#   line-by-line. The producer calls ggrepel::geom_text_repel() with no
#   set.seed(), so that variation is generation noise, not science.
#
# Assertions are written against relationships rather than a future row count,
# so they survive the eventual exporter rerun.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "export_helpers.R"))

ADJ <- testthat::test_path("..", "..", "audits",
                           "phase6h_proposed_figure_export_adjudication.csv")

testthat::test_that("the adjudication covers exactly the 16 proposed-scope rows", {
  testthat::skip_if_not(file.exists(ADJ), "adjudication audit not present")
  d <- utils::read.csv(ADJ, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(d), 16L)
  testthat::expect_identical(sum(d$byte_identical), 15L)
  testthat::expect_identical(sum(!d$byte_identical), 1L)
  testthat::expect_true(all(d$canonical_exists))
  testthat::expect_true(all(d$canonical_currently_exported))
  testthat::expect_setequal(d$scope, c("compareGO_spatial_atlas_validation_proposed",
                                       "microglia_validation_proposed"))
  testthat::expect_identical(sum(d$scope == "compareGO_spatial_atlas_validation_proposed"), 10L)
  testthat::expect_identical(sum(d$scope == "microglia_validation_proposed"), 6L)
})

testthat::test_that("the one non-identical row is the divergence figure, and only its placement differs", {
  testthat::skip_if_not(file.exists(ADJ), "adjudication audit not present")
  d <- utils::read.csv(ADJ, stringsAsFactors = FALSE)
  distinct <- d[!d$byte_identical, ]
  testthat::expect_identical(nrow(distinct), 1L)
  testthat::expect_identical(basename(distinct$proposed_source),
                             "Fig_RES_SUS_divergence_publication.svg")
  # same size, different bytes: a placement difference, not a content difference
  testthat::expect_identical(distinct$proposed_bytes, distinct$canonical_bytes)

  a <- file.path(repo_path(), distinct$canonical_counterpart)
  b <- file.path(repo_path(), distinct$proposed_source)
  testthat::skip_if_not(file.exists(a) && file.exists(b), "divergence figures not present")
  la <- readLines(a, warn = FALSE); lb <- readLines(b, warn = FALSE)

  # the label text is the same set
  txt <- function(l) gsub("^>|</text>$", "",
                          unlist(regmatches(l, gregexpr(">[^<>]+</text>", l))))
  testthat::expect_setequal(txt(la), txt(lb))
  testthat::expect_identical(length(txt(la)), length(txt(lb)))

  # and with every number removed the two files are identical, so nothing but
  # coordinates differs
  strip <- function(l) gsub("[-0-9]+\\.?[0-9]*", "N", l)
  testthat::expect_identical(strip(la), strip(lb))
})

testthat::test_that("the upstream source data behind both variants is byte-identical", {
  a <- repo_path("results", "source_data", "04_differential_expression_enrichment",
                 "compareGO_spatial_atlas",
                 "source_data_RES_SUS_divergence_publication.csv")
  b <- repo_path("results", "source_data", "04_differential_expression_enrichment",
                 "compareGO_spatial_atlas_validation_proposed",
                 "source_data_RES_SUS_divergence_publication.csv")
  testthat::skip_if_not(file.exists(a) && file.exists(b), "source data not present")
  testthat::expect_identical(unname(tools::sha256sum(a)), unname(tools::sha256sum(b)))
  da <- utils::read.csv(a, stringsAsFactors = FALSE)
  db <- utils::read.csv(b, stringsAsFactors = FALSE)
  testthat::expect_identical(dim(da), dim(db))
  testthat::expect_identical(da, db)
})

testthat::test_that("the producer's label placement is non-deterministic, which is why the two differ", {
  src <- readLines(testthat::test_path("..", "..", "analysis", "differential_abundance",
                                       "build_go_program_atlas.R"), warn = FALSE)
  # the same script writes both variants; only the VALIDATION_ONLY flag differs
  testthat::expect_true(any(grepl("SUBSTEP_ID <- if (VALIDATION_ONLY)", src, fixed = TRUE)))
  testthat::expect_true(any(grepl("geom_text_repel", src, fixed = TRUE)))
  # no seed is set, so repel placement varies between runs. This test documents
  # the cause; it is not a requirement that the producer stay unseeded.
  testthat::expect_identical(length(grep("set\\.seed", src)), 0L)
})

testthat::test_that("no proposed-scope figure is selected under the current policy", {
  d <- if (file.exists(ADJ)) utils::read.csv(ADJ, stringsAsFactors = FALSE) else NULL
  testthat::skip_if(is.null(d), "adjudication audit not present")
  prop <- file.path(repo_path(), d$proposed_source)
  testthat::skip_if(!all(file.exists(prop)), "proposed figures not present")
  kept <- drop_noncanonical_proposed_scopes(prop)
  testthat::expect_length(kept, 0L)

  # and every canonical counterpart survives the same filter
  canon <- file.path(repo_path(), d$canonical_counterpart)
  testthat::expect_setequal(drop_noncanonical_proposed_scopes(canon), canon)
})

testthat::test_that("the accepted export still holds the 16 rows this turn did not touch", {
  mp <- path_results("manuscript", "figure_export_manifest.csv")
  testthat::skip_if_not(file.exists(mp), "figure export manifest not present")
  m <- utils::read.csv(mp, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(m), 5598L)
  testthat::expect_identical(
    unname(tools::sha256sum(mp)),
    "68850f547e5ed71ce1223051761c3a636bc4922037c67248a7bd25cd0b42ad28")
  # the adjudication is a record, not an export change
  src <- gsub("\\\\", "/", m$source_file)
  testthat::expect_identical(sum(grepl("_validation_proposed/", src)), 16L)
})
