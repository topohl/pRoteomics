# Publication eligibility for proposed/validation figure scopes.
#
# Until Phase 6H.5B the figure exporter had no notion of a proposed scope.
# Eligibility there was decided by accident: five figures under
# microglia_validation_proposed/ had their extension truncated away by a
# MAX_PATH overrun, so the exporter's \.(svg|pdf|png)$ selector skipped them,
# while their six intact siblings in the same scope were selected. A filesystem
# bug was acting as an implicit publication filter.
#
# The contract these tests pin separates the two concerns:
#   source filename correctness -> valid .svg/.png names, always
#   publication eligibility     -> explicit canonical-scope policy
# so a correctly named figure in a proposed scope is ineligible, and no figure
# depends on filename damage to be excluded.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "export_helpers.R"))

testthat::test_that("the proposed-scope pairs name figure trees exactly", {
  pairs <- source_data_proposed_tree_pairs()
  proposed <- vapply(pairs, `[[`, character(1), 1L)
  testthat::expect_true(any(startsWith(proposed, "results/figures/")))
  # every pair is an exact prefix pair, never a pattern
  for (p in pairs) {
    testthat::expect_length(p, 2L)
    testthat::expect_true(endsWith(p[[1]], "/"))
    testthat::expect_false(grepl("[*?|]", p[[1]]))
  }
  # both known proposed figure trees are covered
  testthat::expect_true(any(grepl(
    "^results/figures/.*microglia_targeted_signature_enrichment/microglia_validation_proposed/$",
    proposed)))
  testthat::expect_true(any(grepl(
    "^results/figures/.*compareGO_spatial_atlas_validation_proposed/$", proposed)))
})

testthat::test_that("a valid-extension figure in a proposed scope is excluded, its canonical twin kept", {
  # The decisive case: both files perfectly named, byte-identical content.
  # Only the scope may decide.
  tmp <- withr::local_tempdir()
  base <- file.path(tmp, "results", "figures", "04_differential_expression_enrichment",
                    "microglia_targeted_signature_enrichment")
  canon <- file.path(base, "microglia", "supplementary_qc")
  prop <- file.path(base, "microglia_validation_proposed", "supplementary_qc")
  dir.create(canon, recursive = TRUE, showWarnings = FALSE)
  dir.create(prop, recursive = TRUE, showWarnings = FALSE)
  body <- '<svg xmlns="http://www.w3.org/2000/svg"></svg>'
  cf <- file.path(canon, "microglia_signature_dotplot.svg")
  pf <- file.path(prop, "microglia_signature_dotplot.svg")
  writeLines(body, cf); writeLines(body, pf)
  # identical bytes, identical valid names
  testthat::expect_identical(unname(tools::sha256sum(cf)), unname(tools::sha256sum(pf)))
  testthat::expect_identical(tools::file_ext(cf), tools::file_ext(pf))

  kept <- drop_noncanonical_proposed_scopes(c(cf, pf),
                                            results_root = file.path(tmp, "results"))
  testthat::expect_true(cf %in% kept)
  testthat::expect_false(pf %in% kept)
  testthat::expect_length(kept, 1L)
})

testthat::test_that("exclusion does not depend on filename damage", {
  tmp <- withr::local_tempdir()
  base <- file.path(tmp, "results", "figures", "04_differential_expression_enrichment",
                    "microglia_targeted_signature_enrichment")
  dir.create(file.path(base, "microglia"), recursive = TRUE, showWarnings = FALSE)
  prop <- file.path(base, "microglia_validation_proposed")
  dir.create(prop, recursive = TRUE, showWarnings = FALSE)
  # a damaged name and a perfect name in the same proposed scope: both excluded
  damaged <- file.path(prop, "microglia_signature_scoreboard.sv")
  perfect <- file.path(prop, "microglia_signature_scoreboard.svg")
  writeLines("x", damaged); writeLines("x", perfect)
  kept <- drop_noncanonical_proposed_scopes(c(damaged, perfect),
                                            results_root = file.path(tmp, "results"))
  testthat::expect_length(kept, 0L)
})

testthat::test_that("the exclusion is exact-scope and spares canonical validation analyses", {
  tmp <- withr::local_tempdir()
  figs <- file.path(tmp, "results", "figures", "04_differential_expression_enrichment")
  # a canonical analysis whose NAME contains "validation"
  ok1 <- file.path(figs, "control_spatial_identity_validation", "global", "p.svg")
  # a canonical domain called spatial_validation
  ok2 <- file.path(tmp, "results", "spatial_validation", "build_spatial_data_contract",
                   "global", "plots", "p.svg")
  for (p in c(ok1, ok2)) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); writeLines("x", p)
  }
  kept <- drop_noncanonical_proposed_scopes(c(ok1, ok2),
                                            results_root = file.path(tmp, "results"))
  testthat::expect_true(ok1 %in% kept)
  testthat::expect_true(ok2 %in% kept)
  testthat::expect_length(kept, 2L)
})

testthat::test_that("the exclusion fails closed when the canonical sibling is absent", {
  tmp <- withr::local_tempdir()
  prop <- file.path(tmp, "results", "figures", "04_differential_expression_enrichment",
                    "microglia_targeted_signature_enrichment", "microglia_validation_proposed")
  dir.create(prop, recursive = TRUE, showWarnings = FALSE)
  p <- file.path(prop, "x.svg"); writeLines("x", p)
  # no canonical microglia/ tree exists, so excluding would be unjustified
  testthat::expect_error(
    drop_noncanonical_proposed_scopes(p, results_root = file.path(tmp, "results")),
    "Refusing to exclude proposed")
})

testthat::test_that("the figure exporter applies the proposed-scope policy", {
  code <- readLines(testthat::test_path("..", "..", "analysis", "publication_source_data",
                                        "08_export_manuscript_figures.R"), warn = FALSE)
  testthat::expect_true(any(grepl("drop_noncanonical_proposed_scopes(candidates)",
                                  code, fixed = TRUE)))
  # and it still applies the failed-run guard
  testthat::expect_true(any(grepl("drop_noncanonical_wgcna_dataset_scopes(candidates)",
                                  code, fixed = TRUE)))
})

testthat::test_that("the validation writer budgets its figure targets", {
  code <- readLines(testthat::test_path("..", "..", "analysis", "differential_abundance",
                                        "test_microglia_targeted_signatures.R"), warn = FALSE)
  testthat::expect_true(any(grepl("budgeted_figure_path(", code, fixed = TRUE)))
  # ggsave_publication routes through ggsave_plot, so one budget covers both
  testthat::expect_true(any(grepl("ggsave_plot(paste0(base_path_no_ext, \".svg\")",
                                  code, fixed = TRUE)))
})

testthat::test_that("svg and png survive budgeting in the deep proposed scope", {
  # the real directory depth that truncated the five
  d <- repo_path("results", "figures", "04_differential_expression_enrichment",
                 "microglia_targeted_signature_enrichment",
                 "microglia_validation_proposed", "main_candidate")
  for (ext in c("svg", "png")) {
    out <- budgeted_figure_target(d, paste0("microglia_vs_neuropil_signature_NES_scatter_main.", ext))
    testthat::expect_identical(tools::file_ext(out), ext)
    testthat::expect_false(endsWith(out, "."))
    for (bad in c(".s", ".sv", ".p", ".pn")) testthat::expect_false(endsWith(out, bad))
    testthat::expect_lte(path_length_chars(file.path(d, out)), FIGURE_WRITE_BUDGET)
  }
})

testthat::test_that("no figure in a proposed scope carries a damaged name any more", {
  root <- repo_path("results", "figures", "04_differential_expression_enrichment",
                    "microglia_targeted_signature_enrichment",
                    "microglia_validation_proposed")
  testthat::skip_if_not(dir.exists(root), "proposed figure scope not present")
  f <- list.files(root, recursive = TRUE)
  testthat::skip_if(!length(f), "proposed figure scope empty")
  testthat::expect_true(all(tolower(tools::file_ext(f)) %in% c("svg", "png", "pdf")))
})

testthat::test_that("the repaired validation figures kept their bytes and their canonical twins", {
  a <- testthat::test_path("..", "..", "audits", "phase6h_validation_figure_repair.csv")
  testthat::skip_if_not(file.exists(a), "repair audit not present")
  d <- utils::read.csv(a, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(d), 5L)
  testthat::expect_true(all(d$same_hash))
  testthat::expect_true(all(d$same_size))
  testthat::expect_true(all(d$mtime_preserved))
  testthat::expect_true(all(d$matches_canonical))
  testthat::expect_setequal(d$intended_extension, c(".svg", ".png"))
  testthat::expect_true(all(d$new_abs_chars < PATH_LENGTH_WALL))
})

testthat::test_that("the figure export carries no proposed-scope selection", {
  mp <- path_results("manuscript", "figure_export_manifest.csv")
  testthat::skip_if_not(file.exists(mp), "figure export manifest not present")
  m <- utils::read.csv(mp, stringsAsFactors = FALSE)
  src <- gsub("\\\\", "/", m$source_file)
  ## Phase 6H.5D applied the policy to the package: the 16 proposed-scope rows
  ## that predated it are gone, so this asserts the invariant rather than a row
  ## count, which is what the policy actually guarantees.
  testthat::expect_identical(sum(grepl("_validation_proposed/", src)), 0L)
  testthat::expect_setequal(src, gsub("\\\\", "/", drop_noncanonical_proposed_scopes(src)))

  # and none of the five repaired validation figures entered the export
  a <- testthat::test_path("..", "..", "audits", "phase6h_validation_figure_repair.csv")
  testthat::skip_if_not(file.exists(a), "repair audit not present")
  d <- utils::read.csv(a, stringsAsFactors = FALSE)
  testthat::expect_false(any(gsub("\\\\", "/", file.path(repo_path(), d$new_path)) %in% src))
})
