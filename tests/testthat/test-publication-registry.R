# Guards for the canonical publication registry and the Extended Data
# canonicalisation performed in Phase 6A.
#
# The point of this phase was to remove every remaining ambiguity about WHICH
# artefact is the publication before the repository is physically restructured.
# These tests encode that: one canonical identity per figure, no competing
# generation, and the three withheld Extended Data figures recorded rather than
# quietly dropped.

repo_rel <- function(...) file.path(testthat::test_path("..", ".."), ...)
pub_contract <- function() {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  yaml::read_yaml(repo_rel("figures", "figure_contract.yml"))
}

testthat::test_that("every numbered publication has exactly one canonical identity", {
  testthat::skip_if_not_installed("yaml")
  y <- pub_contract()
  numbered <- Filter(function(f) isTRUE(f$is_numbered_manuscript_figure), y$figures)

  ids <- vapply(numbered, function(f) {
    as.character(f$canonical_publication_id %||%
      if (!is.null(f$extended_data_number))
        sprintf("extended_data_%02d", as.integer(f$extended_data_number))
      else NA_character_)
  }, character(1))

  testthat::expect_false(anyDuplicated(ids) > 0L)
  testthat::expect_setequal(ids, c(
    "figure_01", "figure_02", "figure_03",
    "extended_data_01", "extended_data_03", "extended_data_05",
    "extended_data_06", "extended_data_08", "extended_data_09"))
})

testthat::test_that("no panel id is reused across publications", {
  testthat::skip_if_not_installed("yaml")
  y <- pub_contract()
  all_ids <- unlist(lapply(y$figures, function(f)
    vapply(f$panels, function(p) as.character(p$id), character(1))), use.names = FALSE)
  dup <- unique(all_ids[duplicated(all_ids)])
  # Extended Data ids are zero-padded precisely so ED1 panel a cannot be confused
  # with Figure 1 panel a. Before that padding both were "1a".
  testthat::expect_identical(dup, character(0))
})

testthat::test_that("the registry matches the contract and records the withheld figures", {
  reg <- repo_rel("manuscript", "canonical_publication_registry.csv")
  testthat::skip_if_not(file.exists(reg), "registry not generated")
  r <- utils::read.csv(reg, stringsAsFactors = FALSE)

  testthat::expect_identical(sum(r$status == "CANONICAL"), 9L)
  withheld <- r[r$status != "CANONICAL", ]
  testthat::expect_setequal(withheld$publication_id,
    c("extended_data_02", "extended_data_04", "extended_data_07"))
  # A withheld identity is RESERVED, not free: nothing may quietly take it.
  testthat::expect_true(all(is.na(withheld$hash) | !nzchar(withheld$hash)))
  testthat::expect_true(all(grepl("^WITHHELD_", withheld$status)))

  y <- pub_contract()
  numbered <- Filter(function(f) isTRUE(f$is_numbered_manuscript_figure), y$figures)
  testthat::expect_identical(length(numbered), 9L)
})

testthat::test_that("the Extended Data promotion audit adjudicates all nine identities", {
  aud <- repo_rel("manuscript", "extended_data_promotion_audit.csv")
  testthat::skip_if_not(file.exists(aud), "promotion audit not generated")
  a <- utils::read.csv(aud, stringsAsFactors = FALSE)

  testthat::expect_identical(nrow(a), 9L)
  testthat::expect_true(all(a$promote %in% c("YES", "NO")))
  testthat::expect_identical(sum(a$promote == "YES"), 6L)
  testthat::expect_setequal(a$canonical_publication_id[a$promote == "NO"],
    c("extended_data_02", "extended_data_04", "extended_data_07"))
  # Every withheld figure must carry a stated defect, never a bare NO.
  for (i in which(a$promote == "NO")) {
    testthat::expect_true(nchar(a$defect_or_note[i]) > 120L,
      info = paste("withheld without a stated defect:", a$ed_figure[i]))
  }
})

testthat::test_that("the promoted Extended Data figures resolve to real assets", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  testthat::skip_if_not_installed("yaml")
  y <- pub_contract()
  for (k in c("ED_01", "ED_03", "ED_06", "ED_08",
              "ED_behaviour_coverage", "ED_behaviour_secondary")) {
    f <- y$figures[[k]]
    testthat::expect_false(isTRUE(f$rendering_repository_computes_statistics))
    for (p in f$panels) {
      testthat::expect_true(file.exists(repo_path(as.character(p$figure_source))),
                            info = paste("missing panel asset:", p$id))
      testthat::expect_true(file.exists(repo_path(as.character(p$primary_source))),
                            info = paste("missing source data:", p$id))
    }
  }
})

testthat::test_that("the pre-restructure freeze manifest is complete", {
  man <- repo_rel("manuscript", "prerestructure_freeze_manifest.csv")
  testthat::skip_if_not(file.exists(man), "freeze manifest not generated")
  m <- utils::read.csv(man, stringsAsFactors = FALSE)

  # The manifest is an equivalence oracle for a structural migration. An object
  # with no hash cannot serve as one.
  testthat::expect_true(all(m$exists))
  testthat::expect_false(any(is.na(m$sha256) | !nzchar(m$sha256)))
  for (cls in c("configuration_contract", "manuscript_provenance",
                "frozen_upstream_import", "canonical_figure_panel",
                "publication_source_data", "canonical_analysis_table",
                "protected_scientific_state", "guard_test")) {
    testthat::expect_gt(sum(m$object_class == cls), 0L,
                        label = paste("objects of class", cls))
  }
  # Every canonical publication must be represented.
  reg <- repo_rel("manuscript", "canonical_publication_registry.csv")
  testthat::skip_if_not(file.exists(reg), "registry not generated")
  r <- utils::read.csv(reg, stringsAsFactors = FALSE)
  canon <- r$publication_id[r$status == "CANONICAL"]
  testthat::expect_true(all(canon %in% m$publication_id[!is.na(m$publication_id)]))
})
