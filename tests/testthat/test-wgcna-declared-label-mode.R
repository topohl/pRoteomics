source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::local_edition(3)

# Unit tests for the declared-label-mode helpers in helper-wgcna-current-artifacts.R.
# Fixtures only: manifest and registry paths are passed explicitly.

write_label_manifest <- function(dir, mode, status, source) {
  path <- file.path(dir, "label_review_status.csv")
  utils::write.csv(data.frame(dataset = "microglia", label_mode = mode,
                              label_review_status = status, label_source = source),
                   path, row.names = FALSE)
  path
}

testthat::test_that("a declared provisional mode skips the reviewed-registry requirement", {
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = NA)
  dir <- withr::local_tempdir()
  decl <- wgcna_declared_label_mode("microglia", write_label_manifest(dir, "provisional", "provisional_unreviewed", "p.csv"))
  testthat::expect_identical(decl$mode, "provisional")
  testthat::expect_condition(
    require_reviewed_registry("microglia", decl, file.path(dir, "missing.csv")),
    class = "skip")
})

testthat::test_that("a declared reviewed mode with a missing registry is an error, not a skip", {
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = NA)
  dir <- withr::local_tempdir()
  decl <- wgcna_declared_label_mode("microglia", write_label_manifest(dir, "reviewed", "reviewed", "config/wgcna_labels/microglia.csv"))
  testthat::expect_error(
    require_reviewed_registry("microglia", decl, file.path(dir, "missing.csv")),
    "reviewed registry is missing")
})

testthat::test_that("without a Stage 07 manifest the default is reviewed", {
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = NA)
  dir <- withr::local_tempdir()
  decl <- wgcna_declared_label_mode("microglia", file.path(dir, "absent.csv"))
  testthat::expect_identical(decl$mode, "reviewed")
  registry <- file.path(dir, "microglia.csv"); writeLines("dataset", registry)
  testthat::expect_identical(require_reviewed_registry("microglia", decl, registry), registry)
})

testthat::test_that("an environment value that contradicts the manifest is an error", {
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = "reviewed")
  dir <- withr::local_tempdir()
  testthat::expect_error(
    wgcna_declared_label_mode("microglia", write_label_manifest(dir, "provisional", "provisional_unreviewed", "p.csv")),
    "contradicts")
})

testthat::test_that("a reviewed row needs a completed human review", {
  r <- data.frame(level = "supermodule", entity_id = c("SM01", "SM02", "SM03"),
                  adjudication_status = "reviewed",
                  reviewer = c("A Human", "pending review (mechanical carry-over)", "A Human"),
                  manual_review_required = c(FALSE, FALSE, TRUE))
  testthat::expect_identical(unname(reviewed_rows_are_complete(r)), c(TRUE, FALSE, FALSE))
})
