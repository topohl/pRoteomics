testthat::test_that("renv lock audit reads package records without guessing versions", {
  source(repo_path("R", "renv_lock_audit.R"))

  lockfile <- tempfile(fileext = ".lock")
  writeLines(c(
    "{",
    '  "R": {"Version": "4.5.1"},',
    '  "Packages": {',
    '    "limma": {',
    '      "Package": "limma",',
    '      "Version": "3.0.0"',
    "    },",
    '    "WGCNA": {',
    '      "Package": "WGCNA",',
    '      "Version": "1.0.0"',
    "    }",
    "  }",
    "}"
  ), lockfile)

  testthat::expect_identical(renv_lock_package_names(lockfile), c("limma", "WGCNA"))
  audit <- audit_renv_lock(lockfile, scientific_sentinels = c("limma", "WGCNA"))
  testthat::expect_true(audit$plausibly_full_scientific_lock)
  testthat::expect_length(audit$missing_scientific_sentinels, 0L)
})

testthat::test_that("renv lock audit fails closed when scientific sentinels are absent", {
  source(repo_path("R", "renv_lock_audit.R"))

  lockfile <- tempfile(fileext = ".lock")
  writeLines(c(
    "{",
    '  "Packages": {',
    '    "renv": {',
    '      "Package": "renv",',
    '      "Version": "1.1.5"',
    "    }",
    "  }",
    "}"
  ), lockfile)

  audit <- audit_renv_lock(lockfile, scientific_sentinels = c("limma", "WGCNA"))
  testthat::expect_false(audit$plausibly_full_scientific_lock)
  testthat::expect_identical(audit$missing_scientific_sentinels, c("limma", "WGCNA"))
})

testthat::test_that("current lock status is computed from the current file", {
  source(repo_path("R", "renv_lock_audit.R"))
  lockfile <- testthat::test_path("..", "..", "renv.lock")
  audit <- audit_renv_lock(lockfile)

  testthat::expect_identical(
    audit$plausibly_full_scientific_lock,
    length(audit$missing_scientific_sentinels) == 0L
  )
  testthat::expect_identical(audit$package_count, length(audit$packages))
})

# ===========================================================================
# Phase 12: the lockfile must describe the frozen analysis environment.
#
# The original failure mode was a hand-seeded three-record lockfile in a
# project with no renv infrastructure at all. These tests are dependency-based
# rather than count-based, with only a loose lower bound as a smoke check.
# ===========================================================================

freeze_lock_path <- function() testthat::test_path("..", "..", "renv.lock")

testthat::test_that("the lockfile records the required scientific stack at frozen versions", {
  source(repo_path("R", "renv_lock_audit.R"))
  testthat::skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite required")
  lock <- jsonlite::fromJSON(freeze_lock_path(), simplifyVector = FALSE)

  expected <- c(
    clusterProfiler = "4.18.4", fgsea = "1.36.2", DOSE = "4.4.0",
    BiocParallel = "1.44.0", WGCNA = "1.74", limma = "3.66.0"
  )
  for (pkg in names(expected)) {
    rec <- lock$Packages[[pkg]]
    testthat::expect_false(is.null(rec), info = paste(pkg, "absent from lockfile"))
    testthat::expect_identical(rec$Version, unname(expected[[pkg]]), info = pkg)
    testthat::expect_identical(rec$Package, pkg, info = pkg)
  }
})

testthat::test_that("the lockfile encodes the R 4.5.1 / Bioconductor 3.22 contract", {
  testthat::skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite required")
  lock <- jsonlite::fromJSON(freeze_lock_path(), simplifyVector = FALSE)

  testthat::expect_identical(lock$R$Version, "4.5.1")
  testthat::expect_identical(lock$Bioconductor$Version, "3.22")

  repo_names <- vapply(lock$R$Repositories, function(r) r$Name, character(1))
  repo_urls <- vapply(lock$R$Repositories, function(r) r$URL, character(1))
  testthat::expect_true("CRAN" %in% repo_names)
  testthat::expect_true("BioCsoft" %in% repo_names)
  # Every Bioconductor repository URL must pin the 3.22 release.
  bioc_urls <- repo_urls[grepl("^BioC", repo_names)]
  testthat::expect_gt(length(bioc_urls), 0L)
  testthat::expect_true(all(grepl("/3.22/", bioc_urls, fixed = TRUE)))
})

testthat::test_that("every active direct dependency is recorded in the lockfile", {
  source(repo_path("R", "renv_lock_audit.R"))
  root <- testthat::test_path("..", "..")
  audit <- audit_renv_lock_completeness(freeze_lock_path(), root = root)

  testthat::expect_gt(length(audit$direct_dependencies), 50L)
  testthat::expect_identical(audit$missing_direct, character(0))
  testthat::expect_true(audit$complete)
})

testthat::test_that("the dependency closure has no unresolved package references", {
  source(repo_path("R", "renv_lock_audit.R"))
  audit <- audit_renv_lock_completeness(freeze_lock_path(),
                                        root = testthat::test_path("..", ".."))
  testthat::expect_identical(audit$unresolved_requirements, character(0))
  testthat::expect_identical(audit$duplicate_records, character(0))
})

testthat::test_that("the lockfile has not regressed to the trivial bootstrap state", {
  source(repo_path("R", "renv_lock_audit.R"))
  audit <- audit_renv_lock(freeze_lock_path())

  # The known failure mode: exactly renv/yaml/testthat and nothing else.
  testthat::expect_false(setequal(audit$packages, c("renv", "yaml", "testthat")))
  testthat::expect_true(audit$plausibly_full_scientific_lock)
  testthat::expect_identical(audit$missing_scientific_sentinels, character(0))
  # Loose lower bound only: correctness is asserted dependency-wise above.
  testthat::expect_gt(audit$package_count, 100L)
})

testthat::test_that("recorded versions match the installed library exactly", {
  testthat::skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite required")
  lock <- jsonlite::fromJSON(freeze_lock_path(), simplifyVector = FALSE)
  drift <- character()
  for (nm in names(lock$Packages)) {
    installed <- tryCatch(as.character(utils::packageVersion(nm)),
                          error = function(e) NA_character_)
    if (!identical(lock$Packages[[nm]]$Version, installed)) {
      drift <- c(drift, sprintf("%s lock=%s installed=%s",
                                nm, lock$Packages[[nm]]$Version, installed))
    }
  }
  testthat::expect_identical(drift, character(0),
                             info = paste(utils::head(drift, 5), collapse = " | "))
})

testthat::test_that("source labels are faithful and no remote metadata is invented", {
  testthat::skip_if_not(requireNamespace("jsonlite", quietly = TRUE), "jsonlite required")
  lock <- jsonlite::fromJSON(freeze_lock_path(), simplifyVector = FALSE)
  sources <- vapply(lock$Packages, function(p) p$Source, character(1))
  testthat::expect_true(all(sources %in% c("Repository", "Bioconductor")))
  # Bioconductor records must exist, and none may claim a git remote, because
  # no package in this library was installed from a git source.
  testthat::expect_gt(sum(sources == "Bioconductor"), 0L)
  for (nm in names(lock$Packages)) {
    p <- lock$Packages[[nm]]
    testthat::expect_null(p$RemoteSha, info = nm)
    testthat::expect_null(p$RemoteRef, info = nm)
    testthat::expect_null(p$Hash, info = nm)
  }
})

testthat::test_that("the freeze manifest can inspect the lockfile deterministically", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "publication_freeze_utils.R"))

  a <- freeze_renv_lockfile_state()
  b <- freeze_renv_lockfile_state()
  testthat::expect_identical(a, b)                  # deterministic
  testthat::expect_true(a$complete)
  testthat::expect_identical(a$r_version, "4.5.1")
  testthat::expect_identical(a$bioconductor_version, "3.22")
  # Hash omission is recorded rather than silently implied.
  testthat::expect_false(a$hash_field_recorded)
  testthat::expect_true(nzchar(a$hash_field_note))

  # The accepted gap must disappear only because the lockfile is complete.
  ids <- vapply(freeze_known_gaps(a), function(g) g$id, character(1))
  testthat::expect_false("renv_lock_incomplete" %in% ids)
  incomplete <- a; incomplete$complete <- FALSE
  ids2 <- vapply(freeze_known_gaps(incomplete), function(g) g$id, character(1))
  testthat::expect_true("renv_lock_incomplete" %in% ids2)
})
