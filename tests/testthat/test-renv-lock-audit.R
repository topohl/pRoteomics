testthat::test_that("renv lock audit reads package records without guessing versions", {
  source(testthat::test_path("..", "..", "R", "renv_lock_audit.R"))

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
  source(testthat::test_path("..", "..", "R", "renv_lock_audit.R"))

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
  source(testthat::test_path("..", "..", "R", "renv_lock_audit.R"))
  lockfile <- testthat::test_path("..", "..", "renv.lock")
  audit <- audit_renv_lock(lockfile)

  testthat::expect_identical(
    audit$plausibly_full_scientific_lock,
    length(audit$missing_scientific_sentinels) == 0L
  )
  testthat::expect_identical(audit$package_count, length(audit$packages))
})
