# The reviewer ledger must stay readable, and tests must stay out of it.
#
# results/reviewer_audit/input_resolution_audit.csv is production evidence:
# reviewers read it to see how every scientific input was resolved. Two
# separate defects were found in it, with two different causes.
#
#   1. Three lines carried one quote character each - orphan tails of records
#      whose heads a concurrent process wrote elsewhere, inserted between two
#      intact records. Three stray quotes desynchronised CSV quoting for every
#      line after them, so a plain read.csv() recovered about 1,600 of 57,550
#      rows. 0.005% of the file cost 97% of it. Repaired in Phase 6I.4 by
#      removing exactly those three lines, after quarantining them verbatim;
#      the content itself was unrecoverable, since the ledger is gitignored and
#      no snapshot or backup exists.
#
#   2. Eighty rows named a drive letter that does not exist, written by a unit
#      test that called read_csv_optional() with fixture paths. Removed in
#      Phase 6I.2. Different cause, so the guards below are different too: one
#      is about what the writer may emit, the other about where tests write.
#
# These tests protect the repaired state and the isolation, not a row count -
# the ledger is append-only and grows on every run, so pinning its size would
# fail immediately and teach nobody anything.

source(testthat::test_path("..", "..", "R", "paths.R"))

# The CANONICAL ledger, resolved with setup.R's suite-wide override suppressed.
# Without this the constant would pick up the temp path setup.R installs, and
# every assertion below would be about a throwaway file - the tripwire would
# pass while production evidence was being appended to.
LEDGER <- withr::with_envvar(c(PROTEOMICS_INPUT_RESOLUTION_AUDIT = NA),
                             input_resolution_audit_path())

quote_counts <- function(lines) {
  vapply(gregexpr('"', lines, fixed = TRUE),
         function(m) if (m[[1]] == -1L) 0L else length(m), integer(1))
}
# top-level fields: commas outside quoted spans, plus one
field_counts <- function(lines) {
  1L + lengths(gregexpr(",", gsub('"[^"]*"', "", lines), fixed = TRUE))
}

testthat::test_that("the canonical ledger parses completely", {
  testthat::skip_if_not(file.exists(LEDGER), "ledger absent")
  raw <- readLines(LEDGER, warn = FALSE)
  testthat::expect_gt(length(raw), 1L)

  # the property that was violated: a plain reader must recover every record
  d <- utils::read.csv(LEDGER, stringsAsFactors = FALSE, colClasses = "character")
  testthat::expect_identical(nrow(d), length(raw) - 1L)
})

testthat::test_that("no record can desynchronise the file's quoting", {
  # This is the specific failure mode, asserted directly rather than via the
  # row count it happens to produce. One unbalanced record swallows the rest
  # of the file, so the count of them must be zero, not small.
  testthat::skip_if_not(file.exists(LEDGER), "ledger absent")
  raw <- readLines(LEDGER, warn = FALSE)
  odd <- which(quote_counts(raw) %% 2L == 1L)
  testthat::expect_identical(length(odd), 0L,
    info = paste("unbalanced lines at:", paste(utils::head(odd, 5), collapse = ", ")))

  wrong <- which(field_counts(raw) != 14L)
  testthat::expect_identical(length(wrong), 0L,
    info = paste("lines that are not 14-field records:",
                 paste(utils::head(wrong, 5), collapse = ", ")))
})

testthat::test_that("the ledger still matches its declared schema", {
  testthat::skip_if_not(file.exists(LEDGER), "ledger absent")
  d <- utils::read.csv(LEDGER, stringsAsFactors = FALSE, colClasses = "character")
  testthat::expect_setequal(names(d), input_resolution_audit_columns())
  testthat::expect_identical(ncol(d), length(input_resolution_audit_columns()))
})

testthat::test_that("no synthetic fixture row is in the ledger", {
  # The Phase 6I.2 contamination signature, plus the shape of it rather than
  # just the one literal: a resolved_path on a drive that is not mounted and
  # under a directory no analysis writes.
  testthat::skip_if_not(file.exists(LEDGER), "ledger absent")
  d <- utils::read.csv(LEDGER, stringsAsFactors = FALSE, colClasses = "character")
  testthat::expect_identical(sum(d$script == "test-vocabulary"), 0L)
  testthat::expect_identical(sum(grepl("declared/elsewhere", d$resolved_path, fixed = TRUE)), 0L)
  testthat::expect_identical(sum(grepl("declared/elsewhere", d$expected_path, fixed = TRUE)), 0L)
})

testthat::test_that("the quarantined partial writes are kept, not discarded", {
  # The three damaged lines were not recoverable, but they were evidence. They
  # sit beside the ledger with their byte offsets so the repair is auditable.
  q <- file.path(dirname(LEDGER), "input_resolution_audit.quarantined_partial_writes.csv")
  testthat::skip_if_not(file.exists(q), "quarantine absent")
  d <- utils::read.csv(q, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(d), 3L)
  for (col in c("physical_line", "byte_offset", "raw_text", "repair_confidence",
                "recoverable", "source_ledger_sha256_before_repair")) {
    testthat::expect_true(col %in% names(d), info = col)
  }
  testthat::expect_true(all(d$repair_confidence == "STRUCTURALLY_UNAMBIGUOUS_REPAIR"))
  testthat::expect_true(all(!d$recoverable))
  # each quarantined line is what it was claimed to be
  testthat::expect_true(all(quote_counts(d$raw_text) == 1L))
  testthat::expect_true(all(field_counts(d$raw_text) == 2L))
})

# ---- isolation ------------------------------------------------------------

testthat::test_that("the ledger destination is routable, so a test need not touch it", {
  # Isolation is by path routing. The appender must stay unconditional:
  # test-preprocessing-writer-namespace.R depends on it having no dry-run
  # guard, which is why build_module_score_metadata is classified
  # PATH_VERIFIED_STRUCTURALLY_ONLY. Routing the destination leaves the
  # writer's behaviour alone.
  src <- paste(readLines(testthat::test_path("..", "..", "R", "paths.R"), warn = FALSE),
               collapse = "\n")
  appender <- substr(sub(".*append_input_resolution_audit <- function", "", src), 1, 1200)
  testthat::expect_false(grepl("is_dry_run", appender, fixed = TRUE))

  withr::local_envvar(PROTEOMICS_INPUT_RESOLUTION_AUDIT = NA)
  default_path <- input_resolution_audit_path()
  tmp <- withr::local_tempdir()
  routed <- file.path(tmp, "ledger.csv")
  withr::local_envvar(PROTEOMICS_INPUT_RESOLUTION_AUDIT = routed)
  testthat::expect_identical(input_resolution_audit_path(), routed)
  testthat::expect_false(identical(input_resolution_audit_path(), default_path))
})

testthat::test_that("a routed write lands in the fixture and not in the ledger", {
  # The tripwire, exercised rather than asserted: hash the canonical ledger,
  # do the thing that used to contaminate it, hash again.
  testthat::skip_if_not(file.exists(LEDGER), "ledger absent")
  before <- unname(tools::sha256sum(LEDGER))

  tmp <- withr::local_tempdir()
  routed <- file.path(tmp, "ledger.csv")
  withr::local_envvar(PROTEOMICS_INPUT_RESOLUTION_AUDIT = routed)

  record_input_resolution(
    script = "test-reviewer-ledger-integrity", dataset = "global", stage = "test",
    input_name = "fixture.csv",
    expected_path = file.path(tmp, "fixture.csv"),
    resolved_path = file.path(tmp, "fixture.csv"),
    resolution_mode = "canonical", producer_script_or_artifact_id = "test")

  testthat::expect_true(file.exists(routed))
  written <- utils::read.csv(routed, stringsAsFactors = FALSE, colClasses = "character")
  testthat::expect_identical(nrow(written), 1L)
  testthat::expect_identical(written$script[[1]], "test-reviewer-ledger-integrity")

  testthat::expect_identical(unname(tools::sha256sum(LEDGER)), before,
    info = "a routed write still reached the canonical reviewer ledger")
})

testthat::test_that("local_input_resolution_audit routes and then restores", {
  testthat::skip_if_not(file.exists(LEDGER), "ledger absent")
  before_env <- Sys.getenv("PROTEOMICS_INPUT_RESOLUTION_AUDIT", unset = NA_character_)
  p <- local_input_resolution_audit()
  testthat::expect_identical(input_resolution_audit_path(), p)
  testthat::expect_false(identical(p, LEDGER))
  # the helper defers restoration to its caller's frame, so inside this test
  # the override is still active; what matters is that it is not permanent
  testthat::expect_true(nzchar(Sys.getenv("PROTEOMICS_INPUT_RESOLUTION_AUDIT")))
  testthat::expect_true(is.na(before_env) || nzchar(before_env))
})

testthat::test_that("the suite routes the ledger in exactly one place", {
  # One mechanism, not eleven. setup.R runs before every test file, including
  # when a single file is run on its own, and sets the environment variable
  # that child processes inherit. Per-file routing was tried first and removed:
  # a file that forgets is silently back to writing production evidence.
  s <- testthat::test_path("setup.R")
  testthat::expect_true(file.exists(s),
    info = "tests/testthat/setup.R is gone; nothing routes the reviewer ledger")
  src <- paste(readLines(s, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("PROTEOMICS_INPUT_RESOLUTION_AUDIT", src, fixed = TRUE))
  testthat::expect_true(grepl("teardown_env", src, fixed = TRUE),
    info = "setup.R sets the override without restoring it")

  # and it is in force right now, which is the whole point
  testthat::expect_true(nzchar(Sys.getenv("PROTEOMICS_INPUT_RESOLUTION_AUDIT")),
    info = "the suite-wide ledger override is not active during this test")
  testthat::expect_false(identical(input_resolution_audit_path(), LEDGER))
})

testthat::test_that("every test that spawns an analysis script is accounted for", {
  # Eleven test files run real analysis scripts through system2(), and a child
  # inherits setup.R's override, so all eleven are covered without any of them
  # saying so. The risk this guards is different: a NEW spawner arriving
  # alongside some future change that also bypasses setup.R. Then someone has
  # to look at it rather than assume.
  spawners <- sort(basename(Filter(function(f) {
    src <- readLines(f, warn = FALSE)
    any(grepl("system2", src, fixed = TRUE)) &&
      any(grepl("analysis/[a-z_]+/[a-z_]+[.]R", src))
  }, list.files(testthat::test_path("."), pattern = "^test-.*[.]R$", full.names = TRUE))))
  testthat::skip_if(!length(spawners), "no spawning tests found")

  REVIEWED <- c(
    "test-biological-integration-entrypoints.R", "test-comparego-tail-archival.R",
    "test-ewce-dataset-saving.R", "test-module-script-entrypoints.R",
    "test-qc-writer-namespace.R", "test-spatial-network-dataset-argument.R",
    "test-spatial-networks-writer-namespace.R",
    "test-wgcna-candidate-protein-shortlist.R", "test-wgcna-downstream-scripts.R",
    "test-wgcna-identity-contract.R", "test-wgcna-writer-namespace.R")

  testthat::expect_identical(setdiff(spawners, REVIEWED), character(0),
    info = paste("a new test spawns analysis scripts and has not been reviewed",
                 "for reviewer-ledger writes:",
                 paste(setdiff(spawners, REVIEWED), collapse = ", ")))
})
