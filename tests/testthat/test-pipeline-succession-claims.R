# pipeline.yml succession claims must be substantively true, not merely
# resolvable.
#
# The `legacy:` block lists 19 scripts excluded from canonical runnable stages.
# Its `replacement:` field is DUAL-PURPOSE by design: sometimes a successor
# script path, sometimes a prose description of what the excluded thing is
# ("repository audit utility"). A prose value is correct usage, not a broken
# path, and a test that flags it as one is wrong.
#
# What went undetected for far longer is the opposite failure. Two entries
# named analysis/preprocessing/extract_protigy_contrasts.R as the replacement
# for archive/01_preprocessing/01_impute.r and 02_excel_convert.r. That target
# exists, so any existence check passed - but it is three steps DOWNSTREAM of
# the first and on the far side of the external ProTigy boundary from the
# second. It consumes their descendants. The dependency arrow was inverted, and
# both output families are still REQUIRED active inputs with no active
# producer.
#
# So existence is not the property worth testing. These tests pin the
# adjudicated RELATIONSHIP, recorded in audits/phase6i8_succession_adjudication.csv,
# and check that it still matches what pipeline.yml says.

source(testthat::test_path("..", "..", "R", "paths.R"))

ADJUDICATION <- repo_path("audits", "phase6i8_succession_adjudication.csv")
PIPELINE <- repo_path("pipeline.yml")

legacy_entries <- function() {
  y <- yaml::read_yaml(PIPELINE)
  lapply(y$legacy, function(e) list(script = e$script, replacement = e$replacement,
                                    status = e$status))
}
is_real_path <- function(v) grepl("[.][Rr]$", v) && file.exists(repo_path(v))

testthat::test_that("every legacy entry is adjudicated", {
  testthat::skip_if_not(file.exists(PIPELINE) && file.exists(ADJUDICATION),
                        "pipeline or adjudication absent")
  e <- legacy_entries()
  a <- utils::read.csv(ADJUDICATION, stringsAsFactors = FALSE)
  testthat::expect_identical(length(e), nrow(a))
  testthat::expect_setequal(vapply(e, function(x) x$script, character(1)),
                            a$superseded_script)
  for (col in c("superseded_script", "value_kind", "adjudication", "basis"))
    testthat::expect_true(col %in% names(a), info = col)
  testthat::expect_true(all(nzchar(a$basis)),
    info = "an adjudication carries no stated basis")
})

testthat::test_that("the adjudication still describes what pipeline.yml says", {
  # The audit is only worth anything while it tracks the file. If someone edits
  # a replacement value from prose to a path, or the reverse, this fails.
  testthat::skip_if_not(file.exists(PIPELINE) && file.exists(ADJUDICATION),
                        "pipeline or adjudication absent")
  e <- legacy_entries()
  a <- utils::read.csv(ADJUDICATION, stringsAsFactors = FALSE)
  for (x in e) {
    row <- a[a$superseded_script == x$script, , drop = FALSE]
    testthat::expect_identical(nrow(row), 1L, info = x$script)
    expected <- if (is_real_path(x$replacement)) "REAL_SCRIPT_PATH" else "PROSE_NOT_PATH"
    testthat::expect_identical(row$value_kind[1], expected,
      info = paste(x$script, "-> value_kind drifted from the actual value"))
  }
})

testthat::test_that("prose replacement values are not treated as broken paths", {
  # The dual-purpose field is deliberate. These entries describe what an
  # excluded script IS; there is no file to find and none should be demanded.
  testthat::skip_if_not(file.exists(ADJUDICATION), "adjudication absent")
  a <- utils::read.csv(ADJUDICATION, stringsAsFactors = FALSE)
  prose <- a[a$value_kind == "PROSE_NOT_PATH", , drop = FALSE]
  testthat::expect_gt(nrow(prose), 0L)

  # Prose comes in two kinds and they must not be flattened together. Most of
  # it DESCRIBES an excluded script ("repository audit utility") and gets the
  # PROSE_NOT_PATH verdict. But the two entries corrected in Phase 6I.8 are
  # prose that STATES a relationship - "no direct successor; outputs preserved
  # and still consumed" - and those keep their substantive verdict, because
  # that verdict is the finding.
  #
  # What no prose value may carry is a POSITIVE successor verdict: prose names
  # no successor, so nothing can have been shown to assume the responsibility.
  positive <- c("EXACT_SUCCESSOR", "FUNCTIONAL_SUCCESSOR", "PARTIAL_SUCCESSOR")
  testthat::expect_false(any(prose$adjudication %in% positive),
    info = paste("a prose value claims a successor:",
                 paste(prose$superseded_script[prose$adjudication %in% positive],
                       collapse = ", ")))
  testthat::expect_true(all(prose$adjudication %in%
                              c("PROSE_NOT_PATH", "NO_SUCCESSOR", "RELATED_BUT_NOT_SUCCESSOR")))
})

testthat::test_that("a real-path replacement carries a real successor verdict", {
  testthat::skip_if_not(file.exists(ADJUDICATION), "adjudication absent")
  a <- utils::read.csv(ADJUDICATION, stringsAsFactors = FALSE)
  paths <- a[a$value_kind == "REAL_SCRIPT_PATH", , drop = FALSE]
  testthat::expect_gt(nrow(paths), 0L)
  allowed <- c("EXACT_SUCCESSOR", "FUNCTIONAL_SUCCESSOR", "PARTIAL_SUCCESSOR",
               "RELATED_BUT_NOT_SUCCESSOR", "NO_SUCCESSOR")
  testthat::expect_true(all(paths$adjudication %in% allowed))
  # and nothing may sit at UNKNOWN or FALSE
  testthat::expect_false(any(paths$adjudication %in% c("UNKNOWN", "FALSE")))
})

testthat::test_that("the two corrected claims cannot silently return", {
  # The specific regression. extract_protigy_contrasts.R must never again be
  # named as the replacement for either preprocessing script: it consumes
  # their descendants rather than replacing them.
  testthat::skip_if_not(file.exists(PIPELINE), "pipeline absent")
  e <- legacy_entries()
  for (s in c("archive/01_preprocessing/01_impute.r",
              "archive/01_preprocessing/02_excel_convert.r")) {
    hit <- Filter(function(x) identical(x$script, s), e)
    testthat::expect_identical(length(hit), 1L, info = s)
    testthat::expect_false(
      grepl("extract_protigy_contrasts", hit[[1]]$replacement, fixed = TRUE),
      info = paste(s, "again names a downstream consumer as its successor"))
    testthat::expect_true(grepl("no direct successor", hit[[1]]$replacement, fixed = TRUE),
      info = paste(s, "no longer states that nothing replaced it"))
    # and the status must still say the outputs are live inputs, because that
    # is the fact that makes "no successor" a finding rather than a shrug
    testthat::expect_true(grepl("REQUIRED active input", hit[[1]]$status, fixed = TRUE),
      info = paste(s, "stopped recording that its outputs are still required"))
  }
})

testthat::test_that("a historical producer is not declared as an active step", {
  # Correcting the succession semantics must not imply archive scripts became
  # runnable. The legacy block is an exclusion list and stays one.
  testthat::skip_if_not(file.exists(PIPELINE), "pipeline absent")
  y <- yaml::read_yaml(PIPELINE)
  legacy_scripts <- vapply(y$legacy, function(e) e$script, character(1))
  archived <- legacy_scripts[startsWith(legacy_scripts, "archive/")]
  testthat::expect_gt(length(archived), 0L)

  # none of them may appear as a registered stage script elsewhere in the file
  stages <- y[setdiff(names(y), "legacy")]
  flat <- unlist(stages, use.names = FALSE)
  flat <- flat[vapply(flat, is.character, logical(1))]
  for (s in archived)
    testthat::expect_false(any(flat == s),
      info = paste(s, "appears as a registered pipeline step"))
})
