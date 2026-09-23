# The archive boundary, and the stochastic code behind it.
#
# Phase 6H.8 established that active statistical RNG is fully seeded. What was
# left open was ARCHIVE_ONLY_RNG_DEBT: whether stochastic code preserved under
# archive/ needs anything doing about it. Phase 6I.6 audited it and closed it
# as A1 - CLOSED_BY_ARCHIVAL_BOUNDARY.
#
# Two properties have to keep holding for that closure to stay true, and they
# are different in kind:
#
#   1. The boundary. archive/ is declared "superseded code kept for provenance"
#      and explicitly not "executable by the registry"
#      (docs/REPOSITORY_ARCHITECTURE.md). If active code ever sources an
#      archived script, its RNG stops being historical and becomes live, and
#      the audit's conclusion is void.
#
#   2. The inventory. There are exactly two stochastic calls under archive/.
#      One is seeded and its seed is recorded three ways; the other is the
#      superseded compareGO bootstrap already adjudicated in Phase 6H. A third
#      appearing means something was archived without being audited.
#
# Deliberately NOT tested: the content of archived files. Archive bytes are
# provenance. Adding a seed to historical source to satisfy a modern standard
# would falsify the record of how the analysis actually ran.

source(testthat::test_path("..", "..", "R", "paths.R"))

r_files_under <- function(prefix) {
  all <- system2("git", c("-C", shQuote(repo_path()), "ls-files"), stdout = TRUE)
  all[grepl("[.][Rr]$", all) & startsWith(all, prefix)]
}

# Parse-aware: a path mentioned in a comment or a message string is not an edge.
exec_calls_naming <- function(path, needle) {
  ex <- tryCatch(parse(repo_path(path), keep.source = TRUE), error = function(e) NULL)
  if (is.null(ex)) return(list(total = 0L, hits = character(0)))
  EXEC <- c("source", "sys.source", "eval", "evalq", "debugSource")
  total <- 0L; hits <- character(0)
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    fn <- e[[1]]
    nm <- if (is.name(fn)) as.character(fn) else
      if (is.call(fn) && identical(as.character(fn[[1]]), "::")) as.character(fn[[3]]) else NA_character_
    if (!is.na(nm) && nm %in% EXEC) {
      total <<- total + 1L
      txt <- paste(deparse(e), collapse = " ")
      if (grepl(needle, txt, fixed = TRUE)) hits <<- c(hits, substr(txt, 1, 160))
    }
    for (i in seq_along(e)) {
      if (identical(e[[i]], quote(expr = ))) next
      walk(e[[i]])
    }
  }
  for (i in seq_along(ex)) walk(ex[[i]])
  list(total = total, hits = hits)
}

testthat::test_that("no non-archive code executes an archived script", {
  files <- setdiff(r_files_under(""), r_files_under("archive/"))
  testthat::skip_if(!length(files), "no tracked R files")

  total <- 0L; offenders <- character(0)
  for (f in files) {
    r <- exec_calls_naming(f, "archive")
    total <- total + r$total
    if (length(r$hits)) offenders <- c(offenders, paste0(f, ": ", r$hits))
  }
  # The denominator matters: a scan that found no source() calls at all would
  # report zero offenders while proving nothing.
  testthat::expect_gt(total, 100L)
  testthat::expect_identical(offenders, character(0),
    info = paste("archive code is executed by:", paste(offenders, collapse = " | ")))
})

testthat::test_that("archive is still declared non-runnable", {
  arch <- repo_path("docs", "REPOSITORY_ARCHITECTURE.md")
  testthat::skip_if_not(file.exists(arch), "architecture doc absent")
  txt <- paste(readLines(arch, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("superseded code kept for provenance", txt, fixed = TRUE),
    info = "archive/ is no longer declared as provenance-only")
})

testthat::test_that("the archive stochastic inventory is exactly what was audited", {
  # Parse-aware census over every tracked archive R file.
  RNG <- c("sample", "sample.int", "runif", "rnorm", "rbinom", "rpois", "rexp",
           "rgamma", "rbeta", "rlnorm", "rmvnorm", "rcauchy", "rchisq",
           "rweibull", "rhyper", "rnbinom", "rgeom", "rmultinom",
           "jitter", "slice_sample", "sample_n", "sample_frac")
  files <- r_files_under("archive/")
  testthat::skip_if(!length(files), "no archive R files in this checkout")

  found <- list()
  for (f in files) {
    pd <- tryCatch(utils::getParseData(parse(repo_path(f), keep.source = TRUE)),
                   error = function(e) NULL)
    if (is.null(pd) || !nrow(pd)) next
    calls <- pd[pd$token == "SYMBOL_FUNCTION_CALL", c("line1", "text")]
    hit <- calls[calls$text %in% RNG, , drop = FALSE]
    if (nrow(hit)) found[[length(found) + 1L]] <-
      data.frame(file = f, fn = hit$text, stringsAsFactors = FALSE)
  }
  inv <- if (length(found)) do.call(rbind, found) else
    data.frame(file = character(), fn = character())

  expected <- c(
    "archive/01_preprocessing/01_impute.r" = "rnorm",
    "archive/04_differential_expression_enrichment/legacy/02_compareGO_superseded_tail.r" = "slice_sample")
  testthat::expect_identical(nrow(inv), 2L,
    info = paste("archive stochastic calls changed:",
                 paste(paste0(inv$file, ":", inv$fn), collapse = ", ")))
  testthat::expect_setequal(inv$file, names(expected))
  for (f in names(expected)) {
    testthat::expect_identical(inv$fn[inv$file == f], unname(expected[f]), info = f)
  }
})

testthat::test_that("the seeded archive draw still records its seed", {
  # The one archive RNG whose output is publication-supporting: the imputation
  # that feeds WGCNA. It is exactly regenerable, and that rests on the seed
  # being recorded rather than remembered. Source is checked, not rewritten.
  f <- repo_path("archive", "01_preprocessing", "01_impute.r")
  testthat::skip_if_not(file.exists(f), "archived imputation absent")
  src <- readLines(f, warn = FALSE)
  testthat::expect_true(any(grepl("IMPUTATION_SEED <- 42L", src, fixed = TRUE)),
    info = "the recorded imputation seed changed in archived source")
  testthat::expect_true(any(grepl("subset_seed <- IMPUTATION_SEED + idx - 1L", src, fixed = TRUE)))
  # and the seed reaches the draw rather than defaulting to NULL
  testthat::expect_true(any(grepl("impute_normal(df_filtered, numeric_cols, seed = subset_seed)",
                                  src, fixed = TRUE)))

  # the executed run wrote the seed out, per subset
  qc <- repo_path("data", "processed", "01_preprocessing", "impute", "imputation_qc.csv")
  testthat::skip_if_not(file.exists(qc), "imputation QC absent")
  d <- utils::read.csv(qc, stringsAsFactors = FALSE)
  for (col in c("base_seed", "subset_seed", "celltype_layer", "output_path")) {
    testthat::expect_true(col %in% names(d), info = col)
  }
  testthat::expect_true(all(d$base_seed == 42L))
  testthat::expect_identical(anyDuplicated(d$subset_seed), 0L,
    info = "two subsets shared a seed; the per-subset draws are not independent")
})

testthat::test_that("active statistical RNG remains fully seeded", {
  # A boundary re-check, not a new sweep. Phase 6H.8 closed active RNG; this
  # only guards against an archive reclassification quietly moving an unseeded
  # call into the active tree.
  RNG <- c("sample", "sample.int", "runif", "rnorm", "rbinom", "rpois", "rexp",
           "rgamma", "rbeta", "rlnorm", "slice_sample", "sample_n", "sample_frac")
  SEED <- c("set.seed", "with_seed", "local_seed", "with_preserve_seed", "RNGkind")
  files <- Filter(function(f) !grepl("^(archive|tests|audits)/", f), r_files_under(""))
  testthat::skip_if(!length(files), "no active R files")

  unseeded <- character(0); n_rng <- 0L
  for (f in files) {
    pd <- tryCatch(utils::getParseData(parse(repo_path(f), keep.source = TRUE)),
                   error = function(e) NULL)
    if (is.null(pd) || !nrow(pd)) next
    calls <- pd[pd$token == "SYMBOL_FUNCTION_CALL", c("line1", "text")]
    r <- calls[calls$text %in% RNG, , drop = FALSE]
    if (!nrow(r)) next
    n_rng <- n_rng + nrow(r)
    # A file carrying an RNG draw must also carry seed control. Callers seed
    # helpers defined above them, so proximity is not the test - presence is.
    if (!any(calls$text %in% SEED)) unseeded <- c(unseeded, paste0(f, ":", r$line1[1]))
  }
  testthat::expect_gt(n_rng, 0L)
  testthat::expect_identical(unseeded, character(0),
    info = paste("active RNG with no seed control in file:",
                 paste(unseeded, collapse = ", ")))
})
