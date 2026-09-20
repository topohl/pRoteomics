# Regression cover for the generalized historical-WGCNA construction scanner.
#
# Two classes of defect are pinned here.
#
# 1. SPAN ARITHMETIC. A parse node spanning exactly two lines has no interior
#    line, and (line1+1):(line2-1) descends instead of being empty. In Batch 3A
#    that silently spliced the wrong lines and made five real call sites look
#    absent, which would have left them unmigrated while the family was
#    reported complete. Single-line, 2-line and 3+-line spans are all fixtured.
#
# 2. DETECTION BY STRUCTURE, NOT SUBSTRING. "06_modules_WGCNA" appears in
#    comments, in user-facing messages and in frozen manifests. Only a call
#    that actually builds a path may be reported.

source(testthat::test_path("..", "..", "R", "paths.R"))
repo_file <- function(...) repo_path(...)
source(repo_path("R", "utilities", "wgcna_construction_scan.R"))
source(repo_path("R", "wgcna_paths.R"))

write_fixture <- function(lines) {
  f <- tempfile(fileext = ".R")
  writeLines(lines, f)
  f
}

# ------------------------------------------------------------------- spans

test_that("wcs_span_text extracts a single-line span", {
  lines <- c("x <- path_results(\"tables\", \"06_modules_WGCNA\")")
  expect_equal(wcs_span_text(lines, 1, 6, 1, nchar(lines[1])),
               "path_results(\"tables\", \"06_modules_WGCNA\")")
})

test_that("wcs_span_text extracts a two-line span with no interior line", {
  lines <- c("p <- path_results(\"tables\", \"06_modules_WGCNA\", ds,",
             "                  \"file.csv\")")
  got <- wcs_span_text(lines, 1, 6, 2, nchar(lines[2]))
  expect_true(grepl("^path_results", got))
  expect_true(grepl("file[.]csv", got))
  expect_equal(length(strsplit(got, "\n")[[1]]), 2L)
  # the decisive property: the interior is empty, not a descending range
  expect_false(grepl("p <- ", got))
})

test_that("wcs_span_text extracts a three-line span including its interior", {
  lines <- c("p <- path_results(",
             "  \"tables\", \"06_modules_WGCNA\", ds,",
             "  \"file.csv\")")
  got <- wcs_span_text(lines, 1, 6, 3, nchar(lines[3]))
  expect_equal(length(strsplit(got, "\n")[[1]]), 3L)
  expect_true(grepl("06_modules_WGCNA", got))
})

test_that("wcs_span_text extracts a five-line span including all interior lines", {
  lines <- c("p <- path_results(", "  \"tables\",", "  \"06_modules_WGCNA\",",
             "  ds,", "  \"file.csv\")")
  got <- wcs_span_text(lines, 1, 6, 5, nchar(lines[5]))
  expect_equal(length(strsplit(got, "\n")[[1]]), 5L)
  expect_true(grepl("ds", got))
})

test_that("wcs_span_text refuses a malformed or out-of-range span", {
  lines <- c("a", "b")
  expect_error(wcs_span_text(lines, 2, 1, 1, 1), "malformed span")
  expect_error(wcs_span_text(lines, 1, 1, 9, 1), "outside file")
  expect_error(wcs_span_text(lines, NA, 1, 2, 1), "malformed span")
})

test_that("every multi-line span round-trips to parseable code", {
  # the property the converter depends on: a span's text must re-parse
  lines <- c("q <- list(", "  a = path_results(", "    \"tables\", \"06_modules_WGCNA\",",
             "    ds, \"f.csv\"", "  )", ")")
  f <- write_fixture(lines)
  got <- wcs_scan_file(f)
  expect_equal(nrow(got), 1L)
  expect_silent(str2lang(gsub("\n", " ", wcs_span_text(
    readLines(f), got$line1, got$col1, got$line2, got$col2))))
})

# -------------------------------------------------------------- detection

test_that("a plain path_results construction is detected with its identity", {
  f <- write_fixture(c(
    "p <- path_results(\"tables\", \"06_modules_WGCNA\", \"interpretable_summary\",",
    "                  \"microglia\", \"WGCNA_final_label_lookup.csv\")"))
  r <- wcs_scan_file(f)
  expect_equal(nrow(r), 1L)
  expect_equal(r$artifact_family, "interpretable_summary")
  expect_equal(r$artifact, "WGCNA_final_label_lookup.csv")
  expect_equal(r$kind, "tables")
})

test_that("a file.path construction spelled without path_results is detected", {
  f <- write_fixture(
    "p <- file.path(\"results\", \"tables\", \"06_modules_WGCNA\", \"group_effects\", ds, \"m.csv\")")
  r <- wcs_scan_file(f)
  expect_equal(nrow(r), 1L)
  expect_equal(r$artifact_family, "group_effects")
})

test_that("a construction through an alias VARIABLE is detected", {
  # the form the original path_results-only oracle could not see
  f <- write_fixture(c(
    "TABLE_DIR <- path_results(\"tables\", \"06_modules_WGCNA\")",
    "p <- file.path(TABLE_DIR, \"interpretable_summary\", DATASET, \"WGCNA_inferential_handoff.csv\")"))
  r <- wcs_scan_file(f)
  aliased <- r[!is.na(r$artifact_family), ]
  expect_equal(nrow(aliased), 1L)
  expect_equal(aliased$artifact_family, "interpretable_summary")
  expect_equal(aliased$artifact, "WGCNA_inferential_handoff.csv")
  expect_equal(aliased$form, "alias")
})

test_that("a construction through an alias FUNCTION is detected", {
  f <- write_fixture(c(
    "base <- function(...) path_results(\"tables\", ...)",
    "p <- base(\"06_modules_WGCNA\", \"module_annotation\", ds, \"a.csv\")"))
  r <- wcs_scan_file(f)
  hit <- r[!is.na(r$artifact_family), ]
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$artifact_family, "module_annotation")
  expect_equal(hit$kind, "tables")
})

test_that("a Sys.glob over the historical tree is detected and flagged", {
  f <- write_fixture(c(
    "g <- Sys.glob(file.path(repo_root(), \"results\", \"tables\",",
    "  \"06_modules_WGCNA\", \"group_effects\", \"*\", \"*group_effects.csv\"))"))
  r <- wcs_scan_file(f)
  expect_true(any(r$has_wildcard))
  expect_true(any(r$artifact_family == "group_effects", na.rm = TRUE))
})

test_that("nested builder calls are counted once, not twice", {
  f <- write_fixture(
    "g <- Sys.glob(file.path(\"results\", \"tables\", \"06_modules_WGCNA\", \"x\", \"*.csv\"))")
  r <- wcs_scan_file(f)
  expect_equal(nrow(r), 1L)
})

test_that("a comment or message string is NOT a construction", {
  f <- write_fixture(c(
    "# reads results/tables/06_modules_WGCNA/interpretable_summary/...",
    "msg <- \"see results/tables/06_modules_WGCNA for details\"",
    "stop(\"missing 06_modules_WGCNA input\")"))
  expect_null(wcs_scan_file(f))
})

test_that("a normalized resolver call is NOT reported as historical", {
  f <- write_fixture(
    "p <- wgcna_interpretable_artifact(\"WGCNA_final_label_lookup.csv\", ds)")
  expect_null(wcs_scan_file(f))
})

test_that("data/processed historical carriers are detected with kind processed", {
  f <- write_fixture(
    "p <- path_processed(\"06_modules_WGCNA\", \"01_WGCNA\", ds, \"net.rds\")")
  r <- wcs_scan_file(f)
  expect_equal(nrow(r), 1L)
  expect_equal(r$kind, "processed")
  expect_equal(r$artifact_family, "01_WGCNA")
})

test_that("a runtime-determined family segment is reported as dynamic, not guessed", {
  f <- write_fixture(
    "p <- path_results(\"tables\", \"06_modules_WGCNA\", fam, ds, \"x.csv\")")
  r <- wcs_scan_file(f)
  expect_equal(nrow(r), 1L)
  expect_true(is.na(r$artifact_family))
  expect_true(r$dynamic_segments >= 1)
})

# ------------------------------------------------------------ classifier

test_that("classification separates fixtures, oracles and tooling from live reads", {
  df <- data.frame(file = c("analysis/wgcna/x.R",
                            "tests/testthat/test-y.R",
                            "R/data_contracts/publication_freeze_utils.R",
                            "tools/audit_wgcna_support_status.R"),
                   stringsAsFactors = FALSE)
  expect_equal(wcs_classify(df),
               c("LIVE_READER", "TEST_FIXTURE", "FREEZE_ORACLE", "TOOLING_REFERENCE"))
})

test_that("the scanner is the single source of construction truth for the audit", {
  # guards against a second, drifting copy of the detection rule
  scan_src <- readLines(repo_file("R", "utilities", "wgcna_construction_scan.R"), warn = FALSE)
  expect_true(any(grepl("^wcs_scan_file <- function", scan_src)))
  expect_equal(sum(grepl("^wcs_span_text <- function", scan_src)), 1L)
})

# ---------------------------------------------- directory selection semantics
#
# Batch 3C: wgcna_dir_any() picks the populated directory. A Stage-01 dataset
# directory holds only sub-directories (modules/, supermodules/), so a
# top-level-only "has files" test judged it empty and returned the normalized
# directory instead -- a path that does not exist -- while the populated
# historical one was ignored. The 6G.5 guard in the other direction still has
# to hold: a normalized tree of empty directories must never shadow it.

test_that("a directory whose content is sub-directories counts as populated", {
  root <- file.path(tempdir(), "wg_pop")
  unlink(root, recursive = TRUE)
  dir.create(file.path(root, "modules"), recursive = TRUE)
  expect_false(wg_has_files(root))            # no file anywhere yet
  writeLines("x", file.path(root, "modules", "a.csv"))
  expect_true(wg_has_files(root))             # file nested one level down
})

test_that("an empty normalized tree never shadows a populated historical one", {
  empty <- file.path(tempdir(), "wg_empty")
  unlink(empty, recursive = TRUE)
  dir.create(file.path(empty, "tables", "modules"), recursive = TRUE)
  expect_false(wg_has_files(empty))
})

test_that("a directory holding only a loose file is populated", {
  d <- file.path(tempdir(), "wg_loose")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE)
  writeLines("x", file.path(d, "f.csv"))
  expect_true(wg_has_files(d))
})

test_that("a missing or empty directory is not populated", {
  expect_false(wg_has_files(file.path(tempdir(), "wg_absent_dir")))
  d <- file.path(tempdir(), "wg_bare"); unlink(d, recursive = TRUE); dir.create(d)
  expect_false(wg_has_files(d))
  expect_false(wg_has_files(NA_character_))
  expect_false(wg_has_files(character(0)))
})
