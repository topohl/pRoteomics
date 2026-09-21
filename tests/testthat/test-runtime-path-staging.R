# Runtime manifest resolution with selective staged reads.
#
# The stored enrichment manifests are immutable provenance recorded under a
# substituted P:/ root. Re-anchoring their short stored suffixes on this
# 69-character repository root is necessary but adds ~66 characters, which
# pushes part of the set past the 260-character wall R can open. So a declared
# path keeps its recorded value and a runtime path is derived beside it, with a
# byte-identical copy staged only where the re-anchored path actually crosses
# the wall and the field is actually read.
#
# These tests pin: the two dimensions never overwrite each other, staging is
# addressability rather than transformation, and nothing unread is staged.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "enrichment", "enrichment_io.R"))

has_pwsh <- function() nzchar(Sys.which("pwsh"))

# A directory whose own path is long but still creatable, so that appending a
# basename crosses the wall. R cannot create a directory at or past the wall
# either, so this stops short of it deliberately.
deep_dir <- function(parent, target = 245L) {
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  d <- parent
  repeat {
    remaining <- target - path_length_chars(d)
    if (remaining <= 2L) break
    nxt <- file.path(d, strrep("d", min(30L, remaining - 1L)))
    if (!dir.create(nxt, recursive = FALSE, showWarnings = FALSE)) break
    d <- nxt
  }
  d
}

testthat::test_that("the field map classifies every path column and no prose column", {
  req <- clusterprofiler_runtime_required_fields()
  prov <- clusterprofiler_provenance_only_fields()
  testthat::expect_setequal(req, c("output_table", "collapsed_gene_input_file",
                                   "term_gene_provenance_file", "input_gene_file"))
  testthat::expect_true("collapsed_gene_provenance_file" %in% prov)
  testthat::expect_length(intersect(req, prov), 0L)
  # gene_mapping_policy holds prose, not a path, and must not be resolved
  testthat::expect_false("gene_mapping_policy" %in% names(CLUSTERPROFILER_MANIFEST_PATH_FIELDS))
})

testthat::test_that("runtime resolution states are a closed vocabulary", {
  testthat::expect_setequal(RUNTIME_RESOLUTION_LEVELS,
    c("direct", "rebased", "staged", "unresolved", "not_consumed"))
  testthat::expect_length(intersect(RUNTIME_RESOLUTION_LEVELS, INPUT_STATUS_LEVELS), 0L)
})

testthat::test_that("a short rebase is used directly and nothing is staged", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "sub"), recursive = TRUE, showWarnings = FALSE)
  src <- file.path(tmp, "sub", "short_input.csv")
  writeLines("a,b", src)

  res <- resolve_runtime_paths(paste0("P://sub/short_input.csv"),
                               repository_root = tmp,
                               stage_root = file.path(tmp, "stage"))
  testthat::expect_identical(res$runtime_resolution, RUNTIME_RESOLUTION_REBASED)
  testthat::expect_identical(res$declared_status, INPUT_STATUS_ROOT_UNMOUNTED)
  testthat::expect_identical(normalizePath(res$runtime_path, winslash = "/", mustWork = FALSE),
                             normalizePath(src, winslash = "/", mustWork = FALSE))
  testthat::expect_false(dir.exists(file.path(tmp, "stage")))
})

testthat::test_that("a declared path that already works resolves as direct", {
  tmp <- withr::local_tempdir()
  src <- file.path(tmp, "already_usable.csv")
  writeLines("a,b", src)
  res <- resolve_runtime_paths(src, repository_root = tmp,
                               stage_root = file.path(tmp, "stage"))
  testthat::expect_identical(res$runtime_resolution, RUNTIME_RESOLUTION_DIRECT)
  testthat::expect_identical(res$declared_status, INPUT_STATUS_PRESENT)
  testthat::expect_identical(res$runtime_path, src)
})

testthat::test_that("a genuinely absent input stays unresolved and is never staged", {
  tmp <- withr::local_tempdir()
  res <- resolve_runtime_paths("P://sub/not_there.csv", repository_root = tmp,
                               stage_root = file.path(tmp, "stage"))
  testthat::expect_identical(res$runtime_resolution, RUNTIME_RESOLUTION_UNRESOLVED)
  testthat::expect_false(dir.exists(file.path(tmp, "stage")))
})

testthat::test_that("the declared status and the runtime resolution are independent", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "sub"), recursive = TRUE, showWarnings = FALSE)
  writeLines("a,b", file.path(tmp, "sub", "x.csv"))
  res <- resolve_runtime_paths("P://sub/x.csv", repository_root = tmp,
                               stage_root = file.path(tmp, "stage"))
  # provenance truth is preserved even though the runtime path works
  testthat::expect_identical(res$declared_status, INPUT_STATUS_ROOT_UNMOUNTED)
  testthat::expect_identical(res$runtime_resolution, RUNTIME_RESOLUTION_REBASED)
  # the declared path itself is untouched
  testthat::expect_identical(res$declared_path, "P://sub/x.csv")
  testthat::expect_false(identical(res$declared_path, res$runtime_path))
})

testthat::test_that("staged destinations are deterministic, collision-free and keep the extension", {
  root <- "S:/stage_root"
  a <- "P://one/deep/place/gsea_kegg_term_gene_provenance.csv"
  b <- "P://another/deep/place/gsea_kegg_term_gene_provenance.csv"
  da <- staged_destination(a, root)
  db <- staged_destination(b, root)

  testthat::expect_identical(da, staged_destination(a, root))   # deterministic
  testthat::expect_false(identical(da, db))                     # same basename, different source
  testthat::expect_identical(basename(da), basename(a))          # extension verbatim
  testthat::expect_match(da, "[.]csv$")
  for (ext in c("csv", "rds", "xlsx", "tsv", "svg")) {
    p <- staged_destination(paste0("P://x/y/object.", ext), root)
    testthat::expect_identical(tools::file_ext(p), ext)
  }
  # no deep hierarchy reproduced: one digest directory, then the basename
  testthat::expect_identical(length(strsplit(sub(paste0("^", root, "/"), "", da), "/")[[1]]), 2L)
})

testthat::test_that("an over-wall source is staged byte-identically and becomes usable", {
  testthat::skip_if_not(has_pwsh(), "pwsh not available for extended-length staging")
  tmp <- withr::local_tempdir()
  deep <- deep_dir(file.path(tmp, "repo"))
  testthat::skip_if_not(dir.exists(deep), "could not build a deep fixture directory")
  base <- "a_source_table_beyond_the_wall.csv"
  src <- file.path(deep, base)
  testthat::skip_if(path_length_chars(src) < PATH_LENGTH_WALL, "fixture did not cross the wall")

  # staged via the OS, because R cannot create this path
  seed <- file.path(tmp, "seed.csv")
  writeLines(c("gene,value", "Aqp4,1.5"), seed)
  testthat::expect_true(all(os_copy_no_clobber(seed, src)))
  testthat::expect_false(file.exists(src))             # the pathology
  facts <- os_path_facts(src)
  testthat::expect_true(facts$exists)

  stage_root <- file.path(tmp, "stg")
  staged <- stage_addressable_copies(src, staged_destination(src, stage_root),
                                     expected_sha256 = facts$sha256)
  testthat::expect_identical(staged$action, "copied")
  testthat::expect_true(staged$same_hash)
  testthat::expect_true(staged$same_size)
  testthat::expect_true(staged$source_unchanged)
  testthat::expect_true(staged$ok)
  testthat::expect_lt(path_length_chars(staged$staged_path), PATH_LENGTH_WALL)
  testthat::expect_true(file.exists(staged$staged_path))

  # scientific equivalence: the staged copy parses to the same data
  direct <- utils::read.csv(seed, stringsAsFactors = FALSE)
  through_stage <- utils::read.csv(staged$staged_path, stringsAsFactors = FALSE)
  testthat::expect_identical(direct, through_stage)
})

testthat::test_that("staging is idempotent and refuses to accept stale content", {
  testthat::skip_if_not(has_pwsh(), "pwsh not available")
  tmp <- withr::local_tempdir()
  src <- file.path(tmp, "src.csv")
  writeLines(c("a,b", "1,2"), src)
  dst <- staged_destination("P://x/src.csv", file.path(tmp, "stg"))
  facts <- os_path_facts(src)

  first <- stage_addressable_copies(src, dst, expected_sha256 = facts$sha256)
  testthat::expect_identical(first$action, "copied")
  testthat::expect_true(first$ok)

  second <- stage_addressable_copies(src, dst, expected_sha256 = facts$sha256)
  testthat::expect_identical(second$action, "reused")
  testthat::expect_true(second$ok)
  testthat::expect_identical(second$staged_sha256, first$staged_sha256)

  # a destination holding different bytes is a mismatch, not a silent overwrite
  writeLines(c("a,b", "9,9"), dst)
  stale <- stage_addressable_copies(src, dst, expected_sha256 = facts$sha256)
  testthat::expect_identical(stale$action, "mismatch")
  testthat::expect_false(stale$ok)
  testthat::expect_false(stale$same_hash)
  testthat::expect_identical(readLines(dst)[2], "9,9")   # not overwritten
})

testthat::test_that("the 259/260 boundary decides whether staging is attempted", {
  tmp <- withr::local_tempdir()
  # below the wall R answers for itself and no OS probe or staging is needed
  dir.create(file.path(tmp, "s"), recursive = TRUE, showWarnings = FALSE)
  writeLines("x", file.path(tmp, "s", "f.csv"))
  res <- resolve_runtime_paths("P://s/f.csv", repository_root = tmp,
                               stage_root = file.path(tmp, "stg"))
  testthat::expect_lt(res$candidate_length, PATH_LENGTH_WALL)
  testthat::expect_identical(res$runtime_resolution, RUNTIME_RESOLUTION_REBASED)
  testthat::expect_true(is.na(res$staged_sha256))
  testthat::expect_identical(PATH_LENGTH_WALL, 260L)
})

# --- manifest level ---------------------------------------------------------

make_runtime_manifest <- function(paths) {
  data.frame(
    dataset = "microglia", comparison = "c1", result_type = "GSEA_GO",
    ontology = "BP", analysis_status = "success_with_terms", n_terms = 1L,
    output_table = paths$output_table,
    collapsed_gene_input_file = paths$collapsed_gene_input_file,
    collapsed_gene_provenance_file = paths$collapsed_gene_provenance_file,
    term_gene_provenance_file = paths$term_gene_provenance_file,
    enrichment_contract_version = canonical_clusterprofiler_manifest_contract_version(),
    gene_annotation_contract_version = "v1",
    stringsAsFactors = FALSE)
}

testthat::test_that("the stored manifest is never mutated by resolution", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "d"), recursive = TRUE, showWarnings = FALSE)
  for (n in c("t", "c", "p", "g")) writeLines("a,b", file.path(tmp, "d", paste0(n, ".csv")))
  declared <- list(output_table = "P://d/t.csv",
                   collapsed_gene_input_file = "P://d/c.csv",
                   collapsed_gene_provenance_file = "P://d/p.csv",
                   term_gene_provenance_file = "P://d/g.csv")
  m <- make_runtime_manifest(declared)
  rr <- clusterprofiler_manifest_runtime_resolution(m, dataset = "microglia",
                                                    repository_root = tmp, stage = FALSE)

  # declared values survive verbatim in the *_declared columns
  for (f in clusterprofiler_runtime_required_fields()) {
    if (!f %in% names(m)) next
    testthat::expect_identical(rr$manifest[[paste0(f, "_declared")]], m[[f]])
  }
  # and the provenance-only column keeps its declared value in place
  testthat::expect_identical(rr$manifest$collapsed_gene_provenance_file,
                             declared$collapsed_gene_provenance_file)
})

testthat::test_that("a provenance-only field is not consumed, not gated and not staged", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "d"), recursive = TRUE, showWarnings = FALSE)
  for (n in c("t", "c", "g")) writeLines("a,b", file.path(tmp, "d", paste0(n, ".csv")))
  declared <- list(output_table = "P://d/t.csv",
                   collapsed_gene_input_file = "P://d/c.csv",
                   # deliberately absent on disk AND unmounted: must not matter
                   collapsed_gene_provenance_file = "P://d/never_read.csv",
                   term_gene_provenance_file = "P://d/g.csv")
  m <- make_runtime_manifest(declared)
  rr <- clusterprofiler_manifest_runtime_resolution(m, dataset = "microglia",
                                                    repository_root = tmp, stage = FALSE)
  prov <- rr$resolution[rr$resolution$field == "collapsed_gene_provenance_file", ]
  testthat::expect_identical(prov$runtime_resolution, RUNTIME_RESOLUTION_NOT_CONSUMED)
  testthat::expect_true(is.na(prov$staged_sha256))

  # the validator passes: runtime-required fields are all usable
  testthat::expect_true(validate_clusterprofiler_manifest_contract(
    rr$manifest, require_files = TRUE))
})

testthat::test_that("the validator stays all-or-nothing over runtime-required fields", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "d"), recursive = TRUE, showWarnings = FALSE)
  for (n in c("t", "c", "p", "g")) writeLines("a,b", file.path(tmp, "d", paste0(n, ".csv")))
  base <- list(output_table = "P://d/t.csv",
               collapsed_gene_input_file = "P://d/c.csv",
               collapsed_gene_provenance_file = "P://d/p.csv",
               term_gene_provenance_file = "P://d/g.csv")

  # one unusable runtime-required field invalidates the whole manifest
  for (f in c("output_table", "collapsed_gene_input_file", "term_gene_provenance_file")) {
    broken <- base
    broken[[f]] <- "P://d/missing.csv"
    rr <- clusterprofiler_manifest_runtime_resolution(
      make_runtime_manifest(broken), dataset = "microglia",
      repository_root = tmp, stage = FALSE)
    testthat::expect_error(
      validate_clusterprofiler_manifest_contract(rr$manifest, require_files = TRUE),
      f)
  }
})

testthat::test_that("the real manifests resolve and validate for every dataset", {
  mans <- Sys.glob(repo_path("data", "processed",
                             "04_differential_expression_enrichment",
                             "clusterProfiler", "*", "clusterProfiler_manifest.csv"))
  testthat::skip_if(!length(mans), "no on-disk clusterProfiler manifests")
  testthat::skip_if_not(has_pwsh(), "pwsh not available for extended-length staging")

  for (m in mans) {
    ds <- basename(dirname(m))
    before <- unname(tools::sha256sum(m))
    mm <- read_canonical_clusterprofiler_manifest(m, dataset = ds, strict = TRUE,
                                                  require_files = TRUE)
    testthat::expect_gt(nrow(mm), 0L)

    # the stored manifest file is byte-identical afterwards
    testthat::expect_identical(unname(tools::sha256sum(m)), before, info = ds)

    res <- attr(mm, "runtime_resolution")
    testthat::expect_false(is.null(res), info = ds)
    req <- res[res$field_kind == "runtime_required", ]
    testthat::expect_identical(sum(req$runtime_resolution == RUNTIME_RESOLUTION_UNRESOLVED),
                               0L, info = ds)
    # every runtime path a consumer will open is under the wall
    rt <- req$runtime_path[req$runtime_resolution != RUNTIME_RESOLUTION_UNRESOLVED]
    testthat::expect_true(all(path_length_chars(rt) < PATH_LENGTH_WALL), info = ds)
    # declared provenance is untouched: still recorded under the unmounted root
    testthat::expect_true(all(grepl("^[Pp]:", req$declared_path)), info = ds)
    testthat::expect_true(all(req$declared_status == INPUT_STATUS_ROOT_UNMOUNTED), info = ds)
    # nothing unread was staged
    prov <- res[res$field_kind == "provenance_only", ]
    testthat::expect_identical(sum(prov$runtime_resolution == RUNTIME_RESOLUTION_STAGED),
                               0L, info = ds)
  }
})

testthat::test_that("staged files are readable from a multisession worker", {
  testthat::skip_if_not(has_pwsh(), "pwsh not available")
  testthat::skip_if_not(requireNamespace("future", quietly = TRUE) &&
                        requireNamespace("future.apply", quietly = TRUE),
                        "future not available")
  tmp <- withr::local_tempdir()
  src <- file.path(tmp, "worker_src.csv")
  writeLines(c("gene,value", "Gfap,2"), src)
  dst <- staged_destination("P://w/worker_src.csv", file.path(tmp, "stg"))
  staged <- stage_addressable_copies(src, dst,
                                     expected_sha256 = os_path_facts(src)$sha256)
  testthat::expect_true(staged$ok)

  suppressWarnings(suppressPackageStartupMessages({
    library(future); library(future.apply)
  }))
  old <- suppressWarnings(future::plan(future::multisession, workers = 2))
  on.exit(suppressWarnings(future::plan(old)), add = TRUE)
  seen <- suppressWarnings(future.apply::future_lapply(seq_len(2), function(i) {
    list(exists = file.exists(dst),
         sha = if (file.exists(dst)) unname(tools::sha256sum(dst)) else NA_character_,
         rows = if (file.exists(dst)) nrow(utils::read.csv(dst)) else NA_integer_)
  }))
  for (s in seen) {
    testthat::expect_true(s$exists)
    testthat::expect_identical(tolower(s$sha), tolower(staged$staged_sha256))
    testthat::expect_identical(s$rows, 1L)
  }
  # the staging root is reached by the same absolute path in both, so this does
  # not depend on any mapped-drive inheritance
  testthat::expect_true(file.exists(dst))
})
