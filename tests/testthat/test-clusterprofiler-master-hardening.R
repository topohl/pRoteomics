make_master_test_manifest <- function(status = "success_with_terms", n_terms = 1L) {
  data.frame(
    result_type = "GSEA_GO", analysis_status = status, n_terms = as.integer(n_terms),
    stringsAsFactors = FALSE
  )
}

make_master_test_result <- function(comparison, manifest_status = "success_with_terms",
                                    n_terms = 1L, worker_status = "SUCCESS", error = NA_character_) {
  list(
    status = worker_status,
    comparison = comparison,
    error = error,
    manifest = make_master_test_manifest(manifest_status, n_terms)
  )
}

testthat::test_that("master assessment distinguishes successful term and zero-term workers", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
  results <- list(
    make_master_test_result("with_terms"),
    make_master_test_result("zero_terms", "success_zero_terms", 0L)
  )
  assessment <- assess_clusterprofiler_worker_results(results, c("with_terms", "zero_terms"))
  counts <- clusterprofiler_master_status_counts(assessment)
  testthat::expect_equal(unname(counts), c(1L, 1L, 0L))
  testthat::expect_equal(clusterprofiler_master_exit_status(assessment), 0L)
  testthat::expect_equal(
    tail(clusterprofiler_master_summary_lines(assessment), 1),
    "ALL COMPARISONS COMPLETED SUCCESSFULLY."
  )
})

testthat::test_that("one failed worker produces a nonzero truthful master status", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
  results <- list(
    make_master_test_result("ok"),
    make_master_test_result("bad", "failed", NA_integer_, "ERROR", "worker failed")
  )
  assessment <- assess_clusterprofiler_worker_results(results, c("ok", "bad"))
  counts <- clusterprofiler_master_status_counts(assessment)
  lines <- clusterprofiler_master_summary_lines(assessment)
  testthat::expect_equal(unname(counts), c(1L, 0L, 1L))
  testthat::expect_equal(clusterprofiler_master_exit_status(assessment), 1L)
  testthat::expect_true("RUN FAILED: 1 comparison(s) failed." %in% lines)
  testthat::expect_false(any(grepl("ALL COMPARISONS COMPLETED SUCCESSFULLY", lines, fixed = TRUE)))
})

testthat::test_that("all failed and malformed worker returns are counted as failed", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
  failed <- list(
    make_master_test_result("a", "failed", NA_integer_, "FAILED", "a failed"),
    make_master_test_result("b", "failed", NA_integer_, "ERROR", "b failed")
  )
  all_failed <- assess_clusterprofiler_worker_results(failed, c("a", "b"))
  testthat::expect_equal(clusterprofiler_master_status_counts(all_failed)[["failed"]], 2L)
  testthat::expect_equal(clusterprofiler_master_exit_status(all_failed), 1L)

  malformed <- assess_clusterprofiler_worker_results(
    list(make_master_test_result("a"), list(status = "SUCCESS")),
    c("a", "b")
  )
  testthat::expect_equal(malformed$analysis_status, c("success_with_terms", "failed"))
  testthat::expect_match(malformed$error[[2]], "missing required fields")

  missing <- assess_clusterprofiler_worker_results(list(make_master_test_result("a")), c("a", "b"))
  testthat::expect_equal(missing$analysis_status, c("success_with_terms", "failed"))
  testthat::expect_match(missing$error[[2]], "missing or malformed")
})

testthat::test_that("output path preflight rejects long roots and accepts P drive paths", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
  excessive <- paste0("C:/", strrep("nested/", 40), "output.csv")
  testthat::expect_error(
    validate_clusterprofiler_output_path_lengths(excessive, safe_limit = 240L),
    "shorter project root such as P:\\\\"
  )
  safe <- validate_clusterprofiler_output_path_lengths(
    "P:/results/tables/04_differential_expression_enrichment/clusterProfiler/microglia/output.csv",
    safe_limit = 240L
  )
  testthat::expect_true(safe$within_safe_limit)
  testthat::expect_lt(safe$path_length, safe$safe_limit)
})

testthat::test_that("strict CSV writing fails loudly and verifies successful writes", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
  root <- tempfile("strict-manifest-")
  dir.create(root)
  path <- file.path(root, "manifest.csv")
  testthat::expect_silent(write_csv_strict(data.frame(a = 1L), path, "test manifest"))
  testthat::expect_true(file.exists(path))
  testthat::expect_error(
    write_csv_strict(data.frame(a = 1L), file.path(root, "missing", "manifest.csv"), "test manifest"),
    "parent directory does not exist"
  )
})

testthat::test_that("child R process exits nonzero for failed assessment", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  helper <- normalizePath(testthat::test_path("..", "..", "R", "enrichment_io.R"), winslash = "/", mustWork = TRUE)
  script <- tempfile("clusterprofiler-exit-", fileext = ".R")
  writeLines(c(
    paste0("source(", dQuote(helper), ")"),
    "assessment <- data.frame(analysis_status = 'failed', stringsAsFactors = FALSE)",
    "quit(status = clusterprofiler_master_exit_status(assessment), save = 'no')"
  ), script)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  output <- suppressWarnings(system2(rscript, script, stdout = TRUE, stderr = TRUE))
  testthat::expect_equal(attr(output, "status"), 1L)
})

testthat::test_that("active master flow preserves failed manifests and uses computed exit status", {
  script <- paste(readLines(
    testthat::test_path("..", "..", "04_differential_expression_enrichment", "01_clusterProfiler.r"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_match(script, "make_failed_manifest_row\\(")
  testthat::expect_match(script, "write_csv_strict\\(manifest, manifest_file")
  testthat::expect_match(script, "quit\\(status = master_exit_status")
  testthat::expect_false(grepl("cat\\(\"ALL COMPARISONS COMPLETED!", script))
})
