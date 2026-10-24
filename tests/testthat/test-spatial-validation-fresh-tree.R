source(testthat::test_path("..", "..", "R", "paths.R"))

# A regeneration tree without the pre-restructure results exposed two defects
# the main tree hid: tied CA2-SLM ranks printed with %d, and a declared output
# written only when the legacy 07_spatial_networks tree existed. Static checks
# only; neither script is run.

sv_file <- function(aid) {
  repo_path("analysis", "spatial_validation", paste0(aid, ".R"))
}
code_of <- function(f) {
  l <- sub("#.*$", "", readLines(f, warn = FALSE))
  paste(l[nzchar(trimws(l))], collapse = "\n")
}
# every data.frame(...) call in a script, unevaluated
df_calls <- function(f) {
  out <- list()
  walk <- function(e) {
    ok <- tryCatch(is.call(e), error = function(...) FALSE)
    if (!ok) return(invisible(NULL))
    if (identical(e[[1]], as.name("data.frame"))) out[[length(out) + 1L]] <<- e
    for (i in seq_along(e)[-1]) tryCatch(walk(e[[i]]), error = function(...) NULL)
  }
  for (x in parse(f, keep.source = FALSE)) walk(x)
  out
}
arg_names <- function(cl) setdiff(names(as.list(cl))[-1], c("", "stringsAsFactors"))

testthat::test_that("CA2-SLM console ranks are integers when the SEs tie", {
  code <- code_of(sv_file("quantify_neuropil_precision"))
  # rank()'s default 'average' ties give 3.5, which sprintf('%d') rejects
  testthat::expect_false(grepl("rank(ctx$", code, fixed = TRUE))
  testthat::expect_true(grepl('ties.method = "min"', code, fixed = TRUE))
  testthat::expect_true(grepl("rk(signif(ctx$median_model_SE, 8))", code, fixed = TRUE))
  testthat::expect_true(grepl("rk(ctx$missing_fraction_preimputation)", code, fixed = TRUE))
  # P1 compares SEs that are equal by construction; noise must not decide it
  testthat::expect_true(grepl("signif(ca2$median_model_SE, 8) <", code, fixed = TRUE))
  se <- 0.15994523864897 + c(-6, 0, 1, 3, -1, 0, 2, 1, -2, 0) * 1e-16
  testthat::expect_gt(length(unique(se)), 1L)
  r <- rank(signif(se, 8), ties.method = "min")
  testthat::expect_type(r, "integer")
  testthat::expect_identical(sprintf("%d of %d", r[1], length(se)), "1 of 10")
  testthat::expect_false(signif(se[1], 8) < stats::median(signif(se[-1], 8)))
})

testthat::test_that("the legacy network comparison is always written with its fixed schema", {
  f <- sv_file("validate_network_workbook")
  code <- code_of(f)
  # declared output: written unconditionally, header-only without a legacy tree
  testthat::expect_true(grepl(
    '\nwrite_csv_safe(legacy, NET_OUT("legacy_spatial_network_comparison.csv"))',
    code, fixed = TRUE))
  testthat::expect_true(grepl("legacy 07_spatial_networks tree absent", code, fixed = TRUE))
  # the zero-row schema and the row builder carry the same columns in order
  calls <- Filter(function(cl) "legacy_file" %in% arg_names(cl), df_calls(f))
  testthat::expect_length(calls, 2L)
  testthat::expect_identical(arg_names(calls[[1]]), arg_names(calls[[2]]))
  schema <- Filter(Negate(is.null), lapply(calls, function(cl)
    tryCatch(eval(cl, envir = baseenv()), error = function(...) NULL)))
  testthat::expect_length(schema, 1L)
  schema <- schema[[1]]
  testthat::expect_equal(nrow(schema), 0L)
  testthat::expect_identical(names(schema), c(
    "legacy_file", "legacy_edge", "legacy_estimate", "legacy_permutation_p",
    "legacy_fdr", "new_animal_level_estimate", "comparable",
    "legacy_inference_status", "legacy_defect", "interpretation"))
  testthat::skip_if_not_installed("readr")
  p <- tempfile(fileext = ".csv")
  on.exit(unlink(p), add = TRUE)
  readr::write_csv(schema, p, na = "")
  back <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  testthat::expect_equal(nrow(back), 0L)
  testthat::expect_identical(names(back), names(schema))
})
