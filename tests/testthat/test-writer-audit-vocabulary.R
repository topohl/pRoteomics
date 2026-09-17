source(testthat::test_path("..", "..", "R", "paths.R"))

# Phase 6G.3: the writer audit and the consumer enumerator are only as good as
# the vocabularies they recognise, and both were found short.
#
#   * tools/audit_writer_namespaces.R knew base R's write calls but not this
#     repository's own wrappers. Three spatial_validation scripts write
#     workbooks through xlsx_save_valid_workbook(wb, wb_path) where wb_path is
#     path_results("tables", "11_spatial_systems", ...). The legacy path was
#     detected, the write was not, so the "active legacy write sites = 0" gate
#     would have passed with those writes still in place.
#
#   * tools/enumerate_output_consumers.R had no branch for an executable
#     script under audits/, so 212 hits from audits/verify_scientific_contracts.R
#     landed in UNKNOWN and blocked the preflight.
#
# A fixed list of names would rot exactly the way the original did, so the
# first test derives the expectation from the code instead: any function in R/
# that writes, and that analysis/ actually calls, must be known to the audit.

AUDIT <- repo_path("tools", "audit_writer_namespaces.R")
ENUM <- repo_path("tools", "enumerate_output_consumers.R")

is_skippable <- function(x) {
  tryCatch(is.null(x) || (is.symbol(x) && !nzchar(as.character(x))),
           error = function(...) TRUE)
}
walk <- function(e, fn) {
  if (!tryCatch({ e; TRUE }, error = function(...) FALSE)) return(invisible(NULL))
  fn(e)
  if (is.call(e) || is.expression(e) || is.list(e)) {
    for (i in seq_along(e)) {
      el <- tryCatch(e[[i]], error = function(...) NULL)
      if (is_skippable(el)) next
      walk(el, fn)
    }
  }
  invisible(NULL)
}
cname <- function(e) {
  if (!is.call(e)) return(NA_character_)
  fn <- e[[1]]
  if (is.name(fn)) return(as.character(fn))
  if (is.call(fn) && length(fn) == 3L && as.character(fn[[1]]) %in% c("::", ":::")) {
    return(as.character(fn[[3]]))
  }
  NA_character_
}

# the vocabulary the audit actually uses, read from its source
audit_write_calls <- function() {
  exprs <- parse(AUDIT)
  found <- NULL
  for (e in exprs) {
    if (is.call(e) && length(e) == 3L && identical(as.character(e[[1]])[1], "<-") &&
        is.name(e[[2]]) && identical(as.character(e[[2]]), "WRITE_CALLS")) {
      found <- e
    }
  }
  testthat::expect_false(is.null(found))
  ## WRITE_CALLS is assembled from two vectors, so evaluate it in a small env
  env <- new.env(parent = baseenv())
  for (e in parse(AUDIT)) {
    if (is.call(e) && length(e) == 3L && identical(as.character(e[[1]])[1], "<-") &&
        is.name(e[[2]]) &&
        as.character(e[[2]]) %in% c("BASE_WRITE_CALLS", "WRAPPER_WRITE_CALLS", "WRITE_CALLS")) {
      eval(e, envir = env)
    }
  }
  get("WRITE_CALLS", envir = env)
}

testthat::test_that("the writer audit knows every write wrapper analysis/ calls", {
  ## a call that actually puts bytes on disk. unlink() is excluded on purpose:
  ## deleting is not writing to a destination.
  PRIMITIVE <- c("write.csv", "write.table", "writeLines", "saveRDS", "ggsave",
                 "saveWorkbook", "write.xlsx", "file.copy", "file.rename",
                 "write_yaml", "png", "pdf", "svg", "jpeg", "tiff", "cairo_pdf",
                 "dir.create")

  rfiles <- list.files(repo_path("R"), pattern = "[.][Rr]$",
                       recursive = TRUE, full.names = TRUE)
  writers <- character(0)
  for (f in rfiles) {
    exprs <- tryCatch(parse(f), error = function(e) NULL)
    if (is.null(exprs)) next
    for (e in exprs) {
      if (!(is.call(e) && length(e) == 3L &&
            as.character(e[[1]])[1] %in% c("<-", "=") && is.name(e[[2]]))) next
      rhs <- e[[3]]
      if (!(is.call(rhs) && identical(as.character(rhs[[1]])[1], "function"))) next
      calls <- character(0)
      walk(rhs[[3]], function(x) {
        n <- cname(x); if (!is.na(n)) calls <<- c(calls, n)
      })
      if (any(calls %in% PRIMITIVE)) writers <- c(writers, as.character(e[[2]]))
    }
  }
  writers <- unique(writers)
  testthat::expect_gt(length(writers), 0L)

  afiles <- list.files(repo_path("analysis"), pattern = "[.][Rr]$",
                       recursive = TRUE, full.names = TRUE)
  used <- character(0)
  for (f in afiles) {
    exprs <- tryCatch(parse(f), error = function(e) NULL)
    if (is.null(exprs)) next
    for (e in exprs) walk(e, function(x) {
      n <- cname(x)
      if (!is.na(n) && n %in% writers) used <<- c(used, n)
    })
  }
  used <- unique(used)
  testthat::expect_gt(length(used), 0L)

  known <- audit_write_calls()
  missing <- setdiff(used, known)
  testthat::expect_identical(
    missing, character(0),
    info = paste0("analysis/ calls these write wrappers but the writer audit ",
                  "does not recognise them, so writes made through them are ",
                  "invisible to the legacy-write gate: ",
                  paste(missing, collapse = ", ")))
})

testthat::test_that("the wrapper that hid spatial_validation's workbook writes is known", {
  known <- audit_write_calls()
  ## the concrete regression: three spatial_validation scripts write through it
  testthat::expect_true("xlsx_save_valid_workbook" %in% known)
})

testthat::test_that("a baseline-keyed audit is provenance, not a live consumer", {
  ## audits/verify_scientific_contracts.R names frozen carriers by their
  ## pre-restructure path and resolves them through the migration map, so it
  ## must never be repointed by a writer migration.
  f <- repo_path("audits", "verify_scientific_contracts.R")
  testthat::expect_true(file.exists(f))
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("restructure_migration_map", txt, fixed = TRUE))

  ## and the enumerator classifies on that derived property rather than on the
  ## directory, so a genuine new reader under audits/ is not excused
  etxt <- paste(readLines(ENUM, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("baseline_keyed_surfaces", etxt, fixed = TRUE))
  testthat::expect_true(grepl("AUDIT_READER", etxt, fixed = TRUE))
})

testthat::test_that("the spatial_validation inventories classify nothing as UNKNOWN", {
  dir <- repo_path("audits", "consumer_inventory")
  testthat::skip_if_not(dir.exists(dir))
  sv <- c("build_spatial_data_contract", "quantify_bilateral_spatial_identity",
          "quantify_empirical_compartments", "quantify_module_bilateral_identity",
          "annotate_module_celltypes", "decompose_bilateral_variance",
          "validate_spatial_foundations", "build_module_spatial_atlas",
          "build_protein_spatial_atlas", "quantify_neuropil_detection_context",
          "summarize_spatial_atlas", "quantify_neuropil_precision",
          "build_animal_spatial_networks", "test_network_group_organization",
          "validate_network_workbook", "audit_ca2_slm_robustness",
          "audit_stress_identity_robustness", "summarize_ca2_slm_robustness")
  checked <- 0L
  for (a in sv) {
    f <- file.path(dir, paste0(a, ".csv"))
    if (!file.exists(f)) next
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    testthat::expect_identical(sum(d$dependency_kind == "UNKNOWN"), 0L,
                               info = paste(a, "has unclassified dependencies"))
    checked <- checked + 1L
  }
  testthat::expect_gt(checked, 0L)
})
