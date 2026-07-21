testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))

stage05_label_test_environment <- function(global_dataset = "wrong_dataset") {
  env <- new.env(parent = globalenv())
  env$DATASET <- global_dataset
  expressions <- parse(testthat::test_path("..", "..", "06_modules_WGCNA", "05_module_supermodule_group_effects.r"))
  for (expr in expressions) {
    is_function_assignment <- is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      length(expr) >= 3L && is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))
    if (is_function_assignment) eval(expr, envir = env)
  }
  env
}

stage01_label_fixture <- function() {
  list(
    supermodule_contents = data.frame(
      dataset = "microglia", supermodule_id = "SM01",
      supermodule_label = paste("SM01", intToUtf8(183), "Stage 01 contents program"),
      stringsAsFactors = FALSE
    ),
    comp = data.frame(
      dataset = "microglia", supermodule_id = "SM03",
      supermodule_label = paste("SM03", intToUtf8(183), "Stage 01 composition program"),
      stringsAsFactors = FALSE
    ),
    super_ann = data.frame(
      dataset = "microglia",
      SupermoduleID = c("SM01", "SM01", "SM03", "SM03", "SM03"),
      Supermodule_DisplayLabel = c(rep("Stage 01 annotation one", 2), rep("Stage 01 annotation three", 3)),
      stringsAsFactors = FALSE
    ),
    super_endpoint_map = data.frame(
      dataset = "microglia", endpoint_id = c("SM01", "SM03"),
      endpoint_label = c("Endpoint one", "Endpoint three"), stringsAsFactors = FALSE
    )
  )
}

build_stage05_lookup <- function(env, fixture = stage01_label_fixture(), dataset = "microglia") {
  env$canonical_supermodule_label_lookup(
    fixture$supermodule_contents, fixture$comp, fixture$super_ann,
    fixture$super_endpoint_map, dataset = dataset
  )
}

testthat::test_that("Stage 05 has no Stage 07 or interpretable-summary dependency", {
  script <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "05_module_supermodule_group_effects.r"), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("interpretable_summary", script, fixed = TRUE))
  testthat::expect_false(grepl("WGCNA_final_label_lookup", script, fixed = TRUE))
  testthat::expect_false(grepl("read_final_label_lookup", script, fixed = TRUE))
})

testthat::test_that("Stage 05 diagnostic labels retain explicit compound keys", {
  env <- stage05_label_test_environment()
  lookup <- build_stage05_lookup(env)
  testthat::expect_identical(names(lookup), names(env$empty_canonical_supermodule_label_lookup()))
  testthat::expect_equal(nrow(lookup), 2L)
  testthat::expect_true(all(lookup$dataset == "microglia"))
  testthat::expect_false(anyDuplicated(lookup[c("dataset", "supermodule_id")]) > 0L)
  testthat::expect_match(lookup$canonical_supermodule_label[lookup$supermodule_id == "SM01"], "Stage 01 contents program")
  testthat::expect_match(lookup$canonical_supermodule_label[lookup$supermodule_id == "SM03"], "Stage 01 composition program")
  testthat::expect_true(all(!lookup$lookup_present))
})

testthat::test_that("Stage 05 label application is many-to-one and order preserving", {
  env <- stage05_label_test_environment()
  lookup <- build_stage05_lookup(env)
  result <- data.frame(
    dataset = "microglia", supermodule_id = c("SM03", "SM01", "SM03", "SM01"),
    row_marker = letters[1:4], supermodule_label = "old", stringsAsFactors = FALSE
  )
  applied <- env$apply_canonical_supermodule_labels(result, lookup)
  testthat::expect_equal(nrow(applied), nrow(result))
  testthat::expect_identical(applied$row_marker, result$row_marker)
  testthat::expect_true(all(nzchar(applied$canonical_supermodule_label)))
  testthat::expect_false(anyDuplicated(names(applied)) > 0L)
})

testthat::test_that("Stage 05 rejects missing or duplicate stable keys", {
  env <- stage05_label_test_environment()
  lookup <- build_stage05_lookup(env)
  result <- data.frame(dataset = "microglia", supermodule_id = "SM01")
  testthat::expect_error(
    env$apply_canonical_supermodule_labels(result, lookup[, setdiff(names(lookup), "dataset")]),
    "missing: dataset"
  )
  testthat::expect_error(
    env$apply_canonical_supermodule_labels(result, dplyr::bind_rows(lookup, lookup[1, ])),
    "duplicated dataset + supermodule_id keys", fixed = TRUE
  )
  testthat::expect_error(
    env$apply_canonical_supermodule_labels(data.frame(dataset = "microglia", supermodule_id = "SM09"), lookup),
    "missing result-table key"
  )
})

testthat::test_that("Stage 05 keeps the same supermodule ID separate by dataset", {
  env <- stage05_label_test_environment()
  one_lookup <- function(dataset, label) {
    env$canonical_supermodule_label_lookup(
      data.frame(dataset = dataset, supermodule_id = "SM01", supermodule_label = label),
      data.frame(), data.frame(), data.frame(), dataset = dataset
    )
  }
  lookup <- dplyr::bind_rows(one_lookup("microglia", "Microglia diagnostic"), one_lookup("neuron_soma", "Soma diagnostic"))
  result <- data.frame(dataset = c("neuron_soma", "microglia"), supermodule_id = "SM01", row_marker = c("soma", "microglia"))
  applied <- env$apply_canonical_supermodule_labels(result, lookup)
  testthat::expect_match(applied$canonical_supermodule_label[[1]], "Soma diagnostic")
  testthat::expect_match(applied$canonical_supermodule_label[[2]], "Microglia diagnostic")
  testthat::expect_identical(applied$row_marker, result$row_marker)
})

testthat::test_that("nine Stage 05 supermodules apply without row multiplication", {
  env <- stage05_label_test_environment()
  ids <- sprintf("SM%02d", 1:9)
  lookup <- env$canonical_supermodule_label_lookup(
    data.frame(dataset = "microglia", supermodule_id = ids, supermodule_label = paste("Diagnostic program", ids)),
    data.frame(), data.frame(), data.frame(), dataset = "microglia"
  )
  result <- data.frame(dataset = "microglia", supermodule_id = rep(ids, each = 3L), row_marker = seq_len(27))
  applied <- env$apply_canonical_supermodule_labels(result, lookup)
  testthat::expect_equal(nrow(lookup), 9L)
  testthat::expect_equal(nrow(applied), 27L)
  testthat::expect_identical(applied$row_marker, result$row_marker)
})
