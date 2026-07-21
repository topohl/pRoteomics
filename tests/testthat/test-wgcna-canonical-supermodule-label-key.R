testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")
testthat::skip_if_not_installed("stringr")

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))

canonical_label_test_environment <- function(global_dataset = "wrong_dataset") {
  env <- new.env(parent = globalenv())
  env$DATASET <- global_dataset
  expressions <- parse(testthat::test_path(
    "..", "..", "06_modules_WGCNA", "05_module_supermodule_group_effects.r"
  ))
  for (expr in expressions) {
    is_function_assignment <- is.call(expr) &&
      identical(expr[[1]], as.name("<-")) &&
      length(expr) >= 3L && is.call(expr[[3]]) &&
      identical(expr[[3]][[1]], as.name("function"))
    if (is_function_assignment) eval(expr, envir = env)
  }
  env
}

canonical_label_fixture <- function() {
  final_labels <- data.frame(
    dataset = "microglia",
    supermodule_id = "SM01",
    final_plot_label = paste("SM01", intToUtf8(183), "Final program one"),
    final_plot_label_short = "SM01\nFinal program one",
    final_label_source = "test_final_lookup",
    stringsAsFactors = FALSE
  )
  attr(final_labels, "lookup_present") <- TRUE
  list(
    supermodule_contents = data.frame(
      dataset = "microglia",
      supermodule_id = "SM01",
      supermodule_label = paste("SM01", intToUtf8(183), "Contents program one"),
      stringsAsFactors = FALSE
    ),
    comp = data.frame(
      dataset = "microglia",
      supermodule_id = "SM03",
      supermodule_label = paste("SM03", intToUtf8(183), "Composition program three"),
      stringsAsFactors = FALSE
    ),
    super_ann = data.frame(
      dataset = "microglia",
      SupermoduleID = c("SM01", "SM01", "SM03", "SM03", "SM03"),
      Supermodule_DisplayLabel = c(
        rep(paste("SM01", intToUtf8(183), "Annotation program one"), 2),
        rep(paste("SM03", intToUtf8(183), "Annotation program three"), 3)
      ),
      stringsAsFactors = FALSE
    ),
    super_endpoint_map = data.frame(
      dataset = "microglia",
      endpoint_id = c("SM01", "SM03"),
      endpoint_label = c("Endpoint one", "Endpoint three"),
      stringsAsFactors = FALSE
    ),
    final_labels = final_labels
  )
}

build_canonical_fixture_lookup <- function(env, fixture = canonical_label_fixture(), dataset = "microglia") {
  env$canonical_supermodule_label_lookup(
    fixture$supermodule_contents,
    fixture$comp,
    fixture$super_ann,
    fixture$super_endpoint_map,
    dataset = dataset,
    final_labels = fixture$final_labels
  )
}

testthat::test_that("nonempty canonical lookup retains explicit dataset keys and labels", {
  env <- canonical_label_test_environment(global_dataset = "wrong_dataset")
  lookup <- build_canonical_fixture_lookup(env)

  testthat::expect_identical(names(lookup), names(env$empty_canonical_supermodule_label_lookup()))
  testthat::expect_equal(nrow(lookup), 2L)
  testthat::expect_true(all(lookup$dataset == "microglia"))
  testthat::expect_false(anyDuplicated(lookup[c("dataset", "supermodule_id")]) > 0L)
  testthat::expect_equal(
    lookup$canonical_supermodule_label[lookup$supermodule_id == "SM01"],
    paste("SM01", intToUtf8(183), "Final program one")
  )
  testthat::expect_equal(
    lookup$canonical_supermodule_label[lookup$supermodule_id == "SM03"],
    paste("SM03", intToUtf8(183), "Composition program three")
  )
  testthat::expect_true(all(nzchar(lookup$canonical_supermodule_plot_label)))
})

testthat::test_that("empty and partially populated sources retain the typed dataset schema", {
  env <- canonical_label_test_environment()
  empty_final <- data.frame(
    dataset = character(), supermodule_id = character(), final_plot_label = character(),
    final_plot_label_short = character(), final_label_source = character()
  )
  attr(empty_final, "lookup_present") <- FALSE
  empty <- env$canonical_supermodule_label_lookup(
    data.frame(), data.frame(), data.frame(), data.frame(),
    dataset = "microglia", final_labels = empty_final
  )
  testthat::expect_equal(nrow(empty), 0L)
  testthat::expect_identical(names(empty), names(env$empty_canonical_supermodule_label_lookup()))
  testthat::expect_type(empty$dataset, "character")

  partial_final <- data.frame(supermodule_id = "SM01", stringsAsFactors = FALSE)
  attr(partial_final, "lookup_present") <- TRUE
  partial <- env$canonical_supermodule_label_lookup(
    data.frame(supermodule_id = "SM01", supermodule_label = "Partial source label"),
    data.frame(), data.frame(), data.frame(), dataset = "microglia",
    final_labels = partial_final
  )
  testthat::expect_equal(nrow(partial), 1L)
  testthat::expect_equal(partial$dataset, "microglia")
  testthat::expect_match(partial$canonical_supermodule_label, "Partial source label")
})

testthat::test_that("canonical application is many-to-one and preserves result order", {
  env <- canonical_label_test_environment()
  lookup <- build_canonical_fixture_lookup(env)
  result <- data.frame(
    dataset = "microglia",
    supermodule_id = c("SM03", "SM01", "SM03", "SM01"),
    row_marker = c("a", "b", "c", "d"),
    supermodule_label = "old label",
    stringsAsFactors = FALSE
  )
  applied <- env$apply_canonical_supermodule_labels(result, lookup)
  testthat::expect_equal(nrow(applied), nrow(result))
  testthat::expect_identical(applied$row_marker, result$row_marker)
  testthat::expect_identical(applied$supermodule_id, result$supermodule_id)
  testthat::expect_true(all(nzchar(applied$canonical_supermodule_label)))
  testthat::expect_true(all(applied$supermodule_label == applied$canonical_supermodule_label))
  testthat::expect_false(anyDuplicated(names(applied)) > 0L)
})

testthat::test_that("invalid canonical lookup key contracts are rejected", {
  env <- canonical_label_test_environment()
  lookup <- build_canonical_fixture_lookup(env)
  result <- data.frame(dataset = "microglia", supermodule_id = "SM01")

  testthat::expect_error(
    env$apply_canonical_supermodule_labels(result, lookup[, setdiff(names(lookup), "dataset")]),
    "missing: dataset"
  )
  duplicated_lookup <- dplyr::bind_rows(lookup, lookup[1, ])
  testthat::expect_error(
    env$apply_canonical_supermodule_labels(result, duplicated_lookup),
    "duplicated dataset + supermodule_id keys",
    fixed = TRUE
  )
  testthat::expect_error(
    env$apply_canonical_supermodule_labels(
      data.frame(dataset = "microglia", supermodule_id = "SM09"),
      lookup
    ),
    "missing result-table key"
  )
})

testthat::test_that("the same supermodule ID remains separate across datasets", {
  env <- canonical_label_test_environment(global_dataset = "wrong_dataset")
  empty_final <- data.frame(
    dataset = character(), supermodule_id = character(), final_plot_label = character(),
    final_plot_label_short = character(), final_label_source = character()
  )
  attr(empty_final, "lookup_present") <- FALSE
  one_lookup <- function(dataset, label) {
    env$canonical_supermodule_label_lookup(
      data.frame(supermodule_id = "SM01", supermodule_label = label),
      data.frame(), data.frame(), data.frame(), dataset = dataset,
      final_labels = empty_final
    )
  }
  lookup <- dplyr::bind_rows(
    one_lookup("microglia", "Microglia program"),
    one_lookup("neuron_soma", "Soma program")
  )
  result <- data.frame(
    dataset = c("neuron_soma", "microglia"),
    supermodule_id = "SM01",
    row_marker = c("soma", "microglia"),
    stringsAsFactors = FALSE
  )
  applied <- env$apply_canonical_supermodule_labels(result, lookup)
  testthat::expect_identical(applied$row_marker, result$row_marker)
  testthat::expect_match(applied$canonical_supermodule_label[[1]], "Soma program")
  testthat::expect_match(applied$canonical_supermodule_label[[2]], "Microglia program")
  testthat::expect_identical(applied$dataset, result$dataset)
})

testthat::test_that("nine-supermodule result shape applies without multiplication", {
  env <- canonical_label_test_environment()
  ids <- sprintf("SM%02d", 1:9)
  empty_final <- data.frame(
    dataset = character(), supermodule_id = character(), final_plot_label = character(),
    final_plot_label_short = character(), final_label_source = character()
  )
  attr(empty_final, "lookup_present") <- FALSE
  lookup <- env$canonical_supermodule_label_lookup(
    data.frame(
      dataset = "microglia", supermodule_id = ids,
      supermodule_label = paste("Microglia program", seq_along(ids))
    ),
    data.frame(), data.frame(), data.frame(), dataset = "microglia",
    final_labels = empty_final
  )
  result <- data.frame(
    dataset = "microglia",
    supermodule_id = rep(ids, each = 3L),
    row_marker = seq_len(27),
    stringsAsFactors = FALSE
  )
  applied <- env$apply_canonical_supermodule_labels(result, lookup)
  testthat::expect_equal(nrow(lookup), 9L)
  testthat::expect_equal(nrow(applied), 27L)
  testthat::expect_identical(applied$row_marker, result$row_marker)
  testthat::expect_equal(length(unique(applied$supermodule_id)), 9L)
  testthat::expect_true(all(nzchar(applied$canonical_supermodule_label)))
})
