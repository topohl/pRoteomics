testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")
testthat::skip_if_not_installed("stringr")

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))

endpoint_map_test_environment <- function() {
  env <- new.env(parent = globalenv())
  env$DATASET <- "microglia"
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

endpoint_map_fixture <- function() {
  module_eigengenes <- c("MEturquoise", "MEblue", "MEbrown", "MEyellow", "MEgreen", "MEred")
  module_eig <- data.frame(
    Sample = paste0("sample", 1:4),
    matrix(seq_len(24), nrow = 4, dimnames = list(NULL, module_eigengenes)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  definitions <- data.frame(
    module_eigengene = module_eigengenes,
    ModuleID = sprintf("WGCNA_m%02d", seq_along(module_eigengenes)),
    ModuleColor = c("#486A8A", "#7FA6C1", "#5D7894", "#8AA0AF", "#2F6F73", "#5B9992"),
    ModuleLabel_Final = paste("Module", seq_along(module_eigengenes)),
    stringsAsFactors = FALSE
  )
  supermodule_id <- c("SM01", "SM01", "SM02", "SM03", "SM03", "SM03")
  program <- c("Program alpha", "Program alpha", "Program beta", rep("Program gamma", 3))
  super_ann <- data.frame(
    dataset = "microglia",
    module_eigengene = module_eigengenes,
    ModuleID = definitions$ModuleID,
    Supermodule_DataDrivenID = supermodule_id,
    SupermoduleID = supermodule_id,
    Supermodule_FinalLabel = program,
    Supermodule_DisplayLabel = paste(supermodule_id, intToUtf8(183), program),
    Supermodule_LabelConfidence = c("high", "high", "low", "medium", "medium", "medium"),
    Supermodule_NameSource = c(rep("recurring_significant_GO", 2), "singleton", rep("recurring_significant_GO", 3)),
    ManualReviewRequired = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
    SupermoduleCutHeight = 0.40,
    supermodule_merge_cut_height = 0.40,
    supermodule_merge_rule = "hclust_average_on_1_minus_module_eigengene_correlation",
    label_source = "preexisting_annotation_source",
    present_in_dataset = TRUE,
    stringsAsFactors = FALSE
  )
  list(module_eig = module_eig, definitions = definitions, super_ann = super_ann)
}

testthat::test_that("endpoint maps resolve labels once per supermodule without multiplying modules", {
  env <- endpoint_map_test_environment()
  fixture <- endpoint_map_fixture()
  maps <- env$make_endpoint_maps(fixture$module_eig, fixture$definitions, fixture$super_ann)

  testthat::expect_equal(nrow(maps$super_map), nrow(fixture$super_ann))
  testthat::expect_equal(nrow(maps$super_label_lookup), 3L)
  testthat::expect_false(anyDuplicated(maps$super_label_lookup[c("dataset", "supermodule_id")]) > 0L)
  testthat::expect_false(anyDuplicated(maps$super_map$module_eigengene) > 0L)
  testthat::expect_false(anyDuplicated(names(maps$super_map)) > 0L)
  testthat::expect_true(all(c(
    "Supermodule_LabelSource", "resolved_label_source", "label_source",
    "Supermodule_DisplayLabel", "Supermodule_LabelConfidence",
    "Supermodule_NameSource", "ManualReviewRequired"
  ) %in% names(maps$super_map)))
  testthat::expect_equal(unique(maps$super_map$label_source), "preexisting_annotation_source")
  testthat::expect_equal(unique(maps$super_map$Supermodule_LabelSource), "Supermodule_FinalLabel")

  expected_labels <- c(
    SM01 = paste("SM01", intToUtf8(183), "Program alpha"),
    SM02 = paste("SM02", intToUtf8(183), "Program beta"),
    SM03 = paste("SM03", intToUtf8(183), "Program gamma")
  )
  testthat::expect_equal(
    unname(maps$super_map$SupermoduleLabel),
    unname(expected_labels[maps$super_map$SupermoduleID])
  )
  testthat::expect_equal(maps$super_map$Supermodule_LabelConfidence, fixture$super_ann$Supermodule_LabelConfidence)
  testthat::expect_equal(maps$super_map$Supermodule_NameSource, fixture$super_ann$Supermodule_NameSource)
  testthat::expect_equal(maps$super_map$ManualReviewRequired, fixture$super_ann$ManualReviewRequired)
  testthat::expect_equal(maps$super_map$SupermoduleCutHeight, fixture$super_ann$SupermoduleCutHeight)

  expected_membership <- fixture$super_ann[, c("ModuleID", "module_eigengene", "SupermoduleID")]
  observed_membership <- maps$super_map[, c("ModuleID", "module_eigengene", "SupermoduleID")]
  testthat::expect_identical(observed_membership, expected_membership)
})

testthat::test_that("endpoint label collapse rejects conflicting supermodule candidates", {
  env <- endpoint_map_test_environment()
  fixture <- endpoint_map_fixture()
  fixture$super_ann$Supermodule_FinalLabel[[2]] <- "Conflicting program"
  testthat::expect_error(
    env$make_endpoint_maps(fixture$module_eig, fixture$definitions, fixture$super_ann),
    "conflicting nonmissing Supermodule_FinalLabel"
  )
})

testthat::test_that("endpoint label join rejects duplicated resolved lookup keys", {
  env <- endpoint_map_test_environment()
  fixture <- endpoint_map_fixture()
  maps <- env$make_endpoint_maps(fixture$module_eig, fixture$definitions, fixture$super_ann)
  duplicated_lookup <- dplyr::bind_rows(maps$super_label_lookup, maps$super_label_lookup[1, ])
  testthat::expect_error(
    env$join_supermodule_endpoint_labels(
      maps$super_map,
      duplicated_lookup,
      expected_modules = setdiff(names(fixture$module_eig), "Sample")
    ),
    "duplicated dataset + SupermoduleID keys",
    fixed = TRUE
  )
})
