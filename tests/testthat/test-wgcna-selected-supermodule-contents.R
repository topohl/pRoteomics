testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")
testthat::skip_if_not_installed("stringr")

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))

selected_contents_test_environment <- function() {
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

selected_contents_fixture <- function() {
  list(
    selected = data.frame(
      dataset = "input_dataset_must_not_win",
      supermodule_id = "SM01",
      supermodule_label = NA_character_,
      contrast = "SUS - RES",
      estimate = 0.42,
      p_value = 0.01,
      FDR_within_dataset_level = 0.04,
      FDR_global = 0.06,
      evidence_status = "supported",
      selection_support = "FDR <= 0.05",
      selection_message = "selected",
      stringsAsFactors = FALSE
    ),
    pca_eigenvalues = data.frame(
      dataset = "input_dataset_must_not_win",
      supermodule_id = "SM01",
      pc = 1L,
      variance_explained = 0.72,
      stringsAsFactors = FALSE
    ),
    comp = data.frame(
      dataset = "input_dataset_must_not_win",
      supermodule_id = "SM01",
      supermodule_eigengene = "SM__SM01",
      n_member_modules = 2L,
      member_modules = "MEblue;MEbrown",
      stringsAsFactors = FALSE
    ),
    definitions = data.frame(
      module_eigengene = c("MEblue", "MEbrown"),
      ModuleID = c("WGCNA_m01", "WGCNA_m02"),
      ModuleLabel_Final = c("Blue module", "Brown module"),
      stringsAsFactors = FALSE
    ),
    super_ann = data.frame(
      dataset = "input_dataset_must_not_win",
      SupermoduleID = rep("SM01", 2),
      module_eigengene = c("MEblue", "MEbrown"),
      ModuleID = c("WGCNA_m01", "WGCNA_m02"),
      ModuleLabel_Final = c("Blue module", "Brown module"),
      stringsAsFactors = FALSE
    ),
    super_summary = data.frame(
      dataset = "input_dataset_must_not_win",
      SupermoduleID = "SM01",
      Supermodule_DisplayLabel = paste("SM01", intToUtf8(183), "Summary program"),
      member_modules = "WGCNA_m01;WGCNA_m02",
      n_modules = 2L,
      top_GO_BP_terms = "vesicle organization",
      top_GO_MF_terms = "RNA binding",
      top_GO_CC_terms = "cytoplasmic vesicle membrane",
      top_hub_symbols = "GeneA;GeneB",
      top_hub_proteins = "P1;P2",
      stringsAsFactors = FALSE
    ),
    modules_long = data.frame(
      ModuleID = c("WGCNA_m01", "WGCNA_m02"),
      abs_kME = c(0.91, 0.84),
      is_core_kME_0.6 = TRUE,
      is_top_hub_25 = c(TRUE, FALSE),
      GeneSymbol = c("GeneA", "GeneB"),
      stringsAsFactors = FALSE
    ),
    module_summary = data.frame(
      ModuleID = c("WGCNA_m01", "WGCNA_m02"),
      top_hub_proteins = c("P1", "P2"),
      stringsAsFactors = FALSE
    ),
    go_enrichment = data.frame(
      ModuleID = c("WGCNA_m01", "WGCNA_m02"),
      ModuleProteinSetType = "all",
      Ontology = "BP",
      Description = c("fallback BP one", "fallback BP two"),
      p.adjust = c(0.01, 0.02),
      stringsAsFactors = FALSE
    ),
    module_name_map = data.frame(
      ModuleID = c("WGCNA_m01", "WGCNA_m02"),
      ModuleLabel_Final = c("Blue module", "Brown module"),
      stringsAsFactors = FALSE
    )
  )
}

call_selected_contents <- function(env, fixture, dataset = "microglia") {
  env$selected_sus_res_supermodule_contents(
    selected = fixture$selected,
    pca_eigenvalues = fixture$pca_eigenvalues,
    comp = fixture$comp,
    definitions = fixture$definitions,
    super_ann = fixture$super_ann,
    super_summary = fixture$super_summary,
    modules_long = fixture$modules_long,
    module_summary = fixture$module_summary,
    go_enrichment = fixture$go_enrichment,
    module_name_map = fixture$module_name_map,
    dataset = dataset
  )
}

testthat::test_that("selected supermodule contents use the explicit dataset and retain evidence", {
  env <- selected_contents_test_environment()
  fixture <- selected_contents_fixture()
  had_global_ds <- exists("ds", envir = .GlobalEnv, inherits = FALSE)
  old_global_ds <- if (had_global_ds) get("ds", envir = .GlobalEnv, inherits = FALSE) else NULL
  if (had_global_ds) rm("ds", envir = .GlobalEnv)
  on.exit({
    if (had_global_ds) assign("ds", old_global_ds, envir = .GlobalEnv) else if (exists("ds", envir = .GlobalEnv, inherits = FALSE)) rm("ds", envir = .GlobalEnv)
  }, add = TRUE)

  testthat::expect_false(exists("ds", envir = .GlobalEnv, inherits = FALSE))
  out <- call_selected_contents(env, fixture, dataset = "microglia")

  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$dataset, "microglia")
  testthat::expect_equal(out$supermodule_id, "SM01")
  testthat::expect_equal(out$supermodule_label, paste("SM01", intToUtf8(183), "Summary program"))
  testthat::expect_equal(out$pca_PC1_variance_explained, 0.72)
  testthat::expect_equal(out$n_member_modules, 2L)
  testthat::expect_equal(out$member_modules, "WGCNA_m01;WGCNA_m02")
  testthat::expect_match(out$member_module_labels, "Blue module")
  testthat::expect_match(out$member_module_labels, "Brown module")
  testthat::expect_equal(out$top_GO_BP_terms, "vesicle organization")
  testthat::expect_equal(out$top_GO_MF_terms, "RNA binding")
  testthat::expect_equal(out$top_GO_CC_terms, "cytoplasmic vesicle membrane")
})

testthat::test_that("a misleading global ds cannot override the dataset argument", {
  env <- selected_contents_test_environment()
  fixture <- selected_contents_fixture()
  had_global_ds <- exists("ds", envir = .GlobalEnv, inherits = FALSE)
  old_global_ds <- if (had_global_ds) get("ds", envir = .GlobalEnv, inherits = FALSE) else NULL
  assign("ds", "wrong_dataset", envir = .GlobalEnv)
  on.exit({
    if (had_global_ds) assign("ds", old_global_ds, envir = .GlobalEnv) else rm("ds", envir = .GlobalEnv)
  }, add = TRUE)

  out <- call_selected_contents(env, fixture, dataset = "microglia")
  testthat::expect_equal(out$dataset, "microglia")
})

testthat::test_that("empty selected contents return the typed contract", {
  env <- selected_contents_test_environment()
  empty <- env$selected_sus_res_supermodule_contents(
    selected = data.frame(), pca_eigenvalues = data.frame(), comp = data.frame(),
    definitions = data.frame(), super_ann = data.frame(), super_summary = data.frame(),
    modules_long = data.frame(), module_summary = data.frame(), go_enrichment = data.frame(),
    module_name_map = data.frame(), dataset = "microglia"
  )
  testthat::expect_equal(nrow(empty), 0L)
  testthat::expect_type(empty$dataset, "character")
  testthat::expect_type(empty$pca_PC1_variance_explained, "double")
  testthat::expect_type(empty$n_member_modules, "integer")
})

testthat::test_that("singular mixed-model rows remain diagnostic and non-claimable", {
  env <- selected_contents_test_environment()
  row <- data.frame(
    model_warning = "boundary (singular) fit: see help('isSingular')",
    model_type = "lmerTest_lmer",
    formula_requested = "eigengene ~ StressGroup + (1 | AnimalID)",
    formula_used = "eigengene ~ StressGroup + (1 | AnimalID)",
    rank_deficient_model = FALSE,
    p_value = 0.01,
    stringsAsFactors = FALSE
  )
  audited <- env$add_claim_model_fields(row)
  testthat::expect_equal(nrow(audited), 1L)
  testthat::expect_true(audited$singular_model)
  testthat::expect_false(audited$primary_model_stable)
  testthat::expect_false(audited$claim_allowed_model)
  testthat::expect_match(audited$model_downgrade_reason, "singular_model")
})
