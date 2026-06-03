testthat::test_that("WGCNA downstream entrypoints exist", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  scripts <- c(
    "03_qc_exploration/06_wgcna_marker_trait_export.r",
    "06_modules_WGCNA/05_module_supermodule_group_effects.r",
    "06_modules_WGCNA/06_annotate_module_microenvironment.r",
    "06_modules_WGCNA/07_wgcna_interpretable_summary.r"
  )
  testthat::expect_true(all(file.exists(repo_path(scripts))))
})

testthat::test_that("WGCNA downstream dry-runs report contracts", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x
  cmd <- file.path(R.home("bin"), "Rscript")
  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  cases <- list(
    c("03_qc_exploration/06_wgcna_marker_trait_export.r", "--dataset", "microglia", "--dry-run"),
    c("06_modules_WGCNA/05_module_supermodule_group_effects.r", "--dataset", "microglia", "--dry-run"),
    c("06_modules_WGCNA/06_annotate_module_microenvironment.r", "--dataset", "microglia", "--dry-run"),
    c("06_modules_WGCNA/07_wgcna_interpretable_summary.r", "--dataset", "all", "--dry-run")
  )
  for (args in cases) {
    out <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
    testthat::expect_equal(attr(out, "status") %||% 0L, 0L, info = paste(args, collapse = " "))
    testthat::expect_true(any(grepl("\\[DRY-RUN", out)), info = paste(args, collapse = " "))
  }
})

testthat::test_that("WGCNA downstream schemas expose required columns", {
  source(testthat::test_path("..", "..", "R", "module_contracts.R"))
  group_cols <- c(
    "dataset", "level", "module_id", "supermodule_id", "module_label",
    "supermodule_label", "spatial_unit", "contrast", "estimate", "SE",
    "statistic", "p_value", "FDR_within_dataset_level", "FDR_global",
    "direction", "n_samples", "formula_used", "dropped_covariates",
    "rank_deficient_model", "model_warning"
  )
  group_df <- as.data.frame(setNames(rep(list(logical()), length(group_cols)), group_cols))
  testthat::expect_silent(validate_wgcna_group_effects(group_df))

  annot_df <- data.frame(
    dataset = character(), ModuleID = character(), ModuleColor = character(),
    n_proteins = integer(), microenvironment_class = character(),
    interpretation_note = character()
  )
  testthat::expect_silent(validate_wgcna_module_annotation(annot_df))

  interp_df <- data.frame(
    dataset = character(), level = character(), contrast = character(),
    estimate = numeric(), p_value = numeric(), FDR_global = numeric(),
    interpretation_sentence = character()
  )
  testthat::expect_silent(validate_wgcna_interpretable_summary(interp_df))
})
