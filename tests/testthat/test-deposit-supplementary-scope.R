source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "export_helpers.R"))

# The deposit's supplementary tables leave out the EWCE cell-type tables and the
# WGCNA module evidence rankings (Tobias, 2026-11-22). The paper uses neither;
# the module table comes back only after the label review.

supplementary_globs <- function() {
  cfg <- repo_path("analysis", "publication_source_data", "config", "export_config.yml")
  unlist(yaml::read_yaml(cfg)$supplementary_table_globs, use.names = FALSE)
}

testthat::test_that("no supplementary glob names EWCE or the module evidence rank", {
  globs <- supplementary_globs()
  testthat::expect_gt(length(globs), 0L)
  bad <- globs[grepl("ewce|module_evidence_rank", globs, ignore.case = TRUE)]
  testthat::expect_identical(bad, character(0), info = paste(bad, collapse = "\n"))
})

testthat::test_that("the remaining globs do not pick up EWCE or module tables either", {
  # a broad glob can reach these files without naming them, so the globs are
  # run against a fixture tree that holds both kinds plus one table they should find
  root <- file.path(tempdir(), paste0("deposit_scope_", as.integer(stats::runif(1, 1, 1e9))))
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  excluded <- c(
    "results/enrichment/run_ewce_celltype_enrichment/microglia/tables/Supplementary_Table_EWCE.xlsx",
    "results/enrichment/run_ewce_celltype_enrichment/microglia/tables/High_Confidence_EWCE_Findings.xlsx",
    "results/tables/05_celltype_enrichment_EWCE/EWCE_E9/microglia/Supplementary_Table_EWCE.xlsx",
    "results/wgcna/build_wgcna_modules/microglia/tables/modules/WGCNA_module_evidence_rank.csv",
    "results/tables/06_modules_WGCNA/01_WGCNA/microglia/modules/WGCNA_module_evidence_rank.csv"
  )
  kept <- "results/differential_abundance/build_sus_res_dap_atlas/global/tables/source_data/sus_res_dap_membership.csv"
  for (f in file.path(root, c(excluded, kept))) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    file.create(f)
  }

  config <- list(supplementary_table_globs = file.path(root, supplementary_globs()))
  selected <- normalize_export_path(supplementary_candidate_files(config, datasets = "microglia"))

  testthat::expect_true(normalize_export_path(file.path(root, kept)) %in% selected)
  testthat::expect_false(any(normalize_export_path(file.path(root, excluded)) %in% selected))
})
