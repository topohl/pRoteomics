source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::test_that("canonical module entrypoints exist", {
  active <- c(
    "analysis/05_wgcna/build_wgcna_modules.R",
    "analysis/05_wgcna/render_module_go_heatmaps.R",
    "analysis/05_wgcna/build_curated_overlap_programs.R",
    "analysis/05_wgcna/score_module_activity.R",
    "analysis/05_wgcna/compare_module_enrichment_overlap.R"
  )
  testthat::expect_true(all(file.exists(repo_path(active))))
})

testthat::test_that("module score implementation lives in 03_score_module_activity", {
  txt <- paste(readLines(repo_path("analysis/05_wgcna/score_module_activity.R"), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("source(repo_path(\"06_modules_WGCNA\", \"05_module_score.r\"))", txt, fixed = TRUE))
  testthat::expect_false(grepl("source(repo_path(\"06_modules_WGCNA\", \"91_module_score.r\"))", txt, fixed = TRUE))
  for (needle in c(
    "default_module_definition_source",
    "module_definition_source",
    "validate_module_score_output",
    "module_feature_mapping_trace.csv",
    "module_gene_coverage.csv",
    "write_run_manifest",
    "module_score"
  )) {
    testthat::expect_true(grepl(needle, txt, fixed = TRUE), info = needle)
  }
})

testthat::test_that("legacy module wrapper scripts have been removed", {
  # Named by basename rather than by a stage path. Checking a pre-migration
  # path would pass vacuously after Phase 6B, because nothing can exist
  # there; checking the basename across the whole active namespace keeps the
  # assertion biting under any directory layout.
  removed <- c(
    "05_module_score.r",
    "91_module_score.r",
    "03_overlap_modules.r",
    "04_overlap_modules.r",
    "05_wgcna_de_gsea_overlap.r"
  )
  active <- c(list.files(repo_path("analysis"), recursive = TRUE),
              list.files(repo_path("R"), recursive = TRUE))
  testthat::expect_false(any(basename(active) %in% removed))
})

testthat::test_that("pipeline module stages use canonical scripts and contracts", {
  testthat::skip_if_not_installed("yaml")
  registry <- yaml::read_yaml(repo_path("pipeline.yml"))
  modules_wgcna <- registry$stages$modules_wgcna$scripts
  modules_downstream <- registry$stages$modules_downstream$scripts
  scripts_wgcna <- vapply(modules_wgcna, function(x) x$script, character(1))
  scripts_downstream <- vapply(modules_downstream, function(x) x$script, character(1))
  testthat::expect_equal(scripts_wgcna, "analysis/05_wgcna/build_wgcna_modules.R")
  required_downstream <- c(
    "analysis/05_wgcna/render_module_go_heatmaps.R",
    "analysis/05_wgcna/compare_recurrent_module_proteins.R",
    "analysis/05_wgcna/build_curated_overlap_programs.R",
    "analysis/05_wgcna/score_module_activity.R",
    "analysis/05_wgcna/compare_module_enrichment_overlap.R",
    "analysis/05_wgcna/build_module_identity_contract.R",
    "analysis/05_wgcna/test_module_phenotypes.R",
    "analysis/05_wgcna/annotate_module_microenvironment.R",
    "analysis/05_wgcna/summarize_module_interpretation.R",
    "analysis/05_wgcna/render_module_figures.R",
    "analysis/05_wgcna/summarize_module_scores.R",
    "analysis/05_wgcna/render_microglia_module_figures.R",
    "analysis/05_wgcna/test_microglia_neuropil_independence.R",
    "analysis/05_wgcna/render_microglia_independence_figures.R",
    "analysis/05_wgcna/summarize_microglia_roi_specificity.R",
    "analysis/05_wgcna/summarize_module_complex_architecture.R",
    "analysis/05_wgcna/audit_module_robustness.R",
    "analysis/05_wgcna/audit_microglia_module_claims.R",
    "analysis/05_wgcna/summarize_microglia_module_claims.R",
    "analysis/05_wgcna/audit_module_claim_readiness.R"
  )
  testthat::expect_true(all(required_downstream %in% scripts_downstream))
  testthat::expect_equal(
    sum(scripts_downstream == "analysis/05_wgcna/score_module_activity.R"),
    2L
  )
  ordered_contract <- c(
    "analysis/05_wgcna/build_module_identity_contract.R",
    "analysis/05_wgcna/test_module_phenotypes.R",
    "analysis/05_wgcna/annotate_module_microenvironment.R",
    "analysis/05_wgcna/summarize_module_interpretation.R",
    "analysis/05_wgcna/audit_module_claim_readiness.R"
  )
  testthat::expect_true(all(diff(match(ordered_contract, scripts_downstream)) > 0L))
  testthat::expect_false(any(c(
    "analysis/05_wgcna/03_overlap_modules.r",
    "analysis/05_wgcna/04_overlap_modules.r",
    "analysis/05_wgcna/05_module_score.r",
    "analysis/05_wgcna/05_wgcna_de_gsea_overlap.r",
    "analysis/05_wgcna/91_module_score.r"
  ) %in% c(scripts_wgcna, scripts_downstream)))
  pipeline_txt <- paste(readLines(repo_path("pipeline.yml"), warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("module_score/<dataset>/", pipeline_txt, fixed = TRUE))
  active_txt <- paste(vapply(c(modules_wgcna, modules_downstream), function(x) paste(unlist(x), collapse = "\n"), character(1)), collapse = "\n")
  testthat::expect_false(grepl("module_score_v0.0.2", active_txt, fixed = TRUE))
})

testthat::test_that("pipeline legacy block does not retain removed wrapper names", {
  pipeline_txt <- paste(readLines(repo_path("pipeline.yml"), warn = FALSE), collapse = "\n")
  for (pair in c("03_overlap_modules.r", "04_overlap_modules.r", "05_module_score.r", "05_wgcna_de_gsea_overlap.r", "91_module_score.r")) {
    testthat::expect_false(grepl(pair, pipeline_txt, fixed = TRUE), info = pair)
  }
  for (target in c(
    "build_curated_overlap_programs.R",
    "score_module_activity.R",
    "compare_module_enrichment_overlap.R"
  )) {
    testthat::expect_true(grepl(target, pipeline_txt, fixed = TRUE), info = target)
  }
})

testthat::test_that("module score dry-run reports dataset-aware source defaults", {
  run <- function(dataset) {
    cmd <- file.path(R.home("bin"), "Rscript")
    old_wd <- setwd(repo_path())
    on.exit(setwd(old_wd), add = TRUE)
    out <- suppressWarnings(system2(
      cmd,
      c("analysis/05_wgcna/score_module_activity.R", "--dataset", dataset, "--dry-run"),
      stdout = TRUE,
      stderr = TRUE
    ))
    paste(out, collapse = "\n")
  }
  microglia <- run("microglia")
  testthat::expect_true(grepl("Module definition source: wgcna", microglia, fixed = TRUE))
  testthat::expect_true(grepl("Will score: WGCNA modules", microglia, fixed = TRUE))
  testthat::expect_true(grepl("Expected dataset/source-scoped table root", microglia, fixed = TRUE))

  soma <- run("neuron_soma")
  testthat::expect_true(grepl("Module definition source: wgcna", soma, fixed = TRUE))

  neuropil <- run("neuron_neuropil")
  testthat::expect_true(grepl("Module definition source: overlap", neuropil, fixed = TRUE))
  testthat::expect_true(grepl("Will score: curated overlap programs", neuropil, fixed = TRUE))
})
