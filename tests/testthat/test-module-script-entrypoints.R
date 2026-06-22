testthat::test_that("canonical module entrypoints exist", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  active <- c(
    "06_modules_WGCNA/01_WGCNA.r",
    "06_modules_WGCNA/02_curated_overlap_programs.r",
    "06_modules_WGCNA/03_score_module_activity.R",
    "06_modules_WGCNA/04_wgcna_de_gsea_overlap.r"
  )
  testthat::expect_true(all(file.exists(repo_path(active))))
})

testthat::test_that("module score implementation lives in 03_score_module_activity", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  txt <- paste(readLines(repo_path("06_modules_WGCNA/03_score_module_activity.R"), warn = FALSE), collapse = "\n")
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
  source(testthat::test_path("..", "..", "R", "paths.R"))
  removed <- c(
    "06_modules_WGCNA/legacy/05_module_score.r",
    "06_modules_WGCNA/legacy/91_module_score.r",
    "06_modules_WGCNA/legacy/03_overlap_modules.r",
    "06_modules_WGCNA/legacy/04_overlap_modules.r",
    "06_modules_WGCNA/legacy/05_wgcna_de_gsea_overlap.r"
  )
  testthat::expect_false(any(file.exists(repo_path(removed))))
})

testthat::test_that("pipeline module stages use canonical scripts and contracts", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  testthat::skip_if_not_installed("yaml")
  registry <- yaml::read_yaml(repo_path("pipeline.yml"))
  modules_wgcna <- registry$stages$modules_wgcna$scripts
  modules_downstream <- registry$stages$modules_downstream$scripts
  scripts_wgcna <- vapply(modules_wgcna, function(x) x$script, character(1))
  scripts_downstream <- vapply(modules_downstream, function(x) x$script, character(1))
  testthat::expect_equal(scripts_wgcna, "06_modules_WGCNA/01_WGCNA.r")
  expected_downstream <- c(
    "06_modules_WGCNA/01a_compare_GO_recurrent_proteins.r",
    "06_modules_WGCNA/02_curated_overlap_programs.r",
    "06_modules_WGCNA/03_score_module_activity.R",
    "06_modules_WGCNA/03_score_module_activity.R",
    "06_modules_WGCNA/04_wgcna_de_gsea_overlap.r",
    "06_modules_WGCNA/05_module_supermodule_group_effects.r",
    "06_modules_WGCNA/06_annotate_module_microenvironment.r",
    "06_modules_WGCNA/07_wgcna_interpretable_summary.r",
    "06_modules_WGCNA/08_wgcna_score_publication_summary.R",
    "06_modules_WGCNA/09_microglia_neuropil_independence.R",
    "06_modules_WGCNA/10_module_complex_architecture.r",
    "06_modules_WGCNA/11_module_robustness_sensitivity.r"
  )
  testthat::expect_equal(scripts_downstream, expected_downstream)
  testthat::expect_false(any(c(
    "06_modules_WGCNA/03_overlap_modules.r",
    "06_modules_WGCNA/04_overlap_modules.r",
    "06_modules_WGCNA/05_module_score.r",
    "06_modules_WGCNA/05_wgcna_de_gsea_overlap.r",
    "06_modules_WGCNA/91_module_score.r"
  ) %in% c(scripts_wgcna, scripts_downstream)))
  pipeline_txt <- paste(readLines(repo_path("pipeline.yml"), warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("module_score/<dataset>/", pipeline_txt, fixed = TRUE))
  active_txt <- paste(vapply(c(modules_wgcna, modules_downstream), function(x) paste(unlist(x), collapse = "\n"), character(1)), collapse = "\n")
  testthat::expect_false(grepl("module_score_v0.0.2", active_txt, fixed = TRUE))
})

testthat::test_that("pipeline legacy block does not retain removed wrapper names", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  pipeline_txt <- paste(readLines(repo_path("pipeline.yml"), warn = FALSE), collapse = "\n")
  for (pair in c("03_overlap_modules.r", "04_overlap_modules.r", "05_module_score.r", "05_wgcna_de_gsea_overlap.r", "91_module_score.r")) {
    testthat::expect_false(grepl(pair, pipeline_txt, fixed = TRUE), info = pair)
  }
  for (target in c(
    "02_curated_overlap_programs.r",
    "03_score_module_activity.R",
    "04_wgcna_de_gsea_overlap.r"
  )) {
    testthat::expect_true(grepl(target, pipeline_txt, fixed = TRUE), info = target)
  }
})

testthat::test_that("module score dry-run reports dataset-aware source defaults", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  run <- function(dataset) {
    cmd <- file.path(R.home("bin"), "Rscript")
    old_wd <- setwd(repo_path())
    on.exit(setwd(old_wd), add = TRUE)
    out <- suppressWarnings(system2(
      cmd,
      c("06_modules_WGCNA/03_score_module_activity.R", "--dataset", dataset, "--dry-run"),
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
