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

testthat::test_that("legacy module scripts are wrappers", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  cases <- list(
    list(path = "06_modules_WGCNA/05_module_score.r", target = "03_score_module_activity.R", word = "Deprecated"),
    list(path = "06_modules_WGCNA/91_module_score.r", target = "03_score_module_activity.R", word = "Legacy"),
    list(path = "06_modules_WGCNA/03_overlap_modules.r", target = "02_curated_overlap_programs.r", word = "deprecated"),
    list(path = "06_modules_WGCNA/04_overlap_modules.r", target = "02_curated_overlap_programs.r", word = "Deprecated"),
    list(path = "06_modules_WGCNA/05_wgcna_de_gsea_overlap.r", target = "04_wgcna_de_gsea_overlap.r", word = "Deprecated")
  )
  for (case in cases) {
    txt <- paste(readLines(repo_path(case$path), warn = FALSE), collapse = "\n")
    testthat::expect_true(grepl(case$target, txt, fixed = TRUE), info = case$path)
    testthat::expect_true(grepl(case$word, txt, ignore.case = TRUE), info = case$path)
  }
})

testthat::test_that("pipeline modules stage uses canonical scripts and contracts", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  testthat::skip_if_not_installed("yaml")
  registry <- yaml::read_yaml(repo_path("pipeline.yml"))
  modules <- registry$stages$modules$scripts
  scripts <- vapply(modules, function(x) x$script, character(1))
  expected <- c(
    "06_modules_WGCNA/01_WGCNA.r",
    "06_modules_WGCNA/02_curated_overlap_programs.r",
    "06_modules_WGCNA/03_score_module_activity.R",
    "06_modules_WGCNA/04_wgcna_de_gsea_overlap.r"
  )
  testthat::expect_equal(scripts, expected)
  testthat::expect_false(any(c(
    "06_modules_WGCNA/03_overlap_modules.r",
    "06_modules_WGCNA/04_overlap_modules.r",
    "06_modules_WGCNA/05_module_score.r",
    "06_modules_WGCNA/05_wgcna_de_gsea_overlap.r",
    "06_modules_WGCNA/91_module_score.r"
  ) %in% scripts))
  pipeline_txt <- paste(readLines(repo_path("pipeline.yml"), warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("module_score/<dataset>/<module_definition_source>", pipeline_txt, fixed = TRUE))
  active_txt <- paste(vapply(modules, function(x) paste(unlist(x), collapse = "\n"), character(1)), collapse = "\n")
  testthat::expect_false(grepl("module_score_v0.0.2", active_txt, fixed = TRUE))
})

testthat::test_that("pipeline legacy block maps old module names", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  pipeline_txt <- paste(readLines(repo_path("pipeline.yml"), warn = FALSE), collapse = "\n")
  for (pair in c(
    "03_overlap_modules.r",
    "04_overlap_modules.r",
    "05_module_score.r",
    "05_wgcna_de_gsea_overlap.r",
    "91_module_score.r"
  )) {
    testthat::expect_true(grepl(pair, pipeline_txt, fixed = TRUE), info = pair)
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
