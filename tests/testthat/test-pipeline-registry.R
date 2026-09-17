source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::test_that("pipeline.yml is valid and references existing active scripts", {
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "pipeline_registry.R"))

  testthat::skip_if_not_installed("yaml")
  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
  testthat::expect_silent(validate_pipeline_scripts_exist(registry))
  testthat::expect_silent(validate_run_order_against_registry(registry))
  testthat::expect_identical(
    run_order_registry_index_references(repo_path("RUN_ORDER.md")),
    pipeline_registry_entries(registry)$script
  )
  steps <- pipeline_steps(registry, "enrichment", dataset = "microglia")
  testthat::expect_true("analysis/04_differential_abundance/01_clusterProfiler.r" %in% steps$script)

  network_steps <- pipeline_steps(registry, "networks", dataset = "microglia")
  testthat::expect_true("analysis/07_spatial_networks/01_network_spatial_relations.r" %in% network_steps$script)
  downstream_networks <- setdiff(network_steps$script, "analysis/07_spatial_networks/01_network_spatial_relations.r")
  testthat::expect_false(any(grepl("02_differential|03_bootstrap|04_bootstrap|05_bootstrap|06_chord", downstream_networks)))

  coupling_steps <- pipeline_steps(registry, "coupling", dataset = "microglia")
  testthat::expect_true("analysis/08_integration/03_module_behavior_coupling.r" %in% coupling_steps$script[coupling_steps$supported])
  testthat::expect_false("analysis/08_integration/02_network_behavior_coupling.r" %in% coupling_steps$script[coupling_steps$supported])

  testthat::expect_true("integration" %in% pipeline_stage_names(registry))
  stage_names <- pipeline_stage_names(registry)
  testthat::expect_lt(match("coupling", stage_names), match("integration", stage_names))
  testthat::expect_lt(match("integration", stage_names), match("export", stage_names))
  integration_steps <- pipeline_steps(
    registry, "integration", dataset = "all", include_unsupported = TRUE
  )$script
  # The Figure 2 and Figure 3 entry points moved to Exp9_manuscript in
  # Phase 6C. What this repository must now guarantee is the inverse:
  # no figure renderer is registered as a pipeline step at all.
  testthat::expect_false(any(grepl("^figures/", integration_steps)))
})

testthat::test_that("deprecated 04d stays excluded and documented as legacy", {
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "pipeline_registry.R"))

  testthat::skip_if_not_installed("yaml")
  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
  active <- pipeline_registry_entries(registry)$script
  legacy <- vapply(registry$legacy, function(x) as.character(x$script), character(1))
  script_04d <- "archive/02_qc/04d_compartment_marker_fidelity.r"
  testthat::expect_false(script_04d %in% active)
  testthat::expect_true(script_04d %in% legacy)

  audit <- utils::read.delim(
    repo_path("docs", "active_script_io_audit.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  row_04d <- audit[audit$script_path == script_04d, , drop = FALSE]
  testthat::expect_equal(nrow(row_04d), 1L)
  testthat::expect_identical(row_04d$status, "legacy_deprecated")
  testthat::expect_match(row_04d$remaining_TODOs, "04e", fixed = TRUE)
})

testthat::test_that("analysis discovery is repository-wide and fail closed", {
  source(repo_path("R", "pipeline_registry.R"))

  candidates <- c(
    "analysis/11_new_analysis/01_new_result.R",
    "tests/testthat/test-new-result.R",
    "archive/deprecated/01_old_result.r",
    "audits/part29/01_new_audit.R",
    "R/utilities/new_helper.R",
    "tools/new_audit.R",
    "analysis/12_future/legacy/01_archived.R",
    "run_dataset_pipeline.R"
  )
  observed <- active_analysis_scripts(candidate_files = candidates)
  testthat::expect_identical(
    observed, "analysis/11_new_analysis/01_new_result.R"
  )
})

testthat::test_that("current optional and superseded blind-spot scripts are classified", {
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "pipeline_registry.R"))

  testthat::skip_if_not_installed("yaml")
  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
  active <- pipeline_registry_entries(registry)$script
  legacy <- vapply(
    registry$legacy, function(x) as.character(x$script), character(1)
  )
  testthat::expect_true(
    "analysis/08_integration/05_manuscript_figure3_wgcna_protein_zoom.R" %in%
      active
  )
  testthat::expect_false(
    "archive/08_integration/01_compartment_fidelity_summary.R" %in%
      active
  )
  testthat::expect_true(
    "archive/08_integration/01_compartment_fidelity_summary.R" %in%
      legacy
  )
  audit <- write_pipeline_validation_tables(registry)
  testthat::expect_equal(nrow(audit$unregistered), 0L)
})

testthat::test_that("documented dry-run stages are real registry stages", {
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "pipeline_registry.R"))

  testthat::skip_if_not_installed("yaml")
  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
  stage_names <- pipeline_stage_names(registry)
  documented <- c("qc", "enrichment", "modules_wgcna", "modules_downstream", "integration")
  testthat::expect_setequal(intersect(documented, stage_names), documented)

  readme <- paste(readLines(repo_path("README.md"), warn = FALSE), collapse = "\n")
  run_order <- paste(readLines(repo_path("RUN_ORDER.md"), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("--stage modules\\b", paste(readme, run_order), perl = TRUE))
  testthat::expect_false(grepl("--stage behavior\\b", paste(readme, run_order), perl = TRUE))
})

testthat::test_that("README and RUN_ORDER do not present legacy scripts as active", {

  readme <- paste(readLines(repo_path("README.md"), warn = FALSE), collapse = "\n")
  run_order <- paste(readLines(repo_path("RUN_ORDER.md"), warn = FALSE), collapse = "\n")

  active_blocks <- paste(readme, run_order, sep = "\n")
  testthat::expect_false(grepl("Backward-compatible retained names", active_blocks, fixed = TRUE))
  testthat::expect_false(grepl("04_neuropil_contamination_annotation.r", run_order, fixed = TRUE))
  testthat::expect_true(grepl("analysis/05_wgcna/03_score_module_activity.R", run_order, fixed = TRUE))
  testthat::expect_false(grepl("analysis/05_wgcna/91_module_score.r", run_order, fixed = TRUE))
  testthat::expect_true(grepl("analysis/09_publication_exports/", readme, fixed = TRUE))
  testthat::expect_true(grepl("legacy", readme, ignore.case = TRUE))
})

testthat::test_that("bespoke enrichment legacy scripts stay out of the active folder", {
  source(repo_path("R", "pipeline_registry.R"))

  testthat::skip_if_not_installed("yaml")
  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
  steps <- pipeline_steps(registry, pipeline_stage_names(registry), dataset = "all", include_unsupported = TRUE)
  moved <- c(
    "analysis/04_differential_abundance/04_compare_pathways.r",
    "analysis/04_differential_abundance/05_compare_sig_expr.r",
    "analysis/04_differential_abundance/07_control_strata_enrichment_figures.r"
  )
  legacy <- c(
    "archive/04_differential_expression_enrichment/legacy/04_compare_pathways.r",
    "archive/04_differential_expression_enrichment/legacy/05_compare_sig_expr.r",
    "archive/04_differential_expression_enrichment/legacy/07_control_strata_enrichment_figures.r"
  )
  testthat::expect_false(any(file.exists(repo_path(moved))))
  testthat::expect_true(all(file.exists(repo_path(legacy))))
  testthat::expect_false(any(legacy %in% steps$script))
})
