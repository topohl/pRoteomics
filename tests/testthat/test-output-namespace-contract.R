testthat::test_that("output namespace contract is complete and classifies roles", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "output_namespace_utils.R"))

  contract <- read_output_namespace_contract()
  testthat::expect_silent(validate_output_namespace_contract(contract))
  audit <- classify_output_namespace(c(
    "results/figures/03_qc_exploration/global/panel.svg",
    "results/figures/manuscript/figure_02/panels/figure_02b.svg",
    "results/figures/manuscript_panels/figure_3/legacy.svg",
    "results/manuscript/figure_2/panels/figure_02b.svg",
    "results/manuscript/_failed_20260901_maxpath/file.svg",
    "results/reviewer_audit/check.csv",
    "results/figures/04_differential_expression_enrichment_comparison/x.svg",
    "data/processed/01_preprocessing/input.rds",
    "pride_submission/metadata/sdrf_like_metadata.tsv",
    "config/marker_panels/wgcna_reference_marker_sets.csv"
  ))
  testthat::expect_identical(audit$namespace, c(
    "canonical_stage_output",
    "manuscript_authoring",
    "legacy_manuscript_authoring",
    "manuscript_export",
    "historical_or_failed_export",
    "diagnostic_or_audit",
    "comparison_or_candidate",
    "processed_data",
    "pride_export",
    "source_control_configuration"
  ))
})

testthat::test_that("manuscript figure paths keep authoring separate from export", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "output_namespace_utils.R"))

  root <- tempfile("output-root-")
  paths <- output_namespace_manuscript_figure_paths(root, "03")
  testthat::expect_identical(
    paths$panels,
    file.path(root, "figures", "manuscript", "figure_03", "panels")
  )
  testthat::expect_identical(
    paths$source_data,
    file.path(root, "source_data", "manuscript", "figure_03")
  )
  testthat::expect_error(
    output_namespace_manuscript_figure_paths(root, "04"),
    "must be 02 or 03"
  )
  testthat::expect_identical(
    classify_output_namespace(
      "results/figures/manuscript/figure_03/panels"
    )$namespace,
    "manuscript_authoring"
  )
  testthat::expect_identical(
    relative_to(output_namespace_manuscript_export_root()),
    "results/manuscript"
  )
})

testthat::test_that("pipeline manuscript entry points declare authoring outputs only", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "pipeline_registry.R"))
  source(repo_path("R", "output_namespace_utils.R"))

  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
  integration <- registry$stages$integration$scripts
  figure_steps <- Filter(
    function(step) startsWith(as.character(step$script), "figures/"),
    integration
  )
  # Inventory guard: figure_02, figure_03, and the two manuscript-supporting
  # immunostaining renderers (three-candidate comparison and ten-candidate
  # panel). Bump deliberately when a figure entry point is added, so an
  # accidental one is still caught.
  testthat::expect_length(figure_steps, 4L)
  for (step in figure_steps) {
    outputs <- as.character(unlist(step$produces, use.names = FALSE))
    testthat::expect_true(all(
      classify_output_namespace(outputs)$namespace == "manuscript_authoring"
    ))
  }
})
