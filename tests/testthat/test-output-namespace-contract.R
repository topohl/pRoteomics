source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::test_that("output namespace contract is complete and classifies roles", {
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
    "must be one of 01, 02, 03"
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
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "pipeline_registry.R"))
  source(repo_path("R", "output_namespace_utils.R"))

  registry <- read_pipeline_registry(repo_path("pipeline.yml"))
  integration <- registry$stages$integration$scripts
  figure_steps <- Filter(
    function(step) startsWith(as.character(step$script), "figures/"),
    integration
  )
  # Phase 6C moved all nine figure entry points to Exp9_manuscript. The
  # inventory guard becomes its inverse, which is stronger: no figure entry
  # point may be registered here at all, so an accidental one is still
  # caught. The namespace classification below is retained so that if one
  # ever reappears it must still declare manuscript_authoring outputs only.
  testthat::expect_length(figure_steps, 0L)
  for (step in figure_steps) {
    outputs <- as.character(unlist(step$produces, use.names = FALSE))
    testthat::expect_true(all(
      classify_output_namespace(outputs)$namespace == "manuscript_authoring"
    ))
  }
})

testthat::test_that("figure 1 is a first-class manuscript figure namespace", {
  source(repo_path("R", "output_namespace_utils.R"))

  root <- tempfile("fig1-root-")
  # Figure 1 was declared as an export destination in pipeline.yml while the
  # validator rejected every ID but 02 and 03, so the slot was unreachable. All
  # three must now resolve, and the authoring convention stays zero-padded.
  for (id in c("01", "02", "03")) {
    p <- output_namespace_manuscript_figure_paths(root, id)
    testthat::expect_identical(
      p$panels,
      file.path(root, "figures", "manuscript", paste0("figure_", id), "panels")
    )
    testthat::expect_identical(
      p$source_data,
      file.path(root, "source_data", "manuscript", paste0("figure_", id))
    )
    testthat::expect_identical(
      classify_output_namespace(
        file.path("results", "figures", "manuscript", paste0("figure_", id))
      )$namespace,
      "manuscript_authoring"
    )
  }
  # unpadded input still normalises to the padded stub
  testthat::expect_identical(
    output_namespace_manuscript_figure_paths(root, 1)$panels,
    output_namespace_manuscript_figure_paths(root, "01")$panels
  )

  # Still fails closed. Widening the allow-list must not become a wildcard: an
  # unrecognised figure number silently creating figure_07 is the failure this
  # validation exists to prevent.
  for (bad in list("04", "00", "10", 7L, "", "abc", NA_character_,
                   "../../etc", "02; rm -rf", c("02", "03"))) {
    testthat::expect_error(
      output_namespace_manuscript_figure_paths(root, bad),
      "Manuscript figure ID must be one of"
    )
  }
})

testthat::test_that("the export router recognises the same figure set", {
  source(repo_path("R", "output_namespace_utils.R"))
  source(repo_path("R", "export_helpers.R"))

  # The router and the validator must not disagree about which figures exist.
  testthat::expect_identical(MANUSCRIPT_FIGURE_IDS, c("01", "02", "03"))

  mroot <- file.path(tempfile("ms-"), "manuscript")
  rel <- c("manuscript/figure_01/panels/a.svg",
           "manuscript/figure_02/panels/b.svg",
           "manuscript/figure_03/panels/c.svg",
           "manuscript/figure_07/panels/d.svg")
  out <- manuscript_curated_figure_target_paths(rel, mroot)
  # authoring stubs are zero-padded, export destinations are not
  testthat::expect_identical(out[[1]], file.path(mroot, "figure_1", "panels", "a.svg"))
  testthat::expect_identical(out[[2]], file.path(mroot, "figure_2", "panels", "b.svg"))
  testthat::expect_identical(out[[3]], file.path(mroot, "figure_3", "panels", "c.svg"))
  # an unrecognised figure is not given a destination of its own
  testthat::expect_false(grepl("figure_7", out[[4]], fixed = TRUE))
  testthat::expect_true(grepl("extended_data", out[[4]], fixed = TRUE))
})
