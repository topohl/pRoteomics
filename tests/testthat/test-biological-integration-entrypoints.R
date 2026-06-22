testthat::test_that("biological integration entrypoints exist and dry-run", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x
  scripts <- c(
    "04_differential_expression_enrichment/08_external_stress_disease_signature_overlap.r",
    "08_behavior_physio_coupling/03_module_behavior_coupling.r",
    "10_biological_integration/01_cross_compartment_program_atlas.r",
    "10_biological_integration/02_manuscript_program_summary.r",
    "10_biological_integration/03_evidence_priority_matrix.r"
  )
  testthat::expect_true(all(file.exists(repo_path(scripts))))

  cmd <- file.path(R.home("bin"), "Rscript")
  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  cases <- list(
    c("04_differential_expression_enrichment/08_external_stress_disease_signature_overlap.r", "--dry-run"),
    c("08_behavior_physio_coupling/03_module_behavior_coupling.r", "--dataset", "microglia", "--dry-run"),
    c("10_biological_integration/01_cross_compartment_program_atlas.r", "--dry-run"),
    c("10_biological_integration/02_manuscript_program_summary.r", "--dry-run"),
    c("10_biological_integration/03_evidence_priority_matrix.r", "--dry-run")
  )
  for (args in cases) {
    out <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
    testthat::expect_equal(attr(out, "status") %||% 0L, 0L, info = paste(args, collapse = " "))
    testthat::expect_true(any(grepl("\\[DRY-RUN", out)), info = paste(args, collapse = " "))
  }
})

testthat::test_that("integration helpers standardize evidence contracts", {
  source(testthat::test_path("..", "..", "R", "integration_utils.R"))
  ev <- standardize_evidence(data.frame(dataset = "microglia", evidence_domain = "test", evidence_id = "x", program_label = "immune", stringsAsFactors = FALSE))
  testthat::expect_true(all(names(empty_evidence()) %in% names(ev)))
  testthat::expect_silent(validate_cross_compartment_program_atlas(ev))
})
