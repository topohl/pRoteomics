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
    status <- attr(out, "status") %||% 0L
    testthat::expect_true(status %in% c(0L, 1L), info = paste(args, collapse = " "))
    if (identical(status, 1L)) {
      testthat::expect_true(
        any(grepl("missing_required_input", out, fixed = TRUE)),
        info = paste(args, collapse = " ")
      )
    }
    testthat::expect_true(any(grepl("\\[DRY-RUN", out)), info = paste(args, collapse = " "))
  }
})

testthat::test_that("integration helpers standardize evidence contracts", {
  source(testthat::test_path("..", "..", "R", "integration_utils.R"))
  ev <- standardize_evidence(data.frame(
    dataset = "microglia", evidence_domain = "test", evidence_id = "x",
    program_label = "immune", tier_specific_family_id = "primary_family",
    tier_specific_family_size = 13L, inferential_claim_gate = "eligible",
    inferential_source_file = "handoff.csv", inferential_source_key = "key",
    stringsAsFactors = FALSE
  ))
  testthat::expect_true(all(names(empty_evidence()) %in% names(ev)))
  testthat::expect_identical(ev$tier_specific_family_size, 13L)
  testthat::expect_identical(ev$inferential_source_key, "key")
  testthat::expect_identical(ev$evidence_source_family, "test")
  testthat::expect_silent(validate_cross_compartment_program_atlas(ev))
})

testthat::test_that("integration counts shared ranked-GSEA summaries as one evidence lineage", {
  source(testthat::test_path("..", "..", "R", "integration_utils.R"))
  evidence <- standardize_evidence(data.frame(
    dataset = c("neuron_soma", "neuron_soma"),
    evidence_domain = c("enrichment_program", "spatial_architecture"),
    evidence_id = c("program", "spatial"), program_label = "RNA/translation",
    entity_type = "program", evidence_status = "available",
    stringsAsFactors = FALSE
  ))
  testthat::expect_identical(evidence$evidence_source_family, c("ranked_GSEA", "ranked_GSEA"))
  testthat::expect_identical(independent_evidence_families(evidence), "ranked_GSEA")

  architecture <- assign_downstream_evidence_roles(data.frame(
    dataset = "microglia", evidence_domain = "wgcna_architecture",
    evidence_id = "WGCNA_m01", program_label = "OXPHOS", entity_type = "module",
    evidence_status = "descriptive_architecture_only",
    evidence_semantic_class = "wgcna_architecture",
    wgcna_group_effect_status = "not_FDR_supported",
    counts_toward_convergence = TRUE,
    stringsAsFactors = FALSE
  ))
  testthat::expect_identical(architecture$evidence_role, "validated_wgcna_architecture")
  testthat::expect_false(architecture$counts_toward_convergence)
})

testthat::test_that("integration CSV reads inspect sparse provenance columns fully", {
  helper_text <- paste(
    readLines(testthat::test_path("..", "..", "R", "integration_utils.R"), warn = FALSE),
    collapse = "\n"
  )
  testthat::expect_match(helper_text, "guess_max = Inf", fixed = TRUE)

  path <- tempfile(fileext = ".csv")
  readr::write_csv(
    data.frame(
      row_id = seq_len(1002L),
      sparse_provenance = c(rep(NA_character_, 1001L), "module")
    ),
    path, na = ""
  )
  parsed <- readr::read_csv(
    path, show_col_types = FALSE, progress = FALSE,
    guess_max = Inf
  )
  testthat::expect_identical(parsed$sparse_provenance[[1002L]], "module")
  testthat::expect_equal(nrow(readr::problems(parsed)), 0L)

  bundle_text <- paste(
    readLines(testthat::test_path("..", "..", "R", "final_evidence_bundle_utils.R"), warn = FALSE),
    collapse = "\n"
  )
  testthat::expect_match(bundle_text, "guess_max = Inf", fixed = TRUE)

  overview_text <- paste(
    readLines(testthat::test_path(
      "..", "..", "10_biological_integration",
      "04_wgcna_cross_compartment_overview.R"
    ), warn = FALSE),
    collapse = "\n"
  )
  testthat::expect_match(overview_text, "guess_max = Inf", fixed = TRUE)
})
