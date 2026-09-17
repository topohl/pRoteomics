testthat::test_that("biological integration entrypoints exist and dry-run", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  scripts <- c(
    "analysis/04_differential_abundance/compare_external_stress_signatures.R",
    "analysis/08_integration/test_module_behaviour_coupling.R",
    "analysis/08_integration/build_cross_compartment_atlas.R",
    "analysis/08_integration/summarize_programs_for_manuscript.R",
    "analysis/08_integration/build_evidence_priority_matrix.R",
    "analysis/08_integration/test_enrichment_module_concordance.R",
    "analysis/08_integration/summarize_enrichment_module_concordance.R",
    "analysis/08_integration/build_candidate_protein_shortlist.R",
    "analysis/08_integration/quantify_candidate_network_position.R",
    "analysis/05_wgcna/audit_module_label_coherence.R",
    "analysis/05_wgcna/adjudicate_module_labels.R",
    "analysis/05_wgcna/build_module_label_registry.R"
  )
  testthat::expect_true(all(file.exists(repo_path(scripts))))

  cmd <- file.path(R.home("bin"), "Rscript")
  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  cases <- list(
    c("analysis/04_differential_abundance/compare_external_stress_signatures.R", "--dry-run"),
    c("analysis/08_integration/test_module_behaviour_coupling.R", "--dataset", "microglia", "--dry-run"),
    c("analysis/08_integration/build_cross_compartment_atlas.R", "--dry-run"),
    c("analysis/08_integration/summarize_programs_for_manuscript.R", "--dry-run"),
    c("analysis/08_integration/build_evidence_priority_matrix.R", "--dry-run"),
    c("analysis/08_integration/test_enrichment_module_concordance.R", "--dataset", "all", "--dry-run"),
    c("analysis/08_integration/summarize_enrichment_module_concordance.R", "--dataset", "all", "--dry-run"),
    c("analysis/08_integration/build_candidate_protein_shortlist.R", "--dataset", "all", "--dry-run"),
    c("analysis/08_integration/quantify_candidate_network_position.R", "--dataset", "all", "--dry-run"),
    c("analysis/05_wgcna/audit_module_label_coherence.R", "--dataset", "all", "--dry-run"),
    c("analysis/05_wgcna/adjudicate_module_labels.R", "--dataset", "all", "--dry-run"),
    c("analysis/05_wgcna/build_module_label_registry.R", "--dataset", "all", "--dry-run")
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
  source(repo_path("R", "integration_utils.R"))
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
  source(repo_path("R", "integration_utils.R"))
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
    readLines(repo_path("R", "integration_utils.R"), warn = FALSE),
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
    readLines(repo_path("R", "evidence_bundle_utils.R"), warn = FALSE),
    collapse = "\n"
  )
  testthat::expect_match(bundle_text, "guess_max = Inf", fixed = TRUE)

  overview_text <- paste(
    readLines(repo_path("analysis/08_integration",
      "summarize_module_cross_compartment.R"
    ), warn = FALSE),
    collapse = "\n"
  )
  testthat::expect_match(overview_text, "guess_max = Inf", fixed = TRUE)
})
