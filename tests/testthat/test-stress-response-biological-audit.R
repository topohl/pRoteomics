source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
source(testthat::test_path("..", "..", "R", "sus_res_spatial_dap_atlas_utils.R"))
source(testthat::test_path("..", "..", "R", "manuscript_go_theme_utils.R"))
source(testthat::test_path("..", "..", "R", "stress_response_biological_audit_utils.R"))

testthat::test_that("three formal contrast signs are deterministic", {
  testthat::expect_identical(stress_response_parse_comparison("unitres_unitcon")$contrast, "RES_vs_CON")
  testthat::expect_identical(stress_response_parse_comparison("unitsus_unitcon")$contrast, "SUS_vs_CON")
  testthat::expect_identical(stress_response_parse_comparison("unitsus_unitres")$contrast, "SUS_vs_RES")
  reversed <- stress_response_parse_comparison("unitcon_unitres")
  testthat::expect_identical(reversed$contrast, "RES_vs_CON")
  testthat::expect_identical(reversed$formal_effect_multiplier, -1)
  testthat::expect_true(reversed$sign_was_flipped)
})

make_protein_contrast <- function(contrast, effect, p) {
  data.frame(
    dataset = "neuron_soma", route_unit = "CA1_sp", spatial_unit = "CA1_sp",
    contrast = contrast, comparison = contrast, ProteinGroupID = c("PG1", "PG2"),
    log2FC = effect, p_value = p, canonical_contrast_FDR = stats::p.adjust(p, "BH"),
    canonical_FDR_support = stats::p.adjust(p, "BH") < 0.05,
    raw_p_input_field = "pval", canonical_FDR_field = "padj", effect_input_field = "log2fc",
    formal_effect_multiplier = 1, formal_effect_definition = contrast,
    source_file = "fixture.csv", source_hash = "fixture", stringsAsFactors = FALSE
  )
}

testthat::test_that("protein joint BH uses the matched 2m raw-p family and algebra is QC only", {
  res <- make_protein_contrast("RES_vs_CON", c(1, -2), c(0.001, 0.4))
  sus <- make_protein_contrast("SUS_vs_CON", c(3, -1), c(0.002, 0.5))
  sr <- make_protein_contrast("SUS_vs_RES", sus$log2FC - res$log2FC, c(0.003, 0.6))
  out <- stress_response_build_protein_outputs(rbind(res, sus, sr))
  expected <- stats::p.adjust(c(res$p_value, sus$p_value), "BH")
  testthat::expect_equal(out$geometry$control_pair_joint_FDR_RES_vs_CON, expected[1:2])
  testthat::expect_equal(out$geometry$control_pair_joint_FDR_SUS_vs_CON, expected[3:4])
  testthat::expect_identical(out$geometry$control_pair_family_size, rep(4L, 2))
  testthat::expect_true(out$algebra_audit$tolerance_pass)
  testthat::expect_equal(out$algebra_audit$max_absolute_residual, 0)
  testthat::expect_identical(out$geometry$derived_geometry_inference_role, rep("descriptive_only_no_p_value", 2))
  testthat::expect_false(out$joint_family_audit$adjusted_p_used_as_input)
})

testthat::test_that("protein joint BH refuses an empty matched universe", {
  res <- make_protein_contrast("RES_vs_CON", c(1, -2), c(0.001, 0.4))
  sus <- make_protein_contrast("SUS_vs_CON", c(3, -1), c(0.002, 0.5))
  sus$ProteinGroupID <- c("PG3", "PG4")
  sr <- make_protein_contrast("SUS_vs_RES", c(2, 1), c(0.003, 0.6))
  testthat::expect_error(stress_response_build_protein_outputs(rbind(res, sus, sr)), "No common valid-p ProteinGroupID universe")
})

testthat::test_that("direct DAP control geometry uses only canonical SUS-vs-RES support and retains zero-DAP contexts", {
  fixture <- data.frame(
    dataset = c("neuron_soma", "neuron_soma", "neuron_soma", "microglia"),
    spatial_unit = c("CA1_sp", "CA1_sp", "CA1_sp", "CA2"),
    ProteinGroupID = c("PG1", "PG2", "PG3", "PG4"),
    log2FC_RES_vs_CON = c(1, -2, 10, 0.5),
    log2FC_SUS_vs_CON = c(-1, 1, -10, 0.2),
    abs_log2FC_RES_vs_CON = c(1, 2, 10, 0.5),
    abs_log2FC_SUS_vs_CON = c(1, 1, 10, 0.2),
    delta_abs_log2FC = c(0, 1, 0, 0.3),
    canonical_FDR_support_SUS_vs_RES = c(TRUE, TRUE, FALSE, FALSE),
    canonical_FDR_support_RES_vs_CON = c(FALSE, FALSE, TRUE, TRUE),
    canonical_FDR_support_SUS_vs_CON = c(FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  out <- stress_response_direct_dap_control_geometry(fixture)
  testthat::expect_identical(out$detail$ProteinGroupID, c("PG1", "PG2"))
  testthat::expect_identical(out$detail$sign_geometry, c("RES_up_SUS_down", "RES_down_SUS_up"))
  testthat::expect_identical(out$detail$relative_displacement, c("equal_distance", "RES_farther"))
  testthat::expect_true(all(out$detail$opposite_side_of_CON))
  testthat::expect_true(all(out$detail$direct_DAP_selection_basis == "canonical SUS_vs_RES FDR < 0.05"))
  testthat::expect_false(any(grepl("^p_|p[_\\.]?value|pvalue", names(out$detail), ignore.case = TRUE)))
  zero <- out$context_summary[out$context_summary$dataset == "microglia" & out$context_summary$spatial_unit == "CA2", , drop = FALSE]
  testthat::expect_identical(zero$n_direct_DAP, 0L)
  testthat::expect_true(is.na(zero$fraction_opposite_side_of_CON))
  testthat::expect_true(is.na(zero$median_delta_abs_log2FC))
  testthat::expect_false(any(grepl("p[_\\.]?value|pvalue", names(out$context_summary), ignore.case = TRUE)))
})

testthat::test_that("GO joint BH uses nominal pvalue and excludes SUS-vs-RES", {
  fixture <- do.call(rbind, lapply(STRESS_RESPONSE_CONTRASTS, function(contrast) {
    p <- if (contrast == "RES_vs_CON") c(0.001, 0.3) else if (contrast == "SUS_vs_CON") c(0.004, 0.2) else c(0.01, 0.1)
    data.frame(
      dataset = "microglia", spatial_unit = "CA1", phenotype_contrast = contrast,
      ID = c("GO:0000001", "GO:0000002"), pvalue = p, p.adjust = stats::p.adjust(p, "BH"),
      stringsAsFactors = FALSE
    )
  }))
  out <- stress_response_add_joint_go_fdr(fixture)
  expected <- stats::p.adjust(c(0.001, 0.3, 0.004, 0.2), "BH")
  testthat::expect_equal(out$joint$control_pair_joint_GO_FDR[out$joint$phenotype_contrast == "RES_vs_CON"], expected[1:2])
  testthat::expect_equal(out$joint$control_pair_joint_GO_FDR[out$joint$phenotype_contrast == "SUS_vs_CON"], expected[3:4])
  testthat::expect_false(any(out$joint$phenotype_contrast == "SUS_vs_RES"))
  testthat::expect_false(out$audit$adjusted_p_used_as_input)
  testthat::expect_identical(out$audit$family_size, 4L)
})

testthat::test_that("theme trajectory keeps continuous geometry separate from direct support", {
  detail <- do.call(rbind, lapply(STRESS_RESPONSE_CONTRASTS, function(contrast) {
    data.frame(
      dataset = "neuron_neuropil", dataset_label = "Neuron neuropil", compartment = "neuropil",
      region = "CA1", layer = "slm", spatial_unit = "CA1_slm", contrast = contrast,
      theme_id = "mitochondrial_respiration_oxphos", manuscript_theme = "Mitochondrial respiration / OXPHOS",
      theme_role = "primary", display_order = 4L,
      median_NES_all_theme_terms = c(1.2, 0.5, -0.7)[match(contrast, STRESS_RESPONSE_CONTRASTS)],
      canonical_GO_FDR_support = contrast == "SUS_vs_RES",
      min_canonical_GO_FDR = if (contrast == "SUS_vs_RES") 0.01 else 0.2,
      n_canonical_FDR05_GO_terms = if (contrast == "SUS_vs_RES") 1L else 0L,
      control_pair_joint_GO_FDR_support = FALSE, min_control_pair_joint_GO_FDR = NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  out <- stress_response_build_theme_trajectory(detail)
  testthat::expect_equal(out$delta_abs_theme_NES, 0.7)
  testthat::expect_true(out$same_direction_RES_and_SUS_vs_CON)
  testthat::expect_true(out$RES_farther_from_CON_descriptively)
  testthat::expect_true(out$canonical_GO_FDR_support_SUS_vs_RES)
  testthat::expect_match(out$trajectory_inference_role, "descriptive_geometry", fixed = TRUE)
  testthat::expect_match(out$NES_algebra_warning, "not algebraically additive", fixed = TRUE)
})

testthat::test_that("theme overview reconstructs delta absolute NES magnitude summaries", {
  detail <- do.call(rbind, lapply(c("CA1_sp", "CA2_sp"), function(unit) {
    do.call(rbind, lapply(STRESS_RESPONSE_CONTRASTS, function(contrast) {
      data.frame(
        dataset = "neuron_soma", spatial_unit = unit, contrast = contrast,
        theme_id = "theme_a", manuscript_theme = "Theme A", theme_role = "primary", display_order = 1L,
        canonical_GO_FDR_support = FALSE, control_pair_joint_GO_FDR_support = FALSE,
        stringsAsFactors = FALSE
      )
    }))
  }))
  trajectory <- data.frame(
    dataset = "neuron_soma", spatial_unit = c("CA1_sp", "CA2_sp"), theme_id = "theme_a",
    abs_median_NES_RES_vs_CON = c(1.2, 0.5), abs_median_NES_SUS_vs_CON = c(0.5, 2.0),
    delta_abs_theme_NES = c(0.7, -1.5), same_direction_RES_and_SUS_vs_CON = c(TRUE, FALSE),
    canonical_GO_FDR_support_SUS_vs_RES = c(FALSE, TRUE), stringsAsFactors = FALSE
  )
  out <- stress_response_theme_overview(detail, trajectory)
  testthat::expect_equal(out$median_delta_abs_theme_NES, median(trajectory$delta_abs_theme_NES))
  testthat::expect_equal(out$q25_delta_abs_theme_NES, unname(quantile(trajectory$delta_abs_theme_NES, 0.25)))
  testthat::expect_equal(out$q75_delta_abs_theme_NES, unname(quantile(trajectory$delta_abs_theme_NES, 0.75)))
  testthat::expect_equal(out$median_abs_delta_theme_NES, median(abs(trajectory$delta_abs_theme_NES)))
  testthat::expect_equal(out$min_delta_abs_theme_NES, min(trajectory$delta_abs_theme_NES))
  testthat::expect_equal(out$max_delta_abs_theme_NES, max(trajectory$delta_abs_theme_NES))
})

testthat::test_that("supported GO occurrence grain stays unique through multi-theme expansion", {
  go <- data.frame(
    dataset = "neuron_soma", spatial_unit = "CA1_sp", phenotype_contrast = "SUS_vs_RES",
    ID = c("GO:0006364", "GO:9999999"), Description = c("rRNA processing", "unclassified fixture"),
    pvalue = c(0.001, 0.002), p.adjust = c(0.01, 0.02), NES = c(2, -2), stringsAsFactors = FALSE
  )
  assignment <- data.frame(
    GO_ID = "GO:0006364", GO_description = "rRNA processing",
    theme_id = c("rna_processing_splicing_rnp", "ribosome_translation"),
    manuscript_theme = c("RNA", "Translation"), theme_role = "primary", display_order = c(1L, 2L),
    anchor_GO_ID = c("GO:0006396", "GO:0042254"), anchor_label = c("RNA processing", "ribosome biogenesis"),
    match_scope = "anchor_and_descendants", match_type = "descendant", path_length = 1L,
    relationship_path = "is_a", GO_path = "fixture", relationship_types_approved = "is_a;part_of",
    mapping_method = "go_id_ontology", registry_version = "fixture_v1", GO_db_package_version = "fixture",
    GO_source_name = "fixture", GO_source_url = "fixture", GO_source_date = "2026-08-11",
    stringsAsFactors = FALSE
  )
  status <- data.frame(
    GO_ID = c("GO:0006364", "GO:9999999"), assignment_status = c("multi_theme", "unclassified"),
    manuscript_themes = c("RNA;Translation", NA), theme_roles = c("primary", NA),
    registry_version = "fixture_v1", mapping_method = c("go_id_ontology", "go_id_unclassified"),
    stringsAsFactors = FALSE
  )
  out <- stress_response_build_supported_go_audits(go, list(assignments = assignment, term_status = status))
  testthat::expect_identical(anyDuplicated(out$occurrence$supported_occurrence_id), 0L)
  testthat::expect_equal(nrow(out$occurrence), 2L)
  testthat::expect_equal(nrow(out$exploded), 3L)
  testthat::expect_equal(length(unique(out$exploded$supported_occurrence_id)), nrow(out$occurrence))
  testthat::expect_identical(anyDuplicated(out$exploded$theme_assignment_row_id), 0L)
  testthat::expect_true(any(grepl("unclassified$", out$exploded$theme_assignment_row_id)))
  testthat::expect_true(all(out$exploded$row_grain == "FDR_supported_GO_occurrence_x_theme_assignment"))
})

testthat::test_that("protected workbook is a deterministic reference artifact, not a Stage-11 input", {
  tmp <- tempfile(fileext = ".xlsx")
  writeLines("protected fixture", tmp, useBytes = TRUE)
  on.exit(unlink(tmp), add = TRUE)
  first <- stress_response_protected_reference_artifacts(
    tmp, producer = "stage10", note = "Stage 11 does not read or modify it.", repository_root = dirname(tmp)
  )
  second <- stress_response_protected_reference_artifacts(
    tmp, producer = "stage10", note = "Stage 11 does not read or modify it.", repository_root = dirname(tmp)
  )
  testthat::expect_identical(first, second)
  testthat::expect_true(first$exists)
  testthat::expect_identical(first$role, "protected_reference_not_consumed")
  testthat::expect_identical(first$sha256_at_stage11_run, file_hash_sha256(tmp))
  script <- readLines(testthat::test_path("..", "..", "04_differential_expression_enrichment", "11_stress_response_biological_audit.r"), warn = FALSE)
  testthat::expect_false(any(grepl("sus_res_workbook_protected =", script, fixed = TRUE)))
  testthat::expect_true(any(grepl("protected_reference_artifacts.csv", script, fixed = TRUE)))
})

testthat::test_that("script 07 marks legacy heuristic behavior classes ineligible for manuscript evidence", {
  script <- readLines(testthat::test_path("..", "..", "04_differential_expression_enrichment", "07_compareGO_spatial_program_atlas.r"), warn = FALSE)
  testthat::expect_true(any(grepl("legacy_heuristic_classification = TRUE", script, fixed = TRUE)))
  testthat::expect_true(any(grepl("manuscript_evidence_eligible = FALSE", script, fixed = TRUE)))
})
