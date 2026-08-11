source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
source(testthat::test_path("..", "..", "R", "sus_res_spatial_dap_atlas_utils.R"))
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

testthat::test_that("script 07 marks legacy heuristic behavior classes ineligible for manuscript evidence", {
  script <- readLines(testthat::test_path("..", "..", "04_differential_expression_enrichment", "07_compareGO_spatial_program_atlas.r"), warn = FALSE)
  testthat::expect_true(any(grepl("legacy_heuristic_classification = TRUE", script, fixed = TRUE)))
  testthat::expect_true(any(grepl("manuscript_evidence_eligible = FALSE", script, fixed = TRUE)))
})
