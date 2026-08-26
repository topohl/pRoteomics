testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "manuscript_go_theme_utils.R"))
source(testthat::test_path("..", "..", "R", "gsea_wgcna_concordance_utils.R"))

testthat::test_that("canonical stage-07 export (comparison + legacy_Comparison) normalizes on comparison", {
  x <- data.frame(
    comparison = "CA1sus_CA1res", legacy_Comparison = "legacy_CA1sus_CA1res",
    ID = "GO:1", stringsAsFactors = FALSE
  )
  out <- normalize_stage07_comparison_column(x, "fixture")
  testthat::expect_identical(out$comparison, "CA1sus_CA1res")
  testthat::expect_false("legacy_Comparison" %in% names(out))
  testthat::expect_false("Comparison" %in% names(out))
})

testthat::test_that("historical pre-fix files carrying only Comparison still normalize", {
  x <- data.frame(Comparison = "CA1sus_CA1res", ID = "GO:1", stringsAsFactors = FALSE)
  out <- normalize_stage07_comparison_column(x, "fixture")
  testthat::expect_identical(out$comparison, "CA1sus_CA1res")
  testthat::expect_false("Comparison" %in% names(out))
})

testthat::test_that("a comparison identity column is required", {
  x <- data.frame(ID = "GO:1", stringsAsFactors = FALSE)
  testthat::expect_error(
    normalize_stage07_comparison_column(x, "fixture"),
    "missing a comparison identity column"
  )
})

gww_term_fixture_with <- function(comparison_column_name) {
  fixture <- data.frame(
    dataset = "microglia",
    phenotype_contrast = "SUS_vs_RES",
    spatial_unit = c("CA1", "CA1", "CA2", "CA2", "CA3"),
    program_class = "Synapse_Vesicle_Organization",
    theme_id = "synaptic_signaling_vesicle",
    manuscript_theme = "Synaptic signaling and vesicle-mediated transport",
    theme_role = "primary",
    theme_claim_eligible = TRUE,
    anchor_GO_IDs = "GO:0007268",
    mapping_method = "go_id_ontology",
    registry_version = "fixture_v1",
    theme_assignment_id = paste0("fixture_", 1:5),
    original_comparison_value = c("ca1_a", "ca1_b", "ca2_a", "ca2_b", "ca3_a"),
    ID = c("GO:2", "GO:1", "GO:3", "GO:4", "GO:5"),
    Description = c("second", "first", "third", "negative", "nonsig"),
    GO_ID = c("GO:2", "GO:1", "GO:3", "GO:4", "GO:5"),
    GO_description = c("second", "first", "third", "negative", "nonsig"),
    NES = c(1.5, 1.8, 1.4, -1.2, 2),
    pvalue = c(0.01, 0.01, 0.02, 0.03, 0.1),
    p.adjust = c(0.02, 0.01, 0.03, 0.04, 0.2),
    core_enrichment = c("P1/P2", "P1", "P2", "P3", "P4"),
    evidence_source_family = "ranked_GSEA",
    source_supplementary_file = "source.xlsx",
    stringsAsFactors = FALSE
  )
  names(fixture)[names(fixture) == "original_comparison_value"] <- comparison_column_name
  fixture
}

testthat::test_that("consumer outputs are unchanged whether the source used legacy Comparison or canonical comparison/legacy_Comparison", {
  legacy_only <- gww_term_fixture_with("Comparison")
  canonical_new <- gww_term_fixture_with("comparison")
  canonical_new$legacy_Comparison <- legacy_only$Comparison

  legacy_only <- normalize_stage07_comparison_column(legacy_only, "legacy fixture")
  canonical_new <- normalize_stage07_comparison_column(canonical_new, "canonical fixture")
  testthat::expect_identical(legacy_only$comparison, canonical_new$comparison)

  evidence_legacy <- gww_build_local_gsea_evidence(legacy_only)
  evidence_canonical <- gww_build_local_gsea_evidence(canonical_new)
  evidence_legacy <- evidence_legacy[order(evidence_legacy$gsea_evidence_id), ]
  evidence_canonical <- evidence_canonical[order(evidence_canonical$gsea_evidence_id), ]
  row.names(evidence_legacy) <- row.names(evidence_canonical) <- NULL
  testthat::expect_identical(evidence_legacy, evidence_canonical)
})
