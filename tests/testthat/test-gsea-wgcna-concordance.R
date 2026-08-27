testthat::local_edition(3)

source(testthat::test_path(
  "..", "..", "R", "gsea_wgcna_concordance_utils.R"
))

gww_term_fixture <- function() {
  data.frame(
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
    comparison = c("ca1_a", "ca1_b", "ca2_a", "ca2_b", "ca3_a"),
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
}

testthat::test_that("formal contrast signs and labels are fixed", {
  testthat::expect_identical(
    gww_formal_contrast(names(gww_formal_contrast_map())),
    unname(gww_formal_contrast_map())
  )
  testthat::expect_identical(
    gww_formal_contrast(unname(gww_formal_contrast_map())),
    unname(gww_formal_contrast_map())
  )
  testthat::expect_true(is.na(gww_formal_contrast("RES_vs_SUS")))
})

testthat::test_that("GSEA representatives and recurrent-cross-spatial direction are deterministic", {
  terms <- gww_term_fixture()
  local_a <- gww_build_local_gsea_evidence(terms)
  local_b <- gww_build_local_gsea_evidence(terms[sample(nrow(terms)), ])
  local_a <- local_a[order(local_a$gsea_evidence_id), ]
  local_b <- local_b[order(local_b$gsea_evidence_id), ]
  row.names(local_a) <- row.names(local_b) <- NULL
  testthat::expect_identical(local_a, local_b)
  ca1 <- local_a[local_a$gsea_spatial_unit == "ca1", ]
  testthat::expect_identical(ca1$GO_ID, "GO:1")
  testthat::expect_identical(ca1$positive_GO_ID, "GO:1")
  testthat::expect_identical(ca1$gsea_direction_status, "positive_only")
  ca2 <- local_a[local_a$gsea_spatial_unit == "ca2", ]
  testthat::expect_identical(ca2$gsea_direction_status, "mixed_direction")
  testthat::expect_true(is.na(ca2$gsea_direction_sign))

  recurrent <- gww_build_recurrent_cross_spatial_gsea_evidence(terms, local_a, 2L)
  testthat::expect_equal(nrow(recurrent), 0L)

  terms$NES[terms$comparison == "ca2_b"] <- 1.2
  local <- gww_build_local_gsea_evidence(terms)
  recurrent <- gww_build_recurrent_cross_spatial_gsea_evidence(terms, local, 2L)
  testthat::expect_equal(nrow(recurrent), 1L)
  testthat::expect_identical(
    recurrent$gsea_direction_status, "positive_recurrent_cross_spatial"
  )
  testthat::expect_identical(recurrent$gsea_direction_sign, 1L)
  testthat::expect_identical(recurrent$comparison_scope, "recurrent_cross_spatial")
  testthat::expect_identical(recurrent$gsea_spatial_unit, "recurrent_cross_spatial")
  testthat::expect_equal(recurrent$n_units_supporting_same_direction, 2L)
})

testthat::test_that("program matching is exact and excludes aliases", {
  registry <- data.frame(
    annotation_token_normalized = "translation / proteostasis",
    registry_annotation_token = "translation / proteostasis",
    biological_program = "Ribosome_Translation",
    mapping_rationale = "fixture",
    stringsAsFactors = FALSE
  )
  handoff <- data.frame(
    dataset = rep("microglia", 2),
    entity_level = rep("module", 2),
    entity_id = rep("WGCNA_m01", 2),
    display_label = rep("Translation", 2),
    module_program_primary = rep("translation / proteostasis", 2),
    module_program_secondary = NA_character_,
    cleaned_biological_label = rep("Translation", 2),
    independent_hypothesis = TRUE,
    claim_entity_role = "canonical_module",
    contrast = c("RES - CON", "SUS - RES"),
    stringsAsFactors = FALSE
  )
  match <- gww_build_entity_program_matches(handoff, registry)
  testthat::expect_equal(nrow(match), 1L)
  testthat::expect_identical(match$annotation_field, "module_program_primary")
  testthat::expect_identical(match$biological_program, "Ribosome_Translation")

  alias <- handoff
  alias$independent_hypothesis <- FALSE
  alias$claim_entity_role <- "compatibility_alias"
  testthat::expect_error(
    gww_build_entity_program_matches(alias, registry),
    "independent non-alias"
  )
})

testthat::test_that("overlap uses canonical ProteinGroupID universe and exact families", {
  universe <- data.frame(
    ProteinGroupID = paste0("PG:microglia:PG", 1:4),
    included_in_wgcna = TRUE,
    ModuleID = c("M1", "M1", "M2", "M2"),
    member_accessions = paste0("P", 1:4),
    RepresentativeUniProt = paste0("P", 1:4),
    MemberUniProts = paste0("P", 1:4),
    representative_accession = paste0("P", 1:4),
    stringsAsFactors = FALSE
  )
  bundle <- gww_build_universe_bundle(universe, "microglia")
  evidence <- data.frame(
    gsea_evidence_id = c("E1", "E2"),
    dataset = "microglia",
    phenotype_contrast = "SUS_vs_RES",
    contrast = "SUS - RES",
    comparison_scope = "local_local",
    gsea_spatial_unit = "ca1",
    biological_program = c("A", "B"),
    GO_ID = c("GO:1", "GO:2"),
    GO_description = c("A", "B"),
    NES = c(1, -1),
    GSEA_FDR = c(0.01, 0.02),
    gsea_source_key = c("K1", "K2"),
    leading_edge_proteins = c("PG:microglia:PG1;PG:microglia:PG3", "PG:microglia:PG3;PG:microglia:PG4"),
    stringsAsFactors = FALSE
  )
  matches <- data.frame(
    dataset = "microglia",
    entity_level = "module",
    entity_id = c("M1", "M2"),
    biological_program = c("A", "B"),
    annotation_field = "module_program_primary",
    annotation_token = c("A", "B"),
    program_module_match_rule = "fixture",
    stringsAsFactors = FALSE
  )
  observed <- gww_program_specific_overlap(
    evidence, matches, list(microglia = bundle)
  )
  testthat::expect_equal(nrow(observed), 2L)
  testthat::expect_true(all(observed$n_universe == 4L))
  testthat::expect_true(all(observed$overlap_family_size == 2L))
  testthat::expect_equal(length(unique(observed$overlap_family_id)), 1L)
  testthat::expect_equal(observed$n_overlap[observed$entity_id == "M1"], 1L)
  testthat::expect_true(all(observed$overlap_universe ==
    "canonical_dataset_WGCNA_feature_universe"))
})

testthat::test_that("low-n concordance classes use effect CI and stability", {
  fixture <- data.frame(
    gsea_direction_sign = c(1L, 1L, 1L, 1L, 1L, NA_integer_),
    wgcna_estimate = c(0.8, 0.6, 0.1, -0.7, 0.6, 0.8),
    wgcna_CI_low = c(0.2, -0.3, -0.4, -1.0, -0.3, 0.2),
    wgcna_CI_high = c(1.4, 1.5, 0.6, -0.2, 1.5, 1.4),
    wgcna_tier_specific_fdr = c(0.01, 0.4, 0.8, 0.2, 0.4, 0.01),
    wgcna_model_valid = TRUE,
    wgcna_model_stability_status = "stable_mixed_model",
    animal_instability_flag = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
    effect_scale_q25_abs = 0.2,
    effect_scale_median_abs = 0.5,
    stringsAsFactors = FALSE
  )
  first <- gww_classify_concordance(fixture)
  second <- gww_classify_concordance(fixture)
  testthat::expect_identical(first$concordance_class, second$concordance_class)
  testthat::expect_identical(
    first$concordance_class,
    c(
      "FDR_supported_concordance", "concordant_imprecise",
      "weak_or_near_zero_module_support", "discordant",
      "animal_sensitive", "unresolved"
    )
  )
})

testthat::test_that("local and recurrent-cross-spatial scopes fail closed when mismatched", {
  good <- data.frame(
    comparison_scope = c("local_local", "recurrent_cross_spatial"),
    gsea_spatial_unit = c("ca1", "recurrent_cross_spatial"),
    wgcna_spatial_unit = c("ca1", "global_spatial_adjusted"),
    wgcna_effect_scope = c("within_spatial_unit", "spatial_adjusted_global"),
    wgcna_analysis_tier = c(
      "exploratory_spatial_localization", "primary_wgcna_global"
    ),
    stringsAsFactors = FALSE
  )
  testthat::expect_silent(gww_validate_local_recurrent_cross_spatial_semantics(good))
  bad <- good
  bad$wgcna_spatial_unit[[1]] <- "ca2"
  testthat::expect_error(
    gww_validate_local_recurrent_cross_spatial_semantics(bad),
    "Local/recurrent-cross-spatial concordance semantics"
  )
  bad <- good
  bad$comparison_scope[[2]] <- "local_local"
  testthat::expect_error(
    gww_validate_local_recurrent_cross_spatial_semantics(bad),
    "Local/recurrent-cross-spatial concordance semantics"
  )
})

testthat::test_that("adaptive pattern requires effects and direct contrast", {
  long <- data.frame(
    dataset = "neuron_neuropil",
    biological_program = "Synapse_Vesicle_Organization",
    contrast = c("RES - CON", "SUS - CON", "SUS - RES"),
    comparison_scope = "recurrent_cross_spatial",
    wgcna_entity_level = "module",
    wgcna_entity_id = "WGCNA_m01",
    gsea_direction_sign = c(1L, -1L, -1L),
    gsea_direction_status = c(
      "positive_recurrent_cross_spatial", "negative_recurrent_cross_spatial",
      "negative_recurrent_cross_spatial"
    ),
    wgcna_estimate = c(0.5, -0.4, -0.9),
    wgcna_CI_low = c(-0.1, -1.0, -1.5),
    wgcna_CI_high = c(1.1, 0.2, -0.3),
    concordance_class = c(
      "concordant_imprecise", "concordant_imprecise",
      "FDR_supported_concordance"
    ),
    stringsAsFactors = FALSE
  )
  observed <- gww_adaptive_pattern_summary(long)
  testthat::expect_identical(
    observed$adaptive_resilience_pattern,
    "RES_supported_shift_candidate"
  )
  testthat::expect_true(observed$RES_gt_CON_gt_SUS_supported)

  no_direct <- long[long$contrast != "SUS - RES", ]
  observed_no_direct <- gww_adaptive_pattern_summary(no_direct)
  testthat::expect_identical(
    observed_no_direct$adaptive_resilience_pattern, "unresolved"
  )
})

testthat::test_that("utility adjusts only the new overlap family", {
  text <- paste(readLines(testthat::test_path(
    "..", "..", "R", "gsea_wgcna_concordance_utils.R"
  ), warn = FALSE), collapse = "\n")
  hits <- gregexpr("p\\.adjust\\s*\\(", text)[[1]]
  testthat::expect_equal(sum(hits > 0), 1L)
  testthat::expect_match(text, "overlap_FDR = stats::p.adjust", fixed = TRUE)
  testthat::expect_false(grepl("combined_p|combined_fdr", text, ignore.case = TRUE))
})
