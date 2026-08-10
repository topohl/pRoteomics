source(testthat::test_path("..", "..", "R", "manuscript_go_theme_utils.R"))

registry_path <- testthat::test_path("..", "..", "config", "manuscript_go_theme_registry.tsv")

go_label <- function(id) AnnotationDbi::Term(GO.db::GOTERM[[id]])

testthat::test_that("manuscript GO-theme registry is versioned and validates live BP anchors", {
  testthat::skip_if_not_installed("GO.db")
  testthat::skip_if_not_installed("AnnotationDbi")
  registry <- read_manuscript_go_theme_registry(registry_path)
  testthat::expect_identical(unique(registry$registry_version), "manuscript_go_themes_v1")
  testthat::expect_true(all(registry$theme_role %in% c("primary", "supporting", "qc_review")))
  testthat::expect_true(all(registry$match_scope %in% c("anchor_and_descendants", "exact_go_id")))
  testthat::expect_true(all(vapply(registry$anchor_go_id, function(id) {
    identical(AnnotationDbi::Ontology(GO.db::GOTERM[[id]]), "BP")
  }, logical(1))))
})

testthat::test_that("ontology traversal admits only is_a and part_of paths", {
  testthat::skip_if_not_installed("GO.db")
  testthat::expect_true(go_bp_shortest_allowed_path("GO:0006119", "GO:0045333")$matched)
  part_of_path <- go_bp_shortest_allowed_path("GO:0030216", "GO:0043588")
  testthat::expect_true(part_of_path$matched)
  testthat::expect_match(part_of_path$relationship_path, "part_of", fixed = TRUE)

  regulatory_only <- go_bp_shortest_allowed_path("GO:0010605", "GO:0043170")
  testthat::expect_false(regulatory_only$matched)
  raw_parents <- GO.db::GOBPPARENTS[["GO:0010605"]]
  testthat::expect_true("negatively regulates" %in% names(raw_parents))
  testthat::expect_identical(manuscript_go_allowed_relationships(), c("is_a", "part_of"))
})

testthat::test_that("known false positives are unclassified and true anchors use GO ancestry", {
  testthat::skip_if_not_installed("GO.db")
  ids <- c(
    "GO:0010605", "GO:0019219", "GO:0009892", "GO:0140694", "GO:1902600",
    "GO:0006119", "GO:0042773", "GO:0042775", "GO:0022900", "GO:0022904",
    "GO:0019646", "GO:0006120", "GO:0045333", "GO:0006396", "GO:0002181", "GO:0006914"
  )
  terms <- data.frame(ID = ids, Description = vapply(ids, go_label, character(1)), stringsAsFactors = FALSE)
  mapping <- map_go_terms_to_manuscript_themes(terms, read_manuscript_go_theme_registry(registry_path))
  status <- mapping$term_status
  false_ids <- ids[seq_len(5)]
  testthat::expect_true(all(status$assignment_status[match(false_ids, status$GO_ID)] == "unclassified"))
  testthat::expect_true(all(is.na(status$manuscript_themes[match(false_ids, status$GO_ID)])))

  mito_ids <- ids[6:13]
  testthat::expect_true(all(grepl(
    "Mitochondrial respiration / oxidative phosphorylation",
    status$manuscript_themes[match(mito_ids, status$GO_ID)], fixed = TRUE
  )))
  testthat::expect_identical(status$manuscript_themes[match("GO:0006396", status$GO_ID)], "RNA processing / splicing / RNP biology")
  testthat::expect_identical(status$manuscript_themes[match("GO:0002181", status$GO_ID)], "Ribosome / translation")
  testthat::expect_identical(status$manuscript_themes[match("GO:0006914", status$GO_ID)], "Autophagy / lysosome / endosome")
  testthat::expect_true(all(mapping$assignments$mapping_method == "go_id_ontology"))
})

testthat::test_that("epidermal terms are QC review and multi-parent theme status remains explicit", {
  testthat::skip_if_not_installed("GO.db")
  ids <- c("GO:0043588", "GO:0008544", "GO:0030216", "GO:0009913", "GO:0006364")
  mapping <- map_go_terms_to_manuscript_themes(
    data.frame(ID = ids, Description = vapply(ids, go_label, character(1)), stringsAsFactors = FALSE),
    read_manuscript_go_theme_registry(registry_path)
  )
  status <- mapping$term_status
  epidermal <- status[status$GO_ID %in% ids[1:4], , drop = FALSE]
  testthat::expect_true(all(epidermal$assignment_status == "qc_review"))
  testthat::expect_true(all(grepl("Epithelial / epidermal", epidermal$manuscript_themes, fixed = TRUE)))
  testthat::expect_false(any(grepl("neural", epidermal$manuscript_themes, ignore.case = TRUE)))
  rrna <- status[status$GO_ID == "GO:0006364", , drop = FALSE]
  testthat::expect_identical(rrna$assignment_status, "multi_theme")
  testthat::expect_equal(rrna$n_manuscript_themes, 2L)
  long <- collapse_go_theme_assignment_audit(mapping)
  testthat::expect_identical(anyDuplicated(paste(long$GO_ID, long$theme_id, sep = "|")), 0L)
})

theme_fixture <- function() {
  ids <- c("GO:0045333", "GO:0006119", "GO:0042773", "GO:0010605")
  data.frame(
    dataset = "neuron_neuropil", dataset_label = "Neuron neuropil",
    compartment = "neuropil", region = c("CA1", "CA1", "CA3", "CA2"),
    layer = c("slm", "slm", "sr", "sr"), spatial_unit = c("CA1_slm", "CA1_slm", "CA3_sr", "CA2_sr"),
    phenotype_contrast = "SUS_vs_RES", Comparison = c("CA1slmsus_CA1slmres", "CA1slmsus_CA1slmres", "CA3srsus_CA3srres", "CA2srsus_CA2srres"),
    ID = ids, Description = vapply(ids, go_label, character(1)),
    NES = c(-1.8, -2.1, 2.2, 4.0), p.adjust = c(0.02, 0.01, 0.001, 0.0001),
    core_enrichment = c("P1/P2", "P2/P3", "P4", "P5"), core_enrichment_gene = c("G1/G2", "G2/G3", "G4", "G5"),
    program_class = c("Other", "Mitochondria_OXPHOS_Metabolism", "Mitochondria_OXPHOS_Metabolism", "Mitochondria_OXPHOS_Metabolism"),
    source_supplementary_file = "ranked.xlsx", source_manifest_file = "manifest.csv",
    evidence_source_family = "ranked_GSEA", stringsAsFactors = FALSE
  )
}

testthat::test_that("supported unclassified terms are retained and representative/directional recurrence are descriptive", {
  testthat::skip_if_not_installed("GO.db")
  fixture <- theme_fixture()
  mapping <- map_go_terms_to_manuscript_themes(
    unique(fixture[c("ID", "Description")]), read_manuscript_go_theme_registry(registry_path)
  )
  audit <- build_sus_res_supported_go_theme_audit(fixture, mapping)
  testthat::expect_equal(nrow(audit), nrow(fixture))
  false_row <- audit[audit$GO_ID == "GO:0010605", , drop = FALSE]
  testthat::expect_identical(false_row$assignment_status, "unclassified")
  testthat::expect_identical(false_row$mapping_method, "go_id_unclassified")

  summary <- build_sus_res_manuscript_theme_summary(audit, mapping)
  ca1 <- summary[summary$spatial_unit == "CA1_slm", , drop = FALSE]
  testthat::expect_identical(ca1$representative_term_id, "GO:0006119")
  testthat::expect_equal(ca1$representative_NES, -2.1)
  testthat::expect_equal(ca1$representative_FDR, 0.01)
  testthat::expect_identical(ca1$direction_consistency, "negative_only")
  testthat::expect_identical(ca1$directional_recurrence, "recurrent_mixed_direction")
  testthat::expect_equal(ca1$n_higher_RES_units, 1L)
  testthat::expect_equal(ca1$n_higher_SUS_units, 1L)
  testthat::expect_true(all(summary$evidence_source_family == "ranked_GSEA"))
})

testthat::test_that("Wang semantic redundancy QA is descriptive and preserves original representatives", {
  testthat::skip_if_not_installed("GOSemSim")
  testthat::skip_if_not_installed("org.Mm.eg.db")
  fixture <- theme_fixture()[1:3, , drop = FALSE]
  mapping <- map_go_terms_to_manuscript_themes(
    unique(fixture[c("ID", "Description")]), read_manuscript_go_theme_registry(registry_path)
  )
  audit <- build_sus_res_supported_go_theme_audit(fixture, mapping)
  semantic <- go_semantic_redundancy_qa(audit, cutoff = 0.70)
  testthat::expect_equal(nrow(semantic$term_audit), length(unique(audit$GO_ID)))
  testthat::expect_true(all(semantic$term_audit$semantic_similarity_method == "Wang_BP"))
  testthat::expect_true(all(semantic$term_audit$semantic_inference_role == "descriptive_redundancy_QA_only"))
  testthat::expect_true(all(semantic$pair_audit$semantic_similarity_cutoff == 0.70))
})
