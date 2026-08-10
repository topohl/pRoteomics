source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "sus_res_spatial_dap_atlas_utils.R"))

make_sus_res_fixture <- function(ids = c("PG1", "PG2"), genes = c("GeneA", "GeneB"),
                                 fdr = c(0.01, 0.2), effect = c(1, -1), eligible = TRUE) {
  n <- length(ids)
  data.frame(
    ProteinGroupID = ids,
    original_identifier = paste0("original_", ids),
    member_accessions = paste0("P", seq_len(n)),
    member_gene_symbols = genes,
    representative_accession = paste0("P", seq_len(n)),
    representative_gene_symbol = genes,
    representative_selection_rule = "reviewed_uniprot_status",
    protein_group_ambiguity_class = "single_accession_single_gene",
    gene_level_claim_allowed = rep(eligible, length.out = n),
    protein_level_claim_allowed = TRUE,
    mapping_status = "mapped",
    official_gene_symbol = genes,
    official_entrez_id = as.character(seq_len(n)),
    protein_group_gene_annotation_status = "concordant_official_gene",
    gene_annotation_contract_version = "mouse_gene_annotation_v1",
    uniprot_mapping_file_hash = "hash",
    orgdb_package_version = "test",
    padj = fdr,
    pval = pmin(fdr, 0.04),
    log2fc = effect,
    t = effect,
    gene_symbol = genes,
    stringsAsFactors = FALSE
  )
}

testthat::test_that("formal SUS-RES orientation is deterministic in both serializations", {
  forward <- sus_res_resolve_orientation("CA1sus_CA1res")
  reverse <- sus_res_resolve_orientation("CA1res_CA1sus")
  testthat::expect_identical(forward$formal_effect_multiplier, 1)
  testthat::expect_identical(reverse$formal_effect_multiplier, -1)
  testthat::expect_false(forward$sign_was_flipped)
  testthat::expect_true(reverse$sign_was_flipped)
})

testthat::test_that("ProteinGroupID defines counts and duplicate identifiers fail", {
  x <- make_sus_res_fixture(ids = c("PG1", "PG2"), genes = c("SameGene", "SameGene"), fdr = c(0.01, 0.01), effect = c(1, 2))
  normalized <- sus_res_normalize_contrast(x, "unitsus_unitres")
  testthat::expect_equal(sum(normalized$is_DAP_FDR05), 2L)
  duplicate <- rbind(x[1, ], x[1, ])
  testthat::expect_error(sus_res_normalize_contrast(duplicate, "unitsus_unitres"), "Duplicate ProteinGroupID")
})

testthat::test_that("BH FDR alone drives primary DAP membership and formal sign drives direction", {
  x <- make_sus_res_fixture(
    ids = c("PG1", "PG2", "PG3", "PG4"), genes = paste0("G", 1:4),
    fdr = c(0.05, 0.05001, 0.01, 0.01), effect = c(0.001, 10, -0.2, 0)
  )
  out <- sus_res_normalize_contrast(x, "unitsus_unitres")
  testthat::expect_identical(out$is_DAP_FDR05, c(TRUE, FALSE, TRUE, FALSE))
  testthat::expect_identical(out$DAP_direction[c(1, 3)], c("Higher in SUS", "Higher in RES"))
  reverse <- sus_res_normalize_contrast(x[1, ], "unitres_unitsus")
  testthat::expect_identical(reverse$DAP_direction, "Higher in RES")
})

make_similarity_unit <- function(ids, tested, dap, direction) {
  data.frame(
    ProteinGroupID = ids, is_tested = tested, is_DAP_FDR05 = dap,
    DAP_direction = direction, stringsAsFactors = FALSE
  )
}

testthat::test_that("pairwise similarity uses only the common tested universe", {
  a <- make_similarity_unit(c("PG1", "PG2"), c(TRUE, TRUE), c(TRUE, TRUE), c("Higher in SUS", "Higher in SUS"))
  b <- make_similarity_unit(c("PG1", "PG2", "PG3"), c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), c("Higher in SUS", NA, "Higher in SUS"))
  result <- sus_res_pairwise_similarity(a, b, "A", "B")
  testthat::expect_equal(result$n_tested_common, 1L)
  testthat::expect_equal(result$n_DAP_A_in_common_universe, 1L)
  testthat::expect_equal(result$n_DAP_B_in_common_universe, 1L)
  testthat::expect_equal(result$direction_aware_jaccard, 1)
})

testthat::test_that("empty DAP sets have non-estimable or zero Jaccard as prespecified", {
  empty_a <- make_similarity_unit(c("PG1", "PG2"), TRUE, FALSE, NA_character_)
  empty_b <- make_similarity_unit(c("PG1", "PG2"), TRUE, FALSE, NA_character_)
  nonempty <- make_similarity_unit(c("PG1", "PG2"), TRUE, c(TRUE, FALSE), c("Higher in SUS", NA_character_))

  empty_empty <- sus_res_pairwise_similarity(empty_a, empty_b, "A", "B")
  empty_nonempty <- sus_res_pairwise_similarity(empty_a, nonempty, "A", "B")

  testthat::expect_true(is.na(empty_empty$direction_aware_jaccard))
  testthat::expect_true(is.na(empty_empty$unsigned_protein_jaccard))
  testthat::expect_equal(empty_nonempty$direction_aware_jaccard, 0)
  testthat::expect_equal(empty_nonempty$unsigned_protein_jaccard, 0)
})

testthat::test_that("same and opposite directions remain distinct in pairwise overlap", {
  a <- make_similarity_unit(c("PG1", "PG2"), TRUE, TRUE, c("Higher in SUS", "Higher in SUS"))
  b <- make_similarity_unit(c("PG1", "PG2"), TRUE, TRUE, c("Higher in SUS", "Higher in RES"))
  result <- sus_res_pairwise_similarity(a, b, "A", "B")
  testthat::expect_equal(result$n_shared_same_direction, 1L)
  testthat::expect_equal(result$n_shared_opposite_direction, 1L)
  testthat::expect_equal(result$direction_aware_jaccard, 1 / 3)
  testthat::expect_equal(result$unsigned_protein_jaccard, 1)
  testthat::expect_equal(result$direction_concordance_among_shared, 0.5)
})

testthat::test_that("ORA queries use only eligible FDR-supported DAP genes and the universe uses all eligible tested genes", {
  x <- make_sus_res_fixture(
    ids = paste0("PG", 1:12), genes = paste0("Gene", 1:12),
    fdr = c(rep(0.01, 10), 0.5, 0.01), effect = c(rep(1, 10), 1, -1),
    eligible = c(rep(TRUE, 11), FALSE)
  )
  normalized <- sus_res_attach_gene_mapping(sus_res_normalize_contrast(x, "unitsus_unitres"))$data
  normalized$dataset <- "neuron_soma"
  normalized$dataset_label <- "Soma"
  normalized$route_unit <- "CA1_sp"
  normalized$spatial_unit <- "CA1"
  plan <- sus_res_build_ora_plan(normalized, min_ora_genes = 10L)
  key <- "neuron_soma|CA1|Higher in SUS"
  testthat::expect_setequal(plan$queries[[key]], paste0("Gene", 1:10))
  testthat::expect_setequal(plan$universes[[key]], paste0("Gene", 1:11))
  status <- plan$summary[plan$summary$direction == "Higher in SUS", ]
  testthat::expect_equal(status$n_DAP_ProteinGroupID, 10L)
  testthat::expect_equal(status$n_unique_query_genes, 10L)
  testthat::expect_identical(status$ORA_status, "eligible_for_ORA")
})

testthat::test_that("fewer than ten eligible DAP genes records the prespecified insufficient status", {
  x <- make_sus_res_fixture(ids = paste0("PG", 1:9), genes = paste0("Gene", 1:9), fdr = rep(0.01, 9), effect = rep(1, 9))
  normalized <- sus_res_attach_gene_mapping(sus_res_normalize_contrast(x, "unitsus_unitres"))$data
  normalized$dataset <- "microglia"
  normalized$dataset_label <- "Microglia-enriched ROI"
  normalized$route_unit <- "CA1_microglia"
  normalized$spatial_unit <- "CA1"
  plan <- sus_res_build_ora_plan(normalized, min_ora_genes = 10L)
  status <- plan$summary[plan$summary$direction == "Higher in SUS", "ORA_status"]
  testthat::expect_identical(status, "insufficient_DAP_genes_for_ORA")
})

testthat::test_that("GO display never fabricates unsupported terms", {
  ora <- data.frame(
    dataset = "neuron_soma", dataset_label = "Soma", spatial_unit = "CA1",
    direction = "Higher in SUS", ID = "GO:1", Description = "term",
    GeneRatio = "1/10", BgRatio = "2/100", pvalue = 0.1, p.adjust = 0.2,
    qvalue = 0.2, geneID = "Gene1", Count = 1L, stringsAsFactors = FALSE
  )
  display <- sus_res_select_go_display(ora, max_terms = 2L, fdr_threshold = 0.05)
  testthat::expect_equal(nrow(display), 0L)
  testthat::expect_true(all(c("ID", "p.adjust", "selected_for_set") %in% names(display)))
})

testthat::test_that("dataset labels and spatial units preserve anatomical and ROI semantics", {
  testthat::expect_identical(sus_res_publication_dataset_label("neuron_neuropil"), "Neuropil")
  testthat::expect_identical(sus_res_publication_dataset_label("neuron_soma"), "Soma")
  testthat::expect_identical(sus_res_publication_dataset_label("microglia"), "Microglia-enriched ROI")
  testthat::expect_identical(sus_res_spatial_unit("neuron_soma", "CA1_sp"), "CA1")
  testthat::expect_identical(sus_res_spatial_unit("microglia", "DG_microglia"), "DG")
  testthat::expect_match(dataset_interpretation("microglia"), "not purified microglia", fixed = TRUE)
})

testthat::test_that("Panel C filters canonical SUS-RES ranked-GSEA rows without changing values", {
  source <- data.frame(
    dataset_compartment = c("neuron_neuropil", "microglia", "neuron_neuropil"),
    dataset_compartment_label = c("Neuron neuropil", "Microglia", "Neuron neuropil"),
    comparison = c("SUS_vs_RES", "SUS_vs_RES", "RES_vs_CON"),
    compartment = c("neuropil", "microglia", "neuropil"),
    region = c("CA1", "CA1", "CA1"), layer = c("slm", NA, "slm"),
    theme_id = c("rna_processing_splicing_rnp", "synaptic_signaling_vesicle", "rna_processing_splicing_rnp"),
    manuscript_theme = c("RNA processing / splicing / RNP biology", "Synaptic signaling / vesicle trafficking", "RNA processing / splicing / RNP biology"),
    theme_role = c("primary", "primary", "primary"), display_order = c(1L, 4L, 1L),
    spatial_unit = c("CA1_slm", "CA1", "CA1_slm"),
    spatial_unit_label = c("CA1 SLM", "CA1", "CA1 SLM"),
    representative_term_id = c("GO:1", "GO:2", "GO:3"),
    representative_term = c("RNA term", "synapse term", "other contrast term"),
    representative_NES = c(-2.375, 1.625, 9),
    representative_FDR = c(0.01, 0.02, 0.01),
    representative_leading_edge_proteins = c("P1/P2", "P3", "P4"),
    representative_leading_edge_genes = c("G1/G2", "G3", "G4"),
    significant_term_count = c(2L, 3L, 99L),
    n_positive_sig_terms = c(0L, 3L, 99L),
    n_negative_sig_terms = c(2L, 0L, 0L),
    direction_consistency = c("negative_only", "positive_only", "positive_only"),
    n_semantic_clusters = c(1L, 2L, 50L),
    n_supported_units = c(1L, 2L, 99L), n_supported_datasets = c(1L, 2L, 3L),
    n_higher_SUS_units = c(0L, 2L, 99L), n_higher_RES_units = c(1L, 0L, 0L),
    n_higher_SUS_datasets = c(0L, 2L, 3L), n_higher_RES_datasets = c(1L, 0L, 0L),
    directional_recurrence = c("spatially_specific", "recurrent_consistent_higher_SUS", "recurrent_consistent_higher_SUS"),
    sus_res_recurrent_units = c(1L, 2L, 99L),
    sus_res_recurrent_datasets = c(1L, 2L, 3L),
    recurrence_annotation = c("spatially_specific", "recurrent_across_datasets", "recurrent_across_datasets"),
    publication_include = c(TRUE, TRUE, TRUE),
    direction = c("higher in RES", "higher in SUS", "higher in SUS"),
    source_supplementary_file = c("neuropil.xlsx", "microglia.xlsx", "neuropil.xlsx"),
    source_manifest_file = c("neuropil_manifest.csv", "microglia_manifest.csv", "neuropil_manifest.csv"),
    representative_source_comparison = c("CA1slmsus_CA1slmres", "CA1sus_CA1res", "CA1slmres_CA1slmcon"),
    representative_source_key = c("n|c1|GO:1", "m|c2|GO:2", "n|c3|GO:3"),
    representative_anchor_GO_IDs = c("GO:0006396", "GO:0099003", "GO:0006396"),
    representative_anchor_labels = c("RNA processing", "vesicle-mediated transport in synapse", "RNA processing"),
    representative_selection_rule = "p.adjust < 0.05; lowest p.adjust; largest abs(NES); GO ID; Description",
    evidence_source_family = "ranked_GSEA",
    panel_data_basis = "ranked_proteome_wide_GO_GSEA",
    recurrence_inference_role = "descriptive_only_not_a_significance_gate",
    mapping_method = "go_id_ontology", registry_version = "manuscript_go_themes_v1",
    GO_db_package_version = "3.22.0", GO_source_name = "Gene Ontology",
    GO_source_url = "http://current.geneontology.org/ontology/go-basic.obo", GO_source_date = "2025-07-22",
    relationship_types_approved = "is_a;part_of",
    stringsAsFactors = FALSE
  )
  out <- sus_res_prepare_ranked_program_panel(source)

  testthat::expect_true(all(out$comparison == "SUS_vs_RES"))
  testthat::expect_identical(out$representative_NES, c(-2.375, 1.625))
  testthat::expect_identical(out$representative_FDR, c(0.01, 0.02))
  testthat::expect_identical(out$significant_term_count, c(2L, 3L))
  testthat::expect_identical(out$dataset_compartment_label[out$dataset_compartment == "microglia"], "Microglia-enriched ROI")
  testthat::expect_true(all(out$panel_data_basis == "ranked_proteome_wide_GO_GSEA"))
  testthat::expect_true(all(!out$panel_c_derived_from_DAP_sets))
  testthat::expect_equal(out$representative_minus_log10_FDR, -log10(c(0.01, 0.02)))

  counts <- data.frame(
    dataset = c("neuron_neuropil", "microglia"), spatial_unit = c("CA1_slm", "CA1"),
    n_DAP_FDR05 = c(4L, 7L), n_higher_in_SUS = c(1L, 5L), n_higher_in_RES = c(3L, 2L),
    stringsAsFactors = FALSE
  )
  joined <- sus_res_attach_dap_counts_to_ranked_programs(out, counts)
  compact <- sus_res_biological_inspection_summary(joined)
  testthat::expect_identical(compact$n_FDR_supported_SUS_RES_DAPs, c(4L, 7L))
  testthat::expect_true(all(c(
    "representative_GO_ID", "representative_GO_term", "representative_NES",
    "representative_FDR", "direction_consistency", "leading_edge_genes",
    "sus_res_recurrent_units", "sus_res_recurrent_datasets", "source_term_key"
  ) %in% names(compact)))

  soma <- out[1, , drop = FALSE]
  soma$dataset_compartment <- "neuron_soma"
  soma$region <- "CA1"
  soma$spatial_unit <- "CA1_sp"
  soma_counts <- data.frame(
    dataset = "neuron_soma", spatial_unit = "CA1", n_DAP_FDR05 = 3L,
    n_higher_in_SUS = 1L, n_higher_in_RES = 2L, stringsAsFactors = FALSE
  )
  testthat::expect_identical(
    sus_res_attach_dap_counts_to_ranked_programs(soma, soma_counts)$n_DAP_FDR05,
    3L
  )
})

testthat::test_that("script 10 consumes script 07 source rather than executing or copying its analysis", {
  script_lines <- readLines(repo_path("04_differential_expression_enrichment", "10_sus_res_spatial_dap_atlas.r"), warn = FALSE)
  script <- paste(script_lines, collapse = "\n")
  testthat::expect_match(script, "source_data_SpatialProgramAtlas_SUS_vs_RES_publication.csv", fixed = TRUE)
  testthat::expect_false(any(grepl("^\\s*source\\(.*07_compareGO_spatial_program_atlas", script_lines)))
  testthat::expect_match(script, "Panel C is ranked proteome-wide GO/GSEA evidence", fixed = TRUE)
  testthat::expect_match(script, "theme_role", fixed = TRUE)
  testthat::expect_match(script, "GO-ID is_a/part_of ancestry", fixed = TRUE)
})

testthat::test_that("pipeline registers one global downstream producer with canonical inputs", {
  testthat::skip_if_not_installed("yaml")
  registry <- yaml::read_yaml(repo_path("pipeline.yml"))
  scripts <- registry$stages$enrichment$scripts
  hits <- which(vapply(scripts, function(x) identical(x$script, "04_differential_expression_enrichment/10_sus_res_spatial_dap_atlas.r"), logical(1)))
  testthat::expect_length(hits, 1L)
  producer <- scripts[[hits]]
  testthat::expect_identical(producer$stage, "enrichment")
  testthat::expect_identical(producer$scope, "global")
  testthat::expect_identical(unlist(producer$datasets), "global")
  testthat::expect_false(producer$recomputes_core_state)
  testthat::expect_true(producer$safe_downstream_rerun)
  testthat::expect_length(producer$consumes_required, 4L)
  testthat::expect_true(any(grepl("source_data_SpatialProgramAtlas_SUS_vs_RES_publication.csv", producer$consumes_required, fixed = TRUE)))
})
