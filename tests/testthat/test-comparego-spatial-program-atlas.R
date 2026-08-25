source(testthat::test_path("..", "..", "R", "paths.R"))

load_spatial_atlas_helpers <- function() {
  env <- new.env(parent = globalenv())
  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  sys.source(
    repo_path("04_differential_expression_enrichment", "07_compareGO_spatial_program_atlas.r"),
    envir = env
  )
  env
}

testthat::test_that("spatial atlas classification reuses the canonical biological-program vocabulary", {
  helpers <- load_spatial_atlas_helpers()
  descriptions <- c(
    "mitochondrial oxidative phosphorylation",
    "autophagy of mitochondrion",
    "unmapped example term"
  )
  canonical <- map_terms_to_programs(
    data.frame(Description = descriptions, stringsAsFactors = FALSE),
    "Description"
  )$biological_program
  observed <- as.character(helpers$classify_program(descriptions))

  testthat::expect_identical(observed[1:2], canonical[1:2])
  testthat::expect_identical(observed[[3]], "Other")
  testthat::expect_identical(helpers$canonical_program_levels, biological_program_patterns()$biological_program)
  testthat::expect_false(exists("program_dictionary", envir = helpers, inherits = FALSE))
})

testthat::test_that("program summary selects an existing FDR-supported representative term deterministically", {
  helpers <- load_spatial_atlas_helpers()
  fixture <- data.frame(
    dataset = "neuron_soma", dataset_label = "Neuron soma", region = "CA1", layer = NA_character_,
    comparison = "CA1spsus_CA1spres", route_category = "phenotype_within_unit", route_unit = "CA1_sp",
    spatial_unit = "CA1", compartment = "soma", phenotype_contrast = "SUS_vs_RES",
    program_class = factor(rep("RNA_RNP_processing", 4), levels = helpers$program_levels),
    ID = c("GO:1", "GO:2", "GO:3", "GO:4"),
    Description = c("term one", "term two", "term three", "term four"),
    NES = c(1.2, 1.0, -3.0, 8.0), p.adjust = c(0.04, 0.01, 0.01, 0.20),
    core_enrichment = c("P1/P2", "P3", "P4/P5", "P6"),
    core_enrichment_gene = c("G1/G2", "G3", "G4/G5", "G6"),
    Comparison = "CA1sus_CA1res",
    source_term_comparison_file = "terms.csv", source_term_gene_provenance_file = "provenance.csv",
    source_analysis_status_file = "status.csv", source_supplementary_file = "terms.csv",
    source_driver_file = "provenance.csv", source_manifest_file = "manifest.csv",
    evidence_source_family = "canonical_compareGO_ranked_GSEA_GO_BP",
    core_enrichment_contract = "gene_level_claim_allowed == TRUE && core_enrichment_member == TRUE",
    input_contract_version = "fixture_contract",
    stringsAsFactors = FALSE
  )
  summary <- helpers$calculate_program_summary(fixture)

  testthat::expect_identical(summary$representative_term_id, "GO:3")
  testthat::expect_identical(summary$representative_term, "term three")
  testthat::expect_equal(summary$representative_NES, -3)
  testthat::expect_equal(summary$representative_FDR, 0.01)
  testthat::expect_equal(summary$n_sig_terms, 3L)
  testthat::expect_equal(summary$n_positive_sig_terms, 2L)
  testthat::expect_equal(summary$n_negative_sig_terms, 1L)
  testthat::expect_identical(summary$direction_consistency, "mixed_direction")
  testthat::expect_identical(summary$representative_leading_edge_genes, "G4/G5")
})

testthat::test_that("SUS-RES recurrence is contrast-specific and does not gate spatially specific supported rows", {
  helpers <- load_spatial_atlas_helpers()
  fixture <- data.frame(
    dataset = c("neuron_soma", "neuron_soma", "neuron_soma", "microglia"),
    dataset_label = c("Neuron soma", "Neuron soma", "Neuron soma", "Microglia"),
    region = c("CA1", "CA2", "CA2", "CA3"), layer = NA_character_,
    spatial_unit = c("CA1", "CA2", "CA2", "CA3"), compartment = c("soma", "soma", "soma", "microglia"),
    phenotype_contrast = c("SUS_vs_RES", "SUS_vs_RES", "RES_vs_CON", "SUS_vs_CON"),
    program_class = factor(rep("RNA_RNP_processing", 4), levels = helpers$program_levels),
    n_sig_terms = c(1L, 0L, 4L, 5L), representative_NES = c(2, NA, 3, -2),
    stringsAsFactors = FALSE
  )
  out <- helpers$prepare_sus_res_publication_summary(fixture, min_recurrent_units = 2L)
  supported <- out[out$spatial_unit == "CA1", , drop = FALSE]

  testthat::expect_true(all(out$comparison == "SUS_vs_RES"))
  testthat::expect_equal(supported$sus_res_recurrent_units, 1L)
  testthat::expect_equal(supported$sus_res_recurrent_datasets, 1L)
  testthat::expect_false(supported$sus_res_is_recurrent)
  testthat::expect_identical(supported$recurrence_annotation, "spatially_specific")
  testthat::expect_true(supported$publication_include)
  testthat::expect_identical(supported$filter_reason, "included_spatially_specific")
})

testthat::test_that("programs without supported terms have no representative effect", {
  helpers <- load_spatial_atlas_helpers()
  fixture <- data.frame(
    dataset = "microglia", dataset_label = "Microglia", region = "CA1", layer = NA_character_,
    comparison = "CA1microgliasus_CA1microgliares", route_category = "phenotype_within_unit", route_unit = "CA1_microglia",
    spatial_unit = "CA1", compartment = "microglia", phenotype_contrast = "SUS_vs_RES",
    program_class = factor("RNA_RNP_processing", levels = helpers$program_levels),
    ID = "GO:1", Description = "RNA processing", NES = 2, p.adjust = 0.2,
    core_enrichment = "P1", core_enrichment_gene = "G1", Comparison = "CA1sus_CA1res",
    source_term_comparison_file = "terms.csv", source_term_gene_provenance_file = "provenance.csv",
    source_analysis_status_file = "status.csv", source_supplementary_file = "terms.csv",
    source_driver_file = "provenance.csv", source_manifest_file = "manifest.csv",
    evidence_source_family = "canonical_compareGO_ranked_GSEA_GO_BP",
    core_enrichment_contract = "gene_level_claim_allowed == TRUE && core_enrichment_member == TRUE",
    input_contract_version = "fixture_contract",
    stringsAsFactors = FALSE
  )
  summary <- helpers$calculate_program_summary(fixture)
  testthat::expect_equal(summary$n_sig_terms, 0L)
  testthat::expect_true(is.na(summary$representative_term_id))
  testthat::expect_true(is.na(summary$representative_NES))
  testthat::expect_true(is.na(summary$representative_FDR))
  testthat::expect_identical(summary$direction_consistency, "no_FDR_supported_term")
})

testthat::test_that("claim-safe core enrichment is deterministic and fails closed", {
  helpers <- load_spatial_atlas_helpers()
  provenance <- data.frame(
    dataset = "microglia", comparison = "CA1microgliasus_CA1microgliares",
    result_type = c("GSEA_GO", "GSEA_GO", "GSEA_GO", "GSEA_GO", "GSEA_KEGG"),
    ontology = c("BP", "BP", "BP", "BP", "KEGG"), term_id = c("GO:1", "GO:1", "GO:1", "GO:1", "path:1"),
    official_gene_symbol = c("GeneB", "GeneA", "Unsafe", "NonCore", "WrongFamily"),
    ProteinGroupID = c("P2", "P1", "PX", "PN", "PK"),
    gene_level_claim_allowed = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    core_enrichment_member = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  observed <- helpers$summarise_claim_safe_core(provenance)
  reversed <- helpers$summarise_claim_safe_core(provenance[nrow(provenance):1, , drop = FALSE])

  testthat::expect_equal(observed$summary, reversed$summary)
  testthat::expect_identical(observed$summary$core_enrichment_gene, "GeneA;GeneB")
  testthat::expect_identical(observed$summary$core_enrichment, "P1;P2")
  testthat::expect_equal(observed$n_rows, 2L)
  testthat::expect_equal(observed$n_unique_genes, 2L)
})

testthat::test_that("driver recurrence uses claim-safe official genes rather than protein aliases", {
  helpers <- load_spatial_atlas_helpers()
  fixture <- data.frame(
    dataset = "microglia", dataset_label = "Microglia", spatial_unit = "CA1",
    program_class = factor("RNA_RNP_processing", levels = helpers$program_levels),
    Description = "RNA processing", p.adjust = 0.01,
    core_enrichment = "ProteinGroup1;ProteinGroup2",
    core_enrichment_gene = "GeneA;GeneB",
    stringsAsFactors = FALSE
  )
  observed <- helpers$make_driver_recurrence(fixture)
  testthat::expect_setequal(observed$Gene, c("GeneA", "GeneB"))
  testthat::expect_false(any(grepl("ProteinGroup", observed$Gene, fixed = TRUE)))
})

testthat::test_that("canonical execution has no legacy workbook dependency", {
  script <- readLines(
    testthat::test_path("..", "..", "04_differential_expression_enrichment", "07_compareGO_spatial_program_atlas.r"),
    warn = FALSE
  )
  text <- paste(script, collapse = "\n")
  testthat::expect_match(text, "canonical_comparego_manifest_path")
  testthat::expect_match(text, "read_canonical_comparego_bundle")
  testthat::expect_match(text, "gene_level_claim_allowed")
  testthat::expect_match(text, "core_enrichment_member")
  testthat::expect_match(text, ".data$dataset == .env$dataset", fixed = TRUE)
  testthat::expect_false(grepl("Supplementary_Data[.]xlsx", text))
  testthat::expect_false(grepl("07_Top_Genes_Driving_TopTerms[.]xlsx", text))
  testthat::expect_false(grepl("readxl::", text, fixed = TRUE))
})

testthat::test_that("XLSX display truncation leaves analytical CSV values untouched", {
  helpers <- load_spatial_atlas_helpers()
  original <- tibble::tibble(id = c("short", "long"), value = c("abc", paste(rep("x", 40000), collapse = "")))
  workbook <- helpers$xlsx_safe_table(original)
  testthat::expect_identical(original$value[[1]], workbook$value[[1]])
  testthat::expect_gt(nchar(original$value[[2]]), 32767L)
  testthat::expect_lte(nchar(workbook$value[[2]]), 32767L)
  testthat::expect_match(workbook$value[[2]], "truncated for XLSX cell limit", fixed = TRUE)
  testthat::expect_gt(nchar(original$value[[2]]), nchar(workbook$value[[2]]))
})
