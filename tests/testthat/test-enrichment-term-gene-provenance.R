testthat::test_that("term core genes expand to every eligible ProteinGroupID deterministically", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  gsea <- data.frame(
    ID = "GO:0001", Description = "example", NES = 1.5, p.adjust = 0.01,
    setSize = 2L, core_enrichment = "GeneB/GeneA", stringsAsFactors = FALSE
  )
  collapsed <- data.frame(
    GeneSymbol = c("GeneA", "GeneB"), EntrezID = c("001", "002"),
    collapsed_statistic = c(2, 1), contributing_ProteinGroupIDs = c("PG2;PG1", "PG3"),
    stringsAsFactors = FALSE
  )
  transform <- data.frame(
    ProteinGroupID = c("PG3", "PG2", "PG1", "PGX"),
    GeneSymbol = c("GeneB", "GeneA", "GeneA", NA), EntrezID = c("002", "001", "001", NA),
    member_accessions = c("P3", "P2", "P1", "PX;PY"),
    protein_group_gene_annotation_status = c(rep("concordant_official_gene", 3), "conflicting_official_genes"),
    gene_level_claim_allowed = c(TRUE, TRUE, TRUE, FALSE),
    gene_annotation_contract_version = "mouse_gene_annotation_v1",
    eligibility_status = c("eligible", "eligible", "eligible", "excluded"),
    stringsAsFactors = FALSE
  )
  a <- build_enrichment_term_gene_provenance(
    gsea, collapsed, transform, "microglia", "comparison_from_manifest", "GSEA_GO", "BP"
  )
  b <- build_enrichment_term_gene_provenance(
    gsea, collapsed[2:1, ], transform[4:1, ], "microglia", "comparison_from_manifest", "GSEA_GO", "BP"
  )
  testthat::expect_identical(a, b)
  testthat::expect_equal(a$ProteinGroupID[a$official_gene_symbol == "GeneA"], c("PG1", "PG2"))
  testthat::expect_false("PGX" %in% a$ProteinGroupID)
  testthat::expect_true(all(a$core_enrichment_member))
  testthat::expect_type(a$official_entrez_id, "character")
})

testthat::test_that("term provenance handles zero terms and rejects missing core provenance", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  empty_gsea <- data.frame(
    ID = character(), Description = character(), NES = numeric(), p.adjust = numeric(),
    setSize = integer(), core_enrichment = character(), stringsAsFactors = FALSE
  )
  collapsed <- data.frame(
    GeneSymbol = "GeneA", EntrezID = "001", collapsed_statistic = 2,
    contributing_ProteinGroupIDs = "PG1", stringsAsFactors = FALSE
  )
  transform <- data.frame(
    ProteinGroupID = "PG1", GeneSymbol = "GeneA", EntrezID = "001", member_accessions = "P1",
    protein_group_gene_annotation_status = "concordant_official_gene", gene_level_claim_allowed = TRUE,
    gene_annotation_contract_version = "mouse_gene_annotation_v1", eligibility_status = "eligible",
    stringsAsFactors = FALSE
  )
  empty <- build_enrichment_term_gene_provenance(
    empty_gsea, collapsed, transform, "microglia", "cmp", "GSEA_GO", "BP"
  )
  testthat::expect_equal(nrow(empty), 0L)
  testthat::expect_identical(names(empty), names(empty_enrichment_term_gene_provenance()))

  bad <- data.frame(
    ID = "GO:1", Description = "missing", NES = 1, p.adjust = 0.1,
    setSize = 1L, core_enrichment = "GeneMissing", stringsAsFactors = FALSE
  )
  testthat::expect_error(
    build_enrichment_term_gene_provenance(bad, collapsed, transform, "microglia", "cmp", "GSEA_GO", "BP"),
    "absent from collapsed official-gene provenance"
  )
})

testthat::test_that("ineligible groups cannot re-enter GO or KEGG term provenance", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  gsea_go <- data.frame(
    ID = "GO:1", Description = "ambiguous", NES = 1, p.adjust = 0.1,
    setSize = 1L, core_enrichment = "GeneX", stringsAsFactors = FALSE
  )
  collapsed <- data.frame(
    GeneSymbol = "GeneX", EntrezID = "009", collapsed_statistic = 1,
    contributing_ProteinGroupIDs = "PGX", stringsAsFactors = FALSE
  )
  transform <- data.frame(
    ProteinGroupID = "PGX", GeneSymbol = "GeneX", EntrezID = "009", member_accessions = "PX;PY",
    protein_group_gene_annotation_status = "conflicting_official_genes", gene_level_claim_allowed = FALSE,
    gene_annotation_contract_version = "mouse_gene_annotation_v1", eligibility_status = "excluded",
    stringsAsFactors = FALSE
  )
  testthat::expect_error(
    build_enrichment_term_gene_provenance(gsea_go, collapsed, transform, "microglia", "cmp", "GSEA_GO", "BP"),
    "ineligible or missing ProteinGroupID"
  )

  transform$protein_group_gene_annotation_status <- "concordant_official_gene"
  transform$gene_level_claim_allowed <- TRUE
  transform$eligibility_status <- "eligible"
  gsea_kegg <- gsea_go
  gsea_kegg$ID <- "mmu0001"
  gsea_kegg$core_enrichment <- "009"
  kegg <- build_enrichment_term_gene_provenance(
    gsea_kegg, collapsed, transform, "microglia", "cmp", "GSEA_KEGG", "KEGG",
    core_identifier = "ENTREZID"
  )
  testthat::expect_equal(kegg$official_gene_symbol, "GeneX")
  testthat::expect_equal(kegg$official_entrez_id, "009")
})
