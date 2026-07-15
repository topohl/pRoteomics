testthat::test_that("canonical enrichment transformation preserves one group observation", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  d <- data.frame(ProteinGroupID = c("PG1", "PG2", "PG3", "PG4", "PG5"), original_identifier = c("A", "B", "C", "D", "E"),
    member_accessions = c("P1", "P2;P3", "P4;P5", "P6", "P7;H1"), member_gene_symbols = c("A", "B", "C;D", "", "E"),
    representative_accession = "", representative_gene_symbol = "", representative_selection_rule = "display_only",
    protein_group_ambiguity_class = c("single_accession_single_gene", "multi_accession_same_gene", "multi_gene_indistinguishable", "unresolved_group", "mixed_species_or_contaminant"),
    gene_level_claim_allowed = c(TRUE, TRUE, FALSE, FALSE, FALSE), protein_level_claim_allowed = FALSE, mapping_status = "mapped",
    official_gene_symbol = c("A", "B", NA, NA, NA), official_entrez_id = c("1", "2", NA, NA, NA),
    protein_group_gene_annotation_status = c("concordant_official_gene", "concordant_official_gene", "conflicting_member_annotations", "no_mapped_accessions", "incomplete_or_ambiguous_member_annotation"),
    gene_annotation_contract_version = "mouse_gene_annotation_v1", uniprot_mapping_file_hash = "abc", orgdb_package_version = "test",
    t = c(2, 4, 1, 3, -2), stringsAsFactors = FALSE)
  x <- build_enrichment_gene_inputs(d)
  testthat::expect_equal(nrow(x$transformation), nrow(d))
  testthat::expect_equal(x$transformation$GeneSymbol[x$transformation$ProteinGroupID == "PG2"], "B")
  testthat::expect_true(all(x$transformation$eligibility_status[3:5] == "excluded"))
  testthat::expect_equal(x$statistic$column, "t")
})

testthat::test_that("duplicate gene collapse is median, provenance-preserving, and order invariant", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  base <- data.frame(ProteinGroupID = c("PG2", "PG1", "PG3"), GeneSymbol = c("A", "A", "A"), source_statistic = c(-1, 3, 2), eligibility_status = "eligible", original_identifier = c("b", "a", "c"), member_accessions = "P", member_gene_symbols = "A", stringsAsFactors = FALSE)
  a <- collapse_protein_group_genes(base); b <- collapse_protein_group_genes(base[3:1, ])
  testthat::expect_equal(a$collapsed_statistic, 2)
  testthat::expect_true(a$discordant_direction)
  testthat::expect_identical(a, b)
})

testthat::test_that("ORA inputs are directional and never use group identifiers as genes", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  c <- data.frame(GeneSymbol = c("A", "B", "C"), collapsed_statistic = c(2, -3, 0.2), stringsAsFactors = FALSE)
  ora <- build_ora_inputs(c, fc_threshold = 1)
  testthat::expect_equal(ora$up, "A"); testthat::expect_equal(ora$down, "B")
})

testthat::test_that("clusterProfiler cannot reintroduce effect sorting duplicate selection", {
  script <- paste(readLines(testthat::test_path("..", "..", "04_differential_expression_enrichment", "01_clusterProfiler.r"), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("sort\\(na.omit\\(original_gene_list\\).*duplicated\\(names\\(gene_list\\)\\)", script))
})

testthat::test_that("member order does not alter a same-gene enrichment mapping", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  make_row <- function(members) data.frame(ProteinGroupID = "PG1", original_identifier = members, member_accessions = members,
    member_gene_symbols = "A", representative_accession = "", representative_gene_symbol = "", representative_selection_rule = "display_only",
    protein_group_ambiguity_class = "multi_accession_same_gene", gene_level_claim_allowed = TRUE, protein_level_claim_allowed = FALSE,
    mapping_status = "mapped", official_gene_symbol = "A", official_entrez_id = "1",
    protein_group_gene_annotation_status = "concordant_official_gene", gene_annotation_contract_version = "mouse_gene_annotation_v1",
    uniprot_mapping_file_hash = "abc", orgdb_package_version = "test", log2fc = 2, stringsAsFactors = FALSE)
  testthat::expect_identical(build_enrichment_gene_inputs(make_row("P1;P2"))$ranked, build_enrichment_gene_inputs(make_row("P2;P1"))$ranked)
})

testthat::test_that("compareGO prefers manifest-provided collapsed gene inputs", {
  script <- paste(readLines(testthat::test_path("..", "..", "04_differential_expression_enrichment", "02_compareGO.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "comparison_input_file")
  testthat::expect_match(script, "GeneSymbol")
})

testthat::test_that("strict mode rejects legacy gene-only input", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  testthat::expect_error(build_enrichment_gene_inputs(data.frame(gene_symbol = "A", log2fc = 1), strict = TRUE), "Canonical ProteinGroupID")
})

testthat::test_that("ORA significance uses collapsed FDR without effect-dependent duplicate selection", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  d <- data.frame(GeneSymbol = c("A", "B"), collapsed_statistic = c(2, -2), collapsed_fdr = c(.01, .2), collapsed_logfc = c(2, -2))
  x <- build_ora_inputs(d, fdr_threshold = .05, fc_threshold = 1)
  testthat::expect_equal(x$up, "A"); testthat::expect_length(x$down, 0)
})

testthat::test_that("ORA input audits preserve direction for empty and non-empty inputs", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  empty <- make_ora_input_audit(character(), "upregulated")
  testthat::expect_equal(names(empty), c("GeneSymbol", "ORA_direction"))
  testthat::expect_equal(nrow(empty), 0L)

  one <- make_ora_input_audit("A", "downregulated")
  testthat::expect_equal(one$GeneSymbol, "A")
  testthat::expect_equal(one$ORA_direction, "downregulated")

  multiple <- make_ora_input_audit(c("A", "B"), "universe")
  testthat::expect_equal(multiple$GeneSymbol, c("A", "B"))
  testthat::expect_equal(multiple$ORA_direction, c("universe", "universe"))
})

testthat::test_that("GSEA diagnostics and result guards handle null, malformed, and empty results", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  ranks <- c(A = 2, B = -1, C = .5)
  diagnostics <- gsea_input_diagnostics(ranks, c("A", "C", "X"))
  testthat::expect_equal(diagnostics$orgdb_symbol_matches, 2)
  testthat::expect_equal(diagnostics$finite_rank_statistics, 3)
  testthat::expect_equal(diagnostics$unique_rank_statistic_values, 3)
  testthat::expect_error(gsea_result_table(NULL, diagnostics), "GO GSEA returned NULL: ranked_genes=3, orgdb_symbol_matches=2")
  testthat::expect_error(gsea_result_table(structure(list(), class = "not_gsea"), diagnostics), "invalid result class")
  empty <- structure(list(result = data.frame()), class = "gseaResult")
  non_empty <- structure(list(result = data.frame(ID = "GO:1")), class = "gseaResult")
  testthat::expect_equal(nrow(gsea_result_table(empty, diagnostics)), 0L)
  testthat::expect_equal(nrow(gsea_result_table(non_empty, diagnostics)), 1L)
})

testthat::test_that("canonical SYMBOL resolution is deterministic and provenance-preserving", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  collapsed <- data.frame(
    GeneSymbol = c("Aaas", "bRaF", "DUP", "missing"),
    collapsed_statistic = c(1, 2, 3, 4),
    contributing_ProteinGroupIDs = c("PG1", "PG2", "PG3", "PG4"),
    stringsAsFactors = FALSE
  )
  resolved <- resolve_enrichment_symbols(collapsed, c("Aaas", "Braf", "Dup", "DUP"))
  audit <- resolved$audit
  testthat::expect_equal(audit$symbol_resolution_status,
    c("exact_symbol_match", "case_normalized_symbol_match", "exact_symbol_match", "unmatched_symbol"))
  testthat::expect_equal(audit$resolved_GeneSymbol, c("Aaas", "Braf", "DUP", NA_character_))
  testthat::expect_equal(audit$contributing_ProteinGroupIDs, collapsed$contributing_ProteinGroupIDs)
  testthat::expect_equal(names(resolved$ranked), c("DUP", "Braf", "Aaas"))

  ambiguous <- resolve_enrichment_symbols(
    data.frame(GeneSymbol = "dup", collapsed_statistic = 1, contributing_ProteinGroupIDs = "PG5"),
    c("Dup", "DUP"))
  testthat::expect_true(ambiguous$audit$ambiguous_case_match)
  testthat::expect_false(ambiguous$audit$included_in_enrichment)
})

testthat::test_that("SYMBOL-to-ENTREZ KEGG preparation is median and order invariant", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  collapsed <- data.frame(GeneSymbol = c("A", "B", "C"), collapsed_statistic = c(3, -1, 5),
    contributing_ProteinGroupIDs = c("PG1", "PG2", "PG3"), stringsAsFactors = FALSE)
  mapping <- data.frame(SYMBOL = c("A", "B", "C", "C"), ENTREZID = c("1", "1", "2", "3"), stringsAsFactors = FALSE)
  a <- prepare_kegg_symbol_ranks(collapsed, mapping)
  b <- prepare_kegg_symbol_ranks(collapsed[3:1, ], mapping[4:1, ])
  testthat::expect_identical(a$ranked, b$ranked)
  testthat::expect_equal(unname(a$ranked), 1)
  testthat::expect_equal(names(a$ranked), "1")
  testthat::expect_equal(a$collapse$gene_collapse_rule, "median_finite_statistics_for_duplicate_entrez")
  testthat::expect_equal(a$audit$entrez_mapping_status[a$audit$GeneSymbol == "C"], "ambiguous_entrez_match")
})

testthat::test_that("canonical symbol vectors are never submitted as UniProt identifiers", {
  script <- paste(readLines(testthat::test_path("..", "..", "04_differential_expression_enrichment", "01_clusterProfiler.r"), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("bitr\\(names\\((original_)?gene_list\\).*fromType\\s*=\\s*['\"]UNIPROT", script, perl = TRUE))
  testthat::expect_false(grepl("enrichGO\\([\\s\\S]{0,500}keyType\\s*=\\s*['\"]UNIPROT", script, perl = TRUE))
  testthat::expect_false(grepl("merge\\(df,\\s*[^,]+,\\s*by.x\\s*=\\s*['\"]gene_symbol", script, perl = TRUE))
})

testthat::test_that("generic enrichment result guard accepts empty results and rejects NULL", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  testthat::expect_error(enrichment_result_table(NULL, "KEGG GSEA"), "KEGG GSEA returned NULL")
  empty <- structure(list(result = data.frame()), class = "enrichResult")
  testthat::expect_equal(nrow(enrichment_result_table(empty, "general GO enrichment")), 0L)
})

testthat::test_that("aliases converge before duplicate-gene median collapse", {
  source(testthat::test_path("..", "..", "R", "protein_mapping_utils.R"))
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  maps <- list(
    uniprot_map = data.frame(UNIPROT = character(), SYMBOL = character(), ENTREZID = character()),
    symbol_map = data.frame(SYMBOL = "GeneA", ENTREZID = "1"),
    alias_map = data.frame(ALIAS = c("OLD1", "OLD2"), SYMBOL = "GeneA", ENTREZID = "1"),
    orgdb_package_version = "test", annotation_contract_version = "mouse_gene_annotation_v1"
  )
  a <- resolve_mouse_gene_annotation("", "OLD1", maps)
  b <- resolve_mouse_gene_annotation("", "OLD2", maps)
  testthat::expect_equal(c(a$official_gene_symbol, b$official_gene_symbol), c("GeneA", "GeneA"))
  d <- data.frame(ProteinGroupID = c("PG1", "PG2"), original_identifier = c("OLD1", "OLD2"),
    member_accessions = c("P1", "P2"), member_gene_symbols = c("OLD1", "OLD2"),
    representative_accession = c("P1", "P2"), representative_gene_symbol = c("GeneA", "GeneA"),
    representative_selection_rule = "display_only", protein_group_ambiguity_class = "single_accession_single_gene",
    gene_level_claim_allowed = TRUE, protein_level_claim_allowed = TRUE, mapping_status = "mapped",
    official_gene_symbol = c(a$official_gene_symbol, b$official_gene_symbol), official_entrez_id = "1",
    protein_group_gene_annotation_status = "concordant_official_gene", gene_annotation_contract_version = "mouse_gene_annotation_v1",
    uniprot_mapping_file_hash = "abc", orgdb_package_version = "test", t = c(1, 3), stringsAsFactors = FALSE)
  collapsed <- build_enrichment_gene_inputs(d)$collapse
  testthat::expect_equal(nrow(collapsed), 1L)
  testthat::expect_equal(collapsed$collapsed_statistic, 2)
  testthat::expect_equal(collapsed$contributing_ProteinGroupIDs, "PG1;PG2")
})

testthat::test_that("strict enrichment validates precomputed symbols without repairing them", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  d <- data.frame(GeneSymbol = c("GeneA", "genea"), contributing_ProteinGroupIDs = c("PG1", "PG2"))
  audit <- validate_precomputed_enrichment_symbols(d, "GeneA")
  testthat::expect_equal(audit$symbol_resolution_status,
    c("precomputed_exact_symbol_validated", "invalid_precomputed_official_symbol"))
  testthat::expect_false(audit$case_normalized[[2]])
  testthat::expect_false(audit$included_in_enrichment[[2]])
})
