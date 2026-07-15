testthat::test_that("canonical enrichment transformation preserves one group observation", {
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R"))
  d <- data.frame(ProteinGroupID = c("PG1", "PG2", "PG3", "PG4", "PG5"), original_identifier = c("A", "B", "C", "D", "E"),
    member_accessions = c("P1", "P2;P3", "P4;P5", "P6", "P7;H1"), member_gene_symbols = c("A", "B", "C;D", "", "E"),
    representative_accession = "", representative_gene_symbol = "", representative_selection_rule = "display_only",
    protein_group_ambiguity_class = c("single_accession_single_gene", "multi_accession_same_gene", "multi_gene_indistinguishable", "unresolved_group", "mixed_species_or_contaminant"),
    gene_level_claim_allowed = c(TRUE, TRUE, FALSE, FALSE, FALSE), protein_level_claim_allowed = FALSE, mapping_status = "mapped", t = c(2, 4, 1, 3, -2), stringsAsFactors = FALSE)
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
    mapping_status = "mapped", log2fc = 2, stringsAsFactors = FALSE)
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
