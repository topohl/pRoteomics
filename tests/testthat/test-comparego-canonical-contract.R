make_test_provenance <- function() {
  data.frame(
    dataset = "microglia", comparison = "manifest_identity", result_type = "GSEA_GO",
    ontology = "BP", term_id = "GO:1", term_description = "term",
    official_gene_symbol = "GeneA", official_entrez_id = "001", ProteinGroupID = "PG1",
    member_accessions = "P1", protein_group_gene_annotation_status = "concordant_official_gene",
    gene_level_claim_allowed = TRUE, rank_statistic = 2, core_enrichment_member = TRUE,
    enrichment_contract_version = "protein_group_enrichment_v3_term_gene_provenance",
    gene_annotation_contract_version = "mouse_gene_annotation_v1", stringsAsFactors = FALSE
  )
}

make_test_cluster_manifest <- function(root, status = "success_with_terms", malformed = FALSE,
                                       missing_provenance = FALSE, weird_filename = FALSE) {
  terms <- if (malformed) data.frame(ID = "GO:1", Description = "term") else data.frame(
    ID = "GO:1", Description = "term", NES = 1.2, p.adjust = 0.02,
    setSize = 1L, core_enrichment = "GeneA", stringsAsFactors = FALSE
  )
  if (status == "success_zero_terms") {
    terms <- data.frame(ID = character(), Description = character(), NES = numeric(),
      p.adjust = numeric(), setSize = integer(), core_enrichment = character(), stringsAsFactors = FALSE)
  }
  output <- file.path(root, if (weird_filename) "does-not-encode-the-comparison.csv" else "terms.csv")
  provenance <- file.path(root, "term_provenance.csv")
  collapsed <- file.path(root, "collapsed.csv")
  collapsed_provenance <- file.path(root, "collapsed_provenance.csv")
  if (status != "failed") {
    utils::write.csv(terms, output, row.names = FALSE)
    utils::write.csv(if (status == "success_zero_terms") make_test_provenance()[0, ] else make_test_provenance(), provenance, row.names = FALSE)
    utils::write.csv(data.frame(official_gene_symbol = "GeneA", official_entrez_id = "001"), collapsed, row.names = FALSE)
    utils::write.csv(data.frame(official_gene_symbol = "GeneA", ProteinGroupID = "PG1"), collapsed_provenance, row.names = FALSE)
  }
  if (missing_provenance && file.exists(provenance)) unlink(provenance)
  data.frame(
    dataset = "microglia", comparison = "manifest_identity", result_type = "GSEA_GO", ontology = "BP",
    analysis_status = status, n_terms = if (status == "success_with_terms") 1L else if (status == "success_zero_terms") 0L else NA_integer_,
    output_table = if (status == "failed") NA_character_ else output,
    collapsed_gene_input_file = if (status == "failed") NA_character_ else collapsed,
    collapsed_gene_provenance_file = if (status == "failed") NA_character_ else collapsed_provenance,
    term_gene_provenance_file = if (status == "failed") NA_character_ else provenance,
    enrichment_contract_version = "clusterProfiler_manifest_v3_term_gene_provenance",
    gene_annotation_contract_version = if (status == "failed") NA_character_ else "mouse_gene_annotation_v1",
    stringsAsFactors = FALSE
  )
}

testthat::test_that("canonical compareGO accepts terms, zero terms, and explicit failures", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
  root <- tempfile("comparego-contract-")
  dir.create(root)
  with_terms <- collect_canonical_comparego_outputs(make_test_cluster_manifest(root), strict = TRUE)
  testthat::expect_equal(with_terms$terms$comparison, "manifest_identity")
  testthat::expect_equal(with_terms$terms$term_id, "GO:1")
  testthat::expect_equal(with_terms$provenance$ProteinGroupID, "PG1")
  testthat::expect_type(with_terms$provenance$official_entrez_id, "character")
  testthat::expect_equal(with_terms$provenance$official_entrez_id, "001")

  zero_root <- tempfile("comparego-zero-")
  dir.create(zero_root)
  zero <- collect_canonical_comparego_outputs(make_test_cluster_manifest(zero_root, "success_zero_terms"), strict = TRUE)
  testthat::expect_equal(nrow(zero$terms), 0L)
  testthat::expect_equal(nrow(zero$provenance), 0L)
  testthat::expect_equal(zero$status$comparego_action, "completed_zero_terms")

  failed_manifest <- make_test_cluster_manifest(root, "failed")
  failed <- collect_canonical_comparego_outputs(failed_manifest, strict = TRUE)
  testthat::expect_equal(failed$status$comparego_action, "recorded_failed")
  testthat::expect_equal(nrow(failed$terms), 0L)

  kegg <- failed_manifest
  kegg$result_type <- "GSEA_KEGG"
  kegg$ontology <- "KEGG"
  mixed <- collect_canonical_comparego_outputs(rbind(make_test_cluster_manifest(root), kegg), strict = TRUE)
  testthat::expect_equal(
    mixed$status$comparego_action[mixed$status$result_type == "GSEA_KEGG"],
    "skipped_unsupported_result_type"
  )
})

testthat::test_that("canonical compareGO rejects stale, malformed, and missing provenance inputs", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
  root <- tempfile("comparego-invalid-")
  dir.create(root)
  stale <- make_test_cluster_manifest(root)
  stale$enrichment_contract_version <- "old_contract"
  testthat::expect_error(collect_canonical_comparego_outputs(stale), "Stale clusterProfiler manifest contract")

  malformed_root <- tempfile("comparego-malformed-")
  dir.create(malformed_root)
  malformed <- make_test_cluster_manifest(malformed_root, malformed = TRUE)
  testthat::expect_error(collect_canonical_comparego_outputs(malformed), "missing required columns")

  missing_root <- tempfile("comparego-missing-")
  dir.create(missing_root)
  missing <- make_test_cluster_manifest(missing_root, missing_provenance = TRUE)
  testthat::expect_error(collect_canonical_comparego_outputs(missing), "missing term_gene_provenance_file")
})

testthat::test_that("comparison identity and output are invariant to filenames and manifest order", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(testthat::test_path("..", "..", "R", "enrichment_io.R"))
  root_a <- tempfile("comparego-order-a-"); dir.create(root_a)
  root_b <- tempfile("comparego-order-b-"); dir.create(root_b)
  a <- make_test_cluster_manifest(root_a, weird_filename = TRUE)
  b <- make_test_cluster_manifest(root_b)
  b$comparison <- "second_identity"
  provenance_b <- make_test_provenance(); provenance_b$comparison <- "second_identity"
  utils::write.csv(provenance_b, b$term_gene_provenance_file, row.names = FALSE)
  forward <- collect_canonical_comparego_outputs(rbind(a, b))
  reverse <- collect_canonical_comparego_outputs(rbind(b, a))
  testthat::expect_identical(forward$terms, reverse$terms)
  testthat::expect_identical(forward$provenance, reverse$provenance)
  testthat::expect_true("manifest_identity" %in% forward$terms$comparison)
})

testthat::test_that("canonical compareGO path terminates before legacy UniProt logic", {
  script <- readLines(testthat::test_path("..", "..", "04_differential_expression_enrichment", "02_compareGO.r"), warn = FALSE)
  marker <- grep("LEGACY_COMPAREGO_TAIL_DISABLED_BY_CANONICAL_EXIT", script, fixed = TRUE)
  testthat::expect_length(marker, 1L)
  canonical <- paste(script[seq_len(marker - 1L)], collapse = "\n")
  testthat::expect_match(canonical, "quit\\(status = 0, save = \\\"no\\\"\\)")
  testthat::expect_false(grepl("keyType\\s*=\\s*['\"]UNIPROT|fromType\\s*=\\s*['\"]UNIPROT", canonical, perl = TRUE))
  testthat::expect_false(grepl("PROTEOMICS_COMPAREGO_QUICK_PLOTS", canonical, fixed = TRUE))
})
