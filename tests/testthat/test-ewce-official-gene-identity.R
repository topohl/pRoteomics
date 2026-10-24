source(testthat::test_path("..", "..", "R", "paths.R"))

# EWCE takes gene identity from the canonical protein-group mapping (official mouse symbols), as the
# differential-abundance GSEA/ORA does. The DIA-NN Genes token must not be re-mapped by alias.

ewce_script_lines <- function() {
  readLines(repo_path("analysis/enrichment", "run_ewce_celltype_enrichment.R"), warn = FALSE)
}

testthat::test_that("EWCE maps genes through the canonical loader, not by alias", {
  code <- paste(sub("^\\s*#.*$", "", ewce_script_lines()), collapse = "\n")
  testthat::expect_false(grepl("ALIAS", code, fixed = TRUE))
  testthat::expect_false(grepl("mapIds(", code, fixed = TRUE))
  testthat::expect_true(grepl("qc_load_canonical_expression(", code, fixed = TRUE))
  testthat::expect_true(grepl("Gene = dplyr::if_else(gene_level_official, official_gene_symbol, NA_character_)",
                              code, fixed = TRUE))
})

testthat::test_that("only gene-level groups with one concordant official gene get a gene", {
  lines <- ewce_script_lines()
  start <- grep("^gene_level_official <- ", lines)
  testthat::expect_length(start, 1L)
  end <- start
  while (grepl("&\\s*$", lines[end])) end <- end + 1L
  canonical_features <- data.frame(
    enrichment_eligible = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
    protein_group_ambiguity_class = c("single_accession_single_gene", "multi_gene_indistinguishable",
      "single_accession_single_gene", "single_accession_single_gene", "multi_accession_same_gene",
      "single_accession_single_gene"),
    protein_group_gene_annotation_status = c("concordant_official_gene", "conflicting_member_annotations",
      "concordant_official_gene", "incomplete_or_ambiguous_member_annotation", "concordant_official_gene",
      "concordant_official_gene"),
    official_gene_symbol = c("Syp", NA, "Lmo7", NA, "Ckb", ""),
    stringsAsFactors = FALSE
  )
  eval(parse(text = lines[start:end]))
  testthat::expect_identical(gene_level_official, c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))
})

testthat::test_that("pipeline.yml declares the inputs EWCE gene identity reads", {
  source(repo_path("R", "pipeline_registry.R"))
  steps <- pipeline_steps(read_pipeline_registry(), "enrichment", dataset = "neuron_neuropil")
  step <- steps[steps$script == "analysis/enrichment/run_ewce_celltype_enrichment.R", , drop = FALSE]
  testthat::expect_equal(nrow(step), 1L)
  required <- split_registry_paths(step$consumes_required)
  testthat::expect_true(all(c(
    "data/processed/01_preprocessing/impute/*_pgmatrix_imputed_<dataset>_*_missing70pct.xlsx",
    "data/metadata/TPE9_sample_metadata_males.xlsx",
    "data/external/MOUSE_10090_idmapping.dat"
  ) %in% required))
  testthat::expect_true(all(c("data/metadata/manual_mapping.xlsx", "data/metadata/manual_gene_annotation_overrides.csv") %in%
    split_registry_paths(step$consumes_optional)))
})
