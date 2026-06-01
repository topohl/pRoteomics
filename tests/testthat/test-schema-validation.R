testthat::test_that("schema validation catches missing columns and accepts valid claims", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "schema_validation.R"))
  testthat::skip_if_not_installed("yaml")

  good <- data.frame(
    claim_id = "CLAIM_0001",
    dataset = "microglia",
    region = NA_character_,
    layer_cell_compartment = NA_character_,
    contrast = "SUS_vs_CON",
    biological_program = "microglia_signature",
    direction = "positive_NES",
    key_proteins_genes = "AIF1",
    evidence_type = "microglia_signature_enrichment",
    effect_size_NES = 1.2,
    raw_p = 0.01,
    FDR = 0.04,
    robustness_stability_metric = "n=1",
    source_file = "example.csv",
    figure_table_target = "supplement",
    interpretation_note = "example",
    claim_grade = "C",
    primary_evidence = "Enrichment/program summary",
    orthogonal_support = "Review companion outputs",
    major_limitation = "ROI, not purified microglia",
    safe_interpretation = "ROI-associated signal",
    unsafe_overinterpretation = "Purified microglia claim",
    stringsAsFactors = FALSE
  )
  testthat::expect_silent(validate_table_schema(good, "biological_claims_table", strict = TRUE))
  bad <- good[, setdiff(names(good), "claim_grade")]
  testthat::expect_error(validate_table_schema(bad, "biological_claims_table", strict = TRUE), "Missing required")

  bad_grade <- good
  bad_grade$claim_grade <- "Z"
  testthat::expect_error(validate_table_schema(bad_grade, "biological_claims_table", strict = TRUE), "claim_grade has invalid")

  bad_p <- good
  bad_p$raw_p <- 1.2
  testthat::expect_error(validate_table_schema(bad_p, "biological_claims_table", strict = TRUE), "raw_p outside allowed range")

  bad_fdr <- good
  bad_fdr$FDR <- -0.01
  testthat::expect_error(validate_table_schema(bad_fdr, "biological_claims_table", strict = TRUE), "FDR outside allowed range")

  bad_dataset <- good
  bad_dataset$dataset <- "microglia_layer"
  testthat::expect_error(validate_table_schema(bad_dataset, "biological_claims_table", strict = TRUE), "dataset has invalid")
})

testthat::test_that("mapped contrast schema validates p-value ranges", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "schema_validation.R"))
  testthat::skip_if_not_installed("yaml")

  good <- data.frame(
    gene_id = "Aif1",
    log2FoldChange = 0.5,
    pvalue = 0.2,
    padj = 0.8,
    dataset = "microglia",
    stringsAsFactors = FALSE
  )
  testthat::expect_silent(validate_table_schema(good, "mapped_contrast", strict = FALSE))

  bad <- good
  bad$padj <- 1.01
  testthat::expect_error(validate_table_schema(bad, "mapped_contrast", strict = FALSE), "padj outside allowed range")
})
