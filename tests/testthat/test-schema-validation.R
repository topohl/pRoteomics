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
})
