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
    claim_type = "microglia_signature",
    claim_use_class = "blocked",
    raw_top_GO_term = NA_character_,
    representative_GO_terms = NA_character_,
    semantic_parent_label = NA_character_,
    safe_program_label = "microglia_signature",
    term_label_risk = "not_applicable",
    label_confidence = "not_applicable",
    label_basis = "not_applicable",
    label_downgrade_reason = "not_applicable",
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
    claim_allowed = FALSE,
    claim_gate_status = "disallowed",
    claim_downgrade_reason = "animal_level_gate=fail; marker_contamination_gate=missing_required",
    model_fit_status = "pass",
    statistical_evidence_status = "pass",
    claim_gate_model_status = "pass",
    primary_model_status = "pass",
    animal_level_gate = "fail",
    qc_gate = "missing_optional",
    missingness_gate = "pass",
    batch_confound_gate = "pass",
    marker_contamination_gate = "missing_required",
    microglia_roi_gate = "pass",
    neuropil_independence_gate = "not_applicable",
    robustness_gate = "pass",
    evidence_independence_gate = "fail",
    missingness_confounded = "not_detected",
    plate_or_batch_confounded = "not_detected",
    region_layer_imbalance_risk = "not_applicable",
    animal_pseudoreplication_risk = "possible",
    early_pc_association = "not_available",
    marker_contamination_risk = "not_available",
    qc_interpretation_flag = "not_available",
    animal_level_status = "sample_level_or_unclear",
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

testthat::test_that("WGCNA group-effect output validation checks required columns and ranges", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "validation_utils.R"))

  good <- data.frame(
    dataset = "microglia",
    level = "module",
    endpoint_id = "M1",
    endpoint_label = "M1",
    contrast = "SUS - CON",
    estimate = 0.2,
    SE = 0.1,
    p_value = 0.03,
    FDR_within_dataset_level = 0.05,
    FDR_global = 0.08,
    evidence_status = "nominal_only",
    n_samples = 8L,
    n_animals = 6L,
    n_animals_total = 6L,
    n_animals_per_group = "CON=3;SUS=3",
    min_animals_per_group = 3L,
    n_samples_total = 8L,
    n_samples_per_group = "CON=4;SUS=4",
    animal_level_status = "animal_level",
    pseudoreplication_guard = "pass",
    model_type = "lm",
    model_family = "linear_model",
    model_formula = "eigengene ~ StressGroup",
    primary_model_stable = TRUE,
    claim_allowed_model = TRUE,
    model_downgrade_reason = "none",
    fallback_used = FALSE,
    fallback_type = NA_character_,
    formula_used = "eigengene ~ StressGroup",
    rank_deficient_model = FALSE,
    singular_model = FALSE,
    emmeans_success = TRUE,
    animal_random_effect_used = FALSE,
    biological_replicate_unit = "animal",
    model_warning = "",
    stringsAsFactors = FALSE
  )
  path <- tempfile("module_group_effects-", fileext = ".csv")
  path <- file.path(dirname(path), "module_group_effects.csv")
  utils::write.csv(good, path, row.names = FALSE)
  ok <- validate_known_pipeline_output(path, dataset = "microglia")
  testthat::expect_equal(ok$validation_status, "ok")

  bad <- good
  bad$p_value <- 1.2
  bad$n_animals <- -1L
  utils::write.csv(bad, path, row.names = FALSE)
  warn <- validate_known_pipeline_output(path, dataset = "microglia")
  testthat::expect_equal(warn$validation_status, "warning")
  testthat::expect_match(warn$validation_message, "p_value")
  testthat::expect_match(warn$validation_message, "n_animals")

  missing_col <- good[, setdiff(names(good), "formula_used")]
  utils::write.csv(missing_col, path, row.names = FALSE)
  warn2 <- validate_known_pipeline_output(path, dataset = "microglia")
  testthat::expect_equal(warn2$validation_status, "warning")
  testthat::expect_match(warn2$validation_message, "formula_used")
})

testthat::test_that("microglia neuropil independence audit schemas validate", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "schema_validation.R"))
  testthat::skip_if_not_installed("yaml")

  claim_gate <- data.frame(
    module_or_supermodule_id = "WGCNA_#000000",
    contrast = "SUS - CON",
    biological_program = "microglia-enriched ROI program",
    microenvironment_class = "microglia_supported",
    adjustment_mode = "predeclared_primary",
    covariate_family = "global_neuropil_score",
    primary_effect_status = "FDR_pass",
    primary_effect_claim_relevant = TRUE,
    primary_effect_threshold = "FDR<=0.05; nominal_p<=0.05 diagnostic_only",
    independence_classification = "neuropil_independent",
    claim_gate_eligible = TRUE,
    downgrade_reason = "none",
    n_matched_animals = 6L,
    min_animals_per_group = 2L,
    effect_before = 0.4,
    effect_after = 0.35,
    effect_before_abs = 0.4,
    effect_before_near_zero = FALSE,
    percent_attenuation = 12.5,
    percent_attenuation_reliable = TRUE,
    direction_preserved = TRUE,
    audit_group_n = 1L,
    eligible_without_primary_effect_count = 0L,
    stringsAsFactors = FALSE
  )
  testthat::expect_silent(validate_table_schema(claim_gate, "microglia_neuropil_independence_claim_gate", strict = TRUE))

  selection <- data.frame(
    endpoint_id = "WGCNA_#000000",
    module_or_supermodule_id = "WGCNA_#000000",
    candidate_covariates = "z_reference_cortical_excitatory_neuron_score",
    predeclared_covariates_available = "global_neuropil_score=z_reference_cortical_excitatory_neuron_score",
    selected_primary_covariate = "z_reference_cortical_excitatory_neuron_score",
    selected_secondary_covariate = "",
    selected_exploratory_covariate = "WGCNA_#111111",
    selection_rule = "predeclared=fixed family priority; exploratory=strongest absolute Spearman",
    primary_claim_gate_eligible = TRUE,
    secondary_claim_gate_eligible = FALSE,
    exploratory_claim_gate_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  testthat::expect_silent(validate_table_schema(selection, "microglia_neuropil_covariate_selection_audit", strict = FALSE))
})
