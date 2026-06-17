testthat::test_that("WGCNA downstream entrypoints exist", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  scripts <- c(
    "03_qc_exploration/04b_import_reference_marker_sources.r",
    "03_qc_exploration/05_empirical_roi_marker_discovery.r",
    "03_qc_exploration/06_wgcna_marker_trait_export.r",
    "06_modules_WGCNA/05_module_supermodule_group_effects.r",
    "06_modules_WGCNA/06_annotate_module_microenvironment.r",
    "06_modules_WGCNA/07_wgcna_interpretable_summary.r",
    "06_modules_WGCNA/08_microglia_neuropil_independence.R",
    "06_modules_WGCNA/08_module_complex_architecture.r",
    "06_modules_WGCNA/09_module_robustness_sensitivity.r"
  )
  testthat::expect_true(all(file.exists(repo_path(scripts))))
})

testthat::test_that("WGCNA downstream dry-runs report contracts", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x
  cmd <- file.path(R.home("bin"), "Rscript")
  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  cases <- list(
    c("03_qc_exploration/04b_import_reference_marker_sources.r", "--dry-run"),
    c("03_qc_exploration/05_empirical_roi_marker_discovery.r", "--dry-run"),
    c("03_qc_exploration/06_wgcna_marker_trait_export.r", "--dataset", "microglia", "--dry-run"),
    c("06_modules_WGCNA/05_module_supermodule_group_effects.r", "--dataset", "microglia", "--dry-run"),
    c("06_modules_WGCNA/06_annotate_module_microenvironment.r", "--dataset", "microglia", "--dry-run"),
    c("06_modules_WGCNA/07_wgcna_interpretable_summary.r", "--dataset", "all", "--dry-run"),
    c("06_modules_WGCNA/08_microglia_neuropil_independence.R", "--dataset", "microglia", "--dry-run"),
    c("06_modules_WGCNA/08_module_complex_architecture.r", "--dataset", "all", "--dry-run"),
    c("06_modules_WGCNA/09_module_robustness_sensitivity.r", "--dataset", "all", "--dry-run")
  )
  for (args in cases) {
    out <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
    testthat::expect_equal(attr(out, "status") %||% 0L, 0L, info = paste(args, collapse = " "))
    testthat::expect_true(any(grepl("\\[DRY-RUN", out)), info = paste(args, collapse = " "))
  }
})

testthat::test_that("marker source manifest is parseable", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  manifest <- repo_path("data", "external", "reference_markers", "reference_marker_sources.yml")
  testthat::expect_true(file.exists(manifest))
  testthat::skip_if_not_installed("yaml")
  parsed <- yaml::read_yaml(manifest)
  testthat::expect_true(length(parsed$sources) >= 5)
  testthat::expect_true(all(c("ewce", "zeisel_linnarsson", "dropviz_saunders", "allen_celltypes", "brainrnaseq_barres") %in% vapply(parsed$sources, `[[`, character(1), "source_name")))
  allen <- parsed$sources[[which(vapply(parsed$sources, `[[`, character(1), "source_name") == "allen_celltypes")]]
  testthat::expect_equal(allen$local_dir, "data/external/reference_markers/derived")
  testthat::expect_true("allen_mouse_ctx_hip_10x_marker_panels_top.csv" %in% unlist(allen$file_patterns, use.names = FALSE))
})

testthat::test_that("marker registry helpers load registry, empirical sets, and legacy fallback", {
  source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
  registry_file <- tempfile(fileext = ".csv")
  empirical_file <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    marker_set = "canonical_microglia_homeostatic",
    cell_class = "microglia",
    cell_state = "homeostatic_identity",
    gene_symbol = c("P2ry12", "Tmem119"),
    source_type = "test",
    source_name = "test_registry",
    source_reference = "test",
    selection_rule = "test",
    confidence = "test",
    use_for = "annotation",
    notes = "",
    stringsAsFactors = FALSE
  ), registry_file, row.names = FALSE)
  utils::write.csv(data.frame(
    marker_set = "empirical_microglia_roi_enriched",
    ProteinID = "Aif1",
    GeneSymbol = "Aif1",
    marker_source = "test_empirical",
    stringsAsFactors = FALSE
  ), empirical_file, row.names = FALSE)
  old_registry <- Sys.getenv("PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE", unset = NA_character_)
  old_empirical <- Sys.getenv("PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE", unset = NA_character_)
  on.exit({
    if (is.na(old_registry)) Sys.unsetenv("PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE") else Sys.setenv(PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE = old_registry)
    if (is.na(old_empirical)) Sys.unsetenv("PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE") else Sys.setenv(PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE = old_empirical)
  }, add = TRUE)
  Sys.setenv(PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE = registry_file, PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE = empirical_file)
  registry <- read_wgcna_marker_registry()
  testthat::expect_true(all(c("marker_set", "gene_symbol", "gene_token") %in% names(registry)))
  sets <- load_wgcna_marker_sets(quiet = TRUE)
  testthat::expect_true("canonical_microglia_homeostatic" %in% names(sets))
  testthat::expect_true("empirical_microglia_roi_enriched" %in% names(sets))

  Sys.setenv(PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE = tempfile(), PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE = tempfile())
  fallback <- suppressWarnings(load_wgcna_marker_sets(quiet = TRUE))
  testthat::expect_true("microglia" %in% names(fallback))
  testthat::expect_equal(attr(fallback, "marker_source_metadata")$marker_source[[1]], "legacy_hardcoded_fallback")
})

testthat::test_that("supermodule display labels keep immutable IDs and cap singleton confidence", {
  source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
  testthat::expect_equal(
    classify_supermodule_label_confidence(
      n_modules = 1L,
      go_class = "GO_supported",
      has_coherent_hubs = TRUE,
      microenvironment_class = "vascular_basement_membrane_ecm"
    ),
    "low"
  )
  lbl <- compose_supermodule_display_label("SM04", "Perivascular ECM")
  testthat::expect_match(lbl, "^SM04\\s+")
  testthat::expect_match(lbl, "Perivascular ECM")
})

testthat::test_that("semantic classifier treats synaptic adhesion scaffold as non-ECM", {
  source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
  assigned <- wgcna_assign_semantic_program(
    bp_terms = "regulation of monoatomic ion transport",
    mf_terms = c("PDZ domain binding", "cell-cell adhesion mediator activity"),
    cc_terms = c("presynaptic membrane", "postsynaptic density", "synaptic membrane")
  )
  testthat::expect_equal(
    assigned$module_specific_label,
    "synaptic membrane / ion-transport regulatory scaffold"
  )
  testthat::expect_equal(assigned$module_program_primary, "synaptic/cytoskeletal trafficking")
  testthat::expect_false(grepl("ECM/adhesion|Perivascular ECM / adhesion", paste(unlist(assigned), collapse = " "), ignore.case = TRUE))
})

testthat::test_that("legacy static supermodule seeds stay opt-in", {
  script <- readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE)
  txt <- paste(script, collapse = "\n")
  testthat::expect_true(grepl('PROTEOMICS_ALLOW_LEGACY_SUPERMODULE_SEED", unset = "false"', txt, fixed = TRUE))
  testthat::expect_match(txt, "legacy_static_seed")
})

testthat::test_that("supermodule sensitivity export uses fixed primary cut height", {
  script <- readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE)
  txt <- paste(script, collapse = "\n")
  testthat::expect_match(txt, "primary_supermodule_cut_height <- cut_height")
  testthat::expect_match(txt, "primary_cut_height = primary_supermodule_cut_height")
  testthat::expect_false(grepl("primary_cut_height\\s*=\\s*cut_height", txt))
  testthat::expect_match(txt, "primary_cut_height must contain exactly one finite selected primary cut height")
  testthat::expect_match(txt, "dataset_profile_resolved")
})

testthat::test_that("downstream supermodule labels prefer display label consistently", {
  scripts <- c(
    testthat::test_path("..", "..", "06_modules_WGCNA", "05_module_supermodule_group_effects.r"),
    testthat::test_path("..", "..", "06_modules_WGCNA", "07_wgcna_interpretable_summary.r"),
    testthat::test_path("..", "..", "09_export_pride_journal", "07_make_biological_claims_table.R")
  )
  txt <- paste(vapply(scripts, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
  testthat::expect_match(txt, "Supermodule_DisplayLabel")
  testthat::expect_true(grepl("SupermoduleLabel = dplyr::coalesce\\(as.character\\(\\.data\\$Supermodule_DisplayLabel\\), as.character\\(\\.data\\$Supermodule_FinalLabel\\), as.character\\(\\.data\\$Macroprogram_Display\\)", txt))
  testthat::expect_true(grepl("biological_program = dplyr::coalesce\\(", txt))
  testthat::expect_true(grepl("blank_to_na\\(\\.data\\$Supermodule_DisplayLabel\\)", txt))
  testthat::expect_true(grepl("blank_to_na\\(\\.data\\$Supermodule_FinalLabel\\)", txt))
  testthat::expect_true(grepl("blank_to_na\\(\\.data\\$Macroprogram_Display\\)", txt))
})

testthat::test_that("WGCNA downstream schemas expose required columns", {
  source(testthat::test_path("..", "..", "R", "module_contracts.R"))
  group_cols <- c(
    "dataset", "level", "module_id", "supermodule_id", "module_label",
    "supermodule_label", "spatial_unit", "contrast", "estimate", "SE",
    "statistic", "p_value", "FDR_within_dataset_level", "FDR_global",
    "direction", "effect_scope", "SpatialUnitType", "model_type",
    "has_repeated_animals", "n_animals", "n_samples", "formula_requested",
    "formula_used", "dropped_covariates",
    "rank_deficient_model", "model_warning"
  )
  group_df <- as.data.frame(setNames(rep(list(logical()), length(group_cols)), group_cols))
  testthat::expect_silent(validate_wgcna_group_effects(group_df))

  annot_df <- data.frame(
    dataset = character(), ModuleID = character(), ModuleColor = character(),
    n_proteins = integer(), microenvironment_class = character(),
    microglia_evidence = logical(), neuropil_evidence = logical(),
    other_cellular_evidence = logical(), classification_threshold = numeric(),
    canonical_microglia_evidence = logical(),
    empirical_microglia_roi_evidence = logical(),
    canonical_neuropil_evidence = logical(),
    empirical_neuropil_evidence = logical(),
    empirical_shared_microenvironment_evidence = logical(),
    microglia_state_or_activation_evidence = logical(),
    peripheral_myeloid_caution_evidence = logical(),
    marker_registry_version = character(),
    empirical_marker_set_version = character(),
    classification_rationale = character(),
    interpretation_note = character()
  )
  testthat::expect_silent(validate_wgcna_module_annotation(annot_df))
  testthat::expect_true(all(c(
    "microglia_supported", "microglia_state_or_activation_supported", "shared_microenvironment", "neuropil_sensitive",
    "other_cellular_or_vascular_sensitive", "ambiguous"
  ) %in% c(
    "microglia_supported", "microglia_state_or_activation_supported", "shared_microenvironment", "neuropil_sensitive",
    "other_cellular_or_vascular_sensitive", "ambiguous"
  )))

  interp_df <- data.frame(
    dataset = character(), level = character(), contrast = character(),
    estimate = numeric(), p_value = numeric(), FDR_global = numeric(),
    interpretation_sentence = character()
  )
  testthat::expect_silent(validate_wgcna_interpretable_summary(interp_df))

  atlas_df <- data.frame(
    dataset = character(), evidence_domain = character(), evidence_id = character(),
    program_label = character(), entity_type = character(), entity_id = character(),
    source_file = character(), evidence_status = character(),
    interpretation_note = character(), qc_flag = character()
  )
  testthat::expect_silent(validate_cross_compartment_program_atlas(atlas_df))

  summary_df <- data.frame(
    program_key = character(), manuscript_claim_scope = character(),
    datasets_supported = character(), evidence_domains = character(),
    strongest_evidence = character(), safe_manuscript_sentence = character(),
    main_limitation = character(), qc_flag = character()
  )
  testthat::expect_silent(validate_manuscript_program_summary(summary_df))

  priority_df <- data.frame(
    priority_id = character(), program_key = character(), dataset = character(),
    priority_tier = character(), evidence_domain_count = integer(),
    strongest_fdr = numeric(), robustness_flag = character(),
    behavior_flag = character(), qc_flag = character(), recommended_use = character()
  )
  testthat::expect_silent(validate_evidence_priority_matrix(priority_df))
})
