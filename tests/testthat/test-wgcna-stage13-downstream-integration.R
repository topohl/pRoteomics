stage13_repo_root <- normalizePath(
  testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE
)
stage13_test_path <- function(...) file.path(stage13_repo_root, ...)

read_required_csv <- function(path) {
  testthat::skip_if_not(file.exists(path), paste("Required generated artifact is unavailable:", path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

require_stage13_generated_inputs <- function(datasets = "microglia") {
  paths <- c(
    stage13_test_path(
      "results", "tables", "06_modules_WGCNA", "claim_readiness",
      "microglia", "WGCNA_entity_claim_readiness.csv"
    ),
    vapply(datasets, function(dataset) {
      stage13_test_path(
        "results", "tables", "06_modules_WGCNA", "interpretable_summary",
        dataset, "WGCNA_inferential_handoff.csv"
      )
    }, character(1))
  )
  missing <- paths[!file.exists(paths)]
  testthat::skip_if(
    length(missing) > 0L,
    paste(
      "Required generated Stage 07/13 artifact(s) are unavailable:",
      paste(missing, collapse = "; ")
    )
  )
}

testthat::test_that("Stage 13 helper enforces the finalized microglia identity contract", {
  require_stage13_generated_inputs("microglia")
  generated_stage13 <- read_required_csv(stage13_test_path(
    "results", "tables", "06_modules_WGCNA", "claim_readiness",
    "microglia", "WGCNA_entity_claim_readiness.csv"
  ))
  testthat::skip_if_not(
    all(c(
      "group_effect_handoff_file", "group_effect_source_artifact",
      "group_effect_source_key_contract"
    ) %in% names(generated_stage13)),
    "Generated Stage 13 output predates direct Stage 07 provenance carry-through"
  )
  source(stage13_test_path("R", "paths.R"))
  source(stage13_test_path("R", "wgcna_claim_readiness_utils.R"))
  contract <- load_microglia_wgcna_claim_readiness()
  x <- contract$all

  expected_aliases <- c(
    SM02 = "WGCNA_m03", SM04 = "WGCNA_m09", SM05 = "WGCNA_m13",
    SM06 = "WGCNA_m01", SM07 = "WGCNA_m02", SM08 = "WGCNA_m06"
  )
  aliases <- x[x$claim_entity_role == "compatibility_alias", ]
  blocks <- x[x$claim_entity_role == "higher_order_block", ]
  modules <- x[x$claim_entity_role == "canonical_module", ]

  testthat::expect_equal(nrow(x), 22L)
  testthat::expect_equal(nrow(modules), 13L)
  testthat::expect_equal(nrow(blocks), 3L)
  testthat::expect_equal(nrow(aliases), 6L)
  testthat::expect_equal(stats::setNames(aliases$canonical_claim_entity_id, aliases$entity_id)[names(expected_aliases)], expected_aliases)
  testthat::expect_false(any(aliases$separate_manuscript_claim_allowed))
  testthat::expect_setequal(blocks$entity_id[blocks$separate_manuscript_claim_allowed], "SM09")
  testthat::expect_equal(sum(modules$separate_manuscript_claim_allowed), 10L)
  testthat::expect_false(any(x$group_effect_status == "FDR_supported"))
  testthat::expect_identical(contract$source_path, normalizePath(microglia_wgcna_claim_readiness_path(), winslash = "/", mustWork = TRUE))
})

testthat::test_that("integration separates validated, diagnostic, annotation, unavailable, and neuronal WGCNA roles", {
  atlas <- read_required_csv(stage13_test_path(
    "results", "tables", "10_biological_integration", "cross_compartment_program_atlas",
    "global", "cross_compartment_program_atlas_long.csv"
  ))
  micro <- atlas[atlas$dataset == "microglia" & nzchar(atlas$claim_entity_role), ]
  aliases <- micro[micro$claim_entity_role == "compatibility_alias", ]
  blocks <- micro[micro$claim_entity_role == "higher_order_block", ]

  testthat::expect_equal(nrow(micro), 22L)
  testthat::expect_equal(nrow(aliases), 6L)
  testthat::expect_false(any(aliases$counts_toward_convergence))
  testthat::expect_false(any(blocks$counts_toward_convergence[blocks$entity_id %in% c("SM01", "SM03")]))
  testthat::expect_true(blocks$counts_toward_convergence[blocks$entity_id == "SM09"])
  testthat::expect_false(any(micro$evidence_semantic_class == "wgcna_stress_group_effect"))
  testthat::expect_equal(sum(micro$counts_toward_convergence), 11L)
  testthat::expect_equal(sum(micro$evidence_role == "validated_wgcna_architecture"), 16L)
  micro_context <- atlas[
    atlas$dataset == "microglia" & atlas$evidence_domain %in% c(
      "microenvironment_marker", "complex_organelle", "robustness_sensitivity",
      "behavior_coupling", "qc_confounding"
    ),
  ]
  testthat::expect_false(any(micro_context$counts_toward_convergence))

  neuronal_wgcna <- atlas[
    atlas$dataset %in% c("neuron_neuropil", "neuron_soma") &
      atlas$evidence_domain == "wgcna_supermodule",
  ]
  testthat::expect_gt(sum(neuronal_wgcna$dataset == "neuron_soma"), 0L)
  testthat::expect_gt(sum(neuronal_wgcna$dataset == "neuron_neuropil"), 0L)
  testthat::expect_false(any(neuronal_wgcna$counts_toward_convergence))
  testthat::expect_true(all(neuronal_wgcna$evidence_role == "unvalidated_neuronal_wgcna"))
  testthat::expect_true(all(neuronal_wgcna$evidence_role_reason == "neuronal_wgcna_readiness_contract_not_available"))

  neuronal_enrichment <- atlas[
    atlas$dataset %in% c("neuron_neuropil", "neuron_soma") & atlas$evidence_domain == "enrichment_program",
  ]
  neuronal_spatial <- atlas[
    atlas$dataset %in% c("neuron_neuropil", "neuron_soma") & atlas$evidence_domain == "spatial_architecture",
  ]
  testthat::expect_true(any(neuronal_enrichment$counts_toward_convergence))
  testthat::expect_true(any(neuronal_spatial$counts_toward_convergence))
  for (domain in c("microenvironment_marker", "complex_organelle", "robustness_sensitivity", "behavior_coupling", "qc_confounding")) {
    rows <- atlas[atlas$evidence_domain == domain, ]
    testthat::expect_false(any(rows$counts_toward_convergence))
  }
  old_wd <- setwd(stage13_repo_root)
  on.exit(setwd(old_wd), add = TRUE)
  source(stage13_test_path("R", "integration_utils.R"))
  unavailable <- assign_downstream_evidence_roles(availability_evidence(
    "neuron_soma", "external_signature_overlap", "missing.csv", "Optional input unavailable."
  ))
  testthat::expect_identical(unavailable$evidence_role, "unavailable_evidence")
  testthat::expect_false(unavailable$counts_toward_convergence)
  placeholder <- assign_downstream_evidence_roles(data.frame(
    dataset = "neuron_soma", evidence_domain = "external_signature_overlap",
    program_label = "Unavailable optional evidence", entity_type = "program",
    evidence_status = "external_signature_overlap", interpretation_note = "External signature evidence imported for atlas.",
    stringsAsFactors = FALSE
  ))
  testthat::expect_identical(placeholder$evidence_role, "unavailable_evidence")
  testthat::expect_false(placeholder$counts_toward_convergence)

  summary <- read_required_csv(stage13_test_path(
    "results", "tables", "10_biological_integration", "manuscript_program_summary",
    "global", "manuscript_program_summary.csv"
  ))
  testthat::expect_true(all(summary$n_evidence_rows_counting_toward_convergence <= summary$n_evidence_rows_total))
  expected_counting <- aggregate(
    atlas$counts_toward_convergence,
    list(program_key = atlas$program_key),
    function(z) sum(z %in% TRUE)
  )
  observed <- summary$n_evidence_rows_counting_toward_convergence[match(expected_counting$program_key, summary$program_key)]
  testthat::expect_equal(observed, expected_counting$x)
})

testthat::test_that("neuronal WGCNA architecture remains descriptive and disallowed", {
  require_stage13_generated_inputs(c(
    "microglia", "neuron_soma", "neuron_neuropil"
  ))
  claims <- read_required_csv(stage13_test_path("results", "tables", "biological_claims_table.csv"))
  architecture <- claims[
    claims$dataset %in% c("neuron_soma", "neuron_neuropil") & claims$claim_type == "wgcna_architecture",
  ]
  testthat::expect_equal(sum(architecture$dataset == "neuron_soma"), 7L)
  testthat::expect_equal(sum(architecture$dataset == "neuron_neuropil"), 15L)
  testthat::expect_false(any(architecture$claim_allowed))
  testthat::expect_true(all(architecture$claim_gate_status == "disallowed"))
  testthat::expect_true(all(architecture$claim_use_class == "annotation_only"))
  testthat::expect_true(all(grepl("neuronal_wgcna_architecture_readiness_not_available", architecture$claim_downgrade_reason, fixed = TRUE)))
  testthat::expect_true(all(grepl("not a stress-effect, network-remodelling, mechanistic, or independently validated architecture claim", architecture$safe_interpretation, fixed = TRUE)))

  group_effects <- claims[
    claims$dataset %in% c("neuron_soma", "neuron_neuropil") & claims$claim_type == "wgcna_group_effect",
  ]
  testthat::expect_true(nrow(group_effects) > 0L)
  testthat::expect_false(any(group_effects$claim_allowed))
})

testthat::test_that("microglia biological claims use Stage 13 stable identities without aliases", {
  require_stage13_generated_inputs("microglia")
  claims <- read_required_csv(stage13_test_path("results", "tables", "biological_claims_table.csv"))
  w <- claims[
    claims$dataset == "microglia" &
      claims$claim_type %in% c("wgcna_architecture", "wgcna_group_effect", "wgcna_convergent_overlap", "wgcna_compatibility_provenance"),
  ]

  testthat::expect_false(any(w$wgcna_claim_entity_role == "compatibility_alias", na.rm = TRUE))
  testthat::expect_false(any(w$wgcna_entity_id %in% c("SM01", "SM03"), na.rm = TRUE))
  group_effects <- w[w$claim_type == "wgcna_group_effect", ]
  testthat::expect_true(nrow(group_effects) > 0L)
  testthat::expect_false(any(group_effects$claim_allowed))
  testthat::expect_true(all(group_effects$wgcna_group_effect_status == "not_FDR_supported"))

  architecture <- w[w$claim_type == "wgcna_architecture", ]
  testthat::expect_setequal(unique(architecture$wgcna_entity_id[architecture$wgcna_level == "supermodule"]), "SM09")
  testthat::expect_equal(length(unique(architecture$wgcna_entity_id[architecture$wgcna_level == "module"])), 10L)
  testthat::expect_false(all(architecture$claim_use_class == "supporting_claim"))
  testthat::expect_false(any(grepl("supports? (stress-dependent regulation|broad network remodelling|cell-intrinsic|caus)", architecture$safe_interpretation, ignore.case = TRUE)))
  testthat::expect_true(any(grepl("covariance|spatial|architecture", architecture$safe_interpretation, ignore.case = TRUE)))

  cardinality <- read_required_csv(stage13_test_path("results", "reviewer_audit", "wgcna_stage13_claim_cardinality_audit.csv"))
  testthat::expect_true(all(cardinality$duplicate_key_is_explicit_endpoints))
})

testthat::test_that("descriptive overviews preserve complete identities while claim selection is canonical", {
  segments <- read_required_csv(stage13_test_path(
    "results", "source_data", "10_biological_integration", "wgcna_circular_atlas",
    "global", "wgcna_circular_atlas_segments.csv"
  ))
  micro <- segments[segments$dataset == "microglia", ]
  expected_segments <- sum(vapply(
    c("neuron_neuropil", "neuron_soma", "microglia"),
    function(dataset) {
      lookup <- read_required_csv(stage13_test_path(
        "results", "tables", "06_modules_WGCNA", "interpretable_summary",
        dataset, "WGCNA_final_label_lookup.csv"
      ))
      sum(lookup$level == "supermodule")
    },
    integer(1)
  ))
  testthat::expect_equal(nrow(segments), expected_segments)
  testthat::expect_equal(nrow(micro), 9L)
  testthat::expect_equal(length(unique(segments$dataset)), 3L)
  testthat::expect_setequal(micro$supermodule_id[micro$selected_for_manuscript_claim], "SM09")
  testthat::expect_false(any(micro$claim_entity_role[micro$selected_for_manuscript_claim] == "compatibility_alias"))
  testthat::expect_true(all(!is.na(micro$claim_entity_role)))
  neuronal <- segments[segments$dataset %in% c("neuron_soma", "neuron_neuropil"), ]
  testthat::expect_false(any(neuronal$selected_for_manuscript_claim[grepl("(^|;)\\s*disallowed\\s*(;|$)", neuronal$claim_display_status)]))

  metrics <- read_required_csv(stage13_test_path(
    "results", "tables", "10_biological_integration", "wgcna_cross_compartment_overview",
    "global", "wgcna_cross_compartment_metrics.csv"
  ))
  m <- metrics[metrics$dataset == "microglia", ]
  testthat::expect_equal(m$n_stage13_technical_identities, 22L)
  testthat::expect_equal(m$n_canonical_modules, 13L)
  testthat::expect_equal(m$n_higher_order_blocks, 3L)
  testthat::expect_equal(m$n_compatibility_aliases, 6L)
  testthat::expect_equal(m$n_independently_eligible_modules, 10L)
  testthat::expect_equal(m$n_independently_eligible_blocks, 1L)
  testthat::expect_equal(m$n_FDR_supported_stage13_group_effects, 0L)
})

testthat::test_that("final bundle preserves all Stage 13 rows and filters manuscript candidates", {
  bundle <- stage13_test_path("results", "tables", "10_biological_integration", "final_evidence_bundle", "global")
  readiness <- read_required_csv(file.path(bundle, "microglia_wgcna_claim_readiness.csv"))
  modules <- read_required_csv(file.path(bundle, "wgcna_key_modules.csv"))
  supermodules <- read_required_csv(file.path(bundle, "wgcna_key_supermodules.csv"))
  micro_modules <- modules[modules$dataset == "microglia", ]
  micro_supermodules <- supermodules[supermodules$dataset == "microglia", ]

  testthat::expect_equal(nrow(readiness), 22L)
  testthat::expect_equal(sum(readiness$claim_entity_role == "compatibility_alias"), 6L)
  testthat::expect_equal(length(unique(micro_modules$module_id)), 10L)
  testthat::expect_false(any(micro_modules$claim_entity_role == "compatibility_alias"))
  testthat::expect_setequal(unique(micro_supermodules$supermodule_id), "SM09")
  testthat::expect_false(any(micro_supermodules$supermodule_id %in% c("SM01", "SM03")))

  neuronal_modules <- modules[modules$dataset %in% c("neuron_soma", "neuron_neuropil"), ]
  neuronal_supermodules <- supermodules[supermodules$dataset %in% c("neuron_soma", "neuron_neuropil"), ]
  testthat::skip_if(
    any(grepl("Stage 05 models", c(
      neuronal_modules$reason, neuronal_supermodules$reason
    ), fixed = TRUE)),
    "Generated final bundle predates the Stage 07 inferential-handoff wording"
  )
  testthat::expect_equal(nrow(neuronal_modules), 2L)
  testthat::expect_equal(nrow(neuronal_supermodules), 2L)
  testthat::expect_true(all(neuronal_modules$status == "no_validated_neuronal_wgcna_key_rows"))
  testthat::expect_true(all(neuronal_supermodules$status == "no_validated_neuronal_wgcna_key_rows"))
  testthat::expect_true(all(neuronal_modules$reason == "neuronal readiness contract unavailable and all current Stage 07 inferential rows are claim-ineligible"))
  testthat::expect_true(all(neuronal_supermodules$reason == "neuronal readiness contract unavailable and all current Stage 07 inferential rows are claim-ineligible"))
})

testthat::test_that("circular claim parsing uses exact status tokens and frozen neuronal states are unchanged", {
  circular_script <- paste(readLines(stage13_test_path("10_biological_integration", "04_wgcna_circular_atlas.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(circular_script, 'tokens == "disallowed"', fixed = TRUE)
  testthat::expect_match(circular_script, 'tokens %in% c("allowed", "downgraded")', fixed = TRUE)
  testthat::expect_false(grepl('grepl("allowed|downgraded", .data$claim_display_status)', circular_script, fixed = TRUE))

  sha256_file <- function(path) {
    raw <- readBin(path, "raw", n = file.info(path)$size)
    paste(format(openssl::sha256(raw)), collapse = "")
  }
  soma_state <- stage13_test_path("data", "processed", "06_modules_WGCNA", "01_WGCNA", "neuron_soma", "wgcna_final_model_state.rds")
  neuropil_state <- stage13_test_path("data", "processed", "06_modules_WGCNA", "01_WGCNA", "neuron_neuropil", "wgcna_final_model_state.rds")
  testthat::skip_if_not(file.exists(soma_state) && file.exists(neuropil_state), "Frozen neuronal WGCNA states are unavailable")
  testthat::expect_identical(sha256_file(soma_state), "f822be9cbf8589be1aef2d6d018b5ee554a08235c699f3dc06bcb45110c0ca73")
  testthat::expect_identical(sha256_file(neuropil_state), "620f8f09ab10e34195e9482ad711609ac3f5f0d6d34ef38bc09f72adb6ed8d4e")
})
