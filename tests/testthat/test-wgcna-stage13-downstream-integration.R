stage13_test_path <- function(...) testthat::test_path("..", "..", ...)

read_required_csv <- function(path) {
  testthat::skip_if_not(file.exists(path), paste("Required generated artifact is unavailable:", path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

testthat::test_that("Stage 13 helper enforces the finalized microglia identity contract", {
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

testthat::test_that("integration convergence uses eligible Stage 13 architecture only", {
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
  neurons <- atlas[atlas$dataset %in% c("neuron_neuropil", "neuron_soma"), ]
  testthat::expect_true(all(neurons$counts_toward_convergence))

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

testthat::test_that("microglia biological claims use Stage 13 stable identities without aliases", {
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
  testthat::expect_equal(nrow(segments), 27L)
  testthat::expect_equal(nrow(micro), 9L)
  testthat::expect_equal(length(unique(segments$dataset)), 3L)
  testthat::expect_setequal(micro$supermodule_id[micro$selected_for_manuscript_claim], "SM09")
  testthat::expect_false(any(micro$claim_entity_role[micro$selected_for_manuscript_claim] == "compatibility_alias"))
  testthat::expect_true(all(!is.na(micro$claim_entity_role)))

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
})
