source(testthat::test_path("..", "..", "R", "paths.R"))

stage13_repo_root <- normalizePath(repo_root(), winslash = "/", mustWork = TRUE)
stage13_test_path <- function(...) file.path(stage13_repo_root, ...)

read_required_csv <- function(path) {
  testthat::skip_if_not(file.exists(path), paste("Required generated artifact is unavailable:", path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# Current canonical producers:
#   Stage 13 readiness  analysis/wgcna/audit_module_claim_readiness.R
#   Stage 07 handoff    analysis/wgcna/summarize_module_interpretation.R
stage13_readiness_file <- function() {
  current_wgcna("audit_module_claim_readiness", "microglia", "WGCNA_entity_claim_readiness.csv")
}
stage07_file <- function(dataset, filename) {
  current_wgcna("summarize_module_interpretation", dataset, filename)
}

require_stage13_generated_inputs <- function(datasets = "microglia") {
  paths <- c(
    stage13_readiness_file(),
    vapply(datasets, stage07_file, character(1), filename = "WGCNA_inferential_handoff.csv")
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

# The Stage 13 roles the current network implies, and the claim decisions the
# current readiness table made. Downstream consumers must reproduce these; they
# are never pinned.
stage13_expected <- function() {
  shape <- current_network_shape("microglia")
  x <- read_required_csv(stage13_readiness_file())
  claimable <- x$separate_manuscript_claim_allowed %in% TRUE
  list(
    shape = shape,
    readiness = x,
    claimable_modules = sort(x$entity_id[claimable & x$claim_entity_role == "canonical_module"]),
    claimable_blocks = sort(x$entity_id[claimable & x$claim_entity_role == "higher_order_block"]),
    unclaimable_blocks = sort(x$entity_id[!claimable & x$claim_entity_role == "higher_order_block"]),
    n_fdr_supported = sum(x$group_effect_status == "FDR_supported", na.rm = TRUE)
  )
}

testthat::test_that("Stage 13 helper enforces the finalized microglia identity contract", {
  require_stage13_generated_inputs("microglia")
  generated_stage13 <- read_required_csv(stage13_readiness_file())
  testthat::skip_if_not(
    all(c(
      "group_effect_handoff_file", "group_effect_source_artifact",
      "group_effect_source_key_contract"
    ) %in% names(generated_stage13)),
    "Generated Stage 13 output predates direct Stage 07 provenance carry-through"
  )
  source(repo_path("R", "paths.R"))
  source(repo_path("R", "wgcna_claim_readiness_utils.R"))
  contract <- load_microglia_wgcna_claim_readiness()
  x <- contract$all
  shape <- current_network_shape("microglia")

  aliases <- x[x$claim_entity_role == "compatibility_alias", ]
  blocks <- x[x$claim_entity_role == "higher_order_block", ]
  modules <- x[x$claim_entity_role == "canonical_module", ]

  testthat::expect_equal(nrow(x), shape$n_identities)
  testthat::expect_setequal(modules$entity_id, shape$modules)
  testthat::expect_setequal(blocks$entity_id, shape$blocks)
  testthat::expect_setequal(aliases$entity_id, names(shape$singletons))
  testthat::expect_equal(
    stats::setNames(aliases$canonical_claim_entity_id, aliases$entity_id)[names(shape$singletons)],
    shape$singletons
  )
  testthat::expect_false(any(aliases$separate_manuscript_claim_allowed))
  # support status follows the tier-specific FDR of the independent endpoint
  independent <- x$claim_entity_role != "compatibility_alias" & !is.na(x$group_effect_tier_specific_fdr)
  testthat::expect_false(any(x$group_effect_status[independent] == "FDR_supported" &
                               x$group_effect_tier_specific_fdr[independent] > 0.05))
  # the helper resolved the CURRENT canonical table, not a historical one
  testthat::expect_identical(contract$source_path, normalizePath(microglia_wgcna_claim_readiness_path(), winslash = "/", mustWork = TRUE))
  testthat::expect_identical(contract$source_path, normalizePath(stage13_readiness_file(), winslash = "/", mustWork = TRUE))
})

testthat::test_that("integration separates validated, diagnostic, annotation, unavailable, and neuronal WGCNA roles", {
  # current canonical atlas (analysis/integration/build_cross_compartment_atlas.R)
  atlas <- read_required_csv(current_integration(
    "build_cross_compartment_atlas", "cross_compartment_program_atlas_long.csv"
  ))
  exp <- stage13_expected()
  micro <- atlas[atlas$dataset == "microglia" & nzchar(atlas$claim_entity_role), ]
  aliases <- micro[micro$claim_entity_role == "compatibility_alias", ]
  blocks <- micro[micro$claim_entity_role == "higher_order_block", ]

  testthat::expect_equal(nrow(micro), exp$shape$n_identities)
  testthat::expect_equal(nrow(aliases), length(exp$shape$singletons))
  testthat::expect_false(any(aliases$counts_toward_convergence))
  testthat::expect_false(any(blocks$counts_toward_convergence))
  testthat::expect_false(any(micro$evidence_semantic_class == "wgcna_stress_group_effect"))
  testthat::expect_equal(sum(micro$counts_toward_convergence), 0L)
  testthat::expect_equal(sum(micro$evidence_role == "validated_wgcna_architecture"),
                         length(exp$shape$modules) + length(exp$shape$blocks))
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
  testthat::expect_true(all(neuronal_enrichment$evidence_source_family == "ranked_GSEA"))
  testthat::expect_true(all(neuronal_spatial$evidence_source_family == "ranked_GSEA"))
  for (domain in c("microenvironment_marker", "complex_organelle", "robustness_sensitivity", "behavior_coupling", "qc_confounding")) {
    rows <- atlas[atlas$evidence_domain == domain, ]
    testthat::expect_false(any(rows$counts_toward_convergence))
  }
  old_wd <- setwd(stage13_repo_root)
  on.exit(setwd(old_wd), add = TRUE)
  source(repo_path("R", "integration_utils.R"))
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

  # current canonical summary (analysis/integration/summarize_programs_for_manuscript.R)
  summary <- read_required_csv(current_integration(
    "summarize_programs_for_manuscript", "manuscript_program_summary.csv"
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
  claims <- read_required_csv(current_claims_table())
  architecture <- claims[
    claims$dataset %in% c("neuron_soma", "neuron_neuropil") & claims$claim_type == "wgcna_architecture",
  ]
  # one architecture row per current module
  testthat::expect_equal(sum(architecture$dataset == "neuron_soma"),
                         length(current_network_shape("neuron_soma")$modules))
  testthat::expect_equal(sum(architecture$dataset == "neuron_neuropil"),
                         length(current_network_shape("neuron_neuropil")$modules))
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
  exp <- stage13_expected()
  claims <- read_required_csv(current_claims_table())
  w <- claims[
    claims$dataset == "microglia" &
      claims$claim_type %in% c("wgcna_architecture", "wgcna_group_effect", "wgcna_convergent_overlap", "wgcna_compatibility_provenance"),
  ]

  testthat::expect_false(any(w$wgcna_claim_entity_role == "compatibility_alias", na.rm = TRUE))
  testthat::expect_false(any(w$wgcna_entity_id %in% exp$unclaimable_blocks, na.rm = TRUE))
  group_effects <- w[w$claim_type == "wgcna_group_effect", ]
  testthat::expect_true(nrow(group_effects) > 0L)
  testthat::expect_false(any(group_effects$claim_allowed))
  # no claims-table row may report support that Stage 13 does not: each row's status is its
  # entity's Stage 13 status (build_biological_claims_table.R joins it by dataset, level, entity)
  readiness <- read_required_csv(current_wgcna("audit_module_claim_readiness", "microglia",
                                               "WGCNA_entity_claim_readiness.csv"))
  status_col <- intersect(c("group_effect_status", "wgcna_group_effect_status"), names(readiness))[1]
  testthat::expect_false(is.na(status_col))
  k <- match(group_effects$wgcna_entity_id, readiness$entity_id)
  testthat::expect_false(anyNA(k))
  supported_here <- group_effects$wgcna_group_effect_status == "FDR_supported"
  testthat::expect_true(all(readiness[[status_col]][k][supported_here] == "FDR_supported"))

  architecture <- w[w$claim_type == "wgcna_architecture", ]
  testthat::expect_setequal(unique(architecture$wgcna_entity_id[architecture$wgcna_level == "supermodule"]), exp$claimable_blocks)
  testthat::expect_setequal(unique(architecture$wgcna_entity_id[architecture$wgcna_level == "module"]), exp$claimable_modules)
  testthat::expect_false(all(architecture$claim_use_class == "supporting_claim"))
  testthat::expect_false(any(grepl("supports? (stress-dependent regulation|broad network remodelling|cell-intrinsic|caus)", architecture$safe_interpretation, ignore.case = TRUE)))
  testthat::expect_true(any(grepl("covariance|spatial|architecture", architecture$safe_interpretation, ignore.case = TRUE)))

  cardinality <- read_required_csv(current_claims_table("wgcna_stage13_claim_cardinality_audit.csv"))
  testthat::expect_true(all(cardinality$duplicate_key_is_explicit_endpoints))
})

testthat::test_that("descriptive overviews preserve complete identities while claim selection is canonical", {
  # current canonical circular atlas (analysis/integration/render_module_circular_atlas.R)
  segments <- read_required_csv(current_integration(
    "render_module_circular_atlas", "source_data", "wgcna_circular_atlas_segments.csv"
  ))
  exp <- stage13_expected()
  micro <- segments[segments$dataset == "microglia", ]
  expected_segments <- sum(vapply(
    wgcna_test_datasets,
    function(dataset) {
      lookup <- read_required_csv(stage07_file(dataset, "WGCNA_final_label_lookup.csv"))
      sum(lookup$level == "supermodule")
    },
    integer(1)
  ))
  testthat::expect_equal(nrow(segments), expected_segments)
  testthat::expect_equal(nrow(micro), length(exp$shape$supermodules))
  testthat::expect_equal(length(unique(segments$dataset)), 3L)
  testthat::expect_setequal(micro$supermodule_id[micro$selected_for_manuscript_claim], exp$claimable_blocks)
  testthat::expect_false(any(micro$claim_entity_role[micro$selected_for_manuscript_claim] == "compatibility_alias"))
  testthat::expect_true(all(!is.na(micro$claim_entity_role)))
  neuronal <- segments[segments$dataset %in% c("neuron_soma", "neuron_neuropil"), ]
  testthat::expect_false(any(neuronal$selected_for_manuscript_claim[grepl("(^|;)\\s*disallowed\\s*(;|$)", neuronal$claim_display_status)]))

  # current canonical overview (analysis/integration/summarize_module_cross_compartment.R)
  metrics <- read_required_csv(current_integration(
    "summarize_module_cross_compartment", "wgcna_cross_compartment_metrics.csv"
  ))
  m <- metrics[metrics$dataset == "microglia", ]
  testthat::expect_equal(m$n_stage13_technical_identities, exp$shape$n_identities)
  testthat::expect_equal(m$n_canonical_modules, length(exp$shape$modules))
  testthat::expect_equal(m$n_higher_order_blocks, length(exp$shape$blocks))
  testthat::expect_equal(m$n_compatibility_aliases, length(exp$shape$singletons))
  testthat::expect_equal(m$n_independently_eligible_modules, length(exp$claimable_modules))
  testthat::expect_equal(m$n_independently_eligible_blocks, length(exp$claimable_blocks))
  testthat::expect_equal(m$n_FDR_supported_stage13_group_effects, exp$n_fdr_supported)
})

testthat::test_that("final bundle preserves all Stage 13 rows and filters manuscript candidates", {
  # current canonical bundle (analysis/integration/build_evidence_priority_matrix.R)
  bundle <- current_integration("build_evidence_priority_matrix", "final_evidence_bundle")
  readiness <- read_required_csv(file.path(bundle, "microglia_wgcna_claim_readiness.csv"))
  modules <- read_required_csv(file.path(bundle, "wgcna_key_modules.csv"))
  supermodules <- read_required_csv(file.path(bundle, "wgcna_key_supermodules.csv"))
  exp <- stage13_expected()
  micro_modules <- modules[modules$dataset == "microglia", ]
  micro_supermodules <- supermodules[supermodules$dataset == "microglia", ]

  testthat::expect_equal(nrow(readiness), exp$shape$n_identities)
  testthat::expect_equal(sum(readiness$claim_entity_role == "compatibility_alias"), length(exp$shape$singletons))
  testthat::expect_setequal(unique(micro_modules$module_id), exp$claimable_modules)
  testthat::expect_false(any(micro_modules$claim_entity_role == "compatibility_alias"))
  testthat::expect_setequal(stats::na.omit(unique(micro_supermodules$supermodule_id)), exp$claimable_blocks)
  testthat::expect_false(any(micro_supermodules$supermodule_id %in% exp$unclaimable_blocks))

  neuronal_modules <- modules[modules$dataset %in% c("neuron_soma", "neuron_neuropil"), ]
  neuronal_supermodules <- supermodules[supermodules$dataset %in% c("neuron_soma", "neuron_neuropil"), ]
  testthat::skip_if(
    any(grepl("Stage 05 models", c(
      neuronal_modules$reason, neuronal_supermodules$reason
    ), fixed = TRUE)),
    "Generated final bundle predates the Stage 07 inferential-handoff wording"
  )
  no_key <- "neuronal readiness contract unavailable and all current Stage 07 inferential rows are claim-ineligible"
  # each neuronal dataset carries either its validated key rows (current
  # supermodule IDs, empty status) or exactly one explicit no-key status row
  for (ds in c("neuron_soma", "neuron_neuropil")) {
    for (tab in list(neuronal_modules, neuronal_supermodules)) {
      rows <- tab[tab$dataset == ds, , drop = FALSE]
      testthat::expect_gt(nrow(rows), 0L)
      placeholder <- rows$status %in% "no_validated_neuronal_wgcna_key_rows"
      if (any(placeholder)) {
        testthat::expect_equal(nrow(rows), 1L, info = ds)
        testthat::expect_identical(rows$reason, no_key)
      }
    }
    key_sm <- neuronal_supermodules$supermodule_id[neuronal_supermodules$dataset == ds &
      !neuronal_supermodules$status %in% "no_validated_neuronal_wgcna_key_rows"]
    testthat::expect_true(all(key_sm %in% current_network_shape(ds)$supermodules), info = ds)
    # key rows are exactly the validated primary SUS-RES rows of the current handoff
    handoff <- read_required_csv(stage07_file(ds, "WGCNA_inferential_handoff.csv"))
    primary <- handoff[handoff$entity_level == "supermodule" &
                         handoff$analysis_tier == "primary_wgcna_global" &
                         handoff$contrast == "SUS - RES" &
                         handoff$result_scope == "primary_global_vulnerability", , drop = FALSE]
    fdr <- suppressWarnings(as.numeric(primary$tier_specific_fdr))
    model_text <- tolower(paste(primary$model_type, primary$support_class, primary$model_warning))
    validated <- !is.na(fdr) & fdr <= 0.10 &
      primary$claim_gate == "eligible_for_readiness_assessment" &
      !grepl("diagnostic|unstable|failed", primary$model_stability_status) &
      !grepl("fallback|diagnostic|model_unstable", model_text)
    testthat::expect_setequal(key_sm, unique(primary$supermodule_id[validated]))
  }
})

testthat::test_that("circular claim parsing uses exact status tokens", {
  circular_script <- paste(readLines(stage13_test_path("analysis/integration", "render_module_circular_atlas.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(circular_script, 'tokens == "disallowed"', fixed = TRUE)
  testthat::expect_match(circular_script, 'tokens %in% c("allowed", "downgraded")', fixed = TRUE)
  testthat::expect_false(grepl('grepl("allowed|downgraded", .data$claim_display_status)', circular_script, fixed = TRUE))
})

testthat::test_that("current network states match the hashes their identity contracts recorded", {
  # Current model-contract consistency (the literal historical neuronal SHA-256
  # snapshot is kept in tests/historical/). The identity contract
  # (analysis/wgcna/build_module_identity_contract.R) records the sha256 of every
  # source it accepted; the current Stage 01 state must still be that file.
  for (ds in wgcna_test_datasets) {
    hashes <- read_required_csv(current_wgcna("build_module_identity_contract", ds,
                                              "WGCNA_identity_source_hashes.csv"))
    state <- hashes[hashes$source_role == "frozen_state", , drop = FALSE]
    testthat::expect_equal(nrow(state), 1L, info = ds)
    testthat::expect_identical(
      state$source_path,
      paste0("results/wgcna/build_wgcna_modules/", ds, "/models/wgcna_final_model_state.rds"),
      info = ds)
    path <- current_wgcna("build_wgcna_modules", ds, "wgcna_final_model_state.rds", child = "models")
    testthat::expect_true(file.exists(path), info = ds)
    testthat::expect_identical(digest::digest(file = path, algo = "sha256"), state$sha256, info = ds)
  }
})
