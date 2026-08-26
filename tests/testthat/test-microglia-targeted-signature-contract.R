testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")
testthat::skip_if_not_installed("readr")

targeted_repo_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
}

source(file.path(
  targeted_repo_root(), "R", "clusterprofiler_reproducibility.R"
))
source(file.path(targeted_repo_root(), "R", "microglia_targeted_signature_utils.R"))

testthat::test_that("stochastic fallback seeds are stable and semantic", {
  first <- targeted_enrichment_reproducibility(
    "microglia", "CA1sus_CA1con", "fgsea", 20260824L, 100000L
  )
  second <- targeted_enrichment_reproducibility(
    "microglia", "CA1sus_CA1con", "fgsea", 20260824L, 100000L
  )
  distinct <- targeted_enrichment_reproducibility(
    "microglia", "CA3sus_CA3con", "fgsea", 20260824L, 100000L
  )
  deterministic <- targeted_enrichment_reproducibility(
    "microglia", "CA1sus_CA1con", "limma_ranked_geneSetTest",
    20260824L, 100000L
  )
  testthat::expect_identical(first, second)
  testthat::expect_false(identical(first$gsea_seed, distinct$gsea_seed))
  testthat::expect_identical(first$fgsea_simple_nperm, 1000L)
  testthat::expect_true(is.na(deterministic$gsea_seed))
  testthat::expect_identical(
    deterministic$rng_kind, "not_applicable_deterministic_method"
  )
})

testthat::test_that("fallback method selection is explicit in provenance", {
  methods <- data.frame(
    method = c(
      "limma_ranked_geneSetTest", "fgsea", "clusterProfiler_GSEA"
    ),
    available = TRUE, priority = 1:3, stringsAsFactors = FALSE
  )
  reproducibility <- targeted_enrichment_reproducibility(
    "microglia", "CA1sus_CA1con", "fgsea", 20260824L, 100000L
  )
  observed <- targeted_method_selection_provenance(
    methods, "fgsea",
    c("limma_ranked_geneSetTest:failed:test", "fgsea:selected"),
    reproducibility
  )
  testthat::expect_true(observed$stochastic_fallback_used)
  testthat::expect_identical(observed$selected_method_priority, 2L)
  testthat::expect_match(
    observed$method_attempt_log,
    "limma_ranked_geneSetTest:failed:test;fgsea:selected", fixed = TRUE
  )
  testthat::expect_match(
    observed$method_priority_order,
    "limma_ranked_geneSetTest;fgsea;clusterProfiler_GSEA", fixed = TRUE
  )
})

testthat::test_that("seeded fallback scope restores caller RNG", {
  identity <- targeted_enrichment_reproducibility(
    "microglia", "CA1sus_CA1con", "fgsea", 20260824L, 100000L
  )
  set.seed(12345)
  kind_before <- RNGkind()
  state_before <- .Random.seed
  first <- run_with_stable_gsea_rng(
    function() stats::runif(5), identity$gsea_seed
  )
  testthat::expect_identical(RNGkind(), kind_before)
  testthat::expect_identical(.Random.seed, state_before)
  second <- run_with_stable_gsea_rng(
    function() stats::runif(5), identity$gsea_seed
  )
  testthat::expect_identical(first, second)
})

testthat::test_that("effect-statistic provenance is method specific", {
  methods <- c("limma_ranked_geneSetTest", "fgsea", "clusterProfiler_GSEA")
  testthat::expect_equal(
    targeted_effect_statistic_type(methods),
    c("standardized_mean_rank_stat", "NES", "NES")
  )
  testthat::expect_equal(
    targeted_effect_statistic_display_label("standardized_mean_rank_stat"),
    "Standardized mean rank statistic"
  )
  testthat::expect_equal(
    targeted_effect_statistic_panel_label(c("standardized_mean_rank_stat", "standardized_mean_rank_stat")),
    "Standardized mean rank statistic"
  )
})

marker_row <- function(marker_set, gene, evidence = "inferential", valid = TRUE,
                       allowed = TRUE, fdr_pass = TRUE, mapped = gene) {
  tibble::tibble(
    marker_set = marker_set,
    marker_evidence_type = evidence,
    inferential_test_valid = valid,
    claim_allowed = allowed,
    FDR_pass = fdr_pass,
    GeneSymbol = gene,
    mapped_gene_symbol = mapped,
    marker_contract_version = TARGETED_EMPIRICAL_MARKER_CONTRACT
  )
}

testthat::test_that("canonical empirical ROI v2 is the only inferential empirical membership", {
  markers <- dplyr::bind_rows(
    marker_row("empirical_microglia_roi_enriched", "Tmem119"),
    marker_row("empirical_microglia_roi_enriched", "TMEM119"),
    marker_row("empirical_microglia_roi_enriched", "BadFdr", fdr_pass = FALSE),
    marker_row("empirical_microglia_roi_high_confidence", "Nested"),
    marker_row("empirical_microglia_neuropil_shared", NA_character_, "descriptive_only",
               allowed = FALSE, fdr_pass = FALSE, mapped = "Shared1")
  )
  out <- build_canonical_empirical_roi_term2gene(markers)

  testthat::expect_equal(out$microglia$gene, "TMEM119")
  testthat::expect_false("NESTED" %in% out$term2gene$gene)
  testthat::expect_equal(out$shared$gene, "SHARED1")
  testthat::expect_true(all(out$microglia$signature_claim_basis_eligible))
  testthat::expect_false(any(out$shared$signature_claim_basis_eligible))
  testthat::expect_match(out$shared$signature_membership_evidence, "not_equivalence")
  testthat::expect_equal(out$diagnostics$detail[out$diagnostics$check == "stress_rank_derived_inferential_genes"], "0")
})

testthat::test_that("canonical empirical marker contract fails closed", {
  markers <- marker_row("empirical_microglia_roi_enriched", "Tmem119")
  markers$marker_contract_version <- "legacy_invalid"
  testthat::expect_error(build_canonical_empirical_roi_term2gene(markers), "contract mismatch")
})

micro_row <- function(left_region = "CA1", right_region = left_region,
                      left_condition = "sus", right_condition = "con", NES = 1) {
  tibble::tibble(
    signature = "sig", comparison = paste0(left_region, left_condition, "_", right_region, right_condition),
    left_region = left_region, right_region = right_region,
    left_condition = left_condition, right_condition = right_condition,
    NES = NES,
    effect_statistic = NES,
    effect_statistic_type = "standardized_mean_rank_stat",
    enrichment_method = "limma_ranked_geneSetTest",
    fdr_scope = "within_comparison_signature_panel"
  )
}

neuropil_row <- function(comparison, left_region, right_region = left_region,
                         left_condition = "sus", right_condition = "con",
                         layer = "slm", NES = 1, padj = 0.01) {
  tibble::tibble(
    signature = "sig", comparison = comparison,
    left_region = left_region, right_region = right_region,
    left_condition = left_condition, right_condition = right_condition,
    layer = layer, NES = NES, effect_statistic = NES,
    effect_statistic_type = "standardized_mean_rank_stat",
    enrichment_method = "limma_ranked_geneSetTest",
    fdr_scope = "within_comparison_signature_panel",
    padj = padj
  )
}

testthat::test_that("claim reference rejects cross-region and cross-contrast fallback", {
  micro <- micro_row("CA1", left_condition = "sus", right_condition = "con")
  neuropil <- dplyr::bind_rows(
    neuropil_row("CA3slmsus_CA3slmcon", "CA3"),
    neuropil_row("CA1slmres_CA1slmcon", "CA1", left_condition = "res", right_condition = "con")
  )
  out <- targeted_attach_neuropil_reference(micro, neuropil)

  testthat::expect_equal(out$reference_match_type, "missing_exact_neuropil_reference")
  testthat::expect_equal(out$matched_neuropil_n_layers, 0L)
  testthat::expect_true(is.na(out$neuropil_reference_NES))
  testthat::expect_true(out$old_global_fallback_would_have_been_used)
  testthat::expect_true(is.finite(out$global_neuropil_NES))
  testthat::expect_false(out$matched_reference_uses_global_fallback)
})

testthat::test_that("all exact-route layers are retained without minimum-FDR selection", {
  micro <- micro_row()
  neuropil <- dplyr::bind_rows(
    neuropil_row("CA1slmsus_CA1slmcon", "CA1", layer = "slm", NES = 0.5, padj = 0.001),
    neuropil_row("CA1srsus_CA1srcon", "CA1", layer = "sr", NES = 1.5, padj = 0.20),
    neuropil_row("CA3slmsus_CA3slmcon", "CA3", layer = "slm", NES = 9, padj = 1e-10)
  )
  out <- targeted_attach_neuropil_reference(micro, neuropil)

  testthat::expect_equal(out$matched_neuropil_n_layers, 2L)
  testthat::expect_equal(out$matched_neuropil_median_NES, 1)
  testthat::expect_equal(out$matched_neuropil_median_effect_statistic, 1)
  testthat::expect_equal(out$neuropil_reference_effect_statistic, 1)
  testthat::expect_equal(out$neuropil_reference_effect_statistic_type, "standardized_mean_rank_stat")
  testthat::expect_equal(out$neuropil_reference_enrichment_method, "limma_ranked_geneSetTest")
  testthat::expect_equal(out$neuropil_reference_NES, 1)
  testthat::expect_true(is.na(out$neuropil_reference_padj))
  testthat::expect_equal(out$matched_neuropil_min_padj, 0.001)
  testthat::expect_equal(out$matched_neuropil_route_keys, "CA1|CA1|sus|con")
  testthat::expect_true(out$matched_neuropil_any_sig_same_direction)
  testthat::expect_false(out$matched_neuropil_any_sig_opposite_direction)
  testthat::expect_match(out$matched_neuropil_layer_results, "CA1slmsus_CA1slmcon")
  testthat::expect_match(out$matched_neuropil_layer_results, "CA1srsus_CA1srcon")
  testthat::expect_match(out$matched_neuropil_layer_results, "effect_statistic=")
  testthat::expect_false(grepl(";NES=", out$matched_neuropil_layer_results, fixed = TRUE))
  testthat::expect_false(grepl("CA3", out$matched_neuropil_layer_results))

  reversed <- targeted_attach_neuropil_reference(micro, dplyr::slice(neuropil, dplyr::n():1))
  testthat::expect_equal(out, reversed)
})

testthat::test_that("same-direction and opposite-direction exact layers are explicit", {
  neuropil <- dplyr::bind_rows(
    neuropil_row("CA1slmsus_CA1slmcon", "CA1", NES = 0.8, padj = 0.01),
    neuropil_row("CA1srsus_CA1srcon", "CA1", NES = -0.7, padj = 0.02)
  )
  out <- targeted_attach_neuropil_reference(micro_row(NES = 1), neuropil)
  testthat::expect_true(out$matched_neuropil_any_sig_same_direction)
  testthat::expect_true(out$matched_neuropil_any_sig_opposite_direction)
})

testthat::test_that("classification uses significant exact-layer evidence and fails non-significant rows closed", {
  out <- targeted_classify_signature_evidence(
    microglia_significant = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    matched_sig_same_direction = c(TRUE, FALSE, FALSE, TRUE, FALSE),
    matched_sig_opposite_direction = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    empirical_support = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    reference_support = c(FALSE, FALSE, FALSE, FALSE, FALSE),
    curated_support = c(FALSE, FALSE, FALSE, FALSE, FALSE)
  )
  testthat::expect_equal(
    out,
    c("neuropil_shared", "mixed_microenvironment", "microglia_enriched_empirical", "ambiguous", "ambiguous")
  )
})

testthat::test_that("formal stress contrasts match independently", {
  micro <- dplyr::bind_rows(
    micro_row(left_condition = "res", right_condition = "con"),
    micro_row(left_condition = "sus", right_condition = "con"),
    micro_row(left_condition = "sus", right_condition = "res")
  )
  neuropil <- dplyr::bind_rows(
    neuropil_row("CA1slmres_CA1slmcon", "CA1", left_condition = "res", right_condition = "con", NES = 0.1),
    neuropil_row("CA1slmsus_CA1slmcon", "CA1", left_condition = "sus", right_condition = "con", NES = 0.2),
    neuropil_row("CA1slmsus_CA1slmres", "CA1", left_condition = "sus", right_condition = "res", NES = 0.3)
  )
  out <- targeted_attach_neuropil_reference(micro, neuropil)
  testthat::expect_equal(out$matched_neuropil_n_layers, rep(1L, 3))
  testthat::expect_equal(out$matched_neuropil_median_NES, c(0.1, 0.2, 0.3))
})

claim_rows <- function() {
  tibble::tibble(
    padj = c(0.01, 0.01, 0.01, 0.01, 0.01),
    contrast_class = rep("within_region_condition", 5),
    reference_celltype_support = rep("microglia_supported", 5),
    microglia_signature_class = c(
      "microglia_enriched_empirical", "curated_microglia_program",
      "microglia_enriched_empirical", "microglia_enriched_reference_supported",
      "microglia_enriched_reference_supported"
    ),
    signature_source = c(
      "canonical_empirical_roi_v2", "curated", "descriptive_stress_rank_diagnostic_only",
      "reference_atlas_EWCE", "reference_atlas_EWCE"
    ),
    signature_membership_contract = c(
      TARGETED_EMPIRICAL_MARKER_CONTRACT, "curated_manual_v1", "noninferential_stress_rank_diagnostic_v1",
      "ewceData_reference_atlas_v1", "ewceData_reference_atlas_v1"
    ),
    signature_claim_basis_eligible = c(TRUE, FALSE, FALSE, TRUE, TRUE),
    reference_match_type = c(
      "exact_same_region_formal_contrast", "exact_same_region_formal_contrast",
      "exact_same_region_formal_contrast", "global_neuropil", "missing_exact_neuropil_reference"
    ),
    matched_reference_uses_global_fallback = c(FALSE, FALSE, FALSE, TRUE, FALSE)
  )
}

testthat::test_that("one fail-closed claim-ready contract excludes circular, curated, and global evidence", {
  out <- targeted_add_claim_ready(claim_rows())
  testthat::expect_equal(out$claim_ready, c(TRUE, FALSE, FALSE, FALSE, TRUE))
  testthat::expect_match(out$claim_ready_failure_reason[[2]], "signature_class_not_claim_eligible")
  testthat::expect_match(out$claim_ready_failure_reason[[3]], "empirical_membership_not_canonical_roi_v2")
  testthat::expect_match(out$claim_ready_failure_reason[[4]], "global_neuropil_fallback_used")
  testthat::expect_true(all(out$claim_ready_contract_version == TARGETED_CLAIM_READY_CONTRACT))
})

testthat::test_that("non-significant microglia enrichment fails claim-ready", {
  x <- claim_rows()[1, ]
  x$padj <- NA_real_
  out <- targeted_add_claim_ready(x)
  testthat::expect_false(out$claim_ready)
  testthat::expect_match(out$claim_ready_failure_reason, "microglia_not_fdr_significant")
})

testthat::test_that("active script cannot bind descriptive stress-rank membership into TERM2GENE", {
  script <- paste(readLines(file.path(targeted_repo_root(), "04_differential_expression_enrichment", "05_microglia_targeted_signature_enrichment.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "bind_rows\\(curated_term2gene, canonical_empirical_term2gene, reference_term2gene\\)")
  testthat::expect_false(grepl("bind_rows\\([^\\n]*descriptive_rank_diagnostics\\$term2gene", script))
  testthat::expect_match(script, "stress_rank_derived_inferential_terms")
  testthat::expect_match(script, "fdr_scope = \"within_comparison_signature_panel\"")
  testthat::expect_match(script, "padj_global_signature_contrast")
  testthat::expect_match(script, "effect_statistic_type = targeted_effect_statistic_type")
  testthat::expect_false(grepl("Microglia ROI NES", script, fixed = TRUE))
  testthat::expect_false(grepl("Neuropil reference NES", script, fixed = TRUE))
  testthat::expect_match(script, "targeted_add_claim_ready\\(\\)")
  testthat::expect_match(script, "filter\\(\\.data\\$claim_ready\\)")
  testthat::expect_false(grepl("claim_ready = \\.data\\$significant", script))
  testthat::expect_match(
    script,
    'unset = "limma_ranked_geneSetTest,fgsea,clusterProfiler_GSEA"',
    fixed = TRUE
  )
  testthat::expect_match(script, "run_with_stable_gsea_rng\\(")
  testthat::expect_match(script, "run_seeded_clusterprofiler_gsea\\(")
  testthat::expect_match(script, "targeted_method_selection_provenance\\(")
  testthat::expect_false(grepl(
    "tryCatch\\(clusterProfiler::GSEA\\(", script, perl = TRUE
  ))
})
