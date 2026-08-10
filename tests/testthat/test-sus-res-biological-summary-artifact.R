source(testthat::test_path("..", "..", "R", "paths.R"))

sus_res_summary_path <- repo_path(
  "results", "tables", "04_differential_expression_enrichment",
  "sus_res_spatial_dap_atlas", "global", "sus_res_manuscript_theme_summary.csv"
)
sus_res_supported_audit_path <- repo_path(
  "results", "source_data", "04_differential_expression_enrichment",
  "compareGO_spatial_atlas", "sus_res_supported_go_term_theme_audit.csv"
)
sus_res_dap_counts_path <- repo_path(
  "results", "tables", "04_differential_expression_enrichment",
  "sus_res_spatial_dap_atlas", "global", "sus_res_dap_counts.csv"
)
sus_res_panel_path <- repo_path(
  "results", "source_data", "04_differential_expression_enrichment",
  "sus_res_spatial_dap_atlas", "global", "panel_c_sus_res_ranked_GSEA_themes.csv"
)

testthat::test_that("live SUS-RES inspection rows preserve exact ontology-mapped representative GSEA terms", {
  required <- c(sus_res_summary_path, sus_res_supported_audit_path)
  testthat::skip_if(any(!file.exists(required)), "Generated ontology-aware SUS-RES summary artifacts are unavailable")
  summary <- utils::read.csv(sus_res_summary_path, check.names = FALSE, stringsAsFactors = FALSE)
  audit <- utils::read.csv(sus_res_supported_audit_path, check.names = FALSE, stringsAsFactors = FALSE)

  testthat::expect_gt(nrow(summary), 0L)
  testthat::expect_equal(nrow(summary), 144L)
  testthat::expect_true(all(summary$theme_role %in% c("primary", "qc_review")))
  testthat::expect_true(all(summary$mapping_method == "go_id_ontology"))
  testthat::expect_true(all(summary$evidence_source_family == "ranked_GSEA"))
  testthat::expect_false(anyDuplicated(paste(summary$dataset, summary$spatial_unit, summary$manuscript_theme_id, sep = "|")) > 0L)
  testthat::expect_true(all(summary$n_theme_terms_tested > 0L))
  testthat::expect_true(all(is.finite(summary$median_NES_all_theme_terms)))
  supported_rows <- as.logical(summary$FDR_support_present)
  testthat::expect_true(all(is.finite(summary$representative_NES[supported_rows])))
  testthat::expect_true(all(is.finite(summary$representative_FDR[supported_rows])))
  testthat::expect_true(all(summary$representative_FDR[supported_rows] < 0.05))
  testthat::expect_true(all(is.na(summary$representative_NES[!supported_rows])))
  testthat::expect_true(all(is.na(summary$representative_FDR[!supported_rows])))

  for (i in which(supported_rows)) {
    theme_pattern <- paste0("(^|;)", summary$manuscript_theme_id[[i]], "($|;)")
    candidate <- audit[
      audit$dataset == summary$dataset[[i]] &
        audit$spatial_unit == summary$spatial_unit[[i]] &
        grepl(theme_pattern, audit$manuscript_theme_ids),
      , drop = FALSE
    ]
    candidate <- candidate[order(
      candidate$p.adjust, -abs(candidate$NES), candidate$GO_ID, candidate$GO_description,
      method = "radix"
    ), , drop = FALSE]
    testthat::expect_gt(nrow(candidate), 0L)
    testthat::expect_identical(summary$representative_GO_ID[[i]], candidate$GO_ID[[1]])
    testthat::expect_identical(summary$representative_GO_term[[i]], candidate$GO_description[[1]])
    testthat::expect_equal(summary$representative_NES[[i]], candidate$NES[[1]], tolerance = 1e-12)
    testthat::expect_equal(summary$representative_FDR[[i]], candidate$p.adjust[[1]], tolerance = 1e-12)
    testthat::expect_identical(summary$leading_edge_proteins[[i]], candidate$leading_edge_proteins[[1]])
    testthat::expect_identical(summary$leading_edge_genes[[i]], candidate$leading_edge_genes[[1]])
  }
})

testthat::test_that("every supported occurrence is retained with explicit assignment status", {
  testthat::skip_if_not(file.exists(sus_res_supported_audit_path), "Generated supported-term audit is unavailable")
  audit <- utils::read.csv(sus_res_supported_audit_path, check.names = FALSE, stringsAsFactors = FALSE)
  testthat::expect_equal(nrow(audit), 336L)
  testthat::expect_equal(length(unique(audit$GO_ID)), 167L)
  testthat::expect_true(all(audit$assignment_status %in% c("single_theme", "multi_theme", "unclassified", "qc_review")))
  testthat::expect_true(all(c("GO:0010605", "GO:0019219", "GO:0009892", "GO:0140694", "GO:1902600") %in% audit$GO_ID))
  false_rows <- audit[audit$GO_ID %in% c("GO:0010605", "GO:0019219", "GO:0009892", "GO:0140694", "GO:1902600"), , drop = FALSE]
  testthat::expect_true(all(false_rows$assignment_status == "unclassified"))
})

testthat::test_that("live DAP joins and directional recurrence remain descriptive", {
  required <- c(sus_res_summary_path, sus_res_dap_counts_path)
  testthat::skip_if(any(!file.exists(required)), "Generated ontology-aware SUS-RES summary artifacts are unavailable")
  summary <- utils::read.csv(sus_res_summary_path, check.names = FALSE, stringsAsFactors = FALSE)
  dap <- utils::read.csv(sus_res_dap_counts_path, check.names = FALSE, stringsAsFactors = FALSE)
  dap_unit <- ifelse(summary$dataset == "neuron_neuropil", summary$spatial_unit, summary$region)
  idx <- match(paste(summary$dataset, dap_unit, sep = "|"), paste(dap$dataset, dap$spatial_unit, sep = "|"))
  testthat::expect_false(anyNA(idx))
  testthat::expect_identical(as.integer(summary$n_FDR_supported_SUS_RES_DAPs), as.integer(dap$n_DAP_FDR05[idx]))

  for (theme in unique(summary$manuscript_theme_id)) {
    rows <- summary$manuscript_theme_id == theme
    supported <- rows & as.logical(summary$FDR_support_present)
    expected_units <- sum(supported)
    expected_datasets <- length(unique(summary$dataset[supported]))
    testthat::expect_true(all(summary$sus_res_recurrent_units[rows] == expected_units))
    testthat::expect_true(all(summary$sus_res_recurrent_datasets[rows] == expected_datasets))
    if (expected_units > 1L && any(summary$direction_consistency[supported] == "mixed_direction")) {
      testthat::expect_true(all(summary$directional_recurrence[rows] == "recurrent_mixed_direction"))
    }
  }
})

testthat::test_that("Panel C contains primary themes only and retains source family", {
  testthat::skip_if_not(file.exists(sus_res_panel_path), "Generated ontology-aware Panel C source is unavailable")
  panel <- utils::read.csv(sus_res_panel_path, check.names = FALSE, stringsAsFactors = FALSE)
  testthat::expect_equal(nrow(panel), 108L)
  testthat::expect_true(all(panel$theme_role == "primary"))
  testthat::expect_true(all(panel$evidence_source_family == "ranked_GSEA"))
  testthat::expect_false(any(panel$theme_role == "qc_review"))
  testthat::expect_true(any(panel$FDR_support_present))
  testthat::expect_true(any(!panel$FDR_support_present))
  testthat::expect_true(all(is.finite(panel$median_NES_all_theme_terms)))
})
