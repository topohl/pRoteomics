testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")
source(testthat::test_path("..", "..", "R", "protein_mapping_utils.R"))
source(testthat::test_path("..", "..", "R", "module_contracts.R"))

cache_parity_fixture <- function() {
  metadata <- data.frame(
    WGCNAInternalColor = c("turquoise", "blue", "brown", "yellow", "green"),
    ModuleID = sprintf("WGCNA_m%02d", 1:5),
    ModuleLegacyID = paste0("WGCNA_", c("turquoise", "blue", "brown", "yellow", "green")),
    ModuleColor = c("#486A8A", "#7FA6C1", "#5D7894", "#8AA0AF", "#2F6F73"),
    ModuleColorName = c("deep_steel_blue", "muted_sky_blue", "slate_blue", "blue_grey", "deep_teal"),
    ModuleColorLabel = c("deep steel blue", "muted sky blue", "slate blue", "blue grey", "deep teal"),
    stringsAsFactors = FALSE
  )
  member_map <- data.frame(
    dataset = "microglia",
    module_eigengene = paste0("ME", metadata$WGCNAInternalColor),
    WGCNAInternalColor = metadata$WGCNAInternalColor,
    ModuleID = metadata$ModuleID,
    ModuleColor = metadata$ModuleColor,
    SupermoduleID = c("SM01", "SM02", "SM02", "SM03", "SM03"),
    Supermodule_DisplayLabel = c(
      "SM01 - singleton", "SM02 - recurring", "SM02 - recurring",
      "SM03 - unresolved", "SM03 - unresolved"
    ),
    stringsAsFactors = FALSE
  )
  go <- dplyr::bind_rows(
    data.frame(ModuleID = "WGCNA_m01", Description = "singleton process", p.adjust = 0.01),
    data.frame(ModuleID = c("WGCNA_m02", "WGCNA_m03"), Description = "shared process", p.adjust = c(0.01, 0.02)),
    data.frame(ModuleID = "WGCNA_m04", Description = "distinct process A", p.adjust = 0.01),
    data.frame(ModuleID = "WGCNA_m05", Description = "distinct process B", p.adjust = 0.01)
  ) |>
    dplyr::left_join(
      metadata[, c("ModuleID", "WGCNAInternalColor", "ModuleColor")],
      by = "ModuleID",
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      ModuleProteinSetType = "all",
      Ontology = "BP",
      .before = 1
    )
  protein_ids <- paste0("PG:microglia:", seq_len(5))
  expression_data <- matrix(
    seq_len(20), nrow = 4,
    dimnames = list(paste0("sample", seq_len(4)), protein_ids)
  )
  merged_colors <- stats::setNames(metadata$WGCNAInternalColor, protein_ids)
  merged_mes <- matrix(
    seq_len(20) / 20, nrow = 4,
    dimnames = list(rownames(expression_data), member_map$module_eigengene)
  )
  state <- list(
    feature_key_contract_version = wgcna_feature_key_contract_version(),
    feature_key_fingerprint = wgcna_feature_key_fingerprint(protein_ids),
    expression.data = expression_data,
    wgcna_feature_table = data.frame(ProteinGroupID = protein_ids),
    wgcna_member_bridge = data.frame(ProteinGroupID = protein_ids),
    WGCNA_feature_universe = data.frame(ProteinGroupID = protein_ids),
    sample_info = data.frame(condition = rep(c("con", "res"), each = 2), row.names = rownames(expression_data)),
    mergedColors = merged_colors,
    mergedMEs = merged_mes,
    kME = data.frame(ProteinGroupID = protein_ids),
    WGCNA_modules_long = data.frame(
      ProteinGroupID = protein_ids,
      ModuleID = metadata$ModuleID,
      stringsAsFactors = FALSE
    ),
    module_summary = data.frame(ModuleID = metadata$ModuleID),
    GO_enrichment = go,
    GO_enrichment_QC = metadata[, c("ModuleID", "WGCNAInternalColor", "ModuleColor")] |>
      dplyr::mutate(ModuleProteinSetType = "all", Ontology = "BP", status = "ok"),
    module_name_map = stats::setNames(paste("Module", metadata$ModuleID), metadata$WGCNAInternalColor),
    module_label_table = metadata,
    color_to_MEcol = stats::setNames(member_map$module_eigengene, metadata$WGCNAInternalColor),
    ME_names_stable = member_map$module_eigengene,
    module_preservation = data.frame(ModuleID = metadata$ModuleID),
    geneTree = list(method = "synthetic"),
    softPower = 6,
    parameters = list(dataset_profile = "microglia")
  )
  list(state = state, go = go, member_map = member_map, metadata = metadata)
}

testthat::test_that("cached hydration preserves live supermodule GO evidence", {
  fixture <- cache_parity_fixture()
  live <- wgcna_supermodule_go_evidence_summary(fixture$go, fixture$member_map)

  hydrated <- wgcna_hydrate_cached_state(fixture$state)
  resumed <- wgcna_supermodule_go_evidence_summary(
    hydrated$GO_enrichment_long,
    fixture$member_map
  )

  parity_columns <- c(
    "SupermoduleID", "member_ModuleIDs", "n_recurring_significant_GO_terms",
    "n_modules_with_GO_support", "GO_label_confidence_class",
    "Supermodule_NameSource", "ManualReviewRequired"
  )
  testthat::expect_identical(resumed[, parity_columns], live[, parity_columns])
  testthat::expect_identical(hydrated$GO_enrichment_long, fixture$state$GO_enrichment)
  testthat::expect_identical(hydrated$GO_enrichment_QC, fixture$state$GO_enrichment_QC)
  testthat::expect_identical(hydrated$WGCNA_module_summary, fixture$state$module_summary)
  testthat::expect_identical(hydrated$module_preservation_long, fixture$state$module_preservation)
  testthat::expect_identical(hydrated$module_color_metadata, fixture$metadata)

  recurring <- resumed[resumed$SupermoduleID == "SM02", ]
  testthat::expect_equal(recurring$n_recurring_significant_GO_terms, 1L)
  testthat::expect_equal(recurring$n_modules_with_GO_support, 2L)
  testthat::expect_equal(recurring$GO_label_confidence_class, "high")
  testthat::expect_equal(recurring$Supermodule_NameSource, "recurring_significant_GO")
  testthat::expect_false(recurring$ManualReviewRequired)

  singleton <- resumed[resumed$SupermoduleID == "SM01", ]
  unresolved <- resumed[resumed$SupermoduleID == "SM03", ]
  testthat::expect_equal(singleton$GO_label_confidence_class, "singleton")
  testthat::expect_true(singleton$ManualReviewRequired)
  testthat::expect_equal(unresolved$n_recurring_significant_GO_terms, 0L)
  testthat::expect_equal(unresolved$GO_label_confidence_class, "mixed_or_unresolved")
  testthat::expect_true(unresolved$ManualReviewRequired)
})

testthat::test_that("a checkpoint reporting GO evidence cannot hydrate GO to NULL", {
  fixture <- cache_parity_fixture()
  fixture$state$GO_enrichment <- NULL
  testthat::expect_error(
    wgcna_hydrate_cached_state(fixture$state),
    "GO enrichment evidence but hydrated GO_enrichment_long is NULL",
    fixed = TRUE
  )
})

testthat::test_that("hydrated GO must agree with stable module metadata", {
  fixture <- cache_parity_fixture()
  fixture$state$GO_enrichment$ModuleID[[1]] <- "turquoise"
  testthat::expect_error(
    wgcna_hydrate_cached_state(fixture$state),
    "stable ModuleID values"
  )

  fixture <- cache_parity_fixture()
  fixture$state$GO_enrichment$ModuleColor[[1]] <- "#FFFFFF"
  testthat::expect_error(
    wgcna_hydrate_cached_state(fixture$state),
    "ModuleColor disagrees"
  )
})
