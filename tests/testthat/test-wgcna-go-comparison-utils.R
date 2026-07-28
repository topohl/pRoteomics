source(testthat::test_path("..", "..", "R", "wgcna_go_comparison_utils.R"))

testthat::test_that("module GO heatmap data retains every module and scores BH-FDR", {
  enrichment <- data.frame(
    ModuleID = c("WGCNA_m01", "WGCNA_m01", "WGCNA_m02"),
    ModuleProteinSetType = "all", Ontology = "BP",
    ID = c("GO:1", "GO:2", "GO:1"),
    Description = c("synapse organization", "vesicle organization", "synapse organization"),
    p.adjust = c(0.001, 0.2, 0.01),
    GeneRatio = c("4/20", "3/20", "5/25"), Count = c(4L, 3L, 5L),
    MappedModuleSize = c(20L, 20L, 25L), ModuleSize = c(21L, 21L, 26L),
    stringsAsFactors = FALSE
  )
  modules <- data.frame(ModuleID = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"), stringsAsFactors = FALSE)
  out <- make_module_go_heatmap_data(enrichment, modules, max_terms = 5L)
  testthat::expect_equal(out$status, "ok")
  testthat::expect_setequal(unique(out$matrix$ModuleID), modules$ModuleID)
  testthat::expect_true(all(out$matrix$EnrichmentScore[out$matrix$p_adjust >= 0.05] == 0))
  testthat::expect_true(any(out$matrix$ModuleID == "WGCNA_m03" & out$matrix$EnrichmentScore == 0))
  testthat::expect_equal(out$matrix$EnrichmentScore[out$matrix$ModuleID == "WGCNA_m01" & out$matrix$TermID == "GO:1"], 3)
})

testthat::test_that("supermodule values are recurrence-weighted member evidence rather than pooled p-values", {
  module_heatmap <- data.frame(
    ModuleID = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m01", "WGCNA_m02"),
    TermID = c("GO:1", "GO:1", "GO:2", "GO:2"),
    Description = c("synapse", "synapse", "immune", "immune"),
    p_adjust = c(0.01, 1, 0.001, 0.001),
    EnrichmentScore = c(2, 0, 3, 3),
    Significant = c(TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  membership <- data.frame(
    ModuleID = c("WGCNA_m01", "WGCNA_m02"), SupermoduleID = c("SM01", "SM01"),
    SupermoduleDisplayLabel = "SM01 · example", stringsAsFactors = FALSE
  )
  out <- make_supermodule_go_heatmap_data(module_heatmap, membership, max_terms = 5L)
  synapse <- out$matrix[out$matrix$TermID == "GO:1", , drop = FALSE]
  immune <- out$matrix[out$matrix$TermID == "GO:2", , drop = FALSE]
  testthat::expect_equal(synapse$fraction_member_modules_FDR_significant, 0.5)
  testthat::expect_equal(synapse$EnrichmentScore, 1)
  testthat::expect_equal(immune$fraction_member_modules_FDR_significant, 1)
  testthat::expect_equal(immune$EnrichmentScore, 3)
  testthat::expect_match(unique(out$matrix$enrichment_measure), "mean capped")
})

testthat::test_that("supermodule term selection uses the full module GO table before module display trimming", {
  enrichment <- data.frame(
    ModuleID = rep(c("WGCNA_m01", "WGCNA_m02"), each = 2L),
    ModuleProteinSetType = "all", Ontology = "BP",
    ID = c("GO:unique_1", "GO:shared", "GO:unique_2", "GO:shared"),
    Description = c("unique module one", "recurring theme", "unique module two", "recurring theme"),
    p.adjust = c(1e-6, 1e-3, 1e-6, 1e-3), stringsAsFactors = FALSE
  )
  modules <- data.frame(ModuleID = c("WGCNA_m01", "WGCNA_m02"), stringsAsFactors = FALSE)
  membership <- data.frame(ModuleID = modules$ModuleID, SupermoduleID = "SM01", stringsAsFactors = FALSE)
  module_panel <- make_module_go_heatmap_data(enrichment, modules, top_terms_per_module = 1L, max_terms = 5L)
  full_terms <- make_full_module_go_term_matrix(enrichment, modules)
  super_panel <- make_supermodule_go_heatmap_data(full_terms$matrix, membership, top_terms_per_supermodule = 1L, max_terms = 5L)
  testthat::expect_false("GO:shared" %in% module_panel$matrix$TermID)
  testthat::expect_true("GO:shared" %in% super_panel$matrix$TermID)
})

testthat::test_that("representative selection falls back cleanly when no term passes FDR", {
  x <- data.frame(
    EntityID = c("M1", "M2"), TermID = c("GO:1", "GO:2"), Description = c("a", "b"),
    EnrichmentScore = c(0, 0), Significant = c(FALSE, FALSE), stringsAsFactors = FALSE
  )
  out <- select_representative_go_terms(x, max_terms = 5L)
  testthat::expect_equal(unique(out$selection_basis), "lowest_FDR_fallback")
  testthat::expect_equal(nrow(out), 2L)
})
