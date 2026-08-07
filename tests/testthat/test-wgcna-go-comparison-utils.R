source(testthat::test_path("..", "..", "R", "wgcna_go_comparison_utils.R"))

testthat::test_that("01b CLI defaults to three focused terms and supports the canonical all compositor", {
  script <- paste(readLines(
    testthat::test_path("..", "..", "06_modules_WGCNA", "01b_module_supermodule_GO_heatmaps.R"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_true(grepl(
    'arg_value("--focused-terms-per-supermodule", "3")', script, fixed = TRUE
  ))
  testthat::expect_true(grepl("wgcna_cli(allow_all = TRUE)", script, fixed = TRUE))
})

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
  focused_panel <- make_focused_supermodule_go_comparison(full_terms$matrix, membership, focused_terms_per_supermodule = 1L)
  testthat::expect_false("GO:shared" %in% module_panel$matrix$TermID)
  testthat::expect_true("GO:shared" %in% super_panel$matrix$TermID)
  testthat::expect_equal(focused_panel$selection$TermID, "GO:shared")
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

focused_go_fixture <- function() {
  membership <- data.frame(
    ModuleID = paste0("M", 1:5),
    SupermoduleID = c("SM01", "SM01", "SM02", "SM03", "SM04"),
    stringsAsFactors = FALSE
  )
  terms <- data.frame(
    TermID = paste0("GO:", 1:5),
    Description = c("recurrent moderate", "single strong", "recurrent strong", "third representative", "unsupported"),
    stringsAsFactors = FALSE
  )
  x <- merge(
    expand.grid(ModuleID = membership$ModuleID, TermID = terms$TermID, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
    terms, by = "TermID", sort = FALSE
  )
  x$p_adjust <- 1
  x$EnrichmentScore <- 0
  set_evidence <- function(module, term, fdr, score) {
    hit <- x$ModuleID == module & x$TermID == term
    x$p_adjust[hit] <<- fdr
    x$EnrichmentScore[hit] <<- score
  }
  set_evidence("M1", "GO:1", 0.01, 2)
  set_evidence("M2", "GO:1", 0.01, 2)
  set_evidence("M1", "GO:2", 1e-05, 5)
  set_evidence("M1", "GO:3", 0.001, 3)
  set_evidence("M2", "GO:3", 0.001, 3)
  set_evidence("M1", "GO:4", 1e-06, 6)
  set_evidence("M3", "GO:2", 0.01, 2)
  set_evidence("M5", "GO:2", 0.001, 3)
  set_evidence("M5", "GO:1", 0.04, 1.4)
  x$Significant <- x$p_adjust <= 0.05
  list(module_go = x, membership = membership)
}

testthat::test_that("focused default caps at three and never pads 0/1/2-term supermodules", {
  fixture <- focused_go_fixture()
  out <- make_focused_supermodule_go_comparison(fixture$module_go, fixture$membership)
  counts <- table(factor(out$selection$selected_for_supermodule, levels = c("SM01", "SM02", "SM03", "SM04")))

  testthat::expect_equal(out$status, "ok")
  testthat::expect_true(all(counts <= 3L))
  testthat::expect_equal(as.integer(counts), c(3L, 1L, 0L, 2L))
  testthat::expect_equal(out$selection$TermID[out$selection$selected_for_supermodule == "SM01"], c("GO:3", "GO:1", "GO:4"))
  testthat::expect_false("GO:2" %in% out$selection$TermID[out$selection$selected_for_supermodule == "SM01"])
  testthat::expect_equal(out$selection$TermID[out$selection$selected_for_supermodule == "SM02"], "GO:2")
  testthat::expect_false("SM03" %in% out$selection$selected_for_supermodule)
  testthat::expect_true(all(grepl("no_nonsignificant_fallback", out$selection$selection_basis, fixed = TRUE)))
})

testthat::test_that("focused ranking remains recurrence-first and deterministic", {
  fixture <- focused_go_fixture()
  out <- make_focused_supermodule_go_comparison(
    fixture$module_go, fixture$membership,
    focused_terms_per_supermodule = 2L
  )
  testthat::expect_equal(
    out$selection$TermID[out$selection$selected_for_supermodule == "SM01"],
    c("GO:3", "GO:1")
  )
})

testthat::test_that("focused union is crossed with every supermodule and retains selector provenance", {
  fixture <- focused_go_fixture()
  out <- make_focused_supermodule_go_comparison(fixture$module_go, fixture$membership)
  n_terms <- length(unique(out$matrix$TermID))
  n_supermodules <- length(unique(fixture$membership$SupermoduleID))
  pairs <- paste(out$matrix$TermID, out$matrix$SupermoduleID, sep = "\r")
  shared <- out$matrix[out$matrix$TermID == "GO:2", , drop = FALSE]

  testthat::expect_equal(nrow(out$matrix), n_terms * n_supermodules)
  testthat::expect_false(anyDuplicated(pairs) > 0L)
  testthat::expect_equal(unique(shared$n_selecting_supermodules), 2L)
  testthat::expect_equal(unique(shared$selecting_supermodules), "SM02;SM04")
  testthat::expect_equal(sum(out$selection$TermID == "GO:2"), 2L)
  testthat::expect_equal(unique(out$matrix$dot_colour_measure), "mean_member_module_enrichment_score")
  testthat::expect_equal(unique(out$matrix$dot_size_measure), "fraction_member_modules_FDR_significant")
  testthat::expect_equal(out$matrix$display_supported_dot, out$matrix$n_modules_FDR_significant > 0L)
})

testthat::test_that("focused summaries expose member evidence without pooled inference", {
  fixture <- focused_go_fixture()
  out <- make_focused_supermodule_go_comparison(fixture$module_go, fixture$membership)
  singleton <- out$matrix[out$matrix$SupermoduleID == "SM02", , drop = FALSE]

  testthat::expect_true(all(singleton$is_singleton_supermodule))
  testthat::expect_true(all(grepl("singleton reflects constituent module", singleton$supermodule_summary_scope, fixed = TRUE)))
  testthat::expect_true(all(grepl("not pooled supermodule ORA", out$matrix$supermodule_summary_scope, fixed = TRUE)))
  testthat::expect_false(any(grepl("pooled|combined", names(out$matrix), ignore.case = TRUE)))
  testthat::expect_false(any(names(out$matrix) %in% c("SupermodulePValue", "SupermoduleFDR", "pooled_p_value", "pooled_FDR")))
})

diverse_focused_go_fixture <- function(include_genes = TRUE) {
  term_ids <- paste0("GO:", LETTERS[1:7])
  x <- data.frame(
    ModuleID = "M1",
    TermID = term_ids,
    Description = c(
      "alpha process", "alpha child", "alpha distant child",
      "unrelated near-identical set", "alpha processes", "distinct fourth", "unsupported"
    ),
    p_adjust = c(1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 0.2),
    EnrichmentScore = c(6, 5.7, 5.5, 5.2, 5, 4.8, 0),
    Significant = c(rep(TRUE, 6L), FALSE),
    stringsAsFactors = FALSE
  )
  if (include_genes) {
    x$geneID <- c(
      "1/2/3/4", "1/2/3/5", "1/6/7/8", "1/2/3/4/9",
      "10/11/12", "20/21", "1/2/3/4"
    )
  }
  hierarchy <- list(
    available = TRUE,
    ancestors = stats::setNames(
      list(character(), "GO:A", "GO:A", character(), character(), character(), character()),
      term_ids
    ),
    ontology = "BP", source = "synthetic", version = "test", reason = NA_character_
  )
  list(
    module_go = x,
    membership = data.frame(ModuleID = "M1", SupermoduleID = "SM01", stringsAsFactors = FALSE),
    hierarchy = hierarchy
  )
}

testthat::test_that("focused diversity selection is sequential, conservative, and capped at three", {
  fixture <- diverse_focused_go_fixture()
  out <- make_focused_supermodule_go_comparison(
    fixture$module_go, fixture$membership, go_hierarchy = fixture$hierarchy
  )
  audit <- out$audit

  testthat::expect_equal(out$redundancy_mode, "GO_hierarchy_plus_gene_overlap")
  testthat::expect_equal(out$selection$TermID, c("GO:A", "GO:C", "GO:E"))
  testthat::expect_equal(out$selection$initial_evidence_rank, c(1L, 3L, 5L))
  testthat::expect_equal(out$selection$focused_representative_rank, 1:3)
  testthat::expect_equal(nrow(out$selection), 3L)
  testthat::expect_false("GO:G" %in% audit$TermID)
  testthat::expect_equal(
    audit$selection_status[audit$TermID == "GO:F"],
    "not_reached_after_three_selected"
  )
  testthat::expect_true(is.na(audit$redundancy_reason[audit$TermID == "GO:F"]))
})

testthat::test_that("ancestor pruning requires substantial overlap when genes are available", {
  fixture <- diverse_focused_go_fixture()
  out <- make_focused_supermodule_go_comparison(
    fixture$module_go, fixture$membership, go_hierarchy = fixture$hierarchy
  )
  high <- out$audit[out$audit$TermID == "GO:B", , drop = FALSE]
  low <- out$audit[out$audit$TermID == "GO:C", , drop = FALSE]

  testthat::expect_equal(high$selection_status, "skipped_redundant_GO_ancestor_descendant")
  testthat::expect_equal(high$redundant_with_TermID, "GO:A")
  testthat::expect_equal(high$go_relationship, "candidate_descendant_of_selected")
  testthat::expect_equal(high$gene_jaccard, 0.6)
  testthat::expect_true(low$selected)
  testthat::expect_lt(wgcna_gene_jaccard("1/6/7/8", "1/2/3/4"), 0.5)
})

testthat::test_that("near-identical genes prune unrelated terms but text similarity alone does not", {
  fixture <- diverse_focused_go_fixture()
  out <- make_focused_supermodule_go_comparison(
    fixture$module_go, fixture$membership, go_hierarchy = fixture$hierarchy
  )
  near_identical <- out$audit[out$audit$TermID == "GO:D", , drop = FALSE]
  text_similar <- out$audit[out$audit$TermID == "GO:E", , drop = FALSE]

  testthat::expect_equal(near_identical$selection_status, "skipped_near_identical_gene_set")
  testthat::expect_equal(near_identical$go_relationship, "none")
  testthat::expect_equal(near_identical$gene_jaccard, 0.8)
  testthat::expect_true(text_similar$selected)
  testthat::expect_equal(text_similar$Description, "alpha processes")
})

testthat::test_that("supporting genes union only significant member-module rows", {
  module_go <- data.frame(
    ModuleID = c("M1", "M2"), TermID = "GO:A", Description = "process",
    p_adjust = c(0.01, 0.2), EnrichmentScore = c(2, 0), Significant = c(TRUE, FALSE),
    geneID = c("1/2", "2/3"), stringsAsFactors = FALSE
  )
  membership <- data.frame(ModuleID = c("M1", "M2"), SupermoduleID = "SM01", stringsAsFactors = FALSE)
  out <- summarize_supermodule_go_evidence(module_go, membership)

  testthat::expect_equal(out$n_supporting_genes, 2L)
  testthat::expect_equal(out$supporting_gene_ids, "1/2")
  testthat::expect_equal(out$supporting_gene_id_delimiter, "/")
})

testthat::test_that("focused fallbacks are explicit and deterministic", {
  no_genes <- diverse_focused_go_fixture(include_genes = FALSE)
  hierarchy_only <- make_focused_supermodule_go_comparison(
    no_genes$module_go, no_genes$membership, go_hierarchy = no_genes$hierarchy
  )
  genes_but_no_hierarchy <- diverse_focused_go_fixture()
  rank_only_a <- make_focused_supermodule_go_comparison(
    genes_but_no_hierarchy$module_go, genes_but_no_hierarchy$membership,
    go_hierarchy = list(available = FALSE, ancestors = list(), reason = "test")
  )
  rank_only_b <- make_focused_supermodule_go_comparison(
    genes_but_no_hierarchy$module_go, genes_but_no_hierarchy$membership,
    go_hierarchy = list(available = FALSE, ancestors = list(), reason = "test")
  )

  testthat::expect_equal(hierarchy_only$redundancy_mode, "GO_hierarchy_only")
  testthat::expect_equal(
    hierarchy_only$audit$selection_status[hierarchy_only$audit$TermID == "GO:B"],
    "skipped_redundant_GO_ancestor_descendant"
  )
  testthat::expect_equal(rank_only_a$redundancy_mode, "evidence_rank_only")
  testthat::expect_equal(rank_only_a$selection$TermID, c("GO:A", "GO:B", "GO:C"))
  testthat::expect_equal(rank_only_a$audit, rank_only_b$audit)
})

testthat::test_that("focused pruning does not alter broad supermodule selection", {
  fixture <- diverse_focused_go_fixture()
  broad <- make_supermodule_go_heatmap_data(
    fixture$module_go, fixture$membership, top_terms_per_supermodule = 7L, max_terms = 20L
  )

  testthat::expect_true(all(c("GO:A", "GO:B", "GO:C", "GO:D", "GO:E", "GO:F") %in% broad$matrix$TermID))
  testthat::expect_false(any(c(
    "n_supporting_genes", "supporting_gene_ids", "supporting_gene_membership_available"
  ) %in% names(broad$matrix)))
})

combined_go_fixture <- function() {
  make_tables <- function(dataset, term_ids, descriptions, singleton = FALSE) {
    selection <- data.frame(
      dataset = dataset, Ontology = "BP", TermID = term_ids, Description = descriptions,
      SupermoduleID = "SM01", selected_for_supermodule = "SM01",
      representative_rank = seq_along(term_ids),
      selection_basis = "FDR_supported_member_module_evidence; max_3_per_supermodule; no_nonsignificant_fallback",
      stringsAsFactors = FALSE
    )
    source <- data.frame(
      dataset = dataset, Ontology = "BP", TermID = term_ids, Description = descriptions,
      SupermoduleID = "SM01", SupermodulePlotLabel = if (singleton) "SM01\u2020" else "SM01",
      row_order = seq_along(term_ids), n_modules_FDR_significant = 1L,
      fraction_member_modules_FDR_significant = seq(0.5, 1, length.out = length(term_ids)),
      mean_member_module_enrichment_score = seq(2, 4, length.out = length(term_ids)),
      is_singleton_supermodule = singleton,
      supermodule_summary_scope = if (singleton) {
        "member-module evidence; singleton reflects constituent module; not pooled supermodule ORA"
      } else {
        "member-module evidence; not pooled supermodule ORA"
      },
      FDRCutoff = 0.05, ScoreCap = 6, display_supported_dot = TRUE,
      dot_colour_measure = "mean_member_module_enrichment_score",
      dot_size_measure = "fraction_member_modules_FDR_significant",
      stringsAsFactors = FALSE
    )
    list(source = source, selection = selection)
  }
  neuropil <- make_tables(
    "neuron_neuropil", c("GO:shared", "GO:similar_A"),
    c("shared process", "same description")
  )
  soma <- make_tables(
    "neuron_soma", c("GO:shared", "GO:similar_B"),
    c("shared process", "same description")
  )
  microglia <- make_tables(
    "microglia", "GO:unique", "ROI-specific process", singleton = TRUE
  )
  list(
    sources = list(
      neuron_neuropil = neuropil$source,
      neuron_soma = soma$source,
      microglia = microglia$source
    ),
    selections = list(
      neuron_neuropil = neuropil$selection,
      neuron_soma = soma$selection,
      microglia = microglia$selection
    )
  )
}

testthat::test_that("combined focused sources retain all datasets and original selections", {
  fixture <- combined_go_fixture()
  original_pairs <- unlist(lapply(fixture$selections, function(x) paste(x$dataset, x$TermID, x$selected_for_supermodule, sep = "\r")))
  out <- combine_focused_supermodule_go_sources(fixture$sources, fixture$selections)
  combined_pairs <- paste(out$selection$dataset, out$selection$TermID, out$selection$selected_for_supermodule, sep = "\r")

  testthat::expect_equal(unique(out$matrix$dataset), wgcna_focused_dataset_order())
  testthat::expect_setequal(combined_pairs, original_pairs)
  testthat::expect_equal(nrow(out$selection), length(original_pairs))
  testthat::expect_equal(out$colour_limits, c(0, 6))
  testthat::expect_equal(out$size_limits, c(0, 1))
  testthat::expect_true(all(out$matrix$combined_colour_measure == "mean_member_module_enrichment_score"))
  testthat::expect_true(all(out$matrix$combined_colour_normalization == "none; shared unnormalized member-module summary scale"))
  testthat::expect_true(all(out$matrix$combined_size_measure == "fraction_member_modules_FDR_significant"))
  testthat::expect_true(all(out$matrix$combined_size_scale_min == 0 & out$matrix$combined_size_scale_max == 1))
})

testthat::test_that("combined recurrence uses exact TermID without semantic merging or inference", {
  fixture <- combined_go_fixture()
  out <- combine_focused_supermodule_go_sources(fixture$sources, fixture$selections)
  shared <- unique(out$selection[out$selection$TermID == "GO:shared", c("n_datasets_selected", "datasets_selected")])
  similar <- unique(out$selection[out$selection$Description == "same description", c("TermID", "n_datasets_selected")])
  singleton <- out$matrix[out$matrix$dataset == "microglia", , drop = FALSE]

  testthat::expect_equal(shared$n_datasets_selected, 2L)
  testthat::expect_equal(shared$datasets_selected, "neuron_neuropil;neuron_soma")
  testthat::expect_setequal(similar$TermID, c("GO:similar_A", "GO:similar_B"))
  testthat::expect_true(all(similar$n_datasets_selected == 1L))
  testthat::expect_true(all(singleton$is_singleton_supermodule))
  testthat::expect_true(all(grepl("singleton reflects constituent module", singleton$supermodule_summary_scope, fixed = TRUE)))
  testthat::expect_true(all(grepl("no cross-dataset inference", out$matrix$cross_dataset_recurrence_scope, fixed = TRUE)))
  testthat::expect_false(any(names(out$matrix) %in% c(
    "cross_dataset_p_value", "cross_dataset_FDR", "combined_p_value", "combined_FDR", "meta_analysis_p_value"
  )))
})
