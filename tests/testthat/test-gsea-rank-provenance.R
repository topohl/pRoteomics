# The GSEA ranked list, and what can honestly be said about it.
#
# Figure 3 g/h/i shows seven proteins per exemplar term, chosen by ranking that
# term's leading edge on the stored rank statistic. The released inventories
# deliberately do NOT carry that statistic, because the upstream theme table
# supplies leading_edge_genes alphabetically sorted, so a reader of the release
# can confirm eligibility and universe size but not which seven the rule picked.
# That disclosure is pinned by test-display-selection-disclosure.R and stays
# true.
#
# This file pins the other half: internally, the rank IS recoverable, and the
# selection IS reproducible. The ranked vector handed to gseGO() was built by
# R/enrichment/protein_group_enrichment_utils.R as
#
#     setNames(collapsed_statistic, GeneSymbol)[order(collapsed_statistic,
#                                                     decreasing = TRUE)]
#
# from the per-comparison collapsed_gene_input.csv, and that file is retained.
# So the historical ranked list is a deterministic exact reconstruction, not an
# approximation - provided nothing re-sorts the stored rows, and provided the
# file is actually reachable, which for the microglia comparisons it is not
# without the Phase 6H.3 staging contract.

source(testthat::test_path("..", "..", "R", "paths.R"))

CP <- repo_path("data", "processed", "04_differential_expression_enrichment",
                "clusterProfiler")
RULES <- repo_path("docs", "FIGURE_SELECTION_RULES.md")

# The three exemplar programs of Figure 3 d/e/f, and the seven proteins each
# contributes to 3 g/h/i. Both halves are recorded in FIGURE_SELECTION_RULES.md;
# the test below checks the document still says so, so these cannot drift apart.
EXEMPLARS <- list(
  list(key = "synaptic", dataset = "neuron_neuropil", unit = "CA3_sr",
       comparison = "CA3srsus_CA3srres", term = "GO:0099536", universe = 255L,
       shown = c("App", "Cnr1", "Dbi", "Eif4ebp2", "Ly6h", "Plppr4", "Synpo")),
  list(key = "mRNA processing", dataset = "neuron_soma", unit = "CA2_sp",
       comparison = "CA2spsus_CA2spres", term = "GO:0006397", universe = 112L,
       shown = c("Cirbp", "Csdc2", "Dcps", "Ddx23", "Lsm3", "Lsm8", "Rbm8a")),
  list(key = "OXPHOS", dataset = "microglia", unit = "CA1_microglia",
       comparison = "CA1microgliasus_CA1microgliares", term = "GO:0006119", universe = 48L,
       shown = c("Cox5b", "Cox6b1", "Iscu", "Ndufs8", "Ndufv2", "Ndufv3", "Uqcrh"))
)

collapsed_input <- function(e) {
  file.path(CP, e$dataset, "phenotype_within_unit", e$unit, e$comparison,
            "protein_group_audits", "collapsed_gene_input.csv")
}
gsea_result <- function(e) {
  file.path(CP, e$dataset, "phenotype_within_unit", e$unit, e$comparison,
            "GO", "BP", "GSEA_BP_results_full.csv")
}

# Resolution goes through the canonical contract, never around it. A path past
# the wall is staged; it is never re-rooted or read directly.
resolve_for_read <- function(path, stage_root) {
  status <- input_addressability(path)
  if (identical(status, INPUT_STATUS_PRESENT)) return(path)
  if (!identical(status, INPUT_STATUS_OVER_LIMIT)) return(NA_character_)
  dest <- staged_destination(path, stage_root)
  staged <- stage_addressable_copies(path, dest)
  if (!nrow(staged) || !file.exists(staged$staged_path[[1]])) return(NA_character_)
  staged$staged_path[[1]]
}

# The one line of arithmetic the whole figure rests on, quoted from the
# producer rather than paraphrased.
rebuild_ranked <- function(collapsed) {
  stats::setNames(collapsed$collapsed_statistic,
                  collapsed$GeneSymbol)[order(collapsed$collapsed_statistic,
                                              decreasing = TRUE)]
}

testthat::test_that("the producer still builds the ranked vector the way this file assumes", {
  # If the construction changes, every claim below is about a list that is no
  # longer the one gseGO() receives, so read it out of the source.
  src <- paste(readLines(repo_path("R", "enrichment", "protein_group_enrichment_utils.R"),
                         warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("order(collapsed$collapsed_statistic, decreasing = TRUE)",
                              src, fixed = TRUE))
  testthat::expect_true(grepl("gene_collapse_rule = \"median_finite_statistics\"",
                              src, fixed = TRUE))
  # and the collapse is over finite source statistics of eligible rows only
  testthat::expect_true(grepl("eligibility_status == \"eligible\" & is.finite(transform$source_statistic)",
                              src, fixed = TRUE))
})

testthat::test_that("every retained ranked input reproduces the gene count GSEA recorded", {
  # The independent check. n_genes was written into the manifest at GSEA time
  # by a different code path than the one that wrote collapsed_gene_input.csv,
  # so agreement is evidence the retained file is the input that was used, not
  # merely a file of the right shape.
  stage <- withr::local_tempdir()
  checked <- 0L; staged <- 0L
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    mf <- file.path(CP, ds, "clusterProfiler_manifest.csv")
    testthat::skip_if_not(file.exists(mf), paste("manifest absent:", ds))
    m <- utils::read.csv(mf, stringsAsFactors = FALSE)
    m <- m[m$result_type == "GSEA_GO" & m$ontology == "BP", , drop = FALSE]
    for (i in seq_len(nrow(m))) {
      declared <- sub("^P://", "", m$collapsed_gene_input_file[i])
      path <- repo_path(declared)
      readable <- resolve_for_read(path, stage)
      if (is.na(readable)) next
      if (!identical(readable, path)) staged <- staged + 1L
      d <- utils::read.csv(readable, stringsAsFactors = FALSE)
      ranked <- rebuild_ranked(d)
      testthat::expect_identical(length(ranked), as.integer(m$n_genes[i]),
        info = paste(ds, m$comparison[i]))
      checked <- checked + 1L
    }
  }
  testthat::skip_if(checked == 0L, "no clusterProfiler inputs retained in this checkout")
  # 54 GSEA_GO/BP comparisons, of which the 9 microglia ones are over the wall
  testthat::expect_gte(checked, 54L)
  testthat::expect_gte(staged, 9L)
})

testthat::test_that("the reconstruction is exact, not merely well-ordered", {
  # order() is stable, so a tie would be resolved by the stored row order and
  # the reconstruction would depend on that order being the one split() made.
  # It is - the stored order is levels(as.factor(GeneSymbol)) - but the point
  # is moot: these statistics carry no ties at all, so the ranked list is
  # unique and no collation rule can reach it.
  stage <- withr::local_tempdir()
  for (e in EXEMPLARS) {
    f <- resolve_for_read(collapsed_input(e), stage)
    testthat::skip_if(is.na(f), paste("ranked input unreachable:", e$comparison))
    d <- utils::read.csv(f, stringsAsFactors = FALSE)

    testthat::expect_identical(anyDuplicated(d$GeneSymbol), 0L, info = e$key)
    testthat::expect_identical(anyDuplicated(d$collapsed_statistic), 0L,
      info = paste(e$key, "gained a tie; the ranked order is no longer unique"))
    testthat::expect_identical(d$GeneSymbol, levels(as.factor(d$GeneSymbol)),
      info = paste(e$key, "stored rows are no longer in split() order"))

    # a permutation of the rows must not move the reconstruction
    set.seed(1L)
    shuffled <- d[sample(nrow(d)), , drop = FALSE]
    testthat::expect_identical(names(rebuild_ranked(shuffled)),
                               names(rebuild_ranked(d)), info = e$key)
    testthat::expect_true(all(diff(unname(rebuild_ranked(d))) <= 0), info = e$key)
  }
})

testthat::test_that("the committed rule reproduces the committed seven, for all three exemplars", {
  # The falsifiable one. The rule is taken from FIGURE_SELECTION_RULES.md -
  # rank the exemplar's leading edge by |stored rank statistic| and keep the
  # top 7 - and is applied to the recovered historical input. Nothing here is
  # tuned to make the answer come out yes: a change to the selection, to the
  # ranked input, or to the leading edge fails this.
  stage <- withr::local_tempdir()
  for (e in EXEMPLARS) {
    gf <- gsea_result(e)
    testthat::skip_if_not(file.exists(gf), paste("GSEA result absent:", e$comparison))
    g <- utils::read.csv(gf, stringsAsFactors = FALSE)
    row <- g[g$ID == e$term, , drop = FALSE]
    testthat::expect_identical(nrow(row), 1L, info = e$key)

    leading_edge <- unlist(strsplit(row$core_enrichment[1], "/", fixed = TRUE),
                           use.names = FALSE)
    testthat::expect_identical(length(leading_edge), e$universe,
      info = paste(e$key, "leading-edge size no longer matches the documented universe"))

    f <- resolve_for_read(collapsed_input(e), stage)
    testthat::skip_if(is.na(f), paste("ranked input unreachable:", e$comparison))
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    statistic <- stats::setNames(d$collapsed_statistic, d$GeneSymbol)
    testthat::expect_true(all(leading_edge %in% names(statistic)), info = e$key)

    v <- statistic[leading_edge]
    top7 <- sort(leading_edge[order(abs(v), decreasing = TRUE)][seq_len(7L)])
    testthat::expect_identical(top7, sort(e$shown), info = e$key)

    # and the cut is not knife-edge: rank 7 is strictly ahead of rank 8
    o <- order(abs(v), decreasing = TRUE)
    testthat::expect_gt(abs(v[[o[7L]]]), abs(v[[o[8L]]]))
  }
})

testthat::test_that("the reproduction discriminates - other plausible rules do not match", {
  # A reproduction only means something if a WRONG rule gives a WRONG answer.
  # Six candidate rules were applied to the same recovered inputs. Four fail,
  # including the two most tempting: the signed ranked-list order, and the
  # order clusterProfiler itself emits core_enrichment in. The committed rule
  # matches 3/3.
  #
  # One alternative also matches 3/3 - smallest p-value - and that is not a
  # second independent rule. The rank statistic is the moderated t and its
  # p-value is a monotone function of |t| at fixed df, so ordering by |t|
  # descending and by p ascending are the same ordering. They agree because
  # they are one rule written two ways, which is asserted below rather than
  # assumed.
  stage <- withr::local_tempdir()
  survives <- c(signed = 0L, core_order = 0L, abs_logfc = 0L, alphabetical = 0L)
  matched <- 0L
  for (e in EXEMPLARS) {
    gf <- gsea_result(e)
    testthat::skip_if_not(file.exists(gf), paste("GSEA result absent:", e$comparison))
    g <- utils::read.csv(gf, stringsAsFactors = FALSE)
    row <- g[g$ID == e$term, , drop = FALSE]
    le <- unlist(strsplit(row$core_enrichment[1], "/", fixed = TRUE), use.names = FALSE)
    f <- resolve_for_read(collapsed_input(e), stage)
    testthat::skip_if(is.na(f), paste("ranked input unreachable:", e$comparison))
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    v <- stats::setNames(d$collapsed_statistic, d$GeneSymbol)[le]
    fc <- stats::setNames(d$collapsed_logfc, d$GeneSymbol)[le]
    p <- stats::setNames(d$collapsed_p_value, d$GeneSymbol)[le]
    cmt <- sort(e$shown)

    if (identical(sort(le[order(abs(v), decreasing = TRUE)][1:7]), cmt)) matched <- matched + 1L
    if (identical(sort(le[order(v, decreasing = TRUE)][1:7]), cmt))
      survives[["signed"]] <- survives[["signed"]] + 1L
    if (identical(sort(le[1:7]), cmt))
      survives[["core_order"]] <- survives[["core_order"]] + 1L
    if (identical(sort(le[order(abs(fc), decreasing = TRUE)][1:7]), cmt))
      survives[["abs_logfc"]] <- survives[["abs_logfc"]] + 1L
    if (identical(sort(sort(le)[1:7]), cmt))
      survives[["alphabetical"]] <- survives[["alphabetical"]] + 1L

    # |statistic| descending and p ascending are the same order, so their
    # agreement is an identity rather than corroboration
    testthat::expect_identical(order(abs(v), decreasing = TRUE), order(p), info = e$key)
  }
  testthat::expect_identical(matched, 3L)
  testthat::expect_lt(survives[["signed"]], 3L)
  testthat::expect_lt(survives[["core_order"]], 3L)
  testthat::expect_identical(unname(survives[["abs_logfc"]]), 0L)
  testthat::expect_identical(unname(survives[["alphabetical"]]), 0L)
})

testthat::test_that("the rules document still states both the rule and its answer", {
  # The seven names above are not a magic constant: they are what the document
  # publishes. If either side is edited alone this fails, which is the point.
  testthat::skip_if_not(file.exists(RULES), "rules document absent")
  txt <- paste(readLines(RULES, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("rank by\\s*\\*\\*\\|stored rank statistic\\|\\*\\*", txt))
  testthat::expect_true(grepl("keep the \\*\\*top 7\\*\\*", txt))
  for (e in EXEMPLARS) {
    testthat::expect_true(grepl(paste(e$shown, collapse = " "), txt, fixed = TRUE),
      info = paste(e$key, "displayed proteins no longer listed in the rules document"))
    testthat::expect_true(grepl(e$term, txt, fixed = TRUE), info = e$key)
  }
})

testthat::test_that("the executed ranked order is STORED, not merely reconstructible", {
  # Stronger than everything above. rank_statistic_sensitivity_audit.csv is
  # written from names(gene_inputs$sensitivity$median), and `median` IS the
  # ranked vector handed to gseGO() - so this file records the executed rank
  # ORDER, row by row, and not just the values it was derived from.
  #
  # That matters for how the provenance may be described. The ranked list is
  # BYTE_EXACT_STORED_ORDER, and the rebuild from collapsed_gene_input.csv is
  # independent corroboration of it rather than the primary evidence. Two files
  # written by different code paths, plus the manifest's n_genes, all agree.
  #
  # It also matters for retention. This file currently reads as a disposable
  # diagnostic - audits/phase6h_over_maxpath_reader_check.csv records one
  # active code reference for it, and that reference is the WRITE site, so it
  # has no readers at all. A file with no readers looks deletable. This one is
  # the primary record of an executed ranking that feeds a published figure.
  stage <- withr::local_tempdir()
  checked <- 0L
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    mf <- file.path(CP, ds, "clusterProfiler_manifest.csv")
    testthat::skip_if_not(file.exists(mf), paste("manifest absent:", ds))
    m <- utils::read.csv(mf, stringsAsFactors = FALSE)
    m <- m[m$result_type == "GSEA_GO" & m$ontology == "BP", , drop = FALSE]
    for (i in seq_len(nrow(m))) {
      dir <- dirname(repo_path(sub("^P://", "", m$collapsed_gene_input_file[i])))
      sens <- resolve_for_read(file.path(dir, "rank_statistic_sensitivity_audit.csv"), stage)
      coll <- resolve_for_read(file.path(dir, "collapsed_gene_input.csv"), stage)
      if (is.na(sens) || is.na(coll)) next
      s <- utils::read.csv(sens, stringsAsFactors = FALSE)
      ranked <- rebuild_ranked(utils::read.csv(coll, stringsAsFactors = FALSE))

      testthat::expect_identical(s$GeneSymbol, names(ranked),
        info = paste(ds, m$comparison[i], "stored order is not the executed order"))
      testthat::expect_true(all(diff(s$median_statistic) <= 0),
        info = paste(ds, m$comparison[i], "stored statistic is not monotone"))
      testthat::expect_equal(unname(s$median_statistic), unname(ranked),
        info = paste(ds, m$comparison[i]))
      testthat::expect_identical(nrow(s), as.integer(m$n_genes[i]),
        info = paste(ds, m$comparison[i]))
      checked <- checked + 1L
    }
  }
  testthat::skip_if(checked == 0L, "no sensitivity audits retained in this checkout")
  testthat::expect_gte(checked, 54L)
})

testthat::test_that("the stored ranked order is reachable only through the staging contract", {
  # 18 of the 54 sensitivity audits are past the 260-character wall, against 9
  # of the collapsed inputs, because the basename is longer. Re-rooting them by
  # hand is what the Phase 6H.3 contract exists to prevent, so the count is
  # pinned: if it reaches zero, either the repository moved or someone started
  # reading these a different way, and both deserve a deliberate look.
  stage <- withr::local_tempdir()
  over <- 0L; total <- 0L
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    mf <- file.path(CP, ds, "clusterProfiler_manifest.csv")
    testthat::skip_if_not(file.exists(mf), paste("manifest absent:", ds))
    m <- utils::read.csv(mf, stringsAsFactors = FALSE)
    m <- m[m$result_type == "GSEA_GO" & m$ontology == "BP", , drop = FALSE]
    for (i in seq_len(nrow(m))) {
      p <- file.path(dirname(repo_path(sub("^P://", "", m$collapsed_gene_input_file[i]))),
                     "rank_statistic_sensitivity_audit.csv")
      total <- total + 1L
      if (identical(input_addressability(p), INPUT_STATUS_OVER_LIMIT)) over <- over + 1L
    }
  }
  testthat::skip_if(total == 0L, "no clusterProfiler outputs in this checkout")
  testthat::expect_identical(total, 54L)
  testthat::expect_identical(over, 18L)
  # and every one of them is recoverable, so "over the wall" never means lost
  recovered <- 0L
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    m <- utils::read.csv(file.path(CP, ds, "clusterProfiler_manifest.csv"),
                         stringsAsFactors = FALSE)
    m <- m[m$result_type == "GSEA_GO" & m$ontology == "BP", , drop = FALSE]
    for (i in seq_len(nrow(m))) {
      p <- file.path(dirname(repo_path(sub("^P://", "", m$collapsed_gene_input_file[i]))),
                     "rank_statistic_sensitivity_audit.csv")
      if (!is.na(resolve_for_read(p, stage))) recovered <- recovered + 1L
    }
  }
  testthat::expect_identical(recovered, 54L)
})

testthat::test_that("recovering the rank does not make it releasable", {
  # Recovering provenance internally is not the same as publishing it, and the
  # release must not quietly acquire the statistic as a side effect of this
  # work. The disclosure contract stays exactly as Phase 6H.11A left it.
  bundle <- repo_path("exports", "publication_source_data",
                      "supplementary_selection_inventories")
  testthat::skip_if_not(dir.exists(bundle), "release bundle absent")
  le <- file.path(bundle, "leading_edge_protein_inventory.csv")
  testthat::skip_if_not(file.exists(le), "leading-edge inventory absent")
  header <- names(utils::read.csv(le, nrows = 1L, stringsAsFactors = FALSE))
  for (forbidden in c("collapsed_statistic", "rank_statistic", "ranked_list_position",
                      "source_statistic", "gsea_rank"))
    testthat::expect_false(forbidden %in% header, info = forbidden)
})
