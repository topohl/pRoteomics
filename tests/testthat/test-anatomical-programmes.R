source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "control_spatial_identity_utils.R"))
source(repo_path("R", "anatomical_programme_utils.R"))

NEUROPIL_UNITS <- c("CA1_SLM", "CA1_SO", "CA1_SR", "CA2_SLM", "CA2_SO", "CA2_SR",
                    "CA3_SO", "CA3_SR", "DG_MO", "DG_PO")

testthat::test_that("locked families name exactly the registry's locked contrasts", {
  fam <- anatomical_programme_families()
  locked <- fam$contrast[fam$source == "locked"]
  reg <- c(names(control_spatial_contrast_registry("neuron_soma", c("CA1", "CA2", "CA3", "DG"))),
           names(control_spatial_contrast_registry("neuron_neuropil", NEUROPIL_UNITS)))
  testthat::expect_true(all(locked %in% reg))
  testthat::expect_false(any(duplicated(fam$contrast)))
})

testthat::test_that("added neuropil contrasts are target versus the mean of the other strata", {
  specs <- anatomical_programme_added_specs("neuron_neuropil", NEUROPIL_UNITS)
  testthat::expect_setequal(names(specs), c(
    "CA2_SO_vs_mean_other_CA2_strata", "CA2_SR_vs_mean_other_CA2_strata",
    "CA2_SLM_vs_mean_other_CA2_strata", "CA3_SO_vs_mean_other_CA3_strata",
    "CA3_SR_vs_mean_other_CA3_strata"))
  testthat::expect_false(any(vapply(specs, function(s) s$manuscript_locked, logical(1))))
  w <- specs[["CA2_SR_vs_mean_other_CA2_strata"]]$weights
  testthat::expect_equal(unname(w[c("CA2_SR", "CA2_SO", "CA2_SLM")]), c(1, -0.5, -0.5))
  w3 <- specs[["CA3_SO_vs_mean_other_CA3_strata"]]$weights
  testthat::expect_equal(unname(w3[c("CA3_SO", "CA3_SR")]), c(1, -1))
  # nothing added collides with, or redefines, a locked contrast
  locked <- names(control_spatial_contrast_registry("neuron_neuropil", NEUROPIL_UNITS))
  testthat::expect_length(intersect(names(specs), locked), 0L)
  # a missing stratum is an error, not a silently different contrast
  testthat::expect_error(anatomical_programme_added_specs("neuron_neuropil", setdiff(NEUROPIL_UNITS, "CA2_SLM")))
})

testthat::test_that("microglia contrasts are the registry's own, unlocked", {
  specs <- anatomical_programme_added_specs("microglia", c("CA1", "CA2", "CA3", "DG"))
  reg <- control_spatial_contrast_registry("microglia", c("CA1", "CA2", "CA3", "DG"))
  testthat::expect_identical(specs, reg[names(specs)])
  testthat::expect_length(specs, 4L)
})

fake_go <- function() {
  fam <- data.frame(family = "f", dataset = "d", source = "added", unit = c("A", "B", "C"),
                    contrast = c("cA", "cB", "cC"), family_label = "F", unit_order = 1:3,
                    stringsAsFactors = FALSE)
  row <- function(ct, id, nes, q, size = 50, core = "g1/g2/g3/g4") {
    data.frame(ID = id, Description = paste("term", id), setSize = size, NES = nes,
               p.adjust = q, core_enrichment = core, contrast = ct, stringsAsFactors = FALSE)
  }
  go <- rbind(
    row("cA", "T1", 2.5, 0.001, core = "a1/a2/a3/a4"),     # specific to A, top
    row("cB", "T1", -1, 0.2), row("cC", "T1", -1, 0.3),
    row("cA", "T2", 2.4, 0.001, core = "a1/a2/a3/x9"),     # nested in T1 -> redundant
    row("cB", "T2", -1, 0.2), row("cC", "T2", 0.5, 0.9),
    row("cA", "T3", 2.3, 0.01, core = "b1/b2"),            # supported in A and B -> not specific
    row("cB", "T3", 2.0, 0.01, core = "b1/b2"), row("cC", "T3", -2, 0.01),
    row("cA", "T4", 2.2, 0.01, size = 500, core = "c1"),   # too large
    row("cB", "T4", 0.1, 0.9), row("cC", "T4", 0.1, 0.9),
    row("cA", "T5", 2.1, 0.02, core = "d1/d2"),            # artefact branch
    row("cB", "T5", 0.1, 0.9), row("cC", "T5", 0.1, 0.9),
    row("cA", "T6", 2.0, 0.03, core = "e1/e2"),            # kept
    row("cB", "T6", 0.1, 0.9), row("cC", "T6", 0.1, 0.9),
    row("cA", "T7", 1.9, 0.03, core = "f1/f2"),            # kept (third)
    row("cB", "T7", 0.1, 0.9), row("cC", "T7", 0.1, 0.9),
    row("cA", "T8", 1.8, 0.04, core = "h1/h2"),            # beyond k
    row("cB", "T8", 0.1, 0.9), row("cC", "T8", 0.1, 0.9),
    row("cA", "T9", 2.6, 0.06, core = "i1"),               # q too high
    row("cB", "T9", 0.1, 0.9), row("cC", "T9", 0.1, 0.9))
  list(go = go, fam = fam, artefact_ids = "T5")
}

testthat::test_that("the unit-specific rule selects, skips and records as stated", {
  x <- fake_go()
  cand <- anatomical_programme_candidates(x$go, x$fam, anatomical_programme_rule(k = 3L), artefact_ids = x$artefact_ids)
  a <- cand[cand$unit == "A", , drop = FALSE]
  testthat::expect_identical(a$ID[a$selected][order(a$display_rank[a$selected])], c("T1", "T6", "T7"))
  testthat::expect_identical(a$redundant_with[a$ID == "T2"], "T1")
  testthat::expect_false(a$unit_specific[a$ID == "T3"])
  testthat::expect_false(a$candidate[a$ID == "T4"])
  testthat::expect_false(a$candidate[a$ID == "T5"])
  testthat::expect_false(a$selected[a$ID == "T8"])
  testthat::expect_false("T9" %in% a$ID)                  # not supported, not listed
  testthat::expect_true(cand$selected[cand$unit == "B" & cand$ID == "T3"] == FALSE)
  d <- anatomical_programme_display(x$go, cand, x$fam)
  testthat::expect_setequal(unique(d$ID), c("T1", "T6", "T7"))
  testthat::expect_equal(nrow(d), 9L)                       # three terms x three units
  testthat::expect_true(all(d$selected_for_unit == "A"))
})

testthat::test_that("GO artefact branches are recognised", {
  testthat::skip_if_not_installed("GO.db")
  testthat::expect_identical(
    anatomical_programme_artefact_flag(c("GO:0035036", "GO:0002027", "GO:0006119")),
    c(TRUE, TRUE, FALSE))  # sperm-egg recognition, regulation of heart rate, oxidative phosphorylation
})

testthat::test_that("sibling pairs cover every pair once; two-unit families reuse their contrast", {
  p <- anatomical_programme_pairs()
  n <- table(p$family)
  testthat::expect_equal(as.integer(n[c("soma_regions", "CA1_strata", "CA2_strata", "CA3_strata", "DG_layers", "microglia_regions")]),
                         c(6L, 3L, 3L, 1L, 1L, 6L))
  testthat::expect_true(all(p$reuses_one_vs_rest[p$family %in% c("CA3_strata", "DG_layers")]))
  testthat::expect_identical(p$contrast[p$family == "DG_layers"], "DG_MO_vs_mean_other_DG_layers")
  specs <- anatomical_programme_pairwise_specs("neuron_neuropil", NEUROPIL_UNITS,
                                               which_families = c("CA1_strata", "CA2_strata", "DG_layers"))
  testthat::expect_length(specs, 6L)
  w <- specs[["pairwise_CA2_strata_SO_vs_SLM"]]$weights
  testthat::expect_equal(unname(w[c("CA2_SO", "CA2_SLM")]), c(1, -1))
  testthat::expect_equal(sum(abs(w)), 2)
})

testthat::test_that("pairwise support is oriented towards the unit and needs every sibling", {
  fam <- data.frame(family = "f", dataset = "d", source = "added", unit = c("A", "B", "C"),
                    anatomical_unit = c("A", "B", "C"), contrast = c("cA", "cB", "cC"), family_label = "F",
                    unit_order = 1:3, shown_in_figure = TRUE, stringsAsFactors = FALSE)
  cand <- data.frame(family = "f", unit = c("B", "C"), contrast = c("cB", "cC"), ID = c("T1", "T2"),
                     core_enrichment = c("g1/g2/g3/g4/g5/g6", "h1/h2"), stringsAsFactors = FALSE)
  pw <- data.frame(contrast = c("pairwise_f_A_vs_B", "pairwise_f_B_vs_C", "pairwise_f_A_vs_B", "pairwise_f_B_vs_C",
                                "pairwise_f_A_vs_C", "pairwise_f_A_vs_C"),
                   ID = c("T1", "T1", "T2", "T2", "T1", "T2"),
                   NES = c(-2.0, 1.5, 0.1, -0.2, 0.3, -1.8), p.adjust = c(0.001, 0.01, 0.9, 0.9, 0.5, 0.002),
                   stringsAsFactors = FALSE)
  side <- data.frame(contrast = c("cB", "cB", "cC", "cC"), side = c("L", "R", "L", "R"), ID = c("T1", "T1", "T2", "T2"),
                     NES = c(1.2, 0.8, 1.0, -0.3), p.adjust = c(0.1, 0.3, 0.2, 0.9), stringsAsFactors = FALSE)
  gs <- data.frame(contrast = "cB", GeneSymbol = paste0("g", 1:6), collapsed_statistic = c(10, 1, 1, 1, 1, 1),
                   collapsed_fdr = c(0.001, 0.01, 0.2, 0.2, 0.2, 0.2), collapsed_logfc = 1, stringsAsFactors = FALSE)
  m <- anatomical_programme_term_metrics(cand, gs, pw, side, fam)
  # T1 in B: B vs A = -(-2.0) = 2.0 (q .001); B vs C = 1.5 (q .01) -> supported against both
  testthat::expect_true(m$pairwise_supported_all[1])
  testthat::expect_equal(m$pairwise_min_NES[1], 1.5)
  testthat::expect_true(m$both_hemispheres_positive[1])
  testthat::expect_equal(m$top_protein_share[1], 10 / 15)
  testthat::expect_identical(m$top_leading_protein[1], "g1")
  testthat::expect_equal(m$n_leading_edge_supported[1], 2L)
  # T2 in C: C vs A = 1.8 (q .002); C vs B = -(-0.2) = 0.2 (q .9) -> not supported against B
  testthat::expect_false(m$pairwise_supported_all[2])
  testthat::expect_false(m$both_hemispheres_positive[2])
})
