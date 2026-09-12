# =====================================================================
# Story-v5 panels: Figure-2 spatial-validation redesign + Extended Data.
#
# WHY 2f AND 2g ARE BEING RE-ENCODED
#
# The audit (figure2_spatial_validation_hierarchy_audit.csv) shows the sparsity
# in both panels is anatomical hierarchy, not missing data:
#
#   2f  occupancy is BLOCK-DIAGONAL. Regional contrasts pair only with regional
#       reference signatures (18 cells) and CA1-strata contrasts only with
#       laminar ones (12 cells); the off-diagonal is empty by construction.
#       52% of a full contrast x signature matrix is structurally meaningless.
#   2g  is a PAIRED LIST, not a matrix: exactly two selected GO terms per
#       contrast, 10 of 12 terms unique to one contrast. A term x contrast
#       grid would be 83% blank.
#
# So neither is drawn as a grid. 2f becomes a hierarchical forest whose blocks
# ARE the anatomical levels; 2g becomes contrast-centric small multiples.
#
# Downstream only: no model fit, no enrichment, no FDR. Values are read from
# the canonical validation source data.
# =====================================================================

# Canonical anatomical ordering, shared by every panel in this layer so the
# external and internal validation panels can be read as one ladder.
s5_domain_levels <- function() {
  c("Soma region — tissue reference",
    "Neuropil region — synaptosome reference",
    "CA1 strata — synaptosome reference")
}
s5_contrast_order <- function() {
  c("CA1_vs_mean_other_soma_regions", "CA2_vs_mean_other_soma_regions",
    "CA3_vs_mean_other_soma_regions", "DG_vs_mean_other_soma_regions",
    "DG_neuropil_vs_mean_non_DG_regions", "CA1_SO_vs_CA3_SO",
    "CA1_SLM_vs_mean_other_CA1_strata", "CA1_SO_vs_mean_other_CA1_strata",
    "CA1_SR_vs_mean_other_CA1_strata")
}
s5_pretty_contrast <- function(x) {
  y <- gsub("_", " ", x)
  y <- sub("vs mean other soma regions", "vs other soma regions", y)
  y <- sub("vs mean other CA1 strata", "vs other CA1 strata", y)
  y <- sub("vs mean non DG regions", "vs non-DG regions", y)
  y
}
s5_block <- function(domain) {
  ifelse(grepl("strata", domain), "CA1 laminar identity", "Regional identity")
}

# ---------------------------------------------------------------- 2f: F2F_A
#
# Hierarchical forest. Rows are the real prespecified contrasts; no cell is
# created for a pairing that cannot exist. Expected pairings are solid, the
# off-target specificity comparisons that make the test meaningful are open.
s5_kaulich_forest <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  k <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family

  z <- data.frame(
    dataset = k$dataset, contrast = k$internal_contrast,
    domain = k$validation_domain, signature = k$external_signature,
    NES = k$NES, FDR = k$p_adjust, expected = k$expected_match %in% TRUE,
    stringsAsFactors = FALSE)
  z <- z[is.finite(z$NES), , drop = FALSE]
  z$supported <- is.finite(z$FDR) & z$FDR < 0.05
  z$block <- factor(s5_block(z$domain),
                    levels = c("Regional identity", "CA1 laminar identity"))
  ord <- intersect(s5_contrast_order(), unique(z$contrast))
  z$row <- factor(s5_pretty_contrast(z$contrast),
                  levels = rev(s5_pretty_contrast(ord)))

  p <- ggplot2::ggplot(z, ggplot2::aes(NES, row)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey65") +
    ggplot2::geom_point(ggplot2::aes(colour = expected, shape = expected,
                                     alpha = supported), size = 1.35,
                        stroke = 0.35) +
    ggplot2::geom_text(data = z[z$expected, , drop = FALSE],
                       ggplot2::aes(label = signature), family = fam,
                       size = nv_size(5), nudge_y = 0.34, colour = "#1F3D52") +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "#1F3D52", "FALSE" = "#9E9A92"),
      labels = c("TRUE" = "expected match", "FALSE" = "specificity comparison"),
      name = NULL) +
    ggplot2::scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = "expected match", "FALSE" = "specificity comparison"),
      name = NULL) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35),
                                guide = "none") +
    ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y",
                        switch = "y") +
    ggplot2::labs(x = "external signature NES", y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5),
                   strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 90, size = 5,
                                                             face = "bold"),
                   legend.position = "bottom",
                   legend.key.height = ggplot2::unit(1.8, "mm"))

  z$structural_note <- paste0(
    "only pairings that exist are plotted. Regional contrasts pair with ",
    "regional reference signatures and CA1-strata contrasts with laminar ones; ",
    "the off-diagonal of a contrast x signature grid is structurally empty, ",
    "not missing data")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---------------------------------------------------------------- 2g: F2G_B
#
# Contrast-centric small multiples. Each contrast carries its own two canonical
# terms, so a shared term axis would be 83% empty. Point size is the stored
# setSize; nothing is invented.
s5_internal_multiples <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  g <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  z <- data.frame(
    contrast = g$contrast, term = g$Description, NES = g$NES,
    FDR = g$p_adjust, setSize = g$setSize, stringsAsFactors = FALSE)
  z <- z[is.finite(z$NES), , drop = FALSE]
  z$supported <- is.finite(z$FDR) & z$FDR < 0.05
  z$block <- factor(ifelse(grepl("CA1_strata", z$contrast),
                           "CA1 laminar identity", "Regional identity"),
                    levels = c("Regional identity", "CA1 laminar identity"))
  ord <- intersect(s5_contrast_order(), unique(z$contrast))
  z$facet <- factor(s5_pretty_contrast(z$contrast),
                    levels = s5_pretty_contrast(ord))
  # term order within its own facet, strongest first
  z <- z[order(z$facet, -z$NES), ]
  z$term_id <- paste(z$facet, z$term)
  z$term_f <- factor(z$term_id, levels = rev(unique(z$term_id)))

  p <- ggplot2::ggplot(z, ggplot2::aes(NES, term_f)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = term_f),
                          colour = "grey85", linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(size = setSize, colour = supported)) +
    ggplot2::geom_text(ggplot2::aes(x = 0.08, label = term), hjust = 0,
                       family = fam, size = nv_size(5), colour = "black") +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#1F3D52", "FALSE" = "#B9B9B4"),
                                 guide = "none") +
    ggplot2::scale_size_area(max_size = 1.8, name = "genes",
                             breaks = c(50, 200, 450)) +
    ggplot2::scale_x_continuous(limits = c(0, max(z$NES) * 1.06),
                                expand = c(0, 0)) +
    ggplot2::facet_grid(facet ~ ., scales = "free_y", space = "free_y",
                        switch = "y") +
    ggplot2::labs(x = "internal anatomical GSEA NES", y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0, size = 4.9,
                                                             hjust = 1),
                   panel.spacing.y = ggplot2::unit(0.4, "mm"),
                   legend.position = "bottom",
                   legend.key.height = ggplot2::unit(1.8, "mm"))

  z$structural_note <- paste0(
    "each contrast carries its own two canonical terms; 10 of 12 terms occur ",
    "in a single contrast, so a term x contrast grid would be 83% blank. ",
    "Point size is the stored setSize. All 14 are FDR-supported")
  write_csv_safe(z[, setdiff(names(z), c("term_id", "term_f"))], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ============================================== EXTENDED DATA helpers

# ED3: CA2-SLM QC. Missingness by sample with QC flags, and the normalisation
# displacement that explains the burden.
s5_ed_ca2_missingness <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  s <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  s$lab <- paste0(s$AnimalID, "_", s$hemisphere)
  s$flag <- ifelse(s$qc_flag == "PASS", "", s$qc_flag)
  p <- ggplot2::ggplot(s, ggplot2::aes(stats::reorder(lab, -fraction_missing_preimputation),
                                       fraction_missing_preimputation,
                                       fill = StressGroup)) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = flag), vjust = -0.3,
                       size = nv_size(5), family = fam) +
    ggplot2::scale_fill_manual(values = nv_group_colours(), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = NULL, y = "fraction missing\nbefore imputation") +
    nv_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1,
                                                       size = 5),
                   legend.position = "right")
  write_csv_safe(s, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

s5_ed_ca2_displacement <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  b <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  b$flag <- ifelse(b$qc_flag == "PASS", "", b$qc_flag)
  p <- ggplot2::ggplot(b, ggplot2::aes(fraction_missing_preimputation,
                                       mean_centred_value_of_always_observed_proteins)) +
    ggplot2::geom_point(ggplot2::aes(colour = StressGroup), size = 1.1) +
    ggplot2::geom_text(ggplot2::aes(label = flag), family = fam, size = nv_size(5),
                       vjust = -0.9, colour = "grey25") +
    ggplot2::scale_colour_manual(values = nv_group_colours(), name = NULL) +
    ggplot2::labs(x = "fraction missing before imputation",
                  y = "centred value of proteins\nobserved in every sample") +
    nv_theme(grid = "y") +
    ggplot2::theme(legend.position = "right")
  b$mechanism <- paste0(
    "per-sample median centring uses observed values only, so a heavily ",
    "missing sample gets an inflated median and every one of its proteins is ",
    "displaced downwards - including proteins with no imputed value anywhere")
  write_csv_safe(b, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

s5_ed_ca2_classes <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  r <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  z <- as.data.frame(table(cls = r$CA2_SLM_robustness_class),
                     stringsAsFactors = FALSE)
  z$lab <- gsub("_", " ", z$cls)
  z$lab <- factor(z$lab, levels = z$lab[order(z$Freq)])
  z$col <- c(robust_to_missingness_and_QC = "#1F3D52",
             not_claimable_due_to_QC = "#D1543A",
             insufficient_observed_data = "#C9C6BF")[z$cls]
  p <- ggplot2::ggplot(z, ggplot2::aes(Freq, lab)) +
    ggplot2::geom_col(fill = z$col[order(z$Freq)], width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = Freq), hjust = -0.35,
                       family = fam, size = nv_size(5.4)) +
    ggplot2::scale_x_continuous(limits = c(0, max(z$Freq) * 1.22),
                                expand = c(0, 0)) +
    ggplot2::labs(x = "proteins of the canonical 28", y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
  write_csv_safe(z[, c("cls", "Freq")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

s5_ed_ca2_sensitivity <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  r <- nv_read_csv(repo_path(panel$primary_source))
  z <- data.frame(
    gene = r$gene_symbol, canonical = r$canonical_log2FC_SUS_minus_RES,
    hemi = r$effect_dropping_QC_failed_hemispheres,
    cls = r$CA2_SLM_robustness_class, stringsAsFactors = FALSE)
  z <- z[is.finite(z$canonical) & is.finite(z$hemi), , drop = FALSE]
  p <- ggplot2::ggplot(z, ggplot2::aes(canonical, hemi)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22",
                         linewidth = nv_lw("reference_pt"), colour = "grey65") +
    ggplot2::geom_hline(yintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey80") +
    ggplot2::geom_point(ggplot2::aes(colour = cls), size = 1.1) +
    ggplot2::scale_colour_manual(values = c(
      robust_to_missingness_and_QC = "#1F3D52",
      not_claimable_due_to_QC = "#D1543A",
      insufficient_observed_data = "#C9C6BF"),
      labels = function(x) gsub("_", " ", x), name = NULL) +
    ggplot2::labs(x = "canonical log2FC (SUS − RES)",
                  y = "with QC-failed\nhemispheres dropped") +
    nv_theme(grid = "y") +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.height = ggplot2::unit(1.8, "mm"))
  z$reading <- paste0("points below the identity line lose magnitude when the ",
                      "two QC-failed SUS acquisitions are removed")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ED5: WGCNA phenotype context, with the null stated on the panel.
s5_ed_wgcna_phenotype <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  e <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  lab <- intersect(c("module_display_label", "ModuleLabel_Final"), names(e))[1]
  z <- data.frame(module = if (is.na(lab)) e$module_id else e[[lab]],
                  module_id = e$module_id, contrast = e$contrast,
                  estimate = e$estimate, fdr = e$tier_specific_fdr,
                  stringsAsFactors = FALSE)
  z$module <- factor(z$module, levels = rev(unique(z$module[order(z$module_id)])))
  lim <- max(abs(z$estimate), na.rm = TRUE)
  n_sig <- sum(z$fdr < 0.05, na.rm = TRUE)
  p <- ggplot2::ggplot(z, ggplot2::aes(contrast, module, fill = estimate)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    nv_diverging(limits = c(-lim, lim), name = "model effect") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = sprintf(paste0("%d of %d module x contrast cells reach ",
                    "tier-specific FDR < 0.05 (minimum %.3f): the colour field is ",
                    "DESCRIPTIVE and no cell is marked as significant"),
                    n_sig, nrow(z), min(z$fdr, na.rm = TRUE))) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 5),
                   legend.position = "right",
                   plot.caption = ggplot2::element_text(size = 5, hjust = 0,
                                                        colour = "grey20",
                                                        lineheight = 1.2))
  z$fdr_supported_cells <- n_sig
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ED7: stress vs baseline spatial identity across every robustness subset.
s5_ed_identity_subsets <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  c1 <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  c1$lab <- gsub("_", " ", c1$subset)
  c1$lab <- factor(c1$lab, levels = rev(c1$lab))
  c1$pct <- 100 * c1$fraction_outside_baseline_affinity
  c1$primary <- c1$subset == "CA2_SLM_robustness_qualified"
  p <- ggplot2::ggplot(c1, ggplot2::aes(pct, lab)) +
    ggplot2::geom_col(ggplot2::aes(fill = primary), width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d/%d",
                                                    effect_outside_baseline_affinity,
                                                    n_hits)),
                       hjust = -0.15, family = fam, size = nv_size(5)) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#1F3D52", "FALSE" = "#9E9A92"),
                               guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 122), expand = c(0, 0)) +
    ggplot2::labs(x = "% of hits outside dominant baseline spatial affinity",
                  y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
  c1$note <- paste0("the result survives every robustness restriction; the old ",
                    "rank-10 tail did not and is deliberately not shown")
  write_csv_safe(c1, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ED8: network and coupling nulls, stated explicitly.
s5_ed_nulls <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  n <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  z <- data.frame(dataset = n$dataset, p = n$exact_p,
                  floor = n$min_attainable_p, n_edges = n$n_edges,
                  stringsAsFactors = FALSE)
  z$lab <- nv_dataset_label(z$dataset)
  z$lab <- factor(z$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))
  p <- ggplot2::ggplot(z, ggplot2::aes(p, lab)) +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "22",
                        linewidth = nv_lw("reference_pt"), colour = "#D1543A") +
    ggplot2::geom_segment(ggplot2::aes(x = floor, xend = p, yend = lab),
                          colour = "grey85", linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(x = floor), colour = "grey60", size = 1,
                        shape = 1) +
    ggplot2::geom_point(size = 1.4, colour = "#1F3D52") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("p = %.2f", p)), hjust = -0.3,
                       family = fam, size = nv_size(5)) +
    ggplot2::scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0),
                                breaks = c(0, 0.05, 0.5, 1)) +
    ggplot2::labs(x = "exact whole-network p (open circle = attainable floor)",
                  y = NULL) +
    nv_theme(grid = "x")
  z$reading <- paste0("an informative null: the enumeration had resolution to ",
                      "0.0036 and found no whole-network group difference")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ED8b: the edge-behaviour coupling null.
#
# The stored BH family has 96 rows but only 48 DISTINCT tests: each estimate is
# duplicated across the Phase/Analysis bookkeeping columns. Plotting 96 points
# would double-count, so the panel deduplicates and states both counts. With
# n = 9 animals a single correlation has almost no resolution, which is the
# point: exactly one of 48 intervals excludes zero, and nothing approaches FDR
# support. This is the panel that prevents the reader from mining the edge list.
s5_ed_coupling <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  x <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  keep <- c("Edge", "Outcome", "Change", "window", "estimate",
            "conf.low", "conf.high", "p.value", "n",
            "p.adj_BH_within_outcome", "p.adj_BH_all_edge_phenotype_tests")
  # The stored BH family has 96 rows keyed on (Analysis, Phase, Change, window,
  # Edge, Outcome), but Phase is structurally inapplicable to the three
  # per-animal outcomes (composite z, corticosterone change, sucrose
  # preference): those are not phase-resolved, so Active and Inactive carry the
  # identical test. Only Movement AUC is genuinely phase-resolved. Collapsing on
  # the statistical content keeps every distinct test and drops only the
  # duplicated ones. Note this means BH was applied over 96 slots for 48 tests,
  # which is conservative and so cannot manufacture the null reported here.
  u <- unique(x[, keep])
  dupkey <- paste(x$Edge, x$Outcome, x$estimate, x$p.value, sep = "\r")
  u$n_stored_rows <- as.integer(table(dupkey)[
    paste(u$Edge, u$Outcome, u$estimate, u$p.value, sep = "\r")])
  u$excludes_zero <- u$conf.low > 0 | u$conf.high < 0
  u$edge_lab <- gsub("_", " ", u$Edge)
  u$out_lab <- c(CombZ = "behavioural composite z",
                 delta_cort = "corticosterone change",
                 Movement_AUC_z_vs_CON = "movement AUC vs CON",
                 sucrose_pref = "sucrose preference")[u$Outcome]
  # deterministic within-cell ordering so repeated builds are byte-identical
  u <- u[order(u$out_lab, u$edge_lab, u$window, u$Change, u$estimate), ]
  u$row <- stats::ave(seq_len(nrow(u)),
                      paste(u$out_lab, u$edge_lab), FUN = seq_along)
  u$ypos <- as.integer(factor(u$edge_lab, levels = rev(sort(unique(u$edge_lab))))) +
    (u$row - mean(unique(u$row))) * 0.22
  p <- ggplot2::ggplot(u, ggplot2::aes(estimate, ypos)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey55") +
    ggplot2::geom_segment(ggplot2::aes(x = conf.low, xend = conf.high,
                                       yend = ypos, colour = excludes_zero),
                          linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(colour = excludes_zero), size = 0.7) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#C0442C", "FALSE" = "#9E9A92"),
                                 guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = seq_along(unique(u$edge_lab)),
      labels = rev(sort(unique(u$edge_lab))), expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
    ggplot2::facet_wrap(~out_lab, nrow = 1) +
    ggplot2::labs(x = "Pearson r (95% CI), n = 9 animals", y = NULL,
                  subtitle = sprintf(paste0(
                    "%d distinct tests; %d interval excludes zero; smallest BH ",
                    "= %.2f within outcome and %.2f across the family \u2014 no ",
                    "edge\u2013behaviour association is claimable"),
                    nrow(u), sum(u$excludes_zero),
                    min(u$p.adj_BH_within_outcome),
                    min(u$p.adj_BH_all_edge_phenotype_tests))) +
    nv_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5),
                   panel.spacing.x = ggplot2::unit(1.6, "mm"),
                   plot.subtitle = ggplot2::element_text(size = 5, colour = "grey20"))
  u$reading <- sprintf(
    paste0("%d stored rows collapse to %d distinct tests; %d of %d intervals ",
           "exclude zero and the smallest BH value is %.2f within-outcome ",
           "(%.2f across the whole family), so no edge-behaviour association ",
           "is claimable"),
    nrow(x), nrow(u), sum(u$excludes_zero), nrow(u),
    min(u$p.adj_BH_within_outcome), min(u$p.adj_BH_all_edge_phenotype_tests))
  write_csv_safe(u, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ED8a: animal-level network distance from the CON consensus.
#
# One AnimalID is one independent network replicate, so this shows the actual
# replicate-level spread behind the whole-network permutation in ED8c.
s5_ed_network_distance <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  fam <- nv_palette()$typography$family
  grp <- nv_group_colours()
  dcol <- names(d)[grepl("^distance", names(d))][1]
  if (is.na(dcol)) stop("no distance column in ", panel$primary_source)
  d$value <- d[[dcol]]
  d$lab <- nv_dataset_label(d$dataset)
  d$lab <- factor(d$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))
  d$Group <- factor(d$StressGroup, levels = c("CON", "RES", "SUS"))
  # deterministic horizontal offset instead of geom_jitter
  d <- d[order(d$lab, d$Group, d$AnimalID), ]
  d$off <- (stats::ave(seq_len(nrow(d)), paste(d$lab, d$Group),
                       FUN = seq_along) - 2) * 0.16
  d$xpos <- as.integer(d$Group) + d$off
  p <- ggplot2::ggplot(d, ggplot2::aes(xpos, value)) +
    ggplot2::stat_summary(ggplot2::aes(x = as.integer(Group)), fun = mean,
                          geom = "crossbar", width = 0.5, linewidth = 0.2,
                          colour = "grey55") +
    ggplot2::geom_point(ggplot2::aes(colour = Group), size = 1.1) +
    ggplot2::scale_colour_manual(values = grp, guide = "none") +
    ggplot2::scale_x_continuous(breaks = 1:3, labels = c("CON", "RES", "SUS"),
                                limits = c(0.5, 3.5)) +
    ggplot2::facet_wrap(~lab, nrow = 1, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "network distance from CON consensus") +
    nv_theme(grid = "y")
  d$reading <- paste0("each point is one animal-level network replicate; CON ",
                      "animals use a leave-one-CON-out centroid so no animal ",
                      "is compared against a centroid containing itself; ",
                      "group means are shown without an inferential claim")
  write_csv_safe(d, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# PCA with a DETERMINISTIC legend order.
#
# svp_pca (Part-17 infrastructure, reused unchanged by Parts 18 and 19) gives
# both the colour and the shape scale name = NULL. ggplot2 then has no
# tie-breaker for guide order, so the dataset legend and the region legend swap
# places at random between runs: repeated renders of the identical input
# alternate between exactly two SVGs. That makes the panel, and any figure
# containing it, non-reproducible.
#
# Parts 16-19 are frozen for this task, so this is NOT a patch to svp_pca. It is
# a v5-local renderer that draws the same data from the same sources with the
# guide order pinned. The frozen layers keep their existing behaviour; the
# defect is reported rather than silently propagated.
s5_pca <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  ve <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]),
                    required = FALSE)
  pc <- function(k) {
    if (is.null(ve)) return(sprintf("PC%d", k))
    v <- ve$variance_explained[ve$PC == paste0("PC", k)]
    if (!length(v)) return(sprintf("PC%d", k))
    sprintf("PC%d (%.0f%%)", k, 100 * v[1])
  }
  ds <- c("neuron_neuropil", "neuron_soma", "microglia")
  d <- d[d$dataset %in% names(nv_dataset_colours()), , drop = FALSE]
  d$lab <- factor(nv_dataset_label(d$dataset), levels = nv_dataset_label(ds))
  d$region <- factor(d$region, levels = c("CA1", "CA2", "CA3", "DG"))
  # stable row order so the drawing order of overlapping points is fixed too
  d <- d[order(d$lab, d$region, d$Sample), , drop = FALSE]

  p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_point(ggplot2::aes(colour = lab, shape = region),
                        size = 0.6, stroke = 0.25, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      unname(nv_dataset_colours()[ds]), nv_dataset_label(ds)), name = NULL) +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 3), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
                    shape = ggplot2::guide_legend(order = 2)) +
    ggplot2::labs(x = pc(1), y = pc(2)) +
    nv_theme() +
    ggplot2::theme(legend.position = "right",
                   legend.spacing.y = ggplot2::unit(0.4, "mm"),
                   legend.key.height = ggplot2::unit(2.2, "mm"))

  write_csv_safe(d[, c("Sample", "PC1", "PC2", "dataset", "region", "layer",
                       "celltype_layer")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
