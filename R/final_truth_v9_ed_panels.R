# Nature-final v7 Extended Data panels.
#
# Extended Data DEEPENS the main figures. A main panel and its ED counterpart
# must never be the same graphic at two sizes, so each renderer here carries
# information the main figure deliberately omits.
#
# DOWNSTREAM ONLY. No model, no test, no p-value, no new inference. 5 pt floor.

# ==========================================================================
# ED1: the FULL bilateral audit
# ==========================================================================
#
# Main Figure 2f shows one number per contrast (the left-right correlation) and
# the regional-versus-laminar split. ED1 shows the complete stored metric
# inventory: Pearson, Spearman, sign agreement and the pair count, for every
# prespecified contrast. Different information, not a bigger copy.
f9_ed_bilateral_full <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  su <- nv_read_csv(repo_path(panel$primary_source))
  su$level <- ifelse(grepl("_strata$|_DG_layers$", su$contrast),
                     "Fine / laminar identity", "Regional identity")
  su$level <- factor(su$level,
                     levels = c("Regional identity", "Fine / laminar identity"))
  su <- su[order(su$level, -su$pearson_r), , drop = FALSE]
  su$lab <- sub(" vs mean other", " vs rest", gsub("_", " ", su$contrast))
  su$ypos <- rev(seq_len(nrow(su)))

  metrics <- c(pearson_r = "Pearson r", spearman_rho = "Spearman rho",
               sign_agreement_fraction = "sign agreement")
  metrics <- metrics[names(metrics) %in% names(su)]
  long <- do.call(rbind, lapply(names(metrics), function(m) data.frame(
    ypos = su$ypos, lab = su$lab, level = su$level,
    metric = factor(unname(metrics[m]), levels = unname(metrics)),
    value = as.numeric(su[[m]]), stringsAsFactors = FALSE)))

  p <- ggplot2::ggplot(long, ggplot2::aes(value, ypos)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = value, yend = ypos,
                                       colour = level),
                          linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(colour = level), size = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)),
                       hjust = -0.25, family = fam, size = nf_sz(5.0),
                       colour = "grey25") +
    ggplot2::scale_colour_manual(
      values = c("Regional identity" = "#1F3D52",
                 "Fine / laminar identity" = "#C2A878"), name = NULL) +
    ggplot2::scale_y_continuous(breaks = su$ypos, labels = su$lab,
                                expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::scale_x_continuous(limits = c(0, 1.28), breaks = c(0, 0.5, 1)) +
    ggplot2::facet_wrap(~metric, nrow = 1) +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    "Complete stored bilateral metric inventory for every ",
                    "prespecified anatomical contrast. Main Figure 2f shows ",
                    "only the left-right correlation and the ",
                    "regional-versus-laminar split.")) +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "bottom",
                   legend.key.size = ggplot2::unit(2.4, "mm"),
                   panel.spacing.x = ggplot2::unit(2.2, "mm"))
  nreg <- sum(su$level == "Regional identity")
  if (nreg > 0 && nreg < nrow(su)) {
    yb <- su$ypos[nreg] - 0.5
    p <- p + ggplot2::annotate("segment", x = -Inf, xend = Inf, y = yb, yend = yb,
                               linewidth = 0.3, colour = "grey55")
  }
  write_csv_safe(su, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# ED2: the COMPLETE external validation inventory
# ==========================================================================
#
# Main Figure 2g shows only the expected pairings. ED2 shows every tested pair,
# expected and specificity alike, so a reader can see that each contrast
# recovers its own signature and NOT the others.
f9_ed_external_full <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  k <- nv_read_csv(repo_path(panel$primary_source))
  k$level <- ifelse(grepl("strata", k$validation_domain),
                    "CA1 laminar identity", "Regional identity")
  k$kind <- ifelse(k$expected_match %in% TRUE, "expected pairing",
                   "specificity comparison")
  k$sig <- is.finite(k$p_adjust) & k$p_adjust < 0.05
  k$ic <- gsub("_", " ", k$internal_contrast)
  ord <- unique(k$ic[order(k$level, k$ic)])
  k$ypos <- match(k$ic, rev(ord))
  k$xpos <- match(k$external_signature, sort(unique(k$external_signature)))
  lim <- max(abs(k$NES), na.rm = TRUE) * c(-1, 1)

  p <- ggplot2::ggplot(k, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_point(ggplot2::aes(fill = NES, size = kind, shape = kind,
                                     colour = sig)) +
    ggplot2::scale_fill_gradient2(low = nv_palette()$diverging$low,
                                  mid = nv_palette()$diverging$mid,
                                  high = nv_palette()$diverging$high,
                                  midpoint = 0, limits = lim, name = "Normalised\nenrichment score",
                                  guide = ggplot2::guide_colourbar(order = 2)) +
    ggplot2::scale_shape_manual(values = c("expected pairing" = 21,
                                           "specificity comparison" = 22),
                                name = NULL,
                                guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_size_manual(values = c("expected pairing" = 2.4,
                                          "specificity comparison" = 1.5),
                               name = NULL,
                               guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "grey10", "FALSE" = "grey70"),
                                 guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(sort(unique(k$external_signature))),
      labels = sort(unique(k$external_signature)), expand = c(0, 0.6)) +
    ggplot2::scale_y_continuous(breaks = seq_along(ord), labels = rev(ord),
                                expand = c(0, 0.6)) +
    ggplot2::labs(x = "external hippocampal signature", y = NULL,
                  caption = paste0(
                    "Complete validation inventory: every tested internal ",
                    "contrast against every external signature. Circles are ",
                    "expected pairings, squares are specificity comparisons; ",
                    "dark outline = FDR < 0.05.\nOnly structurally applicable ",
                    "pairs are drawn. Main Figure 2g shows the expected ",
                    "pairings alone.")) +
    nf_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "right",
                   legend.key.size = ggplot2::unit(2.6, "mm"))
  write_csv_safe(k[, c("internal_contrast", "external_signature", "level",
                       "kind", "NES", "p_adjust")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# ED2: the COMPLETE internal anatomical program inventory
# ==========================================================================
f9_ed_internal_full <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  g <- nv_read_csv(repo_path(panel$primary_source))
  g$level <- ifelse(grepl("CA1_strata", g$contrast),
                    "CA1 laminar identity", "Regional identity")
  g$cl <- sub(" vs mean other", " vs rest", gsub("_", " ", g$contrast))
  g$sig <- is.finite(g$p_adjust) & g$p_adjust < 0.05

  # Part-25 fidelity pass. The panel faceted 14 points into 7 boxes of 57 mm
  # each, so every NES sat on its own axis at its own offset and only the bottom
  # facet of a column carried tick labels: comparing CA1-SLM against DG meant
  # eye-jumping across two columns. One shared axis puts all 14 on the same
  # scale and raises quantitative resolution 2.2x (17.3 -> 38.1 mm per NES
  # unit). It also matches ED1a, which draws the identical data shape - one
  # scalar per labelled row, grouped by a categorical level - as a single-axis
  # multi-row lollipop, so this Extended Data set now teaches one idiom.
  #
  # The shape aesthetic is dropped: max p_adjust over all 14 rows is 8.37e-6, so
  # it never varied and the FDR ceiling is stated once in the caption instead.
  g <- g[order(g$level, g$cl, -g$NES), , drop = FALSE]
  g$ypos <- rev(seq_len(nrow(g)))
  # contrast blocks stay contiguous, so a hairline can divide them at zero
  # row cost
  r <- rle(g$cl)
  ends <- cumsum(r$lengths)
  divs <- g$ypos[utils::head(ends, -1)] - 0.5
  ctag <- data.frame(y = vapply(split(g$ypos, factor(g$cl, levels = r$values)),
                                mean, numeric(1)),
                     l = r$values, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(g, ggplot2::aes(NES, ypos)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey65") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = ypos,
                                       colour = level),
                          linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(size = setSize, colour = level)) +
    ggplot2::geom_text(ggplot2::aes(label = Description), x = -0.12, hjust = 1,
                       family = fam, size = nf_sz(5.0), colour = "grey25") +
    ggplot2::geom_text(data = ctag, ggplot2::aes(x = 3.62, y = y, label = l),
                       inherit.aes = FALSE, hjust = 0, family = fam,
                       size = nf_sz(5.0), colour = "grey45") +
    ggplot2::scale_colour_manual(
      values = c("Regional identity" = "#1F3D52",
                 "CA1 laminar identity" = "#C2A878"), name = NULL,
      guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_size_continuous(range = c(0.6, 2.2), name = "Gene set size",
                                   guide = ggplot2::guide_legend(order = 2)) +
    ggplot2::scale_x_continuous(limits = c(-2.05, 4.5),
                                breaks = c(0, 1, 2, 3)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::labs(x = "NES", y = NULL,
                  caption = paste0(
                    "Every canonical GO term retained for every anatomical ",
                    "contrast, on one shared axis; the right-hand label names ",
                    "the contrast. All 14 terms are FDR < 1e-5. Main Figure 2h ",
                    "shows one term per contrast.")) +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "bottom",
                   legend.key.size = ggplot2::unit(2.4, "mm"))
  if (length(divs)) p <- p + ggplot2::annotate(
    "segment", x = -2.05, xend = 4.5, y = divs, yend = divs,
    linewidth = 0.12, colour = "grey86")
  out <- g[, c("contrast", "Description", "level", "NES", "p_adjust", "setSize")]
  out$axis_note <- paste0(
    "all terms share one NES axis so any two can be compared directly; the ",
    "previous facet grid gave each contrast its own axis")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# ED_WGCNA_COMBINED: module identity, spatial fingerprint and phenotype
# ==========================================================================
#
# ED4 and ED5 are merged. The WGCNA circle is dropped: the module x spatial-unit
# fingerprint in the same figure shows the actual profiles, so the radial layout
# encoded nothing the fingerprint does not encode better.
f9_ed_wgcna_phenotype <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  z <- nv_read_csv(repo_path(panel$primary_source))
  pick <- function(cands, what) {
    for (nm in cands) {
      if (!nm %in% names(z)) next
      v <- z[[nm]]
      if (any(!is.na(v) & nzchar(as.character(v)))) return(nm)
    }
    stop("f9_ed_wgcna_phenotype: no populated ", what, " column", call. = FALSE)
  }
  mcol <- pick(c("module_id", "ModuleID", "endpoint_id"), "module")
  ccol <- pick(c("contrast", "Contrast"), "contrast")
  vcol <- pick(c("estimate", "effect", "logFC", "value"), "effect")
  z$mid <- sub("^WGCNA_", "", as.character(z[[mcol]]))
  z$con <- as.character(z[[ccol]])
  z$val <- as.numeric(z[[vcol]])
  z <- z[!is.na(z$val), , drop = FALSE]

  wa <- nv_read_csv(repo_path("results", "tables", "11_spatial_systems", "atlas",
                              "WGCNA_module_spatial_cell_affinity.csv"))
  ann <- data.frame(mid = sub("^WGCNA_", "", as.character(wa$ModuleID)),
                    ds = as.character(wa$dataset),
                    peak = as.character(wa$peak_unit),
                    cell = as.character(wa$external_celltype_all),
                    stringsAsFactors = FALSE)
  ds <- unique(as.character(z$dataset))
  if (length(ds) == 1L && ds %in% ann$ds) ann <- ann[ann$ds == ds, , drop = FALSE]
  ann <- ann[!duplicated(ann$mid), , drop = FALSE]
  z <- merge(z, ann, by = "mid", all.x = TRUE, sort = FALSE)
  z$peak_label <- ""
  ok <- !is.na(z$peak) & nzchar(z$peak) & !is.na(z$ds)
  if (any(ok)) z$peak_label[ok] <- sg_unit_label(z$peak[ok], z$ds[ok])

  # the descriptive geometry, stated numerically rather than implied
  w <- reshape(unique(z[, c("mid", "con", "val")]), idvar = "mid",
               timevar = "con", direction = "wide")
  names(w) <- sub("^val[.]", "", names(w))
  geom_n <- sum(w[["RES - CON"]] > 0 & w[["SUS - CON"]] < 0 & w[["SUS - RES"]] < 0,
                na.rm = TRUE)

  mods <- unique(z[, c("mid", "peak_label", "cell")])
  mods <- mods[order(mods$mid), , drop = FALSE]
  mods$row_label <- ifelse(nzchar(mods$peak_label),
                           paste0(mods$mid, "  peak ", mods$peak_label), mods$mid)
  z$ypos <- match(z$mid, rev(mods$mid))
  cl <- c("RES - CON", "SUS - CON", "SUS - RES")
  cl <- cl[cl %in% z$con]
  z <- z[z$con %in% cl, , drop = FALSE]
  z$xpos <- match(z$con, cl)
  ny <- nrow(mods)
  lim <- max(abs(z$val), na.rm = TRUE) * c(-1, 1)

  # Brief section 24. The descriptive count and the inferential null were
  # previously two clauses of one wrapped caption sentence, so a reader could
  # take in "13 of 15" without ever reaching "0 of 45". They are now two chips
  # side by side at the top of the panel, read as one statement.
  fcol <- intersect(c("tier_specific_fdr", "q_value", "FDR_global"), names(z))
  n_cells <- nrow(unique(z[, c("mid", "con")]))
  n_sig <- if (length(fcol)) sum(suppressWarnings(as.numeric(z[[fcol[1]]])) <
                                   0.05, na.rm = TRUE) else NA_integer_
  if (is.na(n_sig)) stop("f9_ed_wgcna_phenotype: no FDR column", call. = FALSE)
  # Part-25. The 0/45 chip above is derived from a SPATIALLY ADJUSTED model, and
  # the obvious objection is that a layer-confined module effect would be
  # averaged away by that adjustment. The pipeline already answers this: stage
  # 05 fits a stress-by-spatial-unit interaction omnibus per module. Those rows
  # carry display_allowed = TRUE and model_valid_for_inference = TRUE, so
  # reporting them is a read, not a new test. Nothing is recomputed here.
  #
  # The 690 within-spatial-unit contrasts in the same file are deliberately NOT
  # reported: they are three-animals-per-group-per-unit slices whose smallest
  # FDR is 0.49, and calling that evidence of absence would need an attainable-
  # floor treatment of the kind ED8c carries - which would be new inference.
  omni <- do.call(rbind, lapply(sg_compartments()$id, function(d) {
    f <- repo_path("results", "tables", "06_modules_WGCNA", "group_effects", d,
                   "module_group_effects.csv")
    if (!file.exists(f)) return(NULL)
    e <- nv_read_csv(f)
    e <- e[e$effect_scope == "stress_by_spatial_interaction", , drop = FALSE]
    if (!nrow(e)) return(NULL)
    q <- suppressWarnings(as.numeric(e$FDR_interaction_omnibus))
    data.frame(dataset = d, n_tests = nrow(e),
               n_fdr = sum(q < 0.05, na.rm = TRUE),
               min_fdr = min(q, na.rm = TRUE), stringsAsFactors = FALSE)
  }))
  if (is.null(omni)) stop("f9_ed_wgcna_phenotype: no interaction omnibus rows",
                          call. = FALSE)
  om_n <- sum(omni$n_tests); om_sig <- sum(omni$n_fdr)
  om_min <- min(omni$min_fdr)
  # Part-26. The three counts stay adjacent, as asked in Part 24 section 24,
  # but the badge chrome goes: no filled rectangles and no 6.4 pt bold ratios.
  # Three restrained lines carry the descriptive count and both nulls, divided
  # from the field by a hairline. The two inferential lines are set in the
  # accent colour so they cannot be skipped, which was the point of the chips.
  stat_line <- data.frame(
    y = c(ny + 3.35, ny + 2.55, ny + 1.75),
    l = c(sprintf("%d of %d modules show the descriptive pattern",
                  geom_n, nrow(w)),
          sprintf("%d of %d module × contrast cells are FDR-supported",
                  n_sig, n_cells),
          sprintf(paste0("%d of %d stress × spatial-unit interaction tests ",
                         "are FDR-supported (smallest FDR %.2f)"),
                  om_sig, om_n, om_min)),
    ink = c("grey35", "#B03A24", "#B03A24"),
    stringsAsFactors = FALSE)
  p <- ggplot2::ggplot(z, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = val), colour = "white",
                       linewidth = 0.16) +
    nv_diverging(limits = lim,
                 name = "Module eigengene
difference") +
    ggplot2::scale_x_continuous(breaks = seq_along(cl), labels = cl,
                                expand = c(0, 0),
                                limits = c(0.5, length(cl) + 3.4)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny),
                                labels = rev(mods$row_label),
                                expand = c(0, 0),
                                limits = c(0.5, ny + 3.95)) +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    "DESCRIPTIVE ONLY. The colour scale must not be read",
                    " as evidence of a group difference.")) +
    nf_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
                   axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "right",
                   legend.key.width = ggplot2::unit(1.6, "mm"),
                   legend.key.height = ggplot2::unit(3.4, "mm"))
  p <- p +
    ggplot2::annotate("segment", x = 0.5, xend = length(cl) + 3.38,
                      y = ny + 1.25, yend = ny + 1.25, linewidth = 0.18,
                      colour = "grey80") +
    ggplot2::geom_text(data = stat_line,
                       ggplot2::aes(x = 0.55, y = y, label = l),
                       inherit.aes = FALSE, hjust = 0, family = fam,
                       size = nf_sz(5.0), colour = stat_line$ink)
  tr <- unique(z[, c("mid", "ypos", "cell")])
  p <- p +
    ggplot2::geom_text(data = tr,
                       ggplot2::aes(x = length(cl) + 0.8, y = ypos,
                                    label = ifelse(is.na(cell), "", cell)),
                       inherit.aes = FALSE, hjust = 0, family = fam,
                       size = nf_sz(5.0), colour = "grey35") +
    ggplot2::annotate("text", x = length(cl) + 0.8, y = ny + 0.72, hjust = 0,
                      family = fam, size = nf_sz(5.0), fontface = "bold",
                      colour = "grey20", label = "external cell type")
  z$inferential_status <- sprintf(
    paste0("%d of %d modules show the descriptive geometry; %d of %d ",
           "spatially adjusted cells FDR-supported; %d of %d stress x ",
           "spatial-unit interaction omnibus tests FDR-supported, smallest ",
           "FDR %.2f"),
    geom_n, nrow(w), n_sig, n_cells, om_sig, om_n, om_min)
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# m11 worked example: spatial profile, external affinity, bilateral class.
# The historical label is NOT renamed and the proposed label is NOT activated.
f9_ed_m11 <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  z <- s6_module_fingerprint_table()
  z <- z[z$dataset == "neuron_neuropil" & sub("^WGCNA_", "", z$ModuleID) == "m11" &
           !is.na(z$mean_con_z), , drop = FALSE]
  if (!nrow(z)) stop("f9_ed_m11: no m11 rows", call. = FALSE)
  z$sg_unit <- sg_resolve_unit(z$spatial_unit, z$dataset)
  u <- sg_units(); u <- u[u$dataset == "neuron_neuropil", , drop = FALSE]
  z$xpos <- match(z$sg_unit, u$analysis_key)
  z <- z[order(z$xpos), , drop = FALSE]
  lab <- sub("^WGCNA_", "", z$ModuleID[1])
  peak <- u$display[which.max(z$mean_con_z)]
  cell <- as.character(z$external_celltype_all[1])
  bil <- as.character(z$bilateral_reproducibility_class[1])

  # Brief section 25. Restrained annotation: the peak is marked on the bar that
  # actually is the peak instead of being described in prose, and the label
  # status moves to the title, where it cannot be skipped. The active historical
  # label is unchanged and the proposed oligodendrocyte / myelin label is NOT
  # activated - the external evidence is shown as context only.
  ipk <- which.max(z$mean_con_z)
  pk <- data.frame(x = z$xpos[ipk], y = z$mean_con_z[ipk], stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(z, ggplot2::aes(xpos, mean_con_z)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey70") +
    ggplot2::geom_col(fill = "#C9D6DF", width = 0.72) +
    ggplot2::geom_col(data = z[ipk, , drop = FALSE], fill = "#4E7288",
                      width = 0.72) +
    ggplot2::geom_text(data = pk, ggplot2::aes(x, y, label = "peak"),
                       vjust = -0.45, family = fam, size = nf_sz(5.0),
                       colour = "#4E7288", fontface = "bold") +
    ggplot2::scale_x_continuous(breaks = seq_len(nrow(u)), labels = u$display,
                                expand = c(0, 0.6)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
    ggplot2::labs(
      x = NULL, y = "mean CON z",
      title = paste0("neuropil ", lab,
                     " worked example \u2014 active historical label, unchanged"),
      caption = paste0(
        "Peak ", peak, "  \u00b7  external cell-type affinity ",
        ifelse(is.na(cell), "not assigned", cell), "  \u00b7  bilateral ",
        ifelse(is.na(bil), "not assigned", bil),
        ". The proposed oligodendrocyte / myelin label is NOT activated; the ",
        "external evidence is context only.")) +
    nf_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = NF_MIN_PT,
                                                       angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(size = nf_pt(5.4),
                                                      face = "bold"))
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---- ED8 a. three square-celled compartment blocks -----------------------
#
# Brief section 23. The previous geometry used facet_grid(. ~ comp), which
# shares ONE y scale across the three facets: the soma and microglia matrices
# were therefore 4-column strips with 14 blank rows each, and the three blocks
# could not be compared at a glance.
#
# The replacement draws all three blocks in ONE fixed-aspect coordinate space,
# each with its own local row and column indices. coord_fixed then guarantees
# that a cell is the same physical square in all three compartments, so the
# 10-unit neuropil block and the two 4-region blocks are directly comparable.
# The unit count is stated in each block header, and the "not connectivity"
# disclaimer is the panel title.
f9_ed_similarity <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  e <- nv_read_csv(repo_path(panel$primary_source))
  e$a <- sg_resolve_unit(e$node_a, e$dataset)
  e$b <- sg_resolve_unit(e$node_b, e$dataset)
  sym <- rbind(
    data.frame(dataset = e$dataset, r = e$a, c = e$b, v = e$median_similarity,
               stringsAsFactors = FALSE),
    data.frame(dataset = e$dataset, r = e$b, c = e$a, v = e$median_similarity,
               stringsAsFactors = FALSE))
  u <- sg_units()
  cs <- sg_compartments()

  GUT <- 2.7   # row-label gutter, in cell widths
  GAP <- 1.7   # gap between blocks, in cell widths
  blocks <- list(); xoff <- 0
  for (k in seq_len(nrow(cs))) {
    d <- cs$id[k]
    uk <- u[u$dataset == d, , drop = FALSE]
    blocks[[d]] <- list(dataset = d, units = uk$analysis_key,
                        display = uk$display, n = nrow(uk), x0 = xoff + GUT,
                        title = cs$display[k],
                        count = sprintf("%d %s", nrow(uk),
                                        if (identical(cs$resolution[k],
                                                      "region_x_layer"))
                                          "spatial units" else "regions"))
    xoff <- xoff + GUT + nrow(uk) + GAP
  }
  nmax <- max(vapply(blocks, function(b) b$n, numeric(1)))

  tile <- do.call(rbind, lapply(blocks, function(b) {
    z <- sym[sym$dataset == b$dataset, , drop = FALSE]
    data.frame(x = b$x0 + match(z$c, b$units), y = -match(z$r, b$units),
               v = z$v, dataset = b$dataset, unit_row = z$r, unit_col = z$c,
               stringsAsFactors = FALSE)
  }))
  rlab <- do.call(rbind, lapply(blocks, function(b) data.frame(
    x = b$x0 + 0.35, y = -seq_len(b$n), l = b$display, stringsAsFactors = FALSE)))
  clab <- do.call(rbind, lapply(blocks, function(b) data.frame(
    x = b$x0 + seq_len(b$n), y = -b$n - 0.65, l = b$display,
    stringsAsFactors = FALSE)))
  hdr <- do.call(rbind, lapply(blocks, function(b) data.frame(
    x = b$x0 - GUT + 0.2, t = b$title, c = b$count, stringsAsFactors = FALSE)))

  lim <- max(abs(tile$v), na.rm = TRUE) * c(-1, 1)
  x_end <- 30.0
  ybot <- -(nmax + 0.65 + 1.85)
  # the notes column: line pitch is set in cell units, so it stays tight
  # whatever the block size
  # The metric and the zero reference are common to all three blocks, so ONE
  # global zero-centred diverging scale is retained deliberately. Per-block
  # rescaling would make an identical colour mean a different correlation in
  # each block and would destroy exactly the cross-compartment comparison the
  # panel exists to support. The consequence - that the two region-level blocks
  # occupy a narrower part of the range than neuropil - is a property of the
  # data and is stated here rather than engineered away.
  notes <- data.frame(
    l = c("A value is the correlation between",
          "two units' protein profiles, not a",
          "projection.",
          "CON animals only (n = 3); median",
          "across animals. Diagonal omitted.",
          "Rows and columns of a block carry",
          "the same units in the same order.",
          "One scale centred on zero for all",
          "three blocks: soma and microglia",
          "span a narrower range than",
          "neuropil, which is a property of",
          "the data, not of the scale."),
    y = c(0, 0.62, 1.24, 2.20, 2.82, 3.78, 4.40, 5.36, 5.98, 6.60, 7.22, 7.84),
    stringsAsFactors = FALSE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = tile, ggplot2::aes(x, y, fill = v),
                       colour = "white", linewidth = 0.12) +
    nv_diverging(limits = lim, name = "Median profile\ncorrelation (CON)") +
    ggplot2::geom_text(data = rlab, ggplot2::aes(x, y, label = l), hjust = 1,
                       family = fam, size = nf_sz(5.0), colour = "grey30") +
    ggplot2::geom_text(data = clab, ggplot2::aes(x, y, label = l), hjust = 1,
                       angle = 90, family = fam, size = nf_sz(5.0),
                       colour = "grey30") +
    ggplot2::geom_text(data = hdr, ggplot2::aes(x, 1.45, label = t), hjust = 0,
                       family = fam, size = nf_sz(5.4), fontface = "bold",
                       colour = "grey12") +
    ggplot2::geom_text(data = hdr, ggplot2::aes(x, 0.52, label = c), hjust = 0,
                       family = fam, size = nf_sz(5.0), colour = "grey40") +
    ggplot2::geom_text(data = notes, ggplot2::aes(x = x_end + 1.2,
                                                y = 0.55 - y, label = l),
                       hjust = 0, family = fam, size = nf_sz(5.0),
                       colour = "grey35") +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 41.5), ylim = c(ybot, 2.1),
                         expand = FALSE) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0("Molecular similarity between spatial proteomic profiles",
                     " \u2014 NOT anatomical connectivity")) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(
      text = ggplot2::element_text(family = fam, size = NF_MIN_PT),
      legend.text = ggplot2::element_text(size = NF_MIN_PT),
      legend.title = ggplot2::element_text(size = nf_pt(5.4)),
      plot.title = ggplot2::element_text(size = nf_pt(5.8), face = "bold",
                                         colour = "grey10", hjust = 0),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.6, "mm"),
      legend.key.height = ggplot2::unit(3.6, "mm"),
      plot.margin = ggplot2::margin(0.5, 1, 0.5, 1, "mm"))
  out <- tile[, c("dataset", "unit_row", "unit_col", "v")]
  names(out)[4] <- "median_similarity"
  out$reading <- paste0(
    "similarity between spatial proteomic profiles, NOT anatomical ",
    "connectivity; one cell is the same physical square in all three blocks")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# v7 re-renders of four panels the QA gate caught below the 5 pt floor
# ==========================================================================
#
# The Part-21 originals sit at 4.60-4.82 pt. Parts 16-22 are frozen, so these
# are new v7 renderers over the SAME data helpers, with nf_theme applied and
# the row/column counts trimmed where 5 pt genuinely does not fit.

# ---- full baseline spatial fingerprint (external signature scores) -------
f9_ed_fingerprint_full <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  z <- nv_read_csv(repo_path(panel$primary_source))
  keep <- c("Kaulich CA1", "Kaulich CA2/3", "Kaulich DG", "Kaulich SO",
            "Kaulich SP", "Kaulich SR", "Kaulich SLM", "Excitatory neuron",
            "Inhibitory interneuron", "Synaptic neuropil", "Astrocyte",
            "Oligodendrocyte", "Microglia / PVM", "Vascular")
  z <- z[z$signature %in% keep, , drop = FALSE]
  z$row_label <- sub("^Kaulich ", "", z$signature)
  lev <- sub("^Kaulich ", "", keep)
  b <- sg_blocks(z$sg_unit, z$dataset)
  o <- b$order
  z$xpos <- match(paste(z$dataset, z$sg_unit), paste(o$dataset, o$unit))
  z$ypos <- match(z$row_label, rev(lev))
  n <- nrow(o); ny <- length(lev)
  y_comp <- ny + 2.3; y_reg <- ny + 1.05
  lim <- max(abs(z$score), na.rm = TRUE) * c(-1, 1)

  p <- ggplot2::ggplot(z, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = score), colour = "white",
                       linewidth = 0.1) +
    nv_diverging(limits = lim, name = "Signature score\n(CON z-score)") +
    ggplot2::scale_x_continuous(breaks = seq_len(n), labels = sg_axis_labels(b),
                                limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny), labels = rev(lev),
                                limits = c(0.5, y_comp + 0.9), expand = c(0, 0)) +
    ggplot2::annotate("segment", x = b$compartment$start - 0.5,
                      xend = b$compartment$end + 0.5, y = y_comp - 0.32,
                      yend = y_comp - 0.32, linewidth = 0.45, colour = "grey25") +
    ggplot2::annotate("text", x = b$compartment$mid, y = y_comp,
                      label = b$compartment$label, family = fam,
                      size = nf_sz(5.4), fontface = "bold", colour = "grey12") +
    ggplot2::annotate("text", x = b$region$mid, y = y_reg, label = b$region$label,
                      family = fam, size = nf_sz(5.0), colour = "grey30") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    "CON only. Rows are external prespecified signatures; each ",
                    "protein is standardised within its own compartment. Main ",
                    "Figure 2d shows named proteins from CON-only contrasts.")) +
    nf_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT, colour = "grey25"),
      axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.6, "mm"),
      legend.key.height = ggplot2::unit(3.4, "mm"))
  cend <- utils::head(b$compartment$end, -1)
  if (length(cend)) p <- p + ggplot2::annotate(
    "segment", x = cend + 0.5, xend = cend + 0.5, y = 0.5, yend = y_comp - 0.32,
    linewidth = 0.42, colour = "grey25")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---- module x spatial-unit fingerprint -----------------------------------
f9_ed_module_fingerprint <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  z <- s6_module_fingerprint_table()
  z <- z[!is.na(z$mean_con_z), , drop = FALSE]
  z$sg_unit <- sg_resolve_unit(z$spatial_unit, z$dataset)
  z$mod <- sub("^WGCNA_", "", z$ModuleID)
  u <- sg_units()
  cs <- sg_compartments()
  disp <- stats::setNames(u$display, u$analysis_key)

  # Part-25 fidelity pass. facet_wrap(ncol = 1, scales = "free") gave all three
  # compartment blocks the SAME panel height and the SAME width whatever their
  # row and column counts, so tile aspect was 11.8:1 (neuropil), 13.8:1 (soma)
  # and 25.5:1 (microglia) - a soma module row was 2.14x taller than a neuropil
  # one for no reason, and 160 mm of width was spent on 4 columns in two of the
  # three blocks. Drawing all three blocks in ONE coordinate space with local
  # row and column indices makes a tile a constant 1 x 1 data unit, so it is the
  # same physical rectangle everywhere.
  #
  # Rows were also ordered alphabetically by module id, which is a random
  # permutation of the spatial structure the panel exists to show. They are now
  # seriated by peak_unit, which is already a column of the same table, so each
  # module's spatial preference reads as a block instead of needing a lookup.
  GUT <- 0.85   # row-label gutter, in tile widths
  GAP <- 0.90   # gap between blocks, in tile widths
  blocks <- list(); xoff <- 0
  for (k in seq_len(nrow(cs))) {
    d <- cs$id[k]
    uk <- u[u$dataset == d, , drop = FALSE]
    zk <- z[z$dataset == d, , drop = FALSE]
    if (!nrow(zk)) next
    pk <- stats::setNames(zk$peak_unit[!duplicated(zk$mod)],
                          zk$mod[!duplicated(zk$mod)])
    ord <- names(pk)[order(match(sg_resolve_unit(unname(pk), d),
                                 uk$analysis_key), names(pk))]
    blocks[[d]] <- list(dataset = d, units = uk$analysis_key,
                        display = uk$display, mods = ord, x0 = xoff + GUT,
                        title = cs$display[k])
    xoff <- xoff + GUT + nrow(uk) + GAP
  }
  nrow_max <- max(vapply(blocks, function(b) length(b$mods), numeric(1)))

  tile <- do.call(rbind, lapply(blocks, function(b) {
    zk <- z[z$dataset == b$dataset, , drop = FALSE]
    data.frame(x = b$x0 + match(zk$sg_unit, b$units),
               y = -match(zk$mod, b$mods), v = zk$mean_con_z,
               dataset = b$dataset, module = zk$mod, unit = zk$sg_unit,
               peak_unit = zk$peak_unit, stringsAsFactors = FALSE)
  }))
  tile <- tile[!is.na(tile$x) & !is.na(tile$y), , drop = FALSE]
  rlab <- do.call(rbind, lapply(blocks, function(b) data.frame(
    x = b$x0 + 0.35, y = -seq_along(b$mods), l = b$mods,
    stringsAsFactors = FALSE)))
  clab <- do.call(rbind, lapply(blocks, function(b) data.frame(
    x = b$x0 + seq_along(b$units), y = -length(b$mods) - 0.7,
    l = b$display, stringsAsFactors = FALSE)))
  hdr <- do.call(rbind, lapply(blocks, function(b) data.frame(
    x = b$x0 - GUT + 0.15, l = b$title, stringsAsFactors = FALSE)))

  lim <- max(abs(tile$v), na.rm = TRUE) * c(-1, 1)
  ybot <- -(nrow_max + 0.7 + 1.9)

  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = tile, ggplot2::aes(x, y, fill = v),
                       colour = "white", linewidth = 0.1) +
    nv_diverging(limits = lim,
                 name = "Mean module-member\nabundance (CON z-score)") +
    ggplot2::geom_text(data = rlab, ggplot2::aes(x, y, label = l), hjust = 1,
                       family = fam, size = nf_sz(5.0), colour = "grey30") +
    ggplot2::geom_text(data = clab, ggplot2::aes(x, y, label = l), hjust = 1,
                       angle = 90, family = fam, size = nf_sz(5.0),
                       colour = "grey30") +
    ggplot2::geom_text(data = hdr, ggplot2::aes(x, 0.75, label = l), hjust = 0,
                       family = fam, size = nf_sz(5.4), fontface = "bold",
                       colour = "grey12") +
    ggplot2::coord_cartesian(xlim = c(0, xoff - GAP + 0.5),
                             ylim = c(ybot, 1.35), expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    "CON only. Each cell is the mean within-protein CON z of ",
                    "that module's member proteins. A module exists only in ",
                    "its own compartment, so each block carries only that ",
                    "compartment's units.
A tile is the same size in all ",
                    "three. Modules are ordered by their peak unit.")) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(
      text = ggplot2::element_text(family = fam, size = NF_MIN_PT),
      legend.text = ggplot2::element_text(size = NF_MIN_PT),
      legend.title = ggplot2::element_text(size = nf_pt(5.4)),
      plot.caption = ggplot2::element_text(size = NF_MIN_PT, colour = "grey30",
                                           hjust = 0, lineheight = 1.2),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.6, "mm"),
      legend.key.height = ggplot2::unit(3.4, "mm"),
      plot.margin = ggplot2::margin(0.5, 1, 0.5, 1, "mm"))
  out <- tile[, c("dataset", "module", "unit", "peak_unit", "v")]
  names(out)[5] <- "mean_con_z"
  out$layout_note <- paste0(
    "one coordinate space for all three compartments, so a tile is a constant ",
    "1 x 1 data unit everywhere; modules are seriated by peak_unit")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---- baseline location -> stress-effect location -------------------------
f9_ed_locations <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  h <- s6_ed7_qualified()
  n_all <- nrow(h); n_drop <- sum(!h$qualified)
  z <- h[h$qualified, , drop = FALSE]
  z$comp <- factor(sg_compartment_label(z$dataset),
                   levels = sg_compartments()$short)
  z$y0 <- NA_real_; z$y1 <- NA_real_
  for (d in unique(z$dataset)) {
    m <- s6_unit_positions(d); i <- z$dataset == d
    z$y0[i] <- m[z$baseline_unit[i]]; z$y1[i] <- m[z$effect_unit[i]]
  }
  z <- z[order(z$comp, z$y1, z$y0, z$GeneSymbol), , drop = FALSE]
  z$elsewhere <- z$baseline_unit != z$effect_unit
  z$lab_y <- z$y1
  for (d in unique(z$dataset)) {
    ii <- which(z$dataset == d)
    for (g in split(ii, z$y1[ii])) {
      k <- length(g)
      if (k > 1L) z$lab_y[g] <- z$y1[g] + (seq_len(k) - (k + 1) / 2) * 0.72
    }
  }
  lab <- do.call(rbind, lapply(unique(z$dataset), function(d) {
    u <- sg_units(); u <- u[u$dataset == d, , drop = FALSE]
    data.frame(comp = factor(sg_compartment_label(d),
                             levels = sg_compartments()$short),
               ypos = seq_len(nrow(u)), lab = u$display, stringsAsFactors = FALSE)
  }))

  # Brief section 22, corrected in the Part-25 fidelity pass. TWO different
  # spatial statements were being conflated, and ED7a and ED7b disagreed on the
  # same page as a result:
  #   "different dominant unit"          effect_unit != baseline_unit        15/15
  #   "outside baseline affinity"        the canonical top-2 affinity set    14/15
  # SNU13 (microglia, baseline CA1, effect CA2) is the single protein that
  # separates them: its effect unit is not its baseline peak, but it lies inside
  # that peak's top-2 affinity set, so ED7a counts it as inside. The canonical
  # effect_identity_relationship column is authoritative for affinity; the
  # cruder unit-inequality test is reported alongside it, never instead of it.
  # Neither statement is movement: both are two independent measurements of one
  # protein.
  z$outside_affinity <- z$effect_identity_relationship ==
    "effect_outside_baseline_affinity"
  n_unit <- sum(z$elsewhere)
  n_aff <- sum(z$outside_affinity)
  exc <- z$GeneSymbol[!z$outside_affinity]
  headline <- sprintf(paste0(
    "%d of %d robustness-qualified proteins show their strongest ",
    "phenotype-associated effect outside their dominant baseline spatial unit",
    "\n%d of %d are also outside the canonical baseline affinity set (a); ",
    "%s is the exception — its effect unit lies inside its baseline top-2 set"),
    n_unit, nrow(z), n_aff, nrow(z),
    if (length(exc)) paste(exc, collapse = ", ") else "none")

  p <- ggplot2::ggplot(z) +
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = 2, y = y0, yend = y1,
                                       colour = outside_affinity),
                          linewidth = nv_lw("reference_pt"), alpha = 0.85,
                          arrow = grid::arrow(length = grid::unit(0.9, "mm"),
                                              type = "closed", angle = 22)) +
    ggplot2::geom_point(ggplot2::aes(x = 1, y = y0), size = 0.85,
                        colour = "grey45") +
    ggplot2::geom_segment(ggplot2::aes(x = 2, xend = 2.06, y = y1, yend = lab_y),
                          colour = "grey72", linewidth = 0.12) +
    ggplot2::geom_text(ggplot2::aes(x = 2.09, y = lab_y, label = GeneSymbol,
                                    colour = outside_affinity),
                       hjust = 0, family = fam, size = nf_sz(5.0)) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = 0.92, y = ypos, label = lab),
                       hjust = 1, family = fam, size = nf_sz(5.0),
                       colour = "grey35") +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#C0442C",
                                            "FALSE" = "#6E87A0"), guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2),
      labels = c("most abundant\nHERE at baseline",
                 "strongest SUS\u2212RES\ndifference HERE"),
      limits = c(0.16, 2.95), expand = c(0, 0)) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(add = 0.9)) +
    ggplot2::facet_wrap(~comp, nrow = 1, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, title = headline,
                  caption = sprintf(paste0(
                    "%d of %d FDR-supported SUS-RES proteins: the %d that ",
                    "survived the CA2-SLM robustness audit plus the %d never ",
                    "in CA2-SLM. The %d the audit could not clear are ",
                    "excluded.\nEach arrow joins two independent ",
                    "measurements of one protein; nothing travels between the ",
                    "two units."),
                    nrow(z), n_all,
                    sum(z$qc_class == "robust_to_missingness_and_QC"),
                    sum(z$qc_class == "not_in_CA2_SLM_never_at_risk"), n_drop)) +
    nf_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT, lineheight = 1.05,
                                          face = "bold", colour = "grey15"),
      plot.title = ggplot2::element_text(
        size = nf_pt(5.6), face = "bold", colour = "grey10", hjust = 0,
        lineheight = 1.15,
        margin = ggplot2::margin(3.0, 0, 1.2, 0, "mm")),
      strip.text = ggplot2::element_text(size = nf_pt(5.4), face = "bold"),
      panel.spacing.x = ggplot2::unit(3.4, "mm"))
  out <- z[, c("GeneSymbol", "dataset", "baseline_unit", "effect_unit",
               "effect_identity_relationship", "qc_class", "elsewhere",
               "outside_affinity")]
  out$reading <- paste0(
    "two independent spatial statements about one protein. elsewhere = the ",
    "strongest effect is in a DIFFERENT unit from the baseline peak; ",
    "outside_affinity = it is also outside that peak's canonical top-2 ",
    "affinity set, which is the classification panel a counts. Nothing ",
    "moves: these are two separate measurements of the same protein.")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---- CA2-SLM locator ------------------------------------------------------
f9_ed_ca2_locator <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  accent <- unname(nv_dataset_colours()[["neuron_neuropil"]])
  arcdf <- function(cx, cy, r, t0, t1, n = 120) {
    t <- seq(t0, t1, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  t0 <- pi * 1.06; t1 <- pi * 2.32; span <- t1 - t0
  ca <- arcdf(0, 0, 1.00, t0, t1)
  slm <- arcdf(0, 0, 1.24, pi * 1.08, t0 + 0.66 * span)
  sr <- arcdf(0, 0, 1.12, pi * 1.08, pi * 2.30)
  so <- arcdf(0, 0, 0.88, pi * 1.08, pi * 2.30)
  seg <- function(d, frac, half = 0.11) {
    n <- nrow(d)
    d[max(1L, floor((frac - half) * n)):min(n, ceiling((frac + half) * n)), ,
      drop = FALSE]
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = sr, ggplot2::aes(x, y), colour = "grey82",
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = so, ggplot2::aes(x, y), colour = "grey82",
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = slm, ggplot2::aes(x, y), colour = "grey82",
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = ca, ggplot2::aes(x, y), colour = "grey82",
                       linewidth = 1.0) +
    ggplot2::geom_path(data = seg(slm, 0.62), ggplot2::aes(x, y),
                       colour = accent, linewidth = 1.8, lineend = "round") +
    ggplot2::annotate("text", x = 0, y = -1.28, label = "neuropil",
                      family = fam, size = nf_sz(5.0), colour = "grey35") +
    ggplot2::annotate("text", x = 0, y = -1.62, label = "CA2  SLM",
                      family = fam, size = nf_pt(6.4) * 0.3527777,
                      fontface = "bold", colour = accent) +
    ggplot2::coord_equal(xlim = c(-1.6, 1.6), ylim = c(-1.85, 1.45),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam)
  write_csv_safe(data.frame(compartment = "neuron_neuropil", region = "CA2",
                            layer = "slm", analysis_key = "CA2_slm",
                            note = "locator only; carries no quantity",
                            stringsAsFactors = FALSE), csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
