# Nature-final v7 Figure-3 panels.
#
# DOWNSTREAM ONLY. Reshapes and labels canonical values. No model, no test, no
# p-value, no new inference. Hard 5 pt floor via nf_pt()/nf_sz().

# ==========================================================================
# a. Sparse single-protein differential abundance
# ==========================================================================
#
# The premise Figure 3 rests on: coordinated program-level effects occur DESPITE
# sparse individual-protein differential abundance. Two numeric rows aligned
# column-for-column with the atlas beneath.
#
# Numeric, not bars. A linear bar scale would let the CA2-SLM canonical count of
# 28 - which is the unit the robustness audit could largely not clear - dominate
# the panel and re-emphasise exactly the compromised burden.
#
# Both coordinates come from the canonical atlas: is_sus_res_fdr_supported marks
# the canonical hits, and CA2_SLM_robustness_class marks which survived. Nothing
# is recomputed.
nf_dap_counts <- function() {
  a <- nv_read_csv(repo_path("results", "tables", "11_spatial_systems", "atlas",
                             "protein_spatial_cell_affinity.csv"))
  h <- a[a$is_sus_res_fdr_supported %in% TRUE, , drop = FALSE]
  cls <- as.character(h$CA2_SLM_robustness_class)
  cls[is.na(cls) | !nzchar(cls)] <- "not_in_CA2_SLM_never_at_risk"
  h$claimable <- cls %in% c("robust_to_missingness_and_QC",
                            "not_in_CA2_SLM_never_at_risk")
  h$unit <- sg_resolve_unit(h$sus_res_strongest_spatial_unit_canonical, h$dataset)
  u <- sg_units()
  data.frame(
    dataset = u$dataset, unit = u$analysis_key, display = u$display,
    canonical = vapply(seq_len(nrow(u)), function(i)
      sum(h$dataset == u$dataset[i] & h$unit == u$analysis_key[i]), integer(1)),
    claimable = vapply(seq_len(nrow(u)), function(i)
      sum(h$dataset == u$dataset[i] & h$unit == u$analysis_key[i] & h$claimable),
      integer(1)),
    stringsAsFactors = FALSE)
}

# Panels a and b must behave as ONE block: identical column positions, the same
# left gutter for row labels and the same right gutter for the atlas legend.
# Both therefore use an x range that reserves NF_LAB units at the left, draw
# their row labels inside the panel at a fixed x, and carry the same right
# margin. Relying on ggplot to size two different axis-label sets identically
# would not align them.
NF_LAB <- 7.2   # x units reserved at the left for row labels, in a AND b
NF_RGT <- 7.6   # mm reserved at the right for the atlas legend, in a AND b

nf_dap_track <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  z <- nf_dap_counts()
  b <- sg_blocks(z$unit, z$dataset)
  o <- b$order
  z$xpos <- match(paste(z$dataset, z$unit), paste(o$dataset, o$unit))
  n <- nrow(o)

  long <- rbind(
    data.frame(xpos = z$xpos, row = 2L, n = z$canonical,
               lab = "canonical", stringsAsFactors = FALSE),
    data.frame(xpos = z$xpos, row = 1L, n = z$claimable,
               lab = "claimable", stringsAsFactors = FALSE))
  long$txt <- ifelse(long$n == 0L, "·", as.character(long$n))
  long$col <- ifelse(long$n == 0L, "grey70",
                     ifelse(long$row == 2L, "grey15", "#C0442C"))

  p <- ggplot2::ggplot(long, ggplot2::aes(xpos, row)) +
    ggplot2::geom_tile(fill = "grey97", colour = "white", linewidth = 0.35,
                       width = 1, height = 1) +
    ggplot2::geom_text(ggplot2::aes(label = txt), colour = long$col,
                       family = fam, size = nf_sz(5.2),
                       fontface = ifelse(long$n > 0, "bold", "plain")) +
    ggplot2::geom_text(
      data = data.frame(y = c(2, 1),
                        l = c("canonical DAP", "claimable DAP")),
      ggplot2::aes(x = 0.5 - NF_LAB + 0.2, y = y, label = l),
      inherit.aes = FALSE, hjust = 0, family = fam, size = nf_sz(5.0),
      colour = "grey20") +
    ggplot2::scale_x_continuous(limits = c(0.5 - NF_LAB, n + 0.5),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = c(2, 1), labels = NULL,
                                limits = c(0.5, 2.5), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    nf_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, NF_RGT, 0, 1, "mm"))
  cend <- utils::head(b$compartment$end, -1)
  if (length(cend)) {
    p <- p + ggplot2::annotate("segment", x = cend + 0.5, xend = cend + 0.5,
                               y = 0.5, yend = 2.5, linewidth = 0.42,
                               colour = "grey25")
  }
  out <- z
  out$reading <- sprintf(paste0(
    "FDR-supported SUS-RES proteins per spatial unit. %d of 18 units have none. ",
    "Total canonical %d, robustness-qualified %d. CA2 SLM falls from %d to %d."),
    sum(z$canonical == 0L), sum(z$canonical), sum(z$claimable),
    z$canonical[z$unit == "CA2_slm"], z$claimable[z$unit == "CA2_slm"])
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# b. GSEA spatial atlas
# ==========================================================================
nf_gsea_atlas <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  th <- nv_read_csv(repo_path(panel$primary_source))
  contrast <- as.character(panel$contrast %||% "SUS - RES")
  z <- th[th$contrast == contrast &
            th$theme_claim_eligible %in% TRUE &
            nzchar(as.character(th$theme_id)), , drop = FALSE]
  # seven GO terms are dual-assigned; deduplicate on the true key so they cannot
  # be counted twice
  z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID,
                           z$theme_id)), , drop = FALSE]
  key <- paste(z$dataset, z$spatial_unit, z$theme_id, sep = "\r")
  cells <- do.call(rbind, lapply(split(seq_len(nrow(z)), key), function(ix) {
    w <- z[ix, , drop = FALSE]
    data.frame(dataset = w$dataset[1], spatial_unit = w$spatial_unit[1],
               theme_id = w$theme_id[1], theme = w$manuscript_theme[1],
               n_fdr = sum(is.finite(w$GSEA_FDR) & w$GSEA_FDR < 0.05),
               median_NES = stats::median(w$NES, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(cells) <- NULL
  cells$sg_unit <- sg_resolve_unit(cells$spatial_unit, cells$dataset)
  b <- sg_blocks(cells$sg_unit, cells$dataset)
  o <- b$order
  cells$xpos <- match(paste(cells$dataset, cells$sg_unit),
                      paste(o$dataset, o$unit))
  # short display names so the row labels do not eat the plotting width
  SHORT <- c(synaptic_signaling_vesicle = "Synaptic signalling",
             rna_processing_splicing_rnp = "RNA processing",
             ribosome_translation = "Translation",
             mitochondrial_respiration_oxphos = "OXPHOS",
             autophagy_lysosome_endosome = "Autophagy",
             chromatin_organization = "Chromatin")
  ord <- names(SHORT)[names(SHORT) %in% cells$theme_id]
  cells <- cells[cells$theme_id %in% ord, , drop = FALSE]
  cells$ypos <- match(cells$theme_id, rev(ord))
  n <- nrow(o); ny <- length(ord)
  y_comp <- ny + 2.3; y_reg <- ny + 1.05
  lim <- max(abs(cells$median_NES), na.rm = TRUE) * c(-1, 1)

  p <- ggplot2::ggplot(cells, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = median_NES), colour = "white",
                       linewidth = 0.1) +
    ggplot2::geom_point(data = cells[cells$n_fdr > 0, , drop = FALSE],
                        ggplot2::aes(xpos, ypos), size = 0.45, colour = "grey10") +
    nv_diverging(limits = lim, name = "median\nNES") +
    ggplot2::geom_text(
      data = data.frame(y = seq_len(ny), l = unname(SHORT[rev(ord)])),
      ggplot2::aes(x = 0.5 - NF_LAB + 0.2, y = y, label = l),
      inherit.aes = FALSE, hjust = 0, family = fam, size = nf_sz(5.0),
      colour = "grey20") +
    ggplot2::scale_x_continuous(breaks = seq_len(n), labels = sg_axis_labels(b),
                                limits = c(0.5 - NF_LAB, n + 0.5),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny), labels = NULL,
                                limits = c(0.5, y_comp + 0.9), expand = c(0, 0)) +
    ggplot2::annotate("segment", x = b$compartment$start - 0.5,
                      xend = b$compartment$end + 0.5,
                      y = y_comp - 0.32, yend = y_comp - 0.32,
                      linewidth = 0.45, colour = "grey25") +
    ggplot2::annotate("text", x = b$compartment$mid, y = y_comp,
                      label = b$compartment$label, family = fam,
                      size = nf_sz(5.4), fontface = "bold", colour = "grey12") +
    ggplot2::annotate("text", x = b$region$mid, y = y_reg, label = b$region$label,
                      family = fam, size = nf_sz(5.0), colour = "grey30") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    "Contrast ", contrast,
                    ". Dot = at least one constituent GO term is FDR-supported. ",
                    "Claim-eligible themes only; qc_review themes excluded.")) +
    nf_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT, colour = "grey25"),
      axis.text.y = ggplot2::element_blank(),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.6, "mm"),
      legend.key.height = ggplot2::unit(3.4, "mm"))
  cend <- utils::head(b$compartment$end, -1)
  rend <- setdiff(utils::head(b$region$end, -1), cend)
  if (length(cend)) p <- p + ggplot2::annotate(
    "segment", x = cend + 0.5, xend = cend + 0.5, y = 0.5, yend = y_comp - 0.32,
    linewidth = 0.42, colour = "grey25")
  if (length(rend)) p <- p + ggplot2::annotate(
    "segment", x = rend + 0.5, xend = rend + 0.5, y = 0.5, yend = y_reg + 0.4,
    linewidth = 0.18, colour = "grey72")
  write_csv_safe(cells, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# c. ONE shared anatomical bridge with three anchored callouts
# ==========================================================================
#
# Replaces three repeated miniature hippocampi, which Part 22 measured at 6.2%
# ink over 29.6% of the page with region labels too small to read. One larger
# hippocampus carries all three anchors, and each callout keeps its own
# three-contrast strip.
#
# Wording is derived from the stored FDR, never asserted: an exemplar whose
# RES-CON arm is not FDR-supported is called susceptibility-associated, never
# divergent.
nf_bridge_rows <- function() {
  th <- nv_read_csv(repo_path("results", "tables", "10_biological_integration",
                              "gsea_wgcna_concordance", "global",
                              "ontology_aware_gsea_theme_assignments_all_contrasts.csv"))
  ex <- s6_exemplars()
  rows <- do.call(rbind, lapply(seq_len(nrow(ex)), function(i) {
    e <- ex[i, ]
    z <- th[th$dataset == e$dataset & th$spatial_unit == e$unit &
              th$GO_ID == e$go_id, , drop = FALSE]
    z <- z[!duplicated(z$contrast), , drop = FALSE]
    want <- c("RES - CON", "SUS - CON", "SUS - RES")
    if (!all(want %in% z$contrast)) {
      stop("bridge: exemplar ", e$key, " missing contrast(s)", call. = FALSE)
    }
    z <- z[match(want, z$contrast), , drop = FALSE]
    data.frame(key = e$key, dataset = e$dataset, unit = e$unit,
               program = e$program, contrast = want, NES = z$NES,
               FDR = z$GSEA_FDR, stringsAsFactors = FALSE)
  }))
  rows$supported <- rows$FDR < 0.05
  rows
}

nf_anatomy_bridge <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  dcol <- nv_dataset_colours()
  rows <- nf_bridge_rows()
  ex <- s6_exemplars()
  ex$accent <- unname(dcol[ex$dataset])
  ex$comp <- sg_compartment_label(ex$dataset)
  ex$unit_label <- sg_unit_label(ex$unit, ex$dataset)

  shape_of <- function(k) {
    z <- rows[rows$key == k, ]
    rc <- z[z$contrast == "RES - CON", ]
    sc <- z[z$contrast == "SUS - CON", ]
    if (rc$supported && sc$supported && sign(rc$NES) == sign(sc$NES)) {
      "graded stress-associated"
    } else if (!rc$supported && sc$supported) {
      "susceptibility-associated"
    } else {
      "not FDR-supported vs CON"
    }
  }
  ex$shape <- vapply(ex$key, shape_of, character(1))

  # ---- one hippocampus -------------------------------------------------
  arcdf <- function(cx, cy, r, t0, t1, n = 150) {
    t <- seq(t0, t1, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  t0 <- pi * 1.06; t1 <- pi * 2.32; span <- t1 - t0
  slm_end <- t0 + 0.66 * span
  gap <- data.frame(x = NA_real_, y = NA_real_)
  ca <- arcdf(0, 0, 1.00, t0, t1)
  slm <- arcdf(0, 0, 1.24, pi * 1.08, slm_end)
  sr <- arcdf(0, 0, 1.12, pi * 1.08, pi * 2.30)
  so <- arcdf(0, 0, 0.88, pi * 1.08, pi * 2.30)
  dgg <- rbind(arcdf(0.44, -0.24, 0.26, pi * 1.22, pi * 1.98), gap,
               arcdf(0.44, -0.24, 0.26, pi * 0.30, pi * 0.86))
  pale <- "grey80"
  at <- function(frac, r) {
    t <- t0 + frac * span
    c(x = r * cos(t), y = r * sin(t))
  }
  seg <- function(d, frac, half = 0.10) {
    n <- nrow(d)
    d[max(1L, floor((frac - half) * n)):min(n, ceiling((frac + half) * n)), ,
      drop = FALSE]
  }
  # CA1 ~ 0.10, CA2 ~ 0.46, CA3 ~ 0.84 along the arc
  anch <- c(synaptic = 0.84, rna = 0.46, oxphos = 0.10)

  g <- ggplot2::ggplot() +
    ggplot2::geom_path(data = slm, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = sr, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = so, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = ca, ggplot2::aes(x, y), colour = pale,
                       linewidth = 1.2) +
    ggplot2::geom_path(data = dgg, ggplot2::aes(x, y), colour = pale,
                       linewidth = 1.0)
  # highlight each exemplar on its own band
  g <- g +
    ggplot2::geom_path(data = seg(sr, anch[["synaptic"]]), ggplot2::aes(x, y),
                       colour = ex$accent[ex$key == "synaptic"],
                       linewidth = 1.5, lineend = "round") +
    ggplot2::geom_path(data = seg(ca, anch[["rna"]]), ggplot2::aes(x, y),
                       colour = ex$accent[ex$key == "rna"],
                       linewidth = 2.0, lineend = "round")
  v <- at(anch[["oxphos"]], 0.70)
  g <- g + ggplot2::geom_point(
    data = data.frame(x = v[["x"]], y = v[["y"]]), ggplot2::aes(x, y),
    shape = 18, size = 2.0, colour = ex$accent[ex$key == "oxphos"])
  # numbered anchors
  num <- do.call(rbind, lapply(seq_len(nrow(ex)), function(i) {
    w <- at(anch[[ex$key[i]]], 1.52)
    data.frame(x = w[["x"]], y = w[["y"]], l = as.character(i),
               col = ex$accent[i], stringsAsFactors = FALSE)
  }))
  rl <- do.call(rbind, lapply(list(c(0.10, "CA1"), c(0.46, "CA2"), c(0.84, "CA3")),
                              function(z) {
    w <- at(as.numeric(z[1]), 1.24)
    data.frame(x = w[["x"]], y = w[["y"]], l = z[2], stringsAsFactors = FALSE)
  }))
  g <- g +
    ggplot2::geom_text(data = rl, ggplot2::aes(x, y, label = l), family = fam,
                       size = nf_sz(5.0), colour = "grey45") +
    ggplot2::geom_point(data = num, ggplot2::aes(x, y), size = 2.6,
                        shape = 21, fill = "white", colour = num$col,
                        stroke = 0.5) +
    ggplot2::geom_text(data = num, ggplot2::aes(x, y, label = l), family = fam,
                       size = nf_sz(5.0), fontface = "bold", colour = num$col) +
    ggplot2::annotate("text", x = 0.44, y = -0.72, label = "DG", family = fam,
                      size = nf_sz(5.0), colour = "grey45") +
    ggplot2::coord_equal(xlim = c(-1.66, 1.66), ylim = c(-1.15, 1.72),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam)

  # ---- three callouts --------------------------------------------------
  lim <- max(abs(rows$NES)) * 1.30
  callout <- function(i) {
    z <- rows[rows$key == ex$key[i], ]
    z$contrast <- factor(z$contrast,
                         levels = c("RES - CON", "SUS - CON", "SUS - RES"))
    hdr <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 1.0, hjust = 0, vjust = 1,
                        family = fam, size = nf_sz(5.2), fontface = "bold",
                        colour = ex$accent[i],
                        label = sprintf("%d  %s  %s", i, ex$comp[i],
                                        ex$unit_label[i])) +
      ggplot2::annotate("text", x = 0, y = 0.55, hjust = 0, vjust = 1,
                        family = fam, size = nf_sz(5.0), colour = "grey15",
                        label = ex$program[i]) +
      ggplot2::annotate("text", x = 0, y = 0.14, hjust = 0, vjust = 1,
                        family = fam, size = nf_sz(5.0), colour = "grey40",
                        label = ex$shape[i]) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.25, 1.05),
                               expand = FALSE) +
      ggplot2::theme_void(base_family = fam)
    strip <- ggplot2::ggplot(z, ggplot2::aes(NES, contrast)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                          colour = "grey65") +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = contrast),
                            colour = ex$accent[i],
                            linewidth = nv_lw("reference_pt")) +
      ggplot2::geom_point(ggplot2::aes(shape = supported), size = 1.2,
                          fill = "white", colour = ex$accent[i], stroke = 0.4) +
      ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21),
                                  guide = "none") +
      ggplot2::scale_x_continuous(limits = c(-lim, lim), breaks = c(-2, 0, 2)) +
      ggplot2::scale_y_discrete(limits = rev(levels(z$contrast))) +
      ggplot2::labs(x = NULL, y = NULL) +
      nf_theme(grid = "x") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                     axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
                     plot.margin = ggplot2::margin(0, 1, 0, 1, "mm"))
    patchwork::wrap_plots(hdr, strip, ncol = 1, heights = c(0.42, 0.58))
  }
  cos <- lapply(seq_len(nrow(ex)), callout)
  pl <- patchwork::wrap_plots(
    c(list(g), cos), nrow = 1, widths = c(0.31, 0.23, 0.23, 0.23))
    patchwork::plot_annotation(
      caption = paste0("Filled point = FDR < 0.05; open = not FDR-supported. ",
                       "NES and FDR are read from the canonical enrichment ",
                       "table. The highlighted band marks where the sample was ",
                       "taken, not a quantitative anatomical map."),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(
        family = fam, size = NF_MIN_PT, colour = "grey35", hjust = 0)))

  out <- merge(rows, ex[, c("key", "shape")], by = "key")
  write_csv_safe(out, csv_path)
  nv_save_panel(pl, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# d/e/f. GSEA curve - MAIN version
# ==========================================================================
#
# Calls the SAME verified s4_gsea_scores(), which hard-stops unless it
# reproduces the stored setSize and enrichmentScore. What differs from ED6 is
# editorial, not numerical: the main curve carries compartment, region, layer,
# program, NES, FDR, the member-tick rug and a compact three-contrast strip, and
# drops the audit annotations (leading-edge counts, set sizes, reconstruction
# provenance) that belong in the Extended Data version.
nf_gsea_curve_main <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  pr <- s4_programs()
  prog <- pr[pr$key == as.character(panel$program_key), , drop = FALSE]
  if (!nrow(prog)) stop("nf_gsea_curve_main: unknown program key: ",
                        panel$program_key, call. = FALSE)
  ev <- s4_gsea_scores(prog)
  acc <- prog$accent[1]
  N <- ev$N

  curve <- data.frame(rank = seq_len(N), es = as.numeric(ev$runes))
  curve <- curve[unique(c(seq(1, N, by = 3), ev$peak, N)), ]
  ticks <- data.frame(rank = which(ev$hits))

  top <- ggplot2::ggplot(curve, ggplot2::aes(rank, es)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey78") +
    ggplot2::geom_line(colour = acc, linewidth = 0.45) +
    ggplot2::annotate("segment", x = ev$peak, xend = ev$peak, y = 0,
                      yend = ev$ES, linewidth = nv_lw("reference_pt"),
                      colour = "grey55", linetype = "22") +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, N),
                                breaks = c(1, N), labels = c("SUS", "RES")) +
    ggplot2::labs(x = NULL, y = "ES") +
    nf_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(0.5, 1, 0, 1, "mm"),
                   axis.title.y = ggplot2::element_text(size = NF_MIN_PT))

  mid <- ggplot2::ggplot(ticks, ggplot2::aes(rank)) +
    ggplot2::geom_segment(ggplot2::aes(x = rank, xend = rank, y = 0, yend = 1),
                          colour = acc, linewidth = 0.12) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, N)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    nf_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0, 1, 0, 1, "mm"))

  th <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]))
  tr <- th[th$dataset == prog$dataset[1] & th$spatial_unit == prog$unit[1] &
             th$GO_ID == prog$term[1], , drop = FALSE]
  tr <- tr[match(c("RES - CON", "SUS - CON", "SUS - RES"), tr$contrast), ]
  tr <- tr[!is.na(tr$contrast), ]
  tr$short <- c("R−C", "S−C", "S−R")[seq_len(nrow(tr))]
  tr$short <- factor(tr$short, levels = tr$short)
  tr$sig <- is.finite(tr$GSEA_FDR) & tr$GSEA_FDR < 0.05
  lim <- max(abs(tr$NES), na.rm = TRUE)

  bot <- ggplot2::ggplot(tr, ggplot2::aes(short, 1, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", NES)), family = fam,
                       size = nf_sz(5.0)) +
    ggplot2::geom_point(data = tr[tr$sig, , drop = FALSE],
                        ggplot2::aes(short, 1.44), size = 0.4, colour = "black",
                        inherit.aes = FALSE) +
    nv_diverging(limits = c(-lim, lim), name = NULL, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0.5, 1.62), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    nf_theme_tile() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
                   plot.margin = ggplot2::margin(0, 1, 0.5, 1, "mm"))

  comp <- sg_compartment_label(prog$dataset[1])
  unit <- sg_unit_label(prog$unit[1], prog$dataset[1])
  hdr <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 1.0, hjust = 0, vjust = 1, family = fam,
                      size = nf_sz(5.2), fontface = "bold", colour = "grey12",
                      label = paste0(comp, "  ", unit)) +
    ggplot2::annotate("text", x = 0, y = 0.52, hjust = 0, vjust = 1, family = fam,
                      size = nf_sz(5.0), colour = acc, label = prog$label[1]) +
    ggplot2::annotate("text", x = 0, y = 0.06, hjust = 0, vjust = 1, family = fam,
                      size = nf_sz(5.0), colour = "grey35",
                      label = sprintf("NES %.2f   FDR %.0e", ev$NES, ev$FDR)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.35, 1.05),
                             expand = FALSE) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(plot.margin = ggplot2::margin(0.5, 1, 0, 1, "mm"))

  p <- patchwork::wrap_plots(hdr, top, mid, bot, ncol = 1,
                             heights = c(0.36, 1, 0.09, 0.24))
  out <- data.frame(
    program = prog$label[1], term_id = prog$term[1], dataset = prog$dataset[1],
    compartment_label = comp, spatial_unit = prog$unit[1],
    spatial_unit_label = unit, NES = ev$NES, FDR = ev$FDR,
    RES_CON_NES = tr$NES[1], SUS_CON_NES = tr$NES[2], SUS_RES_NES = tr$NES[3],
    editorial_role = paste0(
      "MAIN: curve, member ticks, NES/FDR, spatial context and a compact ",
      "three-contrast strip. Audit detail (set size, leading-edge counts, ",
      "reconstruction provenance) is carried by the Extended Data version."),
    stringsAsFactors = FALSE)
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# g/h/i. Leading-edge protein zoom, ONE program per panel
# ==========================================================================
#
# The Part-22 audit found the single combined protein panel structurally broken
# at 35 mm: tiles collapsed to hairlines, contrast labels overprinted, titles
# clipped. Splitting into one panel per program gives each block ~56 mm, where
# the same renderer logic is fully legible.
#
# Selection is the established transparent rule, applied identically to each
# program: the stored leading-edge proteins of that program's canonical term,
# ranked by absolute stored rank statistic, top N. No gene is chosen by name.
nf_protein_zoom <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  top_n <- as.integer(panel$top_n %||% 7L)
  pr <- s4_programs()
  prog <- pr[pr$key == as.character(panel$program_key), , drop = FALSE]
  if (!nrow(prog)) stop("nf_protein_zoom: unknown program key: ",
                        panel$program_key, call. = FALSE)
  ev <- s4_gsea_scores(prog)
  st <- stats::setNames(as.numeric(ev$ranked[ev$leading]), ev$leading)
  st <- st[order(-abs(st))]
  keep <- names(st)[seq_len(min(top_n, length(st)))]

  tok <- gsub("_", "", prog$unit[1])
  if (identical(prog$dataset[1], "microglia")) tok <- paste0(tok, "microglia")
  combos <- list(c("res", "con"), c("sus", "con"), c("sus", "res"))
  rows <- list()
  for (cb in combos) {
    f <- repo_path("data", "processed", "02_id_mapping", "mapped",
                   prog$dataset[1], "forward", "per_file",
                   sprintf("%s%s_%s%s.csv", tok, cb[1], tok, cb[2]))
    if (!file.exists(f)) next
    da <- nv_read_csv(f)
    sym <- intersect(c("official_gene_symbol", "gene_symbol"), names(da))[1]
    m <- match(keep, da[[sym]])
    rows[[length(rows) + 1L]] <- data.frame(
      gene = keep,
      contrast = sprintf("%s−%s", toupper(cb[1]), toupper(cb[2])),
      log2FC = da$log2fc[m], BH_FDR = da$padj[m],
      rank_statistic = unname(st[keep]), stringsAsFactors = FALSE)
  }
  z <- do.call(rbind, rows)
  z <- z[!is.na(z$log2FC), , drop = FALSE]
  if (!nrow(z)) stop("nf_protein_zoom: no per-file DA rows for ",
                     panel$program_key, call. = FALSE)
  z$contrast <- factor(z$contrast,
                       levels = c("RES−CON", "SUS−CON", "SUS−RES"))
  ordg <- unique(z$gene[order(-abs(z$rank_statistic))])
  z$gene <- factor(z$gene, levels = rev(ordg))
  z$sig <- is.finite(z$BH_FDR) & z$BH_FDR < 0.05
  lim <- stats::quantile(abs(z$log2FC), 0.96, na.rm = TRUE)

  comp <- sg_compartment_label(prog$dataset[1])
  unit <- sg_unit_label(prog$unit[1], prog$dataset[1])
  p <- ggplot2::ggplot(z, ggplot2::aes(contrast, gene, fill = log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_point(data = z[z$sig, , drop = FALSE], size = 0.45,
                        colour = "black") +
    nv_diverging(limits = c(-lim, lim), name = "log2FC", oob = scales::squish) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0(comp, "  ", unit, "  -  ", prog$label[1])) +
    nf_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
      axis.text.y = ggplot2::element_text(size = NF_MIN_PT, face = "italic"),
      plot.title = ggplot2::element_text(size = NF_MIN_PT, face = "bold",
                                         colour = "grey12",
                                         margin = ggplot2::margin(b = 1)),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.4, "mm"),
      legend.key.height = ggplot2::unit(3, "mm"))
  z$selection_rule <- sprintf(paste0(
    "stored leading-edge proteins of the canonical term for this program, ",
    "ranked by |stored rank statistic|, top %d; identical rule for all three ",
    "programs; no gene chosen by name"), top_n)
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
