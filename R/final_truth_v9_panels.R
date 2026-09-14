# Editorial-v8 panel refinements.
#
# POLISH ONLY. No panel membership changes, no re-selection, no new inference.
# Every renderer here reuses the v7 data path and changes only presentation.

# ==========================================================================
# F2 a. Schematic: less prose, bigger anatomy
# ==========================================================================
#
# The v7 key spelled out each compartment in a sentence. The anatomy should do
# that work. The key is reduced to three two-line entries and the drawing is
# given the larger share of the panel.
f9_schematic <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  dcol <- nv_dataset_colours()
  np <- unname(dcol[["neuron_neuropil"]])
  so_col <- unname(dcol[["neuron_soma"]])
  mg_col <- unname(dcol[["microglia"]])

  arcdf <- function(cx, cy, r, t0, t1, n = 170) {
    t <- seq(t0, t1, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  t0 <- pi * 1.06; t1 <- pi * 2.32; span <- t1 - t0
  slm_end <- t0 + 0.66 * span          # SLM exists in CA1 and CA2 only
  gap <- data.frame(x = NA_real_, y = NA_real_)
  ca  <- arcdf(0, 0, 1.00, t0, t1)
  slm <- arcdf(0, 0, 1.26, pi * 1.08, slm_end)
  sr  <- arcdf(0, 0, 1.13, pi * 1.08, pi * 2.30)
  so  <- arcdf(0, 0, 0.87, pi * 1.08, pi * 2.30)
  dgg <- rbind(arcdf(0.44, -0.24, 0.26, pi * 1.22, pi * 1.98), gap,
               arcdf(0.44, -0.24, 0.26, pi * 0.30, pi * 0.86))
  dgm <- rbind(arcdf(0.44, -0.24, 0.39, pi * 1.22, pi * 1.98), gap,
               arcdf(0.44, -0.24, 0.39, pi * 0.30, pi * 0.86))
  dgp <- arcdf(0.44, -0.24, 0.12, pi * 1.30, pi * 1.90)
  at <- function(frac, r) {
    t <- t0 + frac * span
    c(x = r * cos(t), y = r * sin(t))
  }
  reg <- data.frame(region = c("CA1", "CA2", "CA3"),
                    f0 = c(0.00, 0.30, 0.64), f1 = c(0.28, 0.62, 1.00),
                    stringsAsFactors = FALSE)
  brk <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
    a <- arcdf(0, 0, 1.44, t0 + reg$f0[i] * span, t0 + reg$f1[i] * span, 40)
    a$region <- reg$region[i]; a
  }))
  rlab <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
    v <- at(mean(c(reg$f0[i], reg$f1[i])), 1.62)
    data.frame(x = v[["x"]], y = v[["y"]], l = reg$region[i],
               stringsAsFactors = FALSE)
  }))
  mg <- do.call(rbind, lapply(c(0.12, 0.46, 0.84), function(f) {
    v <- at(f, 0.70); data.frame(x = v[["x"]], y = v[["y"]])
  }))
  mg <- rbind(mg, data.frame(x = 0.44, y = -0.58))
  dgb <- arcdf(0.44, -0.24, 0.58, pi * 1.18, pi * 0.30, 40)

  # laminar ticks with short leaders, so the bands are named on the drawing
  lt <- t0 + 0.015 * span
  lyr <- data.frame(r = c(1.26, 1.13, 1.00, 0.87),
                    l = c("SLM", "SR", "SP", "SO"),
                    col = c(np, np, so_col, np), stringsAsFactors = FALSE)
  lyr$x <- lyr$r * cos(lt); lyr$y <- lyr$r * sin(lt)
  lyr$xe <- lyr$x - 0.30; lyr$ye <- lyr$y - seq(0.00, 0.33, length.out = 4)

  drawing <- ggplot2::ggplot() +
    ggplot2::geom_path(data = slm, ggplot2::aes(x, y), colour = np,
                       linewidth = 0.42, linetype = "22", lineend = "round") +
    ggplot2::geom_path(data = sr, ggplot2::aes(x, y), colour = np,
                       linewidth = 0.42, linetype = "22", lineend = "round") +
    ggplot2::geom_path(data = so, ggplot2::aes(x, y), colour = np,
                       linewidth = 0.42, linetype = "22", lineend = "round") +
    ggplot2::geom_path(data = dgm, ggplot2::aes(x, y), colour = np,
                       linewidth = 0.42, linetype = "22", lineend = "round") +
    ggplot2::geom_path(data = dgp, ggplot2::aes(x, y), colour = np,
                       linewidth = 0.42, linetype = "22", lineend = "round") +
    ggplot2::geom_path(data = ca, ggplot2::aes(x, y), colour = so_col,
                       linewidth = 1.6, lineend = "round") +
    ggplot2::geom_path(data = dgg, ggplot2::aes(x, y), colour = so_col,
                       linewidth = 1.3, lineend = "round") +
    ggplot2::geom_path(data = brk, ggplot2::aes(x, y, group = region),
                       colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_path(data = dgb, ggplot2::aes(x, y), colour = "grey55",
                       linewidth = 0.3) +
    ggplot2::geom_point(data = mg, ggplot2::aes(x, y), shape = 18, size = 1.2,
                        colour = mg_col) +
    ggplot2::geom_segment(data = lyr,
                          ggplot2::aes(x = x, y = y, xend = xe, yend = ye),
                          colour = "grey62", linewidth = 0.2) +
    ggplot2::geom_text(data = lyr,
                       ggplot2::aes(x = xe - 0.05, y = ye, label = l),
                       colour = lyr$col, family = fam, size = nf_sz(5.0),
                       hjust = 1) +
    ggplot2::geom_text(data = rlab, ggplot2::aes(x, y, label = l), family = fam,
                       size = nf_sz(5.8), fontface = "bold", colour = "grey12") +
    ggplot2::annotate("text", x = 0.44, y = -0.98, label = "DG", family = fam,
                      size = nf_sz(5.8), fontface = "bold", colour = "grey12") +
    ggplot2::annotate("text", x = 0, y = 1.92, label = "both hemispheres",
                      family = fam, size = nf_sz(5.0), colour = "grey35") +
    ggplot2::coord_equal(xlim = c(-2.25, 1.95), ylim = c(-1.30, 2.05),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam)

  # compact key: glyph, compartment, resolution - no sentences
  kd <- data.frame(row = 3:1,
                   style = c("neuropil", "soma", "microglia"),
                   name = c("Neuropil", "Neuronal soma", "Microglia-enriched ROI"),
                   res = c("region x layer", "region", "region"),
                   stringsAsFactors = FALSE)
  key <- ggplot2::ggplot(kd) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 0.085, y = row, yend = row,
                                       colour = style),
                          linewidth = ifelse(kd$style == "soma", 1.5, 0.5),
                          linetype = ifelse(kd$style == "neuropil", "22", "solid"),
                          lineend = "round") +
    ggplot2::geom_point(data = kd[kd$style == "microglia", ],
                        ggplot2::aes(x = 0.042, y = row), shape = 18,
                        size = 1.4, colour = mg_col) +
    ggplot2::geom_text(ggplot2::aes(0.12, row + 0.16, label = name), hjust = 0,
                       family = fam, size = nf_sz(5.2), fontface = "bold",
                       colour = "grey12") +
    ggplot2::geom_text(ggplot2::aes(0.12, row - 0.20, label = res), hjust = 0,
                       family = fam, size = nf_sz(5.0), colour = "grey40") +
    ggplot2::scale_colour_manual(
      values = c(neuropil = np, soma = so_col, microglia = mg_col),
      guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-0.02, 1.05), ylim = c(0.3, 3.8),
                             expand = FALSE) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 0, 0, "mm"))

  pl <- patchwork::wrap_plots(drawing, key, ncol = 2, widths = c(0.58, 0.42))
  out <- sg_order_for()
  out$note <- paste0(
    "declared vector schematic, no histology fabricated. DG MO and DG PO are ",
    "neuropil. SLM is measured in CA1 and CA2 only; there is no CA3 SLM unit. ",
    "Sampling and replicate counts are stated in the figure legend, not on the ",
    "artwork.")
  write_csv_safe(out, csv_path)
  nv_save_panel(pl, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# F2 g/h. Validation panels with the method prose removed
# ==========================================================================
#
# Brief section 11: the explanatory sentences move to the figure legend. The
# selection rule is unchanged - only the sentence describing it leaves the
# artwork. Contrast labels are shortened so the lollipops get the width.
f9_short_contrast <- function(x) {
  x <- gsub("_", " ", x)
  x <- sub(" vs mean other ", " vs rest ", x)
  x <- sub(" strata$", "", x)
  x <- sub(" regions$", "", x)
  x <- sub(" layers$", "", x)
  x
}

f9_external_main <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  k <- nv_read_csv(repo_path(panel$primary_source))
  k <- k[k$expected_match %in% TRUE, , drop = FALSE]
  k$level <- ifelse(grepl("strata", k$validation_domain),
                    "CA1 laminar identity", "Regional identity")
  k$level <- factor(k$level, levels = c("Regional identity", "CA1 laminar identity"))
  k$sig <- is.finite(k$p_adjust) & k$p_adjust < 0.05
  k$lab <- paste0(f9_short_contrast(k$internal_contrast), "  →  ",
                  k$external_signature)
  k <- k[order(k$level, k$NES), , drop = FALSE]
  k$ypos <- seq_len(nrow(k))

  p <- ggplot2::ggplot(k, ggplot2::aes(NES, ypos)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey65") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = ypos,
                                       colour = level),
                          linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(colour = level, shape = sig), size = 1.1,
                        fill = "white", stroke = 0.4) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21),
                                guide = "none") +
    ggplot2::scale_colour_manual(
      values = c("Regional identity" = "#1F3D52",
                 "CA1 laminar identity" = "#C2A878"), name = NULL) +
    ggplot2::scale_y_continuous(breaks = k$ypos, labels = k$lab,
                                expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::labs(x = "NES vs external hippocampal signature", y = NULL) +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "none")
  nreg <- sum(k$level == "Regional identity")
  if (nreg > 0 && nreg < nrow(k)) {
    p <- p + ggplot2::annotate("segment", x = -Inf, xend = Inf,
                               y = nreg + 0.5, yend = nreg + 0.5,
                               linewidth = 0.3, colour = "grey55")
  }
  k$legend_text_moved <- paste0(
    "expected pairings only; filled point = FDR < 0.05; the complete ",
    "specificity inventory is in ED2 - stated in the figure legend, not on the ",
    "artwork")
  write_csv_safe(k[, c("internal_contrast", "external_signature", "level",
                       "NES", "p_adjust", "legend_text_moved")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

f9_internal_main <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  g <- nv_read_csv(repo_path(panel$primary_source))
  g$level <- ifelse(grepl("CA1_strata", g$contrast),
                    "CA1 laminar identity", "Regional identity")
  g$level <- factor(g$level, levels = c("Regional identity", "CA1 laminar identity"))
  g <- g[order(g$contrast, -abs(g$NES), g$Description), , drop = FALSE]
  g <- g[!duplicated(g$contrast), , drop = FALSE]
  g <- g[order(g$level, g$NES), , drop = FALSE]
  g$ypos <- seq_len(nrow(g))
  g$lab <- paste0(f9_short_contrast(g$contrast), "  →  ", g$Description)
  g$sig <- is.finite(g$p_adjust) & g$p_adjust < 0.05

  p <- ggplot2::ggplot(g, ggplot2::aes(NES, ypos)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey65") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = ypos,
                                       colour = level),
                          linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(colour = level, shape = sig), size = 1.2,
                        fill = "white", stroke = 0.4) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21),
                                guide = "none") +
    ggplot2::scale_colour_manual(
      values = c("Regional identity" = "#1F3D52",
                 "CA1 laminar identity" = "#C2A878"), name = NULL) +
    ggplot2::scale_y_continuous(breaks = g$ypos, labels = g$lab,
                                expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::labs(x = "NES, internal anatomical GO program", y = NULL) +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "none")
  nreg <- sum(g$level == "Regional identity")
  if (nreg > 0 && nreg < nrow(g)) {
    p <- p + ggplot2::annotate("segment", x = -Inf, xend = Inf,
                               y = nreg + 0.5, yend = nreg + 0.5,
                               linewidth = 0.3, colour = "grey55")
  }
  g$legend_text_moved <- paste0(
    "strongest canonical term per contrast; filled point = FDR < 0.05; the ",
    "complete term inventory is in ED2 - stated in the figure legend, not on ",
    "the artwork. Set-size encoding was dropped: it competed with the effect ",
    "for the same visual channel and the sizes are in the source data.")
  write_csv_safe(g[, c("contrast", "Description", "level", "NES", "p_adjust",
                       "setSize", "legend_text_moved")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# F3 a. DAP track with manuscript-facing labels and visible totals
# ==========================================================================
#
# "canonical" and "claimable" are pipeline words. The panel now says
# FDR-supported and robustness-qualified, and carries the two totals inline so
# the reader gets 37 / 15 without counting cells.
f9_dap_track <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  z <- nf_dap_counts()
  b <- sg_blocks(z$unit, z$dataset)
  o <- b$order
  z$xpos <- match(paste(z$dataset, z$unit), paste(o$dataset, o$unit))
  n <- nrow(o)
  tot_c <- sum(z$canonical); tot_q <- sum(z$claimable)
  n_zero <- sum(z$canonical == 0L)

  long <- rbind(
    data.frame(xpos = z$xpos, row = 2L, n = z$canonical, stringsAsFactors = FALSE),
    data.frame(xpos = z$xpos, row = 1L, n = z$claimable, stringsAsFactors = FALSE))
  long$txt <- as.character(long$n)
  long$col <- ifelse(long$n == 0L, "grey76",
                     ifelse(long$row == 2L, "grey12", "#B23A28"))
  lab <- data.frame(
    y = c(2, 1),
    l = c(sprintf("FDR-supported DAPs   (%d total)", tot_c),
          sprintf("Robustness-qualified   (%d total)", tot_q)),
    stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(long, ggplot2::aes(xpos, row)) +
    ggplot2::geom_tile(fill = "grey97", colour = "white", linewidth = 0.35,
                       width = 1, height = 1) +
    ggplot2::geom_text(ggplot2::aes(label = txt), colour = long$col,
                       family = fam, size = nf_sz(5.4),
                       fontface = ifelse(long$n > 0, "bold", "plain")) +
    ggplot2::geom_text(data = lab,
                       ggplot2::aes(x = 0.5 - NF_LAB + 0.2, y = y, label = l),
                       inherit.aes = FALSE, hjust = 0, family = fam,
                       size = nf_sz(5.0), colour = "grey15") +
    ggplot2::scale_x_continuous(limits = c(0.5 - NF_LAB, n + 0.5),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = c(2, 1), labels = NULL,
                                limits = c(0.4, 2.6), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    nf_theme_tile() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(1, NF_RGT, 0, 1, "mm"))
  cend <- utils::head(b$compartment$end, -1)
  if (length(cend)) {
    p <- p + ggplot2::annotate("segment", x = cend + 0.5, xend = cend + 0.5,
                               y = 0.4, yend = 2.6, linewidth = 0.42,
                               colour = "grey25")
  }
  out <- z
  out$reading <- sprintf(paste0(
    "%d of 18 spatial units have no FDR-supported SUS-RES protein. %d ",
    "FDR-supported in total, %d robustness-qualified. CA2 SLM falls from %d ",
    "to %d, which is why an apparent protein-level concentration there is ",
    "weaker after qualification."),
    n_zero, tot_c, tot_q, z$canonical[z$unit == "CA2_slm"],
    z$claimable[z$unit == "CA2_slm"])
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# F3 d/e/f + g/h/i: three vertical biological columns
# ==========================================================================
#
# Brief section 15. Each column is one program: the GSEA curve on top carries
# the column header (compartment, location, program), and the protein panel
# underneath drops its own title so the pairing reads as one column rather than
# six unrelated plots.
#
# Brief section 17. g/h/i plot the SAME quantity (log2 fold change) over the
# SAME three contrasts, so they are directly comparable and were previously
# rescaled independently (limits 0.513 / 0.566 / 0.513). They now share one
# symmetric limit derived from the combined selected-protein values, and ONE
# colour bar: the legend is drawn on the last column only.

# the shared limit, derived once from all three programs
f9_prot_limit <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    v <- unlist(lapply(c("synaptic", "rna", "oxphos"), function(k) {
      z <- f9_prot_values(k)
      z$log2FC
    }))
    # round the combined maximum up to the next 0.05 so no value is squished
    cache <<- ceiling(max(abs(v), na.rm = TRUE) / 0.05) * 0.05
    cache
  }
})

# the per-program protein table, using the UNCHANGED v7 selection rule
f9_prot_values <- function(program_key, top_n = 7L) {
  pr <- s4_programs()
  prog <- pr[pr$key == program_key, , drop = FALSE]
  if (!nrow(prog)) stop("f9_prot_values: unknown program key: ", program_key,
                        call. = FALSE)
  ev <- s4_gsea_scores(prog)
  st <- stats::setNames(as.numeric(ev$ranked[ev$leading]), ev$leading)
  st <- st[order(-abs(st))]
  keep <- names(st)[seq_len(min(top_n, length(st)))]
  tok <- gsub("_", "", prog$unit[1])
  if (identical(prog$dataset[1], "microglia")) tok <- paste0(tok, "microglia")
  rows <- list()
  for (cb in list(c("res", "con"), c("sus", "con"), c("sus", "res"))) {
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
  z[!is.na(z$log2FC), , drop = FALSE]
}

f9_protein_zoom <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  key <- as.character(panel$program_key)
  show_legend <- isTRUE(as.logical(panel$show_legend %||% FALSE))
  z <- f9_prot_values(key)
  z$contrast <- factor(z$contrast,
                       levels = c("RES−CON", "SUS−CON", "SUS−RES"))
  ordg <- unique(z$gene[order(-abs(z$rank_statistic))])
  z$gene <- factor(z$gene, levels = rev(ordg))
  z$sig <- is.finite(z$BH_FDR) & z$BH_FDR < 0.05
  lim <- f9_prot_limit()

  # Part-26. The tile grid encoded a signed magnitude in colour for only 21
  # cells: the modal |log2FC| is 0.12-0.44 against a +/-0.95 ramp, so a 0.10
  # difference was 5% of the ramp and unreadable, and printing the value inside
  # each tile turned the panel into a table with a colour wash behind it.
  # Position is the accurate channel for a signed magnitude, so the same three
  # numbers per gene are now points on one shared log2FC axis. Contrast is
  # carried redundantly by colour AND shape, so the panel survives greyscale
  # printing and the common colourblindness forms.
  #
  # The FDR overlay is removed rather than restyled: every BH_FDR in all three
  # programs is >= 0.52, so geom_point(data = z[z$sig, ]) drew zero marks and
  # the sentence explaining it described something that was never on the page.
  # The absence of FDR support is stated once in the figure legend instead.
  COL <- c("RES−CON" = "#C0442C", "SUS−CON" = "#7C8A93",
           "SUS−RES" = "#2C6E9B")
  SHP <- c("RES−CON" = 16, "SUS−CON" = 17, "SUS−RES" = 15)

  p <- ggplot2::ggplot(z, ggplot2::aes(log2FC, gene)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey70") +
    ggplot2::geom_line(ggplot2::aes(group = gene), colour = "grey86",
                       linewidth = 0.26) +
    ggplot2::geom_point(ggplot2::aes(colour = contrast, shape = contrast),
                        size = 1.05) +
    ggplot2::scale_colour_manual(values = COL, name = "Contrast",
                                 guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_shape_manual(values = SHP, name = "Contrast",
                                guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_x_continuous(limits = c(-lim, lim), breaks = c(-0.8, 0, 0.8)) +
    ggplot2::labs(x = "log2 fold change", y = NULL) +
    nf_theme(grid = "y") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
      axis.title.x = ggplot2::element_text(size = nf_pt(5.4)),
      axis.text.y = ggplot2::element_text(size = NF_MIN_PT, face = "italic"),
      legend.position = if (show_legend) "inside" else "none",
      legend.position.inside = c(0.80, 0.74),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "mm"),
      legend.title = ggplot2::element_text(size = nf_pt(5.4)),
      legend.text = ggplot2::element_text(size = NF_MIN_PT),
      legend.key.size = ggplot2::unit(2.6, "mm"),
      plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"))
  z$shared_scale_note <- sprintf(paste0(
    "one symmetric log2 fold-change axis shared by all three protein panels, ",
    "limit +/-%.2f derived from the combined selected-protein values; contrast ",
    "is encoded by colour and shape together"), lim)
  z$fdr_note <- sprintf(paste0(
    "no protein in any of the three programs is FDR-supported (smallest BH ",
    "FDR = %.2f), so no significance marking is drawn"),
    min(z$BH_FDR, na.rm = TRUE))
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---- GSEA curve carrying the COLUMN HEADER -------------------------------
f9_gsea_curve <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  pr <- s4_programs()
  prog <- pr[pr$key == as.character(panel$program_key), , drop = FALSE]
  if (!nrow(prog)) stop("f9_gsea_curve: unknown program key", call. = FALSE)
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
  # ONE limit for every three-contrast NES strip in the family (F3 d/e/f and
  # ED6 c/d/e), so the same colour means the same NES everywhere
  lim <- f9_nes_strip_limit(th)
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

  # the COLUMN HEADER: this identifies the whole d/g, e/h or f/i column, so the
  # protein panel below carries no title of its own
  comp <- sg_compartment_label(prog$dataset[1])
  unit <- sg_unit_label(prog$unit[1], prog$dataset[1])
  idx <- as.character(panel$column_index %||% "")
  hdr <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = 0, xend = 1, y = 1.05, yend = 1.05,
                      colour = acc, linewidth = 0.5) +
    ggplot2::annotate("text", x = 0.012, y = 0.955, hjust = 0, vjust = 1,
                      family = fam, size = nf_sz(5.4), fontface = "bold",
                      colour = "grey10",
                      label = paste0(idx, "  ", comp, "  ", unit)) +
    ggplot2::annotate("text", x = 0.012, y = 0.63, hjust = 0, vjust = 1,
                      family = fam, size = nf_sz(5.2), colour = acc,
                      fontface = "bold", label = prog$label[1]) +
    ggplot2::annotate("text", x = 0.012, y = 0.28, hjust = 0, vjust = 1,
                      family = fam, size = nf_sz(5.0), colour = "grey35",
                      label = sprintf("NES %.2f   FDR %.0e", ev$NES, ev$FDR)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.08, 1.02),
                             expand = FALSE) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(plot.margin = ggplot2::margin(0.5, 1, 0, 1, "mm"))

  p <- patchwork::wrap_plots(hdr, top, mid, bot, ncol = 1,
                             heights = c(0.46, 1, 0.09, 0.24))
  out <- data.frame(
    program = prog$label[1], term_id = prog$term[1], dataset = prog$dataset[1],
    compartment_label = comp, spatial_unit_label = unit,
    NES = ev$NES, FDR = ev$FDR,
    RES_CON_NES = tr$NES[1], SUS_CON_NES = tr$NES[2], SUS_RES_NES = tr$NES[3],
    shared_NES_strip_limit = lim,
    column_role = paste0(
      "carries the shared column header for this program; the protein panel ",
      "directly below shows the leading-edge proteins of the same program and ",
      "repeats no context"),
    stringsAsFactors = FALSE)
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# F3 c. Anatomical bridge, strengthened
# ==========================================================================
#
# Brief section 14. The one-hippocampus solution is kept. What changes: the
# drawing gets more of the panel, the numbered anchors move OFF the arc onto
# leader lines so they no longer sit on top of the CA1/CA2/CA3 labels, and the
# same numerals 1/2/3 are reused as the column identifiers in d-i, so the panel
# says "these three programs occur HERE" rather than reading as one anatomy icon
# beside three unrelated mini-plots.
f9_anatomy_bridge <- function(panel, svg_path, csv_path, w_mm, h_mm) {
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

  arcdf <- function(cx, cy, r, t0, t1, n = 170) {
    t <- seq(t0, t1, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  t0 <- pi * 1.06; t1 <- pi * 2.32; span <- t1 - t0
  slm_end <- t0 + 0.66 * span
  gap <- data.frame(x = NA_real_, y = NA_real_)
  ca <- arcdf(0, 0, 1.00, t0, t1)
  slm <- arcdf(0, 0, 1.26, pi * 1.08, slm_end)
  sr <- arcdf(0, 0, 1.13, pi * 1.08, pi * 2.30)
  so <- arcdf(0, 0, 0.87, pi * 1.08, pi * 2.30)
  dgg <- rbind(arcdf(0.44, -0.24, 0.26, pi * 1.22, pi * 1.98), gap,
               arcdf(0.44, -0.24, 0.26, pi * 0.30, pi * 0.86))
  pale <- "grey82"
  at <- function(frac, r) {
    t <- t0 + frac * span
    c(x = r * cos(t), y = r * sin(t))
  }
  seg <- function(d, frac, half = 0.10) {
    n <- nrow(d)
    d[max(1L, floor((frac - half) * n)):min(n, ceiling((frac + half) * n)), ,
      drop = FALSE]
  }
  anch <- c(synaptic = 0.84, rna = 0.46, oxphos = 0.10)

  g <- ggplot2::ggplot() +
    ggplot2::geom_path(data = slm, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = sr, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = so, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = ca, ggplot2::aes(x, y), colour = pale,
                       linewidth = 1.3) +
    ggplot2::geom_path(data = dgg, ggplot2::aes(x, y), colour = pale,
                       linewidth = 1.1) +
    ggplot2::geom_path(data = seg(sr, anch[["synaptic"]]), ggplot2::aes(x, y),
                       colour = ex$accent[ex$key == "synaptic"],
                       linewidth = 1.7, lineend = "round") +
    ggplot2::geom_path(data = seg(ca, anch[["rna"]]), ggplot2::aes(x, y),
                       colour = ex$accent[ex$key == "rna"],
                       linewidth = 2.2, lineend = "round")
  v <- at(anch[["oxphos"]], 0.70)
  g <- g + ggplot2::geom_point(
    data = data.frame(x = v[["x"]], y = v[["y"]]), ggplot2::aes(x, y),
    shape = 18, size = 2.2, colour = ex$accent[ex$key == "oxphos"])

  # numbered anchors OUTSIDE the drawing on short leaders, so they cannot sit on
  # the region labels
  rad_on <- c(synaptic = 1.13, rna = 1.00, oxphos = 0.70)
  num <- do.call(rbind, lapply(seq_len(nrow(ex)), function(i) {
    k <- ex$key[i]
    a <- at(anch[[k]], rad_on[[k]])
    b <- at(anch[[k]], 1.74)
    data.frame(x0 = a[["x"]], y0 = a[["y"]], x = b[["x"]], y = b[["y"]],
               l = as.character(i), col = ex$accent[i], stringsAsFactors = FALSE)
  }))
  # region labels tucked just inside the arc, away from the anchor radius
  rl <- do.call(rbind, lapply(list(c(0.10, "CA1"), c(0.46, "CA2"),
                                   c(0.84, "CA3")), function(z) {
    w <- at(as.numeric(z[1]), 0.50)
    data.frame(x = w[["x"]], y = w[["y"]], l = z[2], stringsAsFactors = FALSE)
  }))
  g <- g +
    ggplot2::geom_segment(data = num,
                          ggplot2::aes(x = x0, y = y0, xend = x, yend = y),
                          colour = num$col, linewidth = 0.22) +
    ggplot2::geom_point(data = num, ggplot2::aes(x, y), size = 2.9,
                        shape = 21, fill = "white", colour = num$col,
                        stroke = 0.55) +
    ggplot2::geom_text(data = num, ggplot2::aes(x, y, label = l), family = fam,
                       size = nf_sz(5.2), fontface = "bold", colour = num$col) +
    ggplot2::geom_text(data = rl, ggplot2::aes(x, y, label = l), family = fam,
                       size = nf_sz(5.0), colour = "grey50") +
    ggplot2::coord_equal(xlim = c(-2.05, 2.05), ylim = c(-1.95, 2.00),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam)

  lim <- max(abs(rows$NES)) * 1.30
  callout <- function(i) {
    z <- rows[rows$key == ex$key[i], ]
    z$contrast <- factor(z$contrast,
                         levels = c("RES - CON", "SUS - CON", "SUS - RES"))
    hdr <- ggplot2::ggplot() +
      ggplot2::annotate("point", x = 0.022, y = 0.90, size = 2.9, shape = 21,
                        fill = "white", colour = ex$accent[i], stroke = 0.55) +
      ggplot2::annotate("text", x = 0.022, y = 0.90, label = as.character(i),
                        family = fam, size = nf_sz(5.2), fontface = "bold",
                        colour = ex$accent[i]) +
      ggplot2::annotate("text", x = 0.075, y = 1.0, hjust = 0, vjust = 1,
                        family = fam, size = nf_sz(5.2), fontface = "bold",
                        colour = "grey10",
                        label = paste0(ex$comp[i], "  ", ex$unit_label[i])) +
      ggplot2::annotate("text", x = 0.075, y = 0.60, hjust = 0, vjust = 1,
                        family = fam, size = nf_sz(5.0), colour = ex$accent[i],
                        label = ex$program[i]) +
      ggplot2::annotate("text", x = 0.075, y = 0.22, hjust = 0, vjust = 1,
                        family = fam, size = nf_sz(5.0), colour = "grey40",
                        label = ex$shape[i]) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.18, 1.06),
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
      ggplot2::theme(axis.text = ggplot2::element_text(size = NF_MIN_PT),
                     plot.margin = ggplot2::margin(0, 1, 0, 1, "mm"))
    patchwork::wrap_plots(hdr, strip, ncol = 1, heights = c(0.44, 0.56))
  }
  cos <- lapply(seq_len(nrow(ex)), callout)
  pl <- patchwork::wrap_plots(c(list(g), cos), nrow = 1,
                              widths = c(0.265, 0.245, 0.245, 0.245))
  out <- merge(rows, ex[, c("key", "shape")], by = "key")
  out$anchor_note <- paste0(
    "the numerals 1-3 identify the same three programs again as the column ",
    "headers of d-i; filled point = FDR < 0.05, open = not FDR-supported")
  write_csv_safe(out, csv_path)
  nv_save_panel(pl, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# The GSEA theme atlas on ONE common NES colour scale
# ==========================================================================
#
# Brief section 20. F3b, ED6a and ED6b are the same plot of the same quantity
# (median NES per theme x spatial unit) for the three contrasts of one
# three-group trajectory. Each previously rescaled to its own maximum
# (+/-2.207 SUS-RES, +/-1.757 RES-CON, +/-2.479 SUS-CON), so the same colour
# meant a different NES in each panel and the trajectory could not be read
# across them.
#
# Audit of the trade-off: under the common limit the weakest atlas (RES-CON)
# still spans 70.9% of the colour range, and the interquartile |NES| of the
# three atlases is nearly identical (0.378 / 0.386 / 0.451), so no atlas is
# compressed into a flat field. The common scale is adopted. The underlying NES
# is unchanged - only the mapping from NES to colour.

# the limit, computed once across ALL cells displayed in any of the three
f9_atlas_limit <- local({
  cache <- NULL
  function(th) {
    if (!is.null(cache)) return(cache)
    v <- unlist(lapply(c("RES - CON", "SUS - CON", "SUS - RES"), function(ct) {
      f9_atlas_cells(th, ct)$median_NES
    }))
    cache <<- max(abs(v), na.rm = TRUE)
    cache
  }
})

# the cell aggregation, identical to the frozen v7 rule
f9_atlas_cells <- function(th, contrast) {
  z <- th[th$contrast == contrast &
            th$theme_claim_eligible %in% TRUE &
            nzchar(as.character(th$theme_id)), , drop = FALSE]
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
  cells
}

f9_gsea_atlas <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  th <- nv_read_csv(repo_path(panel$primary_source))
  contrast <- as.character(panel$contrast %||% "SUS - RES")
  cells <- f9_atlas_cells(th, contrast)
  cells$sg_unit <- sg_resolve_unit(cells$spatial_unit, cells$dataset)
  b <- sg_blocks(cells$sg_unit, cells$dataset)
  o <- b$order
  cells$xpos <- match(paste(cells$dataset, cells$sg_unit),
                      paste(o$dataset, o$unit))
  SHORT <- c(synaptic_signaling_vesicle = "Synaptic signalling",
             rna_processing_splicing_rnp = "RNA processing",
             ribosome_translation = "Translation",
             mitochondrial_respiration_oxphos = "Mitochondrial respiration",
             autophagy_lysosome_endosome = "Autophagy",
             chromatin_organization = "Chromatin")
  ord <- names(SHORT)[names(SHORT) %in% cells$theme_id]
  cells <- cells[cells$theme_id %in% ord, , drop = FALSE]
  cells$ypos <- match(cells$theme_id, rev(ord))
  n <- nrow(o); ny <- length(ord)
  y_comp <- ny + 2.3; y_reg <- ny + 1.05
  lab_w <- as.numeric(panel$label_units %||% NF_LAB)
  shared <- f9_atlas_limit(th)
  lim <- shared * c(-1, 1)

  p <- ggplot2::ggplot(cells, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = median_NES), colour = "white",
                       linewidth = 0.1) +
    ggplot2::geom_point(data = cells[cells$n_fdr > 0, , drop = FALSE],
                        ggplot2::aes(xpos, ypos), size = 0.45,
                        colour = "grey10") +
    nv_diverging(limits = lim, name = "Median normalised\nenrichment score",
                 breaks = c(-2, 0, 2)) +
    ggplot2::geom_text(
      data = data.frame(y = seq_len(ny), l = unname(SHORT[rev(ord)])),
      ggplot2::aes(x = 0.5 - lab_w + 0.2, y = y, label = l),
      inherit.aes = FALSE, hjust = 0, family = fam, size = nf_sz(5.0),
      colour = "grey20") +
    ggplot2::scale_x_continuous(breaks = seq_len(n), labels = sg_axis_labels(b),
                                limits = c(0.5 - lab_w, n + 0.5),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny), labels = NULL,
                                limits = c(0.5, y_comp + 0.9),
                                expand = c(0, 0)) +
    ggplot2::annotate("segment", x = b$compartment$start - 0.5,
                      xend = b$compartment$end + 0.5,
                      y = y_comp - 0.32, yend = y_comp - 0.32,
                      linewidth = 0.45, colour = "grey25") +
    ggplot2::annotate("text", x = b$compartment$mid, y = y_comp,
                      label = b$compartment$label, family = fam,
                      size = nf_sz(5.4), fontface = "bold", colour = "grey12") +
    ggplot2::annotate("text", x = b$region$mid, y = y_reg,
                      label = b$region$label, family = fam, size = nf_sz(5.0),
                      colour = "grey30") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    contrast, ". Dot = at least one constituent GO term is ",
                    "FDR-supported. Colour scale is shared by all three ",
                    "contrast atlases (", sprintf("\u00b1%.2f", shared),
                    "), so they can be compared directly.")) +
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
  cells$shared_NES_scale_limit <- shared
  cells$shared_scale_note <- paste0(
    "one symmetric NES colour scale shared by RES-CON, SUS-CON and SUS-RES, ",
    "derived from the maximum absolute median NES across all cells displayed ",
    "in any of the three atlases")
  write_csv_safe(cells, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
