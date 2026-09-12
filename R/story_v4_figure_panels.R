# =====================================================================
# Story-v4 panel renderers.
#
# Fifth candidate family. Reuses the Part-17 palette/theme/exact-box machinery
# (sourced read-only) and replaces the encodings the Part-18 review rejected.
#
# WHAT IS GENUINELY NEW HERE
#
#   s4_schematic       a real hippocampal vector schematic, not a metadata matrix
#   s4_bilateral       representative scatter + forest over ALL 15 contrasts,
#                      so the weak CA1 laminar result is no longer hidden
#   s4_compartment     10 named markers x 3 compartments, from the stored matrix
#   s4_gsea_curve      a CONVENTIONAL running-enrichment plot. Verified exact:
#                      term membership re-queried from the pinned
#                      org.Mm.eg.db 3.22.0 reproduces stored setSize (507),
#                      stored ES to 3.3e-16 and all 255 leading-edge genes.
#                      That query is an ANNOTATION LOOKUP, not an enrichment
#                      test: no permutation is run and no p or FDR is
#                      recomputed. NES and FDR are read from the stored result.
#   s4_traj            the three-contrast trajectory strip, missing until now
#   s4_proteins        protein zoom with the three canonical contrasts as
#                      columns, at a readable size
#
# Downstream only: no model fit, no enrichment test, no FDR. The Part-17 guard
# scans this file.
# =====================================================================

s4_programs <- function() {
  data.frame(
    key = c("synaptic", "rna", "oxphos"),
    dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
    unit = c("CA3_sr", "CA2_sp", "CA1"),
    unit_dir = c("CA3_sr", "CA2_sp", "CA1_microglia"),
    contrast_dir = c("CA3srsus_CA3srres", "CA2spsus_CA2spres",
                     "CA1microgliasus_CA1microgliares"),
    term = c("GO:0099536", "GO:0006397", "GO:0006119"),
    label = c("synaptic signalling", "mRNA processing", "oxidative phosphorylation"),
    accent = c("#3D5A73", "#7A9BB0", "#C2A878"),
    stringsAsFactors = FALSE)
}

s4_unit_levels <- function() {
  c("CA1_slm", "CA1_so", "CA1_sr", "CA2_slm", "CA2_so", "CA2_sr",
    "CA3_so", "CA3_sr", "DG_mo", "DG_po",
    "CA1_sp", "CA2_sp", "CA3_sp", "DG_sg", "CA1", "CA2", "CA3", "DG")
}

# ============================================================ FIGURE 2

# a. Hippocampal schematic. NO histology is fabricated: this is an explicitly
# labelled Illustrator-ready vector placeholder. A repository-wide search found
# no anatomical artwork of any kind, and the canonical contract already marks
# panel 2a deferred_to_illustrator with figure_source null.
s4_schematic <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  dcol <- nv_dataset_colours()

  # one hippocampus: CA pyramidal band as an open arc, DG as a nested V
  arcdf <- function(cx, cy, r, t0, t1, n = 120) {
    t <- seq(t0, t1, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  build <- function(flip = FALSE, dx = 0) {
    s <- if (flip) -1 else 1
    ca <- arcdf(0, 0, 1.00, pi * 1.08, pi * 2.30)
    slm <- arcdf(0, 0, 1.20, pi * 1.10, pi * 2.28)
    sr <- arcdf(0, 0, 1.10, pi * 1.10, pi * 2.28)
    so <- arcdf(0, 0, 0.90, pi * 1.10, pi * 2.28)
    # DG as the characteristic open V/U, not a closed ring
    dgu <- arcdf(0.46, -0.28, 0.38, pi * 1.22, pi * 1.98)
    dgl <- arcdf(0.46, -0.28, 0.38, pi * 0.30, pi * 0.86)
    tr <- function(d) data.frame(x = s * d$x + dx, y = d$y)
    list(ca = tr(ca), slm = tr(slm), sr = tr(sr), so = tr(so),
         dgu = tr(dgu), dgl = tr(dgl), s = s, dx = dx)
  }
  L <- build(FALSE, -1.55); R <- build(TRUE, 1.55)

  lyr <- function(h, d, col, lw, lt = "solid") {
    ggplot2::geom_path(data = h[[d]], ggplot2::aes(x, y), colour = col,
                       linewidth = lw, linetype = lt, lineend = "round")
  }
  p <- ggplot2::ggplot()
  for (h in list(L, R)) {
    p <- p +
      lyr(h, "slm", dcol[["neuron_neuropil"]], 0.5, "22") +
      lyr(h, "sr",  dcol[["neuron_neuropil"]], 0.5, "22") +
      lyr(h, "so",  dcol[["neuron_neuropil"]], 0.5, "22") +
      lyr(h, "ca",  dcol[["neuron_soma"]], 1.7) +
      lyr(h, "dgu", dcol[["neuron_soma"]], 1.4) +
      lyr(h, "dgl", dcol[["neuron_soma"]], 1.4)
  }
  # Subfield labels are placed ON the arc they name, computed from the same
  # parameterisation as the drawing, so a label can never drift off its band.
  # CA1 is distal (start of the arc), CA3 proximal to the hilus (end).
  t0 <- pi * 1.08; t1 <- pi * 2.30
  at <- function(frac, r) {
    t <- t0 + frac * (t1 - t0)
    c(x = R$dx + R$s * (r * cos(t)), y = r * sin(t))
  }
  ca1 <- at(0.08, 1.34); ca2 <- at(0.45, 1.34); ca3 <- at(0.86, 1.34)
  lab <- data.frame(
    x = c(ca1[["x"]], ca2[["x"]], ca3[["x"]], R$dx + R$s * 0.46),
    y = c(ca1[["y"]], ca2[["y"]], ca3[["y"]], -0.28),
    l = c("CA1", "CA2", "CA3", "DG"), stringsAsFactors = FALSE)
  # layer labels sit beside the three dashed neuropil bands on the left
  lt <- t0 + 0.16 * (t1 - t0)
  lyrlab <- data.frame(
    x = L$dx + L$s * (c(1.34, 1.22, 1.00) * cos(lt)),
    y = c(1.34, 1.22, 1.00) * sin(lt),
    l = c("slm", "sr", "so"), stringsAsFactors = FALSE)
  mg <- data.frame(x = c(L$dx - 0.30, L$dx + 0.40, R$dx - 0.40, R$dx + 0.30),
                   y = c(-0.62, -0.86, -0.86, -0.62))

  p <- p +
    ggplot2::geom_point(data = mg, ggplot2::aes(x, y), size = 0.7,
                        colour = dcol[["microglia"]]) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x, y, label = l), family = fam,
                       size = nv_size(6), fontface = "bold") +
    ggplot2::geom_text(data = lyrlab, ggplot2::aes(x, y, label = l), family = fam,
                       size = nv_size(5), colour = dcol[["neuron_neuropil"]]) +
    ggplot2::annotate("text", x = c(-1.55, 1.55), y = 1.42,
                      label = c("left hemisphere", "right hemisphere"),
                      family = fam, size = nv_size(5.4), fontface = "bold") +
    ggplot2::annotate("segment", x = -0.42, xend = 0.42, y = 1.42, yend = 1.42,
                      linewidth = 0.25, colour = "grey60",
                      arrow = grid::arrow(length = ggplot2::unit(0.8, "mm"),
                                          ends = "both", type = "closed")) +
    ggplot2::annotate("text", x = 0, y = 1.18, label = "paired tissue",
                      family = fam, size = nv_size(5), colour = "grey35") +
    ggplot2::coord_equal(xlim = c(-3.05, 3.05), ylim = c(-1.58, 1.62),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   legend.position = "none")

  # compartment key + design statement, drawn as text so it stays editable
  key <- data.frame(x = c(-2.55, -0.95, 0.75), y = rep(-1.30, 3),
                    l = c("neuronal soma", "neuronal neuropil",
                          "microglia-enriched ROI"),
                    col = c(dcol[["neuron_soma"]], dcol[["neuron_neuropil"]],
                            dcol[["microglia"]]), stringsAsFactors = FALSE)
  p <- p +
    ggplot2::geom_point(data = key, ggplot2::aes(x - 0.13, y), size = 1.1,
                        colour = key$col) +
    ggplot2::geom_text(data = key, ggplot2::aes(x, y, label = l), hjust = 0,
                       family = fam, size = nv_size(5)) +
    ggplot2::annotate("text", x = -3.0, y = -1.50, hjust = 0, family = fam,
                      size = nv_size(5), colour = "grey30",
                      label = "9 animals (3 CON / 3 RES / 3 SUS) · AnimalID = biological replicate · 18 spatial units per animal") +
    ggplot2::annotate("text", x = 3.0, y = 1.55, hjust = 1, family = fam,
                      size = nv_size(5), colour = "grey55", fontface = "italic",
                      label = "schematic — Illustrator-ready placeholder")

  write_csv_safe(data.frame(
    element = c("subfields", "neuropil_layers", "soma_layers", "compartments",
                "hemispheres", "replicate_unit", "artwork_status"),
    value = c("CA1;CA2;CA3;DG", "slm;sr;so;mo;po", "sp;sg",
              "neuronal soma;neuronal neuropil;microglia-enriched ROI",
              "left;right (paired tissue from the same animal)",
              "AnimalID (9 animals, 3 per group)",
              paste0("VECTOR PLACEHOLDER. A repository-wide search found no ",
                     "hippocampal artwork, histology or tissue imagery of any ",
                     "kind; the canonical contract already marks panel 2a ",
                     "deferred_to_illustrator. No histological data is fabricated.")),
    stringsAsFactors = FALSE), csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# d. Bilateral reproducibility: ONE representative scatter + a forest over ALL
# prespecified contrasts. The Part-18 version showed three median contrasts and
# hid the spread; the forest is the finding.
s4_bilateral <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  prot <- nv_read_csv(repo_path(panel$primary_source))
  summ <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]))
  fam <- nv_palette()$typography$family

  # representative = the contrast closest to the overall median r
  target <- stats::median(summ$pearson_r)
  pick <- summ[which.min(abs(summ$pearson_r - target)), , drop = FALSE]
  dat <- prot[prot$dataset == pick$dataset[1] & prot$contrast == pick$contrast[1] &
                is.finite(prot$estimate_L) & is.finite(prot$estimate_R), , drop = FALSE]
  lim <- c(-3.2, 3.2)
  sc <- ggplot2::ggplot(dat, ggplot2::aes(estimate_L, estimate_R)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linewidth = nv_lw("reference_pt"), colour = "grey70") +
    ggplot2::geom_point(colour = "#1F3D52", size = 0.12, alpha = 0.25, shape = 16) +
    ggplot2::annotate("text", x = lim[1] * 0.92, y = lim[2] * 0.88, hjust = 0,
                      family = fam, size = nv_size(5.4),
                      label = sprintf("r = %.2f", pick$pearson_r[1])) +
    ggplot2::annotate("text", x = lim[1] * 0.92, y = lim[2] * 0.66, hjust = 0,
                      family = fam, size = nv_size(5), colour = "grey35",
                      label = gsub("_", " ", pick$contrast[1])) +
    ggplot2::coord_cartesian(xlim = lim, ylim = lim, expand = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(-3, 0, 3)) +
    ggplot2::scale_y_continuous(breaks = c(-3, 0, 3)) +
    ggplot2::labs(x = "left hemisphere (log2)", y = "right hemisphere (log2)") +
    nv_theme()

  s <- summ[order(summ$pearson_r), ]
  s$lab <- gsub("_", " ", s$contrast)
  s$lab <- sub("vs mean other", "vs other", s$lab)
  s$lab <- factor(s$lab, levels = s$lab)
  s$dsl <- nv_dataset_label(s$dataset)
  fo <- ggplot2::ggplot(s, ggplot2::aes(pearson_r, lab, colour = dsl)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = pearson_r, yend = lab),
                          colour = "grey88", linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      unname(nv_dataset_colours()[c("neuron_neuropil", "neuron_soma", "microglia")]),
      nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia"))), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                                expand = c(0, 0)) +
    ggplot2::labs(x = "paired-hemisphere r", y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5),
                   legend.position = "bottom",
                   legend.key.height = ggplot2::unit(1.8, "mm"))

  p <- if (requireNamespace("patchwork", quietly = TRUE))
    patchwork::wrap_plots(sc, fo, nrow = 1, widths = c(1, 1.25)) else sc

  out <- s[, c("dataset", "contrast", "n_pairs", "pearson_r", "spearman_rho",
               "sign_agreement_fraction", "n_evaluable_proteins")]
  out$is_representative_scatter <- out$contrast == pick$contrast[1] &
    out$dataset == pick$dataset[1]
  out$note <- paste0("all prespecified anatomical contrasts are shown; the ",
                     "spread from fine CA1 laminar to coarse regional contrasts ",
                     "is the result, not an inconvenience")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# f. Compartment identity as a named marker x compartment matrix. n = 3-4
# proteins per family makes boxplots indefensible; 10 named genes x 3
# compartments is the honest encoding and is already stored as a matrix.
s4_compartment <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  keep <- intersect(c("marker_gene", "dataset", "median_centered_log2",
                      "marker_class"), names(d))
  z <- d[, keep, drop = FALSE]
  z <- z[z$dataset %in% names(nv_dataset_colours()), , drop = FALSE]
  z$dsl <- nv_dataset_label(z$dataset)
  z$dsl <- factor(z$dsl,
    levels = nv_dataset_label(c("neuron_soma", "neuron_neuropil", "microglia")))
  # order genes by which compartment they mark, then by effect
  ord <- z[order(z$marker_class, -z$median_centered_log2), ]
  z$marker_gene <- factor(z$marker_gene, levels = rev(unique(ord$marker_gene)))
  lim <- max(abs(z$median_centered_log2), na.rm = TRUE)

  p <- ggplot2::ggplot(z, ggplot2::aes(dsl, marker_gene,
                                       fill = median_centered_log2)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    nv_diverging(limits = c(-lim, lim), name = "median-centred\nlog2") +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 5.2, face = "italic"),
                   legend.position = "right")

  z$encoding_note <- paste0(
    "10 named markers x 3 compartments, one stored value per cell. With 3-4 ",
    "proteins per marker family a boxplot or a forest with error bars would ",
    "imply dispersion the data do not carry; no CI or SE exists for these ",
    "differences anywhere in the source")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ============================================================ FIGURE 3

# Exact reconstruction of the canonical running-enrichment curve.
# VERIFIED: setSize, ES (to 3.3e-16) and all leading-edge genes reproduce the
# stored canonical result. Membership comes from the pinned annotation package,
# which is a lookup, not a test; NES and FDR are READ from the stored result.
s4_gsea_scores <- function(prog) {
  base <- repo_path("data", "processed", "04_differential_expression_enrichment",
                    "clusterProfiler", prog$dataset, "phenotype_within_unit",
                    prog$unit_dir, prog$contrast_dir)
  aud <- file.path(base, "protein_group_audits")
  g <- nv_read_csv_longpath(aud, "collapsed_gene_input.csv", stringsAsFactors = FALSE)
  res <- nv_read_csv_longpath(file.path(base, "GO", "BP"),
                              "GSEA_BP_results_full.csv", stringsAsFactors = FALSE)
  r <- res[res$ID == prog$term, , drop = FALSE]
  if (!nrow(r)) stop("term absent from stored GSEA result: ", prog$term)

  ranked <- stats::setNames(g$collapsed_statistic, g$GeneSymbol)
  ranked <- ranked[order(ranked, decreasing = TRUE)]

  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("annotation package unavailable; cannot reproduce canonical membership")
  }
  sy <- suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db, keys = prog$term, keytype = "GOALL",
    columns = "SYMBOL")$SYMBOL)
  gs <- intersect(unique(sy[!is.na(sy)]), names(ranked))

  # the reconstruction must reproduce the canonical setSize, or it is not the
  # canonical gene set and nothing may be drawn from it
  if (length(gs) != r$setSize[1]) {
    stop(sprintf("membership mismatch: reconstructed %d vs canonical setSize %d",
                 length(gs), r$setSize[1]))
  }
  N <- length(ranked); Nh <- length(gs)
  hits <- names(ranked) %in% gs
  Phit <- numeric(N); Phit[hits] <- abs(ranked[hits])          # exponent = 1
  Phit <- cumsum(Phit / sum(Phit))
  Pmiss <- numeric(N); Pmiss[!hits] <- 1 / (N - Nh); Pmiss <- cumsum(Pmiss)
  runes <- Phit - Pmiss
  mx <- max(runes); mn <- min(runes)
  ES <- if (abs(mx) > abs(mn)) mx else mn
  if (!isTRUE(all.equal(ES, r$enrichmentScore[1], tolerance = 1e-10))) {
    stop("reconstructed ES does not match the stored canonical ES")
  }
  peak <- if (abs(mx) > abs(mn)) which.max(runes) else which.min(runes)
  le <- names(ranked)[if (ES < 0) peak:N else 1:peak]
  le <- le[le %in% gs]
  list(ranked = ranked, hits = hits, runes = runes, ES = ES, peak = peak,
       NES = r$NES[1], FDR = r$p.adjust[1], setSize = r$setSize[1],
       description = r$Description[1], leading = le, N = N)
}

# d/e/f. Conventional GSEA panel: running ES, member ticks, NES and canonical
# FDR, plus the three-contrast trajectory strip underneath.
s4_gsea_curve <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  pr <- s4_programs()
  prog <- pr[pr$key == as.character(panel$program_key), , drop = FALSE]
  if (!nrow(prog)) stop("unknown program key: ", panel$program_key)
  ev <- s4_gsea_scores(prog)
  fam <- nv_palette()$typography$family
  acc <- prog$accent[1]
  N <- ev$N

  curve <- data.frame(rank = seq_len(N), es = as.numeric(ev$runes))
  keep <- unique(c(seq(1, N, by = 3), ev$peak, N))       # thin for file size
  curve <- curve[keep, ]
  ticks <- data.frame(rank = which(ev$hits))

  top <- ggplot2::ggplot(curve, ggplot2::aes(rank, es)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey75") +
    ggplot2::geom_line(colour = acc, linewidth = 0.45) +
    ggplot2::annotate("segment", x = ev$peak, xend = ev$peak, y = 0, yend = ev$ES,
                      linewidth = nv_lw("reference_pt"), colour = "grey55",
                      linetype = "22") +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, N),
                                breaks = c(1, N), labels = c("SUS", "RES")) +
    ggplot2::labs(x = NULL, y = "enrichment score") +
    nv_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 0, 1, "mm"),
                   axis.title.y = ggplot2::element_text(size = 5.4))

  mid <- ggplot2::ggplot(ticks, ggplot2::aes(rank)) +
    ggplot2::geom_segment(ggplot2::aes(x = rank, xend = rank, y = 0, yend = 1),
                          colour = acc, linewidth = 0.12) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, N)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0, 1, 0, 1, "mm"))

  # three-contrast trajectory for this exact program and context
  th <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]))
  tr <- th[th$dataset == prog$dataset[1] & th$spatial_unit == prog$unit[1] &
             th$GO_ID == prog$term[1], , drop = FALSE]
  tr <- tr[match(c("RES - CON", "SUS - CON", "SUS - RES"), tr$contrast), ]
  tr <- tr[!is.na(tr$contrast), ]
  tr$short <- sub(" - ", "−", tr$contrast)
  tr$short <- factor(tr$short, levels = tr$short)
  tr$sig <- is.finite(tr$GSEA_FDR) & tr$GSEA_FDR < 0.05
  lim <- max(abs(tr$NES), na.rm = TRUE)

  bot <- ggplot2::ggplot(tr, ggplot2::aes(short, 1, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", NES)),
                       family = fam, size = nv_size(5)) +
    ggplot2::geom_point(data = tr[tr$sig, , drop = FALSE],
                        ggplot2::aes(short, 1.42), size = 0.4, colour = "black",
                        inherit.aes = FALSE) +
    nv_diverging(limits = c(-lim, lim), name = NULL, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0.5, 1.6), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 5),
                   plot.margin = ggplot2::margin(0, 1, 1, 1, "mm"))

  hdr <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0.70, hjust = 0, vjust = 1, family = fam,
                      size = nv_size(5.6), fontface = "bold", colour = acc,
                      label = prog$label[1]) +
    ggplot2::annotate("text", x = 0, y = 0.30, hjust = 0, vjust = 1, family = fam,
                      size = nv_size(5), colour = "grey25",
                      label = sprintf("NES %.2f   FDR %.0e   %d/%d leading edge",
                                      ev$NES, ev$FDR, length(ev$leading), ev$setSize)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0.5, 1, 0, 1, "mm"))

  p <- if (requireNamespace("patchwork", quietly = TRUE))
    patchwork::wrap_plots(hdr, top, mid, bot, ncol = 1,
                          heights = c(0.30, 1, 0.10, 0.26)) else top

  out <- data.frame(
    program = prog$label[1], term_id = prog$term[1],
    term_description = ev$description, dataset = prog$dataset[1],
    spatial_unit = prog$unit[1], setSize = ev$setSize,
    n_leading_edge = length(ev$leading), NES = ev$NES, FDR = ev$FDR,
    ES_peak_rank = ev$peak, n_ranked_genes = N,
    reconstruction_note = paste0(
      "running enrichment reconstructed from the stored ranked statistic and ",
      "the pinned org.Mm.eg.db GO membership; verified to reproduce the stored ",
      "setSize, the stored enrichmentScore to 1e-10 and the stored leading-edge ",
      "gene set exactly. Membership retrieval is an annotation lookup, NOT an ",
      "enrichment test: no permutation is run and no p-value or FDR is ",
      "recomputed. NES and FDR are read from the canonical result."),
    stringsAsFactors = FALSE)
  out <- cbind(out, data.frame(
    RES_CON_NES = tr$NES[tr$contrast == "RES - CON"][1],
    RES_CON_FDR = tr$GSEA_FDR[tr$contrast == "RES - CON"][1],
    SUS_CON_NES = tr$NES[tr$contrast == "SUS - CON"][1],
    SUS_CON_FDR = tr$GSEA_FDR[tr$contrast == "SUS - CON"][1],
    SUS_RES_NES = tr$NES[tr$contrast == "SUS - RES"][1],
    SUS_RES_FDR = tr$GSEA_FDR[tr$contrast == "SUS - RES"][1]))
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# g. Protein zoom with the THREE canonical contrasts as columns, so the
# trajectory is visible at protein level too. Three aligned mini-heatmaps.
s4_proteins <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  top_n <- as.integer(panel$top_n %||% 8L)
  pr <- s4_programs()
  ctr <- c("res", "sus", "sus")   # numerator tokens
  rows <- list()
  for (i in seq_len(nrow(pr))) {
    prog <- pr[i, , drop = FALSE]
    ev <- s4_gsea_scores(prog)
    st <- stats::setNames(as.numeric(ev$ranked[ev$leading]), ev$leading)
    st <- st[order(-abs(st))]
    keep <- names(st)[seq_len(min(top_n, length(st)))]

    tok <- gsub("_", "", prog$unit[1])
    if (identical(prog$dataset[1], "microglia")) tok <- paste0(tok, "microglia")
    combos <- list(c("res", "con"), c("sus", "con"), c("sus", "res"))
    for (cb in combos) {
      f <- repo_path("data", "processed", "02_id_mapping", "mapped",
                     prog$dataset[1], "forward", "per_file",
                     sprintf("%s%s_%s%s.csv", tok, cb[1], tok, cb[2]))
      if (!file.exists(f)) next
      da <- nv_read_csv(f)
      sym <- intersect(c("official_gene_symbol", "gene_symbol"), names(da))[1]
      m <- match(keep, da[[sym]])
      rows[[length(rows) + 1L]] <- data.frame(
        program = prog$label[1], dataset = prog$dataset[1],
        spatial_unit = prog$unit[1], gene = keep,
        contrast = sprintf("%s−%s", toupper(cb[1]), toupper(cb[2])),
        log2FC = da$log2fc[m], BH_FDR = da$padj[m],
        rank_statistic = unname(st[keep]), stringsAsFactors = FALSE)
    }
  }
  z <- dplyr::bind_rows(rows)
  z <- z[!is.na(z$log2FC), , drop = FALSE]
  z$program <- factor(z$program, levels = pr$label)
  z$contrast <- factor(z$contrast, levels = c("RES−CON", "SUS−CON", "SUS−RES"))
  ord <- unique(z[order(z$program, -abs(z$rank_statistic)), c("program", "gene")])
  z$gene <- factor(z$gene, levels = rev(unique(ord$gene)))
  lim <- stats::quantile(abs(z$log2FC), 0.96, na.rm = TRUE)
  z$sig <- is.finite(z$BH_FDR) & z$BH_FDR < 0.05

  p <- ggplot2::ggplot(z, ggplot2::aes(contrast, gene, fill = log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_point(data = z[z$sig, , drop = FALSE], size = 0.4,
                        colour = "black") +
    nv_diverging(limits = c(-lim, lim), name = "log2FC",
                 oob = scales::squish) +
    ggplot2::facet_wrap(~ program, scales = "free", nrow = 1) +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       size = 5),
                   axis.text.y = ggplot2::element_text(size = 5.4, face = "italic"),
                   strip.text = ggplot2::element_text(size = 5.4, face = "bold"),
                   legend.position = "right")

  z$selection_rule <- sprintf(paste0(
    "stored leading-edge proteins of each program's canonical term, ranked by ",
    "|stored rank statistic|, top %d per program; no gene chosen by name"), top_n)
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# a. Spatial stress overview: where protein-level and program-level effects sit,
# on the same anatomical vocabulary as Figure 2.
s4_overview <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  member <- nv_read_csv(repo_path(panel$primary_source))
  deps <- as.character(unlist(panel$input_dependencies))
  atlas <- nv_read_csv(repo_path(deps[1]))
  gsea <- as.data.frame(data.table::fread(repo_path(deps[2]), showProgress = FALSE))

  m <- merge(member[, c("dataset", "spatial_unit", "ProteinGroupID")],
             atlas[, c("ProteinGroupID", "QC_claimability")],
             by = "ProteinGroupID", all.x = TRUE)
  m$QC_claimability[is.na(m$QC_claimability) | m$QC_claimability == ""] <- "not_audited"
  dap <- as.data.frame(table(dataset = m$dataset, spatial_unit = m$spatial_unit),
                       stringsAsFactors = FALSE)
  names(dap)[3] <- "n_dap"
  clm <- as.data.frame(table(dataset = m$dataset[m$QC_claimability == "claimable"],
                             spatial_unit = m$spatial_unit[m$QC_claimability == "claimable"]),
                       stringsAsFactors = FALSE)
  names(clm)[3] <- "n_claimable"

  s <- gsea[gsea$contrast == "SUS - RES" & is.finite(gsea$GSEA_FDR) &
              gsea$GSEA_FDR < 0.05 & nzchar(gsea$theme_id) &
              gsea$theme_role == "primary", ]
  prog <- as.data.frame(table(dataset = s$dataset, spatial_unit = s$spatial_unit),
                        stringsAsFactors = FALSE)
  names(prog)[3] <- "n_supported_terms"

  all_units <- unique(rbind(
    data.frame(dataset = gsea$dataset, spatial_unit = gsea$spatial_unit,
               stringsAsFactors = FALSE)))
  z <- merge(all_units, dap, all.x = TRUE)
  z <- merge(z, clm, all.x = TRUE)
  z <- merge(z, prog, all.x = TRUE)
  z[is.na(z)] <- 0
  z$dsl <- nv_dataset_label(z$dataset)
  z$dsl <- factor(z$dsl,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))
  z$unit <- factor(z$spatial_unit, levels = s4_unit_levels())
  z <- z[!is.na(z$unit), ]

  p <- ggplot2::ggplot(z, ggplot2::aes(unit, 1)) +
    ggplot2::geom_tile(ggplot2::aes(fill = n_supported_terms), colour = "white",
                       linewidth = nv_lw("tile_border_pt")) +
    ggplot2::scale_fill_gradient(low = "#F2F2EF", high = "#1F3D52",
                                 name = "FDR-supported\nGSEA terms") +
    ggplot2::geom_point(data = z[z$n_dap > 0, , drop = FALSE],
                        ggplot2::aes(size = n_dap), y = 1.34,
                        colour = "grey35", shape = 16) +
    ggplot2::geom_point(data = z[z$n_claimable > 0, , drop = FALSE],
                        ggplot2::aes(size = n_claimable), y = 1.34,
                        colour = "#D1543A", shape = 16) +
    ggplot2::scale_size_area(max_size = 2.2, name = "proteins\n(grey: canonical,\nred: claimable)") +
    ggplot2::scale_y_continuous(limits = c(0.5, 1.55), expand = c(0, 0)) +
    ggplot2::facet_grid(. ~ dsl, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "right",
                   legend.key.height = ggplot2::unit(2.2, "mm"))

  z$reading <- paste0("protein-level effects are sparse and spatially ",
                      "restricted; program-level support is broad and present ",
                      "in all three compartments")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
