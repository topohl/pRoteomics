# Editorial-v8 fidelity pass.
#
# This is NOT a redesign. Every renderer here fixes a factual, quantitative or
# print-size defect found by the Part-25 encoding audit, or re-encodes data the
# panel already carries. No analysis is reopened, no new statistic is computed,
# and no panel changes membership.
#
# The defects fixed here, with the measurement that established each:
#
#  1. NES STRIP COLOUR LIMIT. e8_gsea_curve and s4_gsea_curve each computed
#     lim <- max(abs(NES)) INSIDE one panel, so F3 d/e/f and ED6 c/d/e carried
#     five different NES-to-colour mappings - on pages whose atlases state
#     "colour scale is shared by all three contrast atlases". Measured: -1.716
#     rendered at 99.8% of the ramp in one panel while the larger -1.775
#     rendered at 64.3% in another, so the weaker effect was drawn darker.
#     e8_nes_strip_limit() computes ONE limit across all three programs.
#
#  2. ED8c TICK OVERPRINT. With a 31.6 mm plot width and a linear 0-1.05 axis,
#     the "0" tick sits at 0.0 mm and "0.05" at 1.5 mm while their combined
#     half-widths are 2.2 mm, so the two labels overprinted. The attainable
#     floor sat at 0.34% of the axis and was indistinguishable from the axis
#     line - the panel's own stated reading was the one fact it could not show.
#
#  3. ED7a PERCENTAGE AXIS. A proportion was drawn on an axis running to 122%,
#     widened only to fit the k/n labels, and all six values lay in 88.9-100%
#     so ~89% of every bar was ink common to all six rows.
#
#  4. F2e SILENT CENSORING. The plotted column is capped at 3: Ptbp2 (5.100),
#     Npm1 (4.171) and Anp32a (3.105) were all drawn as the identical darkest
#     red, and the sidecar carried only the capped 3.0, so the true values were
#     unrecoverable from the released source data.
#
#  5. F2b RANK SAWTOOTH. The jitter offset was (rank %% 5 - 2) * 0.07 applied
#     after sorting by value, so horizontal position was a deterministic
#     function of rank. Extracting all 323 circle positions from the rendered
#     SVG confirmed exactly 5 discrete lanes and a period-5 staircase in both
#     tails - visible diagonal structure with no data behind it.

# ==========================================================================
# One NES colour limit for every three-contrast strip in the family
# ==========================================================================
#
# F3 d/e/f and ED6 c/d/e draw the same nine numbers: three programs x three
# contrasts. They must therefore share one mapping. The limit is the maximum
# absolute NES over all nine, computed once and cached.
e8_nes_strip_limit <- local({
  cache <- NULL
  function(th) {
    if (!is.null(cache)) return(cache)
    pr <- s4_programs()
    v <- unlist(lapply(seq_len(nrow(pr)), function(i) {
      z <- th[th$dataset == pr$dataset[i] & th$spatial_unit == pr$unit[i] &
                th$GO_ID == pr$term[i], , drop = FALSE]
      z$NES[z$contrast %in% c("RES - CON", "SUS - CON", "SUS - RES")]
    }))
    cache <<- max(abs(v), na.rm = TRUE)
    cache
  }
})

# ==========================================================================
# ED6 c/d/e. The detailed GSEA curve, on the shared strip limit
# ==========================================================================
#
# Derived from the frozen s4_gsea_curve. Two changes only: the NES strip uses
# e8_nes_strip_limit() instead of a per-panel maximum, and the theme is nf_theme
# so every label clears the 5 pt floor.
e8_ed_gsea_curve <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  pr <- s4_programs()
  prog <- pr[pr$key == as.character(panel$program_key), , drop = FALSE]
  if (!nrow(prog)) stop("e8_ed_gsea_curve: unknown program key", call. = FALSE)
  ev <- s4_gsea_scores(prog)
  acc <- prog$accent[1]
  N <- ev$N
  curve <- data.frame(rank = seq_len(N), es = as.numeric(ev$runes))
  curve <- curve[unique(c(seq(1, N, by = 3), ev$peak, N)), ]
  ticks <- data.frame(rank = which(ev$hits))

  top <- ggplot2::ggplot(curve, ggplot2::aes(rank, es)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey75") +
    ggplot2::geom_line(colour = acc, linewidth = 0.45) +
    ggplot2::annotate("segment", x = ev$peak, xend = ev$peak, y = 0,
                      yend = ev$ES, linewidth = nv_lw("reference_pt"),
                      colour = "grey55", linetype = "22") +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, N),
                                breaks = c(1, N), labels = c("SUS", "RES")) +
    ggplot2::labs(x = NULL, y = "enrichment score") +
    nf_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 1, 0, 1, "mm"),
                   axis.title.y = ggplot2::element_text(size = nf_pt(5.4)))
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
  lim <- e8_nes_strip_limit(th)

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

  hdr <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, family = fam,
                      size = nf_sz(5.4), fontface = "bold", colour = acc,
                      label = prog$label[1]) +
    ggplot2::annotate("text", x = 0, y = 0.25, hjust = 0, vjust = 1,
                      family = fam, size = nf_sz(5.0), colour = "grey35",
                      label = sprintf("NES %.2f  FDR %.0e  %d/%d leading edge",
                                      ev$NES, ev$FDR, length(ev$leading),
                                      sum(ev$hits))) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.2, 1.05),
                             expand = FALSE) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(plot.margin = ggplot2::margin(0.5, 1, 0, 1, "mm"))

  p <- patchwork::wrap_plots(hdr, top, mid, bot, ncol = 1,
                             heights = c(0.30, 1, 0.08, 0.22))
  out <- data.frame(
    program = prog$label[1], term_id = prog$term[1], dataset = prog$dataset[1],
    NES = ev$NES, FDR = ev$FDR,
    RES_CON_NES = tr$NES[1], SUS_CON_NES = tr$NES[2], SUS_RES_NES = tr$NES[3],
    shared_NES_strip_limit = lim,
    shared_scale_note = paste0(
      "the three-contrast NES strip uses ONE symmetric colour limit shared by ",
      "Figure 3 d/e/f and ED6 c/d/e, so the same colour means the same NES in ",
      "every strip in the family"),
    stringsAsFactors = FALSE)
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# ED8 c. The whole-network null, on a log axis
# ==========================================================================
#
# Defect 2. Chart type, marks and data are unchanged; only the scale changes.
# On the linear axis the attainable floor (0.0036) and the 0.05 rule were both
# pinned to the y axis and their tick labels overprinted. A log axis separates
# them and makes the panel's stated reading - that the enumeration had
# resolution far below 0.05 and still found nothing - directly readable.
e8_ed_nulls <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  n <- nv_read_csv(repo_path(panel$primary_source))
  z <- data.frame(dataset = n$dataset, p = n$exact_p,
                  floor = n$min_attainable_p, n_edges = n$n_edges,
                  stringsAsFactors = FALSE)
  z$lab <- nv_dataset_label(z$dataset)
  z$lab <- factor(z$lab, levels = nv_dataset_label(
    c("neuron_neuropil", "neuron_soma", "microglia")))
  lo <- 10^floor(log10(min(z$floor, na.rm = TRUE)))

  p <- ggplot2::ggplot(z, ggplot2::aes(p, lab)) +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "22",
                        linewidth = nv_lw("reference_pt"), colour = "#D1543A") +
    ggplot2::geom_segment(ggplot2::aes(x = floor, xend = p, yend = lab),
                          colour = "grey85", linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(x = floor), colour = "grey60", size = 1,
                        shape = 1) +
    ggplot2::geom_point(size = 1.4, colour = "#1F3D52") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("p = %.2f", p)),
                       hjust = -0.3, family = fam, size = nf_sz(5.0)) +
    ggplot2::scale_x_log10(limits = c(lo, 22),
                           breaks = c(0.001, 0.01, 0.05, 1),
                           labels = c("0.001", "0.01", "0.05", "1")) +
    ggplot2::labs(x = "exact whole-network p",
                  y = NULL,
                  caption = paste0(
                    "Open circle = smallest
attainable p. Log axis.")) +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
                   axis.text.y = ggplot2::element_text(size = NF_MIN_PT))
  z$reading <- paste0(
    "an informative null: the enumeration had resolution to ",
    sprintf("%.4f", min(z$floor, na.rm = TRUE)),
    " and found no whole-network group difference. The axis is log10 because ",
    "on a linear axis the attainable floor and the 0.05 rule both collapse ",
    "onto the y axis")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# ED7 a. Counts, not percentages
# ==========================================================================
#
# Defect 3. The quantity is k of n with n between 7 and 37, so the honest
# encoding is the count: a bar running 0 -> n_hits with the "outside" portion
# filled and the remainder left open. Bar length now carries n (which differs
# five-fold across rows and was previously invisible), and the discriminating
# information is the open tail rather than the last 11% of six near-identical
# bars. The two subsets that are the SAME 15 proteins under two definitions are
# labelled as such instead of reading as independent agreement.
e8_ed_identity_subsets <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  c1 <- nv_read_csv(repo_path(panel$primary_source))
  c1$lab <- gsub("_", " ", c1$subset)
  c1$primary <- c1$subset == "CA2_SLM_robustness_qualified"
  c1$pct <- 100 * c1$fraction_outside_baseline_affinity

  # rows that are the identical protein set under two subset definitions
  key <- paste(c1$n_hits, c1$effect_outside_baseline_affinity, c1$n_CA2_SLM)
  dup <- key %in% key[duplicated(key)]
  c1$lab[dup] <- paste0(c1$lab[dup], " †")
  c1$lab <- factor(c1$lab, levels = rev(c1$lab))

  p <- ggplot2::ggplot(c1, ggplot2::aes(y = lab)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = n_hits, yend = lab),
                          colour = "grey82", linewidth = 1.9,
                          lineend = "butt") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = effect_outside_baseline_affinity,
                                       yend = lab, colour = primary),
                          linewidth = 1.9, lineend = "butt") +
    ggplot2::geom_text(ggplot2::aes(x = n_hits,
                                    label = sprintf("%d/%d",
                                                    effect_outside_baseline_affinity,
                                                    n_hits)),
                       hjust = -0.22, family = fam, size = nf_sz(5.0),
                       colour = "grey25") +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#1F3D52",
                                            "FALSE" = "#7C8A93"),
                                 guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 42), expand = c(0, 0),
                                breaks = c(0, 10, 20, 30, 37)) +
    ggplot2::labs(
      x = "proteins whose strongest effect is outside the canonical baseline affinity set (of n)",
      y = NULL,
      caption = paste0(
        "Bar length is n, the size of the subset; the filled portion is the ",
        "count outside the baseline affinity set. † marks two rows that ",
        "are the SAME 15 proteins under two subset definitions.")) +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   axis.text.x = ggplot2::element_text(size = NF_MIN_PT))
  c1$note <- paste0(
    "the result survives every robustness restriction; bar length encodes the ",
    "subset size n, which ranges from 7 to 37 and was not visible when the ",
    "panel plotted percentages on a 0-122 axis")
  write_csv_safe(c1, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# F2 b. Sequencing depth, with an honest horizontal offset
# ==========================================================================
#
# Defect 5. The frozen renderer sorted by value and then set the offset to
# (rank %% 5 - 2) * 0.07, so horizontal position was a deterministic function
# of rank: five discrete lanes and a period-5 diagonal staircase in both tails
# that a reader can mistake for structure. The offset is now proportional to
# LOCAL DENSITY - points are binned on the value axis and spread symmetrically
# across a width proportional to the square root of the bin count - so
# horizontal extent means "many runs here" and nothing else. Quartile and
# median crossbars state the spread instead of leaving it to ink density.
# Fully deterministic: no random jitter is used anywhere.
e8_sina_offset <- function(v, nbin = 34L, wmax = 0.40) {
  b <- cut(v, breaks = nbin, labels = FALSE, include.lowest = TRUE)
  n <- table(b)
  off <- numeric(length(v))
  for (k in unique(b)) {
    i <- which(b == k)
    m <- length(i)
    if (m == 1L) { off[i] <- 0; next }
    w <- wmax * sqrt(m / max(n))
    off[i] <- seq(-w, w, length.out = m)
  }
  off
}

e8_depth_compact <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  d <- nv_read_csv(repo_path(panel$primary_source))
  cn <- intersect(c("Proteins.Identified", "n_proteins", "proteins_identified"),
                  names(d))[1]
  if (is.na(cn)) stop("e8_depth_compact: no protein-count column", call. = FALSE)
  d$val <- as.numeric(d[[cn]])
  ccol <- intersect(c("celltype_layer", "qc_compartment", "dataset"), names(d))[1]
  if (is.na(ccol)) stop("e8_depth_compact: no compartment column", call. = FALSE)
  map <- c(neuron_neuropil = "neuron_neuropil", neuropil = "neuron_neuropil",
           neuron_soma = "neuron_soma", soma = "neuron_soma",
           microglia = "microglia")
  d$dataset <- unname(map[as.character(d[[ccol]])])
  if ("exclude" %in% names(d)) d <- d[!(d$exclude %in% TRUE), , drop = FALSE]
  d <- d[!is.na(d$dataset) & is.finite(d$val), , drop = FALSE]
  d$lab <- factor(sg_compartment_label(d$dataset),
                  levels = sg_compartments()$short)
  d <- d[!is.na(d$lab), , drop = FALSE]
  d <- d[order(d$lab, d$val), , drop = FALSE]
  d$off <- unlist(lapply(split(d$val, d$lab), e8_sina_offset), use.names = FALSE)

  # Part-26. DESCRIPTIVE violin, no test. The frozen renderer drew 323 raw
  # points whose horizontal position was a function of rank; the points are
  # kept but the distribution shape is now stated by a density outline, which
  # is what the panel is actually claiming. A violin is appropriate here only
  # because every group clears the small-sample threshold at which smoothed
  # densities become meaningless (n = 180 / 71 / 72 acquisitions), and the
  # neuropil distribution is multimodal, which a box plot would have hidden.
  #
  # The two n are NOT the same number and are labelled separately: the violin
  # summarises ACQUISITIONS, while the biological replicate is the animal
  # (9 animals, 3 per group). Nothing here is a between-group comparison.
  acq <- table(d$lab)
  n_animals <- if ("AnimalID" %in% names(d))
    length(unique(d$AnimalID[!is.na(d$AnimalID)])) else NA_integer_

  q <- function(p) function(v) stats::quantile(v, p, names = FALSE)
  p <- ggplot2::ggplot(d, ggplot2::aes(as.integer(lab) + off, val)) +
    ggplot2::geom_violin(ggplot2::aes(x = as.integer(lab), group = lab),
                         width = 0.86, colour = "grey45", fill = NA,
                         linewidth = 0.22, bw = "nrd0", trim = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = lab), size = 0.38, alpha = 0.55) +
    ggplot2::stat_summary(ggplot2::aes(x = as.integer(lab)), fun.min = q(0.25),
                          fun.max = q(0.75), fun = stats::median,
                          geom = "crossbar", width = 0.34, linewidth = 0.16,
                          colour = "grey20", fill = NA) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        unname(nv_dataset_colours()[sg_compartment_levels()]),
        sg_compartments()$short), guide = "none") +
    ggplot2::scale_x_continuous(breaks = 1:3, labels = sg_compartments()$short,
                                limits = c(0.5, 3.5)) +
    ggplot2::labs(x = NULL, y = "proteins identified",
                  caption = sprintf(paste0(
                    "%s acquisitions\nfrom %d animals (3 per group), which ",
                    "are the biological replicates"),
                    paste(as.integer(acq[sg_compartments()$short]),
                          collapse = " / "), n_animals)) +
    nf_theme(grid = "y") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT, angle = 30,
                                          hjust = 1),
      plot.caption = ggplot2::element_text(size = NF_MIN_PT, colour = "grey35",
                                           hjust = 0, lineheight = 1.15))
  out <- d[, c("dataset", "val")]
  out$n_acquisitions_in_compartment <- as.integer(acq[as.character(d$lab)])
  out$n_biological_replicates <- n_animals
  out$encoding_note <- paste0(
    "descriptive violin of the acquisition-level distribution with median and ",
    "interquartile crossbar and all acquisitions overplotted; the biological ",
    "replicate is the animal, not the acquisition, and no group comparison is ",
    "made in this panel")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# F2 e. Compartment identity, with censoring made visible
# ==========================================================================
#
# Defect 4. The plotted column is capped for display, so the three strongest
# soma cells were drawn as the identical darkest red and the released sidecar
# carried only the capped value. The heatmap stays - the audit confirmed a dot
# plot would overlap in 6 of 10 rows - but censored tiles are now marked and
# BOTH the displayed and the true value are written to the sidecar.
e8_compartment <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  d <- nv_read_csv(repo_path(panel$primary_source))
  mcol <- intersect(c("marker_label", "marker_gene", "marker", "GeneSymbol",
                      "gene"), names(d))[1]
  ccol <- intersect(c("compartment", "dataset", "celltype_layer",
                      "qc_compartment"), names(d))[1]
  vcol <- intersect(c("displayed_centered_log2", "median_centered_log2",
                      "z", "value"), names(d))[1]
  if (any(is.na(c(mcol, ccol, vcol))))
    stop("e8_compartment: cannot locate marker/compartment/value columns",
         call. = FALSE)
  tcol <- intersect(c("median_centered_log2", "z", "value"), names(d))[1]
  map <- c(neuron_neuropil = "neuron_neuropil", neuropil = "neuron_neuropil",
           neuron_soma = "neuron_soma", soma = "neuron_soma",
           microglia = "microglia")
  d$ds <- unname(map[as.character(d[[ccol]])])
  d <- d[!is.na(d$ds), , drop = FALSE]
  d$marker <- as.character(d[[mcol]])
  d$val <- as.numeric(d[[vcol]])
  d$true_val <- if (!is.na(tcol)) as.numeric(d[[tcol]]) else d$val
  d$comp <- factor(sg_compartment_label(d$ds), levels = sg_compartments()$short)
  mk <- unique(d$marker)
  d$ypos <- match(d$marker, rev(mk))
  d$xpos <- as.integer(d$comp)
  lim <- max(abs(d$val), na.rm = TRUE) * c(-1, 1)
  d$censored <- is.finite(d$true_val) & abs(d$true_val) > abs(lim[2]) + 1e-9
  cen <- d[d$censored, , drop = FALSE]

  p <- ggplot2::ggplot(d, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = val), colour = "white",
                       linewidth = 0.15) +
    nv_diverging(limits = lim, name = NULL) +
    ggplot2::scale_x_continuous(breaks = 1:3, labels = sg_compartments()$short,
                                limits = c(0.5, 3.5), expand = c(0, 0),
                                position = "top") +
    ggplot2::scale_y_continuous(breaks = seq_along(mk), labels = rev(mk),
                                limits = c(0.5, length(mk) + 0.5),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    nf_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT, colour = "grey15",
                                          face = "bold"),
      axis.text.y = ggplot2::element_text(size = NF_MIN_PT, face = "italic"),
      legend.position = "bottom",
      legend.key.height = ggplot2::unit(1.6, "mm"),
      legend.key.width = ggplot2::unit(5, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.caption = ggplot2::element_text(size = NF_MIN_PT, colour = "grey30",
                                           hjust = 0))
  if (nrow(cen)) {
    # a censored tile is marked and its true value printed, so the colour ramp
    # never silently equates values it cannot separate
    p <- p +
      ggplot2::geom_text(data = cen,
                         ggplot2::aes(xpos, ypos,
                                      label = sprintf("%.1f", true_val)),
                         family = fam, size = nf_sz(5.0), colour = "white",
                         fontface = "bold") +
      ggplot2::labs(caption = sprintf(paste0(
        "Colour saturates at %.1f; %d cells exceed it and are printed with ",
        "their true value.
Both values are in the source data."),
        abs(lim[2]), nrow(cen)))
  }
  out <- d[, c("marker", "ds", "comp", "val", "true_val", "censored")]
  names(out) <- c("marker", "dataset", "compartment", "displayed_value",
                  "true_value", "colour_scale_censored")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# ED3 b. CA2-SLM missingness, paired by animal
# ==========================================================================
#
# Sorting 18 acquisitions by value destroyed the two structural facts the page
# depends on. (1) The hemisphere pairing: 755_L was rank 1 and 755_R rank 5, so
# the reader could not see that BOTH QC failures are one hemisphere of an animal
# whose other hemisphere is fine (755: 0.248 vs 0.111; 764: 0.160 vs 0.101) -
# which is exactly the structure panel e's hemisphere-dropping analysis acts on.
# (2) The 9 animals are the biological replicates at n = 3 per group, and group
# was relegated to a fill colour. Animals now sit in three labelled group blocks
# with their two hemispheres joined, so pairing and group are both positional.
e8_ed_ca2_missingness <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  s <- nv_read_csv(repo_path(panel$primary_source))
  s$AnimalID <- as.character(s$AnimalID)
  s$flag <- ifelse(s$qc_flag == "PASS", "", as.character(s$qc_flag))
  grp <- c("CON", "RES", "SUS")
  s$StressGroup <- factor(as.character(s$StressGroup), levels = grp)
  s <- s[order(s$StressGroup, s$AnimalID), , drop = FALSE]

  an <- unique(s[, c("AnimalID", "StressGroup")])
  an <- an[order(an$StressGroup, an$AnimalID), , drop = FALSE]
  # two blank slots separate the three group blocks
  an$slot <- seq_len(nrow(an)) + (as.integer(an$StressGroup) - 1L) * 1.0
  s$slot <- an$slot[match(s$AnimalID, an$AnimalID)]
  s$x <- s$slot + ifelse(s$hemisphere == "L", -0.19, 0.19)

  hdr <- do.call(rbind, lapply(split(an, an$StressGroup), function(z)
    data.frame(g = as.character(z$StressGroup[1]), mid = mean(z$slot),
               lo = min(z$slot) - 0.42, hi = max(z$slot) + 0.42,
               stringsAsFactors = FALSE)))
  ytop <- max(s$fraction_missing_preimputation) * 1.26
  gc <- nv_group_colours()

  p <- ggplot2::ggplot(s) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = fraction_missing_preimputation,
                                    group = AnimalID),
                       colour = "grey72", linewidth = 0.22) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = fraction_missing_preimputation,
                                     colour = StressGroup,
                                     shape = hemisphere), size = 1.15,
                        stroke = 0.4) +
    ggplot2::geom_text(data = s[nzchar(s$flag), , drop = FALSE],
                       ggplot2::aes(x = x, y = fraction_missing_preimputation,
                                    label = flag),
                       vjust = -0.85, family = fam, size = nf_sz(5.0),
                       colour = "grey20") +
    ggplot2::geom_segment(data = hdr,
                          ggplot2::aes(x = lo, xend = hi, y = ytop, yend = ytop),
                          inherit.aes = FALSE, colour = "grey40",
                          linewidth = 0.3) +
    ggplot2::geom_text(data = hdr,
                       ggplot2::aes(x = mid, y = ytop * 1.06, label = g),
                       inherit.aes = FALSE, family = fam, size = nf_sz(5.4),
                       fontface = "bold", colour = "grey15") +
    ggplot2::scale_colour_manual(values = gc, guide = "none") +
    ggplot2::scale_shape_manual(values = c("L" = 1, "R" = 16), name = NULL) +
    ggplot2::scale_x_continuous(breaks = an$slot, labels = an$AnimalID,
                                expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(limits = c(0, ytop * 1.13),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "fraction missing\nbefore imputation",
                  caption = paste0(
                    "One point per acquisition, the two hemispheres of an ",
                    "animal joined. Both QC failures are a single hemisphere ",
                    "of an animal whose other hemisphere passes.")) +
    nf_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
                   axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "right",
                   legend.key.size = ggplot2::unit(2.6, "mm"))
  s$pairing_note <- paste0(
    "animals are the biological replicates (n = 3 per group); hemispheres are ",
    "paired within animal so the unilateral nature of both QC failures is ",
    "visible without reading the caption")
  write_csv_safe(s, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# ED3 e. CA2-SLM sensitivity, as a paired dumbbell on one axis
# ==========================================================================
#
# The scatter put the same quantity on both axes and asked the reader to judge
# distance from the identity line, but the 175 x 54 mm box has no fixed aspect:
# measured from the rendered SVG the x scale was 32.62 mm per log2FC and the y
# scale 11.13 mm per log2FC, a 2.93x anisotropy, so the identity line was drawn
# at 18.8 degrees and the quantity the panel is about - how much effect is lost -
# was compressed threefold against the nuisance axis. One shared axis removes
# the anisotropy entirely and raises the resolution of the loss by ~38%. It also
# names all 28 proteins and blocks the 6 that qualify, so the panel that
# underwrites "28 canonical -> 6 qualified" finally identifies the survivors.
e8_ed_ca2_sensitivity <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  r <- nv_read_csv(repo_path(panel$primary_source))
  z <- data.frame(
    gene = r$gene_symbol, canonical = r$canonical_log2FC_SUS_minus_RES,
    hemi = r$effect_dropping_QC_failed_hemispheres,
    cls = r$CA2_SLM_robustness_class, stringsAsFactors = FALSE)
  z <- z[is.finite(z$canonical) & is.finite(z$hemi), , drop = FALSE]
  z$shrink <- abs(z$canonical) - abs(z$hemi)
  ORD <- c("robust_to_missingness_and_QC", "not_claimable_due_to_QC",
           "insufficient_observed_data")
  z$cls <- factor(z$cls, levels = ORD)
  z <- z[order(z$cls, -abs(z$canonical)), , drop = FALSE]
  n <- nrow(z)
  half <- ceiling(n / 2)
  z$col_id <- rep(c(1L, 2L), times = c(half, n - half))
  z$row <- c(seq_len(half), seq_len(n - half))
  z$ypos <- -z$row
  COL <- c(robust_to_missingness_and_QC = "#1F3D52",
           not_claimable_due_to_QC = "#D1543A",
           insufficient_observed_data = "#B9B6AF")

  p <- ggplot2::ggplot(z) +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey78") +
    ggplot2::geom_segment(ggplot2::aes(x = canonical, xend = hemi, y = ypos,
                                       yend = ypos, colour = cls),
                          linewidth = 0.42, lineend = "round", alpha = 0.75) +
    ggplot2::geom_point(ggplot2::aes(x = canonical, y = ypos, colour = cls),
                        size = 1.05) +
    ggplot2::geom_point(ggplot2::aes(x = hemi, y = ypos, colour = cls),
                        size = 1.05, shape = 21, fill = "white", stroke = 0.4) +
    ggplot2::geom_text(ggplot2::aes(x = -3.15, y = ypos, label = gene),
                       hjust = 0, family = fam, size = nf_sz(5.0),
                       colour = "grey20", fontface = "italic") +
    ggplot2::scale_colour_manual(values = COL, name = NULL,
                                 labels = function(x) gsub("_", " ", x),
                                 guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_x_continuous(limits = c(-3.2, 2.1), breaks = c(-2, -1, 0, 1, 2),
                                expand = c(0, 0)) +
    ggplot2::facet_wrap(~col_id, nrow = 1) +
    ggplot2::labs(
      x = "log2FC (SUS − RES):  filled = canonical,  open = QC-failed hemispheres dropped",
      y = NULL,
      caption = paste0(
        "Segment length is the magnitude lost when the two QC-failed SUS ",
        "acquisitions are removed. Rows are grouped by robustness class, so ",
        "the ", sum(z$cls == ORD[1]), " qualifying proteins read as one block.")) +
    nf_theme(grid = "x") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT),
      strip.text = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.height = ggplot2::unit(1.8, "mm"),
      panel.spacing.x = ggplot2::unit(3.0, "mm"))
  z$reading <- paste0(
    "an open circle to the right of its filled circle means the effect shrank ",
    "towards zero when the QC-failed hemispheres were dropped; both endpoints ",
    "are log2FC on one shared axis, so no aspect-ratio distortion is possible")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# F2 d. The baseline spatial fingerprint, ordered by what it is about
# ==========================================================================
#
# Part-26. The centrepiece of Figure 2 - its largest panel at 175 x 50 mm -
# presented its 19 genes in ALPHABETICAL order (ADCY9, ANKRD63, ATP1A1,
# ATP8A1, CAMK1D, ...), which is a random permutation with respect to the
# spatial structure the panel exists to show. Rows are now seriated by each
# gene's baseline peak spatial unit, so the fingerprint reads as a block
# structure instead of having to be reconstructed cell by cell.
#
# The ordering is PHENOTYPE-BLIND by construction: the peak is taken from
# con_z, which is CON-only baseline abundance, and no stress contrast, group
# label or phenotype statistic enters the sort. Ties break on the gene symbol,
# so the order is deterministic.
e8_fingerprint <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  z <- nv_read_csv(repo_path(panel$primary_source))
  z$sg_unit <- sg_resolve_unit(z$spatial_unit, z$dataset)
  b <- sg_blocks(z$sg_unit, z$dataset)
  o <- b$order
  z$xpos <- match(paste(z$dataset, z$sg_unit), paste(o$dataset, o$unit))

  # each gene's own compartment, and its baseline peak unit within it
  gk <- unique(z[, c("gene", "contrast")])
  gk <- gk[order(gk$contrast, gk$gene), , drop = FALSE]
  gk <- gk[!duplicated(gk$gene), , drop = FALSE]
  pk <- do.call(rbind, lapply(split(z, z$gene), function(g) {
    g <- g[is.finite(g$con_z), , drop = FALSE]
    if (!nrow(g)) return(NULL)
    i <- which.max(g$con_z)
    data.frame(gene = g$gene[1], peak_xpos = g$xpos[i],
               peak_unit = g$sg_unit[i], peak_dataset = g$dataset[i],
               peak_con_z = g$con_z[i], stringsAsFactors = FALSE)
  }))
  ord <- merge(gk, pk, by = "gene", all.x = TRUE, sort = FALSE)
  ord <- ord[order(ord$peak_xpos, ord$gene), , drop = FALSE]
  lev <- ord$gene
  z$ypos <- match(z$gene, rev(lev))
  n <- nrow(o); ny <- length(lev)
  y_comp <- ny + 2.3; y_reg <- ny + 1.05
  lim <- max(abs(z$con_z), na.rm = TRUE) * c(-1, 1)

  p <- ggplot2::ggplot(z, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = con_z), colour = "white",
                       linewidth = 0.1) +
    nv_diverging(limits = lim, name = "Baseline abundance\n(CON z-score)") +
    ggplot2::scale_x_continuous(breaks = seq_len(n),
                                labels = sg_axis_labels(b),
                                limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny), labels = rev(lev),
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
                    "CON animals only (n = 3). Genes are the top hits of each ",
                    "prespecified CON-only anatomical contrast; no stress ",
                    "information enters selection or row order.
Each gene ",
                    "is standardised across its own compartment; rows are ",
                    "ordered by baseline peak unit.")) +
    nf_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT, colour = "grey25"),
      axis.text.y = ggplot2::element_text(size = NF_MIN_PT, face = "italic"),
      legend.title = ggplot2::element_text(size = nf_pt(5.4)),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.6, "mm"),
      legend.key.height = ggplot2::unit(3.6, "mm"))
  cend <- utils::head(b$compartment$end, -1)
  rend <- setdiff(utils::head(b$region$end, -1), cend)
  if (length(cend)) p <- p + ggplot2::annotate(
    "segment", x = cend + 0.5, xend = cend + 0.5, y = 0.5, yend = y_comp - 0.32,
    linewidth = 0.42, colour = "grey25")
  if (length(rend)) p <- p + ggplot2::annotate(
    "segment", x = rend + 0.5, xend = rend + 0.5, y = 0.5, yend = y_reg + 0.4,
    linewidth = 0.18, colour = "grey72")
  out <- merge(z, pk[, c("gene", "peak_unit", "peak_dataset")], by = "gene",
               all.x = TRUE, sort = FALSE)
  out <- out[order(-out$ypos, out$xpos), , drop = FALSE]
  out$row_order_rule <- paste0(
    "rows are seriated by baseline peak spatial unit, taken from CON-only ",
    "con_z; no phenotype information enters the ordering. Ties break on gene ",
    "symbol, so the order is deterministic")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
