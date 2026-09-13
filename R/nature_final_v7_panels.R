# Nature-final v7 Figure-2 panels.
#
# Every renderer here obeys a hard 5 pt floor at final size: nf_pt() clamps, and
# the build re-measures the emitted SVG and fails if anything slipped below.
#
# DOWNSTREAM ONLY. These reshape and label canonical values. No model, no test,
# no p-value, no new inference.

NF_MIN_PT <- 5.0

# clamp any requested point size to the floor
nf_pt <- function(x) pmax(as.numeric(x), NF_MIN_PT)
# ggplot size units from a clamped point size
nf_sz <- function(x) nv_size(nf_pt(x))

nf_fam <- function() nv_palette()$typography$family

# base theme with the floor applied to every text element
nf_theme <- function(grid = "none") {
  fam <- nf_fam()
  nv_theme(grid = grid) +
    ggplot2::theme(
      text = ggplot2::element_text(family = fam, size = NF_MIN_PT),
      axis.text = ggplot2::element_text(size = NF_MIN_PT, colour = "grey20"),
      axis.title = ggplot2::element_text(size = nf_pt(5.6)),
      legend.text = ggplot2::element_text(size = NF_MIN_PT),
      legend.title = ggplot2::element_text(size = nf_pt(5.4)),
      strip.text = ggplot2::element_text(size = nf_pt(5.4)),
      plot.caption = ggplot2::element_text(size = NF_MIN_PT, colour = "grey30",
                                           hjust = 0, lineheight = 1.2))
}

nf_theme_tile <- function() {
  nf_theme() + ggplot2::theme(axis.line = ggplot2::element_blank(),
                              axis.ticks = ggplot2::element_blank())
}

# ==========================================================================
# a. Anatomy and sampling, redesigned to TEACH the spatial grammar
# ==========================================================================
#
# Part 22 found the previous schematic encoded COMPARTMENT (by colour) but left
# REGION and LAYER purely textual, so it could not act as the key to the
# fingerprint in panel d. This version pairs the anatomical drawing with the
# 18-unit grid IN FINGERPRINT COLUMN ORDER, using the same region headers, the
# same layer abbreviations and the same separator weights. A reader can read a
# location off the drawing and find the matching column block in panel d.
#
# Compartment = glyph style (dashed band / solid band / diamond marker).
# Region      = shared header-and-separator system, not a rainbow.
# Layer       = anatomical position plus an explicit short label.
nf_schematic <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  dcol <- nv_dataset_colours()
  np <- unname(dcol[["neuron_neuropil"]])
  so_col <- unname(dcol[["neuron_soma"]])
  mg_col <- unname(dcol[["microglia"]])

  arcdf <- function(cx, cy, r, t0, t1, n = 160) {
    t <- seq(t0, t1, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  t0 <- pi * 1.06; t1 <- pi * 2.32; span <- t1 - t0
  slm_end <- t0 + 0.66 * span   # SLM exists in CA1 and CA2 only

  gap <- data.frame(x = NA_real_, y = NA_real_)
  ca  <- arcdf(0, 0, 1.00, t0, t1)
  slm <- arcdf(0, 0, 1.24, pi * 1.08, slm_end)
  sr  <- arcdf(0, 0, 1.12, pi * 1.08, pi * 2.30)
  so  <- arcdf(0, 0, 0.88, pi * 1.08, pi * 2.30)
  dgg <- rbind(arcdf(0.44, -0.24, 0.26, pi * 1.22, pi * 1.98), gap,
               arcdf(0.44, -0.24, 0.26, pi * 0.30, pi * 0.86))
  dgm <- rbind(arcdf(0.44, -0.24, 0.38, pi * 1.22, pi * 1.98), gap,
               arcdf(0.44, -0.24, 0.38, pi * 0.30, pi * 0.86))
  dgp <- arcdf(0.44, -0.24, 0.12, pi * 1.30, pi * 1.90)

  at <- function(frac, r) {
    t <- t0 + frac * span
    c(x = r * cos(t), y = r * sin(t))
  }
  # region brackets drawn ON the arc, the same device as the fingerprint headers
  reg <- data.frame(
    region = c("CA1", "CA2", "CA3"),
    f0 = c(0.00, 0.30, 0.64), f1 = c(0.28, 0.62, 1.00),
    stringsAsFactors = FALSE)
  brk <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
    a <- arcdf(0, 0, 1.42, t0 + reg$f0[i] * span, t0 + reg$f1[i] * span, 40)
    a$region <- reg$region[i]; a
  }))
  rlab <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
    v <- at(mean(c(reg$f0[i], reg$f1[i])), 1.52)
    data.frame(x = v[["x"]], y = v[["y"]], l = reg$region[i], stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot() +
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
                       linewidth = 1.5, lineend = "round") +
    ggplot2::geom_path(data = dgg, ggplot2::aes(x, y), colour = so_col,
                       linewidth = 1.2, lineend = "round") +
    ggplot2::geom_path(data = brk, ggplot2::aes(x, y, group = region),
                       colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_text(data = rlab, ggplot2::aes(x, y, label = l), family = fam,
                       size = nf_sz(5.8), fontface = "bold", colour = "grey12")

  # DG bracket and label
  dgb <- arcdf(0.44, -0.24, 0.56, pi * 1.18, pi * 0.30, 40)
  p <- p +
    ggplot2::geom_path(data = dgb, ggplot2::aes(x, y), colour = "grey55",
                       linewidth = 0.3) +
    ggplot2::annotate("text", x = 0.44, y = -0.95, label = "DG", family = fam,
                      size = nf_sz(5.8), fontface = "bold", colour = "grey12")

  # microglia-enriched ROI: one diamond per region
  mg <- do.call(rbind, lapply(c(0.12, 0.46, 0.84), function(f) {
    v <- at(f, 0.70); data.frame(x = v[["x"]], y = v[["y"]])
  }))
  mg <- rbind(mg, data.frame(x = 0.44, y = -0.58))
  p <- p + ggplot2::geom_point(data = mg, ggplot2::aes(x, y), shape = 18,
                               size = 1.15, colour = mg_col)

  # Layer abbreviations live in the key beside this drawing, so repeating them
  # on the arcs only creates collisions at 86 mm. The drawing carries REGION,
  # which is the anchor a reader needs to find a block in panel d.
  p <- p +
    ggplot2::annotate("text", x = 0, y = 1.88, label = "both hemispheres",
                      family = fam, size = nf_sz(5.0), colour = "grey30")

  drawing <- p +
    ggplot2::coord_equal(xlim = c(-1.95, 1.95), ylim = c(-1.70, 2.02),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam)

  # ---- compartment / region / layer key --------------------------------
  #
  # An 18-cell grid cannot hold 5 pt tokens in 46 mm, and it would duplicate the
  # header strip panel d already carries. The key states the grammar in words
  # instead: one line per compartment, giving its glyph, its spatial resolution
  # and the exact abbreviations used everywhere else in the figure.
  kd <- data.frame(
    row = 3:1,
    style = c("neuropil", "soma", "microglia"),
    name = c("Neuropil", "Neuronal soma", "Microglia-enriched ROI"),
    res = c("region x layer", "region", "region"),
    detail = c("SO, SR, SLM in CA1 and CA2; SO, SR in CA3; MO, PO in DG",
               "SP in CA1, CA2, CA3; SG in DG",
               "CA1, CA2, CA3, DG"),
    n = c(10L, 4L, 4L), stringsAsFactors = FALSE)
  kcol <- c(neuropil = np, soma = so_col, microglia = mg_col)
  key <- ggplot2::ggplot(kd) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 0.075, y = row, yend = row,
                                       colour = style),
                          linewidth = ifelse(kd$style == "soma", 1.4, 0.5),
                          linetype = ifelse(kd$style == "neuropil", "22", "solid"),
                          lineend = "round") +
    ggplot2::geom_point(data = kd[kd$style == "microglia", ],
                        ggplot2::aes(x = 0.037, y = row), shape = 18,
                        size = 1.3, colour = mg_col) +
    ggplot2::geom_text(ggplot2::aes(0.105, row + 0.24,
                                    label = sprintf("%s  -  %d units, %s", name, n, res)),
                       hjust = 0, family = fam, size = nf_sz(5.2),
                       fontface = "bold", colour = "grey12") +
    ggplot2::geom_text(ggplot2::aes(0.105, row - 0.18, label = detail),
                       hjust = 0, family = fam, size = nf_sz(5.0),
                       colour = "grey35") +
    ggplot2::scale_colour_manual(values = kcol, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-0.02, 1.62), ylim = c(0.25, 3.95),
                             expand = FALSE) +
    ggplot2::labs(caption = "abbreviations and order as in d") +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = NF_MIN_PT,
                                                        colour = "grey45",
                                                        hjust = 0),
                   plot.margin = ggplot2::margin(1, 1, 0, 1, "mm"))

  pl <- patchwork::wrap_plots(drawing, key, ncol = 2, widths = c(0.40, 0.60))

  out <- sg_order_for()
  out$drawn_as <- ifelse(out$dataset == "neuron_neuropil", "dashed neuropil band",
                  ifelse(out$dataset == "neuron_soma", "solid somatic band",
                         "ROI diamond"))
  out$note <- paste0(
    "declared vector schematic, no histology fabricated. Compartment is glyph ",
    "style, region is the shared bracket/header system, layer is anatomical ",
    "position plus a short label. DG MO and DG PO are neuropil. SLM exists in ",
    "CA1 and CA2 only; there is no CA3 SLM unit.")
  write_csv_safe(out, csv_path)
  nv_save_panel(pl, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# b. Proteome depth - the least dominant main panel
# ==========================================================================
nf_depth_compact <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  d <- nv_read_csv(repo_path(panel$primary_source))
  cn <- intersect(c("Proteins.Identified", "n_proteins", "proteins_identified"),
                  names(d))[1]
  if (is.na(cn)) stop("nf_depth_compact: no protein-count column", call. = FALSE)
  d$val <- as.numeric(d[[cn]])
  # this source names the compartment in celltype_layer / qc_compartment, not
  # in a `dataset` column; excluded runs must not count toward depth
  ccol <- intersect(c("celltype_layer", "qc_compartment", "dataset"), names(d))[1]
  if (is.na(ccol)) stop("nf_depth_compact: no compartment column", call. = FALSE)
  map <- c(neuron_neuropil = "neuron_neuropil", neuropil = "neuron_neuropil",
           neuron_soma = "neuron_soma", soma = "neuron_soma",
           microglia = "microglia")
  d$dataset <- unname(map[as.character(d[[ccol]])])
  if ("exclude" %in% names(d)) d <- d[!(d$exclude %in% TRUE), , drop = FALSE]
  d <- d[!is.na(d$dataset), , drop = FALSE]
  d <- d[is.finite(d$val), , drop = FALSE]
  d$lab <- factor(sg_compartment_label(d$dataset),
                  levels = sg_compartments()$short)
  d <- d[!is.na(d$lab), , drop = FALSE]
  d <- d[order(d$lab, d$val), , drop = FALSE]
  d$off <- (stats::ave(seq_len(nrow(d)), d$lab, FUN = seq_along) %% 5 - 2) * 0.07

  p <- ggplot2::ggplot(d, ggplot2::aes(as.integer(lab) + off, val)) +
    ggplot2::stat_summary(ggplot2::aes(x = as.integer(lab)), fun = stats::median,
                          geom = "crossbar", width = 0.55, linewidth = 0.22,
                          colour = "grey45") +
    ggplot2::geom_point(ggplot2::aes(colour = lab), size = 0.5, alpha = 0.8) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(unname(nv_dataset_colours()[sg_compartment_levels()]),
                               sg_compartments()$short), guide = "none") +
    ggplot2::scale_x_continuous(breaks = 1:3, labels = sg_compartments()$short,
                                limits = c(0.5, 3.5)) +
    ggplot2::labs(x = NULL, y = "proteins identified") +
    nf_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  write_csv_safe(d[, c("dataset", "val")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# c. PCA, simplified for the main figure
# ==========================================================================
#
# Part 22 measured two defects in the previous main PCA: a 20 x 15 mm block of
# exactly 0.0% ink, and a legend line terminating exactly on the panel boundary.
# Both came from a right-hand legend competing with the data for a 42 mm box.
# Here the legend moves to the bottom as a single deterministic row, region is
# NOT encoded (it would not be readable at this size and implying it would be
# dishonest), and the caption says so.
nf_pca_compact <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  d <- nv_read_csv(repo_path(panel$primary_source))
  ve <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]),
                    required = FALSE)
  pc <- function(k) {
    if (is.null(ve)) return(sprintf("PC%d", k))
    v <- ve$variance_explained[ve$PC == paste0("PC", k)]
    if (!length(v)) return(sprintf("PC%d", k))
    sprintf("PC%d (%.0f%%)", k, 100 * v[1])
  }
  ds <- sg_compartment_levels()
  d <- d[d$dataset %in% ds, , drop = FALSE]
  d$lab <- factor(sg_compartment_label(d$dataset), levels = sg_compartments()$short)
  d <- d[order(d$lab, d$Sample), , drop = FALSE]

  p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_point(ggplot2::aes(colour = lab), size = 0.5, stroke = 0,
                        alpha = 0.85) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(unname(nv_dataset_colours()[ds]),
                               sg_compartments()$short), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      nrow = 1, override.aes = list(size = 1.1))) +
    ggplot2::labs(x = pc(1), y = pc(2),
                  caption = "compartment only; region is not encoded at this size") +
    nf_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(-2, 0, 0, 0),
      legend.key.size = ggplot2::unit(2.4, "mm"),
      plot.margin = ggplot2::margin(1, 2, 1, 1, "mm"))
  write_csv_safe(d[, c("Sample", "PC1", "PC2", "dataset")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# d. Baseline spatial molecular fingerprint - the dominant Figure-2 panel
# ==========================================================================
#
# Same phenotype-blind selection as Part 21: rows are the top genes of each
# prespecified CON-only anatomical contrast. Values are within-protein CON z.
# Gene labels are held at the 5 pt floor.
nf_fingerprint <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  z <- nv_read_csv(repo_path(panel$primary_source))
  ord <- unique(z[, c("gene", "contrast")])
  ord <- ord[order(ord$contrast, ord$gene), , drop = FALSE]
  ord <- ord[!duplicated(ord$gene), , drop = FALSE]
  lev <- ord$gene
  z$sg_unit <- sg_resolve_unit(z$spatial_unit, z$dataset)
  b <- sg_blocks(z$sg_unit, z$dataset)
  o <- b$order
  z$xpos <- match(paste(z$dataset, z$sg_unit), paste(o$dataset, o$unit))
  z$ypos <- match(z$gene, rev(lev))
  n <- nrow(o); ny <- length(lev)
  y_comp <- ny + 2.3; y_reg <- ny + 1.05
  lim <- max(abs(z$con_z), na.rm = TRUE) * c(-1, 1)

  p <- ggplot2::ggplot(z, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = con_z), colour = "white",
                       linewidth = 0.1) +
    nv_diverging(limits = lim, name = "CON z") +
    ggplot2::scale_x_continuous(breaks = seq_len(n),
                                labels = sg_axis_labels(b),
                                limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny), labels = rev(lev),
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
                    "CON animals only (n = 3). Genes are the top hits of each ",
                    "prespecified CON-only anatomical contrast; no stress ",
                    "information enters selection. Each gene is standardised ",
                    "across its own compartment.")) +
    nf_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = NF_MIN_PT, colour = "grey25"),
      axis.text.y = ggplot2::element_text(size = NF_MIN_PT, face = "italic"),
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
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# e. Compartment identity - given real width this time
# ==========================================================================
#
# Part 22 found the previous version had rotated axis labels escaping the panel
# box and a legend consuming a quarter of 39 mm. Here the compartment labels are
# horizontal, the legend sits at the bottom, and the box is 78 mm.
nf_compartment <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  d <- nv_read_csv(repo_path(panel$primary_source))
  mcol <- intersect(c("marker_label", "marker_gene", "marker", "GeneSymbol",
                      "gene"), names(d))[1]
  ccol <- intersect(c("compartment", "dataset", "celltype_layer",
                      "qc_compartment"), names(d))[1]
  vcol <- intersect(c("displayed_centered_log2", "median_centered_log2",
                      "z", "value"), names(d))[1]
  if (any(is.na(c(mcol, ccol, vcol)))) {
    stop("nf_compartment: cannot locate marker/compartment/value columns; have: ",
         paste(names(d), collapse = ", "), call. = FALSE)
  }
  map <- c(neuron_neuropil = "neuron_neuropil", neuropil = "neuron_neuropil",
           neuron_soma = "neuron_soma", soma = "neuron_soma",
           microglia = "microglia")
  d$ds <- unname(map[as.character(d[[ccol]])])
  d <- d[!is.na(d$ds), , drop = FALSE]
  d$marker <- as.character(d[[mcol]])
  d$val <- as.numeric(d[[vcol]])
  d$comp <- factor(sg_compartment_label(d$ds), levels = sg_compartments()$short)
  mk <- unique(d$marker)
  d$ypos <- match(d$marker, rev(mk))
  d$xpos <- as.integer(d$comp)
  lim <- max(abs(d$val), na.rm = TRUE) * c(-1, 1)

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
      legend.box.margin = ggplot2::margin(-2, 0, 0, 0))
  write_csv_safe(d[, c("marker", "ds", "val")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# f. Bilateral reproducibility - MAIN version
# ==========================================================================
#
# Main carries the story: coarse regional identity reproduces strongly, fine
# laminar identity less so, with an explicit block separator between the two
# anatomical levels. The complete metric inventory (Pearson, Spearman, sign
# agreement, ICC, precision) stays in ED1 - this panel is NOT that panel.
nf_bilateral_main <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  su <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]))
  rc <- intersect(c("pearson_r", "r", "correlation", "bilateral_pearson"),
                  names(su))[1]
  cc <- intersect(c("contrast", "anatomical_contrast", "comparison"), names(su))[1]
  if (any(is.na(c(rc, cc)))) {
    stop("nf_bilateral_main: cannot locate r/contrast columns; have: ",
         paste(names(su), collapse = ", "), call. = FALSE)
  }
  su$r <- as.numeric(su[[rc]])
  su$contrast <- as.character(su[[cc]])
  su <- su[is.finite(su$r), , drop = FALSE]
  # anatomical level from the stored domain where available, else the contrast
  dcol <- intersect(c("validation_domain", "anatomical_level", "level"),
                    names(su))[1]
  su$level <- if (!is.na(dcol)) {
    ifelse(grepl("strat|lamin|layer", su[[dcol]], ignore.case = TRUE),
           "Fine / laminar identity", "Regional identity")
  } else {
    ifelse(grepl("_strata$|_DG_layers$", su$contrast),
           "Fine / laminar identity", "Regional identity")
  }
  su$level <- factor(su$level,
                     levels = c("Regional identity", "Fine / laminar identity"))
  su <- su[order(su$level, -su$r), , drop = FALSE]
  su$lab <- gsub("_", " ", su$contrast)
  su$lab <- sub(" vs mean other", " vs rest", su$lab)
  su$ypos <- rev(seq_len(nrow(su)))

  p <- ggplot2::ggplot(su, ggplot2::aes(r, ypos)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = r, yend = ypos,
                                       colour = level),
                          linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(colour = level), size = 1) +
    ggplot2::scale_colour_manual(
      values = c("Regional identity" = "#1F3D52",
                 "Fine / laminar identity" = "#C2A878"), name = NULL) +
    ggplot2::scale_y_continuous(breaks = su$ypos, labels = su$lab,
                                expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = "left-right correlation of the anatomical effect",
                  y = NULL) +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "bottom",
                   legend.key.size = ggplot2::unit(2.4, "mm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.box.margin = ggplot2::margin(-2, 0, 0, 0))
  # explicit rule between the two anatomical levels
  nreg <- sum(su$level == "Regional identity")
  if (nreg > 0 && nreg < nrow(su)) {
    yb <- rev(seq_len(nrow(su)))[nreg] - 0.5
    p <- p + ggplot2::annotate("segment", x = -Inf, xend = Inf, y = yb, yend = yb,
                               linewidth = 0.3, colour = "grey55")
  }
  write_csv_safe(su[, c("contrast", "level", "r")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# g. External validation - MAIN version (expected pairings only)
# ==========================================================================
#
# Main shows only the story-driving expected pairings, blocked by anatomical
# level. The complete specificity inventory - every off-target comparison -
# stays in ED2, so the two are not the same graphic.
nf_external_main <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  k <- nv_read_csv(repo_path(panel$primary_source))
  k <- k[k$expected_match %in% TRUE, , drop = FALSE]
  k$level <- ifelse(grepl("strata", k$validation_domain),
                    "CA1 laminar identity", "Regional identity")
  k$level <- factor(k$level, levels = c("Regional identity", "CA1 laminar identity"))
  k$sig <- is.finite(k$p_adjust) & k$p_adjust < 0.05
  k$lab <- paste0(gsub("_", " ", k$internal_contrast), "  ->  ", k$external_signature)
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
    ggplot2::labs(x = "NES vs external hippocampal signature", y = NULL,
                  caption = "expected pairings only; filled = FDR < 0.05. Full specificity inventory in ED2.") +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "none")
  nreg <- sum(k$level == "Regional identity")
  if (nreg > 0 && nreg < nrow(k)) {
    p <- p + ggplot2::annotate("segment", x = -Inf, xend = Inf,
                               y = nreg + 0.5, yend = nreg + 0.5,
                               linewidth = 0.3, colour = "grey55")
  }
  write_csv_safe(k[, c("internal_contrast", "external_signature", "level",
                       "NES", "p_adjust")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ==========================================================================
# h. Internal anatomical program validation - MAIN version
# ==========================================================================
#
# Main shows the single strongest canonical term per contrast so the panel fits
# at >=5 pt. The complete term inventory stays in ED2.
nf_internal_main <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nf_fam()
  g <- nv_read_csv(repo_path(panel$primary_source))
  g$level <- ifelse(grepl("CA1_strata", g$contrast),
                    "CA1 laminar identity", "Regional identity")
  g$level <- factor(g$level, levels = c("Regional identity", "CA1 laminar identity"))
  # one term per contrast: the largest |NES|, deterministic on ties
  g <- g[order(g$contrast, -abs(g$NES), g$Description), , drop = FALSE]
  g <- g[!duplicated(g$contrast), , drop = FALSE]
  g <- g[order(g$level, g$NES), , drop = FALSE]
  g$ypos <- seq_len(nrow(g))
  g$lab <- paste0(gsub("_", " ", g$contrast), "  ->  ", g$Description)
  g$sig <- is.finite(g$p_adjust) & g$p_adjust < 0.05

  p <- ggplot2::ggplot(g, ggplot2::aes(NES, ypos)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey65") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = ypos,
                                       colour = level),
                          linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(colour = level, size = setSize,
                                     shape = sig), fill = "white", stroke = 0.4) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21),
                                guide = "none") +
    ggplot2::scale_size_continuous(range = c(0.6, 1.8), guide = "none") +
    ggplot2::scale_colour_manual(
      values = c("Regional identity" = "#1F3D52",
                 "CA1 laminar identity" = "#C2A878"), name = NULL) +
    ggplot2::scale_y_continuous(breaks = g$ypos, labels = g$lab,
                                expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::labs(x = "NES, internal anatomical GO program", y = NULL,
                  caption = "strongest canonical term per contrast; point size = set size. Full term inventory in ED2.") +
    nf_theme(grid = "x") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = NF_MIN_PT),
                   legend.position = "none")
  nreg <- sum(g$level == "Regional identity")
  if (nreg > 0 && nreg < nrow(g)) {
    p <- p + ggplot2::annotate("segment", x = -Inf, xend = Inf,
                               y = nreg + 0.5, yend = nreg + 0.5,
                               linewidth = 0.3, colour = "grey55")
  }
  write_csv_safe(g[, c("contrast", "Description", "level", "NES", "p_adjust",
                       "setSize")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
