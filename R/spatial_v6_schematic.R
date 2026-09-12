# Figure 2a: a spatially CORRECT sampling schematic (Part 21).
#
# Adapted from s4_schematic (Part 19, frozen). The Part-21 story audit found two
# ways in which the Part-19 drawing misstates the design. Because this panel is
# the only map of the 18 units that a reader gets, both errors propagate into
# every later figure:
#
#   1. DG_mo and DG_po are NEUROPIL units, but were drawn in the soma colour, so
#      2 of the 10 neuropil units were visually assigned to the wrong
#      compartment.
#   2. The dashed SO / SR / SLM bands ran continuously from CA1 through CA3,
#      which implies a CA3_slm unit. No CA3_slm exists in the analysis: CA3
#      carries only SO and SR.
#
# This version terminates the SLM band before CA3, draws the DG molecular and
# polymorph layers in the neuropil colour, anchors one microglia-enriched ROI
# marker per region, labels the somatic layers, and states the 10 + 4 + 4 split.
#
# No histology is fabricated: this remains an explicitly declared vector
# schematic. DOWNSTREAM ONLY - no model, no test, no p-value.

s6_schematic <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  dcol <- nv_dataset_colours()
  np <- unname(dcol[["neuron_neuropil"]])
  so_col <- unname(dcol[["neuron_soma"]])
  mg_col <- unname(dcol[["microglia"]])

  arcdf <- function(cx, cy, r, t0, t1, n = 140) {
    t <- seq(t0, t1, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  t0 <- pi * 1.08
  t1 <- pi * 2.30
  span <- t1 - t0
  # CA1 sits near 0.08 of the arc, CA2 near 0.45, CA3 near 0.86. SLM is measured
  # in CA1 and CA2 only, so its band must stop before CA3 starts.
  slm_end <- t0 + 0.66 * span

  build <- function(flip, dx) {
    s <- if (flip) -1 else 1
    tr <- function(d) data.frame(x = s * d$x + dx, y = d$y)
    gap <- data.frame(x = NA_real_, y = NA_real_)
    list(
      ca  = tr(arcdf(0, 0, 1.00, t0, t1)),
      slm = tr(arcdf(0, 0, 1.20, pi * 1.10, slm_end)),
      sr  = tr(arcdf(0, 0, 1.10, pi * 1.10, pi * 2.28)),
      so  = tr(arcdf(0, 0, 0.90, pi * 1.10, pi * 2.28)),
      dgg = tr(rbind(arcdf(0.42, -0.22, 0.24, pi * 1.22, pi * 1.98), gap,
                     arcdf(0.42, -0.22, 0.24, pi * 0.30, pi * 0.86))),
      dgm = tr(rbind(arcdf(0.42, -0.22, 0.34, pi * 1.22, pi * 1.98), gap,
                     arcdf(0.42, -0.22, 0.34, pi * 0.30, pi * 0.86))),
      dgp = tr(arcdf(0.42, -0.22, 0.11, pi * 1.30, pi * 1.90)),
      s = s, dx = dx)
  }
  L <- build(FALSE, -1.55)
  R <- build(TRUE, 1.55)

  lyr <- function(h, d, col, lw, lt = "solid") {
    ggplot2::geom_path(data = h[[d]], ggplot2::aes(x, y), colour = col,
                       linewidth = lw, linetype = lt, lineend = "round")
  }
  p <- ggplot2::ggplot()
  for (h in list(L, R)) {
    p <- p +
      lyr(h, "slm", np, 0.5, "22") +
      lyr(h, "sr", np, 0.5, "22") +
      lyr(h, "so", np, 0.5, "22") +
      # DG molecular and polymorph layers are NEUROPIL units
      lyr(h, "dgm", np, 0.5, "22") +
      lyr(h, "dgp", np, 0.5, "22") +
      # the pyramidal band and the granule blade carry the SOMA units
      lyr(h, "ca", so_col, 1.7) +
      lyr(h, "dgg", so_col, 1.4)
  }

  at <- function(frac, r, h) {
    t <- t0 + frac * span
    c(x = h$dx + h$s * (r * cos(t)), y = r * sin(t))
  }
  ca1 <- at(0.08, 1.38, R)
  ca2 <- at(0.45, 0.66, R)
  ca3 <- at(0.86, 1.38, R)
  lab <- data.frame(
    x = c(ca1[["x"]], ca2[["x"]], ca3[["x"]], R$dx + R$s * 0.42),
    y = c(ca1[["y"]], ca2[["y"]], ca3[["y"]], -0.22),
    l = c("CA1", "CA2", "CA3", "DG"), stringsAsFactors = FALSE)

  # Four laminar labels at one angle on four radii bunch together. Draw them
  # instead as a small keyed column beside the left hippocampus, each with a
  # colour-matched swatch in its own compartment colour, so the
  # neuropil-versus-soma assignment of each layer is unambiguous.
  key_x <- L$dx - 1.42
  lyrlab <- data.frame(
    x = key_x + 0.13,
    y = c(0.62, 0.40, 0.18, -0.04),
    l = c("SLM", "SR", "SO", "SP"),
    col = c(np, np, np, so_col), stringsAsFactors = FALSE)
  keyseg <- data.frame(
    x = key_x - 0.16, xend = key_x + 0.02,
    y = lyrlab$y, col = lyrlab$col,
    lt = c("22", "22", "22", "solid"), stringsAsFactors = FALSE)
  dglab <- data.frame(
    x = L$dx + 0.42 + c(-0.46, -0.30, 0.04),
    y = -0.22 + c(0.20, -0.20, 0.02),
    l = c("MO", "SG", "PO"), col = c(np, so_col, np),
    stringsAsFactors = FALSE)

  # one microglia-enriched ROI marker per region per hemisphere, anchored to the
  # region it samples rather than floating free
  mg <- do.call(rbind, lapply(list(L, R), function(h) {
    d <- do.call(rbind, lapply(c(0.08, 0.45, 0.86), function(f) {
      v <- at(f, 0.74, h)
      data.frame(x = v[["x"]], y = v[["y"]])
    }))
    rbind(d, data.frame(x = h$dx + h$s * 0.42, y = -0.52))
  }))

  p <- p +
    ggplot2::geom_point(data = mg, ggplot2::aes(x, y), size = 0.9,
                        shape = 18, colour = mg_col) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x, y, label = l), family = fam,
                       size = nv_size(6), fontface = "bold") +
    ggplot2::geom_segment(data = keyseg,
                          ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          colour = keyseg$col, linetype = keyseg$lt,
                          linewidth = 0.5, lineend = "round") +
    ggplot2::geom_text(data = lyrlab, ggplot2::aes(x, y, label = l),
                       colour = lyrlab$col, family = fam, size = nv_size(4.8),
                       hjust = 0) +
    ggplot2::geom_text(data = dglab, ggplot2::aes(x, y, label = l),
                       colour = dglab$col, family = fam, size = nv_size(4.4)) +
    ggplot2::annotate("text", x = c(-1.55, 1.55), y = 1.50,
                      label = c("left hemisphere", "right hemisphere"),
                      family = fam, size = nv_size(5.4), fontface = "bold") +
    ggplot2::annotate("segment", x = -0.42, xend = 0.42, y = 1.50, yend = 1.50,
                      linewidth = 0.25, colour = "grey60",
                      arrow = grid::arrow(length = ggplot2::unit(0.8, "mm"),
                                          ends = "both", type = "closed")) +
    ggplot2::annotate("text", x = 0, y = 1.26, label = "paired tissue",
                      family = fam, size = nv_size(5), colour = "grey35") +
    ggplot2::annotate("text", x = 0, y = -1.34, family = fam,
                      size = nv_size(4.9), colour = "grey25",
                      label = paste0("10 neuropil (region x layer)  +  ",
                                     "4 soma (region)  +  ",
                                     "4 microglia ROI (region)  =  ",
                                     "18 spatial units")) +
    ggplot2::annotate("text", x = 0, y = -1.60, family = fam,
                      size = nv_size(4.6), colour = "grey45",
                      label = paste0("9 animals (3 CON / 3 RES / 3 SUS). ",
                                     "SLM is measured in CA1 and CA2 only:",
                                     " there is no CA3 SLM unit.")) +

    ggplot2::coord_equal(xlim = c(-3.05, 3.05), ylim = c(-1.78, 1.66),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam)

  out <- sg_order_for()
  out$compartment_colour <- ifelse(
    out$dataset == "neuron_neuropil", np,
    ifelse(out$dataset == "neuron_soma", so_col, mg_col))
  out$drawn_as <- ifelse(
    out$dataset == "neuron_neuropil", "dashed neuropil band",
    ifelse(out$dataset == "neuron_soma", "solid somatic band", "ROI marker"))
  out$note <- paste0(
    "declared vector schematic; no histology is fabricated. DG MO and DG PO ",
    "are neuropil units and are drawn in the neuropil colour. The SLM band ",
    "terminates before CA3 because no CA3 SLM unit exists in the analysis.")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
