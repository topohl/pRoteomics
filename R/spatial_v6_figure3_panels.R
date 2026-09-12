# Part-21 Figure-3 panels: the anatomical program bridge and the hierarchical
# GSEA atlas.
#
# The bridge exists to translate the text labels "CA3-SR", "CA2-SP" and
# "CA1 microglia-enriched ROI" into a spatial hippocampal concept, using the
# SAME drawing language as the Figure-2 schematic. It carries no new quantity:
# every NES and FDR is read from the canonical theme table.
#
# DOWNSTREAM ONLY. No model, no test, no p-value.

# The three prespecified representative findings. Nothing here is chosen from
# the data; this is the exemplar set the manuscript already commits to.
s6_exemplars <- function() {
  data.frame(
    key = c("synaptic", "rna", "oxphos"),
    dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
    unit = c("CA3_sr", "CA2_sp", "CA1"),
    go_id = c("GO:0099536", "GO:0006397", "GO:0006119"),
    program = c("synaptic signalling", "mRNA processing",
                "oxidative phosphorylation"),
    # where the highlight sits on the schematic: arc fraction along CA1->CA3
    arc_frac = c(0.86, 0.45, 0.08),
    band = c("sr", "ca", "microglia"),
    stringsAsFactors = FALSE)
}

# Read the three stored contrasts for one exemplar straight from the canonical
# theme table. Fails closed if a contrast is missing.
s6_exemplar_contrasts <- function(ex, themes) {
  z <- themes[themes$dataset == ex$dataset &
                themes$spatial_unit == ex$unit &
                themes$GO_ID == ex$go_id, , drop = FALSE]
  z <- z[!duplicated(z$contrast), , drop = FALSE]
  want <- c("RES - CON", "SUS - CON", "SUS - RES")
  if (!all(want %in% z$contrast)) {
    stop("exemplar ", ex$key, " is missing contrast(s): ",
         paste(setdiff(want, z$contrast), collapse = ", "), call. = FALSE)
  }
  z <- z[match(want, z$contrast), , drop = FALSE]
  # the theme table carries no setSize column; only NES and GSEA_FDR are needed
  data.frame(key = ex$key, contrast = want, NES = z$NES, FDR = z$GSEA_FDR,
             stringsAsFactors = FALSE)
}

s6_themes <- function() {
  nv_read_csv(repo_path("results", "tables", "10_biological_integration",
                        "gsea_wgcna_concordance", "global",
                        "ontology_aware_gsea_theme_assignments_all_contrasts.csv"))
}

# ------------------------------------------------------ the locator glyph

# One small hippocampus with a single spatial unit highlighted. Same arc
# parameterisation as s4_schematic, so the two figures read as one drawing
# language rather than two different cartoons.
s6_hippo_glyph <- function(ex, accent, fam) {
  arcdf <- function(cx, cy, r, t0, t1, n = 90) {
    t <- seq(t0, t1, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  t0 <- pi * 1.08; t1 <- pi * 2.30
  ca <- arcdf(0, 0, 1.00, t0, t1)
  slm <- arcdf(0, 0, 1.20, pi * 1.10, pi * 2.28)
  sr <- arcdf(0, 0, 1.10, pi * 1.10, pi * 2.28)
  so <- arcdf(0, 0, 0.90, pi * 1.10, pi * 2.28)
  dgu <- arcdf(0.46, -0.28, 0.38, pi * 1.22, pi * 1.98)
  dgl <- arcdf(0.46, -0.28, 0.38, pi * 0.30, pi * 0.86)

  pale <- "grey78"
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = slm, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = sr, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = so, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = ca, ggplot2::aes(x, y), colour = pale,
                       linewidth = 1.1) +
    ggplot2::geom_path(data = dgu, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.9) +
    ggplot2::geom_path(data = dgl, ggplot2::aes(x, y), colour = pale,
                       linewidth = 0.9)

  # highlight the arc segment belonging to this exemplar's region
  seg <- function(d, frac, half = 0.11) {
    n <- nrow(d)
    i0 <- max(1L, floor((frac - half) * n))
    i1 <- min(n, ceiling((frac + half) * n))
    d[i0:i1, , drop = FALSE]
  }
  if (ex$band == "sr") {
    p <- p + ggplot2::geom_path(data = seg(sr, ex$arc_frac),
                                ggplot2::aes(x, y), colour = accent,
                                linewidth = 1.5, lineend = "round")
  } else if (ex$band == "ca") {
    p <- p + ggplot2::geom_path(data = seg(ca, ex$arc_frac),
                                ggplot2::aes(x, y), colour = accent,
                                linewidth = 2.0, lineend = "round")
  } else {
    # microglia-enriched ROI: a local microenvironment sample, drawn as points
    # beside the CA1 band rather than as a band of its own
    t <- t0 + ex$arc_frac * (t1 - t0)
    mg <- data.frame(x = c(1.04, 0.94, 1.12) * cos(t) + c(0, 0.06, -0.05),
                     y = c(1.04, 0.94, 1.12) * sin(t) + c(0.05, -0.06, 0))
    p <- p + ggplot2::geom_point(data = mg, ggplot2::aes(x, y),
                                 colour = accent, size = 1.15)
  }
  # faint region labels, placed on the same arc parameterisation that draws the
  # bands, so the glyph is readable without the header. CA1 is distal (start of
  # the arc), CA3 proximal to the hilus (end).
  at <- function(frac, r) {
    t <- t0 + frac * (t1 - t0)
    c(x = r * cos(t), y = r * sin(t))
  }
  a1 <- at(0.08, 1.40); a2 <- at(0.45, 1.40); a3 <- at(0.86, 1.40)
  rl <- data.frame(
    x = c(a1[["x"]], a2[["x"]], a3[["x"]], 0.46),
    y = c(a1[["y"]], a2[["y"]], a3[["y"]], -0.28),
    l = c("CA1", "CA2", "CA3", "DG"),
    hit = c(ex$arc_frac < 0.25, abs(ex$arc_frac - 0.45) < 0.2,
            ex$arc_frac > 0.7, FALSE),
    stringsAsFactors = FALSE)
  p +
    ggplot2::geom_text(data = rl, ggplot2::aes(x, y, label = l,
                                               colour = hit, fontface = ifelse(hit, 2, 1)),
                       family = fam, size = nv_size(4.6), show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c("TRUE" = accent, "FALSE" = "grey62"),
                                 guide = "none") +
    ggplot2::coord_equal(xlim = c(-1.62, 1.62), ylim = c(-1.46, 1.46),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam)
}

# ------------------------------------------------- the bridge panel itself

s6_f3_anatomy_bridge <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  dcol <- nv_dataset_colours()
  th <- s6_themes()
  ex <- s6_exemplars()

  rows <- do.call(rbind, lapply(seq_len(nrow(ex)), function(i)
    s6_exemplar_contrasts(ex[i, ], th)))
  rows$supported <- rows$FDR < 0.05
  rows$contrast <- factor(rows$contrast,
                          levels = c("RES - CON", "SUS - CON", "SUS - RES"))

  # Trajectory wording is derived from the STORED numbers, and it refuses to
  # call a pattern divergent when the RES arm is not FDR-supported. For the
  # neuropil and soma exemplars the RES-CON arm sits at FDR 0.188 and 0.776,
  # so those read as SUS-specific, not as bidirectional regulation.
  shape_of <- function(k) {
    z <- rows[rows$key == k, ]
    rc <- z[z$contrast == "RES - CON", ]
    sc <- z[z$contrast == "SUS - CON", ]
    if (rc$supported && sc$supported && sign(rc$NES) == sign(sc$NES)) {
      "graded: RES and SUS shift the same way, both supported"
    } else if (!rc$supported && sc$supported) {
      "susceptibility-associated: RES vs CON not FDR-supported"
    } else if (rc$supported && sc$supported) {
      "divergent: RES and SUS shift in opposite directions, both supported"
    } else {
      "not FDR-supported in either group-vs-CON contrast"
    }
  }
  ex$shape <- vapply(ex$key, shape_of, character(1))
  ex$accent <- unname(dcol[ex$dataset])
  ex$unit_label <- sg_unit_label(ex$unit, ex$dataset)
  ex$comp_label <- sg_compartment_label(ex$dataset)

  lim <- max(abs(rows$NES)) * 1.28
  strips <- lapply(seq_len(nrow(ex)), function(i) {
    z <- rows[rows$key == ex$key[i], ]
    ggplot2::ggplot(z, ggplot2::aes(NES, contrast)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = nv_lw("reference_pt"),
                          colour = "grey55") +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = contrast),
                            colour = ex$accent[i],
                            linewidth = nv_lw("reference_pt")) +
      ggplot2::geom_point(ggplot2::aes(shape = supported), size = 1.5,
                          fill = "white", colour = ex$accent[i], stroke = 0.45) +
      ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21),
                                  guide = "none") +
      ggplot2::scale_x_continuous(limits = c(-lim, lim),
                                  breaks = c(-2, 0, 2)) +
      ggplot2::scale_y_discrete(limits = rev(levels(z$contrast))) +
      ggplot2::labs(x = NULL, y = NULL) +
      nv_theme(grid = "x") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4.9),
                     axis.text.x = ggplot2::element_text(size = 4.7))
  })

  glyphs <- lapply(seq_len(nrow(ex)), function(i)
    s6_hippo_glyph(ex[i, ], ex$accent[i], fam))

  heads <- lapply(seq_len(nrow(ex)), function(i) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 1, hjust = 0.5, vjust = 1,
                        family = fam, size = nv_size(5.5), fontface = "bold",
                        colour = "grey12", label = ex$comp_label[i]) +
      ggplot2::annotate("text", x = 0, y = 0.52, hjust = 0.5, vjust = 1,
                        family = fam, size = nv_size(5.2), colour = ex$accent[i],
                        label = ex$unit_label[i]) +
      ggplot2::annotate("text", x = 0, y = 0.06, hjust = 0.5, vjust = 1,
                        family = fam, size = nv_size(4.9), colour = "grey25",
                        label = ex$program[i]) +
      ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-0.5, 1.05),
                               expand = FALSE) +
      ggplot2::theme_void(base_family = fam)
  })

  notes <- lapply(seq_len(nrow(ex)), function(i) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = -1, y = 0.5, hjust = 0, vjust = 0.5,
                        family = fam, size = nv_size(4.5), colour = "grey35",
                        lineheight = 1.15,
                        label = paste(strwrap(ex$shape[i], 34), collapse = "\n")) +
      ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::theme_void(base_family = fam)
  })

  cols <- lapply(seq_len(nrow(ex)), function(i)
    patchwork::wrap_plots(heads[[i]], glyphs[[i]], strips[[i]], notes[[i]],
                          ncol = 1, heights = c(0.26, 0.40, 0.24, 0.13)))
  p <- patchwork::wrap_plots(cols, nrow = 1) +
    patchwork::plot_annotation(
      caption = paste0(
        "Filled point = FDR < 0.05; open point = not FDR-supported. NES and ",
        "FDR are read from the canonical enrichment table and are not ",
        "recomputed.\nThe highlighted band marks where the sample was taken; ",
        "it is not a quantitative anatomical map."),
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(family = fam, size = 4.6,
                                             colour = "grey35", hjust = 0,
                                             lineheight = 1.15)))

  out <- merge(rows, ex[, c("key", "dataset", "unit", "program", "shape")],
               by = "key")
  out$evidence <- "canonical ontology-aware GSEA theme table; values read, not recomputed"
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---------------------------------------------------------- hierarchical atlas

# The cross-compartment GSEA atlas with nested compartment / region / layer
# headers instead of compound strings on the x-axis.
s6_f3_gsea_atlas <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  th <- s6_themes()
  contrast <- as.character(panel$contrast %||% "SUS - RES")
  z <- th[th$contrast == contrast &
            th$assignment_status %in% c("single_theme", "multi_theme", "qc_review") &
            nzchar(as.character(th$theme_id)), , drop = FALSE]
  # seven GO terms are dual-assigned to two themes; deduplicate on the true key
  # before any count so those terms cannot be counted twice
  z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID,
                           z$theme_id)), , drop = FALSE]

  key <- paste(z$dataset, z$spatial_unit, z$theme_id, sep = "\r")
  cells <- do.call(rbind, lapply(split(seq_len(nrow(z)), key), function(ix) {
    w <- z[ix, , drop = FALSE]
    data.frame(dataset = w$dataset[1], spatial_unit = w$spatial_unit[1],
               theme_id = w$theme_id[1],
               manuscript_theme = w$manuscript_theme[1],
               n_terms = nrow(w),
               n_fdr = sum(is.finite(w$GSEA_FDR) & w$GSEA_FDR < 0.05),
               median_NES = stats::median(w$NES, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(cells) <- NULL
  cells$has_support <- cells$n_fdr > 0L
  cells$sg_unit <- sg_resolve_unit(cells$spatial_unit, cells$dataset)

  blocks <- sg_blocks(cells$sg_unit, cells$dataset)
  ord <- blocks$order
  cells$xpos <- match(paste(cells$dataset, cells$sg_unit),
                      paste(ord$dataset, ord$unit))
  # cytoskeleton_structure and epithelial_epidermal_qc carry
  # theme_claim_eligible = FALSE and theme_role = "qc_review": they are
  # contamination/QC review categories, not claimable biology. Drawing them
  # beside the six primary themes would present them as programs.
  th_order <- c("synaptic_signaling_vesicle", "rna_processing_splicing_rnp",
                "ribosome_translation", "mitochondrial_respiration_oxphos",
                "autophagy_lysosome_endosome", "chromatin_organization")
  dropped <- setdiff(unique(cells$theme_id), th_order)
  th_order <- th_order[th_order %in% cells$theme_id]
  cells <- cells[cells$theme_id %in% th_order, , drop = FALSE]
  cells$ypos <- match(cells$theme_id, rev(th_order))
  lab <- vapply(rev(th_order), function(t)
    cells$manuscript_theme[match(t, cells$theme_id)], character(1))

  n <- nrow(ord); ny <- length(th_order)
  y_comp <- ny + 2.2; y_reg <- ny + 1.05
  lim <- max(abs(cells$median_NES), na.rm = TRUE) * c(-1, 1)

  p <- ggplot2::ggplot(cells, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = median_NES), colour = "white",
                       linewidth = 0.12) +
    ggplot2::geom_point(data = cells[cells$has_support, , drop = FALSE],
                        ggplot2::aes(xpos, ypos), size = 0.42,
                        colour = "grey10") +
    nv_diverging(limits = lim, name = "median\nNES") +
    ggplot2::scale_x_continuous(breaks = seq_len(n),
                                labels = sg_axis_labels(blocks),
                                limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny), labels = lab,
                                limits = c(0.5, y_comp + 0.9), expand = c(0, 0)) +
    ggplot2::annotate("segment", x = blocks$compartment$start - 0.5,
                      xend = blocks$compartment$end + 0.5,
                      y = y_comp - 0.32, yend = y_comp - 0.32,
                      linewidth = 0.45, colour = "grey25") +
    ggplot2::annotate("text", x = blocks$compartment$mid, y = y_comp,
                      label = blocks$compartment$label, family = fam,
                      size = nv_size(5.4), fontface = "bold", colour = "grey15") +
    ggplot2::annotate("text", x = blocks$region$mid, y = y_reg,
                      label = blocks$region$label, family = fam,
                      size = nv_size(5.0), colour = "grey30") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    "Contrast: ", contrast,
                    ". Dot = at least one constituent GO term is FDR-supported ",
                    "in that cell. Colour is the median NES of the theme's\n",
                    "canonical GO-BP terms; theme aggregation is an ",
                    "interpretation layer and is NOT a new FDR family. ",
                    length(dropped), " qc_review theme(s) not claim-eligible ",
                    "are excluded.")) +
    nv_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 4.8, colour = "grey25"),
      axis.text.y = ggplot2::element_text(size = 5),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.8, "mm"),
      legend.key.height = ggplot2::unit(4.2, "mm"),
      plot.caption = ggplot2::element_text(size = 4.6, colour = "grey35",
                                           hjust = 0, lineheight = 1.15))
  cend <- utils::head(blocks$compartment$end, -1)
  rend <- setdiff(utils::head(blocks$region$end, -1), cend)
  if (length(cend)) {
    p <- p + ggplot2::annotate("segment", x = cend + 0.5, xend = cend + 0.5,
                               y = 0.5, yend = y_comp - 0.32,
                               linewidth = 0.42, colour = "grey25")
  }
  if (length(rend)) {
    p <- p + ggplot2::annotate("segment", x = rend + 0.5, xend = rend + 0.5,
                               y = 0.5, yend = y_reg + 0.45,
                               linewidth = 0.18, colour = "grey72")
  }
  write_csv_safe(cells, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---------------------------------------------------------- GSEA curve, v6
#
# Adapted from s4_gsea_curve (Part 19, frozen). The curve, the tick rug, the
# three-contrast strip and the exact reconstruction are unchanged - this calls
# the SAME verified s4_gsea_scores(), which hard-stops unless it reproduces the
# stored setSize and enrichmentScore.
#
# The one substantive change: the Part-19 header names only the PROGRAM
# ("synaptic signalling"), so a reader of that panel alone cannot tell which
# compartment, region or layer it belongs to. The Part-21 story audit flagged
# exactly this. The header now carries compartment and spatial unit, resolved
# through the shared spatial grammar, so a region-level compartment never gains
# a fabricated layer.
s6_gsea_curve <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  pr <- s4_programs()
  prog <- pr[pr$key == as.character(panel$program_key), , drop = FALSE]
  if (!nrow(prog)) stop("unknown program key: ", panel$program_key, call. = FALSE)
  ev <- s4_gsea_scores(prog)
  fam <- nv_palette()$typography$family
  acc <- prog$accent[1]
  N <- ev$N

  curve <- data.frame(rank = seq_len(N), es = as.numeric(ev$runes))
  curve <- curve[unique(c(seq(1, N, by = 3), ev$peak, N)), ]
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

  th <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]))
  tr <- th[th$dataset == prog$dataset[1] & th$spatial_unit == prog$unit[1] &
             th$GO_ID == prog$term[1], , drop = FALSE]
  tr <- tr[match(c("RES - CON", "SUS - CON", "SUS - RES"), tr$contrast), ]
  tr <- tr[!is.na(tr$contrast), ]
  tr$short <- sub(" - ", "\u2212", tr$contrast)
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

  comp <- sg_compartment_label(prog$dataset[1])
  unit <- sg_unit_label(prog$unit[1], prog$dataset[1])
  hdr <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0.98, hjust = 0, vjust = 1, family = fam,
                      size = nv_size(5.2), fontface = "bold", colour = "grey12",
                      label = paste0(comp, "   ", unit)) +
    ggplot2::annotate("text", x = 0, y = 0.62, hjust = 0, vjust = 1, family = fam,
                      size = nv_size(5.4), fontface = "bold", colour = acc,
                      label = prog$label[1]) +
    ggplot2::annotate("text", x = 0, y = 0.26, hjust = 0, vjust = 1, family = fam,
                      size = nv_size(4.8), colour = "grey25",
                      label = sprintf("NES %.2f   FDR %.0e   %d/%d leading edge",
                                      ev$NES, ev$FDR, length(ev$leading),
                                      ev$setSize)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0.5, 1, 0, 1, "mm"))

  p <- patchwork::wrap_plots(hdr, top, mid, bot, ncol = 1,
                             heights = c(0.40, 1, 0.10, 0.26))

  out <- data.frame(
    program = prog$label[1], term_id = prog$term[1],
    term_description = ev$description, dataset = prog$dataset[1],
    compartment_label = comp, spatial_unit = prog$unit[1],
    spatial_unit_label = unit,
    has_layer_resolution = sg_has_layer_resolution(prog$unit[1], prog$dataset[1]),
    setSize = ev$setSize, n_leading_edge = length(ev$leading),
    NES = ev$NES, FDR = ev$FDR, ES_peak_rank = ev$peak, n_ranked_genes = N,
    RES_CON_NES = tr$NES[tr$contrast == "RES - CON"][1],
    RES_CON_FDR = tr$GSEA_FDR[tr$contrast == "RES - CON"][1],
    SUS_CON_NES = tr$NES[tr$contrast == "SUS - CON"][1],
    SUS_CON_FDR = tr$GSEA_FDR[tr$contrast == "SUS - CON"][1],
    SUS_RES_NES = tr$NES[tr$contrast == "SUS - RES"][1],
    SUS_RES_FDR = tr$GSEA_FDR[tr$contrast == "SUS - RES"][1],
    reconstruction_note = paste0(
      "running enrichment reconstructed from the stored ranked statistic and ",
      "the pinned org.Mm.eg.db GO membership; s4_gsea_scores() hard-stops ",
      "unless it reproduces the stored setSize and enrichmentScore. Membership ",
      "retrieval is an annotation lookup, NOT an enrichment test: no ",
      "permutation is run and no p-value or FDR is recomputed."),
    stringsAsFactors = FALSE)
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
