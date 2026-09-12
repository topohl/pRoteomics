# Part-21 spatial_v6 Extended Data renderers.
#
# Same discipline as R/spatial_v6_figure_panels.R: downstream only. These
# reshape and label canonical values. No model, no test, no p-value.

# =====================================================================
# ED7: baseline spatial location vs location of strongest stress effect
# =====================================================================
#
# Both coordinates are ALREADY STORED per protein in the canonical atlas:
#   peak_unit                                = baseline dominant location
#   sus_res_strongest_spatial_unit_canonical = strongest SUS-RES effect location
# Both are written in the same lowercase atlas vocabulary by
# sat_canonical_spatial_unit(), which correctly declines to invent a layer for
# soma or microglia. Nothing is recomputed here.
#
# The displayed set is the QC-QUALIFIED one: the 6 CA2-SLM hits that survived
# the robustness audit, plus the 9 hits that were never in CA2-SLM and so were
# never exposed to that QC problem. The 22 hits the audit could not clear are
# excluded, and the exclusion is stated on the panel.
#
# Deliberately NOT described as redistribution, relocation or migration. A
# protein does not move. The panel makes two separate spatial statements about
# the same protein: where it is most abundant at baseline, and where its
# strongest phenotype-associated difference is.

s6_ed7_qualified <- function() {
  a <- nv_read_csv(repo_path("results", "tables", "11_spatial_systems", "atlas",
                             "protein_spatial_cell_affinity.csv"))
  h <- a[a$is_sus_res_fdr_supported %in% TRUE, , drop = FALSE]
  cls <- as.character(h$CA2_SLM_robustness_class)
  cls[is.na(cls) | !nzchar(cls)] <- "not_in_CA2_SLM_never_at_risk"
  h$qc_class <- cls
  h$qualified <- h$qc_class %in% c("robust_to_missingness_and_QC",
                                   "not_in_CA2_SLM_never_at_risk")
  h$baseline_unit <- sg_resolve_unit(h$peak_unit, h$dataset)
  h$effect_unit <- sg_resolve_unit(h$sus_res_strongest_spatial_unit_canonical,
                                   h$dataset)
  h
}

s6_unit_positions <- function(dataset) {
  u <- sg_units()
  u <- u[u$dataset == dataset, , drop = FALSE]
  stats::setNames(seq_len(nrow(u)), u$analysis_key)
}

s6_ed7_origin_destination <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  h <- s6_ed7_qualified()
  n_all <- nrow(h)
  n_drop <- sum(!h$qualified)
  z <- h[h$qualified, , drop = FALSE]

  z$comp <- factor(sg_compartment_label(z$dataset),
                   levels = sg_compartments()$short)
  z$y0 <- NA_real_
  z$y1 <- NA_real_
  for (d in unique(z$dataset)) {
    m <- s6_unit_positions(d)
    i <- z$dataset == d
    z$y0[i] <- m[z$baseline_unit[i]]
    z$y1[i] <- m[z$effect_unit[i]]
  }
  z <- z[order(z$comp, z$y1, z$y0, z$GeneSymbol), , drop = FALSE]
  z$elsewhere <- z$baseline_unit != z$effect_unit
  # Several proteins can share one destination: 28 of the 31 neuropil hits land
  # in CA2 SLM. Stack their labels deterministically with a short leader rather
  # than overprinting them all on the same point.
  z$lab_y <- z$y1
  for (d in unique(z$dataset)) {
    ii <- which(z$dataset == d)
    for (g in split(ii, z$y1[ii])) {
      k <- length(g)
      if (k > 1L) z$lab_y[g] <- z$y1[g] + (seq_len(k) - (k + 1) / 2) * 0.66
    }
  }

  # y-axis tick labels: each compartment shows only its OWN units
  lab <- do.call(rbind, lapply(unique(z$dataset), function(d) {
    u <- sg_units()
    u <- u[u$dataset == d, , drop = FALSE]
    data.frame(comp = factor(sg_compartment_label(d),
                             levels = sg_compartments()$short),
               ypos = seq_len(nrow(u)), lab = u$display,
               stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot(z) +
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = 2, y = y0, yend = y1,
                                       colour = elsewhere),
                          linewidth = nv_lw("reference_pt"), alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(x = 1, y = y0), size = 0.85,
                        colour = "grey45") +
    ggplot2::geom_point(ggplot2::aes(x = 2, y = y1, colour = elsewhere),
                        size = 1.05) +
    ggplot2::geom_segment(ggplot2::aes(x = 2, xend = 2.055, y = y1,
                                       yend = lab_y),
                          colour = "grey72", linewidth = 0.12) +
    ggplot2::geom_text(ggplot2::aes(x = 2.075, y = lab_y, label = GeneSymbol),
                       hjust = 0, family = fam, size = nv_size(4.6),
                       colour = "grey20") +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = 0.93, y = ypos, label = lab),
                       hjust = 1, family = fam, size = nv_size(4.6),
                       colour = "grey35") +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "#C0442C", "FALSE" = "#9E9A92"), guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2),
      labels = c("baseline\ndominant unit", "strongest\nSUS-RES effect"),
      limits = c(0.24, 2.9), expand = c(0, 0)) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::facet_wrap(~comp, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      x = NULL, y = NULL,
      caption = sprintf(paste0(
        "%d of %d FDR-supported SUS-RES proteins: the %d that survived the ",
        "CA2-SLM robustness audit plus the %d never in CA2-SLM.\nThe %d the ",
        "audit could not clear are excluded. Red = strongest effect outside ",
        "the baseline dominant unit."),
        nrow(z), n_all,
        sum(z$qc_class == "robust_to_missingness_and_QC"),
        sum(z$qc_class == "not_in_CA2_SLM_never_at_risk"), n_drop)) +
    nv_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 4.9, lineheight = 1.05),
      strip.text = ggplot2::element_text(size = 5.4, face = "bold"),
      plot.caption = ggplot2::element_text(size = 4.6, colour = "grey35",
                                           hjust = 0, lineheight = 1.15),
      panel.spacing.x = ggplot2::unit(3.4, "mm"))

  out <- z[, c("GeneSymbol", "dataset", "baseline_unit", "effect_unit",
               "baseline_rank_of_strongest_effect_unit",
               "effect_identity_relationship", "qc_class", "elsewhere")]
  out$reading <- paste0(
    "the strongest phenotype-associated difference is usually NOT in the unit ",
    "where the protein is most abundant at baseline; nothing moves, these are ",
    "two separate spatial statements about the same protein")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# Alternative encoding (brief section 31): a compact marked matrix instead of an
# origin-destination plot, built so the two can be compared directly.
s6_ed7_location_matrix <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  h <- s6_ed7_qualified()
  z <- h[h$qualified, , drop = FALSE]
  z <- z[order(z$dataset, z$effect_unit, z$baseline_unit, z$GeneSymbol), ,
         drop = FALSE]
  blocks <- sg_blocks(c(z$baseline_unit, z$effect_unit),
                      c(z$dataset, z$dataset))
  ord <- blocks$order
  genes <- unique(z$GeneSymbol)
  pos <- function(unit, ds) match(paste(ds, unit), paste(ord$dataset, ord$unit))
  m <- rbind(
    data.frame(gene = z$GeneSymbol, dataset = z$dataset,
               x = pos(z$baseline_unit, z$dataset),
               kind = "baseline dominant", stringsAsFactors = FALSE),
    data.frame(gene = z$GeneSymbol, dataset = z$dataset,
               x = pos(z$effect_unit, z$dataset),
               kind = "strongest SUS-RES effect", stringsAsFactors = FALSE))
  m$y <- match(m$gene, rev(genes))
  n <- nrow(ord)
  ny <- length(genes)
  y_comp <- ny + 2.2
  y_reg <- ny + 1.05

  p <- ggplot2::ggplot(m, ggplot2::aes(x, y)) +
    ggplot2::geom_point(ggplot2::aes(shape = kind, colour = kind), size = 1.3) +
    ggplot2::scale_shape_manual(values = c("baseline dominant" = 1,
                                           "strongest SUS-RES effect" = 4),
                                name = NULL) +
    ggplot2::scale_colour_manual(values = c("baseline dominant" = "grey40",
                                            "strongest SUS-RES effect" = "#C0442C"),
                                 name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq_len(n),
                                labels = sg_axis_labels(blocks),
                                limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny), labels = rev(genes),
                                limits = c(0.5, y_comp + 0.9), expand = c(0, 0)) +
    ggplot2::annotate("segment", x = blocks$compartment$start - 0.5,
                      xend = blocks$compartment$end + 0.5,
                      y = y_comp - 0.32, yend = y_comp - 0.32,
                      linewidth = 0.45, colour = "grey25") +
    ggplot2::annotate("text", x = blocks$compartment$mid, y = y_comp,
                      label = blocks$compartment$label, family = fam,
                      size = nv_size(5.2), fontface = "bold", colour = "grey15") +
    ggplot2::annotate("text", x = blocks$region$mid, y = y_reg,
                      label = blocks$region$label, family = fam,
                      size = nv_size(4.8), colour = "grey30") +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 4.7),
                   axis.text.y = ggplot2::element_text(size = 4.7),
                   legend.position = "bottom",
                   legend.key.size = ggplot2::unit(2.6, "mm"),
                   legend.text = ggplot2::element_text(size = 4.9))
  cend <- utils::head(blocks$compartment$end, -1)
  if (length(cend)) {
    p <- p + ggplot2::annotate("segment", x = cend + 0.5, xend = cend + 0.5,
                               y = 0.5, yend = y_comp - 0.32,
                               linewidth = 0.42, colour = "grey25")
  }
  write_csv_safe(m, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# =====================================================================
# ED8a: what a spatial molecular network actually IS
# =====================================================================
#
# Baseline CON molecular-similarity matrices, anatomically ordered. Reshaped
# from the stored CON edge list; no similarity is recomputed. The title states
# explicitly that this is similarity between spatial proteomic profiles and NOT
# anatomical connectivity, so network edges are never read as projections.
s6_ed8_similarity <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  e <- nv_read_csv(repo_path("results", "tables", "11_spatial_systems",
                             "networks",
                             "CON_spatial_molecular_similarity_matrix.csv"))
  e$a <- sg_resolve_unit(e$node_a, e$dataset)
  e$b <- sg_resolve_unit(e$node_b, e$dataset)
  # symmetric long form, both triangles, so the matrix reads as a matrix
  sym <- rbind(
    data.frame(dataset = e$dataset, r = e$a, c = e$b,
               v = e$median_similarity, stringsAsFactors = FALSE),
    data.frame(dataset = e$dataset, r = e$b, c = e$a,
               v = e$median_similarity, stringsAsFactors = FALSE))
  # explicit self-similarity, drawn as a neutral diagonal rather than left blank
  dg <- do.call(rbind, lapply(unique(e$dataset), function(d) {
    u <- sg_units(); u <- u[u$dataset == d, , drop = FALSE]
    data.frame(dataset = d, r = u$analysis_key, c = u$analysis_key,
               v = NA_real_, stringsAsFactors = FALSE)
  }))
  sym <- rbind(sym, dg)
  sym$comp <- factor(sg_compartment_label(sym$dataset),
                     levels = sg_compartments()$short)
  # All 18 analysis keys are unique across compartments (CA1_so / CA1_sp / CA1),
  # so one factor with free facet scales gives each compartment its own axis
  # showing only its own units - no manual text placement needed.
  lev <- sg_units()$analysis_key
  sym$xf <- factor(sym$c, levels = lev)
  sym$yf <- factor(sym$r, levels = rev(lev))
  disp <- stats::setNames(sg_units()$display, sg_units()$analysis_key)
  relab <- function(x) unname(disp[as.character(x)])
  lim <- max(abs(sym$v), na.rm = TRUE) * c(-1, 1)

  p <- ggplot2::ggplot(sym[!is.na(sym$v), , drop = FALSE],
                       ggplot2::aes(xf, yf)) +
    ggplot2::geom_tile(ggplot2::aes(fill = v), colour = "white",
                       linewidth = 0.12) +
    nv_diverging(limits = lim, name = "median\nsimilarity") +
    ggplot2::scale_x_discrete(labels = relab, drop = TRUE) +
    ggplot2::scale_y_discrete(labels = relab, drop = TRUE) +
    ggplot2::facet_wrap(~comp, nrow = 1, scales = "free") +
    ggplot2::labs(
      x = NULL, y = NULL,
      caption = paste0(
        "Molecular similarity between spatial proteomic profiles. ",
        "NOT anatomical connectivity: a value is the correlation between two ",
        "units' protein profiles,\nnot a projection. CON animals only ",
        "(n = 3); the stored median across CON animals. The diagonal is ",
        "self-similarity and is left empty.")) +
    nv_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 4.4, angle = 90, hjust = 1,
                                          vjust = 0.5, colour = "grey30"),
      axis.text.y = ggplot2::element_text(size = 4.4, colour = "grey30"),
      strip.text = ggplot2::element_text(size = 5.4, face = "bold"),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.8, "mm"),
      legend.key.height = ggplot2::unit(4, "mm"),
      plot.caption = ggplot2::element_text(size = 4.6, colour = "grey35",
                                           hjust = 0, lineheight = 1.15),
      panel.spacing.x = ggplot2::unit(3, "mm"))

  out <- sym[!is.na(sym$v), c("dataset", "r", "c", "v")]
  names(out) <- c("dataset", "unit_row", "unit_col", "median_similarity")
  out$reading <- paste0(
    "similarity between spatial proteomic PROFILES, not anatomical ",
    "connectivity; CON animals only")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
