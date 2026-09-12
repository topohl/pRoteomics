# Part-21 spatial_v6 panel renderers.
#
# Every panel here uses the shared spatial grammar in R/spatial_grammar_utils.R:
# one ordering, one header style, one set of abbreviations. The x-axis of a
# spatial heatmap carries only the LAYER token (or the region, for region-level
# compartments) because the compartment and region are already carried by the
# header strips above it - so the reader never has to decode a compound string
# like "neuron_neuropil_CA1_SLM".
#
# DOWNSTREAM ONLY. These renderers reshape and label canonical values. They fit
# no model, run no test and compute no p-value.

# --------------------------------------------------------------- helpers

s6_tbl <- function(name) {
  p <- path_results("tables", "manuscript_candidates", "spatial_v6", name)
  if (!file.exists(p)) stop("missing_required_input: ", p, call. = FALSE)
  nv_read_csv(p)
}

# A spatial heatmap with compartment/region headers. `df` must already carry
# sg_unit / dataset columns. `row_col` names the y variable; `value_col` the
# fill. Row order is supplied by the caller so it is never data-dependent.
s6_spatial_heatmap <- function(df, row_col, value_col, row_levels,
                               fill_name, fill_limits = NULL,
                               na_note = NULL, y_text_pt = 5) {
  fam <- nv_palette()$typography$family
  blocks <- sg_blocks(df$sg_unit, df$dataset)
  ord <- blocks$order
  df$xpos <- match(paste(df$dataset, df$sg_unit),
                   paste(ord$dataset, ord$unit))
  df$yrow <- factor(as.character(df[[row_col]]), levels = rev(row_levels))
  df$val <- df[[value_col]]

  n <- nrow(ord)
  # headers live above the tiles, in tile-row units
  ny <- length(row_levels)
  y_comp <- ny + 2.15
  y_reg <- ny + 1.05

  p <- ggplot2::ggplot(df, ggplot2::aes(xpos, as.integer(yrow))) +
    ggplot2::geom_tile(ggplot2::aes(fill = val), colour = "white",
                       linewidth = 0.12) +
    nv_diverging(limits = fill_limits, name = fill_name) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(n), labels = sg_axis_labels(blocks),
      expand = c(0, 0), limits = c(0.5, n + 0.5)) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(ny), labels = rev(row_levels),
      expand = c(0, 0), limits = c(0.5, y_comp + 0.9)) +
    # compartment strip
    ggplot2::annotate("segment", x = blocks$compartment$start - 0.5,
                      xend = blocks$compartment$end + 0.5,
                      y = y_comp - 0.30, yend = y_comp - 0.30,
                      linewidth = 0.45, colour = "grey25") +
    ggplot2::annotate("text", x = blocks$compartment$mid, y = y_comp,
                      label = blocks$compartment$label, family = fam,
                      size = nv_size(5.4), fontface = "bold", colour = "grey15",
                      vjust = 0.4) +
    # region strip
    ggplot2::annotate("text", x = blocks$region$mid, y = y_reg,
                      label = blocks$region$label, family = fam,
                      size = nv_size(5.0), colour = "grey30", vjust = 0.4) +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 4.8, colour = "grey25"),
      axis.text.y = ggplot2::element_text(size = y_text_pt),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.8, "mm"),
      legend.key.height = ggplot2::unit(4.2, "mm"))

  # heavy separators between compartments, light between regions
  cend <- utils::head(blocks$compartment$end, -1)
  rend <- setdiff(utils::head(blocks$region$end, -1), cend)
  if (length(cend)) {
    p <- p + ggplot2::annotate("segment", x = cend + 0.5, xend = cend + 0.5,
                               y = 0.5, yend = y_comp - 0.30,
                               linewidth = 0.42, colour = "grey25")
  }
  if (length(rend)) {
    p <- p + ggplot2::annotate("segment", x = rend + 0.5, xend = rend + 0.5,
                               y = 0.5, yend = y_reg + 0.45,
                               linewidth = 0.18, colour = "grey72")
  }
  if (!is.null(na_note)) {
    p <- p + ggplot2::labs(caption = na_note) +
      ggplot2::theme(plot.caption = ggplot2::element_text(
        size = 4.6, colour = "grey35", hjust = 0, lineheight = 1.15))
  }
  p
}

# =====================================================================
# Figure 2 / ED2: the direct baseline spatial molecular fingerprint
# =====================================================================

# Form B - external prespecified signature scores. This is the form that gives
# the microglia-enriched compartment a prespecified spatial row, because no
# CON-only anatomical contrast exists for microglia anywhere in the repository.
s6_fingerprint_scores <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  z <- s6_tbl("figure2_spatial_fingerprint_scores.csv")
  compact <- isTRUE(as.logical(panel$compact %||% FALSE))
  keep <- if (compact) {
    c("Kaulich CA1", "Kaulich CA2/3", "Kaulich DG", "Kaulich SLM", "Kaulich SR",
      "Synaptic neuropil", "Oligodendrocyte", "Microglia / PVM")
  } else {
    c("Kaulich CA1", "Kaulich CA2/3", "Kaulich DG",
      "Kaulich SO", "Kaulich SP", "Kaulich SR", "Kaulich SLM",
      "Excitatory neuron", "Inhibitory interneuron", "Synaptic neuropil",
      "Astrocyte", "Oligodendrocyte", "Microglia / PVM", "Vascular")
  }
  z <- z[z$signature %in% keep, , drop = FALSE]
  z$row_label <- sub("^Kaulich ", "", z$signature)
  lev <- sub("^Kaulich ", "", keep)
  lim <- max(abs(z$score), na.rm = TRUE) * c(-1, 1)

  p <- s6_spatial_heatmap(
    z, row_col = "row_label", value_col = "score", row_levels = lev,
    fill_name = "mean z", fill_limits = lim,
    na_note = paste0(
      "CON animals only (n = 3). Each protein is standardised across the ",
      "spatial units of its OWN compartment, so colour shows spatial pattern ",
      "within a compartment and\nis not an abundance comparison between ",
      "compartments. Rows are external prespecified signatures; no stress ",
      "group, contrast or FDR enters selection or scoring."),
    y_text_pt = if (compact) 5.2 else 4.9)

  # name the two row families on the panel
  fam <- nv_palette()$typography$family
  n_spatial <- sum(grepl("^Kaulich ", keep))
  if (n_spatial > 0 && n_spatial < length(keep)) {
    cut_at <- length(keep) - n_spatial + 0.5
    p <- p + ggplot2::annotate("segment", x = 0.5, xend = nrow(sg_blocks(
      z$sg_unit, z$dataset)$order) + 0.5, y = cut_at, yend = cut_at,
      linewidth = 0.3, colour = "grey40")
  }
  z$value_definition <- paste0(
    "mean within-protein CON z across the signature's measured members; ",
    "descriptive only, no test")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# Form A - individual named proteins selected by prespecified CON-only
# anatomical contrasts. Concrete and nameable; neuropil and soma only, because
# no CON-only anatomical contrast exists for microglia.
s6_fingerprint_proteins <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  z <- s6_tbl("figure2_spatial_fingerprint_proteins.csv")
  # deterministic row order: by the contrast that selected the gene, then gene
  ord <- unique(z[, c("gene", "contrast")])
  ord <- ord[order(ord$contrast, ord$gene), , drop = FALSE]
  ord <- ord[!duplicated(ord$gene), , drop = FALSE]
  z$row_label <- z$gene
  lev <- ord$gene
  lim <- max(abs(z$con_z), na.rm = TRUE) * c(-1, 1)

  p <- s6_spatial_heatmap(
    z, row_col = "row_label", value_col = "con_z", row_levels = lev,
    fill_name = "CON z", fill_limits = lim,
    na_note = paste0(
      "CON animals only (n = 3). Rows are the top genes of each prespecified ",
      "CON-only anatomical contrast; no stress information enters selection.\n",
      "Each gene is standardised across the spatial units of its own ",
      "compartment. No CON-only anatomical contrast exists for the ",
      "microglia-enriched compartment."),
    y_text_pt = 4.8)
  z$value_definition <- paste0(
    "within-protein CON z across that protein's own compartment units; ",
    "descriptive only, no test")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
