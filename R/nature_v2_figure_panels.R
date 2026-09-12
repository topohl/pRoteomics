# =====================================================================
# Nature-style v2 panel renderers.
#
# Every renderer has signature (panel, svg_path, csv_path, w_mm, h_mm) and
# draws at EXACTLY w_mm x h_mm so no text is ever rescaled.
#
# House rules enforced by construction: no plot title, no subtitle, no
# gridlines unless the reading demands them, at most one legend per panel,
# one semantic palette, descriptive information greyed relative to supported
# information. Explanatory text belongs in the figure legend, not the artwork.
#
# Downstream only: nothing here fits a model, runs an enrichment or adjusts a
# p-value. nv_assert_no_model_fitting() scans this file.
# =====================================================================

# ============================================================ FIGURE 2

# a. Spatial / experimental anchor.
# No histology is fabricated. If no suitable source image exists this draws a
# clean schematic of the sampling design and is explicitly labelled as needing
# final Illustrator artwork.
nvp_spatial_anchor <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  dcol <- nv_dataset_colours()

  # A DATA-BACKED SAMPLING MAP, not invented anatomy. Every filled cell is a
  # spatial unit that was actually acquired, read from the canonical sample
  # metadata. The hippocampal drawing itself is deferred to Illustrator.
  meta <- as.data.frame(readxl::read_excel(
    repo_path("data", "metadata", "TPE9_sample_metadata_males.xlsx")))
  meta <- meta[!(meta$exclude %in% TRUE), , drop = FALSE]
  meta$region <- toupper(as.character(meta$region))
  g <- unique(meta[, c("region", "layer", "celltype_layer")])
  g <- g[g$region %in% c("CA1", "CA2", "CA3", "DG"), , drop = FALSE]
  g$row <- ifelse(g$celltype_layer == "microglia", "microglia",
                  ifelse(g$celltype_layer == "neuron_soma",
                         paste0("soma ", g$layer), paste0("neuropil ", g$layer)))
  n_an <- length(unique(meta$AnimalID))
  n_hemi <- length(unique(meta$ReplicateGroup))

  row_order <- c(paste0("neuropil ", c("slm", "sr", "so", "mo", "po")),
                 paste0("soma ", c("sp", "sg")), "microglia")
  g$row <- factor(g$row, levels = rev(intersect(row_order, unique(g$row))))
  g$region <- factor(g$region, levels = c("CA1", "CA2", "CA3", "DG"))
  g$compartment <- ifelse(grepl("^neuropil", as.character(g$row)), "neuron_neuropil",
                          ifelse(grepl("^soma", as.character(g$row)), "neuron_soma",
                                 "microglia"))

  p <- ggplot2::ggplot(g, ggplot2::aes(region, row, fill = compartment)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4, width = 0.9, height = 0.85) +
    ggplot2::scale_fill_manual(values = dcol, labels = nv_dataset_label, name = NULL) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = sprintf(
                    "%d spatial units · left + right hemisphere · %d animals (3 CON / 3 RES / 3 SUS) · AnimalID = biological replicate\nsampling map from canonical metadata; hippocampal artwork to be finalised in Illustrator",
                    nrow(g), n_an)) +
    nv_theme_tile() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.height = ggplot2::unit(2, "mm"),
      axis.text.x = ggplot2::element_text(face = "bold", size = 6),
      axis.text.y = ggplot2::element_text(size = 5.2),
      plot.caption = ggplot2::element_text(size = 5, hjust = 0, colour = "grey30",
                                           lineheight = 1.15,
                                           margin = ggplot2::margin(1.5, 0, 0, 0, "mm")))

  out <- g[, c("region", "layer", "celltype_layer", "row", "compartment")]
  out$n_animals <- n_an
  out$n_hemispheres <- n_hemi
  out$artwork_status <- paste0(
    "sampling map derived from canonical sample metadata; no histology or ",
    "anatomical imagery is fabricated; final anatomical artwork to be drawn ",
    "in Illustrator")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# b. Compact proteome depth.
nvp_depth_compact <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  prot <- intersect(c("Proteins.Identified", "n_proteins_detected", "Proteins"), names(d))[1]
  ds <- intersect(c("celltype_layer", "dataset"), names(d))[1]
  if (is.na(prot) || is.na(ds)) stop("depth source lacks protein-count or dataset column")
  z <- data.frame(dataset = as.character(d[[ds]]),
                  n = suppressWarnings(as.numeric(d[[prot]])), stringsAsFactors = FALSE)
  z <- z[is.finite(z$n) & z$dataset %in% names(nv_dataset_colours()), , drop = FALSE]
  z$lab <- nv_dataset_label(z$dataset)
  z$lab <- factor(z$lab, levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))

  p <- ggplot2::ggplot(z, ggplot2::aes(lab, n, fill = dataset)) +
    ggplot2::geom_boxplot(width = 0.5, linewidth = nv_lw("axis_pt"),
                          outlier.size = 0.25, outlier.colour = "grey55") +
    ggplot2::scale_fill_manual(values = nv_dataset_colours(), guide = "none") +
    ggplot2::scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
    ggplot2::labs(x = NULL, y = "proteins per sample") +
    nv_theme(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# c. Paired hemisphere reproducibility. Minimal annotation: r and n only.
nvp_bilateral <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  prot <- nv_read_csv(repo_path(panel$primary_source))
  summ <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]))
  pick <- do.call(rbind, lapply(split(summ, summ$dataset), function(d) {
    d <- d[order(d$pearson_r), , drop = FALSE]
    d[ceiling(nrow(d) / 2), , drop = FALSE]
  }))
  pick$lab <- nv_dataset_label(pick$dataset)
  dat <- merge(prot, pick[, c("dataset", "contrast", "lab", "pearson_r",
                              "sign_agreement_fraction", "n_evaluable_proteins")],
               by = c("dataset", "contrast"))
  dat <- dat[is.finite(dat$estimate_L) & is.finite(dat$estimate_R), , drop = FALSE]
  dat$lab <- factor(dat$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))

  ann <- unique(dat[, c("lab", "pearson_r", "n_evaluable_proteins")])
  ann$t <- sprintf("r = %.2f", ann$pearson_r)
  lim <- c(-3.2, 3.2)
  ann$x <- lim[1] * 0.92; ann$y <- lim[2] * 0.88

  p <- ggplot2::ggplot(dat, ggplot2::aes(estimate_L, estimate_R)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = nv_lw("reference_pt"),
                         colour = "grey70") +
    ggplot2::geom_point(colour = "#1F3D52", size = 0.12, alpha = 0.25, shape = 16) +
    ggplot2::geom_text(data = ann, ggplot2::aes(x, y, label = t), hjust = 0,
                       size = nv_size(5.4), family = nv_palette()$typography$family,
                       inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ lab, nrow = 1) +
    ggplot2::coord_cartesian(xlim = lim, ylim = lim, expand = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(-3, 0, 3)) +
    ggplot2::scale_y_continuous(breaks = c(-3, 0, 3)) +
    ggplot2::labs(x = "left hemisphere (log2)", y = "right hemisphere (log2)") +
    nv_theme()

  write_csv_safe(dat[, c("dataset", "contrast", "ProteinGroupID", "estimate_L",
                         "estimate_R", "sign_agreement", "pearson_r",
                         "sign_agreement_fraction", "n_evaluable_proteins")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# Precision gain, small enough to sit as an inset beside c.
nvp_precision_inset <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  d <- d[is.finite(d$ICC_single_side) & is.finite(d$ICC_bilateral_mean), , drop = FALSE]
  long <- rbind(
    data.frame(id = seq_len(nrow(d)), sampling = "1 side", ICC = d$ICC_single_side),
    data.frame(id = seq_len(nrow(d)), sampling = "L+R", ICC = d$ICC_bilateral_mean))
  long$sampling <- factor(long$sampling, levels = c("1 side", "L+R"))

  p <- ggplot2::ggplot(long, ggplot2::aes(sampling, ICC)) +
    ggplot2::geom_line(ggplot2::aes(group = id), colour = "grey85",
                       linewidth = nv_lw("reference_pt")) +
    ggplot2::stat_summary(fun = stats::median, geom = "crossbar", width = 0.5,
                          linewidth = nv_lw("data_pt"), colour = "#1F3D52") +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.05, 1.05)) +
    ggplot2::labs(x = NULL, y = "ICC") +
    nv_theme()

  write_csv_safe(d[, c("dataset", "endpoint_class", "endpoint_id",
                       "ICC_single_side", "ICC_bilateral_mean")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# d. Compartment identity: the simplest direct demonstration of enrichment.
nvp_compartment <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  d$marker_class <- as.character(d$marker_class)
  d$intended_dataset <- as.character(d$intended_dataset)
  keep <- d$intended_dataset %in% names(nv_dataset_colours())
  d <- d[keep, , drop = FALSE]
  d$lab <- nv_dataset_label(d$intended_dataset)

  p <- ggplot2::ggplot(d, ggplot2::aes(stats::reorder(marker_class,
                                                      intended_minus_comparator_log2),
                                       intended_minus_comparator_log2)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey70") +
    ggplot2::geom_boxplot(ggplot2::aes(fill = intended_dataset), width = 0.6,
                          linewidth = nv_lw("axis_pt"), outlier.size = 0.2,
                          outlier.colour = "grey60") +
    ggplot2::scale_fill_manual(values = nv_dataset_colours(),
                               labels = nv_dataset_label, name = NULL) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "intended − comparator (log2)") +
    nv_theme(grid = "x") +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.height = ggplot2::unit(2, "mm"))

  write_csv_safe(d[, c("marker_class", "intended_dataset", "ProteinGroupID",
                       "intended_minus_comparator_log2")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# e. External spatial validation (Kaulich), simplified.
nvp_external_validation <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  nes <- intersect(c("NES"), names(d))[1]
  fdr <- intersect(c("p_adjust", "p.adjust", "FDR"), names(d))[1]
  sig <- intersect(c("external_signature", "signature", "ID"), names(d))[1]
  ctr <- intersect(c("internal_contrast", "contrast"), names(d))[1]
  if (any(is.na(c(nes, fdr, sig, ctr)))) stop("kaulich source lacks expected columns")
  z <- data.frame(sig = as.character(d[[sig]]), contrast = as.character(d[[ctr]]),
                  NES = suppressWarnings(as.numeric(d[[nes]])),
                  FDR = suppressWarnings(as.numeric(d[[fdr]])),
                  dataset = if ("dataset" %in% names(d)) as.character(d$dataset) else NA_character_,
                  stringsAsFactors = FALSE)
  z <- z[is.finite(z$NES), , drop = FALSE]
  z$sig <- gsub("_", " ", z$sig)
  z$supported <- is.finite(z$FDR) & z$FDR < 0.05
  lim <- max(abs(z$NES), na.rm = TRUE)

  p <- ggplot2::ggplot(z, ggplot2::aes(contrast, sig, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_point(data = z[z$supported, , drop = FALSE], size = 0.45,
                        colour = "black") +
    nv_diverging(limits = c(-lim, lim), name = "NES") +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "right")

  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# f. Internal anatomical program validation.
nvp_internal_validation <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  nes <- intersect(c("NES"), names(d))[1]
  fdr <- intersect(c("p_adjust", "p.adjust"), names(d))[1]
  if (any(is.na(c(nes, fdr)))) stop("internal GO source lacks NES/FDR")
  z <- data.frame(term = as.character(d$Description),
                  contrast = as.character(d$contrast),
                  NES = suppressWarnings(as.numeric(d[[nes]])),
                  FDR = suppressWarnings(as.numeric(d[[fdr]])), stringsAsFactors = FALSE)
  z <- z[is.finite(z$NES) & is.finite(z$FDR), , drop = FALSE]
  # keep the strongest supported terms only; the rest belongs in a table
  z$supported <- z$FDR < 0.05
  top <- unique(z$term[z$supported][order(z$FDR[z$supported])])
  top <- utils::head(top, 8)
  z <- z[z$term %in% top, , drop = FALSE]
  z$term <- factor(z$term, levels = rev(top))
  lim <- max(abs(z$NES), na.rm = TRUE)

  p <- ggplot2::ggplot(z, ggplot2::aes(contrast, term, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_point(data = z[z$supported, , drop = FALSE], size = 0.4,
                        colour = "black") +
    nv_diverging(limits = c(-lim, lim), name = "NES") +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = nv_pt("minimum_pt")),
                   legend.position = "right")

  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ============================================================ FIGURE 3

nvp_unit_levels <- function() {
  c("CA1_slm", "CA1_so", "CA1_sr", "CA2_slm", "CA2_so", "CA2_sr",
    "CA3_so", "CA3_sr", "DG_mo", "DG_po",
    "CA1_sp", "CA2_sp", "CA3_sp", "DG_sg", "CA1", "CA2", "CA3", "DG")
}
nvp_short_unit <- function(u) sub("^(CA[123]|DG)_?", "", u)

# a. Compact differential-proteome landscape across all three compartments.
# Deliberately NOT led by a giant CA2-SLM bar: counts are on a square-root
# axis and robustness status is encoded compactly.
nvp_da_landscape <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  # built from CANONICAL sources so this layer does not depend on the Part-16
  # candidate layer: the DAP membership table gives the unit of every
  # FDR-supported protein, the atlas gives its audited claimability
  member <- nv_read_csv(repo_path(panel$primary_source))
  atlas <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]))
  m <- merge(member[, c("dataset", "spatial_unit", "ProteinGroupID")],
             atlas[, c("ProteinGroupID", "QC_claimability")],
             by = "ProteinGroupID", all.x = TRUE)
  m$QC_claimability[is.na(m$QC_claimability) | m$QC_claimability == ""] <- "not_audited"
  d <- as.data.frame(table(dataset = m$dataset, spatial_unit = m$spatial_unit,
                           status = m$QC_claimability), stringsAsFactors = FALSE)
  names(d)[names(d) == "Freq"] <- "n_proteins"
  d <- d[d$n_proteins > 0, , drop = FALSE]
  d$lab <- nv_dataset_label(d$dataset)
  d$lab <- factor(d$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))
  d$unit <- factor(d$spatial_unit, levels = nvp_unit_levels())
  d$status <- factor(d$status, levels = c("claimable", "claimable_with_caveat",
                                          "not_claimable", "not_evaluable", "not_audited"))
  # A sqrt axis distorts stacked segments (ggplot stacks before transforming),
  # which made the pale not_evaluable block swallow the bar. Linear axis in a
  # short wide strip is honest and keeps the atlas below as the centrepiece.
  p <- ggplot2::ggplot(d, ggplot2::aes(n_proteins, unit, fill = status)) +
    ggplot2::geom_col(width = 0.66, linewidth = 0) +
    ggplot2::scale_fill_manual(values = nv_claim_colours(), name = NULL, drop = FALSE,
                               labels = function(x) gsub("_", " ", x)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 30),
                                breaks = c(0, 10, 20, 28)) +
    ggplot2::facet_grid(lab ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "FDR-supported proteins", y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(legend.position = "right",
                   strip.text.y = ggplot2::element_text(angle = 0, face = "bold",
                                                        hjust = 0),
                   legend.key.height = ggplot2::unit(1.9, "mm"))
  write_csv_safe(d, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# Shared theme-cell builder (descriptive summary of canonical GSEA terms).
nvp_theme_cells <- function(panel, fdr_cut = 0.05) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  d <- d[d$assignment_status %in% c("single_theme", "multi_theme", "qc_review") &
           nzchar(as.character(d$theme_id)), , drop = FALSE]
  key <- paste(d$dataset, d$spatial_unit, d$contrast, d$theme_id, sep = "\r")
  sp <- split(seq_len(nrow(d)), key)
  out <- do.call(rbind, lapply(sp, function(ix) {
    z <- d[ix, , drop = FALSE]
    data.frame(dataset = z$dataset[1], spatial_unit = z$spatial_unit[1],
               contrast = z$contrast[1], theme_id = z$theme_id[1],
               manuscript_theme = z$manuscript_theme[1], theme_role = z$theme_role[1],
               n_terms = nrow(z),
               n_terms_FDR_supported = sum(is.finite(z$GSEA_FDR) & z$GSEA_FDR < fdr_cut),
               median_NES = stats::median(z$NES, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  out$has_FDR_support <- out$n_terms_FDR_supported > 0L
  out$summary_basis <- paste0(
    "median NES of the theme's constituent canonical GO-BP terms; FDR support ",
    "counted from those same terms at FDR < ", fdr_cut,
    "; theme aggregation is an interpretation layer, NOT a new FDR family")
  rownames(out) <- NULL
  out
}

# b. Cross-dataset ranked-GSEA atlas. The visual centrepiece.
nvp_gsea_atlas <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  cells <- nvp_theme_cells(panel)
  z <- cells[cells$contrast == "SUS - RES", , drop = FALSE]
  z$lab <- nv_dataset_label(z$dataset)
  z$lab <- factor(z$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))
  z$unit <- factor(z$spatial_unit, levels = nvp_unit_levels())
  # qc_review themes are QC context and are pushed to the bottom, greyed label
  # Registry display labels are long; the axis is the wrong place for a full
  # ontology phrase. Short forms here, full names in the source data.
  short <- c(
    "Mitochondrial respiration / OXPHOS" = "Mitochondrial / OXPHOS",
    "Synaptic signaling / vesicle-mediated transport" = "Synaptic / vesicle",
    "RNA processing / splicing / RNP organization" = "RNA processing / RNP",
    "Translation / ribosome biogenesis" = "Translation / ribosome",
    "Chromatin organization / remodeling" = "Chromatin organization",
    "Autophagy / endolysosomal trafficking" = "Autophagy / lysosome",
    "Intermediate filament organization" = "Intermediate filament",
    "Epidermal / keratinocyte differentiation" = "Epidermal / keratinocyte")
  nm <- unname(short[z$manuscript_theme])
  z$theme_lab <- ifelse(is.na(nm), z$manuscript_theme, nm)
  z$theme_lab <- ifelse(z$theme_role == "qc_review",
                        paste0(z$theme_lab, " †"), z$theme_lab)
  ord <- unique(z[order(z$theme_role != "primary", z$theme_lab), "theme_lab"])
  z$theme_lab <- factor(z$theme_lab, levels = rev(ord))
  lim <- max(abs(z$median_NES), na.rm = TRUE)

  p <- ggplot2::ggplot(z, ggplot2::aes(unit, theme_lab, fill = median_NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_point(data = z[z$has_FDR_support, , drop = FALSE],
                        size = 0.5, colour = "black") +
    nv_diverging(limits = c(-lim, lim), name = "NES") +
    ggplot2::facet_grid(. ~ lab, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 5.4),
                   legend.position = "right",
                   strip.text = ggplot2::element_text(face = "bold"))

  write_csv_safe(cells, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# c-e. Representative ranked evidence, one per compartment.
#
# A direct representation of the ranked evidence rather than another summary
# heatmap: every measured gene is placed at its stored rank statistic and the
# term's stored leading-edge genes are marked. Both quantities are READ from
# the canonical clusterProfiler audit files. No enrichment is recomputed and no
# running-enrichment score is derived here.
nvp_ranked_example <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  base <- repo_path("data", "processed", "04_differential_expression_enrichment",
                    "clusterProfiler", as.character(panel$example_dataset),
                    "phenotype_within_unit", as.character(panel$example_unit_dir),
                    as.character(panel$example_contrast_dir))
  audits <- file.path(base, "protein_group_audits")
  ranked <- nv_read_csv_longpath(audits, "collapsed_gene_input.csv",
                                 stringsAsFactors = FALSE)
  res <- nv_read_csv_longpath(file.path(base, "GO", "BP"), "GSEA_BP_results_full.csv",
                              stringsAsFactors = FALSE)
  term <- as.character(panel$example_term_id)
  r <- res[res$ID == term, , drop = FALSE]
  if (!nrow(r)) stop("example term not present in the stored GSEA result: ", term)

  prov <- as.data.frame(nv_read_csv_longpath(
    audits, "gsea_go_term_gene_provenance.csv",
    reader = function(f, ...) data.table::fread(f, showProgress = FALSE)))
  le <- prov[prov$term_id == term, , drop = FALSE]

  stat_col <- intersect(c("collapsed_statistic", "log2fc"), names(ranked))[1]
  rk <- data.frame(gene = as.character(ranked$official_gene_symbol),
                   stat = suppressWarnings(as.numeric(ranked[[stat_col]])),
                   stringsAsFactors = FALSE)
  rk <- rk[is.finite(rk$stat), , drop = FALSE]
  rk <- rk[order(-rk$stat), , drop = FALSE]
  rk$rank <- seq_len(nrow(rk))
  rk$leading_edge <- rk$gene %in% as.character(le$official_gene_symbol)

  dcol <- nv_dataset_colours()[[as.character(panel$example_dataset)]]
  fam <- nv_palette()$typography$family
  n <- nrow(rk)

  top <- ggplot2::ggplot(rk, ggplot2::aes(rank, stat)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey75") +
    ggplot2::geom_area(fill = "grey88") +
    ggplot2::geom_segment(data = rk[rk$leading_edge, , drop = FALSE],
                          ggplot2::aes(x = rank, xend = rank,
                                       y = min(rk$stat), yend = min(rk$stat) * 0.72),
                          colour = dcol, linewidth = nv_lw("reference_pt")) +
    ggplot2::annotate("text", x = n * 0.97, y = max(rk$stat) * 0.92, hjust = 1,
                      family = fam, size = nv_size(5.2),
                      label = sprintf("NES %.2f\nFDR %.0e\n%d leading-edge genes",
                                      r$NES[1], r$p.adjust[1], sum(rk$leading_edge))) +
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                breaks = c(1, round(n / 2), n)) +
    ggplot2::labs(x = "rank in SUS − RES ordered list",
                  y = "moderated t") +
    nv_theme()

  out <- rk
  out$term_id <- term
  out$term_description <- r$Description[1]
  out$NES <- r$NES[1]; out$FDR <- r$p.adjust[1]
  out$dataset <- as.character(panel$example_dataset)
  out$spatial_unit <- as.character(panel$example_unit_dir)
  out$evidence_note <- paste0(
    "gene positions and leading-edge membership are READ from the canonical ",
    "clusterProfiler audit; no enrichment statistic is recomputed here")
  write_csv_safe(out, csv_path)
  nv_save_panel(top, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# f. Stress effect versus baseline spatial affinity, redesigned.
nvp_stress_identity <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  atlas <- nv_read_csv(repo_path(panel$primary_source))
  a <- atlas[atlas$is_sus_res_fdr_supported %in% TRUE, , drop = FALSE]
  a$claim <- ifelse(is.na(a$QC_claimability) | a$QC_claimability == "",
                    "not_audited", a$QC_claimability)
  a$dsl <- nv_dataset_label(a$dataset)
  a$dsl <- factor(a$dsl,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))
  # rank is expressed as a fraction of that protein's own unit count, so
  # 10-unit neuropil and 4-unit soma/microglia are directly comparable
  a$rel <- (a$baseline_rank_of_strongest_effect_unit - 1) /
    pmax(a$n_spatial_units - 1, 1)
  a$claim <- factor(a$claim, levels = c("claimable", "claimable_with_caveat",
                                        "not_claimable", "not_evaluable", "not_audited"))

  # DETERMINISTIC vertical spread. geom_jitter draws from the RNG, so the same
  # inputs would produce a different SVG on every run; overlapping points are
  # fanned out by their stable rank instead.
  a <- a[order(a$dsl, a$rel, a$ProteinGroupID), , drop = FALSE]
  a$offset <- unlist(lapply(split(seq_len(nrow(a)), a$dsl), function(ix) {
    k <- seq_along(ix)
    c(0, 0.2, -0.2, 0.34, -0.34)[((k - 1L) %% 5L) + 1L]
  }), use.names = FALSE)
  a$ypos <- as.integer(a$dsl) + a$offset

  p <- ggplot2::ggplot(a, ggplot2::aes(rel, ypos)) +
    ggplot2::annotate("rect", xmin = -0.03, xmax = 1 / 3, ymin = -Inf, ymax = Inf,
                      fill = "grey96") +
    ggplot2::annotate("segment", x = 1 / 3, xend = 1 / 3, y = -Inf, yend = Inf,
                      linewidth = nv_lw("reference_pt"), colour = "grey70",
                      linetype = "22") +
    ggplot2::geom_point(ggplot2::aes(colour = claim), size = 0.8, alpha = 0.95) +
    ggplot2::scale_colour_manual(values = nv_claim_colours(), name = NULL, drop = FALSE,
                                 labels = function(x) gsub("_", " ", x)) +
    ggplot2::scale_x_continuous(breaks = c(0, 1 / 3, 1),
                                labels = c("own peak\nunit", "high-affinity\nboundary",
                                           "lowest-affinity\nunit"),
                                limits = c(-0.03, 1.03), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(levels(a$dsl)),
                                labels = levels(a$dsl),
                                limits = c(0.4, length(levels(a$dsl)) + 0.6)) +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme() +
    ggplot2::theme(legend.position = "right",
                   axis.text.x = ggplot2::element_text(size = nv_pt("minimum_pt")),
                   legend.key.height = ggplot2::unit(2, "mm"))

  keep <- c("dataset", "ProteinGroupID", "GeneSymbol",
            "baseline_rank_of_strongest_effect_unit", "n_spatial_units",
            "effect_identity_relationship", "QC_claimability")
  o <- a[, intersect(keep, names(a)), drop = FALSE]
  o$relative_baseline_rank <- a$rel
  o$interpretation_note <- paste0(
    "position of the strongest stress effect within the protein's OWN baseline ",
    "spatial ranking, scaled to its number of units. Grey band is the ",
    "prespecified high-affinity third. This is where the effect sits, not ",
    "movement of protein between compartments")
  write_csv_safe(o, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# g. Compact cross-method synthesis: strongest recurring programs only.
nvp_synthesis <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  cells <- nvp_theme_cells(panel)
  s <- cells[cells$contrast == "SUS - RES" & cells$theme_role == "primary", , drop = FALSE]
  agg <- stats::aggregate(has_FDR_support ~ manuscript_theme + dataset, s, sum)
  names(agg)[3] <- "n_supported_contexts"
  agg <- agg[agg$n_supported_contexts > 0, , drop = FALSE]
  tot <- stats::aggregate(n_supported_contexts ~ manuscript_theme, agg, sum)
  keep <- tot$manuscript_theme[order(-tot$n_supported_contexts)]
  agg <- agg[agg$manuscript_theme %in% keep, , drop = FALSE]
  agg$manuscript_theme <- factor(agg$manuscript_theme, levels = rev(keep))
  agg$lab <- nv_dataset_label(agg$dataset)
  agg$lab <- factor(agg$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))

  p <- ggplot2::ggplot(agg, ggplot2::aes(n_supported_contexts, manuscript_theme,
                                         fill = lab)) +
    ggplot2::geom_col(width = 0.62, linewidth = 0) +
    ggplot2::scale_fill_manual(values = stats::setNames(
      unname(nv_dataset_colours()[c("neuron_neuropil", "neuron_soma", "microglia")]),
      nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia"))), name = NULL) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "spatial contexts with FDR-supported terms", y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(size = 5.2),
                   legend.key.height = ggplot2::unit(2, "mm"))

  agg$evidence_role <- paste0(
    "contextual convergence across compartments; the columns are NOT five ",
    "statistically independent validations")
  write_csv_safe(agg, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# Small WGCNA supporting strip, for the variant that keeps WGCNA in the main.
nvp_wgcna_small <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  e <- nv_read_csv(repo_path(panel$primary_source))
  lab <- intersect(c("module_display_label", "ModuleLabel_Final"), names(e))[1]
  z <- data.frame(module = if (is.na(lab)) e$module_id else e[[lab]],
                  module_id = e$module_id, contrast = e$contrast,
                  estimate = e$estimate, fdr = e$tier_specific_fdr,
                  stringsAsFactors = FALSE)
  z$module <- factor(z$module, levels = rev(unique(z$module[order(z$module_id)])))
  lim <- max(abs(z$estimate), na.rm = TRUE)
  n_sig <- sum(z$fdr < 0.05, na.rm = TRUE)

  p <- ggplot2::ggplot(z, ggplot2::aes(contrast, module, fill = estimate)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    nv_diverging(limits = c(-lim, lim), name = "effect") +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 5),
                   legend.position = "right")

  z$fdr_supported_cells_in_panel <- n_sig
  z$status_note <- paste0(
    "descriptive only: ", n_sig, " of ", nrow(z),
    " module x contrast cells reach tier-specific FDR < 0.05")
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ================================================ WGCNA EXTENDED DATA

nvp_ed_module_identity <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  a <- nv_read_csv(repo_path(panel$primary_source))
  z <- a[, intersect(c("dataset", "ModuleID", "canonical_display_label",
                       "module_size", "external_celltype_all",
                       "bilateral_reproducibility_class", "spatial_tau",
                       "peak_unit"), names(a)), drop = FALSE]
  z$lab <- nv_dataset_label(z$dataset)
  z$lab <- factor(z$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))
  z$reproducible <- z$bilateral_reproducibility_class == "reproducible_level_and_pattern"
  z$mod <- paste0(sub("WGCNA_", "", z$ModuleID))

  p <- ggplot2::ggplot(z, ggplot2::aes(spatial_tau, stats::reorder(mod, spatial_tau))) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = spatial_tau,
                                       yend = stats::reorder(mod, spatial_tau)),
                          colour = "grey85", linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(colour = reproducible), size = 1) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = unname(nv_evidence_colours()["supported"]),
                 "FALSE" = unname(nv_evidence_colours()["descriptive"])),
      labels = c("TRUE" = "bilaterally reproducible", "FALSE" = "poor reproducibility"),
      name = NULL) +
    ggplot2::facet_grid(lab ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(x = "spatial specificity (tau)", y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(size = 5),
                   strip.text.y = ggplot2::element_text(angle = 0, face = "bold"))

  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

nvp_ed_wgcna_heatmap <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  nvp_wgcna_small(panel, svg_path, csv_path, w_mm, h_mm)
}

nvp_ed_celltype <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  a <- nv_read_csv(repo_path(panel$primary_source))
  col <- intersect(c("external_celltype_all", "external_celltype_core_kME06"), names(a))[1]
  if (is.na(col)) stop("no external cell-type column in the atlas")
  z <- a[, c("dataset", "ModuleID", col), drop = FALSE]
  names(z)[3] <- "celltype"
  z$celltype[is.na(z$celltype) | z$celltype == ""] <- "none"
  z$lab <- nv_dataset_label(z$dataset)
  tab <- as.data.frame(table(lab = z$lab, celltype = z$celltype), stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  tab$lab <- factor(tab$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))

  p <- ggplot2::ggplot(tab, ggplot2::aes(lab, stats::reorder(celltype, Freq), fill = Freq)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_text(ggplot2::aes(label = Freq), size = nv_size(5),
                       family = nv_palette()$typography$family) +
    ggplot2::scale_fill_gradient(low = "#EFEFEC", high = "#3D5A73", guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 5))
  write_csv_safe(tab, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
