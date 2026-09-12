# =====================================================================
# Story-v3 panel renderers.
#
# A fourth candidate family. It REUSES the Part-17 typography, palette and
# exact-box machinery by sourcing R/nature_v2_figure_utils.R read-only, and
# reuses the Part-17 renderers that the audit judged already correct. What is
# new here are the panels the audit said were missing or wrongly encoded:
#
#   sv2_pca          restored from the canonical layer, de-cluttered
#   sv2_precision    dumbbell instead of a 3x3 line-plot grid
#   sv3_anatomy      the missing BRIDGE panel: where the programs actually are
#   sv3_proteins     the old protein-heatmap DESIGN with QC-clean biology
#   sv3_wgcna_circle / sv3_wgcna_strip   two encodings of module STRUCTURE
#
# Downstream only: no model fit, no enrichment, no FDR. The Part-17 guard
# scans this file too.
# =====================================================================

# ---------------------------------------------------------- FIGURE 2

# c. Global molecular structure. Restored because it answers a question no
# other panel answers: do the measured proteomes occupy structured molecular
# spaces at all? De-cluttered - no per-sample labels, compartment colour,
# region shape, variance explained on the axes.
svp_pca <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  ve <- nv_read_csv(repo_path(as.character(unlist(panel$input_dependencies))[1]),
                    required = FALSE)
  pc <- function(k) {
    if (is.null(ve)) return(sprintf("PC%d", k))
    v <- ve$variance_explained[ve$PC == paste0("PC", k)]
    if (!length(v)) return(sprintf("PC%d", k))
    sprintf("PC%d (%.0f%%)", k, 100 * v[1])
  }
  d <- d[d$dataset %in% names(nv_dataset_colours()), , drop = FALSE]
  d$lab <- nv_dataset_label(d$dataset)
  d$lab <- factor(d$lab,
    levels = nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia")))
  d$region <- factor(d$region, levels = c("CA1", "CA2", "CA3", "DG"))

  p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_point(ggplot2::aes(colour = lab, shape = region),
                        size = 0.6, stroke = 0.25, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      unname(nv_dataset_colours()[c("neuron_neuropil", "neuron_soma", "microglia")]),
      nv_dataset_label(c("neuron_neuropil", "neuron_soma", "microglia"))),
      name = NULL) +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 3), name = NULL) +
    ggplot2::labs(x = pc(1), y = pc(2)) +
    nv_theme() +
    ggplot2::theme(legend.position = "right",
                   legend.spacing.y = ggplot2::unit(0.4, "mm"),
                   legend.key.height = ggplot2::unit(2.2, "mm"))

  write_csv_safe(d[, c("Sample", "PC1", "PC2", "dataset", "region", "layer",
                       "celltype_layer")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# e. Precision as a DUMBBELL. The audit found the 3x3 line grid abstract: the
# reader had to infer improvement from nine small panels. One row per endpoint
# class, faint per-endpoint pairs behind a bold class summary, makes the single
# point immediate: averaging hemispheres raises reliability.
svp_precision_dumbbell <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  d <- nv_read_csv(repo_path(panel$primary_source))
  d <- d[is.finite(d$ICC_single_side) & is.finite(d$ICC_bilateral_mean), , drop = FALSE]
  lab <- c(wgcna_module_eigengene = "WGCNA modules",
           reference_marker_score = "reference markers",
           empirical_compartment_score = "compartment scores")
  d$cls <- unname(lab[d$endpoint_class])
  d$cls <- factor(d$cls, levels = rev(unname(lab)))
  d$ymid <- as.numeric(d$cls)
  # deterministic vertical spread of the faint individual pairs
  d <- d[order(d$cls, d$endpoint_id, d$dataset), , drop = FALSE]
  d$yoff <- unlist(lapply(split(seq_len(nrow(d)), d$cls), function(ix) {
    n <- length(ix); if (n == 1L) 0 else seq(-0.3, 0.3, length.out = n)
  }), use.names = FALSE)
  d$ypos <- d$ymid + d$yoff

  sm <- do.call(rbind, lapply(split(d, d$cls), function(z) data.frame(
    cls = z$cls[1], ymid = z$ymid[1],
    single = stats::median(z$ICC_single_side),
    bilateral = stats::median(z$ICC_bilateral_mean), stringsAsFactors = FALSE)))

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = d,
      ggplot2::aes(x = ICC_single_side, xend = ICC_bilateral_mean,
                   y = ypos, yend = ypos),
      colour = "grey88", linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_segment(data = sm,
      ggplot2::aes(x = single, xend = bilateral, y = ymid, yend = ymid),
      colour = "grey45", linewidth = nv_lw("data_pt") * 1.6,
      arrow = grid::arrow(length = ggplot2::unit(1.1, "mm"), type = "closed")) +
    ggplot2::geom_point(data = sm, ggplot2::aes(single, ymid),
                        colour = unname(nv_evidence_colours()["descriptive"]),
                        size = 1.5) +
    ggplot2::geom_point(data = sm, ggplot2::aes(bilateral, ymid),
                        colour = unname(nv_evidence_colours()["supported"]),
                        size = 1.5) +
    ggplot2::scale_y_continuous(breaks = sm$ymid, labels = as.character(sm$cls),
                                limits = c(0.4, nrow(sm) + 0.6)) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                                expand = c(0, 0)) +
    ggplot2::labs(x = "ICC:  one hemisphere → bilateral mean", y = NULL) +
    nv_theme(grid = "x")

  write_csv_safe(d[, c("dataset", "endpoint_class", "endpoint_id",
                       "ICC_single_side", "ICC_bilateral_mean")], csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ---------------------------------------------------------- FIGURE 3

# The three representative programs, defined ONCE so every downstream panel
# refers to the same biology with the same label and the same accent colour.
svp_programs <- function() {
  data.frame(
    key = c("synaptic", "rna", "oxphos"),
    dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
    unit = c("CA3_sr", "CA2_sp", "CA1"),
    unit_dir = c("CA3_sr", "CA2_sp", "CA1_microglia"),
    contrast_dir = c("CA3srsus_CA3srres", "CA2spsus_CA2spres",
                     "CA1microgliasus_CA1microgliares"),
    term = c("GO:0099536", "GO:0006397", "GO:0006119"),
    label = c("synaptic signalling", "mRNA processing",
              "oxidative phosphorylation"),
    accent = c("#3D5A73", "#7A9BB0", "#C2A878"),
    stringsAsFactors = FALSE)
}

# c. THE BRIDGE PANEL. Panel b says which programs differ; this says WHERE.
# It reuses the Figure-2 sampling map as the spatial substrate - the reader has
# already learned that grid - and marks only the three prespecified programs
# with their direction. No anatomical geometry is invented.
svp_anatomy_bridge <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  meta <- as.data.frame(readxl::read_excel(
    repo_path("data", "metadata", "TPE9_sample_metadata_males.xlsx")))
  meta <- meta[!(meta$exclude %in% TRUE), , drop = FALSE]
  meta$region <- toupper(as.character(meta$region))
  g <- unique(meta[, c("region", "layer", "celltype_layer")])
  g <- g[g$region %in% c("CA1", "CA2", "CA3", "DG"), , drop = FALSE]
  g$row <- ifelse(g$celltype_layer == "microglia", "microglia",
                  ifelse(g$celltype_layer == "neuron_soma",
                         paste0("soma ", g$layer), paste0("neuropil ", g$layer)))
  row_order <- c(paste0("neuropil ", c("slm", "sr", "so", "mo", "po")),
                 paste0("soma ", c("sp", "sg")), "microglia")
  lv <- rev(intersect(row_order, unique(g$row)))
  g$row <- factor(g$row, levels = lv)
  g$region <- factor(g$region, levels = c("CA1", "CA2", "CA3", "DG"))

  pr <- svp_programs()
  # place each program on its own sampled cell
  pr$region <- factor(sub("_.*$", "", pr$unit), levels = c("CA1", "CA2", "CA3", "DG"))
  pr$row <- factor(c("neuropil sr", "soma sp", "microglia"), levels = lv)

  # direction of the program from the canonical theme table
  cells <- nv_read_csv(repo_path(panel$primary_source))
  pr$NES <- vapply(seq_len(nrow(pr)), function(i) {
    z <- cells[cells$dataset == pr$dataset[i] & cells$spatial_unit == pr$unit[i] &
                 cells$contrast == "SUS - RES" &
                 cells$GO_ID == pr$term[i], , drop = FALSE]
    if (!nrow(z)) NA_real_ else z$NES[1]
  }, numeric(1))
  pr$dir <- ifelse(is.finite(pr$NES) & pr$NES > 0, "higher in SUS", "higher in RES")
  pr$num <- as.character(seq_len(nrow(pr)))
  # numbered markers on the map with a key underneath: callout labels placed on
  # a 4-column grid collide with the panel edges at this width
  key_txt <- paste(sprintf("%s  %s (%s), %s", pr$num, pr$label,
                           nv_dataset_label(pr$dataset), pr$dir), collapse = "\n")

  p <- ggplot2::ggplot(g, ggplot2::aes(region, row)) +
    ggplot2::geom_tile(fill = "grey94", colour = "white", linewidth = 0.4,
                       width = 0.9, height = 0.85) +
    ggplot2::geom_tile(data = pr, ggplot2::aes(fill = key), colour = "white",
                       linewidth = 0.5, width = 0.9, height = 0.85,
                       show.legend = FALSE) +
    ggplot2::geom_text(data = pr, ggplot2::aes(label = num), size = nv_size(6),
                       family = fam, fontface = "bold", colour = "white") +
    ggplot2::scale_fill_manual(values = stats::setNames(pr$accent, pr$key)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL, caption = key_txt) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", size = 6),
                   axis.text.y = ggplot2::element_text(size = 5.2),
                   plot.caption = ggplot2::element_text(size = 5, hjust = 0,
                                                        colour = "grey20",
                                                        lineheight = 1.25))

  out <- pr[, c("key", "label", "dataset", "unit", "term", "NES", "dir")]
  out$substrate <- "Figure 2 sampling map; no anatomical geometry is invented"
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# Shared loader for one program's stored ranked evidence.
svp_ranked <- function(prog) {
  base <- repo_path("data", "processed", "04_differential_expression_enrichment",
                    "clusterProfiler", prog$dataset, "phenotype_within_unit",
                    prog$unit_dir, prog$contrast_dir)
  audits <- file.path(base, "protein_group_audits")
  ranked <- nv_read_csv_longpath(audits, "collapsed_gene_input.csv",
                                 stringsAsFactors = FALSE)
  res <- nv_read_csv_longpath(file.path(base, "GO", "BP"),
                              "GSEA_BP_results_full.csv", stringsAsFactors = FALSE)
  prov <- as.data.frame(nv_read_csv_longpath(
    audits, "gsea_go_term_gene_provenance.csv",
    reader = function(f, ...) data.table::fread(f, showProgress = FALSE)))
  r <- res[res$ID == prog$term, , drop = FALSE]
  if (!nrow(r)) stop("term absent from stored GSEA result: ", prog$term)
  le <- prov[prov$term_id == prog$term, , drop = FALSE]
  stat_col <- intersect(c("collapsed_statistic", "log2fc"), names(ranked))[1]
  rk <- data.frame(gene = as.character(ranked$official_gene_symbol),
                   stat = suppressWarnings(as.numeric(ranked[[stat_col]])),
                   stringsAsFactors = FALSE)
  rk <- rk[is.finite(rk$stat), , drop = FALSE]
  rk <- rk[order(-rk$stat), , drop = FALSE]
  rk$rank <- seq_len(nrow(rk))
  rk$leading_edge <- rk$gene %in% as.character(le$official_gene_symbol)
  list(rank = rk, NES = r$NES[1], FDR = r$p.adjust[1],
       description = r$Description[1],
       le_stat = stats::setNames(le$rank_statistic, le$official_gene_symbol))
}

# d/e/f. Direct ranked evidence, matched grammar across the three compartments.
# The audit asked for the leading-edge marks to be unmissable and the numbers
# out of the data area.
svp_rank_example <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  pr <- svp_programs()
  prog <- pr[pr$key == as.character(panel$program_key), , drop = FALSE]
  if (!nrow(prog)) stop("unknown program key: ", panel$program_key)
  ev <- svp_ranked(prog)
  rk <- ev$rank
  n <- nrow(rk)
  fam <- nv_palette()$typography$family
  acc <- prog$accent[1]
  ymin <- min(rk$stat); ymax <- max(rk$stat)
  band <- ymin - (ymax - ymin) * 0.16

  p <- ggplot2::ggplot(rk, ggplot2::aes(rank, stat)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = nv_lw("reference_pt"),
                        colour = "grey75") +
    ggplot2::geom_area(fill = "grey90") +
    ggplot2::geom_segment(data = rk[rk$leading_edge, , drop = FALSE],
                          ggplot2::aes(x = rank, xend = rank,
                                       y = band, yend = ymin * 0.97),
                          colour = acc, linewidth = 0.22) +
    ggplot2::annotate("text", x = n * 0.5, y = ymax * 0.96, hjust = 0.5,
                      family = fam, size = nv_size(5.4), colour = acc,
                      fontface = "bold", label = prog$label[1]) +
    ggplot2::annotate("text", x = n * 0.5, y = ymax * 0.62, hjust = 0.5,
                      family = fam, size = nv_size(5),
                      label = sprintf("NES %.2f   FDR %.0e", ev$NES, ev$FDR)) +
    ggplot2::annotate("text", x = n * 0.5, y = band * 1.06, hjust = 0.5,
                      family = fam, size = nv_size(5), colour = "grey35",
                      label = sprintf("%d leading-edge proteins", sum(rk$leading_edge))) +
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                breaks = c(1, n), labels = c("SUS", "RES")) +
    ggplot2::coord_cartesian(ylim = c(band * 1.12, ymax * 1.04), expand = FALSE) +
    ggplot2::labs(x = NULL, y = "moderated t") +
    nv_theme()

  out <- rk
  out$program <- prog$label[1]; out$term_id <- prog$term[1]
  out$term_description <- ev$description
  out$NES <- ev$NES; out$FDR <- ev$FDR
  out$dataset <- prog$dataset[1]; out$spatial_unit <- prog$unit[1]
  out$evidence_note <- paste0(
    "gene positions and leading-edge membership READ from the canonical ",
    "clusterProfiler audit; no enrichment statistic is recomputed")
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# g. PROTEIN-LEVEL ZOOM. The old Figure-3 heatmap design, with the QC-clean
# biology the audit asked for: the leading-edge proteins of the three
# supported programs, shown across their own compartment's spatial units.
#
# SELECTION RULE, applied identically to all three programs: take that
# program's stored leading-edge proteins, rank by |stored rank statistic|,
# keep the top N. No gene is chosen by name.
svp_protein_zoom <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  top_n <- as.integer(panel$top_n %||% 7L)
  pr <- svp_programs()
  units_of <- list(
    neuron_neuropil = c("CA1_slm", "CA1_so", "CA1_sr", "CA2_slm", "CA2_so",
                        "CA2_sr", "CA3_so", "CA3_sr", "DG_mo", "DG_po"),
    neuron_soma = c("CA1_sp", "CA2_sp", "CA3_sp", "DG_sg"),
    microglia = c("CA1", "CA2", "CA3", "DG"))

  rows <- list()
  for (i in seq_len(nrow(pr))) {
    prog <- pr[i, , drop = FALSE]
    ev <- svp_ranked(prog)
    le <- ev$le_stat
    le <- le[order(-abs(le))]
    keep <- names(le)[seq_len(min(top_n, length(le)))]

    for (u in units_of[[prog$dataset[1]]]) {
      # microglia contrast files carry the compartment token in the unit name
      # (CA1microgliasus_...), unlike neuropil and soma
      tok <- gsub("_", "", u)
      if (identical(prog$dataset[1], "microglia")) tok <- paste0(tok, "microglia")
      f <- repo_path("data", "processed", "02_id_mapping", "mapped",
                     prog$dataset[1], "forward", "per_file",
                     sprintf("%ssus_%sres.csv", tok, tok))
      if (!file.exists(f)) next
      da <- nv_read_csv(f)
      sym <- intersect(c("official_gene_symbol", "gene_symbol"), names(da))[1]
      m <- match(keep, da[[sym]])
      rows[[length(rows) + 1L]] <- data.frame(
        program = prog$label[1], program_key = prog$key[1],
        dataset = prog$dataset[1], gene = keep, spatial_unit = u,
        log2FC = da$log2fc[m], BH_FDR = da$padj[m],
        rank_statistic = unname(le[keep]),
        is_example_unit = u == prog$unit[1], stringsAsFactors = FALSE)
    }
  }
  z <- dplyr::bind_rows(rows)
  z <- z[!is.na(z$log2FC), , drop = FALSE]
  z$program <- factor(z$program, levels = pr$label)
  z$unit <- factor(z$spatial_unit,
                   levels = unlist(units_of[c("neuron_neuropil", "neuron_soma",
                                              "microglia")], use.names = FALSE))
  # genes ordered within program by their stored rank statistic
  ord <- unique(z[order(z$program, -abs(z$rank_statistic)), c("program", "gene")])
  z$gene <- factor(z$gene, levels = rev(unique(ord$gene)))
  lim <- stats::quantile(abs(z$log2FC), 0.97, na.rm = TRUE)

  # Three program-specific heatmaps side by side (OPTION 1). A single shared
  # x-axis of all 18 units would leave each program's block empty across the
  # other two compartments - more than a third of the panel as whitespace.
  z$unit <- droplevels(z$unit)
  p <- ggplot2::ggplot(z, ggplot2::aes(unit, gene, fill = log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = nv_lw("tile_border_pt")) +
    ggplot2::geom_point(data = z[z$is_example_unit, , drop = FALSE],
                        size = 0.35, colour = "black") +
    nv_diverging(limits = c(-lim, lim), name = "log2FC\nSUS − RES",
                 oob = scales::squish) +
    # facet_wrap frees BOTH axes per facet; facet_grid(. ~ program) shares the
    # y axis, so every program's genes would be listed in every facet
    ggplot2::facet_wrap(~ program, scales = "free", nrow = 1) +
    ggplot2::labs(x = NULL, y = NULL) +
    nv_theme_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 5),
                   strip.text = ggplot2::element_text(size = 5, face = "bold"),
                   legend.position = "right")

  z$selection_rule <- sprintf(paste0(
    "stored leading-edge proteins of each program's example term, ranked by ",
    "|stored rank statistic|, top %d per program; no gene chosen by name"), top_n)
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# h. WGCNA STRUCTURE, encoding 1: simplified circle. The old circular atlas was
# visually the most "systems-level" panel in the library; its weakness was the
# phenotype interpretation, not the geometry. NO phenotype significance is
# encoded here - only module identity, spatial preference, cell-type affinity
# and paired-hemisphere reproducibility.
svp_wgcna_circle <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  a <- nv_read_csv(repo_path(panel$primary_source))
  z <- a[a$dataset == "neuron_neuropil", , drop = FALSE]
  z <- z[order(z$ModuleID), , drop = FALSE]
  n <- nrow(z)
  z$idx <- seq_len(n)
  z$ang <- (z$idx - 0.5) / n * 2 * pi
  z$mod <- sub("WGCNA_", "", z$ModuleID)
  z$repro <- z$bilateral_reproducibility_class == "reproducible_level_and_pattern"
  ext <- intersect(c("external_celltype_all", "external_celltype_core_kME06"), names(z))[1]
  z$cell <- if (is.na(ext)) NA_character_ else z[[ext]]
  z$cell[is.na(z$cell) | z$cell == ""] <- "none"
  z$peak <- toupper(sub("_.*$", "", as.character(z$peak_unit)))

  ring <- function(r0, r1, val, lab) {
    do.call(rbind, lapply(seq_len(n), function(i) {
      a0 <- (i - 1) / n * 2 * pi; a1 <- i / n * 2 * pi
      t <- seq(a0, a1, length.out = 8)
      data.frame(x = c(r0 * sin(t), rev(r1 * sin(t))),
                 y = c(r0 * cos(t), rev(r1 * cos(t))),
                 grp = paste0(lab, i), v = val[i], ring = lab,
                 stringsAsFactors = FALSE)
    }))
  }
  r_tau <- ring(0.62, 0.78, z$spatial_tau, "tau")
  r_peak <- ring(0.80, 0.90, z$peak, "peak")
  r_rep <- ring(0.92, 1.00, ifelse(z$repro, "yes", "no"), "repro")

  fam <- nv_palette()$typography$family
  labs <- data.frame(x = 1.13 * sin(z$ang), y = 1.13 * cos(z$ang),
                     l = z$mod, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = r_tau,
      ggplot2::aes(x, y, group = grp, fill = as.numeric(v)),
      colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient(low = "#EFEFEC", high = "#1F3D52",
                                 name = "spatial\nspecificity", limits = c(0, 1)) +
    ggnewscale_or_identity() +
    ggplot2::geom_polygon(data = r_peak,
      ggplot2::aes(x, y, group = grp), fill = "grey93",
      colour = "white", linewidth = 0.15) +
    ggplot2::geom_text(data = data.frame(
        x = 0.85 * sin(z$ang), y = 0.85 * cos(z$ang), l = z$peak),
      ggplot2::aes(x, y, label = l), size = nv_size(5), family = fam) +
    ggplot2::geom_polygon(data = r_rep[r_rep$v == "yes", , drop = FALSE],
      ggplot2::aes(x, y, group = grp), fill = unname(nv_evidence_colours()["supported"]),
      colour = "white", linewidth = 0.15) +
    ggplot2::geom_polygon(data = r_rep[r_rep$v == "no", , drop = FALSE],
      ggplot2::aes(x, y, group = grp), fill = unname(nv_evidence_colours()["descriptive"]),
      colour = "white", linewidth = 0.15) +
    ggplot2::geom_text(data = labs, ggplot2::aes(x, y, label = l),
                       size = nv_size(5), family = fam) +
    ggplot2::annotate("text", x = 0, y = 0.12, label = "WGCNA\nmodule structure",
                      family = fam, size = nv_size(5.4), fontface = "bold",
                      lineheight = 1.05) +
    ggplot2::annotate("text", x = 0, y = -0.16,
                      label = "outer: bilateral reproducibility\nmiddle: spatial peak\ninner: spatial specificity",
                      family = fam, size = nv_size(5), colour = "grey35",
                      lineheight = 1.15) +
    ggplot2::coord_equal(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), expand = FALSE) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(legend.position = "right",
                   legend.title = ggplot2::element_text(size = nv_pt("legend_title_pt")),
                   legend.text = ggplot2::element_text(size = nv_pt("legend_text_pt")),
                   legend.key.size = ggplot2::unit(2.4, "mm"),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))

  out <- z[, c("dataset", "ModuleID", "canonical_display_label", "spatial_tau",
               "peak_unit", "bilateral_reproducibility_class")]
  out$external_celltype <- z$cell
  out$encodes_phenotype_significance <- FALSE
  out$note <- "module STRUCTURE only; no phenotype effect or significance is encoded"
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# ggnewscale is optional; without it the later rings simply use fixed fills.
ggnewscale_or_identity <- function() {
  if (requireNamespace("ggnewscale", quietly = TRUE)) ggnewscale::new_scale_fill()
  else ggplot2::theme()
}

# h. WGCNA STRUCTURE, encoding 2: linear strips carrying exactly the same
# information, for a like-for-like legibility comparison against the circle.
svp_wgcna_strip <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  a <- nv_read_csv(repo_path(panel$primary_source))
  z <- a[a$dataset == "neuron_neuropil", , drop = FALSE]
  ext <- intersect(c("external_celltype_all", "external_celltype_core_kME06"), names(z))[1]
  z$cell <- if (is.na(ext)) NA_character_ else z[[ext]]
  z$cell[is.na(z$cell) | z$cell == ""] <- "none"
  z$mod <- sub("WGCNA_", "", z$ModuleID)
  z$repro <- z$bilateral_reproducibility_class == "reproducible_level_and_pattern"
  z$peak <- toupper(sub("_.*$", "", as.character(z$peak_unit)))
  z <- z[order(-z$spatial_tau), , drop = FALSE]
  z$mod <- factor(z$mod, levels = rev(z$mod))

  fam <- nv_palette()$typography$family
  p <- ggplot2::ggplot(z, ggplot2::aes(spatial_tau, mod)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = spatial_tau, yend = mod),
                          colour = "grey88", linewidth = nv_lw("reference_pt")) +
    ggplot2::geom_point(ggplot2::aes(colour = repro), size = 1.2) +
    ggplot2::geom_text(ggplot2::aes(x = 1.02, label = peak), hjust = 0,
                       size = nv_size(5), family = fam) +
    ggplot2::geom_text(ggplot2::aes(x = 1.22, label = substr(cell, 1, 14)),
                       hjust = 0, size = nv_size(5), family = fam,
                       colour = "grey30") +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = unname(nv_evidence_colours()["supported"]),
                 "FALSE" = unname(nv_evidence_colours()["descriptive"])),
      labels = c("TRUE" = "bilaterally reproducible", "FALSE" = "poor"),
      name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1.9), breaks = c(0, 0.5, 1),
                                expand = c(0, 0)) +
    ggplot2::labs(x = "spatial specificity (tau)", y = NULL) +
    nv_theme(grid = "x") +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(size = 5),
                   legend.key.height = ggplot2::unit(2, "mm"))

  out <- z[, c("dataset", "ModuleID", "canonical_display_label", "spatial_tau",
               "peak_unit", "bilateral_reproducibility_class")]
  out$external_celltype <- z$cell
  out$encodes_phenotype_significance <- FALSE
  write_csv_safe(out, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
