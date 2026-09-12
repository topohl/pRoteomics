# Part-21 WGCNA spatial panels and small locators.
#
# DOWNSTREAM ONLY. WGCNA is never re-run: module membership is read from the
# canonical shortlist, and the spatial profile is a descriptive average of
# already-canonical CON baseline values.

# =====================================================================
# ED4: the module x spatial-unit fingerprint
# =====================================================================
#
# ED4 previously reported only a module's PEAK unit, a specificity number and an
# annotation. That states a spatial fact without showing spatial structure. This
# panel shows the actual profile: for each canonical module, the mean
# within-protein CON z of its member proteins in every spatial unit of its own
# compartment.
#
# Membership comes from the canonical WGCNA shortlist (ModuleID), read as a
# lookup. Values come from the phenotype-blind CON baseline. Nothing is refitted
# and no phenotype effect enters.
s6_module_fingerprint_table <- function(min_members = 5L) {
  base <- nv_read_csv(path_results("tables", "manuscript_candidates",
                                   "spatial_v6",
                                   "spatial_v6_con_baseline_profile_long.csv"))
  mem <- nv_read_csv(repo_path("results", "tables", "10_biological_integration",
                               "wgcna_candidate_protein_shortlist", "global",
                               "wgcna_candidate_proteins_all.csv"))
  mem <- unique(mem[, c("dataset", "ProteinGroupID", "ModuleID")])
  mem <- mem[!is.na(mem$ModuleID) & nzchar(mem$ModuleID), , drop = FALSE]

  j <- merge(base[, c("dataset", "ProteinGroupID", "spatial_unit", "con_z")],
             mem, by = c("dataset", "ProteinGroupID"))
  key <- paste(j$dataset, j$ModuleID, j$spatial_unit, sep = "\r")
  out <- do.call(rbind, lapply(split(seq_len(nrow(j)), key), function(ix) {
    z <- j[ix, , drop = FALSE]
    data.frame(dataset = z$dataset[1], ModuleID = z$ModuleID[1],
               spatial_unit = z$spatial_unit[1],
               n_members = length(unique(z$ProteinGroupID)),
               mean_con_z = mean(z$con_z, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out$mean_con_z[out$n_members < min_members] <- NA_real_

  # annotation tracks from the canonical module table
  wa <- nv_read_csv(repo_path("results", "tables", "11_spatial_systems", "atlas",
                              "WGCNA_module_spatial_cell_affinity.csv"))
  keep <- c("dataset", "ModuleID", "canonical_display_label",
            "bilateral_reproducibility_class", "external_celltype_all",
            "external_FDR_all", "peak_unit", "spatial_tau", "module_size")
  keep <- intersect(keep, names(wa))
  out <- merge(out, unique(wa[, keep]), by = c("dataset", "ModuleID"),
               all.x = TRUE)
  out$value_definition <- paste0(
    "mean within-protein CON z across the module's member proteins in that ",
    "spatial unit; descriptive, WGCNA not recomputed, no phenotype effect")
  out
}

s6_ed_module_fingerprint <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  z <- s6_module_fingerprint_table()
  z <- z[!is.na(z$mean_con_z), , drop = FALSE]
  z$sg_unit <- sg_resolve_unit(z$spatial_unit, z$dataset)

  # A module exists only in its own compartment, so a single 18-column grid
  # would be about two thirds structurally blank - the same defect Part 20
  # removed from Figure 2f. One block per compartment, each with only its own
  # units, so no structurally inapplicable cell is ever drawn.
  z$comp <- factor(sg_compartment_label(z$dataset),
                   levels = sg_compartments()$short)
  z$unit_f <- factor(z$sg_unit, levels = sg_units()$analysis_key)
  z$mod_f <- factor(z$ModuleID, levels = rev(sort(unique(z$ModuleID))))
  disp <- stats::setNames(sg_units()$display, sg_units()$analysis_key)
  lyr <- stats::setNames(
    ifelse(sg_units()$layer_is_resolution & !is.na(sg_units()$layer),
           sg_layer_display(sg_units()$layer), sg_units()$region),
    sg_units()$analysis_key)
  lim <- max(abs(z$mean_con_z), na.rm = TRUE) * c(-1, 1)

  p <- ggplot2::ggplot(z, ggplot2::aes(unit_f, mod_f)) +
    ggplot2::geom_tile(ggplot2::aes(fill = mean_con_z), colour = "white",
                       linewidth = 0.1) +
    nv_diverging(limits = lim, name = "mean\nCON z") +
    ggplot2::scale_x_discrete(labels = function(x) unname(disp[as.character(x)]),
                              drop = TRUE) +
    ggplot2::scale_y_discrete(drop = TRUE) +
    ggplot2::facet_wrap(~comp, ncol = 1, scales = "free",
                        strip.position = "left") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    "CON animals only. Each cell is the mean within-protein CON ",
                    "z of that module's member proteins in that spatial unit. ",
                    "A module exists only in\nits own compartment, so each ",
                    "block carries only that compartment's units and no ",
                    "structurally inapplicable cell is drawn. Descriptive: ",
                    "WGCNA is not\nrecomputed and no phenotype effect enters ",
                    "this panel.")) +
    nv_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 4.8, colour = "grey25"),
      axis.text.y = ggplot2::element_text(size = 4.6),
      strip.text.y.left = ggplot2::element_text(size = 5.4, face = "bold",
                                                angle = 90, colour = "grey15"),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      panel.spacing.y = ggplot2::unit(1.6, "mm"),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.7, "mm"),
      legend.key.height = ggplot2::unit(4, "mm"),
      plot.caption = ggplot2::element_text(size = 4.6, colour = "grey35",
                                           hjust = 0, lineheight = 1.15))
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# =====================================================================
# ED5: WGCNA phenotype field with spatial row annotation
# =====================================================================
#
# Keeps the canonical module x contrast effect field and adds the spatial
# context tracks the brief asks for: module spatial peak, bilateral
# reproducibility class and external cell-type affinity. The colouring must not
# imply inferential support, so the null is stated on the panel and no cell is
# starred.
s6_ed_wgcna_phenotype <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  z <- nv_read_csv(repo_path(panel$primary_source))
  wa <- nv_read_csv(repo_path("results", "tables", "11_spatial_systems", "atlas",
                              "WGCNA_module_spatial_cell_affinity.csv"))

  # This source carries BOTH an all-empty `ModuleID` column and the populated
  # `module_id`. Taking the first name that merely EXISTS silently selects the
  # empty one and collapses 15 modules to a single NA row, so a candidate only
  # counts if it actually holds values.
  pick <- function(cands, what) {
    for (nm in cands) {
      if (!nm %in% names(z)) next
      v <- z[[nm]]
      if (any(!is.na(v) & nzchar(as.character(v)))) return(nm)
    }
    stop("s6_ed_wgcna_phenotype: no populated ", what, " column among ",
         paste(cands, collapse = ", "), call. = FALSE)
  }
  mcol <- pick(c("module_id", "ModuleID", "module", "endpoint_id"), "module")
  ccol <- pick(c("contrast", "Contrast"), "contrast")
  vcol <- pick(c("effect", "estimate", "logFC", "value", "median_effect"),
               "effect")
  z$mid <- as.character(z[[mcol]])
  z$con <- as.character(z[[ccol]])
  z$val <- as.numeric(z[[vcol]])
  z <- z[!is.na(z$val), , drop = FALSE]

  # The effect source spells modules "WGCNA_m01" while the atlas spells them
  # "m01"; and both tables carry a `dataset` column, so a naive merge would
  # silently rename it to dataset.x/dataset.y and break every downstream
  # reference. Normalise the key and prefix the annotation columns.
  norm_mid <- function(x) sub("^WGCNA_", "", as.character(x))
  z$mid_key <- norm_mid(z$mid)
  ann <- data.frame(
    mid_key = norm_mid(wa$ModuleID),
    ann_dataset = as.character(wa$dataset),
    ann_peak_unit = as.character(wa$peak_unit),
    ann_bilateral = as.character(wa$bilateral_reproducibility_class),
    ann_celltype = as.character(wa$external_celltype_all),
    stringsAsFactors = FALSE)
  # the effect field is neuropil-level, so resolve the annotation within that
  # compartment rather than taking whichever module happens to sort first
  ds <- unique(as.character(z$dataset))
  if (length(ds) == 1L && ds %in% ann$ann_dataset) {
    ann <- ann[ann$ann_dataset == ds, , drop = FALSE]
  }
  ann <- ann[!duplicated(ann$mid_key), , drop = FALSE]
  z <- merge(z, ann, by = "mid_key", all.x = TRUE, sort = FALSE)
  z$peak_label <- ""
  ok <- !is.na(z$ann_peak_unit) & nzchar(z$ann_peak_unit) & !is.na(z$ann_dataset)
  if (any(ok)) {
    z$peak_label[ok] <- sg_unit_label(z$ann_peak_unit[ok], z$ann_dataset[ok])
  }
  z$bilateral_reproducibility_class <- z$ann_bilateral
  z$external_celltype_all <- z$ann_celltype

  mods <- unique(z[, c("mid", "peak_label", "bilateral_reproducibility_class",
                       "external_celltype_all")])
  mods <- mods[order(mods$mid), , drop = FALSE]
  mods$row_label <- ifelse(nzchar(mods$peak_label),
                           paste0(mods$mid, "   peak ", mods$peak_label),
                           mods$mid)
  z$ypos <- match(z$mid, rev(mods$mid))
  cl <- c("RES - CON", "SUS - CON", "SUS - RES")
  cl <- cl[cl %in% z$con]
  if (!length(cl)) cl <- sort(unique(z$con))
  z <- z[z$con %in% cl, , drop = FALSE]
  z$xpos <- match(z$con, cl)

  ny <- nrow(mods)
  lim <- max(abs(z$val), na.rm = TRUE) * c(-1, 1)
  p <- ggplot2::ggplot(z, ggplot2::aes(xpos, ypos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = val), colour = "white",
                       linewidth = 0.14) +
    nv_diverging(limits = lim, name = "effect") +
    # the scale limit must leave room for the annotation tracks: a scale limit
    # narrower than the annotation x positions silently DROPS them, which is
    # not the same as clipping them
    ggplot2::scale_x_continuous(breaks = seq_along(cl), labels = cl,
                                expand = c(0, 0),
                                limits = c(0.5, length(cl) + 4.2)) +
    ggplot2::scale_y_continuous(breaks = seq_len(ny), labels = rev(mods$row_label),
                                expand = c(0, 0), limits = c(0.5, ny + 0.5)) +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0(
                    "DESCRIPTIVE ONLY. No module x contrast cell reaches ",
                    "tier-specific FDR support anywhere in this panel, so the ",
                    "colour scale must not be read\nas evidence of a group ",
                    "difference. Row labels carry each module's baseline ",
                    "spatial peak. WGCNA is not recomputed.")) +
    nv_theme_tile() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 5),
      axis.text.y = ggplot2::element_text(size = 4.8),
      legend.position = "right",
      legend.key.width = ggplot2::unit(1.8, "mm"),
      legend.key.height = ggplot2::unit(4, "mm"),
      plot.caption = ggplot2::element_text(size = 4.6, colour = "grey35",
                                           hjust = 0, lineheight = 1.15))

  # annotation tracks to the right of the field
  tr <- unique(z[, c("mid", "ypos", "bilateral_reproducibility_class",
                     "external_celltype_all")])
  tr$x1 <- length(cl) + 0.9
  tr$x2 <- length(cl) + 1.9
  p <- p +
    ggplot2::geom_text(data = tr,
                       ggplot2::aes(x = x1, y = ypos,
                                    label = ifelse(is.na(bilateral_reproducibility_class),
                                                   "", bilateral_reproducibility_class)),
                       inherit.aes = FALSE, hjust = 0, family = fam,
                       size = nv_size(4.3), colour = "grey35") +
    ggplot2::geom_text(data = tr,
                       ggplot2::aes(x = x2, y = ypos,
                                    label = ifelse(is.na(external_celltype_all),
                                                   "", external_celltype_all)),
                       inherit.aes = FALSE, hjust = 0, family = fam,
                       size = nv_size(4.3), colour = "grey35") +
    ggplot2::annotate("text", x = c(length(cl) + 0.9, length(cl) + 1.9),
                      y = ny + 0.35, hjust = 0, family = fam,
                      size = nv_size(4.5), fontface = "bold", colour = "grey20",
                      label = c("bilateral class", "external cell type")) +
    ggplot2::coord_cartesian(xlim = c(0.5, length(cl) + 3.6),
                             ylim = c(0.5, ny + 0.8), clip = "off")

  z$inferential_status <- "0 of 45 module x contrast cells FDR-supported"
  write_csv_safe(z, csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}

# =====================================================================
# ED3: a tiny CA2-SLM locator
# =====================================================================
#
# So the reader immediately knows the QC problem sits in neuropil, region CA2,
# layer SLM - rather than having to decode the string "CA2_slm".
s6_ed_ca2_locator <- function(panel, svg_path, csv_path, w_mm, h_mm) {
  fam <- nv_palette()$typography$family
  accent <- unname(nv_dataset_colours()[["neuron_neuropil"]])
  ex <- data.frame(arc_frac = 0.45, band = "slm", stringsAsFactors = FALSE)

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
  seg <- function(d, frac, half = 0.11) {
    n <- nrow(d)
    d[max(1L, floor((frac - half) * n)):min(n, ceiling((frac + half) * n)), ,
      drop = FALSE]
  }
  at <- function(frac, r) {
    t <- t0 + frac * (t1 - t0)
    c(x = r * cos(t), y = r * sin(t))
  }
  a1 <- at(0.08, 1.42); a2 <- at(0.45, 1.46); a3 <- at(0.86, 1.42)
  rl <- data.frame(x = c(a1[["x"]], a2[["x"]], a3[["x"]], 0.46),
                   y = c(a1[["y"]], a2[["y"]], a3[["y"]], -0.28),
                   l = c("CA1", "CA2", "CA3", "DG"),
                   hit = c(FALSE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = sr, ggplot2::aes(x, y), colour = "grey80",
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = so, ggplot2::aes(x, y), colour = "grey80",
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = slm, ggplot2::aes(x, y), colour = "grey80",
                       linewidth = 0.3, linetype = "22") +
    ggplot2::geom_path(data = ca, ggplot2::aes(x, y), colour = "grey80",
                       linewidth = 1.0) +
    ggplot2::geom_path(data = dgu, ggplot2::aes(x, y), colour = "grey80",
                       linewidth = 0.8) +
    ggplot2::geom_path(data = dgl, ggplot2::aes(x, y), colour = "grey80",
                       linewidth = 0.8) +
    ggplot2::geom_path(data = seg(slm, 0.45), ggplot2::aes(x, y),
                       colour = accent, linewidth = 1.7, lineend = "round") +
    ggplot2::geom_text(data = rl,
                       ggplot2::aes(x, y, label = l, colour = hit,
                                    fontface = ifelse(hit, 2, 1)),
                       family = fam, size = nv_size(4.8), show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c("TRUE" = accent, "FALSE" = "grey62"),
                                 guide = "none") +
    ggplot2::annotate("text", x = 0, y = -1.32, label = "neuropil",
                      family = fam, size = nv_size(5.0), colour = "grey30") +
    ggplot2::annotate("text", x = 0, y = -1.62, label = "CA2  SLM",
                      family = fam, size = nv_size(6.0), fontface = "bold",
                      colour = accent) +
    ggplot2::coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.85, 1.62),
                         expand = FALSE) +
    ggplot2::theme_void(base_family = fam)

  write_csv_safe(data.frame(
    compartment = "neuron_neuropil", region = "CA2", layer = "slm",
    analysis_key = "CA2_slm",
    note = "locator only; carries no quantity", stringsAsFactors = FALSE),
    csv_path)
  nv_save_panel(p, svg_path, w_mm, h_mm)
  invisible(list(status = "ok"))
}
