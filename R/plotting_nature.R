# Compact plotting helpers for manuscript-scale enrichment figures.

mm_to_in <- function(mm) {
  as.numeric(mm) / 25.4
}

theme_nature_base <- function(base_size = 7, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      line = ggplot2::element_line(linewidth = 0.25),
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "black", margin = ggplot2::margin(b = 2)),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
      legend.key.size = grid::unit(3.0, "mm"),
      legend.spacing.y = grid::unit(0.5, "mm"),
      legend.box.spacing = grid::unit(1.0, "mm"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(hjust = 0),
      panel.spacing = grid::unit(1.0, "mm")
    )
}

theme_nature_heatmap <- function(base_size = 7, base_family = "Arial") {
  theme_nature_base(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right"
    )
}

theme_nature_dotplot <- function(base_size = 7, base_family = "Arial") {
  theme_nature_base(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = "grey90", linewidth = 0.2),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right"
    )
}

save_nature_svg <- function(plot, filename, width_mm, height_mm) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  device <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = mm_to_in(width_mm),
    height = mm_to_in(height_mm),
    units = "in",
    device = device,
    limitsize = FALSE
  )
  if (!requireNamespace("svglite", quietly = TRUE)) {
    pdf_file <- sub("\\.svg$", ".pdf", filename)
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      width = mm_to_in(width_mm),
      height = mm_to_in(height_mm),
      units = "in",
      device = grDevices::pdf,
      limitsize = FALSE
    )
  }
  invisible(filename)
}

clean_program_label <- function(x) {
  x <- as.character(x)
  recode <- c(
    RNA_RNP_processing = "RNA/RNP processing",
    Ribosome_Translation = "Ribosome/translation",
    Translation_Ribosome = "Ribosome/translation",
    Mitochondria_OXPHOS_Metabolism = "Mitochondria/OXPHOS/metabolism",
    Mitochondria_OXPHOS = "Mitochondria/OXPHOS/metabolism",
    Proteostasis_Ubiquitin_Folding = "Proteostasis/ubiquitin/protein folding",
    Proteostasis_Lysosome = "Proteostasis/ubiquitin/protein folding",
    Synapse_Vesicle_Organization = "Synapse/vesicle organization",
    Synapse_Plasticity = "Synapse/vesicle organization",
    Cytoskeleton_Motility = "Cytoskeleton/motility",
    Cytoskeleton_Transport = "Cytoskeleton/motility",
    Development_Patterning = "Development/patterning",
    Immune_Microglia = "Immune/microglia",
    Other = "Other"
  )
  out <- unname(recode[x])
  out[is.na(out)] <- gsub("_", " ", x[is.na(out)])
  out
}

clean_spatial_unit_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- gsub("\\bso\\b", "SO", x, ignore.case = TRUE)
  x <- gsub("\\bsr\\b", "SR", x, ignore.case = TRUE)
  x <- gsub("\\bslm\\b", "SLM", x, ignore.case = TRUE)
  x <- gsub("\\bsp\\b", "SP", x, ignore.case = TRUE)
  x <- gsub("\\bmo\\b", "MO", x, ignore.case = TRUE)
  x <- gsub("\\bpo\\b", "PO", x, ignore.case = TRUE)
  x <- gsub("\\bsg\\b", "SG", x, ignore.case = TRUE)
  x <- gsub("\\bgranule\\b", "granule", x, ignore.case = TRUE)
  x
}

clean_comparison_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_vs_", " vs ", x, fixed = TRUE)
  x <- gsub("_", " ", x, fixed = TRUE)
  x
}

anatomical_spatial_unit_levels <- function(units) {
  units <- unique(as.character(units))
  preferred <- c(
    as.vector(outer(c("CA1", "CA2", "CA3"), c("so", "sr", "slm", "sp"), paste, sep = "_")),
    "DG_mo", "DG_po", "DG_sg", "DG_granule",
    "CA1", "CA2", "CA3", "DG"
  )
  c(preferred[preferred %in% units], sort(setdiff(units, preferred)))
}
