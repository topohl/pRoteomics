# Lightweight enrichment plotting helpers.

save_plot_dual <- function(plot, path_svg, width = 7, height = 5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NA_character_))
  dir.create(dirname(path_svg), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path_svg, plot, width = width, height = height, units = "in", limitsize = FALSE)
  pdf_path <- sub("\\.svg$", ".pdf", path_svg)
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in", limitsize = FALSE)
  invisible(path_svg)
}
