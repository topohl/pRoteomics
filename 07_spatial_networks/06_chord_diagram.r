# =========================================================
# Proteomics overlaps: plots + Excel shared-protein exports
# =========================================================

# Packages
pkgs <- c("dplyr","tidyr","circlize","RColorBrewer","tools","magick","grid","openxlsx","digest")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# Paths
in_dir  <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/neuron-phenotypeWithinUnit"
out_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/prot_sign"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------- Read CSVs and build wide tables --------
csv_files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(csv_files) > 0)

read_both <- function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("gene_symbol","padj","log2fc")
  if (!all(need %in% names(df))) {
    stop(sprintf("File %s missing columns: %s", basename(f),
                 paste(setdiff(need, names(df)), collapse = ", ")))
  }
  base <- tools::file_path_sans_ext(basename(f))
  list(
    padj = setNames(df[, c("gene_symbol","padj")],   c("gene_symbol", base)),
    fc   = setNames(df[, c("gene_symbol","log2fc")], c("gene_symbol", base))
  )
}

lst <- lapply(csv_files, read_both)
wide_padj <- Reduce(function(x, y) dplyr::full_join(x, y, by = "gene_symbol"), lapply(lst, `[[`, "padj"))
wide_fc   <- Reduce(function(x, y) dplyr::full_join(x, y, by = "gene_symbol"), lapply(lst, `[[`, "fc"))

# -------- Build significance mask --------
alpha <- 0.05
sig <- wide_padj %>% mutate(across(-gene_symbol, ~ .x < alpha))
comp_cols <- setdiff(names(sig), "gene_symbol")

# -------- Pair overlaps by direction (inclusive) --------
pair_edges_dir <- function(dir_sign = +1, lfc_min = 0) {
  pairs <- combn(comp_cols, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    a <- p[1]; b <- p[2]
    in_both <- sig[[a]] & sig[[b]]
    if (!any(in_both, na.rm = TRUE)) return(tibble(from=a,to=b,weight=0))
    fa <- wide_fc[[a]][in_both]
    fb <- wide_fc[[b]][in_both]
    ok <- (sign(fa) == dir_sign) & (sign(fb) == dir_sign) &
          (abs(fa) >= lfc_min) & (abs(fb) >= lfc_min)
    tibble(from = a, to = b, weight = sum(ok, na.rm = TRUE))
  })
  bind_rows(rows) %>% filter(weight > 0)
}

edges_up   <- pair_edges_dir(+1)   # both significant and up
edges_down <- pair_edges_dir(-1)   # both significant and down

# -------- Plotting helpers --------
plot_chord <- function(edge_df, title_text) {
  if (nrow(edge_df) == 0) {
    grid::grid.newpage()
    grid::grid.text(paste0("No edges: ", title_text))
    return(invisible(NULL))
  }
  nodes <- sort(unique(c(edge_df$from, edge_df$to)))
  pal <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(nodes))
  grid.col <- setNames(pal, nodes)

  circlize::circos.clear()
  circlize::circos.par(start.degree = 90,
                       gap.after = c(rep(4, length(nodes)-1), 12),
                       track.margin = c(0.01, 0.01),
                       points.overflow.warning = FALSE)

  circlize::chordDiagram(edge_df,
                         grid.col = grid.col,
                         transparency = 0.45,
                         link.sort = TRUE,
                         link.largest.ontop = TRUE,
                         link.border = NA,
                         directional = 0,
                         annotationTrack = "grid",
                         preAllocateTracks = 1)

  circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    si <- circlize::get.cell.meta.data("sector.index")
    circlize::circos.text(circlize::CELL_META$xcenter,
                          circlize::CELL_META$ylim[1] + circlize::mm_y(3),
                          si, facing = "clockwise", niceFacing = TRUE,
                          adj = c(0, 0.5), cex = 0.9, col = "#222")
  }, bg.border = NA)

  grid::grid.text(title_text, x = grid::unit(0.5, "npc"), y = grid::unit(0.98, "npc"),
                  gp = grid::gpar(fontsize = 12, col = "#222"))
}

# -------- Save plots --------
f_up   <- file.path(out_dir, "chord_upregulated.png")
f_down <- file.path(out_dir, "chord_downregulated.png")
f_side <- file.path(out_dir, "chord_up_down_side_by_side.png")

# produce both SVG (requested) and PNG (kept for downstream magick usage)
f_up_svg   <- sub("\\.png$", ".svg", f_up)
f_down_svg <- sub("\\.png$", ".svg", f_down)

# sizes: convert px/res -> inches (svg() expects inches)
inch_w <- 2200 / 250
inch_h <- 2200 / 250

# SVG outputs
svg(filename = f_up_svg, width = inch_w, height = inch_h)
plot_chord(edges_up, "Shared significant upregulated proteins")
dev.off()

svg(filename = f_down_svg, width = inch_w, height = inch_h)
plot_chord(edges_down, "Shared significant downregulated proteins")
dev.off()

# Also keep PNGs so subsequent magick::image_read(...) (which expects PNG in original code) works
png(f_up, width = 2200, height = 2200, res = 250)
plot_chord(edges_up, "Shared significant upregulated proteins")
dev.off()

png(f_down, width = 2200, height = 2200, res = 250)
plot_chord(edges_down, "Shared significant downregulated proteins")
dev.off()

# Combine side-by-side
img1 <- magick::image_read(f_up)
img2 <- magick::image_read(f_down)
combo <- magick::image_append(c(img1, img2))
magick::image_write(combo, f_side)

# -------- Excel exports: shared proteins per pair --------
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# Build long summary for a given direction
summary_by_dir <- function(dir_sign = +1, lfc_min = 0) {
  # Direction-specific significant mask
  is_sig <- sig
  for (cc in comp_cols) {
    vv <- wide_fc[[cc]]
    is_sig[[cc]] <- is_sig[[cc]] & (sign(vv) == dir_sign) & (abs(vv) >= lfc_min)
  }

  # For each pair, collect shared genes and return long rows
  pairs <- combn(comp_cols, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    a <- p[1]; b <- p[2]
    in_both <- is_sig[[a]] & is_sig[[b]]
    if (!any(in_both, na.rm = TRUE)) return(NULL)
    idx <- which(in_both)
    tibble::tibble(
      gene_symbol = wide_padj$gene_symbol[idx],
      comparison_a = a,
      padj_a = wide_padj[[a]][idx],
      log2fc_a = wide_fc[[a]][idx],
      comparison_b = b,
      padj_b = wide_padj[[b]][idx],
      log2fc_b = wide_fc[[b]][idx]
    )
  })
  dplyr::bind_rows(rows)
}

# Create a single-sheet workbook per direction
write_summary_wb <- function(df_long, label = "Up") {
  if (is.null(df_long) || nrow(df_long) == 0) {
    df_long <- tibble::tibble(
      gene_symbol = character(), comparison_a = character(), padj_a = numeric(), log2fc_a = numeric(),
      comparison_b = character(), padj_b = numeric(), log2fc_b = numeric()
    )
  }

  # Add per-protein stats: how many pairs it is shared in and the list of pairs
  pair_str <- paste(df_long$comparison_a, df_long$comparison_b, sep = " vs ")
  df_with_pair <- dplyr::mutate(df_long, pair = pair_str)

  per_gene <- df_with_pair %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      pairs = paste(sort(unique(pair)), collapse = " | "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_pairs), gene_symbol)

  # Final summary sheet: one row per gene with ranking, plus optional exemplar stats
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "summary")

  # Left part: per-gene ranking
  openxlsx::writeData(wb, "summary", per_gene, startRow = 1, startCol = 1)

  # Right part: detailed long table (optional)
  openxlsx::addWorksheet(wb, "details")
  # Order details by gene then by decreasing evidence
  details <- df_with_pair %>%
    dplyr::mutate(max_padj = pmax(padj_a, padj_b, na.rm = TRUE),
                  mean_abs_lfc = rowMeans(cbind(abs(log2fc_a), abs(log2fc_b)), na.rm = TRUE)) %>%
    dplyr::arrange(gene_symbol, max_padj, dplyr::desc(mean_abs_lfc)) %>%
    dplyr::select(gene_symbol, comparison_a, padj_a, log2fc_a, comparison_b, padj_b, log2fc_b)
  openxlsx::writeData(wb, "details", details, startRow = 1, startCol = 1)

  out_file <- file.path(out_dir, paste0("shared_proteins_", tolower(label), "_summary.xlsx"))
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  message("Saved ", label, " summary workbook: ", out_file)
  invisible(out_file)
}

# Build and save Up and Down summaries
up_long   <- summary_by_dir(+1, lfc_min = 0)
down_long <- summary_by_dir(-1, lfc_min = 0)

up_xlsx_summary   <- write_summary_wb(up_long,   label = "Up")
down_xlsx_summary <- write_summary_wb(down_long, label = "Down")

message("All outputs saved to: ", out_dir)