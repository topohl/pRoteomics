# ================================================================
# Script: analysis/spatial_networks/render_network_chord_diagram.R
# Stage: networks
# Scope: dataset_specific
# Consumes: required results/preprocessing/map_protein_identifiers/<dataset>/tables/mapped/forward/per_file/*.csv; data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv; optional none declared in pipeline.yml
# Produces: results/spatial_networks/render_network_chord_diagram/<dataset>/plots; results/spatial_networks/render_network_chord_diagram/<dataset>/tables; results/spatial_networks/render_network_chord_diagram/<dataset>/manifests
# Dataset behavior: runs for neuron_neuropil according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Chord diagrams of significant proteins shared between phenotype-within-unit contrasts.
# ================================================================
# Reads the mapped contrast CSVs of map_protein_identifiers.R (forward
# direction) and keeps the phenotype-within-unit contrasts, whose two sides
# name the same spatial unit. Proteins are matched across contrasts by
# ProteinGroupID and labelled with official_gene_symbol (gene_symbol is a
# display-only field in the mapping contract and repeats across protein
# groups); only protein groups with protein_level_claim_allowed are counted.
# Shared-protein workbooks go to tables/, figures to plots/. A direction with no
# shared protein gets a note instead of an empty figure.
# =========================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "preprocessing_paths.R"))
# Phase 6G: destinations resolve through the normalized output contract in
# config/output_layout.yml, addressed by this analysis's own identity.
#
# Its historical namespace was 07_spatial_networks/<substep>. Outputs already
# written there stay exactly where they are, registered LEGACY_READ_ONLY in
# config/legacy_output_registry.csv; readers still resolve them and nothing
# writes there again. The name appears only in comments, which a test
# enforces, so it cannot return as a destination.
ANALYSIS_ID <- "render_network_chord_diagram"
SUBSTEP_ID <- "chord_diagram"
NETWORK_DATASET <- current_dataset_from_cli()
assert_dataset_capability(NETWORK_DATASET, "layer", analysis = "spatial chord overlap analysis")
CANONICAL_PATHS <- canonical_module_dirs("spatial_networks", ANALYSIS_ID,
                                         scope = NETWORK_DATASET, create = FALSE)

in_dir  <- preprocessing_mapped_contrast_dir(NETWORK_DATASET, "forward")
out_dir <- CANONICAL_PATHS$figures
table_dir <- CANONICAL_PATHS$tables

# Phenotype-within-unit contrast: both sides are <unit><phenotype> with the
# same unit, e.g. CA1slmsus_CA1slmres.
contrast_unit <- function(name) {
  parts <- strsplit(name, "_", fixed = TRUE)[[1]]
  if (length(parts) != 2L) return(NA_character_)
  side <- regmatches(tolower(parts), regexec("^(.+?)(con|res|sus)$", tolower(parts)))
  units <- vapply(side, function(m) if (length(m) == 3L) m[[2]] else NA_character_, character(1))
  if (anyNA(units) || units[[1]] != units[[2]]) return(NA_character_)
  units[[1]]
}

all_csv <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
contrast_units <- vapply(tools::file_path_sans_ext(basename(all_csv)), contrast_unit, character(1))
csv_files <- all_csv[!is.na(contrast_units)]

if (is_dry_run()) {
  dry_run_line("Script", "analysis/spatial_networks/render_network_chord_diagram.R")
  dry_run_line("Mapped input directory", in_dir, if (dir.exists(in_dir)) "PASS" else "FAIL")
  dry_run_line("Phenotype-within-unit CSV count", length(csv_files), if (length(csv_files) > 0) "PASS" else "FAIL")
  dry_run_line("Other contrast CSVs skipped", length(all_csv) - length(csv_files), "INFO")
  dry_run_line("Output folders", paste(unlist(CANONICAL_PATHS), collapse = "; "))
  quit(status = if (length(csv_files) > 0) 0 else 1, save = "no")
}
if (length(csv_files) == 0) stop("No phenotype-within-unit mapped CSV files found in: ", in_dir, call. = FALSE)
for (d in c(out_dir, table_dir, CANONICAL_PATHS$logs)) dir_create(d)
write_session_info(file.path(CANONICAL_PATHS$logs, "sessionInfo.txt"))

# Packages
pkgs <- c("dplyr", "tidyr", "circlize", "RColorBrewer", "tools", "magick", "grid", "openxlsx", "digest")
missing_packages <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "),
       ". Install them explicitly before running this script.", call. = FALSE)
}
invisible(lapply(pkgs, library, character.only = TRUE))

read_both <- function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("ProteinGroupID", "official_gene_symbol", "protein_level_claim_allowed", "padj", "log2fc")
  if (!all(need %in% names(df))) {
    stop(sprintf("File %s missing columns: %s", basename(f),
                 paste(setdiff(need, names(df)), collapse = ", ")))
  }
  df <- df[as.logical(df$protein_level_claim_allowed) %in% TRUE &
             !is.na(df$ProteinGroupID) & nzchar(df$ProteinGroupID), , drop = FALSE]
  if (anyDuplicated(df$ProteinGroupID)) {
    stop(sprintf("File %s has duplicate ProteinGroupID rows.", basename(f)), call. = FALSE)
  }
  base <- tools::file_path_sans_ext(basename(f))
  list(
    labels = df[, c("ProteinGroupID", "official_gene_symbol")],
    padj = setNames(df[, c("ProteinGroupID", "padj")],   c("ProteinGroupID", base)),
    fc   = setNames(df[, c("ProteinGroupID", "log2fc")], c("ProteinGroupID", base))
  )
}

lst <- lapply(csv_files, read_both)
wide_padj <- Reduce(function(x, y) dplyr::full_join(x, y, by = "ProteinGroupID"), lapply(lst, `[[`, "padj"))
wide_fc   <- Reduce(function(x, y) dplyr::full_join(x, y, by = "ProteinGroupID"), lapply(lst, `[[`, "fc"))
wide_fc   <- wide_fc[match(wide_padj$ProteinGroupID, wide_fc$ProteinGroupID), , drop = FALSE]

# one label per protein group: its official gene symbol, else the group ID
labels <- do.call(rbind, lapply(lst, `[[`, "labels"))
labels <- labels[!is.na(labels$official_gene_symbol) & nzchar(labels$official_gene_symbol), , drop = FALSE]
labels <- labels[!duplicated(labels$ProteinGroupID), , drop = FALSE]
protein_label <- labels$official_gene_symbol[match(wide_padj$ProteinGroupID, labels$ProteinGroupID)]
protein_label[is.na(protein_label)] <- wide_padj$ProteinGroupID[is.na(protein_label)]

# -------- Build significance mask --------
alpha <- 0.05
sig <- wide_padj %>% mutate(across(-ProteinGroupID, ~ .x < alpha))
comp_cols <- setdiff(names(sig), "ProteinGroupID")

# -------- Pair overlaps by direction (inclusive) --------
pair_edges_dir <- function(dir_sign = +1, lfc_min = 0) {
  if (length(comp_cols) < 2L) return(tibble::tibble(from = character(), to = character(), weight = integer()))
  pairs <- combn(comp_cols, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    a <- p[1]; b <- p[2]
    in_both <- sig[[a]] & sig[[b]]
    if (!any(in_both, na.rm = TRUE)) return(tibble::tibble(from = a, to = b, weight = 0L))
    fa <- wide_fc[[a]][in_both]
    fb <- wide_fc[[b]][in_both]
    ok <- (sign(fa) == dir_sign) & (sign(fb) == dir_sign) &
          (abs(fa) >= lfc_min) & (abs(fb) >= lfc_min)
    tibble::tibble(from = a, to = b, weight = sum(ok, na.rm = TRUE))
  })
  bind_rows(rows) %>% filter(weight > 0)
}

edges_up   <- pair_edges_dir(+1)   # both significant and up
edges_down <- pair_edges_dir(-1)   # both significant and down

# -------- Plotting helpers --------
plot_chord <- function(edge_df, title_text) {
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
  circlize::circos.clear()
}

# sizes: convert px/res -> inches (svg() expects inches)
inch_w <- 2200 / 250
inch_h <- 2200 / 250

# One direction: SVG and PNG (the PNG feeds the side-by-side panel), or a note
# when no pair of contrasts shares a significant protein in that direction.
save_chord <- function(edge_df, label, title_text) {
  f_png <- file.path(out_dir, paste0("chord_", label, ".png"))
  f_svg <- sub("\\.png$", ".svg", f_png)
  f_note <- sub("\\.png$", "_note.txt", f_png)
  if (nrow(edge_df) == 0) {
    unlink(c(f_png, f_svg))
    writeLines(c(
      paste0("No chord diagram drawn: no pair of the ", length(comp_cols),
             " phenotype-within-unit contrasts shares a significant ", label,
             " protein (padj < ", alpha, ", protein_level_claim_allowed)."),
      paste0("Dataset: ", NETWORK_DATASET, "; input: ", in_dir)
    ), f_note)
    message("No shared ", label, " proteins; wrote ", f_note)
    return(NA_character_)
  }
  unlink(f_note)
  svg(filename = f_svg, width = inch_w, height = inch_h)
  plot_chord(edge_df, title_text)
  dev.off()
  png(f_png, width = 2200, height = 2200, res = 250)
  plot_chord(edge_df, title_text)
  dev.off()
  f_png
}

f_up   <- save_chord(edges_up, "upregulated", "Shared significant upregulated proteins")
f_down <- save_chord(edges_down, "downregulated", "Shared significant downregulated proteins")

# Combine side-by-side when both directions were drawn
f_side <- file.path(out_dir, "chord_up_down_side_by_side.png")
if (!is.na(f_up) && !is.na(f_down)) {
  img1 <- magick::image_read(f_up)
  img2 <- magick::image_read(f_down)
  combo <- magick::image_append(c(img1, img2))
  magick::image_write(combo, f_side)
} else {
  unlink(f_side)
}

# -------- Excel exports: shared proteins per pair --------
library(openxlsx)

# Build long summary for a given direction
summary_by_dir <- function(dir_sign = +1, lfc_min = 0) {
  # Direction-specific significant mask
  is_sig <- sig
  for (cc in comp_cols) {
    vv <- wide_fc[[cc]]
    is_sig[[cc]] <- is_sig[[cc]] & (sign(vv) == dir_sign) & (abs(vv) >= lfc_min)
  }

  # For each pair, collect shared proteins and return long rows
  if (length(comp_cols) < 2L) return(NULL)
  pairs <- combn(comp_cols, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    a <- p[1]; b <- p[2]
    in_both <- is_sig[[a]] & is_sig[[b]]
    if (!any(in_both, na.rm = TRUE)) return(NULL)
    idx <- which(in_both)
    tibble::tibble(
      ProteinGroupID = wide_padj$ProteinGroupID[idx],
      official_gene_symbol = protein_label[idx],
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

# Create a two-sheet workbook per direction
write_summary_wb <- function(df_long, label = "Up") {
  if (is.null(df_long) || nrow(df_long) == 0) {
    df_long <- tibble::tibble(
      ProteinGroupID = character(), official_gene_symbol = character(),
      comparison_a = character(), padj_a = numeric(), log2fc_a = numeric(),
      comparison_b = character(), padj_b = numeric(), log2fc_b = numeric()
    )
  }

  # Add per-protein stats: how many pairs it is shared in and the list of pairs
  pair_str <- paste(df_long$comparison_a, df_long$comparison_b, sep = " vs ")
  df_with_pair <- dplyr::mutate(df_long, pair = pair_str)

  per_protein <- df_with_pair %>%
    dplyr::group_by(ProteinGroupID, official_gene_symbol) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      pairs = paste(sort(unique(pair)), collapse = " | "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_pairs), official_gene_symbol, ProteinGroupID)

  # Final summary sheet: one row per protein group with ranking
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "summary")
  openxlsx::writeData(wb, "summary", per_protein, startRow = 1, startCol = 1)

  # Detailed long table
  openxlsx::addWorksheet(wb, "details")
  # Order details by protein then by decreasing evidence
  details <- df_with_pair %>%
    dplyr::mutate(max_padj = pmax(padj_a, padj_b, na.rm = TRUE),
                  mean_abs_lfc = rowMeans(cbind(abs(log2fc_a), abs(log2fc_b)), na.rm = TRUE)) %>%
    dplyr::arrange(official_gene_symbol, ProteinGroupID, max_padj, dplyr::desc(mean_abs_lfc)) %>%
    dplyr::select(ProteinGroupID, official_gene_symbol, comparison_a, padj_a, log2fc_a,
                  comparison_b, padj_b, log2fc_b)
  openxlsx::writeData(wb, "details", details, startRow = 1, startCol = 1)

  out_file <- file.path(table_dir, paste0("shared_proteins_", tolower(label), "_summary.xlsx"))
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  message("Saved ", label, " summary workbook: ", out_file, " (", nrow(per_protein), " proteins)")
  invisible(out_file)
}

# Build and save Up and Down summaries
up_long   <- summary_by_dir(+1, lfc_min = 0)
down_long <- summary_by_dir(-1, lfc_min = 0)

up_xlsx_summary   <- write_summary_wb(up_long,   label = "Up")
down_xlsx_summary <- write_summary_wb(down_long, label = "Down")

message("Contrasts: ", length(comp_cols), "; contrast pairs sharing an up protein: ", nrow(edges_up),
        "; sharing a down protein: ", nrow(edges_down))
message("Figures saved to: ", out_dir)
message("Tables saved to: ", table_dir)
