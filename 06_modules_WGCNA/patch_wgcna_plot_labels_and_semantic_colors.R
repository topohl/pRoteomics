#!/usr/bin/env Rscript
# ============================================================
# Patch WGCNA interpretable-summary plotting labels and colors
# ============================================================
# This script edits 06_modules_WGCNA/07_wgcna_interpretable_summary.r
# in a local checkout. It keeps raw WGCNA IDs unchanged, but adds:
#   - ModulePlotLabel: m01 | dark grey, etc.
#   - Supermodule_PlotLabel: canonical concise supermodule labels.
#   - MacroprogramColorKey / SemanticProgramColor for semantic QC colors.
# It also fixes the vectorized dataset_label() bug and the microglia
# composition reorder() warning.
#
# Usage from repository root:
#   Rscript 06_modules_WGCNA/patch_wgcna_plot_labels_and_semantic_colors.R
# ============================================================

options(stringsAsFactors = FALSE)

target <- file.path("06_modules_WGCNA", "07_wgcna_interpretable_summary.r")
if (!file.exists(target)) stop("Run from repository root; missing: ", target, call. = FALSE)

x <- readLines(target, warn = FALSE)
backup <- paste0(target, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
writeLines(x, backup)
message("Backup written: ", backup)

find_one <- function(pattern, fixed = TRUE, from = 1L) {
  idx <- if (fixed) grep(pattern, x[from:length(x)], fixed = TRUE) else grep(pattern, x[from:length(x)], perl = TRUE)
  if (!length(idx)) return(integer())
  from + idx[[1]] - 1L
}

replace_block <- function(start_pattern, end_regex, replacement) {
  s <- find_one(start_pattern, fixed = TRUE)
  if (!length(s)) stop("Start not found: ", start_pattern, call. = FALSE)
  e_rel <- grep(end_regex, x[s:length(x)], perl = TRUE)
  if (!length(e_rel)) stop("End not found after: ", start_pattern, call. = FALSE)
  e <- s + e_rel[[1]] - 1L
  x <<- c(x[seq_len(s - 1L)], replacement, x[(e + 1L):length(x)])
}

# ------------------------------------------------------------------
# 1) Vectorize dataset_label()
# ------------------------------------------------------------------
replace_block(
  "dataset_label <- function(ds) {",
  "^}$",
  c(
    "dataset_label <- function(ds) {",
    "  ds_chr <- as.character(ds)",
    "  dplyr::case_when(",
    "    ds_chr == \"neuron_neuropil\" ~ \"Neuron neuropil\",",
    "    ds_chr == \"neuron_soma\" ~ \"Neuron soma\",",
    "    ds_chr == \"microglia\" ~ \"Microglia ROI\",",
    "    TRUE ~ ds_chr",
    "  )",
    "}"
  )
)

# ------------------------------------------------------------------
# 2) Insert helper functions after ensure_columns()
# ------------------------------------------------------------------
if (!any(grepl("make_supermodule_plot_label <- function", x, fixed = TRUE))) {
  s <- find_one("ensure_columns <- function(df, cols) {", fixed = TRUE)
  if (!length(s)) stop("ensure_columns() not found", call. = FALSE)
  e_rel <- grep("^}$", x[s:length(x)], perl = TRUE)
  if (!length(e_rel)) stop("ensure_columns() end not found", call. = FALSE)
  e <- s + e_rel[[1]] - 1L

  helpers <- c(
    "",
    "clean_label_value <- function(x) {",
    "  x <- as.character(x)",
    "  x <- stringr::str_squish(x)",
    "  x[x %in% c(\"\", \"NA\", \"NaN\", \"Unlabelled\", \"Unassigned\")] <- NA_character_",
    "  x",
    "}",
    "",
    "hex_to_rgb01 <- function(hex) {",
    "  hex <- as.character(hex)",
    "  hex <- trimws(hex)",
    "  hex <- sub(\"^WGCNA_\", \"\", hex, ignore.case = TRUE)",
    "  hex <- sub(\"^ME\", \"\", hex)",
    "  hex <- toupper(hex)",
    "  if (!grepl(\"^#[0-9A-F]{6}$\", hex)) return(c(NA_real_, NA_real_, NA_real_))",
    "  as.numeric(grDevices::col2rgb(hex)) / 255",
    "}",
    "",
    "nearest_named_color <- function(hex) {",
    "  rgb <- hex_to_rgb01(hex)",
    "  if (any(!is.finite(rgb))) {",
    "    fallback <- as.character(hex)",
    "    fallback <- sub(\"^WGCNA_\", \"\", fallback, ignore.case = TRUE)",
    "    fallback <- sub(\"^ME\", \"\", fallback)",
    "    fallback <- gsub(\"^#\", \"\", fallback)",
    "    fallback <- gsub(\"[_-]\", \" \", fallback)",
    "    fallback <- stringr::str_squish(fallback)",
    "    if (!nzchar(fallback)) fallback <- \"unknown colour\"",
    "    return(fallback)",
    "  }",
    "  palette <- c(",
    "    black = \"#000000\", charcoal = \"#252525\", dark_grey = \"#4D4D4D\",",
    "    grey = \"#737373\", mid_grey = \"#969696\", light_grey = \"#BDBDBD\",",
    "    very_light_grey = \"#D9D9D9\", white = \"#FFFFFF\",",
    "    dark_brown = \"#7F2704\", brown = \"#8C510A\", copper_brown = \"#A6611A\",",
    "    tan = \"#BF812D\", beige = \"#DFC27D\",",
    "    dark_green = \"#006D2C\", forest_green = \"#238B45\", green = \"#41AB5D\", light_green = \"#A1D99B\",",
    "    navy = \"#08306B\", dark_blue = \"#08519C\", blue = \"#2B8CBE\", sky_blue = \"#3182BD\", light_blue = \"#9ECAE1\",",
    "    dark_purple = \"#3F007D\", purple = \"#756BB1\", lavender_purple = \"#9E9AC8\",",
    "    red = \"#CB181D\", orange = \"#E6550D\", yellow = \"#FDD049\", pink = \"#F768A1\"",
    "  )",
    "  pal_rgb <- grDevices::col2rgb(palette) / 255",
    "  d <- colSums((pal_rgb - rgb)^2)",
    "  gsub(\"_\", \" \", names(which.min(d)))",
    "}",
    "",
    "semantic_program_key <- function(x) {",
    "  z <- tolower(as.character(x))",
    "  dplyr::case_when(",
    "    grepl(\"mitochond|respirat|oxidative|\\\\batp\\\\b|acetyl|tca|electron transport\", z) ~ \"mitochondrial_metabolism\",",
    "    grepl(\"\\\\brna\\\\b|translation|ribosom|splic|mrna|ncrna|rnp\", z) ~ \"rna_translation\",",
    "    grepl(\"synap|vesicle|postsynap|presynap|cytoskeleton|actin|microtubule\", z) ~ \"synaptic_cytoskeletal\",",
    "    grepl(\"ecm|adhesion|collagen|laminin|integrin|basement membrane|extracellular matrix\", z) ~ \"ecm_adhesion\",",
    "    grepl(\"microglia|immune|phago|lysosom|complement|inflamm\", z) ~ \"microglia_state\",",
    "    grepl(\"vascular|bbb|endothelial|pericyte|blood vessel\", z) ~ \"vascular_bbb\",",
    "    grepl(\"myelin|oligodendro\", z) ~ \"myelin_oligodendrocyte\",",
    "    grepl(\"astrocyte|endfoot\", z) ~ \"astrocyte_endfoot\",",
    "    TRUE ~ \"mixed_unresolved\"",
    "  )",
    "}",
    "",
    "semantic_program_palette <- function() {",
    "  c(mitochondrial_metabolism = \"#6A51A3\", rna_translation = \"#3182BD\",",
    "    synaptic_cytoskeletal = \"#238B45\", ecm_adhesion = \"#8C510A\",",
    "    microglia_state = \"#C51B7D\", vascular_bbb = \"#E6550D\",",
    "    myelin_oligodendrocyte = \"#41B6C4\", astrocyte_endfoot = \"#756BB1\",",
    "    mixed_unresolved = \"#969696\")",
    "}",
    "",
    "semantic_program_color <- function(key) {",
    "  pal <- semantic_program_palette()",
    "  out <- unname(pal[as.character(key)])",
    "  out[is.na(out)] <- unname(pal[\"mixed_unresolved\"])",
    "  out",
    "}",
    "",
    "make_module_plot_label <- function(df) {",
    "  n <- nrow(df)",
    "  raw <- dplyr::coalesce(clean_label_value(col_or_na(df, \"module_eigengene\")),",
    "                         clean_label_value(col_or_na(df, \"module_id\")),",
    "                         clean_label_value(col_or_na(df, \"ModuleID\")),",
    "                         clean_label_value(col_or_na(df, \"ModuleColor\")),",
    "                         rep(NA_character_, n))",
    "  raw_color <- sub(\"^WGCNA_\", \"\", raw, ignore.case = TRUE)",
    "  raw_color <- sub(\"^ME\", \"\", raw_color)",
    "  display_id <- clean_label_value(col_or_na(df, \"ModuleDisplayID\"))",
    "  miss <- is.na(display_id)",
    "  if (any(miss)) {",
    "    levels <- unique(raw_color[!is.na(raw_color) & nzchar(raw_color)])",
    "    id_lookup <- stats::setNames(sprintf(\"m%02d\", seq_along(levels)), levels)",
    "    display_id[miss] <- unname(id_lookup[raw_color[miss]])",
    "  }",
    "  display_id[is.na(display_id)] <- \"m??\"",
    "  color_name <- clean_label_value(col_or_na(df, \"ModuleColorName\"))",
    "  miss_color <- is.na(color_name)",
    "  if (any(miss_color)) color_name[miss_color] <- vapply(raw_color[miss_color], nearest_named_color, character(1))",
    "  paste0(display_id, \" | \", color_name)",
    "}",
    "",
    "make_supermodule_plot_label <- function(df) {",
    "  n <- nrow(df)",
    "  id <- dplyr::coalesce(clean_label_value(col_or_na(df, \"SupermoduleID\")),",
    "                        clean_label_value(col_or_na(df, \"supermodule_id\")),",
    "                        clean_label_value(col_or_na(df, \"Supermodule_DataDrivenID\")),",
    "                        clean_label_value(col_or_na(df, \"supermodule_id_for_module\")),",
    "                        clean_label_value(col_or_na(df, \"module_supermodule_id\")),",
    "                        rep(\"SM??\", n))",
    "  display <- dplyr::coalesce(clean_label_value(col_or_na(df, \"Supermodule_DisplayLabel\")),",
    "                             clean_label_value(col_or_na(df, \"Macroprogram_Display\")),",
    "                             clean_label_value(col_or_na(df, \"Supermodule_ShortLabel\")),",
    "                             clean_label_value(col_or_na(df, \"supermodule_label_for_module\")),",
    "                             clean_label_value(col_or_na(df, \"module_supermodule_label\")),",
    "                             clean_label_value(col_or_na(df, \"Supermodule_FinalLabel\")),",
    "                             clean_label_value(col_or_na(df, \"supermodule_label\")),",
    "                             rep(NA_character_, n))",
    "  bad <- is.na(display) | display %in% c(\"Mixed / unresolved\", \"Unresolved / mixed\") | grepl(\"^Unresolved module cluster\", display, ignore.case = TRUE)",
    "  display[bad] <- paste0(id[bad], \" · Mixed / unresolved\")",
    "  needs_id <- !grepl(\"^SM[0-9]+\\\\s*[·:-]\", display)",
    "  display[needs_id] <- paste0(id[needs_id], \" · \", display[needs_id])",
    "  shorten_supermodule_label(display, max_chars = 45)",
    "}",
    "",
    "add_semantic_columns <- function(df) {",
    "  label_input <- coalesce_chr(col_or_na(df, \"Macroprogram_Display\"), col_or_na(df, \"Supermodule_PlotLabel\"),",
    "                              col_or_na(df, \"Supermodule_DisplayLabel\"), col_or_na(df, \"ModuleLabel_Final\"),",
    "                              col_or_na(df, \"best_GO_BP\"), col_or_na(df, \"best_GO_MF\"), col_or_na(df, \"best_GO_CC\"),",
    "                              col_or_na(df, \"supermodule_label_for_module\"), col_or_na(df, \"supermodule_label\"))",
    "  df$MacroprogramColorKey <- semantic_program_key(label_input)",
    "  df$SemanticProgramColor <- semantic_program_color(df$MacroprogramColorKey)",
    "  df",
    "}"
  )
  x <- c(x[seq_len(e)], helpers, x[(e + 1L):length(x)])
}

# ------------------------------------------------------------------
# 3) Add canonical labels after joins
# ------------------------------------------------------------------
needle <- "    ensure_columns(c(\"Supermodule_DisplayLabel\", \"Supermodule_LongLabel\", \"Macroprogram_Display\", \"Supermodule_FinalLabel\", \"supermodule_label\", \"supermodule_id\", \"dominant_microenvironment_class\", \"p_value\", \"estimate\", \"FDR_global\", \"contrast\", \"spatial_unit\", \"effect_scope\", \"evidence_status\", \"direction\"))"
i <- grep(needle, x, fixed = TRUE)
if (length(i) && !any(grepl("super_join <- add_semantic_columns", x, fixed = TRUE))) {
  i <- i[[1]]
  x <- c(x[seq_len(i)],
         "  super_join$Supermodule_PlotLabel <- make_supermodule_plot_label(super_join)",
         "  super_join <- add_semantic_columns(super_join)",
         x[(i + 1L):length(x)])
}

needle <- "  module_join <- attach_module_supermodules(module_join, module_super_map)"
i <- grep(needle, x, fixed = TRUE)
if (length(i) && !any(grepl("module_join$ModulePlotLabel", x, fixed = TRUE))) {
  i <- i[[1]]
  x <- c(x[seq_len(i)],
         "  module_join$ModulePlotLabel <- make_module_plot_label(module_join)",
         "  module_join$Supermodule_PlotLabel <- make_supermodule_plot_label(module_join)",
         "  module_join <- add_semantic_columns(module_join)",
         x[(i + 1L):length(x)])
}

# ------------------------------------------------------------------
# 4) Plot label substitutions
# ------------------------------------------------------------------
x <- gsub("program_label_col\\(plot_df, \\\"supermodule\\\"\\)",
          "coalesce_chr(col_or_na(plot_df, \"Supermodule_PlotLabel\"), program_label_col(plot_df, \"supermodule\"))",
          x, perl = TRUE)
x <- gsub("pretty_program_label\\(plot_df\\$supermodule_label_for_module\\)",
          "pretty_program_label(coalesce_chr(col_or_na(plot_df, \"Supermodule_PlotLabel\"), col_or_na(plot_df, \"supermodule_label_for_module\")))",
          x, perl = TRUE)
x <- gsub("pretty_program_label\\(df\\$supermodule_label_for_module\\)",
          "pretty_program_label(coalesce_chr(col_or_na(df, \"Supermodule_PlotLabel\"), col_or_na(df, \"supermodule_label_for_module\")))",
          x, perl = TRUE)
x <- gsub("as.character\\(plot_df\\$supermodule_label\\)",
          "coalesce_chr(col_or_na(plot_df, \"Supermodule_PlotLabel\"), col_or_na(plot_df, \"supermodule_label\"))",
          x, perl = TRUE)

raw_pat <- "  raw_id <- pretty_module_label(coalesce_chr(col_or_na(plot_df, \"module_id\"), col_or_na(plot_df, \"ModuleID\"), col_or_na(plot_df, \"ModuleColor\")))"
raw_rep <- "  raw_id <- coalesce_chr(col_or_na(plot_df, \"ModulePlotLabel\"), make_module_plot_label(plot_df))"
hits <- grep(raw_pat, x, fixed = TRUE)
if (length(hits)) x[hits] <- raw_rep
x <- gsub("plot_df\\$module_label_plot <- program_label_col\\(plot_df, \\\"module\\\"\\)",
          "plot_df$module_label_plot <- coalesce_chr(col_or_na(plot_df, \"ModulePlotLabel\"), program_label_col(plot_df, \"module\"))",
          x, perl = TRUE)

# ------------------------------------------------------------------
# 5) QC outputs
# ------------------------------------------------------------------
qc_needle <- "  write_table_and_source(super_join, paths$tables, paths$source_data, \"WGCNA_supermodule_group_effects_interpretable.csv\")"
i <- grep(qc_needle, x, fixed = TRUE)
if (length(i) && !any(grepl("WGCNA_supermodule_plot_label_qc.csv", x, fixed = TRUE))) {
  i <- i[[1]]
  qc <- c(
    "  supermodule_plot_label_qc <- super_join |>",
    "    dplyr::distinct(dplyr::across(dplyr::any_of(c(",
    "      \"dataset\", \"supermodule_id\", \"SupermoduleID\", \"Supermodule_DisplayLabel\",",
    "      \"Macroprogram_Display\", \"Supermodule_ShortLabel\", \"Supermodule_FinalLabel\",",
    "      \"supermodule_label\", \"Supermodule_PlotLabel\", \"MacroprogramColorKey\",",
    "      \"SemanticProgramColor\", \"Supermodule_LabelConfidence\", \"DataDrivenClusterSize\"",
    "    ))), .keep_all = FALSE)",
    "  module_plot_label_qc <- module_join |>",
    "    dplyr::distinct(dplyr::across(dplyr::any_of(c(",
    "      \"dataset\", \"module_id\", \"ModuleID\", \"ModuleColor\", \"module_eigengene\",",
    "      \"ModuleDisplayID\", \"ModuleColorName\", \"ModulePlotLabel\", \"ModuleLabel_Final\",",
    "      \"microenvironment_label\", \"supermodule_id_for_module\", \"Supermodule_PlotLabel\",",
    "      \"MacroprogramColorKey\", \"SemanticProgramColor\"",
    "    ))), .keep_all = FALSE)",
    "  write_table_and_source(supermodule_plot_label_qc, paths$tables, paths$source_data, \"WGCNA_supermodule_plot_label_qc.csv\")",
    "  write_table_and_source(module_plot_label_qc, paths$tables, paths$source_data, \"WGCNA_module_plot_label_qc.csv\")",
    ""
  )
  x <- c(x[seq_len(i - 1L)], qc, x[i:length(x)])
}

# ------------------------------------------------------------------
# 6) Cross-dataset label priority and microglia reorder warning
# ------------------------------------------------------------------
# Keep this conservative: make_cross_dataset_summary will use Supermodule_PlotLabel
# because program_label_col() now prioritizes it when present.
old_line <- "  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$fraction, y = stats::reorder(.data$program_label, .data$SupermoduleID), fill = .data$class)) +"
new_line <- "  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$fraction, y = .data$program_label, fill = .data$class)) +"
hits <- grep(old_line, x, fixed = TRUE)
if (length(hits)) x[hits] <- new_line

required <- c("make_supermodule_plot_label <- function", "make_module_plot_label <- function", "semantic_program_key <- function")
missing <- required[!vapply(required, function(tok) any(grepl(tok, x, fixed = TRUE)), logical(1))]
if (length(missing)) stop("Patch failed; missing: ", paste(missing, collapse = ", "), call. = FALSE)

writeLines(x, target)
message("Patched: ", target)
message("Validate with: Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset microglia")
