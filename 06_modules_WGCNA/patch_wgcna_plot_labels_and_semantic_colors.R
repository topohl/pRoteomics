#!/usr/bin/env Rscript
# ============================================================
# Patch WGCNA interpretable-summary plotting labels and colors
# ============================================================
# Purpose:
# - Preserve raw WGCNA identifiers such as ME#737373 / #737373.
# - Add stable display labels such as m01 | grey.
# - Add canonical supermodule plot labels.
# - Add semantic macroprogram color keys/colors for plotting/QC.
# - Fix the vectorized dataset_label() bug.
# - Fix the microglia composition reorder() warning.
#
# Usage from repo root:
#   Rscript 06_modules_WGCNA/patch_wgcna_plot_labels_and_semantic_colors.R
#
# Then rerun:
#   Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset microglia
#   Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset neuron_soma
#   Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset neuron_neuropil
#   Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all
# ============================================================

options(stringsAsFactors = FALSE)

script_path <- file.path("06_modules_WGCNA", "07_wgcna_interpretable_summary.r")
if (!file.exists(script_path)) {
  stop("Cannot find ", script_path, ". Run this script from the repository root.", call. = FALSE)
}

x <- readLines(script_path, warn = FALSE)
original <- x
backup_path <- paste0(script_path, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
writeLines(original, backup_path)
message("Backup written: ", backup_path)

replace_once <- function(text, pattern, replacement, fixed = TRUE, label = pattern) {
  hit <- if (fixed) grep(pattern, text, fixed = TRUE) else grep(pattern, text, perl = TRUE)
  if (!length(hit)) stop("Pattern not found for replacement: ", label, call. = FALSE)
  if (length(hit) > 1L) message("Pattern matched ", length(hit), " times; replacing first occurrence for: ", label)
  i <- hit[[1]]
  if (fixed) {
    text[[i]] <- sub(pattern, replacement, text[[i]], fixed = TRUE)
  } else {
    text[[i]] <- sub(pattern, replacement, text[[i]], perl = TRUE)
  }
  text
}

replace_block <- function(text, start_pattern, end_pattern, replacement_lines, fixed = TRUE, label = start_pattern) {
  s <- if (fixed) grep(start_pattern, text, fixed = TRUE) else grep(start_pattern, text, perl = TRUE)
  if (!length(s)) stop("Start pattern not found: ", label, call. = FALSE)
  s <- s[[1]]
  e_rel <- if (fixed) grep(end_pattern, text[s:length(text)], fixed = TRUE) else grep(end_pattern, text[s:length(text)], perl = TRUE)
  if (!length(e_rel)) stop("End pattern not found after start for: ", label, call. = FALSE)
  e <- s + e_rel[[1]] - 1L
  c(text[seq_len(s - 1L)], replacement_lines, text[(e + 1L):length(text)])
}

insert_after <- function(text, pattern, insert_lines, fixed = TRUE, label = pattern) {
  hit <- if (fixed) grep(pattern, text, fixed = TRUE) else grep(pattern, text, perl = TRUE)
  if (!length(hit)) stop("Insert pattern not found: ", label, call. = FALSE)
  i <- hit[[1]]
  c(text[seq_len(i)], insert_lines, text[(i + 1L):length(text)])
}

# 1. Vectorize dataset_label().
x <- replace_block(
  x,
  "dataset_label <- function(ds) {",
  "}",
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
  ),
  label = "dataset_label function"
)

# 2. Insert label/color helpers after ensure_columns().
helper_block <- c(
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
  "  c(",
  "    mitochondrial_metabolism = \"#6A51A3\",",
  "    rna_translation = \"#3182BD\",",
  "    synaptic_cytoskeletal = \"#238B45\",",
  "    ecm_adhesion = \"#8C510A\",",
  "    microglia_state = \"#C51B7D\",",
  "    vascular_bbb = \"#E6550D\",",
  "    myelin_oligodendrocyte = \"#41B6C4\",",
  "    astrocyte_endfoot = \"#756BB1\",",
  "    mixed_unresolved = \"#969696\"",
  "  )",
  "}",
  "",
  "semantic_program_color <- function(key) {",
  "  pal <- semantic_program_palette()",
  "  key <- as.character(key)",
  "  out <- unname(pal[key])",
  "  out[is.na(out)] <- unname(pal[\"mixed_unresolved\"])",
  "  out",
  "}",
  "",
  "make_module_plot_label <- function(df) {",
  "  n <- nrow(df)",
  "  raw <- dplyr::coalesce(",
  "    clean_label_value(col_or_na(df, \"module_eigengene\")),",
  "    clean_label_value(col_or_na(df, \"module_id\")),",
  "    clean_label_value(col_or_na(df, \"ModuleID\")),",
  "    clean_label_value(col_or_na(df, \"ModuleColor\")),",
  "    rep(NA_character_, n)",
  "  )",
  "  raw_color <- sub(\"^WGCNA_\", \"\", raw, ignore.case = TRUE)",
  "  raw_color <- sub(\"^ME\", \"\", raw_color)",
  "  display_id <- clean_label_value(col_or_na(df, \"ModuleDisplayID\"))",
  "  missing_id <- is.na(display_id)",
  "  if (any(missing_id)) {",
  "    stable_levels <- unique(raw_color[!is.na(raw_color) & nzchar(raw_color)])",
  "    id_lookup <- stats::setNames(sprintf(\"m%02d\", seq_along(stable_levels)), stable_levels)",
  "    display_id[missing_id] <- unname(id_lookup[raw_color[missing_id]])",
  "  }",
  "  display_id[is.na(display_id)] <- \"m??\"",
  "  color_name <- clean_label_value(col_or_na(df, \"ModuleColorName\"))",
  "  missing_color <- is.na(color_name)",
  "  if (any(missing_color)) color_name[missing_color] <- vapply(raw_color[missing_color], nearest_named_color, character(1))",
  "  paste0(display_id, \" | \<\", color_name, \"\") |> gsub(\" \\| <\", \" | \", x = _)",
  "}",
  "",
  "make_supermodule_plot_label <- function(df) {",
  "  n <- nrow(df)",
  "  id <- dplyr::coalesce(",
  "    clean_label_value(col_or_na(df, \"SupermoduleID\")),",
  "    clean_label_value(col_or_na(df, \"supermodule_id\")),",
  "    clean_label_value(col_or_na(df, \"Supermodule_DataDrivenID\")),",
  "    clean_label_value(col_or_na(df, \"supermodule_id_for_module\")),",
  "    clean_label_value(col_or_na(df, \"module_supermodule_id\")),",
  "    rep(\"SM??\", n)",
  "  )",
  "  display <- dplyr::coalesce(",
  "    clean_label_value(col_or_na(df, \"Supermodule_DisplayLabel\")),",
  "    clean_label_value(col_or_na(df, \"Macroprogram_Display\")),",
  "    clean_label_value(col_or_na(df, \"Supermodule_ShortLabel\")),",
  "    clean_label_value(col_or_na(df, \"supermodule_label_for_module\")),",
  "    clean_label_value(col_or_na(df, \"module_supermodule_label\")),",
  "    clean_label_value(col_or_na(df, \"Supermodule_FinalLabel\")),",
  "    clean_label_value(col_or_na(df, \"supermodule_label\")),",
  "    rep(NA_character_, n)",
  "  )",
  "  bad <- is.na(display) | display %in% c(\"Mixed / unresolved\", \"Unresolved / mixed\") | grepl(\"^Unresolved module cluster\", display, ignore.case = TRUE)",
  "  display[bad] <- paste0(id[bad], \" · Mixed / unresolved\")",
  "  needs_id <- !grepl(\"^SM[0-9]+\\\\s*[·:-]\", display)",
  "  display[needs_id] <- paste0(id[needs_id], \" · \", display[needs_id])",
  "  shorten_supermodule_label(display, max_chars = 45)",
  "}",
  "",
  "add_semantic_columns <- function(df) {",
  "  label_input <- coalesce_chr(",
  "    col_or_na(df, \"Macroprogram_Display\"),",
  "    col_or_na(df, \"Supermodule_PlotLabel\"),",
  "    col_or_na(df, \"Supermodule_DisplayLabel\"),",
  "    col_or_na(df, \"ModuleLabel_Final\"),",
  "    col_or_na(df, \"best_GO_BP\"),",
  "    col_or_na(df, \"best_GO_MF\"),",
  "    col_or_na(df, \"best_GO_CC\"),",
  "    col_or_na(df, \"supermodule_label_for_module\"),",
  "    col_or_na(df, \"supermodule_label\")",
  "  )",
  "  df$MacroprogramColorKey <- semantic_program_key(label_input)",
  "  df$SemanticProgramColor <- semantic_program_color(df$MacroprogramColorKey)",
  "  df",
  "}"
)

if (!any(grepl("make_supermodule_plot_label <- function", x, fixed = TRUE))) {
  x <- insert_after(x, "ensure_columns <- function(df, cols) {", helper_block, label = "after ensure_columns")
}

# The helper insertion above inserts too early if placed inside ensure_columns; correct by moving only if needed.
# Safer fallback: if helper block landed inside ensure_columns, abort before writing.
if (any(grepl("for \(nm in cols\)", x, perl = TRUE)) && any(grepl("clean_label_value <- function", x, fixed = TRUE))) {
  clean_idx <- grep("clean_label_value <- function", x, fixed = TRUE)[[1]]
  ensure_idx <- grep("ensure_columns <- function", x, fixed = TRUE)[[1]]
  effect_idx <- grep("effect_sentence <- function", x, fixed = TRUE)[[1]]
  if (clean_idx > ensure_idx && clean_idx < effect_idx) {
    # Rebuild the early helper area deterministically.
    x <- original
    x <- replace_block(
      x,
      "dataset_label <- function(ds) {",
      "}",
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
      ),
      label = "dataset_label function"
    )
    end_ensure <- grep("^}$", x)
    ensure_start <- grep("ensure_columns <- function", x, fixed = TRUE)[[1]]
    ensure_end <- end_ensure[end_ensure > ensure_start][[1]]
    x <- c(x[seq_len(ensure_end)], helper_block, x[(ensure_end + 1L):length(x)])
  }
}

# 3. Prefer Macroprogram_Display before verbose Supermodule_FinalLabel in program_label_col().
x <- gsub(
  "col_or_na\\(df, \\\"Supermodule_DisplayLabel\\\"\\),\\n      col_or_na\\(df, \\\"Supermodule_FinalLabel\\\"\\),\\n      col_or_na\\(df, \\\"Supermodule_FinalLabel.y\\\"\\),\\n      col_or_na\\(df, \\\"Supermodule_FinalLabel.x\\\"\\),\\n      col_or_na\\(df, \\\"Macroprogram_Display\\\"\\)",
  "col_or_na(df, \"Supermodule_PlotLabel\"),\n      col_or_na(df, \"Supermodule_DisplayLabel\"),\n      col_or_na(df, \"Macroprogram_Display\"),\n      col_or_na(df, \"Supermodule_ShortLabel\"),\n      col_or_na(df, \"Supermodule_FinalLabel\"),\n      col_or_na(df, \"Supermodule_FinalLabel.y\"),\n      col_or_na(df, \"Supermodule_FinalLabel.x\")",
  x,
  perl = TRUE
)

# 4. Add plot labels after super_join/module_join construction.
x <- replace_once(
  x,
  "    ensure_columns(c(\"Supermodule_DisplayLabel\", \"Supermodule_LongLabel\", \"Macroprogram_Display\", \"Supermodule_FinalLabel\", \"supermodule_label\", \"supermodule_id\", \"dominant_microenvironment_class\", \"p_value\", \"estimate\", \"FDR_global\", \"contrast\", \"spatial_unit\", \"effect_scope\", \"evidence_status\", \"direction\"))",
  "    ensure_columns(c(\"Supermodule_DisplayLabel\", \"Supermodule_LongLabel\", \"Macroprogram_Display\", \"Supermodule_FinalLabel\", \"supermodule_label\", \"supermodule_id\", \"dominant_microenvironment_class\", \"p_value\", \"estimate\", \"FDR_global\", \"contrast\", \"spatial_unit\", \"effect_scope\", \"evidence_status\", \"direction\"))\n  super_join$Supermodule_PlotLabel <- make_supermodule_plot_label(super_join)\n  super_join <- add_semantic_columns(super_join)",
  label = "super_join ensure columns"
)

x <- replace_once(
  x,
  "  module_join <- attach_module_supermodules(module_join, module_super_map)",
  "  module_join <- attach_module_supermodules(module_join, module_super_map)\n  module_join$ModulePlotLabel <- make_module_plot_label(module_join)\n  module_join$Supermodule_PlotLabel <- make_supermodule_plot_label(module_join)\n  module_join <- add_semantic_columns(module_join)",
  label = "attach_module_supermodules"
)

# 5. Use canonical labels in selected plots.
x <- gsub("program_label_col\\(plot_df, \\\"supermodule\\\"\\)", "coalesce_chr(col_or_na(plot_df, \"Supermodule_PlotLabel\"), program_label_col(plot_df, \"supermodule\"))", x, perl = TRUE)
x <- gsub("pretty_program_label\\(plot_df\\$supermodule_label_for_module\\)", "pretty_program_label(coalesce_chr(col_or_na(plot_df, \"Supermodule_PlotLabel\"), col_or_na(plot_df, \"supermodule_label_for_module\")))", x, perl = TRUE)
x <- gsub("pretty_program_label\\(df\\$supermodule_label_for_module\\)", "pretty_program_label(coalesce_chr(col_or_na(df, \"Supermodule_PlotLabel\"), col_or_na(df, \"supermodule_label_for_module\")))", x, perl = TRUE)
x <- gsub("as.character\\(plot_df\\$supermodule_label\\)", "coalesce_chr(col_or_na(plot_df, \"Supermodule_PlotLabel\"), col_or_na(plot_df, \"supermodule_label\"))", x, perl = TRUE)

# 6. Improve module labels where raw ME#/WGCNA_# labels were used.
x <- replace_once(
  x,
  "  raw_id <- pretty_module_label(coalesce_chr(col_or_na(plot_df, \"module_id\"), col_or_na(plot_df, \"ModuleID\"), col_or_na(plot_df, \"ModuleColor\")))",
  "  raw_id <- coalesce_chr(col_or_na(plot_df, \"ModulePlotLabel\"), make_module_plot_label(plot_df))",
  label = "module main raw_id"
)
x <- replace_once(
  x,
  "  raw_id <- pretty_module_label(coalesce_chr(col_or_na(plot_df, \"module_id\"), col_or_na(plot_df, \"ModuleID\"), col_or_na(plot_df, \"ModuleColor\")))",
  "  raw_id <- coalesce_chr(col_or_na(plot_df, \"ModulePlotLabel\"), make_module_plot_label(plot_df))",
  label = "module spatial raw_id"
)
x <- gsub("plot_df\\$module_label_plot <- program_label_col\\(plot_df, \\\"module\\\"\\)", "plot_df$module_label_plot <- coalesce_chr(col_or_na(plot_df, \"ModulePlotLabel\"), program_label_col(plot_df, \"module\"))", x, perl = TRUE)

# 7. Fix microglia composition label priority and reorder warning.
x <- gsub(
  "col_or_na\\(super_annot, \\\"Supermodule_DisplayLabel\\\"\\),\\n        col_or_na\\(super_annot, \\\"Supermodule_FinalLabel\\\"\\),\\n        col_or_na\\(super_annot, \\\"Macroprogram_Display\\\"\\)",
  "col_or_na(super_annot, \"Supermodule_DisplayLabel\"),\n        col_or_na(super_annot, \"Macroprogram_Display\"),\n        col_or_na(super_annot, \"Supermodule_ShortLabel\"),\n        col_or_na(super_annot, \"Supermodule_FinalLabel\")",
  x,
  perl = TRUE
)
x <- gsub(
  "ggplot2::aes\\(x = \\.data\\$fraction, y = stats::reorder\\(\\.data\\$program_label, \\.data\\$SupermoduleID\\), fill = \\.data\\$class\\)",
  "ggplot2::aes(x = .data$fraction, y = .data$program_label, fill = .data$class)",
  x,
  perl = TRUE
)

# Insert deterministic ordering before microglia composition ggplot if not present.
needle <- "  if (!nrow(plot_df)) return(invisible(NULL))\n  write_csv_safe2(plot_df, file.path(paths$source_data, \"microglia_supermodule_annotation_composition_source.csv\"))"
replacement <- "  if (!nrow(plot_df)) return(invisible(NULL))\n  plot_df <- plot_df |>\n    dplyr::mutate(\n      supermodule_order = suppressWarnings(as.integer(gsub(\"[^0-9]\", \"\", .data$SupermoduleID))),\n      supermodule_order = dplyr::coalesce(.data$supermodule_order, dplyr::row_number())\n    ) |>\n    dplyr::arrange(.data$supermodule_order, .data$program_label) |>\n    dplyr::mutate(program_label = factor(.data$program_label, levels = rev(unique(.data$program_label))))\n  write_csv_safe2(plot_df, file.path(paths$source_data, \"microglia_supermodule_annotation_composition_source.csv\"))"
if (any(grepl(needle, x, fixed = TRUE))) x <- replace_once(x, needle, replacement, label = "microglia composition ordering")

# 8. Add QC tables before main outputs are written.
qc_needle <- "  write_table_and_source(super_join, paths$tables, paths$source_data, \"WGCNA_supermodule_group_effects_interpretable.csv\")"
qc_block <- paste0(
  "  supermodule_plot_label_qc <- super_join |>\n",
  "    dplyr::distinct(dplyr::across(dplyr::any_of(c(\n",
  "      \"dataset\", \"supermodule_id\", \"SupermoduleID\", \"Supermodule_DisplayLabel\",\n",
  "      \"Macroprogram_Display\", \"Supermodule_ShortLabel\", \"Supermodule_FinalLabel\",\n",
  "      \"supermodule_label\", \"Supermodule_PlotLabel\", \"MacroprogramColorKey\",\n",
  "      \"SemanticProgramColor\", \"Supermodule_LabelConfidence\", \"DataDrivenClusterSize\"\n",
  "    ))), .keep_all = FALSE)\n",
  "  module_plot_label_qc <- module_join |>\n",
  "    dplyr::distinct(dplyr::across(dplyr::any_of(c(\n",
  "      \"dataset\", \"module_id\", \"ModuleID\", \"ModuleColor\", \"module_eigengene\",\n",
  "      \"ModuleDisplayID\", \"ModuleColorName\", \"ModulePlotLabel\", \"ModuleLabel_Final\",\n",
  "      \"microenvironment_label\", \"supermodule_id_for_module\", \"Supermodule_PlotLabel\",\n",
  "      \"MacroprogramColorKey\", \"SemanticProgramColor\"\n",
  "    ))), .keep_all = FALSE)\n",
  "  write_table_and_source(supermodule_plot_label_qc, paths$tables, paths$source_data, \"WGCNA_supermodule_plot_label_qc.csv\")\n",
  "  write_table_and_source(module_plot_label_qc, paths$tables, paths$source_data, \"WGCNA_module_plot_label_qc.csv\")\n\n",
  qc_needle
)
if (any(grepl(qc_needle, x, fixed = TRUE)) && !any(grepl("WGCNA_supermodule_plot_label_qc.csv", x, fixed = TRUE))) {
  x <- replace_once(x, qc_needle, qc_block, label = "QC table insertion")
}

# 9. Cross-dataset: prefer canonical labels and fix vectorized dataset_label handled above.
x <- gsub(
  "col_or_na\\(all_super, \\\"Supermodule_DisplayLabel\\\"\\),\\n        col_or_na\\(all_super, \\\"Supermodule_FinalLabel\\\"\\),\\n        col_or_na\\(all_super, \\\"Macroprogram_Display\\\"\\)",
  "col_or_na(all_super, \"Supermodule_PlotLabel\"),\n        col_or_na(all_super, \"Supermodule_DisplayLabel\"),\n        col_or_na(all_super, \"Macroprogram_Display\"),\n        col_or_na(all_super, \"Supermodule_FinalLabel\")",
  x,
  perl = TRUE
)

# 10. Basic sanity checks before writing.
text_all <- paste(x, collapse = "\n")
required_tokens <- c(
  "make_supermodule_plot_label <- function",
  "make_module_plot_label <- function",
  "semantic_program_key <- function",
  "Supermodule_PlotLabel",
  "ModulePlotLabel",
  "WGCNA_supermodule_plot_label_qc.csv",
  "WGCNA_module_plot_label_qc.csv"
)
missing_tokens <- required_tokens[!vapply(required_tokens, grepl, logical(1), x = text_all, fixed = TRUE)]
if (length(missing_tokens)) {
  stop("Patch sanity check failed; missing token(s): ", paste(missing_tokens, collapse = ", "), call. = FALSE)
}

# Avoid writing known malformed placeholder from early drafts.
if (grepl("\\| <", text_all, fixed = TRUE)) {
  text_all <- gsub("paste0\\(display_id, \\\" \\\\| \\\".*?\\)", "paste0(display_id, \" | \", color_name)", text_all, perl = TRUE)
  x <- strsplit(text_all, "\n", fixed = TRUE)[[1]]
}

writeLines(x, script_path)
message("Patched: ", script_path)
message("Next: run Rscript ", script_path, " --dataset microglia to validate before full rerun.")
