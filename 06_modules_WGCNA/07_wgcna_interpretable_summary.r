#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/07_wgcna_interpretable_summary.r
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: required results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv; results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv; +2 more; optional results/tables/06_modules_WGCNA/04_wgcna_de_gsea_overlap/<dataset>/.
# Produces: results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_interpretable_summary.xlsx.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Runs after module/supermodule group effects and microenvironment annotation.
# ================================================================

#
# Combine WGCNA group effects and biological annotation into interpretable summaries.
#
# Design intent:
#   - keep the primary inference from 05_module_supermodule_group_effects.r
#   - join biological annotation from 06_annotate_module_microenvironment.r
#   - produce readable, publication-ready summary figures
#   - avoid overloading one figure with spatial-unit, contrast, dataset, and annotation dimensions at once

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_labeling_utils.R"))
source(repo_path("R", "schema_validation.R"))

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "readr", "stringr", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
if (!length(missing_pkgs)) {
  suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))
}

run <- wgcna_cli(allow_all = TRUE)
DATASET_ARG <- run$dataset

dataset_label <- function(ds) {
  ds_chr <- as.character(ds)
  dplyr::case_when(
    ds_chr == "neuron_neuropil" ~ "Neuron neuropil",
    ds_chr == "neuron_soma" ~ "Neuron soma",
    ds_chr == "microglia" ~ "Microglia ROI",
    TRUE ~ ds_chr
  )
}

spatial_unit_label <- function(ds) {
  if (identical(ds, "neuron_neuropil")) "region-layer" else "region"
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

first_existing_col <- function(df, candidates, fallback = NULL) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else fallback
}

coalesce_chr <- function(...) {
  x <- dplyr::coalesce(...)
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "Unlabelled"
  x
}

label_wrap <- function(x, width = 36) {
  vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

clip_p <- function(p) pmax(safe_num(p), 1e-300)

col_or_na <- function(df, nm) {
  if (nm %in% names(df)) return(df[[nm]])
  rep(NA_character_, nrow(df))
}

ensure_columns <- function(df, cols) {
  for (nm in cols) {
    if (!nm %in% names(df)) df[[nm]] <- rep(NA, nrow(df))
  }
  df
}

clean_label_value <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[x %in% c("", "NA", "NaN", "Unlabelled", "Unassigned")] <- NA_character_
  x
}

supermodule_label_sep <- function() paste0(" ", intToUtf8(183), " ")

is_bad_supermodule_label <- function(x) {
  x <- stringr::str_squish(as.character(x))
  is.na(x) | !nzchar(trimws(x)) |
    x %in% c("Mixed / unresolved", "Unresolved / mixed", "Unlabelled", "Unassigned") |
    grepl(paste0("^SM[0-9]+\\s*(", intToUtf8(183), "|:|-)\\s*Mixed / unresolved$"), x)
}

clean_supermodule_label_value <- function(x) {
  x <- clean_label_value(x)
  x[is_bad_supermodule_label(x)] <- NA_character_
  x
}

displayed_supermodule_theme_count <- function(label) {
  vapply(as.character(label), function(z) {
    z <- stringr::str_squish(z)
    if (is.na(z) || !nzchar(z)) return(0L)
    z <- sub(paste0("^SM[0-9]+\\s*(", intToUtf8(183), "|:|-)\\s*"), "", z, ignore.case = TRUE)
    z <- sub("^[^:]{1,80}:\\s*(?=(mostly\\s+|mixed:|mixed\\s+multi-program|mixed /))", "", z, perl = TRUE, ignore.case = TRUE)
    if (grepl("^mixed\\s+multi-program$", z, ignore.case = TRUE)) return(0L)
    if (grepl("^mixed:\\s*", z, ignore.case = TRUE)) {
      parts <- trimws(unlist(strsplit(sub("^mixed:\\s*", "", z, ignore.case = TRUE), "\\s*;\\s*", perl = TRUE), use.names = FALSE))
      parts <- parts[nzchar(parts) & !grepl("mixed / low-specificity|mixed / unresolved|unresolved / mixed", parts, ignore.case = TRUE)]
      return(max(2L, length(unique(parts))))
    }
    z <- sub("^mostly\\s+", "", z, ignore.case = TRUE)
    z <- sub("^mixed:\\s*", "", z, ignore.case = TRUE)
    parts <- trimws(unlist(strsplit(z, "\\s*;\\s*", perl = TRUE), use.names = FALSE))
    parts <- parts[nzchar(parts) & !grepl("mixed / low-specificity|mixed / unresolved|unresolved / mixed", parts, ignore.case = TRUE)]
    length(unique(parts))
  }, integer(1))
}

member_theme_fraction_values <- function(x) {
  parts <- trimws(unlist(strsplit(as.character(x), "\\s*;\\s*", perl = TRUE), use.names = FALSE))
  parts <- parts[nzchar(parts) & grepl("=", parts, fixed = TRUE)]
  vals <- suppressWarnings(as.numeric(sub("^.*=", "", parts)))
  names(vals) <- trimws(sub("=.*$", "", parts))
  vals[is.finite(vals)]
}

supermodule_theme_qc_warnings <- function(df) {
  n <- nrow(df)
  if (!n) return(character())
  out <- character(n)
  plot_label <- as.character(col_or_na(df, "Supermodule_PlotLabel"))
  comp_label <- as.character(col_or_na(df, "Supermodule_CompositionLabel"))
  fractions <- as.character(col_or_na(df, "MemberThemeFractions"))
  omitted <- as.character(col_or_na(df, "themes_omitted_from_display_label"))
  n_distinct <- suppressWarnings(as.integer(col_or_na(df, "n_distinct_member_themes")))
  dominant_frac <- suppressWarnings(as.numeric(col_or_na(df, "DominantMemberThemeFraction")))
  shown_n <- displayed_supermodule_theme_count(plot_label)

  for (i in seq_len(n)) {
    warn <- character()
    vals <- member_theme_fraction_values(fractions[[i]])
    informative_vals <- vals[names(vals) != "mixed / low-specificity"]
    if (!length(vals)) {
      warn <- c(warn, "MemberThemeFractions is empty")
    }
    if (!is.na(omitted[[i]]) && nzchar(trimws(omitted[[i]]))) {
      warn <- c(warn, "themes_omitted_from_display_label is non-empty")
    }
    if (is.finite(n_distinct[[i]]) && n_distinct[[i]] >= 3L &&
        length(informative_vals) >= 3L && all(informative_vals >= 0.20) &&
        shown_n[[i]] < 3L) {
      warn <- c(warn, "plot label shows fewer than 3 displayed themes despite >=3 member themes all at fraction >=0.20")
    }
    single_theme_display <- shown_n[[i]] <= 1L &&
      !grepl("^SM[0-9]+\\s*(·|:|-)\\s*(mixed:|mixed multi-program|mixed / low-specificity|mixed / unresolved)", plot_label[[i]], ignore.case = TRUE) &&
      !grepl("^(mixed:|mixed multi-program|mixed / low-specificity|mixed / unresolved)", comp_label[[i]], ignore.case = TRUE)
    if ((!is.finite(dominant_frac[[i]]) || dominant_frac[[i]] < 0.60) && single_theme_display) {
      warn <- c(warn, "supermodule without dominant theme >=0.60 is displayed as a single-theme program")
    }
    out[[i]] <- paste(unique(warn), collapse = " | ")
  }
  out
}

add_supermodule_id_prefix <- function(id, label) {
  id <- clean_label_value(id)
  label <- clean_label_value(label)
  out <- label
  needs_id <- !is.na(id) & !is.na(out) & !grepl(paste0("^SM[0-9]+\\s*(", intToUtf8(183), "|:|-)"), out)
  out[needs_id] <- paste0(id[needs_id], supermodule_label_sep(), out[needs_id])
  out
}

hex_to_rgb01 <- function(hex) {
  hex <- as.character(hex)
  hex <- trimws(hex)
  hex <- sub("^WGCNA_", "", hex, ignore.case = TRUE)
  hex <- sub("^ME", "", hex)
  hex <- toupper(hex)
  if (!grepl("^#[0-9A-F]{6}$", hex)) return(c(NA_real_, NA_real_, NA_real_))
  as.numeric(grDevices::col2rgb(hex)) / 255
}

nearest_named_color <- function(hex) {
  rgb <- hex_to_rgb01(hex)
  if (any(!is.finite(rgb))) {
    fallback <- as.character(hex)
    fallback <- sub("^WGCNA_", "", fallback, ignore.case = TRUE)
    fallback <- sub("^ME", "", fallback)
    fallback <- gsub("^#", "", fallback)
    fallback <- gsub("[_-]", " ", fallback)
    fallback <- stringr::str_squish(fallback)
    if (!nzchar(fallback)) fallback <- "unknown colour"
    return(fallback)
  }
  palette <- c(
    black = "#000000", charcoal = "#252525", dark_grey = "#4D4D4D",
    grey = "#737373", mid_grey = "#969696", light_grey = "#BDBDBD",
    very_light_grey = "#D9D9D9", white = "#FFFFFF",
    dark_brown = "#7F2704", brown = "#8C510A", copper_brown = "#A6611A",
    tan = "#BF812D", beige = "#DFC27D",
    dark_green = "#006D2C", forest_green = "#238B45", green = "#41AB5D", light_green = "#A1D99B",
    navy = "#08306B", dark_blue = "#08519C", blue = "#2B8CBE", sky_blue = "#3182BD", light_blue = "#9ECAE1",
    dark_purple = "#3F007D", purple = "#756BB1", lavender_purple = "#9E9AC8",
    red = "#CB181D", orange = "#E6550D", yellow = "#FDD049", pink = "#F768A1"
  )
  pal_rgb <- grDevices::col2rgb(palette) / 255
  d <- colSums((pal_rgb - rgb)^2)
  gsub("_", " ", names(which.min(d)))
}

semantic_program_key <- function(x) {
  z <- tolower(as.character(x))
  dplyr::case_when(
    grepl("mitochond|respirat|oxidative|\\batp\\b|acetyl|tca|electron transport", z) ~ "mitochondrial_metabolism",
    grepl("\\brna\\b|translation|ribosom|splic|mrna|ncrna|rnp", z) ~ "rna_translation",
    grepl("synap|vesicle|postsynap|presynap|cytoskeleton|actin|microtubule", z) ~ "synaptic_cytoskeletal",
    grepl("ecm|extracellular matrix|basement membrane|collagen|laminin|nidogen|\\bagrn\\b|hspg2|perlecan|\\bcol4a[12]?\\b|\\blama[0-9]?\\b|\\blamb[0-9]?\\b|\\blamc[0-9]?\\b|\\bnid[12]\\b|\\bbcam\\b|serpinh1", z) ~ "ecm_adhesion",
    grepl("microglia|immune|phago|lysosom|complement|inflamm", z) ~ "microglia_state",
    grepl("vascular|bbb|endothelial|pericyte|blood vessel", z) ~ "vascular_bbb",
    grepl("myelin|oligodendro", z) ~ "myelin_oligodendrocyte",
    grepl("astrocyte|endfoot", z) ~ "astrocyte_endfoot",
    TRUE ~ "mixed_unresolved"
  )
}

semantic_program_palette <- function() {
  c(mitochondrial_metabolism = "#6A51A3", rna_translation = "#3182BD",
    synaptic_cytoskeletal = "#238B45", ecm_adhesion = "#8C510A",
    microglia_state = "#C51B7D", vascular_bbb = "#E6550D",
    myelin_oligodendrocyte = "#41B6C4", astrocyte_endfoot = "#756BB1",
    mixed_unresolved = "#969696")
}

semantic_program_color <- function(key) {
  pal <- semantic_program_palette()
  out <- unname(pal[as.character(key)])
  out[is.na(out)] <- unname(pal["mixed_unresolved"])
  out
}

make_module_plot_label <- function(df) {
  n <- nrow(df)
  raw <- dplyr::coalesce(clean_label_value(col_or_na(df, "module_eigengene")),
                         clean_label_value(col_or_na(df, "module_id")),
                         clean_label_value(col_or_na(df, "ModuleID")),
                         clean_label_value(col_or_na(df, "ModuleColor")),
                         rep(NA_character_, n))
  raw_color <- sub("^WGCNA_", "", raw, ignore.case = TRUE)
  raw_color <- sub("^ME", "", raw_color)
  display_id <- clean_label_value(col_or_na(df, "ModuleDisplayID"))
  color_name <- clean_label_value(col_or_na(df, "ModuleColorName"))
  miss_color <- is.na(color_name)
  if (any(miss_color)) color_name[miss_color] <- vapply(raw_color[miss_color], nearest_named_color, character(1))
  out <- color_name
  has_display_id <- !is.na(display_id)
  out[has_display_id] <- paste0(display_id[has_display_id], " | ", color_name[has_display_id])
  out[is.na(out)] <- "unlabelled module"
  out
}

make_supermodule_plot_label <- function(df) {
  n <- nrow(df)
  id <- dplyr::coalesce(clean_label_value(col_or_na(df, "SupermoduleID")),
                        clean_label_value(col_or_na(df, "supermodule_id")),
                        clean_label_value(col_or_na(df, "Supermodule_DataDrivenID")),
                        clean_label_value(col_or_na(df, "supermodule_id_for_module")),
                        clean_label_value(col_or_na(df, "module_supermodule_id")),
                        rep("SM??", n))
  display <- dplyr::coalesce(clean_supermodule_label_value(col_or_na(df, "Supermodule_DisplayLabel")),
                             clean_supermodule_label_value(col_or_na(df, "Macroprogram_Display")),
                             clean_supermodule_label_value(col_or_na(df, "Supermodule_FinalLabel")),
                             clean_supermodule_label_value(col_or_na(df, "Supermodule_ShortLabel")),
                             clean_supermodule_label_value(col_or_na(df, "supermodule_label_for_module")),
                             clean_supermodule_label_value(col_or_na(df, "module_supermodule_label")),
                             clean_supermodule_label_value(col_or_na(df, "supermodule_label")),
                             rep(NA_character_, n))
  bad <- is.na(display) | grepl("^Unresolved module cluster", display, ignore.case = TRUE)
  display[bad] <- paste0(id[bad], supermodule_label_sep(), "Mixed / unresolved")
  display <- add_supermodule_id_prefix(id, display)
  shorten_supermodule_label(display, max_chars = 45)
}

add_semantic_columns <- function(df) {
  label_input <- coalesce_chr(col_or_na(df, "Macroprogram_Display"), col_or_na(df, "Supermodule_PlotLabel"),
                              col_or_na(df, "Supermodule_DisplayLabel"), col_or_na(df, "ModuleLabel_Final"),
                              col_or_na(df, "best_GO_BP"), col_or_na(df, "best_GO_MF"), col_or_na(df, "best_GO_CC"),
                              col_or_na(df, "supermodule_label_for_module"), col_or_na(df, "supermodule_label"))
  df$MacroprogramColorKey <- semantic_program_key(label_input)
  df$SemanticProgramColor <- semantic_program_color(df$MacroprogramColorKey)
  df
}

roman_label <- function(i) {
  out <- as.character(utils::as.roman(i))
  out[is.na(out)] <- as.character(i[is.na(out)])
  out
}

supermodule_broad_label <- function(x) {
  z <- tolower(as.character(x))
  dplyr::case_when(
    grepl("targeted microglia signature overlap", z) ~ "microglia-signature overlap",
    grepl("microglia[- ]enriched roi|microglia[- ]associated roi", z) ~ "microglia-associated ROI",
    grepl("shared microglia[- ]neuropil microenvironment|shared local microenvironment", z) ~ "shared microglia-neuropil ROI",
    grepl("phagolysosom|phago|lysosom|complement|inflamm|immune activation|\\bdam\\b|dam[- ]like|disease[- ]associated microglia|\\bmhc\\b|major histocompatibility|antigen|interferon|cytokine|chemokine", z) ~ "microglia/immune",
    grepl("mitochond|respirat|oxidative|\\batp\\b|acetyl|tca|electron transport", z) ~ "mitochondrial",
    grepl("\\brna\\b|translation|ribosom|splic|mrna|ncrna|rnp", z) ~ "RNA/translation",
    grepl("synap|vesicle|postsynap|presynap|cytoskeleton|actin|microtubule", z) ~ "synaptic/cytoskeletal",
    grepl("ecm|extracellular matrix|basement membrane|collagen|laminin|nidogen|\\bagrn\\b|hspg2|perlecan|\\bcol4a[12]?\\b|\\blama[0-9]?\\b|\\blamb[0-9]?\\b|\\blamc[0-9]?\\b|\\bnid[12]\\b|\\bbcam\\b|serpinh1|perivascular", z) ~ "perivascular ECM",
    grepl("skin|keratin|epithelial|epiderm", z) ~ "epithelial/keratinocyte-like",
    grepl("stimulus|ion|cation|homeostasis|endosomal transport", z) ~ "stimulus/ion homeostasis",
    grepl("vascular|bbb|endothelial|pericyte|blood vessel", z) ~ "vascular/BBB",
    TRUE ~ "mixed"
  )
}

add_short_supermodule_labels <- function(df) {
  if (!nrow(df)) {
    df$Supermodule_FullAnnotationLabel <- character()
    df$Supermodule_DisplayShort <- character()
    df$Supermodule_CleanPlotLabel <- character()
    df$Supermodule_PlotLabel <- character()
    return(df)
  }
  id <- dplyr::coalesce(
    clean_label_value(col_or_na(df, "supermodule_id")),
    clean_label_value(col_or_na(df, "SupermoduleID")),
    clean_label_value(col_or_na(df, "supermodule_id_for_module")),
    clean_label_value(col_or_na(df, "module_supermodule_id")),
    rep("SM??", nrow(df))
  )
  full <- dplyr::coalesce(
    clean_supermodule_label_value(col_or_na(df, "Supermodule_CompositionDisplayLabel")),
    clean_supermodule_label_value(col_or_na(df, "Supermodule_CompositionLabel")),
    clean_supermodule_label_value(col_or_na(df, "cleaned_biological_label")),
    clean_supermodule_label_value(col_or_na(df, "Supermodule_FinalLabel")),
    clean_supermodule_label_value(col_or_na(df, "Supermodule_LongLabel")),
    clean_supermodule_label_value(col_or_na(df, "Macroprogram_Display")),
    clean_supermodule_label_value(col_or_na(df, "Supermodule_DisplayLabel")),
    clean_supermodule_label_value(col_or_na(df, "supermodule_label")),
    id
  )
  broad <- supermodule_broad_label(full)
  rank_df <- data.frame(supermodule_id = id, broad = broad, stringsAsFactors = FALSE) |>
    dplyr::distinct(.data$supermodule_id, .data$broad) |>
    dplyr::arrange(.data$broad, .data$supermodule_id) |>
    dplyr::group_by(.data$broad) |>
    dplyr::mutate(
      broad_n = dplyr::n(),
      broad_rank = dplyr::row_number(),
      suffix = dplyr::if_else(.data$broad_n > 1L, paste0(" ", roman_label(.data$broad_rank)), "")
    ) |>
    dplyr::ungroup()
  suffix <- rank_df$suffix[match(id, rank_df$supermodule_id)]
  clean_display <- dplyr::coalesce(
    clean_supermodule_label_value(col_or_na(df, "Supermodule_CompositionDisplayLabel")),
    clean_supermodule_label_value(col_or_na(df, "Supermodule_CompositionLabel")),
    clean_supermodule_label_value(col_or_na(df, "Supermodule_CleanPlotLabel"))
  )
  display <- dplyr::coalesce(add_supermodule_id_prefix(id, clean_display), paste0(id, supermodule_label_sep(), broad, suffix))
  df$Supermodule_FullAnnotationLabel <- full
  df$Supermodule_DisplayShort <- paste0(broad, suffix)
  df$Supermodule_CleanPlotLabel <- display
  df$Supermodule_PlotLabel <- display
  df
}

effect_sentence <- function(row, level = "supermodule") {
  row_get <- function(name, default = NA) {
    if (name %in% names(row)) row[[name]][[1]] else default
  }

  ds_label <- dataset_label(row$dataset)
  id <- if (level == "supermodule") row_get("supermodule_id") else row_get("module_id")
  label <- if (level == "supermodule") {
    row_get("Supermodule_DisplayLabel", row_get("Supermodule_FinalLabel", row_get("Macroprogram_Display", row_get("supermodule_label", row_get("supermodule_id")))))
  } else {
    row_get("module_label", row_get("ModuleID", row_get("module_id")))
  }
  cls <- row_get("dominant_microenvironment_class", row_get("microenvironment_class"))
  estimate <- safe_num(row_get("estimate"))
  direction <- ifelse(is.na(estimate), "is altered", ifelse(estimate > 0, "is higher", "is lower"))
  spatial_unit <- row_get("spatial_unit")
  spatial <- ifelse(
    is.na(spatial_unit) || spatial_unit == "global_spatial_adjusted",
    "after spatial adjustment",
    paste0("in ", spatial_unit)
  )
  effect_scope <- row_get("effect_scope")
  scope <- ifelse(!is.na(effect_scope), paste0(" (", effect_scope, ")"), "")
  base <- paste(ds_label, level, id, "annotated as", label, direction, "for", row_get("contrast"), spatial, scope)
  base <- gsub("[[:space:]]+\\(", " (", base)

  if (row$dataset != "microglia" || is.na(cls)) return(paste0(base, "."))
  if (cls == "microglia_supported") return(paste0(base, "; this ROI signal has microglia-supported ROI evidence."))
  if (cls == "microglia_state_or_activation_supported") return(paste0(base, "; this ROI signal has microglia state/phagolysosomal support."))
  if (cls == "shared_microenvironment") return(paste0(base, "; this ROI signal is shared local microenvironment signal, not purified microglial regulation."))
  if (cls == "neuropil_sensitive") return(paste0(base, "; this ROI signal is neuropil-sensitive; do not interpret as purified microglial regulation."))
  if (cls == "other_cellular_or_vascular_sensitive") return(paste0(base, "; this ROI signal shows other cellular/vascular marker support; interpret cautiously."))
  paste0(base, "; microenvironment support is ambiguous.")
}

add_interpretation_sentences <- function(df, level) {
  if (!nrow(df)) {
    df$interpretation_sentence <- character()
    return(df)
  }
  df$interpretation_sentence <- vapply(seq_len(nrow(df)), function(i) {
    effect_sentence(df[i, , drop = FALSE], level)
  }, character(1))
  df
}

theme_clean <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black", size = base_size),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = grid::unit(1.2, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain", color = "black", size = base_size * 0.9),
      legend.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 1, margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = base_size * 0.9, margin = ggplot2::margin(b = 4)),
      plot.margin = ggplot2::margin(4, 5, 4, 4),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(3.5, "mm")
    )
}

contrast_plot_levels <- function() c("RES - CON", "SUS - CON", "SUS - RES")

orient_contrasts_for_plot <- function(df) {
  if (!nrow(df) || !"contrast" %in% names(df)) return(df)
  df$contrast[!is.na(df$contrast) & df$contrast == "RES - SUS"] <- "SUS - RES"
  df
}

region_plot_levels <- function() c("DG", "CA3", "CA2", "CA1")

neuropil_region_layer_levels <- function() {
  c(
    "CA1_so", "CA1_sr", "CA1_slm",
    "CA2_so", "CA2_sr", "CA2_slm",
    "CA3_so", "CA3_sr",
    "DG_mo", "DG_po"
  )
}

spatial_unit_plot_label <- function(x, ds = NULL) {
  x_raw <- as.character(x)

  if (identical(ds, "neuron_neuropil")) {
    out <- toupper(x_raw)
    out <- sub("_", "_", out, fixed = TRUE)
    out <- gsub("^CA([123])_", "CA\\1_", out)
    out <- gsub("^DG_", "DG_", out)
    return(out)
  }

  x_low <- tolower(x_raw)
  dplyr::case_when(
    grepl("^dg", x_low) ~ "DG",
    grepl("^ca3", x_low) ~ "CA3",
    grepl("^ca2", x_low) ~ "CA2",
    grepl("^ca1", x_low) ~ "CA1",
    TRUE ~ x_raw
  )
}

apply_spatial_unit_order <- function(df, ds = NULL) {
  if (!nrow(df) || !"spatial_unit" %in% names(df)) return(df)
  raw <- if ("spatial_unit_raw" %in% names(df)) as.character(df$spatial_unit_raw) else as.character(df$spatial_unit)
  labels <- spatial_unit_plot_label(raw, ds = ds)

  spatial_levels <- if (identical(ds, "neuron_neuropil")) {
    c(
      "DG_MO", "DG_PO",
      "CA3_SO", "CA3_SR",
      "CA2_SO", "CA2_SR", "CA2_SLM",
      "CA1_SO", "CA1_SR", "CA1_SLM"
    )
  } else {
    c(region_plot_levels(), sort(setdiff(unique(labels), region_plot_levels())))
  }

  spatial_levels <- c(spatial_levels, sort(setdiff(unique(labels), spatial_levels)))
  df$spatial_unit_raw <- raw
  df$spatial_unit <- factor(labels, levels = spatial_levels)
  df
}

effect_limits <- function(x) {
  z <- abs(safe_num(x))
  lim <- suppressWarnings(max(z[is.finite(z)], na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) lim <- 1
  c(-lim, lim)
}

scale_effect_fill <- function(limits = NULL, name = "Estimate") {
  ggplot2::scale_fill_gradient2(
    low = "#3B6EA8", mid = "#F7F7F7", high = "#B24A4A",
    midpoint = 0, limits = limits, oob = scales::squish,
    name = name, na.value = "grey92"
  )
}

q_value_col <- function(df) {
  n <- nrow(df)
  dplyr::coalesce(
    safe_num(col_or_na(df, "FDR_global")),
    safe_num(col_or_na(df, "FDR_within_dataset_level")),
    safe_num(col_or_na(df, "q_value")),
    rep(NA_real_, n)
  )
}

add_plot_metrics <- function(df, fdr_col = "FDR_global", p_col = "p_value") {
  if (!fdr_col %in% names(df)) df[[fdr_col]] <- NA_real_
  if (!p_col %in% names(df)) df[[p_col]] <- NA_real_
  if (!"estimate" %in% names(df)) df$estimate <- NA_real_
  df$q_value <- q_value_col(df)

  df |>
    dplyr::mutate(
      estimate = safe_num(.data$estimate),
      p_value = safe_num(.data[[p_col]]),
      FDR_global = safe_num(.data[[fdr_col]]),
      neg_log10_FDR = -log10(clip_p(.data$FDR_global)),
      neg_log10_P = -log10(clip_p(.data$p_value)),
      neg_log10_q = dplyr::if_else(is.na(.data$q_value), NA_real_, -log10(clip_p(.data$q_value))),
      sig_label = dplyr::case_when(
        !is.na(.data$q_value) & .data$q_value < 0.01 ~ "***",
        !is.na(.data$q_value) & .data$q_value < 0.05 ~ "**",
        !is.na(.data$q_value) & .data$q_value < 0.10 ~ "*",
        TRUE ~ ""
      ),
      evidence_rank = dplyr::case_when(
        .data$evidence_status == "robust_FDR" ~ 1L,
        .data$evidence_status == "suggestive_FDR10" ~ 2L,
        .data$evidence_status == "nominal_only" ~ 3L,
        .data$evidence_status == "model_unstable" ~ 4L,
        TRUE ~ 5L
      )
    )
}

main_effect_rows <- function(df) {
  if (!nrow(df)) return(df)
  out <- df |> dplyr::filter(!is.na(.data$p_value))
  if ("spatial_unit" %in% names(out) && any(out$spatial_unit == "global_spatial_adjusted", na.rm = TRUE)) {
    out <- out |> dplyr::filter(.data$spatial_unit == "global_spatial_adjusted")
  }
  out
}

spatial_effect_rows <- function(df) {
  if (!nrow(df) || !"spatial_unit" %in% names(df)) return(df[0, , drop = FALSE])
  df |>
    dplyr::filter(!is.na(.data$p_value), !is.na(.data$spatial_unit), .data$spatial_unit != "global_spatial_adjusted")
}

program_label_col <- function(df, level = "supermodule") {
  if (level == "supermodule") {
    coalesce_chr(
      col_or_na(df, "Supermodule_PlotLabel"),
      col_or_na(df, "Supermodule_CleanPlotLabel"),
      col_or_na(df, "Supermodule_CompositionDisplayLabel"),
      col_or_na(df, "Supermodule_CompositionLabel"),
      col_or_na(df, "Supermodule_DisplayLabel"),
      col_or_na(df, "Supermodule_FinalLabel"),
      col_or_na(df, "Supermodule_FinalLabel.y"),
      col_or_na(df, "Supermodule_FinalLabel.x"),
      col_or_na(df, "Macroprogram_Display"),
      col_or_na(df, "supermodule_label_for_module"),
      col_or_na(df, "module_supermodule_label"),
      col_or_na(df, "supermodule_label"),
      col_or_na(df, "supermodule_label.y"),
      col_or_na(df, "supermodule_label.x"),
      col_or_na(df, "SupermoduleID"),
      col_or_na(df, "supermodule_id_for_module"),
      col_or_na(df, "module_supermodule_id"),
      col_or_na(df, "supermodule_id")
    )
  } else {
    # Display labels should prefer the annotation table from 06_annotate_module_microenvironment.r.
    # The group-effect table from 05 usually only knows the raw WGCNA module label/color.
    coalesce_chr(
      col_or_na(df, "ModulePlotLabel"),
      col_or_na(df, "Module_CleanPlotLabel"),
      col_or_na(df, "cleaned_biological_label"),
      col_or_na(df, "module_display_label"),
      col_or_na(df, "ModuleLabel_Final"),
      col_or_na(df, "module_label.y"),
      col_or_na(df, "module_label"),
      col_or_na(df, "endpoint_label"),
      col_or_na(df, "ModuleID"),
      col_or_na(df, "module_id"),
      col_or_na(df, "ModuleColor")
    )
  }
}

module_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\s*[â€”|-].*$", "", x)  # strip display suffix if present
  x <- sub("^ME", "", x, ignore.case = FALSE)
  x <- sub("^WGCNA[_-]?", "", x, ignore.case = TRUE)
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

module_join_key <- function(df) {
  # Raw IDs/colors only. Do not use module_display_label here because it may contain biology text.
  candidates <- c(
    "module_id", "ModuleID", "ModuleColor", "module_eigengene",
    "endpoint_id", "module_label.x", "module_label", "module_label.y"
  )
  vals <- rep(NA_character_, nrow(df))
  for (nm in candidates) {
    if (nm %in% names(df)) vals <- dplyr::coalesce(vals, as.character(df[[nm]]))
  }
  module_key(vals)
}

pretty_module_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^WGCNA[_-]?", "", x, ignore.case = TRUE)
  x <- sub("^ME", "", x, ignore.case = FALSE)
  x[is.na(x) | !nzchar(x)] <- "unlabelled"
  x
}

pretty_program_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | !nzchar(x)] <- "Unassigned"
  x
}

candidate_cols <- function(df, candidates) intersect(candidates, names(df))

build_module_supermodule_map <- function(module_effects, module_to_supermodule_map = NULL, supermodule_composition = NULL, super_annot = NULL) {
  out <- data.frame(
    module_key = character(),
    supermodule_id = character(),
    supermodule_label = character(),
    map_source = character(),
    stringsAsFactors = FALSE
  )

  add_map_rows <- function(keys, sid, slabel, source) {
    keys <- unique(module_key(keys))
    keys <- keys[nzchar(keys) & !is.na(keys)]
    if (!length(keys)) return(NULL)
    data.frame(
      module_key = keys,
      supermodule_id = as.character(sid %||% NA_character_),
      supermodule_label = as.character(slabel %||% sid %||% NA_character_),
      map_source = source,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(module_to_supermodule_map) && nrow(module_to_supermodule_map)) {
    map_df <- as.data.frame(module_to_supermodule_map, stringsAsFactors = FALSE)

    module_cols <- candidate_cols(map_df, c(
      "ModuleID", "module_id", "module", "Module", "ModuleColor", "module_color",
      "module_eigengene", "moduleEigengene", "module_eigengene_col"
    ))
    sid_cols <- candidate_cols(map_df, c(
      "Supermodule_DataDrivenID", "Supermodule_DataDriven", "SupermoduleID",
      "supermodule_id", "Supermodule", "Supermodule_DataDrivenLabel"
    ))
    slabel_cols <- candidate_cols(map_df, c(
      "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "SupermoduleLabel", "Supermodule_DataDrivenLabel",
      "supermodule_label", "Supermodule", "SupermoduleID"
    ))

    if (length(module_cols) && length(sid_cols)) {
      rows <- lapply(seq_len(nrow(map_df)), function(i) {
        keys <- unlist(map_df[i, module_cols, drop = TRUE], use.names = FALSE)
        sid_vals <- unlist(map_df[i, sid_cols, drop = TRUE], use.names = FALSE)
        sid_vals <- as.character(sid_vals[!is.na(sid_vals) & nzchar(as.character(sid_vals))])
        sid <- if (length(sid_vals)) sid_vals[[1]] else NA_character_
        if (length(slabel_cols)) {
          lab_vals <- unlist(map_df[i, slabel_cols, drop = TRUE], use.names = FALSE)
          lab_vals <- as.character(lab_vals[!is.na(lab_vals) & nzchar(as.character(lab_vals))])
          slabel <- if (length(lab_vals)) lab_vals[[1]] else sid
        } else {
          slabel <- sid
        }
        add_map_rows(keys, sid, slabel, "module_to_supermodule_map")
      })
      out <- dplyr::bind_rows(out, dplyr::bind_rows(rows))
    }
  }

  if (!is.null(supermodule_composition) && nrow(supermodule_composition) && "member_modules" %in% names(supermodule_composition)) {
    comp <- as.data.frame(supermodule_composition, stringsAsFactors = FALSE)
    for (nm in c("supermodule_id", "supermodule_label", "SupermoduleID", "SupermoduleLabel")) if (!nm %in% names(comp)) comp[[nm]] <- NA_character_
    rows <- lapply(seq_len(nrow(comp)), function(i) {
      members <- unlist(strsplit(as.character(comp$member_modules[[i]]), "[;|, ]+"), use.names = FALSE)
      sid <- dplyr::coalesce(as.character(comp$supermodule_id[[i]]), as.character(comp$SupermoduleID[[i]]))
      slabel <- dplyr::coalesce(as.character(comp$supermodule_label[[i]]), as.character(comp$SupermoduleLabel[[i]]), sid)
      add_map_rows(members, sid, slabel, "supermodule_composition")
    })
    out <- dplyr::bind_rows(out, dplyr::bind_rows(rows))
  }

  if (!nrow(out)) {
    return(data.frame(module_key = character(), module_supermodule_id = character(), module_supermodule_label = character(), module_supermodule_map_source = character()))
  }

  if (!is.null(super_annot) && nrow(super_annot)) {
    ann <- as.data.frame(super_annot, stringsAsFactors = FALSE)
    sid_cols <- candidate_cols(ann, c("SupermoduleID", "supermodule_id", "Supermodule_DataDrivenID", "Supermodule_DataDriven", "Supermodule"))
    label_cols <- candidate_cols(ann, c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "SupermoduleLabel", "Supermodule_DataDrivenLabel", "supermodule_label"))
    if (length(sid_cols) && length(label_cols)) {
      ann2 <- dplyr::bind_rows(lapply(seq_len(nrow(ann)), function(i) {
        sid_vals <- unlist(ann[i, sid_cols, drop = TRUE], use.names = FALSE)
        lab_vals <- unlist(ann[i, label_cols, drop = TRUE], use.names = FALSE)
        sid_vals <- as.character(sid_vals[!is.na(sid_vals) & nzchar(as.character(sid_vals))])
        lab_vals <- as.character(lab_vals[!is.na(lab_vals) & nzchar(as.character(lab_vals))])
        data.frame(
          supermodule_id = if (length(sid_vals)) sid_vals[[1]] else NA_character_,
          curated_supermodule_label = if (length(lab_vals)) lab_vals[[1]] else NA_character_,
          stringsAsFactors = FALSE
        )
      })) |>
        dplyr::filter(!is.na(.data$supermodule_id), nzchar(.data$supermodule_id)) |>
        dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)
      out <- out |>
        dplyr::left_join(ann2, by = "supermodule_id") |>
        dplyr::mutate(supermodule_label = dplyr::coalesce(.data$curated_supermodule_label, .data$supermodule_label)) |>
        dplyr::select(-dplyr::any_of("curated_supermodule_label"))
    }
  }

  out |>
    dplyr::filter(!is.na(.data$module_key), nzchar(.data$module_key)) |>
    dplyr::filter(!is.na(.data$supermodule_id), nzchar(.data$supermodule_id)) |>
    dplyr::distinct(.data$module_key, .keep_all = TRUE) |>
    dplyr::rename(
      module_supermodule_id = "supermodule_id",
      module_supermodule_label = "supermodule_label",
      module_supermodule_map_source = "map_source"
    )
}

attach_module_supermodules <- function(module_join, module_super_map) {
  if (!nrow(module_join)) return(module_join)
  module_join$module_key <- module_join_key(module_join)
  if (!is.null(module_super_map) && nrow(module_super_map)) {
    module_join <- module_join |>
      dplyr::left_join(module_super_map, by = "module_key")
    module_join$module_supermodule_id <- dplyr::coalesce(
      clean_label_value(col_or_na(module_join, "module_supermodule_id")),
      clean_label_value(col_or_na(module_join, "module_supermodule_id.y")),
      clean_label_value(col_or_na(module_join, "module_supermodule_id.x"))
    )
    module_join$module_supermodule_label <- dplyr::coalesce(
      clean_label_value(col_or_na(module_join, "module_supermodule_label")),
      clean_label_value(col_or_na(module_join, "module_supermodule_label.y")),
      clean_label_value(col_or_na(module_join, "module_supermodule_label.x"))
    )
    module_join$module_supermodule_map_source <- dplyr::coalesce(
      clean_label_value(col_or_na(module_join, "module_supermodule_map_source")),
      clean_label_value(col_or_na(module_join, "module_supermodule_map_source.y")),
      clean_label_value(col_or_na(module_join, "module_supermodule_map_source.x"))
    )
  } else {
    module_join$module_supermodule_id <- NA_character_
    module_join$module_supermodule_label <- NA_character_
    module_join$module_supermodule_map_source <- NA_character_
  }
  module_join |>
    dplyr::mutate(
      supermodule_id_for_module = dplyr::coalesce(.data$module_supermodule_id, .data$supermodule_id),
      supermodule_label_for_module = dplyr::coalesce(.data$module_supermodule_label, .data$supermodule_label, "Unassigned"),
      supermodule_map_source_for_module = dplyr::coalesce(.data$module_supermodule_map_source, "not_mapped")
    )
}

save_svg <- function(plot, filename, width = 170, height = 110) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot, width = width, height = height, units = "mm", device = svglite::svglite, limitsize = FALSE)
  pdf_file <- sub("\\.svg$", ".pdf", filename)
  tryCatch(
    ggplot2::ggsave(pdf_file, plot, width = width, height = height, units = "mm", device = grDevices::cairo_pdf, limitsize = FALSE),
    error = function(e) warning("PDF export skipped for ", basename(filename), ": ", conditionMessage(e), call. = FALSE)
  )
}

plot_supermodule_main_heatmap <- function(super_join, paths, ds) {
  plot_df <- main_effect_rows(super_join) |>
    orient_contrasts_for_plot() |>
    add_plot_metrics()

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df$program_label <- label_wrap(coalesce_chr(col_or_na(plot_df, "Supermodule_PlotLabel"), program_label_col(plot_df, "supermodule")), 34)
  plot_df <- plot_df |>
    dplyr::mutate(
      contrast = factor(.data$contrast, levels = contrast_plot_levels()),
      program_label = stats::reorder(.data$program_label, .data$estimate, FUN = function(z) max(abs(z), na.rm = TRUE))
    )
  fill_limits <- effect_limits(plot_df$estimate)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data$program_label, fill = .data$estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.18) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 2.0, color = "black") +
    scale_effect_fill(fill_limits) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": supermodule effects"),
      x = NULL, y = NULL
    ) +
    theme_clean(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "supermodule_group_effects_main_heatmap.svg"), width = 90, height = 90)
}

plot_supermodule_spatial_heatmap <- function(super_join, paths, ds) {
  plot_df <- spatial_effect_rows(super_join) |>
    orient_contrasts_for_plot() |>
    add_plot_metrics()

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df$program_label <- label_wrap(coalesce_chr(col_or_na(plot_df, "Supermodule_PlotLabel"), program_label_col(plot_df, "supermodule")), 30)
  plot_df <- plot_df |>
    apply_spatial_unit_order(ds = ds) |>
    dplyr::mutate(
      contrast = factor(.data$contrast, levels = contrast_plot_levels()),
      panel = paste(.data$contrast, .data$spatial_unit, sep = " | ")
    )

  # Keep the spatial plot readable by showing strongest rows when the matrix is too large.
  max_rows <- 45L
  keep_labels <- plot_df |>
    dplyr::group_by(.data$program_label) |>
    dplyr::summarise(rank_metric = min(.data$FDR_global, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$rank_metric) |>
    dplyr::slice_head(n = max_rows) |>
    dplyr::pull(.data$program_label)

  plot_df <- plot_df |>
    dplyr::filter(.data$program_label %in% keep_labels) |>
    dplyr::arrange(.data$contrast, .data$spatial_unit, .data$program_label)

  fill_limits <- effect_limits(plot_df$estimate)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$spatial_unit, y = .data$program_label, fill = .data$estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.16) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 1.8, color = "black") +
    ggplot2::facet_wrap(~ contrast, nrow = 1) +
    scale_effect_fill(fill_limits) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": spatial supermodule effects"),
      x = spatial_unit_label(ds), y = NULL
    ) +
    theme_clean(base_size = 7.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "supermodule_group_effects_spatial_heatmap.svg"), width = 160, height = 100)
}


plot_module_main_heatmap <- function(module_join, paths, ds) {
  plot_df <- main_effect_rows(module_join) |>
    orient_contrasts_for_plot() |>
    add_plot_metrics() |>
    dplyr::filter(!is.na(.data$p_value))

  if (!nrow(plot_df)) return(invisible(NULL))

  raw_id <- coalesce_chr(col_or_na(plot_df, "ModulePlotLabel"), make_module_plot_label(plot_df))
  ann_lab <- shorten_supermodule_label(coalesce_chr(col_or_na(plot_df, "microenvironment_label"), program_label_col(plot_df, "module")), max_chars = 28)
  plot_df$module_label_plot <- paste0(raw_id, " | ", ann_lab)
  plot_df$module_label_plot <- label_wrap(plot_df$module_label_plot, 26)
  plot_df$supermodule_plot <- label_wrap(shorten_supermodule_label(pretty_program_label(coalesce_chr(col_or_na(plot_df, "Supermodule_PlotLabel"), col_or_na(plot_df, "supermodule_label_for_module"))), max_chars = 45), 30)

  max_modules <- 45L
  keep_modules <- plot_df |>
    dplyr::group_by(.data$module_label_plot) |>
    dplyr::summarise(best_p = min(.data$p_value, na.rm = TRUE), best_abs_estimate = max(abs(.data$estimate), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$best_p, dplyr::desc(.data$best_abs_estimate)) |>
    dplyr::slice_head(n = max_modules) |>
    dplyr::pull(.data$module_label_plot)

  plot_df <- plot_df |>
    dplyr::filter(.data$module_label_plot %in% keep_modules) |>
    dplyr::mutate(
      contrast = factor(.data$contrast, levels = contrast_plot_levels()),
      module_label_plot = stats::reorder(.data$module_label_plot, .data$estimate, FUN = function(z) max(abs(z), na.rm = TRUE))
    )
  write_csv_safe2(plot_df, file.path(paths$source_data, "module_group_effects_main_heatmap_source.csv"))
  fig_h <- min(180, max(130, 34 + 4.4 * length(unique(plot_df$module_label_plot)) + 3.0 * length(unique(plot_df$supermodule_plot))))
  fill_limits <- effect_limits(plot_df$estimate)

  p_heat <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data$module_label_plot, fill = .data$estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.16) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 1.7, color = "black") +
    ggplot2::facet_grid(supermodule_plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_effect_fill(fill_limits) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": module effects"),
      x = NULL, y = NULL
    ) +
    theme_clean(base_size = 6.8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 5.8),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
      legend.position = "right"
    )

  save_svg(p_heat, file.path(paths$figures, "module_group_effects_main_heatmap.svg"), width = 175, height = fig_h)
}

plot_module_spatial_heatmap <- function(module_join, paths, ds) {
  plot_df <- spatial_effect_rows(module_join) |>
    orient_contrasts_for_plot() |>
    add_plot_metrics() |>
    dplyr::filter(!is.na(.data$p_value))

  if (!nrow(plot_df)) return(invisible(NULL))

  raw_id <- coalesce_chr(col_or_na(plot_df, "ModulePlotLabel"), make_module_plot_label(plot_df))
  ann_lab <- shorten_supermodule_label(coalesce_chr(col_or_na(plot_df, "microenvironment_label"), program_label_col(plot_df, "module")), max_chars = 26)
  plot_df$module_label_plot <- label_wrap(paste0(raw_id, " | ", ann_lab), 25)
  plot_df$supermodule_plot <- label_wrap(shorten_supermodule_label(pretty_program_label(coalesce_chr(col_or_na(plot_df, "Supermodule_PlotLabel"), col_or_na(plot_df, "supermodule_label_for_module"))), max_chars = 45), 30)

  max_modules <- 35L
  keep_modules <- plot_df |>
    dplyr::group_by(.data$module_label_plot) |>
    dplyr::summarise(best_p = min(.data$p_value, na.rm = TRUE), best_abs_estimate = max(abs(.data$estimate), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$best_p, dplyr::desc(.data$best_abs_estimate)) |>
    dplyr::slice_head(n = max_modules) |>
    dplyr::pull(.data$module_label_plot)

  plot_df <- plot_df |>
    apply_spatial_unit_order(ds = ds) |>
    dplyr::filter(.data$module_label_plot %in% keep_modules) |>
    dplyr::mutate(
      contrast = factor(.data$contrast, levels = contrast_plot_levels())
    ) |>
    dplyr::arrange(.data$contrast, .data$spatial_unit, .data$supermodule_plot, .data$module_label_plot)
  write_csv_safe2(plot_df, file.path(paths$source_data, "module_group_effects_spatial_heatmap_source.csv"))
  fig_h <- min(180, max(130, 36 + 4.6 * length(unique(plot_df$module_label_plot)) + 3.2 * length(unique(plot_df$supermodule_plot))))
  fill_limits <- effect_limits(plot_df$estimate)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$spatial_unit, y = .data$module_label_plot, fill = .data$estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.12) +
    ggplot2::facet_grid(supermodule_plot ~ contrast, scales = "free_y", space = "free_y", switch = "y") +
    scale_effect_fill(fill_limits) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": spatial module effects"),
      x = spatial_unit_label(ds), y = NULL
    ) +
    theme_clean(base_size = 6.2) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 5.2),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "module_group_effects_spatial_heatmap.svg"), width = 180, height = fig_h)
}

plot_supermodule_membership_overview <- function(module_join, super_join, paths, ds) {
  if (!nrow(module_join)) return(invisible(NULL))

  df <- module_join |>
    dplyr::distinct(.data$module_id, .data$supermodule_id_for_module, .data$supermodule_label_for_module, .keep_all = TRUE)

  df$supermodule_plot <- label_wrap(pretty_program_label(coalesce_chr(col_or_na(df, "Supermodule_PlotLabel"), col_or_na(df, "supermodule_label_for_module"))), 30)
  protein_col <- first_existing_col(df, c("n_proteins", "nGenes", "module_size", "ModuleSize"), fallback = NULL)

  overview <- df |>
    dplyr::group_by(.data$supermodule_plot) |>
    dplyr::summarise(
      n_modules = dplyr::n_distinct(.data$module_id),
      n_proteins = if (!is.null(protein_col)) sum(safe_num(.data[[protein_col]]), na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      membership_type = dplyr::if_else(.data$n_modules <= 1L, "singleton module group", "multi-module group"),
      supermodule_plot_display = dplyr::if_else(.data$n_modules <= 1L, paste0(.data$supermodule_plot, "\n(singleton)"), .data$supermodule_plot)
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_modules))

  if (!nrow(overview)) return(invisible(NULL))
  write_csv_safe2(overview, file.path(paths$source_data, "supermodule_membership_overview_source.csv"))

  p <- ggplot2::ggplot(overview, ggplot2::aes(x = .data$n_modules, y = stats::reorder(.data$supermodule_plot_display, .data$n_modules), fill = .data$membership_type)) +
    ggplot2::geom_col(width = 0.64, color = "white", linewidth = 0.15) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_modules), hjust = -0.18, size = 2.5) +
    ggplot2::scale_fill_manual(values = c("multi-module group" = "grey45", "singleton module group" = "grey75"), drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": module groups"),
      x = "WGCNA modules", y = NULL, fill = NULL
    ) +
    theme_clean(base_size = 7) +
    ggplot2::theme(legend.position = "bottom")

  save_svg(p, file.path(paths$figures, "supermodule_membership_overview.svg"), width = 120, height = 95)
}

plot_top_supermodules <- function(top_super, paths, ds, n_modules = 10) {

  plot_df <- top_super |>
    orient_contrasts_for_plot() |>
    add_plot_metrics() |>
    dplyr::filter(
      !is.na(.data$p_value),
      !is.na(.data$estimate)
    )

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df <- plot_df |>
    dplyr::mutate(
      program_label_raw = coalesce_chr(
        col_or_na(plot_df, "Supermodule_PlotLabel"),
        program_label_col(plot_df, "supermodule")
      ),
      contrast = factor(.data$contrast, levels = contrast_plot_levels())
    )

  top_modules <- plot_df |>
    dplyr::group_by(.data$program_label_raw) |>
    dplyr::summarise(
      max_abs_estimate = max(abs(.data$estimate), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$max_abs_estimate)) |>
    dplyr::slice_head(n = n_modules)

  plot_df <- plot_df |>
    dplyr::semi_join(top_modules, by = "program_label_raw") |>
    dplyr::left_join(top_modules, by = "program_label_raw") |>
    dplyr::mutate(
      program_label = label_wrap(.data$program_label_raw, 34),
      program_label = stats::reorder(.data$program_label, .data$max_abs_estimate),
      neg_log10_q_plot = pmin(.data$neg_log10_q, 2.5)
    )

  x_max <- max(abs(plot_df$estimate), na.rm = TRUE)
  x_lim <- c(-1, 1) * x_max * 1.15

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$estimate,
      y = .data$program_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.25,
      color = "grey60"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = .data$contrast,
        size = .data$neg_log10_q_plot
      ),
      shape = 21,
      color = "white",
      stroke = 0.20,
      alpha = 0.7,
      position = ggplot2::position_jitter(
        width = 0,
        height = 0.08,
        seed = 1
      )
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "RES - CON" = "#3E3C6F",
        "SUS - CON" = "#E63A48",
        "SUS - RES" = "#8D8982"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
      range = c(1.8, 4.8),
      breaks = c(0.5, 1, 1.5, 2, 2.5),
      limits = c(0, 2.5),
      name = expression(-log[10]("FDR"))
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        override.aes = list(
          fill = "grey45",
          color = "grey20",
          alpha = 1,
          stroke = 0.25
        )
      ),
      fill = ggplot2::guide_legend(
        override.aes = list(
          size = 3,
          color = "grey20",
          alpha = 1,
          stroke = 0.25
        )
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = x_lim,
      expand = ggplot2::expansion(mult = c(0.03, 0.06))
    ) +
    ggplot2::labs(
      title = NULL,
      x = "Model estimate",
      y = NULL,
      fill = "Contrast"
    ) +
    theme_clean(base_size = 7) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial"),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 6.5),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(3.2, "mm"),
      axis.text.y = ggplot2::element_text(size = 6.2),
      axis.text.x = ggplot2::element_text(size = 6.2),
      axis.title.x = ggplot2::element_text(size = 7),
      #panel.grid.major.y = ggplot2::element_line(
      #  linewidth = 0.15,
      #  color = "grey90"
      #),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 3, 2, 2)
    )

  save_svg(
    p,
    file.path(paths$figures, "top_supermodule_effects_dotplot.svg"),
    width = 90,
    height = 62
  )

  invisible(p)
}

plot_module_effects_by_supermodule <- function(module_join, paths, ds) {
  plot_df <- module_join |>
    dplyr::filter(!is.na(.data$p_value), !is.na(.data$Supermodule_PlotLabel)) |>
    orient_contrasts_for_plot() |>
    add_plot_metrics()

  if (!nrow(plot_df)) return(invisible(NULL))

  raw_id <- coalesce_chr(col_or_na(plot_df, "ModulePlotLabel"), make_module_plot_label(plot_df))
  ann_lab <- shorten_supermodule_label(coalesce_chr(col_or_na(plot_df, "microenvironment_label"), program_label_col(plot_df, "module")), max_chars = 24)
  plot_df$module_label_plot <- label_wrap(paste0(raw_id, " | ", ann_lab), 23)
  plot_df$supermodule_label_plot <- label_wrap(coalesce_chr(col_or_na(plot_df, "Supermodule_PlotLabel"), col_or_na(plot_df, "supermodule_label")), 24)

  # This is a supplementary diagnostic figure. Limit each facet to the strongest module rows.
  plot_df <- plot_df |>
    dplyr::group_by(.data$supermodule_label_plot, .data$module_id) |>
    dplyr::mutate(module_rank = min(.data$p_value, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$supermodule_label_plot) |>
    dplyr::arrange(.data$module_rank, .by_group = TRUE) |>
    dplyr::filter(dplyr::row_number() <= 30L) |>
    dplyr::ungroup()

  plot_df <- plot_df |>
    dplyr::mutate(contrast = factor(.data$contrast, levels = contrast_plot_levels()))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$estimate, y = .data$module_label_plot, fill = .data$contrast, size = .data$neg_log10_q)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.22, color = "grey70") +
    ggplot2::geom_point(shape = 21, color = "black", stroke = 0.12, alpha = 0.86) +
    ggplot2::facet_grid(supermodule_label_plot ~ contrast, scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_manual(values = c("RES - CON" = "#4D9221", "SUS - CON" = "#C51B7D", "SUS - RES" = "#2166AC"), drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(0.8, 2.8), name = "-log10 q") +
    ggplot2::labs(
      title = paste0(dataset_label(ds), ": module effects by supermodule"),
      x = "Model estimate", y = NULL, size = "-log10 q", fill = "Contrast"
    ) +
    theme_clean(base_size = 6.3) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 5.2),
      strip.placement = "outside",
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "module_effects_by_supermodule_dotplot.svg"), width = 170, height = 140)
}

plot_microglia_effect_heatmap <- function(super_join, paths) {
  plot_df <- main_effect_rows(super_join) |>
    orient_contrasts_for_plot() |>
    add_plot_metrics()

  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df$program_label <- label_wrap(coalesce_chr(col_or_na(plot_df, "Supermodule_PlotLabel"), program_label_col(plot_df, "supermodule")), 34)
  if (!"dominant_microenvironment_class" %in% names(plot_df)) plot_df$dominant_microenvironment_class <- NA_character_
  plot_df <- plot_df |>
    dplyr::mutate(
      microenvironment_class = dplyr::coalesce(.data$dominant_microenvironment_class, "unclassified"),
      contrast = factor(.data$contrast, levels = contrast_plot_levels())
    )
  fill_limits <- effect_limits(plot_df$estimate)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data$program_label, fill = .data$estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.16) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 1.9, color = "black") +
    ggplot2::facet_grid(microenvironment_class ~ ., scales = "free_y", space = "free_y") +
    scale_effect_fill(fill_limits) +
    ggplot2::labs(
      title = "Microglia ROI: supermodule effects",
      x = NULL, y = NULL
    ) +
    theme_clean(base_size = 6.8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 5.8),
      legend.position = "right"
    )

  save_svg(p, file.path(paths$figures, "microglia_supermodule_effects_by_microenvironment.svg"), width = 120, height = 95)
}

plot_microglia_composition_heatmap <- function(super_annot, paths) {
  frac_cols <- grep("^fraction_modules_", names(super_annot), value = TRUE)
  if (!length(frac_cols) || !nrow(super_annot)) return(invisible(NULL))

  plot_df <- super_annot |>
    dplyr::mutate(
      program_label = make_supermodule_plot_label(super_annot),
      program_label = label_wrap(shorten_supermodule_label(.data$program_label, max_chars = 45), 34)
    ) |>
    dplyr::select(dplyr::any_of(c("SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_LongLabel", "Macroprogram_Display")), "program_label", dplyr::all_of(frac_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(frac_cols), names_to = "class", values_to = "fraction") |>
    dplyr::mutate(
      fraction = safe_num(.data$fraction),
      class = gsub("^fraction_modules_", "", .data$class),
      class = gsub("_", " ", .data$class)
    ) |>
    dplyr::filter(is.finite(.data$fraction), .data$fraction > 0)

  if (!nrow(plot_df)) return(invisible(NULL))
  write_csv_safe2(plot_df, file.path(paths$source_data, "microglia_supermodule_annotation_composition_source.csv"))
  fig_h <- max(92, 9.5 * length(unique(plot_df$program_label)) + 30)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$fraction, y = .data$program_label, fill = .data$class)) +
    ggplot2::geom_col(width = 0.66, color = "white", linewidth = 0.12) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::scale_fill_manual(
      values = c(
        ambiguous = "#BDBDBD",
        "microglia supported" = "#C51B7D",
        "microglia state or activation supported" = "#E6550D",
        "neuropil sensitive" = "#3182BD",
        "other cellular or vascular sensitive" = "#756BB1",
        "shared microenvironment" = "#238B45"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Microglia ROI: supermodule annotation composition",
      subtitle = "Fractions summarize the module classes contributing to each supermodule",
      x = "Fraction of member modules", y = NULL, fill = NULL
    ) +
    theme_clean(base_size = 6.7) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box.spacing = grid::unit(1, "mm"),
      axis.text.y = ggplot2::element_text(size = 5.8)
    )

  save_svg(p, file.path(paths$figures, "microglia_supermodule_annotation_composition.svg"), width = 166, height = fig_h)
}

plot_microglia_marker_evidence <- function(module_annot, paths) {
  if (!nrow(module_annot)) return(invisible(NULL))

  x_col <- first_existing_col(
    module_annot,
    c("empirical_neuropil_enriched_fraction", "canonical_neuronal_synaptic_neuropil_fraction", "neuropil_synaptic_neuronal_marker_fraction")
  )
  y_col <- first_existing_col(
    module_annot,
    c("empirical_microglia_roi_enriched_fraction", "canonical_microglia_homeostatic_fraction", "microglia_marker_fraction")
  )
  if (is.null(x_col) || is.null(y_col)) return(invisible(NULL))

  class_col <- first_existing_col(module_annot, c("microenvironment_class", "dominant_microenvironment_class"), fallback = NULL)

  plot_df <- module_annot |>
    dplyr::mutate(
      x_value = safe_num(.data[[x_col]]),
      y_value = safe_num(.data[[y_col]])
    )
  plot_df$module_label_plot <- coalesce_chr(col_or_na(plot_df, "ModulePlotLabel"), program_label_col(plot_df, "module"))

  if (!is.null(class_col)) {
    plot_df$marker_class <- as.character(plot_df[[class_col]])
  } else {
    plot_df$marker_class <- "unclassified"
  }

  plot_df <- plot_df |> dplyr::filter(!is.na(.data$x_value), !is.na(.data$y_value))
  if (!nrow(plot_df)) return(invisible(NULL))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x_value, y = .data$y_value)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey82", linewidth = 0.22) +
    ggplot2::geom_vline(xintercept = 0, color = "grey82", linewidth = 0.22) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.22, color = "grey55", linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$marker_class), shape = 21, color = "black", stroke = 0.14, size = 1.55, alpha = 0.82) +
    ggplot2::scale_fill_manual(
      values = c(
        ambiguous = "#737373",
        microglia_supported = "#C51B7D",
        microglia_state_or_activation_supported = "#E6550D",
        neuropil_sensitive = "#3182BD",
        other_cellular_or_vascular_sensitive = "#756BB1",
        shared_microenvironment = "#238B45",
        unclassified = "#969696"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Microglia ROI: module marker evidence",
      subtitle = "Numeric scatter; diagonal separates neuropil-skewed from microglia-skewed marker evidence",
      x = x_col,
      y = y_col,
      fill = "Class"
    ) +
    theme_clean(base_size = 6.8) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.18),
      panel.grid.minor = ggplot2::element_blank()
    )

  save_svg(p, file.path(paths$figures, "microglia_module_marker_evidence_scatter.svg"), width = 118, height = 96)
}

make_dataset_summary <- function(ds) {
  paths <- wgcna_downstream_paths("interpretable_summary", ds)

  # Correct upstream inputs:
  #   - group_effects/<dataset>/module_group_effects.csv
  #   - group_effects/<dataset>/supermodule_group_effects.csv
  #   - module_annotation/<dataset>/WGCNA_module_biological_annotation.csv
  #   - module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv
  #   - 04_wgcna_de_gsea_overlap/<dataset>/WGCNA_vs_DE_GSEA_overlap.csv
  module_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"))
  super_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"))
  module_annot <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"))
  super_annot <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"))
  overlap <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "04_wgcna_de_gsea_overlap", ds, "WGCNA_vs_DE_GSEA_overlap.csv"))
  super_comp <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_composition.csv"))
  module_super_map_file <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_to_supermodule_map_with_annotations.csv"))

  if (is.null(module_effects)) module_effects <- empty_group_effects(ds, "module", "missing module_group_effects.csv")
  if (is.null(super_effects)) super_effects <- empty_group_effects(ds, "supermodule", "missing supermodule_group_effects.csv")
  if (is.null(module_annot)) module_annot <- data.frame(dataset = character(), ModuleID = character(), microenvironment_class = character())
  if (is.null(super_annot)) super_annot <- data.frame(dataset = character(), SupermoduleID = character(), Supermodule_FinalLabel = character(), dominant_microenvironment_class = character())

  needed_super_cols <- c(
    "Supermodule_DisplayLabel",
    "Supermodule_LongLabel",
    "Macroprogram_Display",
    "dominant_microenvironment_class",
    "fraction_modules_microglia_supported",
    "fraction_modules_microglia_state_or_activation_supported",
    "fraction_modules_shared_microenvironment",
    "fraction_modules_neuropil_sensitive",
    "fraction_modules_other_cellular_or_vascular_sensitive",
    "fraction_modules_ambiguous",
    "marker_registry_version",
    "empirical_marker_set_version",
    "classification_rationale",
    "dominant_module_labels",
    "Supermodule_LabelSource",
    "Supermodule_LabelConfidence",
    "Supermodule_LabelRationale",
    "GO_label_confidence_class",
    "annotation_scope",
    "manual_label_status",
    "ManualReviewRequired",
    "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms",
    "raw_top_GO_label", "raw_module_label", "raw_hub_proteins",
    "raw_marker_or_signature_label",
    "cleaned_biological_label", "cleaned_biological_label_short",
    "cleaned_biological_label_source", "cleaned_biological_label_confidence",
    "cleaned_biological_label_rationale", "GO_label_relevance_flag",
    "GO_label_relevance_rationale",
    "microenvironment_caution_label", "microenvironment_caution_class",
    "microenvironment_caution_rationale",
    "Supermodule_ConservativeLabel",
    "Supermodule_CompositionLabel",
    "Supermodule_CompositionShortLabel",
    "Supermodule_CompositionDisplayLabel",
    "Supermodule_CompositionLabelSource",
    "Supermodule_CompositionConfidence",
    "Supermodule_CompositionRationale",
    "Supermodule_CleanPlotLabel",
    "DominantMemberTheme",
    "DominantMemberThemeFraction",
    "SecondMemberTheme",
    "SecondMemberThemeFraction",
    "MemberThemeCounts",
    "MemberThemeFractions",
    "n_distinct_member_themes",
    "is_multi_theme_supermodule",
    "themes_above_display_threshold",
    "themes_omitted_from_display_label",
    "TopMemberModuleLabels",
    "TopMemberGOTerms",
    "n_member_modules_with_informative_labels",
    "fraction_member_modules_with_informative_labels",
    "annotation_confidence", "annotation_basis", "annotation_downgrade_reason",
    "annotation_stable_across_thresholds", "unsafe_interpretation",
    "raw_annotation_label", "cleaned_annotation_label", "safe_display_label",
    "label_confidence", "label_basis", "label_downgrade_reason",
    "marker_fraction_primary", "marker_panels_supporting",
    "class_at_0.05", "class_at_0.10", "class_at_0.20"
  )
  for (nm in needed_super_cols) {
    if (!nm %in% names(super_annot)) super_annot[[nm]] <- rep(NA, nrow(super_annot))
  }

  super_annot <- super_annot |>
    dplyr::arrange(.data$SupermoduleID) |>
    dplyr::distinct(.data$dataset, .data$SupermoduleID, .keep_all = TRUE)

  module_join <- module_effects |>
    dplyr::left_join(module_annot, by = c("dataset" = "dataset", "module_id" = "ModuleID")) |>
    add_interpretation_sentences("module") |>
    ensure_columns(c(
      "module_label", "module_label.x", "module_label.y", "ModuleID", "ModuleColor", "module_id",
      "module_display_label", "ModuleLabel_Final", "endpoint_label", "Module_CleanPlotLabel",
      "module_biological_label", "module_biological_label_short", "module_label_display",
      "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms", "raw_top_GO_label",
      "raw_module_label", "raw_hub_proteins", "raw_marker_or_signature_label",
      "cleaned_biological_label", "cleaned_biological_label_short",
      "cleaned_biological_label_source", "cleaned_biological_label_confidence",
      "cleaned_biological_label_rationale", "GO_label_relevance_flag",
      "GO_label_relevance_rationale", "microenvironment_caution_label",
      "microenvironment_caution_class", "microenvironment_caution_rationale",
      "microenvironment_label", "microenvironment_class", "microenvironment_confidence",
      "supermodule_label", "supermodule_id", "supermodule_id_for_module", "module_supermodule_id",
      "p_value", "estimate", "FDR_global", "contrast"
    ))
  module_join$module_biological_label <- dplyr::coalesce(
    clean_label_value(module_join$module_biological_label),
    clean_label_value(module_join$cleaned_biological_label),
    clean_label_value(module_join$ModuleLabel_Final),
    clean_label_value(module_join$module_label.y),
    clean_label_value(module_join$module_label.x),
    clean_label_value(module_join$raw_top_GO_label),
    clean_label_value(module_join$endpoint_label)
  )
  module_join$module_biological_label_short <- dplyr::coalesce(
    clean_label_value(module_join$module_biological_label_short),
    clean_label_value(module_join$cleaned_biological_label_short),
    clean_label_value(module_join$module_biological_label)
  )
  module_join$module_label_display <- dplyr::coalesce(
    clean_label_value(module_join$module_label_display),
    clean_label_value(module_join$Module_CleanPlotLabel),
    paste(module_join$module_id, module_join$module_biological_label_short, sep = " | ")
  )
  module_join$module_label <- dplyr::coalesce(
    clean_label_value(module_join$module_label_display),
    clean_label_value(module_join$Module_CleanPlotLabel),
    clean_label_value(module_join$module_biological_label),
    clean_label_value(module_join$module_display_label),
    clean_label_value(module_join$module_id)
  )

  super_join <- super_effects |>
    dplyr::left_join(super_annot, by = c("dataset" = "dataset", "supermodule_id" = "SupermoduleID")) |>
    add_interpretation_sentences("supermodule") |>
    ensure_columns(c("Supermodule_DisplayLabel", "Supermodule_LongLabel", "Macroprogram_Display", "Supermodule_FinalLabel", "supermodule_label", "supermodule_id", "SupermoduleID", "Supermodule_CleanPlotLabel", "Supermodule_CompositionDisplayLabel", "Supermodule_CompositionLabel", "dominant_microenvironment_class", "p_value", "estimate", "FDR_global", "contrast", "spatial_unit", "effect_scope", "evidence_status", "direction", "MemberThemeCounts", "MemberThemeFractions", "n_distinct_member_themes", "is_multi_theme_supermodule", "themes_above_display_threshold", "themes_omitted_from_display_label"))
  super_join <- add_short_supermodule_labels(super_join)
  super_join$SupermoduleID <- dplyr::coalesce(
    clean_label_value(super_join$SupermoduleID),
    clean_label_value(super_join$supermodule_id),
    clean_label_value(col_or_na(super_join, "Supermodule_DataDrivenID")),
    clean_label_value(col_or_na(super_join, "Supermodule_DataDriven"))
  )
  super_join$supermodule_id <- dplyr::coalesce(clean_label_value(super_join$supermodule_id), clean_label_value(super_join$SupermoduleID))
  super_join$supermodule_label <- dplyr::coalesce(
    clean_label_value(super_join$supermodule_label),
    clean_label_value(super_join$Supermodule_CleanPlotLabel),
    clean_label_value(super_join$Supermodule_CompositionDisplayLabel),
    clean_label_value(super_join$Supermodule_CompositionLabel),
    clean_label_value(super_join$Supermodule_PlotLabel),
    clean_label_value(super_join$Supermodule_DisplayLabel)
  )
  super_join <- add_semantic_columns(super_join)
  super_join$supermodule_theme_label_qc_warning <- supermodule_theme_qc_warnings(super_join)

  module_super_map <- build_module_supermodule_map(
    module_effects = module_join,
    module_to_supermodule_map = module_super_map_file,
    supermodule_composition = super_comp,
    super_annot = super_annot
  )
  module_join <- attach_module_supermodules(module_join, module_super_map)
  module_join$ModulePlotLabel <- dplyr::coalesce(clean_label_value(module_join$module_label_display), clean_label_value(module_join$Module_CleanPlotLabel), make_module_plot_label(module_join))

  # Use the richer supermodule annotation/effect table as authoritative label source
  # for module-level plots. The module-to-supermodule map may only carry IDs or
  # weak fallback labels such as "Mixed / unresolved".
  supermodule_label_lookup <- super_join |>
    dplyr::mutate(
      Supermodule_PlotLabel_from_super = clean_label_value(.data$Supermodule_PlotLabel)
    ) |>
    dplyr::arrange(.data$supermodule_id) |>
    dplyr::distinct(.data$supermodule_id, .keep_all = TRUE) |>
    dplyr::select(
      supermodule_id_for_module = "supermodule_id",
      Supermodule_PlotLabel_from_super,
      Supermodule_DisplayLabel_from_super = "Supermodule_DisplayLabel",
      Macroprogram_Display_from_super = "Macroprogram_Display",
      Supermodule_FinalLabel_from_super = "Supermodule_FinalLabel",
      Supermodule_FullAnnotationLabel_from_super = "Supermodule_FullAnnotationLabel",
      Supermodule_DisplayShort_from_super = "Supermodule_DisplayShort",
      Supermodule_ConservativeLabel_from_super = "Supermodule_ConservativeLabel",
      Supermodule_CompositionLabel_from_super = "Supermodule_CompositionLabel",
      Supermodule_CompositionShortLabel_from_super = "Supermodule_CompositionShortLabel",
      Supermodule_CompositionDisplayLabel_from_super = "Supermodule_CompositionDisplayLabel",
      Supermodule_CompositionLabelSource_from_super = "Supermodule_CompositionLabelSource",
      Supermodule_CompositionConfidence_from_super = "Supermodule_CompositionConfidence",
      Supermodule_CleanPlotLabel_from_super = "Supermodule_CleanPlotLabel",
      MemberThemeCounts_from_super = "MemberThemeCounts",
      MemberThemeFractions_from_super = "MemberThemeFractions",
      n_distinct_member_themes_from_super = "n_distinct_member_themes",
      is_multi_theme_supermodule_from_super = "is_multi_theme_supermodule",
      themes_above_display_threshold_from_super = "themes_above_display_threshold",
      themes_omitted_from_display_label_from_super = "themes_omitted_from_display_label",
      supermodule_theme_label_qc_warning_from_super = "supermodule_theme_label_qc_warning",
      microenvironment_caution_label_from_super = "microenvironment_caution_label",
      microenvironment_caution_class_from_super = "microenvironment_caution_class",
      supermodule_label_from_super = "supermodule_label",
      SemanticProgramColor_from_super = "SemanticProgramColor",
      MacroprogramColorKey_from_super = "MacroprogramColorKey"
    ) |>
    dplyr::mutate(supermodule_id_for_module = clean_label_value(.data$supermodule_id_for_module))

  module_join <- module_join |>
    dplyr::mutate(supermodule_id_for_module = clean_label_value(.data$supermodule_id_for_module)) |>
    dplyr::left_join(supermodule_label_lookup, by = "supermodule_id_for_module")

  existing_supermodule_plot_label <- add_supermodule_id_prefix(
    module_join$supermodule_id_for_module,
    clean_supermodule_label_value(module_join$supermodule_label_for_module)
  )
  module_join$Supermodule_PlotLabel <- dplyr::coalesce(
    clean_supermodule_label_value(module_join$Supermodule_PlotLabel_from_super),
    existing_supermodule_plot_label,
    paste0(dplyr::coalesce(clean_label_value(module_join$supermodule_id_for_module), "SM??"), supermodule_label_sep(), "Mixed / unresolved")
  )

  module_join$Macroprogram_Display <- dplyr::coalesce(
    clean_label_value(col_or_na(module_join, "Macroprogram_Display")),
    clean_label_value(module_join$Macroprogram_Display_from_super)
  )

  module_join$MacroprogramColorKey <- dplyr::coalesce(
    clean_label_value(module_join$MacroprogramColorKey_from_super),
    clean_label_value(col_or_na(module_join, "MacroprogramColorKey"))
  )

  module_join$SemanticProgramColor <- dplyr::coalesce(
    clean_label_value(module_join$SemanticProgramColor_from_super),
    clean_label_value(col_or_na(module_join, "SemanticProgramColor"))
  )
  module_join$Supermodule_FullAnnotationLabel <- clean_label_value(module_join$Supermodule_FullAnnotationLabel_from_super)
  module_join$Supermodule_DisplayShort <- clean_label_value(module_join$Supermodule_DisplayShort_from_super)
  module_join$Supermodule_ConservativeLabel <- clean_label_value(module_join$Supermodule_ConservativeLabel_from_super)
  module_join$Supermodule_CompositionLabel <- clean_label_value(module_join$Supermodule_CompositionLabel_from_super)
  module_join$Supermodule_CompositionShortLabel <- clean_label_value(module_join$Supermodule_CompositionShortLabel_from_super)
  module_join$Supermodule_CompositionDisplayLabel <- clean_label_value(module_join$Supermodule_CompositionDisplayLabel_from_super)
  module_join$Supermodule_CompositionLabelSource <- clean_label_value(module_join$Supermodule_CompositionLabelSource_from_super)
  module_join$Supermodule_CompositionConfidence <- clean_label_value(module_join$Supermodule_CompositionConfidence_from_super)
  module_join$Supermodule_CleanPlotLabel <- clean_label_value(module_join$Supermodule_CleanPlotLabel_from_super)
  module_join$MemberThemeCounts <- clean_label_value(module_join$MemberThemeCounts_from_super)
  module_join$MemberThemeFractions <- clean_label_value(module_join$MemberThemeFractions_from_super)
  module_join$n_distinct_member_themes <- suppressWarnings(as.integer(module_join$n_distinct_member_themes_from_super))
  module_join$is_multi_theme_supermodule <- suppressWarnings(as.logical(module_join$is_multi_theme_supermodule_from_super))
  module_join$themes_above_display_threshold <- clean_label_value(module_join$themes_above_display_threshold_from_super)
  module_join$themes_omitted_from_display_label <- clean_label_value(module_join$themes_omitted_from_display_label_from_super)
  module_join$supermodule_theme_label_qc_warning <- clean_label_value(module_join$supermodule_theme_label_qc_warning_from_super)
  module_join$microenvironment_caution_label_for_supermodule <- clean_label_value(module_join$microenvironment_caution_label_from_super)
  module_join$microenvironment_caution_class_for_supermodule <- clean_label_value(module_join$microenvironment_caution_class_from_super)
  module_join$supermodule_id <- dplyr::coalesce(
    clean_label_value(module_join$supermodule_id),
    clean_label_value(module_join$supermodule_id_for_module),
    clean_label_value(module_join$module_supermodule_id),
    clean_label_value(col_or_na(module_join, "SupermoduleID"))
  )
  module_join$supermodule_label <- dplyr::coalesce(
    clean_label_value(module_join$supermodule_label),
    clean_label_value(module_join$Supermodule_CleanPlotLabel),
    clean_label_value(module_join$Supermodule_CompositionDisplayLabel),
    clean_label_value(module_join$Supermodule_CompositionLabel),
    clean_label_value(module_join$supermodule_label_for_module),
    clean_label_value(module_join$Supermodule_PlotLabel)
  )

  module_supermodule_label_qc <- module_join |>
    dplyr::transmute(
      dataset,
      module_id,
      module_key,
      supermodule_id_for_module,
      supermodule_label_for_module,
      Supermodule_PlotLabel,
      Supermodule_DisplayLabel = .data$Supermodule_DisplayLabel_from_super,
      Macroprogram_Display = .data$Macroprogram_Display_from_super,
      Supermodule_FinalLabel = .data$Supermodule_FinalLabel_from_super,
      Supermodule_FullAnnotationLabel = .data$Supermodule_FullAnnotationLabel_from_super,
      Supermodule_DisplayShort = .data$Supermodule_DisplayShort_from_super,
      Supermodule_ConservativeLabel = .data$Supermodule_ConservativeLabel_from_super,
      Supermodule_CompositionLabel = .data$Supermodule_CompositionLabel_from_super,
      Supermodule_CompositionShortLabel = .data$Supermodule_CompositionShortLabel_from_super,
      Supermodule_CompositionDisplayLabel = .data$Supermodule_CompositionDisplayLabel_from_super,
      Supermodule_CompositionLabelSource = .data$Supermodule_CompositionLabelSource_from_super,
      Supermodule_CompositionConfidence = .data$Supermodule_CompositionConfidence_from_super,
      Supermodule_CleanPlotLabel = .data$Supermodule_CleanPlotLabel_from_super,
      MemberThemeCounts = .data$MemberThemeCounts_from_super,
      MemberThemeFractions = .data$MemberThemeFractions_from_super,
      n_distinct_member_themes = .data$n_distinct_member_themes_from_super,
      is_multi_theme_supermodule = .data$is_multi_theme_supermodule_from_super,
      themes_above_display_threshold = .data$themes_above_display_threshold_from_super,
      themes_omitted_from_display_label = .data$themes_omitted_from_display_label_from_super,
      supermodule_theme_label_qc_warning = .data$supermodule_theme_label_qc_warning_from_super,
      microenvironment_caution_label = .data$microenvironment_caution_label_from_super,
      microenvironment_caution_class = .data$microenvironment_caution_class_from_super,
      supermodule_label = .data$supermodule_label_from_super,
      supermodule_map_source_for_module
    ) |>
    dplyr::distinct()

  mixed_plot_fraction <- mean(is_bad_supermodule_label(module_supermodule_label_qc$Supermodule_PlotLabel), na.rm = TRUE)
  if (is.finite(mixed_plot_fraction) && mixed_plot_fraction > 0.50) {
    warning(
      sprintf(
        "%.1f%% of module-level Supermodule_PlotLabel values are still Mixed / unresolved; upstream supermodule annotation may be weak.",
        100 * mixed_plot_fraction
      ),
      call. = FALSE
    )
  }

  module_join <- module_join |>
    dplyr::select(-dplyr::any_of(c(
      "Supermodule_PlotLabel_from_super",
      "Supermodule_DisplayLabel_from_super",
      "Macroprogram_Display_from_super",
      "Supermodule_FinalLabel_from_super",
      "Supermodule_FullAnnotationLabel_from_super",
      "Supermodule_DisplayShort_from_super",
      "Supermodule_ConservativeLabel_from_super",
      "Supermodule_CompositionLabel_from_super",
      "Supermodule_CompositionShortLabel_from_super",
      "Supermodule_CompositionDisplayLabel_from_super",
      "Supermodule_CompositionLabelSource_from_super",
      "Supermodule_CompositionConfidence_from_super",
      "Supermodule_CleanPlotLabel_from_super",
      "MemberThemeCounts_from_super",
      "MemberThemeFractions_from_super",
      "n_distinct_member_themes_from_super",
      "is_multi_theme_supermodule_from_super",
      "themes_above_display_threshold_from_super",
      "themes_omitted_from_display_label_from_super",
      "supermodule_theme_label_qc_warning_from_super",
      "microenvironment_caution_label_from_super",
      "microenvironment_caution_class_from_super",
      "supermodule_label_from_super",
      "SemanticProgramColor_from_super",
      "MacroprogramColorKey_from_super"
    )))

  module_join <- add_semantic_columns(module_join)

  label_candidates <- wgcna_make_label_candidates(module_join, super_join, dataset = ds) |>
    wgcna_score_label_candidates() |>
    wgcna_select_final_labels()
  final_label_lookup <- wgcna_build_final_label_lookup(label_candidates, module_join, super_join, dataset = ds)
  wgcna_validate_label_lookup(final_label_lookup)
  validate_table_schema(label_candidates, "wgcna_label_candidates", strict = TRUE)
  validate_table_schema(final_label_lookup, "wgcna_final_label_lookup", strict = TRUE)

  module_final_labels <- final_label_lookup |>
    dplyr::filter(.data$level == "module") |>
    dplyr::select(module_id = "entity_id", canonical_module_plot_label = "final_plot_label")
  super_final_labels <- final_label_lookup |>
    dplyr::filter(.data$level == "supermodule") |>
    dplyr::select(supermodule_id = "entity_id", canonical_supermodule_plot_label = "final_plot_label")

  module_join <- module_join |>
    dplyr::left_join(module_final_labels, by = "module_id") |>
    dplyr::left_join(super_final_labels, by = "supermodule_id") |>
    dplyr::mutate(
      ModulePlotLabel = dplyr::coalesce(.data$canonical_module_plot_label, .data$ModulePlotLabel),
      Supermodule_PlotLabel = dplyr::coalesce(.data$canonical_supermodule_plot_label, .data$Supermodule_PlotLabel)
    ) |>
    dplyr::select(-"canonical_module_plot_label", -"canonical_supermodule_plot_label")
  super_join <- super_join |>
    dplyr::left_join(super_final_labels, by = "supermodule_id") |>
    dplyr::mutate(
      Supermodule_PlotLabel = dplyr::coalesce(.data$canonical_supermodule_plot_label, .data$Supermodule_PlotLabel),
      supermodule_label = .data$Supermodule_PlotLabel
    ) |>
    dplyr::select(-"canonical_supermodule_plot_label")

  module_supermodule_join_qc <- module_join |>
    dplyr::distinct(.data$module_id, .data$module_key, .data$module_label, .data$module_biological_label, .data$module_biological_label_short, .data$module_label_display, .data$module_display_label, .data$Module_CleanPlotLabel, .data$microenvironment_label, .data$microenvironment_caution_label, .data$supermodule_id, .data$supermodule_id_for_module, .data$supermodule_label, .data$supermodule_label_for_module, .data$Supermodule_CleanPlotLabel, .data$supermodule_map_source_for_module)

  write_table_and_source(
    module_supermodule_join_qc,
    paths$tables,
    paths$source_data,
    "WGCNA_module_supermodule_join_qc.csv"
  )
  write_table_and_source(
    module_supermodule_label_qc,
    paths$tables,
    paths$source_data,
    "WGCNA_module_supermodule_label_qc.csv"
  )

  if (!is.null(overlap) && nrow(overlap) && "ModuleID" %in% names(overlap)) {
    overlap_summary <- overlap |>
      dplyr::group_by(.data$ModuleID) |>
      dplyr::summarise(
        best_wgcna_de_gsea_overlap = paste(utils::head(unique(.data$contrast), 5), collapse = ";"),
        .groups = "drop"
      )
    module_join <- module_join |> dplyr::left_join(overlap_summary, by = c("module_id" = "ModuleID"))
  }

  top_super <- super_join |>
    main_effect_rows() |>
    dplyr::filter(!is.na(.data$p_value)) |>
    add_plot_metrics() |>
    dplyr::group_by(.data$dataset, .data$supermodule_id, .data$contrast) |>
    dplyr::arrange(
      .data$evidence_rank,
      .data$FDR_global,
      .data$p_value,
      dplyr::desc(abs(.data$estimate)),
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      .data$evidence_rank,
      .data$FDR_global,
      .data$p_value,
      dplyr::desc(abs(.data$estimate))
    ) |>
    dplyr::select(-dplyr::any_of(c(
      "neg_log10_FDR",
      "neg_log10_P",
      "signed_FDR_score",
      "signed_P_score",
      "sig_label",
      "evidence_rank"
    ))) |>
    dplyr::slice_head(n = 50)

  supermodule_plot_label_qc <- super_join |>
    dplyr::distinct(dplyr::across(dplyr::any_of(c(
      "dataset", "supermodule_id", "SupermoduleID", "Supermodule_DisplayLabel",
      "Macroprogram_Display", "Supermodule_ShortLabel", "Supermodule_FinalLabel",
      "Supermodule_FullAnnotationLabel", "Supermodule_DisplayShort",
      "supermodule_label", "Supermodule_PlotLabel", "MacroprogramColorKey",
      "SemanticProgramColor", "Supermodule_LabelConfidence", "DataDrivenClusterSize",
      "Supermodule_CleanPlotLabel", "Supermodule_CompositionDisplayLabel",
      "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms",
      "raw_top_GO_label", "raw_module_label", "raw_hub_proteins",
      "raw_marker_or_signature_label", "cleaned_biological_label",
      "cleaned_biological_label_short", "cleaned_biological_label_source",
      "cleaned_biological_label_confidence", "GO_label_relevance_flag",
      "GO_label_relevance_rationale", "microenvironment_caution_label",
      "microenvironment_caution_class", "microenvironment_caution_rationale",
      "dominant_microenvironment_class", "dominant_module_labels", "Supermodule_LabelRationale",
      "Supermodule_ConservativeLabel", "Supermodule_CompositionLabel",
      "Supermodule_CompositionShortLabel", "Supermodule_CompositionLabelSource",
      "Supermodule_CompositionConfidence", "Supermodule_CompositionRationale",
      "DominantMemberTheme", "DominantMemberThemeFraction", "SecondMemberTheme",
      "SecondMemberThemeFraction", "TopMemberModuleLabels", "TopMemberGOTerms",
      "MemberThemeCounts", "MemberThemeFractions", "n_distinct_member_themes",
      "is_multi_theme_supermodule", "themes_above_display_threshold",
      "themes_omitted_from_display_label", "supermodule_theme_label_qc_warning",
      "n_member_modules_with_informative_labels", "fraction_member_modules_with_informative_labels",
      "annotation_confidence", "annotation_basis", "annotation_downgrade_reason",
      "annotation_stable_across_thresholds", "unsafe_interpretation",
      "raw_annotation_label", "cleaned_annotation_label", "safe_display_label",
      "label_confidence", "label_basis", "label_downgrade_reason",
      "marker_fraction_primary", "marker_panels_supporting",
      "ManualReviewRequired"
    ))), .keep_all = FALSE)
  supermodule_label_audit <- super_join |>
    dplyr::arrange(.data$supermodule_id) |>
    dplyr::distinct(.data$supermodule_id, .keep_all = TRUE) |>
    dplyr::transmute(
      SupermoduleID = .data$supermodule_id,
      Supermodule_CleanPlotLabel,
      Supermodule_PlotLabel,
      Supermodule_FullAnnotationLabel,
      Supermodule_DisplayShort,
      raw_GO_BP_terms,
      raw_GO_MF_terms,
      raw_GO_CC_terms,
      raw_top_GO_label,
      raw_module_label,
      raw_hub_proteins,
      raw_marker_or_signature_label,
      cleaned_biological_label,
      cleaned_biological_label_short,
      cleaned_biological_label_source,
      cleaned_biological_label_confidence,
      cleaned_biological_label_rationale,
      GO_label_relevance_flag,
      GO_label_relevance_rationale,
      microenvironment_caution_label,
      microenvironment_caution_class,
      microenvironment_caution_rationale,
      dominant_microenvironment_class,
      dominant_module_labels,
      Supermodule_ConservativeLabel,
      Supermodule_CompositionLabel,
      Supermodule_CompositionShortLabel,
      Supermodule_CompositionLabelSource,
      Supermodule_CompositionConfidence,
      Supermodule_CompositionRationale,
      DominantMemberTheme,
      DominantMemberThemeFraction,
      SecondMemberTheme,
      SecondMemberThemeFraction,
      MemberThemeCounts,
      MemberThemeFractions,
      n_distinct_member_themes,
      is_multi_theme_supermodule,
      themes_above_display_threshold,
      themes_omitted_from_display_label,
      supermodule_theme_label_qc_warning,
      TopMemberModuleLabels,
      TopMemberGOTerms,
      n_member_modules_with_informative_labels,
      fraction_member_modules_with_informative_labels,
      annotation_confidence,
      annotation_basis,
      annotation_downgrade_reason,
      annotation_stable_across_thresholds,
      unsafe_interpretation,
      raw_annotation_label,
      cleaned_annotation_label,
      safe_display_label,
      label_confidence,
      label_basis,
      label_downgrade_reason,
      marker_fraction_primary,
      marker_panels_supporting,
      ManualReviewRequired,
      Supermodule_LabelRationale
    )
  module_plot_label_qc <- module_join |>
    dplyr::distinct(dplyr::across(dplyr::any_of(c(
      "dataset", "module_id", "ModuleID", "ModuleColor", "module_eigengene",
      "ModuleDisplayID", "ModuleColorName", "ModulePlotLabel", "ModuleLabel_Final",
      "Module_CleanPlotLabel", "module_biological_label", "module_biological_label_short",
      "module_label_display", "module_label", "raw_GO_BP_terms", "raw_GO_MF_terms",
      "raw_GO_CC_terms", "raw_top_GO_label", "raw_module_label", "raw_hub_proteins",
      "raw_marker_or_signature_label", "cleaned_biological_label",
      "cleaned_biological_label_short", "cleaned_biological_label_source",
      "cleaned_biological_label_confidence", "GO_label_relevance_flag",
      "GO_label_relevance_rationale", "microenvironment_label",
      "microenvironment_caution_label", "microenvironment_caution_class",
      "annotation_confidence", "annotation_basis", "annotation_downgrade_reason",
      "annotation_stable_across_thresholds", "unsafe_interpretation",
      "raw_annotation_label", "cleaned_annotation_label", "safe_display_label",
      "label_confidence", "label_basis", "label_downgrade_reason",
      "marker_fraction_primary", "marker_panels_supporting",
      "supermodule_id", "supermodule_id_for_module", "Supermodule_CleanPlotLabel",
      "Supermodule_CompositionDisplayLabel", "Supermodule_PlotLabel",
      "MemberThemeCounts", "MemberThemeFractions", "n_distinct_member_themes",
      "is_multi_theme_supermodule", "themes_above_display_threshold",
      "themes_omitted_from_display_label", "supermodule_theme_label_qc_warning",
      "MacroprogramColorKey", "SemanticProgramColor"
    ))), .keep_all = FALSE)
  write_table_and_source(supermodule_plot_label_qc, paths$tables, paths$source_data, "WGCNA_supermodule_plot_label_qc.csv")
  write_table_and_source(supermodule_label_audit, paths$tables, paths$source_data, "WGCNA_supermodule_label_audit.csv")
  write_table_and_source(module_plot_label_qc, paths$tables, paths$source_data, "WGCNA_module_plot_label_qc.csv")
  write_table_and_source(label_candidates, paths$tables, paths$source_data, "WGCNA_label_candidates.csv")
  write_table_and_source(final_label_lookup, paths$tables, paths$source_data, "WGCNA_final_label_lookup.csv")

  write_table_and_source(super_join, paths$tables, paths$source_data, "WGCNA_supermodule_group_effects_interpretable.csv")
  write_table_and_source(module_join, paths$tables, paths$source_data, "WGCNA_module_group_effects_interpretable.csv")
  write_table_and_source(top_super, paths$tables, paths$source_data, "WGCNA_top_changed_supermodules.csv")

  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list(supermodules = super_join, modules = module_join, top_supermodules = top_super, label_candidates = label_candidates, final_label_lookup = final_label_lookup),
      file.path(paths$tables, "WGCNA_interpretable_summary.xlsx")
    )
  }

  plot_supermodule_main_heatmap(super_join, paths, ds)
  plot_supermodule_spatial_heatmap(super_join, paths, ds)
  plot_supermodule_membership_overview(module_join, super_join, paths, ds)
  plot_module_main_heatmap(module_join, paths, ds)
  plot_module_spatial_heatmap(module_join, paths, ds)
  plot_top_supermodules(top_super, paths, ds)
  plot_module_effects_by_supermodule(module_join, paths, ds)

  if (identical(ds, "microglia")) {
    plot_microglia_effect_heatmap(super_join, paths)
    plot_microglia_composition_heatmap(super_annot, paths)
    plot_microglia_marker_evidence(module_annot, paths)
  }

  write_run_manifest(
    file.path(paths$logs, "run_manifest.yml"),
    inputs = list(
      module_effects = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"),
      supermodule_effects = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"),
      module_annotation = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"),
      supermodule_annotation = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"),
      de_gsea_overlap = path_results("tables", "06_modules_WGCNA", "04_wgcna_de_gsea_overlap", ds, "WGCNA_vs_DE_GSEA_overlap.csv"),
      supermodule_composition = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_composition.csv"),
      module_to_supermodule_map = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_to_supermodule_map_with_annotations.csv")
    ),
    outputs = list(
      tables = paths$tables,
      source_data = paths$source_data,
      figures = paths$figures,
      label_candidates = file.path(paths$tables, "WGCNA_label_candidates.csv"),
      final_label_lookup = file.path(paths$tables, "WGCNA_final_label_lookup.csv")
    ),
    parameters = list(dataset = ds),
    notes = "Interpretable layer. Supermodules are treated as compressed overview units; module-level dotplots/heatmaps are exported as the biological-resolution view after explicitly joining the module-to-supermodule map; WGCNA_module_supermodule_join_qc.csv records whether labels came from 06 annotations and whether modules mapped to supermodules. Main heatmaps use global spatial-adjusted rows where available; spatial rows are plotted separately."
  )

  list(super = super_join, module = module_join, top = top_super)
}

make_cross_dataset_summary <- function(summaries) {
  paths_all <- wgcna_downstream_paths("interpretable_summary", "all")

  all_super <- dplyr::bind_rows(lapply(summaries, `[[`, "super"))
  if (!nrow(all_super)) {
    cross <- data.frame()
    write_table_and_source(cross, paths_all$tables, paths_all$source_data, "WGCNA_cross_dataset_supermodule_program_summary.csv")
    return(invisible(cross))
  }

  cross_base <- all_super |>
    orient_contrasts_for_plot() |>
    add_plot_metrics() |>
    dplyr::mutate(
      program_label = coalesce_chr(
        col_or_na(all_super, "Supermodule_PlotLabel"),
        col_or_na(all_super, "Supermodule_DisplayShort"),
        col_or_na(all_super, "Macroprogram_Display"),
        col_or_na(all_super, "Supermodule_DisplayLabel"),
        col_or_na(all_super, "Supermodule_FinalLabel"),
        col_or_na(all_super, "supermodule_label"),
        col_or_na(all_super, "supermodule_id")
      ),
      is_global = !is.na(.data$spatial_unit) & .data$spatial_unit == "global_spatial_adjusted",
      row_priority = dplyr::case_when(
        .data$is_global ~ 1L,
        TRUE ~ 2L
      )
    ) |>
    dplyr::filter(!is.na(.data$p_value))

  # One row per dataset x broad program x contrast.
  # Prefer global spatial-adjusted rows; if absent, use the strongest available spatial row.
  cross <- cross_base |>
    dplyr::group_by(.data$dataset, .data$program_label, .data$contrast) |>
    dplyr::arrange(.data$row_priority, .data$FDR_global, .data$p_value, dplyr::desc(abs(.data$estimate)), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      dataset,
      program_label,
      contrast,
      selected_spatial_unit = .data$spatial_unit,
      selected_effect_scope = .data$effect_scope,
      estimate,
      p_value,
      FDR_global,
      neg_log10_FDR,
      q_value,
      neg_log10_q,
      evidence_status,
      dominant_microenvironment_class,
      interpretation_sentence
    )

  write_table_and_source(cross, paths_all$tables, paths_all$source_data, "WGCNA_cross_dataset_supermodule_program_summary.csv")

  if (nrow(cross)) {
    plot_df <- cross |>
      dplyr::mutate(
        dataset = factor(.data$dataset, levels = c("microglia", "neuron_soma", "neuron_neuropil")),
        dataset_label = factor(dataset_label(as.character(.data$dataset)), levels = dataset_label(c("microglia", "neuron_soma", "neuron_neuropil"))),
        contrast = factor(.data$contrast, levels = contrast_plot_levels()),
        program_label_plot = label_wrap(.data$program_label, 34),
        sig_label = dplyr::case_when(
          !is.na(.data$q_value) & .data$q_value < 0.01 ~ "***",
          !is.na(.data$q_value) & .data$q_value < 0.05 ~ "**",
          !is.na(.data$q_value) & .data$q_value < 0.10 ~ "*",
          TRUE ~ ""
        )
      )

    # Limit crowded cross-dataset plot to strongest adjusted-evidence programs.
    keep_programs <- plot_df |>
      dplyr::group_by(.data$program_label_plot) |>
      dplyr::summarise(best_q = min(.data$q_value, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.data$best_q) |>
      dplyr::slice_head(n = 40) |>
      dplyr::pull(.data$program_label_plot)

    plot_df <- plot_df |> dplyr::filter(.data$program_label_plot %in% keep_programs)
    fill_limits <- effect_limits(plot_df$estimate)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$dataset_label, y = .data$program_label_plot, fill = .data$estimate)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.16) +
      ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), size = 1.9, color = "black") +
      ggplot2::facet_wrap(~ contrast, nrow = 1) +
      scale_effect_fill(fill_limits) +
      ggplot2::labs(
        title = "Cross-dataset supermodule effects",
        x = NULL, y = NULL
      ) +
      theme_clean(base_size = 6.8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

    save_svg(p, file.path(paths_all$figures, "cross_dataset_supermodule_program_effect_heatmap.svg"), width = 178, height = 118)
  }

  invisible(cross)
}

if (run$dry_run) {
  datasets <- if (DATASET_ARG == "all") valid_datasets() else DATASET_ARG
  for (ds in datasets) {
    paths <- wgcna_downstream_paths("interpretable_summary", ds)
    invisible(lapply(unlist(paths), dir_create))
    dry_run_line("Dataset", ds)
    dry_run_line(
      "Module effects",
      path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"),
      if (file.exists(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"))) "PASS" else "WARN"
    )
    dry_run_line(
      "Supermodule effects",
      path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"),
      if (file.exists(path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"))) "PASS" else "WARN"
    )
    dry_run_line(
      "Module annotation",
      path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"),
      if (file.exists(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"))) "PASS" else "WARN"
    )
    dry_run_line(
      "Supermodule annotation",
      path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"),
      if (file.exists(path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"))) "PASS" else "WARN"
    )
  }
  if (DATASET_ARG == "all") {
    dry_run_line("Cross-dataset output", path_results("tables", "06_modules_WGCNA", "interpretable_summary", "all"))
  }
  quit(status = 0, save = "no")
}

datasets <- if (DATASET_ARG == "all") valid_datasets() else DATASET_ARG
summaries <- lapply(datasets, make_dataset_summary)
names(summaries) <- datasets

if (DATASET_ARG == "all") {
  make_cross_dataset_summary(summaries)
}

message("WGCNA interpretable summary complete for dataset argument: ", DATASET_ARG)

