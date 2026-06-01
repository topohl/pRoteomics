#!/usr/bin/env Rscript

# Microglia-targeted signature enrichment for ranked proteomics contrasts.
#
# The microglia dataset represents microglia-enriched ROI / local
# microenvironment samples, not purified microglia. Neuropil is therefore used
# as a reference annotation layer only. This script never subtracts neuropil
# intensities or logFC values.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}
dataset_cli <- arg_value("--dataset", default = "")
if (nzchar(dataset_cli)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset_cli, source = "--dataset"))

MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- "microglia_targeted_signature_enrichment"
DATASET <- current_dataset()
REFERENCE_DATASET <- validate_dataset(Sys.getenv("PROTEOMICS_NEUROPIL_REFERENCE_DATASET", unset = "neuron_neuropil"), source = "PROTEOMICS_NEUROPIL_REFERENCE_DATASET")
DRY_RUN <- is_dry_run()
RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
SIGNATURE_METHOD_PRIORITY <- strsplit(Sys.getenv("PROTEOMICS_MICROGLIA_SIGNATURE_METHOD", unset = "limma_ranked_geneSetTest,fgsea,clusterProfiler_GSEA"), "[,;[:space:]]+")[[1]]
SIGNATURE_METHOD_PRIORITY <- SIGNATURE_METHOD_PRIORITY[nzchar(SIGNATURE_METHOD_PRIORITY)]

PATHS <- create_module_dirs(MODULE_ID, file.path(SUBSTEP_ID, DATASET))
invisible(lapply(PATHS, dir_create))

required_pkgs <- c("dplyr", "readr", "tidyr", "stringr", "purrr", "tibble", "ggplot2")
missing_required <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) && !isTRUE(DRY_RUN)) {
  stop("Missing required package(s): ", paste(missing_required, collapse = ", "), call. = FALSE)
}
if (!length(missing_required)) {
  suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

FIG_WIDTH_SINGLE <- 3.5
FIG_WIDTH_DOUBLE <- 7.2
FIG_HEIGHT_SHORT <- 2.6
FIG_HEIGHT_STANDARD <- 4.2
FIG_HEIGHT_TALL <- 5.2
FONT_SIZE_PANEL <- 7.0
FONT_SIZE_AXIS <- 7.0
FONT_SIZE_TICKS <- 6.5
FONT_SIZE_STRIP <- 7.0
FONT_SIZE_LEGEND <- 6.5

palette_diverging <- c(low = "#3C5488", mid = "#F7F7F7", high = "#B23A48")
palette_signature_class <- c(
  microglia_enriched_empirical = "#3C5488",
  microglia_enriched_reference_supported = "#4DBBD5",
  curated_microglia_program = "#91D1C2",
  mixed_microenvironment = "#7E6148",
  neuropil_shared = "#E64B35",
  ambiguous = "#8F8F8F"
)
palette_contrast_class <- c(
  within_region_condition = "#3C5488",
  cross_region_same_condition = "#00A087",
  cross_region_cross_condition = "#E64B35",
  same_region_same_condition = "#7E6148",
  unclassified = "#8F8F8F"
)
palette_direction <- c(up = "#B23A48", down = "#3C5488", neutral = "#8F8F8F")

theme_manuscript <- function(base_size = FONT_SIZE_PANEL, grid = c("none", "x", "y", "xy")) {
  grid <- match.arg(grid)
  theme <- ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = FONT_SIZE_PANEL, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = FONT_SIZE_LEGEND, color = "#4D4D4D"),
      axis.title = ggplot2::element_text(face = "plain", color = "black", size = FONT_SIZE_AXIS),
      axis.text = ggplot2::element_text(color = "black", size = FONT_SIZE_TICKS),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.4, "mm"),
      legend.title = ggplot2::element_text(face = "plain", size = FONT_SIZE_LEGEND),
      legend.text = ggplot2::element_text(color = "black", size = FONT_SIZE_LEGEND),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.box.spacing = grid::unit(0.8, "mm"),
      legend.spacing.x = grid::unit(0.8, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", color = "black", size = FONT_SIZE_STRIP),
      panel.border = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 4, 2, 3)
    )

  if (grid == "none") {
    theme <- theme + ggplot2::theme(panel.grid = ggplot2::element_blank())
  } else if (grid == "x") {
    theme <- theme + ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.25),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  } else if (grid == "y") {
    theme <- theme + ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.25),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  } else {
    theme <- theme + ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank()
    )
  }

  theme
}

scale_fill_diverging <- function(name = "NES") {
  ggplot2::scale_fill_gradient2(
    low = palette_diverging[["low"]],
    mid = palette_diverging[["mid"]],
    high = palette_diverging[["high"]],
    midpoint = 0,
    name = name
  )
}

scale_color_diverging <- function(name = "NES") {
  ggplot2::scale_color_gradient2(
    low = palette_diverging[["low"]],
    mid = palette_diverging[["mid"]],
    high = palette_diverging[["high"]],
    midpoint = 0,
    name = name
  )
}

pick_col <- function(df, candidates) {
  detect_column(df, candidates, required = FALSE)
}

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^ +| +$", "", x)
  toupper(x[nzchar(x) & !is.na(x)])
}

normalize_id_all <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^ +| +$", "", x)
  x <- toupper(x)
  x[is.na(x)] <- ""
  x
}

signature_sets <- list(
  homeostatic_microglia = c("P2RY12", "TMEM119", "CX3CR1", "SALL1", "GPR34", "HEXB", "CSF1R", "OLFML3", "SELPLG", "SIGLECH", "FCRLS", "MAF"),
  disease_associated_microglia_DAM = c("TREM2", "APOE", "TYROBP", "AXL", "LPL", "CST7", "CTSD", "CTSB", "LGALS3", "ITGAX", "CLEC7A", "SPP1"),
  interferon_response_microglia = c("STAT1", "STAT2", "IRF7", "IRF9", "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1A", "OASL2"),
  phagolysosomal_microglia = c("LAMP1", "LAMP2", "CTSB", "CTSD", "CTSS", "HEXB", "LAPTM5", "ATP6V0D1", "ATP6V1B2", "LYZ2", "CD68", "RAB7A"),
  complement_phagocytosis = c("C1QA", "C1QB", "C1QC", "C3", "C4B", "CFH", "CR1L", "ITGAM", "ITGB2", "TYROBP", "FCGR3", "VSIG4"),
  antigen_presentation_MHC = c("H2-AA", "H2-AB1", "H2-EB1", "H2-D1", "H2-K1", "B2M", "CD74", "CTSS", "TAP1", "TAP2", "PSMB8", "PSMB9"),
  chemokine_cytokine_microglia = c("CCL2", "CCL3", "CCL4", "CCL5", "CXCL10", "CXCL16", "IL1B", "TNF", "TGFBR1", "CSF1R", "CX3CR1", "NFKB1"),
  lipid_metabolism_microglia = c("APOE", "LPL", "ABCA1", "ABCG1", "SOAT1", "LIPA", "PLIN2", "CST7", "TREM2", "NPC2", "LITAF", "FABP5"),
  synapse_pruning_microglia = c("C1QA", "C1QB", "C1QC", "C3", "ITGAM", "ITGB2", "TREM2", "TYROBP", "CX3CR1", "P2RY12", "GRN", "MERTK"),
  oxidative_stress_mitochondrial_microglia = c("SOD1", "SOD2", "PRDX1", "PRDX3", "GPX1", "TXN", "TXNRD1", "NFE2L2", "PARK7", "VDAC1", "NDUFA9", "UQCRC2")
)

marker_sets <- list(
  microglia = unique(unlist(signature_sets[c(
    "homeostatic_microglia", "disease_associated_microglia_DAM",
    "phagolysosomal_microglia", "complement_phagocytosis"
  )], use.names = FALSE)),
  synaptic_neuronal = c("SYN1", "SYP", "SNAP25", "STX1A", "STXBP1", "DLG4", "CAMK2A", "CAMK2B", "MAP2", "RBFOX3", "TUBB3", "GRIN1", "GRIA1", "VAMP2")
)

make_term2gene <- function(signatures) {
  tibble::tibble(
    term = rep(names(signatures), lengths(signatures)),
    gene = normalize_id(unlist(signatures, use.names = FALSE))
  ) %>%
    dplyr::mutate(signature_source = "curated") %>%
    dplyr::distinct()
}

read_mouse_id_map <- function() {
  candidates <- c(
    path_external("MOUSE_10090_idmapping.dat"),
    path_external("MOUSE_10090_idmapping.dat.gz")
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) return(tibble::tibble(UNIPROT = character(), SYMBOL = character()))
  df <- tryCatch(readr::read_tsv(path, col_names = c("UNIPROT", "type", "value"), show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (is.null(df)) return(tibble::tibble(UNIPROT = character(), SYMBOL = character()))
  df %>%
    dplyr::filter(.data$type %in% c("Gene_Name", "Gene_ORFName", "Gene_Synonym")) %>%
    dplyr::transmute(UNIPROT = normalize_id_all(.data$UNIPROT), SYMBOL = normalize_id_all(.data$value)) %>%
    dplyr::filter(nzchar(.data$UNIPROT), nzchar(.data$SYMBOL)) %>%
    dplyr::distinct()
}

expand_term2gene_for_uniprot <- function(term2gene, id_map) {
  if (!nrow(id_map)) return(term2gene)
  mapped <- suppressWarnings(
    term2gene %>%
      dplyr::inner_join(id_map, by = c("gene" = "SYMBOL")) %>%
      dplyr::transmute(term, gene = UNIPROT, signature_source)
  )
  dplyr::bind_rows(term2gene, mapped) %>% dplyr::distinct()
}

add_symbol_aliases <- function(tbl, id_map) {
  if (!nrow(tbl) || !"gene" %in% names(tbl) || !nrow(id_map)) return(tbl)
  symbols <- tbl %>%
    dplyr::inner_join(id_map, by = c("gene" = "UNIPROT")) %>%
    dplyr::transmute(gene, gene_symbol = SYMBOL)
  tbl %>%
    dplyr::left_join(symbols %>% dplyr::group_by(.data$gene) %>% dplyr::summarise(gene_symbol = dplyr::first(.data$gene_symbol), .groups = "drop"), by = "gene")
}

contrast_dir <- function(dataset) {
  path_processed("02_id_mapping", "mapped", dataset, "forward", "per_file")
}

safe_read_contrast <- function(path) {
  tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
}

extract_unit_condition <- function(unit) {
  x <- tolower(as.character(unit))
  condition <- stringr::str_extract(x, "(con|res|sus)$")
  base <- ifelse(is.na(condition), x, sub("(con|res|sus)$", "", x))
  list(unit = base, condition = condition %||% NA_character_)
}

parse_region_from_unit <- function(unit) {
  x <- tolower(as.character(unit))
  x <- sub("\\.csv$", "", basename(x))
  x <- sub("(con|res|sus)$", "", x)
  x <- sub("microglia$", "", x)
  dplyr::case_when(
    grepl("^ca1", x) ~ "CA1",
    grepl("^ca2", x) ~ "CA2",
    grepl("^ca3", x) ~ "CA3",
    grepl("^dg|^hilus", x) ~ "DG",
    TRUE ~ toupper(x)
  )
}

collapse_neuropil_layer_to_region <- function(unit) {
  parse_region_from_unit(unit)
}

classify_contrast_class <- function(left_region, right_region, left_condition, right_condition) {
  if (is.na(left_region) || is.na(right_region) || is.na(left_condition) || is.na(right_condition)) {
    return("unclassified")
  }
  if (left_region == right_region && left_condition != right_condition) return("within_region_condition")
  if (left_region != right_region && left_condition == right_condition) return("cross_region_same_condition")
  if (left_region != right_region && left_condition != right_condition) return("cross_region_cross_condition")
  if (left_region == right_region && left_condition == right_condition) return("same_region_same_condition")
  "unclassified"
}

parse_comparison <- function(comparison) {
  parts <- strsplit(sub("\\.csv$", "", basename(comparison)), "_", fixed = TRUE)[[1]]
  left <- parts[[1]]
  right <- if (length(parts) >= 2) parts[[2]] else NA_character_
  left_parsed <- extract_unit_condition(left)
  right_parsed <- extract_unit_condition(right)
  left_region <- parse_region_from_unit(left_parsed$unit)
  right_region <- parse_region_from_unit(right_parsed$unit)
  contrast_class <- classify_contrast_class(left_region, right_region, left_parsed$condition, right_parsed$condition)
  cond_pair <- paste(na.omit(c(left_parsed$condition, right_parsed$condition)), collapse = "_vs_")
  if (!nzchar(cond_pair)) cond_pair <- NA_character_
  tibble::tibble(
    comparison = sub("\\.csv$", "", basename(comparison)),
    left_unit = left_parsed$unit,
    right_unit = right_parsed$unit,
    left_region = left_region,
    right_region = right_region,
    left_condition = left_parsed$condition,
    right_condition = right_parsed$condition,
    contrast_class = contrast_class,
    condition_pair = cond_pair,
    region = left_region,
    layer = ifelse(grepl("microglia", tolower(left)), NA_character_, sub("^[Cc][Aa][123]|^[Dd][Gg]", "", left_parsed$unit)),
    region_level_unit = collapse_neuropil_layer_to_region(left_parsed$unit)
  )
}

map_comparison_to_region_level <- function(comparison) {
  parse_comparison(comparison) %>%
    dplyr::select(
      .data$comparison, .data$region, .data$region_level_unit, .data$condition_pair,
      .data$left_unit, .data$right_unit, .data$left_region, .data$right_region,
      .data$left_condition, .data$right_condition, .data$contrast_class
    )
}

prepare_ranked_contrast <- function(path, dataset) {
  df <- safe_read_contrast(path)
  if (is.null(df) || !nrow(df)) return(NULL)
  id_col <- pick_col(df, c("gene_symbol", "Gene", "SYMBOL", "UNIPROT", "protein"))
  stat_col <- pick_col(df, c("log2fc", "logFC"))
  p_col <- pick_col(df, c("padj", "p.adjust", "FDR", "pval", "pvalue", "P.Value"))
  if (is.na(id_col) || is.na(stat_col)) return(NULL)
  ranks <- suppressWarnings(as.numeric(df[[stat_col]]))
  ids <- normalize_id_all(df[[id_col]])
  keep <- !is.na(ranks) & nzchar(ids) & !duplicated(ids)
  if (!any(keep)) return(NULL)
  meta <- parse_comparison(basename(path))
  tibble::tibble(
    dataset = dataset,
    comparison = meta$comparison[[1]],
    gene = ids[keep],
    rank_stat = ranks[keep],
    p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]][keep])) else NA_real_,
    input_file = normalizePath(path, winslash = "/", mustWork = FALSE),
    region = meta$region[[1]],
    layer = meta$layer[[1]],
    region_level_unit = meta$region_level_unit[[1]],
    condition_pair = meta$condition_pair[[1]],
    left_unit = meta$left_unit[[1]],
    right_unit = meta$right_unit[[1]],
    left_region = meta$left_region[[1]],
    right_region = meta$right_region[[1]],
    left_condition = meta$left_condition[[1]],
    right_condition = meta$right_condition[[1]],
    contrast_class = meta$contrast_class[[1]]
  )
}

load_ranked_contrasts <- function(dataset) {
  dir <- contrast_dir(dataset)
  files <- if (dir.exists(dir)) list.files(dir, pattern = "\\.csv$", full.names = TRUE) else character()
  rows <- purrr::map(files, prepare_ranked_contrast, dataset = dataset)
  dplyr::bind_rows(rows)
}

rank_support_summary <- function(ranked) {
  if (!nrow(ranked)) {
    return(tibble::tibble(
      gene = character(), n_contrasts = integer(), detection_fraction = numeric(),
      mean_abs_rank_stat = numeric(), mean_signed_rank_stat = numeric(),
      top10_fraction = numeric(), top20_fraction = numeric()
    ))
  }
  ranked %>%
    dplyr::group_by(.data$comparison) %>%
    dplyr::mutate(
      abs_rank = abs(.data$rank_stat),
      abs_percentile = dplyr::percent_rank(.data$abs_rank),
      signed_percentile = dplyr::percent_rank(.data$rank_stat)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::summarise(
      n_contrasts = dplyr::n_distinct(.data$comparison),
      detection_fraction = dplyr::n_distinct(.data$comparison) / dplyr::n_distinct(ranked$comparison),
      mean_abs_rank_stat = mean(abs(.data$rank_stat), na.rm = TRUE),
      mean_signed_rank_stat = mean(.data$rank_stat, na.rm = TRUE),
      mean_abs_percentile = mean(.data$abs_percentile, na.rm = TRUE),
      mean_signed_percentile = mean(.data$signed_percentile, na.rm = TRUE),
      top10_fraction = mean(.data$abs_percentile >= 0.90, na.rm = TRUE),
      top20_fraction = mean(.data$abs_percentile >= 0.80, na.rm = TRUE),
      .groups = "drop"
    )
}

derive_empirical_signatures <- function(micro_ranked, neuropil_ranked, id_map) {
  micro <- rank_support_summary(micro_ranked) %>%
    dplyr::rename_with(~paste0("microglia_", .x), -gene)
  neuropil <- rank_support_summary(neuropil_ranked) %>%
    dplyr::rename_with(~paste0("neuropil_", .x), -gene)
  empirical <- dplyr::full_join(micro, neuropil, by = "gene") %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~dplyr::coalesce(.x, 0)),
      rank_delta = .data$microglia_mean_abs_percentile - .data$neuropil_mean_abs_percentile,
      top10_delta = .data$microglia_top10_fraction - .data$neuropil_top10_fraction,
      empirical_class = dplyr::case_when(
        .data$microglia_detection_fraction >= 0.50 &
          .data$microglia_top20_fraction >= 0.20 &
          .data$rank_delta >= 0.08 &
          .data$top10_delta >= 0.03 ~ "empirical_microglia_enriched",
        .data$microglia_detection_fraction >= 0.50 &
          .data$neuropil_detection_fraction >= 0.50 &
          .data$microglia_top20_fraction >= 0.20 &
          .data$neuropil_top20_fraction >= 0.20 &
          abs(.data$rank_delta) < 0.08 ~ "empirical_neuropil_shared",
        TRUE ~ "not_empirical_signature"
      )
    ) %>%
    add_symbol_aliases(id_map) %>%
    dplyr::arrange(.data$empirical_class, dplyr::desc(.data$rank_delta), dplyr::desc(.data$microglia_mean_abs_percentile))

  micro_sig <- empirical %>%
    dplyr::filter(.data$empirical_class == "empirical_microglia_enriched") %>%
    dplyr::slice_head(n = 250)
  shared_sig <- empirical %>%
    dplyr::filter(.data$empirical_class == "empirical_neuropil_shared") %>%
    dplyr::arrange(dplyr::desc(pmin(.data$microglia_mean_abs_percentile, .data$neuropil_mean_abs_percentile))) %>%
    dplyr::slice_head(n = 250)

  term2gene_empirical <- dplyr::bind_rows(
    tibble::tibble(term = "empirical_microglia_enriched", gene = micro_sig$gene, signature_source = "empirical_microglia_vs_neuropil"),
    tibble::tibble(term = "empirical_neuropil_shared", gene = shared_sig$gene, signature_source = "empirical_microglia_vs_neuropil")
  ) %>%
    dplyr::filter(nzchar(.data$gene)) %>%
    dplyr::distinct()

  diagnostics <- tibble::tibble(
    check = c("empirical_gene_rows", "empirical_microglia_enriched_genes", "empirical_neuropil_shared_genes"),
    status = c(ifelse(nrow(empirical) > 0, "PASS", "WARN"), ifelse(nrow(micro_sig) >= 10, "PASS", "WARN"), ifelse(nrow(shared_sig) >= 10, "PASS", "WARN")),
    detail = as.character(c(nrow(empirical), nrow(micro_sig), nrow(shared_sig)))
  )

  list(all = empirical, microglia = micro_sig, shared = shared_sig, term2gene = term2gene_empirical, diagnostics = diagnostics)
}

ewce_results_path <- function(dataset = DATASET) {
  candidates <- c(
    path_processed("05_celltype_enrichment_EWCE", "EWCE_E9", dataset, "EWCE_results_full.rds"),
    path_processed("05_celltype_enrichment_EWCE", "EWCE_E9", "EWCE_results_full.rds"),
    path_processed("05_celltype_enrichment_EWCE", "EWCE_E9", "neuron_neuropil", "EWCE_results_full.rds")
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) NA_character_ else hit
}

reference_support_from_ewce_results <- function(dataset = DATASET) {
  path <- ewce_results_path(dataset)
  if (is.na(path)) {
    return(list(
      gene_support = tibble::tibble(),
      diagnostics = tibble::tibble(check = "ewce_results_full_rds", status = "WARN", detail = "EWCE_results_full.rds not found")
    ))
  }
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || is.null(obj$target_gene_tbl)) {
    return(list(
      gene_support = tibble::tibble(),
      diagnostics = tibble::tibble(check = "ewce_results_full_rds", status = "WARN", detail = paste("Unreadable or missing target_gene_tbl:", path))
    ))
  }
  # This is not cell-type specificity, but records whether EWCE target lists already
  # contain each protein/gene in this project-specific reference workflow.
  support <- obj$target_gene_tbl %>%
    dplyr::mutate(gene_symbol = normalize_id_all(.data$Gene)) %>%
    dplyr::group_by(.data$gene_symbol) %>%
    dplyr::summarise(
      ewce_target_list_fraction = dplyr::n_distinct(.data$Target) / max(1, dplyr::n_distinct(obj$target_gene_tbl$Target)),
      ewce_target_top100_fraction = mean(.data$TopN <= 100, na.rm = TRUE),
      .groups = "drop"
    )
  list(
    gene_support = support,
    diagnostics = tibble::tibble(check = "ewce_results_full_rds", status = "PASS", detail = path)
  )
}

read_reference_specificity <- function(id_map) {
  if (!requireNamespace("ewceData", quietly = TRUE)) {
    return(list(
      specificity = tibble::tibble(),
      diagnostics = tibble::tibble(check = "ewceData_ctd", status = "WARN", detail = "ewceData package not available")
    ))
  }
  ctd <- tryCatch(ewceData::ctd(), error = function(e) NULL)
  if (is.null(ctd) || !length(ctd)) {
    return(list(
      specificity = tibble::tibble(),
      diagnostics = tibble::tibble(check = "ewceData_ctd", status = "WARN", detail = "Could not load ewceData::ctd()")
    ))
  }

  specificity_rows <- lapply(seq_along(ctd), function(level) {
    mat <- ctd[[level]]$specificity
    if (is.null(mat) || !nrow(mat)) return(NULL)
    tibble::as_tibble(mat, rownames = "gene_symbol_raw") %>%
      tidyr::pivot_longer(-gene_symbol_raw, names_to = "celltype", values_to = "specificity") %>%
      dplyr::mutate(
        gene_symbol = normalize_id_all(.data$gene_symbol_raw),
        celltype_group = dplyr::case_when(
          grepl("micro|mgl|pvm", .data$celltype, ignore.case = TRUE) ~ "microglia",
          grepl("pyr|int|neuron|glut|gaba|sst|vip|pvalb", .data$celltype, ignore.case = TRUE) ~ "neuron",
          grepl("astro", .data$celltype, ignore.case = TRUE) ~ "astrocyte",
          grepl("oligo|opc", .data$celltype, ignore.case = TRUE) ~ "oligodendrocyte",
          TRUE ~ "other"
        )
      ) %>%
      dplyr::group_by(.data$gene_symbol, .data$celltype_group) %>%
      dplyr::summarise(specificity = max(.data$specificity, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = "celltype_group", values_from = "specificity", values_fill = 0) %>%
      dplyr::mutate(reference_level = level)
  })

  specificity <- dplyr::bind_rows(specificity_rows)
  if (!nrow(specificity)) {
    return(list(
      specificity = tibble::tibble(),
      diagnostics = tibble::tibble(check = "ewceData_specificity_rows", status = "WARN", detail = "No specificity rows parsed")
    ))
  }
  for (col in c("microglia", "neuron", "astrocyte", "oligodendrocyte")) {
    if (!col %in% names(specificity)) specificity[[col]] <- 0
  }
  specificity <- specificity %>%
    dplyr::group_by(.data$gene_symbol) %>%
    dplyr::summarise(
      reference_microglia_specificity = max(.data$microglia, na.rm = TRUE),
      reference_neuron_specificity = max(.data$neuron, na.rm = TRUE),
      reference_astrocyte_specificity = max(.data$astrocyte, na.rm = TRUE),
      reference_oligodendrocyte_specificity = max(.data$oligodendrocyte, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      reference_celltype_support = dplyr::case_when(
        .data$reference_microglia_specificity >= pmax(.data$reference_neuron_specificity, .data$reference_astrocyte_specificity, .data$reference_oligodendrocyte_specificity) &
          .data$reference_microglia_specificity >= 0.20 ~ "microglia_supported",
        .data$reference_neuron_specificity >= 0.20 ~ "neuron_supported",
        .data$reference_astrocyte_specificity >= 0.20 ~ "astrocyte_supported",
        .data$reference_oligodendrocyte_specificity >= 0.20 ~ "oligodendrocyte_supported",
        TRUE ~ "low_specificity_or_unclassified"
      )
    )

  if (nrow(id_map)) {
    id_map_symbol_unique <- id_map %>%
      dplyr::group_by(.data$SYMBOL) %>%
      dplyr::summarise(UNIPROT = dplyr::first(.data$UNIPROT), .groups = "drop")
    specificity_uniprot <- specificity %>%
      dplyr::inner_join(id_map_symbol_unique, by = c("gene_symbol" = "SYMBOL")) %>%
      dplyr::transmute(
        gene = UNIPROT,
        gene_symbol,
        reference_microglia_specificity,
        reference_neuron_specificity,
        reference_astrocyte_specificity,
        reference_oligodendrocyte_specificity,
        reference_celltype_support
      )
    specificity_symbol <- specificity %>%
      dplyr::transmute(
        gene = gene_symbol,
        gene_symbol,
        reference_microglia_specificity,
        reference_neuron_specificity,
        reference_astrocyte_specificity,
        reference_oligodendrocyte_specificity,
        reference_celltype_support
      )
    specificity <- dplyr::bind_rows(specificity_symbol, specificity_uniprot) %>% dplyr::distinct()
  } else {
    specificity <- specificity %>% dplyr::mutate(gene = .data$gene_symbol)
  }

  list(
    specificity = specificity,
    diagnostics = tibble::tibble(check = "ewceData_specificity_rows", status = "PASS", detail = as.character(nrow(specificity)))
  )
}

make_reference_term2gene <- function(reference_specificity) {
  if (!nrow(reference_specificity)) {
    return(tibble::tibble(term = character(), gene = character(), signature_source = character()))
  }
  reference_specificity %>%
    dplyr::filter(
      .data$reference_celltype_support == "microglia_supported",
      .data$reference_microglia_specificity >= 0.20,
      .data$reference_microglia_specificity >= .data$reference_neuron_specificity,
      .data$reference_microglia_specificity >= .data$reference_astrocyte_specificity,
      .data$reference_microglia_specificity >= .data$reference_oligodendrocyte_specificity
    ) %>%
    dplyr::arrange(dplyr::desc(.data$reference_microglia_specificity)) %>%
    dplyr::slice_head(n = 300) %>%
    dplyr::transmute(term = "reference_atlas_EWCE_microglia_specific", gene = .data$gene, signature_source = "reference_atlas_EWCE") %>%
    dplyr::distinct()
}

method_status <- function() {
  tibble::tibble(
    method = c("fgsea", "limma_ranked_geneSetTest", "clusterProfiler_GSEA"),
    available = c(
      requireNamespace("fgsea", quietly = TRUE),
      requireNamespace("limma", quietly = TRUE),
      requireNamespace("clusterProfiler", quietly = TRUE)
    ),
    note = c(
      "fgsea::fgsea with custom TERM2GENE pathways",
      "limma::geneSetTest on ranked contrast statistics; camera/roast require expression matrix/design and are not used from mapped tables",
      "clusterProfiler::GSEA with custom TERM2GENE fallback"
    )
  ) %>%
    dplyr::mutate(priority = match(.data$method, SIGNATURE_METHOD_PRIORITY)) %>%
    dplyr::arrange(is.na(.data$priority), .data$priority)
}

run_one_gsea <- function(stats, term2gene, methods) {
  stats <- stats[!is.na(stats)]
  stats <- sort(stats, decreasing = TRUE)
  stats <- stats[!duplicated(names(stats))]
  pathways <- split(term2gene$gene, term2gene$term)
  pathways <- lapply(pathways, unique)
  pathways <- pathways[vapply(pathways, function(g) sum(g %in% names(stats)) >= 3, logical(1))]
  if (!length(pathways) || length(stats) < 10) return(tibble())

  for (method_i in methods$method[methods$available]) {
    if (identical(method_i, "fgsea")) {
    res <- tryCatch({
      if ("fgseaSimple" %in% getNamespaceExports("fgsea")) {
        fgsea::fgseaSimple(pathways = pathways, stats = stats, minSize = 3, maxSize = 500, nperm = 1000, nproc = 1)
      } else if (requireNamespace("BiocParallel", quietly = TRUE)) {
        fgsea::fgsea(pathways = pathways, stats = stats, minSize = 3, maxSize = 500, BPPARAM = BiocParallel::SerialParam())
      } else {
        fgsea::fgsea(pathways = pathways, stats = stats, minSize = 3, maxSize = 500, nproc = 1)
      }
    }, error = function(e) NULL)
    if (!is.null(res) && nrow(res)) {
      return(as.data.frame(res) %>%
        tibble::as_tibble() %>%
        dplyr::transmute(signature = .data$pathway, method = "fgsea", NES = .data$NES, pvalue = .data$pval, padj = .data$padj, set_size = .data$size, leading_edge = vapply(.data$leadingEdge, paste, character(1), collapse = "/")))
    }
    }

    if (identical(method_i, "limma_ranked_geneSetTest")) {
    out <- lapply(names(pathways), function(term) {
      idx <- which(names(stats) %in% pathways[[term]])
      if (length(idx) < 3) return(NULL)
      p <- tryCatch(limma::geneSetTest(idx, stats, alternative = "either"), error = function(e) NA_real_)
      leading <- idx[order(abs(stats[idx]), decreasing = TRUE)]
      leading <- leading[seq_len(min(25, length(leading)))]
      tibble::tibble(
        signature = term,
        method = "limma_ranked_geneSetTest",
        NES = mean(stats[idx], na.rm = TRUE) / stats::sd(stats, na.rm = TRUE),
        pvalue = p,
        set_size = length(idx),
        leading_edge = paste(names(stats)[leading], collapse = "/")
      )
    })
    res <- dplyr::bind_rows(out)
    if (nrow(res)) return(res %>% dplyr::mutate(padj = p.adjust(.data$pvalue, method = "BH")))
    }

    if (identical(method_i, "clusterProfiler_GSEA")) {
    res <- tryCatch(clusterProfiler::GSEA(stats, TERM2GENE = term2gene, minGSSize = 3, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE), error = function(e) NULL)
    if (!is.null(res) && nrow(as.data.frame(res))) {
      return(as.data.frame(res) %>%
        tibble::as_tibble() %>%
        dplyr::transmute(signature = .data$ID, method = "clusterProfiler_GSEA", NES = .data$NES, pvalue = .data$pvalue, padj = .data$p.adjust, set_size = .data$setSize, leading_edge = .data$core_enrichment))
    }
    }
  }

  tibble()
}

summarise_reference_for_genes <- function(gene_string, reference_specificity) {
  empty <- tibble::tibble(
    reference_microglia_specificity = NA_real_,
    reference_neuron_specificity = NA_real_,
    reference_astrocyte_specificity = NA_real_,
    reference_oligodendrocyte_specificity = NA_real_,
    reference_celltype_support = NA_character_
  )
  if (!nrow(reference_specificity) || is.na(gene_string) || !nzchar(gene_string)) return(empty)
  genes <- normalize_id(unlist(strsplit(gene_string, "[/;,|[:space:]]+")))
  ref <- reference_specificity %>% dplyr::filter(.data$gene %in% genes)
  if (!nrow(ref)) return(empty)
  support <- ref$reference_celltype_support
  support <- support[!is.na(support) & nzchar(support)]
  tibble::tibble(
    reference_microglia_specificity = mean(ref$reference_microglia_specificity, na.rm = TRUE),
    reference_neuron_specificity = mean(ref$reference_neuron_specificity, na.rm = TRUE),
    reference_astrocyte_specificity = mean(ref$reference_astrocyte_specificity, na.rm = TRUE),
    reference_oligodendrocyte_specificity = mean(ref$reference_oligodendrocyte_specificity, na.rm = TRUE),
    reference_celltype_support = if (length(support)) names(sort(table(support), decreasing = TRUE))[[1]] else NA_character_
  )
}

run_signature_enrichment <- function(ranked, term2gene, reference_specificity = tibble::tibble()) {
  methods <- method_status()
  signature_meta <- term2gene %>%
    dplyr::transmute(signature = .data$term, signature_source = .data$signature_source) %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$signature) %>%
    dplyr::summarise(signature_source = dplyr::first(.data$signature_source), .groups = "drop")
  ranked %>%
    dplyr::group_split(.data$dataset, .data$comparison) %>%
    purrr::map_dfr(function(df) {
      stats <- df$rank_stat
      names(stats) <- df$gene
      res <- run_one_gsea(stats, term2gene, methods)
      if (!nrow(res)) return(tibble())
      meta <- df[1, c(
        "dataset", "comparison", "region", "layer", "region_level_unit", "condition_pair", "input_file",
        "left_unit", "right_unit", "left_region", "right_region", "left_condition", "right_condition", "contrast_class"
      )]
      dplyr::bind_cols(meta[rep(1, nrow(res)), ], res) %>%
        dplyr::mutate(
          matched_genes = purrr::map_chr(.data$signature, ~paste(intersect(term2gene$gene[term2gene$term == .x], names(stats)), collapse = ";")),
          signature_gene_count = lengths(strsplit(.data$matched_genes, ";", fixed = TRUE))
        ) %>%
        dplyr::left_join(signature_meta, by = "signature") %>%
        dplyr::bind_cols(purrr::map_dfr(.$matched_genes, summarise_reference_for_genes, reference_specificity = reference_specificity))
    })
}

marker_fraction_from_string <- function(x, markers) {
  genes <- normalize_id(unlist(strsplit(as.character(x), "[/;,|[:space:]]+")))
  if (!length(genes)) return(NA_real_)
  mean(genes %in% normalize_id(markers))
}

attach_neuropil_reference <- function(micro, neuropil) {
  if (!nrow(micro)) return(micro)
  if (!nrow(neuropil)) {
    return(micro %>% dplyr::mutate(
      reference_dataset = REFERENCE_DATASET,
      neuropil_reference_NES = NA_real_,
      neuropil_reference_padj = NA_real_,
      neuropil_reference_comparison = NA_character_,
      reference_match_type = "missing_neuropil_reference"
    ))
  }

  neuropil_summary <- neuropil %>%
    dplyr::group_by(.data$signature, .data$condition_pair, .data$region_level_unit) %>%
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)), .by_group = TRUE) %>%
    dplyr::summarise(
      neuropil_reference_NES = dplyr::first(.data$NES),
      neuropil_reference_padj = dplyr::first(.data$padj),
      neuropil_reference_comparison = dplyr::first(.data$comparison),
      .groups = "drop"
    )
  global_condition_summary <- neuropil %>%
    dplyr::group_by(.data$signature, .data$condition_pair) %>%
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)), .by_group = TRUE) %>%
    dplyr::summarise(
      global_neuropil_NES = dplyr::first(.data$NES),
      global_neuropil_padj = dplyr::first(.data$padj),
      global_neuropil_comparison = dplyr::first(.data$comparison),
      .groups = "drop"
    )
  global_signature_summary <- neuropil %>%
    dplyr::group_by(.data$signature) %>%
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)), .by_group = TRUE) %>%
    dplyr::summarise(
      signature_global_neuropil_NES = dplyr::first(.data$NES),
      signature_global_neuropil_padj = dplyr::first(.data$padj),
      signature_global_neuropil_comparison = dplyr::first(.data$comparison),
      .groups = "drop"
    )

  micro %>%
    dplyr::left_join(neuropil_summary, by = c("signature", "condition_pair", "region_level_unit")) %>%
    dplyr::left_join(global_condition_summary, by = c("signature", "condition_pair")) %>%
    dplyr::left_join(global_signature_summary, by = "signature") %>%
    dplyr::mutate(
      reference_dataset = REFERENCE_DATASET,
      reference_match_type = dplyr::case_when(
        !is.na(.data$neuropil_reference_NES) ~ "region_matched_neuropil",
        !is.na(.data$global_neuropil_NES) | !is.na(.data$signature_global_neuropil_NES) ~ "global_neuropil",
        TRUE ~ "missing_neuropil_reference"
      ),
      neuropil_reference_NES = dplyr::coalesce(.data$neuropil_reference_NES, .data$global_neuropil_NES, .data$signature_global_neuropil_NES),
      neuropil_reference_padj = dplyr::coalesce(.data$neuropil_reference_padj, .data$global_neuropil_padj, .data$signature_global_neuropil_padj),
      neuropil_reference_comparison = dplyr::coalesce(.data$neuropil_reference_comparison, .data$global_neuropil_comparison, .data$signature_global_neuropil_comparison)
    ) %>%
    dplyr::select(-dplyr::any_of(c(
      "global_neuropil_NES", "global_neuropil_padj", "global_neuropil_comparison",
      "signature_global_neuropil_NES", "signature_global_neuropil_padj", "signature_global_neuropil_comparison"
    )))
}

classify_signature_rows <- function(df) {
  if (!nrow(df)) return(df)
  df %>%
    dplyr::mutate(
      microglia_marker_fraction = vapply(.data$matched_genes, marker_fraction_from_string, numeric(1), markers = marker_sets$microglia),
      neuropil_marker_fraction = vapply(.data$matched_genes, marker_fraction_from_string, numeric(1), markers = marker_sets$synaptic_neuronal),
      microglia_sig = !is.na(.data$padj) & .data$padj < 0.05,
      neuropil_sig = !is.na(.data$neuropil_reference_padj) & .data$neuropil_reference_padj < 0.05,
      same_direction_as_neuropil = !is.na(.data$NES) & !is.na(.data$neuropil_reference_NES) & sign(.data$NES) == sign(.data$neuropil_reference_NES),
      stronger_in_microglia = is.na(.data$neuropil_reference_NES) | abs(.data$NES) >= abs(.data$neuropil_reference_NES) + 0.25,
      rank_based_flag = is.na(.data$padj) & !is.na(.data$pvalue) & .data$pvalue < 0.05,
      empirical_support = .data$signature_source == "empirical_microglia_vs_neuropil" | .data$signature == "empirical_microglia_enriched",
      reference_support = .data$signature_source == "reference_atlas_EWCE" |
        (!is.na(.data$reference_microglia_specificity) & .data$reference_microglia_specificity >= 0.20 & .data$reference_celltype_support == "microglia_supported"),
      signature_confidence = dplyr::case_when(
        .data$empirical_support & .data$reference_support ~ "high",
        .data$empirical_support | .data$reference_support ~ "medium",
        .data$signature_source == "curated" ~ "exploratory",
        TRUE ~ "exploratory"
      ),
      microglia_signature_class = dplyr::case_when(
        (.data$microglia_sig | .data$rank_based_flag) & .data$neuropil_sig & .data$same_direction_as_neuropil & .data$neuropil_marker_fraction >= 0.10 ~ "neuropil_shared",
        (.data$microglia_sig | .data$rank_based_flag) & .data$empirical_support & (!.data$neuropil_sig | .data$stronger_in_microglia) ~ "microglia_enriched_empirical",
        (.data$microglia_sig | .data$rank_based_flag) & .data$reference_support & (!.data$neuropil_sig | .data$stronger_in_microglia) ~ "microglia_enriched_reference_supported",
        (.data$microglia_sig | .data$rank_based_flag) & .data$signature_source == "curated" & (!.data$neuropil_sig | .data$stronger_in_microglia) ~ "curated_microglia_program",
        (.data$microglia_sig | .data$rank_based_flag) & .data$neuropil_sig & .data$same_direction_as_neuropil ~ "neuropil_shared",
        TRUE ~ "ambiguous"
      )
    )
}

split_token_field <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(character())
  out <- unlist(strsplit(as.character(x), "[/;,|[:space:]]+"), use.names = FALSE)
  out <- normalize_id(out)
  unique(out[nzchar(out)])
}

build_leading_edge_recurrence <- function(df) {
  if (!nrow(df)) return(tibble::tibble())
  sig <- df %>% dplyr::filter(!is.na(.data$padj) & .data$padj < 0.05)
  if (!nrow(sig)) return(tibble::tibble())

  long <- sig %>%
    dplyr::mutate(
      direction_sign = dplyr::case_when(
        .data$NES > 0 ~ "up",
        .data$NES < 0 ~ "down",
        TRUE ~ "neutral"
      ),
      region_pair = ifelse(!is.na(.data$left_region) & !is.na(.data$right_region), paste(.data$left_region, .data$right_region, sep = "_vs_"), NA_character_),
      proteins = purrr::map2(.data$leading_edge, .data$matched_genes, ~unique(c(split_token_field(.x), split_token_field(.y))))
    ) %>%
    dplyr::filter(lengths(.data$proteins) > 0) %>%
    tidyr::unnest_longer("proteins", values_to = "protein") %>%
    dplyr::filter(!is.na(.data$protein), nzchar(.data$protein))

  long %>%
    dplyr::group_by(
      .data$signature, .data$protein, .data$contrast_class, .data$direction_sign,
      .data$region_pair, .data$condition_pair
    ) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_comparisons = dplyr::n_distinct(.data$comparison),
      min_fdr = min(.data$padj, na.rm = TRUE),
      median_fdr = stats::median(.data$padj, na.rm = TRUE),
      mean_nes = mean(.data$NES, na.rm = TRUE),
      median_nes = stats::median(.data$NES, na.rm = TRUE),
      max_abs_nes = max(abs(.data$NES), na.rm = TRUE),
      example_comparison = dplyr::first(.data$comparison),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$n_comparisons), .data$min_fdr, dplyr::desc(abs(.data$median_nes)))
}

build_claims_ready <- function(df) {
  if (!nrow(df)) return(tibble::tibble())
  df %>%
    dplyr::filter(
      !is.na(.data$padj) & .data$padj < 0.05,
      .data$reference_celltype_support == "microglia_supported",
      .data$contrast_class %in% c("within_region_condition", "cross_region_same_condition")
    ) %>%
    dplyr::mutate(
      claim_type = dplyr::case_when(
        .data$contrast_class == "within_region_condition" ~ "within_region_condition_microglia_program",
        .data$contrast_class == "cross_region_same_condition" ~ "regional_microglia_program",
        TRUE ~ NA_character_
      ),
      claim_interpretation_note = dplyr::case_when(
        .data$contrast_class == "within_region_condition" ~ "Condition-sensitive microglia program within matched region.",
        .data$contrast_class == "cross_region_same_condition" ~ "Regional microglia specialization within matched condition.",
        TRUE ~ "Exploratory only."
      )
    ) %>%
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)))
}

write_review_xlsx <- function(all_tbl, within_tbl, cross_same_tbl, cross_cross_tbl, recurrence_tbl, claims_tbl, diagnostics_tbl) {
  if (!requireNamespace("writexl", quietly = TRUE)) return(invisible(FALSE))
  out <- file.path(PATHS$tables, "microglia_signature_enrichment_review.xlsx")
  writexl::write_xlsx(
    list(
      all = all_tbl,
      within_region_condition = within_tbl,
      cross_region_same_condition = cross_same_tbl,
      cross_region_cross_condition = cross_cross_tbl,
      leading_edge_recurrence = recurrence_tbl,
      claims_ready = claims_tbl,
      diagnostics = diagnostics_tbl
    ),
    path = out
  )
  invisible(TRUE)
}

write_empty_outputs <- function(reason, diagnostics) {
  empty <- tibble::tibble(dataset = DATASET, status = reason, run_id = RUN_ID)
  files <- c(
    "microglia_signature_enrichment.csv",
    "microglia_signature_enrichment_with_neuropil_reference.csv",
    "microglia_signature_enrichment_with_contrast_class.csv",
    "microglia_signature_enrichment_summary.csv",
    "microglia_signature_within_region_condition.csv",
    "microglia_signature_cross_region_same_condition.csv",
    "microglia_signature_cross_region_cross_condition.csv",
    "microglia_signature_unclassified_technical.csv",
    "microglia_signature_leading_edge_recurrence.csv",
    "microglia_signature_claims_ready.csv",
    "empirical_microglia_enriched_signature.csv",
    "empirical_neuropil_shared_signature.csv",
    "empirical_microglia_signature_diagnostics.csv",
    "robust_microglia_signatures.csv",
    "mixed_microenvironment_signatures.csv",
    "neuropil_shared_signatures.csv",
    "ambiguous_signatures.csv"
  )
  invisible(lapply(file.path(PATHS$tables, files), function(path) readr::write_csv(empty, path, na = "")))
  readr::write_csv(diagnostics, file.path(PATHS$tables, "microglia_signature_diagnostics.csv"))
}

diagnostics <- method_status() %>%
  dplyr::transmute(check = paste0("method_", .data$method), status = ifelse(.data$available, "PASS", "WARN"), detail = .data$note)

if (isTRUE(DRY_RUN)) {
  input_dir <- contrast_dir(DATASET)
  reference_dir <- contrast_dir(REFERENCE_DATASET)
  dry_diag <- dplyr::bind_rows(
    diagnostics,
    tibble::tibble(
      check = c("dataset", "reference_dataset", "input_dir", "reference_input_dir", "output_tables", "output_figures"),
      status = c(ifelse(DATASET == "microglia", "PASS", "WARN"), "PASS", ifelse(dir.exists(input_dir), "PASS", "WARN"), ifelse(dir.exists(reference_dir), "PASS", "WARN"), "PASS", "PASS"),
      detail = c(DATASET, REFERENCE_DATASET, input_dir, reference_dir, PATHS$tables, PATHS$figures)
    )
  )
  dry_run_line("Script", "04_differential_expression_enrichment/05_microglia_targeted_signature_enrichment.r")
  dry_run_line("Dataset", DATASET, ifelse(DATASET == "microglia", "PASS", "WARN"))
  dry_run_line("Input contrasts", input_dir, ifelse(dir.exists(input_dir), "PASS", "WARN"))
  dry_run_line("Neuropil reference contrasts", reference_dir, ifelse(dir.exists(reference_dir), "PASS", "WARN"))
  dry_run_line("Output tables", PATHS$tables)
  dry_run_line("With contrast class", file.path(PATHS$tables, "microglia_signature_enrichment_with_contrast_class.csv"))
  dry_run_line("Within-region condition", file.path(PATHS$tables, "microglia_signature_within_region_condition.csv"))
  dry_run_line("Cross-region same-condition", file.path(PATHS$tables, "microglia_signature_cross_region_same_condition.csv"))
  dry_run_line("Cross-region cross-condition", file.path(PATHS$tables, "microglia_signature_cross_region_cross_condition.csv"))
  dry_run_line("Leading-edge recurrence", file.path(PATHS$tables, "microglia_signature_leading_edge_recurrence.csv"))
  dry_run_line("Claims-ready", file.path(PATHS$tables, "microglia_signature_claims_ready.csv"))
  dry_run_line("Review workbook (optional writexl)", file.path(PATHS$tables, "microglia_signature_enrichment_review.xlsx"))
  dry_run_line("Empirical microglia signature", file.path(PATHS$tables, "empirical_microglia_enriched_signature.csv"))
  dry_run_line("Empirical neuropil-shared signature", file.path(PATHS$tables, "empirical_neuropil_shared_signature.csv"))
  readr::write_csv(dry_diag, file.path(PATHS$tables, "microglia_signature_diagnostics.csv"))
  quit(status = 0, save = "no")
}

if (DATASET != "microglia") {
  diagnostics <- dplyr::bind_rows(diagnostics, tibble::tibble(check = "dataset", status = "WARN", detail = "Skipped because dataset is not microglia."))
  write_empty_outputs("skipped_non_microglia_dataset", diagnostics)
  quit(status = 0, save = "no")
}

id_map <- read_mouse_id_map()
curated_term2gene <- expand_term2gene_for_uniprot(make_term2gene(signature_sets), id_map)

micro_ranked <- load_ranked_contrasts(DATASET)
neuropil_ranked <- load_ranked_contrasts(REFERENCE_DATASET)

empirical_signatures <- derive_empirical_signatures(micro_ranked, neuropil_ranked, id_map)
reference_specificity_result <- read_reference_specificity(id_map)
ewce_workflow_support <- reference_support_from_ewce_results(DATASET)
reference_specificity <- reference_specificity_result$specificity
if (nrow(reference_specificity) && nrow(ewce_workflow_support$gene_support)) {
  reference_specificity <- reference_specificity %>%
    dplyr::left_join(ewce_workflow_support$gene_support, by = "gene_symbol")
} else if (nrow(reference_specificity)) {
  reference_specificity$ewce_target_list_fraction <- NA_real_
  reference_specificity$ewce_target_top100_fraction <- NA_real_
}
reference_term2gene <- make_reference_term2gene(reference_specificity)

term2gene <- dplyr::bind_rows(curated_term2gene, empirical_signatures$term2gene, reference_term2gene) %>%
  dplyr::distinct()

diagnostics <- dplyr::bind_rows(
  diagnostics,
  empirical_signatures$diagnostics,
  reference_specificity_result$diagnostics,
  ewce_workflow_support$diagnostics,
  tibble::tibble(
    check = c("signature_terms", "curated_signature_terms", "empirical_signature_terms", "reference_signature_terms", "id_map_rows", "microglia_ranked_rows", "neuropil_ranked_rows"),
    status = c("PASS", "PASS", ifelse(nrow(empirical_signatures$term2gene) > 0, "PASS", "WARN"), ifelse(nrow(reference_term2gene) > 0, "PASS", "WARN"), ifelse(nrow(id_map) > 0, "PASS", "WARN"), ifelse(nrow(micro_ranked) > 0, "PASS", "WARN"), ifelse(nrow(neuropil_ranked) > 0, "PASS", "WARN")),
    detail = as.character(c(nrow(term2gene), nrow(curated_term2gene), nrow(empirical_signatures$term2gene), nrow(reference_term2gene), nrow(id_map), nrow(micro_ranked), nrow(neuropil_ranked)))
  )
)

if (!nrow(micro_ranked)) {
  write_empty_outputs("missing_microglia_ranked_contrasts", diagnostics)
  warning("No readable microglia mapped contrast tables found under: ", contrast_dir(DATASET))
  quit(status = 0, save = "no")
}

micro_enrichment <- run_signature_enrichment(micro_ranked, term2gene, reference_specificity = reference_specificity)
neuropil_enrichment <- run_signature_enrichment(neuropil_ranked, term2gene, reference_specificity = reference_specificity)
if (nrow(neuropil_enrichment)) {
  neuropil_enrichment <- neuropil_enrichment %>%
    dplyr::mutate(region = .data$region_level_unit, layer = .data$layer, region_level_unit = collapse_neuropil_layer_to_region(.data$region_level_unit))
}

if (!nrow(micro_enrichment)) {
  diagnostics <- dplyr::bind_rows(diagnostics, tibble::tibble(check = "enrichment_rows", status = "WARN", detail = "No methods produced signature enrichment rows."))
  write_empty_outputs("no_signature_enrichment_rows", diagnostics)
  quit(status = 0, save = "no")
}

with_reference <- attach_neuropil_reference(micro_enrichment, neuropil_enrichment) %>%
  classify_signature_rows() %>%
  dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)))

with_contrast <- with_reference %>%
  dplyr::mutate(
    contrast_class = dplyr::coalesce(
      .data$contrast_class,
      purrr::pmap_chr(
        list(.data$left_region, .data$right_region, .data$left_condition, .data$right_condition),
        classify_contrast_class
      )
    ),
    region_pair = ifelse(!is.na(.data$left_region) & !is.na(.data$right_region), paste(.data$left_region, .data$right_region, sep = "_vs_"), NA_character_),
    condition_pair = dplyr::coalesce(
      .data$condition_pair,
      ifelse(!is.na(.data$left_condition) & !is.na(.data$right_condition), paste(.data$left_condition, .data$right_condition, sep = "_vs_"), NA_character_)
    ),
    direction_sign = dplyr::case_when(
      .data$NES > 0 ~ "up",
      .data$NES < 0 ~ "down",
      TRUE ~ "neutral"
    )
  )

within_region_tbl <- with_contrast %>% dplyr::filter(.data$contrast_class == "within_region_condition")
cross_region_same_tbl <- with_contrast %>% dplyr::filter(.data$contrast_class == "cross_region_same_condition")
cross_region_cross_tbl <- with_contrast %>% dplyr::filter(.data$contrast_class == "cross_region_cross_condition")
unclassified_technical_tbl <- with_contrast %>% dplyr::filter(.data$contrast_class %in% c("same_region_same_condition", "unclassified"))
leading_edge_recurrence <- build_leading_edge_recurrence(with_contrast)
claims_ready <- build_claims_ready(with_contrast)

summary_tbl <- with_contrast %>%
  dplyr::count(.data$contrast_class, .data$microglia_signature_class, .data$signature, .data$reference_match_type, name = "n_contrasts") %>%
  dplyr::arrange(.data$contrast_class, .data$microglia_signature_class, dplyr::desc(.data$n_contrasts))

readr::write_csv(micro_enrichment, file.path(PATHS$tables, "microglia_signature_enrichment.csv"), na = "")
readr::write_csv(with_reference, file.path(PATHS$tables, "microglia_signature_enrichment_with_neuropil_reference.csv"), na = "")
readr::write_csv(with_contrast, file.path(PATHS$tables, "microglia_signature_enrichment_with_contrast_class.csv"), na = "")
readr::write_csv(summary_tbl, file.path(PATHS$tables, "microglia_signature_enrichment_summary.csv"), na = "")
readr::write_csv(within_region_tbl, file.path(PATHS$tables, "microglia_signature_within_region_condition.csv"), na = "")
readr::write_csv(cross_region_same_tbl, file.path(PATHS$tables, "microglia_signature_cross_region_same_condition.csv"), na = "")
readr::write_csv(cross_region_cross_tbl, file.path(PATHS$tables, "microglia_signature_cross_region_cross_condition.csv"), na = "")
readr::write_csv(unclassified_technical_tbl, file.path(PATHS$tables, "microglia_signature_unclassified_technical.csv"), na = "")
readr::write_csv(leading_edge_recurrence, file.path(PATHS$tables, "microglia_signature_leading_edge_recurrence.csv"), na = "")
readr::write_csv(claims_ready, file.path(PATHS$tables, "microglia_signature_claims_ready.csv"), na = "")
readr::write_csv(diagnostics, file.path(PATHS$tables, "microglia_signature_diagnostics.csv"), na = "")
readr::write_csv(empirical_signatures$microglia, file.path(PATHS$tables, "empirical_microglia_enriched_signature.csv"), na = "")
readr::write_csv(empirical_signatures$shared, file.path(PATHS$tables, "empirical_neuropil_shared_signature.csv"), na = "")
readr::write_csv(empirical_signatures$diagnostics, file.path(PATHS$tables, "empirical_microglia_signature_diagnostics.csv"), na = "")
readr::write_csv(with_contrast %>% dplyr::filter(.data$microglia_signature_class %in% c("microglia_enriched_empirical", "microglia_enriched_reference_supported")), file.path(PATHS$tables, "robust_microglia_signatures.csv"), na = "")
readr::write_csv(with_contrast %>% dplyr::filter(.data$microglia_signature_class == "mixed_microenvironment"), file.path(PATHS$tables, "mixed_microenvironment_signatures.csv"), na = "")
readr::write_csv(with_contrast %>% dplyr::filter(.data$microglia_signature_class == "neuropil_shared"), file.path(PATHS$tables, "neuropil_shared_signatures.csv"), na = "")
readr::write_csv(with_contrast %>% dplyr::filter(.data$microglia_signature_class == "ambiguous"), file.path(PATHS$tables, "ambiguous_signatures.csv"), na = "")

write_review_xlsx(
  all_tbl = with_contrast,
  within_tbl = within_region_tbl,
  cross_same_tbl = cross_region_same_tbl,
  cross_cross_tbl = cross_region_cross_tbl,
  recurrence_tbl = leading_edge_recurrence,
  claims_tbl = claims_ready,
  diagnostics_tbl = diagnostics
)

if (nrow(with_contrast)) {
  dot_df <- with_contrast %>%
    dplyr::mutate(
      sig_label = ifelse(!is.na(.data$padj) & .data$padj < 0.05, "FDR<0.05", "nominal/unthresholded"),
      comparison = stats::reorder(.data$comparison, .data$NES, FUN = median, na.rm = TRUE),
      signature = stats::reorder(.data$signature, abs(.data$NES), FUN = median, na.rm = TRUE),
      plot_size = pmin(6, pmax(1.5, -log10(pmax(.data$padj, .Machine$double.xmin)))),
      contrast_class = factor(
        .data$contrast_class,
        levels = c("within_region_condition", "cross_region_same_condition", "cross_region_cross_condition", "same_region_same_condition", "unclassified")
      )
    )
  p1 <- ggplot2::ggplot(dot_df, ggplot2::aes(x = .data$comparison, y = .data$signature)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$NES, size = .data$plot_size),
      shape = 21, color = "black", stroke = 0.2, alpha = 0.92
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$contrast_class), scales = "free_x", space = "free_x") +
    scale_fill_diverging(name = "NES") +
    ggplot2::scale_size_continuous(name = expression(-log[10]~FDR), range = c(1.2, 5.0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = FONT_SIZE_PANEL, grid = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      strip.placement = "outside",
      panel.spacing.x = grid::unit(1.2, "mm"),
      panel.spacing.y = grid::unit(1.0, "mm")
    )
  ggplot2::ggsave(
    file.path(PATHS$figures, "microglia_signature_dotplot.svg"),
    p1,
    width = max(FIG_WIDTH_DOUBLE, length(unique(dot_df$comparison)) * 0.22),
    height = max(FIG_HEIGHT_STANDARD, length(unique(dot_df$signature)) * 0.18),
    units = "in",
    limitsize = FALSE
  )

  scatter_df <- dot_df %>% dplyr::filter(!is.na(.data$neuropil_reference_NES))
  if (nrow(scatter_df)) {
    p2 <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = .data$neuropil_reference_NES, y = .data$NES, color = .data$microglia_signature_class)) +
      ggplot2::geom_hline(yintercept = 0, color = "#BDBDBD", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = 0, color = "#BDBDBD", linewidth = 0.3) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#7F7F7F", linewidth = 0.35) +
      ggplot2::geom_point(size = 2.2, alpha = 0.9) +
      ggplot2::stat_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.35, linetype = "solid") +
      ggplot2::scale_color_manual(values = palette_signature_class, drop = FALSE, name = "Signature class") +
      ggplot2::labs(x = "Neuropil reference NES", y = "Microglia ROI NES", color = "Class") +
      theme_manuscript(base_size = FONT_SIZE_PANEL, grid = "none") +
      ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal"
      ) +
      ggplot2::coord_equal()
    ggplot2::ggsave(file.path(PATHS$figures, "microglia_vs_neuropil_signature_NES_scatter.svg"), p2, width = FIG_WIDTH_SINGLE + 0.4, height = FIG_WIDTH_SINGLE + 0.2, units = "in")
  }

  p3_df <- with_contrast %>%
    dplyr::count(.data$microglia_signature_class, name = "n") %>%
    dplyr::mutate(microglia_signature_class = stats::reorder(.data$microglia_signature_class, .data$n))
  p3 <- ggplot2::ggplot(p3_df, ggplot2::aes(y = .data$microglia_signature_class, x = .data$n, color = .data$microglia_signature_class)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$n, yend = .data$microglia_signature_class), linewidth = 0.6, alpha = 0.85, show.legend = FALSE) +
    ggplot2::geom_point(size = 2.8, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = palette_signature_class, drop = FALSE) +
    ggplot2::labs(x = "Signature-contrast rows", y = NULL) +
    theme_manuscript(base_size = FONT_SIZE_PANEL, grid = "x")
  ggplot2::ggsave(file.path(PATHS$figures, "microglia_signature_classification_barplot.svg"), p3, width = FIG_WIDTH_SINGLE + 0.7, height = FIG_HEIGHT_SHORT + 0.8, units = "in")

  p4_df <- with_contrast %>%
    dplyr::count(.data$contrast_class, name = "n_rows") %>%
    dplyr::mutate(contrast_class = stats::reorder(.data$contrast_class, .data$n_rows))
  p4 <- ggplot2::ggplot(p4_df, ggplot2::aes(y = .data$contrast_class, x = .data$n_rows, color = .data$contrast_class)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$n_rows, yend = .data$contrast_class), linewidth = 0.7, alpha = 0.9, show.legend = FALSE) +
    ggplot2::geom_point(size = 3.0, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = palette_contrast_class, drop = FALSE) +
    ggplot2::labs(x = "Signature rows", y = NULL) +
    theme_manuscript(base_size = FONT_SIZE_PANEL, grid = "x")
  ggplot2::ggsave(file.path(PATHS$figures, "microglia_signature_contrast_class_overview.svg"), p4, width = FIG_WIDTH_SINGLE + 0.9, height = FIG_HEIGHT_SHORT + 0.9, units = "in")

  if (nrow(cross_region_same_tbl)) {
    heat_df <- cross_region_same_tbl %>%
      dplyr::group_by(.data$signature, .data$region_pair) %>%
      dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)), .by_group = TRUE) %>%
      dplyr::summarise(nes_rep = dplyr::first(.data$NES), fdr_rep = dplyr::first(.data$padj), .groups = "drop") %>%
      dplyr::mutate(
        region_pair = stats::reorder(.data$region_pair, .data$nes_rep, FUN = median, na.rm = TRUE),
        signature = stats::reorder(.data$signature, abs(.data$nes_rep), FUN = max, na.rm = TRUE),
        fdr_flag = !is.na(.data$fdr_rep) & .data$fdr_rep < 0.05
      )
    p5 <- ggplot2::ggplot(heat_df, ggplot2::aes(x = .data$region_pair, y = .data$signature, fill = .data$nes_rep)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::geom_point(data = dplyr::filter(heat_df, .data$fdr_flag), shape = 21, size = 1.4, stroke = 0.2, fill = "black", color = "black") +
      scale_fill_diverging(name = "NES") +
      ggplot2::labs(x = "Region pair", y = "Signature") +
      theme_manuscript(base_size = FONT_SIZE_PANEL, grid = "none") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom",
        legend.direction = "horizontal"
      )
    ggplot2::ggsave(file.path(PATHS$figures, "microglia_signature_region_program_heatmap.svg"), p5, width = max(FIG_WIDTH_DOUBLE - 0.2, length(unique(heat_df$region_pair)) * 0.42), height = max(FIG_HEIGHT_STANDARD, length(unique(heat_df$signature)) * 0.19), units = "in", limitsize = FALSE)
  }

  if (nrow(within_region_tbl)) {
    within_plot_df <- within_region_tbl %>%
      dplyr::mutate(
        comparison = stats::reorder(.data$comparison, .data$NES, FUN = median, na.rm = TRUE),
        signature = stats::reorder(.data$signature, abs(.data$NES), FUN = median, na.rm = TRUE),
        plot_size = pmin(6, pmax(1.5, -log10(pmax(.data$padj, .Machine$double.xmin)))),
        left_region = factor(.data$left_region, levels = c("CA1", "CA2", "CA3", "DG"))
      )
    p6 <- ggplot2::ggplot(within_plot_df, ggplot2::aes(x = .data$comparison, y = .data$signature)) +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data$NES, size = .data$plot_size),
        shape = 21, color = "black", stroke = 0.22, alpha = 0.95
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$left_region), scales = "free_x", space = "free_x") +
      scale_fill_diverging(name = "NES") +
      ggplot2::scale_size_continuous(name = expression(-log[10]~FDR), range = c(1.3, 5.2)) +
      ggplot2::labs(x = NULL, y = NULL) +
      theme_manuscript(base_size = FONT_SIZE_PANEL, grid = "none") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.placement = "outside",
        panel.spacing.x = grid::unit(1.0, "mm"),
        panel.spacing.y = grid::unit(1.0, "mm")
      )
    ggplot2::ggsave(file.path(PATHS$figures, "microglia_signature_within_region_condition_dotplot.svg"), p6, width = max(FIG_WIDTH_DOUBLE - 0.4, length(unique(within_plot_df$comparison)) * 0.26), height = max(FIG_HEIGHT_STANDARD, length(unique(within_plot_df$signature)) * 0.18), units = "in", limitsize = FALSE)
  }

  if (nrow(leading_edge_recurrence)) {
    rec_plot <- leading_edge_recurrence %>%
      dplyr::group_by(.data$protein) %>%
      dplyr::summarise(
        total_recurrence = sum(.data$n_comparisons, na.rm = TRUE),
        dominant_direction = names(sort(table(.data$direction_sign), decreasing = TRUE))[1],
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(.data$total_recurrence)) %>%
      dplyr::slice_head(n = 18) %>%
      dplyr::mutate(protein = stats::reorder(.data$protein, .data$total_recurrence))
    p7 <- ggplot2::ggplot(rec_plot, ggplot2::aes(y = .data$protein, x = .data$total_recurrence, color = .data$dominant_direction)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$total_recurrence, yend = .data$protein), linewidth = 0.7, alpha = 0.9, show.legend = FALSE) +
      ggplot2::geom_point(size = 2.8, show.legend = TRUE) +
      ggplot2::scale_color_manual(values = palette_direction, drop = FALSE, name = "Direction") +
      ggplot2::labs(x = "Recurrence across significant rows", y = "Leading-edge protein") +
      theme_manuscript(base_size = FONT_SIZE_PANEL, grid = "x") +
      ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal"
      )
    ggplot2::ggsave(file.path(PATHS$figures, "microglia_signature_leading_edge_recurrence.svg"), p7, width = FIG_WIDTH_SINGLE + 1.1, height = FIG_HEIGHT_STANDARD, units = "in")
  }
}

writeLines(c(
  "Microglia-targeted signature enrichment",
  "",
  paste0("Dataset: ", DATASET),
  paste0("Neuropil reference dataset: ", REFERENCE_DATASET),
  paste0("Run ID: ", RUN_ID),
  "",
  "Method note:",
  "All detected proteins in each mapped contrast table are used as the ranked universe.",
  "Neuropil is used as a reference annotation layer; no intensities or logFC values are subtracted.",
  "Microglia region-only contrasts are compared to neuropil references after collapsing neuropil layer units to parent regions where possible.",
  "Contrast metadata is parsed explicitly (left/right unit, region, condition, and contrast_class).",
  "Condition_pair alone is not interpreted as a stress/condition effect without contrast_class context.",
  "",
  "Signature evidence:",
  "Curated microglia programs are retained as exploratory/manual biological programs.",
  "Empirical signatures are derived from proteins repeatedly high-ranking in microglia ROI contrasts relative to neuropil reference contrasts.",
  "Reference-atlas signatures and per-row specificity columns are imported from EWCE/ewceData when available and skipped softly otherwise.",
  "Classification prioritizes empirical and reference-backed evidence over curated-only programs.",
  "Claims-ready rows are restricted to within_region_condition and cross_region_same_condition contrasts with FDR<0.05 and microglia-supported reference specificity.",
  "Cross_region_cross_condition contrasts are retained as exploratory/confounded outputs and excluded from claims-ready exports."
), file.path(PATHS$reports, "microglia_signature_methods_note.txt"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(microglia_contrasts = contrast_dir(DATASET), neuropil_contrasts = contrast_dir(REFERENCE_DATASET)),
  outputs = list(tables = PATHS$tables, figures = PATHS$figures),
  parameters = list(dataset = DATASET, reference_dataset = REFERENCE_DATASET, run_id = RUN_ID),
  notes = "Microglia signatures are interpreted for ROI-local microenvironment samples. Empirical and EWCE/reference-backed evidence is prioritized over curated-only programs; neuropil remains an annotation reference, not a subtraction baseline."
)

message("Microglia-targeted signature enrichment completed.")
message("Annotated table: ", file.path(PATHS$tables, "microglia_signature_enrichment_with_neuropil_reference.csv"))
message("Contrast-class table: ", file.path(PATHS$tables, "microglia_signature_enrichment_with_contrast_class.csv"))
message("Claims-ready table: ", file.path(PATHS$tables, "microglia_signature_claims_ready.csv"))
