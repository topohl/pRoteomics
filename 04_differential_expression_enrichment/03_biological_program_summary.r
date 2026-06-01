#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "enrichment_plots.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}
dataset_cli <- arg_value("--dataset", default = "")
if (nzchar(dataset_cli)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset_cli, source = "--dataset"))
DATASET <- current_dataset()

MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- file.path("biological_program_summary", DATASET)
PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

required_pkgs <- c("dplyr", "tidyr", "readr", "stringr", "tibble", "ggplot2")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

latest_manifest <- function() {
  first_existing_path(c(
    path_processed("04_differential_expression_enrichment", "compareGO", DATASET, "compareGO_input_manifest.csv"),
    latest_file(path_processed("04_differential_expression_enrichment", "compareGO", DATASET), "^compareGO_input_manifest.*\\.csv$"),
    latest_file(path_processed("04_differential_expression_enrichment", "clusterProfiler", DATASET), "^clusterProfiler_manifest.*\\.csv$"),
    latest_file(path_results("reports", "04_differential_expression_enrichment", "clusterProfiler", DATASET), "^clusterProfiler_manifest.*\\.csv$")
  ))
}

manifest_file <- latest_manifest()

if (is_dry_run()) {
  dry_run_line("Script", "04_differential_expression_enrichment/03_biological_program_summary.r")
  dry_run_line("Dataset", DATASET)
  dry_run_line("Manifest", manifest_file, if (file.exists(manifest_file)) "PASS" else "FAIL")
  dry_run_line("Program summary", file.path(PATHS$tables, "program_summary.csv"))
  dry_run_line("Program summary wide", file.path(PATHS$tables, "program_summary_wide.csv"))
  dry_run_line("Program heatmap-ready table", file.path(PATHS$tables, "program_summary_heatmap_ready.csv"))
  dry_run_line("Program evidence", file.path(PATHS$source_data, "program_term_gene_evidence.csv"))
  dry_run_line("Neuropil-annotated summary", file.path(PATHS$tables, "program_summary_neuropil_annotated.csv"))
  dry_run_line("Microglia signature-annotated summary", file.path(PATHS$tables, "program_summary_microglia_signature_annotated.csv"))
  dry_run_line("Integrated interpretation", file.path(PATHS$tables, "program_summary_integrated_interpretation.csv"))
  dry_run_line("Program atlas heatmap", file.path(PATHS$figures, "program_atlas_heatmap.svg"))
  quit(status = if (file.exists(manifest_file)) 0 else 1, save = "no")
}

if (!file.exists(manifest_file)) {
  stop("No compareGO/clusterProfiler manifest found for dataset: ", DATASET, call. = FALSE)
}

manifest <- readr::read_csv(manifest_file, show_col_types = FALSE)
if ("dataset" %in% names(manifest)) manifest <- dplyr::filter(manifest, .data$dataset == DATASET)
if (!"output_table" %in% names(manifest)) stop("Manifest lacks output_table column: ", manifest_file, call. = FALSE)
path_qc <- validate_manifest_paths(manifest, path_cols = "output_table", allow_missing = TRUE)
readr::write_csv(path_qc, file.path(PATHS$logs, "program_summary_manifest_path_qc.csv"), na = "")

term_tables <- lapply(seq_len(nrow(manifest)), function(i) {
  table_path <- as.character(manifest$output_table[[i]])
  if (is.na(table_path) || !file.exists(table_path)) return(NULL)
  tbl <- tryCatch(readr::read_csv(table_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(tbl) || !"Description" %in% names(tbl)) return(NULL)
  tbl$dataset <- DATASET
  tbl$comparison <- if ("comparison" %in% names(manifest)) manifest$comparison[[i]] else if ("contrast" %in% names(manifest)) manifest$contrast[[i]] else NA_character_
  tbl$route_category <- if ("route_category" %in% names(manifest)) manifest$route_category[[i]] else NA_character_
  tbl$route_unit <- if ("route_unit" %in% names(manifest)) manifest$route_unit[[i]] else NA_character_
  tbl$result_type <- if ("result_type" %in% names(manifest)) manifest$result_type[[i]] else NA_character_
  tbl$source_file <- table_path
  tbl
})
terms <- dplyr::bind_rows(term_tables)

if (!nrow(terms)) {
  status <- tibble::tibble(dataset = DATASET, status = "skipped", reason = "No readable enrichment output tables from manifest.")
  readr::write_csv(status, file.path(PATHS$tables, "program_summary.csv"), na = "")
  readr::write_csv(status, file.path(PATHS$tables, "program_summary_wide.csv"), na = "")
  readr::write_csv(status, file.path(PATHS$tables, "program_summary_heatmap_ready.csv"), na = "")
  readr::write_csv(status, file.path(PATHS$source_data, "program_term_gene_evidence.csv"), na = "")
  write_run_manifest(
    file.path(PATHS$logs, "run_manifest.yml"),
    inputs = list(manifest = manifest_file),
    outputs = list(
      program_summary = file.path(PATHS$tables, "program_summary.csv"),
      program_summary_wide = file.path(PATHS$tables, "program_summary_wide.csv"),
      heatmap_ready = file.path(PATHS$tables, "program_summary_heatmap_ready.csv"),
      evidence = file.path(PATHS$source_data, "program_term_gene_evidence.csv")
    ),
    parameters = list(dataset = DATASET, program_patterns = biological_program_patterns()),
    notes = "Skipped: no readable enrichment output tables from manifest."
  )
  quit(status = 0, save = "no")
}

if (!"pvalue" %in% names(terms)) terms$pvalue <- NA_real_
if (!"p.value" %in% names(terms)) terms$`p.value` <- NA_real_

terms <- map_terms_to_programs(terms, "Description") %>%
  dplyr::filter(!is.na(.data$biological_program)) %>%
  dplyr::mutate(
    NES = if ("NES" %in% names(.)) suppressWarnings(as.numeric(.data$NES)) else NA_real_,
    p.adjust = if ("p.adjust" %in% names(.)) suppressWarnings(as.numeric(.data$p.adjust)) else NA_real_,
    FDR = .data$p.adjust,
    core_genes = if ("core_enrichment" %in% names(.)) as.character(.data$core_enrichment) else NA_character_,
    raw_p = dplyr::coalesce(
      suppressWarnings(as.numeric(.data$pvalue)),
      suppressWarnings(as.numeric(.data$p.value))
    )
  )

split_gene_tokens <- function(x) {
  x <- stats::na.omit(as.character(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(character())
  genes <- unlist(strsplit(paste(x, collapse = "/"), "[/;,[:space:]]+"))
  genes <- unique(genes[nzchar(genes)])
  genes
}

frequent_genes <- function(x, n = 12L) {
  genes <- unlist(lapply(x, split_gene_tokens), use.names = FALSE)
  genes <- genes[nzchar(genes)]
  if (!length(genes)) return(NA_character_)
  tab <- sort(table(genes), decreasing = TRUE)
  paste(names(tab)[seq_len(min(n, length(tab)))], collapse = ";")
}

evidence <- terms %>%
  dplyr::transmute(
    dataset,
    comparison,
    route_category,
    route_unit,
    result_type,
    biological_program,
    ID = if ("ID" %in% names(terms)) as.character(.data$ID) else NA_character_,
    Description = as.character(.data$Description),
    NES,
    effect_direction = dplyr::case_when(
      is.na(.data$NES) ~ "undirected",
      .data$NES > 0 ~ "positive_NES",
      .data$NES < 0 ~ "negative_NES",
      TRUE ~ "neutral"
    ),
    raw_p,
    p.adjust,
    FDR,
    core_genes,
    source_file
  )

program_summary <- evidence %>%
  dplyr::group_by(.data$dataset, .data$comparison, .data$route_category, .data$route_unit, .data$biological_program) %>%
  dplyr::arrange(.data$FDR, dplyr::desc(abs(.data$NES)), .by_group = TRUE) %>%
  dplyr::summarise(
    n_terms = dplyr::n(),
    min_fdr = suppressWarnings(min(.data$FDR, na.rm = TRUE)),
    min_raw_p = suppressWarnings(min(.data$raw_p, na.rm = TRUE)),
    representative_NES = dplyr::first(.data$NES),
    effect_direction = dplyr::first(.data$effect_direction),
    top_term = dplyr::first(.data$Description),
    key_genes = frequent_genes(.data$core_genes),
    source_file = dplyr::first(.data$source_file),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    min_fdr = ifelse(is.infinite(.data$min_fdr), NA_real_, .data$min_fdr),
    min_raw_p = ifelse(is.infinite(.data$min_raw_p), NA_real_, .data$min_raw_p),
    direction = dplyr::case_when(
      !is.na(.data$effect_direction) ~ .data$effect_direction,
      is.na(.data$representative_NES) ~ "undirected",
      .data$representative_NES > 0 ~ "positive_NES",
      .data$representative_NES < 0 ~ "negative_NES",
      TRUE ~ "neutral"
    )
  ) %>%
  dplyr::arrange(.data$min_fdr, .data$dataset, .data$comparison, .data$biological_program)

readr::write_csv(program_summary, file.path(PATHS$tables, "program_summary.csv"), na = "")
readr::write_csv(evidence, file.path(PATHS$source_data, "program_term_gene_evidence.csv"), na = "")

heatmap_ready <- program_summary %>%
  dplyr::mutate(
    comparison_label = paste(.data$comparison, .data$route_category, .data$route_unit, sep = " | "),
    signed_neg_log10_fdr = sign(dplyr::coalesce(.data$representative_NES, 0)) *
      -log10(pmax(dplyr::coalesce(.data$min_fdr, 1), .Machine$double.xmin))
  ) %>%
  dplyr::select(
    dataset, comparison, route_category, route_unit, comparison_label,
    biological_program, direction, representative_NES, min_raw_p, min_fdr,
    signed_neg_log10_fdr, n_terms, top_term, key_genes, source_file
  )
readr::write_csv(heatmap_ready, file.path(PATHS$tables, "program_summary_heatmap_ready.csv"), na = "")

program_wide <- heatmap_ready %>%
  dplyr::mutate(cell_value = paste0(
    direction,
    "; NES=", ifelse(is.na(representative_NES), "NA", signif(representative_NES, 3)),
    "; FDR=", ifelse(is.na(min_fdr), "NA", signif(min_fdr, 3)),
    "; genes=", dplyr::coalesce(key_genes, "NA")
  )) %>%
  dplyr::select(dataset, comparison_label, biological_program, cell_value) %>%
  tidyr::pivot_wider(names_from = biological_program, values_from = cell_value)
readr::write_csv(program_wide, file.path(PATHS$tables, "program_summary_wide.csv"), na = "")

optional_read_csv <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
}

mode_value <- function(x) {
  x <- stats::na.omit(as.character(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[[1]]
}

latest_neuropil_annotation <- file.path(
  path_results("tables", MODULE_ID, "neuropil_contamination_annotation", DATASET),
  "microglia_neuropil_annotation_latest.csv"
)
latest_microglia_signature <- file.path(
  path_results("tables", MODULE_ID, "microglia_targeted_signature_enrichment", DATASET),
  "microglia_signature_enrichment_with_neuropil_reference.csv"
)

neuropil_annotation <- optional_read_csv(latest_neuropil_annotation)
signature_annotation <- optional_read_csv(latest_microglia_signature)

program_summary_neuropil <- program_summary
if (!is.null(neuropil_annotation) && nrow(neuropil_annotation) && "comparison" %in% names(neuropil_annotation)) {
  neuropil_by_comparison <- neuropil_annotation %>%
    dplyr::group_by(.data$comparison) %>%
    dplyr::summarise(
      interpretation_class = mode_value(.data$interpretation_class),
      gene_overlap_fraction = suppressWarnings(mean(.data$gene_overlap_fraction, na.rm = TRUE)),
      gene_jaccard = suppressWarnings(mean(.data$gene_jaccard, na.rm = TRUE)),
      neuropil_marker_fraction = suppressWarnings(mean(.data$neuropil_marker_fraction, na.rm = TRUE)),
      microglia_marker_fraction = suppressWarnings(mean(.data$microglia_marker_fraction, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ifelse(is.nan(.x), NA_real_, .x)))
  program_summary_neuropil <- program_summary %>%
    dplyr::left_join(neuropil_by_comparison, by = "comparison")
} else {
  program_summary_neuropil <- program_summary_neuropil %>%
    dplyr::mutate(
      interpretation_class = NA_character_,
      gene_overlap_fraction = NA_real_,
      gene_jaccard = NA_real_,
      neuropil_marker_fraction = NA_real_,
      microglia_marker_fraction = NA_real_
    )
}

program_summary_signature <- program_summary
if (!is.null(signature_annotation) && nrow(signature_annotation) && "comparison" %in% names(signature_annotation)) {
  signature_by_comparison <- signature_annotation %>%
    dplyr::group_by(.data$comparison) %>%
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)), .by_group = TRUE) %>%
    dplyr::summarise(
      microglia_signature_class = mode_value(.data$microglia_signature_class),
      top_microglia_signature = dplyr::first(.data$signature),
      microglia_signature_NES = dplyr::first(.data$NES),
      neuropil_reference_NES = dplyr::first(.data$neuropil_reference_NES),
      reference_match_type = mode_value(.data$reference_match_type),
      .groups = "drop"
    )
  program_summary_signature <- program_summary %>%
    dplyr::left_join(signature_by_comparison, by = "comparison")
} else {
  program_summary_signature <- program_summary_signature %>%
    dplyr::mutate(
      microglia_signature_class = NA_character_,
      top_microglia_signature = NA_character_,
      microglia_signature_NES = NA_real_,
      neuropil_reference_NES = NA_real_,
      reference_match_type = NA_character_
    )
}

program_summary_integrated <- program_summary_neuropil %>%
  dplyr::left_join(
    program_summary_signature %>%
      dplyr::select(
        dataset, comparison, route_category, route_unit, biological_program,
        microglia_signature_class, top_microglia_signature,
        microglia_signature_NES, neuropil_reference_NES, reference_match_type
      ),
    by = c("dataset", "comparison", "route_category", "route_unit", "biological_program")
  ) %>%
  dplyr::mutate(
    integrated_interpretation = dplyr::case_when(
      .data$microglia_signature_class == "microglia_enriched" & !.data$interpretation_class %in% c("neuropil_sensitive", "neuropil_marker_enriched") ~ "microglia_supported_program",
      .data$microglia_signature_class == "neuropil_shared" | .data$interpretation_class %in% c("neuropil_sensitive", "neuropil_marker_enriched") ~ "neuropil_shared_or_sensitive_program",
      .data$microglia_signature_class == "mixed_microenvironment" | .data$interpretation_class == "mixed_microenvironment" ~ "mixed_microenvironment_program",
      is.na(.data$microglia_signature_class) & is.na(.data$interpretation_class) ~ "unannotated_program",
      TRUE ~ "ambiguous_program"
    ),
    interpretation_note = dplyr::case_when(
      .data$integrated_interpretation == "microglia_supported_program" ~ "Supported by microglia-targeted signatures with weak/absent or weaker neuropil reference signal.",
      .data$integrated_interpretation == "neuropil_shared_or_sensitive_program" ~ "Shared with or sensitive to neuropil reference; do not interpret as microglia-intrinsic without orthogonal support.",
      .data$integrated_interpretation == "mixed_microenvironment_program" ~ "Present in microglia ROI and neuropil reference in a pattern consistent with local microenvironment biology.",
      .data$integrated_interpretation == "unannotated_program" ~ "No optional neuropil or microglia signature annotation was available.",
      TRUE ~ "Weak, discordant, broad, or partially missing annotation support."
    )
  )

readr::write_csv(program_summary_neuropil, file.path(PATHS$tables, "program_summary_neuropil_annotated.csv"), na = "")
readr::write_csv(program_summary_signature, file.path(PATHS$tables, "program_summary_microglia_signature_annotated.csv"), na = "")
readr::write_csv(program_summary_integrated, file.path(PATHS$tables, "program_summary_integrated_interpretation.csv"), na = "")

plot_df <- heatmap_ready %>% dplyr::rename(score = "signed_neg_log10_fdr")
if (nrow(plot_df)) {
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$comparison_label, y = .data$biological_program, fill = .data$score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", na.value = "grey90", name = "signed -log10 FDR") +
    ggplot2::labs(x = NULL, y = NULL, title = paste("Biological program atlas:", DATASET)) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
  save_plot_dual(p, file.path(PATHS$figures, "program_atlas_heatmap.svg"), width = max(6, length(unique(plot_df$comparison_label)) * 0.35), height = 4.8)
}

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(manifest = manifest_file, enrichment_tables = unique(evidence$source_file)),
  outputs = list(
    program_summary = file.path(PATHS$tables, "program_summary.csv"),
    program_summary_wide = file.path(PATHS$tables, "program_summary_wide.csv"),
    heatmap_ready = file.path(PATHS$tables, "program_summary_heatmap_ready.csv"),
    evidence = file.path(PATHS$source_data, "program_term_gene_evidence.csv"),
    heatmap = file.path(PATHS$figures, "program_atlas_heatmap.svg")
  ),
  parameters = list(dataset = DATASET, program_patterns = biological_program_patterns()),
  notes = "Program mapping uses conservative regex on enrichment Description; claims remain exploratory unless independently supported."
)
