#!/usr/bin/env Rscript
# ================================================================
# Script: 04_differential_expression_enrichment/04_neuropil_reference_annotation.r
# Stage: enrichment
# Scope: dataset_specific
# Consumes: required data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv; data/processed/04_differential_expression_enrichment/clusterProfiler/neuron_neuropil/clusterProfiler_manifest.csv; optional config/marker_panels/wgcna_reference_marker_sets.csv.
# Produces: results/tables/04_differential_expression_enrichment/neuropil_reference_annotation/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Runs after neuron_neuropil enrichment exists.
# ================================================================

# Neuropil reference annotation workflow
#
# This script annotates microglia-enriched GO/GSEA results using a separately
# processed neuropil dataset as a reference layer. It does not subtract raw
# intensities or logFC values, because microglia and neuropil matrices may have
# been normalized and imputed independently. The intended output is an
# interpretation/sensitivity layer for manuscript use:
#   - microglia_robust
#   - mixed_microenvironment
#   - neuropil_sensitive
#   - ambiguous

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))

MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- "neuropil_reference_annotation"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

DATASET <- current_dataset_from_cli()
REFERENCE_DATASET <- Sys.getenv("PROTEOMICS_NEUROPIL_REFERENCE_DATASET", unset = "neuron_neuropil")
REFERENCE_DATASET <- validate_dataset(REFERENCE_DATASET, source = "PROTEOMICS_NEUROPIL_REFERENCE_DATASET")
RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
DRY_RUN <- is_dry_run()

CANONICAL_PATHS <- lapply(CANONICAL_PATHS, function(path) file.path(path, DATASET))
invisible(lapply(CANONICAL_PATHS, dir_create))

message("Neuropil reference annotation")
message("Dataset: ", DATASET)
message("Neuropil reference dataset: ", REFERENCE_DATASET)
message("Dry run: ", DRY_RUN)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

required_pkgs <- c("dplyr", "readr", "tidyr", "stringr", "purrr", "tibble", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0 && !isTRUE(DRY_RUN)) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}
if (length(missing_pkgs) == 0) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
    library(purrr)
    library(tibble)
    library(ggplot2)
  })
}

marker_sets <- list(
  microglia = c(
    "AIF1", "P2RY12", "TMEM119", "CX3CR1", "CSF1R", "TYROBP", "HEXB",
    "C1QA", "C1QB", "C1QC", "ITGAM", "SPI1", "TREM2", "LAPTM5"
  ),
  neuropil_synaptic_neuronal = c(
    "SYN1", "SYP", "SNAP25", "STX1A", "STXBP1", "DLG4", "CAMK2A", "CAMK2B",
    "MAP2", "NEFL", "NEFM", "RBFOX3", "TUBB3", "GRIN1", "GRIA1", "VAMP2"
  ),
  astrocyte = c("GFAP", "ALDH1L1", "AQP4", "SLC1A3", "SLC1A2", "GLUL", "GJA1"),
  oligodendrocyte_myelin = c("MBP", "PLP1", "MAG", "MOG", "MOBP", "CNP", "CLDN11"),
  endothelial_pericyte = c("PECAM1", "CLDN5", "KDR", "FLT1", "PDGFRB", "RGS5", "ACTA2")
)

normalize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^ +| +$", "", x)
  toupper(x)
}

split_gene_field <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(character(0))
  genes <- unlist(strsplit(as.character(x), "[/;,|[:space:]]+"), use.names = FALSE)
  genes <- normalize_id(genes)
  unique(genes[nzchar(genes)])
}

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) NA_character_ else hit
}

safe_read_csv <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) {
      warning("Failed to read CSV: ", path, " | ", e$message)
      NULL
    }
  )
}

manifest_path <- function(dataset) {
  path_processed(
    MODULE_ID,
    "clusterProfiler",
    dataset,
    "clusterProfiler_manifest.csv"
  )
}

load_manifest <- function(dataset) {
  path <- manifest_path(dataset)
  manifest <- safe_read_csv(path)
  if (is.null(manifest)) {
    return(list(path = path, data = tibble(), status = "missing"))
  }
  if (!"dataset" %in% colnames(manifest)) manifest$dataset <- dataset
  list(path = path, data = manifest, status = "loaded")
}

resolve_table_path <- function(path) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(path)) return(NA_character_)
  path <- as.character(path)
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = FALSE))
  candidate <- repo_path(path)
  if (file.exists(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

make_empty_annotation <- function(reason) {
  tibble(
    dataset = DATASET,
    reference_dataset = REFERENCE_DATASET,
    status = reason,
    run_id = RUN_ID,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}

extract_enrichment_tables <- function(manifest, dataset_label) {
  if (nrow(manifest) == 0) return(tibble())

  table_col <- pick_col(manifest, c("output_table", "result_table", "table_path", "path"))
  if (is.na(table_col)) {
    warning("Manifest has no output table column for dataset: ", dataset_label)
    return(tibble())
  }

  comparison_col <- pick_col(manifest, c("comparison", "comparison_name", "Comparison"))
  result_type_col <- pick_col(manifest, c("result_type", "type", "analysis_type"))
  ontology_col <- pick_col(manifest, c("ontology", "ONTOLOGY"))
  route_category_col <- pick_col(manifest, c("route_category", "category"))
  route_unit_col <- pick_col(manifest, c("route_unit", "unit_folder", "unit"))

  rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    table_path_i <- resolve_table_path(manifest[[table_col]][[i]])
    tab <- safe_read_csv(table_path_i)
    if (is.null(tab) || nrow(tab) == 0) next

    if (!"ID" %in% colnames(tab) && "id" %in% colnames(tab)) tab$ID <- tab$id
    if (!"Description" %in% colnames(tab) && "description" %in% colnames(tab)) tab$Description <- tab$description

    if (!"ID" %in% colnames(tab) && !"Description" %in% colnames(tab)) next

    tab$.source_dataset <- dataset_label
    tab$.source_table <- table_path_i
    tab$.manifest_row <- i
    tab$.comparison <- if (!is.na(comparison_col)) as.character(manifest[[comparison_col]][[i]]) else NA_character_
    tab$.result_type <- if (!is.na(result_type_col)) as.character(manifest[[result_type_col]][[i]]) else NA_character_
    tab$.ontology_manifest <- if (!is.na(ontology_col)) as.character(manifest[[ontology_col]][[i]]) else NA_character_
    tab$.route_category <- if (!is.na(route_category_col)) as.character(manifest[[route_category_col]][[i]]) else NA_character_
    tab$.route_unit <- if (!is.na(route_unit_col)) as.character(manifest[[route_unit_col]][[i]]) else NA_character_
    rows[[i]] <- tab
  }

  dplyr::bind_rows(rows)
}

prepare_terms <- function(df) {
  if (nrow(df) == 0) return(tibble())

  gene_col <- pick_col(df, c("core_enrichment", "geneID", "gene_id", "genes", "leading_edge"))
  padj_col <- pick_col(df, c("p.adjust", "p_adj", "padj", "qvalue"))
  nes_col <- pick_col(df, c("NES", "nes", "enrichmentScore", "score"))
  ontology_col <- pick_col(df, c("ONTOLOGY", "ontology", ".ontology_manifest"))

  df %>%
    mutate(
      term_id = if ("ID" %in% colnames(.)) as.character(.data$ID) else NA_character_,
      term_description = if ("Description" %in% colnames(.)) as.character(.data$Description) else term_id,
      term_key = ifelse(!is.na(term_id) & nzchar(term_id), term_id, term_description),
      comparison = .data$.comparison,
      result_type = .data$.result_type,
      ontology = if (!is.na(ontology_col)) as.character(.data[[ontology_col]]) else NA_character_,
      route_category = .data$.route_category,
      route_unit = .data$.route_unit,
      p_adjust = if (!is.na(padj_col)) suppressWarnings(as.numeric(.data[[padj_col]])) else NA_real_,
      score = if (!is.na(nes_col)) suppressWarnings(as.numeric(.data[[nes_col]])) else NA_real_,
      gene_string = if (!is.na(gene_col)) as.character(.data[[gene_col]]) else NA_character_
    ) %>%
    rowwise() %>%
    mutate(
      genes = list(split_gene_field(gene_string)),
      n_genes = length(genes)
    ) %>%
    ungroup() %>%
    select(
      source_dataset = .data$.source_dataset,
      source_table = .data$.source_table,
      comparison, result_type, ontology, route_category, route_unit,
      term_id, term_description, term_key, p_adjust, score, gene_string, genes, n_genes
    )
}

marker_fraction <- function(genes, marker_vector) {
  genes <- normalize_id(genes)
  marker_vector <- normalize_id(marker_vector)
  if (length(genes) == 0) return(NA_real_)
  sum(genes %in% marker_vector) / length(genes)
}

marker_hits <- function(genes, marker_vector) {
  genes <- normalize_id(genes)
  marker_vector <- normalize_id(marker_vector)
  paste(sort(intersect(genes, marker_vector)), collapse = ";")
}

term_overlap_fraction <- function(a, b) {
  a <- normalize_id(a)
  b <- normalize_id(b)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(unique(a))
}

jaccard <- function(a, b) {
  a <- normalize_id(a)
  b <- normalize_id(b)
  u <- union(a, b)
  if (length(u) == 0) return(0)
  length(intersect(a, b)) / length(u)
}

classify_term <- function(microglia_marker_fraction,
                          neuropil_marker_fraction,
                          neuropil_term_present,
                          gene_overlap_fraction,
                          same_direction,
                          microglia_padj,
                          neuropil_padj) {
  microglia_sig <- !is.na(microglia_padj) && microglia_padj < 0.05
  neuropil_sig <- !is.na(neuropil_padj) && neuropil_padj < 0.05

  if (isTRUE(microglia_sig) &&
      !isTRUE(neuropil_sig) &&
      !is.na(microglia_marker_fraction) && microglia_marker_fraction >= 0.15 &&
      (is.na(neuropil_marker_fraction) || neuropil_marker_fraction < 0.15) &&
      gene_overlap_fraction < 0.20) {
    return("microglia_robust")
  }

  if ((isTRUE(neuropil_sig) || isTRUE(neuropil_term_present)) &&
      !is.na(neuropil_marker_fraction) && neuropil_marker_fraction >= 0.15 &&
      gene_overlap_fraction >= 0.30 &&
      isTRUE(same_direction)) {
    return("neuropil_sensitive")
  }

  if ((isTRUE(neuropil_sig) || isTRUE(neuropil_term_present)) && gene_overlap_fraction >= 0.20) {
    return("mixed_microenvironment")
  }

  if (!is.na(neuropil_marker_fraction) && neuropil_marker_fraction >= 0.25) {
    return("neuropil_marker_enriched")
  }

  "ambiguous"
}

find_best_neuropil_match <- function(m_row, neuropil_terms) {
  if (nrow(neuropil_terms) == 0) {
    return(tibble(
      neuropil_term_present = FALSE,
      neuropil_comparison = NA_character_,
      neuropil_p_adjust = NA_real_,
      neuropil_score = NA_real_,
      neuropil_source_table = NA_character_,
      gene_overlap_fraction = 0,
      gene_jaccard = 0,
      overlapping_genes = ""
    ))
  }

  same_term <- neuropil_terms %>%
    filter(.data$term_key == m_row$term_key)

  if (nrow(same_term) == 0) {
    return(tibble(
      neuropil_term_present = FALSE,
      neuropil_comparison = NA_character_,
      neuropil_p_adjust = NA_real_,
      neuropil_score = NA_real_,
      neuropil_source_table = NA_character_,
      gene_overlap_fraction = 0,
      gene_jaccard = 0,
      overlapping_genes = ""
    ))
  }

  genes_m <- m_row$genes[[1]]
  scored <- same_term %>%
    rowwise() %>%
    mutate(
      gene_overlap_fraction = term_overlap_fraction(genes_m, genes),
      gene_jaccard = jaccard(genes_m, genes),
      overlapping_genes = paste(sort(intersect(normalize_id(genes_m), normalize_id(genes))), collapse = ";")
    ) %>%
    ungroup() %>%
    arrange(
      desc(!is.na(.data$p_adjust) & .data$p_adjust < 0.05),
      desc(.data$gene_jaccard),
      .data$p_adjust
    )

  best <- scored[1, , drop = FALSE]
  tibble(
    neuropil_term_present = TRUE,
    neuropil_comparison = best$comparison[[1]],
    neuropil_p_adjust = best$p_adjust[[1]],
    neuropil_score = best$score[[1]],
    neuropil_source_table = best$source_table[[1]],
    gene_overlap_fraction = best$gene_overlap_fraction[[1]],
    gene_jaccard = best$gene_jaccard[[1]],
    overlapping_genes = best$overlapping_genes[[1]]
  )
}

annotate_microglia_terms <- function(microglia_terms, neuropil_terms) {
  if (nrow(microglia_terms) == 0) return(tibble())

  annotations <- vector("list", nrow(microglia_terms))
  for (i in seq_len(nrow(microglia_terms))) {
    m <- microglia_terms[i, , drop = FALSE]
    best <- find_best_neuropil_match(m, neuropil_terms)
    genes <- m$genes[[1]]

    microglia_fraction <- marker_fraction(genes, marker_sets$microglia)
    neuropil_fraction <- marker_fraction(genes, marker_sets$neuropil_synaptic_neuronal)
    astrocyte_fraction <- marker_fraction(genes, marker_sets$astrocyte)
    oligodendrocyte_fraction <- marker_fraction(genes, marker_sets$oligodendrocyte_myelin)
    vascular_fraction <- marker_fraction(genes, marker_sets$endothelial_pericyte)

    same_direction <- !is.na(m$score[[1]]) && !is.na(best$neuropil_score[[1]]) &&
      sign(m$score[[1]]) == sign(best$neuropil_score[[1]])

    cls <- classify_term(
      microglia_marker_fraction = microglia_fraction,
      neuropil_marker_fraction = neuropil_fraction,
      neuropil_term_present = best$neuropil_term_present[[1]],
      gene_overlap_fraction = best$gene_overlap_fraction[[1]],
      same_direction = same_direction,
      microglia_padj = m$p_adjust[[1]],
      neuropil_padj = best$neuropil_p_adjust[[1]]
    )

    annotations[[i]] <- bind_cols(
      m %>% select(-genes),
      tibble(
        reference_dataset = REFERENCE_DATASET,
        microglia_marker_fraction = microglia_fraction,
        neuropil_marker_fraction = neuropil_fraction,
        astrocyte_marker_fraction = astrocyte_fraction,
        oligodendrocyte_marker_fraction = oligodendrocyte_fraction,
        vascular_marker_fraction = vascular_fraction,
        microglia_marker_hits = marker_hits(genes, marker_sets$microglia),
        neuropil_marker_hits = marker_hits(genes, marker_sets$neuropil_synaptic_neuronal),
        astrocyte_marker_hits = marker_hits(genes, marker_sets$astrocyte),
        oligodendrocyte_marker_hits = marker_hits(genes, marker_sets$oligodendrocyte_myelin),
        vascular_marker_hits = marker_hits(genes, marker_sets$endothelial_pericyte),
        same_direction_as_neuropil = same_direction,
        interpretation_class = cls
      ),
      best
    )
  }

  bind_rows(annotations)
}

write_outputs <- function(annotated, diagnostics) {
  table_dir <- CANONICAL_PATHS$tables
  figure_dir <- CANONICAL_PATHS$figures
  report_dir <- CANONICAL_PATHS$reports
  log_dir <- CANONICAL_PATHS$logs
  invisible(lapply(c(table_dir, figure_dir, report_dir, log_dir), dir_create))

  annotation_csv <- file.path(table_dir, paste0("microglia_neuropil_annotation_", RUN_ID, ".csv"))
  latest_csv <- file.path(table_dir, "microglia_neuropil_annotation_latest.csv")
  summary_csv <- file.path(table_dir, paste0("microglia_neuropil_annotation_summary_", RUN_ID, ".csv"))
  diagnostics_csv <- file.path(report_dir, paste0("neuropil_annotation_diagnostics_", RUN_ID, ".csv"))

  if (nrow(annotated) > 0) {
    readr::write_csv(annotated, annotation_csv)
    readr::write_csv(annotated, latest_csv)

    summary <- annotated %>%
      count(.data$interpretation_class, .data$result_type, .data$ontology, name = "n_terms") %>%
      arrange(.data$result_type, .data$ontology, desc(.data$n_terms))
    readr::write_csv(summary, summary_csv)

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- annotated %>%
        count(.data$interpretation_class, name = "n_terms") %>%
        ggplot(aes(x = reorder(.data$interpretation_class, .data$n_terms), y = .data$n_terms)) +
        geom_col() +
        coord_flip() +
        labs(
          x = "Interpretation class",
          y = "GO/GSEA terms",
          title = "Neuropil annotation of microglia-enriched terms"
        ) +
        theme_minimal(base_size = 10)
      ggsave(file.path(figure_dir, paste0("neuropil_annotation_class_counts_", RUN_ID, ".svg")), p, width = 7, height = 4)
    }
  } else {
    readr::write_csv(tibble(), annotation_csv)
    readr::write_csv(tibble(), latest_csv)
    readr::write_csv(tibble(), summary_csv)
  }

  readr::write_csv(diagnostics, diagnostics_csv)

  methods_note <- c(
    "Neuropil reference annotation",
    "",
    paste0("Dataset: ", DATASET),
    paste0("Reference dataset: ", REFERENCE_DATASET),
    paste0("Run ID: ", RUN_ID),
    "",
    "Interpretation:",
    "This workflow uses the neuropil dataset as an annotation/reference layer.",
    "It does not subtract logFC or protein intensities because separately normalized and imputed datasets are not on a guaranteed common quantitative scale.",
    "Terms are classified using neuropil term presence, core-gene overlap, direction concordance, adjusted P values, and marker composition.",
    "",
    "Main output:",
    annotation_csv
  )
  writeLines(methods_note, file.path(report_dir, paste0("neuropil_annotation_methods_note_", RUN_ID, ".txt")))

  invisible(list(annotation = annotation_csv, summary = summary_csv, diagnostics = diagnostics_csv))
}

microglia_manifest <- load_manifest(DATASET)
neuropil_manifest <- load_manifest(REFERENCE_DATASET)

if (isTRUE(DRY_RUN)) {
  diagnostics <- tibble(
    check = c(
      "dataset",
      "reference_dataset",
      "microglia_manifest_exists",
      "neuropil_manifest_exists",
      "tables_output_dir",
      "figures_output_dir",
      "reports_output_dir"
    ),
    status = c(
      if (DATASET == "microglia") "PASS" else "WARN",
      "PASS",
      if (file.exists(microglia_manifest$path)) "PASS" else "WARN",
      if (file.exists(neuropil_manifest$path)) "PASS" else "WARN",
      "PASS",
      "PASS",
      "PASS"
    ),
    detail = c(
      DATASET,
      REFERENCE_DATASET,
      microglia_manifest$path,
      neuropil_manifest$path,
      CANONICAL_PATHS$tables,
      CANONICAL_PATHS$figures,
      CANONICAL_PATHS$reports
    )
  )
  dry_run_line("Script", "04_neuropil_reference_annotation.r")
  dry_run_line("Dataset", DATASET, diagnostics$status[1])
  dry_run_line("Reference dataset", REFERENCE_DATASET)
  dry_run_line("Microglia/current manifest", microglia_manifest$path, diagnostics$status[3])
  dry_run_line("Neuropil reference manifest", neuropil_manifest$path, diagnostics$status[4])
  dry_run_line("Output tables", CANONICAL_PATHS$tables)
  write_outputs(tibble(), diagnostics)
  quit(status = 0, save = "no")
}

if (DATASET != "microglia" && tolower(Sys.getenv("PROTEOMICS_FORCE_NEUROPIL_ANNOTATION", unset = "false")) != "true") {
  diagnostics <- make_empty_annotation("skipped_non_microglia_dataset")
  write_outputs(tibble(), diagnostics)
  message("Skipping neuropil annotation because DATASET != microglia. Set PROTEOMICS_FORCE_NEUROPIL_ANNOTATION=true to force.")
  quit(status = 0, save = "no")
}

if (microglia_manifest$status != "loaded" || nrow(microglia_manifest$data) == 0) {
  diagnostics <- make_empty_annotation("missing_current_dataset_clusterProfiler_manifest")
  write_outputs(tibble(), diagnostics)
  warning("Current dataset clusterProfiler manifest missing or empty: ", microglia_manifest$path)
  quit(status = 0, save = "no")
}

if (neuropil_manifest$status != "loaded" || nrow(neuropil_manifest$data) == 0) {
  diagnostics <- make_empty_annotation("missing_neuropil_reference_clusterProfiler_manifest")
  write_outputs(tibble(), diagnostics)
  warning("Neuropil reference clusterProfiler manifest missing or empty: ", neuropil_manifest$path)
  quit(status = 0, save = "no")
}

microglia_tables <- extract_enrichment_tables(microglia_manifest$data, DATASET)
neuropil_tables <- extract_enrichment_tables(neuropil_manifest$data, REFERENCE_DATASET)

microglia_terms <- prepare_terms(microglia_tables)
neuropil_terms <- prepare_terms(neuropil_tables)

annotated <- annotate_microglia_terms(microglia_terms, neuropil_terms)

diagnostics <- tibble(
  check = c(
    "current_dataset",
    "reference_dataset",
    "current_manifest_rows",
    "reference_manifest_rows",
    "current_enrichment_rows",
    "reference_enrichment_rows",
    "annotated_terms"
  ),
  status = c("PASS", "PASS", "PASS", "PASS", if (nrow(microglia_terms) > 0) "PASS" else "WARN", if (nrow(neuropil_terms) > 0) "PASS" else "WARN", if (nrow(annotated) > 0) "PASS" else "WARN"),
  detail = c(
    DATASET,
    REFERENCE_DATASET,
    as.character(nrow(microglia_manifest$data)),
    as.character(nrow(neuropil_manifest$data)),
    as.character(nrow(microglia_terms)),
    as.character(nrow(neuropil_terms)),
    as.character(nrow(annotated))
  )
)

outputs <- write_outputs(annotated, diagnostics)
message("Neuropil annotation completed.")
message("Annotation table: ", outputs$annotation)
message("Summary table: ", outputs$summary)
message("Diagnostics: ", outputs$diagnostics)
