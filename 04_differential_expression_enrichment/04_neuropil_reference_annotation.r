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
source(repo_path("R", "enrichment_io.R"))

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

REFERENCE_MATCHING_CONTRACT_VERSION <- "microglia_neuropil_reference_v2_region_contrast_fdr_aggregate"
CLASSIFICATION_CONTRACT_VERSION <- "microglia_neuropil_interpretation_v2_microglia_and_reference_fdr"

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
    "Aif1", "P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Tyrobp", "Hexb",
    "C1qa", "C1qb", "C1qc", "Itgam", "Spi1", "Trem2", "Laptm5"
  ),
  neuropil_synaptic_neuronal = c(
    "Syn1", "Syp", "Snap25", "Stx1a", "Stxbp1", "Dlg4", "Camk2a", "Camk2b",
    "Map2", "Nefl", "Nefm", "Rbfox3", "Tubb3", "Grin1", "Gria1", "Vamp2"
  ),
  astrocyte = c("Gfap", "Aldh1l1", "Aqp4", "Slc1a3", "Slc1a2", "Glul", "Gja1"),
  oligodendrocyte_myelin = c("Mbp", "Plp1", "Mag", "Mog", "Mobp", "Cnp", "Cldn11"),
  endothelial_pericyte = c("Pecam1", "Cldn5", "Kdr", "Flt1", "Pdgfrb", "Rgs5", "Acta2")
)

normalize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[nzchar(x)]
}

parse_reference_identity <- function(comparison, route_category, route_unit) {
  n <- max(length(comparison), length(route_category), length(route_unit))
  lengths_ok <- vapply(
    list(comparison, route_category, route_unit),
    function(x) length(x) %in% c(1L, n),
    logical(1)
  )
  if (!all(lengths_ok)) {
    stop("comparison, route_category, and route_unit must have compatible lengths.", call. = FALSE)
  }

  comparison <- rep(as.character(comparison), length.out = n)
  route_category <- rep(as.character(route_category), length.out = n)
  route_unit <- rep(as.character(route_unit), length.out = n)

  region_match <- regexec("^(CA1|CA2|CA3|DG)(?:_|$)", route_unit, perl = TRUE)
  region_parts <- regmatches(route_unit, region_match)
  anatomical_region <- vapply(
    region_parts,
    function(x) if (length(x) >= 2L) x[[2L]] else NA_character_,
    character(1)
  )

  parse_contrast <- function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    matched <- regmatches(tolower(x), regexec("^(.*)(con|res|sus)_(.*)(con|res|sus)$", tolower(x), perl = TRUE))[[1]]
    if (length(matched) != 5L) return(NA_character_)
    paste0(toupper(matched[[3L]]), "-vs-", toupper(matched[[5L]]))
  }
  formal_group_contrast <- unname(vapply(comparison, parse_contrast, character(1)))
  identity_valid <- route_category == "phenotype_within_unit" &
    !is.na(anatomical_region) & !is.na(formal_group_contrast)

  tibble::tibble(
    anatomical_region = anatomical_region,
    formal_group_contrast = formal_group_contrast,
    reference_identity_status = ifelse(identity_valid, "valid", "invalid"),
    reference_identity_key = ifelse(
      identity_valid,
      paste(anatomical_region, formal_group_contrast, sep = "::"),
      NA_character_
    )
  )
}

add_reference_identity <- function(terms, source_label = "canonical enrichment") {
  if (!nrow(terms)) {
    terms$anatomical_region <- character()
    terms$formal_group_contrast <- character()
    terms$reference_identity_status <- character()
    terms$reference_identity_key <- character()
    terms$neuropil_match_key <- character()
    return(terms)
  }

  identity <- parse_reference_identity(terms$comparison, terms$route_category, terms$route_unit)
  invalid <- which(identity$reference_identity_status != "valid")
  if (length(invalid)) {
    examples <- unique(paste(
      terms$comparison[invalid], terms$route_category[invalid], terms$route_unit[invalid],
      sep = " | "
    ))
    stop(
      "Unable to derive a valid region/formal-contrast identity for ", source_label,
      ": ", paste(utils::head(examples, 5L), collapse = "; "),
      call. = FALSE
    )
  }

  dplyr::bind_cols(terms, identity) %>%
    dplyr::mutate(
      neuropil_match_key = paste(.data$term_key, .data$reference_identity_key, sep = "::")
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
  if (!file.exists(path)) return(list(path = path, data = data.frame(), bundle = NULL, status = "missing"))
  manifest <- read_canonical_clusterprofiler_manifest(
    path, dataset, strict = TRUE, require_files = !isTRUE(DRY_RUN)
  )
  bundle <- if (isTRUE(DRY_RUN)) NULL else read_canonical_clusterprofiler_bundle(
    path, dataset, result_types = "GSEA_GO", strict = TRUE
  )
  list(path = path, data = manifest, bundle = bundle, status = "loaded")
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

prepare_terms <- function(bundle) {
  manifest_identity <- bundle$manifest %>%
    select(dataset, comparison, result_type, ontology, route_category, route_unit,
      output_table, term_gene_provenance_file, analysis_status)
  provenance <- bundle$provenance %>%
    filter(.data$gene_level_claim_allowed, .data$core_enrichment_member) %>%
    group_by(.data$dataset, .data$comparison, .data$result_type, .data$ontology,
      .data$term_id, .data$term_description) %>%
    summarise(
      genes = list(sort(unique(as.character(.data$official_gene_symbol)), method = "radix")),
      official_gene_symbol = paste(sort(unique(as.character(.data$official_gene_symbol)), method = "radix"), collapse = ";"),
      ProteinGroupID = paste(sort(unique(as.character(.data$ProteinGroupID)), method = "radix"), collapse = ";"),
      member_accessions = paste(sort(unique(as.character(.data$member_accessions)), method = "radix"), collapse = ";"),
      protein_group_gene_annotation_status = paste(sort(unique(as.character(.data$protein_group_gene_annotation_status)), method = "radix"), collapse = ";"),
      gene_level_claim_allowed = all(.data$gene_level_claim_allowed),
      .groups = "drop"
    )
  terms <- bundle$terms %>%
    transmute(
      dataset, comparison, result_type, ontology, term_id, term_description,
      p_adjust = suppressWarnings(as.numeric(.data$`p.adjust`)),
      score = suppressWarnings(as.numeric(.data$NES))
    ) %>%
    left_join(provenance, by = c("dataset", "comparison", "result_type", "ontology", "term_id", "term_description")) %>%
    left_join(manifest_identity, by = c("dataset", "comparison", "result_type", "ontology"))
  if (nrow(terms) && any(vapply(terms$genes, is.null, logical(1)))) {
    stop("Canonical enrichment term is missing required term-gene provenance.", call. = FALSE)
  }
  terms <- terms %>%
    mutate(
      source_dataset = .data$dataset,
      source_table = .data$output_table,
      source_manifest = bundle$manifest_source,
      source_term_provenance = .data$term_gene_provenance_file,
      term_key = paste(.data$result_type, .data$ontology, .data$term_id, sep = "::"),
      gene_string = .data$official_gene_symbol,
      n_genes = lengths(.data$genes)
    )
  add_reference_identity(terms, paste0("dataset ", unique(terms$dataset))) %>%
    arrange(.data$dataset, .data$comparison, .data$result_type, .data$ontology, .data$term_id)
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
                          microglia_padj,
                          any_significant_neuropil_match,
                          max_matched_gene_overlap_fraction,
                          max_significant_gene_overlap_fraction,
                          max_significant_same_direction_gene_overlap_fraction) {
  microglia_sig <- !is.na(microglia_padj) && microglia_padj < 0.05

  if (!isTRUE(microglia_sig)) {
    return("not_evaluated_non_significant")
  }

  if (!isTRUE(any_significant_neuropil_match) &&
      !is.na(microglia_marker_fraction) && microglia_marker_fraction >= 0.15 &&
      (is.na(neuropil_marker_fraction) || neuropil_marker_fraction < 0.15) &&
      !is.na(max_matched_gene_overlap_fraction) && max_matched_gene_overlap_fraction < 0.20) {
    return("microglia_robust")
  }

  if (isTRUE(any_significant_neuropil_match) &&
      !is.na(neuropil_marker_fraction) && neuropil_marker_fraction >= 0.15 &&
      !is.na(max_significant_same_direction_gene_overlap_fraction) &&
      max_significant_same_direction_gene_overlap_fraction >= 0.30) {
    return("neuropil_sensitive")
  }

  if (isTRUE(any_significant_neuropil_match) &&
      !is.na(max_significant_gene_overlap_fraction) &&
      max_significant_gene_overlap_fraction >= 0.20) {
    return("mixed_microenvironment")
  }

  if (!is.na(neuropil_marker_fraction) && neuropil_marker_fraction >= 0.25) {
    return("neuropil_marker_enriched")
  }

  "ambiguous"
}

empty_neuropil_evidence <- function() {
  tibble::tibble(
    neuropil_term_present = FALSE,
    neuropil_comparison = NA_character_,
    neuropil_p_adjust = NA_real_,
    neuropil_score = NA_real_,
    neuropil_source_table = NA_character_,
    gene_overlap_fraction = 0,
    gene_jaccard = 0,
    overlapping_genes = "",
    same_direction_as_neuropil = FALSE,
    matched_neuropil_comparisons = "",
    matched_neuropil_layers = "",
    significant_neuropil_comparisons = "",
    n_matched_neuropil_comparisons = 0L,
    n_matched_neuropil_layers = 0L,
    n_significant_neuropil_matches = 0L,
    any_significant_neuropil_match = FALSE,
    min_matched_neuropil_p_adjust = NA_real_,
    any_significant_neuropil_same_direction = FALSE,
    max_matched_gene_overlap_fraction = 0,
    max_matched_gene_jaccard = 0,
    max_significant_gene_overlap_fraction = NA_real_,
    max_significant_gene_jaccard = NA_real_,
    max_significant_same_direction_gene_overlap_fraction = NA_real_,
    max_significant_same_direction_gene_jaccard = NA_real_
  )
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

safe_max <- function(x, empty = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) empty else max(x)
}

find_best_neuropil_match <- function(m_row, neuropil_terms, candidate_rows = NULL) {
  if (nrow(neuropil_terms) == 0) return(empty_neuropil_evidence())

  if (!is.null(candidate_rows)) {
    candidate_rows <- as.integer(candidate_rows)
    candidate_rows <- candidate_rows[
      !is.na(candidate_rows) & candidate_rows >= 1L & candidate_rows <= nrow(neuropil_terms)
    ]
    neuropil_terms <- neuropil_terms[candidate_rows, , drop = FALSE]
  }
  if (nrow(neuropil_terms) == 0) return(empty_neuropil_evidence())

  same_term <- neuropil_terms %>%
    dplyr::filter(
      .data$term_key == m_row$term_key[[1]],
      .data$anatomical_region == m_row$anatomical_region[[1]],
      .data$formal_group_contrast == m_row$formal_group_contrast[[1]]
    )
  if (nrow(same_term) == 0) return(empty_neuropil_evidence())

  genes_m <- m_row$genes[[1]]
  microglia_score <- suppressWarnings(as.numeric(m_row$score[[1]]))
  scored <- same_term %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      gene_overlap_fraction = term_overlap_fraction(genes_m, genes),
      gene_jaccard = jaccard(genes_m, genes),
      overlapping_genes = paste(sort(intersect(normalize_id(genes_m), normalize_id(genes))), collapse = ";"),
      neuropil_significant = !is.na(.data$p_adjust) & .data$p_adjust < 0.05,
      same_direction = is.finite(microglia_score) & is.finite(.data$score) &
        sign(microglia_score) == sign(.data$score)
    ) %>%
    dplyr::ungroup()

  representative <- scored %>%
    dplyr::arrange(
      dplyr::desc(.data$neuropil_significant),
      dplyr::desc(.data$gene_overlap_fraction),
      dplyr::desc(.data$gene_jaccard),
      .data$p_adjust,
      .data$comparison,
      .data$route_unit,
      .data$source_table
    ) %>%
    dplyr::slice(1L)
  significant <- scored %>% dplyr::filter(.data$neuropil_significant)
  significant_same_direction <- significant %>% dplyr::filter(.data$same_direction)

  tibble::tibble(
    neuropil_term_present = TRUE,
    neuropil_comparison = representative$comparison[[1]],
    neuropil_p_adjust = representative$p_adjust[[1]],
    neuropil_score = representative$score[[1]],
    neuropil_source_table = representative$source_table[[1]],
    gene_overlap_fraction = representative$gene_overlap_fraction[[1]],
    gene_jaccard = representative$gene_jaccard[[1]],
    overlapping_genes = representative$overlapping_genes[[1]],
    same_direction_as_neuropil = representative$same_direction[[1]],
    matched_neuropil_comparisons = paste(sort(unique(scored$comparison), method = "radix"), collapse = ";"),
    matched_neuropil_layers = paste(sort(unique(scored$route_unit), method = "radix"), collapse = ";"),
    significant_neuropil_comparisons = paste(sort(unique(significant$comparison), method = "radix"), collapse = ";"),
    n_matched_neuropil_comparisons = dplyr::n_distinct(scored$comparison),
    n_matched_neuropil_layers = dplyr::n_distinct(scored$route_unit),
    n_significant_neuropil_matches = nrow(significant),
    any_significant_neuropil_match = nrow(significant) > 0L,
    min_matched_neuropil_p_adjust = safe_min(scored$p_adjust),
    any_significant_neuropil_same_direction = any(significant$same_direction),
    max_matched_gene_overlap_fraction = safe_max(scored$gene_overlap_fraction, empty = 0),
    max_matched_gene_jaccard = safe_max(scored$gene_jaccard, empty = 0),
    max_significant_gene_overlap_fraction = safe_max(significant$gene_overlap_fraction),
    max_significant_gene_jaccard = safe_max(significant$gene_jaccard),
    max_significant_same_direction_gene_overlap_fraction = safe_max(significant_same_direction$gene_overlap_fraction),
    max_significant_same_direction_gene_jaccard = safe_max(significant_same_direction$gene_jaccard)
  )
}

annotate_microglia_terms <- function(microglia_terms, neuropil_terms) {
  if (nrow(microglia_terms) == 0) {
    return(tibble::tibble(
      dataset = character(), comparison = character(), result_type = character(), ontology = character(),
      term_id = character(), term_description = character(), p_adjust = numeric(), score = numeric(),
      official_gene_symbol = character(), ProteinGroupID = character(), member_accessions = character(),
      protein_group_gene_annotation_status = character(), gene_level_claim_allowed = logical(),
      route_category = character(), route_unit = character(), output_table = character(),
      term_gene_provenance_file = character(), analysis_status = character(),
      source_dataset = character(), source_table = character(), source_manifest = character(),
      source_term_provenance = character(), term_key = character(), gene_string = character(),
      n_genes = integer(), anatomical_region = character(), formal_group_contrast = character(),
      reference_identity_status = character(), reference_identity_key = character(),
      neuropil_match_key = character(), reference_dataset = character(), reference_region = character(),
      reference_contrast = character(), reference_matching_contract_version = character(),
      classification_contract_version = character(), microglia_significant = logical(),
      microglia_marker_fraction = numeric(),
      neuropil_marker_fraction = numeric(), astrocyte_marker_fraction = numeric(),
      oligodendrocyte_marker_fraction = numeric(), vascular_marker_fraction = numeric(),
      microglia_marker_hits = character(), neuropil_marker_hits = character(),
      astrocyte_marker_hits = character(), oligodendrocyte_marker_hits = character(),
      vascular_marker_hits = character(), same_direction_as_neuropil = logical(),
      interpretation_class = character(), neuropil_term_present = logical(),
      neuropil_comparison = character(), neuropil_p_adjust = numeric(), neuropil_score = numeric(),
      neuropil_source_table = character(), gene_overlap_fraction = numeric(),
      gene_jaccard = numeric(), overlapping_genes = character(),
      matched_neuropil_comparisons = character(), matched_neuropil_layers = character(),
      significant_neuropil_comparisons = character(), n_matched_neuropil_comparisons = integer(),
      n_matched_neuropil_layers = integer(), n_significant_neuropil_matches = integer(),
      any_significant_neuropil_match = logical(), min_matched_neuropil_p_adjust = numeric(),
      any_significant_neuropil_same_direction = logical(),
      max_matched_gene_overlap_fraction = numeric(), max_matched_gene_jaccard = numeric(),
      max_significant_gene_overlap_fraction = numeric(), max_significant_gene_jaccard = numeric(),
      max_significant_same_direction_gene_overlap_fraction = numeric(),
      max_significant_same_direction_gene_jaccard = numeric()
    ))
  }

  neuropil_index <- if (nrow(neuropil_terms)) {
    split(seq_len(nrow(neuropil_terms)), neuropil_terms$neuropil_match_key, drop = TRUE)
  } else {
    list()
  }
  annotations <- vector("list", nrow(microglia_terms))
  for (i in seq_len(nrow(microglia_terms))) {
    m <- microglia_terms[i, , drop = FALSE]
    candidate_rows <- neuropil_index[[m$neuropil_match_key[[1]]]]
    if (is.null(candidate_rows)) candidate_rows <- integer()
    reference_evidence <- find_best_neuropil_match(m, neuropil_terms, candidate_rows)
    genes <- m$genes[[1]]

    microglia_fraction <- marker_fraction(genes, marker_sets$microglia)
    neuropil_fraction <- marker_fraction(genes, marker_sets$neuropil_synaptic_neuronal)
    astrocyte_fraction <- marker_fraction(genes, marker_sets$astrocyte)
    oligodendrocyte_fraction <- marker_fraction(genes, marker_sets$oligodendrocyte_myelin)
    vascular_fraction <- marker_fraction(genes, marker_sets$endothelial_pericyte)

    cls <- classify_term(
      microglia_marker_fraction = microglia_fraction,
      neuropil_marker_fraction = neuropil_fraction,
      microglia_padj = m$p_adjust[[1]],
      any_significant_neuropil_match = reference_evidence$any_significant_neuropil_match[[1]],
      max_matched_gene_overlap_fraction = reference_evidence$max_matched_gene_overlap_fraction[[1]],
      max_significant_gene_overlap_fraction = reference_evidence$max_significant_gene_overlap_fraction[[1]],
      max_significant_same_direction_gene_overlap_fraction = reference_evidence$max_significant_same_direction_gene_overlap_fraction[[1]]
    )

    annotations[[i]] <- dplyr::bind_cols(
      m %>% dplyr::select(-genes),
      tibble::tibble(
        reference_dataset = REFERENCE_DATASET,
        reference_region = m$anatomical_region[[1]],
        reference_contrast = m$formal_group_contrast[[1]],
        reference_matching_contract_version = REFERENCE_MATCHING_CONTRACT_VERSION,
        classification_contract_version = CLASSIFICATION_CONTRACT_VERSION,
        microglia_significant = !is.na(m$p_adjust[[1]]) && m$p_adjust[[1]] < 0.05,
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
        interpretation_class = cls
      ),
      reference_evidence
    )
  }

  dplyr::bind_rows(annotations)
}

write_outputs <- function(annotated, diagnostics, analysis_status = tibble()) {
  table_dir <- CANONICAL_PATHS$tables
  figure_dir <- CANONICAL_PATHS$figures
  report_dir <- CANONICAL_PATHS$reports
  log_dir <- CANONICAL_PATHS$logs
  invisible(lapply(c(table_dir, figure_dir, report_dir, log_dir), dir_create))

  annotation_csv <- file.path(table_dir, paste0("microglia_neuropil_annotation_", RUN_ID, ".csv"))
  latest_csv <- file.path(table_dir, "microglia_neuropil_annotation_latest.csv")
  summary_csv <- file.path(table_dir, paste0("microglia_neuropil_annotation_summary_", RUN_ID, ".csv"))
  status_csv <- file.path(table_dir, "microglia_neuropil_analysis_status.csv")
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
    readr::write_csv(annotated, annotation_csv)
    readr::write_csv(annotated, latest_csv)
    readr::write_csv(tibble(interpretation_class = character(), result_type = character(), ontology = character(), n_terms = integer()), summary_csv)
  }

  readr::write_csv(analysis_status, status_csv)
  readr::write_csv(diagnostics, diagnostics_csv)

  methods_note <- c(
    "Neuropil reference annotation",
    "",
    paste0("Dataset: ", DATASET),
    paste0("Reference dataset: ", REFERENCE_DATASET),
    paste0("Run ID: ", RUN_ID),
    paste0("Reference matching contract: ", REFERENCE_MATCHING_CONTRACT_VERSION),
    paste0("Classification contract: ", CLASSIFICATION_CONTRACT_VERSION),
    "",
    "Interpretation:",
    "This workflow uses the neuropil dataset as an annotation/reference layer.",
    "It does not subtract logFC or protein intensities because separately normalized and imputed datasets are not on a guaranteed common quantitative scale.",
    "Neuropil evidence is restricted to identical GO terms from the same anatomical region and formal condition contrast; all corresponding neuropil layers are retained as aggregate evidence.",
    "Term presence is descriptive only. Biological interpretation classes require microglia FDR < 0.05, and neuropil-supported classes require matched neuropil FDR < 0.05.",
    "The retained representative neuropil comparison is selected deterministically by significance, overlap, Jaccard, adjusted P value, and lexical tie-breaks; classification uses aggregate matched-reference evidence.",
    "",
    "Main output:",
    annotation_csv
  )
  writeLines(methods_note, file.path(report_dir, paste0("neuropil_annotation_methods_note_", RUN_ID, ".txt")))

  invisible(list(annotation = annotation_csv, summary = summary_csv, status = status_csv, diagnostics = diagnostics_csv))
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
  quit(status = 0, save = "no")
}

if (DATASET != "microglia" && tolower(Sys.getenv("PROTEOMICS_FORCE_NEUROPIL_ANNOTATION", unset = "false")) != "true") {
  diagnostics <- make_empty_annotation("skipped_non_microglia_dataset")
  write_outputs(annotate_microglia_terms(tibble(), tibble()), diagnostics, tibble())
  message("Skipping neuropil annotation because DATASET != microglia. Set PROTEOMICS_FORCE_NEUROPIL_ANNOTATION=true to force.")
  quit(status = 0, save = "no")
}

if (microglia_manifest$status != "loaded" || nrow(microglia_manifest$data) == 0) {
  diagnostics <- make_empty_annotation("missing_current_dataset_clusterProfiler_manifest")
  write_outputs(annotate_microglia_terms(tibble(), tibble()), diagnostics, tibble())
  warning("Current dataset clusterProfiler manifest missing or empty: ", microglia_manifest$path)
  quit(status = 0, save = "no")
}

if (neuropil_manifest$status != "loaded" || nrow(neuropil_manifest$data) == 0) {
  diagnostics <- make_empty_annotation("missing_neuropil_reference_clusterProfiler_manifest")
  write_outputs(annotate_microglia_terms(tibble(), tibble()), diagnostics, tibble())
  warning("Neuropil reference clusterProfiler manifest missing or empty: ", neuropil_manifest$path)
  quit(status = 0, save = "no")
}

microglia_terms <- prepare_terms(microglia_manifest$bundle)
neuropil_terms <- prepare_terms(neuropil_manifest$bundle)
analysis_status <- dplyr::bind_rows(
  microglia_manifest$bundle$status %>% dplyr::mutate(source_role = "target", source_manifest = microglia_manifest$path),
  neuropil_manifest$bundle$status %>% dplyr::mutate(source_role = "neuropil_reference", source_manifest = neuropil_manifest$path)
) %>%
  dplyr::arrange(.data$source_role, .data$dataset, .data$comparison, .data$result_type, .data$ontology)

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

outputs <- write_outputs(annotated, diagnostics, analysis_status)
message("Neuropil annotation completed.")
message("Annotation table: ", outputs$annotation)
message("Summary table: ", outputs$summary)
message("Diagnostics: ", outputs$diagnostics)
