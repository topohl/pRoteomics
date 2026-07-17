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
  terms %>%
    mutate(
      source_dataset = .data$dataset,
      source_table = .data$output_table,
      source_manifest = bundle$manifest_source,
      source_term_provenance = .data$term_gene_provenance_file,
      term_key = paste(.data$result_type, .data$ontology, .data$term_id, sep = "::"),
      gene_string = .data$official_gene_symbol,
      n_genes = lengths(.data$genes)
    ) %>%
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
  if (nrow(microglia_terms) == 0) {
    return(tibble(
      dataset = character(), comparison = character(), result_type = character(), ontology = character(),
      term_id = character(), term_description = character(), p_adjust = numeric(), score = numeric(),
      official_gene_symbol = character(), ProteinGroupID = character(), member_accessions = character(),
      protein_group_gene_annotation_status = character(), gene_level_claim_allowed = logical(),
      route_category = character(), route_unit = character(), output_table = character(),
      term_gene_provenance_file = character(), analysis_status = character(),
      source_dataset = character(), source_table = character(), source_manifest = character(),
      source_term_provenance = character(), term_key = character(), gene_string = character(),
      n_genes = integer(), reference_dataset = character(), microglia_marker_fraction = numeric(),
      neuropil_marker_fraction = numeric(), astrocyte_marker_fraction = numeric(),
      oligodendrocyte_marker_fraction = numeric(), vascular_marker_fraction = numeric(),
      microglia_marker_hits = character(), neuropil_marker_hits = character(),
      astrocyte_marker_hits = character(), oligodendrocyte_marker_hits = character(),
      vascular_marker_hits = character(), same_direction_as_neuropil = logical(),
      interpretation_class = character(), neuropil_term_present = logical(),
      neuropil_comparison = character(), neuropil_p_adjust = numeric(), neuropil_score = numeric(),
      neuropil_source_table = character(), gene_overlap_fraction = numeric(),
      gene_jaccard = numeric(), overlapping_genes = character()
    ))
  }

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
  write_outputs(annotate_microglia_terms(tibble(), tibble()), diagnostics, tibble())
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
