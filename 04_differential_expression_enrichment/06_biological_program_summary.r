#!/usr/bin/env Rscript

# ================================================================
# Script: 04_differential_expression_enrichment/06_biological_program_summary.r
# Stage: enrichment
# Scope: dataset_specific
# Consumes: canonical compareGO manifest and declared term/provenance/status tables; optional canonical neuropil-reference and targeted-signature annotations.
# Produces: results/tables/04_differential_expression_enrichment/biological_program_summary/<dataset>/program_summary.csv.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Runs after clusterProfiler and compareGO; uses annotations/signatures where supported.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "enrichment_plots.R"))
source(repo_path("R", "plotting_nature.R"))

SCRIPT_ID <- "04_differential_expression_enrichment/06_biological_program_summary.r"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
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

manifest_file <- canonical_comparego_manifest_path(DATASET)

if (is_dry_run()) {
  if (file.exists(manifest_file)) {
    read_canonical_comparego_manifest(manifest_file, DATASET, require_files = FALSE)
  }
  dry_run_line("Script", "04_differential_expression_enrichment/06_biological_program_summary.r")
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
  stop("Canonical compareGO manifest not found for dataset ", DATASET, ": ", manifest_file, call. = FALSE)
}

bundle <- read_canonical_comparego_bundle(manifest_file, DATASET)
manifest <- bundle$manifest
path_qc <- validate_manifest_paths(
  manifest,
  path_cols = c("term_comparison_file", "term_gene_provenance_output_file", "analysis_status_summary_file"),
  allow_missing = FALSE
)
readr::write_csv(path_qc, file.path(PATHS$logs, "program_summary_manifest_path_qc.csv"), na = "")
route_identity <- manifest %>%
  dplyr::select(.data$dataset, .data$comparison, .data$result_type, .data$ontology,
    .data$route_category, .data$route_unit) %>%
  dplyr::distinct()
terms <- bundle$terms %>%
  dplyr::left_join(route_identity, by = c("dataset", "comparison", "result_type", "ontology"))

if (!"pvalue" %in% names(terms)) terms$pvalue <- rep(NA_real_, nrow(terms))
if (!"p.value" %in% names(terms)) terms$`p.value` <- rep(NA_real_, nrow(terms))

terms <- map_terms_to_programs(terms, "term_description") %>%
  dplyr::filter(!is.na(.data$biological_program)) %>%
  dplyr::mutate(
    ID = as.character(.data$term_id),
    Description = as.character(.data$term_description),
    NES = if ("NES" %in% names(.)) suppressWarnings(as.numeric(.data$NES)) else NA_real_,
    p.adjust = if ("p.adjust" %in% names(.)) suppressWarnings(as.numeric(.data$p.adjust)) else NA_real_,
    FDR = .data$p.adjust,
    source_file = bundle$term_source,
    term_provenance_source = bundle$provenance_source,
    comparego_manifest_source = bundle$manifest_source,
    raw_p = dplyr::coalesce(
      suppressWarnings(as.numeric(.data$pvalue)),
      suppressWarnings(as.numeric(.data$p.value))
    )
  )

frequent_genes <- function(x, n = 12L) {
  genes <- stats::na.omit(as.character(x))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) return(NA_character_)
  tab <- as.data.frame(table(genes), stringsAsFactors = FALSE)
  tab <- tab[order(-tab$Freq, tab$genes, method = "radix"), , drop = FALSE]
  paste(tab$genes[seq_len(min(n, nrow(tab)))], collapse = ";")
}

strongest_term <- function(description, nes, direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  nes <- suppressWarnings(as.numeric(nes))
  ok <- is.finite(nes)
  if (!any(ok)) return(NA_character_)
  target <- if (direction == "positive") max(nes[ok], na.rm = TRUE) else min(nes[ok], na.rm = TRUE)
  hit <- which(ok & nes == target)
  as.character(description[hit[[1]]])
}

term_evidence <- terms %>%
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
    FDR, source_file, term_provenance_source, comparego_manifest_source
  )

evidence <- bundle$provenance %>%
  dplyr::filter(.data$gene_level_claim_allowed, .data$core_enrichment_member) %>%
  dplyr::inner_join(
    terms %>%
      dplyr::select(.data$dataset, .data$comparison, .data$result_type, .data$ontology,
        .data$term_id, .data$term_description, .data$route_category, .data$route_unit,
        .data$biological_program, .data$NES, raw_p = .data$raw_p,
        p.adjust = .data$p.adjust, FDR = .data$FDR),
    by = c("dataset", "comparison", "result_type", "ontology", "term_id", "term_description")
  ) %>%
  dplyr::mutate(
    effect_direction = dplyr::case_when(
      is.na(.data$NES) ~ "undirected",
      .data$NES > 0 ~ "positive_NES",
      .data$NES < 0 ~ "negative_NES",
      TRUE ~ "neutral"
    ),
    term_provenance_source = bundle$provenance_source,
    comparego_manifest_source = bundle$manifest_source
  ) %>%
  dplyr::select(
    dataset, comparison, route_category, route_unit, result_type, ontology,
    biological_program, term_id, term_description, official_gene_symbol,
    official_entrez_id, ProteinGroupID, member_accessions,
    protein_group_gene_annotation_status, gene_level_claim_allowed,
    rank_statistic, core_enrichment_member, NES, effect_direction, raw_p,
    p.adjust, FDR, term_provenance_source, comparego_manifest_source
  ) %>%
  dplyr::arrange(.data$dataset, .data$comparison, .data$ontology, .data$term_id,
    .data$official_gene_symbol, .data$ProteinGroupID)

key_gene_summary <- evidence %>%
  dplyr::distinct(.data$dataset, .data$comparison, .data$route_category, .data$route_unit,
    .data$biological_program, .data$term_id, .data$official_gene_symbol) %>%
  dplyr::group_by(.data$dataset, .data$comparison, .data$route_category, .data$route_unit,
    .data$biological_program) %>%
  dplyr::summarise(key_genes = frequent_genes(.data$official_gene_symbol), .groups = "drop")

program_summary <- term_evidence %>%
  dplyr::group_by(.data$dataset, .data$comparison, .data$route_category, .data$route_unit, .data$biological_program) %>%
  dplyr::arrange(.data$FDR, dplyr::desc(abs(.data$NES)), .by_group = TRUE) %>%
  dplyr::summarise(
    n_terms = dplyr::n(),
    min_fdr = suppressWarnings(min(.data$FDR, na.rm = TRUE)),
    min_raw_p = suppressWarnings(min(.data$raw_p, na.rm = TRUE)),
    representative_NES = dplyr::first(.data$NES),
    effect_direction = dplyr::first(.data$effect_direction),
    n_positive_terms = sum(.data$NES > 0, na.rm = TRUE),
    n_negative_terms = sum(.data$NES < 0, na.rm = TRUE),
    median_NES = suppressWarnings(stats::median(.data$NES, na.rm = TRUE)),
    strongest_positive_term = strongest_term(.data$Description, .data$NES, "positive"),
    strongest_negative_term = strongest_term(.data$Description, .data$NES, "negative"),
    top_term = dplyr::first(.data$Description),
    source_file = dplyr::first(.data$source_file),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    min_fdr = ifelse(is.infinite(.data$min_fdr), NA_real_, .data$min_fdr),
    min_raw_p = ifelse(is.infinite(.data$min_raw_p), NA_real_, .data$min_raw_p),
    median_NES = ifelse(is.nan(.data$median_NES), NA_real_, .data$median_NES),
    direction_consistency = dplyr::case_when(
      .data$n_positive_terms > 0 & .data$n_negative_terms > 0 ~ "mixed_direction",
      .data$n_positive_terms > 0 ~ "positive_only",
      .data$n_negative_terms > 0 ~ "negative_only",
      TRUE ~ "undirected"
    ),
    direction = dplyr::case_when(
      !is.na(.data$effect_direction) ~ .data$effect_direction,
      is.na(.data$representative_NES) ~ "undirected",
      .data$representative_NES > 0 ~ "positive_NES",
      .data$representative_NES < 0 ~ "negative_NES",
      TRUE ~ "neutral"
    )
  ) %>%
  dplyr::left_join(
    key_gene_summary,
    by = c("dataset", "comparison", "route_category", "route_unit", "biological_program")
  ) %>%
  dplyr::arrange(.data$min_fdr, .data$dataset, .data$comparison, .data$biological_program)

readr::write_csv(program_summary, file.path(PATHS$tables, "program_summary.csv"), na = "")
readr::write_csv(evidence, file.path(PATHS$source_data, "program_term_gene_evidence.csv"), na = "")
readr::write_csv(bundle$status, file.path(PATHS$tables, "program_analysis_status.csv"), na = "")

heatmap_ready <- program_summary %>%
  dplyr::mutate(
    comparison_label = paste(.data$comparison, .data$route_category, .data$route_unit, sep = " | "),
    signed_neg_log10_fdr = sign(dplyr::coalesce(.data$representative_NES, 0)) *
      -log10(pmax(dplyr::coalesce(.data$min_fdr, 1), .Machine$double.xmin))
  ) %>%
  dplyr::select(
    dataset, comparison, route_category, route_unit, comparison_label,
    biological_program, direction, representative_NES, min_raw_p, min_fdr,
    signed_neg_log10_fdr, n_terms, n_positive_terms, n_negative_terms,
    median_NES, direction_consistency, strongest_positive_term,
    strongest_negative_term, top_term, key_genes, source_file
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
  counts <- table(x)
  sort(names(counts)[counts == max(counts)], method = "radix")[[1]]
}

NEUROPIL_PROGRAM_ANNOTATION_CONTRACT_VERSION <- "microglia_neuropil_program_v2_significant_terms_by_program"

empty_neuropil_program_annotation <- function() {
  tibble::tibble(
    dataset = character(),
    comparison = character(),
    route_category = character(),
    route_unit = character(),
    biological_program = character(),
    interpretation_class = character(),
    dominant_interpretation_class = character(),
    neuropil_program_annotation_status = character(),
    neuropil_program_annotation_contract_version = character(),
    source_classification_contract_version = character(),
    n_significant_microglia_terms = integer(),
    n_mixed_microenvironment = integer(),
    n_neuropil_sensitive = integer(),
    n_microglia_robust = integer(),
    n_neuropil_marker_enriched = integer(),
    n_ambiguous = integer(),
    fraction_mixed_microenvironment = numeric(),
    fraction_neuropil_sensitive = numeric(),
    fraction_microglia_robust = numeric(),
    fraction_neuropil_marker_enriched = numeric(),
    fraction_ambiguous = numeric(),
    gene_overlap_fraction = numeric(),
    gene_jaccard = numeric(),
    neuropil_marker_fraction = numeric(),
    microglia_marker_fraction = numeric()
  )
}

safe_mean <- function(x) {
  value <- suppressWarnings(mean(as.numeric(x), na.rm = TRUE))
  if (is.nan(value)) NA_real_ else value
}

summarise_neuropil_annotation_by_program <- function(neuropil_annotation) {
  required <- c(
    "dataset", "comparison", "route_category", "route_unit", "term_description",
    "p_adjust", "interpretation_class"
  )
  missing <- setdiff(required, names(neuropil_annotation))
  if (length(missing)) {
    stop(
      "Neuropil annotation is missing required program-integration columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!nrow(neuropil_annotation)) return(empty_neuropil_program_annotation())

  optional_numeric <- c(
    "gene_overlap_fraction", "gene_jaccard", "neuropil_marker_fraction",
    "microglia_marker_fraction"
  )
  for (column in setdiff(optional_numeric, names(neuropil_annotation))) {
    neuropil_annotation[[column]] <- NA_real_
  }
  if (!"classification_contract_version" %in% names(neuropil_annotation)) {
    neuropil_annotation$classification_contract_version <- NA_character_
  }

  significant <- map_terms_to_programs(neuropil_annotation, "term_description") %>%
    dplyr::filter(
      is.finite(suppressWarnings(as.numeric(.data$p_adjust))),
      suppressWarnings(as.numeric(.data$p_adjust)) < 0.05,
      !is.na(.data$biological_program)
    )
  if (!nrow(significant)) return(empty_neuropil_program_annotation())

  significant %>%
    dplyr::group_by(
      .data$dataset, .data$comparison, .data$route_category, .data$route_unit,
      .data$biological_program
    ) %>%
    dplyr::summarise(
      source_classification_contract_version = mode_value(.data$classification_contract_version),
      n_significant_microglia_terms = dplyr::n(),
      n_mixed_microenvironment = sum(.data$interpretation_class == "mixed_microenvironment", na.rm = TRUE),
      n_neuropil_sensitive = sum(.data$interpretation_class == "neuropil_sensitive", na.rm = TRUE),
      n_microglia_robust = sum(.data$interpretation_class == "microglia_robust", na.rm = TRUE),
      n_neuropil_marker_enriched = sum(.data$interpretation_class == "neuropil_marker_enriched", na.rm = TRUE),
      n_ambiguous = sum(.data$interpretation_class == "ambiguous", na.rm = TRUE),
      fraction_mixed_microenvironment = .data$n_mixed_microenvironment / .data$n_significant_microglia_terms,
      fraction_neuropil_sensitive = .data$n_neuropil_sensitive / .data$n_significant_microglia_terms,
      fraction_microglia_robust = .data$n_microglia_robust / .data$n_significant_microglia_terms,
      fraction_neuropil_marker_enriched = .data$n_neuropil_marker_enriched / .data$n_significant_microglia_terms,
      fraction_ambiguous = .data$n_ambiguous / .data$n_significant_microglia_terms,
      gene_overlap_fraction = safe_mean(.data$gene_overlap_fraction),
      gene_jaccard = safe_mean(.data$gene_jaccard),
      neuropil_marker_fraction = safe_mean(.data$neuropil_marker_fraction),
      microglia_marker_fraction = safe_mean(.data$microglia_marker_fraction),
      interpretation_class = mode_value(.data$interpretation_class),
      dominant_interpretation_class = .data$interpretation_class,
      neuropil_program_annotation_status = "significant_microglia_terms",
      neuropil_program_annotation_contract_version = NEUROPIL_PROGRAM_ANNOTATION_CONTRACT_VERSION,
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      .data$dataset, .data$comparison, .data$route_category, .data$route_unit,
      .data$biological_program
    )
}

attach_neuropil_program_annotation <- function(program_summary, neuropil_annotation) {
  join_keys <- c("dataset", "comparison", "route_category", "route_unit", "biological_program")
  missing <- setdiff(join_keys, names(program_summary))
  if (length(missing)) {
    stop("Program summary is missing join columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  annotation_summary <- summarise_neuropil_annotation_by_program(neuropil_annotation)
  count_columns <- c(
    "n_significant_microglia_terms", "n_mixed_microenvironment", "n_neuropil_sensitive",
    "n_microglia_robust", "n_neuropil_marker_enriched", "n_ambiguous"
  )
  program_summary %>%
    dplyr::left_join(annotation_summary, by = join_keys) %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(count_columns), ~dplyr::coalesce(as.integer(.x), 0L)),
      interpretation_class = dplyr::if_else(
        .data$n_significant_microglia_terms == 0L,
        "no_significant_microglia_terms",
        .data$interpretation_class
      ),
      dominant_interpretation_class = dplyr::if_else(
        .data$n_significant_microglia_terms == 0L,
        "no_significant_microglia_terms",
        .data$dominant_interpretation_class
      ),
      neuropil_program_annotation_status = dplyr::if_else(
        .data$n_significant_microglia_terms == 0L,
        "no_significant_microglia_terms",
        .data$neuropil_program_annotation_status
      ),
      neuropil_program_annotation_contract_version = NEUROPIL_PROGRAM_ANNOTATION_CONTRACT_VERSION
    )
}

add_unavailable_neuropil_program_annotation <- function(program_summary) {
  program_summary %>%
    dplyr::mutate(
      interpretation_class = NA_character_,
      dominant_interpretation_class = NA_character_,
      neuropil_program_annotation_status = "neuropil_annotation_unavailable",
      neuropil_program_annotation_contract_version = NEUROPIL_PROGRAM_ANNOTATION_CONTRACT_VERSION,
      source_classification_contract_version = NA_character_,
      n_significant_microglia_terms = NA_integer_,
      n_mixed_microenvironment = NA_integer_,
      n_neuropil_sensitive = NA_integer_,
      n_microglia_robust = NA_integer_,
      n_neuropil_marker_enriched = NA_integer_,
      n_ambiguous = NA_integer_,
      fraction_mixed_microenvironment = NA_real_,
      fraction_neuropil_sensitive = NA_real_,
      fraction_microglia_robust = NA_real_,
      fraction_neuropil_marker_enriched = NA_real_,
      fraction_ambiguous = NA_real_,
      gene_overlap_fraction = NA_real_,
      gene_jaccard = NA_real_,
      neuropil_marker_fraction = NA_real_,
      microglia_marker_fraction = NA_real_
    )
}

neuropil_annotation_path <- file.path(
  path_results("tables", MODULE_ID, "neuropil_reference_annotation", DATASET),
  "microglia_neuropil_annotation_latest.csv"
)
microglia_signature_path <- file.path(
  path_results("tables", MODULE_ID, "microglia_targeted_signature_enrichment", DATASET),
  "microglia_signature_enrichment_with_contrast_class.csv"
)

record_input_resolution(
  script = SCRIPT_ID,
  dataset = DATASET,
  stage = "enrichment",
  input_name = "neuropil_reference_annotation",
  expected_path = file.path(path_results("tables", MODULE_ID, "neuropil_reference_annotation", DATASET), "microglia_neuropil_annotation_latest.csv"),
  resolved_path = neuropil_annotation_path,
  resolution_mode = if (file.exists(neuropil_annotation_path)) "canonical" else "missing",
  strict_mode = strict_inputs_enabled(),
  allowed_in_strict_mode = TRUE,
  producer_script_or_artifact_id = "04_differential_expression_enrichment/neuropil_reference_annotation"
)
record_input_resolution(
  script = SCRIPT_ID,
  dataset = DATASET,
  stage = "enrichment",
  input_name = "microglia_targeted_signature_annotation",
  expected_path = file.path(path_results("tables", MODULE_ID, "microglia_targeted_signature_enrichment", DATASET), "microglia_signature_enrichment_with_contrast_class.csv"),
  resolved_path = microglia_signature_path,
  resolution_mode = if (file.exists(microglia_signature_path)) "canonical" else "missing",
  strict_mode = strict_inputs_enabled(),
  allowed_in_strict_mode = TRUE,
  producer_script_or_artifact_id = "04_differential_expression_enrichment/05_microglia_targeted_signature_enrichment.r"
)

neuropil_annotation <- optional_read_csv(neuropil_annotation_path)
signature_annotation <- optional_read_csv(microglia_signature_path)

program_summary_neuropil <- program_summary
if (!is.null(neuropil_annotation) && nrow(neuropil_annotation) && "comparison" %in% names(neuropil_annotation)) {
  program_summary_neuropil <- attach_neuropil_program_annotation(program_summary, neuropil_annotation)
} else {
  program_summary_neuropil <- add_unavailable_neuropil_program_annotation(program_summary_neuropil)
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
      .data$microglia_signature_class %in% c("microglia_enriched_empirical", "microglia_enriched_reference_supported") & !.data$interpretation_class %in% c("neuropil_sensitive", "neuropil_marker_enriched") ~ "microglia_supported_program",
      .data$microglia_signature_class == "curated_microglia_program" & !.data$interpretation_class %in% c("neuropil_sensitive", "neuropil_marker_enriched") ~ "curated_microglia_relevant_program",
      .data$microglia_signature_class == "neuropil_shared" | .data$interpretation_class %in% c("neuropil_sensitive", "neuropil_marker_enriched") ~ "neuropil_shared_or_sensitive_program",
      .data$microglia_signature_class == "mixed_microenvironment" | .data$interpretation_class == "mixed_microenvironment" ~ "mixed_microenvironment_program",
      is.na(.data$microglia_signature_class) & .data$neuropil_program_annotation_status == "no_significant_microglia_terms" ~ "no_significant_microglia_terms",
      is.na(.data$microglia_signature_class) & is.na(.data$interpretation_class) ~ "unannotated_program",
      TRUE ~ "ambiguous_program"
    ),
    interpretation_note = dplyr::case_when(
      .data$integrated_interpretation == "microglia_supported_program" ~ "Supported by microglia-targeted signatures with weak/absent or weaker neuropil reference signal.",
      .data$integrated_interpretation == "curated_microglia_relevant_program" ~ "Curated microglia-relevant gene set; not microglia-specific or claim-ready without empirical/reference support.",
      .data$integrated_interpretation == "neuropil_shared_or_sensitive_program" ~ "Shared with or sensitive to neuropil reference; do not interpret as microglia-intrinsic without orthogonal support.",
      .data$integrated_interpretation == "mixed_microenvironment_program" ~ "Present in microglia ROI and neuropil reference in a pattern consistent with local microenvironment biology.",
      .data$integrated_interpretation == "no_significant_microglia_terms" ~ "No FDR-significant microglia GO terms mapped to this biological program.",
      .data$integrated_interpretation == "unannotated_program" ~ "No optional neuropil or microglia signature annotation was available.",
      TRUE ~ "Weak, discordant, broad, or partially missing annotation support."
    )
  )

readr::write_csv(program_summary_neuropil, file.path(PATHS$tables, "program_summary_neuropil_annotated.csv"), na = "")
readr::write_csv(program_summary_signature, file.path(PATHS$tables, "program_summary_microglia_signature_annotated.csv"), na = "")
readr::write_csv(program_summary_integrated, file.path(PATHS$tables, "program_summary_integrated_interpretation.csv"), na = "")

plot_df <- heatmap_ready %>% dplyr::rename(score = "signed_neg_log10_fdr")
if (nrow(plot_df)) {
  p <- ggplot2::ggplot(
    plot_df %>%
      dplyr::mutate(
        comparison_label = clean_comparison_label(.data$comparison_label),
        biological_program = clean_program_label(.data$biological_program)
      ),
    ggplot2::aes(x = .data$comparison_label, y = .data$biological_program, fill = .data$score)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", na.value = "grey90", name = "signed -log10 FDR") +
    ggplot2::labs(x = NULL, y = NULL, title = paste("Biological program atlas:", DATASET)) +
    theme_nature_heatmap(7, base_family = "sans")
  save_plot_dual(p, file.path(PATHS$figures, "program_atlas_heatmap.svg"), width = max(6, length(unique(plot_df$comparison_label)) * 0.35), height = 4.8)
}

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(
    manifest = manifest_file,
    term_comparison = bundle$term_source,
    term_gene_provenance = bundle$provenance_source,
    analysis_status = bundle$status_source
  ),
  outputs = list(
    program_summary = file.path(PATHS$tables, "program_summary.csv"),
    program_summary_wide = file.path(PATHS$tables, "program_summary_wide.csv"),
    heatmap_ready = file.path(PATHS$tables, "program_summary_heatmap_ready.csv"),
    evidence = file.path(PATHS$source_data, "program_term_gene_evidence.csv"),
    analysis_status = file.path(PATHS$tables, "program_analysis_status.csv"),
    heatmap = file.path(PATHS$figures, "program_atlas_heatmap.svg")
  ),
  parameters = list(dataset = DATASET, program_patterns = biological_program_patterns()),
  notes = "Program mapping uses conservative regex on enrichment Description; claims remain exploratory unless independently supported."
)
