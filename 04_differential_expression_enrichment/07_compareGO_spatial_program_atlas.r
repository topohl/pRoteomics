# ================================================================
# Script: 04_differential_expression_enrichment/07_compareGO_spatial_program_atlas.r
# Stage: enrichment
# Scope: dataset_specific
# Consumes: canonical compareGO manifest-declared term comparison, term-gene provenance, and analysis-status tables for each dataset; config/manuscript_go_theme_registry.tsv; optional finalized biological_program_summary/<dataset>/program_summary.csv for validation cross-checks.
# Produces: legacy broad heuristic spatial-program outputs plus ontology-aware SUS-RES manuscript-theme tables/audits under compareGO_spatial_atlas.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Spatial program atlas from compareGO outputs.
# ================================================================

#' Spatial Program Atlas from compareGO outputs
#'
#' Consumes existing compareGO outputs across dataset families and synthesizes
#' generic spatial program summaries and ontology-aware SUS-RES manuscript
#' theme outputs. This script does not rerun clusterProfiler or compareGO.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "plotting_nature.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "manuscript_go_theme_utils.R"))
dataset_config_file <- repo_path("R", "dataset_config.R")
if (!file.exists(dataset_config_file)) {
  stop("Missing required config: ", dataset_config_file, call. = FALSE)
}
source(dataset_config_file)

if (!exists("proteomics_dataset_order", inherits = TRUE)) {
  proteomics_dataset_order <- if (exists("valid_datasets", inherits = TRUE)) valid_datasets() else c("neuron_neuropil", "neuron_soma", "microglia")
}
if (!exists("dataset_label", inherits = TRUE)) {
  dataset_label <- function(dataset) {
    labels <- c(neuron_neuropil = "Neuron neuropil", neuron_soma = "Neuron soma", microglia = "Microglia")
    out <- unname(labels[dataset])
    ifelse(is.na(out), dataset, out)
  }
}
if (!exists("dataset_compartment", inherits = TRUE)) {
  dataset_compartment <- function(dataset) {
    compartments <- c(neuron_neuropil = "neuropil", neuron_soma = "soma", microglia = "microglia")
    out <- unname(compartments[dataset])
    ifelse(is.na(out), dataset, out)
  }
}

MODULE_ID <- "04_differential_expression_enrichment"
args <- commandArgs(trailingOnly = TRUE)
VALIDATION_ONLY <- "--validation-only" %in% args || tolower(Sys.getenv("PROTEOMICS_SPATIAL_ATLAS_VALIDATION_ONLY", unset = "")) %in% c("1", "true", "yes", "y")
SUBSTEP_ID <- if (VALIDATION_ONLY) "compareGO_spatial_atlas_validation_proposed" else "compareGO_spatial_atlas"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

required_pkgs <- c("dplyr", "tidyr", "purrr", "readr", "writexl", "ggplot2", "stringr", "tibble", "scales")
MANUSCRIPT_GO_THEME_REGISTRY <- repo_path("config", "manuscript_go_theme_registry.tsv")
MANUSCRIPT_GO_SEMANTIC_CUTOFF <- 0.70
ATLAS_RESULT_TYPE <- "GSEA_GO"
ATLAS_ONTOLOGY <- "BP"
ATLAS_ROUTE_CATEGORY <- "phenotype_within_unit"
ATLAS_INPUT_CONTRACT_VERSION <- "compareGO_spatial_atlas_v3_canonical_bundle_claim_safe_core"

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

load_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required R package(s): ", paste(missing, collapse = ", "),
         ". Install them explicitly before running this script.", call. = FALSE)
  }
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}

canonical_program_levels <- biological_program_patterns()$biological_program
program_levels <- c(canonical_program_levels, "Other")
contrast_levels <- c("RES_vs_CON", "SUS_vs_CON", "SUS_vs_RES")
region_levels <- c("CA1", "CA2", "CA3", "DG")
publication_program_labels <- clean_program_label(canonical_program_levels)
atlas_min_recurrent_units <- as.integer(Sys.getenv("PROTEOMICS_PUBLICATION_ATLAS_MIN_RECURRENT_UNITS", unset = "2"))
divergence_delta_threshold <- as.numeric(Sys.getenv("PROTEOMICS_PUBLICATION_DIVERGENCE_DELTA_NES", unset = "0.25"))

is_truthy_env <- function(name) {
  tolower(Sys.getenv(name, unset = "")) %in% c("1", "true", "yes", "y")
}

is_script_dry_run <- function() {
  is_dry_run() || is_truthy_env("PROTEOMICS_SPATIAL_ATLAS_DRY_RUN")
}

atlas_theme <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70", linewidth = 0.3),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      panel.spacing = grid::unit(0.8, "lines")
    )
}

get_requested_datasets <- function() {
  raw <- Sys.getenv("PROTEOMICS_SPATIAL_ATLAS_DATASETS", unset = "neuron_neuropil,neuron_soma,microglia")
  datasets <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  datasets[nzchar(datasets)]
}

discover_dataset_files <- function(dataset) {
  manifest_path <- canonical_comparego_manifest_path(dataset)
  manifest <- if (file.exists(manifest_path)) {
    read_canonical_comparego_manifest(
      manifest_path, dataset, require_files = FALSE, repository_root = repo_path()
    )
  } else {
    NULL
  }
  declared_path <- function(column) {
    if (is.null(manifest) || !column %in% names(manifest)) return(NA_character_)
    values <- sort(unique(as.character(manifest[[column]][!is.na(manifest[[column]]) & nzchar(manifest[[column]])])), method = "radix")
    if (length(values) != 1L) return(NA_character_)
    values[[1]]
  }
  list(
    dataset = dataset,
    manifest = manifest_path,
    term_comparison = declared_path("term_comparison_file"),
    term_gene_provenance = declared_path("term_gene_provenance_output_file"),
    analysis_status = declared_path("analysis_status_summary_file")
  )
}

diagnose_inputs <- function(files) {
  paths <- c(files$manifest, files$term_comparison, files$term_gene_provenance, files$analysis_status)
  tibble::tibble(
    dataset = files$dataset,
    input_type = c("canonical_comparego_manifest", "term_comparison", "term_gene_provenance", "analysis_status"),
    path = paths,
    exists = !is.na(paths) & nzchar(paths) & file.exists(paths),
    input_contract_version = ATLAS_INPUT_CONTRACT_VERSION
  )
}

collapse_sorted_unique <- function(x) {
  values <- sort(unique(trimws(as.character(x))), method = "radix")
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values)) paste(values, collapse = ";") else NA_character_
}

canonical_term_key <- function(x) {
  paste(x$dataset, x$comparison, x$result_type, x$ontology, x$term_id, sep = "\r")
}

summarise_claim_safe_core <- function(provenance) {
  required <- c(
    "dataset", "comparison", "result_type", "ontology", "term_id",
    "official_gene_symbol", "ProteinGroupID", "gene_level_claim_allowed",
    "core_enrichment_member"
  )
  missing <- setdiff(required, names(provenance))
  if (length(missing)) {
    stop("Canonical term-gene provenance is missing stage-07 columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  safe <- provenance %>%
    dplyr::filter(
      .data$result_type == ATLAS_RESULT_TYPE,
      toupper(.data$ontology) == ATLAS_ONTOLOGY,
      .data$gene_level_claim_allowed %in% TRUE,
      .data$core_enrichment_member %in% TRUE,
      !is.na(.data$official_gene_symbol), nzchar(trimws(.data$official_gene_symbol))
    )
  summary <- safe %>%
    dplyr::group_by(.data$dataset, .data$comparison, .data$result_type, .data$ontology, .data$term_id) %>%
    dplyr::summarise(
      core_enrichment_gene = collapse_sorted_unique(.data$official_gene_symbol),
      core_enrichment = collapse_sorted_unique(.data$ProteinGroupID),
      n_claim_safe_core_rows = dplyr::n(),
      n_claim_safe_core_genes = dplyr::n_distinct(.data$official_gene_symbol),
      n_claim_safe_core_protein_groups = dplyr::n_distinct(.data$ProteinGroupID[!is.na(.data$ProteinGroupID) & nzchar(.data$ProteinGroupID)]),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$dataset, .data$comparison, .data$result_type, .data$ontology, .data$term_id)
  list(
    summary = summary,
    n_rows = nrow(safe),
    n_unique_genes = dplyr::n_distinct(safe$official_gene_symbol),
    n_unique_protein_groups = dplyr::n_distinct(safe$ProteinGroupID[!is.na(safe$ProteinGroupID) & nzchar(safe$ProteinGroupID)])
  )
}

parse_unit <- function(unit, dataset) {
  unit <- as.character(unit)
  unit <- gsub("_", "", unit)
  m <- stringr::str_match(unit, "^(CA1|CA2|CA3|DG)(.*)$")
  region <- ifelse(is.na(m[, 2]), NA_character_, m[, 2])
  layer <- ifelse(is.na(m[, 3]) | !nzchar(m[, 3]), NA_character_, m[, 3])
  if (identical(dataset, "microglia")) layer <- NA_character_
  tibble::tibble(region = region, layer = layer)
}

parse_comparison_name <- function(comparison, dataset, route_unit = NA_character_) {
  parts <- strsplit(as.character(comparison), "_", fixed = TRUE)[[1]]
  left <- parts[[1]]
  right <- if (length(parts) >= 2) parts[[2]] else NA_character_

  parse_side <- function(x) {
    m <- stringr::str_match(tolower(x), "^(.*?)(con|res|sus)$")
    tibble::tibble(unit = m[, 2], phenotype = toupper(m[, 3]))
  }

  left_parsed <- parse_side(left)
  right_parsed <- parse_side(right)
  unit <- if (!is.na(route_unit) && nzchar(as.character(route_unit))) as.character(route_unit) else left_parsed$unit
  unit_parsed <- parse_unit(unit, dataset)
  contrast <- paste0(left_parsed$phenotype, "_vs_", right_parsed$phenotype)

  tibble::tibble(
    Comparison = comparison,
    phenotype_contrast = contrast,
    region = unit_parsed$region,
    layer = unit_parsed$layer,
    compartment = dataset_compartment(dataset),
    spatial_unit = ifelse(is.na(unit_parsed$layer), unit_parsed$region, paste(unit_parsed$region, unit_parsed$layer, sep = "_"))
  )
}

classify_program <- function(description) {
  mapped <- map_terms_to_programs(
    data.frame(Description = as.character(description), stringsAsFactors = FALSE),
    description_col = "Description"
  )
  out <- as.character(mapped$biological_program)
  out[is.na(out) | !nzchar(out)] <- "Other"
  factor(out, levels = program_levels)
}

split_genes <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  unique(trimws(unlist(strsplit(x, "/|;|,|\\s+"))))
}

collapse_top_terms <- function(df, n = 5) {
  df %>%
    dplyr::arrange(dplyr::desc(!is.na(.data$p.adjust) & .data$p.adjust < 0.05), .data$p.adjust, dplyr::desc(abs(.data$NES))) %>%
    dplyr::pull(.data$Description) %>%
    unique() %>%
    utils::head(n) %>%
    paste(collapse = "; ")
}

top_driver_genes_for_group <- function(core_gene_sets, n = 10L) {
  genes <- unlist(lapply(core_gene_sets, split_genes), use.names = FALSE)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) return(NA_character_)
  counts <- as.data.frame(table(genes), stringsAsFactors = FALSE)
  counts <- counts[order(-counts$Freq, counts$genes, method = "radix"), , drop = FALSE]
  paste(utils::head(as.character(counts$genes), n), collapse = "; ")
}

load_dataset_atlas <- function(files) {
  bundle <- read_canonical_comparego_bundle(files$manifest, files$dataset)
  admitted_manifest <- bundle$manifest %>%
    dplyr::filter(
      .data$route_category == ATLAS_ROUTE_CATEGORY,
      .data$result_type == ATLAS_RESULT_TYPE,
      toupper(.data$ontology) == ATLAS_ONTOLOGY,
      .data$comparego_analysis_status == "included"
    )
  if (!nrow(admitted_manifest)) {
    stop("Canonical compareGO manifest has no admitted GSEA_GO/BP phenotype-within-unit rows for ", files$dataset, ".", call. = FALSE)
  }
  route_identity <- admitted_manifest %>%
    dplyr::select(.data$dataset, .data$comparison, .data$result_type, .data$ontology,
      .data$route_category, .data$route_unit) %>%
    dplyr::distinct()
  route_keys <- paste(route_identity$dataset, route_identity$comparison, route_identity$result_type, route_identity$ontology, sep = "\r")
  if (anyDuplicated(route_keys)) {
    stop("Canonical compareGO manifest has non-unique structured route identity for ", files$dataset, ".", call. = FALSE)
  }

  terms <- bundle$terms %>%
    dplyr::filter(.data$result_type == ATLAS_RESULT_TYPE, toupper(.data$ontology) == ATLAS_ONTOLOGY)
  term_comparisons <- sort(unique(as.character(terms$comparison)), method = "radix")
  manifest_comparisons <- sort(unique(as.character(admitted_manifest$comparison)), method = "radix")
  if (!identical(term_comparisons, manifest_comparisons)) {
    stop("Canonical compareGO term comparisons do not exactly match the admitted manifest comparisons for ", files$dataset, ".", call. = FALSE)
  }
  term_keys <- canonical_term_key(terms)
  if (anyDuplicated(term_keys)) {
    stop("Canonical compareGO term comparison contains duplicate exact term keys for ", files$dataset, ".", call. = FALSE)
  }

  core <- summarise_claim_safe_core(bundle$provenance)
  parsed <- terms %>%
    dplyr::rename(canonical_term_core_enrichment_raw = .data$core_enrichment) %>%
    dplyr::left_join(
      route_identity,
      by = c("dataset", "comparison", "result_type", "ontology")
    ) %>%
    dplyr::left_join(
      core$summary,
      by = c("dataset", "comparison", "result_type", "ontology", "term_id")
    ) %>%
    dplyr::mutate(
      dataset = files$dataset,
      dataset_label = dataset_label(files$dataset),
      ID = as.character(.data$term_id),
      NES = suppressWarnings(as.numeric(.data$NES)),
      p.adjust = suppressWarnings(as.numeric(.data$p.adjust)),
      Comparison = as.character(.data$comparison),
      Description = as.character(.data$term_description),
      core_enrichment = as.character(.data$core_enrichment),
      core_enrichment_gene = as.character(.data$core_enrichment_gene),
      source_term_comparison_file = bundle$term_source,
      source_term_gene_provenance_file = bundle$provenance_source,
      source_analysis_status_file = bundle$status_source,
      source_manifest_file = bundle$manifest_source,
      source_supplementary_file = .data$source_term_comparison_file,
      source_driver_file = .data$source_term_gene_provenance_file,
      evidence_source_family = "canonical_compareGO_ranked_GSEA_GO_BP",
      core_enrichment_contract = "gene_level_claim_allowed == TRUE && core_enrichment_member == TRUE",
      input_contract_version = ATLAS_INPUT_CONTRACT_VERSION
    )
  if (any(is.na(parsed$route_unit) | !nzchar(as.character(parsed$route_unit)))) {
    stop("Canonical compareGO terms are missing structured route_unit after manifest join for ", files$dataset, ".", call. = FALSE)
  }

  comparison_meta <- purrr::pmap_dfr(
    list(parsed$Comparison, parsed$dataset, parsed$route_unit),
    parse_comparison_name
  ) %>%
    dplyr::distinct(.data$Comparison, .keep_all = TRUE)

  parsed <- parsed %>%
    dplyr::left_join(comparison_meta, by = "Comparison") %>%
    dplyr::filter(.data$phenotype_contrast %in% contrast_levels) %>%
    dplyr::mutate(
      program_class = classify_program(.data$Description),
      spatial_unit_order = factor(.data$spatial_unit, levels = order_spatial_units(unique(.data$spatial_unit))),
      phenotype_contrast = factor(.data$phenotype_contrast, levels = contrast_levels),
      phenotype_contrast_source = "standardized_comparison_name_no_structured_condition_pair_available",
      route_identity_source = "canonical_compareGO_manifest"
    )

  list(
    enrichment = parsed,
    manifest = admitted_manifest,
    status = bundle$status,
    claim_safe_core = core,
    sources = tibble::tibble(
      dataset = files$dataset,
      manifest = bundle$manifest_source,
      term_comparison = bundle$term_source,
      term_gene_provenance = bundle$provenance_source,
      analysis_status = bundle$status_source
    )
  )
}

order_spatial_units <- function(units) {
  anatomical_spatial_unit_levels(units)
}

calculate_program_summary <- function(enrichment_df) {
  group_columns <- c(
    "dataset", "dataset_label", "comparison", "route_category", "route_unit",
    "region", "layer", "spatial_unit", "compartment", "phenotype_contrast", "program_class"
  )
  collapse_core <- function(x) {
    values <- unique(trimws(unlist(lapply(x, split_genes), use.names = FALSE)))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values)) paste(values, collapse = "; ") else NA_character_
  }

  base_summary <- enrichment_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) %>%
    dplyr::summarise(
      n_terms = dplyr::n_distinct(.data$ID),
      n_sig_terms = dplyr::n_distinct(.data$ID[!is.na(.data$p.adjust) & .data$p.adjust < 0.05]),
      n_positive_sig_terms = dplyr::n_distinct(.data$ID[!is.na(.data$p.adjust) & .data$p.adjust < 0.05 & .data$NES > 0]),
      n_negative_sig_terms = dplyr::n_distinct(.data$ID[!is.na(.data$p.adjust) & .data$p.adjust < 0.05 & .data$NES < 0]),
      min_fdr = suppressWarnings(min(.data$p.adjust, na.rm = TRUE)),
      mean_NES = mean(.data$NES, na.rm = TRUE),
      median_NES = median(.data$NES, na.rm = TRUE),
      mean_signed_log10FDR = mean(sign(.data$NES) * -log10(pmax(.data$p.adjust, .Machine$double.xmin)), na.rm = TRUE),
      leading_edge_union_size = length(unique(unlist(purrr::map(.data$core_enrichment, split_genes)))),
      significant_leading_edge_proteins = collapse_core(.data$core_enrichment[!is.na(.data$p.adjust) & .data$p.adjust < 0.05]),
      significant_leading_edge_genes = collapse_core(.data$core_enrichment_gene[!is.na(.data$p.adjust) & .data$p.adjust < 0.05]),
      top_driver_genes = top_driver_genes_for_group(.data$core_enrichment_gene[!is.na(.data$p.adjust) & .data$p.adjust < 0.05]),
      top_GO_terms = collapse_top_terms(dplyr::pick(dplyr::everything())),
      comparisons = paste(unique(.data$Comparison), collapse = ";"),
      source_term_comparison_file = dplyr::first(.data$source_term_comparison_file),
      source_term_gene_provenance_file = dplyr::first(.data$source_term_gene_provenance_file),
      source_analysis_status_file = dplyr::first(.data$source_analysis_status_file),
      source_supplementary_file = dplyr::first(.data$source_supplementary_file),
      source_driver_file = dplyr::first(.data$source_driver_file),
      source_manifest_file = dplyr::first(.data$source_manifest_file),
      evidence_source_family = dplyr::first(.data$evidence_source_family),
      core_enrichment_contract = dplyr::first(.data$core_enrichment_contract),
      input_contract_version = dplyr::first(.data$input_contract_version),
      .groups = "drop"
    )

  representative <- enrichment_df %>%
    dplyr::filter(!is.na(.data$p.adjust), .data$p.adjust < 0.05) %>%
    dplyr::mutate(representative_abs_NES = abs(.data$NES)) %>%
    dplyr::arrange(
      dplyr::across(dplyr::all_of(group_columns)),
      .data$p.adjust, dplyr::desc(.data$representative_abs_NES), .data$ID, .data$Description
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) %>%
    dplyr::slice_head(n = 1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      dplyr::across(dplyr::all_of(group_columns)),
      representative_term_id = as.character(.data$ID),
      representative_term = as.character(.data$Description),
      representative_NES = as.numeric(.data$NES),
      representative_FDR = as.numeric(.data$p.adjust),
      representative_leading_edge_proteins = as.character(.data$core_enrichment),
      representative_leading_edge_genes = as.character(.data$core_enrichment_gene),
      representative_source_comparison = as.character(.data$Comparison),
      representative_source_key = paste(.data$dataset, .data$Comparison, .data$ID, sep = "|")
    )

  base_summary %>%
    dplyr::left_join(representative, by = group_columns) %>%
    dplyr::mutate(
      min_fdr = ifelse(is.infinite(.data$min_fdr), NA_real_, .data$min_fdr),
      direction_consistency = dplyr::case_when(
        .data$n_sig_terms == 0 ~ "no_FDR_supported_term",
        .data$n_positive_sig_terms > 0 & .data$n_negative_sig_terms > 0 ~ "mixed_direction",
        .data$n_positive_sig_terms > 0 ~ "positive_only",
        .data$n_negative_sig_terms > 0 ~ "negative_only",
        TRUE ~ "neutral_or_undirected"
      ),
      representative_direction = dplyr::case_when(
        is.na(.data$representative_NES) ~ "no_FDR_supported_term",
        .data$representative_NES > 0 ~ "positive_NES",
        .data$representative_NES < 0 ~ "negative_NES",
        TRUE ~ "neutral"
      ),
      representative_selection_rule = "p.adjust < 0.05; lowest p.adjust; largest abs(NES); GO ID; Description"
    )
}

# Backward-compatible legacy heuristic only. These classes mix significance
# status and arbitrary mean-NES cutoffs and must not be used as manuscript
# evidence for a difference between RES-vs-CON and SUS-vs-CON.
classify_program_behavior <- function(summary_df, min_abs_nes = 0.15) {
  wide <- summary_df %>%
    dplyr::group_by(.data$dataset, .data$dataset_label, .data$region, .data$layer, .data$spatial_unit,
                    .data$compartment, .data$program_class, .data$phenotype_contrast) %>%
    dplyr::summarise(mean_NES = mean(.data$mean_NES, na.rm = TRUE), n_sig_terms = sum(.data$n_sig_terms), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = "phenotype_contrast",
      values_from = c("mean_NES", "n_sig_terms"),
      values_fill = list(mean_NES = NA_real_, n_sig_terms = 0)
    )

  get_col <- function(df, nm, default = NA_real_) {
    if (nm %in% names(df)) df[[nm]] else rep(default, nrow(df))
  }

  res <- get_col(wide, "mean_NES_RES_vs_CON")
  sus <- get_col(wide, "mean_NES_SUS_vs_CON")
  sus_res <- get_col(wide, "mean_NES_SUS_vs_RES")
  res_sig <- get_col(wide, "n_sig_terms_RES_vs_CON", 0)
  sus_sig <- get_col(wide, "n_sig_terms_SUS_vs_CON", 0)
  sr_sig <- get_col(wide, "n_sig_terms_SUS_vs_RES", 0)

  wide %>%
    dplyr::mutate(
      RES_vs_CON_mean_NES = res,
      SUS_vs_CON_mean_NES = sus,
      SUS_vs_RES_mean_NES = sus_res,
      behavior_class = dplyr::case_when(
        !is.na(res) & !is.na(sus) & abs(res) >= min_abs_nes & abs(sus) >= min_abs_nes & sign(res) != sign(sus) ~ "divergent",
        (sr_sig > 0 | (!is.na(sus_res) & abs(sus_res) >= min_abs_nes)) & !is.na(res) & !is.na(sus) ~ "phenotype_separating",
        res_sig > 0 & sus_sig > 0 & !is.na(res) & !is.na(sus) & sign(res) == sign(sus) ~ "shared_stress",
        res_sig > 0 & (sus_sig == 0 | is.na(sus) | abs(sus) < min_abs_nes) ~ "resilience_specific",
        sus_sig > 0 & (res_sig == 0 | is.na(res) | abs(res) < min_abs_nes) ~ "susceptibility_specific",
        TRUE ~ "phenotype_separating"
      ),
      legacy_heuristic_classification = TRUE,
      manuscript_evidence_eligible = FALSE,
      legacy_classification_warning = "Descriptive backward-compatibility field; significance in one contrast and nonsignificance in another does not establish a between-contrast difference."
    )
}

make_driver_recurrence <- function(enrichment_df, max_genes_per_class = 20) {
  long_genes <- enrichment_df %>%
    dplyr::filter(!is.na(.data$p.adjust), .data$p.adjust < 0.05) %>%
    dplyr::mutate(Gene = purrr::map(.data$core_enrichment_gene, split_genes)) %>%
    tidyr::unnest("Gene") %>%
    dplyr::filter(!is.na(.data$Gene), nzchar(.data$Gene)) %>%
    dplyr::distinct(.data$dataset, .data$dataset_label, .data$spatial_unit, .data$program_class, .data$Gene, .data$Description)

  top_genes <- long_genes %>%
    dplyr::group_by(.data$program_class, .data$Gene) %>%
    dplyr::summarise(n_spatial_units = dplyr::n_distinct(paste(.data$dataset, .data$spatial_unit)), n_terms = dplyr::n_distinct(.data$Description), .groups = "drop") %>%
    dplyr::arrange(.data$program_class, dplyr::desc(.data$n_spatial_units), dplyr::desc(.data$n_terms)) %>%
    dplyr::group_by(.data$program_class) %>%
    dplyr::slice_head(n = max_genes_per_class) %>%
    dplyr::ungroup()

  long_genes %>%
    dplyr::semi_join(top_genes, by = c("program_class", "Gene")) %>%
    dplyr::group_by(.data$program_class, .data$Gene, .data$dataset_label, .data$spatial_unit) %>%
    dplyr::summarise(n_terms = dplyr::n_distinct(.data$Description), .groups = "drop") %>%
    dplyr::mutate(unit_label = paste(.data$dataset_label, .data$spatial_unit, sep = " | "))
}

plot_spatial_program_atlas <- function(summary_df, output_file) {
  plot_df <- summary_df %>%
    dplyr::mutate(
      spatial_unit = factor(.data$spatial_unit, levels = order_spatial_units(unique(.data$spatial_unit))),
      program_class = factor(.data$program_class, levels = rev(program_levels)),
      phenotype_contrast = factor(.data$phenotype_contrast, levels = contrast_levels)
    )
  max_abs <- max(abs(plot_df$mean_NES), na.rm = TRUE)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$spatial_unit, y = .data$program_class, color = .data$mean_NES, size = .data$n_sig_terms)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(.data$dataset_label ~ .data$phenotype_contrast, scales = "free_x", space = "free_x") +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-max_abs, max_abs), name = "Mean NES") +
    ggplot2::scale_size_continuous(range = c(0.6, 5), name = "Sig. terms") +
    ggplot2::labs(x = NULL, y = NULL, title = "Spatial GO Program Atlas") +
    atlas_theme(7)
  ggplot2::ggsave(output_file, p, width = 12, height = 7, device = "svg", limitsize = FALSE)
}

plot_res_sus_divergence <- function(behavior_df, output_file) {
  plot_df <- behavior_df %>%
    dplyr::filter(!is.na(.data$RES_vs_CON_mean_NES), !is.na(.data$SUS_vs_CON_mean_NES)) %>%
    dplyr::mutate(program_class = factor(.data$program_class, levels = program_levels))
  if (nrow(plot_df) == 0) return(invisible(FALSE))
  lim <- max(abs(c(plot_df$RES_vs_CON_mean_NES, plot_df$SUS_vs_CON_mean_NES)), na.rm = TRUE)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$RES_vs_CON_mean_NES, .data$SUS_vs_CON_mean_NES, label = .data$spatial_unit, color = .data$behavior_class)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey75", linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.85, size = 1.8) +
    ggplot2::geom_text(size = 2, check_overlap = TRUE, vjust = -0.5) +
    ggplot2::facet_grid(.data$dataset_label ~ .data$program_class) +
    ggplot2::coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = "RES vs CON mean NES", y = "SUS vs CON mean NES", color = "Behavior", title = "Resilience and Susceptibility Program Divergence") +
    atlas_theme(7)
  ggplot2::ggsave(output_file, p, width = 13, height = 7, device = "svg", limitsize = FALSE)
}

plot_driver_recurrence <- function(driver_df, output_file) {
  if (nrow(driver_df) == 0) return(invisible(FALSE))
  driver_df <- driver_df %>%
    dplyr::mutate(
      gene_label = paste(.data$program_class, .data$Gene, sep = " | "),
      gene_label = factor(.data$gene_label, levels = rev(unique(.data$gene_label[order(as.character(.data$program_class), .data$Gene)]))),
      unit_label = factor(.data$unit_label, levels = unique(.data$unit_label))
    )
  p <- ggplot2::ggplot(driver_df, ggplot2::aes(.data$unit_label, .data$gene_label, fill = .data$n_terms)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "white", high = "#3B7EA1", name = "Terms") +
    ggplot2::labs(x = NULL, y = "Program | leading-edge protein", title = "Recurrent Leading-Edge Driver Proteins") +
    atlas_theme(6) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
  ggplot2::ggsave(output_file, p, width = 13, height = max(6, length(levels(driver_df$gene_label)) * 0.12), device = "svg", limitsize = FALSE)
}

plot_compartment_comparison <- function(summary_df, output_file) {
  plot_df <- summary_df %>%
    dplyr::group_by(.data$dataset_label, .data$compartment, .data$region, .data$phenotype_contrast, .data$program_class) %>%
    dplyr::summarise(mean_NES = mean(.data$mean_NES, na.rm = TRUE), n_sig_terms = sum(.data$n_sig_terms), .groups = "drop") %>%
    dplyr::mutate(
      region = factor(.data$region, levels = region_levels),
      program_class = factor(.data$program_class, levels = rev(program_levels))
    )
  max_abs <- max(abs(plot_df$mean_NES), na.rm = TRUE)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$region, .data$program_class, color = .data$mean_NES, size = .data$n_sig_terms)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(.data$dataset_label ~ .data$phenotype_contrast) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-max_abs, max_abs), name = "Mean NES") +
    ggplot2::scale_size_continuous(range = c(0.6, 5), name = "Sig. terms") +
    ggplot2::labs(x = NULL, y = NULL, title = "Compartment-Level Program Comparison by Region") +
    atlas_theme(7)
  ggplot2::ggsave(output_file, p, width = 10, height = 7, device = "svg", limitsize = FALSE)
}

publication_color_limits <- function(x, cap = 2.5) {
  lim <- suppressWarnings(stats::quantile(abs(x), probs = 0.98, na.rm = TRUE, names = FALSE))
  if (!is.finite(lim) || lim <= 0) lim <- suppressWarnings(max(abs(x), na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) lim <- 1
  lim <- min(lim, cap)
  c(-lim, lim)
}

sig_size_breaks <- function(x) {
  x <- sort(unique(as.numeric(x[is.finite(x) & x > 0])))
  if (!length(x)) return(c(1, 2, 3))
  unique(round(stats::quantile(x, probs = c(0, 0.5, 1), names = FALSE)))
}

prepare_publication_summary <- function(summary_df, min_recurrent_units = atlas_min_recurrent_units) {
  summary_df %>%
    dplyr::mutate(
      program = clean_program_label(.data$program_class),
      dataset_compartment = as.character(.data$dataset),
      dataset_compartment_label = clean_spatial_unit_label(.data$dataset_label),
      comparison = as.character(.data$phenotype_contrast),
      comparison_label = clean_comparison_label(.data$comparison),
      spatial_unit_label = clean_spatial_unit_label(.data$spatial_unit),
      region = as.character(.data$region),
      layer = as.character(.data$layer),
      significant_term_count = as.integer(.data$n_sig_terms),
      mean_NES = as.numeric(.data$mean_NES),
      min_recurrent_units = min_recurrent_units,
      interpretable_program = .data$program %in% publication_program_labels,
      nonzero_signal = .data$significant_term_count > 0
    ) %>%
    dplyr::group_by(.data$program) %>%
    dplyr::mutate(
      recurrent_units = dplyr::n_distinct(paste(.data$dataset_compartment, .data$spatial_unit)[.data$nonzero_signal]),
      recurrent_comparisons = dplyr::n_distinct(.data$comparison[.data$nonzero_signal]),
      passes_recurrence_filter = .data$recurrent_units >= .data$min_recurrent_units | .data$recurrent_comparisons >= .data$min_recurrent_units,
      publication_include = .data$interpretable_program & .data$nonzero_signal & .data$passes_recurrence_filter,
      filter_reason = dplyr::case_when(
        !.data$interpretable_program ~ "excluded_non_manuscript_program",
        !.data$nonzero_signal ~ "excluded_zero_significant_terms",
        !.data$passes_recurrence_filter ~ "excluded_below_recurrence_threshold",
        TRUE ~ "included"
      )
    ) %>%
    dplyr::ungroup()
}

prepare_sus_res_publication_summary <- function(summary_df, min_recurrent_units = atlas_min_recurrent_units) {
  summary_df %>%
    dplyr::filter(as.character(.data$phenotype_contrast) == "SUS_vs_RES") %>%
    dplyr::mutate(
      program = clean_program_label(.data$program_class),
      dataset_compartment = as.character(.data$dataset),
      dataset_compartment_label = clean_spatial_unit_label(.data$dataset_label),
      comparison = as.character(.data$phenotype_contrast),
      comparison_label = clean_comparison_label(.data$comparison),
      spatial_unit_label = clean_spatial_unit_label(.data$spatial_unit),
      region = as.character(.data$region),
      layer = as.character(.data$layer),
      significant_term_count = as.integer(.data$n_sig_terms),
      min_recurrent_units = min_recurrent_units,
      interpretable_program = .data$program %in% publication_program_labels,
      nonzero_signal = .data$significant_term_count > 0,
      dataset_compartment_label = dplyr::if_else(
        .data$dataset_compartment == "microglia",
        "Microglia-enriched ROI",
        .data$dataset_compartment_label
      )
    ) %>%
    dplyr::group_by(.data$program) %>%
    dplyr::mutate(
      sus_res_recurrent_units = dplyr::n_distinct(paste(.data$dataset_compartment, .data$spatial_unit)[.data$nonzero_signal]),
      sus_res_recurrent_datasets = dplyr::n_distinct(.data$dataset_compartment[.data$nonzero_signal]),
      sus_res_recurrent_compartments = .data$sus_res_recurrent_datasets,
      sus_res_is_recurrent = .data$sus_res_recurrent_units >= .data$min_recurrent_units,
      recurrence_annotation = dplyr::case_when(
        !.data$nonzero_signal ~ "no_FDR_supported_term",
        .data$sus_res_recurrent_datasets >= 2L ~ "recurrent_across_datasets",
        .data$sus_res_is_recurrent ~ "recurrent_across_spatial_units",
        TRUE ~ "spatially_specific"
      ),
      recurrent_units = .data$sus_res_recurrent_units,
      recurrent_comparisons = dplyr::n_distinct(.data$comparison[.data$nonzero_signal]),
      passes_recurrence_filter = .data$sus_res_is_recurrent,
      publication_include = .data$interpretable_program & .data$nonzero_signal,
      filter_reason = dplyr::case_when(
        !.data$interpretable_program ~ "excluded_noncanonical_program",
        !.data$nonzero_signal ~ "excluded_no_FDR_supported_term",
        .data$sus_res_is_recurrent ~ "included_recurrent",
        TRUE ~ "included_spatially_specific"
      ),
      direction = dplyr::case_when(
        is.na(.data$representative_NES) ~ "no FDR-supported term",
        .data$representative_NES > 0 ~ "higher in SUS",
        .data$representative_NES < 0 ~ "higher in RES",
        TRUE ~ "neutral"
      ),
      panel_data_basis = "ranked_proteome_wide_GO_GSEA",
      NES_interpretation = "positive = higher in SUS; negative = higher in RES",
      recurrence_inference_role = "descriptive_only_not_a_significance_gate"
    ) %>%
    dplyr::ungroup()
}

plot_spatial_program_atlas_publication <- function(summary_df, figure_file, source_file, min_recurrent_units = atlas_min_recurrent_units) {
  source_df <- prepare_publication_summary(summary_df, min_recurrent_units)
  readr::write_csv(source_df, source_file, na = "")
  plot_df <- source_df %>%
    dplyr::filter(.data$publication_include) %>%
    dplyr::mutate(
      spatial_unit_label = factor(.data$spatial_unit_label, levels = clean_spatial_unit_label(anatomical_spatial_unit_levels(unique(.data$spatial_unit)))),
      program = factor(.data$program, levels = rev(publication_program_labels)),
      comparison_label = factor(.data$comparison_label, levels = clean_comparison_label(contrast_levels)),
      dataset_compartment_label = factor(.data$dataset_compartment_label, levels = clean_spatial_unit_label(dataset_label(proteomics_dataset_order)))
    )
  if (!nrow(plot_df)) return(invisible(FALSE))
  lims <- publication_color_limits(plot_df$mean_NES)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$spatial_unit_label, .data$program, color = .data$mean_NES, size = .data$significant_term_count)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(.data$dataset_compartment_label ~ .data$comparison_label, scales = "free_x", space = "free_x") +
    ggplot2::scale_color_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = lims, oob = scales::squish, name = "Mean NES") +
    ggplot2::scale_size_continuous(range = c(1.0, 4.2), breaks = sig_size_breaks(plot_df$significant_term_count), name = "Sig. terms") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature_dotplot(7)
  save_nature_svg(p, figure_file, width_mm = 180, height_mm = 105)
}

plot_sus_res_spatial_program_atlas_publication <- function(summary_df, figure_file, source_file, min_recurrent_units = atlas_min_recurrent_units) {
  source_df <- prepare_sus_res_publication_summary(summary_df, min_recurrent_units)
  readr::write_csv(source_df, source_file, na = "")
  axis_df <- source_df %>%
    dplyr::filter(.data$interpretable_program) %>%
    dplyr::mutate(
      spatial_unit_label = factor(
        .data$spatial_unit_label,
        levels = clean_spatial_unit_label(anatomical_spatial_unit_levels(unique(.data$spatial_unit)))
      ),
      program = factor(.data$program, levels = rev(publication_program_labels)),
      dataset_compartment_label = factor(
        .data$dataset_compartment_label,
        levels = c("Neuron neuropil", "Neuron soma", "Microglia-enriched ROI")
      )
    )
  plot_df <- axis_df %>% dplyr::filter(.data$publication_include)
  if (!nrow(plot_df)) return(invisible(FALSE))
  lims <- publication_color_limits(plot_df$representative_NES)
  p <- ggplot2::ggplot(axis_df, ggplot2::aes(.data$spatial_unit_label, .data$program)) +
    ggplot2::geom_blank() +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(color = .data$representative_NES, size = .data$significant_term_count),
      alpha = 0.92
    ) +
    ggplot2::facet_wrap(~dataset_compartment_label, ncol = 1, scales = "free_x", strip.position = "right") +
    ggplot2::scale_color_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = lims, oob = scales::squish, name = "Representative\nterm NES") +
    ggplot2::scale_size_continuous(range = c(1.0, 4.2), breaks = sig_size_breaks(plot_df$significant_term_count), name = "Sig. terms") +
    ggplot2::labs(
      x = NULL, y = NULL,
      subtitle = "Positive NES = higher in SUS; negative NES = higher in RES",
      caption = paste(
        "Point = lowest-FDR supported GO term per program/unit (ties: largest |NES|); recurrence is descriptive.",
        "Proteome-wide ranked GSEA; not DAP-set ORA.",
        sep = "\n"
      )
    ) +
    theme_nature_dotplot(7)
  save_nature_svg(p, figure_file, width_mm = 183, height_mm = 110)
}

plot_sus_res_manuscript_theme_atlas <- function(theme_summary, figure_file, source_file) {
  readr::write_csv(theme_summary, source_file, na = "")
  plot_df <- theme_summary %>%
    dplyr::filter(.data$theme_role == "primary", .data$n_theme_terms_tested > 0) %>%
    dplyr::mutate(
      manuscript_theme = factor(
        .data$manuscript_theme,
        levels = rev(unique(.data$manuscript_theme[order(.data$display_order)]))
      ),
      spatial_unit_label = factor(
        .data$spatial_unit_label,
        levels = clean_spatial_unit_label(anatomical_spatial_unit_levels(unique(.data$spatial_unit)))
      ),
      dataset_compartment_label = factor(
        .data$dataset_compartment_label,
        levels = c("Neuron neuropil", "Neuron soma", "Microglia-enriched ROI")
      ),
      unit_direction_status = ifelse(.data$direction_consistency == "mixed_direction", "mixed supported directions", "supported direction not mixed"),
      representative_minus_log10_FDR = ifelse(
        .data$FDR_support_present,
        -log10(pmax(.data$representative_supported_FDR, .Machine$double.xmin)),
        NA_real_
      )
    )
  if (!nrow(plot_df)) return(invisible(FALSE))
  lims <- publication_color_limits(plot_df$median_NES_all_theme_terms)
  descriptive_only <- plot_df %>% dplyr::filter(!.data$FDR_support_present)
  supported <- plot_df %>% dplyr::filter(.data$FDR_support_present)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$spatial_unit_label, .data$manuscript_theme)) +
    ggplot2::geom_point(data = descriptive_only, ggplot2::aes(fill = .data$median_NES_all_theme_terms), shape = 21, size = 1.4, colour = "transparent", alpha = 0.45) +
    ggplot2::geom_point(data = supported, ggplot2::aes(fill = .data$median_NES_all_theme_terms, size = .data$representative_minus_log10_FDR, shape = .data$unit_direction_status), colour = "black", stroke = 0.35, alpha = 0.96) +
    ggplot2::facet_wrap(~dataset_compartment_label, ncol = 1, scales = "free_x", strip.position = "right") +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
      limits = lims, oob = scales::squish, name = "Median NES\n(all mapped terms)"
    ) +
    ggplot2::scale_size_continuous(range = c(2.2, 4.5), name = expression(-log[10]("representative GO BH FDR"))) +
    ggplot2::scale_shape_manual(values = c("supported direction not mixed" = 21, "mixed supported directions" = 23), name = "FDR-supported\nterm directions") +
    ggplot2::labs(
      x = NULL, y = NULL,
      subtitle = "Positive NES = higher in SUS; negative NES = higher in RES",
      caption = paste(
        "Color represents the median NES across all tested GO-BP terms assigned to each ontology theme and spatial unit.",
        "Outlined/emphasized points indicate at least one constituent GO term with original BH FDR < 0.05.",
        "Non-emphasized values are descriptive directional summaries and are not statistically significant theme-level effects.",
        sep = "\n"
      )
    ) +
    theme_nature_dotplot(7)
  save_nature_svg(p, figure_file, width_mm = 183, height_mm = 110)
}

plot_neuropil_focus_publication <- function(summary_df, figure_file, source_file, min_recurrent_units = 1L) {
  source_df <- prepare_publication_summary(summary_df, min_recurrent_units) %>%
    dplyr::filter(.data$dataset_compartment == "neuron_neuropil")
  if (!nrow(source_df)) {
    readr::write_csv(source_df, source_file, na = "")
    return(invisible(FALSE))
  }
  available_focus <- intersect(clean_comparison_label(contrast_levels), unique(source_df$comparison_label))
  if (length(available_focus)) source_df <- dplyr::filter(source_df, .data$comparison_label %in% available_focus)
  readr::write_csv(source_df, source_file, na = "")
  plot_df <- source_df %>%
    dplyr::filter(.data$interpretable_program, .data$nonzero_signal) %>%
    dplyr::mutate(
      spatial_unit_label = factor(.data$spatial_unit_label, levels = clean_spatial_unit_label(anatomical_spatial_unit_levels(unique(.data$spatial_unit)))),
      program = factor(.data$program, levels = rev(publication_program_labels)),
      comparison_label = factor(.data$comparison_label, levels = clean_comparison_label(contrast_levels))
    )
  if (!nrow(plot_df)) return(invisible(FALSE))
  lims <- publication_color_limits(plot_df$mean_NES)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$spatial_unit_label, .data$program, color = .data$mean_NES, size = .data$significant_term_count)) +
    ggplot2::geom_point(alpha = 0.92) +
    ggplot2::facet_wrap(~comparison_label, nrow = 1, scales = "free_x") +
    ggplot2::scale_color_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = lims, oob = scales::squish, name = "Mean NES") +
    ggplot2::scale_size_continuous(range = c(1.0, 4.0), breaks = sig_size_breaks(plot_df$significant_term_count), name = "Sig. terms") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature_dotplot(7)
  save_nature_svg(p, figure_file, width_mm = 180, height_mm = 72)
}

plot_compartment_comparison_publication <- function(summary_df, figure_file, source_file) {
  source_df <- prepare_publication_summary(summary_df, min_recurrent_units = 1L) %>%
    dplyr::filter(.data$interpretable_program, .data$nonzero_signal) %>%
    dplyr::group_by(.data$dataset_compartment, .data$dataset_compartment_label, .data$comparison, .data$comparison_label, .data$region, .data$program) %>%
    dplyr::summarise(
      mean_NES = mean(.data$mean_NES, na.rm = TRUE),
      significant_term_count = sum(.data$significant_term_count, na.rm = TRUE),
      enriched_units = dplyr::n_distinct(.data$spatial_unit[.data$significant_term_count > 0]),
      min_fdr = suppressWarnings(min(.data$min_fdr, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(min_fdr = ifelse(is.infinite(.data$min_fdr), NA_real_, .data$min_fdr))
  readr::write_csv(source_df, source_file, na = "")
  if (!nrow(source_df)) return(invisible(FALSE))
  plot_df <- source_df %>%
    dplyr::mutate(
      program = factor(.data$program, levels = rev(publication_program_labels)),
      region = factor(.data$region, levels = region_levels),
      dataset_compartment_label = factor(.data$dataset_compartment_label, levels = clean_spatial_unit_label(dataset_label(proteomics_dataset_order))),
      comparison_label = factor(.data$comparison_label, levels = clean_comparison_label(contrast_levels))
    )
  lims <- publication_color_limits(plot_df$mean_NES)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$region, .data$program, color = .data$mean_NES, size = .data$enriched_units)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(.data$dataset_compartment_label ~ .data$comparison_label) +
    ggplot2::scale_color_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = lims, oob = scales::squish, name = "Mean NES") +
    ggplot2::scale_size_continuous(range = c(1.0, 4.0), breaks = sig_size_breaks(plot_df$enriched_units), name = "Enriched units") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature_dotplot(7)
  save_nature_svg(p, figure_file, width_mm = 170, height_mm = 105)
}

plot_res_sus_divergence_publication <- function(behavior_df, figure_file, source_file, delta_threshold = divergence_delta_threshold) {
  source_df <- behavior_df %>%
    dplyr::mutate(
      program = clean_program_label(.data$program_class),
      dataset_compartment = as.character(.data$dataset),
      dataset_compartment_label = clean_spatial_unit_label(.data$dataset_label),
      spatial_unit_label = clean_spatial_unit_label(.data$spatial_unit),
      RES_value = as.numeric(.data$RES_vs_CON_mean_NES),
      SUS_value = as.numeric(.data$SUS_vs_CON_mean_NES),
      delta = .data$SUS_value - .data$RES_value,
      abs_delta = abs(.data$delta),
      divergence_threshold = delta_threshold,
      divergence_class = dplyr::case_when(
        is.na(.data$RES_value) | is.na(.data$SUS_value) ~ "weak/ambiguous",
        abs(.data$RES_value) < 0.15 & abs(.data$SUS_value) < 0.15 ~ "weak/ambiguous",
        sign(.data$RES_value) != sign(.data$SUS_value) & abs(.data$RES_value) >= 0.15 & abs(.data$SUS_value) >= 0.15 ~ "opposite direction",
        .data$RES_value >= 0.15 & .data$SUS_value >= 0.15 & .data$abs_delta < delta_threshold ~ "shared up",
        .data$RES_value <= -0.15 & .data$SUS_value <= -0.15 & .data$abs_delta < delta_threshold ~ "shared down",
        .data$RES_value - .data$SUS_value >= delta_threshold ~ "RES-biased",
        .data$SUS_value - .data$RES_value >= delta_threshold ~ "SUS-biased",
        TRUE ~ "weak/ambiguous"
      ),
      label_candidate = .data$program %in% publication_program_labels & is.finite(.data$abs_delta)
    )
  readr::write_csv(source_df, source_file, na = "")
  plot_df <- source_df %>%
    dplyr::filter(.data$program %in% publication_program_labels, !is.na(.data$RES_value), !is.na(.data$SUS_value))
  if (!nrow(plot_df)) return(invisible(FALSE))
  labels_df <- plot_df %>%
    dplyr::arrange(dplyr::desc(.data$abs_delta)) %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::mutate(point_label = paste(.data$spatial_unit_label, .data$program, sep = "\n"))
  lim <- max(abs(c(plot_df$RES_value, plot_df$SUS_value)), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$RES_value, .data$SUS_value, color = .data$divergence_class)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey65", linewidth = 0.25, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.85, size = 1.8) +
    ggplot2::scale_color_manual(
      values = c("shared up" = "#B40426", "shared down" = "#3B4CC0", "RES-biased" = "#009E73", "SUS-biased" = "#D55E00", "opposite direction" = "#7A3E9D", "weak/ambiguous" = "grey65"),
      name = NULL
    ) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = "RES vs CON mean NES", y = "SUS vs CON mean NES") +
    theme_nature_dotplot(7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = labels_df, ggplot2::aes(label = .data$point_label), size = 2, min.segment.length = 0, max.overlaps = 12, show.legend = FALSE)
  } else {
    p <- p + ggplot2::geom_text(data = labels_df, ggplot2::aes(label = .data$point_label), size = 2, check_overlap = TRUE, vjust = -0.6, show.legend = FALSE)
  }
  save_nature_svg(p, figure_file, width_mm = 90, height_mm = 82)
}

plot_driver_recurrence_topdrivers_publication <- function(driver_df, figure_file, source_file, max_drivers = 20L) {
  source_df <- driver_df %>%
    dplyr::mutate(program = clean_program_label(.data$program_class), gene_or_protein = as.character(.data$Gene)) %>%
    dplyr::filter(.data$program %in% publication_program_labels, nzchar(.data$gene_or_protein)) %>%
    dplyr::group_by(.data$program, .data$gene_or_protein) %>%
    dplyr::summarise(
      recurrence_count = dplyr::n_distinct(.data$unit_label),
      significant_term_count = sum(.data$n_terms, na.rm = TRUE),
      spatial_unit_support = paste(sort(unique(.data$unit_label)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$recurrence_count), dplyr::desc(.data$significant_term_count), .data$program, .data$gene_or_protein) %>%
    dplyr::slice_head(n = max_drivers)
  readr::write_csv(source_df, source_file, na = "")
  if (!nrow(source_df)) return(invisible(FALSE))
  plot_df <- source_df %>%
    dplyr::mutate(
      driver_label = factor(.data$gene_or_protein, levels = rev(unique(.data$gene_or_protein))),
      program = factor(.data$program, levels = publication_program_labels)
    )
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$recurrence_count, .data$driver_label, fill = .data$program)) +
    ggplot2::geom_col(width = 0.72, color = "black", linewidth = 0.15) +
    ggplot2::facet_grid(.data$program ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_manual(values = c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#BBBBBB"), guide = "none") +
    ggplot2::labs(x = "Spatial-unit recurrence", y = NULL) +
    theme_nature_base(7) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "grey90", linewidth = 0.2))
  save_nature_svg(p, figure_file, width_mm = 90, height_mm = 95)
}

write_publication_provenance <- function(paths, analyzed_datasets, row_counts) {
  ontology <- go_ontology_provenance()
  registry <- read_manuscript_go_theme_registry(MANUSCRIPT_GO_THEME_REGISTRY)
  writeLines(
    c(
      "compareGO spatial program atlas publication outputs",
      paste0("Datasets analyzed: ", paste(analyzed_datasets, collapse = ", ")),
      paste0("Input contract: ", ATLAS_INPUT_CONTRACT_VERSION),
      paste0("Admitted canonical result family: ", ATLAS_RESULT_TYPE, " / ", ATLAS_ONTOLOGY, " / ", ATLAS_ROUTE_CATEGORY),
      "Term statistics are read from the manifest-declared canonical compareGO term comparison table.",
      "Leading-edge genes and ProteinGroupIDs are derived only from canonical term-gene provenance rows with gene_level_claim_allowed == TRUE and core_enrichment_member == TRUE.",
      paste0("Legacy broad heuristic classes (technical generic outputs only): ", paste(publication_program_labels, collapse = "; ")),
      paste0("SUS-RES manuscript themes: ", paste(unique(registry$display_label[order(registry$display_order)]), collapse = "; ")),
      paste0("Manuscript theme registry: ", MANUSCRIPT_GO_THEME_REGISTRY, " (", unique(registry$registry_version), ")"),
      paste0("GO.db version/source date: ", ontology$go_db_package_version, " / ", ontology$go_source_date),
      paste0("Approved ontology relationships: ", ontology$approved_relationships),
      paste0("Semantic redundancy QA: GOSemSim Wang BP; descriptive cutoff = ", MANUSCRIPT_GO_SEMANTIC_CUTOFF),
      paste0("Descriptive SUS-RES recurrence threshold (not an inferential/display gate): ", atlas_min_recurrent_units),
      paste0("Divergence delta NES threshold: ", divergence_delta_threshold),
      paste0("Supported SUS-RES GO term occurrences: ", row_counts$before),
      paste0("Primary manuscript theme/unit rows: ", row_counts$after),
      "Figures:",
      paste0("  ", paste(paths$figures, collapse = "\n  ")),
      "Source data:",
      paste0("  ", paste(paths$source_data, collapse = "\n  "))
    ),
    file.path(CANONICAL_PATHS$reports, "compareGO_spatial_program_atlas_publication_provenance.txt")
  )
}

assert_no_case_duplicate_columns <- function(df, label) {
  lowered <- tolower(names(df))
  if (anyDuplicated(lowered)) {
    dupes <- unique(lowered[duplicated(lowered)])
    stop(label, " has case-insensitive duplicate column name(s): ", paste(dupes, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# CSV export boundary only: internal stage-07 logic keeps using the legacy
# `Comparison` field (e.g. build_sus_res_supported_go_theme_audit(),
# build_sus_res_theme_spatial_detail()); case-insensitive consumers such as
# PowerShell Import-Csv collide on comparison/Comparison, so the exported
# header renames the legacy field to legacy_Comparison and keeps `comparison`
# canonical.
prepare_enrichment_export <- function(enrichment_df) {
  export_df <- enrichment_df
  if ("Comparison" %in% names(export_df)) {
    names(export_df)[names(export_df) == "Comparison"] <- "legacy_Comparison"
  }
  assert_no_case_duplicate_columns(export_df, "spatial_atlas_enrichment_long.csv")
  export_df
}

xlsx_safe_table <- function(x, limit = 32000L) {
  truncate_cell <- function(value) {
    value <- as.character(value)
    too_long <- !is.na(value) & nchar(value, type = "chars") > limit
    value[too_long] <- paste0(
      substr(value[too_long], 1L, limit - 70L),
      "...[truncated for XLSX cell limit; complete value retained in canonical CSV]"
    )
    value
  }
  x %>% dplyr::mutate(dplyr::across(where(is.character), truncate_cell))
}

write_outputs <- function(enrichment_df, summary_df, behavior_df, driver_df, diagnostics) {
  source_dir <- CANONICAL_PATHS$source_data
  table_dir <- CANONICAL_PATHS$tables
  report_dir <- CANONICAL_PATHS$reports

  readr::write_csv(diagnostics, file.path(report_dir, "input_diagnostics.csv"))
  readr::write_csv(prepare_enrichment_export(enrichment_df), file.path(source_dir, "spatial_atlas_enrichment_long.csv"))
  readr::write_csv(summary_df, file.path(table_dir, "spatial_program_summary.csv"))
  readr::write_csv(behavior_df, file.path(table_dir, "spatial_program_behavior.csv"))
  readr::write_csv(driver_df, file.path(source_dir, "leading_edge_driver_recurrence.csv"))
  writexl::write_xlsx(
    list(
      workbook_contract = tibble::tibble(
        field = c("input_contract_version", "xlsx_character_cell_policy", "complete_provenance_location"),
        value = c(
          ATLAS_INPUT_CONTRACT_VERSION,
          "Character cells longer than 32,000 characters are display-truncated only in this workbook.",
          "Complete untruncated values are retained in the sibling CSV/source-data exports."
        )
      ),
      input_diagnostics = xlsx_safe_table(diagnostics),
      spatial_program_summary = xlsx_safe_table(summary_df),
      spatial_program_behavior = xlsx_safe_table(behavior_df),
      leading_edge_driver_recurrence = xlsx_safe_table(driver_df)
    ),
    file.path(table_dir, "compareGO_spatial_program_atlas_tables.xlsx")
  )
}

build_sus_res_manuscript_theme_outputs <- function(
    enrichment_df,
    registry_path = MANUSCRIPT_GO_THEME_REGISTRY,
    semantic_cutoff = MANUSCRIPT_GO_SEMANTIC_CUTOFF) {
  registry <- read_manuscript_go_theme_registry(registry_path)
  tested_terms <- enrichment_df %>%
    dplyr::filter(
      as.character(.data$phenotype_contrast) == "SUS_vs_RES",
      !is.na(.data$ID), nzchar(as.character(.data$ID))
    ) %>%
    dplyr::distinct(.data$ID, .data$Description)
  mapping <- map_go_terms_to_manuscript_themes(tested_terms, registry)
  supported_audit <- build_sus_res_supported_go_theme_audit(enrichment_df, mapping)
  semantic <- go_semantic_redundancy_qa(supported_audit, cutoff = semantic_cutoff)
  supported_audit <- supported_audit %>%
    dplyr::left_join(
      semantic$term_audit %>% dplyr::select(
        .data$GO_ID, .data$semantic_cluster_id, .data$semantic_cluster_size,
        .data$semantic_cluster_representative_GO_ID,
        .data$semantic_cluster_representative_term
      ),
      by = "GO_ID"
    )
  theme_summary <- build_sus_res_theme_spatial_detail(
    enrichment_df, mapping, semantic$term_audit
  ) %>%
    dplyr::mutate(
      program = .data$manuscript_theme,
      dataset_compartment_label = dplyr::case_when(
        .data$dataset_compartment == "neuron_neuropil" ~ "Neuron neuropil",
        .data$dataset_compartment == "neuron_soma" ~ "Neuron soma",
        .data$dataset_compartment == "microglia" ~ "Microglia-enriched ROI",
        TRUE ~ .data$dataset_compartment
      ),
      spatial_unit_label = clean_spatial_unit_label(.data$spatial_unit),
      interpretable_program = .data$theme_role == "primary"
    )
  list(
    registry = registry,
    mapping = mapping,
    assignment_long = collapse_go_theme_assignment_audit(mapping),
    supported_audit = supported_audit,
    semantic_term_audit = semantic$term_audit,
    semantic_pair_audit = semantic$pair_audit,
    theme_summary = theme_summary,
    assignment_status_summary = summarize_go_theme_assignment_status(supported_audit)
  )
}

build_contract_validation <- function(loaded, enrichment_df, program_summary) {
  bundle_diagnostics <- purrr::map_dfr(loaded, function(x) {
    tibble::tibble(
      dataset = unique(as.character(x$enrichment$dataset))[[1]],
      input_contract_version = ATLAS_INPUT_CONTRACT_VERSION,
      manifest_rows_admitted = nrow(x$manifest),
      manifest_comparisons = dplyr::n_distinct(x$manifest$comparison),
      term_rows_admitted = nrow(x$enrichment),
      term_comparisons = dplyr::n_distinct(x$enrichment$comparison),
      admitted_result_types = collapse_sorted_unique(x$enrichment$result_type),
      admitted_ontologies = collapse_sorted_unique(x$enrichment$ontology),
      admitted_route_categories = collapse_sorted_unique(x$enrichment$route_category),
      claim_safe_core_rows = x$claim_safe_core$n_rows,
      claim_safe_unique_genes = x$claim_safe_core$n_unique_genes,
      claim_safe_unique_protein_groups = x$claim_safe_core$n_unique_protein_groups,
      analysis_status_rows = nrow(x$status)
    )
  })

  admitted_term_counts <- enrichment_df %>%
    dplyr::count(.data$dataset, .data$comparison, .data$result_type, .data$ontology, name = "n_term_rows") %>%
    dplyr::arrange(.data$dataset, .data$comparison, .data$result_type, .data$ontology)

  comparison_identity <- purrr::map_dfr(loaded, function(x) {
    dataset <- unique(as.character(x$enrichment$dataset))[[1]]
    manifest_ids <- sort(unique(as.character(x$manifest$comparison)), method = "radix")
    term_ids <- sort(unique(as.character(x$enrichment$comparison)), method = "radix")
    tibble::tibble(
      dataset = dataset,
      manifest_comparisons = length(manifest_ids),
      admitted_term_comparisons = length(term_ids),
      identities_exact = identical(manifest_ids, term_ids),
      missing_from_terms = paste(setdiff(manifest_ids, term_ids), collapse = ";"),
      extra_in_terms = paste(setdiff(term_ids, manifest_ids), collapse = ";")
    )
  })

  status <- purrr::map_dfr(loaded, "status") %>%
    dplyr::mutate(
      admitted_to_stage07 = .data$result_type == ATLAS_RESULT_TYPE &
        toupper(.data$ontology) == ATLAS_ONTOLOGY &
        .data$comparego_action == "included"
    ) %>%
    dplyr::count(
      .data$dataset, .data$result_type, .data$ontology, .data$analysis_status,
      .data$comparego_action, .data$admitted_to_stage07,
      name = "n_status_rows"
    ) %>%
    dplyr::arrange(.data$dataset, dplyr::desc(.data$admitted_to_stage07), .data$result_type, .data$ontology)

  stage06_detail <- purrr::map_dfr(unique(as.character(program_summary$dataset)), function(dataset) {
    path <- path_results("tables", MODULE_ID, "biological_program_summary", dataset, "program_summary.csv")
    stage07 <- program_summary %>%
      dplyr::filter(.data$dataset == .env$dataset, as.character(.data$program_class) %in% canonical_program_levels) %>%
      dplyr::transmute(
        dataset, comparison, biological_program = as.character(.data$program_class),
        stage07_min_fdr = .data$min_fdr,
        stage07_representative_direction = .data$representative_direction,
        stage07_n_terms = .data$n_terms
      )
    if (!file.exists(path)) {
      return(stage07 %>%
        dplyr::mutate(
          stage06_program_summary_path = path,
          stage06_available = FALSE,
          stage06_min_fdr = NA_real_,
          stage06_representative_direction = NA_character_,
          stage06_n_terms = NA_integer_,
          stage06_key_present = FALSE,
          min_fdr_abs_diff = NA_real_,
          min_fdr_match = NA,
          representative_direction_comparable = FALSE,
          representative_direction_match = NA
        ))
    }
    stage06 <- readr::read_csv(path, show_col_types = FALSE) %>%
      dplyr::transmute(
        dataset = as.character(.data$dataset),
        comparison = as.character(.data$comparison),
        biological_program = as.character(.data$biological_program),
        stage06_min_fdr = suppressWarnings(as.numeric(.data$min_fdr)),
        stage06_representative_direction = as.character(.data$direction),
        stage06_n_terms = suppressWarnings(as.integer(.data$n_terms))
      )
    dplyr::full_join(
      stage07 %>% dplyr::mutate(stage07_key_present = TRUE),
      stage06 %>% dplyr::mutate(stage06_key_present = TRUE),
      by = c("dataset", "comparison", "biological_program")
    ) %>%
      dplyr::mutate(
        stage07_key_present = dplyr::coalesce(.data$stage07_key_present, FALSE),
        stage06_key_present = dplyr::coalesce(.data$stage06_key_present, FALSE),
        stage06_program_summary_path = path,
        stage06_available = TRUE,
        min_fdr_abs_diff = abs(.data$stage07_min_fdr - .data$stage06_min_fdr),
        min_fdr_match = dplyr::if_else(
          .data$stage07_key_present & .data$stage06_key_present,
          dplyr::near(.data$stage07_min_fdr, .data$stage06_min_fdr, tol = 1e-12),
          NA
        ),
        representative_direction_comparable = .data$stage07_key_present & .data$stage06_key_present &
          .data$stage07_representative_direction != "no_FDR_supported_term",
        representative_direction_match = dplyr::if_else(
          .data$representative_direction_comparable,
          .data$stage07_representative_direction == .data$stage06_representative_direction,
          NA
        )
      )
  })

  stage06_summary <- stage06_detail %>%
    dplyr::group_by(.data$dataset) %>%
    dplyr::summarise(
      stage06_available = all(.data$stage06_available),
      stage07_comparisons = dplyr::n_distinct(.data$comparison[.data$stage07_key_present]),
      stage06_comparisons = dplyr::n_distinct(.data$comparison[.data$stage06_key_present]),
      stage07_programs = dplyr::n_distinct(.data$biological_program[.data$stage07_key_present]),
      stage06_programs = dplyr::n_distinct(.data$biological_program[.data$stage06_key_present]),
      overlapping_keys = sum(.data$stage07_key_present & .data$stage06_key_present),
      stage07_only_keys = sum(.data$stage07_key_present & !.data$stage06_key_present),
      stage06_only_keys = sum(!.data$stage07_key_present & .data$stage06_key_present),
      min_fdr_matches = sum(.data$min_fdr_match %in% TRUE),
      min_fdr_mismatches = sum(.data$min_fdr_match %in% FALSE),
      max_min_fdr_abs_diff = suppressWarnings(max(.data$min_fdr_abs_diff, na.rm = TRUE)),
      comparable_representative_directions = sum(.data$representative_direction_comparable),
      representative_direction_matches = sum(.data$representative_direction_match %in% TRUE),
      representative_direction_mismatches = sum(.data$representative_direction_match %in% FALSE),
      direction_scope_note = "Stage 07 exposes no_FDR_supported_term when no term passes FDR<0.05; stage 06 reports the minimum-FDR term direction regardless of significance. Direction is compared only for stage-07 FDR-supported representatives.",
      .groups = "drop"
    ) %>%
    dplyr::mutate(max_min_fdr_abs_diff = ifelse(is.infinite(.data$max_min_fdr_abs_diff), NA_real_, .data$max_min_fdr_abs_diff))

  list(
    bundle_diagnostics = bundle_diagnostics,
    admitted_term_counts = admitted_term_counts,
    comparison_identity = comparison_identity,
    analysis_status = status,
    stage06_detail = stage06_detail,
    stage06_summary = stage06_summary
  )
}

write_contract_validation <- function(validation) {
  readr::write_csv(validation$bundle_diagnostics, file.path(CANONICAL_PATHS$reports, "canonical_bundle_diagnostics.csv"), na = "")
  readr::write_csv(validation$admitted_term_counts, file.path(CANONICAL_PATHS$reports, "admitted_canonical_term_counts.csv"), na = "")
  readr::write_csv(validation$comparison_identity, file.path(CANONICAL_PATHS$reports, "canonical_comparison_identity_audit.csv"), na = "")
  readr::write_csv(validation$analysis_status, file.path(CANONICAL_PATHS$reports, "canonical_analysis_status_admission.csv"), na = "")
  readr::write_csv(validation$stage06_detail, file.path(CANONICAL_PATHS$reports, "stage06_program_summary_crosscheck.csv"), na = "")
  readr::write_csv(validation$stage06_summary, file.path(CANONICAL_PATHS$reports, "stage06_program_summary_crosscheck_summary.csv"), na = "")
}

main <- function() {
  load_required_packages(required_pkgs)
  require_go_ontology_dependencies(include_semantic = TRUE)
  registry_dry_check <- read_manuscript_go_theme_registry(MANUSCRIPT_GO_THEME_REGISTRY)

  datasets <- get_requested_datasets()
  files <- purrr::map(datasets, discover_dataset_files)
  diagnostics <- purrr::map_dfr(files, diagnose_inputs)

  if (is_script_dry_run()) {
    message("[DRY-RUN] compareGO spatial atlas")
    message("[DRY-RUN] input contract: ", ATLAS_INPUT_CONTRACT_VERSION)
    message("[DRY-RUN] admitted result family: ", ATLAS_RESULT_TYPE, " / ", ATLAS_ONTOLOGY, " / ", ATLAS_ROUTE_CATEGORY)
    message("[DRY-RUN] manuscript GO-theme registry: ", MANUSCRIPT_GO_THEME_REGISTRY)
    message("[DRY-RUN] registry version/themes: ", unique(registry_dry_check$registry_version), " / ", length(unique(registry_dry_check$theme_id)))
    message("[DRY-RUN] ontology relationships: ", paste(manuscript_go_allowed_relationships(), collapse = ", "))
    print(diagnostics)
    missing <- diagnostics %>% dplyr::filter(!.data$exists)
    if (nrow(missing) > 0) {
      message("[DRY-RUN] FAIL: missing canonical compareGO contract inputs.")
    } else {
      message("[DRY-RUN] PASS: canonical compareGO manifest, term comparison, term-gene provenance, and analysis status exist for all requested datasets.")
    }
    return(invisible(diagnostics))
  }

  complete_datasets <- diagnostics %>%
    dplyr::group_by(.data$dataset) %>%
    dplyr::summarise(complete = all(.data$exists), .groups = "drop")
  missing_datasets <- complete_datasets$dataset[!complete_datasets$complete]
  if (length(missing_datasets) > 0) {
    stop("Missing canonical compareGO atlas inputs for dataset(s): ", paste(missing_datasets, collapse = ", "), call. = FALSE)
  }

  loaded <- purrr::map(files, load_dataset_atlas)
  enrichment_df <- purrr::map_dfr(loaded, "enrichment")

  program_summary <- calculate_program_summary(enrichment_df)
  behavior_df <- classify_program_behavior(program_summary)
  driver_recurrence <- make_driver_recurrence(enrichment_df)
  manuscript <- build_sus_res_manuscript_theme_outputs(enrichment_df)
  contract_validation <- build_contract_validation(loaded, enrichment_df, program_summary)

  write_outputs(enrichment_df, program_summary, behavior_df, driver_recurrence, diagnostics)
  write_contract_validation(contract_validation)
  readr::write_csv(
    manuscript$assignment_long,
    file.path(CANONICAL_PATHS$source_data, "sus_res_go_term_theme_assignment_long.csv"), na = ""
  )
  readr::write_csv(
    manuscript$mapping$assignments,
    file.path(CANONICAL_PATHS$source_data, "sus_res_go_term_theme_path_audit.csv"), na = ""
  )
  readr::write_csv(
    manuscript$supported_audit,
    file.path(CANONICAL_PATHS$source_data, "sus_res_supported_go_term_theme_audit.csv"), na = ""
  )
  readr::write_csv(
    manuscript$semantic_term_audit,
    file.path(CANONICAL_PATHS$source_data, "sus_res_go_semantic_redundancy_audit.csv"), na = ""
  )
  readr::write_csv(
    manuscript$semantic_pair_audit,
    file.path(CANONICAL_PATHS$source_data, "sus_res_go_semantic_similarity_pairs.csv"), na = ""
  )
  readr::write_csv(
    manuscript$assignment_status_summary,
    file.path(CANONICAL_PATHS$tables, "sus_res_go_theme_assignment_summary.csv"), na = ""
  )
  readr::write_csv(
    manuscript$supported_audit %>% dplyr::select(dplyr::all_of(c(
      "supported_occurrence_id", "dataset", "compartment", "region", "layer",
      "spatial_unit", "GO_ID", "GO_description", "NES", "p.adjust",
      "legacy_program_class", "manuscript_themes", "theme_roles", "assignment_status",
      "mapping_method", "registry_version", "GO_db_package_version", "GO_source_date",
      "relationship_types_approved"
    ))),
    file.path(CANONICAL_PATHS$source_data, "sus_res_legacy_vs_manuscript_classification_audit.csv"), na = ""
  )
  readr::write_csv(
    manuscript$theme_summary,
    file.path(CANONICAL_PATHS$tables, "sus_res_manuscript_theme_summary.csv"), na = ""
  )
  readr::write_csv(
    manuscript$theme_summary,
    file.path(CANONICAL_PATHS$source_data, "sus_res_theme_spatial_detail_all_terms.csv"), na = ""
  )

  plot_spatial_program_atlas(program_summary, file.path(CANONICAL_PATHS$figures, "Fig_SpatialProgramAtlas_dotheatmap.svg"))
  plot_res_sus_divergence(behavior_df, file.path(CANONICAL_PATHS$figures, "Fig_RES_SUS_divergence.svg"))
  plot_driver_recurrence(driver_recurrence, file.path(CANONICAL_PATHS$figures, "Fig_LeadingEdgeDriverRecurrence.svg"))
  plot_compartment_comparison(program_summary, file.path(CANONICAL_PATHS$figures, "Fig_Compartment_Comparison.svg"))

  publication_figures <- c(
    spatial_atlas = file.path(CANONICAL_PATHS$figures, "Fig_SpatialProgramAtlas_dotheatmap_publication.svg"),
    sus_res_spatial_atlas = file.path(CANONICAL_PATHS$figures, "Fig_SpatialProgramAtlas_SUS_vs_RES_publication.svg"),
    neuropil_focus = file.path(CANONICAL_PATHS$figures, "Fig_SpatialProgramAtlas_neuropil_focus_publication.svg"),
    compartment = file.path(CANONICAL_PATHS$figures, "Fig_Compartment_Comparison_publication.svg"),
    divergence = file.path(CANONICAL_PATHS$figures, "Fig_RES_SUS_divergence_publication.svg"),
    topdrivers = file.path(CANONICAL_PATHS$figures, "Fig_LeadingEdgeDriverRecurrence_topdrivers_publication.svg")
  )
  publication_source <- c(
    spatial_atlas = file.path(CANONICAL_PATHS$source_data, "source_data_SpatialProgramAtlas_publication.csv"),
    sus_res_spatial_atlas = file.path(CANONICAL_PATHS$source_data, "source_data_SpatialProgramAtlas_SUS_vs_RES_publication.csv"),
    neuropil_focus = file.path(CANONICAL_PATHS$source_data, "source_data_SpatialProgramAtlas_neuropil_focus_publication.csv"),
    compartment = file.path(CANONICAL_PATHS$source_data, "source_data_Compartment_Comparison_publication.csv"),
    divergence = file.path(CANONICAL_PATHS$source_data, "source_data_RES_SUS_divergence_publication.csv"),
    topdrivers = file.path(CANONICAL_PATHS$source_data, "source_data_LeadingEdgeDriverRecurrence_topdrivers_publication.csv")
  )

  legacy_sus_res_source <- file.path(
    CANONICAL_PATHS$source_data,
    "source_data_SpatialProgramAtlas_SUS_vs_RES_legacy_programs.csv"
  )
  readr::write_csv(
    prepare_sus_res_publication_summary(program_summary, atlas_min_recurrent_units),
    legacy_sus_res_source, na = ""
  )

  publication_source_df <- prepare_publication_summary(program_summary, atlas_min_recurrent_units)
  plot_spatial_program_atlas_publication(program_summary, publication_figures[["spatial_atlas"]], publication_source[["spatial_atlas"]], atlas_min_recurrent_units)
  plot_sus_res_manuscript_theme_atlas(
    manuscript$theme_summary,
    publication_figures[["sus_res_spatial_atlas"]],
    publication_source[["sus_res_spatial_atlas"]]
  )
  plot_neuropil_focus_publication(program_summary, publication_figures[["neuropil_focus"]], publication_source[["neuropil_focus"]])
  plot_compartment_comparison_publication(program_summary, publication_figures[["compartment"]], publication_source[["compartment"]])
  plot_res_sus_divergence_publication(behavior_df, publication_figures[["divergence"]], publication_source[["divergence"]], divergence_delta_threshold)
  plot_driver_recurrence_topdrivers_publication(driver_recurrence, publication_figures[["topdrivers"]], publication_source[["topdrivers"]])

  write_publication_provenance(
    paths = list(figures = unname(publication_figures), source_data = unname(publication_source)),
    analyzed_datasets = vapply(files, `[[`, character(1), "dataset"),
    row_counts = list(
      before = nrow(manuscript$supported_audit),
      after = sum(manuscript$theme_summary$theme_role == "primary", na.rm = TRUE)
    )
  )

  writeLines(
    c(
      "compareGO spatial program atlas complete",
      paste0("Validation-only isolated output: ", VALIDATION_ONLY),
      paste0("Input contract: ", ATLAS_INPUT_CONTRACT_VERSION),
      paste0("Admitted result family: ", ATLAS_RESULT_TYPE, " / ", ATLAS_ONTOLOGY, " / ", ATLAS_ROUTE_CATEGORY),
      paste0("Datasets requested: ", paste(datasets, collapse = ", ")),
      paste0("Datasets analyzed: ", paste(vapply(files, `[[`, character(1), "dataset"), collapse = ", ")),
      paste0("Program summary rows: ", nrow(program_summary)),
      paste0("Supported SUS-RES GO term occurrences: ", nrow(manuscript$supported_audit)),
      paste0("Manuscript theme summary rows: ", nrow(manuscript$theme_summary)),
      paste0("Publication atlas rows: ", sum(publication_source_df$publication_include, na.rm = TRUE)),
      paste0("Figures directory: ", CANONICAL_PATHS$figures)
    ),
    file.path(CANONICAL_PATHS$reports, "compareGO_spatial_program_atlas_summary.txt")
  )

  message("Spatial program atlas complete: ", CANONICAL_PATHS$figures)
  invisible(list(
    summary = program_summary,
    behavior = behavior_df,
    driver_recurrence = driver_recurrence,
    manuscript = manuscript,
    contract_validation = contract_validation
  ))
}

if (sys.nframe() == 0) {
  main()
}
