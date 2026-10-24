#!/usr/bin/env Rscript
# ================================================================
# Script: analysis/wgcna/compare_module_enrichment_overlap.R
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/; data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv; +1 more; optional results/tables/04_differential_expression_enrichment/biological_program_summary/<dataset>/program_summary.csv.
# Produces: results/tables/06_modules_WGCNA/04_wgcna_de_gsea_overlap/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Must run after WGCNA and enrichment.
# ================================================================

#
# File contract:
#   analysis/wgcna/compare_module_enrichment_overlap.R
#   results/tables/06_modules_WGCNA/04_wgcna_de_gsea_overlap/<dataset>/

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "module_contracts.R"))
source(repo_path("R", "wgcna_paths.R"))
source(repo_path("R", "gsea_wgcna_concordance_utils.R"))
source(repo_path("R", "enrichment_io.R"))  # read_canonical_clusterprofiler_manifest(): runtime paths, long-path staging

packages <- c("dplyr", "tidyr", "readr", "stringr", "tibble", "purrr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}

find_latest <- function(root, pattern) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(root)) return(NA_character_)
  candidates <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(NA_character_)
  info <- file.info(candidates)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

write_overlap_xlsx <- function(sheets, path) {
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path)
  }
  invisible(path)
}

# Both sides are compared as dataset-scoped ProteinGroupIDs over the WGCNA universe.
de_set_rule <- "DA table significant (padj < 0.05)"
leading_edge_set_rule <- "core_enrichment_member of p.adjust < 0.05 terms, via term_gene_provenance"

# Enrichment manifests record ABSOLUTE paths from the machine that produced them
# (a drive letter followed by the processed-data path). On any other machine
# that drive does not exist, every referenced set silently reads as empty and
# the overlap table is
# written as "skipped" - which is how it came to be frozen at a superseded run.
# Fall back to resolving the recorded path against the current repository root.
resolve_manifest_path <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return(NA_character_)
  if (file.exists(path)) return(path)
  rel <- sub("^[A-Za-z]:[/\\\\]+", "", gsub("\\\\", "/", path))
  rel <- sub("^/+", "", rel)
  candidate <- repo_path(rel)
  if (file.exists(candidate)) return(candidate)
  NA_character_
}

# A declared file that cannot be read, or lacks the contract columns, fails loudly
# instead of reading as an empty set.
read_contract_columns <- function(path, columns, optional = character()) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return(NULL)
  resolved <- resolve_manifest_path(path)
  if (is.na(resolved)) stop("Declared enrichment input is not readable: ", path, call. = FALSE)
  df <- readr::read_csv(
    resolved, col_select = dplyr::any_of(c(columns, optional)),
    show_col_types = FALSE, progress = FALSE
  )
  missing <- setdiff(columns, names(df))
  if (length(missing)) {
    stop("Declared enrichment input ", resolved, " is missing column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  df
}

# DE set: ProteinGroupIDs the per-contrast DA table flags significant.
read_de_set <- function(path) {
  df <- read_contract_columns(path, c("ProteinGroupID", "significant"))
  if (is.null(df)) return(character())
  unique(as.character(df$ProteinGroupID[as.logical(df$significant) %in% TRUE]))
}

# Leading edge: core members of this comparison's p.adjust < 0.05 terms. The producer's
# term-gene provenance bridges official symbols (KEGG: Entrez) to ProteinGroupIDs; same
# rule as build_go_program_atlas.R.
read_leading_edge_set <- function(term_path, provenance_path, comparison) {
  terms <- read_contract_columns(term_path, c("ID", "p.adjust"), optional = "comparison")
  provenance <- read_contract_columns(provenance_path, c(
    "comparison", "term_id", "ProteinGroupID", "core_enrichment_member", "gene_level_claim_allowed"
  ))
  if (is.null(terms) || is.null(provenance)) return(character())
  # The compareGO term table pools all comparisons.
  if ("comparison" %in% names(terms)) terms <- terms[terms$comparison %in% comparison, , drop = FALSE]
  significant <- terms$ID[is.finite(terms$p.adjust) & terms$p.adjust < 0.05]
  keep <- provenance$comparison %in% comparison & provenance$term_id %in% significant &
    provenance$core_enrichment_member %in% TRUE & provenance$gene_level_claim_allowed %in% TRUE
  unique(as.character(provenance$ProteinGroupID[keep]))
}

# Exact, namespace-checked identity through the GSEA-WGCNA concordance token map; a
# ProteinGroupID from another dataset fails closed there.
map_to_wgcna_universe <- function(ids, bundle, dataset) {
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (!length(ids)) return(character())
  if (all(gww_is_protein_group_id(ids, dataset)) && !any(ids %in% bundle$universe$ProteinGroupID)) {
    return(character())
  }
  evidence <- data.frame(
    dataset = dataset, gsea_evidence_id = "set",
    leading_edge_accessions = paste(ids, collapse = ";"),
    leading_edge_token_type = "protein_group_id", stringsAsFactors = FALSE
  )
  gww_split_proteins(gww_map_leading_edges(evidence, bundle$token_map)$leading_edge_proteins)
}

empty_overlap_status <- function(dataset, reason) {
  tibble::tibble(
    dataset = dataset,
    status = "skipped",
    reason = reason,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}

overlap_stats <- function(module_set, signature_set, universe, labels = NULL) {
  module_set <- intersect(unique(module_set), universe)
  signature_set <- intersect(unique(signature_set), universe)
  overlap <- intersect(module_set, signature_set)
  p <- NA_real_
  if (length(universe) && length(module_set) && length(signature_set)) {
    a <- length(overlap)
    b <- length(module_set) - a
    c <- length(signature_set) - a
    d <- max(length(universe) - a - b - c, 0)
    p <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
  }
  top <- utils::head(overlap, 25)
  list(
    n_signature = length(signature_set),
    n_overlap = length(overlap),
    jaccard = if (length(union(module_set, signature_set))) length(overlap) / length(union(module_set, signature_set)) else NA_real_,
    fisher_p = p,
    overlap_proteins = paste(top, collapse = ";"),
    overlap_genes = paste(if (is.null(labels)) top else dplyr::coalesce(unname(labels[top]), top), collapse = ";")
  )
}

run_wgcna_de_gsea_overlap <- function(dataset = current_dataset(), dry_run = is_dry_run()) {
  dataset <- validate_dataset(dataset)
  out_paths <- wgcna_dirs("compare_module_enrichment_overlap", dataset, create = TRUE)
  out_csv <- file.path(out_paths$tables, "WGCNA_vs_DE_GSEA_overlap.csv")
  out_xlsx <- file.path(out_paths$tables, "WGCNA_vs_DE_GSEA_overlap.xlsx")
  status_csv <- file.path(out_paths$logs, "WGCNA_vs_DE_GSEA_overlap_status.csv")

  wgcna_modules_dir <- wgcna_modules_dir(dataset)
  definitions_file <- wgcna_modules_artifact("WGCNA_module_definitions_for_downstream.csv", dataset, child = "tables", "modules")
  universe_file <- wgcna_modules_artifact("WGCNA_feature_universe.csv", dataset, child = "tables", "modules")
  priority_file <- wgcna_modules_artifact("WGCNA_module_priority_summary.csv", dataset, child = "tables", "modules")
  # The current enrichment producers' manifests; the historical 04_* copies are older runs.
  cluster_manifest <- path_results("differential_abundance", "run_clusterprofiler_enrichment", dataset, "models", "clusterProfiler_manifest.csv")
  compare_manifest <- path_results("differential_abundance", "compare_go_enrichment", dataset, "models", "compareGO_input_manifest.csv")

  if (isTRUE(dry_run)) {
    dry_run_line("Script", "analysis/wgcna/compare_module_enrichment_overlap.R")
    dry_run_line("Dataset", dataset)
    dry_run_line("WGCNA downstream definitions", definitions_file, if (file.exists(definitions_file)) "PASS" else "WARN")
    dry_run_line("WGCNA Fisher background universe", universe_file, if (file.exists(universe_file)) "PASS" else "FAIL")
    dry_run_line("clusterProfiler manifest", cluster_manifest, if (file.exists(cluster_manifest)) "PASS" else "WARN")
    dry_run_line("compareGO input manifest", compare_manifest, if (file.exists(compare_manifest)) "PASS" else "WARN")
    dry_run_line("Output overlap table", out_csv)
    return(invisible(file.exists(definitions_file) && file.exists(universe_file) && (file.exists(cluster_manifest) || file.exists(compare_manifest))))
  }

  if (!file.exists(definitions_file)) {
    status <- empty_overlap_status(dataset, paste("Missing WGCNA downstream definitions:", definitions_file))
    readr::write_csv(status, status_csv, na = "")
    return(invisible(status))
  }
  if (!file.exists(universe_file)) {
    status <- empty_overlap_status(dataset, paste("Missing WGCNA Fisher-test universe:", universe_file))
    readr::write_csv(status, status_csv, na = "")
    readr::write_csv(status, out_csv, na = "")
    write_overlap_xlsx(list(status = status), out_xlsx)
    return(invisible(status))
  }
  if (!file.exists(cluster_manifest) && !file.exists(compare_manifest)) {
    status <- empty_overlap_status(dataset, "No clusterProfiler_manifest.csv or compareGO_input_manifest.csv found for dataset.")
    readr::write_csv(status, status_csv, na = "")
    readr::write_csv(status, out_csv, na = "")
    write_overlap_xlsx(list(status = status), out_xlsx)
    return(invisible(status))
  }

  modules_raw <- readr::read_csv(definitions_file, show_col_types = FALSE)
  validate_wgcna_module_definitions(modules_raw, "WGCNA downstream definitions")
  universe_raw <- readr::read_csv(
    universe_file,
    col_select = dplyr::any_of(c(
      "ProteinGroupID", "included_in_wgcna", "ModuleID", "member_accessions", "RepresentativeUniProt",
      "MemberUniProts", "representative_accession", "official_gene_symbol"
    )),
    show_col_types = FALSE
  )
  # The ProteinGroupID universe and token map the GSEA-WGCNA concordance uses.
  bundle <- gww_build_universe_bundle(universe_raw, dataset)
  universe <- bundle$universe$ProteinGroupID
  if (!length(universe)) {
    status <- empty_overlap_status(dataset, paste("WGCNA feature universe was empty:", universe_file))
    readr::write_csv(status, status_csv, na = "")
    readr::write_csv(status, out_csv, na = "")
    write_overlap_xlsx(list(status = status), out_xlsx)
    return(invisible(status))
  }
  gene_labels <- if ("official_gene_symbol" %in% names(universe_raw)) {
    stats::setNames(as.character(universe_raw$official_gene_symbol), as.character(universe_raw$ProteinGroupID))
  }
  modules <- modules_raw %>%
    dplyr::mutate(
      ModuleID = as.character(.data$ModuleID),
      ModuleColor = as.character(.data$ModuleColor),
      ProteinGroupID = as.character(.data$ProteinGroupID)
    ) %>%
    dplyr::filter(!is.na(.data$ProteinGroupID), nzchar(.data$ProteinGroupID)) %>%
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) %>%
    dplyr::summarise(module_proteins = list(unique(.data$ProteinGroupID)), .groups = "drop")

  manifests <- list()
  if (file.exists(cluster_manifest)) {
    manifests$clusterProfiler <- tibble::as_tibble(read_canonical_clusterprofiler_manifest(cluster_manifest, dataset)) %>%
      dplyr::mutate(route = "clusterProfiler")
  }
  if (file.exists(compare_manifest)) {
    manifests$compareGO <- readr::read_csv(compare_manifest, show_col_types = FALSE) %>%
      dplyr::mutate(route = "compareGO")
  }
  # compareGO re-collates the clusterProfiler GSEA_GO runs under the same analysis_id; read it
  # only for runs clusterProfiler lacks, so no test enters the BH family twice.
  if (length(manifests) == 2L && all(vapply(manifests, function(x) "analysis_id" %in% names(x), logical(1)))) {
    manifests$compareGO <- dplyr::filter(manifests$compareGO, !.data$analysis_id %in% manifests$clusterProfiler$analysis_id)
  }
  manifest <- dplyr::bind_rows(manifests)
  if ("dataset" %in% names(manifest)) manifest <- dplyr::filter(manifest, .data$dataset == !!dataset)
  # A failed enrichment run declares outputs it never wrote.
  if ("analysis_status" %in% names(manifest)) manifest <- dplyr::filter(manifest, !.data$analysis_status %in% "failed")
  if (!"input_gene_file" %in% names(manifest)) manifest$input_gene_file <- NA_character_
  if (!"output_table" %in% names(manifest)) manifest$output_table <- NA_character_
  if (!"contrast" %in% names(manifest)) {
    contrast_candidates <- intersect(c("comparison", "analysis_id", "condition"), names(manifest))
    manifest$contrast <- if (length(contrast_candidates)) as.character(manifest[[contrast_candidates[[1]]]]) else NA_character_
  }
  if (!"route_category" %in% names(manifest)) manifest$route_category <- NA_character_
  if (!"comparison" %in% names(manifest)) manifest$comparison <- manifest$contrast
  for (column in c("result_type", "term_gene_provenance_file", "term_gene_provenance_output_file")) {
    if (!column %in% names(manifest)) manifest[[column]] <- NA_character_
  }
  # A compareGO row pairs its pooled term table with the pooled provenance it declares.
  manifest$provenance_file <- dplyr::coalesce(
    as.character(manifest$term_gene_provenance_output_file),
    as.character(manifest$term_gene_provenance_file)
  )
  if (!nrow(manifest)) {
    status <- empty_overlap_status(dataset, "Manifests were found but no rows matched the requested dataset.")
    readr::write_csv(status, status_csv, na = "")
    readr::write_csv(status, out_csv, na = "")
    write_overlap_xlsx(list(status = status), out_xlsx)
    return(invisible(status))
  }

  rows <- purrr::pmap_dfr(
    list(manifest$contrast, manifest$comparison, manifest$route, manifest$route_category, manifest$result_type,
         manifest$input_gene_file, manifest$output_table, manifest$provenance_file),
    function(contrast, comparison, route, route_category, result_type, input_gene_file, output_table, provenance_file) {
      de_set <- map_to_wgcna_universe(read_de_set(input_gene_file), bundle, dataset)
      le_set <- map_to_wgcna_universe(read_leading_edge_set(output_table, provenance_file, comparison), bundle, dataset)
      if (!length(de_set) && !length(le_set)) return(tibble::tibble())
      purrr::pmap_dfr(list(modules$ModuleID, modules$ModuleColor, modules$module_proteins), function(ModuleID, ModuleColor, module_set) {
        module_set <- intersect(module_set, universe)
        de_set <- intersect(de_set, universe)
        le_set <- intersect(le_set, universe)
        union_set <- union(de_set, le_set)
        de <- overlap_stats(module_set, de_set, universe)
        le <- overlap_stats(module_set, le_set, universe)
        both <- overlap_stats(module_set, union_set, universe, labels = gene_labels)
        tibble::tibble(
          dataset = dataset,
          ModuleID = ModuleID,
          ModuleColor = ModuleColor,
          contrast = as.character(contrast),
          route = as.character(route),
          route_category = as.character(route_category),
          result_type = as.character(result_type),
          universe_type = "WGCNA_feature_universe",
          DE_set_rule = de_set_rule,
          leading_edge_set_rule = leading_edge_set_rule,
          n_universe = length(universe),
          n_module_in_universe = length(module_set),
          n_signature_in_universe = length(union_set),
          n_module_proteins = length(module_set),
          input_status = dplyr::case_when(
            length(de_set) > 0L & length(le_set) > 0L ~ "DE_and_leading_edge_available",
            length(de_set) > 0L ~ "DE_only_available",
            length(le_set) > 0L ~ "leading_edge_only_available",
            TRUE ~ "no_readable_input"
          ),
          n_DE_proteins = de$n_signature,
          n_leading_edge_proteins = le$n_signature,
          n_union_DE_leading_edge_proteins = both$n_signature,
          n_DE_overlap = de$n_overlap,
          n_leading_edge_overlap = le$n_overlap,
          n_union_DE_leading_edge_overlap = both$n_overlap,
          jaccard_DE = de$jaccard,
          jaccard_leading_edge = le$jaccard,
          jaccard_union_DE_leading_edge = both$jaccard,
          fisher_p_DE = de$fisher_p,
          fisher_p_leading_edge = le$fisher_p,
          fisher_p_union_DE_leading_edge = both$fisher_p,
          fisher_p = both$fisher_p,
          top_DE_overlap_proteins = de$overlap_proteins,
          top_leading_edge_overlap_proteins = le$overlap_proteins,
          top_union_overlap_proteins = both$overlap_proteins,
          top_overlap_proteins = both$overlap_proteins,
          top_overlap_genes = both$overlap_genes
        )
      })
    }
  )

  if (!nrow(rows)) {
    status <- empty_overlap_status(dataset, "No manifest entry had a significant DE protein or a significant-term leading edge in the WGCNA universe.")
    readr::write_csv(status, status_csv, na = "")
    readr::write_csv(status, out_csv, na = "")
    write_overlap_xlsx(list(status = status), out_xlsx)
    return(invisible(status))
  }

  rows <- rows %>%
    dplyr::mutate(
      fisher_fdr_DE = stats::p.adjust(.data$fisher_p_DE, method = "BH"),
      fisher_fdr_leading_edge = stats::p.adjust(.data$fisher_p_leading_edge, method = "BH"),
      fisher_fdr_union_DE_leading_edge = stats::p.adjust(.data$fisher_p_union_DE_leading_edge, method = "BH"),
      fisher_fdr = .data$fisher_fdr_union_DE_leading_edge,
      best_overlap_fdr = pmin(.data$fisher_fdr_DE, .data$fisher_fdr_leading_edge, .data$fisher_fdr_union_DE_leading_edge, na.rm = TRUE),
      best_overlap_fdr = ifelse(is.infinite(.data$best_overlap_fdr), NA_real_, .data$best_overlap_fdr)
    ) %>%
    dplyr::arrange(.data$best_overlap_fdr, dplyr::desc(.data$n_union_DE_leading_edge_overlap), dplyr::desc(.data$n_DE_overlap), dplyr::desc(.data$n_leading_edge_overlap))
  summary <- rows %>%
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) %>%
    dplyr::arrange(.data$best_overlap_fdr, dplyr::desc(.data$n_union_DE_leading_edge_overlap), dplyr::desc(.data$n_DE_overlap), dplyr::desc(.data$n_leading_edge_overlap), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      ModuleID,
      ModuleColor,
      strongest_DE_GSEA_contrast = .data$contrast,
      strongest_DE_GSEA_route = .data$route,
      strongest_DE_GSEA_result_type = .data$result_type,
      input_status,
      n_DE_overlap,
      n_leading_edge_overlap,
      n_union_DE_leading_edge_overlap,
      jaccard_DE,
      jaccard_leading_edge,
      jaccard_union_DE_leading_edge,
      fisher_p_DE,
      fisher_fdr_DE,
      fisher_p_leading_edge,
      fisher_fdr_leading_edge,
      fisher_p_union_DE_leading_edge,
      fisher_fdr_union_DE_leading_edge,
      fisher_p,
      fisher_fdr,
      best_overlap_fdr,
      top_overlap_proteins,
      top_overlap_genes
    )

  readr::write_csv(rows, out_csv, na = "")
  write_overlap_xlsx(list(overlap_long = rows, strongest_overlap_per_module = summary), out_xlsx)
  availability <- rows %>%
    dplyr::distinct(.data$contrast, .data$route, .data$route_category, .data$result_type, .data$input_status) %>%
    dplyr::count(.data$input_status, name = "n_manifest_entries")
  n_without_sets <- nrow(manifest) - sum(availability$n_manifest_entries)
  status <- empty_overlap_status(dataset, "completed") %>%
    dplyr::mutate(
      status = "completed",
      input_availability = paste(c(
        paste0(availability$input_status, "=", availability$n_manifest_entries),
        paste0("no_significant_set=", n_without_sets)
      ), collapse = ";"),
      test_semantics = paste0(
        "Separate Fisher enrichment tests for DE, GSEA leading edge, and union(DE, leading edge); alternative='greater'. ",
        "ProteinGroupID identity over the WGCNA universe; DE = ", de_set_rule, "; leading edge = ", leading_edge_set_rule, "."
      )
    )
  readr::write_csv(status, status_csv, na = "")

  if (file.exists(priority_file)) {
    priority <- readr::read_csv(priority_file, show_col_types = FALSE) %>%
      dplyr::select(-dplyr::any_of(names(summary)[!names(summary) %in% c("ModuleID", "ModuleLegacyID", "ModuleColor", "ModuleColorName", "ModuleColorLabel")])) %>%
      dplyr::left_join(
        summary %>% dplyr::select(-dplyr::any_of(c("ModuleLegacyID", "ModuleColor", "ModuleColorName", "ModuleColorLabel"))),
        by = "ModuleID"
      )
    readr::write_csv(priority, priority_file, na = "")
  }
  invisible(summary)
}

if (sys.nframe() == 0) {
  dataset <- arg_value("--dataset", default = current_dataset())
  if (nzchar(dataset)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset, source = "--dataset"))
  run_wgcna_de_gsea_overlap(current_dataset(), dry_run = "--dry-run" %in% args || is_dry_run())
}
