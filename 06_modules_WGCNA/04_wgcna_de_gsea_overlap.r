#!/usr/bin/env Rscript
#
# File contract:
#   06_modules_WGCNA/04_wgcna_de_gsea_overlap.r
#   results/tables/06_modules_WGCNA/04_wgcna_de_gsea_overlap/<dataset>/

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "module_contracts.R"))

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

tokenize_proteins <- function(x) {
  x <- normalize_module_identifier(x)
  x <- unlist(strsplit(x, "[;/,|[:space:]]+"))
  unique(x[nzchar(x) & !is.na(x)])
}

read_protein_set <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  df <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(character())
  preferred <- intersect(c("UniProt", "ProteinID", "GeneSymbol", "Gene", "gene", "SYMBOL", "symbol"), names(df))
  col <- if (length(preferred)) preferred[[1]] else names(df)[[1]]
  tokenize_proteins(df[[col]])
}

read_leading_edge_set <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  df <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || !"core_enrichment" %in% names(df)) return(character())
  tokenize_proteins(df$core_enrichment)
}

empty_overlap_status <- function(dataset, reason) {
  tibble::tibble(
    dataset = dataset,
    status = "skipped",
    reason = reason,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}

run_wgcna_de_gsea_overlap <- function(dataset = current_dataset(), dry_run = is_dry_run()) {
  dataset <- validate_dataset(dataset)
  out_paths <- create_module_dirs("06_modules_WGCNA", file.path("04_wgcna_de_gsea_overlap", dataset))
  out_csv <- file.path(out_paths$tables, "WGCNA_vs_DE_GSEA_overlap.csv")
  out_xlsx <- file.path(out_paths$tables, "WGCNA_vs_DE_GSEA_overlap.xlsx")
  status_csv <- file.path(out_paths$logs, "WGCNA_vs_DE_GSEA_overlap_status.csv")

  wgcna_modules_dir <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules")
  definitions_file <- file.path(wgcna_modules_dir, "WGCNA_module_definitions_for_downstream.csv")
  universe_file <- file.path(wgcna_modules_dir, "WGCNA_feature_universe.csv")
  priority_file <- file.path(wgcna_modules_dir, "WGCNA_module_priority_summary.csv")
  cluster_manifest <- find_latest(
    path_results("reports", "04_differential_expression_enrichment", "clusterProfiler", dataset),
    "^clusterProfiler_manifest.*\\.csv$"
  )
  compare_manifest <- first_existing_path(c(
    path_processed("04_differential_expression_enrichment", "compareGO", dataset, "compareGO_input_manifest.csv"),
    find_latest(path_processed("04_differential_expression_enrichment", "compareGO", dataset), "^compareGO_input_manifest.*\\.csv$")
  ))

  if (isTRUE(dry_run)) {
    dry_run_line("Script", "06_modules_WGCNA/04_wgcna_de_gsea_overlap.r")
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
  universe_raw <- readr::read_csv(universe_file, show_col_types = FALSE)
  require_module_contract_columns(universe_raw, c("ProteinID", "UniProt", "GeneSymbol", "included_in_wgcna"), "WGCNA feature universe")
  universe <- unique(c(
    tokenize_proteins(universe_raw$ProteinID),
    tokenize_proteins(universe_raw$UniProt),
    tokenize_proteins(universe_raw$GeneSymbol)
  ))
  universe <- universe[nzchar(universe)]
  if (!length(universe)) {
    status <- empty_overlap_status(dataset, paste("WGCNA feature universe was empty:", universe_file))
    readr::write_csv(status, status_csv, na = "")
    readr::write_csv(status, out_csv, na = "")
    write_overlap_xlsx(list(status = status), out_xlsx)
    return(invisible(status))
  }
  for (needed_col in c("UniProt", "ProteinID", "GeneSymbol")) {
    if (!needed_col %in% names(modules_raw)) modules_raw[[needed_col]] <- NA_character_
  }
  modules <- modules_raw %>%
    dplyr::mutate(
      ModuleID = as.character(.data$ModuleID),
      ModuleColor = as.character(.data$ModuleColor),
      protein_token = dplyr::coalesce(.data$UniProt, .data$ProteinID, .data$GeneSymbol)
    ) %>%
    dplyr::filter(!is.na(.data$protein_token), nzchar(.data$protein_token)) %>%
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) %>%
    dplyr::summarise(module_proteins = list(unique(tokenize_proteins(.data$protein_token))), .groups = "drop")

  manifests <- list()
  if (file.exists(cluster_manifest)) {
    manifests$clusterProfiler <- readr::read_csv(cluster_manifest, show_col_types = FALSE) %>%
      dplyr::mutate(route = "clusterProfiler")
  }
  if (file.exists(compare_manifest)) {
    manifests$compareGO <- readr::read_csv(compare_manifest, show_col_types = FALSE) %>%
      dplyr::mutate(route = "compareGO")
  }
  manifest <- dplyr::bind_rows(manifests)
  if ("dataset" %in% names(manifest)) manifest <- dplyr::filter(manifest, .data$dataset == !!dataset)
  if (!"input_gene_file" %in% names(manifest)) manifest$input_gene_file <- NA_character_
  if (!"output_table" %in% names(manifest)) manifest$output_table <- NA_character_
  if (!"contrast" %in% names(manifest)) {
    contrast_candidates <- intersect(c("comparison", "analysis_id", "condition"), names(manifest))
    manifest$contrast <- if (length(contrast_candidates)) as.character(manifest[[contrast_candidates[[1]]]]) else NA_character_
  }
  if (!"route_category" %in% names(manifest)) manifest$route_category <- NA_character_
  if (!nrow(manifest)) {
    status <- empty_overlap_status(dataset, "Manifests were found but no rows matched the requested dataset.")
    readr::write_csv(status, status_csv, na = "")
    readr::write_csv(status, out_csv, na = "")
    write_overlap_xlsx(list(status = status), out_xlsx)
    return(invisible(status))
  }

  rows <- purrr::pmap_dfr(
    list(manifest$contrast, manifest$route, manifest$route_category, manifest$input_gene_file, manifest$output_table),
    function(contrast, route, route_category, input_gene_file, output_table) {
      de_set <- read_protein_set(input_gene_file)
      le_set <- read_leading_edge_set(output_table)
      if (!length(de_set) && !length(le_set)) return(tibble::tibble())
      purrr::pmap_dfr(list(modules$ModuleID, modules$ModuleColor, modules$module_proteins), function(ModuleID, ModuleColor, module_set) {
        module_set <- intersect(module_set, universe)
        de_set <- intersect(de_set, universe)
        le_set <- intersect(le_set, universe)
        signature_set <- union(de_set, le_set)
        de_overlap <- intersect(module_set, de_set)
        le_overlap <- intersect(module_set, le_set)
        fisher_p <- NA_real_
        if (length(universe) && length(de_set)) {
          a <- length(de_overlap)
          b <- length(module_set) - a
          c <- length(de_set) - a
          d <- max(length(universe) - a - b - c, 0)
          fisher_p <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
        }
        tibble::tibble(
          dataset = dataset,
          ModuleID = ModuleID,
          ModuleColor = ModuleColor,
          contrast = as.character(contrast),
          route = as.character(route),
          route_category = as.character(route_category),
          universe_type = "WGCNA_feature_universe",
          n_universe = length(universe),
          n_module_in_universe = length(module_set),
          n_signature_in_universe = length(signature_set),
          n_module_proteins = length(module_set),
          n_DE_proteins = length(de_set),
          n_leading_edge_proteins = length(le_set),
          n_DE_overlap = length(de_overlap),
          n_leading_edge_overlap = length(le_overlap),
          jaccard_DE = if (length(union(module_set, de_set))) length(de_overlap) / length(union(module_set, de_set)) else NA_real_,
          fisher_p = fisher_p,
          top_overlap_proteins = paste(utils::head(unique(c(de_overlap, le_overlap)), 25), collapse = ";")
        )
      })
    }
  )

  if (!nrow(rows)) {
    status <- empty_overlap_status(dataset, "No readable DE/GSEA protein sets were available from manifests.")
    readr::write_csv(status, status_csv, na = "")
    readr::write_csv(status, out_csv, na = "")
    write_overlap_xlsx(list(status = status), out_xlsx)
    return(invisible(status))
  }

  rows <- rows %>%
    dplyr::mutate(fisher_fdr = stats::p.adjust(.data$fisher_p, method = "BH")) %>%
    dplyr::arrange(.data$fisher_fdr, dplyr::desc(.data$n_DE_overlap), dplyr::desc(.data$n_leading_edge_overlap))
  summary <- rows %>%
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) %>%
    dplyr::arrange(.data$fisher_fdr, dplyr::desc(.data$n_DE_overlap), dplyr::desc(.data$n_leading_edge_overlap), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      ModuleID,
      ModuleColor,
      strongest_DE_GSEA_contrast = .data$contrast,
      strongest_DE_GSEA_route = .data$route,
      n_DE_overlap,
      n_leading_edge_overlap,
      jaccard_DE,
      fisher_p,
      fisher_fdr,
      top_overlap_proteins
    )

  readr::write_csv(rows, out_csv, na = "")
  write_overlap_xlsx(list(overlap_long = rows, strongest_overlap_per_module = summary), out_xlsx)
  readr::write_csv(empty_overlap_status(dataset, "completed") %>% dplyr::mutate(status = "completed"), status_csv, na = "")

  if (file.exists(priority_file)) {
    priority <- readr::read_csv(priority_file, show_col_types = FALSE) %>%
      dplyr::select(-dplyr::any_of(names(summary)[!names(summary) %in% c("ModuleID", "ModuleColor")])) %>%
      dplyr::left_join(summary, by = c("ModuleID", "ModuleColor"))
    readr::write_csv(priority, priority_file, na = "")
  }
  invisible(summary)
}

if (sys.nframe() == 0) {
  dataset <- arg_value("--dataset", default = current_dataset())
  if (nzchar(dataset)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset, source = "--dataset"))
  run_wgcna_de_gsea_overlap(current_dataset(), dry_run = "--dry-run" %in% args || is_dry_run())
}
