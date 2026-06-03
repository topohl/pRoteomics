#!/usr/bin/env Rscript
#
# Import cached/package-accessible reference marker evidence into a normalized WGCNA marker registry.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
dry_run <- is_dry_run()
SUBSTEP_ID <- "reference_marker_import"
PATHS <- create_module_dirs("03_qc_exploration", file.path("reference_marker_import"))
manifest_file <- Sys.getenv(
  "PROTEOMICS_REFERENCE_MARKER_SOURCES_YML",
  unset = path_external("reference_markers", "reference_marker_sources.yml")
)
allow_download <- tolower(Sys.getenv("PROTEOMICS_REFERENCE_MARKERS_ALLOW_DOWNLOAD", unset = "false")) %in% c("1", "true", "yes", "y")
top_n <- suppressWarnings(as.integer(Sys.getenv("PROTEOMICS_REFERENCE_MARKER_TOP_N_PER_CLASS", unset = "50")))
if (!is.finite(top_n) || top_n < 1L) top_n <- 50L
min_specificity <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_REFERENCE_MARKER_MIN_SPECIFICITY", unset = NA_character_)))

if (dry_run) {
  invisible(lapply(unlist(PATHS), dir_create))
  dry_run_line("Script", "03_qc_exploration/04b_import_reference_marker_sources.r")
  dry_run_line("Manifest", manifest_file, if (file.exists(manifest_file)) "PASS" else "FAIL")
  dry_run_line("Live download allowed", allow_download, "INFO")
  dry_run_line("Output registry", repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"), "INFO")
  quit(status = if (file.exists(manifest_file)) 0L else 1L, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "svglite", "readr", "yaml")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(manifest_file)) stop("Missing marker source manifest: ", manifest_file, call. = FALSE)
manifest <- yaml::read_yaml(manifest_file)
sources <- manifest$sources %||% list()

canonical_marker_panels <- list(
  canonical_microglia_homeostatic = c("P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Hexb", "Fcrls", "Sall1", "Siglech", "Gpr34", "Mertk", "Aif1"),
  canonical_microglia_phagolysosomal_state = c("Tyrobp", "Trem2", "Apoe", "Lpl", "Cst7", "Ctsb", "Ctsd", "Lgals3", "Itgax", "Axl", "C1qa", "C1qb", "C1qc"),
  canonical_neuronal_synaptic_neuropil = c("Snap25", "Syp", "Syn1", "Vamp2", "Stx1a", "Stxbp1", "Dlg4", "Camk2a", "Camk2b", "Map2", "Nefl", "Nefm", "Rbfox3", "Grin1", "Gria1"),
  canonical_neuronal_soma_nuclear = c("Rbfox3", "Map2", "Tubb3", "H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3", "Matr3", "Srsf3", "Ddx39b"),
  canonical_astrocyte = c("Aqp4", "Gfap", "Aldh1l1", "Slc1a2", "Slc1a3", "Glul", "Aldoc", "Gja1", "S100b"),
  canonical_oligodendrocyte_myelin = c("Mbp", "Plp1", "Mog", "Cnp", "Mag", "Mobp", "Cldn11", "Myrf", "Olig1", "Olig2"),
  canonical_opc = c("Pdgfra", "Cspg4", "Vcan", "Sox10", "Olig1", "Olig2"),
  canonical_endothelial_vascular = c("Cldn5", "Pecam1", "Kdr", "Flt1", "Slco1a4"),
  canonical_pericyte_vascular = c("Rgs5", "Pdgfrb", "Vtn", "Acta2"),
  canonical_peripheral_myeloid_caution = c("Lyz2", "Cd74", "H2-Ab1", "Fcgr1", "Ccr2", "Ly6c2", "Itgam"),
  canonical_mitochondrial_oxphos = c("Ndufs1", "Ndufa9", "Sdha", "Uqcrc2", "Cox4i1", "Atp5f1a", "Atp5f1b"),
  canonical_ribosomal_translation = c("Rpl3", "Rpl4", "Rpl5", "Rps3", "Rps6", "Eef1a1", "Eef2"),
  canonical_rnp_rna_processing = c("Hnrnpa2b1", "Hnrnpc", "Sfpq", "Snrnp70", "Ddx5", "Ddx17", "Pabpc1")
)

panel_meta <- data.frame(
  marker_set = names(canonical_marker_panels),
  cell_class = c("microglia", "microglia", "neuron", "neuron", "astrocyte", "oligodendrocyte", "opc", "endothelial", "pericyte", "peripheral_myeloid", "mitochondrial", "ribosomal", "rnp"),
  cell_state = c("homeostatic_identity", "phagolysosomal_complement_state", "synaptic_neuropil", "soma_nuclear", "canonical", "myelin", "canonical", "vascular", "vascular", "caution", "oxphos", "translation", "rna_processing"),
  stringsAsFactors = FALSE
)

read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(tryCatch(readRDS(path), error = function(e) NULL))
  if (ext == "xlsx" || ext == "xls") {
    if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
    return(tryCatch(as.data.frame(readxl::read_excel(path), check.names = FALSE), error = function(e) NULL))
  }
  if (ext == "tsv") return(tryCatch(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL))
  if (ext == "csv") return(tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL))
  NULL
}

coerce_marker_table <- function(df, source) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  gene_col <- first_present_col(df, c("gene_symbol", "GeneSymbol", "gene", "Gene", "symbol", "marker", "Marker"))
  class_col <- first_present_col(df, c("cell_class", "cell_type", "celltype", "CellType", "class", "cluster", "reference_cell_label"))
  score_col <- first_present_col(df, c("specificity_score", "specificity", "auc", "score", "logFC", "avg_log2FC", "rank"))
  if (is.na(gene_col) || is.na(class_col)) return(NULL)
  out <- data.frame(
    source_name = source$source_name %||% "unknown",
    source_version = source$source_version %||% "cached_local",
    source_type = source$source_type %||% "local_file",
    gene_symbol = as.character(df[[gene_col]]),
    cell_class = as.character(df[[class_col]]),
    cell_state = NA_character_,
    reference_cell_label = as.character(df[[class_col]]),
    mean_target = NA_real_,
    mean_other = NA_real_,
    logFC_target_vs_other = if (!is.na(score_col) && grepl("log", score_col, ignore.case = TRUE)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_,
    specificity_score = if (!is.na(score_col)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_,
    pct_target = NA_real_,
    pct_other = NA_real_,
    rank_within_source = NA_integer_,
    selection_rule = "cached_marker_table",
    selected = TRUE,
    confidence = "external_cached",
    notes = source$notes %||% "",
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(normalize_gene_token(out$gene_symbol)), , drop = FALSE]
  out <- out |>
    dplyr::group_by(.data$source_name, .data$cell_class) |>
    dplyr::arrange(dplyr::desc(dplyr::coalesce(.data$specificity_score, abs(.data$logFC_target_vs_other), 0)), .by_group = TRUE) |>
    dplyr::mutate(rank_within_source = dplyr::row_number()) |>
    dplyr::ungroup()
  if (is.finite(min_specificity)) out <- out |> dplyr::filter(is.na(.data$specificity_score) | .data$specificity_score >= min_specificity)
  out |> dplyr::filter(.data$rank_within_source <= top_n)
}

source_status <- list()
external_evidence <- list()
for (src in sources) {
  src_name <- src$source_name %||% "unknown"
  src_dir <- normalizePath(repo_path(src$local_dir %||% ""), winslash = "/", mustWork = FALSE)
  status <- "missing_optional"
  n_rows <- 0L
  if (identical(src$source_type, "package") && requireNamespace(src$package %||% "", quietly = TRUE)) {
    status <- "package_available_not_auto_extracted"
  }
  patterns <- unlist(src$file_patterns %||% c("*.csv", "*.tsv", "*.xlsx", "*.rds"), use.names = FALSE)
  files <- character()
  if (dir.exists(src_dir)) {
    regex <- paste0("(", paste(gsub("\\*", ".*", gsub("\\.", "\\\\.", patterns)), collapse = "|"), ")$")
    files <- list.files(src_dir, pattern = regex, full.names = TRUE, ignore.case = TRUE)
  }
  if (length(files)) {
    parsed <- lapply(files, function(f) coerce_marker_table(read_any_table(f), src))
    parsed <- parsed[!vapply(parsed, is.null, logical(1))]
    if (length(parsed)) {
      tab <- dplyr::bind_rows(parsed)
      external_evidence[[src_name]] <- tab
      n_rows <- nrow(tab)
      status <- "imported_cached_file"
    }
  }
  if (allow_download && identical(status, "missing_optional")) status <- "download_allowed_but_not_implemented"
  source_status[[src_name]] <- data.frame(source_name = src_name, status = status, n_rows = n_rows, local_dir = src_dir, stringsAsFactors = FALSE)
}

fallback_registry <- dplyr::bind_rows(lapply(names(canonical_marker_panels), function(ms) {
  meta <- panel_meta[panel_meta$marker_set == ms, , drop = FALSE]
  data.frame(
    marker_set = ms,
    cell_class = meta$cell_class,
    cell_state = meta$cell_state,
    gene_symbol = canonical_marker_panels[[ms]],
    source_type = "curated_fallback",
    source_name = "proteomics_curated_fallback",
    source_reference = "03_qc_exploration/04b_import_reference_marker_sources.r curated conservative fallback",
    selection_rule = "conservative_curated_fallback; annotation_only",
    confidence = "curated_conservative",
    use_for = "annotation_reporting_sensitivity_interpretation",
    notes = "Marker evidence only; not purity correction, not WGCNA filtering, not CON/RES/SUS-derived.",
    stringsAsFactors = FALSE
  )
}))

evidence <- dplyr::bind_rows(external_evidence)
if (!nrow(evidence)) {
  evidence <- fallback_registry |>
    dplyr::transmute(
      source_name, source_version = as.character(Sys.Date()), source_type, gene_symbol, cell_class, cell_state,
      reference_cell_label = cell_class, mean_target = NA_real_, mean_other = NA_real_,
      logFC_target_vs_other = NA_real_, specificity_score = NA_real_, pct_target = NA_real_, pct_other = NA_real_,
      rank_within_source = dplyr::row_number(), selection_rule, selected = TRUE, confidence, notes
    )
} else {
  evidence <- evidence |>
    dplyr::mutate(selected = .data$rank_within_source <= top_n)
}

registry_external <- if (nrow(dplyr::bind_rows(external_evidence))) {
  dplyr::bind_rows(external_evidence) |>
    dplyr::filter(.data$selected %in% TRUE) |>
    dplyr::transmute(
      marker_set = paste0("reference_", safe_filename(tolower(.data$cell_class))),
      cell_class, cell_state, gene_symbol, source_type, source_name,
      source_reference = source_name, selection_rule, confidence, use_for = "annotation_reporting_sensitivity_interpretation", notes
    )
} else {
  data.frame()
}
registry <- dplyr::bind_rows(fallback_registry, registry_external) |>
  dplyr::mutate(gene_token = normalize_gene_token(.data$gene_symbol)) |>
  dplyr::filter(nzchar(.data$gene_token)) |>
  dplyr::distinct(.data$marker_set, .data$gene_token, .keep_all = TRUE) |>
  dplyr::select(-"gene_token")

ranked <- registry |>
  dplyr::group_by(.data$marker_set) |>
  dplyr::mutate(rank = dplyr::row_number(), n_markers = dplyr::n()) |>
  dplyr::ungroup()

write_table_and_source(evidence, PATHS$tables, PATHS$source_data, "reference_marker_evidence_long.csv")
write_table_and_source(ranked, PATHS$tables, PATHS$source_data, "reference_marker_sets_ranked.csv")
write_csv_safe2(registry, repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"))

status_tbl <- dplyr::bind_rows(source_status)
p_status <- ggplot2::ggplot(status_tbl, ggplot2::aes(x = .data$source_name, y = .data$n_rows, fill = .data$status)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Imported evidence rows", fill = "Status") +
  ggplot2::theme_classic(base_size = 8)
ggplot2::ggsave(file.path(PATHS$figures, "reference_marker_source_summary.svg"), p_status, width = 140, height = 80, units = "mm", device = svglite::svglite)

counts <- registry |> dplyr::count(.data$cell_class, name = "n_markers")
p_counts <- ggplot2::ggplot(counts, ggplot2::aes(x = .data$cell_class, y = .data$n_markers)) +
  ggplot2::geom_col(fill = "#2F6F73") +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Registry markers") +
  ggplot2::theme_classic(base_size = 8)
ggplot2::ggsave(file.path(PATHS$figures, "reference_marker_cellclass_counts.svg"), p_counts, width = 120, height = 90, units = "mm", device = svglite::svglite)

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(manifest = manifest_file, source_dirs = vapply(sources, function(x) repo_path(x$local_dir %||% ""), character(1))),
  outputs = list(
    evidence = file.path(PATHS$tables, "reference_marker_evidence_long.csv"),
    ranked = file.path(PATHS$tables, "reference_marker_sets_ranked.csv"),
    registry = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"),
    figures = PATHS$figures
  ),
  parameters = list(allow_download = allow_download, top_n_per_class = top_n, min_specificity = min_specificity),
  notes = "Reference marker evidence is cached/versioned/optional and used for annotation/reporting/sensitivity only."
)

message("Reference marker import complete. Registry: ", repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"))
