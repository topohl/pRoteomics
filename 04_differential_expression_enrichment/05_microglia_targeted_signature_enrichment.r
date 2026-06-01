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
  ) %>% dplyr::distinct()
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
      dplyr::transmute(term, gene = UNIPROT)
  )
  dplyr::bind_rows(term2gene, mapped) %>% dplyr::distinct()
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

parse_comparison <- function(comparison) {
  parts <- strsplit(sub("\\.csv$", "", basename(comparison)), "_", fixed = TRUE)[[1]]
  left <- parts[[1]]
  right <- if (length(parts) >= 2) parts[[2]] else NA_character_
  left_parsed <- extract_unit_condition(left)
  right_parsed <- extract_unit_condition(right)
  cond_pair <- paste(na.omit(c(left_parsed$condition, right_parsed$condition)), collapse = "_vs_")
  if (!nzchar(cond_pair)) cond_pair <- NA_character_
  tibble::tibble(
    comparison = sub("\\.csv$", "", basename(comparison)),
    left_unit = left_parsed$unit,
    right_unit = right_parsed$unit,
    left_condition = left_parsed$condition,
    right_condition = right_parsed$condition,
    condition_pair = cond_pair,
    region = parse_region_from_unit(left_parsed$unit),
    layer = ifelse(grepl("microglia", tolower(left)), NA_character_, sub("^[Cc][Aa][123]|^[Dd][Gg]", "", left_parsed$unit)),
    region_level_unit = collapse_neuropil_layer_to_region(left_parsed$unit)
  )
}

map_comparison_to_region_level <- function(comparison) {
  parse_comparison(comparison) %>% dplyr::select(.data$comparison, .data$region, .data$region_level_unit, .data$condition_pair)
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
    condition_pair = meta$condition_pair[[1]]
  )
}

load_ranked_contrasts <- function(dataset) {
  dir <- contrast_dir(dataset)
  files <- if (dir.exists(dir)) list.files(dir, pattern = "\\.csv$", full.names = TRUE) else character()
  rows <- purrr::map(files, prepare_ranked_contrast, dataset = dataset)
  dplyr::bind_rows(rows)
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
  )
}

run_one_gsea <- function(stats, term2gene, methods) {
  stats <- stats[!is.na(stats)]
  stats <- sort(stats, decreasing = TRUE)
  stats <- stats[!duplicated(names(stats))]
  pathways <- split(term2gene$gene, term2gene$term)
  pathways <- lapply(pathways, unique)
  pathways <- pathways[vapply(pathways, function(g) sum(g %in% names(stats)) >= 3, logical(1))]
  if (!length(pathways) || length(stats) < 10) return(tibble())

  if (methods$available[methods$method == "fgsea"]) {
    res <- tryCatch({
      if ("fgseaSimple" %in% getNamespaceExports("fgsea")) {
        fgsea::fgseaSimple(pathways = pathways, stats = stats, minSize = 3, maxSize = 500, nperm = 10000, nproc = 1)
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

  if (methods$available[methods$method == "limma_ranked_geneSetTest"]) {
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

  if (methods$available[methods$method == "clusterProfiler_GSEA"]) {
    res <- tryCatch(clusterProfiler::GSEA(stats, TERM2GENE = term2gene, minGSSize = 3, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE), error = function(e) NULL)
    if (!is.null(res) && nrow(as.data.frame(res))) {
      return(as.data.frame(res) %>%
        tibble::as_tibble() %>%
        dplyr::transmute(signature = .data$ID, method = "clusterProfiler_GSEA", NES = .data$NES, pvalue = .data$pvalue, padj = .data$p.adjust, set_size = .data$setSize, leading_edge = .data$core_enrichment))
    }
  }

  tibble()
}

run_signature_enrichment <- function(ranked, term2gene) {
  methods <- method_status()
  ranked %>%
    dplyr::group_split(.data$dataset, .data$comparison) %>%
    purrr::map_dfr(function(df) {
      stats <- df$rank_stat
      names(stats) <- df$gene
      res <- run_one_gsea(stats, term2gene, methods)
      if (!nrow(res)) return(tibble())
      meta <- df[1, c("dataset", "comparison", "region", "layer", "region_level_unit", "condition_pair", "input_file")]
      dplyr::bind_cols(meta[rep(1, nrow(res)), ], res) %>%
        dplyr::mutate(
          matched_genes = purrr::map_chr(.data$signature, ~paste(intersect(term2gene$gene[term2gene$term == .x], names(stats)), collapse = ";")),
          signature_gene_count = lengths(strsplit(.data$matched_genes, ";", fixed = TRUE))
        )
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
      microglia_signature_class = dplyr::case_when(
        (.data$microglia_sig | .data$rank_based_flag) & (.data$microglia_marker_fraction >= 0.10) & (!.data$neuropil_sig | .data$stronger_in_microglia) ~ "microglia_enriched",
        (.data$microglia_sig | .data$rank_based_flag) & .data$neuropil_sig & .data$same_direction_as_neuropil & .data$neuropil_marker_fraction >= 0.10 ~ "neuropil_shared",
        (.data$microglia_sig | .data$rank_based_flag) & .data$neuropil_sig & .data$same_direction_as_neuropil ~ "mixed_microenvironment",
        TRUE ~ "ambiguous"
      )
    )
}

write_empty_outputs <- function(reason, diagnostics) {
  empty <- tibble::tibble(dataset = DATASET, status = reason, run_id = RUN_ID)
  files <- c(
    "microglia_signature_enrichment.csv",
    "microglia_signature_enrichment_with_neuropil_reference.csv",
    "microglia_signature_enrichment_summary.csv",
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
  readr::write_csv(dry_diag, file.path(PATHS$tables, "microglia_signature_diagnostics.csv"))
  quit(status = 0, save = "no")
}

if (DATASET != "microglia") {
  diagnostics <- dplyr::bind_rows(diagnostics, tibble::tibble(check = "dataset", status = "WARN", detail = "Skipped because dataset is not microglia."))
  write_empty_outputs("skipped_non_microglia_dataset", diagnostics)
  quit(status = 0, save = "no")
}

term2gene <- make_term2gene(signature_sets)
id_map <- read_mouse_id_map()
term2gene <- expand_term2gene_for_uniprot(term2gene, id_map)

micro_ranked <- load_ranked_contrasts(DATASET)
neuropil_ranked <- load_ranked_contrasts(REFERENCE_DATASET)

diagnostics <- dplyr::bind_rows(
  diagnostics,
  tibble::tibble(
    check = c("signature_terms", "id_map_rows", "microglia_ranked_rows", "neuropil_ranked_rows"),
    status = c("PASS", ifelse(nrow(id_map) > 0, "PASS", "WARN"), ifelse(nrow(micro_ranked) > 0, "PASS", "WARN"), ifelse(nrow(neuropil_ranked) > 0, "PASS", "WARN")),
    detail = as.character(c(nrow(term2gene), nrow(id_map), nrow(micro_ranked), nrow(neuropil_ranked)))
  )
)

if (!nrow(micro_ranked)) {
  write_empty_outputs("missing_microglia_ranked_contrasts", diagnostics)
  warning("No readable microglia mapped contrast tables found under: ", contrast_dir(DATASET))
  quit(status = 0, save = "no")
}

micro_enrichment <- run_signature_enrichment(micro_ranked, term2gene)
neuropil_enrichment <- run_signature_enrichment(neuropil_ranked, term2gene)
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

summary_tbl <- with_reference %>%
  dplyr::count(.data$microglia_signature_class, .data$signature, .data$reference_match_type, name = "n_contrasts") %>%
  dplyr::arrange(.data$microglia_signature_class, dplyr::desc(.data$n_contrasts))

readr::write_csv(micro_enrichment, file.path(PATHS$tables, "microglia_signature_enrichment.csv"), na = "")
readr::write_csv(with_reference, file.path(PATHS$tables, "microglia_signature_enrichment_with_neuropil_reference.csv"), na = "")
readr::write_csv(summary_tbl, file.path(PATHS$tables, "microglia_signature_enrichment_summary.csv"), na = "")
readr::write_csv(diagnostics, file.path(PATHS$tables, "microglia_signature_diagnostics.csv"), na = "")
readr::write_csv(with_reference %>% dplyr::filter(.data$microglia_signature_class == "microglia_enriched"), file.path(PATHS$tables, "robust_microglia_signatures.csv"), na = "")
readr::write_csv(with_reference %>% dplyr::filter(.data$microglia_signature_class == "mixed_microenvironment"), file.path(PATHS$tables, "mixed_microenvironment_signatures.csv"), na = "")
readr::write_csv(with_reference %>% dplyr::filter(.data$microglia_signature_class == "neuropil_shared"), file.path(PATHS$tables, "neuropil_shared_signatures.csv"), na = "")
readr::write_csv(with_reference %>% dplyr::filter(.data$microglia_signature_class == "ambiguous"), file.path(PATHS$tables, "ambiguous_signatures.csv"), na = "")

if (nrow(with_reference)) {
  dot_df <- with_reference %>% dplyr::mutate(sig_label = ifelse(!is.na(.data$padj) & .data$padj < 0.05, "FDR<0.05", "nominal/unthresholded"))
  p1 <- ggplot2::ggplot(dot_df, ggplot2::aes(x = .data$comparison, y = .data$signature, color = .data$NES, size = -log10(pmax(.data$padj, .Machine$double.xmin)))) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "NES") +
    ggplot2::labs(x = NULL, y = NULL, size = "-log10 FDR") +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(file.path(PATHS$figures, "microglia_signature_dotplot.svg"), p1, width = max(7, length(unique(dot_df$comparison)) * 0.35), height = 4.8)

  scatter_df <- dot_df %>% dplyr::filter(!is.na(.data$neuropil_reference_NES))
  if (nrow(scatter_df)) {
    p2 <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = .data$neuropil_reference_NES, y = .data$NES, color = .data$microglia_signature_class)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey75") +
      ggplot2::geom_vline(xintercept = 0, color = "grey75") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
      ggplot2::geom_point(alpha = 0.8) +
      ggplot2::labs(x = "Neuropil reference NES", y = "Microglia ROI NES", color = "Class") +
      ggplot2::theme_minimal(base_size = 9)
    ggplot2::ggsave(file.path(PATHS$figures, "microglia_vs_neuropil_signature_NES_scatter.svg"), p2, width = 6.2, height = 4.6)
  }

  p3 <- with_reference %>%
    dplyr::count(.data$microglia_signature_class, name = "n") %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(.data$microglia_signature_class, .data$n), y = .data$n, fill = .data$microglia_signature_class)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Signature-contrast rows") +
    ggplot2::theme_minimal(base_size = 9)
  ggplot2::ggsave(file.path(PATHS$figures, "microglia_signature_classification_barplot.svg"), p3, width = 5.8, height = 3.8)
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
  "Microglia region-only contrasts are compared to neuropil references after collapsing neuropil layer units to parent regions where possible."
), file.path(PATHS$reports, "microglia_signature_methods_note.txt"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(microglia_contrasts = contrast_dir(DATASET), neuropil_contrasts = contrast_dir(REFERENCE_DATASET)),
  outputs = list(tables = PATHS$tables, figures = PATHS$figures),
  parameters = list(dataset = DATASET, reference_dataset = REFERENCE_DATASET, run_id = RUN_ID),
  notes = "Microglia signatures are interpreted for ROI-local microenvironment samples and compared to region-collapsed neuropil reference results."
)

message("Microglia-targeted signature enrichment completed.")
message("Annotated table: ", file.path(PATHS$tables, "microglia_signature_enrichment_with_neuropil_reference.csv"))
