#!/usr/bin/env Rscript
#
# Discover experiment-specific ROI-enrichment marker sets across neuron neuropil, neuron soma, and microglia.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))

dry_run <- is_dry_run()
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
SUBSTEP_ID <- "05_empirical_roi_marker_discovery"
PATHS <- create_module_dirs("03_qc_exploration", SUBSTEP_ID)
min_detection <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_EMPIRICAL_MARKER_MIN_DETECTION", unset = "0.30")))
min_abs_logfc <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_EMPIRICAL_MARKER_MIN_ABS_LOGFC", unset = "0.50")))
fdr_threshold <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_EMPIRICAL_MARKER_FDR", unset = "0.10")))
if (!is.finite(min_detection)) min_detection <- 0.30
if (!is.finite(min_abs_logfc)) min_abs_logfc <- 0.50
if (!is.finite(fdr_threshold)) fdr_threshold <- 0.10

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

inputs <- setNames(lapply(DATASETS, function(ds) {
  list(matrix = qc_resolve_matrix(ds, env = paste0("PROTEOMICS_EMPIRICAL_MARKER_", toupper(ds), "_MATRIX_FILE")),
       metadata = qc_resolve_metadata(ds, env = paste0("PROTEOMICS_EMPIRICAL_MARKER_", toupper(ds), "_METADATA_FILE")))
}), DATASETS)

if (dry_run) {
  invisible(lapply(unlist(PATHS), dir_create))
  dry_run_line("Script", "03_qc_exploration/05_empirical_roi_marker_discovery.r")
  for (ds in DATASETS) {
    dry_run_line(paste(ds, "matrix"), inputs[[ds]]$matrix, if (file.exists(inputs[[ds]]$matrix)) "PASS" else "WARN")
    dry_run_line(paste(ds, "metadata"), inputs[[ds]]$metadata, if (file.exists(inputs[[ds]]$metadata)) "PASS" else "WARN")
  }
  dry_run_line("Output", file.path(PATHS$tables, "empirical_roi_marker_sets.csv"), "INFO")
  quit(status = 0, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "svglite", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

read_dataset <- function(ds) {
  if (!file.exists(inputs[[ds]]$matrix)) stop("Missing empirical marker matrix for ", ds, ": ", inputs[[ds]]$matrix, call. = FALSE)
  expr <- qc_read_expression(inputs[[ds]]$matrix, inputs[[ds]]$metadata, ds)
  mat <- expr$mat
  gene <- rownames(mat)
  data.frame(
    dataset = ds,
    ProteinID = gene,
    GeneSymbol = gene,
    gene_token = normalize_gene_token(gene),
    detection_rate = rowMeans(is.finite(mat) & !is.na(mat)),
    mean_abundance = rowMeans(mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(nzchar(.data$gene_token)) |>
    dplyr::group_by(.data$gene_token) |>
    dplyr::summarise(
      ProteinID = dplyr::first(.data$ProteinID),
      GeneSymbol = dplyr::first(.data$GeneSymbol),
      detection_rate = safe_max(.data$detection_rate),
      mean_abundance = safe_max(.data$mean_abundance),
      .groups = "drop"
    ) |>
    dplyr::mutate(dataset = ds, .before = "ProteinID")
}

stats_long <- dplyr::bind_rows(lapply(DATASETS, read_dataset))
wide <- stats_long |>
  dplyr::select("gene_token", "ProteinID", "GeneSymbol", "dataset", "detection_rate", "mean_abundance") |>
  tidyr::pivot_wider(
    names_from = "dataset",
    values_from = c("detection_rate", "mean_abundance", "ProteinID", "GeneSymbol"),
    values_fn = list(detection_rate = max, mean_abundance = max, ProteinID = dplyr::first, GeneSymbol = dplyr::first)
  )

first_nonmissing <- function(...) {
  vals <- c(...)
  vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
  if (length(vals)) vals[[1]] else NA_character_
}
safe_diff <- function(a, b) ifelse(is.finite(a) & is.finite(b), a - b, NA_real_)

wide <- wide |>
  dplyr::mutate(
    ProteinID = mapply(first_nonmissing, .data$ProteinID_microglia, .data$ProteinID_neuron_neuropil, .data$ProteinID_neuron_soma),
    GeneSymbol = mapply(first_nonmissing, .data$GeneSymbol_microglia, .data$GeneSymbol_neuron_neuropil, .data$GeneSymbol_neuron_soma),
    detection_rate_microglia = .data$detection_rate_microglia,
    detection_rate_neuropil = .data$detection_rate_neuron_neuropil,
    detection_rate_soma = .data$detection_rate_neuron_soma,
    logFC_microglia_vs_neuropil = safe_diff(.data$mean_abundance_microglia, .data$mean_abundance_neuron_neuropil),
    logFC_microglia_vs_soma = safe_diff(.data$mean_abundance_microglia, .data$mean_abundance_neuron_soma),
    logFC_neuropil_vs_microglia = -.data$logFC_microglia_vs_neuropil,
    logFC_soma_vs_microglia = -.data$logFC_microglia_vs_soma,
    p_value = NA_real_,
    FDR = NA_real_,
    model_type = "detection_logFC_fallback"
  )

pvals <- lapply(split(stats_long, stats_long$gene_token), function(tab) {
  if (nrow(tab) < 3L || length(unique(tab$dataset)) < 2L) return(NA_real_)
  fit <- suppressWarnings(try(stats::lm(mean_abundance ~ dataset, data = tab), silent = TRUE))
  if (inherits(fit, "try-error")) return(NA_real_)
  a <- suppressWarnings(stats::anova(fit))
  if ("Pr(>F)" %in% names(a)) a[1, "Pr(>F)"] else NA_real_
})
p_df <- data.frame(gene_token = names(pvals), p_value = as.numeric(pvals), stringsAsFactors = FALSE)
p_df$FDR <- stats::p.adjust(p_df$p_value, method = "BH")
wide <- wide |>
  dplyr::select(-"p_value", -"FDR") |>
  dplyr::left_join(p_df, by = "gene_token") |>
  dplyr::mutate(
    FDR_pass = is.na(.data$FDR) | .data$FDR <= fdr_threshold,
    dataset_enriched_in = dplyr::case_when(
      .data$logFC_microglia_vs_neuropil >= min_abs_logfc & .data$logFC_microglia_vs_soma >= min_abs_logfc ~ "microglia",
      .data$logFC_neuropil_vs_microglia >= min_abs_logfc ~ "neuron_neuropil",
      .data$logFC_soma_vs_microglia >= min_abs_logfc ~ "neuron_soma",
      TRUE ~ "shared_or_ambiguous"
    )
  )

make_set <- function(marker_set, keep, comparison, confidence, notes) {
  tab <- wide[keep & wide$FDR_pass, , drop = FALSE]
  if (!nrow(tab)) return(data.frame())
  tab |>
    dplyr::transmute(
      marker_set = marker_set,
      ProteinID, GeneSymbol, dataset_enriched_in, comparison,
      logFC_microglia_vs_neuropil, logFC_microglia_vs_soma,
      logFC_neuropil_vs_microglia, logFC_soma_vs_microglia,
      p_value, FDR,
      detection_rate_microglia, detection_rate_neuropil, detection_rate_soma,
      marker_source = "empirical_roi_marker_sets",
      selection_rule = paste0("min_detection=", min_detection, "; min_abs_logFC=", min_abs_logfc, "; FDR<=", fdr_threshold, "; no CON/RES/SUS marker contrast"),
      confidence = confidence,
      notes = notes,
      model_type = .data$model_type
    )
}

strong_detection <- min(0.80, min_detection + 0.20)
strong_logfc <- min_abs_logfc + 0.50
out <- dplyr::bind_rows(
  make_set(
    "empirical_microglia_roi_enriched",
    wide$detection_rate_microglia >= min_detection & wide$logFC_microglia_vs_neuropil >= min_abs_logfc & wide$logFC_microglia_vs_soma >= min_abs_logfc,
    "microglia_vs_neuron_neuropil_and_soma", "empirical_roi_enriched",
    "ROI enrichment only; microglia ROI samples are not purified microglia."
  ),
  make_set(
    "empirical_microglia_roi_high_confidence",
    wide$detection_rate_microglia >= strong_detection & wide$logFC_microglia_vs_neuropil >= strong_logfc & wide$logFC_microglia_vs_soma >= strong_logfc,
    "microglia_vs_neuron_neuropil_and_soma", "empirical_high_confidence",
    "Stronger empirical ROI enrichment; annotation only."
  ),
  make_set(
    "empirical_neuropil_enriched",
    wide$detection_rate_neuropil >= min_detection & wide$logFC_neuropil_vs_microglia >= min_abs_logfc,
    "neuron_neuropil_vs_microglia", "empirical_roi_enriched",
    "Neuropil-sensitive marker evidence; do not subtract from WGCNA."
  ),
  make_set(
    "empirical_neuropil_sensitive_high_confidence",
    wide$detection_rate_neuropil >= strong_detection & wide$logFC_neuropil_vs_microglia >= strong_logfc,
    "neuron_neuropil_vs_microglia", "empirical_high_confidence",
    "High-confidence neuropil-sensitive marker evidence."
  ),
  make_set(
    "empirical_soma_enriched",
    wide$detection_rate_soma >= min_detection & wide$logFC_soma_vs_microglia >= min_abs_logfc,
    "neuron_soma_vs_microglia", "empirical_roi_enriched",
    "Soma-enriched marker evidence; annotation only."
  ),
  make_set(
    "empirical_microglia_neuropil_shared",
    wide$detection_rate_microglia >= min_detection & wide$detection_rate_neuropil >= min_detection &
      abs(wide$logFC_microglia_vs_neuropil) < min_abs_logfc & wide$dataset_enriched_in == "shared_or_ambiguous",
    "microglia_neuropil_shared_detection", "empirical_shared",
    "Shared local microenvironment/carryover interpretation set."
  )
)

if (!nrow(out)) {
  out <- data.frame(
    marker_set = character(), ProteinID = character(), GeneSymbol = character(),
    dataset_enriched_in = character(), comparison = character(),
    logFC_microglia_vs_neuropil = numeric(), logFC_microglia_vs_soma = numeric(),
    logFC_neuropil_vs_microglia = numeric(), logFC_soma_vs_microglia = numeric(),
    p_value = numeric(), FDR = numeric(), detection_rate_microglia = numeric(),
    detection_rate_neuropil = numeric(), detection_rate_soma = numeric(),
    marker_source = character(), selection_rule = character(), confidence = character(),
    notes = character(), model_type = character()
  )
}

write_table_and_source(out, PATHS$tables, PATHS$source_data, "empirical_roi_marker_sets.csv")

heat <- wide |>
  dplyr::select("GeneSymbol", "detection_rate_microglia", "detection_rate_neuropil", "detection_rate_soma") |>
  dplyr::slice_max(order_by = pmax(.data$detection_rate_microglia, .data$detection_rate_neuropil, .data$detection_rate_soma, na.rm = TRUE), n = 75) |>
  tidyr::pivot_longer(-"GeneSymbol", names_to = "dataset", values_to = "detection_rate")
p1 <- ggplot2::ggplot(heat, ggplot2::aes(x = .data$dataset, y = .data$GeneSymbol, fill = .data$detection_rate)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "#2F6F73", na.value = "grey90") +
  ggplot2::labs(x = NULL, y = NULL, fill = "Detection") +
  ggplot2::theme_classic(base_size = 7)
ggplot2::ggsave(file.path(PATHS$figures, "empirical_marker_detection_heatmap.svg"), p1, width = 120, height = 150, units = "mm", device = svglite::svglite)

lfc <- wide |>
  dplyr::select("GeneSymbol", "logFC_microglia_vs_neuropil", "logFC_microglia_vs_soma", "logFC_neuropil_vs_microglia", "logFC_soma_vs_microglia") |>
  tidyr::pivot_longer(-"GeneSymbol", names_to = "comparison", values_to = "logFC")
p2 <- ggplot2::ggplot(lfc, ggplot2::aes(x = .data$comparison, y = .data$logFC)) +
  ggplot2::geom_hline(yintercept = c(-min_abs_logfc, min_abs_logfc), linetype = "dashed", color = "grey55") +
  ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey92", color = "grey35") +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Mean abundance difference") +
  ggplot2::theme_classic(base_size = 8)
ggplot2::ggsave(file.path(PATHS$figures, "empirical_marker_logFC_summary.svg"), p2, width = 140, height = 90, units = "mm", device = svglite::svglite)

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = inputs,
  outputs = list(markers = file.path(PATHS$tables, "empirical_roi_marker_sets.csv"), figures = PATHS$figures),
  parameters = list(min_detection = min_detection, min_abs_logFC = min_abs_logfc, FDR = fdr_threshold, model_type = "detection_logFC_fallback"),
  notes = "Empirical marker discovery uses dataset/ROI contrasts only; StressGroup is not a marker-defining contrast."
)

message("Empirical ROI marker discovery complete: ", file.path(PATHS$tables, "empirical_roi_marker_sets.csv"))
