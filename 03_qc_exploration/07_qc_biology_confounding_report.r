# ================================================================
# Script: 03_qc_exploration/07_qc_biology_confounding_report.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required results/tables/03_qc_exploration/05_pca_confounding_qc/<dataset>/; results/tables/03_qc_exploration/06_variance_partitioning/<dataset>/; optional results/tables/03_qc_exploration/04_marker_rank_abundance_qc/<dataset>/.
# Produces: results/reports/03_qc_exploration/07_qc_biology_confounding_report/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Summarizes QC biology/confounding evidence.
# ================================================================

# Combine QC, missingness, marker, PCA, and metadata into a compact confounding report.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
PATHS <- qc_paths("07_qc_biology_confounding_report", DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_CONFOUNDING_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_CONFOUNDING_METADATA_FILE")
if (run$dry_run) {
  status <- qc_dry_run_contract("03_qc_exploration/07_qc_biology_confounding_report.r", DATASET,
                                matrix_file = matrix_file, metadata_file = metadata_file, paths = PATHS,
                                extra = c("Also reuses upstream QC tables when present; otherwise recomputes minimal metrics."))
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "svglite", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))
if (!file.exists(matrix_file)) stop("Confounding report matrix not found: ", matrix_file, call. = FALSE)

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
raw_mat <- expr$mat
meta <- expr$meta
mat <- qc_impute_for_pca(raw_mat)
meta <- meta[colnames(mat), , drop = FALSE]

qc_metrics <- data.frame(
  Sample = colnames(raw_mat),
  missing_fraction = colMeans(is.na(raw_mat)),
  detected_proteins = colSums(!is.na(raw_mat)),
  median_abundance = apply(raw_mat, 2L, median, na.rm = TRUE),
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(meta, by = "Sample")

pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
pc_scores <- data.frame(Sample = rownames(pca$x), pca$x[, seq_len(min(10L, ncol(pca$x))), drop = FALSE], check.names = FALSE)
combined <- qc_metrics |>
  dplyr::left_join(pc_scores, by = "Sample")

marker_file <- file.path(module_paths("03_qc_exploration", file.path("04_marker_rank_abundance_qc", DATASET))$tables,
                         "marker_scores_by_sample.csv")
if (file.exists(marker_file)) {
  marker_scores <- utils::read.csv(marker_file, check.names = FALSE)
  marker_wide <- marker_scores |>
    dplyr::select(Sample, marker_panel, marker_score) |>
    tidyr::pivot_wider(names_from = marker_panel, values_from = marker_score, names_prefix = "marker_")
  combined <- combined |> dplyr::left_join(marker_wide, by = "Sample")
}

terms <- qc_metadata_terms(meta)
metric_cols <- setdiff(names(combined)[vapply(combined, is.numeric, logical(1))],
                       intersect(c("sampleNumber", "order"), names(combined)))
assoc <- dplyr::bind_rows(lapply(metric_cols, function(metric) {
  dplyr::bind_rows(lapply(terms, function(term) {
    data.frame(metric = metric, term = term,
               p_value = qc_safe_lm_p(combined[[metric]], combined[[term]]),
               eta2 = qc_eta2(combined[[metric]], combined[[term]]),
               stringsAsFactors = FALSE)
  }))
}))
assoc$q_value <- p.adjust(assoc$p_value, method = "BH")
assoc$flag <- ifelse(is.finite(assoc$eta2) & assoc$eta2 >= 0.25 & assoc$q_value <= 0.1, "FAIL",
                     ifelse(is.finite(assoc$eta2) & assoc$eta2 >= 0.15, "WARN", "PASS"))
assoc <- assoc[order(match(assoc$flag, c("FAIL", "WARN", "PASS")), -assoc$eta2), ]

qc_write_csv(combined, file.path(PATHS$tables, "qc_biology_combined_sample_table.csv"))
qc_write_csv(assoc, file.path(PATHS$tables, "qc_biology_confounding_associations.csv"))
qc_write_xlsx(list(sample_table = combined, associations = assoc),
              file.path(PATHS$tables, "qc_biology_confounding_report.xlsx"))

heat <- assoc |> dplyr::filter(is.finite(eta2))
if (nrow(heat)) {
  p <- ggplot(heat, aes(term, metric, fill = eta2)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = "#B2182B", labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = NULL, fill = "Eta2") +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(PATHS$figures, "qc_biology_association_heatmap.svg"), p,
         width = 160, height = 120, units = "mm", device = svglite::svglite)
}

top_flags <- assoc |> dplyr::filter(flag != "PASS")
status <- if (any(top_flags$flag == "FAIL")) "FAIL" else if (nrow(top_flags)) "WARN" else "PASS"
lines <- c("# QC-to-Biology Confounding Report", "", paste0("Dataset: ", DATASET),
           paste0(status, ": strongest QC/metadata associations summarized in qc_biology_confounding_associations.csv."),
           "Interpret Group biology cautiously when Group, plate/batch, missingness, marker scores, or early PCs are associated.")
if (nrow(top_flags)) {
  lines <- c(lines, "", "## Flags",
             paste0("- ", top_flags$flag, ": ", top_flags$metric, " ~ ", top_flags$term,
                    " eta2=", signif(top_flags$eta2, 3), " q=", signif(top_flags$q_value, 3)))
}
writeLines(lines, file.path(PATHS$reports, "qc_biology_confounding_summary.md"))

write_run_manifest(file.path(PATHS$logs, "run_manifest.yml"),
                   inputs = list(matrix = matrix_file, metadata = metadata_file, marker_scores = marker_file),
                   outputs = list(tables = PATHS$tables, figures = PATHS$figures, reports = PATHS$reports),
                   parameters = list(dataset = DATASET),
                   notes = "Compact PASS/WARN/FAIL screen for QC metrics, markers, PCs, and metadata.")

message("QC-to-biology confounding report complete for dataset: ", DATASET)
