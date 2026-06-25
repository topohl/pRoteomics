#!/usr/bin/env Rscript

# ================================================================
# Script: 03_qc_exploration/00_dataset_qc_report.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv; optional data/metadata/*.xlsx.
# Produces: results/tables/03_qc_exploration/00_dataset_qc_report/<dataset>/dataset_qc_report.xlsx; results/reports/03_qc_exploration/00_dataset_qc_report/<dataset>/dataset_qc_summary.md.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Baseline dataset QC.
# ================================================================

# Canonical one-stop dataset QC report.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
PATHS <- qc_paths("00_dataset_qc_report", DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_DATASET_QC_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_DATASET_QC_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/00_dataset_qc_report.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c(
      "Writes missingness, imputation footprint, sample/protein counts, PCA, metadata structure, abundance distributions, and outlier flags.",
      "Use PROTEOMICS_DATASET_QC_MATRIX_FILE to point at a pre-imputation matrix when available."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(matrix_file)) {
  stop("Dataset QC matrix not found for dataset '", DATASET, "': ", matrix_file,
       ". Set PROTEOMICS_DATASET_QC_MATRIX_FILE or PROTEOMICS_QC_MATRIX_FILE.", call. = FALSE)
}

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
raw_mat <- expr$mat
meta <- expr$meta
input_appears_imputed <- grepl("imput", basename(matrix_file), ignore.case = TRUE)

sample_qc <- data.frame(
  Sample = colnames(raw_mat),
  n_proteins_total = nrow(raw_mat),
  n_proteins_detected = colSums(!is.na(raw_mat)),
  missing_fraction = colMeans(is.na(raw_mat)),
  median_abundance = apply(raw_mat, 2L, median, na.rm = TRUE),
  mean_abundance = colMeans(raw_mat, na.rm = TRUE),
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(meta, by = "Sample")

protein_qc <- data.frame(
  Protein = rownames(raw_mat),
  n_samples_total = ncol(raw_mat),
  n_samples_detected = rowSums(!is.na(raw_mat)),
  missing_fraction = rowMeans(is.na(raw_mat)),
  median_abundance = apply(raw_mat, 1L, median, na.rm = TRUE),
  stringsAsFactors = FALSE
)

mat_for_pca <- qc_impute_for_pca(raw_mat)
imputation_footprint <- data.frame(
  dataset = DATASET,
  input_matrix = matrix_file,
  input_appears_imputed = input_appears_imputed,
  n_samples = ncol(raw_mat),
  n_proteins_raw = nrow(raw_mat),
  n_proteins_pca = nrow(mat_for_pca),
  n_missing_raw = sum(is.na(raw_mat)),
  missing_fraction_raw = mean(is.na(raw_mat)),
  n_values_median_imputed_for_pca = sum(is.na(raw_mat[rownames(mat_for_pca), colnames(mat_for_pca), drop = FALSE])),
  imputation_note = if (input_appears_imputed) {
    "Input filename suggests upstream imputation; pre-imputation missingness may be underestimated."
  } else {
    "Median imputation is used only inside this report for PCA calculation."
  },
  stringsAsFactors = FALSE
)

robust_z <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  mad <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(mad) || mad == 0) return(rep(NA_real_, length(x)))
  (x - med) / mad
}
outlier_metrics <- intersect(c("missing_fraction", "n_proteins_detected", "median_abundance", "mean_abundance"), names(sample_qc))
outlier_z <- as.data.frame(lapply(sample_qc[outlier_metrics], robust_z), check.names = FALSE)
names(outlier_z) <- paste0(outlier_metrics, "_robust_z")
outlier_z_abs <- lapply(outlier_z, abs)
sample_outliers <- dplyr::bind_cols(sample_qc["Sample"], outlier_z) |>
  dplyr::mutate(
    n_outlier_metrics = rowSums(abs(dplyr::across(dplyr::ends_with("_robust_z"))) > 3, na.rm = TRUE),
    max_abs_robust_z = do.call(pmax, c(outlier_z_abs, na.rm = TRUE)),
    outlier_flag = dplyr::case_when(
      .data$n_outlier_metrics >= 2 ~ "FAIL",
      .data$max_abs_robust_z > 3 ~ "WARN",
      TRUE ~ "PASS"
    )
  )

terms <- qc_metadata_terms(meta)
sample_structure <- dplyr::bind_rows(lapply(terms, function(term) {
  data.frame(
    dataset = DATASET,
    metadata_term = term,
    level = as.character(names(table(meta[[term]], useNA = "ifany"))),
    n_samples = as.integer(table(meta[[term]], useNA = "ifany")),
    stringsAsFactors = FALSE
  )
}))

pca <- stats::prcomp(t(mat_for_pca), center = TRUE, scale. = TRUE)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
pca_scores <- data.frame(Sample = rownames(pca$x), pca$x[, seq_len(min(10L, ncol(pca$x))), drop = FALSE], check.names = FALSE) |>
  dplyr::left_join(sample_qc, by = "Sample")
pc_terms <- intersect(paste0("PC", 1:min(10L, ncol(pca$x))), names(pca_scores))
pc_assoc <- dplyr::bind_rows(lapply(pc_terms, function(pc) {
  dplyr::bind_rows(lapply(terms, function(term) {
    data.frame(
      PC = pc,
      metadata_term = term,
      p_value = qc_safe_lm_p(pca_scores[[pc]], pca_scores[[term]]),
      eta2 = qc_eta2(pca_scores[[pc]], pca_scores[[term]]),
      stringsAsFactors = FALSE
    )
  }))
}))
if (nrow(pc_assoc)) pc_assoc$q_value <- p.adjust(pc_assoc$p_value, method = "BH")

qc_write_csv(sample_qc, file.path(PATHS$tables, "dataset_qc_sample_metrics.csv"))
qc_write_csv(protein_qc, file.path(PATHS$tables, "dataset_qc_protein_metrics.csv"))
qc_write_csv(imputation_footprint, file.path(PATHS$tables, "dataset_qc_imputation_footprint.csv"))
qc_write_csv(sample_outliers, file.path(PATHS$tables, "dataset_qc_outlier_flags.csv"))
qc_write_csv(sample_structure, file.path(PATHS$tables, "dataset_qc_metadata_structure.csv"))
qc_write_csv(pca_scores, file.path(PATHS$tables, "dataset_qc_pca_scores.csv"))
qc_write_csv(pc_assoc, file.path(PATHS$tables, "dataset_qc_pca_metadata_associations.csv"))
qc_write_xlsx(
  list(
    sample_metrics = sample_qc,
    protein_metrics = protein_qc,
    imputation_footprint = imputation_footprint,
    outlier_flags = sample_outliers,
    metadata_structure = sample_structure,
    pca_metadata = pc_assoc
  ),
  file.path(PATHS$tables, "dataset_qc_report.xlsx")
)

theme_qc <- function() ggplot2::theme_classic(base_size = 8) + ggplot2::theme(legend.position = "bottom")
theme_embedding <- function() qc_embedding_theme(base_size = 8)
color_term <- intersect(c("Group", "group", "ExpGroup", "Region", "region", "Layer", "layer", "batch", "plate", "Sex", "sex"), names(pca_scores))[1]
if (length(color_term) && !is.na(color_term)) {
  pca_plot_df <- pca_scores
  pca_plot_df[[color_term]] <- factor(pca_plot_df[[color_term]])
  p_pca <- ggplot2::ggplot(pca_plot_df, ggplot2::aes(PC1, PC2, color = .data[[color_term]])) +
    ggplot2::geom_point(size = 2, alpha = 0.85) +
    ggplot2::scale_color_manual(values = qc_embedding_palette(nlevels(pca_plot_df[[color_term]]))) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_explained[[1]]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_explained[[2]]),
      color = color_term
    ) +
    ggplot2::coord_equal() +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = qc_legend_rows(pca_plot_df[[color_term]]), byrow = TRUE, override.aes = list(alpha = 1, size = 2))) +
    theme_embedding()
  qc_save_square_svg(file.path(PATHS$figures, "dataset_qc_pca.svg"), p_pca, size_mm = 90)
}

if (run$run_embeddings) {
  message("Dataset QC UMAP/t-SNE are exploratory only and should not be used as primary evidence.")
  embedding_input <- pca$x[, seq_len(min(20L, ncol(pca$x))), drop = FALSE]
  if (requireNamespace("uwot", quietly = TRUE)) {
    set.seed(20260623)
    umap <- uwot::umap(embedding_input, n_threads = 1)
    umap_scores <- data.frame(Sample = rownames(pca$x), UMAP1 = umap[, 1], UMAP2 = umap[, 2], check.names = FALSE) |>
      dplyr::left_join(sample_qc, by = "Sample")
    qc_write_csv(umap_scores, file.path(PATHS$tables, "dataset_qc_exploratory_umap_scores.csv"))
    if (length(color_term) && !is.na(color_term) && color_term %in% names(umap_scores)) {
      keep <- !is.na(umap_scores[[color_term]]) & nzchar(as.character(umap_scores[[color_term]]))
      if (sum(keep) >= 3L && length(unique(umap_scores[[color_term]][keep])) >= 2L) {
        umap_plot_df <- umap_scores[keep, ]
        umap_plot_df[[color_term]] <- factor(umap_plot_df[[color_term]])
        p_umap <- ggplot2::ggplot(umap_plot_df, ggplot2::aes(UMAP1, UMAP2, color = .data[[color_term]])) +
          ggplot2::geom_point(size = 2, alpha = 0.85) +
          ggplot2::scale_color_manual(values = qc_embedding_palette(nlevels(umap_plot_df[[color_term]]))) +
          ggplot2::labs(
            title = "Exploratory UMAP",
            subtitle = "Not primary evidence; computed from the first up to 20 PCs",
            x = "UMAP1",
            y = "UMAP2",
            color = color_term
          ) +
          ggplot2::coord_equal() +
          ggplot2::guides(color = ggplot2::guide_legend(nrow = qc_legend_rows(umap_plot_df[[color_term]]), byrow = TRUE, override.aes = list(alpha = 1, size = 2))) +
          theme_embedding()
        qc_save_square_svg(file.path(PATHS$figures, "dataset_qc_umap.svg"), p_umap, size_mm = 90)
      }
    }
  } else {
    warning("Package 'uwot' is not installed; skipped dataset QC exploratory UMAP outputs.", call. = FALSE)
  }
  if (requireNamespace("Rtsne", quietly = TRUE)) {
    perplexity <- min(30L, floor((nrow(embedding_input) - 1L) / 3L))
    if (perplexity >= 1L) {
      set.seed(20260623)
      tsne <- Rtsne::Rtsne(embedding_input, dims = 2, perplexity = perplexity, pca = FALSE,
                           theta = 0.5, check_duplicates = FALSE, verbose = FALSE)
      tsne_scores <- data.frame(Sample = rownames(pca$x), TSNE1 = tsne$Y[, 1], TSNE2 = tsne$Y[, 2],
                                tsne_perplexity = perplexity, check.names = FALSE) |>
        dplyr::left_join(sample_qc, by = "Sample")
      qc_write_csv(tsne_scores, file.path(PATHS$tables, "dataset_qc_exploratory_tsne_scores.csv"))
      if (length(color_term) && !is.na(color_term) && color_term %in% names(tsne_scores)) {
        keep <- !is.na(tsne_scores[[color_term]]) & nzchar(as.character(tsne_scores[[color_term]]))
        if (sum(keep) >= 3L && length(unique(tsne_scores[[color_term]][keep])) >= 2L) {
          tsne_plot_df <- tsne_scores[keep, ]
          tsne_plot_df[[color_term]] <- factor(tsne_plot_df[[color_term]])
          p_tsne <- ggplot2::ggplot(tsne_plot_df, ggplot2::aes(TSNE1, TSNE2, color = .data[[color_term]])) +
            ggplot2::geom_point(size = 2, alpha = 0.85) +
            ggplot2::scale_color_manual(values = qc_embedding_palette(nlevels(tsne_plot_df[[color_term]]))) +
            ggplot2::labs(
              title = "Exploratory t-SNE",
              subtitle = sprintf("Not primary evidence; computed from the first up to 20 PCs, perplexity %d", perplexity),
              x = "t-SNE1",
              y = "t-SNE2",
              color = color_term
            ) +
            ggplot2::coord_equal() +
            ggplot2::guides(color = ggplot2::guide_legend(nrow = qc_legend_rows(tsne_plot_df[[color_term]]), byrow = TRUE, override.aes = list(alpha = 1, size = 2))) +
            theme_embedding()
          qc_save_square_svg(file.path(PATHS$figures, "dataset_qc_tsne.svg"), p_tsne, size_mm = 90)
        }
      }
    } else {
      warning("Too few samples for dataset QC exploratory t-SNE; skipped t-SNE outputs.", call. = FALSE)
    }
  } else {
    warning("Package 'Rtsne' is not installed; skipped dataset QC exploratory t-SNE outputs.", call. = FALSE)
  }
}

p_missing <- ggplot2::ggplot(sample_qc, ggplot2::aes(stats::reorder(.data$Sample, .data$missing_fraction), .data$missing_fraction)) +
  ggplot2::geom_col(fill = "#4C78A8") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = NULL, y = "Missing fraction") +
  theme_qc() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
ggplot2::ggsave(file.path(PATHS$figures, "dataset_qc_missingness_by_sample.svg"), p_missing,
                width = 160, height = 80, units = "mm", device = svglite::svglite)

abund_df <- as.data.frame(raw_mat) |>
  tibble::rownames_to_column("Protein") |>
  tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "abundance") |>
  dplyr::filter(is.finite(.data$abundance))
p_abund <- ggplot2::ggplot(abund_df, ggplot2::aes(.data$abundance)) +
  ggplot2::geom_histogram(bins = 80, fill = "#777777", color = "white", linewidth = 0.1) +
  ggplot2::labs(x = "Abundance", y = "Values") +
  theme_qc()
ggplot2::ggsave(file.path(PATHS$figures, "dataset_qc_abundance_distribution.svg"), p_abund,
                width = 120, height = 75, units = "mm", device = svglite::svglite)

if (nrow(pc_assoc) && any(is.finite(pc_assoc$eta2))) {
  p_assoc <- ggplot2::ggplot(pc_assoc, ggplot2::aes(.data$metadata_term, .data$PC, fill = .data$eta2)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", labels = scales::percent_format(accuracy = 1), na.value = "grey90") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Eta2") +
    theme_qc() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(PATHS$figures, "dataset_qc_pca_metadata_heatmap.svg"), p_assoc,
                  width = 150, height = 90, units = "mm", device = svglite::svglite)
}

summary_lines <- c(
  "# Dataset QC Report",
  "",
  paste0("Dataset: ", DATASET),
  paste0("Matrix: ", relative_to(matrix_file)),
  paste0("Metadata: ", ifelse(file.exists(metadata_file), relative_to(metadata_file), "not found")),
  paste0("Samples: ", ncol(raw_mat)),
  paste0("Proteins: ", nrow(raw_mat)),
  paste0("Raw missing fraction: ", signif(mean(is.na(raw_mat)), 3)),
  paste0("Input appears imputed: ", input_appears_imputed),
  paste0("Outlier flags: PASS=", sum(sample_outliers$outlier_flag == "PASS"), "; WARN=", sum(sample_outliers$outlier_flag == "WARN"), "; FAIL=", sum(sample_outliers$outlier_flag == "FAIL")),
  paste0("Exploratory UMAP/t-SNE requested: ", run$run_embeddings, " (not primary evidence)"),
  "",
  "No statistics are invented here: unavailable metadata terms simply do not appear in the association tables."
)
writeLines(summary_lines, file.path(PATHS$reports, "dataset_qc_summary.md"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(matrix = matrix_file, metadata = metadata_file),
  outputs = list(tables = PATHS$tables, figures = PATHS$figures, reports = PATHS$reports),
  parameters = list(dataset = DATASET, input_appears_imputed = input_appears_imputed, run_embeddings = run$run_embeddings),
  notes = "Canonical dataset-level QC report combining missingness, imputation footprint, counts, PCA, metadata structure, abundance distributions, and outlier flags."
)

message("Dataset QC report complete for dataset: ", DATASET)
