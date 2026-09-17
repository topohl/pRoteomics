# ================================================================
# Script: 03_qc_exploration/05_pca_confounding_qc.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv; optional data/metadata/*.xlsx.
# Produces: results/tables/03_qc_exploration/05_pca_confounding_qc/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: PCA and confounding QC.
# ================================================================

# Dataset-aware PCA QC for spatial proteomics.
# UMAP/t-SNE/clustering are explicitly optional exploratory outputs.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "05_pca_confounding_qc"
PATHS <- qc_paths(SUBSTEP_ID, DATASET)
matrix_file <- path_or_env("PROTEOMICS_PCA_MATRIX_FILE", qc_resolve_matrix(DATASET), must_exist = FALSE)
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_PCA_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/05_pca_confounding_qc.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c(
      paste0("Optional embeddings enabled: ", run$run_embeddings),
      paste0("Optional clustering enabled: ", run$run_clustering),
      "Writes PCA_confounding_summary.csv for PC1-PC10."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "svglite", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(matrix_file)) {
  stop("PCA matrix not found for dataset '", DATASET, "': ", matrix_file,
       ". Set PROTEOMICS_PCA_MATRIX_FILE or PROTEOMICS_QC_MATRIX_FILE.", call. = FALSE)
}

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
raw_mat <- expr$mat
meta <- expr$meta
missing_sample <- data.frame(Sample = colnames(raw_mat), missing_fraction = colMeans(is.na(raw_mat)))
mat <- qc_impute_for_pca(raw_mat)
if (nrow(mat) < 2L || ncol(mat) < 2L) stop("Matrix too small for PCA after filtering.", call. = FALSE)
meta <- meta[colnames(mat), , drop = FALSE]

pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
scores <- data.frame(Sample = rownames(pca$x), pca$x, check.names = FALSE) |>
  dplyr::left_join(meta, by = "Sample") |>
  dplyr::left_join(missing_sample, by = "Sample")
loadings <- data.frame(Protein = rownames(pca$rotation), pca$rotation, check.names = FALSE)
scree <- data.frame(PC = paste0("PC", seq_along(var_explained)),
                    PC_index = seq_along(var_explained),
                    variance_explained = var_explained,
                    cumulative_variance = cumsum(var_explained))

qc_write_csv(scores, file.path(PATHS$tables, "pca_scores.csv"))
qc_write_csv(loadings, file.path(PATHS$tables, "pca_loadings.csv"))
qc_write_csv(scree, file.path(PATHS$tables, "pca_variance_explained.csv"))

theme_qc <- function() {
  ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(face = "bold"))
}
theme_embedding <- function() qc_embedding_theme(base_size = 8)

plot_pca_by <- function(term) {
  if (!term %in% names(scores)) return(invisible(NULL))
  keep <- !is.na(scores[[term]]) & nzchar(as.character(scores[[term]]))
  if (sum(keep) < 3L || length(unique(scores[[term]][keep])) < 2L) return(invisible(NULL))
  plot_df <- scores[keep, ]
  plot_df[[term]] <- factor(plot_df[[term]])
  p <- ggplot(plot_df, aes(PC1, PC2, color = .data[[term]])) +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_manual(values = qc_embedding_palette(nlevels(plot_df[[term]]))) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * var_explained[[1]]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_explained[[2]]),
      color = term
    ) +
    coord_equal() +
    guides(color = guide_legend(nrow = qc_legend_rows(plot_df[[term]]), byrow = TRUE, override.aes = list(alpha = 1, size = 2))) +
    theme_embedding()
  qc_save_square_svg(file.path(PATHS$figures, paste0("pca_by_", safe_filename(term), ".svg")), p, size_mm = 90)
}

terms <- qc_metadata_terms(meta)
invisible(lapply(unique(c("Group", "group", "ExpGroup", "Region", "region", "Layer", "layer", "ReplicateGroup", "AnimalID", "plate", terms)), plot_pca_by))

p_scree <- ggplot(scree[seq_len(min(20L, nrow(scree))), ], aes(PC_index, variance_explained)) +
  geom_col(fill = "#4C78A8") +
  geom_line(aes(y = cumulative_variance), color = "#D55E00") +
  geom_point(aes(y = cumulative_variance), color = "#D55E00", size = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Principal component", y = "Variance explained / cumulative") +
  theme_qc()
ggsave(file.path(PATHS$figures, "pca_scree_cumulative.svg"), p_scree,
       width = 110, height = 80, units = "mm", device = svglite::svglite)

pc_terms <- intersect(paste0("PC", 1:min(10L, ncol(pca$x))), names(scores))
assoc <- dplyr::bind_rows(lapply(pc_terms, function(pc) {
  dplyr::bind_rows(lapply(terms, function(term) {
    data.frame(
      PC = pc,
      term = term,
      p_value = qc_safe_lm_p(scores[[pc]], scores[[term]]),
      eta2 = qc_eta2(scores[[pc]], scores[[term]]),
      stringsAsFactors = FALSE
    )
  }))
}))
if (nrow(assoc)) {
  assoc$q_value <- p.adjust(assoc$p_value, method = "BH")
  assoc <- assoc[order(assoc$PC, -assoc$eta2), ]
  qc_write_csv(assoc, file.path(PATHS$tables, "pc_metadata_associations.csv"))
  finite_max <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    max(x)
  }
  finite_min <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    min(x)
  }
  conf <- assoc |>
    dplyr::group_by(term) |>
    dplyr::summarise(max_eta2_pc1_pc10 = finite_max(eta2),
                     min_q_value = finite_min(q_value), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_eta2_pc1_pc10))
  qc_write_csv(conf, file.path(PATHS$tables, "PCA_confounding_summary.csv"))
}

pc_scores <- as.matrix(scores[, pc_terms[seq_len(min(2L, length(pc_terms)))], drop = FALSE])
cors <- dplyr::bind_rows(lapply(colnames(pc_scores), function(pc) {
  y <- pc_scores[, pc]
  out <- t(apply(mat, 1L, function(x) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(c(r = NA_real_, p = NA_real_))
    test <- suppressWarnings(stats::cor.test(x[ok], y[ok]))
    c(r = unname(test$estimate), p = test$p.value)
  }))
  data.frame(Protein = rownames(out), PC = pc, r = out[, "r"], p_value = out[, "p"], row.names = NULL)
}))
cors$q_value <- ave(cors$p_value, cors$PC, FUN = function(x) p.adjust(x, "BH"))
qc_write_csv(cors, file.path(PATHS$tables, "protein_pc_correlations.csv"))

for (pc in intersect(c("PC1", "PC2"), names(loadings))) {
  top <- rbind(
    head(loadings[order(-loadings[[pc]]), c("Protein", pc)], 50),
    head(loadings[order(loadings[[pc]]), c("Protein", pc)], 50)
  )
  qc_write_csv(top, file.path(PATHS$tables, paste0("top_loadings_", pc, ".csv")))
}

if (run$run_embeddings) {
  message("UMAP/t-SNE are exploratory only and should not be used as primary evidence in small or structured spatial proteomics datasets.")
  embedding_input <- pca$x[, seq_len(min(20L, ncol(pca$x))), drop = FALSE]
  if (requireNamespace("uwot", quietly = TRUE)) {
    set.seed(20260623)
    emb <- uwot::umap(embedding_input, n_threads = 1)
    emb_df <- data.frame(Sample = rownames(pca$x), UMAP1 = emb[, 1], UMAP2 = emb[, 2]) |>
      dplyr::left_join(meta, by = "Sample")
    qc_write_csv(emb_df, file.path(PATHS$tables, "exploratory_umap_scores.csv"))
    plot_umap_by <- function(term) {
      if (!term %in% names(emb_df)) return(invisible(NULL))
      keep <- !is.na(emb_df[[term]]) & nzchar(as.character(emb_df[[term]]))
      if (sum(keep) < 3L || length(unique(emb_df[[term]][keep])) < 2L) return(invisible(NULL))
      plot_df <- emb_df[keep, ]
      plot_df[[term]] <- factor(plot_df[[term]])
      p <- ggplot(plot_df, aes(UMAP1, UMAP2, color = .data[[term]])) +
        geom_point(size = 2, alpha = 0.85) +
        scale_color_manual(values = qc_embedding_palette(nlevels(plot_df[[term]]))) +
        labs(
          title = "Exploratory UMAP",
          subtitle = "Not primary evidence; computed from the first up to 20 PCs",
          x = "UMAP1",
          y = "UMAP2",
          color = term
        ) +
        coord_equal() +
        guides(color = guide_legend(nrow = qc_legend_rows(plot_df[[term]]), byrow = TRUE, override.aes = list(alpha = 1, size = 2))) +
        theme_embedding()
      qc_save_square_svg(file.path(PATHS$figures, paste0("umap_by_", safe_filename(term), ".svg")), p, size_mm = 90)
    }
    invisible(lapply(unique(c("Group", "group", "ExpGroup", "Region", "region", "Layer", "layer", "ReplicateGroup", "AnimalID", "plate", terms)), plot_umap_by))
  } else {
    warning("Package 'uwot' is not installed; skipped exploratory UMAP outputs.", call. = FALSE)
  }
  if (requireNamespace("Rtsne", quietly = TRUE)) {
    perplexity <- min(30L, floor((nrow(embedding_input) - 1L) / 3L))
    if (perplexity >= 1L) {
      set.seed(20260623)
      tsne <- Rtsne::Rtsne(embedding_input, dims = 2, perplexity = perplexity, pca = FALSE,
                           theta = 0.5, check_duplicates = FALSE, verbose = FALSE)
      tsne_df <- data.frame(Sample = rownames(pca$x), TSNE1 = tsne$Y[, 1], TSNE2 = tsne$Y[, 2],
                            tsne_perplexity = perplexity) |>
        dplyr::left_join(meta, by = "Sample")
      qc_write_csv(tsne_df, file.path(PATHS$tables, "exploratory_tsne_scores.csv"))
      plot_tsne_by <- function(term) {
        if (!term %in% names(tsne_df)) return(invisible(NULL))
        keep <- !is.na(tsne_df[[term]]) & nzchar(as.character(tsne_df[[term]]))
        if (sum(keep) < 3L || length(unique(tsne_df[[term]][keep])) < 2L) return(invisible(NULL))
        plot_df <- tsne_df[keep, ]
        plot_df[[term]] <- factor(plot_df[[term]])
        p <- ggplot(plot_df, aes(TSNE1, TSNE2, color = .data[[term]])) +
          geom_point(size = 2, alpha = 0.85) +
          scale_color_manual(values = qc_embedding_palette(nlevels(plot_df[[term]]))) +
          labs(
            title = "Exploratory t-SNE",
            subtitle = sprintf("Not primary evidence; computed from the first up to 20 PCs, perplexity %d", perplexity),
            x = "t-SNE1",
            y = "t-SNE2",
            color = term
          ) +
          coord_equal() +
          guides(color = guide_legend(nrow = qc_legend_rows(plot_df[[term]]), byrow = TRUE, override.aes = list(alpha = 1, size = 2))) +
          theme_embedding()
        qc_save_square_svg(file.path(PATHS$figures, paste0("tsne_by_", safe_filename(term), ".svg")), p, size_mm = 90)
      }
      invisible(lapply(unique(c("Group", "group", "ExpGroup", "Region", "region", "Layer", "layer", "ReplicateGroup", "AnimalID", "plate", terms)), plot_tsne_by))
    } else {
      warning("Too few samples for exploratory t-SNE; skipped t-SNE outputs.", call. = FALSE)
    }
  } else {
    warning("Package 'Rtsne' is not installed; skipped exploratory t-SNE outputs.", call. = FALSE)
  }
}

if (run$run_clustering) {
  message("Clustering is exploratory only and should not be used as primary evidence in small or structured spatial proteomics datasets.")
  pcs <- scale(pca$x[, seq_len(min(10L, ncol(pca$x))), drop = FALSE])
  km <- stats::kmeans(pcs, centers = min(3L, nrow(pcs) - 1L), nstart = 50)
  qc_write_csv(data.frame(Sample = rownames(pcs), exploratory_cluster = km$cluster),
               file.path(PATHS$tables, "exploratory_kmeans_clusters.csv"))
}

writeLines(c(
  "# PCA QC Notes",
  "",
  "PCA uses features with <80% missingness and non-zero variance.",
  "Remaining missing values are median-imputed only for PCA calculation.",
  "UMAP, t-SNE, and clustering outputs, when explicitly enabled, are exploratory only and not primary evidence."
), file.path(PATHS$reports, "pca_qc_notes.md"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(matrix = matrix_file, metadata = metadata_file),
  outputs = list(figures = PATHS$figures, tables = PATHS$tables, reports = PATHS$reports),
  parameters = list(dataset = DATASET, run_embeddings = run$run_embeddings, run_clustering = run$run_clustering),
  notes = "Dataset-aware PCA with PC1-PC10 metadata association and confounding summary."
)

message("PCA QC complete for dataset: ", DATASET)
