#!/usr/bin/env Rscript
# Global joint-compartment QC consumer. Exploratory embeddings never replace PCA.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "joint_compartment_qc_utils.R"))

global_arg <- tolower(script_arg_value("--dataset", "all"))
if (!global_arg %in% c("all", "global")) stop("This is a global QC script; use --dataset all or global.", call. = FALSE)
runtime <- list(script = "03_qc_exploration/00b_joint_compartment_qc.r", stage = "qc_global", dataset = "global", args = commandArgs(trailingOnly = TRUE), dry_run = is_dry_run(), started_at = Sys.time())
processed_root <- Sys.getenv("PROTEOMICS_JOINT_QC_PROCESSED_DIR", unset = path_processed("01_preprocessing", "joint_compartment_qc", "global"))
bundle_file <- file.path(processed_root, "joint_compartment_qc_matrices.rds")
out_root <- list(tables = path_results("tables", "03_qc_exploration", "00b_joint_compartment_qc", "global"), figures = path_results("figures", "03_qc_exploration", "00b_joint_compartment_qc", "global"), reports = path_results("reports", "03_qc_exploration", "00b_joint_compartment_qc", "global"), logs = path_results("logs", "03_qc_exploration", "00b_joint_compartment_qc", "global"))
seed <- suppressWarnings(as.integer(Sys.getenv("PROTEOMICS_JOINT_QC_SEED", unset = "42"))); if (is.na(seed)) seed <- 42L

if (runtime$dry_run) {
  cat("[DRY-RUN] joint preprocessing bundle: ", normalizePath(bundle_file, winslash = "/", mustWork = FALSE), " [", if (file.exists(bundle_file)) "PASS" else "FAIL", "]\n", sep = "")
  cat("[DRY-RUN] tables: ", out_root$tables, "\n[DRY-RUN] figures: ", out_root$figures, "\n[DRY-RUN] reports: ", out_root$reports, "\n", sep = "")
  quit(status = if (file.exists(bundle_file)) 0 else 1, save = "no")
}
required <- c("ggplot2", "svglite", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
for (p in out_root) dir_create(p)
if (!file.exists(bundle_file)) stop("Joint preprocessing bundle is missing: ", bundle_file, ". Run 01_prepare_joint_protigy_input.r first.", call. = FALSE)
bundle <- readRDS(bundle_file)
if (!identical(bundle$contract_version, "joint_compartment_qc_v1")) stop("Unsupported joint QC input contract.", call. = FALSE)
mat <- bundle$primary$matrix; complete <- bundle$complete_case$matrix; broad <- bundle$broad_union$detected_binary; meta <- bundle$metadata
if (anyNA(mat) || any(!is.finite(mat))) stop("Primary joint matrix must be finite after preprocessing imputation.", call. = FALSE)
if (anyNA(complete) || any(!is.finite(complete))) stop("Complete-case sensitivity matrix must be complete and finite.", call. = FALSE)
if (!identical(colnames(mat), as.character(meta$Sample))) stop("Joint QC metadata alignment failed.", call. = FALSE)

run_pca <- function(x) stats::prcomp(t(x), center = TRUE, scale. = TRUE)
pca <- run_pca(mat); complete_pca <- run_pca(complete)
# `meta` already carries Sample for matrix alignment. Keep one authoritative
# sample column in scores/embeddings so ggplot receives uniquely named data.
meta_for_outputs <- meta[, setdiff(names(meta), "Sample"), drop = FALSE]
npcs <- min(10L, ncol(pca$x)); scores <- data.frame(Sample = rownames(pca$x), pca$x[, seq_len(npcs), drop = FALSE], meta_for_outputs, check.names = FALSE)
loadings <- data.frame(ProteinGroupID = rownames(pca$rotation), pca$rotation[, seq_len(npcs), drop = FALSE], check.names = FALSE)
variance <- data.frame(PC = paste0("PC", seq_along(pca$sdev)), variance_explained = pca$sdev^2 / sum(pca$sdev^2), stringsAsFactors = FALSE)

metadata_term <- function(candidates, label) {
  col <- joint_qc_first_col(meta, candidates)
  if (is.na(col)) return(NULL)
  val <- as.character(meta[[col]]); val[is.na(val) | !nzchar(val)] <- NA_character_
  if (sum(!is.na(val)) < 3L || length(unique(val[!is.na(val)])) < 2L) return(NULL)
  list(label = label, value = val)
}
region_col <- joint_qc_first_col(meta, c("region", "Region")); layer_col <- joint_qc_first_col(meta, c("layer", "Layer"))
terms <- Filter(Negate(is.null), list(metadata_term(c("dataset"), "dataset"), metadata_term(c("compartment", "celltype_layer", "CellTypeLayer"), "compartment"), metadata_term(c("batch", "Batch"), "batch"), metadata_term(c("plate", "Plate"), "plate"), metadata_term(c("run", "Run", "run_order"), "run"), metadata_term(c("region", "Region"), "region"), metadata_term(c("layer", "Layer"), "layer"), metadata_term(c("ExpGroup", "group", "Group"), "experimental_group"), metadata_term(c("AnimalID", "animal_id"), "animal"), metadata_term(c("hemisphere", "Hemisphere"), "hemisphere"), metadata_term(c("sex", "Sex"), "sex"), metadata_term(c("exclusion_status"), "exclusion_status")))
if (!is.na(region_col) && !is.na(layer_col)) terms[[length(terms) + 1L]] <- list(label = "region_x_layer", value = interaction(meta[[region_col]], meta[[layer_col]], drop = TRUE))
assoc_one <- function(y, term) {
  keep <- is.finite(y) & !is.na(term$value)
  if (sum(keep) < 4L || length(unique(term$value[keep])) < 2L) return(data.frame())
  fit <- stats::lm(y[keep] ~ as.factor(term$value[keep])); an <- stats::anova(fit); total <- sum((y[keep] - mean(y[keep]))^2)
  data.frame(metadata_term = term$label, n = sum(keep), n_levels = length(unique(term$value[keep])), eta_squared = if (total > 0) an$`Sum Sq`[1] / total else NA_real_, p_value = an$`Pr(>F)`[1], stringsAsFactors = FALSE)
}
association_parts <- lapply(seq_len(npcs), function(i) {
  pieces <- Filter(function(x) nrow(x) > 0L, lapply(terms, function(t) assoc_one(scores[[paste0("PC", i)]], t)))
  if (!length(pieces)) return(NULL)
  z <- do.call(rbind, pieces); z$PC <- paste0("PC", i); z
})
association_parts <- Filter(Negate(is.null), association_parts)
association <- if (length(association_parts)) do.call(rbind, association_parts) else data.frame(PC = character(), metadata_term = character(), n = integer(), n_levels = integer(), eta_squared = numeric(), p_value = numeric())
if (nrow(association)) association$FDR <- stats::p.adjust(association$p_value, method = "BH")

pc_input <- pca$x[, seq_len(min(20L, ncol(pca$x))), drop = FALSE]
embedding_status <- data.frame(method = c("UMAP", "tSNE"), status = "not_run", reason = "package unavailable or insufficient samples", stringsAsFactors = FALSE)
if (nrow(pc_input) >= 4L && requireNamespace("uwot", quietly = TRUE)) {
  set.seed(seed); umap <- uwot::umap(pc_input, n_neighbors = min(15L, nrow(pc_input) - 1L), n_threads = 1, verbose = FALSE, ret_model = FALSE)
  umap_scores <- cbind(data.frame(Sample = rownames(pc_input), UMAP1 = umap[, 1], UMAP2 = umap[, 2]), meta_for_outputs)
  utils::write.csv(umap_scores, file.path(out_root$tables, "joint_umap_exploratory_scores.csv"), row.names = FALSE); embedding_status[1, 2:3] <- c("completed", "fixed seed; first up to 20 primary PCs")
} else umap_scores <- NULL
if (nrow(pc_input) >= 4L && requireNamespace("Rtsne", quietly = TRUE)) {
  perplexity <- min(30L, floor((nrow(pc_input) - 1L) / 3L))
  if (perplexity >= 1L) { set.seed(seed); tsne <- Rtsne::Rtsne(pc_input, dims = 2, perplexity = perplexity, check_duplicates = FALSE, pca = FALSE, verbose = FALSE, max_iter = 1000); tsne_scores <- cbind(data.frame(Sample = rownames(pc_input), tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2]), meta_for_outputs); utils::write.csv(tsne_scores, file.path(out_root$tables, "joint_tsne_exploratory_scores.csv"), row.names = FALSE); embedding_status[2, 2:3] <- c("completed", paste0("fixed seed; perplexity=", perplexity)) } else tsne_scores <- NULL
} else tsne_scores <- NULL

corr <- stats::cor(mat, method = "pearson"); distance <- as.matrix(stats::dist(t(mat)))
nearest <- apply(distance + diag(Inf, nrow(distance)), 1, which.min)
nearest_table <- data.frame(Sample = rownames(distance), nearest_sample = rownames(distance)[nearest], dataset = meta$dataset, nearest_dataset = meta$dataset[nearest], nearest_compartment_consistent = meta$dataset == meta$dataset[nearest], nearest_distance = distance[cbind(seq_len(nrow(distance)), nearest)], leverage = rowSums(scale(pca$x[, seq_len(min(10L, ncol(pca$x))), drop = FALSE])^2), stringsAsFactors = FALSE)
outlier_cutoff <- stats::quantile(nearest_table$leverage, .99, na.rm = TRUE); nearest_table$flagged_leverage_outlier <- nearest_table$leverage > outlier_cutoff
cluster <- stats::hclust(stats::as.dist(1 - corr), method = "average")
primary_dist <- as.vector(stats::dist(t(mat))); complete_dist <- as.vector(stats::dist(t(complete)))
concordance <- data.frame(metric = c("pairwise_distance_correlation", "PC1_score_correlation"), value = c(stats::cor(primary_dist, complete_dist), abs(stats::cor(pca$x[, 1], complete_pca$x[, 1]))), stringsAsFactors = FALSE)

write.csv(scores, file.path(out_root$tables, "joint_primary_pca_scores.csv"), row.names = FALSE); write.csv(loadings, file.path(out_root$tables, "joint_primary_pca_loadings.csv"), row.names = FALSE); write.csv(variance, file.path(out_root$tables, "joint_primary_pca_variance_explained.csv"), row.names = FALSE); write.csv(association, file.path(out_root$tables, "joint_pc_metadata_associations.csv"), row.names = FALSE); write.csv(nearest_table, file.path(out_root$tables, "joint_sample_nearest_neighbour_and_outlier_metrics.csv"), row.names = FALSE); write.csv(as.data.frame(corr), file.path(out_root$tables, "joint_sample_correlations.csv"), row.names = TRUE); write.csv(data.frame(Sample = cluster$labels[cluster$order], cluster_order = seq_along(cluster$order)), file.path(out_root$tables, "joint_hierarchical_clustering_order.csv"), row.names = FALSE); write.csv(concordance, file.path(out_root$tables, "primary_vs_complete_case_pca_concordance.csv"), row.names = FALSE); write.csv(embedding_status, file.path(out_root$tables, "embedding_status.csv"), row.names = FALSE)

theme_qc <- ggplot2::theme_classic(base_size = 8) + ggplot2::theme(legend.position = "bottom")
dataset_palette <- c(neuron_neuropil = "#0072B2", neuron_soma = "#D55E00", microglia = "#009E73")
save_plot <- function(name, plot, width = 100, height = 85) ggplot2::ggsave(file.path(out_root$figures, name), plot, width = width, height = height, units = "mm", device = svglite::svglite)
score2 <- scores[, c("PC1", "PC2", "dataset"), drop = FALSE]
save_plot("joint_pca_by_dataset.svg", ggplot2::ggplot(score2, ggplot2::aes(PC1, PC2, color = dataset)) + ggplot2::geom_point(size = 1.8) + ggplot2::scale_color_manual(values = dataset_palette) + ggplot2::labs(title = "Joint PCA: biological dataset", subtitle = "Raw-derived primary shared core; no batch correction") + theme_qc)
block_var <- bundle$technical_block$variable; block_values <- bundle$technical_block$block
score2$technical_block <- block_values; save_plot("joint_pca_dataset_by_technical_block.svg", ggplot2::ggplot(score2, ggplot2::aes(PC1, PC2, color = dataset, shape = technical_block)) + ggplot2::geom_point(size = 1.8) + ggplot2::scale_color_manual(values = dataset_palette) + ggplot2::labs(title = "Joint PCA: dataset and technical block", shape = block_var) + theme_qc)
for (item in list(list("region_layer", if (!is.na(region_col) && !is.na(layer_col)) interaction(meta[[region_col]], meta[[layer_col]], drop = TRUE) else NULL), list("experimental_group", if (!is.na(joint_qc_first_col(meta, c("ExpGroup", "group", "Group")))) meta[[joint_qc_first_col(meta, c("ExpGroup", "group", "Group"))]] else NULL))) { if (!is.null(item[[2]])) { score2$term <- item[[2]]; save_plot(paste0("joint_pca_by_", item[[1]], ".svg"), ggplot2::ggplot(score2, ggplot2::aes(PC1, PC2, color = term)) + ggplot2::geom_point(size = 1.7) + ggplot2::labs(title = paste("Joint PCA:", gsub("_", " ", item[[1]])), color = gsub("_", " ", item[[1]])) + theme_qc, 120, 90) } }
save_plot("joint_pca_scree.svg", ggplot2::ggplot(variance[seq_len(min(20, nrow(variance))), ], ggplot2::aes(PC, variance_explained)) + ggplot2::geom_col(fill = "#0072B2") + ggplot2::scale_y_continuous(labels = scales::percent) + ggplot2::labs(title = "Joint PCA variance explained", x = NULL, y = "Variance explained") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
if (!is.null(umap_scores)) save_plot("joint_umap_by_dataset_exploratory.svg", ggplot2::ggplot(umap_scores, ggplot2::aes(UMAP1, UMAP2, color = dataset)) + ggplot2::geom_point(size = 1.8) + ggplot2::scale_color_manual(values = dataset_palette) + ggplot2::labs(title = "Exploratory joint UMAP") + theme_qc)
if (!is.null(tsne_scores)) save_plot("joint_tsne_by_dataset_exploratory.svg", ggplot2::ggplot(tsne_scores, ggplot2::aes(tSNE1, tSNE2, color = dataset)) + ggplot2::geom_point(size = 1.8) + ggplot2::scale_color_manual(values = dataset_palette) + ggplot2::labs(title = "Exploratory joint t-SNE") + theme_qc)
corr_long <- as.data.frame(as.table(corr)); names(corr_long) <- c("sample_x", "sample_y", "correlation"); save_plot("joint_sample_correlation_heatmap.svg", ggplot2::ggplot(corr_long, ggplot2::aes(sample_x, sample_y, fill = correlation)) + ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", limits = c(-1, 1)) + ggplot2::labs(title = "Joint sample correlations", x = NULL, y = NULL) + theme_qc + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()), 115, 105)
missing_long <- as.data.frame(as.table(broad[seq_len(min(200L, nrow(broad))), , drop = FALSE])); names(missing_long) <- c("ProteinGroupID", "Sample", "detected"); save_plot("joint_broad_union_missingness_heatmap.svg", ggplot2::ggplot(missing_long, ggplot2::aes(Sample, ProteinGroupID, fill = factor(detected))) + ggplot2::geom_tile() + ggplot2::scale_fill_manual(values = c("0" = "white", "1" = "#2F6F73")) + ggplot2::labs(title = "Broad union detection/missingness (first 200 proteins)", fill = "Detected") + theme_qc + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()), 115, 120)
norm <- utils::read.csv(file.path(processed_root, "sample_normalization_audit.csv")); norm_long <- rbind(data.frame(Sample = norm$Sample, dataset = norm$dataset, stage = "pre", median = norm$pre_normalization_median_log2), data.frame(Sample = norm$Sample, dataset = norm$dataset, stage = "post", median = norm$post_normalization_median_log2)); save_plot("joint_pre_post_normalization_medians.svg", ggplot2::ggplot(norm_long, ggplot2::aes(dataset, median, color = stage)) + ggplot2::geom_jitter(width = .15) + ggplot2::labs(title = "Joint normalization: sample medians", y = "log2 intensity median") + theme_qc)
distribution <- rbind(data.frame(value = as.vector(bundle$primary$unnormalized_log2), stage = "pre-normalization"), data.frame(value = as.vector(mat), stage = "post-normalization")); distribution <- distribution[is.finite(distribution$value), , drop = FALSE]; save_plot("joint_pre_post_normalization_abundance_distributions.svg", ggplot2::ggplot(distribution, ggplot2::aes(value, color = stage)) + ggplot2::geom_density() + ggplot2::labs(title = "Joint abundance distributions", x = "log2 intensity", y = "Density") + theme_qc)
imp <- utils::read.csv(file.path(processed_root, "imputation_footprint_by_sample.csv")); save_plot("joint_imputation_footprint_by_sample.svg", ggplot2::ggplot(imp, ggplot2::aes(Sample, fraction_imputed, color = dataset)) + ggplot2::geom_point() + ggplot2::scale_color_manual(values = dataset_palette) + ggplot2::labs(title = "Primary imputation footprint", y = "Fraction imputed", x = NULL) + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()))
save_plot("joint_primary_complete_case_concordance.svg", ggplot2::ggplot(data.frame(x = primary_dist, y = complete_dist), ggplot2::aes(x, y)) + ggplot2::geom_point(size = .4, alpha = .4) + ggplot2::labs(title = "Primary vs complete-case distance concordance", x = "Primary distance", y = "Complete-case distance") + theme_qc)

assoc_summary <- if (nrow(association)) aggregate(eta_squared ~ metadata_term, association, max) else data.frame(metadata_term = character(), eta_squared = numeric())
writeLines(c("# Joint compartment QC summary", paste0("Primary matrix: ", nrow(mat), " proteins × ", ncol(mat), " retained samples."), paste0("Technical-block variable: ", bundle$technical_block$variable, " (", bundle$technical_block$mode, ")."), paste0("Complete-case proteins: ", nrow(complete), "; broad-union proteins: ", nrow(broad), "."), paste0("Primary imputation fraction: ", signif(mean(utils::read.csv(file.path(processed_root, "imputation_footprint_by_sample.csv"))$fraction_imputed), 4), "."), paste0("Flagged leverage samples: ", paste(nearest_table$Sample[nearest_table$flagged_leverage_outlier], collapse = ", ")), paste0("Primary/complete distance concordance: ", signif(concordance$value[1], 4), "."), "PCA is uncorrected and primary. UMAP/t-SNE are exploratory only. Dataset separation is expected biology; association tables quantify, but do not prove, technical or biological causes. Microglia observations are microglia-enriched ROIs, not purified microglia."), file.path(out_root$reports, "joint_compartment_qc_summary.md"))
write_run_manifest(file.path(out_root$logs, "run_manifest.yml"), inputs = list(joint_preprocessing_bundle = bundle_file), outputs = out_root, parameters = list(seed = seed, primary_pca = "centered_scaled_proteins", embeddings = "exploratory; first up to 20 PCs; no batch correction"), notes = "Global QC and biological-identity analysis only; no combined DE/WGCNA.")
