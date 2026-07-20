#!/usr/bin/env Rscript
# Plot-only Nature-style rendering of the completed joint-compartment QC.
# PCA/UMAP/t-SNE coordinates and all statistical source tables are consumed
# unchanged from 00b_joint_compartment_qc.r.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "joint_compartment_qc_plotting.R"))

global_arg <- tolower(script_arg_value("--dataset", "global"))
if (!global_arg %in% c("all", "global")) stop("This is a global figure script; use --dataset all or global.", call. = FALSE)
runtime <- list(
  script = "03_qc_exploration/00c_joint_compartment_qc_publication_figures.r",
  stage = "qc_global_publication_figures",
  dataset = "global",
  args = commandArgs(trailingOnly = TRUE),
  dry_run = is_dry_run(),
  started_at = Sys.time()
)

processed_root <- Sys.getenv(
  "PROTEOMICS_JOINT_QC_PROCESSED_DIR",
  unset = path_processed("01_preprocessing", "joint_compartment_qc", "global")
)
qc_table_root <- path_results("tables", "03_qc_exploration", "00b_joint_compartment_qc", "global")
figure_root <- Sys.getenv(
  "PROTEOMICS_JOINT_QC_PUBLICATION_FIGURE_DIR",
  unset = path_results("figures", "03_qc_exploration", "00b_joint_compartment_qc", "publication_style", "global")
)
panel_root <- file.path(figure_root, "panels")
report_root <- path_results("reports", "03_qc_exploration", "00c_joint_compartment_qc_publication_figures", "global")
log_root <- path_results("logs", "03_qc_exploration", "00c_joint_compartment_qc_publication_figures", "global")

inputs <- c(
  bundle = file.path(processed_root, "joint_compartment_qc_matrices.rds"),
  normalization = file.path(processed_root, "sample_normalization_audit.csv"),
  imputation = file.path(processed_root, "imputation_footprint_by_sample.csv"),
  pca_scores = file.path(qc_table_root, "joint_primary_pca_scores.csv"),
  pca_variance = file.path(qc_table_root, "joint_primary_pca_variance_explained.csv"),
  umap_scores = file.path(qc_table_root, "joint_umap_exploratory_scores.csv"),
  tsne_scores = file.path(qc_table_root, "joint_tsne_exploratory_scores.csv"),
  correlations = file.path(qc_table_root, "joint_sample_correlations.csv"),
  cluster_order = file.path(qc_table_root, "joint_hierarchical_clustering_order.csv"),
  concordance = file.path(qc_table_root, "primary_vs_complete_case_pca_concordance.csv")
)

if (runtime$dry_run) {
  for (nm in names(inputs)) {
    cat("[DRY-RUN] ", nm, ": ", normalizePath(inputs[[nm]], winslash = "/", mustWork = FALSE), " [", if (file.exists(inputs[[nm]])) "PASS" else "FAIL", "]\n", sep = "")
  }
  cat("[DRY-RUN] publication figures: ", normalizePath(figure_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
  quit(status = if (all(file.exists(inputs))) 0 else 1, save = "no")
}

missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing completed joint-QC input(s): ", paste(names(missing_inputs), collapse = ", "), ". Run 00b_joint_compartment_qc.r first.", call. = FALSE)
}
required_packages <- c("ggplot2", "svglite", "scales", "patchwork", "cowplot", "pheatmap", "ggrepel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
dir_create(figure_root)
dir_create(panel_root)
dir_create(report_root)
dir_create(log_root)

read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

bundle <- readRDS(inputs[["bundle"]])
if (!identical(bundle$contract_version, "joint_compartment_qc_v1")) stop("Unsupported joint-QC bundle contract.", call. = FALSE)
mat <- bundle$primary$matrix
complete <- bundle$complete_case$matrix
broad <- bundle$broad_union$detected_binary
meta <- bundle$metadata
if (!identical(colnames(mat), as.character(meta$Sample)) || !identical(colnames(broad), as.character(meta$Sample))) {
  stop("Publication plotting inputs are not aligned to bundle metadata.", call. = FALSE)
}

scores <- read_csv(inputs[["pca_scores"]])
variance <- joint_pub_order_pcs(read_csv(inputs[["pca_variance"]]))
umap_scores <- read_csv(inputs[["umap_scores"]])
tsne_scores <- read_csv(inputs[["tsne_scores"]])
normalization <- read_csv(inputs[["normalization"]])
imputation <- read_csv(inputs[["imputation"]])
concordance <- read_csv(inputs[["concordance"]])
require_columns(scores, c("Sample", "PC1", "PC2", "dataset", "plate", "region", "layer", "ExpGroup"), "PCA score table")
require_columns(umap_scores, c("Sample", "UMAP1", "UMAP2", "dataset"), "UMAP score table")
require_columns(tsne_scores, c("Sample", "tSNE1", "tSNE2", "dataset"), "t-SNE score table")
require_columns(normalization, c("Sample", "dataset", "pre_normalization_median_log2", "post_normalization_median_log2"), "Normalization audit")
require_columns(imputation, c("Sample", "dataset", "fraction_imputed"), "Imputation audit")
if (!identical(as.character(scores$Sample), as.character(meta$Sample))) stop("PCA score order differs from the completed bundle.", call. = FALSE)
joint_pub_assert_coordinates(scores, umap_scores, c("dataset"), "UMAP source validation")
joint_pub_assert_coordinates(scores, tsne_scores, c("dataset"), "t-SNE source validation")

dataset_palette <- joint_pub_dataset_palette()
dataset_shapes <- c("Microglia-enriched ROI" = 16, "Neuropil" = 17, "Soma" = 15)
region_palette <- c(CA1 = "#4477AA", CA2 = "#EE6677", CA3 = "#228833", DG = "#CCBB44")
plate_palette <- c(B = "#4C78A8", C = "#D89C4A")
group_palette <- c(CON = "#4C78A8", RES = "#D6A74A", SUS = "#B35C66")
layer_shapes <- c(ROI = 16, SLM = 17, SO = 15, SP = 18, SR = 3, MO = 7, PO = 8, SG = 4)
base_theme <- joint_pub_theme()

dataset_factor <- joint_pub_dataset_factor(scores$dataset)
pc1_variance <- variance$variance_explained[variance$pc_number == 1L]
pc2_variance <- variance$variance_explained[variance$pc_number == 2L]
if (length(pc1_variance) != 1L || length(pc2_variance) != 1L) stop("PC1/PC2 variance values are unavailable.", call. = FALSE)
pc_x_label <- paste0("PC1 (", sprintf("%.1f", 100 * pc1_variance), "%)")
pc_y_label <- paste0("PC2 (", sprintf("%.1f", 100 * pc2_variance), "%)")
pc_limits <- joint_pub_equal_limits(scores$PC1, scores$PC2)

equal_coordinate_scales <- function(limits) {
  list(
    ggplot2::scale_x_continuous(limits = limits$x, expand = ggplot2::expansion(mult = 0)),
    ggplot2::scale_y_continuous(limits = limits$y, expand = ggplot2::expansion(mult = 0)),
    ggplot2::coord_fixed(ratio = 1, clip = "off")
  )
}
dataset_guide <- ggplot2::guide_legend(
  title.position = "top", title.hjust = 0, nrow = 1, byrow = TRUE, order = 1,
  override.aes = list(size = 1.5, alpha = 1)
)
plate_guide <- ggplot2::guide_legend(
  title.position = "top", title.hjust = 0, nrow = 1, byrow = TRUE, order = 2,
  override.aes = list(size = 1.5, alpha = 1)
)

pca_plot_data <- data.frame(
  Sample = scores$Sample,
  PC1 = scores$PC1,
  PC2 = scores$PC2,
  Dataset = dataset_factor,
  Plate = factor(as.character(scores$plate)),
  Region = factor(toupper(as.character(scores$region)), levels = names(region_palette)),
  Layer = toupper(as.character(scores$layer)),
  Group = joint_pub_group_factor(scores$ExpGroup),
  stringsAsFactors = FALSE
)
pca_plot_data$Layer[pca_plot_data$Layer == "MICROGLIA"] <- "ROI"
pca_plot_data$Layer <- factor(pca_plot_data$Layer, levels = names(layer_shapes))

p_pca_dataset <- ggplot2::ggplot(pca_plot_data, ggplot2::aes(PC1, PC2, colour = Dataset)) +
  ggplot2::geom_point(size = 0.72, alpha = 0.72) +
  ggplot2::scale_colour_manual(values = dataset_palette, drop = FALSE) +
  equal_coordinate_scales(pc_limits) +
  ggplot2::labs(x = pc_x_label, y = pc_y_label, colour = "Dataset") +
  ggplot2::guides(colour = dataset_guide) +
  base_theme

p_pca_technical <- ggplot2::ggplot(pca_plot_data, ggplot2::aes(PC1, PC2, colour = Dataset, shape = Plate)) +
  ggplot2::geom_point(size = 0.74, alpha = 0.72) +
  ggplot2::scale_colour_manual(values = dataset_palette, drop = FALSE) +
  ggplot2::scale_shape_manual(values = c(B = 16, C = 17), drop = FALSE) +
  equal_coordinate_scales(pc_limits) +
  ggplot2::labs(x = pc_x_label, y = pc_y_label, colour = "Dataset", shape = "Plate") +
  ggplot2::guides(colour = dataset_guide, shape = plate_guide) +
  base_theme +
  ggplot2::theme(legend.box = "vertical", legend.box.just = "left")

p_pca_anatomy <- ggplot2::ggplot(pca_plot_data, ggplot2::aes(PC1, PC2, colour = Region, shape = Layer)) +
  ggplot2::geom_point(size = 0.62, alpha = 0.68) +
  ggplot2::facet_wrap(~Dataset, nrow = 1, drop = FALSE) +
  ggplot2::scale_colour_manual(values = region_palette, drop = FALSE, na.translate = FALSE) +
  ggplot2::scale_shape_manual(values = layer_shapes, drop = TRUE, na.translate = FALSE) +
  equal_coordinate_scales(pc_limits) +
  ggplot2::labs(x = pc_x_label, y = pc_y_label, colour = "Region", shape = "Layer") +
  ggplot2::guides(
    colour = ggplot2::guide_legend(title.position = "top", title.hjust = 0, nrow = 1, order = 1, override.aes = list(size = 1.4, alpha = 1)),
    shape = ggplot2::guide_legend(title.position = "top", title.hjust = 0, nrow = 2, byrow = TRUE, order = 2, override.aes = list(size = 1.4, alpha = 1))
  ) +
  base_theme +
  ggplot2::theme(legend.box = "vertical", legend.box.just = "left")

p_pca_group <- ggplot2::ggplot(pca_plot_data, ggplot2::aes(PC1, PC2, colour = Group, shape = Dataset)) +
  ggplot2::geom_point(size = 0.68, alpha = 0.68) +
  ggplot2::scale_colour_manual(values = group_palette, drop = FALSE) +
  ggplot2::scale_shape_manual(values = dataset_shapes, drop = FALSE) +
  equal_coordinate_scales(pc_limits) +
  ggplot2::labs(x = pc_x_label, y = pc_y_label, colour = "Experimental group", shape = "Dataset") +
  ggplot2::guides(
    colour = ggplot2::guide_legend(title.position = "top", title.hjust = 0, nrow = 1, order = 1, override.aes = list(size = 1.4, alpha = 1)),
    shape = ggplot2::guide_legend(title.position = "top", title.hjust = 0, nrow = 2, order = 2, override.aes = list(size = 1.4, alpha = 1))
  ) +
  base_theme

scree_data <- variance[variance$pc_number <= 20L, , drop = FALSE]
p_scree <- ggplot2::ggplot(scree_data, ggplot2::aes(PC, variance_explained)) +
  ggplot2::geom_col(width = 0.76, fill = unname(dataset_palette[["Neuropil"]])) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(x = NULL, y = "Variance explained") +
  base_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6.2), legend.position = "none")

primary_dist <- as.vector(stats::dist(t(mat)))
complete_dist <- as.vector(stats::dist(t(complete)))
distance_r <- stats::cor(primary_dist, complete_dist)
reported_r <- concordance$value[match("pairwise_distance_correlation", concordance$metric)]
if (length(reported_r) != 1L || !isTRUE(all.equal(distance_r, reported_r, tolerance = 1e-12))) {
  stop("Distance concordance no longer matches the completed QC source table.", call. = FALSE)
}
distance_limits <- joint_pub_equal_limits(primary_dist, complete_dist, padding = 0.025)
distance_data <- data.frame(primary = primary_dist, complete_case = complete_dist)
p_concordance <- ggplot2::ggplot(distance_data, ggplot2::aes(primary, complete_case)) +
  ggplot2::geom_point(size = 0.14, alpha = 0.055, colour = "#176a5a") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "22", colour = "#3f3f3f") +
  ggplot2::annotate(
    "text",
    x = distance_limits$x[[1]] + 0.04 * distance_limits$span,
    y = distance_limits$y[[2]] - 0.04 * distance_limits$span,
    label = paste0("Pearson r = ", sprintf("%.3f", distance_r)),
    hjust = 0, vjust = 1, family = "sans", size = 2.45, colour = "#1a1a1a"
  ) +
  equal_coordinate_scales(distance_limits) +
  ggplot2::labs(x = "Primary-matrix sample distance", y = "Complete-case sample distance") +
  base_theme +
  ggplot2::theme(legend.position = "none")

normalization$Dataset <- joint_pub_dataset_factor(normalization$dataset)
norm_long <- rbind(
  data.frame(Sample = normalization$Sample, Dataset = normalization$Dataset, Stage = "Pre", median = normalization$pre_normalization_median_log2),
  data.frame(Sample = normalization$Sample, Dataset = normalization$Dataset, Stage = "Post", median = normalization$post_normalization_median_log2)
)
norm_long$Stage <- factor(norm_long$Stage, levels = c("Pre", "Post"))
p_norm_medians <- ggplot2::ggplot(norm_long, ggplot2::aes(Stage, median, group = Sample)) +
  ggplot2::geom_line(linewidth = 0.18, alpha = 0.22, colour = "#9a9a9a") +
  ggplot2::geom_point(ggplot2::aes(colour = Dataset), size = 0.42, alpha = 0.58) +
  ggplot2::scale_colour_manual(values = dataset_palette, drop = FALSE) +
  ggplot2::labs(x = NULL, y = "Sample median (log2)") +
  base_theme +
  ggplot2::theme(legend.position = "none")

density_data <- rbind(
  data.frame(intensity = as.vector(bundle$primary$unnormalized_log2), Stage = "Pre"),
  data.frame(intensity = as.vector(mat), Stage = "Post")
)
density_data <- density_data[is.finite(density_data$intensity), , drop = FALSE]
density_data$Stage <- factor(density_data$Stage, levels = c("Pre", "Post"))
p_norm_density <- ggplot2::ggplot(density_data, ggplot2::aes(intensity, colour = Stage, linetype = Stage)) +
  ggplot2::geom_density(linewidth = 0.46, adjust = 1, show.legend = TRUE) +
  ggplot2::scale_colour_manual(values = c(Pre = "#767577", Post = "#176a5a")) +
  ggplot2::scale_linetype_manual(values = c(Pre = "solid", Post = "22")) +
  ggplot2::labs(x = "Protein intensity (log2)", y = "Density", colour = "", linetype = "") +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 1, override.aes = list(linewidth = 0.7)),
    linetype = "none"
  ) +
  base_theme
p_normalization <- cowplot::plot_grid(p_norm_medians, p_norm_density, nrow = 1, rel_widths = c(0.95, 1.05), align = "h", axis = "tb")

imp_idx <- match(imputation$Sample, meta$Sample)
if (anyNA(imp_idx) || nrow(imputation) != nrow(meta)) stop("Imputation audit and metadata sample sets differ.", call. = FALSE)
imputation$Dataset <- joint_pub_dataset_factor(imputation$dataset)
imputation$AnimalID <- as.character(meta$AnimalID[imp_idx])
imputation$Region <- toupper(as.character(meta$region[imp_idx]))
imputation$Layer <- toupper(as.character(meta$layer[imp_idx]))
imputation$sample_number <- as.character(meta$sample_number[imp_idx])
imp_order <- order(-imputation$fraction_imputed, as.character(imputation$Sample), method = "radix")
imputation <- imputation[imp_order, , drop = FALSE]
imputation$rank <- seq_len(nrow(imputation))
imputation$label <- paste(imputation$AnimalID, paste(imputation$Region, imputation$Layer), imputation$sample_number, sep = " · ")
label_data <- head(imputation, 5L)
p_imputation <- ggplot2::ggplot(imputation, ggplot2::aes(rank, fraction_imputed, colour = Dataset)) +
  ggplot2::geom_hline(yintercept = c(0.05, 0.10), linewidth = 0.35, linetype = "22", colour = "#767577") +
  ggplot2::geom_point(size = 0.62, alpha = 0.74) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = label),
    family = "sans", size = 2.1, colour = "#1a1a1a",
    min.segment.length = 0, segment.size = 0.22, seed = 42,
    box.padding = 0.22, point.padding = 0.15, max.overlaps = Inf,
    direction = "both", hjust = 0
  ) +
  ggplot2::annotate("text", x = nrow(imputation), y = 0.05, label = "5%", hjust = 1, vjust = -0.35, family = "sans", size = 2.1, colour = "#5f5f5f") +
  ggplot2::annotate("text", x = nrow(imputation), y = 0.10, label = "10%", hjust = 1, vjust = -0.35, family = "sans", size = 2.1, colour = "#5f5f5f") +
  ggplot2::scale_colour_manual(values = dataset_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
  ggplot2::labs(x = "Sample rank (highest imputation first)", y = "Imputed values", colour = "Dataset") +
  ggplot2::guides(colour = dataset_guide) +
  base_theme +
  ggplot2::theme(plot.margin = ggplot2::margin(2.2, 6, 2.2, 2.2))

correlation <- as.matrix(utils::read.csv(inputs[["correlations"]], row.names = 1, check.names = FALSE))
storage.mode(correlation) <- "numeric"
if (!identical(rownames(correlation), colnames(correlation)) || !setequal(rownames(correlation), meta$Sample)) {
  stop("Sample-correlation matrix identifiers are invalid.", call. = FALSE)
}
cor_hc <- stats::hclust(stats::as.dist(1 - correlation), method = "average")
reported_cluster <- read_csv(inputs[["cluster_order"]])
clustered_samples <- cor_hc$labels[cor_hc$order]
if (!identical(clustered_samples, as.character(reported_cluster$Sample[order(reported_cluster$cluster_order)]))) {
  stop("Correlation clustering differs from the completed QC cluster-order table.", call. = FALSE)
}
annotation_index <- match(colnames(correlation), meta$Sample)
annotation <- data.frame(
  Dataset = joint_pub_dataset_factor(meta$dataset[annotation_index]),
  Plate = factor(as.character(meta$plate[annotation_index])),
  Anatomy = factor(toupper(as.character(meta$region[annotation_index])), levels = names(region_palette)),
  row.names = colnames(correlation),
  check.names = FALSE
)
annotation_colours <- list(Dataset = dataset_palette, Plate = plate_palette, Anatomy = region_palette)
corr_floor <- floor(min(correlation) * 100) / 100
corr_breaks <- seq(corr_floor, 1, length.out = 101L)
corr_colours <- grDevices::colorRampPalette(c("#f7f7f5", "#a8d8d5", "#176a5a"))(100L)
corr_heatmap <- pheatmap::pheatmap(
  correlation,
  color = corr_colours,
  breaks = corr_breaks,
  cluster_rows = cor_hc,
  cluster_cols = cor_hc,
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  annotation_col = annotation,
  annotation_colors = annotation_colours,
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  legend_breaks = c(0.65, 0.8, 0.9, 1),
  legend_labels = sprintf("%.2f", c(0.65, 0.8, 0.9, 1)),
  treeheight_row = 8,
  treeheight_col = 8,
  fontsize = 6,
  fontsize_row = 5.5,
  fontsize_col = 5.5,
  fontsize_legend = 5.5,
  main = NA,
  silent = TRUE
)
p_correlation <- cowplot::ggdraw() + cowplot::draw_grob(corr_heatmap$gtable)

selection <- joint_pub_select_detection_features(broad, meta$dataset, n = 200L)
selected_detection <- broad[selection$ProteinGroupID, , drop = FALSE]
sample_order <- joint_pub_order_samples(meta)
selected_detection <- selected_detection[, sample_order, drop = FALSE]
protein_hc <- stats::hclust(stats::dist(selected_detection), method = "average")
missing_annotation <- data.frame(
  Dataset = joint_pub_dataset_factor(meta$dataset[sample_order]),
  Plate = factor(as.character(meta$plate[sample_order])),
  Anatomy = factor(toupper(as.character(meta$region[sample_order])), levels = names(region_palette)),
  row.names = meta$Sample[sample_order],
  check.names = FALSE
)
missing_heatmap <- pheatmap::pheatmap(
  selected_detection,
  color = c("#f4f4f2", "#176a5a"),
  breaks = c(-0.5, 0.5, 1.5),
  cluster_rows = protein_hc,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  annotation_col = missing_annotation,
  annotation_colors = annotation_colours,
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  legend_breaks = c(0, 1),
  legend_labels = c("Missing", "Detected"),
  treeheight_row = 8,
  treeheight_col = 0,
  fontsize = 6,
  fontsize_row = 5.5,
  fontsize_col = 5.5,
  fontsize_legend = 5.5,
  main = NA,
  silent = TRUE
)
p_missingness <- cowplot::ggdraw() + cowplot::draw_grob(missing_heatmap$gtable)

umap_limits <- joint_pub_equal_limits(umap_scores$UMAP1, umap_scores$UMAP2)
umap_data <- data.frame(UMAP1 = umap_scores$UMAP1, UMAP2 = umap_scores$UMAP2, Dataset = joint_pub_dataset_factor(umap_scores$dataset))
p_umap <- ggplot2::ggplot(umap_data, ggplot2::aes(UMAP1, UMAP2, colour = Dataset)) +
  ggplot2::geom_point(size = 0.58, alpha = 0.72) +
  ggplot2::scale_colour_manual(values = dataset_palette, drop = FALSE) +
  equal_coordinate_scales(umap_limits) +
  ggplot2::labs(x = "UMAP1", y = "UMAP2", colour = "Dataset") +
  ggplot2::guides(colour = dataset_guide) +
  base_theme

tsne_limits <- joint_pub_equal_limits(tsne_scores$tSNE1, tsne_scores$tSNE2)
tsne_data <- data.frame(tSNE1 = tsne_scores$tSNE1, tSNE2 = tsne_scores$tSNE2, Dataset = joint_pub_dataset_factor(tsne_scores$dataset))
p_tsne <- ggplot2::ggplot(tsne_data, ggplot2::aes(tSNE1, tSNE2, colour = Dataset)) +
  ggplot2::geom_point(size = 0.58, alpha = 0.72) +
  ggplot2::scale_colour_manual(values = dataset_palette, drop = FALSE) +
  equal_coordinate_scales(tsne_limits) +
  ggplot2::labs(x = "t-SNE1", y = "t-SNE2", colour = "Dataset") +
  ggplot2::guides(colour = dataset_guide) +
  base_theme

# Editable title-free manuscript panels. Panel letters are reserved for the
# assembled figures below.
panel_files <- c(
  main_pca = file.path(panel_root, "main_pca_by_dataset_89mm.svg"),
  technical_pca = file.path(panel_root, "technical_block_pca_89mm.svg"),
  anatomy_pca = file.path(panel_root, "anatomical_pca_183mm.svg"),
  scree = file.path(panel_root, "scree_numeric_pc_order_89mm.svg"),
  concordance = file.path(panel_root, "primary_complete_case_concordance_89mm.svg"),
  experimental_group = file.path(panel_root, "experimental_group_pca_extended_data_89mm.svg"),
  normalization = file.path(panel_root, "normalization_diagnostics_183mm.svg"),
  imputation = file.path(panel_root, "ranked_imputation_183mm.svg"),
  correlation = file.path(panel_root, "clustered_sample_correlation_89mm.svg"),
  missingness = file.path(panel_root, "variable_detection_missingness_89mm.svg"),
  umap = file.path(panel_root, "umap_exploratory_89mm.svg"),
  tsne = file.path(panel_root, "tsne_exploratory_89mm.svg")
)
joint_pub_save_svg(p_pca_dataset, panel_files[["main_pca"]], 89, 83)
joint_pub_save_svg(p_pca_technical, panel_files[["technical_pca"]], 89, 91)
joint_pub_save_svg(p_pca_anatomy, panel_files[["anatomy_pca"]], 183, 76)
joint_pub_save_svg(p_scree, panel_files[["scree"]], 89, 62)
joint_pub_save_svg(p_concordance, panel_files[["concordance"]], 89, 70)
joint_pub_save_svg(p_pca_group, panel_files[["experimental_group"]], 89, 91)
joint_pub_save_svg(p_normalization, panel_files[["normalization"]], 183, 58)
joint_pub_save_svg(p_imputation, panel_files[["imputation"]], 183, 62)
joint_pub_save_svg(p_correlation, panel_files[["correlation"]], 89, 89)
joint_pub_save_svg(p_missingness, panel_files[["missingness"]], 89, 92)
joint_pub_save_svg(p_umap, panel_files[["umap"]], 89, 83)
joint_pub_save_svg(p_tsne, panel_files[["tsne"]], 89, 83)

main_row_1 <- patchwork::wrap_plots(
  p_pca_dataset + ggplot2::theme(legend.position = "none"),
  p_pca_technical,
  nrow = 1,
  widths = c(0.92, 1.08)
)
main_row_3 <- patchwork::wrap_plots(p_scree, p_concordance, nrow = 1, widths = c(0.9, 1.1))
main_figure <- patchwork::wrap_plots(
  main_row_1,
  p_pca_anatomy,
  main_row_3,
  ncol = 1,
  heights = c(1.04, 0.92, 0.82)
) +
  patchwork::plot_annotation(tag_levels = "a", theme = joint_pub_panel_tag_theme())

extended_row_1 <- patchwork::wrap_plots(p_normalization, p_pca_group, nrow = 1, widths = c(1.12, 0.88))
extended_row_3 <- patchwork::wrap_plots(p_correlation, p_missingness, nrow = 1)
extended_row_4 <- patchwork::wrap_plots(p_umap, p_tsne + ggplot2::theme(legend.position = "none"), nrow = 1)
extended_figure <- patchwork::wrap_plots(
  extended_row_1,
  p_imputation,
  extended_row_3,
  extended_row_4,
  ncol = 1,
  heights = c(0.78, 0.62, 1.18, 0.92)
) +
  patchwork::plot_annotation(tag_levels = "a", theme = joint_pub_panel_tag_theme())

main_file <- file.path(figure_root, "joint_compartment_qc_main_figure_183mm.svg")
extended_file <- file.path(figure_root, "joint_compartment_qc_extended_data_figure_183mm.svg")
joint_pub_save_svg(main_figure, main_file, 183, 188)
joint_pub_save_svg(extended_figure, extended_file, 183, 235)

report_file <- file.path(report_root, "joint_compartment_qc_publication_figure_manifest.md")
writeLines(c(
  "# Joint-compartment QC publication figures",
  "",
  "This is a plot-only rendering of the completed joint-QC outputs. PCA, UMAP and t-SNE coordinates, sample inclusion, normalization, imputation and statistical results were not recomputed or changed.",
  "",
  paste0("- Main figure: `", normalizePath(main_file, winslash = "/", mustWork = FALSE), "` (183 x 188 mm)"),
  paste0("- Extended Data figure: `", normalizePath(extended_file, winslash = "/", mustWork = FALSE), "` (183 x 235 mm)"),
  paste0("- Editable source panels: `", normalizePath(panel_root, winslash = "/", mustWork = FALSE), "`"),
  paste0("- PCA axes: ", pc_x_label, "; ", pc_y_label, "."),
  paste0("- Primary/complete-case distance concordance: Pearson r = ", sprintf("%.6f", distance_r), "."),
  paste0("- Correlation heatmap scale: ", sprintf("%.2f", corr_floor), " to 1.00; average-linkage clustering on 1 - Pearson correlation."),
  paste0("- Missingness heatmap: ", nrow(selection), " broad-union ProteinGroupIDs with the greatest max-minus-min detection-rate range across biological datasets; samples ordered by dataset, plate, region and layer; proteins average-linkage clustered by binary detection pattern."),
  "- Experimental groups are discrete CON/RES/SUS values mapped from 1/2/3.",
  "- Dataset colours: Microglia-enriched ROI #a8d8d5; Neuropil #176a5a; Soma #767577.",
  "- Standalone panels have no titles, subtitles or panel letters. Lowercase bold panel letters occur only in assembled figures.",
  "- No existing source-data table was overwritten."
), report_file)

write_run_manifest(
  file.path(log_root, "run_manifest.yml"),
  inputs = as.list(inputs),
  outputs = list(
    main_figure = main_file,
    extended_data_figure = extended_file,
    panel_directory = panel_root,
    report = report_file
  ),
  parameters = list(
    rendering_only = TRUE,
    base_font = "sans-serif 7 pt",
    axis_line_width_pt = "0.4",
    main_dimensions_mm = "183x188",
    extended_dimensions_mm = "183x235",
    missingness_feature_rule = "top 200 max-minus-min detection-rate range across datasets",
    correlation_scale = paste0(sprintf("%.2f", corr_floor), " to 1.00")
  ),
  notes = "Publication rendering only; completed PCA/UMAP/t-SNE coordinates and source tables are consumed unchanged."
)

cat("Joint-compartment QC publication figures written to:\n", normalizePath(figure_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
