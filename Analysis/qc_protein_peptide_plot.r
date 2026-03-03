#' @title QC Protein Peptide Plot
#' @description
#' Creates a quality control visualization of protein and peptide counts across different cell types.
#' The function generates a faceted bar plot with dual y-axes, where protein counts are shown on the left axis and peptide counts on the right axis. Each cell type is displayed in a separate panel, with horizontal dashed lines indicating overall averages for proteins and peptides. The plot excludes background samples and uses custom color schemes to distinguish between cell types and measurement types.
#' @param data A data frame containing protein and peptide counts.
#' @return A ggplot object.
#' @export

# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  ggplot2, tidyr, dplyr, readxl, UpSetR, ComplexHeatmap, FactoMineR, factoextra, reshape2, ggrepel, patchwork, viridis
)
library(RColorBrewer)

# Read data
data <- read_excel("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/pg_matrix/raw/quicksearch.stats.annotated.xlsx")
saving_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/QC/"

# Prepare long format for protein/peptide plot
data_long <- data %>%
  pivot_longer(cols = c(Proteins.Identified, Precursors.Identified), names_to = "Type", values_to = "Count") %>%
  mutate(ColorGroup = paste(celltype_layer, Type, sep = "_")) %>%
  arrange(celltype_layer, sample_id) %>%
  mutate(id = as.character(sample_id))

# Custom colors for protein/peptide plot
custom_colors <- c(
  "microglia_Proteins.Identified" = "#418caf",
  "microglia_Precursors.Identified" = "#4aa2c5",
  "neuron_soma_Proteins.Identified" = "#e06565",
  "neuron_soma_Precursors.Identified" = "#f0aca0",
  "neuron_neuropil_Proteins.Identified" = "#aa9f8f",
  "neuron_neuropil_Precursors.Identified" = "#ddceb9"
)

# Axis scaling
protein_range <- range(data_long$Count[data_long$Type == "Proteins.Identified"], na.rm = TRUE)
peptide_range <- range(data_long$Count[data_long$Type == "Precursors.Identified"], na.rm = TRUE)
protein_min <- protein_range[1]
protein_max <- protein_range[2]
peptide_min <- peptide_range[1]
peptide_max <- peptide_range[2]
protein_span <- protein_max - protein_min
peptide_span <- peptide_max - peptide_min

data_long <- data_long %>%
  mutate(Count_scaled = case_when(
    Type == "Proteins.Identified" ~ Count,
    Type == "Precursors.Identified" ~ (Count - peptide_min) / peptide_span * protein_span + protein_min
  ))

# Overall averages
overall_avg <- data %>%
  pivot_longer(cols = c(Proteins.Identified, Precursors.Identified), names_to = "Type", values_to = "Count") %>%
  group_by(Type) %>%
  summarise(AvgCount = mean(Count, na.rm = TRUE)) %>%
  mutate(AvgCount_scaled = case_when(
    Type == "Proteins.Identified" ~ AvgCount,
    Type == "Precursors.Identified" ~ (AvgCount - peptide_min) / peptide_span * protein_span + protein_min
  ))

# Protein/peptide plot
plt <- ggplot(data_long, aes(x = id, y = Count_scaled, fill = ColorGroup)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(data = overall_avg, aes(yintercept = AvgCount_scaled, color = Type), linetype = "dashed", size = 1, inherit.aes = FALSE) +
  facet_grid(. ~ celltype_layer, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = custom_colors, labels = names(custom_colors), name = NULL) +
  scale_color_manual(values = c("Proteins.Identified" = "#457b9d", "Precursors.Identified" = "#e63946"), guide = "none") +
  scale_y_continuous(
    name = "Protein Count",
    expand = c(0, 0),
    sec.axis = sec_axis(
      trans = ~ (.-protein_min) / protein_span * peptide_span + peptide_min,
      name = "Peptide Count"
    )
  ) +
  scale_x_discrete(expand = c(0, 0), breaks = unique(data_long$id), labels = unique(data_long$id)) +
  labs(x = "Sample ID") +
  theme_minimal(base_size = 30) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 35),
    axis.title.y.right = element_text(size = 35),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 25),
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(face = "bold", size = 35),
    panel.spacing = unit(0, 'lines')
  )

ggsave(filename = paste0(saving_dir, "protein_peptide_counts.svg"), plot = plt, width = 60, height = 30, units = "cm")
print(plt)

# Use ColorBrewer Set2 palette for celltype colors
celltype_levels <- unique(data$celltype_layer)
palette_size <- length(celltype_levels)
celltype_colors <- brewer.pal(min(palette_size, 8), "Pastel2")
if (palette_size > 8) celltype_colors <- rep(celltype_colors, length.out = palette_size)
names(celltype_colors) <- celltype_levels

# 9. Missing Value Distribution
missing_per_sample <- rowSums(is.na(data))
plt_missing_sample <- ggplot(data.frame(sample_id = data$sample_id, missing = missing_per_sample), aes(x = sample_id, y = missing)) +
  geom_col(fill = "#e63946") +
  labs(title = "Missing Values per Sample", x = "Sample ID", y = "Missing Value Count") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0(saving_dir, "qc_missing_sample.svg"), plot = plt_missing_sample, width = 15, height = 8, units = "cm")

missing_per_celltype <- data %>% group_by(celltype_layer) %>% summarise(missing = sum(is.na(.)))
plt_missing_celltype <- ggplot(missing_per_celltype, aes(x = celltype_layer, y = missing, fill = celltype_layer)) +
  geom_col() +
  labs(title = "Missing Values per Celltype", x = "Celltype", y = "Missing Value Count") +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal(base_size = 18)
ggsave(filename = paste0(saving_dir, "qc_missing_celltype.svg"), plot = plt_missing_celltype, width = 10, height = 8, units = "cm")

# 10. Protein/Peptide Overlap (UpSet Plot)
protein_sets <- split(data$Proteins.Identified, data$celltype_layer)
protein_sets <- lapply(protein_sets, unique)
protein_list <- list()
for (ct in names(protein_sets)) protein_list[[ct]] <- protein_sets[[ct]]
upset_data <- fromList(protein_list)
png(filename = paste0(saving_dir, "qc_upset_protein.png"), width = 1200, height = 800)
upset(upset_data, nsets = length(protein_list), order.by = "freq", main.bar.color = "#457b9d", sets.bar.color = "#e63946")
dev.off()

# 11. Intensity Distribution
plt_ms1_box <- ggplot(data, aes(x = celltype_layer, y = MS1.Signal, fill = celltype_layer)) +
  geom_boxplot() +
  labs(title = "MS1 Intensity Distribution", x = "Celltype", y = "MS1 Signal") +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal(base_size = 18)
ggsave(filename = paste0(saving_dir, "qc_ms1_intensity_box.svg"), plot = plt_ms1_box, width = 10, height = 8, units = "cm")

plt_ms2_box <- ggplot(data, aes(x = celltype_layer, y = MS2.Signal, fill = celltype_layer)) +
  geom_boxplot() +
  labs(title = "MS2 Intensity Distribution", x = "Celltype", y = "MS2 Signal") +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal(base_size = 18)
ggsave(filename = paste0(saving_dir, "qc_ms2_intensity_box.svg"), plot = plt_ms2_box, width = 10, height = 8, units = "cm")

# 12. PCA
pca_vars <- c("Proteins.Identified", "Precursors.Identified", "MS1.Signal", "MS2.Signal", "Median.Mass.Acc.MS1", "Median.Mass.Acc.MS2", "Normalisation.Instability", "Median.RT.Prediction.Acc")
pca_data <- data[, pca_vars]
pca_data <- na.omit(pca_data)
pca_res <- PCA(pca_data, graph = FALSE)
plt_pca <- fviz_pca_ind(
  pca_res,
  geom = "point",
  pointshape = 21,
  pointsize = 3,
  col.ind = data$celltype_layer[as.numeric(rownames(pca_data))],
  palette = celltype_colors,
  addEllipses = TRUE,
  legend.title = "Celltype"
) + ggtitle("PCA of QC Metrics")
ggsave(filename = paste0(saving_dir, "qc_pca.svg"), plot = plt_pca, width = 12, height = 10, units = "cm")

# 13. Correlation Heatmaps
cor_data <- data[, pca_vars]
cor_data <- cor_data[, sapply(cor_data, is.numeric)]
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")
plt_cor_heat <- Heatmap(cor_matrix, name = "Correlation", col = colorRampPalette(c("#457b9d", "white", "#e63946"))(100))
png(filename = paste0(saving_dir, "qc_correlation_heatmap.png"), width = 1200, height = 1000)
draw(plt_cor_heat)
dev.off()

# 14. Batch Effects
if ("batch" %in% colnames(data)) {
  plt_batch <- ggplot(data, aes(x = batch, y = Normalisation.Instability, fill = batch)) +
    geom_boxplot() +
    labs(title = "Normalisation Instability by Batch", x = "Batch", y = "Normalisation Instability") +
    theme_minimal(base_size = 18)
  ggsave(filename = paste0(saving_dir, "qc_batch_effects.svg"), plot = plt_batch, width = 10, height = 8, units = "cm")
}

# 15. Peptide Modification Frequencies
if ("PTM" %in% colnames(data)) {
  ptm_freq <- data %>% group_by(celltype_layer, PTM) %>% summarise(count = n())
  plt_ptm <- ggplot(ptm_freq, aes(x = PTM, y = count, fill = celltype_layer)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "PTM Frequencies by Celltype", x = "PTM", y = "Count") +
    scale_fill_manual(values = celltype_colors) +
    theme_minimal(base_size = 18)
  ggsave(filename = paste0(saving_dir, "qc_ptm_freq.svg"), plot = plt_ptm, width = 15, height = 8, units = "cm")
}

# 16. Identification Rate
if ("Total.Peptides" %in% colnames(data)) {
  data$id_rate_pep <- data$Precursors.Identified / data$Total.Peptides
  plt_idrate_pep <- ggplot(data, aes(x = sample_id, y = id_rate_pep, fill = celltype_layer)) +
    geom_col() +
    labs(title = "Peptide Identification Rate", x = "Sample ID", y = "ID Rate") +
    scale_fill_manual(values = celltype_colors) +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename = paste0(saving_dir, "qc_idrate_pep.svg"), plot = plt_idrate_pep, width = 15, height = 8, units = "cm")
}
if ("Total.Proteins" %in% colnames(data)) {
  data$id_rate_prot <- data$Proteins.Identified / data$Total.Proteins
  plt_idrate_prot <- ggplot(data, aes(x = sample_id, y = id_rate_prot, fill = celltype_layer)) +
    geom_col() +
    labs(title = "Protein Identification Rate", x = "Sample ID", y = "ID Rate") +
    scale_fill_manual(values = celltype_colors) +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename = paste0(saving_dir, "qc_idrate_prot.svg"), plot = plt_idrate_prot, width = 15, height = 8, units = "cm")
}

# 17. Sample Outlier Detection (Mahalanobis Distance)
mahal_data <- data[, pca_vars]
mahal_data <- na.omit(mahal_data)
mahal_dist <- mahalanobis(mahal_data, colMeans(mahal_data), cov(mahal_data))
mahal_df <- data.frame(
  sample_id = data$sample_id[as.numeric(rownames(mahal_data))],
  mahal_dist = mahal_dist,
  celltype_layer = data$celltype_layer[as.numeric(rownames(mahal_data))]
)
plt_mahal <- ggplot(mahal_df, aes(x = sample_id, y = mahal_dist, fill = celltype_layer)) +
  geom_col() +
  labs(title = "Sample Outlier Detection (Mahalanobis Distance)", x = "Sample ID", y = "Mahalanobis Distance") +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0(saving_dir, "qc_mahalanobis.svg"), plot = plt_mahal, width = 15, height = 8, units = "cm")

# Improved QC Plots and Combined Output

# 1. MS1 and MS2 Signal Distribution
plt_ms_signal <- ggplot(data, aes(x = MS1.Signal, y = MS2.Signal, color = celltype_layer)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "MS1 vs\nMS2 Signal", x = "MS1 Signal", y = "MS2 Signal") +
  scale_color_manual(values = celltype_colors, name = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 2. FWHM Distribution (Scans and RT)
plt_fwhm <- ggplot(data, aes(x = FWHM.Scans, y = FWHM.RT, color = celltype_layer)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "FWHM:\nScans vs RT", x = "FWHM (Scans)", y = "FWHM (RT)") +
  scale_color_manual(values = celltype_colors, name = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# mass accuracy MA plot
# MA plot: log2 ratio (y) vs log2 intensity (x)
if (all(c("MS1.Signal", "MS2.Signal") %in% colnames(data))) {
  ma_data <- data %>%
    mutate(
      log2_intensity = log2((MS1.Signal + MS2.Signal) / 2),
      log2_ratio = log2(MS2.Signal / MS1.Signal)
    )
  plt_ma <- ggplot(ma_data, aes(x = log2_intensity, y = log2_ratio, color = celltype_layer)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "MA Plot: log2 Ratio vs log2 Intensity", x = "log2(Intensity)", y = "log2(MS2/MS1)") +
    scale_color_manual(values = celltype_colors, name = NULL) +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "bottom",
      aspect.ratio = 1,
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(filename = paste0(saving_dir, "qc_ma_plot.svg"), plot = plt_ma, width = 10, height = 10, units = "cm")
}

# 3. Mass Accuracy (MS1 and MS2, raw and corrected)
plt_massacc <- ggplot(data, aes(x = Median.Mass.Acc.MS1, y = Median.Mass.Acc.MS2, color = celltype_layer)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Median Mass Accuracy:\nMS1 vs MS2", x = "MS1 (ppm)", y = "MS2 (ppm)") +
  scale_color_manual(values = celltype_colors, name = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

plt_massacc_corr <- ggplot(data, aes(x = Median.Mass.Acc.MS1.Corrected, y = Median.Mass.Acc.MS2.Corrected, color = celltype_layer)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Corrected Median Mass Accuracy:\nMS1 vs MS2", x = "MS1 Corrected (ppm)", y = "MS2 Corrected (ppm)") +
  scale_color_manual(values = celltype_colors, name = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 4. Normalisation Instability
plt_norminstab <- ggplot(data, aes(x = celltype_layer, y = Normalisation.Instability, fill = celltype_layer)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = celltype_layer), shape = 16, size = 1.5) +
  labs(title = "Normalisation Instability\nby Celltype", x = "Celltype", y = "Normalisation Instability") +
  scale_fill_manual(values = celltype_colors, name = NULL) +
  scale_color_manual(values = celltype_colors, guide = "none") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0, 'pt'),
    legend.margin = margin(0,0,0,0),
    legend.key.size = unit(0.7, "lines"),
    legend.text = element_text(size = 12),
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 5. RT Prediction Accuracy
plt_rtacc <- ggplot(data, aes(x = celltype_layer, y = Median.RT.Prediction.Acc, fill = celltype_layer)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = celltype_layer), shape = 16, size = 1.5) +
  labs(title = "Median RT Prediction Accuracy\nby Celltype", x = "Celltype", y = "Median RT Prediction Accuracy") +
  scale_fill_manual(values = celltype_colors, name = NULL) +
  scale_color_manual(values = celltype_colors, guide = "none") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0, 'pt'),
    legend.margin = margin(0,0,0,0),
    legend.key.size = unit(0.7, "lines"),
    legend.text = element_text(size = 12),
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 6. Peptide Length, Charge, Missed Cleavages
plt_peplen <- ggplot(data, aes(x = celltype_layer, y = Average.Peptide.Length, fill = celltype_layer)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = celltype_layer), shape = 16, size = 1.5) +
  labs(title = "Average Peptide Length\nby Celltype", x = "Celltype", y = "Average Peptide Length") +
  scale_fill_manual(values = celltype_colors, name = NULL) +
  scale_color_manual(values = celltype_colors, guide = "none") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0, 'pt'),
    legend.margin = margin(0,0,0,0),
    legend.key.size = unit(0.7, "lines"),
    legend.text = element_text(size = 12),
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

plt_pepcharge <- ggplot(data, aes(x = celltype_layer, y = Average.Peptide.Charge, fill = celltype_layer)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = celltype_layer), shape = 16, size = 1.5) +
  labs(title = "Average Peptide Charge\nby Celltype", x = "Celltype", y = "Average Peptide Charge") +
  scale_fill_manual(values = celltype_colors, name = NULL) +
  scale_color_manual(values = celltype_colors, guide = "none") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0, 'pt'),
    legend.margin = margin(0,0,0,0),
    legend.key.size = unit(0.7, "lines"),
    legend.text = element_text(size = 12),
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

plt_pepmiss <- ggplot(data, aes(x = celltype_layer, y = Average.Missed.Tryptic.Cleavages, fill = celltype_layer)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = celltype_layer), shape = 16, size = 1.5) +
  labs(title = "Average Missed Tryptic Cleavages\nby Celltype", x = "Celltype", y = "Average Missed Tryptic Cleavages") +
  scale_fill_manual(values = celltype_colors, name = NULL) +
  scale_color_manual(values = celltype_colors, guide = "none") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0, 'pt'),
    legend.margin = margin(0,0,0,0),
    legend.key.size = unit(0.7, "lines"),
    legend.text = element_text(size = 12),
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save each plot individually
ggsave(filename = paste0(saving_dir, "qc_ms_signal.svg"), plot = plt_ms_signal, width = 10, height = 10, units = "cm")
ggsave(filename = paste0(saving_dir, "qc_fwhm.svg"), plot = plt_fwhm, width = 10, height = 10, units = "cm")
ggsave(filename = paste0(saving_dir, "qc_massacc.svg"), plot = plt_massacc, width = 10, height = 10, units = "cm")
ggsave(filename = paste0(saving_dir, "qc_massacc_corr.svg"), plot = plt_massacc_corr, width = 10, height = 10, units = "cm")
ggsave(filename = paste0(saving_dir, "qc_norminstab.svg"), plot = plt_norminstab, width = 10, height = 10, units = "cm")
ggsave(filename = paste0(saving_dir, "qc_rtacc.svg"), plot = plt_rtacc, width = 10, height = 10, units = "cm")
ggsave(filename = paste0(saving_dir, "qc_peplen.svg"), plot = plt_peplen, width = 10, height = 10, units = "cm")
ggsave(filename = paste0(saving_dir, "qc_pepcharge.svg"), plot = plt_pepcharge, width = 10, height = 10, units = "cm")
ggsave(filename = paste0(saving_dir, "qc_pepmiss.svg"), plot = plt_pepmiss, width = 10, height = 10, units = "cm")

# Combine all plots into one file using patchwork
all_plots <- (
  plt_ms_signal + plt_fwhm + plt_massacc + plt_massacc_corr
) /
  (
    plt_norminstab + plt_rtacc + plt_peplen + plt_pepcharge + plt_pepmiss
  )
ggsave(filename = paste0(saving_dir, "qc_all_plots.svg"), plot = all_plots, width = 40, height = 40, units = "cm")

# 7. Summary Table: Mean/SD for key metrics by celltype_layer
qc_summary <- data %>%
  group_by(celltype_layer) %>%
  summarise(
    n = n(),
    Proteins.Identified_mean = mean(Proteins.Identified, na.rm = TRUE),
    Proteins.Identified_sd = sd(Proteins.Identified, na.rm = TRUE),
    Precursors.Identified_mean = mean(Precursors.Identified, na.rm = TRUE),
    Precursors.Identified_sd = sd(Precursors.Identified, na.rm = TRUE),
    MS1.Signal_mean = mean(MS1.Signal, na.rm = TRUE),
    MS1.Signal_sd = sd(MS1.Signal, na.rm = TRUE),
    MS2.Signal_mean = mean(MS2.Signal, na.rm = TRUE),
    MS2.Signal_sd = sd(MS2.Signal, na.rm = TRUE),
    FWHM.Scans_mean = mean(FWHM.Scans, na.rm = TRUE),
    FWHM.Scans_sd = sd(FWHM.Scans, na.rm = TRUE),
    FWHM.RT_mean = mean(FWHM.RT, na.rm = TRUE),
    FWHM.RT_sd = sd(FWHM.RT, na.rm = TRUE),
    Median.Mass.Acc.MS1_mean = mean(Median.Mass.Acc.MS1, na.rm = TRUE),
    Median.Mass.Acc.MS1_sd = sd(Median.Mass.Acc.MS1, na.rm = TRUE),
    Median.Mass.Acc.MS2_mean = mean(Median.Mass.Acc.MS2, na.rm = TRUE),
    Median.Mass.Acc.MS2_sd = sd(Median.Mass.Acc.MS2, na.rm = TRUE),
    Normalisation.Instability_mean = mean(Normalisation.Instability, na.rm = TRUE),
    Normalisation.Instability_sd = sd(Normalisation.Instability, na.rm = TRUE),
    Median.RT.Prediction.Acc_mean = mean(Median.RT.Prediction.Acc, na.rm = TRUE),
    Median.RT.Prediction.Acc_sd = sd(Median.RT.Prediction.Acc, na.rm = TRUE),
    Average.Peptide.Length_mean = mean(Average.Peptide.Length, na.rm = TRUE),
    Average.Peptide.Length_sd = sd(Average.Peptide.Length, na.rm = TRUE),
    Average.Peptide.Charge_mean = mean(Average.Peptide.Charge, na.rm = TRUE),
    Average.Peptide.Charge_sd = sd(Average.Peptide.Charge, na.rm = TRUE),
    Average.Missed.Tryptic.Cleavages_mean = mean(Average.Missed.Tryptic.Cleavages, na.rm = TRUE),
    Average.Missed.Tryptic.Cleavages_sd = sd(Average.Missed.Tryptic.Cleavages, na.rm = TRUE)
  )
print(qc_summary)
write.csv(qc_summary, file = paste0(saving_dir, "qc_summary.csv"), row.names = FALSE)

# 8. Table: Outlier detection for Normalisation Instability
norminstab_outliers <- data %>%
  group_by(celltype_layer) %>%
  mutate(norminstab_z = (Normalisation.Instability - mean(Normalisation.Instability, na.rm = TRUE)) / sd(Normalisation.Instability, na.rm = TRUE)) %>%
  filter(abs(norminstab_z) > 3)
print(norminstab_outliers)
write.csv(norminstab_outliers, file = paste0(saving_dir, "normalisation_instability_outliers.csv"), row.names = FALSE)


