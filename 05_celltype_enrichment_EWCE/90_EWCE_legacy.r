# ==========================================
# 1. SETUP, LIBRARIES & DIRECTORIES
# ==========================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("EWCE", quietly = TRUE)) BiocManager::install("EWCE")

# Added Rtsne for dimensionality reduction
packages <- c("readxl", "dplyr", "tibble", "tidyr", "ggplot2", "pheatmap", 
              "svglite", "ggridges", "ggrepel", "ggsci", "viridis", 
              "UpSetR", "openxlsx", "stringr", "Rtsne")
invisible(lapply(packages, function(x) if (!require(x, character.only = TRUE)) install.packages(x)))

library(EWCE); library(ewceData); library(readxl); library(dplyr)
library(tibble); library(tidyr); library(ggplot2); library(pheatmap)
library(openxlsx); library(stringr); library(ggridges); library(ggrepel)
library(ggsci); library(Rtsne); library(UpSetR); library(viridis)

# Define and create output structure
base_results <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/EWCE"
plot_path    <- file.path(base_results, "Plots")
table_path   <- file.path(base_results, "Tables")
qc_path      <- file.path(base_results, "QC_Logs")

lapply(c(plot_path, table_path, qc_path), function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

# ==========================================
# 2. DATA CLEANING & SPECIES FILTERING
# ==========================================
message("Step 1: Loading Data and Filtering for Mouse Proteins...")
ctd <- ewceData::ctd()
ref_genes <- rownames(ctd[[1]]$specificity)

path_to_file <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/imputed_data.xlsx"
raw_df <- read_excel(path_to_file)

# SPECIES FILTER: Keep only _MOUSE and sanitize symbols (e.g., KRT10 -> Krt10)
clean_df <- raw_df %>%
  filter(!is.na(Genes), Genes != "0") %>%
  filter(grepl("_MOUSE", Gene_Group)) %>% 
  mutate(Genes = str_to_title(Genes)) %>%  
  distinct(Genes, .keep_all = TRUE)

# Track non-mouse removals for QC
contaminants <- raw_df %>% filter(!grepl("_MOUSE", Gene_Group))
mapped_genes <- intersect(clean_df$Genes, ref_genes)
dropped_genes <- setdiff(clean_df$Genes, ref_genes)

# Create Grouped Matrix
exp_matrix_grouped <- clean_df %>%
  filter(Genes %in% mapped_genes) %>%
  select(Genes, matches("neuron|mcherry|cfos|bg")) %>%
  pivot_longer(cols = -Genes, names_to = "Raw_Sample", values_to = "Expression") %>%
  mutate(Grouped_Name = gsub("\\.\\.\\.[0-9]+$", "", Raw_Sample)) %>%
  group_by(Genes, Grouped_Name) %>%
  summarise(Mean_Exp = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Grouped_Name, values_from = Mean_Exp) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

# ==========================================
# 3. EWCE ANALYSIS LOOP
# ==========================================
all_results_lvl1 <- data.frame()
all_results_lvl2 <- data.frame()
gene_specificity_map <- data.frame()
group_names <- colnames(exp_matrix_grouped)

for (gname in group_names) {
  message("Processing Group: ", gname)
  hits <- rownames(exp_matrix_grouped)[order(exp_matrix_grouped[, gname], decreasing = TRUE)[1:200]]
  curr_cond <- gsub("_[0-9]+.*$", "", gname)
  
  # Enrichment Level 1 (Broad)
  res_l1 <- bootstrap_enrichment_test(sct_data = ctd, hits = hits, reps = 1000, 
                                      annotLevel = 1, genelistSpecies = "mouse", sctSpecies = "mouse")
  all_results_lvl1 <- rbind(all_results_lvl1, res_l1$results %>% mutate(Group = gname, Condition = curr_cond))
  
  # Enrichment Level 2 (Granular)
  res_l2 <- bootstrap_enrichment_test(sct_data = ctd, hits = hits, reps = 1000, 
                                      annotLevel = 2, genelistSpecies = "mouse", sctSpecies = "mouse")
  all_results_lvl2 <- rbind(all_results_lvl2, res_l2$results %>% mutate(Group = gname, Condition = curr_cond))
  
  # Extract specificity for gene-level plots
  common_hits <- intersect(hits, ref_genes)
  gene_spec <- as.data.frame(ctd[[1]]$specificity[common_hits, ]) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "CellType", values_to = "SpecificityScore") %>%
    mutate(Group = gname, Condition = curr_cond)
  gene_specificity_map <- rbind(gene_specificity_map, gene_spec)
}

all_results_lvl2 <- all_results_lvl2 %>% group_by(Group) %>% mutate(q = p.adjust(p, method = "BH")) %>% ungroup()

# ==========================================
# 4. DIMENSIONALITY REDUCTION & CONSENSUS
# ==========================================
# Prep Z-matrix for PCA/t-SNE
z_matrix <- all_results_lvl1 %>% 
  select(Group, CellType, sd_from_mean) %>% 
  pivot_wider(names_from = Group, values_from = sd_from_mean) %>% 
  column_to_rownames("CellType")

# PCA
pca_res <- prcomp(t(z_matrix), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x) %>% rownames_to_column("Group") %>% mutate(Condition = gsub("_[0-9]+.*$", "", Group))

# t-SNE (perplexity set low for small sample size)
set.seed(42)
tsne_res <- Rtsne(t(z_matrix), perplexity = 1, check_duplicates = FALSE) 
tsne_df <- as.data.frame(tsne_res$Y) %>% rename(tSNE1=V1, tSNE2=V2) %>% 
  mutate(Group = colnames(z_matrix), Condition = gsub("_[0-9]+.*$", "", Group))

# Consensus Calculation
consensus <- all_results_lvl1 %>% 
  group_by(Condition, CellType) %>% 
  summarise(Mean_Z = mean(sd_from_mean), SE = sd(sd_from_mean)/sqrt(n()), .groups = "drop")

# ==========================================
# 5. PUBLICATION VISUALIZATIONS
# ==========================================
theme_publication <- function() {
  theme_classic(base_size = 8) +
    theme(text = element_text(family = "sans", color = "black"),
          axis.text = element_text(size = 7, color = "black"), axis.title = element_text(size = 8, face = "bold"),
          strip.text = element_text(size = 8, face = "bold"), strip.background = element_blank(),
          legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6))
}

# FIG 1: Broad Identity Dot Plot
p1 <- ggplot(all_results_lvl1 %>% filter(grepl("interneurons|pyramidal", CellType)), 
              aes(x = Group, y = CellType, size = -log10(p), color = sd_from_mean)) +
  geom_point(alpha = 0.8) + scale_color_distiller(palette = "RdBu", direction = -1) +
  theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# FIG 3: Volcano Plot
p3 <- ggplot(all_results_lvl2, aes(x = sd_from_mean, y = -log10(p))) +
  geom_point(aes(color = Condition), alpha = 0.6, size = 1.5) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  geom_text_repel(data = subset(all_results_lvl2, q < 0.05 & sd_from_mean > 3), 
                  aes(label = CellType), size = 2, max.overlaps = 20) +
  scale_color_npg() + theme_publication()

# FIG 4: Lollipop Plot
sig_lvl2 <- all_results_lvl2 %>% filter(q < 0.05) %>% group_by(Group) %>% slice_max(sd_from_mean, n = 5) %>% ungroup()
p4 <- ggplot(sig_lvl2, aes(x = reorder(CellType, sd_from_mean), y = sd_from_mean, color = Condition)) +
  geom_segment(aes(xend = CellType, yend = 0), linewidth = 0.5) + geom_point(size = 2) + coord_flip() + 
  facet_wrap(~Group, scales = "free_y", ncol = 4) + scale_color_npg() + theme_publication()

# FIG 5: Ridge Plot (Specificity)
p5 <- ggplot(gene_specificity_map %>% filter(SpecificityScore > 0.1), 
              aes(x = SpecificityScore, y = Condition, fill = Condition)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2) + scale_fill_npg() + theme_publication()

# FIG 6: Consensus Pointrange
p6 <- ggplot(consensus %>% filter(Mean_Z > 1), aes(x = reorder(CellType, Mean_Z), y = Mean_Z, color = Condition)) +
  geom_pointrange(aes(ymin = Mean_Z - SE, ymax = Mean_Z + SE), position = position_dodge(width = 0.5), size = 0.3) +
  coord_flip() + scale_color_npg() + theme_publication()

# FIG 7: PCA
p7 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) + geom_text_repel(aes(label = Group), size = 2.5) +
  scale_color_npg() + theme_publication()

# FIG 8: Driver Bubble Plot
drivers <- gene_specificity_map %>% group_by(Group, CellType) %>% slice_max(SpecificityScore, n=15) %>% ungroup()
top_genes <- drivers %>% group_by(Condition) %>% slice_max(SpecificityScore, n=3, with_ties=FALSE) %>% pull(Gene) %>% unique()
p8 <- ggplot(gene_specificity_map %>% filter(Gene %in% top_genes, SpecificityScore > 0.05), 
              aes(x = Condition, y = Gene, size = SpecificityScore, color = SpecificityScore)) +
  geom_point(alpha = 0.8) + scale_color_viridis(option = "magma") + theme_publication()

# FIG 10: t-SNE Plot
p10 <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Condition)) +
  geom_point(size = 3) + geom_text_repel(aes(label = Group), size = 2.5) +
  scale_color_npg() + theme_publication() + labs(title="t-SNE of Sample Cell-Type Profiles")

# ==========================================
# 6. EXPORT TABLES & VISUALS
# ==========================================
message("Step 4: Final Exporting...")

# Tables
write.csv(all_results_lvl1, file.path(table_path, "Table_S1_Lvl1_Full.csv"), row.names = FALSE)
write.csv(all_results_lvl2, file.path(table_path, "Table_S2_Lvl2_Full.csv"), row.names = FALSE)
write.csv(sig_lvl2,         file.path(table_path, "Table_S3_Top_Significant_Subtypes.csv"), row.names = FALSE)
write.csv(consensus,        file.path(table_path, "Table_S4_Condition_Consensus.csv"), row.names = FALSE)
write.csv(gene_specificity_map, file.path(table_path, "Table_S5_Full_Specificity_Map.csv"), row.names = FALSE)

# Plots
ggsave(file.path(plot_path, "Fig1_DotPlot.svg"), p1, width = 120, height = 90, units = "mm")
ggsave(file.path(plot_path, "Fig3_Volcano.svg"), p3, width = 150, height = 100, units = "mm")
ggsave(file.path(plot_path, "Fig4_Lollipop.svg"), p4, width = 250, height = 200, units = "mm")
ggsave(file.path(plot_path, "Fig5_RidgePlot.svg"), p5, width = 120, height = 80, units = "mm")
ggsave(file.path(plot_path, "Fig6_Consensus.svg"), p6, width = 120, height = 100, units = "mm")
ggsave(file.path(plot_path, "Fig7_PCA.svg"), p7, width = 100, height = 90, units = "mm")
ggsave(file.path(plot_path, "Fig8_GeneBubble.svg"), p8, width = 100, height = 120, units = "mm")
ggsave(file.path(plot_path, "Fig10_tSNE.svg"), p10, width = 100, height = 90, units = "mm")

# Special PDF Exports
pdf(file.path(plot_path, "Fig2_Clustered_Heatmap.pdf"), width = 7, height = 6)
pheatmap(z_matrix, clustering_distance_cols = "euclidean", color = colorRampPalette(c("#3B4992FF", "white", "#EE0000FF"))(50), border_color = NA, fontsize = 7)
dev.off()

sig_by_cond <- all_results_lvl2 %>% filter(q < 0.05) %>% split(.$Condition) %>% lapply(function(x) unique(x$CellType))
if(length(sig_by_cond) > 1) {
  pdf(file.path(plot_path, "Fig9_UpSet_Intersection.pdf"), width = 7, height = 5)
  print(upset(fromList(sig_by_cond), order.by = "freq", main.bar.color = "#3B4992FF"))
  dev.off()
}

# Final QC Workbook
wb <- createWorkbook()
addWorksheet(wb, "Metrics"); writeData(wb, 1, data.frame(Metric=c("Original Rows", "Mouse Kept", "Contaminants Purged"), Value=c(nrow(raw_df), nrow(clean_df), nrow(contaminants))))
addWorksheet(wb, "Purged_List"); writeData(wb, 2, contaminants)
addWorksheet(wb, "Unmapped_Mouse"); writeData(wb, 3, data.frame(Gene=dropped_genes))
saveWorkbook(wb, file.path(qc_path, "EWCE_Species_QC_Final.xlsx"), overwrite = TRUE)

message("SUCCESS: Analysis complete.")






























# ==========================================
# 1. SETUP, LOGGING & LIBRARIES
# ==========================================
# Install missing core components
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("readxl", "dplyr", "tibble", "tidyr", "ggplot2", "pheatmap", 
              "svglite", "ggridges", "ggrepel", "ggsci", "viridis", 
              "UpSetR", "openxlsx", "stringr", "Rtsne", "patchwork")

invisible(lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) install.packages(x)
}))

library(EWCE); library(ewceData); library(readxl); library(dplyr)
library(tibble); library(tidyr); library(ggplot2); library(pheatmap)
library(openxlsx); library(stringr); library(ggridges); library(ggrepel)
library(ggsci); library(Rtsne); library(UpSetR); library(patchwork)

# Paths
base_results <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/EWCE"
plot_path    <- file.path(base_results, "Plots")
table_path   <- file.path(base_results, "Tables")
qc_path      <- file.path(base_results, "QC_Logs")
lapply(c(plot_path, table_path, qc_path), function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

# Start Logging
log_file <- file(file.path(qc_path, paste0("Run_Log_", Sys.Date(), ".txt")), open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

cat("========================================================\n")
cat("EWCE ANALYSIS PIPELINE - PUBLICATION STYLE\n")
cat("Run Date:", as.character(Sys.Date()), "\n")
cat("========================================================\n\n")

# ==========================================
# 2. DATA PREPROCESSING (SPECIES FIX)
# ==========================================
cat("Step 1: Data Cleaning and Species Filtering...\n")
ctd <- ewceData::ctd()
ref_genes <- rownames(ctd[[1]]$specificity)

raw_df <- read_excel("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/imputed_data.xlsx")

# Filter for Mouse, convert to Sentence Case (e.g., GAD2 -> Gad2)
clean_df <- raw_df %>%
  filter(!is.na(Genes), Genes != "0") %>%
  filter(grepl("_MOUSE", Gene_Group)) %>% 
  mutate(Genes = str_to_title(Genes)) %>%  
  distinct(Genes, .keep_all = TRUE)

mapped_genes <- intersect(clean_df$Genes, ref_genes)
dropped_genes <- setdiff(clean_df$Genes, ref_genes)
contaminants  <- raw_df %>% filter(!grepl("_MOUSE", Gene_Group))

cat(sprintf("Success: %d genes mapped to CTD. %d genes dropped. %d contaminants removed.\n", 
            length(mapped_genes), length(dropped_genes), nrow(contaminants)))

# Build Matrix
exp_matrix <- clean_df %>%
  filter(Genes %in% mapped_genes) %>%
  select(Genes, matches("neuron|mcherry|cfos|bg")) %>%
  pivot_longer(cols = -Genes, names_to = "Raw_Sample", values_to = "Exp") %>%
  mutate(Group = gsub("\\.\\.\\.[0-9]+$", "", Raw_Sample)) %>%
  group_by(Genes, Group) %>%
  summarise(Mean_Exp = mean(Exp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = Mean_Exp) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

# ==========================================
# 3. ANALYSIS LOOP
# ==========================================
cat("\nStep 2: Running EWCE Bootstrap (Levels 1 & 2)...\n")
all_lvl1 <- data.frame(); all_lvl2 <- data.frame(); gene_map <- data.frame()

for (gname in colnames(exp_matrix)) {
  cat("Processing:", gname, "... ")
  # Selection: Top 250 genes (Standard for Proteomics)
  hits <- rownames(exp_matrix)[order(exp_matrix[, gname], decreasing = TRUE)[1:250]]
  cond <- gsub("_[0-9]+.*$", "", gname)
  
  res1 <- bootstrap_enrichment_test(ctd, hits=hits, reps=1000, annotLevel=1, genelistSpecies="mouse", sctSpecies="mouse")
  all_lvl1 <- rbind(all_lvl1, res1$results %>% mutate(Group = gname, Condition = cond))
  
  res2 <- bootstrap_enrichment_test(ctd, hits=hits, reps=1000, annotLevel=2, genelistSpecies="mouse", sctSpecies="mouse")
  all_lvl2 <- rbind(all_lvl2, res2$results %>% mutate(Group = gname, Condition = cond))
  
  # Spec mapping
  spec <- as.data.frame(ctd[[1]]$specificity[intersect(hits, ref_genes), ]) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "CellType", values_to = "Score") %>%
    mutate(Group = gname, Condition = cond)
  gene_map <- rbind(gene_map, spec)
  cat("Done.\n")
}

all_lvl2 <- all_lvl2 %>% group_by(Group) %>% mutate(q = p.adjust(p, method = "BH")) %>% ungroup()

# ==========================================
# 4. PUBLICATION VISUALIZATION ENGINE
# ==========================================
theme_publication <- function() {
  theme_classic(base_size = 7) + # Compact publication text
    theme(text = element_text(family = "sans"),
          axis.title = element_text(face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 8),
          legend.key.size = unit(3, "mm"))
}

# A: PCA of Signatures
z_mat <- all_lvl1 %>% select(Group, CellType, sd_from_mean) %>% 
  pivot_wider(names_from = Group, values_from = sd_from_mean) %>% column_to_rownames("CellType")
pca <- prcomp(t(z_mat), scale. = TRUE)
pca_df <- as.data.frame(pca$x) %>% rownames_to_column("Group") %>% mutate(Cond = gsub("_[0-9]+.*$", "", Group))
p_pca <- ggplot(pca_df, aes(PC1, PC2, color=Cond)) + 
  geom_point(size=2) + geom_text_repel(aes(label=Group), size=1.8) + 
  scale_color_npg() + theme_publication() + labs(title="a Signature PCA")

# B: Broad Enrichment DotPlot
p_dot <- ggplot(all_lvl1, aes(x = Group, y = CellType, size = -log10(p), color = sd_from_mean)) +
  geom_point() + scale_color_distiller(palette = "RdBu") + 
  theme_publication() + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="b Broad Cell Identity", color="Z-score", size="-log10(P)")

# C: Granular Lollipop (Top Hits)
sig_hits <- all_lvl2 %>% filter(q < 0.05) %>% group_by(Group) %>% slice_max(sd_from_mean, n=3)
p_loli <- ggplot(sig_hits, aes(x=reorder(CellType, sd_from_mean), y=sd_from_mean, color=Condition)) +
  geom_point(size=2) + geom_segment(aes(xend=CellType, yend=0)) + coord_flip() +
  facet_wrap(~Group, scales="free_y", ncol=2) + scale_color_npg() + 
  theme_publication() + labs(title="c Top Subtype Hits", x=NULL, y="Enrichment Z-score")

# D: Specificity Ridge
p_ridge <- ggplot(gene_map %>% filter(Score > 0.1), aes(x=Score, y=Condition, fill=Condition)) +
  geom_density_ridges(alpha=0.6, scale=1.2) + scale_fill_npg() + 
  theme_publication() + labs(title="d Marker Specificity Distribution")

# ==========================================
# 5. COMPOSITE ASSEMBLY & EXPORT
# ==========================================
cat("\nStep 3: Assembling and Exporting Figures...\n")

# Assemble Figure 1 (Publication Style Multi-panel)
figure_1 <- (p_pca | p_dot) / (p_loli | p_ridge) + plot_layout(widths = c(1, 1.5))

ggsave(file.path(plot_path, "Figure_1_Composite.pdf"), figure_1, width = 180, height = 200, units = "mm")
ggsave(file.path(plot_path, "Figure_1_Composite.svg"), figure_1, width = 180, height = 200, units = "mm")

# Individual Tables
write.csv(all_lvl2, file.path(table_path, "Full_Results_Lvl2.csv"))
write.csv(gene_map %>% group_by(Group, CellType) %>% slice_max(Score, n=10), 
          file.path(table_path, "Top_Driver_Genes_per_CellType.csv"))

# Save Workspace
saveRDS(list(lvl1=all_lvl1, lvl2=all_lvl2, genes=gene_map), file.path(qc_path, "Final_Analysis_Objects.rds"))

cat("\nAnalysis Finished Successfully.\n")
cat("Session Information:\n")
print(sessionInfo())

# Close log
sink()
sink()



















# ==========================================
# 1. SETUP, LOGGING & LIBRARIES
# ==========================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("readxl", "dplyr", "tibble", "tidyr", "ggplot2", "pheatmap", 
              "svglite", "ggridges", "ggrepel", "ggsci", "viridis", 
              "UpSetR", "openxlsx", "stringr", "Rtsne", "patchwork")

invisible(lapply(packages, function(x) if (!require(x, character.only = TRUE)) install.packages(x)))

library(EWCE); library(ewceData); library(readxl); library(dplyr)
library(tibble); library(tidyr); library(ggplot2); library(pheatmap)
library(openxlsx); library(stringr); library(ggridges); library(ggrepel)
library(ggsci); library(Rtsne); library(UpSetR); library(patchwork)

# Define and create output structure
base_results <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/EWCE"
plot_path    <- file.path(base_results, "Plots")
table_path   <- file.path(base_results, "Tables")
qc_path      <- file.path(base_results, "QC_Logs")
lapply(c(plot_path, table_path, qc_path), function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

# Initialize Log File
log_con <- file(file.path(qc_path, paste0("Run_Log_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt")), open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

cat("========================================================\n")
cat("EWCE PUBLICATION-STYLE ANALYSIS PIPELINE\n")
cat("Run Started:", as.character(Sys.time()), "\n")
cat("========================================================\n\n")

# ==========================================
# 2. DATA CLEANING & SPECIES FILTERING
# ==========================================
message("Step 1: Data Cleaning...")
ctd <- ewceData::ctd()
ref_genes <- rownames(ctd[[1]]$specificity)

path_to_file <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/imputed_data.xlsx"
raw_df <- read_excel(path_to_file)

# SPECIES FILTER: Only _MOUSE + Fix Case (GAD2 -> Gad2)
clean_df <- raw_df %>%
  filter(!is.na(Genes), Genes != "0") %>%
  filter(grepl("_MOUSE", Gene_Group)) %>% 
  mutate(Genes = str_to_title(Genes)) %>%  
  distinct(Genes, .keep_all = TRUE)

mapped_genes <- intersect(clean_df$Genes, ref_genes)
dropped_genes <- setdiff(clean_df$Genes, ref_genes)
contaminants <- raw_df %>% filter(!grepl("_MOUSE", Gene_Group))

cat(sprintf("Total Input Rows: %d\nMouse Genes Cleaned: %d\nSuccessfully Mapped to Atlas: %d\nDropped (Non-Atlas): %d\n", 
            nrow(raw_df), nrow(clean_df), length(mapped_genes), length(dropped_genes)))

# Create Grouped Matrix for top-gene selection
exp_matrix_grouped <- clean_df %>%
  filter(Genes %in% mapped_genes) %>%
  select(Genes, matches("neuron|mcherry|cfos|bg")) %>%
  pivot_longer(cols = -Genes, names_to = "Raw_Sample", values_to = "Expression") %>%
  mutate(Grouped_Name = gsub("\\.\\.\\.[0-9]+$", "", Raw_Sample)) %>%
  group_by(Genes, Grouped_Name) %>%
  summarise(Mean_Exp = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Grouped_Name, values_from = Mean_Exp) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

# ==========================================
# 3. EWCE ANALYSIS LOOP
# ==========================================
message("Step 2: Bootstrap Enrichment...")
all_results_lvl1 <- data.frame(); all_results_lvl2 <- data.frame()
gene_specificity_map <- data.frame()
group_names <- colnames(exp_matrix_grouped)

for (gname in group_names) {
  cat("Processing:", gname, "... ")
  hits <- rownames(exp_matrix_grouped)[order(exp_matrix_grouped[, gname], decreasing = TRUE)[1:250]] # Use top 250 markers for compact reporting
  curr_cond <- gsub("_[0-9]+.*$", "", gname)
  
  # Enrichment (Level 1 and 2)
  res_l1 <- bootstrap_enrichment_test(ctd, hits=hits, reps=1000, annotLevel=1, genelistSpecies="mouse", sctSpecies="mouse")
  all_results_lvl1 <- rbind(all_results_lvl1, res_l1$results %>% mutate(Group = gname, Condition = curr_cond))
  
  res_l2 <- bootstrap_enrichment_test(ctd, hits=hits, reps=1000, annotLevel=2, genelistSpecies="mouse", sctSpecies="mouse")
  all_results_lvl2 <- rbind(all_results_lvl2, res_l2$results %>% mutate(Group = gname, Condition = curr_cond))
  
  # Spec mapping for bubble/ridge
  gene_spec <- as.data.frame(ctd[[1]]$specificity[intersect(hits, ref_genes), ]) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "CellType", values_to = "SpecificityScore") %>%
    mutate(Group = gname, Condition = curr_cond)
  gene_specificity_map <- rbind(gene_specificity_map, gene_spec)
  cat("Done.\n")
}

all_results_lvl2 <- all_results_lvl2 %>% group_by(Group) %>% mutate(q = p.adjust(p, method = "BH")) %>% ungroup()

# ==========================================
# 4. DIM REDUCTION & STATS
# ==========================================
z_matrix <- all_results_lvl1 %>% select(Group, CellType, sd_from_mean) %>% 
  pivot_wider(names_from = Group, values_from = sd_from_mean) %>% column_to_rownames("CellType")

# PCA
pca_df <- as.data.frame(prcomp(t(z_matrix), scale. = TRUE)$x) %>% 
  rownames_to_column("Group") %>% mutate(Condition = gsub("_[0-9]+.*$", "", Group))

# t-SNE
set.seed(42)
tsne_df <- as.data.frame(Rtsne(t(z_matrix), perplexity = 1, check_duplicates = FALSE)$Y) %>% 
  rename(tSNE1=V1, tSNE2=V2) %>% mutate(Group = colnames(z_matrix), Condition = gsub("_[0-9]+.*$", "", Group))

# Consensus
consensus <- all_results_lvl1 %>% group_by(Condition, CellType) %>% 
  summarise(Mean_Z = mean(sd_from_mean), SE = sd(sd_from_mean)/sqrt(n()), .groups = "drop")

# ==========================================
# 5. PUBLICATION VISUALS (Patchwork Assembly)
# ==========================================
theme_publication <- function() {
  theme_classic(base_size = 7) + 
    theme(text = element_text(family = "sans", color = "black"),
          axis.text = element_text(size = 6), axis.title = element_text(size = 7, face = "bold"),
          strip.text = element_text(size = 7, face = "bold"), legend.title = element_text(size = 6, face = "bold"))
}

# Panel A: PCA
p_pca <- ggplot(pca_df, aes(PC1, PC2, color=Condition)) + geom_point(size=2) + 
  geom_text_repel(aes(label=Group), size=1.5) + scale_color_npg() + theme_publication() + labs(title="A: Profile Similarity")

# Panel B: Identity DotPlot
p_dot <- ggplot(all_results_lvl1 %>% filter(grepl("interneurons|pyramidal|astrocyte|microglia", CellType)), 
                aes(x = Group, y = reorder(CellType, sd_from_mean), size = -log10(p), color = sd_from_mean)) +
  geom_point() + scale_color_distiller(palette = "RdBu") + theme_publication() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title="B: Broad Enrichment")

# Panel C: Volcano
p_volc <- ggplot(all_results_lvl2, aes(sd_from_mean, -log10(p), color=Condition)) +
  geom_point(alpha=0.5, size=1) + geom_hline(yintercept = -log10(0.05), linetype="dashed", size=0.2) +
  geom_text_repel(data=subset(all_results_lvl2, q<0.05 & sd_from_mean>4), aes(label=CellType), size=1.5) +
  scale_color_npg() + theme_publication() + labs(title="C: Subtype Significance")

# Panel D: Spec Ridge
p_ridge <- ggplot(gene_specificity_map %>% filter(SpecificityScore > 0.1), aes(SpecificityScore, Condition, fill=Condition)) +
  geom_density_ridges(alpha=0.5) + scale_fill_npg() + theme_publication() + labs(title="D: Marker Specificity")

# COMBINE INTO FIGURE 1
main_fig <- (p_pca | p_dot) / (p_volc | p_ridge) + plot_annotation(tag_levels = 'a')

# ==========================================
# 6. EXPORT & FINAL LOGGING
# ==========================================
message("Step 4: Exporting...")

# Save Figures
ggsave(file.path(plot_path, "Main_Figure_1.pdf"), main_fig, width = 180, height = 180, units = "mm")
ggsave(file.path(plot_path, "Fig10_tSNE.svg"), p_pca, width = 100, height = 90, units = "mm")

# Save Tables
write.csv(all_results_lvl2, file.path(table_path, "Full_Results_Lvl2.csv"))
write.csv(consensus, file.path(table_path, "Condition_Consensus.csv"))

# Final QC Workbook
wb <- createWorkbook()
addWorksheet(wb, "Summary"); writeData(wb, 1, data.frame(Metric=c("Mapped", "Dropped"), Count=c(length(mapped_genes), length(dropped_genes))))
addWorksheet(wb, "DroppedGenes"); writeData(wb, 2, data.frame(Gene=dropped_genes))
saveWorkbook(wb, file.path(qc_path, "Analysis_QC_Final.xlsx"), overwrite = TRUE)

cat("\nAnalysis Complete. Session Info:\n")
print(sessionInfo())
sink(); sink() # Close logging
message("SUCCESS: Check Results folder.")












# ==========================================
# 1. SETUP, LOGGING & LIBRARIES
# ==========================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("readxl", "dplyr", "tibble", "tidyr", "ggplot2", "pheatmap", 
              "svglite", "ggridges", "ggrepel", "ggsci", "viridis", 
              "UpSetR", "openxlsx", "stringr", "Rtsne", "patchwork")

invisible(lapply(packages, function(x) if (!require(x, character.only = TRUE)) install.packages(x)))

library(EWCE); library(ewceData); library(readxl); library(dplyr)
library(tibble); library(tidyr); library(ggplot2); library(pheatmap)
library(openxlsx); library(stringr); library(ggridges); library(ggrepel)
library(ggsci); library(Rtsne); library(UpSetR); library(patchwork); library(svglite)

# Define and create output structure
base_results <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/EWCE"
plot_path    <- file.path(base_results, "Plots")
table_path   <- file.path(base_results, "Tables")
qc_path      <- file.path(base_results, "QC_Logs")
lapply(c(plot_path, table_path, qc_path), function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

# Initialize Log File
log_con <- file(file.path(qc_path, paste0("Run_Log_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt")), open = "wt")
sink(log_con, type = "output"); sink(log_con, type = "message")

cat("========================================================\n")
cat("EWCE PUBLICATION-STYLE COMPREHENSIVE PIPELINE\n")
cat("Run Started:", as.character(Sys.time()), "\n")
cat("========================================================\n\n")

# ==========================================
# 2. DATA CLEANING & SPECIES FILTERING
# ==========================================
message("Step 1: Data Cleaning...")
ctd <- ewceData::ctd()
ref_genes <- rownames(ctd[[1]]$specificity)

path_to_file <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/imputed_data.xlsx"
raw_df <- read_excel(path_to_file)

# SPECIES FILTER: Keep only _MOUSE and sanitize symbols (e.g., KRT10 -> Krt10)
clean_df <- raw_df %>%
  filter(!is.na(Genes), Genes != "0") %>%
  filter(grepl("_MOUSE", Gene_Group)) %>% 
  mutate(Genes = str_to_title(Genes)) %>%  
  distinct(Genes, .keep_all = TRUE)

mapped_genes  <- intersect(clean_df$Genes, ref_genes)
dropped_genes <- setdiff(clean_df$Genes, ref_genes)
contaminants  <- raw_df %>% filter(!grepl("_MOUSE", Gene_Group))

cat(sprintf("Mapping Stats:\n- Raw Input: %d\n- Mouse Cleaned: %d\n- Atlas Mapped: %d\n- Dropped: %d\n", 
            nrow(raw_df), nrow(clean_df), length(mapped_genes), length(dropped_genes)))

# Create Grouped Matrix
exp_matrix_grouped <- clean_df %>%
  filter(Genes %in% mapped_genes) %>%
  select(Genes, matches("neuron|mcherry|cfos|bg")) %>%
  pivot_longer(cols = -Genes, names_to = "Raw_Sample", values_to = "Expression") %>%
  mutate(Grouped_Name = gsub("\\.\\.\\.[0-9]+$", "", Raw_Sample)) %>%
  group_by(Genes, Grouped_Name) %>%
  summarise(Mean_Exp = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Grouped_Name, values_from = Mean_Exp) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

# ==========================================
# 3. EWCE ANALYSIS LOOP
# ==========================================
library(future); library(furrr)
plan(multisession, workers = availableCores() - 1)

message("Step 2: Bootstrap Enrichment (10000 reps, parallelized)...")
all_results_lvl1 <- data.frame(); all_results_lvl2 <- data.frame()
gene_specificity_map <- data.frame()
group_names <- colnames(exp_matrix_grouped)

# Parallel processing with future_map
results_list <- future_map(group_names, function(gname) {
  hits <- rownames(exp_matrix_grouped)[order(exp_matrix_grouped[, gname], decreasing = TRUE)[1:250]]
  curr_cond <- gsub("_[0-9]+.*$", "", gname)
  
  res_l1 <- bootstrap_enrichment_test(ctd, hits=hits, reps=10000, annotLevel=1, genelistSpecies="mouse", sctSpecies="mouse")
  lvl1_data <- res_l1$results %>% mutate(Group = gname, Condition = curr_cond)
  
  res_l2 <- bootstrap_enrichment_test(ctd, hits=hits, reps=10000, annotLevel=2, genelistSpecies="mouse", sctSpecies="mouse")
  lvl2_data <- res_l2$results %>% mutate(Group = gname, Condition = curr_cond)
  
  gene_spec <- as.data.frame(ctd[[1]]$specificity[intersect(hits, ref_genes), ]) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "CellType", values_to = "SpecificityScore") %>%
    mutate(Group = gname, Condition = curr_cond)
  
  list(lvl1 = lvl1_data, lvl2 = lvl2_data, spec = gene_spec)
}, .progress = TRUE)

# Combine results
all_results_lvl1 <- do.call(rbind, lapply(results_list, `[[`, "lvl1"))
all_results_lvl2 <- do.call(rbind, lapply(results_list, `[[`, "lvl2"))
gene_specificity_map <- do.call(rbind, lapply(results_list, `[[`, "spec"))

all_results_lvl2 <- all_results_lvl2 %>% group_by(Group) %>% mutate(q = p.adjust(p, method = "BH")) %>% ungroup()

plan(sequential)  # Reset to sequential

# ==========================================
# 4. DIM REDUCTION & STATS
# ==========================================
z_matrix <- all_results_lvl1 %>% select(Group, CellType, sd_from_mean) %>% 
  pivot_wider(names_from = Group, values_from = sd_from_mean) %>% column_to_rownames("CellType")

# PCA & t-SNE
pca_res <- prcomp(t(z_matrix), scale. = TRUE)
pca_df  <- as.data.frame(pca_res$x) %>% rownames_to_column("Group") %>% mutate(Condition = gsub("_[0-9]+.*$", "", Group))

set.seed(42)
tsne_res <- Rtsne(t(z_matrix), perplexity = 1, check_duplicates = FALSE)
tsne_df  <- as.data.frame(tsne_res$Y) %>% rename(tSNE1=V1, tSNE2=V2) %>% 
  mutate(Group = colnames(z_matrix), Condition = gsub("_[0-9]+.*$", "", Group))

# Consensus
consensus <- all_results_lvl1 %>% group_by(Condition, CellType) %>% 
  summarise(Mean_Z = mean(sd_from_mean), SE = sd(sd_from_mean)/sqrt(n()), .groups = "drop")

# ==========================================
# 5. PUBLICATION-STYLE PLOTTING ENGINE
# ==========================================
theme_publication <- function() {
  theme_classic(base_size = 8) + 
    theme(text = element_text(family = "sans", color = "black"),
          axis.text = element_text(size = 7), axis.title = element_text(size = 8, face = "bold"),
          strip.text = element_text(size = 8, face = "bold"), legend.title = element_text(size = 7, face = "bold"))
}

# Individual Plot Definitions
p_pca  <- ggplot(pca_df, aes(PC1, PC2, color=Condition)) + geom_point(size=3) + geom_text_repel(aes(label=Group), size=2) + scale_color_npg() + theme_publication() + labs(title="PCA: Sample Identity")
p_tsne <- ggplot(tsne_df, aes(tSNE1, tSNE2, color=Condition)) + geom_point(size=3) + geom_text_repel(aes(label=Group), size=2) + scale_color_npg() + theme_publication() + labs(title="t-SNE: Clustering")
p_dot  <- ggplot(all_results_lvl1 %>% filter(grepl("neuron|astrocyte|microglia", CellType, ignore.case = T)), 
                aes(x = Group, y = reorder(CellType, sd_from_mean), size = -log10(p), color = sd_from_mean)) +
          geom_point() + scale_color_distiller(palette = "RdBu") + theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title="Broad Identity")
p_volc <- ggplot(all_results_lvl2, aes(sd_from_mean, -log10(p), color=Condition)) +
          geom_point(alpha=0.5) + geom_hline(yintercept = -log10(0.05), linetype="dashed", linewidth=0.3) +
          geom_text_repel(data=subset(all_results_lvl2, q<0.05 & sd_from_mean>3.5), aes(label=CellType), size=1.8) +
          scale_color_npg() + theme_publication() + labs(title="Subtype Volcano")
p_ridge <- ggplot(gene_specificity_map %>% filter(SpecificityScore > 0.1), aes(SpecificityScore, Condition, fill=Condition)) +
          geom_density_ridges(alpha=0.6) + scale_fill_npg() + theme_publication() + labs(title="Marker Distribution")

# Create Main Figure Composite
main_fig <- (p_pca | p_dot) / (p_volc | p_ridge) + plot_annotation(tag_levels = 'A')

# ==========================================
# 6. COMPREHENSIVE EXPORT
# ==========================================
message("Step 4: Full Export...")

# 1. Save Composite (PDF & SVG)
ggsave(file.path(plot_path, "Main_Figure_1_Composite.pdf"), main_fig, width = 180, height = 180, units = "mm")
ggsave(file.path(plot_path, "Main_Figure_1_Composite.svg"), main_fig, width = 180, height = 180, units = "mm")

# 2. Save Individual SVGs for every plot
ggsave(file.path(plot_path, "Indiv_PCA.svg"), p_pca, width = 100, height = 100, units = "mm")
ggsave(file.path(plot_path, "Indiv_tSNE.svg"), p_tsne, width = 100, height = 100, units = "mm")
ggsave(file.path(plot_path, "Indiv_DotPlot.svg"), p_dot, width = 120, height = 120, units = "mm")
ggsave(file.path(plot_path, "Indiv_Volcano.svg"), p_volc, width = 120, height = 120, units = "mm")
ggsave(file.path(plot_path, "Indiv_Ridge.svg"), p_ridge, width = 120, height = 100, units = "mm")

# 3. Create Expanded Excel Workbook (The "Full Data" File)
wb <- createWorkbook()
addWorksheet(wb, "0_Analysis_Summary"); writeData(wb, 1, data.frame(Metric=c("Date", "Input Rows", "Mouse Genes", "Mapped to Atlas", "Dropped"), Value=c(as.character(Sys.Date()), nrow(raw_df), nrow(clean_df), length(mapped_genes), length(dropped_genes))))
addWorksheet(wb, "1_Level1_Broad"); writeData(wb, 2, all_results_lvl1)
addWorksheet(wb, "2_Level2_Granular"); writeData(wb, 3, all_results_lvl2)
addWorksheet(wb, "3_Condition_Consensus"); writeData(wb, 4, consensus)
addWorksheet(wb, "4_Top_Driver_Genes"); writeData(wb, 5, gene_specificity_map %>% group_by(Group, CellType) %>% slice_max(SpecificityScore, n=20) %>% ungroup())
addWorksheet(wb, "5_Dropped_Genes"); writeData(wb, 6, data.frame(Gene=dropped_genes))
addWorksheet(wb, "6_Contaminants_Purged"); writeData(wb, 7, contaminants)
saveWorkbook(wb, file.path(table_path, "EWCE_Full_Analysis_Results.xlsx"), overwrite = TRUE)

# 4. Save RDS for later R-session loading
saveRDS(list(lvl1=all_results_lvl1, lvl2=all_results_lvl2, pca=pca_df, spec=gene_specificity_map), 
        file.path(qc_path, "EWCE_Final_Objects.rds"))

cat("\nAnalysis Finished. Check the 'Plots', 'Tables', and 'QC_Logs' folders.\n")
print(sessionInfo())
sink(); sink()
