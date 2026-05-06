# ==========================================
# 1. SETUP, LIBRARIES & DIRECTORIES
# ==========================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("readxl", "dplyr", "tibble", "tidyr", "ggplot2", "pheatmap", 
              "svglite", "ggridges", "ggrepel", "ggsci", "viridis", 
              "openxlsx", "stringr", "patchwork", "future", "future.apply",
              "EWCE", "ewceData", "org.Mm.eg.db", "AnnotationDbi" )

invisible(lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    if (x %in% c("EWCE", "ewceData", "org.Mm.eg.db", "AnnotationDbi")) {
      BiocManager::install(x, update = FALSE)
    } else {
      install.packages(x)
    }
  }
  library(x, character.only = TRUE)
}))

# --- Improved Folder Structure ---
data_path    <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/pg_matrix/imputed/20260218_pgmatrix_imputed_neuron_soma_71samples_missing70pct_groups.xlsx"
base_results <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/EWCE/neuron_soma"

# Logical Sub-folders
dirs <- list(
  plots   = file.path(base_results, "01_Figures_Main"),
  svgs    = file.path(base_results, "01_Figures_Main/SVG_Editable"),
  tables  = file.path(base_results, "02_Tables_Supplements"),
  qc      = file.path(base_results, "03_QC_Mapping_Logs"),
  data    = file.path(base_results, "04_Processed_Data_Objects")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# ==========================================
# 2. ROBUST GENE MAPPING & CLEANING
# ==========================================
message("Step 1: Robust Gene Mapping...")
ctd <- ewceData::ctd()
ref_genes <- rownames(ctd[[1]]$specificity)

raw_df <- read_excel(data_path)

clean_df <- raw_df %>%
  filter(!is.na(Genes), Genes != "0") %>%
  mutate(Genes = str_split_i(Genes, ";", 1))

mapped_ids <- mapIds(org.Mm.eg.db, keys = clean_df$Genes, column = "SYMBOL", keytype = "ALIAS", multiVals = "first")
clean_df$Gene <- as.character(mapped_ids)

clean_df <- clean_df %>% filter(!is.na(Gene)) %>% distinct(Gene, .keep_all = TRUE)
mapped_genes <- intersect(clean_df$Gene, ref_genes)
dropped_genes <- setdiff(clean_df$Gene, ref_genes)

# Save Mapping QC
write.table(dropped_genes, file.path(dirs$qc, "dropped_genes_not_in_CTD.txt"), row.names = FALSE, col.names = FALSE)

# ==========================================
# 3. GENERATE TARGET SIGNATURES
# ==========================================
sample_cols <- colnames(clean_df)[grep("^D:", colnames(clean_df))]

long_df <- clean_df %>%
  filter(Gene %in% mapped_genes) %>%
  dplyr::select(Gene, all_of(sample_cols)) %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Exp") %>%
  mutate(Region = str_extract(Sample, "CA1|CA2|CA3|DG"),
         Cond   = str_extract(tolower(Sample), "con|res|sus"))

# Baseline Abundance
baseline_mat <- long_df %>%
  group_by(Gene, Region, Cond) %>%
  summarise(MeanExp = mean(Exp, na.rm = TRUE), .groups = "drop") %>%
  mutate(ID = paste(Region, Cond, sep = "_")) %>%
  dplyr::select(-Region, -Cond) %>%
  pivot_wider(names_from = ID, values_from = MeanExp) %>%
  column_to_rownames("Gene")

# Differential Expression (LogFC style)
diff_list <- long_df %>%
  group_by(Gene, Region, Cond) %>%
  summarise(MeanExp = mean(Exp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cond, values_from = MeanExp) %>%
  mutate(Sus_up = sus - con, Res_up = res - con)

# ==========================================
# 4. PARALLEL BOOTSTRAPPING
# ==========================================
message("Step 2: Bootstrapping 10,000 reps...")
plan(multisession, workers = parallel::detectCores() - 1)

all_targets <- c(colnames(baseline_mat), 
                  paste0(rep(unique(long_df$Region), each=2), c("_Sus_up", "_Res_up")))

run_ewce_core <- function(target) {
  if(grepl("_up", target)) {
    reg <- str_split_i(target, "_", 1)
    col <- paste0(str_split_i(target, "_", 2), "_up")
    hits <- diff_list %>% filter(Region == reg) %>% arrange(desc(.data[[col]])) %>% slice_head(n=250) %>% pull(Gene)
    type <- "Differential"
  } else {
    hits <- rownames(baseline_mat)[order(baseline_mat[,target], decreasing = TRUE)[1:250]]
    type <- "Baseline"
  }
  
  res <- EWCE::bootstrap_enrichment_test(ctd, hits=hits, reps=10000, annotLevel=2, genelistSpecies="mouse", sctSpecies="mouse")
  res$results %>% mutate(Target = target, AnalysisType = type, 
                         Region = str_extract(target, "CA1|CA2|CA3|DG"),
                         Metric = str_remove(target, "^(CA1|CA2|CA3|DG)_"))
}

results_all <- future_lapply(all_targets, run_ewce_core, future.seed = TRUE) %>% bind_rows()
results_all <- results_all %>% group_by(Target) %>% mutate(q = p.adjust(p, method="BH")) %>% ungroup()
results_all <- results_all %>%
  mutate(
    neglog10q = -log10(pmax(q, 1e-300)),
    SignedSig = sign(sd_from_mean) * neglog10q,
    Direction = case_when(
      sd_from_mean > 0 ~ "Positive",
      sd_from_mean < 0 ~ "Negative",
      TRUE ~ "Neutral"
    )
  )
plan(sequential)

# ==========================================
# 5. NATURE STYLE VISUALIZATION
# ==========================================

# Standardized Nature Theme
theme_nature <- function() {
  theme_bw(base_size = 7) + # Nature uses 5-7pt for labels
    theme(
      text = element_text(family = "sans", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      axis.line = element_blank(), 
      axis.ticks = element_line(linewidth = 0.2, color = "black"),
      axis.text = element_text(size = 6, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 6, face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0)
    )
}

col_baseline <- "#7ea6a9" 
col_sus      <- "#e63f4b"
col_res      <- "#3f3d6f"

# P1: Dotplot Baseline
p1 <- ggplot(results_all %>% filter(AnalysisType == "Baseline", q < 0.05), 
             aes(x = Metric, y = CellType, size = sd_from_mean, color = sd_from_mean)) +
  geom_point() + 
  facet_wrap(~Region, nrow = 1) + 
  scale_color_viridis_c(option = "magma", name = "Z-score") +
  scale_size_area(max_size = 3, name = "SD") +
  theme_nature() + 
  labs(x = NULL, y = "Cell Type", title = "Baseline Cell Enrichment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# P2: Heatmap Differential
p2 <- ggplot(results_all %>% filter(AnalysisType == "Differential"), 
             aes(x = Metric, y = CellType, fill = sd_from_mean)) +
  geom_tile(color = "white", linewidth = 0.1) + 
  facet_wrap(~Region, nrow = 1) + 
  scale_fill_gradient2(low = col_res, mid = "white", high = col_sus, 
                       midpoint = 0, name = "ΔZ") +
  theme_nature() + 
  labs(x = NULL, y = "Cell Type", title = "Stress-Induced Changes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# P3: Volcano
p3 <- ggplot(results_all, aes(x = sd_from_mean, y = -log10(p))) +
  geom_point(aes(color = AnalysisType), alpha = 0.6, size = 1) + 
  scale_color_manual(values = c("Baseline" = col_baseline, "Differential" = col_sus)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2) +
  theme_nature() + 
  labs(x = "Effect Size (Z-score)", y = "-log10(p-value)", title = "Enrichment Significance")

# P4: Signed Significance Heatmap
p4 <- ggplot(results_all %>% filter(AnalysisType == "Differential"),
             aes(x = Metric, y = CellType, fill = SignedSig)) +
  geom_tile(color = "white", linewidth = 0.1) +
  facet_wrap(~Region, nrow = 1) +
  scale_fill_gradient2(low = col_res, mid = "white", high = col_sus,
                       midpoint = 0, name = "Signed -log10(q)") +
  theme_nature() +
  labs(x = NULL, y = "Cell Type", title = "Direction + Significance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# P5: Top Cell Types by Max Absolute Effect
top_celltypes <- results_all %>%
  group_by(CellType) %>%
  summarise(max_abs_effect = max(abs(sd_from_mean), na.rm = TRUE),
            best_q = min(q, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_effect)) %>%
  slice_head(n = 20)

p5 <- ggplot(top_celltypes,
             aes(x = reorder(CellType, max_abs_effect), y = max_abs_effect,
                 fill = -log10(pmax(best_q, 1e-300)))) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_viridis_c(option = "plasma", name = "-log10(q)") +
  theme_nature() +
  labs(x = NULL, y = "Max |Z-score|", title = "Top Cell-Type Effects")

# P6: Effect Size Distributions
p6 <- ggplot(results_all %>% filter(AnalysisType == "Differential"),
             aes(x = sd_from_mean, y = Region, fill = Region)) +
  ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9, color = "white", linewidth = 0.2) +
  scale_fill_viridis_d(option = "D", guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2) +
  theme_nature() +
  labs(x = "Z-score", y = NULL, title = "Regional Effect Distributions")

# Assemble
main_fig <- (p1 / p2 | p3) + 
  plot_layout(widths = c(2, 1)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold", size = 10))

supp_fig <- (p4 / p5 / p6) +
  plot_layout(heights = c(1.1, 1, 0.8)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 10))

# Clustered Matrix Heatmap (publication-ready supplement)
diff_heatmap_tbl <- results_all %>%
  filter(AnalysisType == "Differential") %>%
  mutate(TargetLabel = paste(Region, Metric, sep = "_")) %>%
  dplyr::select(CellType, TargetLabel, SignedSig) %>%
  distinct() %>%
  pivot_wider(names_from = TargetLabel, values_from = SignedSig)

diff_heatmap_mat <- diff_heatmap_tbl %>%
  column_to_rownames("CellType") %>%
  as.matrix()

if (nrow(diff_heatmap_mat) > 2 && ncol(diff_heatmap_mat) > 2) {
  pheatmap::pheatmap(
    diff_heatmap_mat,
    color = colorRampPalette(c(col_res, "white", col_sus))(101),
    border_color = NA,
    fontsize = 6,
    fontsize_row = 5,
    fontsize_col = 6,
    angle_col = 45,
    main = "Signed Significance by Target"
  )
}

# ==========================================
# 6. EXPORTING
# ==========================================
message("Step 3: Exporting Files...")

# Save Plots
ggsave(file.path(dirs$plots, "Fig1_EWCE_Summary.pdf"), main_fig, width = 180, height = 200, units = "mm", device = cairo_pdf)
ggsave(file.path(dirs$svgs, "Fig1_EWCE_Summary.svg"), main_fig, width = 180, height = 200, units = "mm")
ggsave(file.path(dirs$svgs, "Volcano_Panel.svg"), p3, width = 80, height = 80, units = "mm")
ggsave(file.path(dirs$plots, "FigS1_EWCE_Additional_Visuals.pdf"), supp_fig, width = 180, height = 240, units = "mm", device = cairo_pdf)
ggsave(file.path(dirs$svgs, "FigS1_EWCE_Additional_Visuals.svg"), supp_fig, width = 180, height = 240, units = "mm")

if (nrow(diff_heatmap_mat) > 2 && ncol(diff_heatmap_mat) > 2) {
  pdf(file.path(dirs$plots, "FigS2_EWCE_ClusteredHeatmap.pdf"), width = 8, height = 10, family = "sans")
  pheatmap::pheatmap(
    diff_heatmap_mat,
    color = colorRampPalette(c(col_res, "white", col_sus))(101),
    border_color = NA,
    fontsize = 6,
    fontsize_row = 5,
    fontsize_col = 6,
    angle_col = 45,
    main = "Signed Significance by Target"
  )
  dev.off()
}

# Complex Table Export
wb <- createWorkbook()
addWorksheet(wb, "Full_Results"); writeDataTable(wb, "Full_Results", results_all)
addWorksheet(wb, "Significant_Results"); writeDataTable(wb, "Significant_Results", results_all %>% filter(q < 0.05))
addWorksheet(wb, "Input_Gene_Stats"); writeDataTable(wb, "Input_Gene_Stats", data.frame(Total=nrow(clean_df), Mapped=length(mapped_genes), Dropped=length(dropped_genes)))

summary_tbl <- results_all %>%
  group_by(AnalysisType, Region, Metric) %>%
  summarise(
    N_Tested = n(),
    N_Significant = sum(q < 0.05, na.rm = TRUE),
    Median_Abs_Effect = median(abs(sd_from_mean), na.rm = TRUE),
    Mean_Abs_Effect = mean(abs(sd_from_mean), na.rm = TRUE),
    Top_CellType = CellType[which.max(abs(sd_from_mean))],
    Top_Effect = max(abs(sd_from_mean), na.rm = TRUE),
    Best_q = min(q, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(AnalysisType, Region, desc(N_Significant), desc(Mean_Abs_Effect))

top_hits_tbl <- results_all %>%
  group_by(Target) %>%
  arrange(q, desc(abs(sd_from_mean)), .by_group = TRUE) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  dplyr::select(Target, AnalysisType, Region, Metric, CellType, sd_from_mean, p, q, SignedSig)

sig_rank_tbl <- results_all %>%
  mutate(Significant = q < 0.05) %>%
  group_by(CellType, AnalysisType) %>%
  summarise(
    Significant_Count = sum(Significant, na.rm = TRUE),
    Mean_SignedSig = mean(SignedSig, na.rm = TRUE),
    Max_Abs_Effect = max(abs(sd_from_mean), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Significant_Count), desc(Max_Abs_Effect))

addWorksheet(wb, "Summary_by_Contrast"); writeDataTable(wb, "Summary_by_Contrast", summary_tbl)
addWorksheet(wb, "Top_Hits_per_Target"); writeDataTable(wb, "Top_Hits_per_Target", top_hits_tbl)
addWorksheet(wb, "CellType_Significance_Rank"); writeDataTable(wb, "CellType_Significance_Rank", sig_rank_tbl)

saveWorkbook(wb, file.path(dirs$tables, "Supplementary_Table_EWCE.xlsx"), overwrite = TRUE)

# Data Object for R
saveRDS(results_all, file.path(dirs$data, "EWCE_results_full.rds"))

cat("\nPipeline Complete. Files organized in:", base_results, "\n")