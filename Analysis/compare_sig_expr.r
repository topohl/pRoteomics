#' @title Compare Proteomics Expression to Learning Signatures
#' 
#' @description
#' This script analyzes proteomics expression data to evaluate the concordance between 
#' an experimental treatment (specifically CNO vs. Vehicle) and established "learning" 
#' molecular signatures. The pipeline processes raw expression data, integrates extenal 
#' differential expression statistics, and visualizes the correlation between the 
#' experiment and the reference signature.
#' 
#' @section Workflow:
#' 1. **Setup**: Loads required libraries (pheatmap, ggplot2, etc.) and defines input/output paths.
#' 2. **Data Ingestion**: 
#'    - Parses a GCT file to extract protein expression matrix and sample metadata.
#'    - Filters samples based on user-defined cell type (`COMP_TYPE`) and standardizes metadata terms.
#' 3. **Signature Integration**: 
#'    - Aggregates external Excel files containing "learning" signature statistics (Up/Down regulated lists).
#'    - Aligns the signature data with the experimental dataset via UniProt IDs.
#' 4. **Statistical Analysis**:
#'    - Calculates raw Log2 Fold Change (Log2FC) of CNO vs. Vehicle in the current experiment.
#'    - computes global Z-scores for these Log2FC values.
#'    - Annotates proteins based on directionality concordance (Learning UP/DOWN vs CNO UP/DOWN).
#' 5. **Visualization**: Generates vector graphics and PDFs, including:
#'    - Ordered heatmaps of Z-scores.
#'    - Per-sample and group-averaged expression heatmaps.
#'    - Side-by-side comparison heatmaps (Reference vs. Experiment).
#'    - Correlation scatter plots with linear regression.
#' 
#' @param COMP_TYPE Input string (e.g., "mcherry") used to filter `celltype` in metadata and define output filenames.
#' @param GCT_FILE_PATH Path to the source .gct file containing unnormalized protein matrix and header metadata.
#' @param SIG_PATH Directory path containing reference Excel files for the learning signature (matching pattern `log2fc_.*_(down|up)regulated_all.xlsx`).
#' 
#' @return 
#' Writes output to `BASE_OUT_DIR` within `/tables` and `/plots`:
#' - **CSV Table**: Merged dataframe of fold-changes, Z-scores, and directional annotations.
#' - **PDF/SVG Plots**: Various heatmaps and a correlation scatter plot illustrating the relationship between the drug effect and the learning signature.
#' 
#' @dependencies pacman, readxl, dplyr, tidyr, pheatmap, svglite, ggplot2, RColorBrewer
#' @author Tobi
## =============================================================================
## TITLE: Compare Proteomics Expression to Learning Signatures
## =============================================================================

## ---------- 0. Libraries and Setup ----------
if (!require("pacman")) install.packages("pacman")
library(pacman)
# Added ggplot2 and RColorBrewer here to ensure they are loaded early
p_load(readxl, dplyr, tidyr, pheatmap, svglite, ggplot2, RColorBrewer)

## ---------- 1. Configuration (USER INPUT) ----------

# 1.1 Comparison Settings
# Define the primary comparison variable (e.g., "mcherry", "cfos").
# Determines cell type filtering and contrast name (pattern: "{comp_type}2_{comp_type}4").
COMP_TYPE <- "mcherry" 
TARGET_CONTRAST <- paste0(COMP_TYPE, "2_", COMP_TYPE, "4")

# 1.2 Input File Paths
# GCT File
GCT_FILE_PATH <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/pg.matrix_filtered_pcaAdjusted_unnormalized.gct"

# Folder containing differential expression results (up/down regulated lists)
SIG_PATH <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/compareGO/BP/learning_signature/memory_ensemble/core_enrichment/significant_proteins/"

# 1.3 Output Settings
BASE_OUT_DIR <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/compare_sig_expr"
DIR_TABLES   <- file.path(BASE_OUT_DIR, "tables")
DIR_PLOTS    <- file.path(BASE_OUT_DIR, "plots")

# Check directories
if (!dir.exists(DIR_TABLES)) dir.create(DIR_TABLES, recursive = TRUE)
if (!dir.exists(DIR_PLOTS))  dir.create(DIR_PLOTS,  recursive = TRUE)

## ---------- 2. Helper Functions ----------

# Parse GCT files: handle skipping metadata lines and standardizing column names.
read_gct <- function(file_path) {
  if (!file.exists(file_path)) stop("GCT file not found: ", file_path)
  
  gct_data <- read.delim(file_path, skip = 2, header = FALSE, check.names = FALSE)
  if ("Name" %in% colnames(gct_data)) {
    colnames(gct_data)[1] <- "Gene"
  }
  gct_data
}

## ---------- 3. Data Loaded: GCT File ----------

gct_data <- read_gct(GCT_FILE_PATH)

# Locate the start of the expression data (first UniProt ID ending in "_MOUSE").
num_start <- which(grepl("_MOUSE$", gct_data[[1]]))[1]

# Separate metadata and expression data.
meta_rows <- gct_data[1:(num_start - 1), ]
expr_raw  <- gct_data[num_start:nrow(gct_data), ]

# Construct Expression Matrix.
prot_ids  <- expr_raw[[1]]
expr_mat  <- as.matrix(expr_raw[,-1])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- prot_ids

# Construct Sample Annotation.
sample_anno <- as.data.frame(t(meta_rows[,-1]), stringsAsFactors = FALSE)
colnames(sample_anno) <- trimws(meta_rows[[1]])
sample_anno[] <- lapply(sample_anno, function(x) if (is.character(x)) trimws(x) else x)
sample_anno$sample_id <- colnames(expr_mat)

# Process 'ExpGroup' factor.
if (!"ExpGroup" %in% colnames(sample_anno)) {
  warning("Metadata 'ExpGroup' not found. Available: ", paste(colnames(sample_anno), collapse = ", "))
} else {
  sample_anno$ExpGroup <- factor(as.numeric(sample_anno$ExpGroup),
                                 levels = c(1, 2, 3, 4),
                                 labels = c("1", "2", "3", "4"))
}

# Filter by 'celltype' based on configuration.
if ("celltype" %in% colnames(sample_anno)) {
  keep_samples <- sample_anno$sample_id[sample_anno$celltype == COMP_TYPE]
  
  if (length(keep_samples) > 0) {
    sample_anno <- sample_anno[sample_anno$sample_id %in% keep_samples, ]
    expr_mat    <- expr_mat[, keep_samples, drop = FALSE]
    message(sprintf("Filtered to %d %s samples.", length(keep_samples), COMP_TYPE))
  } else {
    warning(sprintf("No samples found with celltype == '%s'.", COMP_TYPE))
  }
} else {
  warning("Column 'celltype' not found in metadata.")
}

## ---------- 4. Data Loading: Learning Signature ----------

# List and aggregate Excel files.
file_list <- list.files(path = SIG_PATH, pattern = "log2fc_.*_(down|up)regulated_all\\.xlsx$", full.names = TRUE)
file_names  <- tools::file_path_sans_ext(basename(file_list))
group_names <- gsub("^log2fc_|_(down|up)regulated_all$", "", file_names)

data_list <- list()
unique_groups <- unique(group_names)

for (grp in unique_groups) {
  current_files <- file_list[group_names == grp]
  data_list[[grp]] <- lapply(current_files, read_excel) %>% bind_rows()
}

if (!TARGET_CONTRAST %in% names(data_list)) {
  stop("Contrast '", TARGET_CONTRAST, "' not found. Available: ", paste(names(data_list), collapse = ", "))
}

# Clean and format the learning dataframe.
learn_df_clean <- data_list[[TARGET_CONTRAST]] %>%
  mutate(
    log2fc     = as.numeric(log2fc),
    abs_log2fc = as.numeric(abs_log2fc),
    padj       = as.numeric(padj),
    Rank       = as.numeric(Rank)
  ) %>%
  dplyr::select(Gene_Name, `UniProtKB-ID`, Uniprot_Accession, Gene_Synonym, 
                Comparison, log2fc, abs_log2fc, padj, Direction, Rank)

# Intersect UniProt IDs.
learn_ids <- intersect(unique(learn_df_clean$`UniProtKB-ID`), rownames(expr_mat))

## ---------- 5. Statistics: Calculate Differential Effect & Z-Scores ----------

# Define groups (1 = CNO, 2 = VEH).
cno_samples <- intersect(sample_anno$sample_id[sample_anno$ExpGroup == "1"], colnames(expr_mat))
veh_samples <- intersect(sample_anno$sample_id[sample_anno$ExpGroup == "2"], colnames(expr_mat))

if (length(cno_samples) == 0 || length(veh_samples) == 0) stop("No samples found for CNO or VEH groups.")

expr_cno <- expr_mat[, cno_samples, drop = FALSE]
expr_veh <- expr_mat[, veh_samples, drop = FALSE]

# A. Calculate Raw Fold Change
log2fc_cno_veh <- rowMeans(expr_cno, na.rm = TRUE) - rowMeans(expr_veh, na.rm = TRUE)

# B. Global Z-Score Calculation (All Proteins)
z_all <- (log2fc_cno_veh - mean(log2fc_cno_veh, na.rm = TRUE)) / sd(log2fc_cno_veh, na.rm = TRUE)

# C. Extract Signatures and Integreate
z_learning <- z_all[learn_ids]
learning_z_df <- data.frame(
  `UniProtKB-ID` = names(z_learning),
  z_cno_vs_veh   = as.numeric(z_learning),
  stringsAsFactors = FALSE, check.names = FALSE
)

learn_annot_z <- learn_df_clean %>%
  inner_join(learning_z_df, by = "UniProtKB-ID") %>%
  mutate(
    learning_dir = ifelse(log2fc > 0, "learning_up", "learning_down"),
    cno_dir      = ifelse(z_cno_vs_veh > 0, "cno_up", "cno_down"),
    same_direction = case_when(
      learning_dir == "learning_up"   & cno_dir == "cno_up"   ~ "same",
      learning_dir == "learning_down" & cno_dir == "cno_down" ~ "same",
      TRUE ~ "opposite"
    )
  )

# Save intermediate table
write.csv(learn_annot_z, file.path(DIR_TABLES, paste0("learning_signature_zscores_cno_vs_veh_", COMP_TYPE, ".csv")), row.names = FALSE)

## ---------- 6. Visualization Setup (Colors & Palettes) ----------

# Define colors for pheatmap annotations
ann_colors <- list()

# Helper to assign colors safely if the column exists/is not empty
unique_exp <- if(!all(is.na(sample_anno$ExpGroup))) unique(sample_anno$ExpGroup) else NULL
unique_ld  <- unique(learn_annot_z$learning_dir)
unique_cd  <- unique(learn_annot_z$cno_dir)
unique_sd  <- unique(learn_annot_z$same_direction)

if (!is.null(unique_exp)) ann_colors$ExpGroup       <- setNames(brewer.pal(max(2, length(unique_exp)), "Set1")[seq_along(unique_exp)], unique_exp)
if (!is.null(unique_ld))  ann_colors$learning_dir   <- setNames(c("#e31a1c", "#1f78b4")[seq_along(unique_ld)], unique_ld)
if (!is.null(unique_cd))  ann_colors$cno_dir        <- setNames(c("#fb9a99", "#a6cee3")[seq_along(unique_cd)], unique_cd)
if (!is.null(unique_sd))  ann_colors$same_direction <- setNames(c("#33a02c", "#e31a1c")[seq_along(unique_sd)], unique_sd)

heatmap_cols <- colorRampPalette(c("#6698CC", "white", "#F08C21"))(100)

## ---------- 7. Generate Plots ----------

## --- Plot A: 1-Column Heatmap of Z-Scores ---
ord <- order(learn_annot_z$learning_dir, -learn_annot_z$z_cno_vs_veh)
heatmap_mat1 <- matrix(learn_annot_z$z_cno_vs_veh[ord], ncol = 1, 
                       dimnames = list(learn_annot_z$`UniProtKB-ID`[ord], "Z-score (CNO vs VEH)"))

pheatmap(
  heatmap_mat1,
  cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
  main = paste0("Learning-up vs learning-down (CNO z-scores) - ", COMP_TYPE),
  filename = file.path(DIR_PLOTS, paste0("heatmap_learning_signature_zscores_1col_ordered_", COMP_TYPE, ".pdf")),
  width = 4, height = 8
)

## --- Plot B: Per-Sample Expression Heatmap ---
expr_learning <- expr_mat[learn_ids, c(cno_samples, veh_samples), drop = FALSE]
expr_learning_scaled <- t(scale(t(as.matrix(expr_learning)))) # Z-score scaling

# Annotations
col_anno <- data.frame(row.names = colnames(expr_learning_scaled),
                       ExpGroup = sample_anno$ExpGroup[match(colnames(expr_learning_scaled), sample_anno$sample_id)])

row_anno <- learn_annot_z %>% 
  dplyr::select(`UniProtKB-ID`, learning_dir, cno_dir, same_direction) %>% distinct() %>%
  as.data.frame()
rownames(row_anno) <- row_anno$`UniProtKB-ID`
row_anno <- row_anno[rownames(expr_learning_scaled), c("learning_dir", "cno_dir", "same_direction")]

pheatmap(
  expr_learning_scaled,
  annotation_col = col_anno, annotation_row = row_anno, annotation_colors = ann_colors,
  color = heatmap_cols, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE,
  main = paste0("Learning-signature proteins: paired CNO vs VEH (per sample) - ", COMP_TYPE),
  filename = file.path(DIR_PLOTS, paste0("heatmap_learning_signature_expr_per_sample_annotated_", COMP_TYPE, ".pdf")),
  width = 7, height = 8
)

## --- Plot C: Averaged Expression per ExpGroup ---
expr_learning_all <- expr_mat[learn_ids, , drop = FALSE]
group_vec <- sample_anno$ExpGroup[match(colnames(expr_learning_all), sample_anno$sample_id)]

# Calculate averages
avg_list <- lapply(sort(unique(group_vec)), function(g) {
  rowMeans(expr_learning_all[, colnames(expr_learning_all)[group_vec == g], drop = FALSE], na.rm = TRUE)
})
avg_expr <- do.call(cbind, avg_list)
colnames(avg_expr) <- paste0("Group_", sort(unique(group_vec)))
avg_expr_scaled <- t(scale(t(avg_expr)))

# Match row annotation to averages
row_anno_avg <- row_anno[rownames(avg_expr_scaled), , drop = FALSE]

pheatmap(
  avg_expr_scaled,
  annotation_row = row_anno_avg, annotation_colors = ann_colors,
  color = heatmap_cols, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE,
  main = paste0("Average Expression per ExpGroup - ", COMP_TYPE),
  filename = file.path(DIR_PLOTS, paste0("heatmap_learning_signature_avg_per_group_", COMP_TYPE, ".pdf")),
  width = 5, height = 8
)

## --- Plot D: Side-by-Side Heatmap (Learning vs CNO) ---
# Normalize Learning FC to Z-score
mu_learn  <- mean(learn_annot_z$log2fc, na.rm = TRUE)
sig_learn <- sd(learn_annot_z$log2fc,   na.rm = TRUE)
learn_annot_z$z_learning_fc <- (learn_annot_z$log2fc - mu_learn) / sig_learn

ord_learn <- order(learn_annot_z$z_learning_fc, decreasing = TRUE)
mat_pair <- rbind(
  learning_z   = learn_annot_z$z_learning_fc[ord_learn],
  cno_vs_veh_z = learn_annot_z$z_cno_vs_veh[ord_learn]
)
colnames(mat_pair) <- learn_annot_z$`UniProtKB-ID`[ord_learn]

limit <- max(abs(mat_pair), na.rm = TRUE)
pheatmap(
  mat_pair,
  color = heatmap_cols, breaks = seq(-limit, limit, length.out = 101),
  cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
  main = paste0("Learning z-score vs CNO z-score - ", COMP_TYPE),
  filename = file.path(DIR_PLOTS, paste0("heatmap_learning_vs_cno_signature_", COMP_TYPE, ".pdf")),
  width = 6, height = 3
)

## --- Plot E: Scatter Plot (Learning vs CNO Effect) ---
cor_val <- cor(learn_annot_z$z_learning_fc, learn_annot_z$z_cno_vs_veh, use = "complete.obs")

p_scatter <- ggplot(learn_annot_z, aes(x = z_learning_fc, y = z_cno_vs_veh)) +
  geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
  geom_point(color = "#2c3e50", alpha = 0.3, size = 4, shape = 16) +
  geom_smooth(method = "lm", color = "#e74c3c", fill = "#e74c3c", alpha = 0.2) +
  labs(
    x = paste0("Learning z-score (", TARGET_CONTRAST, ")"),
    y = "CNO vs VEH z-score (paired, group 1 vs 2)",
    title = paste0("Learning signature vs CNO effect - ", COMP_TYPE),
    subtitle = paste0("Pearson correlation: r = ", round(cor_val, 2))
  ) +
  theme_minimal(base_size = 22, base_family = "Arial Nova") +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 20)
  ) +
  coord_fixed(ratio = 1)

ggsave(file.path(DIR_PLOTS, paste0("scatter_learning_z_vs_cno_z_", COMP_TYPE, ".svg")),
       plot = p_scatter, width = 6, height = 6)

message("Analysis complete. Files saved to: ", BASE_OUT_DIR)
