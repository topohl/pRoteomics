#' Generate Z-Scored Expression Plots by Cell Type and Group
#'
#' This script reads an Excel dataset with expression data, filters out rows where
#' the "Celltype" column is labeled as "Background", and computes the z-scores of
#' TKNK_MOUSE expression values for each combination of Celltype and Group.
#'
#' The script performs the following steps:
#'   - Reads the dataset from an Excel file.
#'   - Filters out records with "Background" in the Celltype column.
#'   - Calculates the z-score for TKNK_MOUSE expression within each Celltype/Group grouping.
#'   - Iterates over unique groups to create plots that combine:
#'       * A boxplot to show the distribution of z-scored expression values.
#'       * A jitter plot to overlay individual data points.
#'   - Applies a custom color palette and minimal theme for clear visualization.
#'   - Saves each generated plot as an SVG file in the specified results directory.
#'
#' @note Ensure that the file paths for the Excel input and the output directory are correctly
#' set before running the script. Adjust the color palette and plot aesthetics as needed.


if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, readxl, tibble, tidyr, pheatmap)

# Read the data (adjust path and separator as needed), and remove "Background" celltype data
df <- read_excel("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/rawDataFC.xlsx") %>%
  filter(Celltype != "Background") %>%
  mutate(Group = factor(Group, levels = c(1,2,3,4),
                        labels = c("CS CNO", "CS VEH", "US CNO", "US VEH")))

results_dir <- "C:/Users/topohl/Documents/Rtestground"

# Optional: Check available group values
available_groups <- unique(df$Group)
print(available_groups)

# Generate z-score dataframe for TKNK_MOUSE expression per Celltype and Group
df_z <- df %>%
  group_by(Celltype, Group) %>%
  mutate(TKNK_zscore = scale(as.numeric(TKNK_MOUSE))) %>%
  ungroup()

# Loop over each group and generate a plot
for (grp in available_groups) {

  # Filter the z-scored data for current group
  df_grp_z <- df_z %>% filter(Group == grp)

  # Optional: Print cell type counts for the current group
  print(table(df_grp_z$Celltype))

  # Create the plot
  p <- ggplot(df_grp_z, aes(x = Celltype, y = TKNK_zscore, fill = Celltype)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_boxplot(outlier.shape = NA, color = "#4D4D4D", size = 0.5, width = 0.7, alpha = 0.8) +
    geom_jitter(width = 0.1, size = 5, alpha = 0.8, shape = 21, stroke = 0.5, aes(fill = Celltype), color = "#DFDFDF") +
    scale_fill_manual(values = c("#A8DADC", "#457B9D", "#E63946", "#FFB703", "#8ecae6")) +
    theme_minimal(base_family = "Helvetica", base_size = 14) +
    labs(
      title = paste("TKNK_MOUSE Expression (z-score) Across Cell Types (Group", grp, ")"),
      y = "TKNK_MOUSE (z-score)",
      x = "Cell Type"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "#4A4A4A"),
      axis.title.x = element_text(face = "plain", size = 14, color = "#4A4A4A"),
      axis.title.y = element_text(face = "plain", size = 14, color = "#4A4A4A"),
      axis.text.x = element_text(color = "#4A4A4A", angle = 45, hjust = 1),
      axis.text.y = element_text(color = "#4A4A4A"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "none"
    )

  print(p)

  # Save the plot as SVG with group in file name
  filename <- file.path(results_dir, paste0("TKNK_zscore_plot_group_", grp, ".svg"))
  ggsave(filename = filename, plot = p, width = 3, height = 5, dpi = 300)
}

# Combined plot: All groups in one figure, faceted by Group
p_all <- ggplot(df_z, aes(x = Celltype, y = TKNK_zscore, fill = Celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_boxplot(outlier.shape = NA, color = "#4D4D4D", size = 0.5, width = 0.7, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 5, alpha = 0.7, shape = 21, stroke = 0.4, aes(fill = Celltype), color = "#DFDFDF") +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = c("#A8DADC", "#457B9D", "#E63946", "#FFB703", "#8ecae6")) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  labs(
    title = "TKNK_MOUSE Expression (z-score) Across Cell Types and Groups",
    y = "TKNK_MOUSE (z-score)",
    x = "Cell Type"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "#4A4A4A"),
    axis.title.x = element_text(face = "plain", size = 14, color = "#4A4A4A"),
    axis.title.y = element_text(face = "plain", size = 14, color = "#4A4A4A"),
    axis.text.x = element_text(color = "#4A4A4A", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "#4A4A4A"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

print(p_all)

# Save the combined plot
filename_all <- file.path(results_dir, "TKNK_zscore_plot_all_groups.svg")
ggsave(filename = filename_all, plot = p_all, width = 8, height = 6, dpi = 300)

# Aggregate: mean z-score per Celltype and Group
df_heat <- df_z %>%
  group_by(Celltype, Group) %>%
  summarise(mean_z = mean(TKNK_zscore, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_z)

# Convert to matrix for pheatmap (Celltypes as rownames)
mat <- df_heat %>%
  column_to_rownames("Celltype") %>%
  as.matrix()

heat_colors <- colorRampPalette(c("#457B9D", "white", "#E63946"))(100)

# Determine the range of the matrix to center 0
max_abs <- max(abs(mat), na.rm = TRUE)
breaks_seq <- seq(-max_abs, max_abs, length.out = 101)

# Generate the heatmap with custom breaks
pheatmap(
  mat,
  color = heat_colors,
  breaks = breaks_seq,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 12,
  angle_col = 45,
  main = "TKNK_MOUSE Z-Score Heatmap",
  filename = file.path(results_dir, "TKNK_zscore_heatmap.pdf"),
  width = 6,
  height = 5
)