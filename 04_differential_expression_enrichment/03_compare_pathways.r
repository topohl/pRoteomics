library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(svglite)
library(gridExtra)
library(pheatmap)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

# Define file paths
path_memory_ensemble <- Sys.getenv(
    "PROTEOMICS_COMPARE_PATHWAYS_MEMORY",
    unset = path_results("tables", "04_differential_expression_enrichment", "compareGO", "BP", "learning_signature", "memory_ensemble", "01_Tables_and_Data", "Supplementary_Data.xlsx")
)
path_effects_inhibition <- Sys.getenv(
    "PROTEOMICS_COMPARE_PATHWAYS_EFFECTS",
    unset = path_results("tables", "04_differential_expression_enrichment", "compareGO", "BP", "learning_signature", "effects_inhibition_memory_ensemble", "01_Tables_and_Data", "Supplementary_Data.xlsx")
)
if (!file.exists(path_memory_ensemble)) stop("Memory ensemble pathway workbook not found: ", path_memory_ensemble, call. = FALSE)
if (!file.exists(path_effects_inhibition)) stop("Effects inhibition pathway workbook not found: ", path_effects_inhibition, call. = FALSE)

# Define directories for saving output
output_dir <- Sys.getenv("PROTEOMICS_COMPARE_PATHWAYS_OUTPUT_DIR", unset = path_results("compare_sig_expr"))
tables_dir <- file.path(output_dir, "tables")
plots_dir <- file.path(output_dir, "plots")
ensure_dir(tables_dir)
ensure_dir(plots_dir)

# Read the Excel files into data frames
data_memory_ensemble <- read_excel(path_memory_ensemble, sheet = "GO_Enrichment_Results")
data_effects_inhibition <- read_excel(path_effects_inhibition, sheet = "GO_Enrichment_Results")

# Inspect the first few rows of the datasets
head(data_memory_ensemble)
head(data_effects_inhibition)
library(dplyr)
library(tidyr)
library(ggplot2)

# Add dataset labels
data_memory_ensemble <- data_memory_ensemble %>% mutate(Dataset = "Memory_Ensemble")
data_effects_inhibition <- data_effects_inhibition %>% mutate(Dataset = "Effects_Inhibition")

# Combine the relevant columns from both datasets and filter by p.adjust
combined_data <- bind_rows(data_memory_ensemble, data_effects_inhibition) %>%
    filter(p.adjust <= 0.05) %>%
    dplyr::select(Dataset, Comparison, Description, NES, p.adjust)

# Create an aligned comparison group to match specific corresponding comparisons
combined_data <- combined_data %>%
    mutate(Comparison_Group = case_when(
        Comparison %in% c("bg2_bg4", "bg1_bg2") ~ "bg",
        Comparison %in% c("mcherry2_mcherry4", "mcherry1_mcherry2") ~ "mcherry",
        Comparison %in% c("neuron2_neuron4", "neuron1_neuron2") ~ "neuron",
        Comparison %in% c("cfos2_cfos4", "cfos1_cfos2") ~ "cfos",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Comparison_Group))

# Identify overlapping pathways present in BOTH datasets for matching comparison groups
overlapping_pathways_df <- combined_data %>%
    group_by(Description, Comparison_Group) %>%
    filter(n_distinct(Dataset) == 2) %>%
    ungroup()

# Extract the descriptions and comparisons, adding p.adjust, NES, and Comparison by dataset
overlapping_pathways <- overlapping_pathways_df %>%
    distinct(Dataset, Description, Comparison_Group, .keep_all = TRUE) %>%
    pivot_wider(
        id_cols = c(Description, Comparison_Group),
        names_from = Dataset,
        values_from = c(NES, p.adjust, Comparison)
    )

View(overlapping_pathways)

# Calculate the difference in NES between the two datasets
heatmap_data <- overlapping_pathways %>%
    mutate(NES_Difference = NES_Effects_Inhibition - NES_Memory_Ensemble)

# Define color palette
comparison_colors <- c("bg" = "#4c87c6", "cfos" = "#6ccff6", "mcherry" = "#faa51a", "neuron" = "#fdd700")

# Plot the improved heatmap
ggplot(heatmap_data, aes(x = Comparison_Group, y = Description, fill = NES_Difference)) +
    geom_tile(color = "white", linewidth = 0.4, width = 0.6, height = 0.6) +
    scale_fill_gradient2(
        low = "#2166AC", 
        mid = "#F7F7F7", 
        high = "#B2182B", 
        midpoint = 0, 
        name = "\u0394 NES\n(Inhibition - Memory)",
        na.value = "grey90",
        limits = c(-max(abs(heatmap_data$NES_Difference), na.rm = TRUE), 
                   max(abs(heatmap_data$NES_Difference), na.rm = TRUE))
    ) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 11, color = "#2b2b2b"),
        axis.text.y = element_text(size = 10, color = "#2b2b2b"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 8), color = "#1a1a1a"),
        plot.subtitle = element_text(size = 11, color = "#555555", hjust = 0.5, margin = margin(b = 15)),
        plot.caption = element_text(size = 9, color = "#777777", hjust = 1, margin = margin(t = 15)),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
        title = "Differences in Pathway Regulation (NES)",
        subtitle = "Red: Higher in Effects Inhibition  |  Blue: Higher in Memory Ensemble",
        caption = "Only pathways significantly enriched (p.adj ≤ 0.05) in both datasets are shown."
    ) +
    coord_fixed(ratio = 1)


# Create heatmap matrix based on the overlapping_pathways data frame
library(pheatmap)
library(tibble)

# Create unique row names by combining Description and Comparison_Group, and sort
heatmap_mat_data <- overlapping_pathways %>%
    arrange(Comparison_Group, Description) %>%
    mutate(RowName = paste(Description, Comparison_Group, sep = " | ")) %>%
    dplyr::select(RowName, NES_Memory_Ensemble, NES_Effects_Inhibition) %>%
    column_to_rownames("RowName")

# Extract Comparison_Group for row annotations
annotation_row <- overlapping_pathways %>%
    arrange(Comparison_Group, Description) %>%
    mutate(RowName = paste(Description, Comparison_Group, sep = " | ")) %>%
    dplyr::select(RowName, Comparison_Group) %>%
    column_to_rownames("RowName")

# Convert to matrix
nes_matrix <- as.matrix(heatmap_mat_data)

# Extract descriptions in the exact sorted order for labels
sorted_descriptions <- overlapping_pathways %>%
    arrange(Comparison_Group, Description) %>%
    pull(Description)

# Define annotation colors
annotation_colors <- list(Comparison_Group = comparison_colors)

# Create heatmap using pheatmap
pheatmap(
    nes_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    labels_row = sorted_descriptions,
    color = colorRampPalette(c("#2d8be9", "#F7F7F7", "#d12f42"))(50),
    border_color = NA,
    main = "Heatmap of NES Values Sorted by Comparison Group",
    fontsize_row = 8, 
    cellwidth = 15,
    cellheight = 15
)

# Save the heatmap to a file 
my_pheatmap <-pheatmap(
    nes_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    labels_row = sorted_descriptions,
    color = colorRampPalette(c("#2d8be9", "#F7F7F7", "#d12f42"))(50),
    border_color = "white",
    main = "Heatmap of NES Values Sorted by Comparison Group",
    fontsize_row = 8,
    cellwidth = 15,
    cellheight = 15,
    filename = "NES_Heatmap_Sorted_By_Comparison_Group.png",
    width = 10,
    height = 12
)

my_plot <- grid.arrange(my_pheatmap, nrow = 1, ncol = 1)
ggsave(filename = "NES_Heatmap_Sorted_By_Comparison_Group.svg", plot = my_plot, width = 10, height = 12, bg = "white")

# State where it's saved
cat("Heatmap saved to:", file.path(getwd(), "NES_Heatmap_Sorted_By_Comparison_Group.png"), "\n")

# --- NEW: ggplot version of the absolute NES Heatmap ---
# Prepare data in long format for ggplot
ggplot_nes_data <- heatmap_mat_data %>%
    rownames_to_column("RowName") %>%
    pivot_longer(
        cols = c(NES_Memory_Ensemble, NES_Effects_Inhibition),
        names_to = "Dataset",
        values_to = "NES"
    ) %>%
    mutate(Dataset = gsub("NES_", "", Dataset) %>% gsub("_", " ", .)) %>%
    # Split RowName back to get grouping
    separate(RowName, into = c("Description", "Comparison_Group"), sep = " \\| ", remove = FALSE) %>%
    # Ensure correct ordering of rows (top-to-bottom usually means reversing levels in ggplot)
    arrange(Comparison_Group, Description) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))

# Create ggplot heatmap grouped by Comparison_Group
nes_absolute_ggplot <- ggplot(ggplot_nes_data, aes(x = Dataset, y = Description, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
        low = "#2d8be9", 
        mid = "#F7F7F7", 
        high = "#d12f42", 
        midpoint = 0, 
        name = "NES"
    ) +
    # Mimic the row annotations visually by faceting
    facet_grid(Comparison_Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
        strip.placement = "outside",
        panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 10)),
        panel.grid = element_blank()
    ) +
    labs(title = "Heatmap of NES Values by Comparison Group (ggplot)")
# -----------------------------------------------------


# Create directories if they don't exist
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Save tables as CSV files
write.csv(combined_data, file.path(tables_dir, "combined_data.csv"), row.names = FALSE)
write.csv(overlapping_pathways, file.path(tables_dir, "overlapping_pathways.csv"), row.names = FALSE)
write.csv(heatmap_data, file.path(tables_dir, "heatmap_data.csv"), row.names = FALSE)

# Save the ggplot (Difference Heatmap)
ggsave(
    filename = file.path(plots_dir, "NES_Difference_Heatmap.svg"), 
    plot = last_plot(), # Depending on execution order, you might want to specify `plot = ...`
    width = 10, 
    height = 10, 
    bg = "white"
)

# Save the absolute NES ggplot heatmap
ggsave(
    filename = file.path(plots_dir, "NES_Absolute_Heatmap_ggplot.svg"), 
    plot = nes_absolute_ggplot,
    width = 8, 
    height = 12, 
    bg = "white"
)

# Save the pheatmap to the new plots directory
my_heatmap <- pheatmap(
    t(nes_matrix),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = annotation_row,
    annotation_colors = annotation_colors,
    labels_col = sorted_descriptions,
    color = colorRampPalette(c("#2d8be9", "#F7F7F7", "#d12f42"))(50),
    border_color = "white",
    main = "Heatmap of NES Values Sorted by Comparison Group",
    fontsize_col = 8,
    angle_col = 45,
    cellwidth = 15,
    cellheight = 15,
    filename = file.path(plots_dir, "NES_Heatmap_Sorted_By_Comparison_Group.png"),
    width = 12,
    height = 10
)

my_plot <- grid.arrange(my_heatmap$gtable, nrow = 1, ncol = 1)
ggsave(filename = file.path(plots_dir, "NES_Heatmap_Sorted_By_Comparison_Group.svg"), plot = my_plot, width = 12, height = 10, bg = "white")

cat("Successfully saved tables and plots to:", output_dir, "\n")
