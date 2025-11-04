#' @title QC Protein Peptide Plot
#' @description
#' Creates a quality control visualization of protein and peptide counts across different cell types.
#' The function generates a faceted bar plot with dual y-axes, where protein counts are shown on the 
#' left axis and peptide counts on the right axis. Each cell type is displayed in a separate panel,
#' with horizontal dashed lines indicating overall averages for proteins and peptides. The plot 
#' excludes background samples and uses custom color schemes to distinguish between cell types and 
#' measurement types.
#' @param data A data frame containing protein and peptide counts.
#' @return A ggplot object.
#' @export

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyr, dplyr, readxl)

# Read the data
data <- read_excel("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/protein_count.xlsx")
saving_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/Results/Plots/"

# Reshape data to long format
data_long <- data %>%
    filter(Celltype != "Background") %>%
    pivot_longer(cols = c(proteins, peptides),
                 names_to = "Type",
                 values_to = "Count") %>%
    mutate(ColorGroup = paste(Celltype, Type, sep = "_")) %>%
    arrange(Celltype, id) %>%
    mutate(id = as.character(id))

# Define custom colors
custom_colors <- c(
  "cFosN_proteins" = "#f4a261", "cFosN_peptides" = "#f6bd60",
  "Background_proteins" = "#f9c74f", "Background_peptides" = "#fce38a",
  "mCherryN_proteins" = "#90be6d", "mCherryN_peptides" = "#b9fbc0",
  "Neuron_proteins" = "#43aa8b", "Neuron_peptides" = "#a8e6cf"
)

# Calculate ranges
protein_range <- range(data_long$Count[data_long$Type == "proteins"], na.rm = TRUE)
peptide_range <- range(data_long$Count[data_long$Type == "peptides"], na.rm = TRUE)

# Define min and max for axes
protein_min <- protein_range[1]
protein_max <- protein_range[2]
peptide_min <- peptide_range[1]
peptide_max <- peptide_range[2]

# Calculate spans
protein_span <- protein_max - protein_min
peptide_span <- peptide_max - peptide_min

# Add 'Count_scaled' to data_long
data_long <- data_long %>%
  mutate(Count_scaled = case_when(
    Type == "proteins" ~ Count,
    Type == "peptides" ~ (Count - peptide_min) / peptide_span * protein_span + protein_min
  ))

# Compute overall averages
overall_avg <- data %>%
    filter(Celltype != "Background") %>%
    pivot_longer(cols = c(proteins, peptides),
                 names_to = "Type",
                 values_to = "Count") %>%
    group_by(Type) %>%
    summarise(AvgCount = mean(Count, na.rm = TRUE)) %>%
    mutate(AvgCount_scaled = case_when(
        Type == "proteins" ~ AvgCount,
        Type == "peptides" ~ (AvgCount - peptide_min) / peptide_span * protein_span + protein_min
    ))

# Plot with correct secondary axis transformation
plt <- ggplot(data_long, aes(x = id, y = Count_scaled, fill = ColorGroup)) +
  geom_rect(aes(
    xmin = as.numeric(id) - 0.5,
    xmax = as.numeric(id) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = ColorGroup
  ), alpha = 0.01) +
  geom_bar(stat = "identity", position = "dodge", width = 1) +
  geom_hline(data = overall_avg, aes(yintercept = AvgCount_scaled, color = Type),
             linetype = "dashed", size = 1) +
  facet_grid(. ~ Celltype, scales = "free", space = "free") +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  scale_color_manual(values = c("proteins" = "black", "peptides" = "red")) +
  scale_y_continuous(
    name = "Proteins",
    expand = c(0, 0),
    limits = c(0, protein_max),
    sec.axis = sec_axis(
      trans = ~ (. / protein_max) * peptide_max, 
      name = "Peptides"
    )
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = NULL, fill = NULL) +
  theme_minimal(base_size = 50) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 55),
    axis.title.y.left = element_text(size = 55),
    axis.title.y.right = element_text(size = 55),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 40),
    plot.margin = margin(10, 10, 10, 10),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(face = "bold", size = 55),
    panel.spacing = unit(0, 'lines')
  )

# Save plot
ggsave(filename = paste0(saving_dir, "protein_peptide_counts.svg"),
       plot = plt,
       width = 100,
       height = 50,
       units = "cm")

# Print plot
print(plt)
