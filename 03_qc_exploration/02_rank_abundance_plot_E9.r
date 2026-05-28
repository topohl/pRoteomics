# ==============================================================================
# E9 Social instability stress - Rank abundance QC plots
# Unified Soma-style visual language for soma, neuropil, microglia
# ==============================================================================

library(readxl)
library(tidyverse)
library(scales)
library(ggrepel)
library(writexl)
library(svglite)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

# ==============================================================================
# 0. Paths
# ==============================================================================

saving_dir <- Sys.getenv("PROTEOMICS_RANK_ABUNDANCE_DIR", unset = path_results("QC"))
ensure_dir(saving_dir)

files <- list(
  soma = Sys.getenv("PROTEOMICS_RANK_ABUNDANCE_SOMA", unset = path_processed("pg_matrix", "imputed", "20260218_pgmatrix_imputed_neuron_soma_71samples_missing70pct.xlsx")),
  neuropil = Sys.getenv("PROTEOMICS_RANK_ABUNDANCE_NEUROPIL", unset = path_processed("pg_matrix", "imputed", "20260218_pgmatrix_imputed_neuron_neuropil_180samples_missing70pct.xlsx")),
  microglia = Sys.getenv("PROTEOMICS_RANK_ABUNDANCE_MICROGLIA", unset = path_processed("pg_matrix", "imputed", "20260218_pgmatrix_imputed_microglia_72samples_missing70pct.xlsx"))
)
missing_files <- files[!file.exists(unlist(files))]
if (length(missing_files) > 0) {
  stop("Rank abundance input file(s) not found: ", paste(unlist(missing_files), collapse = ", "), call. = FALSE)
}

# ==============================================================================
# 1. Palettes and markers
# ==============================================================================

region_colors <- c(
  "CA1" = "#2C3E50",
  "CA2" = "#16A085",
  "CA3" = "#E67E22",
  "DG"  = "#8E44AD"
)

marker_colors <- c(
  region_colors,
  "Specific" = "#E67E22",
  "Neuropil" = "#2980B9",
  "Microglia" = "#C0392B",
  "Contamination" = "grey40"
)

soma_markers <- c(
  "H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3",
  "Matr3", "Srsf3", "Ddx39b"
)

neuropil_markers <- c(
  "Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2", "Snap25"
)

microglia_markers <- c(
  "Aif1", "Tmem119", "P2ry12", "Cx3cr1",
  "Csf1r", "C1qa", "Hexb", "Mertk"
)

neuronal_contamination_markers <- c(
  "Syn1", "Syp", "Atp5f1a", "Tuba1a", "Nefl"
)

# ==============================================================================
# 2. Helper functions
# ==============================================================================

extract_region <- function(sample) {
  case_when(
    str_detect(sample, fixed("CA1", ignore_case = TRUE)) ~ "CA1",
    str_detect(sample, fixed("CA2", ignore_case = TRUE)) ~ "CA2",
    str_detect(sample, fixed("CA3", ignore_case = TRUE)) ~ "CA3",
    str_detect(sample, fixed("DG",  ignore_case = TRUE)) ~ "DG",
    TRUE ~ NA_character_
  )
}

extract_layer <- function(sample) {
  case_when(
    str_detect(sample, fixed("SLM", ignore_case = TRUE)) ~ "SLM",
    str_detect(sample, fixed("SO",  ignore_case = TRUE)) ~ "SO",
    str_detect(sample, fixed("SR",  ignore_case = TRUE)) ~ "SR",
    str_detect(sample, fixed("MO",  ignore_case = TRUE)) ~ "MO",
    str_detect(sample, fixed("PO",  ignore_case = TRUE)) ~ "PO",
    TRUE ~ NA_character_
  )
}

make_rank_data <- function(file_path, mode = c("region", "region_layer")) {
  mode <- match.arg(mode)

  read_excel(file_path) %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "Sample",
      values_to = "Log2Intensity"
    ) %>%
    mutate(
      Region = extract_region(Sample),
      Layer  = extract_layer(Sample)
    ) %>%
    {
      if (mode == "region") {
        filter(., !is.na(Region)) %>%
          group_by(Region, Genes) %>%
          summarise(MeanLog2 = mean(Log2Intensity, na.rm = TRUE), .groups = "drop") %>%
          mutate(Group = Region)
      } else {
        filter(., !is.na(Region), !is.na(Layer)) %>%
          group_by(Region, Layer, Genes) %>%
          summarise(MeanLog2 = mean(Log2Intensity, na.rm = TRUE), .groups = "drop") %>%
          mutate(Group = paste(Region, Layer, sep = " - "))
      }
    } %>%
    mutate(LinearValue = 2^MeanLog2) %>%
    group_by(Group) %>%
    arrange(desc(LinearValue), .by_group = TRUE) %>%
    mutate(Rank = row_number()) %>%
    ungroup()
}

save_rank_excel <- function(plot_data, file_name) {
  processed <- plot_data %>%
    select(any_of(c(
      "Region", "Layer", "Group", "Rank", "Genes",
      "MeanLog2", "LinearValue", "MarkerType"
    ))) %>%
    arrange(Group, Rank)

  excel_list <- split(processed, processed$Group)

  sheet_names <- names(excel_list) %>%
    str_replace_all("[:\\\\/?*\\[\\]]", "-") %>%
    str_sub(1, 31)

  names(excel_list) <- make.unique(sheet_names, sep = "_")

  write_xlsx(
    excel_list,
    path = file.path(saving_dir, paste0(file_name, "_processed_ranks.xlsx"))
  )
}

plot_rank_abundance <- function(
    plot_data,
    file_name,
    ncol = 2,
    width_mm = 120,
    height_mm = 130,
    show_legend = FALSE
) {

  label_df <- plot_data %>%
    filter(MarkerType != "None")

  rank_plot <- ggplot(plot_data, aes(x = Rank, y = LinearValue)) +

    geom_point(
      alpha = 0.08,
      size = 0.15,
      color = "grey82"
    ) +

    geom_line(
      aes(color = Region),
      alpha = 0.85,
      linewidth = 0.35
    ) +

    geom_point(
      data = label_df,
      aes(color = MarkerColor),
      size = 1.45,
      alpha = 1
    ) +

    geom_label_repel(
      data = label_df,
      aes(label = Genes, fill = MarkerColor),
      color = "white",
      size = 2.35,
      fontface = "bold",
      family = "sans",
      label.padding = unit(0.15, "lines"),
      label.r = unit(0, "lines"),
      label.size = 0,
      segment.color = "grey20",
      segment.linewidth = 0.25,
      segment.alpha = 0.85,
      box.padding = 0.35,
      point.padding = 0.15,
      force = 10,
      max.overlaps = Inf
    ) +

    scale_y_log10(
      expand = expansion(mult = c(0.05, 0.22)),
      labels = trans_format("log10", math_format(10^.x))
    ) +

    scale_x_continuous(
      expand = expansion(mult = c(0.025, 0.02)),
      labels = label_comma()
    ) +

    scale_color_manual(values = marker_colors) +
    scale_fill_manual(values = marker_colors) +

    facet_wrap(~Group, scales = "free_x", ncol = ncol) +

    labs(
      x = "Protein rank",
      y = "Intensity (log10)"
    ) +

    theme_minimal(base_size = 8) +
    theme(
      aspect.ratio = 1,
      text = element_text(family = "sans", color = "black"),
      strip.text = element_text(
        face = "bold",
        size = 9,
        margin = margin(b = 8)
      ),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 6, color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      legend.position = ifelse(show_legend, "bottom", "none"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  svg_path <- file.path(saving_dir, paste0(file_name, ".svg"))

  ggsave(
    svg_path,
    rank_plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = svglite
  )

  save_rank_excel(plot_data, file_name)

  message("Saved plot: ", svg_path)

  return(rank_plot)
}

# ==============================================================================
# 3. Soma
# ==============================================================================

soma_data <- make_rank_data(files$soma, mode = "region") %>%
  mutate(
    MarkerType = case_when(
      Genes %in% soma_markers ~ "Specific",
      TRUE ~ "None"
    ),
    MarkerColor = case_when(
      MarkerType == "Specific" ~ Region,
      TRUE ~ "None"
    )
  )

plot_soma <- plot_rank_abundance(
  soma_data,
  file_name = "rank_abundance_E9_soma_somaStyle",
  ncol = 2,
  width_mm = 120,
  height_mm = 135,
  show_legend = FALSE
)

# ==============================================================================
# 4. Neuropil
# ==============================================================================

neuropil_data <- make_rank_data(files$neuropil, mode = "region_layer") %>%
  mutate(
    MarkerType = case_when(
      Genes %in% soma_markers     ~ "Specific",
      Genes %in% neuropil_markers ~ "Neuropil",
      TRUE ~ "None"
    ),
    MarkerColor = case_when(
      MarkerType == "Specific" ~ "Specific",
      MarkerType == "Neuropil" ~ "Neuropil",
      TRUE ~ "None"
    )
  )

plot_neuropil <- plot_rank_abundance(
  neuropil_data,
  file_name = "rank_abundance_E9_neuropil_somaStyle",
  ncol = 4,
  width_mm = 180,
  height_mm = 210,
  show_legend = FALSE
)

# ==============================================================================
# 5. Microglia
# ==============================================================================

microglia_data <- make_rank_data(files$microglia, mode = "region") %>%
  mutate(
    MarkerType = case_when(
      Genes %in% microglia_markers              ~ "Microglia",
      Genes %in% neuronal_contamination_markers ~ "Contamination",
      TRUE ~ "None"
    ),
    MarkerColor = case_when(
      MarkerType == "Microglia"      ~ "Microglia",
      MarkerType == "Contamination"  ~ "Contamination",
      TRUE ~ "None"
    )
  )

plot_microglia <- plot_rank_abundance(
  microglia_data,
  file_name = "rank_abundance_E9_microglia_somaStyle",
  ncol = 2,
  width_mm = 120,
  height_mm = 135,
  show_legend = FALSE
)

message("All rank abundance plots completed.")

# ==============================================================================
# 6. Cross-compartment marker comparison:
#    Microglia markers in microglia vs neuropil
#    Neuropil markers in neuropil vs microglia
# ==============================================================================

cross_marker_data <- bind_rows(
  microglia_data %>%
    mutate(Compartment = "Microglia"),

  neuropil_data %>%
    mutate(Compartment = "Neuropil")
) %>%
  mutate(
    CrossMarkerType = case_when(
      Genes %in% microglia_markers ~ "Microglia markers",
      Genes %in% neuropil_markers  ~ "Neuropil markers",
      Genes %in% neuronal_contamination_markers ~ "Neuronal / synaptic markers",
      Genes %in% soma_markers ~ "Nuclear / soma markers",
      TRUE ~ "None"
    ),
    CrossMarkerColor = case_when(
      CrossMarkerType == "Microglia markers" ~ "Microglia",
      CrossMarkerType == "Neuropil markers" ~ "Neuropil",
      CrossMarkerType == "Neuronal / synaptic markers" ~ "Contamination",
      CrossMarkerType == "Nuclear / soma markers" ~ "Specific",
      TRUE ~ "None"
    ),
    PlotGroup = paste(Compartment, Group, sep = " | ")
  )

cross_marker_plot_data <- cross_marker_data %>%
  filter(CrossMarkerType != "None")

# ==============================================================================
# 6A. Rank-abundance plot: same markers shown in both compartments
# ==============================================================================

cross_rank_plot <- ggplot(cross_marker_data, aes(x = Rank, y = LinearValue)) +

  geom_point(
    alpha = 0.06,
    size = 0.12,
    color = "grey82"
  ) +

  geom_line(
    aes(color = Region),
    alpha = 0.75,
    linewidth = 0.32
  ) +

  geom_point(
    data = cross_marker_plot_data,
    aes(color = CrossMarkerColor),
    size = 1.35,
    alpha = 1
  ) +

  geom_label_repel(
    data = cross_marker_plot_data,
    aes(label = Genes, fill = CrossMarkerColor),
    color = "white",
    size = 2.15,
    fontface = "bold",
    family = "sans",
    label.padding = unit(0.13, "lines"),
    label.r = unit(0, "lines"),
    label.size = 0,
    segment.color = "grey20",
    segment.linewidth = 0.22,
    segment.alpha = 0.8,
    box.padding = 0.35,
    point.padding = 0.12,
    force = 8,
    max.overlaps = Inf
  ) +

  scale_y_log10(
    expand = expansion(mult = c(0.05, 0.25)),
    labels = trans_format("log10", math_format(10^.x))
  ) +

  scale_x_continuous(
    expand = expansion(mult = c(0.025, 0.02)),
    labels = label_comma()
  ) +

  scale_color_manual(values = marker_colors) +
  scale_fill_manual(values = marker_colors) +

  facet_wrap(~PlotGroup, scales = "free_x", ncol = 4) +

  labs(
    x = "Protein rank",
    y = "Intensity (log10)"
  ) +

  theme_minimal(base_size = 8) +
  theme(
    aspect.ratio = 1,
    text = element_text(family = "sans", color = "black"),
    strip.text = element_text(
      face = "bold",
      size = 6.5,
      margin = margin(b = 6)
    ),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 5.5, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.7, "lines"),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(saving_dir, "rank_abundance_crossCompartment_markerCheck.svg"),
  cross_rank_plot,
  width = 240,
  height = 260,
  units = "mm",
  device = svglite
)

# ==============================================================================
# 6B. Cleaner summary plot: marker abundance by compartment
# ==============================================================================

marker_summary <- cross_marker_data %>%
  filter(CrossMarkerType != "None") %>%
  group_by(Compartment, Region, Genes, CrossMarkerType) %>%
  summarise(
    MeanLog2 = mean(MeanLog2, na.rm = TRUE),
    MedianRank = median(Rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    CrossMarkerType = factor(
      CrossMarkerType,
      levels = c(
        "Microglia markers",
        "Neuropil markers",
        "Neuronal / synaptic markers",
        "Nuclear / soma markers"
      )
    )
  )

marker_summary_plot <- ggplot(
  marker_summary,
  aes(x = Compartment, y = MeanLog2)
) +

  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    linewidth = 0.3,
    fill = "grey92",
    color = "black"
  ) +

  geom_point(
    aes(color = CrossMarkerType),
    position = position_jitter(width = 0.12, height = 0),
    size = 1.6,
    alpha = 0.9
  ) +

  facet_wrap(~CrossMarkerType, scales = "free_y", nrow = 1) +

  scale_color_manual(
    values = c(
      "Microglia markers" = marker_colors["Microglia"],
      "Neuropil markers" = marker_colors["Neuropil"],
      "Neuronal / synaptic markers" = marker_colors["Contamination"],
      "Nuclear / soma markers" = marker_colors["Specific"]
    )
  ) +

  labs(
    x = NULL,
    y = "Mean log2 intensity"
  ) +

  theme_classic(base_size = 8) +
  theme(
    text = element_text(family = "sans", color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 6, color = "black"),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  )

ggsave(
  file.path(saving_dir, "crossCompartment_marker_abundance_summary.svg"),
  marker_summary_plot,
  width = 160,
  height = 60,
  units = "mm",
  device = svglite
)

# ==============================================================================
# 6C. Export table
# ==============================================================================

cross_marker_table <- marker_summary %>%
  arrange(CrossMarkerType, Genes, Compartment, Region)

write_xlsx(
  list(
    marker_summary = cross_marker_table,
    full_cross_marker_data = cross_marker_data %>%
      filter(CrossMarkerType != "None") %>%
      arrange(CrossMarkerType, Genes, Compartment, Group)
  ),
  path = file.path(saving_dir, "crossCompartment_marker_check.xlsx")
)

message("Cross-compartment marker check completed.")


# ==============================================================================
# 7. CROSS-COMPARTMENT MARKER MODULE SCORES
#    Soma vs Neuropil vs Microglia, with statistics
# ==============================================================================

library(ggpubr)
library(rstatix)

# ------------------------------------------------------------------------------
# 7.1 Marker sets
# ------------------------------------------------------------------------------

marker_sets <- list(
  Soma_Nuclear = soma_markers,
  Neuropil_Synaptic = unique(c(neuropil_markers, neuronal_contamination_markers)),
  Microglia = microglia_markers
)

compartment_colors <- c(
  "Soma"     = "#2C3E50",
  "Neuropil" = "#2980B9",
  "Microglia" = "#C0392B"
)

# ------------------------------------------------------------------------------
# 7.2 Read all raw datasets
# ------------------------------------------------------------------------------

raw_cross_all <- bind_rows(
  read_excel(files$soma) %>%
    mutate(Compartment = "Soma"),

  read_excel(files$neuropil) %>%
    mutate(Compartment = "Neuropil"),

  read_excel(files$microglia) %>%
    mutate(Compartment = "Microglia")
)

long_cross_all <- raw_cross_all %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "Sample",
    values_to = "Log2Intensity"
  ) %>%
  filter(!is.na(Log2Intensity)) %>%
  mutate(
    Region = extract_region(Sample),
    Layer  = extract_layer(Sample),
    Layer  = case_when(
      Compartment == "Soma" & is.na(Layer) ~ "Soma",
      Compartment == "Microglia" & is.na(Layer) ~ "Microglia",
      TRUE ~ Layer
    ),
    Group = paste(Region, Layer, sep = " - ")
  ) %>%
  filter(!is.na(Region))

# ------------------------------------------------------------------------------
# 7.3 Calculate marker module scores per sample
# ------------------------------------------------------------------------------

score_all <- imap_dfr(marker_sets, function(markers, set_name) {

  long_cross_all %>%
    filter(Genes %in% markers) %>%
    group_by(Compartment, Region, Layer, Group, Sample) %>%
    summarise(
      Score = mean(Log2Intensity, na.rm = TRUE),
      DetectedMarkers = sum(!is.na(Log2Intensity)),
      .groups = "drop"
    ) %>%
    mutate(
      MarkerSet = set_name,
      MarkerN = length(markers)
    )
}) %>%
  mutate(
    Compartment = factor(
      Compartment,
      levels = c("Soma", "Neuropil", "Microglia")
    ),
    MarkerSet = factor(
      MarkerSet,
      levels = c("Soma_Nuclear", "Neuropil_Synaptic", "Microglia"),
      labels = c("Soma / nuclear markers", "Neuropil / synaptic markers", "Microglia markers")
    )
  )

# ------------------------------------------------------------------------------
# 7.4 Pairwise statistics: compartments within each marker set
# ------------------------------------------------------------------------------

pairwise_stats <- score_all %>%
  group_by(MarkerSet) %>%
  pairwise_wilcox_test(
    Score ~ Compartment,
    p.adjust.method = "holm"
  ) %>%
  add_xy_position(
    x = "Compartment",
    dodge = 0.8,
    step.increase = 0.12
  ) %>%
  mutate(
    p_label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

global_stats <- score_all %>%
  group_by(MarkerSet) %>%
  kruskal_test(Score ~ Compartment)

write_xlsx(
  list(
    module_scores = score_all,
    pairwise_wilcox = pairwise_stats,
    global_kruskal = global_stats
  ),
  file.path(saving_dir, "crossCompartment_allCellTypes_marker_module_scores_stats.xlsx")
)

# ------------------------------------------------------------------------------
# 7.5 Publication-style module score plot with cleaner statistics
# ------------------------------------------------------------------------------

# Cleaner p-value labels
pairwise_stats_plot <- pairwise_stats %>%
  mutate(
    p_label = case_when(
      p.adj < 0.0001 ~ "P < 0.0001",
      p.adj < 0.001  ~ "P < 0.001",
      p.adj < 0.01   ~ paste0("P = ", signif(p.adj, 2)),
      p.adj < 0.05   ~ paste0("P = ", signif(p.adj, 2)),
      TRUE ~ "n.s."
    )
  )

# Optional: only show significant brackets to reduce clutter
pairwise_stats_sig <- pairwise_stats_plot %>%
  filter(p.adj < 0.05)

module_score_plot_all <- ggplot(
  score_all,
  aes(x = Compartment, y = Score)
) +

  # thin reference structure
  geom_boxplot(
    aes(fill = Compartment),
    width = 0.48,
    outlier.shape = NA,
    linewidth = 0.32,
    color = "black",
    alpha = 0.92
  ) +

  # individual samples
  geom_point(
    aes(color = Compartment),
    position = position_jitter(width = 0.09, height = 0),
    size = 1.05,
    alpha = 0.55,
    stroke = 0
  ) +

  # median emphasis
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.42,
    linewidth = 0.28,
    color = "black"
  ) +

  stat_pvalue_manual(
    pairwise_stats_sig,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.006,
    bracket.size = 0.24,
    size = 2.1,
    hide.ns = TRUE
  ) +

  facet_wrap(
    ~MarkerSet,
    scales = "free_y",
    nrow = 1
  ) +

  scale_fill_manual(values = compartment_colors) +
  scale_color_manual(values = compartment_colors) +

  labs(
    x = NULL,
    y = "Marker module score\n(mean log2 intensity)"
  ) +

  theme_classic(base_size = 7) +
  theme(
    text = element_text(family = "Arial", color = "black"),

    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 7.2,
      margin = margin(b = 6)
    ),

    axis.title.y = element_text(
      size = 7.2,
      face = "bold",
      margin = margin(r = 5)
    ),

    axis.text.x = element_text(
      size = 6.2,
      color = "black",
      angle = 35,
      hjust = 1,
      vjust = 1
    ),

    axis.text.y = element_text(
      size = 6.2,
      color = "black"
    ),

    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    axis.ticks.length = unit(1.5, "mm"),

    panel.spacing = unit(1.1, "lines"),
    legend.position = "none",

    plot.margin = margin(4, 4, 4, 4),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(saving_dir, "publication_panel_allCellTypes_marker_module_scores_withStats_clean.svg"),
  module_score_plot_all,
  width = 155,
  height = 65,
  units = "mm",
  device = svglite
)

module_score_plot_all
