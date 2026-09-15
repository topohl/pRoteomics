#!/usr/bin/env Rscript
# ================================================================
# Script: figures/manuscript_supporting_immunostaining_candidates.R
# Stage: manuscript_candidates
# Scope: global
# Consumes: required results/source_data/10_biological_integration/immunostaining_candidate_comparison/*.csv
# Produces: results/figures/manuscript/immunostaining_candidates/
# Notes: Rendering only. Reads frozen source data; computes no statistics.
# ================================================================
#
# Manuscript-supporting candidate comparison: three nominated immunostaining
# targets across the ten neuropil region x layer units.
#
# This renderer is deliberately thin. It reads only the frozen source-data CSVs
# written by 10_biological_integration/10_prepare_immunostaining_candidate_comparison.R
# and never reaches back to the GCT, the mapped DA files or any model. It does
# no aggregation: the abundance table is already bilateral animal-level, one
# point per animal, and the log2FC and FDR values are canonical copies.
#
# This is a candidate-comparison / validation-planning figure. It is NOT part of
# canonical Figure 2 or Figure 3 and makes no new inferential claim.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

SRC <- path_results("source_data", "10_biological_integration",
                    "immunostaining_candidate_comparison")
OUT <- path_results("figures", "manuscript", "immunostaining_candidates")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

need <- function(f) {
  p <- file.path(SRC, f)
  if (!file.exists(p))
    stop("missing_required_input: ", p,
         "\nRun 10_biological_integration/10_prepare_immunostaining_candidate_comparison.R first.",
         call. = FALSE)
  utils::read.csv(p, stringsAsFactors = FALSE)
}
abundance <- need("immunostaining_candidates_animal_abundance.csv")
effects   <- need("immunostaining_candidates_sus_res_effects.csv")

PAL <- yaml::read_yaml(repo_path("config", "manuscript_palette.yml"))
GROUP_COL <- unlist(PAL$group)[c("CON", "RES", "SUS")]
FAM <- "sans"
BASE_PT <- 6

# One anatomical order, taken from the source data's spatial_order column and
# applied identically to every protein in both panels.
lev <- abundance %>% distinct(spatial_unit_display, spatial_order) %>%
  arrange(spatial_order) %>% pull(spatial_unit_display)
prot_lev <- c("Oga", "Slc22a23", "Anxa2")
prot_lab <- c(Oga = "OGA (Q9EQQ9)", Slc22a23 = "SLC22A23 (Q3UHH2)",
              Anxa2 = "ANXA2 (P07356)")

abundance <- abundance %>%
  mutate(unit = factor(spatial_unit_display, levels = lev),
         ExpGroup = factor(ExpGroup, levels = c("CON", "RES", "SUS")),
         protein = factor(GeneSymbol, levels = prot_lev))
effects <- effects %>%
  mutate(unit = factor(spatial_unit_display, levels = lev),
         protein = factor(GeneSymbol, levels = prot_lev))
stopifnot(identical(levels(abundance$unit), levels(effects$unit)))

group_mean <- abundance %>%
  group_by(protein, unit, ExpGroup) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop")

theme_ms <- function() {
  theme_classic(base_size = BASE_PT, base_family = FAM) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = BASE_PT + 1, face = "bold", hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_PT - 1),
      axis.text.y = element_text(size = BASE_PT - 1),
      axis.title = element_text(size = BASE_PT),
      legend.position = "top",
      legend.key.size = unit(3, "mm"),
      legend.title = element_text(size = BASE_PT),
      legend.text = element_text(size = BASE_PT),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2),
      plot.title = element_text(size = BASE_PT + 1, face = "bold"),
      plot.margin = margin(2, 3, 2, 3, "mm"))
}

# ------------------------------------------------------------------ panel A
# One point per animal. No bars, no hemisphere points, no error bars: with three
# animals per group a dispersion bar would imply more than the design supports.
# The group mean is drawn FIRST so the animal points stay on top of it; drawn
# the other way round the wide mean tick hides the very points it summarises.
pA <- ggplot(abundance, aes(unit, abundance, colour = ExpGroup)) +
  geom_point(data = group_mean,
             aes(unit, mean_abundance, colour = ExpGroup),
             position = position_dodge(width = 0.66), shape = 95, size = 4.2,
             alpha = 0.55, show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.66), size = 0.85,
             alpha = 0.95, stroke = 0) +
  facet_wrap(~ protein, nrow = 1, scales = "free_y",
             labeller = labeller(protein = prot_lab)) +
  scale_colour_manual(values = GROUP_COL, name = NULL) +
  labs(x = NULL, y = "Animal-level abundance (log2)") +
  theme_ms()

# ------------------------------------------------------------------ panel B
# Canonical SUS-RES log2FC. No confidence interval is drawn because the
# canonical differential-abundance output does not carry one, and deriving one
# here would be a new statistic.
eff_sup <- effects %>% filter(fdr_supported)
pB <- ggplot(effects, aes(unit, log2FC)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey45") +
  geom_point(size = 1.1, colour = "grey30", stroke = 0) +
  geom_point(data = eff_sup, size = 1.9, shape = 21, fill = NA,
             colour = "#1F3D52", stroke = 0.5) +
  # A supported unit sitting at the right-hand edge would push its label off the
  # panel, so the label flips to the inside for units in the last third.
  geom_text(data = eff_sup %>%
              mutate(hj = ifelse(as.integer(unit) > length(lev) * 2 / 3, 1,
                                 ifelse(as.integer(unit) <= length(lev) / 3, 0, 0.5)),
                     nudge = ifelse(hj == 1, -0.25, ifelse(hj == 0, 0.25, 0))),
            aes(x = as.integer(unit) + nudge, y = log2FC, hjust = hj,
                label = sprintf("FDR %.3g", fdr_bh)),
            vjust = -1.3, size = (BASE_PT - 2) / .pt, colour = "#1F3D52",
            inherit.aes = FALSE) +
  facet_wrap(~ protein, nrow = 1, scales = "free_y",
             labeller = labeller(protein = prot_lab)) +
  scale_y_continuous(expand = expansion(mult = c(0.16, 0.22))) +
  scale_x_discrete(expand = expansion(add = 0.7)) +
  labs(x = NULL, y = expression("SUS - RES log"[2]*" fold change")) +
  theme_ms() +
  theme(legend.position = "none")

# ------------------------------------------------------------------- output
# svglite is the repository's vector convention; text stays live and nothing is
# rasterised.
save_trio <- function(plot, stem, w_mm, h_mm) {
  svglite::svglite(file.path(OUT, paste0(stem, ".svg")),
                   width = w_mm / 25.4, height = h_mm / 25.4, bg = "white",
                   fix_text_size = FALSE)
  print(plot); grDevices::dev.off()
  grDevices::cairo_pdf(file.path(OUT, paste0(stem, ".pdf")),
                       width = w_mm / 25.4, height = h_mm / 25.4, bg = "white")
  print(plot); grDevices::dev.off()
  grDevices::png(file.path(OUT, paste0(stem, ".png")),
                 width = w_mm, height = h_mm, units = "mm", res = 600,
                 type = "cairo", bg = "white")
  print(plot); grDevices::dev.off()
  invisible(file.path(OUT, paste0(stem, ".svg")))
}

W <- 183; HA <- 62; HB <- 52

# composable panel-level assets, matching how the figure layers keep panels
svglite::svglite(file.path(OUT, "immunostaining_candidates_abundance.svg"),
                 width = W / 25.4, height = HA / 25.4, bg = "white",
                 fix_text_size = FALSE)
print(pA); grDevices::dev.off()
svglite::svglite(file.path(OUT, "immunostaining_candidates_sus_res_log2fc.svg"),
                 width = W / 25.4, height = HB / 25.4, bg = "white",
                 fix_text_size = FALSE)
print(pB); grDevices::dev.off()

combined <- patchwork::wrap_plots(pA, pB, ncol = 1,
                                  heights = c(HA, HB)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(plot.margin = margin(1, 1, 1, 1, "mm")))
save_trio(combined, "immunostaining_candidates_spatial_comparison", W, HA + HB)

svg_main <- file.path(OUT, "immunostaining_candidates_spatial_comparison.svg")
sz <- file.size(svg_main)
cat("immunostaining candidate figure ->", OUT, "\n")
cat("  proteins:", paste(prot_lev, collapse = ", "), "\n")
cat("  spatial units:", length(lev), "in canonical order:",
    paste(lev, collapse = " > "), "\n")
cat("  animals plotted:", length(unique(abundance$AnimalID)),
    "| points in panel A:", nrow(abundance), "\n")
cat("  FDR-supported points marked:", nrow(eff_sup), "\n")
cat("  main SVG:", basename(svg_main), "|", sz, "bytes | non-empty:", sz > 0, "\n")
