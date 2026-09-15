#!/usr/bin/env Rscript
# ================================================================
# Script: figures/manuscript_supporting_immunostaining_panel.R
# Stage: manuscript_candidates
# Scope: global
# Consumes: required results/source_data/10_biological_integration/immunostaining_candidate_panel/*.csv
# Produces: results/figures/manuscript/immunostaining_candidate_panel/
# Notes: Rendering only. Reads frozen source data; computes no statistics.
# ================================================================
#
# The ten screened immunostaining candidates across the ten neuropil
# region x layer units. Same conventions as the three-candidate comparison:
# one point per animal, group mean overlaid, canonical SUS-RES log2FC below,
# identical anatomical order everywhere.
#
# Candidates carry an evidence class, and the figure shows it rather than
# hiding it. Class 1 is FDR-supported outside CA2-SLM; class 2 is FDR-supported
# in CA2-SLM and passed the CA2-SLM robustness audit; class 3 has no
# FDR-supported context at all and was nominated on effect size and directional
# consistency only. A class-3 candidate is a planning suggestion, not a result.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

# One renderer serves both candidate panels. The screen that produced the
# source data sets these; the ten-candidate FDR panel is the default.
PANEL_SET <- Sys.getenv("IMMUNO_PANEL_SET", unset = "fdr")
CFG <- switch(PANEL_SET,
  fdr = list(src = "immunostaining_candidate_panel",
             out = "immunostaining_candidate_panel",
             abund = "candidate_panel_animal_abundance.csv",
             eff = "candidate_panel_sus_res_effects.csv",
             sel = "candidate_panel_selection.csv",
             stem = "immunostaining_panel"),
  separation = list(src = "immunostaining_separation_screen",
                    out = "immunostaining_separation_panel",
                    abund = "separation_panel_animal_abundance.csv",
                    eff = "separation_panel_sus_res_effects.csv",
                    sel = "separation_candidate_panel.csv",
                    stem = "immunostaining_separation"),
  stop("unknown IMMUNO_PANEL_SET: ", PANEL_SET, call. = FALSE))
SRC <- path_results("source_data", "10_biological_integration", CFG$src)
OUT <- path_results("figures", "manuscript", CFG$out)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

need <- function(f) {
  p <- file.path(SRC, f)
  if (!file.exists(p))
    stop("missing_required_input: ", p,
         "\nRun 10_biological_integration/11_screen_immunostaining_candidate_panel.R first.",
         call. = FALSE)
  utils::read.csv(p, stringsAsFactors = FALSE)
}
abundance <- need(CFG$abund)
effects   <- need(CFG$eff)
sel       <- need(CFG$sel)

PAL <- yaml::read_yaml(repo_path("config", "manuscript_palette.yml"))
GROUP_COL <- unlist(PAL$group)[c("CON", "RES", "SUS")]
BASE_PT <- 5.5

lev <- abundance %>% distinct(spatial_unit_display, spatial_order) %>%
  arrange(spatial_order) %>% pull(spatial_unit_display)

# facet label carries the evidence class, so a class-3 candidate can never be
# mistaken for an FDR-supported one at a glance
cls_short <- c("1_FDR_supported_outside_CA2SLM" = "class 1",
               "2_FDR_supported_CA2SLM_robustness_qualified" = "class 2",
               "3_no_FDR_support_effect_and_consistency_only" = "class 3")
sel <- sel %>% rename_with(function(z) ifelse(z == "rank", "panel_rank", z)) %>%
  arrange(panel_rank)
# The FDR panel labels by evidence class; the separation panel has no such
# column and is labelled by its rank and separation strength instead.
if ("evidence_class" %in% names(sel)) {
  tag <- unname(cls_short[sel$evidence_class])
  tag[is.na(tag)] <- sel$evidence_class[is.na(tag)]
  sel <- sel %>% mutate(facet_label = sprintf("%s (%s)  ·  %s",
                                              toupper(GeneSymbol), UniProt, tag))
} else {
  # Keep this short or the strip text is clipped: gene, SSMD, and a one-word
  # caution only where one applies.
  tag <- sprintf("SSMD %.1f", sel$best_ssmd)
  if ("housekeeping_flag" %in% names(sel))
    tag <- ifelse(as.logical(sel$housekeeping_flag),
                  paste0(tag, " · housekeeping?"), tag)
  if ("low_abundance_flag" %in% names(sel))
    tag <- ifelse(as.logical(sel$low_abundance_flag),
                  paste0(tag, " · low abundance"), tag)
  sel <- sel %>% mutate(facet_label = sprintf("%d. %s · %s", panel_rank,
                                              toupper(GeneSymbol), tag))
}
lab_map <- setNames(sel$facet_label, sel$GeneSymbol)
prot_lev <- sel$GeneSymbol

prep <- function(d) d %>%
  mutate(unit = factor(spatial_unit_display, levels = lev),
         protein = factor(GeneSymbol, levels = prot_lev))
abundance <- prep(abundance) %>%
  mutate(ExpGroup = factor(ExpGroup, levels = c("CON", "RES", "SUS")))
effects <- prep(effects)
stopifnot(identical(levels(abundance$unit), levels(effects$unit)))

group_mean <- abundance %>%
  group_by(protein, unit, ExpGroup) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop")

theme_ms <- function() {
  theme_classic(base_size = BASE_PT, base_family = "sans") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = BASE_PT, face = "bold", hjust = 0),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 size = BASE_PT - 1.5),
      axis.text.y = element_text(size = BASE_PT - 1.5),
      axis.title = element_text(size = BASE_PT),
      legend.position = "top", legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = BASE_PT),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2),
      panel.spacing = unit(1.6, "mm"),
      plot.margin = margin(2, 3, 2, 3, "mm"))
}

# mean drawn first so the animal points sit on top of it
pA <- ggplot(abundance, aes(unit, abundance, colour = ExpGroup)) +
  geom_point(data = group_mean, aes(unit, mean_abundance, colour = ExpGroup),
             position = position_dodge(width = 0.68), shape = 95, size = 3.0,
             alpha = 0.55, show.legend = FALSE) +
  geom_point(position = position_dodge(width = 0.68), size = 0.55,
             alpha = 0.95, stroke = 0) +
  facet_wrap(~ protein, nrow = 2, scales = "free_y",
             labeller = labeller(protein = lab_map)) +
  scale_colour_manual(values = GROUP_COL, name = NULL) +
  labs(x = NULL, y = "Animal-level abundance (log2)") +
  theme_ms()

# The separation panel has no FDR-supported point by construction - that is the
# whole point of that screen - so the marking layers are added only when there
# is something to mark.
eff_sup <- if ("fdr_supported" %in% names(effects) &&
               "fdr_bh" %in% names(effects)) {
  effects %>% filter(fdr_supported)
} else effects[0, , drop = FALSE]

pB <- ggplot(effects, aes(unit, log2FC)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey45") +
  geom_point(size = 0.9, colour = "grey30", stroke = 0)
if (nrow(eff_sup)) pB <- pB +
  geom_point(data = eff_sup, size = 1.7, shape = 21, fill = NA,
             colour = "#1F3D52", stroke = 0.45) +
  geom_text(data = eff_sup %>%
              mutate(hj = ifelse(as.integer(unit) > length(lev) * 2 / 3, 1,
                                 ifelse(as.integer(unit) <= length(lev) / 3, 0, 0.5)),
                     nudge = ifelse(hj == 1, -0.25, ifelse(hj == 0, 0.25, 0))),
            aes(x = as.integer(unit) + nudge, y = log2FC, hjust = hj,
                label = sprintf("FDR %.3g", fdr_bh)),
            vjust = -1.3, size = (BASE_PT - 2.2) / .pt, colour = "#1F3D52",
            inherit.aes = FALSE)

pB <- pB +
  facet_wrap(~ protein, nrow = 2, scales = "free_y",
             labeller = labeller(protein = lab_map)) +
  scale_y_continuous(expand = expansion(mult = c(0.16, 0.24))) +
  scale_x_discrete(expand = expansion(add = 0.7)) +
  labs(x = NULL, y = expression("SUS - RES log"[2]*" fold change")) +
  theme_ms() + theme(legend.position = "none")

W <- 183; HA <- 92; HB <- 86
save_svg <- function(p, stem, h) {
  svglite::svglite(file.path(OUT, paste0(stem, ".svg")), width = W / 25.4,
                   height = h / 25.4, bg = "white", fix_text_size = FALSE)
  print(p); grDevices::dev.off()
}
save_svg(pA, paste0(CFG$stem, "_abundance"), HA)
save_svg(pB, paste0(CFG$stem, "_sus_res_log2fc"), HB)

combined <- patchwork::wrap_plots(pA, pB, ncol = 1, heights = c(HA, HB)) +
  patchwork::plot_annotation(tag_levels = "A",
    theme = theme(plot.margin = margin(1, 1, 1, 1, "mm")))
stem <- paste0(CFG$stem, "_spatial_comparison")
save_svg(combined, stem, HA + HB)
grDevices::cairo_pdf(file.path(OUT, paste0(stem, ".pdf")), width = W / 25.4,
                     height = (HA + HB) / 25.4, bg = "white")
print(combined); grDevices::dev.off()
grDevices::png(file.path(OUT, paste0(stem, ".png")), width = W,
               height = HA + HB, units = "mm", res = 600, type = "cairo",
               bg = "white")
print(combined); grDevices::dev.off()

svg_main <- file.path(OUT, paste0(stem, ".svg"))
cat("immunostaining candidate panel ->", OUT, "\n")
cat("  candidates:", nrow(sel), "|",
    paste(sprintf("%s(%s)", sel$GeneSymbol, cls_short[sel$evidence_class]),
          collapse = " "), "\n")
cat("  animals:", length(unique(abundance$AnimalID)),
    "| points in panel A:", nrow(abundance),
    "| FDR-supported marked:", nrow(eff_sup), "\n")
cat("  main SVG:", basename(svg_main), "|", file.size(svg_main), "bytes\n")
