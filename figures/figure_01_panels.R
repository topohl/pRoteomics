#!/usr/bin/env Rscript
# ================================================================
# Script: figures/figure_01_panels.R
# Stage: manuscript_figures
# Scope: global
# Consumes: required manuscript/figure1_bridge_mmmsociability/source_data/*.csv
# Produces: results/figures/manuscript/figure_01_panels/
# Notes: Rendering only. Reads frozen imported source data; computes no statistic.
# ================================================================
#
# Figure 1 panels. The behavioural analysis lives in a separate repository and is
# frozen; this script plots values that were computed there and imported with
# hashes. It fits no correlation, no model and runs no permutation. Every number
# it prints is read from figure1_panel_statistics.csv, which is the contract for
# what the figure is allowed to say - a value absent from that index cannot
# appear on the figure, because there is nowhere for this script to get it.
#
# VISUAL LANGUAGE. This layout descends from the original behaviour main figure,
# whose closest surviving ancestor is
# MMMSociability Analysis/27_candidate_recompose_behavior_main_figure.R at commit
# 4b0f90f. Recovered from it verbatim: the held-out prediction panel with its
# permutation-null companion, the dashed identity line, the circle/triangle sex
# encoding, and the point and histogram marks. The original's green timeline
# schematic and its six-component forest have no source in either repository at
# any commit, so the timeline here is authored fresh in that visual language and
# the forest is deliberately not reproduced: those six components construct the
# outcome and their group separation is true by construction.
#
# Panels:
#   a  design and temporal ordering - the predictor window closes before any
#      outcome component is measured, so the classification did not exist yet
#   b  CombZ and how RES/SUS was defined from it (definition, not validation)
#   c  early mean movement versus later CombZ, the primary association
#   d  held-out prediction with the full-refit permutation null beside it
#
# Four panels, not five: the prediction and its null are one argument and the
# original presented them together. The repeated grouped five-fold result is
# annotated rather than drawn, because upstream persists summary rows only, so no
# per-repeat distribution exists and none may be invented.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

SRC <- repo_path("manuscript", "figure1_bridge_mmmsociability", "source_data")
OUT <- path_results("figures", "manuscript", "figure_01_panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

need <- function(f) {
  p <- file.path(SRC, f)
  if (!file.exists(p))
    stop("missing_required_input: ", p,
         "\nThe frozen Figure 1 panel source data must be imported first.",
         call. = FALSE)
  utils::read.csv(p, stringsAsFactors = FALSE)
}
tl    <- need("figure1a_timeline_source.csv")
cls   <- need("figure1b_combz_classification_source.csv")
assoc <- need("figure1c_movement_combz_source.csv")
loao  <- need("figure1d_loao_predictions_source.csv")
perm  <- need("figure1e_permutation_source.csv")
rcv   <- need("figure1f_repeated_cv_source.csv")
STATS <- need("figure1_panel_statistics.csv")

# Every printed number resolves through here. There is deliberately no fallback:
# if a statistic is not in the frozen index, the render fails rather than
# inventing or recomputing it.
stat <- function(panel, statistic) {
  v <- STATS$value[STATS$figure_panel == panel & STATS$statistic == statistic]
  if (length(v) != 1L)
    stop("figure1_panel_statistics.csv does not resolve ", panel, " / ",
         statistic, call. = FALSE)
  v
}
statn <- function(panel, statistic) as.numeric(stat(panel, statistic))

# ---------------------------------------------------------------- house style
#
# Group colours come from R/plotting_nature.R, which is what Figures 2 and 3
# actually render with and which is byte-identical to MMM_GROUP_COLOURS in the
# behavioural repository. config/manuscript_palette.yml declares a different
# set that no numbered figure has ever used; this figure previously obeyed it
# and was the only file in the repository doing so.
source(repo_path("R", "plotting_nature.R"))
GROUP_COL <- NATURE_SEMANTIC_PALETTES$group[c("CON", "RES", "SUS")]
GROUP_LEV <- c("CON", "RES", "SUS")
# RES is deliberately light, so colour is never the only group encoding: shape
# carries the same information redundantly, as upstream requires.
GROUP_SHAPE <- c(CON = 21L, RES = 24L, SUS = 22L)
SEX_SHAPE <- c(Female = 21L, Male = 24L)

INK <- "#2B2B2B"; MUTED <- "#6E6E6E"; RULE <- "#B5B5B5"
GREEN_DARK <- "#2F6F62"   # paper's own dataset green, used by Figures 2 and 3
GREEN_MID <- "#6FA79B"
GREEN_LIGHT <- "#A8D5CF"
NEUTRAL_BOX <- "#EDEDED"

# Authored at the exact size the assembler places, so scale = 1.0 and the point
# sizes set here are the ones the reader sees. The assembler reserves 4 mm of
# each row for the panel letter, so a 54 mm row places a 50 mm image box.
BASE_PT <- 7      # titles, axis titles, strip labels
BODY_PT <- 6.5    # tick labels and in-panel statistics
NOTE_PT <- 6      # subordinate captions; nothing may print below this

theme_f1 <- function(base = BASE_PT) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.text = element_text(size = BODY_PT, colour = INK),
      axis.title = element_text(size = base, colour = INK),
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(size = base, colour = INK, face = "italic",
                                margin = margin(b = 1.5, t = 1)),
      legend.position = "top",
      legend.justification = "left",
      legend.key.height = unit(7, "pt"),
      legend.key.width = unit(11.2, "pt"),
      legend.text = element_text(size = NOTE_PT),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, -3.5, 0),
      panel.grid = element_blank(),
      plot.title = element_text(size = base, colour = INK, hjust = 0,
                                face = "plain", margin = margin(b = 0.5)),
      plot.subtitle = element_text(size = NOTE_PT, colour = "grey25",
                                   margin = margin(b = 1.5)),
      plot.margin = margin(3.5, 3.5, 3.5, 3.5))
}

# =============================================================== panel a
# Design schematic in the original figure's visual language: a stage band above a
# black developmental timeline with circular milestone markers. The predictor
# stage is the only saturated block, so the eye lands on it first; everything
# downstream of it is desaturated. No source for the original exists, so this is
# authored, but the greens are the paper's own dataset greens rather than new ink.
# Stage content is DERIVED FROM THE FROZEN TIMELINE, not from the original
# artwork. The original showed an RFID implant and a tissue-collection endpoint;
# neither is an event in figure1a_timeline_source.csv, so neither is drawn here.
# Only the wording is presentational, exactly as a display-label map is.
DISPLAY <- c(
  "first cage change (CC1)"        = "Cage change\nCC1",
  "early RFID window"              = "Early RFID window\n18:30-06:30, 12 h",
  "social-instability paradigm"    = "Social instability stress\nrepeated regrouping",
  "later behavioural outcome tests" = "Behavioural\noutcome tests",
  "terminal physiology"            = "Terminal\nphysiology",
  "CombZ computed"                 = "CombZ\ncomputed",
  "RES/SUS assigned"               = "Resilient /\nsusceptible")
missing_events <- setdiff(names(DISPLAY), tl$event)
if (length(missing_events))
  stop("the frozen timeline no longer contains: ",
       paste(missing_events, collapse = ", "), call. = FALSE)
tl <- tl[match(names(DISPLAY), tl$event), , drop = FALSE]

# The predictor stage is the only saturated block; everything downstream of it
# is desaturated, so the eye lands on the recorded window first.
n <- nrow(tl)
STAGE <- data.frame(
  x0 = seq_len(n) - 0.44,
  x1 = seq_len(n) + 0.44,
  label = unname(DISPLAY),
  # Saturation falls with distance from the predictor: the reference event and
  # the recorded window carry the colour, the stressor sits mid, and everything
  # that happens after the window closes is deliberately quiet.
  fill = ifelse(tl$used_as_predictor, GREEN_MID,
                ifelse(tl$used_in_outcome, NEUTRAL_BOX, GREEN_LIGHT)),
  ink = ifelse(tl$used_as_predictor, "white", INK),
  stringsAsFactors = FALSE)
STAGE$fill[1] <- GREEN_DARK
STAGE$ink[1] <- "white"
MILE <- data.frame(
  x = seq_len(n),
  lab = c("P25", "", "", "", "", "", "labels\nassigned"),
  stringsAsFactors = FALSE)
PRED_I <- which(tl$used_as_predictor)
OUT_I <- which(tl$used_in_outcome)

cov_txt <- sprintf(
  "%d expected 10-min slots per animal; %d of %d complete, the rest missing leading slots only (mean coverage %.1f%%)",
  as.integer(statn("1a", "expected 10-min slots per animal")),
  as.integer(statn("1a", "animals with all expected slots")),
  as.integer(tl$n_animals[!is.na(tl$n_animals)][1]),
  100 * statn("1a", "mean coverage fraction"))

pa <- ggplot() +
  geom_rect(data = STAGE, aes(xmin = x0, xmax = x1, ymin = 0.30, ymax = 1.15,
                              fill = I(fill)), colour = NA) +
  geom_text(data = STAGE, aes(x = (x0 + x1) / 2, y = 0.725, label = label,
                              colour = I(ink)),
            size = BODY_PT / .pt, lineheight = 0.95) +
  # developmental timeline
  annotate("segment", x = 0.5, xend = n + 0.5, y = 0, yend = 0,
           linewidth = 0.55, colour = "black") +
  geom_point(data = MILE, aes(x, 0), size = 1.5, colour = "black") +
  geom_text(data = MILE, aes(x, -0.20, label = lab), size = NOTE_PT / .pt,
            colour = INK, lineheight = 0.95, vjust = 1) +
  # The load-bearing statement of the panel: the predictor bracket closes before
  # the outcome bracket opens, and both are derived from the frozen flags.
  annotate("segment", x = min(PRED_I) - 0.44, xend = max(PRED_I) + 0.44,
           y = 1.34, yend = 1.34, linewidth = 0.45, colour = GREEN_DARK) +
  annotate("text", x = mean(PRED_I), y = 1.52, label = "predictor recorded here",
           size = BODY_PT / .pt, colour = GREEN_DARK) +
  annotate("segment", x = min(OUT_I) - 0.44, xend = max(OUT_I) + 0.44,
           y = 1.34, yend = 1.34, linewidth = 0.45, colour = MUTED) +
  annotate("text", x = mean(OUT_I), y = 1.52,
           label = "outcome, CombZ and group labels derived here",
           size = BODY_PT / .pt, colour = MUTED) +
  annotate("text", x = 0.5, y = 1.86, hjust = 0, fontface = "italic",
           label = "the resilient / susceptible label does not exist while the predictor is recorded",
           size = BODY_PT / .pt, colour = INK) +
  scale_x_continuous(limits = c(0.4, n + 0.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.72, 2.02), expand = c(0, 0)) +
  labs(caption = cov_txt) +
  theme_void(base_size = BODY_PT, base_family = "sans") +
  theme(plot.caption = element_text(size = NOTE_PT - 0.5, colour = MUTED,
                                    hjust = 0, margin = margin(t = 0.5)),
        plot.margin = margin(3.5, 3.5, 3.5, 3.5))

# =============================================================== panel b
# Phenotype definition. The threshold is the content; the group separation is
# true by construction and must not read as a discovered difference. Hence a
# title that says so, no brackets, no stars, and points that are quieter than
# the rule that produced them.
cls$Sex <- factor(cls$Sex, levels = c("Female", "Male"))
cls$final_group <- factor(cls$final_group, levels = GROUP_LEV)
thr <- cls %>% distinct(Sex, susceptibility_threshold, control_reference_mean)
set.seed(1)  # jitter is presentational only and affects no value
pb <- ggplot(cls, aes(final_group, CombZ)) +
  geom_hline(data = thr, aes(yintercept = control_reference_mean),
             linewidth = 0.25, colour = RULE) +
  geom_hline(data = thr, aes(yintercept = susceptibility_threshold),
             linewidth = 0.45, colour = INK, linetype = "22") +
  geom_point(aes(fill = final_group, shape = final_group),
             position = position_jitter(width = 0.20, height = 0),
             size = 1.05, stroke = 0.22, colour = "grey20", alpha = 0.85) +
  facet_wrap(~ Sex, nrow = 1) +
  scale_fill_manual(values = GROUP_COL, guide = "none") +
  scale_shape_manual(values = GROUP_SHAPE, guide = "none") +
  geom_text(data = thr, aes(x = 0.40, y = susceptibility_threshold,
                            label = sprintf("%.3f", susceptibility_threshold)),
            inherit.aes = FALSE, hjust = 0, vjust = 1.4,
            size = NOTE_PT / .pt, colour = INK) +
  scale_x_discrete(expand = expansion(add = c(0.62, 0.35))) +
  labs(x = NULL, y = "Later CombZ",
       title = "Later composite outcome and phenotype definition",
       # The separation here is produced by the dashed line, not discovered.
       # Saying so is the cheapest defence against reading this panel as a test,
       # and it is why there are no brackets and no stars.
       subtitle = "the threshold defines these groups; it is not a test of them") +
  theme_f1()

# =============================================================== panel c
# The primary association. Combined, never split by sex: the formal
# feature-by-sex interaction is unsupported at all q = 0.90. Group is shown for
# orientation only and never entered the model, so it is a quiet fill under a
# common outline rather than a second signal.
assoc$Group <- factor(assoc$Group, levels = GROUP_LEV)
c_lab <- sprintf("Spearman ρ = %.2f\n95%% CI [%.2f, %.2f]\nBH q = %.2e, n = %d",
                 statn("1c", "Spearman rho"), statn("1c", "95% CI lower"),
                 statn("1c", "95% CI upper"), statn("1c", "BH q"),
                 as.integer(statn("1c", "n animals")))
# Single colour, deliberately, and this is the one place the original's group
# colouring is not recovered. Susceptible animals sit low on this y axis and
# resilient animals high BY CONSTRUCTION, because the axis is the score the
# groups were cut from. Colouring by group here would let a reader see the
# correlation as group separation - the strongest available misreading of the
# figure's central claim. The association is between two continuous variables
# over 111 animals and group entered nothing.
pc <- ggplot(assoc, aes(Movement_mean, CombZ)) +
  geom_point(size = 1.15, stroke = 0.2, alpha = 0.9, shape = 21,
             colour = "grey20", fill = "#6E8B99") +
  annotate("text", x = Inf, y = Inf, label = c_lab, hjust = 1.04, vjust = 1.18,
           size = NOTE_PT / .pt, colour = INK, lineheight = 1.08) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = "Early mean movement (transitions / 10 min)",
       y = "Later CombZ",
       subtitle = "one point per animal; rank correlation, no model fitted, not split by sex") +
  theme_f1()

# =============================================================== panel d
# Held-out prediction and its null, recovered from the original design: dashed
# identity line, circle/triangle sex encoding, group fill, and the permutation
# null beside it rather than as a separate panel. coord_equal is the original's
# own choice and is what makes the identity line read at 45 degrees.
loao$Group <- factor(loao$Group, levels = GROUP_LEV)
loao$Sex <- factor(loao$Sex, levels = c("Female", "Male"))
lim <- range(c(loao$observed_CombZ, loao$predicted_CombZ), na.rm = TRUE)
d_lab <- sprintf("LOAO R² = %.3f\nbaseline R² = %.3f",
                 statn("1d", "LOAO R2"), statn("1d", "intercept-only baseline R2"))

pd_scatter <- ggplot(loao, aes(observed_CombZ, predicted_CombZ)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.35,
              colour = "grey62") +
  geom_point(aes(fill = Group, shape = Sex), size = 1.0, stroke = 0.24,
             colour = "grey20", alpha = 0.86) +
  scale_fill_manual(values = GROUP_COL, name = NULL) +
  scale_shape_manual(values = SEX_SHAPE, name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  annotate("text", x = -Inf, y = Inf, label = d_lab, hjust = -0.06, vjust = 1.2,
           size = NOTE_PT / .pt, colour = INK, lineheight = 1.08) +
  guides(fill = guide_legend(order = 1, override.aes =
                               list(shape = 21, size = 1.5, alpha = 1)),
         shape = guide_legend(order = 2, override.aes =
                                list(fill = "grey70", size = 1.5))) +
  labs(x = "Observed CombZ", y = "Held-out predicted CombZ",
       title = "Held-out prediction of later CombZ and its permutation null",
       subtitle = "movement-mean model") +
  theme_f1()

nulls <- perm[perm$row_role == "permutation_draw", , drop = FALSE]
obs <- statn("1e", "observed LOAO R2")
e_sub <- sprintf("repeated grouped five-fold: R² = %.3f (%.3f-%.3f across repeats)",
                 statn("1f", "repeated grouped CV mean R2"),
                 statn("1f", "2.5th percentile across repeats"),
                 statn("1f", "97.5th percentile across repeats"))
pd_null <- ggplot(nulls, aes(performance_value)) +
  geom_histogram(bins = 24, fill = "grey78", colour = "white",
                 linewidth = 0.15) +
  geom_vline(xintercept = obs, linewidth = 0.45, colour = GROUP_COL[["SUS"]]) +
  annotate("text", x = obs, y = Inf, vjust = 1.3, hjust = 1.08,
           label = sprintf("observed R² = %.3f\np = %s", obs,
                           stat("1e", "empirical p")),
           size = NOTE_PT / .pt, colour = INK, lineheight = 1.08) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.10, 0.15)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = expression("null R"^2*" (permuted outcomes)"),
       y = sprintf("Permutations (n = %d)", as.integer(statn("1e", "permutation draws"))),
       subtitle = e_sub) +
  theme_f1()

# The scatter is square by coord_equal, so its column is sized to what a square
# plot plus its axes actually needs; giving it half the panel would leave a wide
# internal gap before the null.
pd <- patchwork::wrap_plots(pd_scatter, pd_null, nrow = 1, widths = c(0.62, 1))

# ---------------------------------------------------------------- write out
# Authored at the exact placed size, so scale = 1.0 and native pt == printed pt.
# Three equal rows of 54 mm place 50 mm image boxes; a and d span both columns.
DIMS <- list(a = c(173, 50), b = c(84.5, 50), c = c(84.5, 50), d = c(173, 50))
PLOTS <- list(a = pa, b = pb, c = pc, d = pd)

save_svg <- function(p, id) {
  d <- DIMS[[id]]
  f <- file.path(OUT, paste0("figure_01", id, ".svg"))
  svglite::svglite(f, width = d[1] / 25.4, height = d[2] / 25.4, bg = "white",
                   fix_text_size = FALSE)
  print(p); grDevices::dev.off()
  f
}
written <- vapply(names(PLOTS), function(id) save_svg(PLOTS[[id]], id),
                  character(1))

# Panel e no longer exists: its content is now the right half of panel d. Remove
# a stale export rather than leave a file that still looks canonical.
stale <- file.path(OUT, "figure_01e.svg")
if (file.exists(stale)) {
  file.remove(stale)
  cat("removed superseded panel export: figure_01e.svg\n")
}

cat("Figure 1 panels ->", OUT, "\n")
for (id in names(written))
  cat(sprintf("  figure_01%s.svg  %6d bytes  %.1f x %.0f mm\n", id,
              file.size(written[[id]]), DIMS[[id]][1], DIMS[[id]][2]))
cat("  animals: 1b", nrow(cls), "| 1c", nrow(assoc), "| 1d", nrow(loao),
    "| null draws", nrow(nulls), "\n")
cat("  group palette: CON", GROUP_COL[["CON"]], "RES", GROUP_COL[["RES"]],
    "SUS", GROUP_COL[["SUS"]], "(NATURE_SEMANTIC_PALETTES, as Figures 2 and 3)\n")
cat("  every printed number resolved through figure1_panel_statistics.csv\n")
