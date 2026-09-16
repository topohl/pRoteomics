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
# Panels:
#   a  design and temporal ordering - the predictor window closes before any
#      outcome component is measured, so the classification did not exist yet
#   b  CombZ and how RES/SUS was defined from it (definition, not validation)
#   c  early mean movement versus later CombZ, the primary association
#   d  observed versus leave-one-animal-out predicted CombZ
#   e  full-refit permutation null of the LOAO R2, with the observed value
#
# The repeated grouped five-fold result is annotated on panel e rather than
# given its own panel: the upstream analysis persists summary rows only, so no
# per-repeat distribution exists to draw and none may be invented.

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
tl   <- need("figure1a_timeline_source.csv")
cls  <- need("figure1b_combz_classification_source.csv")
assoc <- need("figure1c_movement_combz_source.csv")
loao <- need("figure1d_loao_predictions_source.csv")
perm <- need("figure1e_permutation_source.csv")
rcv  <- need("figure1f_repeated_cv_source.csv")
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

PAL <- yaml::read_yaml(repo_path("config", "manuscript_palette.yml"))
GROUP_COL <- unlist(PAL$group)[c("CON", "RES", "SUS")]
INK <- "#1A1A1A"; MUTED <- "#6E6E6E"; RULE <- "#B5B5B5"
BASE_PT <- 6

theme_f1 <- function() {
  theme_classic(base_size = BASE_PT, base_family = "sans") +
    theme(
      axis.text = element_text(size = BASE_PT - 1, colour = INK),
      axis.title = element_text(size = BASE_PT, colour = INK),
      axis.line = element_line(linewidth = 0.3, colour = INK),
      axis.ticks = element_line(linewidth = 0.3, colour = INK),
      strip.background = element_blank(),
      strip.text = element_text(size = BASE_PT, colour = INK, face = "plain"),
      legend.position = "top",
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = BASE_PT - 1),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      panel.grid = element_blank(),
      plot.margin = margin(2, 3, 2, 3, "mm"),
      plot.subtitle = element_text(size = BASE_PT - 1, colour = MUTED))
}

# ------------------------------------------------------------------ panel a
# A design schematic, not a data plot. Temporal precedence is the content: the
# recorded window sits entirely to the left of every outcome component.
lane <- data.frame(
  x = c(1, 2, 3, 4, 5, 6, 7),
  label = c("Cage change\n(CC1, P25)", "Early RFID window\n18:30-06:30, 12 h",
            "Social-instability\nparadigm", "Behavioural\noutcome tests",
            "Terminal\nphysiology", "CombZ\ncomputed", "Resilient /\nsusceptible"),
  role = c("reference", "predictor", "paradigm", "outcome", "outcome",
           "derived", "derived"),
  stringsAsFactors = FALSE)
ROLE_FILL <- c(reference = "#FFFFFF", predictor = "#DCE7EC",
               paradigm = "#F2F2F2", outcome = "#F7EDE6", derived = "#EDEDED")
cov_txt <- sprintf(
  "%d expected 10-min slots per animal; %d of %d animals complete, the rest missing leading slots only (mean coverage %.1f%%)",
  as.integer(statn("1a", "expected 10-min slots per animal")),
  as.integer(statn("1a", "animals with all expected slots")),
  as.integer(tl$n_animals[!is.na(tl$n_animals)][1]),
  100 * statn("1a", "mean coverage fraction"))

pa <- ggplot(lane, aes(x, 0)) +
  annotate("segment", x = 0.55, xend = 7.45, y = -0.84, yend = -0.84,
           linewidth = 0.3, colour = MUTED,
           arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_tile(aes(fill = role), width = 0.86, height = 0.72, colour = INK,
            linewidth = 0.25, show.legend = FALSE) +
  scale_fill_manual(values = ROLE_FILL) +
  geom_text(aes(label = label), size = (BASE_PT - 1.4) / .pt, colour = INK,
            lineheight = 0.95) +
  # The outcome bracket sits ABOVE the boxes and the predictor marker below, so
  # neither can collide with the box band, which spans -0.36 to 0.36.
  annotate("segment", x = 4, xend = 7, y = 0.48, yend = 0.48,
           linewidth = 0.25, colour = MUTED) +
  annotate("segment", x = c(4, 7), xend = c(4, 7), y = 0.42, yend = 0.48,
           linewidth = 0.25, colour = MUTED) +
  annotate("text", x = 5.5, y = 0.57, label = "outcome measured and derived here",
           size = (BASE_PT - 1.6) / .pt, colour = MUTED) +
  annotate("segment", x = 2, xend = 2, y = -0.38, yend = -0.54,
           linewidth = 0.25, colour = INK) +
  annotate("text", x = 2, y = -0.64, label = "predictor recorded here",
           size = (BASE_PT - 1.6) / .pt, colour = INK) +
  annotate("text", x = 4, y = 0.82,
           label = "the resilient / susceptible label does not exist during the recorded window",
           size = (BASE_PT - 1.6) / .pt, colour = INK, fontface = "italic") +
  scale_y_continuous(limits = c(-1.02, 0.98), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.4, 7.6), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, caption = cov_txt) +
  theme_void(base_size = BASE_PT, base_family = "sans") +
  theme(plot.caption = element_text(size = BASE_PT - 2, colour = MUTED,
                                    hjust = 0, margin = margin(t = 1)),
        plot.margin = margin(2, 3, 2, 3, "mm"))

# ------------------------------------------------------------------ panel b
# How the phenotype was defined. The threshold is the point of the panel; the
# six components that build CombZ are deliberately not shown, because their
# group differences are true by construction.
cls$Sex <- factor(cls$Sex, levels = c("Female", "Male"))
cls$final_group <- factor(cls$final_group, levels = c("CON", "RES", "SUS"))
thr <- cls %>% distinct(Sex, susceptibility_threshold, control_reference_mean)
set.seed(1)  # jitter is presentational only and affects no value
pb <- ggplot(cls, aes(final_group, CombZ, colour = final_group)) +
  geom_hline(data = thr, aes(yintercept = control_reference_mean),
             linewidth = 0.25, colour = RULE) +
  geom_hline(data = thr, aes(yintercept = susceptibility_threshold),
             linewidth = 0.35, colour = INK, linetype = "22") +
  geom_point(position = position_jitter(width = 0.22, height = 0),
             size = 0.5, alpha = 0.9, stroke = 0) +
  facet_wrap(~ Sex, nrow = 1) +
  scale_colour_manual(values = GROUP_COL, guide = "none") +
  geom_text(data = thr, aes(x = 0.62, y = susceptibility_threshold,
                            label = sprintf("%.3f", susceptibility_threshold)),
            inherit.aes = FALSE, hjust = 0, vjust = 1.45,
            size = (BASE_PT - 2) / .pt, colour = INK) +
  labs(x = NULL, y = "Later composite outcome (CombZ)",
       subtitle = "dashed: same-sex control mean minus one control population SD; grey: control mean") +
  theme_f1()

# ------------------------------------------------------------------ panel c
# The primary association. Deliberately unstratified and deliberately
# single-colour: the formal sex interaction is unsupported, and colouring by
# outcome group would invite the reading that group entered the model.
c_rho <- statn("1c", "Spearman rho")
c_lo <- statn("1c", "95% CI lower"); c_hi <- statn("1c", "95% CI upper")
c_q <- statn("1c", "BH q"); c_n <- as.integer(statn("1c", "n animals"))
c_lab <- sprintf("Spearman ρ = %.2f  [%.2f, %.2f]\nq = %.1e   n = %d",
                 c_rho, c_lo, c_hi, c_q, c_n)
pc <- ggplot(assoc, aes(Movement_mean, CombZ)) +
  geom_point(size = 0.6, alpha = 0.85, stroke = 0, colour = "#3D5A73") +
  annotate("text", x = Inf, y = Inf, label = c_lab, hjust = 1.05, vjust = 1.25,
           size = (BASE_PT - 1.6) / .pt, colour = INK, lineheight = 1.05) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = "Early mean movement (first 12 h active window)",
       y = "Later composite outcome (CombZ)",
       subtitle = "higher CombZ = more resilient-like; rank correlation, no model fitted here") +
  theme_f1()

# ------------------------------------------------------------------ panel d
# Held-out predictions. Each point is an animal predicted by a model refitted
# without it. The reference line is identity, not a fit.
d_r2 <- statn("1d", "LOAO R2"); d_base <- statn("1d", "intercept-only baseline R2")
d_n <- as.integer(statn("1d", "n held-out animals"))
d_lab <- sprintf("LOAO R² = %.3f\nintercept-only baseline R² = %.3f\nn = %d animals",
                 d_r2, d_base, d_n)
lim <- range(c(loao$observed_CombZ, loao$predicted_CombZ))
pd <- ggplot(loao, aes(observed_CombZ, predicted_CombZ)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = RULE) +
  geom_point(size = 0.6, alpha = 0.85, stroke = 0, colour = "#3D5A73") +
  annotate("text", x = -Inf, y = Inf, label = d_lab, hjust = -0.06, vjust = 1.2,
           size = (BASE_PT - 1.6) / .pt, colour = INK, lineheight = 1.05) +
  coord_cartesian(xlim = lim, ylim = lim) +
  labs(x = "Observed CombZ", y = "Predicted CombZ (held out)",
       subtitle = "movement-mean model; every animal predicted by a model fitted without it") +
  theme_f1()

# ------------------------------------------------------------------ panel e
# The real persisted null, not a density reconstructed from quantiles.
nulls <- perm[perm$row_role == "permutation_draw", , drop = FALSE]
obs <- statn("1e", "observed LOAO R2")
e_p <- stat("1e", "empirical p")
e_draws <- as.integer(statn("1e", "permutation draws"))
rcv_h <- rcv[rcv$model_id == "movement_mean", ][1, ]
e_sub <- sprintf(
  "repeated grouped five-fold companion: R² = %.3f
2.5th-97.5th percentile range across %d repeats %.3f-%.3f (not a confidence interval)",
  statn("1f", "repeated grouped CV mean R2"),
  as.integer(rcv_h$repeated_cv_repeats),
  statn("1f", "2.5th percentile across repeats"),
  statn("1f", "97.5th percentile across repeats"))
pe <- ggplot(nulls, aes(performance_value)) +
  geom_histogram(bins = 46, fill = "#C9D3DA", colour = "#FFFFFF",
                 linewidth = 0.12) +
  geom_vline(xintercept = obs, linewidth = 0.4, colour = GROUP_COL[["SUS"]]) +
  annotate("text", x = obs, y = Inf, vjust = 1.4, hjust = 1.08,
           label = sprintf("observed R² = %.3f\np = %s", obs, e_p),
           size = (BASE_PT - 1.6) / .pt, colour = INK, lineheight = 1.05) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression("Leave-one-animal-out R"^2*" under permuted outcomes"),
       y = sprintf("Permutation draws (n = %d)", e_draws),
       subtitle = e_sub) +
  theme_f1()

# ---------------------------------------------------------------- write out
W <- 183
DIMS <- list(a = c(W, 46), b = c(W / 2, 50), c = c(W / 2, 50),
             d = c(W / 2, 50), e = c(W / 2, 50))
PLOTS <- list(a = pa, b = pb, c = pc, d = pd, e = pe)

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

cat("Figure 1 panels ->", OUT, "\n")
for (id in names(written))
  cat(sprintf("  figure_01%s.svg  %6d bytes  %.0f x %.0f mm\n", id,
              file.size(written[[id]]), DIMS[[id]][1], DIMS[[id]][2]))
cat("  animals: 1b", nrow(cls), "| 1c", nrow(assoc), "| 1d", nrow(loao),
    "| 1e draws", nrow(nulls), "\n")
cat("  every printed number resolved through figure1_panel_statistics.csv\n")
