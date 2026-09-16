#!/usr/bin/env Rscript
# ================================================================
# Script: figures/extended_data_behaviour_panels.R
# Stage: manuscript_figures
# Scope: global
# Consumes: required manuscript/figure1_bridge_mmmsociability/behavior_prediction_model_ladder.csv;
#           required manuscript/figure1_bridge_mmmsociability/behavior_sex_effect_contract.csv
# Produces: results/figures/manuscript/extended_data_behaviour_panels/
# Notes: Rendering only. Reads frozen imported contracts; computes no statistic.
# ================================================================
#
# Behavioural Extended Data. Figure 1 shows one model and one association. This
# figure shows the two things a reader needs in order to know that neither was
# chosen after the fact: the complete a-priori model ladder, and the formal
# feature-by-sex interaction tests.
#
# Both tables are frozen upstream in topohl/MMMSociability and imported
# byte-identically with SHA-256 into manuscript/figure1_bridge_mmmsociability/.
# Nothing here is fitted, resampled or corrected. Every value printed is a cell
# of one of those two tables.
#
# Panels:
#   a  the five a-priori models: leave-one-animal-out R2, the repeated grouped
#      five-fold interval, and the permutation P where one was run
#   b  the three feature-by-sex interaction tests with 95% CI and BH q
#   c  the sex-stratified correlations, drawn as descriptive context only
#
# WORDING DISCIPLINE. The sex panels exist to license a negative statement and
# to block a positive one. The upstream contract carries an explicit prohibited
# wording list, and this script refuses to render if every interaction row is
# not classified FORMAL_INTERACTION_NOT_SUPPORTED - because the moment one is,
# panel c stops being descriptive context and the figure has to be rethought.
#
# The ladder panel likewise states no winner. The models differ by amounts far
# smaller than their own intervals, and no test compares them; the panel is
# ordered by the fixed registry, never by performance.

suppressPackageStartupMessages({
  library(ggplot2)
})

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

BRIDGE <- repo_path("manuscript", "figure1_bridge_mmmsociability")
OUT <- path_results("figures", "manuscript", "extended_data_behaviour_panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

need <- function(f) {
  p <- file.path(BRIDGE, f)
  if (!file.exists(p))
    stop("missing_required_input: ", p,
         "\nThe frozen behavioural contracts must be imported first.",
         call. = FALSE)
  utils::read.csv(p, stringsAsFactors = FALSE)
}
ladder <- need("behavior_prediction_model_ladder.csv")
sexc <- need("behavior_sex_effect_contract.csv")

# ------------------------------------------------------------------ contracts
#
# These are not defensive niceties. Each one guards a sentence the manuscript is
# allowed to write only while the guarded condition holds.

REQ_LADDER <- c("model_id", "features", "sex_scope", "cv_design", "loao_r2",
                "repeated_cv_mean_r2", "cv_r2_q025", "cv_r2_q975",
                "permutation_p", "n", "primary_or_exploratory", "current_status")
absent <- setdiff(REQ_LADDER, names(ladder))
if (length(absent))
  stop("behavior_prediction_model_ladder.csv is missing: ",
       paste(absent, collapse = ", "), call. = FALSE)

REQ_SEX <- c("feature", "interaction_estimate", "interaction_ci_low",
             "interaction_ci_high", "interaction_p", "interaction_q_bh",
             "female_rho", "male_rho", "classification", "allowed_wording",
             "prohibited_wording")
absent <- setdiff(REQ_SEX, names(sexc))
if (length(absent))
  stop("behavior_sex_effect_contract.csv is missing: ",
       paste(absent, collapse = ", "), call. = FALSE)

# Panel c is drawn as descriptive context. That is only honest while no formal
# interaction is supported. If upstream ever changes one, this figure must be
# redesigned rather than quietly redrawn.
if (!all(sexc$classification == "FORMAL_INTERACTION_NOT_SUPPORTED"))
  stop("At least one feature-by-sex interaction is no longer classified ",
       "FORMAL_INTERACTION_NOT_SUPPORTED. Panel c presents stratified ",
       "correlations as descriptive context, which is no longer defensible. ",
       "Redesign the panel; do not re-render it.", call. = FALSE)

# The manuscript's headline model must remain the one Figure 1 draws.
headline <- ladder$model_id[ladder$current_status == "HEADLINE PRIMARY"]
if (length(headline) != 1L || headline != "movement_mean")
  stop("The ladder no longer names movement_mean as the single headline ",
       "primary model; Figure 1 and this panel would disagree.", call. = FALSE)

if (length(unique(ladder$n)) != 1L)
  stop("The ladder models are not all fitted on the same animals.", call. = FALSE)
N_ANIMALS <- unique(ladder$n)

# ---------------------------------------------------------------- house style
# Identical to figures/figure_01_panels.R, because this figure is read beside it.
source(repo_path("R", "plotting_nature.R"))

INK <- "#2B2B2B"; MUTED <- "#6E6E6E"; RULE <- "#B5B5B5"
GREEN_DARK <- "#2F6F62"
GREEN_MID <- "#6FA79B"
GREEN_LIGHT <- "#A8D5CF"

BASE_PT <- 7
BODY_PT <- 6.5
NOTE_PT <- 6

theme_ed <- function(base = BASE_PT) {
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
      # Titles align to the figure, not to the plot panel. Every panel here is a
      # horizontal forest whose y-axis labels are long sentences, so a
      # panel-aligned title would start a third of the way across and sit
      # nowhere near the panel letter. Figure 1 keeps the ggplot default because
      # none of its panels has a label column of that width.
      plot.title.position = "plot",
      plot.caption.position = "plot",
      # Indented past the 7 mm panel-letter box the assembler paints over the
      # top-left corner of every panel. Figure-aligned titles start exactly
      # where that box sits, so without this the first word is covered.
      plot.title = element_text(size = base, colour = INK, hjust = 0,
                                face = "plain", margin = margin(b = 0.5, l = 22)),
      plot.subtitle = element_text(size = NOTE_PT, colour = "grey25",
                                   margin = margin(b = 1.5, l = 22)),
      plot.margin = margin(3.5, 3.5, 3.5, 3.5))
}

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

# ================================================================== panel a
#
# The ladder is drawn in registry order, bottom to top, so the reader sees the
# reference baseline first and the fixed order is visible as order. The point is
# the leave-one-animal-out estimate, which is the primary quantity; the bar is
# the repeated grouped five-fold interval, which belongs to the repeated-CV mean
# and NOT to the point. They are drawn in different ink for that reason.

# The model's role is carried in its own axis label rather than in a legend.
# A legend would cost three of the panel's 48 mm to encode five rows that are
# already named on the axis, and the roles are the reason the rows are ordered
# the way they are, so they belong beside the names.
MODEL_LABEL <- c(
  mean_only = "Intercept only  (baseline)",
  movement_mean = "Mean movement",
  primary_behavior_family = "Mean movement + RMSSD + entropy ACF1",
  movement_mean_sex = "Mean movement + sex  (sensitivity)",
  primary_behavior_family_sex = "Mean movement + RMSSD + entropy ACF1 + sex  (sensitivity)")
absent <- setdiff(ladder$model_id, names(MODEL_LABEL))
if (length(absent))
  stop("The ladder contains models this figure has no label for: ",
       paste(absent, collapse = ", "),
       "\nA new model may not be rendered under a guessed name.", call. = FALSE)

ROLE_FILL <- c(reference_baseline = "grey72",
               primary_behavior_only = GREEN_DARK,
               sex_adjusted_sensitivity = GREEN_MID)
absent <- setdiff(ladder$primary_or_exploratory, names(ROLE_FILL))
if (length(absent))
  stop("Unrecognised model role(s): ", paste(absent, collapse = ", "),
       call. = FALSE)

# Registry order as stored, reversed so row 1 of the file prints at the bottom.
la <- ladder
la$label <- factor(MODEL_LABEL[la$model_id],
                   levels = rev(MODEL_LABEL[la$model_id]))
la$role <- factor(la$primary_or_exploratory,
                  levels = c("reference_baseline", "primary_behavior_only",
                             "sex_adjusted_sensitivity"),
                  labels = c("Reference baseline", "Behaviour only",
                             "Sex-adjusted sensitivity"))

# The permutation was run for the behaviour-only models only. Saying so is the
# point of the column: an absent P is a design fact, not a missing value.
la$p_text <- ifelse(
  is.na(la$permutation_p), "not run",
  ifelse(la$permutation_p < 0.001, "P < 0.001",
         paste0("P = ", fmt(la$permutation_p, 3))))

x_lo <- min(la$cv_r2_q025, la$loao_r2)
x_hi <- max(la$cv_r2_q975, la$loao_r2)
pad <- 0.045 * (x_hi - x_lo)
x_text <- x_hi + pad

pa <- ggplot(la, aes(y = .data$label)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = RULE) +
  geom_segment(aes(x = .data$cv_r2_q025, xend = .data$cv_r2_q975,
                   yend = .data$label),
               linewidth = 1.6, colour = GREEN_LIGHT, lineend = "round") +
  geom_point(aes(x = .data$repeated_cv_mean_r2), shape = 124, size = 1.5,
             colour = MUTED) +
  geom_point(aes(x = .data$loao_r2, fill = .data$role), shape = 21, size = 1.9,
             stroke = 0.3, colour = "grey20") +
  geom_text(aes(x = x_text, label = .data$p_text), hjust = 0, size = NOTE_PT / .pt,
            colour = MUTED) +
  scale_fill_manual(values = unname(ROLE_FILL[c("reference_baseline",
                                                "primary_behavior_only",
                                                "sex_adjusted_sensitivity")]),
                    drop = FALSE, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.22))) +
  labs(x = expression("Cross-validated " * italic(R)^2),
       y = NULL,
       title = "The complete a-priori model ladder",
       subtitle = paste0(
         "Point, leave-one-animal-out estimate; bar, 2.5-97.5% across repeated grouped five-fold\n",
         "splits, with its mean as a tick. n = ", N_ANIMALS, " animals. Rows are in registry order, not\n",
         "performance order, and no test compares two models.")) +
  theme_ed() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# ================================================================== panel b
#
# The formal test. Zero sits inside every interval, and the BH q values are
# printed rather than starred, because the result the manuscript relies on is
# the absence of an interaction and a reader is entitled to see how far from
# supported each one is.

FEATURE_LABEL <- c(Movement_mean = "Mean movement",
                   Movement_rmssd = "Movement RMSSD",
                   Entropy_acf1 = "Entropy ACF1")
absent <- setdiff(sexc$feature, names(FEATURE_LABEL))
if (length(absent))
  stop("The sex contract contains features this figure has no label for: ",
       paste(absent, collapse = ", "), call. = FALSE)

sx <- sexc
sx$label <- factor(FEATURE_LABEL[sx$feature],
                   levels = rev(FEATURE_LABEL[sx$feature]))
sx$q_text <- paste0("q = ", fmt(sx$interaction_q_bh, 2))

b_lo <- min(sx$interaction_ci_low)
b_hi <- max(sx$interaction_ci_high)
b_pad <- 0.05 * (b_hi - b_lo)

pb <- ggplot(sx, aes(y = .data$label)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = RULE) +
  geom_segment(aes(x = .data$interaction_ci_low, xend = .data$interaction_ci_high,
                   yend = .data$label),
               linewidth = 0.7, colour = MUTED, lineend = "round") +
  geom_point(aes(x = .data$interaction_estimate), shape = 21, size = 1.9,
             stroke = 0.3, fill = "white", colour = "grey20") +
  geom_text(aes(x = b_hi + b_pad, label = .data$q_text), hjust = 0,
            size = NOTE_PT / .pt, colour = MUTED) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.34))) +
  labs(x = "Feature x sex interaction estimate (95% CI)",
       y = NULL,
       title = "The association did not differ detectably by sex",
       subtitle = paste0(
         "One formal interaction test per early feature,\n",
         "Benjamini-Hochberg across the three. Zero lies\n",
         "inside every interval.")) +
  theme_ed() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# ================================================================== panel c
#
# Descriptive only, and labelled as such on the panel itself rather than only in
# the legend. These two numbers are the ones most likely to be misread as a
# sex difference, so they are drawn on a shared axis with a connecting rule that
# makes their closeness the visible fact.

rho <- data.frame(
  label = rep(FEATURE_LABEL[sexc$feature], 2L),
  sex = rep(c("Female", "Male"), each = nrow(sexc)),
  rho = c(sexc$female_rho, sexc$male_rho),
  stringsAsFactors = FALSE)
rho$label <- factor(rho$label, levels = rev(FEATURE_LABEL[sexc$feature]))

link <- data.frame(
  label = factor(FEATURE_LABEL[sexc$feature],
                 levels = rev(FEATURE_LABEL[sexc$feature])),
  lo = pmin(sexc$female_rho, sexc$male_rho),
  hi = pmax(sexc$female_rho, sexc$male_rho),
  stringsAsFactors = FALSE)

SEX_SHAPE <- c(Female = 21L, Male = 24L)

pc <- ggplot(rho, aes(y = .data$label)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = RULE) +
  geom_segment(data = link, aes(x = .data$lo, xend = .data$hi, yend = .data$label),
               linewidth = 0.5, colour = RULE) +
  geom_point(aes(x = .data$rho, shape = .data$sex, fill = .data$sex), size = 1.9,
             stroke = 0.3, colour = "grey20") +
  scale_shape_manual(values = SEX_SHAPE) +
  scale_fill_manual(values = c(Female = GREEN_DARK, Male = "white")) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
  labs(x = "Spearman rho, early feature versus later CombZ",
       y = NULL,
       title = "Sex-stratified estimates, descriptive",
       subtitle = paste0(
         "Shown because panel b is a negative result. These\n",
         "are not tests, carry no correction, and may not be\n",
         "reported as a sex difference.")) +
  theme_ed() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# ---------------------------------------------------------------- write out
# Authored at the exact placed size, so scale = 1.0 and native pt == printed pt.
# Two equal rows of 52 mm place 48 mm image boxes; a spans both columns.
DIMS <- list(a = c(173, 48), b = c(84.5, 48), c = c(84.5, 48))
PLOTS <- list(a = pa, b = pb, c = pc)

save_svg <- function(p, id) {
  d <- DIMS[[id]]
  f <- file.path(OUT, paste0("extended_data_behaviour_", id, ".svg"))
  svglite::svglite(f, width = d[1] / 25.4, height = d[2] / 25.4, bg = "white",
                   fix_text_size = FALSE)
  print(p); grDevices::dev.off()
  f
}
written <- vapply(names(PLOTS), function(id) save_svg(PLOTS[[id]], id),
                  character(1))

cat("Behavioural Extended Data panels ->", OUT, "\n")
for (id in names(written))
  cat(sprintf("  extended_data_behaviour_%s.svg  %6d bytes  %.1f x %.0f mm\n",
              id, file.size(written[[id]]), DIMS[[id]][1], DIMS[[id]][2]))
cat("  models:", nrow(ladder), "| interaction tests:", nrow(sexc),
    "| animals:", N_ANIMALS, "\n")
cat("  every value read from the frozen bridge contracts; nothing computed\n")
