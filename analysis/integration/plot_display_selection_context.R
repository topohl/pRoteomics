#!/usr/bin/env Rscript
# ================================================================
# Script: analysis/integration/plot_display_selection_context.R
# Stage: integration
# Scope: global
# Consumes: required results/integration/build_display_selection_inventories/global/tables/pathway_enrichment_inventory.csv; results/integration/build_display_selection_inventories/global/tables/leading_edge_protein_recurrence.csv; results/integration/build_display_selection_inventories/global/tables/leading_edge_protein_inventory.csv
# Produces: results/integration/plot_display_selection_context/global/plots/display_selection_context.svg; results/integration/plot_display_selection_context/global/plots/display_selection_context.pdf; results/integration/plot_display_selection_context/global/tables/source_data/display_selection_context_panel_a.csv; results/integration/plot_display_selection_context/global/tables/source_data/display_selection_context_panel_b.csv; results/integration/plot_display_selection_context/global/tables/source_data/display_selection_context_panel_c.csv
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Diagnostic figure showing each displayed selection inside its complete distribution. Not a manuscript panel.
# ================================================================
#
# WHAT THIS IS
#
# The visual form of docs/FIGURE_SELECTION_RULES.md: it shows each selection
# that Figures 2 and 3 display sitting inside the complete set it was drawn
# from. A reader will not open a 203,073-row CSV; they will look at this.
#
# THIS IS A DIAGNOSTIC PLOT, NOT A MANUSCRIPT PANEL. It is deliberately built
# in this repository, under plots/, and not in the manuscript figure layer.
# Promoting it to an Extended Data figure is a separate, deliberate decision
# that would touch the frozen figure contract, the legend layer, the figure
# index and the promotion tests. Nothing here presumes that decision.
#
# No statistic is computed. Every value is read from the released inventories.
#
# -------------------------------------------------------------- palette
#
# Colours come from config/manuscript_palette.yml, which is a frozen contract,
# so the encodings match the paper's other figures. The palette was checked
# against the computable colour rules rather than judged by eye, and the two
# results below changed the design:
#
#   1. The compartment triple (#3D5A73 neuropil, #7A9BB0 soma, #C2A878
#      microglia) passes CVD separation (worst min(protan,deutan) dE = 11.3)
#      but the soma/microglia pair sits at dE = 13.6 for NORMAL vision, under
#      the floor of 15. Full-colour readers would find that pair harder to
#      separate than a colour-blind reader finds the others. The contract is
#      frozen and must not be edited, so compartment is carried as a LABELLED
#      FACET instead of a colour series, and the failing pair never has to be
#      told apart by colour at all.
#
#   2. The evidence pair (#1F3D52 supported, #B9B9B4 descriptive) separates
#      enormously (dE = 44) but both slots fall outside the categorical
#      lightness band. That check assumes series of equal standing; this is a
#      figure/ground pair in which the grey is MEANT to recede, which is the
#      palette contract's own stated rule. It is used as intended, and carries
#      size and direct labels as well as colour, so identity is never
#      colour-alone.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

ANALYSIS_ID <- "plot_display_selection_context"
SCRIPT_ID <- "analysis/integration/plot_display_selection_context.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

run <- integration_cli(default_dataset = "all", allow_all = TRUE)
paths <- integration_dirs(ANALYSIS_ID, "global", create = !isTRUE(run$dry_run))

SRC <- function(f) integration_find(f, owner = "build_display_selection_inventories",
                                    legacy_stage = "10_biological_integration",
                                    legacy_substep = "display_selection_inventories")
inputs <- list(
  pathways = SRC("pathway_enrichment_inventory.csv"),
  recurrence = SRC("leading_edge_protein_recurrence.csv")
)

if (isTRUE(run$dry_run)) {
  dry_run_inputs(SCRIPT_ID, inputs)
  quit(status = 0, save = "no")
}

rd <- function(p) {
  if (!file.exists(p)) {
    stop("missing_required_input: run analysis/integration/",
         "build_display_selection_inventories.R first: ", p, call. = FALSE)
  }
  utils::read.csv(p, stringsAsFactors = FALSE)
}
inv <- rd(inputs$pathways)
rec <- rd(inputs$recurrence)
## leading_edge_protein_inventory.csv is deliberately NOT read here. It was
## loaded into `le` and never used: 46 MB through utils::read.csv for nothing.
## The recurrence table already carries the per-protein classification these
## panels draw, so the dependency was declared but not real, and it has been
## removed from this script's consumes_required in pipeline.yml to match.

CONTRAST <- "SUS - RES"

# ------------------------------------------------------------------ palette
pal <- yaml::read_yaml(repo_path("config", "manuscript_palette.yml"))
SUPPORTED <- pal$evidence$supported        # #1F3D52, the displayed selection
BACKGROUND <- pal$evidence$descriptive     # #B9B9B4, everything it was drawn from
GRID <- "grey92"
INK <- "grey20"; INK2 <- "grey40"
FAM <- pal$typography$family

# Diagnostic sizing: readable on a slide and in a rebuttal letter, rather than
# the 5 pt manuscript floor, which exists for a 183 mm print canvas.
theme_diag <- function(grid = c("x", "y", "none")) {
  grid <- match.arg(grid)
  theme_minimal(base_family = FAM, base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = if (grid == "x") element_line(colour = GRID, linewidth = 0.25) else element_blank(),
      panel.grid.major.y = if (grid == "y") element_line(colour = GRID, linewidth = 0.25) else element_blank(),
      axis.text = element_text(size = 7, colour = INK2),
      axis.title = element_text(size = 7.5, colour = INK),
      plot.title = element_text(size = 9, face = "bold", colour = INK),
      plot.subtitle = element_text(size = 7.5, colour = INK2),
      strip.text = element_text(size = 7.5, face = "bold", colour = INK),
      legend.position = "none",
      panel.spacing = unit(7, "pt"),
      plot.margin = margin(3, 6, 3, 4))
}

EXEMPLARS <- data.frame(
  dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
  spatial_unit = c("CA3_sr", "CA2_sp", "CA1"),
  GO_ID = c("GO:0099536", "GO:0006397", "GO:0006119"),
  program = c("synaptic signalling", "mRNA processing", "oxidative phosphorylation"),
  stringsAsFactors = FALSE)
# Compartment is a labelled facet, not a colour series: see the palette note.
FACET <- c(neuron_neuropil = "Neuropil · CA3-SR",
           neuron_soma = "Neuronal soma · CA2-SP",
           microglia = "Microglia-enriched · CA1")

# ============================================================ PANEL A
# Every FDR-supported claim-eligible term in the exemplar's own spatial unit,
# with the exemplar marked. This is the "is it cherry-picked?" panel, and it
# answers the question against the exemplar in one case out of three.
a_rows <- do.call(rbind, lapply(seq_len(nrow(EXEMPLARS)), function(i) {
  z <- inv[inv$contrast == CONTRAST & inv$dataset == EXEMPLARS$dataset[i] &
             inv$spatial_unit == EXEMPLARS$spatial_unit[i] &
             inv$fdr_supported & inv$theme_claim_eligible, , drop = FALSE]
  z <- z[order(-abs(z$NES)), , drop = FALSE]
  z$rank_by_abs_nes <- seq_len(nrow(z))
  z$facet <- FACET[[EXEMPLARS$dataset[i]]]
  z$is_exemplar <- z$GO_ID == EXEMPLARS$GO_ID[i]
  z$n_in_unit <- nrow(z)
  z[, c("facet", "dataset", "spatial_unit", "GO_ID", "GO_description", "NES",
        "BH_FDR", "rank_by_abs_nes", "n_in_unit", "is_exemplar")]
}))
a_rows$facet <- factor(a_rows$facet, levels = unname(FACET))
ex_a <- a_rows[a_rows$is_exemplar, , drop = FALSE]
# The curve descends left to right, so the top-right of each facet is the one
# reliably empty region. Labelling there instead of beside the point stops the
# text from being drawn across the very distribution it is describing, and a
# leader line keeps the association explicit.
ex_a$label <- sprintf("%s\nrank %d of %d", ex_a$GO_description, ex_a$rank_by_abs_nes,
                      ex_a$n_in_unit)
ex_a$lab_x <- ex_a$n_in_unit * 0.55
ex_a$lab_y <- vapply(seq_len(nrow(ex_a)), function(i)
  max(abs(a_rows$NES[a_rows$facet == ex_a$facet[i]])), numeric(1))

pa <- ggplot(a_rows, aes(rank_by_abs_nes, abs(NES))) +
  geom_point(colour = BACKGROUND, size = 1.1, stroke = 0) +
  geom_segment(data = ex_a, aes(x = lab_x, y = lab_y, xend = rank_by_abs_nes,
                                yend = abs(NES)),
               colour = SUPPORTED, linewidth = 0.2, linetype = "22") +
  geom_point(data = ex_a, colour = SUPPORTED, size = 2.6, stroke = 0) +
  geom_text(data = ex_a, aes(lab_x, lab_y, label = label), hjust = 0.5,
            vjust = -0.15, family = FAM, size = 2.35, colour = SUPPORTED,
            lineheight = 0.95) +
  facet_wrap(~facet, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.14))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.30))) +
  labs(title = "a  The three exemplar programs inside their own spatial unit",
       subtitle = paste0("Every FDR-supported claim-eligible GO term in that unit for ",
                         CONTRAST, ". Dark point = the term drawn in Figure 3 d/e/f."),
       x = "Rank by |NES| within the spatial unit", y = "|NES|") +
  theme_diag("y")

# ============================================================ PANEL B
# The 21 displayed proteins against the recurrence rule they were NOT selected
# by. They were chosen by rank statistic; recurrence is an independent count.
SHOWN <- c("App", "Cnr1", "Dbi", "Eif4ebp2", "Ly6h", "Plppr4", "Synpo",
           "Cirbp", "Csdc2", "Dcps", "Ddx23", "Lsm3", "Lsm8", "Rbm8a",
           "Cox5b", "Cox6b1", "Iscu", "Ndufs8", "Ndufv2", "Ndufv3", "Uqcrh")
r <- rec[rec$contrast == CONTRAST, , drop = FALSE]
CLASSES <- c("SINGLE_APPEARANCE", "INTERMEDIATE", "RECURRENT_CORE")
CLASS_LAB <- c(SINGLE_APPEARANCE = "Single appearance",
               INTERMEDIATE = "Intermediate",
               RECURRENT_CORE = "Recurrent core\n(≥3 terms, ≥3 contexts)")
grp <- function(z, nm) {
  z <- z[z$recurrence_class %in% CLASSES, , drop = FALSE]
  tb <- table(factor(z$recurrence_class, levels = CLASSES))
  data.frame(group = nm, recurrence_class = names(tb), n = as.integer(tb),
             pct = 100 * as.integer(tb) / sum(tb), total = sum(tb),
             stringsAsFactors = FALSE)
}
# Both groups are restricted to proteins that carry a supported term inside a
# claim-eligible theme, which is the only comparison that is like for like:
# the 21 displayed proteins are leading edge of claim-eligible exemplar terms
# by construction, so the wider pool of leading-edge proteins whose terms all
# fall outside the seven atlas rows is not a population they could have been
# drawn from. Naming the restriction in the facet keeps the denominator visible.
b_rows <- rbind(
  grp(r, sprintf("All in a claim-eligible theme (n = %d)",
                 sum(r$recurrence_class %in% CLASSES))),
  grp(r[r$gene %in% SHOWN, , drop = FALSE], "Shown in Figure 3 g/h/i (n = 21)"))
b_rows$recurrence_class <- factor(b_rows$recurrence_class, levels = CLASSES)
b_rows$group <- factor(b_rows$group, levels = unique(b_rows$group))
b_rows$is_shown <- grepl("^Shown", b_rows$group)

pb <- ggplot(b_rows, aes(pct, recurrence_class, fill = is_shown)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f%%  (%d)", pct, n)), hjust = -0.12,
            family = FAM, size = 2.4, colour = INK2) +
  facet_wrap(~group, nrow = 1) +
  scale_fill_manual(values = c(`FALSE` = BACKGROUND, `TRUE` = SUPPORTED)) +
  scale_y_discrete(labels = CLASS_LAB) +
  scale_x_continuous(limits = c(0, 118), breaks = c(0, 50, 100),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "b  The displayed proteins against a rule they were not selected by",
       subtitle = paste0("Selection used the GSEA rank statistic. Recurrence is an ",
                         "independent count of enrichment memberships, ", CONTRAST, "."),
       x = "% of proteins in that group", y = NULL) +
  theme_diag("x")

# ============================================================ PANEL C
# The denominator, made concrete.
z <- inv[inv$contrast == CONTRAST, , drop = FALSE]
c_rows <- data.frame(
  step = c("GO terms tested", "In a claim-eligible atlas theme",
           "and FDR-supported", "Atlas rows drawn in Figure 3b"),
  n = c(nrow(z), sum(z$theme_claim_eligible),
        sum(z$theme_claim_eligible & z$fdr_supported), 7L),
  stringsAsFactors = FALSE)
c_rows$step <- factor(c_rows$step, levels = rev(c_rows$step))
c_rows$is_shown <- c_rows$n == 7L

# A bar encodes magnitude by LENGTH, and on a log axis length is not
# proportional to value: a bar chart here would make 67,691 look about four
# times 7. The counts span four orders of magnitude, so the axis has to be log
# and the mark therefore has to be a position mark, not a length one.
pc <- ggplot(c_rows, aes(n, step, colour = is_shown)) +
  geom_segment(aes(x = 1, xend = n, yend = step), linewidth = 0.3,
               colour = GRID) +
  geom_point(size = 2.6) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.35,
            family = FAM, size = 2.4, colour = INK2) +
  scale_colour_manual(values = c(`FALSE` = BACKGROUND, `TRUE` = SUPPORTED)) +
  scale_x_log10(limits = c(1, 10^6), breaks = c(1, 10, 100, 1000, 10000, 100000),
                labels = c("1", "10", "100", "1,000", "10,000", "100,000"),
                expand = expansion(mult = c(0, 0))) +
  labs(title = "c  What the seven atlas rows were selected from",
       subtitle = paste0(CONTRAST, ". Log scale. A cell carries a support dot when ",
                         "at least one constituent term is FDR-supported."),
       x = "Count (log scale)", y = NULL) +
  theme_diag("x")

# ============================================================ ASSEMBLE
fig <- pa / pb / pc + plot_layout(heights = c(1.15, 1, 0.78))

plots_dir <- paths$plots %||% paths$figures
dir_create(plots_dir)
svg_path <- file.path(plots_dir, "display_selection_context.svg")
pdf_path <- file.path(plots_dir, "display_selection_context.pdf")
W <- 183 / 25.4; H <- 175 / 25.4
ggsave(svg_path, fig, width = W, height = H, units = "in", device = "svg")
ggsave(pdf_path, fig, width = W, height = H, units = "in", device = grDevices::cairo_pdf)

write_csv_safe(a_rows, file.path(paths$source_data, "display_selection_context_panel_a.csv"))
write_csv_safe(b_rows, file.path(paths$source_data, "display_selection_context_panel_b.csv"))
write_csv_safe(c_rows, file.path(paths$source_data, "display_selection_context_panel_c.csv"))

write_integration_manifest(
  paths, inputs, list(plots = plots_dir, source_data = paths$source_data),
  list(dataset = run$dataset, contrast = CONTRAST),
  paste0("Diagnostic figure placing each displayed selection inside its complete ",
         "distribution. Reads the released inventories; computes no statistic. ",
         "Not a manuscript panel."))

message("panel a: ", nrow(a_rows), " terms across 3 units")
message("panel b: ", sum(b_rows$n[b_rows$is_shown]), " displayed proteins classified")
message("panel c: ", format(c_rows$n[c_rows$step == "GO terms tested"], big.mark = ","),
        " terms tested for ", CONTRAST)
message("figure: ", relative_to(svg_path))
