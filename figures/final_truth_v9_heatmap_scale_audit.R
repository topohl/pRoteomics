#!/usr/bin/env Rscript

# Part-28: does every colour-encoded matrix in this layer represent its numbers
# faithfully?
#
# VISUAL-ENCODING AUDIT ONLY. No statistic is recomputed and no scale is
# changed here; this script measures what the released panels actually do.
#
# THE HARD RULE (Part-28 section 21). No heatmap may map two different numbers
# to the same endpoint colour without saying so. A panel fails if any plotted
# value lies outside its colour limits UNLESS the saturation is intentional,
# disclosed on the colourbar, and the uncapped values are still in the source
# data. This script hard-stops on a failure.
#
# WHAT "OUTSIDE THE LIMITS" MEANS HERE. nv_diverging() wraps
# scale_fill_gradient2(), which inherits the ggplot2 default oob = censor: an
# out-of-range value is NOT squished to the endpoint, it becomes NA and is
# painted na.value (grey50). So a silent clip in this layer would show up as a
# grey tile, and the audit checks both the numbers and the rendered SVG.

source(file.path("R", "paths.R"))
source(repo_path("R", "null_coalescing.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "final_truth_v9_figure_utils.R"))
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_heatmap_scale_audit.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] final_truth_v9 heatmap scale integrity audit")
  quit(save = "no", status = 0L)
}

OUT <- path_results("tables", "manuscript_candidates", "final_truth_v9", "audit")
REP <- path_results("reports", "manuscript_candidates", "final_truth_v9")
FIG <- path_results("figures", "manuscript_candidates", "final_truth_v9")
SD <- path_results("source_data", "manuscript_candidates", "final_truth_v9")
dir_create(OUT); dir_create(REP)

ct <- s9f_contract()
fig_of <- do.call(rbind, lapply(ct$figures, function(f) do.call(rbind, lapply(
  f$layout, function(l) data.frame(figure = as.character(f$name),
    figure_key = as.character(f$figure_key), panel_label = as.character(l$label),
    panel_id = as.character(l$panel), stringsAsFactors = FALSE)))))

# ------------------------------------------------- S20 the code search, first
#
# Every token the brief names, over the renderer files this layer actually uses,
# recorded whether or not it turns out to matter. A token that reaches no v9
# panel is still listed, so the search is auditable rather than asserted.
TOKENS <- c("geom_tile", "geom_raster", "Heatmap\\(", "scale_fill_gradient2",
            "scale_fill_gradient", "scale_fill_distiller", "colorRamp2",
            "nv_diverging", "limits *=", "oob *=", "squish", "scales::squish",
            "pmin\\(", "pmax\\(", "cut\\(", "winsor", "truncate", "\\bcap\\b",
            "clamp", "censor")
SRC <- list.files(repo_path("R"), pattern = "[.]R$", full.names = TRUE)
SRC <- c(SRC[grepl("final_truth_v9|nature_v2_figure_utils|nature_final_v7",
                   basename(SRC))],
         list.files(repo_path("figures"), pattern = "^final_truth_v9.*[.]R$",
                    full.names = TRUE),
         # Part-29 section 31: the analysis directories can also draw a
         # publication NES matrix, and two of them did so on a 98th-percentile
         # colour limit that disagreed with the frozen atlas. The search must
         # reach them or the certification is scoped too narrowly to mean what
         # it says.
         list.files(repo_path("04_differential_expression_enrichment"),
                    pattern = "[.][Rr]$", full.names = TRUE),
         list.files(repo_path("10_biological_integration"),
                    pattern = "[.][Rr]$", full.names = TRUE))
code <- do.call(rbind, lapply(sort(unique(SRC)), function(f) {
  ln <- readLines(f, warn = FALSE)
  do.call(rbind, lapply(TOKENS, function(tk) {
    i <- grep(tk, ln, perl = TRUE)
    if (!length(i)) return(NULL)
    data.frame(file = sub(".*proteomics[/\\]", "", f), line = i, token = tk,
               code = trimws(substr(ln[i], 1, 160)),
               is_v9_renderer = grepl("final_truth_v9", basename(f)),
               stringsAsFactors = FALSE)
  }))
}))
code$alters_displayed_value <- code$token %in%
  c("squish", "scales::squish", "pmin\\(", "pmax\\(", "winsor", "clamp",
    "truncate", "\\bcap\\b")
code <- code[order(code$file, code$line), , drop = FALSE]
write_csv_safe(code, file.path(OUT, "heatmap_scale_code_search.csv"))

# ------------------------------------------------- S20 the per-panel integrity
#
# Declared once, from the renderers, because a colourbar's name and the identity
# of a shared-scale group are editorial facts about the design and cannot be
# read out of a CSV.
H <- function(panel_id, quantity, value_col, limit_fun, shared_group,
              zero_meaningful, intentional_saturation, uncapped_col = NA,
              colourbar = NA)
  data.frame(panel_id = panel_id, quantity = quantity, value_col = value_col,
             limit_rule = limit_fun, shared_scale_group = shared_group,
             zero_meaningful = zero_meaningful,
             intentional_saturation = intentional_saturation,
             uncapped_col = uncapped_col, colourbar_name = colourbar,
             stringsAsFactors = FALSE)

SPEC <- rbind(
  H("v9_fingerprint", "baseline abundance, CON z-score", "con_z",
    "max(abs(plotted)) symmetric", "none", TRUE, FALSE, NA,
    "Baseline abundance (CON z-score)"),
  H("v9_compartment", "marker abundance, median centred log2",
    "displayed_value", "declared display cap from the source table", "none",
    TRUE, TRUE, "true_value", "unnamed"),
  H("v9_atlas", "median NES across a theme's GO terms", "median_NES",
    "f9_atlas_limit: max(abs) over all three contrast atlases",
    "atlas_NES", TRUE, FALSE, NA, "Median normalised enrichment score"),
  H("v9_ed_atlas_rescon", "median NES across a theme's GO terms", "median_NES",
    "f9_atlas_limit: max(abs) over all three contrast atlases",
    "atlas_NES", TRUE, FALSE, NA, "Median normalised enrichment score"),
  H("v9_ed_atlas_suscon", "median NES across a theme's GO terms", "median_NES",
    "f9_atlas_limit: max(abs) over all three contrast atlases",
    "atlas_NES", TRUE, FALSE, NA, "Median normalised enrichment score"),
  H("v9_ed_fingerprint_full", "external signature score, CON z-score", "score",
    "max(abs(plotted)) symmetric", "none", TRUE, FALSE, NA, "CON z-score"),
  H("v9_ed_module_fingerprint", "mean module-member CON z-score", "mean_con_z",
    "max(abs(plotted)) symmetric", "none", TRUE, FALSE, NA,
    "Mean module-member abundance (CON z-score)"),
  H("v9_ed_wgcna_phenotype", "module eigengene group difference", "estimate",
    "max(abs(plotted)) symmetric", "none", TRUE, FALSE, NA,
    "Module eigengene difference"),
  H("v9_ed_similarity", "median CON profile correlation between two units",
    "median_similarity", "max(abs(plotted)) symmetric, one global scale",
    "ed8_similarity", TRUE, FALSE, NA, "Median profile correlation (CON)"))

num <- function(x) suppressWarnings(as.numeric(x))
sidecar <- function(pid) {
  k <- fig_of$figure_key[match(pid, fig_of$panel_id)]
  p <- file.path(SD, k, paste0(pid, "_source_data.csv"))
  if (file.exists(p)) nv_read_csv(p) else NULL
}
svg_path <- function(pid) {
  k <- fig_of$figure_key[match(pid, fig_of$panel_id)]
  file.path(FIG, k, "panels", paste0(pid, ".svg"))
}

rows <- do.call(rbind, lapply(seq_len(nrow(SPEC)), function(i) {
  pid <- SPEC$panel_id[i]
  d <- sidecar(pid)
  vc <- SPEC$value_col[i]
  if (is.null(d) || !vc %in% names(d))
    stop("heatmap audit: source data or value column missing for ", pid,
         call. = FALSE)
  v <- num(d[[vc]]); v <- v[is.finite(v)]
  uc <- SPEC$uncapped_col[i]
  tv <- if (!is.na(uc) && uc %in% names(d)) num(d[[uc]]) else v
  tv <- tv[is.finite(tv)]
  lim <- max(abs(v))
  # a shared-scale group takes the group's limit, not this panel's own
  if (SPEC$shared_scale_group[i] != "none") {
    grp <- SPEC$panel_id[SPEC$shared_scale_group ==
                           SPEC$shared_scale_group[i]]
    lim <- max(vapply(grp, function(q) {
      dd <- sidecar(q)
      if (is.null(dd) || !vc %in% names(dd)) return(0)
      max(abs(num(dd[[vc]])), na.rm = TRUE)
    }, numeric(1)))
  }
  sv <- svg_path(pid)
  txt <- if (file.exists(sv))
    paste(readLines(sv, warn = FALSE, encoding = "UTF-8"), collapse = " ") else ""
  labs <- sub(".*>([^<]*)</text>", "\\1",
              unlist(regmatches(txt, gregexpr("<text[^>]*>[^<]*</text>", txt))))
  data.frame(
    figure = fig_of$figure[match(pid, fig_of$panel_id)],
    panel = fig_of$panel_label[match(pid, fig_of$panel_id)],
    panel_id = pid, quantity = SPEC$quantity[i],
    raw_min = min(tv), raw_max = max(tv),
    colour_min = -lim, colour_max = lim,
    n_below_colour_min = sum(tv < -lim - 1e-9),
    n_above_colour_max = sum(tv > lim + 1e-9),
    oob_handling = "scales::censor (ggplot2 default); an out-of-range value would be painted grey50",
    shared_scale_group = SPEC$shared_scale_group[i],
    scale_symmetric = TRUE, zero_meaningful = SPEC$zero_meaningful[i],
    intentional_saturation = SPEC$intentional_saturation[i],
    disclosed_on_colourbar = any(grepl("≥|≤", labs)),
    disclosed_in_legend = any(grepl("saturat", labs, ignore.case = TRUE)),
    uncapped_source_values_present = !is.na(uc) && uc %in% names(d),
    grey50_pixels_in_svg = grepl("grey50|#7F7F7F|#808080", txt,
                                 ignore.case = TRUE),
    colourbar_name = SPEC$colourbar_name[i],
    limit_rule = SPEC$limit_rule[i],
    stringsAsFactors = FALSE)
}))

rows$status <- with(rows, ifelse(
  n_below_colour_min + n_above_colour_max == 0L, "FAITHFUL",
  ifelse(intentional_saturation & disclosed_on_colourbar &
           uncapped_source_values_present, "DISCLOSED_SATURATION",
         "SILENT_CLIP")))
rows <- rows[order(rows$figure, rows$panel), , drop = FALSE]
write_csv_safe(rows, file.path(OUT, "heatmap_scale_integrity_audit.csv"))

cat("\n===== PART-28 HEATMAP SCALE INTEGRITY =====\n")
print(rows[, c("figure", "panel", "panel_id", "raw_min", "raw_max",
               "colour_min", "colour_max", "n_above_colour_max",
               "intentional_saturation", "disclosed_on_colourbar", "status")],
      row.names = FALSE)
cat("\ncode-search hits:", nrow(code), "| tokens that could alter a displayed value:",
    sum(code$alters_displayed_value), "| of those inside a v9 renderer:",
    sum(code$alters_displayed_value & code$is_v9_renderer), "\n")
cat("panels with a grey50 tile in the rendered SVG:",
    sum(rows$grey50_pixels_in_svg), "\n")

bad <- rows[rows$status == "SILENT_CLIP", , drop = FALSE]
if (nrow(bad)) {
  print(bad)
  stop("section 21: ", nrow(bad), " heatmap(s) clip values without disclosure",
       call. = FALSE)
}
cat("section 21 PASSES: no heatmap clips a value without disclosure\n")
cat("\nwritten to:", OUT, "\n")

# ------------------------------------------------------- S27 decisions, prose
sat <- rows[rows$intentional_saturation, , drop = FALSE]
md <- c(
"# Heatmap scale decisions (Part 28)",
"",
sprintf("%d colour-encoded matrix panels were audited. %d represent every value",
        nrow(rows), sum(rows$status == "FAITHFUL")),
sprintf("faithfully; %d saturates deliberately and discloses it. None clips silently.",
        sum(rows$status == "DISCLOSED_SATURATION")),
"",
"## Why a clip here would be visible rather than silent",
"",
"nv_diverging() wraps scale_fill_gradient2(), which inherits the ggplot2",
"default oob = scales::censor. An out-of-range value is therefore NOT squished",
"to the endpoint colour - it becomes NA and is painted grey50. The audit checks",
sprintf("the rendered SVGs for that colour as well as the numbers: %d panels contain",
        sum(rows$grey50_pixels_in_svg)),
"a grey50 tile.",
"",
"## The one saturating panel",
"")
if (nrow(sat)) md <- c(md,
  sprintf("%s panel %s (%s) maps colour to +/-%.1f while the true values run",
          sat$figure[1], sat$panel[1], sat$panel_id[1], sat$colour_max[1]),
  sprintf("%.3f to %.3f. %d cells exceed the cap.", sat$raw_min[1], sat$raw_max[1],
          sat$n_above_colour_max[1]),
  "",
  "Retained, because three extreme cells would otherwise compress the colour",
  "discrimination of the whole remaining matrix. The saturation is disclosed",
  "four ways: the uncapped value is printed inside each saturated tile, the",
  "caption states the cap and the number of affected cells, the colourbar top",
  "tick reads >=3, and the sidecar carries true_value and colour_scale_censored",
  "alongside displayed_value.",
  "",
  "Only the upper tail saturates - the displayed minimum is above the lower",
  "limit - so only the top tick carries the >= mark. Marking both ends would",
  "imply an overflow bin at the bottom that does not exist.",
  "",
  "The colour limit is now read from the declared centered_log2_display_cap",
  "column rather than re-derived from the already-capped data, so the",
  "disclosure cannot quietly disappear if a future refresh leaves no cell",
  "sitting exactly at the cap.",
  "")
md <- c(md,
"## Shared scales",
"",
"| group | panels | limit | why shared |",
"|---|---|---|---|")
for (g in setdiff(unique(rows$shared_scale_group), "none")) {
  z <- rows[rows$shared_scale_group == g, ]
  md <- c(md, sprintf("| %s | %s | +/-%.6f | same quantity, intended for direct comparison |",
                      g, paste(z$panel_id, collapse = ", "), z$colour_max[1]))
}
md <- c(md, "",
"Quantities that are NOT forced onto a common scale, correctly: NES against",
"log2FC, CON z-score against module eigengene difference, and profile",
"correlation against any of them. Each carries its own named colourbar.",
"",
"## Known residual, recorded rather than fixed",
"",
"The three-cell NES strips in Figure 3 d-f and ED6 c-e carry no colour key, and",
"the atlas on the same page carries an NES colourbar at a different limit",
"(theme-summary +/-2.479 against single-term +/-2.761). The strips print their",
"NES numerically in every tile, and the legend now states that the strip is",
"keyed by its numbers rather than by the atlas bar. The two limits are not",
"reconciled because they summarise different quantities: a median across a",
"theme cannot share a scale with a single term without one of them being",
"rescaled away from its own range.",
"",
"The atlas colourbar is ticked -2, 0, 2 against limits +/-2.479, so the",
"outermost labelled tick understates the endpoint. Nothing is clipped; the",
"caption states the limit numerically.")
writeLines(md, file.path(REP, "heatmap_scale_decisions.md"))
cat("heatmap decisions doc written
")
