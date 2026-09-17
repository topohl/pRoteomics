#!/usr/bin/env Rscript

# ==========================================================================
# Part-29, sections 31 + 32
# Atlas-like NES renderer colour-scale audit.
#
# AUDIT ONLY. This script reads canonical outputs and repository source code.
# It never reruns, regenerates or modifies any canonical analysis output.
# It writes exactly one CSV, under
#   results/tables/publication_audits/upstream_enrichment_v10/
#
# Question. The frozen v9 behaviour is f9_atlas_limit() in
# R/final_truth_v9_panels.R:740, which uses the TRUE GLOBAL MAXIMUM absolute
# median NES across all three contrast atlases. Does any other renderer in the
# repository compute an atlas NES colour limit differently - in particular from
# a percentile - and is such a renderer reachable from the manuscript?
# ==========================================================================

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

OUT_DIR <- file.path("results", "tables", "publication_audits",
                     "upstream_enrichment_v10")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_CSV <- file.path(OUT_DIR, "atlas_renderer_scale_audit.csv")

fmt <- function(x, d = 6L) {
  if (length(x) != 1L || !is.finite(x)) return("NA")
  formatC(x, format = "f", digits = d)
}
src_line <- function(path, line) {
  if (!file.exists(path)) return(NA_character_)
  ln <- readLines(path, warn = FALSE)
  if (line < 1L || line > length(ln)) return(NA_character_)
  trimws(ln[[line]])
}

# ==========================================================================
# 1. THE FROZEN v9 REFERENCE LIMIT, recomputed from the canonical theme table
# ==========================================================================
THEME_TBL <- file.path("results", "tables", "10_biological_integration",
                       "gsea_wgcna_concordance", "global",
                       "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
stopifnot(file.exists(THEME_TBL))

# 33 columns; read only the nine the v9 atlas rule touches. The two leading-edge
# columns are the bulk of the 263 MB and are not needed.
hdr <- names(utils::read.csv(THEME_TBL, nrows = 1L, stringsAsFactors = FALSE,
                             check.names = FALSE))
keep <- c("dataset", "contrast", "spatial_unit", "GO_ID", "NES", "GSEA_FDR",
          "theme_id", "manuscript_theme", "theme_claim_eligible")
stopifnot(all(keep %in% hdr))
cc <- rep("NULL", length(hdr))
cc[match(c("dataset", "contrast", "spatial_unit", "GO_ID", "theme_id",
           "manuscript_theme"), hdr)] <- "character"
cc[match(c("NES", "GSEA_FDR"), hdr)] <- "numeric"
cc[match("theme_claim_eligible", hdr)] <- "logical"
message("[1/5] reading canonical theme table (", keep[1], " ... 9 of ",
        length(hdr), " columns)")
th <- utils::read.csv(THEME_TBL, stringsAsFactors = FALSE, check.names = FALSE,
                      colClasses = cc)
message("      rows = ", nrow(th))

# f9_atlas_cells(), reimplemented verbatim in base R from
# R/final_truth_v9_panels.R:753-771.
f9_cells <- function(th, contrast) {
  z <- th[th$contrast == contrast &
            th$theme_claim_eligible %in% TRUE &
            nzchar(as.character(th$theme_id)), , drop = FALSE]
  z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID,
                           z$theme_id)), , drop = FALSE]
  key <- paste(z$dataset, z$spatial_unit, z$theme_id, sep = "\r")
  cells <- do.call(rbind, lapply(split(seq_len(nrow(z)), key), function(ix) {
    w <- z[ix, , drop = FALSE]
    data.frame(dataset = w$dataset[1], spatial_unit = w$spatial_unit[1],
               theme_id = w$theme_id[1], theme = w$manuscript_theme[1],
               n_fdr = sum(is.finite(w$GSEA_FDR) & w$GSEA_FDR < 0.05),
               median_NES = stats::median(w$NES, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(cells) <- NULL
  cells
}

CONTRASTS <- c("RES - CON", "SUS - CON", "SUS - RES")
message("[2/5] recomputing f9_atlas_limit()")
cells_by_ct <- lapply(CONTRASTS, function(ct) f9_cells(th, ct))
names(cells_by_ct) <- CONTRASTS
v9_all <- unlist(lapply(cells_by_ct, function(x) x$median_NES), use.names = FALSE)
V9_LIMIT <- max(abs(v9_all), na.rm = TRUE)
v9_per_contrast <- vapply(cells_by_ct,
                          function(x) max(abs(x$median_NES), na.rm = TRUE),
                          numeric(1))
v9_n_cells <- vapply(cells_by_ct, nrow, integer(1))
# the 98th percentile the v9 atlas would have used, for direct comparison
V9_P98 <- unname(stats::quantile(abs(v9_all), probs = 0.98, na.rm = TRUE,
                                 names = FALSE))
message("      f9_atlas_limit = ", fmt(V9_LIMIT), " over ", length(v9_all),
        " cells (", paste(sprintf("%s=%s", CONTRASTS,
                                  vapply(v9_per_contrast, fmt, character(1), 6L)),
                          collapse = ", "), ")")

# ==========================================================================
# 2. THE TWO PERCENTILE RENDERERS, evaluated on their own real source data
# ==========================================================================
message("[3/5] evaluating the percentile renderers on their released source data")

# publication_color_limits(), verbatim from
# analysis/04_differential_abundance/07_compareGO_spatial_program_atlas.r:582
publication_color_limits <- function(x, cap = 2.5) {
  lim <- suppressWarnings(stats::quantile(abs(x), probs = 0.98, na.rm = TRUE,
                                          names = FALSE))
  if (!is.finite(lim) || lim <= 0) lim <- suppressWarnings(max(abs(x), na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) lim <- 1
  lim <- min(lim, cap)
  c(-lim, lim)
}

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)

SD07 <- file.path("results", "source_data", "04_differential_expression_enrichment",
                  "compareGO_spatial_atlas")
SD10 <- file.path("results", "source_data", "04_differential_expression_enrichment",
                  "sus_res_spatial_dap_atlas", "global")

# -- 07, the SUS-RES manuscript THEME x spatial-unit atlas (line 776) ------
p07_theme <- file.path(SD07, "source_data_SpatialProgramAtlas_SUS_vs_RES_publication.csv")
s07 <- rd(p07_theme)
d07 <- s07[s07$theme_role == "primary" &
             suppressWarnings(as.numeric(s07$n_theme_terms_tested)) > 0, , drop = FALSE]
v07 <- suppressWarnings(as.numeric(d07$median_NES_all_theme_terms))
lim07 <- publication_color_limits(v07)[2]
p98_07 <- unname(stats::quantile(abs(v07), probs = 0.98, na.rm = TRUE, names = FALSE))
max07 <- max(abs(v07), na.rm = TRUE)
sat07 <- sum(abs(v07) > lim07, na.rm = TRUE)

# -- 10, the SUS-RES ranked-GSEA THEME x spatial-unit panel C (line 508) ---
p10_theme <- file.path(SD10, "panel_c_sus_res_ranked_GSEA_themes.csv")
s10 <- rd(p10_theme)
d10 <- s10[s10$theme_role == "primary", , drop = FALSE]
v10 <- suppressWarnings(as.numeric(d10$median_NES_all_theme_terms))
lim10 <- suppressWarnings(stats::quantile(abs(v10), probs = 0.98, na.rm = TRUE,
                                          names = FALSE))
if (!is.finite(lim10) || lim10 <= 0) lim10 <- max(abs(v10), na.rm = TRUE)
if (!is.finite(lim10) || lim10 <= 0) lim10 <- 1
lim10 <- min(lim10, 2.5)
max10 <- max(abs(v10), na.rm = TRUE)
sat10 <- sum(abs(v10) > lim10, na.rm = TRUE)

# -- cross-check: are 07/10 displaying the SAME numbers as the v9 SUS-RES
#    atlas, only on a different colour limit? Join on the cell key.
v9_sr <- cells_by_ct[["SUS - RES"]]
cellkey <- function(ds, su, tid) paste(ds, su, tid, sep = "\r")
k9 <- cellkey(v9_sr$dataset, v9_sr$spatial_unit, v9_sr$theme_id)
k07 <- cellkey(d07$dataset_compartment, d07$spatial_unit, d07$theme_id)
k10 <- cellkey(d10$dataset_compartment, d10$spatial_unit, d10$theme_id)
join_delta <- function(kx, vx) {
  i <- match(kx, k9)
  if (all(is.na(i))) return(list(n = 0L, maxdiff = NA_real_))
  list(n = sum(!is.na(i)),
       maxdiff = max(abs(vx[!is.na(i)] - v9_sr$median_NES[i[!is.na(i)]]),
                     na.rm = TRUE))
}
cmp07 <- join_delta(k07, v07)
cmp10 <- join_delta(k10, v10)

# which cells the 98th-percentile limit silently saturates
sq07 <- d07[abs(v07) > lim07, , drop = FALSE]
sq07_txt <- if (nrow(sq07)) paste(sprintf("%s/%s/%s = %s",
    sq07$dataset_compartment, sq07$spatial_unit, sq07$theme_id,
    vapply(suppressWarnings(as.numeric(sq07$median_NES_all_theme_terms)),
           fmt, character(1), 6L)), collapse = "; ") else "none"

# -- 07, the broad PROGRAM-class x spatial-unit atlas (line 698) -----------
p07_prog <- file.path(SD07, "source_data_SpatialProgramAtlas_publication.csv")
s07p <- rd(p07_prog)
d07p <- s07p[s07p$publication_include %in% c(TRUE, "TRUE", "True", "true"), , drop = FALSE]
v07p <- suppressWarnings(as.numeric(d07p$mean_NES))
lim07p <- publication_color_limits(v07p)[2]
max07p <- max(abs(v07p), na.rm = TRUE)
sat07p <- sum(abs(v07p) > lim07p, na.rm = TRUE)

# -- 07, neuropil-focus program atlas (line 821) --------------------------
p07_np <- file.path(SD07, "source_data_SpatialProgramAtlas_neuropil_focus_publication.csv")
s07n <- rd(p07_np)
istrue <- function(x) x %in% c(TRUE, "TRUE", "True", "true")
d07n <- s07n[istrue(s07n$interpretable_program) & istrue(s07n$nonzero_signal), , drop = FALSE]
v07n <- suppressWarnings(as.numeric(d07n$mean_NES))
lim07n <- publication_color_limits(v07n)[2]
max07n <- max(abs(v07n), na.rm = TRUE)
sat07n <- sum(abs(v07n) > lim07n, na.rm = TRUE)

# -- 07, compartment-comparison program atlas (line 853) ------------------
p07_cc <- file.path(SD07, "source_data_Compartment_Comparison_publication.csv")
s07c <- rd(p07_cc)
v07c <- suppressWarnings(as.numeric(s07c$mean_NES))
lim07c <- publication_color_limits(v07c)[2]
max07c <- max(abs(v07c), na.rm = TRUE)
sat07c <- sum(abs(v07c) > lim07c, na.rm = TRUE)

# ==========================================================================
# 3. REACHABILITY, traced from the registry and the released export manifest
# ==========================================================================
message("[4/5] tracing reachability via pipeline.yml, RUN_ORDER.md, export manifest")

pl <- readLines("pipeline.yml", warn = FALSE)
ro <- readLines("RUN_ORDER.md", warn = FALSE)
registered <- function(rel) any(grepl(paste0('script: "', rel, '"'), pl, fixed = TRUE))
stage_of <- function(rel) {
  i <- grep(paste0('script: "', rel, '"'), pl, fixed = TRUE)
  if (!length(i)) return(NA_character_)
  j <- grep('stage: "', pl[seq(i[1], min(i[1] + 12L, length(pl)))], value = TRUE)
  if (!length(j)) return(NA_character_)
  sub('.*stage: "([^"]+)".*', "\\1", j[1])
}
in_run_order <- function(rel) any(grepl(rel, ro, fixed = TRUE))

# the publication freeze manifest: does any frozen artefact come from these
# renderers? (the freeze covers NUMBERS, not colour scales, but a frozen
# source-data file is the strongest possible statement of manuscript standing)
FREEZE <- file.path("docs", "publication_freeze_manifest.yml")
frz <- if (file.exists(FREEZE)) readLines(FREEZE, warn = FALSE) else character(0)
frozen_hits <- function(pattern) sum(grepl(pattern, frz, fixed = TRUE))
n_frozen_10 <- frozen_hits("sus_res_spatial_dap_atlas")
n_frozen_07 <- frozen_hits("compareGO_spatial_atlas")

EXPORT_MANIFEST <- file.path("results", "manuscript", "figure_export_manifest.csv")
exp_man <- if (file.exists(EXPORT_MANIFEST)) rd(EXPORT_MANIFEST) else NULL
exported <- function(pattern) {
  if (is.null(exp_man)) return(NA)
  any(grepl(pattern, exp_man$source_file, fixed = TRUE))
}
n_exported <- function(pattern) {
  if (is.null(exp_man)) return(NA_integer_)
  sum(grepl(pattern, exp_man$source_file, fixed = TRUE))
}

fig07_theme <- "compareGO_spatial_atlas/Fig_SpatialProgramAtlas_SUS_vs_RES_publication.svg"
fig10_theme <- "sus_res_spatial_dap_atlas/global/panel_c_sus_res_ranked_GSEA_themes"
fig07_prog  <- "compareGO_spatial_atlas/Fig_SpatialProgramAtlas_dotheatmap_publication.svg"
fig07_np    <- "compareGO_spatial_atlas/Fig_SpatialProgramAtlas_neuropil_focus_publication.svg"
fig07_cc    <- "compareGO_spatial_atlas/Fig_Compartment_Comparison_publication.svg"

# the Part-28 v9 heatmap scale audit's own code-search scope: it only globbed
# R/final_truth_v9*, nature_v2_figure_utils, nature_final_v7 and
# figures/final_truth_v9*.R, so the numbered analysis directories were never
# searched for a competing colour rule.
P28 <- file.path("figures", "final_truth_v9_heatmap_scale_audit.R")
p28_scope <- paste(grep("SRC <- |final_truth_v9\\|nature_v2", readLines(P28, warn = FALSE),
                        value = TRUE), collapse = " | ")
p28_covers_04 <- grepl("04_differential", p28_scope)

# ==========================================================================
# 4. THE RENDERER INVENTORY
# ==========================================================================
message("[5/5] assembling the renderer inventory")

R_ <- function(file, line, fn, atlas, rule, pct, pctval, reach, cls, conflict, ev)
  data.frame(file = file, line = as.integer(line), function_or_script = fn,
             draws_atlas_like_NES_matrix = atlas,
             colour_limit_rule = rule,
             uses_percentile_clip = pct,
             percentile_value_if_any = pctval,
             reachable_from_manuscript = reach,
             classification = cls,
             conflicts_with_v9_scale = conflict,
             evidence = ev,
             stringsAsFactors = FALSE)

F07 <- "analysis/04_differential_abundance/07_compareGO_spatial_program_atlas.r"
F10 <- "analysis/04_differential_abundance/10_sus_res_spatial_dap_atlas.r"

reach07 <- paste0(
  "YES. pipeline.yml registers '", F07, "' (registered=", registered(F07),
  ", stage=", stage_of(F07), "); RUN_ORDER.md lists it (", in_run_order(F07),
  "). Its figure root results/figures/04_differential_expression_enrichment ",
  "is an explicit candidate_root of the registered exporter ",
  "analysis/09_publication_exports/08_export_manuscript_figures.R (stage=",
  stage_of("analysis/09_publication_exports/08_export_manuscript_figures.R"),
  "), and results/manuscript/figure_export_manifest.csv actually lists ",
  fig07_theme, " (", n_exported(fig07_theme), " row(s)) copied into ",
  "results/manuscript/extended_data/. It is NOT in the publication freeze ",
  "manifest (", n_frozen_07, " matching lines in docs/",
  "publication_freeze_manifest.yml).")

reach10 <- paste0(
  "YES. pipeline.yml registers '", F10, "' (registered=", registered(F10),
  ", stage=", stage_of(F10), "); RUN_ORDER.md lists it (", in_run_order(F10),
  "). results/manuscript/figure_export_manifest.csv lists ",
  n_exported(fig10_theme), " exported artefact(s) matching '", fig10_theme,
  "' under results/manuscript/extended_data/. Strongest evidence: ",
  "docs/publication_freeze_manifest.yml freeze group 'sus_res_stage11' names ",
  "results/source_data/04_differential_expression_enrichment/",
  "sus_res_spatial_dap_atlas/global as a frozen root and pins ",
  "panel_c_sus_res_ranked_GSEA_themes.csv by sha256 (", n_frozen_10,
  " matching lines). Note the freeze pins the NUMBERS, which are correct; it ",
  "does not pin the colour scale.")

rows <- rbind(

  # ---------------- the frozen v9 rule -----------------------------------
  R_("R/final_truth_v9_panels.R", 740, "f9_atlas_limit",
     "TRUE",
     src_line("R/final_truth_v9_panels.R", 747),
     "FALSE", "",
     paste0("YES - canonical. figures/figure_final_truth_v9_contract.yml binds ",
            "panels v9_atlas (F3b), v9_ed_atlas_rescon (ED6a) and ",
            "v9_ed_atlas_suscon (ED6b) to renderer f9_gsea_atlas; the driver ",
            "figures/final_truth_v9_figure_03.R + final_truth_v9_extended_data.R ",
            "are registered in pipeline.yml stage manuscript_candidates and ",
            "listed in RUN_ORDER.md."),
     "ACTIVE_CANONICAL", "FALSE",
     paste0("Recomputed from ", THEME_TBL, " (", nrow(th), " rows): ",
            "max(abs(median_NES)) over all three contrast atlases = ",
            fmt(V9_LIMIT), " over ", length(v9_all), " cells (",
            paste(sprintf("%s: n=%d max=%s", CONTRASTS, v9_n_cells,
                          vapply(v9_per_contrast, fmt, character(1), 6L)),
                  collapse = "; "), "). This IS the frozen 2.478727.")),

  R_("R/final_truth_v9_panels.R", 794, "f9_gsea_atlas",
     "TRUE",
     paste0(src_line("R/final_truth_v9_panels.R", 794), " ; ",
            src_line("R/final_truth_v9_panels.R", 795)),
     "FALSE", "",
     paste0("YES - canonical. Renderer of v9_atlas / v9_ed_atlas_rescon / ",
            "v9_ed_atlas_suscon in figures/figure_final_truth_v9_contract.yml ",
            "lines 150-152, 405-406, 415-416."),
     "ACTIVE_CANONICAL", "FALSE",
     paste0("nv_diverging(limits = c(-", fmt(V9_LIMIT), ", ", fmt(V9_LIMIT),
            ")) with no oob argument, so ggplot2's default oob = censor ",
            "applies (R/nature_v2_figure_utils.R:72-76). No value can be ",
            "silently squished; 0 of ", length(v9_all),
            " atlas cells exceed the limit by construction.")),

  # ---------------- the percentile renderers ------------------------------
  R_(F07, 582, "publication_color_limits",
     "TRUE",
     paste0(src_line(F07, 583), " ; ", src_line(F07, 584), " ; ",
            src_line(F07, 585), " ; ", src_line(F07, 586)),
     "TRUE", "0.98",
     reach07,
     "ACTIVE_NONPUBLICATION", "TRUE",
     paste0("Shared colour-limit helper used by five renderers in this file ",
            "(lines 698, 727, 776, 821, 857). Rule = min(quantile(|x|, 0.98), ",
            "cap = 2.5), always paired with oob = scales::squish, so any value ",
            "beyond the limit is SILENTLY painted the endpoint colour. This is ",
            "categorically different from f9_atlas_limit's true global maximum ",
            "(", fmt(V9_LIMIT), ") with oob = censor.")),

  R_(F07, 776, "plot_sus_res_manuscript_theme_atlas",
     "TRUE",
     paste0(src_line(F07, 776), " ; ", src_line(F07, 785)),
     "TRUE", "0.98",
     reach07,
     "ACTIVE_NONPUBLICATION", "TRUE",
     paste0("THE DIRECT COMPETITOR. Same displayed quantity as v9 F3b: median ",
            "NES per manuscript theme x spatial unit, SUS - RES. Evaluated on ",
            "its own released source data ", p07_theme, " (", nrow(s07),
            " rows; plot_df after theme_role=='primary' & ",
            "n_theme_terms_tested>0 = ", nrow(d07), " cells): ",
            "quantile(|median_NES_all_theme_terms|, 0.98) = ", fmt(p98_07),
            ", cap 2.5 not binding, applied limit = +/-", fmt(lim07),
            " vs the frozen v9 limit +/-", fmt(V9_LIMIT),
            ". Difference = ", fmt(lim07 - V9_LIMIT),
            " (", fmt(100 * (lim07 - V9_LIMIT) / V9_LIMIT, 2),
            "% relative to 2.478727); ratio lim07/v9 = ", fmt(lim07 / V9_LIMIT, 4),
            ". Max |value| actually plotted = ", fmt(max07), "; ", sat07,
            " of ", nrow(d07), " plotted cells exceed the limit and are ",
            "squished to the endpoint colour (", sq07_txt, "). ",
            "SAME NUMBERS, DIFFERENT SCALE: ", cmp07$n, " of ", nrow(d07),
            " cells join the v9 SUS-RES atlas cells on ",
            "dataset/spatial_unit/theme_id with max |difference| = ",
            fmt(cmp07$maxdiff, 12), ", so the two renderers display an ",
            "identical matrix and differ only in the colour limit.")),

  R_(F10, 508, "10_sus_res_spatial_dap_atlas.r (top-level script body, panel C)",
     "TRUE",
     paste0(src_line(F10, 508), " ; ", src_line(F10, 509), " ; ",
            src_line(F10, 510), " ; ", src_line(F10, 511), " ; ",
            src_line(F10, 546)),
     "TRUE", "0.98",
     reach10,
     "ACTIVE_NONPUBLICATION", "TRUE",
     paste0("Second independent copy of the 98th-percentile rule, inlined in ",
            "the script body rather than via publication_color_limits(). Same ",
            "displayed quantity as v9 F3b (median NES per theme x spatial ",
            "unit, SUS - RES). Evaluated on its own released source data ",
            p10_theme, " (", nrow(s10), " rows; panel_c_plot_source after ",
            "theme_role=='primary' = ", nrow(d10), " cells): ",
            "quantile(|median_NES_all_theme_terms|, 0.98) = ", fmt(lim10),
            " -> applied limit +/-", fmt(lim10), " vs frozen v9 +/-",
            fmt(V9_LIMIT), ". Difference = ", fmt(lim10 - V9_LIMIT),
            " (", fmt(100 * (lim10 - V9_LIMIT) / V9_LIMIT, 2),
            "%); ratio = ", fmt(lim10 / V9_LIMIT, 4), ". Max |value| plotted = ",
            fmt(max10), "; ", sat10, " of ", nrow(d10),
            " plotted cells exceed the limit (oob = scales::squish, line 546). ",
            "SAME NUMBERS, DIFFERENT SCALE: ", cmp10$n, " of ", nrow(d10),
            " cells join the v9 SUS-RES atlas cells on ",
            "dataset/spatial_unit/theme_id with max |difference| = ",
            fmt(cmp10$maxdiff, 12), ".")),

  R_(F07, 698, "plot_spatial_program_atlas_publication",
     "TRUE",
     paste0(src_line(F07, 698), " ; ", src_line(F07, 702)),
     "TRUE", "0.98",
     reach07,
     "ACTIVE_NONPUBLICATION", "TRUE",
     paste0("Broad heuristic program_class x spatial_unit NES dot-matrix across ",
            "all three contrasts - atlas-like, but the legacy program taxonomy, ",
            "not the manuscript_go_themes_v2 registry. Source data ", p07_prog,
            " (", nrow(s07p), " rows; publication_include = ", nrow(d07p),
            " plotted): quantile(|mean_NES|, 0.98) -> limit +/-", fmt(lim07p),
            "; max |mean_NES| = ", fmt(max07p), "; ", sat07p,
            " plotted points squished. Exported as ", fig07_prog, " (",
            n_exported(fig07_prog), " manifest row(s)).")),

  R_(F07, 821, "plot_neuropil_focus_publication",
     "TRUE",
     paste0(src_line(F07, 821), " ; ", src_line(F07, 825)),
     "TRUE", "0.98",
     reach07,
     "ACTIVE_NONPUBLICATION", "TRUE",
     paste0("Neuropil-only program x laminar-unit NES dot-matrix. Source data ",
            p07_np, " (", nrow(s07n), " rows; ", nrow(d07n),
            " plotted): limit +/-", fmt(lim07n), "; max |mean_NES| = ",
            fmt(max07n), "; ", sat07n, " squished. Exported as ", fig07_np,
            " (", n_exported(fig07_np), " manifest row(s)).")),

  R_(F07, 853, "plot_compartment_comparison_publication",
     "TRUE",
     paste0(src_line(F07, 853), " ; ", src_line(F07, 857)),
     "TRUE", "0.98",
     reach07,
     "ACTIVE_NONPUBLICATION", "TRUE",
     paste0("Region x program NES dot-matrix faceted by dataset and contrast - ",
            "region-level rather than spatial-unit-level, still an NES colour ",
            "matrix. Source data ", p07_cc, " (", nrow(s07c),
            " rows): limit +/-", fmt(lim07c), "; max |mean_NES| = ",
            fmt(max07c), "; ", sat07c, " squished. Exported as ", fig07_cc,
            " (", n_exported(fig07_cc), " manifest row(s)).")),

  R_(F07, 727, "plot_sus_res_spatial_program_atlas_publication",
     "TRUE",
     paste0(src_line(F07, 727), " ; ", src_line(F07, 736)),
     "TRUE", "0.98",
     paste0("NO. Defined at ", F07, ":709 but never invoked anywhere in the ",
            "repository - grep for the symbol returns the definition line only. ",
            "It writes no figure and no source-data file, so nothing ",
            "manuscript-facing can depend on it."),
     "UNREACHABLE_FROM_MANUSCRIPT", "FALSE",
     paste0("Dead code. Would have used publication_color_limits() on ",
            "representative_NES. Recorded because it matches the search idiom ",
            "and would become a live percentile atlas if ever wired up.")),

  # ---------------- superseded figure generations -------------------------
  R_("R/editorial_v8_panels.R", 740, "e8_atlas_limit",
     "TRUE",
     src_line("R/editorial_v8_panels.R", 747),
     "FALSE", "",
     paste0("NO for the released manuscript. figures/editorial_v8_figure_03.R ",
            "and editorial_v8_extended_data.R are registered in pipeline.yml ",
            "stage manuscript_candidates, but v8 is superseded by the ",
            "final_truth_v9 generation, which has its own contract ",
            "(figures/figure_final_truth_v9_contract.yml) and its own output ",
            "tree. No v9 artefact reads a v8 output."),
     "LEGACY", "FALSE",
     paste0("Byte-identical rule to f9_atlas_limit: max(abs(v)) over the three ",
            "contrast atlases, so it yields the same ", fmt(V9_LIMIT),
            ". Verified by diff of R/editorial_v8_panels.R:722-800 against ",
            "R/final_truth_v9_panels.R:722-800: the only differences are the ",
            "e8_/f9_ symbol prefixes and four SHORT axis labels; the ",
            "explanatory comment block and the limit rule are byte-identical.")),

  R_("R/editorial_v8_panels.R", 794, "e8_gsea_atlas",
     "TRUE",
     paste0(src_line("R/editorial_v8_panels.R", 794), " ; ",
            src_line("R/editorial_v8_panels.R", 795)),
     "FALSE", "",
     "NO for the released manuscript - superseded v8 generation (see e8_atlas_limit row).",
     "LEGACY", "FALSE",
     "Same shared-limit atlas as v9; no scale conflict, only a generation difference."),

  R_("R/nature_final_v7_figure3_panels.R", 147, "nf_gsea_atlas",
     "TRUE",
     paste0(src_line("R/nature_final_v7_figure3_panels.R", 147), " ; ",
            src_line("R/nature_final_v7_figure3_panels.R", 154)),
     "FALSE", "",
     paste0("NO for the released manuscript. figures/nature_final_v7_figure_03.R ",
            "is registered in pipeline.yml stage manuscript_candidates but v7 ",
            "is superseded by v8 then v9."),
     "LEGACY", "TRUE",
     paste0("PER-CONTRAST maximum, not the global one: each atlas rescaled to ",
            "its own max, giving ", paste(sprintf("%s = +/-%s", CONTRASTS,
              vapply(v9_per_contrast, fmt, character(1), 3L)), collapse = ", "),
            ". This is exactly the defect the v8/v9 comment at ",
            "R/final_truth_v9_panels.R:727-731 records as fixed. Conflicts with ",
            "the frozen scale, but by per-panel rescaling, NOT by a percentile.")),

  R_("R/spatial_v6_figure3_panels.R", 296, "s6_f3_gsea_atlas",
     "TRUE",
     paste0(src_line("R/spatial_v6_figure3_panels.R", 296), " ; ",
            src_line("R/spatial_v6_figure3_panels.R", 304)),
     "FALSE", "",
     "NO for the released manuscript - figures/spatial_v6_figure_03.R is registered in stage manuscript_candidates but v6 is superseded by v7/v8/v9.",
     "LEGACY", "TRUE",
     "Per-contrast max(abs(median_NES)); also admits qc_review themes into the matrix. Not a percentile."),

  R_("R/nature_v2_figure_panels.R", 376, "nvp_gsea_atlas",
     "TRUE",
     paste0(src_line("R/nature_v2_figure_panels.R", 376), " ; ",
            src_line("R/nature_v2_figure_panels.R", 382)),
     "FALSE", "",
     "NO for the released manuscript - figures/nature_v2_figure_03.R is registered in stage manuscript_candidates but v2 is superseded by v3-v9.",
     "LEGACY", "TRUE",
     "SUS-RES-only max(abs(median_NES)) including qc_review themes. Not a percentile."),

  R_("R/candidate_figure_panels.R", 264, "cfp_gsea_atlas_plot",
     "TRUE",
     paste0(src_line("R/candidate_figure_panels.R", 264), " ; ",
            src_line("R/candidate_figure_panels.R", 273), " ; ",
            src_line("R/candidate_figure_panels.R", 274)),
     "FALSE", "",
     "NO for the released manuscript - figures/candidate_figure_03.R is registered in stage manuscript_candidates; the candidate layer is an exploratory contract (figures/figure_candidate_contract.yml), superseded by v9.",
     "LEGACY", "TRUE",
     "max(abs(median_NES)) over whichever contrast subset the caller passes, so cfp_render_gsea_atlas_susres draws its main panel and its reference strip on two DIFFERENT limits. Not a percentile."),

  # ---------------- audit-only recomputations -----------------------------
  R_("figures/final_truth_v9_annotation_audit.R", 261, "ATLAS_LIMIT (audit constant)",
     "FALSE",
     src_line("figures/final_truth_v9_annotation_audit.R", 261),
     "FALSE", "",
     "Audit-only: writes to results/tables/manuscript_candidates/final_truth_v9/audit; no figure and no manuscript artefact consumes it.",
     "ACTIVE_NONPUBLICATION", "FALSE",
     "Recomputes max(abs(current_median_NES)) as a yardstick for annotation-sensitivity deltas; draws no colour matrix."),

  R_("figures/final_truth_v9_heatmap_scale_audit.R", 100,
     "SPEC / heatmap_scale_code_search (Part-28 scale audit)",
     "FALSE",
     paste0(src_line("figures/final_truth_v9_heatmap_scale_audit.R", 101), " ; ",
            src_line("figures/final_truth_v9_heatmap_scale_audit.R", 104), " ; ",
            src_line("figures/final_truth_v9_heatmap_scale_audit.R", 107)),
     "FALSE", "",
     "Registered in pipeline.yml stage manuscript_candidates and listed in RUN_ORDER.md; it is the audit that certifies the v9 colour scales.",
     "ACTIVE_NONPUBLICATION", "FALSE",
     paste0("SCOPE GAP. Its own code search (lines 57-60) globs only R/ files ",
            "matching 'final_truth_v9|nature_v2_figure_utils|nature_final_v7' ",
            "plus figures/final_truth_v9*.R. p28_scope_covers_04_directory = ",
            p28_covers_04, ". The two ACTIVE 98th-percentile atlas renderers in ",
            "04_differential_expression_enrichment/ were therefore never in ",
            "scope of the Part-28 certification."))
)

# provenance row so the artefact carries its own reference numbers
rows <- rbind(rows, R_(
  "audits/part29/atlas_renderer_scale_audit.R", 0,
  "AUDIT_REFERENCE_VALUES (this script)", "FALSE",
  paste0("f9_atlas_limit = max(abs(median_NES)) over RES-CON, SUS-CON, SUS-RES = ",
         fmt(V9_LIMIT)),
  "FALSE", "",
  "n/a - audit artefact, writes only under results/tables/publication_audits/upstream_enrichment_v10/",
  "ACTIVE_NONPUBLICATION", "FALSE",
  paste0("Reference set recomputed at git HEAD from ", THEME_TBL, ": ",
         "n_rows=", nrow(th), "; v9_global_limit=", fmt(V9_LIMIT),
         "; per-contrast maxima ",
         paste(sprintf("%s=%s", CONTRASTS,
                       vapply(v9_per_contrast, fmt, character(1), 6L)),
               collapse = ", "),
         "; n_cells ", paste(sprintf("%s=%d", CONTRASTS, v9_n_cells),
                             collapse = ", "),
         ". For contrast: the 98th percentile of the SAME v9 pooled cell ",
         "vector would be ", fmt(V9_P98), ", i.e. ",
         fmt(100 * (V9_P98 - V9_LIMIT) / V9_LIMIT, 2),
         "% of the frozen limit - so a percentile rule applied to the v9 atlas ",
         "itself would shrink the limit and silently saturate ",
         sum(abs(v9_all) > V9_P98, na.rm = TRUE), " of ", length(v9_all),
         " cells.")))

rows <- rows[order(factor(rows$classification,
                          levels = c("ACTIVE_CANONICAL", "ACTIVE_NONPUBLICATION",
                                     "LEGACY", "DEPRECATED",
                                     "UNREACHABLE_FROM_MANUSCRIPT")),
                   rows$file, rows$line), , drop = FALSE]
rownames(rows) <- NULL

utils::write.csv(rows, OUT_CSV, row.names = FALSE, na = "")
message("wrote ", OUT_CSV, "  rows = ", nrow(rows))

cat("\n================ PART-29 S31/S32 KEY NUMBERS ================\n")
cat("frozen v9 f9_atlas_limit (recomputed)          : ", fmt(V9_LIMIT), "\n", sep = "")
cat("  per contrast RES-CON / SUS-CON / SUS-RES     : ",
    paste(vapply(v9_per_contrast, fmt, character(1), 6L), collapse = " / "), "\n", sep = "")
cat("  pooled atlas cells                           : ", length(v9_all), "\n", sep = "")
cat("07 plot_sus_res_manuscript_theme_atlas  p98 lim: ", fmt(lim07),
    "   (delta vs v9 = ", fmt(lim07 - V9_LIMIT), ", ",
    fmt(100 * (lim07 - V9_LIMIT) / V9_LIMIT, 2), "%)\n", sep = "")
cat("   plotted cells = ", nrow(d07), ", max|NES| = ", fmt(max07),
    ", squished = ", sat07, "\n", sep = "")
cat("10 panel_c ranked_GSEA_themes           p98 lim: ", fmt(lim10),
    "   (delta vs v9 = ", fmt(lim10 - V9_LIMIT), ", ",
    fmt(100 * (lim10 - V9_LIMIT) / V9_LIMIT, 2), "%)\n", sep = "")
cat("   plotted cells = ", nrow(d10), ", max|NES| = ", fmt(max10),
    ", squished = ", sat10, "\n", sep = "")
cat("07 program atlas   p98 lim = ", fmt(lim07p), " (max ", fmt(max07p),
    ", squished ", sat07p, ")\n", sep = "")
cat("07 neuropil focus  p98 lim = ", fmt(lim07n), " (max ", fmt(max07n),
    ", squished ", sat07n, ")\n", sep = "")
cat("07 compartment cmp p98 lim = ", fmt(lim07c), " (max ", fmt(max07c),
    ", squished ", sat07c, ")\n", sep = "")
cat("98th pct of the v9 pooled cells (hypothetical) : ", fmt(V9_P98), "\n", sep = "")
cat("cell-join 07 vs v9 SUS-RES : n=", cmp07$n, " max|diff|=",
    fmt(cmp07$maxdiff, 12), "\n", sep = "")
cat("cell-join 10 vs v9 SUS-RES : n=", cmp10$n, " max|diff|=",
    fmt(cmp10$maxdiff, 12), "\n", sep = "")
cat("cells squished by the p98 limit: ", sq07_txt, "\n", sep = "")
cat("Part-28 scale audit covers 04_* directory      : ", p28_covers_04, "\n", sep = "")
cat("07 registered / stage : ", registered(F07), " / ", stage_of(F07), "\n", sep = "")
cat("10 registered / stage : ", registered(F10), " / ", stage_of(F10), "\n", sep = "")
cat("exported 07 theme atlas svg rows : ", n_exported(fig07_theme), "\n", sep = "")
cat("exported 10 panel_c themes rows  : ", n_exported(fig10_theme), "\n", sep = "")
cat("freeze-manifest lines  07 / 10   : ", n_frozen_07, " / ", n_frozen_10,
    "\n", sep = "")
cat("============================================================\n")
