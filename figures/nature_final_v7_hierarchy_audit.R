#!/usr/bin/env Rscript

# Part-23: the main/ED hierarchy audit, the spatial x phenotype coverage matrix
# and the main-figure story map.
#
# The hierarchy rule this audit enforces: a main panel and an Extended Data
# panel must never be the same graphic at two sizes. Part 22 measured four such
# duplications in the spatial_v6 layer, three of which were LARGER in Extended
# Data than in the main figure. byte_identical must be FALSE for every repeated
# analysis unless an exception is documented here.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "nature_final_v7_figure_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/nature_final_v7_hierarchy_audit.R")

OUT <- function(...) {
  d <- path_results("tables", "manuscript_candidates", "nature_final_v7")
  dir_create(d)
  file.path(d, ...)
}
REP <- function(...) {
  d <- path_results("reports", "manuscript_candidates", "nature_final_v7")
  dir_create(d)
  file.path(d, ...)
}

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] main/ED hierarchy audit over ", relative_to(s7e_contract_path()))
  quit(save = "no", status = 0L)
}

ct <- s7e_contract()
panels <- ct$panels
names(panels) <- vapply(panels, function(p) as.character(p$id), character(1))

use <- do.call(rbind, lapply(ct$figures, function(f) do.call(rbind, lapply(
  f$layout, function(it) data.frame(
    figure = as.character(f$name), key = as.character(f$figure_key),
    label = as.character(it$label), panel = as.character(it$panel),
    w = as.numeric(it$w), h = as.numeric(it$h), stringsAsFactors = FALSE)))))
use$area <- use$w * use$h
use$side <- ifelse(use$key == "extended_data", "ED", "MAIN")

svg_of <- function(key, panel, w, h) {
  d <- path_results("figures", "manuscript_candidates", "nature_final_v7", key,
                    "panels")
  cand <- c(sprintf("%s_%gx%g.svg", panel, w, h), sprintf("%s.svg", panel))
  hit <- cand[file.exists(file.path(d, cand))]
  if (!length(hit)) return(NA_character_)
  file.path(d, hit[1])
}
md5 <- function(p) if (is.na(p) || !file.exists(p)) NA_character_ else
  unname(tools::md5sum(p))

# Explicit main <-> ED pairing by ANALYSIS SUBJECT. Matching on the role string
# alone would report every analysis as one-sided, because main and ED panels
# deliberately carry different roles; the point of this audit is to show the
# pairs and prove they are not the same graphic.
PAIRS <- list(
  list(a = "bilateral reproducibility", m = "v7_bilateral_main",  e = "v7_ed_bilateral_full"),
  list(a = "external spatial validation", m = "v7_external_main", e = "v7_ed_external_full"),
  list(a = "internal anatomical programs", m = "v7_internal_main", e = "v7_ed_internal_full"),
  list(a = "baseline spatial fingerprint", m = "v7_fingerprint",   e = "v7_ed_fingerprint_full"),
  list(a = "GSEA curve, neuropil CA3-SR", m = "v7_curve_syn",      e = "v7_ed_gsea_curve_syn"),
  list(a = "GSEA curve, soma CA2-SP",     m = "v7_curve_rna",      e = "v7_ed_gsea_curve_rna"),
  list(a = "GSEA curve, microglia CA1",   m = "v7_curve_ox",       e = "v7_ed_gsea_curve_ox"),
  list(a = "GSEA spatial atlas",          m = "v7_atlas",          e = "v7_ed_atlas_rescon"),
  list(a = "PCA",                         m = "v7_pca",            e = NA_character_),
  list(a = "compartment identity",        m = "v7_compartment",    e = NA_character_),
  list(a = "proteome depth",              m = "v7_depth",          e = NA_character_),
  list(a = "anatomy and sampling",        m = "v7_schematic",      e = NA_character_),
  list(a = "sparse single-protein DA",    m = "v7_dap_track",      e = NA_character_),
  list(a = "anatomical program bridge",   m = "v7_bridge",         e = NA_character_),
  list(a = "leading-edge proteins, neuropil",  m = "v7_prot_syn",  e = NA_character_),
  list(a = "leading-edge proteins, soma",      m = "v7_prot_rna",  e = NA_character_),
  list(a = "leading-edge proteins, microglia", m = "v7_prot_ox",   e = NA_character_))

EXCEPT <- c(
  "PCA" = "ED PCA deliberately REMOVED (brief section 13): the main PCA already communicates everything needed and an enlarged identical copy would be exactly the defect Part 22 found.",
  "compartment identity" = "the canonical marker source contains 10 markers x 3 compartments and the main panel shows all 30 cells, so Extended Data has nothing further to add.",
  "proteome depth" = "a single coverage statement; the main panel is the complete evidence.",
  "anatomy and sampling" = "a schematic, not an analysis; there is no deeper version.",
  "sparse single-protein DA" = "the track already shows every spatial unit; the per-protein detail lives in ED7 and ED3 rather than as a larger copy.",
  "anatomical program bridge" = "an orientation device unique to the main figure.",
  "leading-edge proteins, neuropil" = "the full leading-edge inventory is a supplementary table, not a larger figure panel.",
  "leading-edge proteins, soma" = "as above.",
  "leading-edge proteins, microglia" = "as above.")

rows <- lapply(PAIRS, function(pp) {
  m <- use[use$panel == pp$m, , drop = FALSE]
  e <- if (is.na(pp$e)) use[0, ] else use[use$panel == pp$e, , drop = FALSE]
  mh <- if (nrow(m)) md5(svg_of(m$key[1], m$panel[1], m$w[1], m$h[1])) else NA
  eh <- if (nrow(e)) md5(svg_of(e$key[1], e$panel[1], e$w[1], e$h[1])) else NA
  ident <- if (is.na(mh) || is.na(eh)) FALSE else identical(mh, eh)
  data.frame(
    analysis = pp$a,
    main_panel = if (nrow(m)) m$panel[1] else NA_character_,
    main_size_mm = if (nrow(m)) sprintf("%g x %g", m$w[1], m$h[1]) else NA_character_,
    ed_panel = if (nrow(e)) e$panel[1] else NA_character_,
    ed_size_mm = if (nrow(e)) sprintf("%g x %g", e$w[1], e$h[1]) else NA_character_,
    byte_identical = ident,
    main_unique_content = if (nrow(m))
      trimws(as.character(panels[[m$panel[1]]]$narrative %||% "")) else "",
    ed_unique_content = if (nrow(e))
      trimws(as.character(panels[[e$panel[1]]]$narrative %||% ""))
      else unname(EXCEPT[pp$a] %||% "main only"),
    hierarchy_valid = if (!nrow(e)) "VALID - main only, exception documented"
      else if (ident) "INVALID - identical graphic on both sides"
      else "VALID - different graphic, different editorial role",
    recommended_fix = if (nrow(e) && ident)
      "split the content: main concise, ED complete" else "none required",
    stringsAsFactors = FALSE)
})
aud <- do.call(rbind, rows)
write_csv_safe(aud, OUT("main_ed_panel_hierarchy_audit.csv"))

# ---------------------------------------------------------------------------
# spatial x phenotype coverage matrix
# ---------------------------------------------------------------------------
has <- function(s, tok) if (is.na(s) || !nzchar(s)) "" else
  if (grepl(tok, s, fixed = TRUE)) "x" else ""
mat <- do.call(rbind, lapply(seq_len(nrow(use)), function(i) {
  u <- use[i, ]
  p <- panels[[u$panel]]
  sax <- as.character(p$spatial_axes %||% "")
  pax <- as.character(p$phenotype_axes %||% "")
  role <- as.character(p$role %||% "")
  blind <- identical(pax, "none_intentional") || identical(pax, "none_applicable")
  data.frame(
    figure = u$figure, label = u$label, panel = u$panel,
    compartment = has(sax, "compartment"), region = has(sax, "region"),
    layer = has(sax, "layer"),
    CON = if (!blind && grepl("CON", pax)) "x" else "",
    RES = if (!blind && grepl("RES", pax)) "x" else "",
    SUS = if (!blind && grepl("SUS", pax)) "x" else "",
    `RES-CON` = if (grepl("RES-CON", pax)) "x" else "",
    `SUS-CON` = if (grepl("SUS-CON", pax)) "x" else "",
    `SUS-RES` = if (grepl("SUS-RES", pax)) "x" else "",
    baseline_spatial_identity = if (blind && nzchar(sax)) "x" else "",
    stress_sensitive_spatial_identity =
      if (!blind && grepl("region|layer", sax)) "x" else "",
    direct_proteins = if (grepl("protein|leading_edge", role)) "x" else "",
    program_level_biology = if (grepl("program|gsea|evidence", role)) "x" else "",
    systems_level_biology = if (grepl("wgcna|network|coupling", role)) "x" else "",
    qc_robustness = if (grepl("qc|robust|reproduc|technical|measurement|audit",
                              role)) "x" else "",
    check.names = FALSE, stringsAsFactors = FALSE)
}))
write_csv_safe(mat, OUT("final_spatial_phenotype_coverage_matrix.csv"))

# ---------------------------------------------------------------------------
# main-figure story map
# ---------------------------------------------------------------------------
story <- c(
  "# Final main-figure story map",
  "",
  "Generated by `figures/nature_final_v7_hierarchy_audit.R`. Each line states",
  "what that panel contributes to the argument. If any transition would need",
  "the words \"separately, we also analysed\", the figure logic is wrong.",
  "", "## Figure 2", "",
  "> Does the spatial proteomics experiment faithfully resolve biologically",
  "> meaningful hippocampal molecular architecture?", "")
f2 <- use[use$figure == "F2_NATURE_FINAL_V7", ]
f2 <- f2[order(f2$label), ]
F2TXT <- c(a = "tells the reader what was sampled, and teaches the compartment / region / layer grammar used by every later panel.",
           b = "shows the coverage is sufficient in all three compartments.",
           c = "establishes that the global proteome is broadly organised before any specific spatial signature is inspected.",
           d = "shows the actual baseline spatial molecular fingerprint - the direct answer to what the hippocampal proteome looks like in space.",
           e = "verifies that each compartment carries its expected molecular identity.",
           f = "tests paired-hemisphere reproducibility, separating coarse regional from fine laminar identity.",
           g = "validates the spatial identities against an external hippocampal reference.",
           h = "validates that the measured anatomical differences correspond to coherent biological programs.")
for (i in seq_len(nrow(f2))) {
  story <- c(story, sprintf("**%s** (`%s`, %g x %g mm) %s", f2$label[i],
                            f2$panel[i], f2$w[i], f2$h[i],
                            F2TXT[[f2$label[i]]]))
}
story <- c(story, "", "## Figure 3", "",
  "> What molecular programs distinguish later resilient and susceptible",
  "> outcomes, where do they occur, and which proteins carry them?", "")
f3 <- use[use$figure == "F3_NATURE_FINAL_V7", ]
f3 <- f3[order(f3$label), ]
F3TXT <- c(a = "establishes the premise: individual-protein differential abundance is sparse.",
           b = "shows that coordinated program-level effects are nevertheless present, across the same spatial units.",
           c = "maps the three representative programs onto hippocampal anatomy and carries their three-group trajectories.",
           d = "gives direct ranked-enrichment evidence in the neuropil.",
           e = "gives direct ranked-enrichment evidence in the neuronal soma.",
           f = "gives direct ranked-enrichment evidence in the microglia-enriched ROI.",
           g = "shows the proteins carrying the neuropil program.",
           h = "shows the proteins carrying the soma program.",
           i = "shows the proteins carrying the microglia-enriched program.")
for (i in seq_len(nrow(f3))) {
  story <- c(story, sprintf("**%s** (`%s`, %g x %g mm) %s", f3$label[i],
                            f3$panel[i], f3$w[i], f3$h[i],
                            F3TXT[[f3$label[i]]]))
}
story <- c(story, "",
  "## Why a, b run as one block in Figure 3", "",
  "Panels a and b share their column geometry exactly - same 18 spatial units,",
  "same compartment and region boundaries, no gutter between them - so the",
  "reader sees few significant single proteins directly above richer",
  "program-level structure, in the same columns.")
writeLines(story, REP("final_main_figure_story_map.md"))

cat("\n===== MAIN / ED HIERARCHY AUDIT =====\n")
cat("analyses audited:", nrow(aud), "\n")
print(table(aud$hierarchy_valid))
bad <- aud[aud$byte_identical, , drop = FALSE]
cat("\nbyte-identical main/ED pairs:", nrow(bad), "(target 0)\n")
if (nrow(bad)) print(bad[, c("analysis", "main_panel", "ed_panel")], row.names = FALSE)
both <- aud[!is.na(aud$main_panel) & !is.na(aud$ed_panel), , drop = FALSE]
cat("\nanalyses appearing on BOTH sides:", nrow(both), "\n")
if (nrow(both)) print(both[, c("analysis", "main_panel", "main_size_mm",
                               "ed_panel", "ed_size_mm", "byte_identical")],
                      row.names = FALSE)
cat("\ncoverage matrix rows:", nrow(mat), " columns:", ncol(mat), "\n")
cat("\nwritten:\n  ", relative_to(OUT("main_ed_panel_hierarchy_audit.csv")),
    "\n  ", relative_to(OUT("final_spatial_phenotype_coverage_matrix.csv")),
    "\n  ", relative_to(REP("final_main_figure_story_map.md")), "\n")
if (nrow(bad)) stop("main/ED hierarchy violated: identical graphic on both sides",
                    call. = FALSE)
