#!/usr/bin/env Rscript

# Part-24: the vector-export audit, and the final story write-up.
#
# Parts 16-23 shipped PDFs that were a single full-page 300 dpi raster with zero
# font objects, because nv_pdf() falls back to ImageMagick when rsvg is absent.
# This audit proves the editorial_v8 PDFs are genuine vector: embedded fonts,
# real text-showing operators, and no page-sized raster.
#
# qpdf is the only PDF tool available here (no poppler, no ghostscript), so the
# probe decompresses with --qdf and reads the object structure directly.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "editorial_v8_export.R"))
source(repo_path("R", "final_truth_v9_figure_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_vector_audit.R")

OUT <- function(...) {
  d <- path_results("tables", "manuscript_candidates", "final_truth_v9")
  dir_create(d)
  file.path(d, ...)
}
REP <- function(...) {
  d <- path_results("reports", "manuscript_candidates", "final_truth_v9")
  dir_create(d)
  file.path(d, ...)
}

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] vector export audit over editorial_v8")
  quit(save = "no", status = 0L)
}

ct <- s9f_contract()
pdfs <- character(0)
for (f in ct$figures) {
  p <- path_results("figures", "manuscript_candidates", "final_truth_v9",
                    as.character(f$figure_key), "assembled",
                    paste0(as.character(f$name), ".pdf"))
  if (file.exists(p)) pdfs <- c(pdfs, p)
}
if (!length(pdfs)) stop("no editorial_v8 PDFs found", call. = FALSE)

aud <- e8_vector_audit(pdfs)
write_csv_safe(aud, OUT("vector_export_audit.csv"))

cat("\n===== VECTOR EXPORT AUDIT =====\n")
print(aud[, c("figure", "page_size_mm", "font_count", "embedded_fonts",
              "page_sized_raster_present", "largest_raster_dimensions",
              "text_show_operators", "vector_export_pass")], row.names = FALSE)

# for contrast, probe the Part-23 PDFs so the fix is demonstrable
old <- character(0)
for (nm in c("figure_02/assembled/F2_NATURE_FINAL_V7.pdf",
             "figure_03/assembled/F3_NATURE_FINAL_V7.pdf")) {
  p <- path_results("figures", "manuscript_candidates", "nature_final_v7", nm)
  if (file.exists(p)) old <- c(old, p)
}
if (length(old)) {
  cat("\n--- the Part-23 PDFs, for comparison ---\n")
  print(e8_vector_audit(old)[, c("figure", "font_count", "embedded_fonts",
                                 "page_sized_raster_present",
                                 "largest_raster_dimensions",
                                 "vector_export_pass")], row.names = FALSE)
}


# ------------------------------------------------------- SVG write integrity
#
# svglite writes incrementally, so an interrupted write on a network share
# leaves a syntactically truncated file that every downstream step accepts
# until something tries to rasterise it. Assert the closing tag, and that the
# file is not implausibly small, for every assembled figure and every panel.
svgs <- sort(list.files(
  path_results("figures", "manuscript_candidates", "final_truth_v9"),
  pattern = "[.]svg$", recursive = TRUE, full.names = TRUE))
svg_ok <- vapply(svgs, function(p) {
  n <- file.info(p)$size
  con <- file(p, "rb"); on.exit(close(con))
  seek(con, max(0, n - 64))
  last <- rawToChar(readBin(con, "raw", 64L))
  isTRUE(n > 2000) && grepl("</svg>", last, fixed = TRUE)
}, logical(1))
if (!all(svg_ok)) {
  bad <- basename(svgs[!svg_ok])
  stop("truncated or empty SVG (incomplete write): ",
       paste(bad, collapse = ", "), call. = FALSE)
}
cat("SVG write integrity:", sum(svg_ok), "of", length(svg_ok),
    "complete (closing tag present)\n")

cat("\nfigures passing:", sum(aud$vector_export_pass), "of", nrow(aud), "\n")

# ---------------------------------------------------------------- story map
story <- c(
"# Final figure story (final truth v9)",
"",
"## Figure 2",
"",
"**The spatial proteomics platform resolves a reproducible multi-resolution",
"hippocampal molecular architecture.**",
"",
"Panel a states what was sampled and teaches the grammar the rest of the figure",
"uses: compartment by glyph, region by bracket, layer by position and a short",
"key. b shows coverage is sufficient in all three compartments. c shows the",
"global proteome is already organised before any spatial signature is inspected.",
"d is the centrepiece: the actual control-animal spatial molecular",
"fingerprint, with",
"genes chosen only by prespecified CON-only anatomical contrasts. e confirms each",
"compartment carries its expected molecular identity. f shows paired hemispheres",
"reproduce the architecture, with coarse regional identity stronger than fine",
"laminar identity - the weaker CA1 laminar result is deliberately left visible.",
"g validates the same spatial identities against an external hippocampal",
"reference; h characterises them functionally against canonical GO",
"programs and is not independent validation.",
"",
"## Figure 3",
"",
"**Later resilient and susceptible outcomes are associated with sparse",
"protein-level effects but stronger spatially selective coordinated molecular",
"programs across distinct hippocampal compartments.**",
"",
"Panel a establishes the premise numerically: most spatial units carry no",
"FDR-supported SUS-RES protein at all, and the apparent protein-level",
"concentration in CA2-SLM falls sharply once robustness qualification is",
"applied. b, sharing a column geometry with a, shows that program-level",
"structure is nevertheless present across the same units. c places the three",
"representative programs on one hippocampal schematic. The figure then reads as",
"three vertical biological columns: each column pairs its direct ranked",
"enrichment evidence (d, e, f) with the leading-edge proteins that carry it",
"(g, h, i), on one shared log2 fold-change scale so the three programs can be",
"compared directly.",
"",
"## Supporting Extended Data",
"",
"**QC and systems analyses constrain the interpretation.** The CA2-SLM",
"protein-level concentration weakens after robustness qualification (ED3).",
"WGCNA modules remain spatially organised, and while most neuropil modules share",
"a descriptive directional pattern, no module x contrast test is FDR-supported.",
"Spatial molecular networks show no detectable whole-network group",
"difference, and no edge-behaviour association survived multiple-testing",
"correction (ED8).",
"",
"## The defensible conclusion",
"",
"Stress outcome is associated with selective local molecular-program differences",
"superimposed on a hippocampal spatial molecular architecture that remained",
"evident across groups: bilaterally reproducible, and with no whole-network",
"group difference detectable at three animals per group.",
"",
"This is deliberately NOT a claim of global reprogramming, redistribution,",
"relocation, network rewiring, or any cell-intrinsic microglial property: the",
"microglia-enriched compartment is a local microenvironment measurement.")
writeLines(story, REP("final_figure_story_v9.md"))
cat("written:", relative_to(REP("final_figure_story_v9.md")), "\n")

if (any(aud$page_sized_raster_present)) {
  stop("VECTOR EXPORT FAILED: page-sized raster present in ",
       paste(aud$figure[aud$page_sized_raster_present], collapse = ", "),
       call. = FALSE)
}
