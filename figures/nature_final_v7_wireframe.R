#!/usr/bin/env Rscript

# Part-23 step 1: TRUE-SIZE WIREFRAMES, built before any scientific panel is
# rendered (brief section 5).
#
# Each box states its letter, mm geometry, intended content, the minimum font
# size the panel is expected to carry, and whether the content is unique to the
# main figure or also represented in Extended Data. The point is to prove the
# composition fits 183 x <=170 mm with a sensible area hierarchy BEFORE any
# expensive rendering, and to expose an area hierarchy that contradicts the
# declared scientific priority.
#
# DOWNSTREAM ONLY: draws rectangles and text. No data, no model, no inference.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "nature_final_v7_figure_utils.R"))
suppressPackageStartupMessages({ library(ggplot2); library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/nature_final_v7_wireframe.R")

OUT <- function(...) {
  d <- path_results("figures", "manuscript_candidates", "nature_final_v7", "wireframe")
  dir_create(d)
  file.path(d, ...)
}
TBL <- function(...) {
  d <- path_results("tables", "manuscript_candidates", "nature_final_v7")
  dir_create(d)
  file.path(d, ...)
}

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] wireframes for ", relative_to(s7e_contract_path()))
  quit(save = "no", status = 0L)
}

ct <- s7e_contract()
panels <- ct$panels
names(panels) <- vapply(panels, function(p) as.character(p$id), character(1))

# what each panel is expected to carry, and its declared main/ED role
CONTENT <- c(
  v7_schematic = "anatomy + sampling; teaches compartment / region / layer",
  v7_depth = "proteome depth per compartment",
  v7_pca = "global PCA, compartment encoded",
  v7_fingerprint = "baseline spatial molecular fingerprint (phenotype-blind)",
  v7_compartment = "compartment identity markers",
  v7_bilateral_main = "L-vs-R scatter + contrast summary, regional then laminar",
  v7_external_main = "external hippocampal validation, expected pairings",
  v7_internal_main = "internal anatomical GO programs",
  v7_dap_track = "sparse single-protein DA: canonical + claimable counts",
  v7_atlas = "GSEA theme x spatial-unit atlas (SUS-RES)",
  v7_bridge = "ONE hippocampus + 3 anchored callouts + 3-contrast strips",
  v7_curve_syn = "GSEA curve: neuropil CA3-SR synaptic",
  v7_curve_rna = "GSEA curve: soma CA2-SP mRNA processing",
  v7_curve_ox = "GSEA curve: microglia CA1 OXPHOS",
  v7_prot_syn = "leading-edge proteins: CA3-SR synaptic",
  v7_prot_rna = "leading-edge proteins: CA2-SP mRNA processing",
  v7_prot_ox = "leading-edge proteins: CA1 OXPHOS")

MINPT <- c(v7_schematic = 5.0, v7_depth = 5.0, v7_pca = 5.0, v7_fingerprint = 5.0,
           v7_compartment = 5.0, v7_bilateral_main = 5.0, v7_external_main = 5.0,
           v7_internal_main = 5.0, v7_dap_track = 5.0, v7_atlas = 5.0,
           v7_bridge = 5.0, v7_curve_syn = 5.0, v7_curve_rna = 5.0,
           v7_curve_ox = 5.0, v7_prot_syn = 5.0, v7_prot_rna = 5.0,
           v7_prot_ox = 5.0)

# declared Figure-2 scientific priority (brief section 21); 1 = most important
PRIORITY <- c(v7_fingerprint = 1, v7_schematic = 2, v7_compartment = 3,
              v7_bilateral_main = 4, v7_external_main = 5, v7_internal_main = 6,
              v7_pca = 7, v7_depth = 8)

fam <- nv_palette()$typography$family
rows <- list()

for (f in ct$figures) {
  W <- as.numeric(f$width_mm); H <- as.numeric(f$height_mm)
  d <- do.call(rbind, lapply(f$layout, function(it) {
    id <- as.character(it$panel)
    data.frame(figure = as.character(f$name), label = as.character(it$label),
               panel = id, x = as.numeric(it$x), y = as.numeric(it$y),
               w = as.numeric(it$w), h = as.numeric(it$h),
               content = unname(CONTENT[id] %||% ""),
               min_pt = unname(MINPT[id] %||% 5),
               role = as.character(panels[[id]]$main_or_ed %||% ""),
               stringsAsFactors = FALSE)
  }))
  d$area <- d$w * d$h
  d$share <- d$area / (W * H)
  d$fits <- (d$x + d$w) <= W + 1e-9 & (d$y + d$h) <= H + 1e-9

  lab <- sprintf("%s\n%g x %g mm at (%g, %g)\n%s\n>= %.1f pt | %s",
                 d$label, d$w, d$h, d$x, d$y,
                 vapply(d$content, function(s)
                   paste(strwrap(s, 34), collapse = "\n"), character(1)),
                 d$min_pt, d$role)

  p <- ggplot2::ggplot(d) +
    ggplot2::geom_rect(ggplot2::aes(xmin = x, xmax = x + w,
                                    ymin = y, ymax = y + h),
                       fill = "grey96", colour = "grey35", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(x = x + w / 2, y = y + h / 2, label = lab),
                       family = fam, size = nv_size(5.2), lineheight = 1.25,
                       colour = "grey15") +
    ggplot2::annotate("rect", xmin = 0, xmax = W, ymin = 0, ymax = H,
                      fill = NA, colour = "grey70", linewidth = 0.4) +
    ggplot2::scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
    ggplot2::scale_y_reverse(limits = c(H, 0), expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = sprintf("%s  -  %g x %g mm  -  %d panels  -  %.1f%% of page covered",
                                  f$name, W, H, nrow(d), 100 * sum(d$share))) +
    ggplot2::theme_void(base_family = fam) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 7, face = "bold",
                                                      hjust = 0,
                                                      margin = ggplot2::margin(b = 2)))

  svg <- OUT(paste0(sub("_V7$", "", f$name), "_WIREFRAME.svg"))
  nv_save_panel(p, svg, W, H)
  nv_pdf(svg, sub("[.]svg$", ".pdf", svg))
  rows[[length(rows) + 1L]] <- d
}

wf <- do.call(rbind, rows)
wf$priority_rank <- unname(PRIORITY[wf$panel])
write_csv_safe(wf, TBL("nature_final_v7_wireframe_geometry.csv"))

cat("\n===== TRUE-SIZE WIREFRAMES =====\n")
for (fg in unique(wf$figure)) {
  z <- wf[wf$figure == fg, ]
  cat("\n---", fg, "---\n")
  cat("page coverage:", sprintf("%.1f%%", 100 * sum(z$share)),
      "  all panels inside the page:", all(z$fits), "\n")
  print(z[order(-z$area), c("label", "panel", "w", "h", "area", "share",
                            "priority_rank")],
        row.names = FALSE, digits = 3)
}

cat("\n===== AREA HIERARCHY vs DECLARED PRIORITY (Figure 2) =====\n")
z <- wf[!is.na(wf$priority_rank), ]
z <- z[order(z$priority_rank), ]
z$area_rank <- rank(-z$area, ties.method = "min")
print(z[, c("label", "panel", "priority_rank", "area_rank", "w", "h", "area")],
      row.names = FALSE)
inv <- z$area_rank > z$priority_rank + 1
cat("\npanels whose area rank trails their priority by more than one place:",
    if (any(inv)) paste(z$panel[inv], collapse = ", ") else "none", "\n")
# the specific inversion Part 22 flagged
pca <- z$area[z$panel == "v7_pca"]; cmp <- z$area[z$panel == "v7_compartment"]
cat("PCA area", pca, "mm2 vs compartment identity", cmp, "mm2 -> ",
    if (pca < cmp) "FIXED: PCA no longer outranks compartment identity"
    else "STILL INVERTED", "\n")

if (!all(wf$fits)) stop("wireframe: panel(s) overflow the page", call. = FALSE)
cat("\nwritten:", relative_to(OUT()), "\n")
