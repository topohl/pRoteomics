#!/usr/bin/env Rscript

# Part-24 sections 31 and 54: the reviewable raster outputs.
#
# The delivered artwork is the vector SVG and the cairo PDF. These previews
# exist only so the composition can be checked at true print size without a PDF
# rasteriser: 300 dpi is what a referee sees on screen, 900 dpi is what settles
# whether 5 pt text is genuinely legible in print. Part 22 mis-read crisp 4.8 pt
# text as illegible because it was judged at 300 dpi only, so BOTH are emitted.
#
# A contact sheet puts all nine pages on one image, which is the only view in
# which main/Extended Data balance is visible at a glance.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "final_truth_v9_figure_utils.R"))
suppressPackageStartupMessages({ library(magick) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_previews.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] editorial_v8 previews and contact sheet")
  quit(save = "no", status = 0L)
}

ROOT <- path_results("figures", "manuscript_candidates", "final_truth_v9")
PREV <- file.path(ROOT, "previews")
dir_create(file.path(PREV, "300dpi"))
dir_create(file.path(PREV, "900dpi"))

ct <- s9f_contract()
rows <- list(); sheet <- list()
for (f in ct$figures) {
  nm <- as.character(f$name)
  svg <- file.path(ROOT, as.character(f$figure_key), "assembled",
                   paste0(nm, ".svg"))
  if (!file.exists(svg)) {
    message("[previews] missing: ", relative_to(svg)); next
  }
  flat <- function(im) image_flatten(image_background(im, "white"))
  p3 <- file.path(PREV, "300dpi", paste0(nm, "_300dpi.png"))
  p9 <- file.path(PREV, "900dpi", paste0(nm, "_900dpi.png"))
  i3 <- flat(image_read(svg, density = 300))
  i9 <- flat(image_read(svg, density = 900))
  image_write(i3, p3, format = "png")
  image_write(i9, p9, format = "png")
  sheet[[nm]] <- image_annotate(
    image_border(image_scale(i3, "620"), "white", "6x6"),
    nm, size = 15, gravity = "northwest", location = "+8+2", color = "grey25")
  rows[[length(rows) + 1L]] <- data.frame(
    figure = nm,
    page_mm = sprintf("%g x %g", as.numeric(f$width_mm),
                      as.numeric(f$height_mm)),
    px_300dpi = paste(image_info(i3)$width, "x", image_info(i3)$height),
    px_900dpi = paste(image_info(i9)$width, "x", image_info(i9)$height),
    preview_300dpi = relative_to(p3), preview_900dpi = relative_to(p9),
    stringsAsFactors = FALSE)
}
if (!length(rows)) stop("no assembled editorial_v8 SVGs found", call. = FALSE)

cs <- file.path(PREV, "final_truth_v9_contact_sheet.png")
image_write(image_montage(image_join(sheet), tile = "3x", geometry = "640x640+8+8",
                          bg = "white"), cs, format = "png")

out <- do.call(rbind, rows)
tab <- path_results("tables", "manuscript_candidates", "final_truth_v9")
dir_create(tab)
write_csv_safe(out, file.path(tab, "preview_render_index.csv"))

cat("\n===== EDITORIAL V8 PREVIEWS =====\n")
print(out[, c("figure", "page_mm", "px_300dpi", "px_900dpi")], row.names = FALSE)
cat("\ncontact sheet:", relative_to(cs), "\n")
