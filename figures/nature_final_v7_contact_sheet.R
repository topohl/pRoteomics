#!/usr/bin/env Rscript

# Part-23: contact sheet plus true-size raster previews.
#
# Emits a 300 dpi and a 900 dpi PNG for every assembled figure, and one contact
# sheet showing the whole family at a glance. The 900 dpi pass exists because a
# 300 dpi raster softens 5 pt type enough to look illegible when it is not -
# Part 22 nearly mis-reported that as a print defect.
#
# DOWNSTREAM ONLY: rasterises existing SVGs. No data, no inference.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "nature_final_v7_figure_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/nature_final_v7_contact_sheet.R")

OUT <- function(...) {
  d <- path_results("figures", "manuscript_candidates", "nature_final_v7",
                    "previews")
  dir_create(d)
  file.path(d, ...)
}

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] contact sheet and previews for nature_final_v7")
  quit(save = "no", status = 0L)
}
if (!requireNamespace("magick", quietly = TRUE)) {
  message("magick unavailable; skipping previews")
  quit(save = "no", status = 0L)
}

ct <- s7e_contract()
rows <- list()
thumbs <- list()
for (f in ct$figures) {
  key <- as.character(f$figure_key)
  svg <- path_results("figures", "manuscript_candidates", "nature_final_v7", key,
                      "assembled", paste0(as.character(f$name), ".svg"))
  if (!file.exists(svg)) next
  for (dpi in c(300, 900)) {
    im <- magick::image_background(magick::image_read(svg, density = dpi),
                                   "white", flatten = TRUE)
    magick::image_write(im, OUT(sprintf("%s_%ddpi.png", f$name, dpi)),
                        format = "png")
    if (dpi == 300) {
      info <- magick::image_info(im)
      thumbs[[as.character(f$name)]] <- magick::image_annotate(
        magick::image_border(magick::image_scale(im, "600"), "grey60", "2x2"),
        as.character(f$name), size = 22, gravity = "northwest",
        location = "+6+6", color = "black")
      rows[[length(rows) + 1L]] <- data.frame(
        figure = as.character(f$name), figure_key = key,
        width_mm = as.numeric(f$width_mm), height_mm = as.numeric(f$height_mm),
        px_300dpi = sprintf("%d x %d", info$width, info$height),
        n_panels = length(f$layout), stringsAsFactors = FALSE)
    }
  }
}

sheet <- magick::image_append(
  magick::image_append(magick::image_join(thumbs[1:3]), stack = FALSE),
  stack = TRUE)
if (length(thumbs) > 3) {
  grp <- split(seq_along(thumbs), ceiling(seq_along(thumbs) / 3))
  bands <- lapply(grp, function(ix)
    magick::image_append(magick::image_join(thumbs[ix]), stack = FALSE))
  sheet <- magick::image_append(magick::image_join(bands), stack = TRUE)
}
sheet <- magick::image_background(sheet, "white", flatten = TRUE)
magick::image_write(sheet, OUT("nature_final_v7_contact_sheet.png"),
                    format = "png")
magick::image_write(sheet, OUT("nature_final_v7_contact_sheet.pdf"),
                    format = "pdf")

inv <- do.call(rbind, rows)
write_csv_safe(inv, path_results("reports", "manuscript_candidates",
                                 "nature_final_v7",
                                 "nature_final_v7_preview_inventory.csv"))
cat("\n===== previews and contact sheet =====\n")
print(inv, row.names = FALSE)
cat("\nwritten:", relative_to(OUT()), "\n")
