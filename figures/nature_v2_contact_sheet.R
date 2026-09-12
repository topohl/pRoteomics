#!/usr/bin/env Rscript

# Nature-style v2 contact sheet, at TRUE RELATIVE PANEL SIZE.
#
# Every variant is drawn at the same scale as every other, so the comparison is
# about composition and information density rather than about how large a
# thumbnail happens to be. Nothing is promoted.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "nature_v2_figure_panels.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr) })

Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/nature_v2_contact_sheet.R")
args <- commandArgs(trailingOnly = TRUE)
contract <- nv_contract()
shared <- nv_shared_paths()

if ("--dry-run" %in% args || is_dry_run()) {
  cat("[DRY-RUN] Nature-style v2 contact sheet and panel inventory.\n")
  dry_run_inputs("figures/nature_v2_contact_sheet.R",
                 list(nature_v2_contract = nv_contract_path(),
                      palette = nv_palette_path()))
  quit(save = "no", status = 0L)
}

# -------------------------------------------------------- inventories

prow <- list()
for (p in contract$panels) {
  prow[[length(prow) + 1L]] <- data.frame(
    panel_id = as.character(p$id),
    renderer = as.character(p$renderer),
    role = as.character(p$role %||% ""),
    scientific_question = gsub("\\s+", " ", trimws(as.character(p$question %||% ""))),
    source_artifact = as.character(p$primary_source %||% "derived from canonical metadata"),
    caveat = gsub("\\s+", " ", trimws(as.character(p$caveat %||% ""))),
    stringsAsFactors = FALSE)
}
panel_inventory <- dplyr::bind_rows(prow)
panel_inventory$promotion_status <- "not_promoted"
dir_create(shared$tables)
write_csv_safe(panel_inventory, file.path(shared$tables, "nature_v2_panel_inventory.csv"))

vrow <- list()
for (f in contract$figures) {
  key <- as.character(f$figure_key)
  svg <- file.path(nv_output_paths(key)$assembled, paste0(as.character(f$name), ".svg"))
  used <- vapply(f$layout, function(x) as.character(x$panel), character(1))
  area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
  vrow[[length(vrow) + 1L]] <- data.frame(
    variant = as.character(f$name), figure_key = key,
    width_mm = as.numeric(f$width_mm), height_mm = as.numeric(f$height_mm),
    n_panels = length(used), panels = paste(used, collapse = ";"),
    largest_panel = used[which.max(area)],
    largest_panel_area_share = round(max(area) / sum(area), 3),
    message = gsub("\\s+", " ", trimws(as.character(f$message %||% ""))),
    svg = relative_to(svg), exists = file.exists(svg), stringsAsFactors = FALSE)
}
variant_inventory <- dplyr::bind_rows(vrow)
write_csv_safe(variant_inventory, file.path(shared$tables, "nature_v2_variant_inventory.csv"))

missing <- variant_inventory$variant[!variant_inventory$exists]
if (length(missing)) {
  message("[WARN] not yet built: ", paste(missing, collapse = ", "),
          " - run the nature_v2 entry points first")
}

# ------------------------------------------------------- contact sheet

built <- variant_inventory[variant_inventory$exists, , drop = FALSE]
sheet <- file.path(shared$figures, "nature_v2_contact_sheet.svg")
if (nrow(built)) {
  # true relative size: every variant keeps its real mm footprint, scaled by
  # one shared factor so the whole sheet fits a sensible page
  scale <- 0.42
  gap <- 6; margin <- 8; header <- 10
  n_col <- 3
  cw <- max(built$width_mm) * scale
  ch <- max(built$height_mm) * scale
  n_row <- ceiling(nrow(built) / n_col)
  W <- margin * 2 + n_col * cw + (n_col - 1) * gap
  H <- margin * 2 + header + n_row * (ch + 6) + (n_row - 1) * gap

  items <- sprintf(
    '<text x="%d" y="%d" font-family="Arial" font-size="5" font-weight="bold" fill="#1F3D52">%s</text>',
    margin, 7, "Nature-style v2 candidates - true relative size - editorial comparison only, nothing promoted")
  for (i in seq_len(nrow(built))) {
    r <- (i - 1L) %/% n_col; c <- (i - 1L) %% n_col
    x <- margin + c * (cw + gap)
    y <- margin + header + r * (ch + 6 + gap)
    w <- built$width_mm[i] * scale; h <- built$height_mm[i] * scale
    uri <- base64enc::dataURI(file = repo_path(built$svg[i]), mime = "image/svg+xml")
    items <- c(items,
      sprintf('<text x="%.2f" y="%.2f" font-family="Arial" font-size="3.6" font-weight="bold" fill="#333">%s</text>',
              x, y + 3, built$variant[i]),
      sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="none" stroke="#CCC" stroke-width="0.25"/>',
              x, y + 5, w, h),
      sprintf('<image x="%.2f" y="%.2f" width="%.2f" height="%.2f" href="%s"/>',
              x, y + 5, w, h, uri))
  }
  dir_create(dirname(sheet))
  writeLines(c('<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%.1fmm" height="%.1fmm" viewBox="0 0 %.1f %.1f">',
            W, H, W, H),
    '<rect width="100%" height="100%" fill="white"/>', items, '</svg>'),
    sheet, useBytes = TRUE)
  nv_pdf(sheet, sub("[.]svg$", ".pdf", sheet))
}

dir_create(shared$reports)
write_csv_safe(data.frame(
  artifact = c("nature_v2_panel_inventory", "nature_v2_variant_inventory",
               "nature_v2_contact_sheet"),
  path = c(relative_to(file.path(shared$tables, "nature_v2_panel_inventory.csv")),
           relative_to(file.path(shared$tables, "nature_v2_variant_inventory.csv")),
           relative_to(sheet)),
  n_rows = c(nrow(panel_inventory), nrow(variant_inventory), nrow(built)),
  stringsAsFactors = FALSE),
  file.path(shared$reports, "nature_v2_contact_sheet_status.csv"))

cat("\n===== Nature-style v2 contact sheet =====\n")
cat("panels  :", nrow(panel_inventory), "\n")
cat("variants:", nrow(built), "of", nrow(variant_inventory), "\n")
print(built[, c("variant", "n_panels", "largest_panel", "largest_panel_area_share")],
      row.names = FALSE)
cat("\nsheet:", relative_to(sheet), "\n")
