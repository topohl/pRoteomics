#!/usr/bin/env Rscript

# Story-v3 contact sheet + figure_story_map.md.
#
# The story map is generated FROM the contract's per-panel narrative lines, so
# the argument recorded in the contract and the argument shown to a reviewer
# cannot drift apart. If a panel has no narrative line, the map says so rather
# than inventing a transition.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "story_v3_figure_utils.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr) })

Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/story_v3_contact_sheet.R")
args <- commandArgs(trailingOnly = TRUE)
contract <- sv_contract()
shared <- sv_shared_paths()

if ("--dry-run" %in% args || is_dry_run()) {
  cat("[DRY-RUN] Story-v3 contact sheet, inventory and story map.\n")
  dry_run_inputs("figures/story_v3_contact_sheet.R",
                 list(story_v3_contract = sv_contract_path()))
  quit(save = "no", status = 0L)
}

narr <- vapply(contract$panels, function(p)
  gsub("\\s+", " ", trimws(as.character(p$narrative %||% ""))), character(1))
names(narr) <- vapply(contract$panels, function(p) as.character(p$id), character(1))
role <- vapply(contract$panels, function(p) as.character(p$role %||% ""), character(1))
names(role) <- names(narr)

# ------------------------------------------------------ story map

lines <- c("# Story-v3 figure story maps", "",
  "Generated from the `narrative` field of every panel in",
  "`figures/figure_story_v3_contract.yml`. Each sentence states what the panel",
  "establishes; the figure is coherent only if the next panel follows from it",
  "without an \"and separately we also found\".", "")

for (f in contract$figures) {
  lines <- c(lines, sprintf("## %s", as.character(f$name)), "",
             sprintf("**Question:** %s",
                     gsub("\\s+", " ", trimws(as.character(f$question %||% "")))), "")
  lay <- f$layout
  ord <- order(vapply(lay, function(x) as.character(x$label), character(1)))
  for (i in ord) {
    it <- lay[[i]]
    id <- as.character(it$panel); lb <- as.character(it$label)
    n <- narr[[id]]
    lines <- c(lines, sprintf("- **Panel %s** (`%s`, %s) %s", lb, id, role[[id]],
                              if (nzchar(n)) n else
                                "HAS NO NARRATIVE LINE - it does not yet earn its place."))
  }
  lines <- c(lines, "")
}

dir_create(shared$reports)
writeLines(lines, file.path(shared$reports, "figure_story_map.md"))

# ------------------------------------------------------ inventories

prow <- list()
for (p in contract$panels) {
  prow[[length(prow) + 1L]] <- data.frame(
    panel_id = as.character(p$id), renderer = as.character(p$renderer),
    role = as.character(p$role %||% ""),
    narrative = gsub("\\s+", " ", trimws(as.character(p$narrative %||% ""))),
    source_artifact = as.character(p$primary_source %||% "derived"),
    promotion_status = "not_promoted", stringsAsFactors = FALSE)
}
panel_inventory <- dplyr::bind_rows(prow)
dir_create(shared$tables)
write_csv_safe(panel_inventory, file.path(shared$tables, "story_v3_panel_inventory.csv"))

vrow <- list()
for (f in contract$figures) {
  key <- as.character(f$figure_key)
  svg <- file.path(sv_output_paths(key)$assembled, paste0(as.character(f$name), ".svg"))
  ids <- vapply(f$layout, function(x) as.character(x$panel), character(1))
  area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
  vrow[[length(vrow) + 1L]] <- data.frame(
    variant = as.character(f$name), figure_key = key,
    width_mm = as.numeric(f$width_mm), height_mm = as.numeric(f$height_mm),
    n_panels = length(ids), panels = paste(ids, collapse = ";"),
    largest_panel = ids[which.max(area)],
    largest_area_share = round(max(area) / sum(area), 3),
    used_area_share = round(sum(area) /
      (as.numeric(f$width_mm) * as.numeric(f$height_mm)), 3),
    question = gsub("\\s+", " ", trimws(as.character(f$question %||% ""))),
    svg = relative_to(svg), exists = file.exists(svg), stringsAsFactors = FALSE)
}
variant_inventory <- dplyr::bind_rows(vrow)
write_csv_safe(variant_inventory, file.path(shared$tables, "story_v3_variant_inventory.csv"))

# ------------------------------------------------------ contact sheet

built <- variant_inventory[variant_inventory$exists, , drop = FALSE]
sheet <- file.path(shared$figures, "story_v3_contact_sheet.svg")
if (nrow(built)) {
  scale <- 0.42; gap <- 6; margin <- 8; header <- 10; n_col <- 3
  cw <- max(built$width_mm) * scale; ch <- max(built$height_mm) * scale
  n_row <- ceiling(nrow(built) / n_col)
  W <- margin * 2 + n_col * cw + (n_col - 1) * gap
  H <- margin * 2 + header + n_row * (ch + 6) + (n_row - 1) * gap
  items <- sprintf(
    '<text x="%d" y="%d" font-family="Arial" font-size="5" font-weight="bold" fill="#1F3D52">%s</text>',
    margin, 7, "Story-v3 candidates - true relative size - editorial comparison only, nothing promoted")
  for (i in seq_len(nrow(built))) {
    r <- (i - 1L) %/% n_col; c <- (i - 1L) %% n_col
    x <- margin + c * (cw + gap); y <- margin + header + r * (ch + 6 + gap)
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

cat("\n===== Story-v3 contact sheet =====\n")
cat("panels  :", nrow(panel_inventory),
    "  with a narrative line:", sum(nzchar(panel_inventory$narrative)), "\n")
cat("variants:", nrow(built), "of", nrow(variant_inventory), "\n")
print(built[, c("variant", "n_panels", "largest_panel", "largest_area_share",
                "used_area_share")], row.names = FALSE)
cat("\nstory map:", relative_to(file.path(shared$reports, "figure_story_map.md")), "\n")
