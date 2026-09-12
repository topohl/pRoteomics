#!/usr/bin/env Rscript

# Candidate comparison contact sheet and panel inventory.
#
# Builds one page of thumbnails covering every declared whole-figure variant,
# plus candidate_panel_inventory.csv. Editorial decisions are NOT encoded here:
# the inventory records a suggested role per panel, taken verbatim from the
# candidate contract, and nothing is promoted or deleted.
#
# USAGE
#   Rscript figures/candidate_figure_contact_sheet.R
#   Rscript figures/candidate_figure_contact_sheet.R --dry-run

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "candidate_figure_utils.R"))
source(repo_path("R", "candidate_figure_panels.R"))
suppressPackageStartupMessages({ library(readr); library(dplyr) })

Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/candidate_figure_contact_sheet.R")
a <- cf_args()
contract <- cf_contract()
shared <- cf_shared_paths()

if (isTRUE(a$dry_run)) {
  cat("[DRY-RUN] Candidate contact sheet and panel inventory.\n")
  dry_run_inputs("figures/candidate_figure_contact_sheet.R",
                 list(candidate_contract = cf_contract_path(),
                      canonical_contract = cf_canonical_contract_path()))
  cat("[DRY-RUN] Reads canonical contract READ-ONLY; writes only under manuscript_candidates.\n")
  quit(save = "no", status = 0L)
}

# ------------------------------------------------- panel inventory

rows <- list()
for (p in contract$panels) {
  fig <- sprintf("%02d", as.integer(p$candidate_figure))
  src <- if (identical(as.character(p$render_mode), "copy_canonical")) {
    paste0("canonical panel ", p$canonical_panel, " (unmodified copy)")
  } else {
    as.character(p$primary_source %||% NA_character_)
  }
  rows[[length(rows) + 1L]] <- data.frame(
    candidate_panel_id = as.character(p$id),
    candidate_figure = fig,
    scientific_question = gsub("\\s+", " ", trimws(as.character(p$scientific_question %||% ""))),
    evidence_type = as.character(p$evidence_type %||% ""),
    inferential_status = as.character(p$inferential_status %||% ""),
    source_artifact = src,
    overlaps_existing_panel = as.character(p$overlaps_existing_panel %||% "none"),
    likely_role = as.character(p$likely_role %||% "undecided"),
    caveat = gsub("\\s+", " ", trimws(as.character(p$caveat %||% ""))),
    render_mode = as.character(p$render_mode),
    stringsAsFactors = FALSE)
}
inventory <- dplyr::bind_rows(rows)
inventory$editorial_decision <- "none_recorded_in_this_task"
inventory$promotion_status <- "not_promoted"

dir_create(shared$tables)
write_csv_safe(inventory, file.path(shared$tables, "candidate_panel_inventory.csv"))

# ------------------------------------------------- variant inventory

vrows <- list()
for (asm in contract$assemblies) {
  fig <- sprintf("%02d", as.integer(asm$candidate_figure))
  svg <- cf_output_paths(fig)$assembled
  svg <- file.path(svg, paste0(as.character(asm$name), ".svg"))
  vrows[[length(vrows) + 1L]] <- data.frame(
    assembly = as.character(asm$name), candidate_figure = fig,
    purpose = gsub("\\s+", " ", trimws(as.character(asm$purpose %||% ""))),
    n_panels = length(asm$layout),
    panels = paste(vapply(asm$layout, function(x) as.character(x$panel), character(1)),
                   collapse = ";"),
    svg = relative_to(svg), exists = file.exists(svg),
    stringsAsFactors = FALSE)
}
variants <- dplyr::bind_rows(vrows)
write_csv_safe(variants, file.path(shared$tables, "candidate_variant_inventory.csv"))

missing <- variants$assembly[!variants$exists]
if (length(missing)) {
  message("[WARN] variant(s) not yet built: ", paste(missing, collapse = ", "),
          " - run figures/candidate_figure_02.R and figures/candidate_figure_03.R first")
}

# ------------------------------------------------- contact sheet

built <- variants[variants$exists, , drop = FALSE]
sheet_svg <- file.path(shared$figures, "candidate_contact_sheet.svg")
if (nrow(built)) {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required for the contact sheet.", call. = FALSE)
  }
  # two rows: Figure 2 variants on top, Figure 3 variants below
  f2 <- built[built$candidate_figure == "02", , drop = FALSE]
  f3 <- built[built$candidate_figure == "03", , drop = FALSE]
  n_col <- max(nrow(f2), nrow(f3), 1L)
  width <- 60 * n_col + 20
  height <- 220
  margin <- 10; gap <- 6; header <- 16
  cell_w <- (width - 2 * margin - (n_col - 1) * gap) / n_col
  cell_h <- (height - 2 * margin - header - gap) / 2

  items <- c(
    sprintf('<text x="%d" y="12" font-family="Arial" font-size="6" font-weight="bold" fill="#23384D">%s</text>',
            margin, "CANDIDATE FIGURE VARIANTS - editorial comparison only, nothing promoted"))
  place <- function(df, row_i) {
    out <- character()
    for (i in seq_len(nrow(df))) {
      x <- margin + (i - 1L) * (cell_w + gap)
      y <- margin + header + (row_i - 1L) * (cell_h + gap)
      uri <- base64enc::dataURI(file = repo_path(df$svg[i]), mime = "image/svg+xml")
      out <- c(out,
        sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="none" stroke="#BBB" stroke-width="0.3"/>',
                x, y + 5, cell_w, cell_h - 5),
        sprintf('<image x="%.2f" y="%.2f" width="%.2f" height="%.2f" preserveAspectRatio="xMidYMid meet" href="%s"/>',
                x, y + 5, cell_w, cell_h - 5, uri),
        sprintf('<text x="%.2f" y="%.2f" font-family="Arial" font-size="3.2" font-weight="bold" fill="#333">%s</text>',
                x, y + 3, df$assembly[i]))
    }
    out
  }
  items <- c(items, place(f2, 1L), place(f3, 2L))
  dir_create(dirname(sheet_svg))
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%smm" height="%smm" viewBox="0 0 %s %s">',
            width, height, width, height),
    '<rect width="100%" height="100%" fill="white"/>', items, '</svg>'),
    sheet_svg, useBytes = TRUE)
  cf_raster_companion(sheet_svg, sub("[.]svg$", ".pdf", sheet_svg))
}

# ------------------------------------------------- report

dir_create(shared$reports)
write_csv_safe(data.frame(
  artifact = c("candidate_panel_inventory", "candidate_variant_inventory",
               "candidate_contact_sheet"),
  path = c(relative_to(file.path(shared$tables, "candidate_panel_inventory.csv")),
           relative_to(file.path(shared$tables, "candidate_variant_inventory.csv")),
           relative_to(sheet_svg)),
  n_rows = c(nrow(inventory), nrow(variants), nrow(built)),
  note = c("one row per candidate panel; likely_role is a suggestion, not a decision",
           "one row per whole-figure variant",
           "thumbnails of every built variant"),
  stringsAsFactors = FALSE),
  file.path(shared$reports, "candidate_contact_sheet_status.csv"))

cat("\n===== Candidate contact sheet =====\n")
cat("panels inventoried :", nrow(inventory), "\n")
cat("variants built     :", nrow(built), "of", nrow(variants), "\n")
cat("contact sheet      :", relative_to(sheet_svg), "\n\n")
print(as.data.frame(table(inventory$likely_role)), row.names = FALSE)
cat("\nNothing is promoted or deleted by this script.\n")
