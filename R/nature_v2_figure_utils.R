# =====================================================================
# Nature-style v2 candidate figure engine.
#
# A THIRD parallel layer. It does not touch the canonical manuscript contract
# and it does not touch the Part-16 candidate layer: different contract,
# different engine, different output root (manuscript_candidates/nature_v2).
#
# THE ONE STRUCTURAL DIFFERENCE FROM THE PART-16 ENGINE
#
# The Part-16 assembler scaled whole panel SVGs into uniform grid cells with
# preserveAspectRatio, so a panel authored at 183 mm that landed in an 85 mm
# cell had all of its text shrunk by more than half. That is why every label
# in the Part-16 PDFs is illegible at print size.
#
# Here each panel declares its EXACT final box in millimetres, is rendered at
# that box, and is placed at 1:1. Scale is never applied, so 6 pt text on the
# canvas really is 6 pt on the page. nv_verify_scale() asserts it.
#
# Renderers remain downstream only: no model fit, no enrichment, no FDR.
# =====================================================================

nv_contract_path <- function() repo_path("figures", "figure_nature_v2_contract.yml")
nv_palette_path <- function() repo_path("config", "manuscript_palette.yml")
nv_contract_version <- function() "manuscript_nature_v2_figures_v1"

nv_output_paths <- function(figure_key, output_root = path_results()) {
  list(
    panels = file.path(output_root, "figures", "manuscript_candidates", "nature_v2",
                       figure_key, "panels"),
    assembled = file.path(output_root, "figures", "manuscript_candidates", "nature_v2",
                          figure_key, "assembled"),
    source_data = file.path(output_root, "source_data", "manuscript_candidates",
                            "nature_v2", figure_key),
    reports = file.path(output_root, "reports", "manuscript_candidates", "nature_v2",
                        figure_key)
  )
}

nv_shared_paths <- function(output_root = path_results()) {
  list(
    figures = file.path(output_root, "figures", "manuscript_candidates", "nature_v2"),
    tables = file.path(output_root, "tables", "manuscript_candidates", "nature_v2"),
    reports = file.path(output_root, "reports", "manuscript_candidates", "nature_v2")
  )
}

# ------------------------------------------------------------- palette

nv_palette <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- yaml::read_yaml(nv_palette_path())
    cache
  }
})

nv_group_colours <- function() unlist(nv_palette()$group)
nv_dataset_colours <- function() unlist(nv_palette()$dataset)
nv_dataset_label <- function(x) {
  l <- unlist(nv_palette()$dataset_label)
  out <- unname(l[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}
nv_claim_colours <- function() unlist(nv_palette()$claimability)
nv_evidence_colours <- function() unlist(nv_palette()$evidence)
nv_pt <- function(key) as.numeric(nv_palette()$typography[[key]])
nv_lw <- function(key) as.numeric(nv_palette()$line[[key]]) * 0.75  # pt -> ggplot linewidth

# ggplot2 sizes are in mm for text; 1 pt = 0.3527 mm
nv_size <- function(pt) pt * 0.3527777

nv_diverging <- function(limits = NULL, name = "NES", ...) {
  d <- nv_palette()$diverging
  ggplot2::scale_fill_gradient2(low = d$low, mid = d$mid, high = d$high,
                                midpoint = 0, limits = limits, name = name, ...)
}

# --------------------------------------------------------------- theme
#
# Nature-style: no plot title, no subtitle, no gridlines by default, hairline
# axes, compact legend. Explanatory text belongs in the figure legend, not
# inside the artwork.
nv_theme <- function(grid = "none", base = nv_pt("axis_text_pt")) {
  fam <- nv_palette()$typography$family
  th <- ggplot2::theme_bw(base_size = base, base_family = fam) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.line = ggplot2::element_line(linewidth = nv_lw("axis_pt"), colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = nv_lw("axis_pt"), colour = "black"),
      axis.ticks.length = ggplot2::unit(0.6, "mm"),
      axis.text = ggplot2::element_text(size = nv_pt("axis_text_pt"), colour = "black"),
      axis.title = ggplot2::element_text(size = nv_pt("axis_title_pt"), colour = "black"),
      legend.title = ggplot2::element_text(size = nv_pt("legend_title_pt")),
      legend.text = ggplot2::element_text(size = nv_pt("legend_text_pt")),
      legend.key.size = ggplot2::unit(2.4, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(1, "mm"),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = nv_pt("axis_text_pt"), colour = "black",
                                         margin = ggplot2::margin(0.4, 0, 0.4, 0, "mm")),
      panel.grid = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.8, "mm")
    )
  if (identical(grid, "y")) {
    th <- th + ggplot2::theme(panel.grid.major.y =
      ggplot2::element_line(linewidth = nv_lw("reference_pt"), colour = "grey92"))
  } else if (identical(grid, "x")) {
    th <- th + ggplot2::theme(panel.grid.major.x =
      ggplot2::element_line(linewidth = nv_lw("reference_pt"), colour = "grey92"))
  }
  th
}

# Heatmap variant: no axis line, tight tiles
nv_theme_tile <- function() {
  nv_theme() + ggplot2::theme(
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank())
}

# ------------------------------------------------------------- guards

# Tokens are CALL forms. `p.adjust` and `p.adjust.method` are also COLUMN names
# in canonical GSEA tables, so matching the bare word would forbid reading the
# stored FDR - which is exactly what a renderer is supposed to do. Requiring the
# opening parenthesis distinguishes "call the function" from "read the column".
nv_forbidden_tokens <- function() {
  c("lmFit(", "eBayes(", "limma::", "lmer(", "lme4::", "glm(", "aov(",
    "t.test(", "wilcox.test(", "cor.test(", "p.adjust(", "fdrtool(", "qvalue(",
    "gseGO(", "gseKEGG(", "GSEA(", "enricher(", "enrichGO(",
    "bootstrap_enrichment_test(", "blockwiseModules(", "TOMsimilarity(", "WGCNA::")
}

nv_renderer_sources <- function() {
  c(repo_path("R", "nature_v2_figure_panels.R"),
    repo_path("R", "nature_v2_figure_utils.R"),
    repo_path("figures", "nature_v2_figure_02.R"),
    repo_path("figures", "nature_v2_figure_03.R"),
    repo_path("figures", "nature_v2_contact_sheet.R"))
}

nv_assert_no_model_fitting <- function(files = nv_renderer_sources()) {
  bad <- character()
  for (f in files) {
    if (!file.exists(f)) next
    lines <- sub("#.*$", "", readLines(f, warn = FALSE))
    start <- grep("^nv_forbidden_tokens <- function", lines)
    if (length(start)) {
      close <- grep("^\\}", lines)
      stop_at <- close[close > start[1]][1]
      if (!is.na(stop_at)) lines <- lines[-seq.int(start[1], stop_at)]
    }
    txt <- paste(lines, collapse = "\n")
    for (tok in nv_forbidden_tokens()) {
      if (grepl(tok, txt, fixed = TRUE)) bad <- c(bad, paste0(basename(f), ": ", tok))
    }
  }
  if (length(bad)) {
    stop("nature_v2 renderers must not create new inference. Found: ",
         paste(bad, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

# ------------------------------------------------------------- contract

nv_contract <- function(path = nv_contract_path()) {
  y <- yaml::read_yaml(path)
  ids <- vapply(y$panels, function(p) as.character(p$id), character(1))
  if (anyDuplicated(ids)) {
    stop("duplicate nature_v2 panel id(s): ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
  }
  if (any(grepl("^[23][a-f]$", ids))) {
    stop("nature_v2 panel id collides with the canonical namespace", call. = FALSE)
  }
  # every figure must fit the Nature size contract
  for (f in y$figures) {
    if (as.numeric(f$width_mm) != 183) {
      stop("figure ", f$name, " is not 183 mm wide", call. = FALSE)
    }
    if (as.numeric(f$height_mm) > 170) {
      stop("figure ", f$name, " exceeds the 170 mm height ceiling", call. = FALSE)
    }
    # boxes must stay inside the canvas
    for (it in f$layout) {
      if (as.numeric(it$x) + as.numeric(it$w) > as.numeric(f$width_mm) + 1e-6 ||
          as.numeric(it$y) + as.numeric(it$h) > as.numeric(f$height_mm) + 1e-6) {
        stop("panel ", it$panel, " overflows figure ", f$name, call. = FALSE)
      }
    }
  }
  y
}

nv_panel_by_id <- function(contract, id) {
  hit <- Filter(function(p) identical(as.character(p$id), as.character(id)),
                contract$panels)
  if (!length(hit)) stop("nature_v2 panel not declared: ", id, call. = FALSE)
  hit[[1]]
}

# ------------------------------------------------------------- helpers

# Some clusterProfiler audit paths exceed the Windows 260-character MAX_PATH,
# so file.exists() reports FALSE and the readers fail on files that are really
# present. Changing directory does not help, because R still resolves the
# relative name against the full working directory. The Windows
# extended-length prefix does work, but only with backslashes and a fully
# qualified path.
nv_long_path <- function(path) {
  if (.Platform$OS.type != "windows") return(path)
  if (nchar(path) < 250L || grepl("^\\\\\\\\[?]", path)) return(path)
  paste0("\\\\?\\", gsub("/", "\\\\", path))
}

nv_read_csv_longpath <- function(dir, file, reader = utils::read.csv, ...) {
  p <- file.path(dir, file)
  if (file.exists(p)) return(reader(p, ...))
  lp <- nv_long_path(p)
  if (!identical(lp, p) && file.exists(lp)) return(reader(lp, ...))
  stop("missing_required_input: ", p, call. = FALSE)
}

nv_read_csv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("missing_required_input: ", path, call. = FALSE)
    return(NULL)
  }
  as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
}

# Render at the EXACT declared box. Nothing downstream rescales it.
nv_save_panel <- function(plot, svg_path, w_mm, h_mm) {
  dir_create(dirname(svg_path))
  ggplot2::ggsave(svg_path, plot, width = w_mm, height = h_mm, units = "mm",
                  device = svglite::svglite, bg = "white")
  invisible(svg_path)
}

nv_placeholder <- function(path, id, reason, w_mm, h_mm) {
  dir_create(dirname(path))
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%smm" height="%smm" viewBox="0 0 %s %s">',
            w_mm, h_mm, w_mm, h_mm),
    '<rect width="100%" height="100%" fill="white"/>',
    sprintf('<rect x="0.3" y="0.3" width="%.2f" height="%.2f" fill="none" stroke="#B9B9B4" stroke-width="0.3" stroke-dasharray="1.5,1.2"/>',
            w_mm - 0.6, h_mm - 0.6),
    sprintf('<text x="%.2f" y="%.2f" font-family="Arial" font-size="6" fill="#555" text-anchor="middle">%s</text>',
            w_mm / 2, h_mm / 2 - 1, substr(reason, 1, 60)),
    sprintf('<text x="%.2f" y="%.2f" font-family="Arial" font-size="4.6" fill="#999" text-anchor="middle">%s</text>',
            w_mm / 2, h_mm / 2 + 4, id),
    '</svg>'), path, useBytes = TRUE)
  invisible(path)
}

# ------------------------------------------------------------ assembly
#
# 1:1 placement. Each panel is embedded at exactly the width and height it was
# authored at, so no text is rescaled. Panel letters are lowercase bold.
nv_assemble <- function(panel_paths, figure, target) {
  W <- as.numeric(figure$width_mm); H <- as.numeric(figure$height_mm)
  lab_pt <- nv_pt("panel_label_pt")
  items <- character()
  for (i in seq_along(figure$layout)) {
    it <- figure$layout[[i]]
    id <- as.character(it$panel)
    x <- as.numeric(it$x); y <- as.numeric(it$y)
    w <- as.numeric(it$w); h <- as.numeric(it$h)
    p <- panel_paths[[sprintf("%s@%gx%g", id, w, h)]]
    if (!is.null(p) && file.exists(p)) {
      uri <- base64enc::dataURI(file = p, mime = "image/svg+xml")
      items <- c(items, sprintf(
        '<image x="%.3f" y="%.3f" width="%.3f" height="%.3f" href="%s"/>', x, y, w, h, uri))
    }
    lab <- as.character(it$label %||% letters[i])
    items <- c(items, sprintf(
      '<text x="%.3f" y="%.3f" font-family="Arial" font-size="%s" font-weight="bold" fill="black">%s</text>',
      x, y + lab_pt * 0.3528 * 0.92, lab_pt, lab))
  }
  dir_create(dirname(target))
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%smm" height="%smm" viewBox="0 0 %s %s">',
            W, H, W, H),
    '<rect width="100%" height="100%" fill="white"/>', items, '</svg>'),
    target, useBytes = TRUE)
  invisible(target)
}

nv_pdf <- function(svg_path, pdf_path) {
  if (!requireNamespace("rsvg", quietly = TRUE)) {
    if (!requireNamespace("magick", quietly = TRUE)) return(NA_character_)
    out <- tryCatch({
      i <- magick::image_read(svg_path, density = 300)
      magick::image_write(magick::image_background(i, "white", flatten = TRUE),
                          path = pdf_path, format = "pdf"); pdf_path
    }, error = function(e) NA_character_)
    return(out)
  }
  out <- tryCatch({ rsvg::rsvg_pdf(svg_path, pdf_path); pdf_path },
                  error = function(e) NA_character_)
  out
}

# A panel authored at its declared box must be embedded at that same box.
nv_verify_scale <- function(figure, panel_paths, tol = 0.01) {
  bad <- character()
  for (it in figure$layout) {
    id <- as.character(it$panel)
    p <- panel_paths[[sprintf("%s@%gx%g", id, as.numeric(it$w), as.numeric(it$h))]]
    if (is.null(p) || !file.exists(p)) {
      stop("no panel rendered for ", id, " at its declared box ",
           it$w, "x", it$h, " mm", call. = FALSE)
    }
    hdr <- paste(readLines(p, n = 4, warn = FALSE), collapse = " ")
    # svglite writes single-quoted attributes; a double-quote-only pattern
    # silently fails to parse and would make this whole check a no-op
    grab <- function(attr) {
      m <- regmatches(hdr, regexpr(paste0(attr, "=['\"]([0-9.]+)pt"), hdr))
      if (!length(m)) return(NA_real_)
      suppressWarnings(as.numeric(sub(paste0(attr, "=['\"]"), "", m)))
    }
    wm <- grab("width"); hm <- grab("height")
    if (!is.finite(wm) || !is.finite(hm)) next        # placeholder svg, mm units
    w_mm <- wm / 72 * 25.4; h_mm <- hm / 72 * 25.4
    if (abs(w_mm - as.numeric(it$w)) > max(tol * as.numeric(it$w), 0.5) ||
        abs(h_mm - as.numeric(it$h)) > max(tol * as.numeric(it$h), 0.5)) {
      bad <- c(bad, sprintf("%s authored %.1fx%.1f mm but placed at %.1fx%.1f mm",
                            id, w_mm, h_mm, as.numeric(it$w), as.numeric(it$h)))
    }
  }
  if (length(bad)) {
    stop("nature_v2 panels would be rescaled, which shrinks text below the ",
         "minimum point size: ", paste(bad, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------- main

nv_build <- function(figure_key, dry_run = FALSE) {
  contract <- nv_contract()
  nv_assert_no_model_fitting()
  figs <- Filter(function(f) identical(as.character(f$figure_key), figure_key),
                 contract$figures)
  if (!length(figs)) stop("no nature_v2 figure with key: ", figure_key, call. = FALSE)
  needed <- unique(unlist(lapply(figs, function(f)
    vapply(f$layout, function(it) as.character(it$panel), character(1)))))
  panels <- Filter(function(p) as.character(p$id) %in% needed, contract$panels)
  paths <- nv_output_paths(figure_key)

  if (dry_run) {
    for (p in panels) {
      src <- as.character(p$primary_source %||% NA_character_)
      ok <- is.na(src) || file.exists(repo_path(src))
      message("[DRY-RUN ", if (ok) "PASS" else "WARN", "] ", p$id, ": ",
              if (is.na(src)) "(no external source)" else src)
    }
    message("[DRY-RUN] nature_v2 candidate layer; canonical and Part-16 outputs untouched.")
    return(invisible(NULL))
  }

  invisible(lapply(paths, dir_create))

  # A panel used by two variants at two different boxes must be rendered TWICE,
  # once per box. Rendering it once and letting the assembler stretch it is
  # exactly the defect that made the Part-16 figures illegible, so panel files
  # are keyed by (id, box) and the assembler looks them up the same way.
  boxes_of <- list()
  for (f in figs) for (it in f$layout) {
    id <- as.character(it$panel)
    b <- sprintf("%gx%g", as.numeric(it$w), as.numeric(it$h))
    boxes_of[[id]] <- unique(c(boxes_of[[id]], b))
  }
  box_key <- function(id, w, h) sprintf("%s@%gx%g", id, w, h)

  panel_paths <- list(); records <- list()
  for (p in panels) {
   id <- as.character(p$id)
   for (bstr in boxes_of[[id]]) {
    box <- as.numeric(strsplit(bstr, "x", fixed = TRUE)[[1]])
    tag <- if (length(boxes_of[[id]]) > 1L) paste0(id, "_", sub("x", "x", bstr)) else id
    svg <- file.path(paths$panels, paste0(tag, ".svg"))
    csv <- file.path(paths$source_data, paste0(id, "_source_data.csv"))
    status <- "ok"; note <- ""
    fn <- tryCatch(get(as.character(p$renderer), mode = "function"),
                   error = function(e) NULL)
    if (is.null(fn)) {
      status <- "renderer_missing"; note <- as.character(p$renderer)
      nv_placeholder(svg, id, note, box[1], box[2])
    } else {
      res <- tryCatch(fn(p, svg, csv, box[1], box[2]), error = function(e) {
        nv_placeholder(svg, id, conditionMessage(e), box[1], box[2])
        write_csv_safe(data.frame(panel = id, status = "render_error",
                                  note = conditionMessage(e), stringsAsFactors = FALSE), csv)
        structure(list(note = conditionMessage(e)), class = "nv_failed")
      })
      if (inherits(res, "nv_failed")) { status <- "render_error"; note <- res$note }
    }
    panel_paths[[box_key(id, box[1], box[2])]] <- svg
    records[[length(records) + 1L]] <- data.frame(
      panel_id = id, figure_key = figure_key,
      box_w_mm = box[1], box_h_mm = box[2],
      role = as.character(p$role %||% ""), status = status, note = note,
      svg = relative_to(svg), source_data = relative_to(csv),
      stringsAsFactors = FALSE)
   }
  }
  panel_records <- dplyr::bind_rows(records)

  asm <- list()
  for (f in figs) {
    nv_verify_scale(f, panel_paths)
    target <- file.path(paths$assembled, paste0(as.character(f$name), ".svg"))
    nv_assemble(panel_paths, f, target)
    pdf <- sub("[.]svg$", ".pdf", target)
    nv_pdf(target, pdf)
    asm[[length(asm) + 1L]] <- data.frame(
      variant = as.character(f$name), figure_key = figure_key,
      width_mm = as.numeric(f$width_mm), height_mm = as.numeric(f$height_mm),
      n_panels = length(f$layout),
      panels = paste(vapply(f$layout, function(x) as.character(x$panel), character(1)),
                     collapse = ";"),
      message = gsub("\\s+", " ", trimws(as.character(f$message %||% ""))),
      svg = relative_to(target),
      pdf = if (file.exists(pdf)) relative_to(pdf) else NA_character_,
      stringsAsFactors = FALSE)
  }
  assembly_records <- dplyr::bind_rows(asm)

  write_csv_safe(panel_records, file.path(paths$reports, "nature_v2_panel_status.csv"))
  write_csv_safe(assembly_records, file.path(paths$reports, "nature_v2_variant_status.csv"))
  invisible(list(panels = panel_records, variants = assembly_records))
}
