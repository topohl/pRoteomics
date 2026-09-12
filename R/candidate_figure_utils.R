# =====================================================================
# Candidate manuscript-figure engine.
#
# A PARALLEL layer beside R/manuscript_figure_utils.R. It deliberately does
# not call, wrap, patch or extend the canonical engine, so the whole candidate
# layer can be deleted without touching the canonical pipeline. The only
# contact point is READ-ONLY: copy_canonical panels resolve their asset path
# out of figures/figure_contract.yml so that a "current reference" variant
# always shows the genuine current canonical panel rather than a stale copy.
#
# Candidate renderers are DOWNSTREAM RENDERERS. They join, filter, reshape,
# annotate and visualize validated canonical outputs. They never fit a
# differential-abundance model, refit WGCNA, run GSEA or EWCE, recompute an
# FDR, or generate a new phenotype test. cf_assert_no_model_fitting() enforces
# that by scanning the renderer sources, and a test pins it.
# =====================================================================

cf_contract_path <- function() repo_path("figures", "figure_candidate_contract.yml")
cf_canonical_contract_path <- function() repo_path("figures", "figure_contract.yml")

cf_contract_version <- function() "manuscript_candidate_figures_v1"

# Output namespace. Isolated from results/figures/manuscript/** and from
# results/figures/manuscript_panels/**, which belong to the canonical layer.
cf_output_paths <- function(figure_id, output_root = path_results()) {
  figure_id <- sprintf("%02d", as.integer(figure_id))
  stub <- paste0("figure_", figure_id)
  list(
    panels = file.path(output_root, "figures", "manuscript_candidates", stub, "panels"),
    assembled = file.path(output_root, "figures", "manuscript_candidates", stub, "assembled"),
    source_data = file.path(output_root, "source_data", "manuscript_candidates", stub),
    tables = file.path(output_root, "tables", "manuscript_candidates", stub),
    reports = file.path(output_root, "reports", "manuscript_candidates", stub)
  )
}

cf_shared_paths <- function(output_root = path_results()) {
  list(
    figures = file.path(output_root, "figures", "manuscript_candidates"),
    tables = file.path(output_root, "tables", "manuscript_candidates"),
    reports = file.path(output_root, "reports", "manuscript_candidates")
  )
}

# ------------------------------------------------------------- contract

cf_contract <- function(path = cf_contract_path()) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for the candidate figure contract.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("missing candidate contract: ", path, call. = FALSE)
  }
  y <- yaml::read_yaml(path)
  ids <- vapply(y$panels, function(p) as.character(p$id), character(1))
  if (anyDuplicated(ids)) {
    stop("duplicate candidate panel id(s): ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
  }
  # A candidate id may never look like a canonical manuscript panel id.
  canon_like <- grepl("^[23][a-f]$", ids)
  if (any(canon_like)) {
    stop("candidate panel id(s) collide with the canonical panel namespace: ",
         paste(ids[canon_like], collapse = ", "),
         ". Candidate ids must be distinct from 2a-2f and 3a-3e.", call. = FALSE)
  }
  y
}

cf_panels_for <- function(contract, figure_id) {
  figure_id <- sprintf("%02d", as.integer(figure_id))
  Filter(function(p) identical(sprintf("%02d", as.integer(p$candidate_figure)), figure_id),
         contract$panels)
}

cf_assemblies_for <- function(contract, figure_id) {
  figure_id <- sprintf("%02d", as.integer(figure_id))
  Filter(function(a) identical(sprintf("%02d", as.integer(a$candidate_figure)), figure_id),
         contract$assemblies)
}

cf_panel_by_id <- function(contract, id) {
  hit <- Filter(function(p) identical(as.character(p$id), as.character(id)), contract$panels)
  if (!length(hit)) stop("candidate panel not declared: ", id, call. = FALSE)
  hit[[1]]
}

# Resolve a canonical panel's SVG asset READ-ONLY.
#
# Preference order:
#   1. the panel the CANONICAL ENGINE actually rendered, under
#      results/figures/manuscript/figure_NN/panels/. This is what the current
#      canonical figure really shows, including panels like 3e that the engine
#      renders inline and that therefore have no static figure_source.
#   2. the contract's declared figure_source, for the copy_svg panels.
# Nothing is ever written to either location; the canonical entry points are
# NOT invoked, so canonical outputs cannot change.
cf_canonical_panel_source <- function(canonical_panel_id) {
  id <- as.character(canonical_panel_id)
  fig <- sprintf("%02d", as.integer(substr(id, 1, 1)))
  rendered <- path_results("figures", "manuscript", paste0("figure_", fig),
                           "panels", paste0("figure_", fig, substring(id, 2), ".svg"))
  if (file.exists(rendered)) return(rendered)

  y <- yaml::read_yaml(cf_canonical_contract_path())
  for (f in y$figures) {
    for (p in f$panels) {
      if (identical(as.character(p$id), id)) {
        src <- p$figure_source
        if (is.null(src) || !nzchar(as.character(src))) return(NA_character_)
        return(repo_path(as.character(src)))
      }
    }
  }
  stop("canonical panel not found in the canonical contract: ", id, call. = FALSE)
}

# ------------------------------------------------------- guard rails

# The candidate layer is a renderer. Any of these tokens in a renderer source
# means a new scientific inference is being created here, which is forbidden.
cf_forbidden_tokens <- function() {
  c("lmFit", "eBayes", "limma::", "lmer(", "lme4::", "glm(", "aov(",
    "t.test", "wilcox.test", "cor.test", "p.adjust", "fdrtool", "qvalue",
    "gseGO", "gseKEGG", "GSEA(", "enricher(", "enrichGO", "bootstrap_enrichment_test",
    "blockwiseModules", "TOMsimilarity", "WGCNA::")
}

cf_assert_no_model_fitting <- function(files = cf_renderer_sources()) {
  offending <- character()
  for (f in files) {
    if (!file.exists(f)) next
    lines <- readLines(f, warn = FALSE)
    lines <- sub("#.*$", "", lines)           # code only, not the prose that names them
    # The token list itself names every forbidden call, so its own definition
    # must be cut out before scanning or the guard always trips on itself.
    start <- grep("^cf_forbidden_tokens <- function", lines)
    if (length(start)) {
      close <- grep("^\\}", lines)
      stop_at <- close[close > start[1]][1]
      if (!is.na(stop_at)) lines <- lines[-seq.int(start[1], stop_at)]
    }
    txt <- paste(lines, collapse = "\n")
    for (tok in cf_forbidden_tokens()) {
      if (grepl(tok, txt, fixed = TRUE)) {
        offending <- c(offending, paste0(basename(f), ": ", tok))
      }
    }
  }
  if (length(offending)) {
    stop("candidate renderers must not create new inference. Found: ",
         paste(offending, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

cf_renderer_sources <- function() {
  c(repo_path("R", "candidate_figure_panels.R"),
    repo_path("R", "candidate_figure_utils.R"),
    repo_path("figures", "candidate_figure_02.R"),
    repo_path("figures", "candidate_figure_03.R"),
    repo_path("figures", "candidate_figure_contact_sheet.R"))
}

# ------------------------------------------------------------- helpers

cf_read_csv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("missing_required_input: ", path, call. = FALSE)
    return(NULL)
  }
  as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
}

cf_require_columns <- function(df, cols, label) {
  absent <- setdiff(cols, names(df))
  if (length(absent)) {
    stop("candidate panel ", label, " source is missing required column(s): ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

cf_theme <- function(base_size = 7) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5,
                                            colour = "grey25"),
      legend.key.size = ggplot2::unit(3, "mm"),
      strip.text = ggplot2::element_text(size = base_size - 0.5, face = "bold")
    )
}

cf_group_colours <- function() c(CON = "#3E3C6F", RES = "#9E9A92", SUS = "#D7303F")

# Panels are narrow; an unwrapped title silently runs off the canvas.
cf_wrap <- function(x, width = 95) paste(strwrap(x, width = width), collapse = "\n")

cf_save_panel <- function(plot, svg_path, width_mm, height_mm) {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Package 'svglite' is required to write candidate panels.", call. = FALSE)
  }
  dir_create(dirname(svg_path))
  ggplot2::ggsave(svg_path, plot, width = width_mm / 25.4, height = height_mm / 25.4,
                  units = "in", device = svglite::svglite, bg = "white")
  invisible(svg_path)
}

cf_placeholder_svg <- function(path, panel_id, reason) {
  dir_create(dirname(path))
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<svg xmlns="http://www.w3.org/2000/svg" width="89mm" height="70mm" viewBox="0 0 89 70">',
    '<rect width="100%" height="100%" fill="white" stroke="#B00020" stroke-width="1"/>',
    sprintf('<text x="4" y="12" font-family="Arial" font-size="5" fill="#B00020">candidate panel %s unavailable</text>',
            panel_id),
    sprintf('<text x="4" y="22" font-family="Arial" font-size="3.4" fill="#444">%s</text>',
            substr(gsub("[<>&]", " ", reason), 1, 90)),
    '</svg>'
  ), path, useBytes = TRUE)
  invisible(path)
}

# ------------------------------------------------------------ assembly

# Grid composition into a single self-contained SVG, mirroring the canonical
# layout conventions (base64 data URIs, panel letters) without calling it.
cf_assemble_svg <- function(panel_paths, layout, width_mm, height_mm, target,
                            title = NULL) {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required for candidate SVG assembly.", call. = FALSE)
  }
  n_rows <- max(vapply(layout, function(x) as.integer(x$row), integer(1)))
  n_cols <- max(vapply(layout, function(x) {
    as.integer(x$col) + as.integer(x$colspan %||% 1L) - 1L
  }, integer(1)))
  margin <- 5; gap <- 4; header <- if (is.null(title)) 0 else 7
  cell_w <- (width_mm - 2 * margin - (n_cols - 1) * gap) / n_cols
  cell_h <- (height_mm - 2 * margin - header - (n_rows - 1) * gap) / n_rows
  items <- character()
  if (!is.null(title)) {
    items <- c(items, sprintf(
      '<text x="%.3f" y="%.3f" font-family="Arial" font-size="4.6" font-weight="bold" fill="#23384D">%s</text>',
      margin, margin + 4, title))
  }
  letters_used <- LETTERS[seq_along(layout)]
  for (i in seq_along(layout)) {
    it <- layout[[i]]
    colspan <- as.integer(it$colspan %||% 1L)
    x <- margin + (as.integer(it$col) - 1L) * (cell_w + gap)
    y <- margin + header + (as.integer(it$row) - 1L) * (cell_h + gap)
    w <- cell_w * colspan + gap * (colspan - 1L)
    h <- cell_h
    p <- panel_paths[[as.character(it$panel)]]
    if (is.null(p) || !file.exists(p)) next
    uri <- base64enc::dataURI(file = p, mime = "image/svg+xml")
    items <- c(items,
      sprintf('<image x="%.3f" y="%.3f" width="%.3f" height="%.3f" preserveAspectRatio="xMidYMid meet" href="%s"/>',
              x, y + 4, w, h - 4, uri),
      sprintf('<rect x="%.3f" y="%.3f" width="7" height="7" fill="white" fill-opacity="0.9"/>', x, y),
      sprintf('<text x="%.3f" y="%.3f" fill="black" font-family="Arial" font-size="5" font-weight="bold">%s</text>',
              x + 0.8, y + 5.2, letters_used[i]))
  }
  dir_create(dirname(target))
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%smm" height="%smm" viewBox="0 0 %s %s">',
            width_mm, height_mm, width_mm, height_mm),
    '<rect width="100%" height="100%" fill="white"/>',
    items,
    '</svg>'
  ), target, useBytes = TRUE)
  invisible(target)
}

cf_raster_companion <- function(svg_path, pdf_path) {
  if (!requireNamespace("magick", quietly = TRUE)) return(character())
  out <- tryCatch({
    img <- magick::image_read(svg_path, density = 200)
    img <- magick::image_background(img, "white", flatten = TRUE)
    magick::image_write(img, path = pdf_path, format = "pdf")
    pdf_path
  }, error = function(e) character())
  invisible(out)
}

# ---------------------------------------------------------------- main

cf_build_figure <- function(figure_id, output_root = path_results(),
                            check_only = FALSE, dry_run = FALSE) {
  contract <- cf_contract()
  cf_assert_no_model_fitting()
  panels <- cf_panels_for(contract, figure_id)
  assemblies <- cf_assemblies_for(contract, figure_id)
  paths <- cf_output_paths(figure_id, output_root)

  # -------------------------------------------------- input validation
  rows <- list()
  for (p in panels) {
    srcs <- c(primary = p$primary_source %||% NA_character_,
              stats::setNames(as.character(unlist(p$input_dependencies %||% list())),
                              rep("dependency", length(unlist(p$input_dependencies %||% list())))))
    if (identical(as.character(p$render_mode), "copy_canonical")) {
      srcs <- c(canonical_asset = cf_canonical_panel_source(p$canonical_panel))
    }
    for (k in seq_along(srcs)) {
      v <- srcs[[k]]
      if (is.na(v) || !nzchar(v)) next
      abs <- if (grepl("^(?:[A-Za-z]:|//|\\\\\\\\)", v)) v else repo_path(v)
      rows[[length(rows) + 1L]] <- data.frame(
        panel = as.character(p$id), role = names(srcs)[k],
        input_relative_path = relative_to(abs), exists = file.exists(abs),
        stringsAsFactors = FALSE)
    }
  }
  manifest <- if (length(rows)) dplyr::bind_rows(rows) else
    data.frame(panel = character(), role = character(),
               input_relative_path = character(), exists = logical(),
               stringsAsFactors = FALSE)

  if (dry_run) {
    for (i in seq_len(nrow(manifest))) {
      message("[DRY-RUN ", if (manifest$exists[i]) "PASS" else "WARN", "] panel ",
              manifest$panel[i], " ", manifest$role[i], ": ",
              manifest$input_relative_path[i])
    }
    message("[DRY-RUN] candidate layer only; canonical Figure 2/3 outputs are never written.")
    return(invisible(manifest))
  }
  if (check_only) {
    message("Candidate Figure ", figure_id, " contract check passed for ",
            length(panels), " panel(s) and ", length(assemblies), " assembly variant(s).")
    return(invisible(manifest))
  }

  invisible(lapply(paths, dir_create))

  # -------------------------------------------------------- panels
  panel_paths <- list()
  records <- list()
  for (p in panels) {
    id <- as.character(p$id)
    svg <- file.path(paths$panels, paste0(id, ".svg"))
    src_csv <- file.path(paths$source_data, paste0(id, "_source_data.csv"))
    status <- "ok"; note <- ""
    if (identical(as.character(p$render_mode), "copy_canonical")) {
      asset <- cf_canonical_panel_source(p$canonical_panel)
      if (!is.na(asset) && file.exists(asset)) {
        dir_create(dirname(svg))
        file.copy(asset, svg, overwrite = TRUE, copy.date = TRUE)
        # A reference copy carries a provenance stub, not a re-derived table.
        write_csv_safe(data.frame(
          candidate_panel_id = id, canonical_panel = as.character(p$canonical_panel),
          canonical_asset = relative_to(asset),
          note = "unmodified byte copy of the canonical panel asset; no analysis performed",
          stringsAsFactors = FALSE), src_csv)
      } else {
        status <- "missing_canonical_asset"
        note <- paste0("canonical asset for panel ", p$canonical_panel, " not on disk")
        cf_placeholder_svg(svg, id, note)
        write_csv_safe(data.frame(candidate_panel_id = id, status = status,
                                  note = note, stringsAsFactors = FALSE), src_csv)
      }
    } else {
      fn <- get(as.character(p$renderer), mode = "function")
      res <- tryCatch(fn(p, svg, src_csv), error = function(e) {
        cf_placeholder_svg(svg, id, conditionMessage(e))
        write_csv_safe(data.frame(candidate_panel_id = id, status = "render_error",
                                  note = conditionMessage(e), stringsAsFactors = FALSE),
                       src_csv)
        structure(list(status = "render_error", note = conditionMessage(e)),
                  class = "cf_failed")
      })
      if (inherits(res, "cf_failed")) { status <- res$status; note <- res$note }
    }
    panel_paths[[id]] <- svg
    records[[length(records) + 1L]] <- data.frame(
      candidate_panel_id = id, candidate_figure = sprintf("%02d", as.integer(figure_id)),
      title = as.character(p$title %||% ""),
      render_mode = as.character(p$render_mode),
      status = status, note = note,
      panel_svg = relative_to(svg), source_data = relative_to(src_csv),
      stringsAsFactors = FALSE)
  }
  panel_records <- dplyr::bind_rows(records)

  # ----------------------------------------------------- assemblies
  asm <- list()
  for (a in assemblies) {
    nm <- as.character(a$name)
    target <- file.path(paths$assembled, paste0(nm, ".svg"))
    cf_assemble_svg(panel_paths, a$layout, as.numeric(a$width_mm),
                    as.numeric(a$height_mm), target,
                    title = paste0("CANDIDATE - NOT A MANUSCRIPT FIGURE - ", nm))
    pdf <- sub("[.]svg$", ".pdf", target)
    cf_raster_companion(target, pdf)
    asm[[length(asm) + 1L]] <- data.frame(
      assembly = nm, candidate_figure = sprintf("%02d", as.integer(figure_id)),
      n_panels = length(a$layout),
      panels = paste(vapply(a$layout, function(x) as.character(x$panel), character(1)),
                     collapse = ";"),
      purpose = gsub("\\s+", " ", trimws(as.character(a$purpose %||% ""))),
      svg = relative_to(target),
      pdf = if (file.exists(pdf)) relative_to(pdf) else NA_character_,
      stringsAsFactors = FALSE)
  }
  assembly_records <- dplyr::bind_rows(asm)

  write_csv_safe(panel_records, file.path(paths$reports, "candidate_panel_status.csv"))
  write_csv_safe(assembly_records, file.path(paths$reports, "candidate_assembly_status.csv"))
  write_csv_safe(manifest, file.path(paths$reports, "candidate_input_manifest.csv"))

  invisible(list(panels = panel_records, assemblies = assembly_records,
                 manifest = manifest, paths = paths))
}

# Minimal CLI shared by both candidate entry points.
cf_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  list(
    check_only = "--check-only" %in% args,
    dry_run = "--dry-run" %in% args || is_dry_run()
  )
}
