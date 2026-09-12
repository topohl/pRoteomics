# =====================================================================
# Story-v3 engine.
#
# A fourth candidate family, isolated from canonical, Part-16 and Part-17.
# It SOURCES R/nature_v2_figure_utils.R read-only to reuse the palette, theme,
# exact-box rendering, long-path reader, assembler and model-fitting guard.
# Those Part-17 files are not modified; only new s4e_* functions are added here.
# =====================================================================

if (!exists("nv_contract", mode = "function")) {
  source(repo_path("R", "nature_v2_figure_utils.R"))
}

s4e_contract_path <- function() repo_path("figures", "figure_story_v4_contract.yml")
s4e_contract_version <- function() "manuscript_story_v4_figures_v1"

s4e_output_paths <- function(figure_key, output_root = path_results()) {
  list(
    panels = file.path(output_root, "figures", "manuscript_candidates", "story_v4",
                       figure_key, "panels"),
    assembled = file.path(output_root, "figures", "manuscript_candidates", "story_v4",
                          figure_key, "assembled"),
    source_data = file.path(output_root, "source_data", "manuscript_candidates",
                            "story_v4", figure_key),
    reports = file.path(output_root, "reports", "manuscript_candidates", "story_v4",
                        figure_key)
  )
}

s4e_shared_paths <- function(output_root = path_results()) {
  list(
    figures = file.path(output_root, "figures", "manuscript_candidates", "story_v4"),
    tables = file.path(output_root, "tables", "manuscript_candidates", "story_v4"),
    reports = file.path(output_root, "reports", "manuscript_candidates", "story_v4")
  )
}

s4e_renderer_sources <- function() {
  c(repo_path("R", "story_v4_figure_panels.R"),
    repo_path("R", "story_v4_figure_utils.R"),
    repo_path("figures", "story_v4_figure_02.R"),
    repo_path("figures", "story_v4_figure_03.R"),
    repo_path("figures", "story_v4_extended_data.R"),
    repo_path("figures", "story_v4_contact_sheet.R"))
}

s4e_contract <- function(path = s4e_contract_path()) {
  y <- yaml::read_yaml(path)
  ids <- vapply(y$panels, function(p) as.character(p$id), character(1))
  if (anyDuplicated(ids)) {
    stop("duplicate story_v4 panel id(s): ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
  }
  if (any(grepl("^[23][a-f]$", ids))) {
    stop("story_v4 panel id collides with the canonical namespace", call. = FALSE)
  }
  for (f in y$figures) {
    if (as.numeric(f$width_mm) != 183) {
      stop("figure ", f$name, " is not 183 mm wide", call. = FALSE)
    }
    if (as.numeric(f$height_mm) > 170) {
      stop("figure ", f$name, " exceeds the 170 mm height ceiling", call. = FALSE)
    }
    for (it in f$layout) {
      if (as.numeric(it$x) + as.numeric(it$w) > as.numeric(f$width_mm) + 1e-6 ||
          as.numeric(it$y) + as.numeric(it$h) > as.numeric(f$height_mm) + 1e-6) {
        stop("panel ", it$panel, " overflows figure ", f$name, call. = FALSE)
      }
    }
  }
  y
}

s4e_panel_by_id <- function(contract, id) {
  hit <- Filter(function(p) identical(as.character(p$id), as.character(id)),
                contract$panels)
  if (!length(hit)) stop("story_v4 panel not declared: ", id, call. = FALSE)
  hit[[1]]
}

s4e_build <- function(figure_key, dry_run = FALSE) {
  contract <- s4e_contract()
  nv_assert_no_model_fitting(s4e_renderer_sources())
  figs <- Filter(function(f) identical(as.character(f$figure_key), figure_key),
                 contract$figures)
  if (!length(figs)) stop("no story_v4 figure with key: ", figure_key, call. = FALSE)
  needed <- unique(unlist(lapply(figs, function(f)
    vapply(f$layout, function(it) as.character(it$panel), character(1)))))
  panels <- Filter(function(p) as.character(p$id) %in% needed, contract$panels)
  paths <- s4e_output_paths(figure_key)

  if (dry_run) {
    for (p in panels) {
      src <- as.character(p$primary_source %||% NA_character_)
      ok <- is.na(src) || file.exists(repo_path(src))
      message("[DRY-RUN ", if (ok) "PASS" else "WARN", "] ", p$id, ": ",
              if (is.na(src)) "(derived)" else src)
    }
    message("[DRY-RUN] story_v4 candidate layer; canonical, Part-16 and Part-17 untouched.")
    return(invisible(NULL))
  }

  invisible(lapply(paths, dir_create))
  boxes_of <- list()
  for (f in figs) for (it in f$layout) {
    id <- as.character(it$panel)
    boxes_of[[id]] <- unique(c(boxes_of[[id]],
                               sprintf("%gx%g", as.numeric(it$w), as.numeric(it$h))))
  }

  panel_paths <- list(); records <- list()
  for (p in panels) {
    id <- as.character(p$id)
    for (bstr in boxes_of[[id]]) {
      box <- as.numeric(strsplit(bstr, "x", fixed = TRUE)[[1]])
      tag <- if (length(boxes_of[[id]]) > 1L) paste0(id, "_", bstr) else id
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
                                    note = conditionMessage(e),
                                    stringsAsFactors = FALSE), csv)
          structure(list(note = conditionMessage(e)), class = "nv_failed")
        })
        if (inherits(res, "nv_failed")) { status <- "render_error"; note <- res$note }
      }
      panel_paths[[sprintf("%s@%gx%g", id, box[1], box[2])]] <- svg
      records[[length(records) + 1L]] <- data.frame(
        panel_id = id, figure_key = figure_key,
        box_w_mm = box[1], box_h_mm = box[2],
        role = as.character(p$role %||% ""),
        narrative = gsub("\\s+", " ", trimws(as.character(p$narrative %||% ""))),
        status = status, note = note,
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
    area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
    ids <- vapply(f$layout, function(x) as.character(x$panel), character(1))
    asm[[length(asm) + 1L]] <- data.frame(
      variant = as.character(f$name), figure_key = figure_key,
      width_mm = as.numeric(f$width_mm), height_mm = as.numeric(f$height_mm),
      n_panels = length(ids), panels = paste(ids, collapse = ";"),
      largest_panel = ids[which.max(area)],
      largest_area_share = round(max(area) / sum(area), 3),
      question = gsub("\\s+", " ", trimws(as.character(f$question %||% ""))),
      svg = relative_to(target),
      pdf = if (file.exists(pdf)) relative_to(pdf) else NA_character_,
      stringsAsFactors = FALSE)
  }
  assembly_records <- dplyr::bind_rows(asm)

  write_csv_safe(panel_records, file.path(paths$reports, "story_v4_panel_status.csv"))
  write_csv_safe(assembly_records, file.path(paths$reports, "story_v4_variant_status.csv"))
  invisible(list(panels = panel_records, variants = assembly_records))
}
