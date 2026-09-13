# =====================================================================
# Story-v3 engine.
#
# A fourth candidate family, isolated from canonical, Part-16 and Part-17.
# It SOURCES R/nature_v2_figure_utils.R read-only to reuse the palette, theme,
# exact-box rendering, long-path reader, assembler and model-fitting guard.
# Those Part-17 files are not modified; only new s7e_* functions are added here.
# =====================================================================

if (!exists("nv_contract", mode = "function")) {
  source(repo_path("R", "nature_v2_figure_utils.R"))
}

s7e_contract_path <- function() repo_path("figures", "figure_nature_final_v7_contract.yml")
s7e_contract_version <- function() "manuscript_nature_final_v7_figures_v1"

s7e_output_paths <- function(figure_key, output_root = path_results()) {
  list(
    panels = file.path(output_root, "figures", "manuscript_candidates", "nature_final_v7",
                       figure_key, "panels"),
    assembled = file.path(output_root, "figures", "manuscript_candidates", "nature_final_v7",
                          figure_key, "assembled"),
    source_data = file.path(output_root, "source_data", "manuscript_candidates",
                            "nature_final_v7", figure_key),
    reports = file.path(output_root, "reports", "manuscript_candidates", "nature_final_v7",
                        figure_key)
  )
}

s7e_shared_paths <- function(output_root = path_results()) {
  list(
    figures = file.path(output_root, "figures", "manuscript_candidates", "nature_final_v7"),
    tables = file.path(output_root, "tables", "manuscript_candidates", "nature_final_v7"),
    reports = file.path(output_root, "reports", "manuscript_candidates", "nature_final_v7")
  )
}

s7e_renderer_sources <- function() {
  c(repo_path("R", "nature_final_v7_panels.R"),
    repo_path("R", "nature_final_v7_figure3_panels.R"),
    repo_path("R", "nature_final_v7_ed_panels.R"),
    repo_path("R", "nature_final_v7_figure_utils.R"),
    repo_path("R", "spatial_grammar_utils.R"),
    repo_path("figures", "nature_final_v7_figure_02.R"),
    repo_path("figures", "nature_final_v7_figure_03.R"),
    repo_path("figures", "nature_final_v7_extended_data.R"))
}

s7e_contract <- function(path = s7e_contract_path()) {
  y <- yaml::read_yaml(path)
  ids <- vapply(y$panels, function(p) as.character(p$id), character(1))
  if (anyDuplicated(ids)) {
    stop("duplicate nature_final_v7 panel id(s): ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
  }
  if (any(grepl("^[23][a-f]$", ids))) {
    stop("nature_final_v7 panel id collides with the canonical namespace", call. = FALSE)
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

s7e_panel_by_id <- function(contract, id) {
  hit <- Filter(function(p) identical(as.character(p$id), as.character(id)),
                contract$panels)
  if (!length(hit)) stop("nature_final_v7 panel not declared: ", id, call. = FALSE)
  hit[[1]]
}

s7e_build <- function(figure_key, dry_run = FALSE) {
  contract <- s7e_contract()
  nv_assert_no_model_fitting(s7e_renderer_sources())
  figs <- Filter(function(f) identical(as.character(f$figure_key), figure_key),
                 contract$figures)
  if (!length(figs)) stop("no nature_final_v7 figure with key: ", figure_key, call. = FALSE)
  needed <- unique(unlist(lapply(figs, function(f)
    vapply(f$layout, function(it) as.character(it$panel), character(1)))))
  panels <- Filter(function(p) as.character(p$id) %in% needed, contract$panels)
  paths <- s7e_output_paths(figure_key)

  if (dry_run) {
    for (p in panels) {
      src <- as.character(p$primary_source %||% NA_character_)
      ok <- is.na(src) || file.exists(repo_path(src))
      message("[DRY-RUN ", if (ok) "PASS" else "WARN", "] ", p$id, ": ",
              if (is.na(src)) "(derived)" else src)
    }
    message("[DRY-RUN] nature_final_v7 candidate layer; canonical, Part-16 and Part-17 untouched.")
    return(invisible(NULL))
  }

  invisible(lapply(paths, dir_create))
  boxes_of <- list()
  for (f in figs) for (it in f$layout) {
    id <- as.character(it$panel)
    boxes_of[[id]] <- unique(c(boxes_of[[id]],
                               sprintf("%gx%g", as.numeric(it$w), as.numeric(it$h))))
  }


  # Purge panel SVGs left over from a previous contract revision. Without this a
  # renamed box silently persists on disk and pollutes both the layout audit and
  # any byte-for-byte determinism check.
  expected_svg <- character(0)
  for (p in panels) for (bstr in boxes_of[[as.character(p$id)]]) {
    tag <- if (length(boxes_of[[as.character(p$id)]]) > 1L)
      paste0(as.character(p$id), "_", bstr) else as.character(p$id)
    expected_svg <- c(expected_svg, paste0(tag, ".svg"))
  }
  stale <- setdiff(list.files(paths$panels, "[.]svg$"), expected_svg)
  if (length(stale)) {
    message("[nature_final_v7] removing ", length(stale), " stale panel file(s): ",
            paste(stale, collapse = ", "))
    file.remove(file.path(paths$panels, stale))
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

  write_csv_safe(panel_records, file.path(paths$reports, "nature_final_v7_panel_status.csv"))
  write_csv_safe(assembly_records, file.path(paths$reports, "nature_final_v7_variant_status.csv"))
  qa <- s7e_qa(panel_records, contract, figure_key)
  qa_path <- file.path(paths$reports, "nature_final_v7_qa.csv")
  write_csv_safe(data.frame(
    figure_key = figure_key,
    check = if (length(qa)) qa else "all hard QA conditions passed",
    status = if (length(qa)) "FAIL" else "PASS", stringsAsFactors = FALSE),
    qa_path)
  if (length(qa)) {
    stop("nature_final_v7 QA FAILED:
  ", paste(qa, collapse = "
  "),
         call. = FALSE)
  }
  message("[nature_final_v7] QA passed for ", figure_key)
  invisible(list(panels = panel_records, variants = assembly_records, qa = qa))
}

# --------------------------------------------------------------------- QA
#
# Hard fail conditions from the brief. These run against the EMITTED SVG, not
# against the contract, so a renderer that silently drops below the floor is
# caught rather than trusted.
s7e_font_sizes <- function(svg) {
  if (!file.exists(svg)) return(numeric(0))
  txt <- paste(readLines(svg, warn = FALSE), collapse = "")
  m <- regmatches(txt, gregexpr("font-size: [0-9.]+px", txt))[[1]]
  as.numeric(sub("px", "", sub("font-size: ", "", m)))
}

s7e_qa <- function(records, contract, figure_key) {
  fails <- character(0)
  floor_pt <- as.numeric(contract$font_floor_pt %||% 5)

  # 1. every panel rendered
  bad <- records$panel_id[records$status != "ok"]
  if (length(bad)) fails <- c(fails, paste0("panel did not render: ",
                                            paste(unique(bad), collapse = ", ")))

  # 2. FONT FLOOR, measured on the emitted SVG
  for (i in seq_len(nrow(records))) {
    fs <- s7e_font_sizes(repo_path(records$svg[i]))
    lo <- fs[fs < floor_pt - 1e-6]
    if (length(lo)) {
      fails <- c(fails, sprintf("%s: %d text element(s) below %.1f pt (min %.2f)",
                                records$panel_id[i], length(lo), floor_pt, min(lo)))
    }
  }

  # 3. page size and overlap
  figs <- Filter(function(f) identical(as.character(f$figure_key), figure_key),
                 contract$figures)
  for (f in figs) {
    W <- as.numeric(f$width_mm); H <- as.numeric(f$height_mm)
    if (W != 183) fails <- c(fails, paste0(f$name, ": not 183 mm wide"))
    if (H > 170) fails <- c(fails, paste0(f$name, ": exceeds 170 mm"))
    n <- length(f$layout)
    for (i in seq_len(n)) {
      a <- f$layout[[i]]
      if (as.numeric(a$x) + as.numeric(a$w) > W + 1e-6 ||
          as.numeric(a$y) + as.numeric(a$h) > H + 1e-6) {
        fails <- c(fails, sprintf("%s: panel %s overflows the page", f$name, a$panel))
      }
      if (i < n) for (j in seq(i + 1L, n)) {
        b <- f$layout[[j]]
        sep <- as.numeric(a$x) + as.numeric(a$w) <= as.numeric(b$x) ||
          as.numeric(b$x) + as.numeric(b$w) <= as.numeric(a$x) ||
          as.numeric(a$y) + as.numeric(a$h) <= as.numeric(b$y) ||
          as.numeric(b$y) + as.numeric(b$h) <= as.numeric(a$y)
        if (!sep) fails <- c(fails, sprintf("%s: %s overlaps %s", f$name,
                                            a$panel, b$panel))
      }
    }
    # 4. a central scientific panel must not be smaller than a minor technical one
    area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h),
                   numeric(1))
    ids <- vapply(f$layout, function(x) as.character(x$panel), character(1))
    minor <- c("v7_depth", "v7_pca")
    central <- c("v7_fingerprint", "v7_compartment", "v7_bilateral_main",
                 "v7_atlas", "v7_prot_syn", "v7_prot_rna", "v7_prot_ox")
    if (any(ids %in% minor) && any(ids %in% central)) {
      if (max(area[ids %in% minor]) > min(area[ids %in% central])) {
        fails <- c(fails, sprintf(
          "%s: a minor technical panel is larger than a central scientific panel",
          f$name))
      }
    }
  }
  fails
}
