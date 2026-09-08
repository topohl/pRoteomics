# Utilities for explicit manuscript Figure 2/3 entry points.
#
# This layer consumes frozen/canonical downstream artifacts. It does not fit
# differential-abundance, enrichment, or WGCNA models.

if (!exists("output_namespace_manuscript_figure_paths", mode = "function")) {
  source(repo_path("R", "output_namespace_utils.R"))
}

manuscript_figure_contract_path <- function() {
  repo_path("figures", "figure_contract.yml")
}

manuscript_figure_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  value_for <- function(flag) {
    exact <- which(args == flag)
    if (length(exact)) {
      if (exact[[1]] == length(args)) {
        stop("Missing value after ", flag, ".", call. = FALSE)
      }
      return(args[[exact[[1]] + 1L]])
    }
    prefixed <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (length(prefixed)) return(sub(paste0("^", flag, "="), "", prefixed[[1]]))
    NULL
  }

  output_root <- value_for("--output-root")
  panel <- value_for("--panel")
  dataset <- value_for("--dataset")
  if (!is.null(dataset) && !identical(dataset, "global")) {
    stop("Manuscript figure entry points support only --dataset global.", call. = FALSE)
  }

  consumed <- c(
    "--check-only", "--dry-run", "--allow-incomplete",
    grep("^--output-root=", args, value = TRUE),
    grep("^--panel=", args, value = TRUE),
    grep("^--dataset=", args, value = TRUE)
  )
  for (flag in c("--output-root", "--panel", "--dataset")) {
    at <- which(args == flag)
    if (length(at)) consumed <- c(consumed, flag, args[at[[1]] + 1L])
  }
  unknown <- setdiff(args, consumed)
  if (length(unknown)) {
    stop("Unknown manuscript-figure argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  output_explicit <- !is.null(output_root)
  if (is.null(output_root)) {
    output_root <- path_results()
  } else if (!grepl("^(?:[A-Za-z]:[/\\\\]|//|\\\\\\\\)", output_root, perl = TRUE)) {
    output_root <- repo_path(output_root)
  }
  output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)

  list(
    check_only = "--check-only" %in% args,
    dry_run = "--dry-run" %in% args || is_dry_run(),
    allow_incomplete = "--allow-incomplete" %in% args,
    output_root = output_root,
    output_explicit = output_explicit,
    panel = if (is.null(panel)) NULL else tolower(panel)
  )
}

manuscript_figure_contract <- function(figure_id) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for manuscript figure contracts.", call. = FALSE)
  }
  contract_path <- manuscript_figure_contract_path()
  if (!file.exists(contract_path)) stop("Missing figure contract: ", contract_path, call. = FALSE)
  contract <- yaml::read_yaml(contract_path)
  figure <- contract$figures[[figure_id]]
  if (is.null(figure)) stop("Unknown manuscript figure: ", figure_id, call. = FALSE)
  panel_ids <- vapply(figure$panels, function(x) as.character(x$id), character(1))
  if (anyDuplicated(panel_ids)) stop("Duplicate panel IDs in Figure ", figure_id, ".", call. = FALSE)
  figure$contract_version <- as.character(contract$contract_version)
  figure
}

manuscript_figure_panel_is_automated <- function(panel) {
  include <- panel$include_in_automated_assembly
  if (is.null(include)) return(TRUE)
  isTRUE(include)
}

manuscript_figure_resolve <- function(path) {
  path <- as.character(path %||% NA_character_)
  if (!length(path) || is.na(path[[1]]) || !nzchar(path[[1]])) return(NA_character_)
  normalizePath(repo_path(path[[1]]), winslash = "/", mustWork = FALSE)
}

manuscript_figure_input_rows <- function(panel) {
  paths <- c(
    figure_source = as.character(panel$figure_source %||% character()),
    primary_source = as.character(panel$primary_source %||% character()),
    input_dependency = as.character(unlist(panel$input_dependencies %||% character(), use.names = FALSE))
  )
  roles <- names(paths)
  roles[roles == ""] <- "input_dependency"
  paths <- unname(paths)
  keep <- !is.na(paths) & nzchar(paths)
  paths <- paths[keep]
  roles <- roles[keep]
  resolved <- vapply(paths, manuscript_figure_resolve, character(1))
  exists <- file.exists(resolved)
  size <- rep(NA_real_, length(resolved))
  mtime <- rep(NA_character_, length(resolved))
  size[exists] <- as.numeric(file.info(resolved[exists])$size)
  mtime[exists] <- format(file.info(resolved[exists])$mtime, "%Y-%m-%d %H:%M:%S %z")
  sha <- vapply(resolved, file_hash_sha256, character(1))
  data.frame(
    panel = as.character(panel$id),
    role = roles,
    input_relative_path = paths,
    input_resolved_path = resolved,
    exists = exists,
    size_bytes = size,
    mtime = mtime,
    sha256 = sha,
    stringsAsFactors = FALSE
  )
}

manuscript_figure_primary_header <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("csv", "tsv")) return(character())
  sep <- if (ext == "tsv") "\t" else ","
  names(utils::read.table(
    path, header = TRUE, sep = sep, nrows = 0L, check.names = FALSE,
    quote = "\"", comment.char = "", stringsAsFactors = FALSE
  ))
}

manuscript_figure_validate_panel <- function(panel) {
  inputs <- manuscript_figure_input_rows(panel)
  missing <- inputs[!inputs$exists, , drop = FALSE]

  primary <- manuscript_figure_resolve(panel$primary_source)
  required_columns <- as.character(unlist(panel$required_columns %||% character(), use.names = FALSE))
  if (file.exists(primary) && length(required_columns)) {
    header <- manuscript_figure_primary_header(primary)
    absent <- setdiff(required_columns, header)
    if (length(absent)) {
      stop(
        "Panel ", panel$id, " primary source is missing required column(s): ",
        paste(absent, collapse = ", "), call. = FALSE
      )
    }
  }

  expected_rows <- suppressWarnings(as.integer(panel$expected_rows %||% NA_integer_))
  if (file.exists(primary) && is.finite(expected_rows)) {
    observed <- nrow(utils::read.csv(primary, check.names = FALSE))
    if (!identical(observed, expected_rows)) {
      stop(
        "Panel ", panel$id, " expected ", expected_rows,
        " rows but observed ", observed, ".", call. = FALSE
      )
    }
  }
  list(inputs = inputs, missing = missing)
}

manuscript_figure_output_paths <- function(output_root, figure_id) {
  output_namespace_manuscript_figure_paths(output_root, figure_id)
}

manuscript_figure_copy <- function(source, target) {
  dir_create(dirname(target))
  if (!file.copy(source, target, overwrite = TRUE, copy.date = TRUE)) {
    stop("Could not copy manuscript figure artifact: ", source, " -> ", target, call. = FALSE)
  }
  invisible(target)
}

manuscript_figure_derive_2b_source <- function(workbook, target) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required.", call. = FALSE)
  if (!requireNamespace("readr", quietly = TRUE)) stop("Package 'readr' is required.", call. = FALSE)
  x <- readxl::read_excel(workbook)
  needed <- c("sample_id", "Precursors.Identified", "Proteins.Identified")
  absent <- setdiff(needed, names(x))
  if (length(absent)) stop("Figure 2b workbook is missing: ", paste(absent, collapse = ", "), call. = FALSE)
  sample_id <- as.character(x$sample_id)
  celltype_layer <- if ("celltype_layer" %in% names(x)) as.character(x$celltype_layer) else rep(NA_character_, nrow(x))
  keep <- !grepl("background|blank|bg", sample_id, ignore.case = TRUE) &
    !grepl("background|blank|bg", celltype_layer, ignore.case = TRUE)
  x <- x[keep, , drop = FALSE]

  candidates <- intersect(c("celltype_layer", "CellTypeLayer", "celltype", "CellType"), names(x))
  compartment <- rep(NA_character_, nrow(x))
  for (column in candidates) {
    key <- tolower(trimws(as.character(x[[column]])))
    key <- gsub("[ -]+", "_", key)
    mapped <- ifelse(
      key %in% c("microglia", "microglial"), "microglia",
      ifelse(key %in% c("neuron_neuropil", "neuropil"), "neuropil",
             ifelse(key %in% c("neuron_soma", "soma"), "soma", NA_character_))
    )
    fill <- is.na(compartment) & !is.na(mapped)
    compartment[fill] <- mapped[fill]
  }
  if (any(is.na(compartment))) {
    bad <- sort(unique(celltype_layer[is.na(compartment)]))
    stop("Figure 2b has unmapped compartment rows: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  x$qc_compartment <- compartment
  x$proteins_plot_order <- ave(
    as.numeric(x$Proteins.Identified), x$qc_compartment,
    FUN = function(v) rank(v, ties.method = "first", na.last = "keep")
  )
  x$precursors_plot_order <- ave(
    as.numeric(x$Precursors.Identified), x$qc_compartment,
    FUN = function(v) rank(v, ties.method = "first", na.last = "keep")
  )
  keep_columns <- intersect(c(
    "sample_id", "File.Name", "AnimalID", "celltype", "celltype_layer",
    "group", "group2", "region", "layer", "exclude", "qc_compartment",
    "Precursors.Identified", "Proteins.Identified", "proteins_plot_order",
    "precursors_plot_order"
  ), names(x))
  out <- as.data.frame(x[, keep_columns, drop = FALSE], stringsAsFactors = FALSE)
  out <- out[order(match(out$qc_compartment, c("microglia", "neuropil", "soma")), out$sample_id), , drop = FALSE]
  dir_create(dirname(target))
  readr::write_csv(out, target, na = "")
  invisible(target)
}

manuscript_figure_render_m12 <- function(source, panel_svg, source_target) {
  required_packages <- c("readr", "dplyr", "ggplot2", "svglite", "scales")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
  if (!exists("nature_palette", mode = "function")) source(repo_path("R", "plotting_nature.R"))
  if (!exists("anatomical_spatial_unit_levels", mode = "function")) {
    source(repo_path("R", "sus_res_spatial_dap_atlas_utils.R"))
  }

  all_rows <- readr::read_csv(source, show_col_types = FALSE, progress = FALSE)
  block <- all_rows[all_rows$ModuleID == "WGCNA_m12", , drop = FALSE]
  if (!nrow(block)) stop("Figure 3e source has no WGCNA_m12 rows.", call. = FALSE)
  if (length(unique(block$ProteinGroupID)) != 15L) {
    stop("Figure 3e requires exactly 15 unique WGCNA_m12 proteins.", call. = FALSE)
  }
  contrast_levels <- c("RES - CON", "SUS - CON", "SUS - RES")
  if (!setequal(unique(block$contrast), contrast_levels)) {
    stop("Figure 3e requires all three canonical contrasts.", call. = FALSE)
  }
  block$spatial_unit <- factor(
    block$spatial_unit,
    levels = anatomical_spatial_unit_levels(unique(block$spatial_unit))
  )
  block$contrast_plot_label <- factor(
    block$contrast, levels = contrast_levels,
    labels = c("RES\n- CON", "SUS\n- CON", "SUS\n- RES")
  )
  protein_order <- unique(block[order(block$display_rank), c("ProteinGroupID", "protein_label")])
  block$protein_plot_key <- factor(
    as.character(block$ProteinGroupID),
    levels = rev(as.character(protein_order$ProteinGroupID))
  )
  protein_labels <- stats::setNames(
    as.character(protein_order$protein_label), as.character(protein_order$ProteinGroupID)
  )
  effect_limit <- max(abs(all_rows$log2FC), na.rm = TRUE)
  manuscript_text <- nature_manuscript_text_sizes_pt()
  compact_label <- function(x) gsub(" ", "-", clean_spatial_unit_label(x), fixed = TRUE)
  supported <- block[block$fdr_supported %in% TRUE, , drop = FALSE]
  p <- ggplot2::ggplot(
    block,
    ggplot2::aes(.data$spatial_unit, .data$protein_plot_key, fill = .data$log2FC)
  ) +
    ggplot2::geom_tile(width = 0.96, height = 0.96, colour = NA) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$contrast_plot_label)) +
    ggplot2::scale_x_discrete(labels = compact_label, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = protein_labels, drop = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = nature_palette("signed")[["low"]],
      mid = nature_palette("signed")[["mid"]],
      high = nature_palette("signed")[["high"]],
      midpoint = 0, limits = c(-effect_limit, effect_limit),
      breaks = c(-effect_limit, 0, effect_limit),
      name = expression("Protein log"[2] * "FC")
    ) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(
      title.position = "top", barwidth = grid::unit(35, "mm"),
      barheight = grid::unit(2.1, "mm")
    )) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(x = NULL, y = "m12") +
    theme_nature_manuscript_panel(
      base_size = manuscript_text[["normal"]], base_family = "Arial",
      axes = FALSE, publication_legible = TRUE
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        size = manuscript_text[["dense"]], lineheight = 0.85
      ),
      axis.text.y = ggplot2::element_text(size = manuscript_text[["dense"]], face = "italic"),
      axis.title.y = ggplot2::element_text(
        angle = 0, size = manuscript_text[["normal"]], face = "bold",
        margin = ggplot2::margin(r = 1.5)
      ),
      strip.text.x = ggplot2::element_text(
        size = manuscript_text[["normal"]], face = "plain", lineheight = 0.85
      ),
      panel.spacing.x = grid::unit(0.65, "mm"),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5)
    )
  if (nrow(supported)) {
    p <- p + ggplot2::geom_point(data = supported, shape = 8, size = 0.9, colour = "black")
  }
  dir_create(dirname(panel_svg))
  ggplot2::ggsave(
    panel_svg, p, width = 132, height = 66, units = "mm",
    device = function(...) svglite::svglite(..., pointsize = 7, fix_text_size = FALSE),
    bg = "white", limitsize = FALSE
  )
  dir_create(dirname(source_target))
  readr::write_csv(block, source_target, na = "")
  invisible(c(panel = panel_svg, source_data = source_target))
}

manuscript_figure_placeholder_svg <- function(path, panel_id) {
  dir_create(dirname(path))
  lines <- c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<svg xmlns="http://www.w3.org/2000/svg" width="89mm" height="70mm" viewBox="0 0 89 70">',
    '<rect width="89" height="70" fill="#fff7f7" stroke="#b2182b" stroke-width="0.5"/>',
    sprintf('<text x="44.5" y="30" text-anchor="middle" font-family="Arial" font-size="4" font-weight="bold">Panel %s asset missing</text>', panel_id),
    '<text x="44.5" y="38" text-anchor="middle" font-family="Arial" font-size="2.8">Candidate render only - not publication ready</text>',
    '</svg>'
  )
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

manuscript_figure_assemble_svg <- function(panel_paths, panels, figure, target) {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required for self-contained SVG assembly.", call. = FALSE)
  }
  width <- as.numeric(figure$width_mm)
  height <- as.numeric(figure$height_mm)
  n_rows <- max(vapply(panels, function(x) as.integer(x$row), integer(1)))
  n_cols <- max(vapply(panels, function(x) as.integer(x$col), integer(1)))
  margin <- 5
  gap <- 4
  cell_width <- (width - 2 * margin - (n_cols - 1) * gap) / n_cols
  cell_height <- (height - 2 * margin - (n_rows - 1) * gap) / n_rows
  items <- character()
  for (i in seq_along(panels)) {
    panel <- panels[[i]]
    colspan <- as.integer(panel$colspan %||% 1L)
    x <- margin + (as.integer(panel$col) - 1L) * (cell_width + gap)
    y <- margin + (as.integer(panel$row) - 1L) * (cell_height + gap)
    w <- cell_width * colspan + gap * (colspan - 1L)
    h <- cell_height
    uri <- base64enc::dataURI(file = panel_paths[[as.character(panel$id)]], mime = "image/svg+xml")
    label <- sub("^[0-9]+", "", as.character(panel$id))
    items <- c(items,
      sprintf('<image x="%.3f" y="%.3f" width="%.3f" height="%.3f" preserveAspectRatio="xMidYMid meet" href="%s"/>', x, y + 4, w, h - 4, uri),
      sprintf('<rect x="%.3f" y="%.3f" width="7" height="7" fill="white" fill-opacity="0.9"/>', x, y),
      sprintf('<text x="%.3f" y="%.3f" fill="black" font-family="Arial" font-size="5" font-weight="bold">%s</text>', x + 0.8, y + 5.2, label)
    )
  }
  dir_create(dirname(target))
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%smm" height="%smm" viewBox="0 0 %s %s">', width, height, width, height),
    '<rect width="100%" height="100%" fill="white"/>',
    items,
    '</svg>'
  ), target, useBytes = TRUE)
  invisible(target)
}

manuscript_figure_raster_companions <- function(svg_path, png_path, pdf_path) {
  if (!requireNamespace("magick", quietly = TRUE)) {
    warning("Package 'magick' is unavailable; assembled PNG/PDF were not written.", call. = FALSE)
    return(character())
  }
  result <- tryCatch({
    image <- magick::image_read(svg_path, density = 300)
    image <- magick::image_background(image, "white", flatten = TRUE)
    magick::image_write(image, path = png_path, format = "png")
    magick::image_write(image, path = pdf_path, format = "pdf")
    c(png_path, pdf_path)
  }, error = function(e) {
    warning("Could not render assembled PNG/PDF companions: ", conditionMessage(e), call. = FALSE)
    character()
  })
  invisible(result)
}

manuscript_figure_records_as_list <- function(x) {
  if (!nrow(x)) return(list())
  lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
}

manuscript_figure_write_manifest <- function(
    figure_id, figure, paths, panel_manifest, input_manifest, outputs,
    incomplete_panels = character(), deferred_panels = character()) {
  dir_create(paths$reports)
  dir_create(paths$logs)
  panel_manifest_path <- file.path(paths$reports, "panel_manifest.csv")
  input_manifest_path <- file.path(paths$reports, "input_manifest.csv")
  utils::write.csv(panel_manifest, panel_manifest_path, row.names = FALSE, na = "")
  utils::write.csv(input_manifest, input_manifest_path, row.names = FALSE, na = "")
  session_path <- file.path(paths$logs, "sessionInfo.txt")
  capture.output(utils::sessionInfo(), file = session_path)
  manifest_path <- file.path(paths$logs, "run_manifest.yml")
  manifest <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    contract_version = figure$contract_version,
    figure = paste0("Figure ", as.integer(figure_id)),
    repository_relative_contract = relative_to(manuscript_figure_contract_path()),
    repository_root_at_render = repo_root(),
    render_git_commit = git_commit_sha(),
    scientific_recomputation = FALSE,
    inputs = manuscript_figure_records_as_list(input_manifest),
    outputs = lapply(outputs[file.exists(outputs)], function(x) list(
      path = normalizePath(x, winslash = "/", mustWork = FALSE),
      sha256 = file_hash_sha256(x)
    )),
    panel_manifest = normalizePath(panel_manifest_path, winslash = "/", mustWork = FALSE),
    input_manifest = normalizePath(input_manifest_path, winslash = "/", mustWork = FALSE),
    incomplete_panels = as.list(incomplete_panels),
    deferred_panels = as.list(deferred_panels),
    assembled_svg_contract = "self_contained_vector_svg_with_embedded_automated_panel_SVGs_deferred_panels_excluded",
    assembled_png_contract = "rasterized_from_assembled_SVG_via_magick_at_300_dpi",
    assembled_pdf_contract = "raster_backed_PDF_from_assembled_SVG_via_magick_at_300_dpi",
    notes = paste(
      "Explicit manuscript rendering/materialization layer.",
      "Primary statistical models, p-values, FDRs, module identities, and enrichment results are not recomputed."
    ),
    session_info = normalizePath(session_path, winslash = "/", mustWork = FALSE)
  )
  writeLines(yaml::as.yaml(manifest), manifest_path, useBytes = TRUE)
  invisible(manifest_path)
}

manuscript_figure_main <- function(figure_id) {
  args <- manuscript_figure_args()
  if (args$allow_incomplete && !args$output_explicit) {
    stop("--allow-incomplete requires an explicit --output-root; incomplete placeholders cannot enter canonical results.", call. = FALSE)
  }
  figure <- manuscript_figure_contract(figure_id)
  declared_panels <- figure$panels
  deferred_panels <- Filter(function(x) !manuscript_figure_panel_is_automated(x), declared_panels)
  panels <- Filter(manuscript_figure_panel_is_automated, declared_panels)
  if (!is.null(args$panel)) {
    requested <- Filter(function(x) identical(tolower(as.character(x$id)), args$panel), declared_panels)
    if (!length(requested)) stop("Panel not declared for Figure ", figure_id, ": ", args$panel, call. = FALSE)
    if (!manuscript_figure_panel_is_automated(requested[[1]])) {
      stop(
        "Panel ", requested[[1]]$id, " is declared as ",
        as.character(requested[[1]]$automation_status %||% "not_automated"),
        " and must be added during final Illustrator composition.", call. = FALSE
      )
    }
    panels <- requested
  }

  validation <- lapply(panels, manuscript_figure_validate_panel)
  input_manifest <- unique(do.call(rbind, lapply(validation, `[[`, "inputs")))
  missing_rows <- input_manifest[!input_manifest$exists, , drop = FALSE]
  missing_panels <- unique(missing_rows$panel)
  externally_incomplete <- vapply(panels, function(panel) {
    as.character(panel$id) %in% missing_panels &&
      identical(as.character(panel$producer_script), "external_tracked_asset")
  }, logical(1))
  unsafe_missing <- setdiff(missing_panels, vapply(panels[externally_incomplete], function(x) as.character(x$id), character(1)))
  if (nrow(missing_rows)) {
    for (i in seq_len(nrow(missing_rows))) {
      message("[MISSING] panel ", missing_rows$panel[[i]], " ", missing_rows$role[[i]], ": ", missing_rows$input_relative_path[[i]])
    }
  }
  if (args$dry_run) {
    for (i in seq_len(nrow(input_manifest))) {
      message(
        "[DRY-RUN ", if (input_manifest$exists[[i]]) "PASS" else "WARN", "] panel ",
        input_manifest$panel[[i]], " ", input_manifest$role[[i]], ": ",
        input_manifest$input_relative_path[[i]]
      )
    }
    return(invisible(input_manifest))
  }
  if (length(unsafe_missing)) {
    stop(
      "Missing required canonical inputs for panel(s): ", paste(unsafe_missing, collapse = ", "),
      ". See --dry-run for the exact paths.", call. = FALSE
    )
  }
  if (nrow(missing_rows) && !args$allow_incomplete) {
    stop(
      "Figure ", as.integer(figure_id), " is incomplete. Supply the tracked external asset or use ",
      "--allow-incomplete together with an isolated --output-root for a clearly marked candidate render.",
      call. = FALSE
    )
  }
  if (args$check_only) {
    message(
      "Figure ", as.integer(figure_id), " contract check passed for ",
      length(panels), " automated panel(s)",
      if (is.null(args$panel) && length(deferred_panels)) {
        paste0("; deferred to Illustrator: ", paste(vapply(deferred_panels, function(x) x$id, character(1)), collapse = ", "))
      } else "",
      "."
    )
    return(invisible(input_manifest))
  }

  paths <- manuscript_figure_output_paths(args$output_root, figure_id)
  invisible(lapply(paths, dir_create))
  panel_paths <- character()
  panel_records <- list()
  incomplete <- character()

  for (panel in panels) {
    panel_id <- as.character(panel$id)
    panel_svg <- file.path(paths$panels, paste0("figure_", figure_id, sub("^[0-9]+", "", panel_id), ".svg"))
    source_target <- file.path(paths$source_data, paste0("figure_", figure_id, sub("^[0-9]+", "", panel_id), "_source_data.csv"))
    primary <- manuscript_figure_resolve(panel$primary_source)
    figure_source <- manuscript_figure_resolve(panel$figure_source)
    is_external_missing <- panel_id %in% missing_panels &&
      identical(as.character(panel$producer_script), "external_tracked_asset")

    if (is_external_missing) {
      manuscript_figure_placeholder_svg(panel_svg, panel_id)
      source_target <- NA_character_
      incomplete <- c(incomplete, panel_id)
    } else if (identical(as.character(panel$render_mode), "copy_svg")) {
      manuscript_figure_copy(figure_source, panel_svg)
      source_mode <- as.character(panel$source_data_mode %||% "none")
      if (identical(source_mode, "copy_csv")) {
        manuscript_figure_copy(primary, source_target)
      } else if (identical(source_mode, "derive_figure2b")) {
        manuscript_figure_derive_2b_source(primary, source_target)
      } else {
        source_target <- NA_character_
      }
    } else if (identical(as.character(panel$render_mode), "render_m12_heatmap")) {
      manuscript_figure_render_m12(primary, panel_svg, source_target)
    } else {
      stop("Unsupported render mode for panel ", panel_id, ": ", panel$render_mode, call. = FALSE)
    }

    panel_paths[[panel_id]] <- panel_svg
    panel_records[[length(panel_records) + 1L]] <- data.frame(
      figure = paste0("Figure ", as.integer(figure_id)),
      panel = panel_id,
      description = as.character(panel$description),
      producer_script = as.character(panel$producer_script),
      render_mode = as.character(panel$render_mode),
      scientific_contract = as.character(panel$scientific_contract),
      biological_unit = as.character(panel$biological_unit),
      hemisphere_handling = as.character(panel$hemisphere_handling),
      primary_input_relative_path = as.character(panel$primary_source),
      primary_input_sha256 = file_hash_sha256(primary),
      source_data_path = if (is.na(source_target)) NA_character_ else normalizePath(source_target, winslash = "/", mustWork = FALSE),
      source_data_sha256 = if (is.na(source_target)) NA_character_ else file_hash_sha256(source_target),
      panel_path = normalizePath(panel_svg, winslash = "/", mustWork = FALSE),
      panel_sha256 = file_hash_sha256(panel_svg),
      status = if (is_external_missing) "incomplete_placeholder_not_for_publication" else "materialized",
      render_git_commit = git_commit_sha(),
      stringsAsFactors = FALSE
    )
  }

  outputs <- unname(panel_paths)
  if (is.null(args$panel)) {
    assembled_svg <- file.path(paths$assembled, paste0("figure_", figure_id, ".svg"))
    assembled_png <- file.path(paths$assembled, paste0("figure_", figure_id, ".png"))
    assembled_pdf <- file.path(paths$assembled, paste0("figure_", figure_id, ".pdf"))
    manuscript_figure_assemble_svg(panel_paths, panels, figure, assembled_svg)
    companions <- manuscript_figure_raster_companions(assembled_svg, assembled_png, assembled_pdf)
    outputs <- c(outputs, assembled_svg, companions)
  }
  if (is.null(args$panel) && length(deferred_panels)) {
    for (panel in deferred_panels) {
      panel_records[[length(panel_records) + 1L]] <- data.frame(
        figure = paste0("Figure ", as.integer(figure_id)),
        panel = as.character(panel$id),
        description = as.character(panel$description),
        producer_script = as.character(panel$producer_script),
        render_mode = as.character(panel$render_mode),
        scientific_contract = as.character(panel$scientific_contract),
        biological_unit = as.character(panel$biological_unit),
        hemisphere_handling = as.character(panel$hemisphere_handling),
        primary_input_relative_path = NA_character_,
        primary_input_sha256 = NA_character_,
        source_data_path = NA_character_,
        source_data_sha256 = NA_character_,
        panel_path = NA_character_,
        panel_sha256 = NA_character_,
        status = as.character(panel$automation_status %||% "not_automated"),
        render_git_commit = git_commit_sha(),
        stringsAsFactors = FALSE
      )
    }
  }
  panel_manifest <- do.call(rbind, panel_records)
  declared_order <- vapply(declared_panels, function(x) as.character(x$id), character(1))
  panel_manifest <- panel_manifest[order(match(panel_manifest$panel, declared_order)), , drop = FALSE]
  manifest_path <- manuscript_figure_write_manifest(
    figure_id, figure, paths, panel_manifest, input_manifest, outputs,
    incomplete_panels = incomplete,
    deferred_panels = if (is.null(args$panel)) {
      vapply(deferred_panels, function(x) as.character(x$id), character(1))
    } else character()
  )
  message(
    "Figure ", as.integer(figure_id), " materialized under ", paths$figures,
    if (length(incomplete)) paste0("; incomplete panel(s): ", paste(incomplete, collapse = ", ")) else ""
  )
  invisible(list(
    panels = panel_paths, outputs = outputs, panel_manifest = panel_manifest,
    input_manifest = input_manifest, run_manifest = manifest_path
  ))
}
