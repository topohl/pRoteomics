# Shared spatial visual grammar for the manuscript figures (Part 21).
#
# One ordering, one set of labels, one header/separator style, used by every
# spatial panel so a reader learns the grammar once and recognises it
# everywhere. Everything is driven by config/manuscript_spatial_order.yml;
# nothing here parses anatomy out of a unit string at plot time.
#
# This module is DOWNSTREAM ONLY. It reshapes and labels canonical values. It
# fits no model, runs no test and computes no p-value.

sg_contract_path <- function() repo_path("config", "manuscript_spatial_order.yml")

sg_contract <- local({
  cache <- NULL
  function(path = sg_contract_path()) {
    if (!is.null(cache)) return(cache)
    y <- yaml::read_yaml(path)
    cache <<- y
    y
  }
})

# ---------------------------------------------------------------- unit table

# The 18 canonical units as a data.frame in display order.
sg_units <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    u <- sg_contract()$units
    out <- do.call(rbind, lapply(u, function(z) data.frame(
      dataset = as.character(z$dataset),
      order = as.integer(z$order),
      region = as.character(z$region),
      layer = if (is.null(z$layer)) NA_character_ else as.character(z$layer),
      layer_is_resolution = !identical(z$layer_is_resolution, FALSE) &&
        !is.null(z$layer),
      analysis_key = as.character(z$analysis_key),
      atlas_key = as.character(z$atlas_key),
      display = as.character(z$display),
      stringsAsFactors = FALSE)))
    out <- out[order(out$order), , drop = FALSE]
    rownames(out) <- NULL
    cache <<- out
    out
  }
})

sg_compartments <- function() {
  cs <- sg_contract()$compartments
  d <- do.call(rbind, lapply(cs, function(z) data.frame(
    id = as.character(z$id), order = as.integer(z$order),
    display = as.character(z$display), short = as.character(z$short),
    resolution = as.character(z$resolution), stringsAsFactors = FALSE)))
  d[order(d$order), , drop = FALSE]
}

sg_compartment_levels <- function() sg_compartments()$id
sg_compartment_label <- function(x) {
  cs <- sg_compartments()
  cs$short[match(as.character(x), cs$id)]
}

sg_region_levels <- function() {
  r <- sg_contract()$regions
  vapply(r[order(vapply(r, function(z) z$order, numeric(1)))],
         function(z) as.character(z$id), character(1))
}

sg_layer_display <- function(x) {
  ls <- sg_contract()$layers
  key <- vapply(ls, function(z) as.character(z$id), character(1))
  val <- vapply(ls, function(z) as.character(z$display), character(1))
  unname(val[match(as.character(x), key)])
}

# ------------------------------------------------------------- alias resolve

# Every alias, resolved WITHIN its dataset. The unit string alone is ambiguous:
# soma and microglia both spell their units ca1/ca2/ca3/dg in the atlas tables.
sg_alias_map <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    u <- sg_contract()$units
    rows <- list()
    for (z in u) {
      keys <- unique(c(as.character(z$analysis_key), as.character(z$atlas_key),
                       as.character(unlist(z$aliases %||% character(0)))))
      keys <- keys[nzchar(keys)]
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = as.character(z$dataset),
        alias = tolower(keys),
        analysis_key = as.character(z$analysis_key),
        stringsAsFactors = FALSE)
    }
    m <- unique(do.call(rbind, rows))
    dup <- duplicated(paste(m$dataset, m$alias))
    if (any(dup)) {
      stop("manuscript_spatial_order.yml: alias collides within a dataset: ",
           paste(unique(paste(m$dataset[dup], m$alias[dup])), collapse = ", "),
           call. = FALSE)
    }
    cache <<- m
    m
  }
})

# Resolve (dataset, unit) to the canonical analysis_key. Fails closed: an
# unrecognised pair is an error, never a silently dropped row.
sg_resolve_unit <- function(unit, dataset, required = TRUE) {
  u <- tolower(trimws(as.character(unit)))
  d <- as.character(dataset)
  if (length(d) == 1L) d <- rep(d, length(u))
  if (length(d) != length(u)) {
    stop("sg_resolve_unit: dataset must be length 1 or length(unit)", call. = FALSE)
  }
  m <- sg_alias_map()
  out <- m$analysis_key[match(paste(d, u), paste(m$dataset, m$alias))]
  bad <- is.na(out) & !is.na(u) & nzchar(u)
  if (any(bad) && required) {
    stop("sg_resolve_unit: unrecognised (dataset, unit) pair(s): ",
         paste(unique(paste0(d[bad], "/", unit[bad])), collapse = ", "),
         call. = FALSE)
  }
  out
}

# Ordered factor over the canonical analysis keys, restricted to what is present.
sg_unit_factor <- function(unit, dataset, drop_unused = TRUE) {
  key <- sg_resolve_unit(unit, dataset)
  lev <- sg_units()$analysis_key
  if (drop_unused) lev <- lev[lev %in% key]
  factor(key, levels = lev)
}

# Short display label ("CA1 SLM", "CA1"). Region-only compartments never gain a
# layer suffix, because they have no laminar resolution.
sg_unit_label <- function(unit, dataset) {
  key <- sg_resolve_unit(unit, dataset)
  u <- sg_units()
  u$display[match(paste(dataset, key), paste(u$dataset, u$analysis_key))]
}

sg_unit_region <- function(unit, dataset) {
  key <- sg_resolve_unit(unit, dataset)
  u <- sg_units()
  u$region[match(paste(dataset, key), paste(u$dataset, u$analysis_key))]
}

sg_unit_layer <- function(unit, dataset) {
  key <- sg_resolve_unit(unit, dataset)
  u <- sg_units()
  u$layer[match(paste(dataset, key), paste(u$dataset, u$analysis_key))]
}

# TRUE only where the layer token reflects real laminar resolution. Soma
# carries an sp/sg token naming the layer that was dissected, but only one
# layer per region was sampled, so it is NOT laminar resolution.
sg_has_layer_resolution <- function(unit, dataset) {
  key <- sg_resolve_unit(unit, dataset)
  u <- sg_units()
  u$layer_is_resolution[match(paste(dataset, key), paste(u$dataset, u$analysis_key))]
}

# Annotate an arbitrary canonical table in one call.
sg_annotate <- function(df, unit_col = "spatial_unit", dataset_col = "dataset") {
  if (!unit_col %in% names(df)) stop("sg_annotate: no column ", unit_col, call. = FALSE)
  if (!dataset_col %in% names(df)) stop("sg_annotate: no column ", dataset_col, call. = FALSE)
  df$sg_unit <- sg_resolve_unit(df[[unit_col]], df[[dataset_col]])
  df$sg_compartment <- factor(as.character(df[[dataset_col]]),
                              levels = sg_compartment_levels())
  df$sg_compartment_label <- factor(sg_compartment_label(df[[dataset_col]]),
                                    levels = sg_compartments()$short)
  df$sg_region <- factor(sg_unit_region(df[[unit_col]], df[[dataset_col]]),
                         levels = sg_region_levels())
  df$sg_layer <- sg_unit_layer(df[[unit_col]], df[[dataset_col]])
  df$sg_layer_display <- sg_layer_display(df$sg_layer)
  df$sg_layer_is_resolution <- sg_has_layer_resolution(df[[unit_col]], df[[dataset_col]])
  df$sg_unit_label <- sg_unit_label(df[[unit_col]], df[[dataset_col]])
  df$sg_unit_f <- sg_unit_factor(df[[unit_col]], df[[dataset_col]])
  df
}

# ------------------------------------------------------------------ ordering

# The ordered unit keys for a set of datasets, region-major with layers nested.
sg_order_for <- function(datasets = sg_compartment_levels()) {
  u <- sg_units()
  u <- u[u$dataset %in% datasets, , drop = FALSE]
  u[order(u$order), c("dataset", "analysis_key", "display", "region", "layer")]
}

# ------------------------------------------------------- headers / separators

# Where compartment and region blocks begin and end along an ordered axis.
# Returns 1-based positions so a caller can draw strips and rules without
# recomputing the grouping.
#
# Input may be the full long vectors of a tidy table; the axis is built from the
# UNIQUE (dataset, unit) pairs present, in canonical order. Forgetting to
# deduplicate here silently produces one axis position per data row.
sg_blocks <- function(unit_keys, datasets) {
  key <- sg_resolve_unit(unit_keys, datasets)
  ds <- as.character(datasets)
  if (length(ds) == 1L) ds <- rep(ds, length(key))
  keep <- !is.na(key) & !is.na(ds)
  pair <- unique(data.frame(dataset = ds[keep], unit = key[keep],
                            stringsAsFactors = FALSE))
  u <- sg_units()
  idx <- match(paste(pair$dataset, pair$unit), paste(u$dataset, u$analysis_key))
  if (anyNA(idx)) {
    stop("sg_blocks: unit not in the spatial contract: ",
         paste(paste0(pair$dataset, "/", pair$unit)[is.na(idx)], collapse = ", "),
         call. = FALSE)
  }
  ord <- order(u$order[idx])
  pair <- pair[ord, , drop = FALSE]; idx <- idx[ord]
  pos <- seq_len(nrow(pair))
  mk <- function(grp, label) {
    r <- rle(grp)
    end <- cumsum(r$lengths); start <- end - r$lengths + 1L
    data.frame(label = label[start], start = start, end = end,
               mid = (start + end) / 2, stringsAsFactors = FALSE)
  }
  list(
    order = data.frame(pos = pos, dataset = pair$dataset, unit = pair$unit,
                       display = u$display[idx],
                       region = u$region[idx],
                       layer_display = sg_layer_display(u$layer[idx]),
                       layer_is_resolution = u$layer_is_resolution[idx],
                       stringsAsFactors = FALSE),
    compartment = mk(pair$dataset, sg_compartment_label(pair$dataset)),
    region = mk(paste(pair$dataset, u$region[idx]), u$region[idx]))
}

# A compartment header strip plus region subheaders, as ggplot layers to add on
# top of an existing x-axis of ordered spatial units.
sg_header_layers <- function(blocks, y_comp = 1.06, y_region = 1.015,
                             size_comp = 5.4, size_region = 5.0,
                             fam = NULL) {
  if (is.null(fam)) fam <- nv_palette()$typography$family
  cmp <- blocks$compartment; reg <- blocks$region
  n <- nrow(blocks$order)
  list(
    # heavy rules bound compartments, light rules bound regions
    ggplot2::annotate("segment", x = cmp$start - 0.5, xend = cmp$end + 0.5,
                      y = y_comp - 0.012, yend = y_comp - 0.012,
                      linewidth = 0.5, colour = "grey25"),
    ggplot2::annotate("text", x = cmp$mid, y = y_comp, label = cmp$label,
                      family = fam, size = nv_size(size_comp), fontface = "bold",
                      colour = "grey15", vjust = 0),
    ggplot2::annotate("segment", x = reg$start - 0.5, xend = reg$end + 0.5,
                      y = y_region - 0.008, yend = y_region - 0.008,
                      linewidth = 0.25, colour = "grey55"),
    ggplot2::annotate("text", x = reg$mid, y = y_region, label = reg$label,
                      family = fam, size = nv_size(size_region), colour = "grey30",
                      vjust = 0),
    # heavy vertical separators BETWEEN compartments only
    if (nrow(cmp) > 1L)
      ggplot2::annotate("segment", x = head(cmp$end, -1) + 0.5,
                        xend = head(cmp$end, -1) + 0.5,
                        y = -Inf, yend = y_comp - 0.02,
                        linewidth = 0.4, colour = "grey25") else NULL,
    # light separators between regions inside a compartment
    if (nrow(reg) > 1L)
      ggplot2::annotate("segment",
                        x = setdiff(head(reg$end, -1), head(cmp$end, -1)) + 0.5,
                        xend = setdiff(head(reg$end, -1), head(cmp$end, -1)) + 0.5,
                        y = -Inf, yend = y_region - 0.015,
                        linewidth = 0.18, colour = "grey78") else NULL)
}

# Axis labels for an ordered unit axis. Neuropil shows only the LAYER token,
# because its region is already carried by the header strip above it. For
# region-level compartments the region header IS the label, so the axis is left
# blank rather than repeating it. This is what stops the axis from carrying
# compound strings like "neuron_neuropil_CA1_SLM".
sg_axis_labels <- function(blocks, style = c("layer_only", "full", "region")) {
  style <- match.arg(style)
  o <- blocks$order
  if (style == "full") return(o$display)
  if (style == "region") return(o$region)
  ifelse(o$layer_is_resolution & !is.na(o$layer_display), o$layer_display, "")
}

# ------------------------------------------------------------------- guards

# Hard-stop if a caller tries to give a region-level compartment laminar
# resolution, or to drop a real neuropil layer distinction.
sg_assert_resolution <- function(df, unit_col = "spatial_unit",
                                 dataset_col = "dataset") {
  d <- sg_annotate(df, unit_col, dataset_col)
  region_only <- d$sg_compartment %in% c("neuron_soma", "microglia")
  bad <- region_only & !is.na(d$sg_layer) & d$sg_layer_is_resolution
  if (any(bad)) {
    stop("sg_assert_resolution: laminar resolution asserted for a region-level ",
         "compartment: ", paste(unique(d$sg_unit[bad]), collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}
