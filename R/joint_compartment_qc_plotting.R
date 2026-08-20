# Plot-only helpers for the publication rendering of joint-compartment QC.
# These functions do not transform expression values or recompute embeddings.

joint_pub_dataset_palette <- function() {
  palette <- nature_palette("dataset")
  c(
    "Microglia-enriched ROI" = unname(palette[["microglia"]]),
    "Neuropil" = unname(palette[["neuropil"]]),
    "Soma" = unname(palette[["soma"]])
  )
}

joint_pub_dataset_labels <- function() {
  c(
    microglia = "Microglia-enriched ROI",
    neuron_neuropil = "Neuropil",
    neuron_soma = "Soma"
  )
}

joint_pub_dataset_factor <- function(x) {
  labels <- joint_pub_dataset_labels()
  raw <- as.character(x)
  unknown <- setdiff(unique(raw[!is.na(raw)]), names(labels))
  if (length(unknown)) {
    stop("Unknown joint-QC dataset label(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  factor(unname(labels[raw]), levels = unname(labels))
}

joint_pub_group_factor <- function(x) {
  token <- toupper(trimws(as.character(x)))
  mapped <- ifelse(
    token %in% c("1", "CON", "CTRL", "CONTROL"), "CON",
    ifelse(
      token %in% c("2", "RES", "RESILIENT"), "RES",
      ifelse(token %in% c("3", "SUS", "SUSCEPTIBLE"), "SUS", NA_character_)
    )
  )
  bad <- unique(token[!is.na(token) & nzchar(token) & is.na(mapped)])
  if (length(bad)) {
    stop("Unsupported experimental-group value(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  factor(mapped, levels = c("CON", "RES", "SUS"))
}

joint_pub_order_pcs <- function(variance_table) {
  required <- c("PC", "variance_explained")
  missing <- setdiff(required, names(variance_table))
  if (length(missing)) stop("PCA variance table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  pc_number <- suppressWarnings(as.integer(sub("^PC", "", as.character(variance_table$PC))))
  if (anyNA(pc_number) || any(pc_number < 1L) || anyDuplicated(pc_number)) {
    stop("PCA labels must be unique PC<number> values.", call. = FALSE)
  }
  out <- variance_table[order(pc_number), , drop = FALSE]
  out$pc_number <- sort(pc_number)
  out$PC <- factor(as.character(out$PC), levels = as.character(out$PC), ordered = TRUE)
  rownames(out) <- NULL
  out
}

joint_pub_equal_limits <- function(x, y, padding = 0.04) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (!length(x) || !length(y) || any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Equal-axis limits require complete finite coordinates.", call. = FALSE)
  }
  if (!is.finite(padding) || padding < 0) stop("padding must be non-negative.", call. = FALSE)
  xr <- range(x)
  yr <- range(y)
  span <- max(diff(xr), diff(yr))
  if (!is.finite(span) || span <= 0) span <- 1
  half <- span * (0.5 + padding)
  list(
    x = mean(xr) + c(-half, half),
    y = mean(yr) + c(-half, half),
    span = 2 * half
  )
}

joint_pub_pca_class_label_positions <- function(pca_data) {
  required <- c("Dataset", "PC1", "PC2")
  missing <- setdiff(required, names(pca_data))
  if (length(missing)) {
    stop("PCA class-label data are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (any(!is.finite(pca_data$PC1)) || any(!is.finite(pca_data$PC2))) {
    stop("PCA class-label positions require finite coordinates.", call. = FALSE)
  }
  datasets <- c("Microglia-enriched ROI", "Neuropil", "Soma")
  if (!setequal(unique(as.character(pca_data$Dataset)), datasets)) {
    stop("PCA class-label positions require all three canonical datasets.", call. = FALSE)
  }
  span <- max(diff(range(pca_data$PC1)), diff(range(pca_data$PC2)))
  rules <- data.frame(
    Dataset = datasets,
    x_quantile = c(0.50, 0.90, 0.50),
    y_quantile = c(0.10, 0.90, 0.75),
    x_offset_fraction = c(0, 0.02, -0.02),
    y_offset_fraction = c(-0.02, 0.015, 0.015),
    class_label = c("Microglia", "Neuropil", "Neuron soma"),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(rules)), function(i) {
    z <- pca_data[as.character(pca_data$Dataset) == rules$Dataset[[i]], , drop = FALSE]
    data.frame(
      Dataset = rules$Dataset[[i]],
      class_label = rules$class_label[[i]],
      PC1 = unname(stats::quantile(z$PC1, rules$x_quantile[[i]], names = FALSE)) +
        rules$x_offset_fraction[[i]] * span,
      PC2 = unname(stats::quantile(z$PC2, rules$y_quantile[[i]], names = FALSE)) +
        rules$y_offset_fraction[[i]] * span,
      label_position_rule = paste0(
        "data_quantiles_q", rules$x_quantile[[i]], "_q", rules$y_quantile[[i]],
        "_plus_span_offsets_", rules$x_offset_fraction[[i]], "_",
        rules$y_offset_fraction[[i]]
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$Dataset <- factor(out$Dataset, levels = datasets)
  rownames(out) <- NULL
  out
}

joint_pub_select_detection_features <- function(detected_binary, dataset, n = 200L) {
  detected_binary <- as.matrix(detected_binary)
  storage.mode(detected_binary) <- "numeric"
  dataset <- as.character(dataset)
  n <- as.integer(n)
  if (is.null(rownames(detected_binary)) || anyNA(rownames(detected_binary)) || any(!nzchar(rownames(detected_binary))) || anyDuplicated(rownames(detected_binary))) {
    stop("Detection matrix requires unique nonblank ProteinGroupID row names.", call. = FALSE)
  }
  if (length(dataset) != ncol(detected_binary) || anyNA(dataset) || any(!nzchar(dataset))) {
    stop("Detection matrix and dataset labels are not aligned.", call. = FALSE)
  }
  if (anyNA(detected_binary) || any(!detected_binary %in% c(0, 1))) {
    stop("Detection matrix must contain only complete 0/1 values.", call. = FALSE)
  }
  if (!is.finite(n) || n < 1L) stop("n must be a positive integer.", call. = FALSE)
  dataset_levels <- intersect(names(joint_pub_dataset_labels()), unique(dataset))
  if (length(dataset_levels) < 2L) stop("At least two biological datasets are required.", call. = FALSE)
  rates <- vapply(dataset_levels, function(ds) rowMeans(detected_binary[, dataset == ds, drop = FALSE]), numeric(nrow(detected_binary)))
  if (is.null(dim(rates))) rates <- matrix(rates, ncol = 1L, dimnames = list(rownames(detected_binary), dataset_levels))
  variation <- apply(rates, 1L, function(z) max(z) - min(z))
  rate_sd <- apply(rates, 1L, stats::sd)
  rate_sd[!is.finite(rate_sd)] <- 0
  ord <- order(-variation, -rate_sd, rownames(detected_binary), method = "radix")
  ord <- head(ord, min(n, length(ord)))
  out <- data.frame(
    ProteinGroupID = rownames(detected_binary)[ord],
    between_dataset_detection_range = unname(variation[ord]),
    between_dataset_detection_sd = unname(rate_sd[ord]),
    stringsAsFactors = FALSE
  )
  for (ds in dataset_levels) out[[paste0("detection_rate_", ds)]] <- rates[ord, ds]
  out
}

joint_pub_order_samples <- function(metadata) {
  required <- c("Sample", "dataset", "plate", "region", "layer")
  missing <- setdiff(required, names(metadata))
  if (length(missing)) stop("Sample metadata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyNA(metadata$Sample) || any(!nzchar(as.character(metadata$Sample))) || anyDuplicated(as.character(metadata$Sample))) {
    stop("Sample metadata requires unique nonblank Sample identifiers.", call. = FALSE)
  }
  dataset_rank <- match(as.character(metadata$dataset), names(joint_pub_dataset_labels()))
  if (anyNA(dataset_rank)) stop("Sample metadata contains unsupported dataset labels.", call. = FALSE)
  order(
    dataset_rank,
    as.character(metadata$plate),
    as.character(metadata$region),
    as.character(metadata$layer),
    as.character(metadata$Sample),
    method = "radix"
  )
}

joint_pub_assert_coordinates <- function(reference, candidate, coordinate_columns, label) {
  if (!"Sample" %in% names(reference) || !"Sample" %in% names(candidate)) {
    stop(label, ": both tables require Sample.", call. = FALSE)
  }
  missing <- setdiff(coordinate_columns, intersect(names(reference), names(candidate)))
  if (length(missing)) stop(label, ": missing coordinate columns: ", paste(missing, collapse = ", "), call. = FALSE)
  idx <- match(as.character(reference$Sample), as.character(candidate$Sample))
  if (anyNA(idx) || nrow(reference) != nrow(candidate)) stop(label, ": sample sets differ.", call. = FALSE)
  for (column in coordinate_columns) {
    if (!identical(reference[[column]], candidate[[column]][idx])) {
      stop(label, ": coordinate values changed for ", column, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

joint_pub_theme <- function(base_size = 7, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "#1a1a1a"),
      line = ggplot2::element_line(linewidth = 0.4, colour = "#1a1a1a"),
      axis.line = ggplot2::element_line(linewidth = 0.4, colour = "#1a1a1a"),
      axis.ticks = ggplot2::element_line(linewidth = 0.4, colour = "#1a1a1a"),
      axis.ticks.length = grid::unit(1.1, "mm"),
      axis.text = ggplot2::element_text(size = base_size, colour = "#1a1a1a"),
      axis.title = ggplot2::element_text(size = base_size, colour = "#1a1a1a"),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "plain", colour = "#1a1a1a"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.title = ggplot2::element_text(size = base_size - 0.4, face = "plain"),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.key.width = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(0.8, "mm"),
      legend.spacing.y = grid::unit(0.2, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2.2, 2.2, 2.2, 2.2),
      panel.spacing = grid::unit(1.4, "mm")
    )
}

joint_pub_panel_tag_theme <- function(base_size = 7, base_family = "sans") {
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = base_family,
      face = "bold",
      size = base_size,
      colour = "#1a1a1a",
      hjust = 0,
      vjust = 1
    ),
    plot.tag.position = c(0, 1)
  )
}

joint_pub_save_svg <- function(plot, filename, width_mm, height_mm) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = function(...) svglite::svglite(..., pointsize = 7, fix_text_size = FALSE),
    limitsize = FALSE,
    bg = "white"
  )
  invisible(filename)
}
