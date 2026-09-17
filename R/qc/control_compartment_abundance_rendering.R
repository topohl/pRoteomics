# Rendering-only helpers for the control compartment-abundance figures.
# These functions consume completed 04e source tables and do not recalculate
# inclusion, normalization, detection, aggregation, ranks, or robust z-scores.

ca_display_dataset_labels <- c(
  neuron_soma = "Soma",
  neuron_neuropil = "Neuropil",
  microglia = "MG/PVM ROI"
)

ca_display_class_order <- c("Soma markers", "Neuropil markers", "Microglia/PVM markers")
ca_display_class_headings <- c(
  "Soma markers" = "Soma",
  "Neuropil markers" = "Neuropil",
  "Microglia/PVM markers" = "MG/PVM"
)
ca_manuscript_dataset_labels <- c(
  neuron_soma = "Neuron soma",
  neuron_neuropil = "Neuron neuropil",
  microglia = "Microglia/PVM-enriched ROI"
)
ca_manuscript_marker_class_headings <- c(
  "Soma markers" = "Neuron soma markers",
  "Neuropil markers" = "Neuron neuropil markers",
  "Microglia/PVM markers" = "Microglia/PVM-enriched ROI markers"
)
ca_manuscript_marker_strip_labels <- c(
  "Neuron soma markers" = "Neuron soma markers",
  "Neuron neuropil markers" = "Neuron neuropil markers",
  "Microglia/PVM-enriched ROI markers" = "Microglia/PVM-\nenriched ROI markers"
)
ca_dataset_palette <- nature_palette("dataset")
ca_display_class_palette <- c(
  "Soma markers" = unname(ca_dataset_palette[["soma"]]),
  "Neuropil markers" = unname(ca_dataset_palette[["neuropil"]]),
  "Microglia/PVM markers" = unname(ca_dataset_palette[["microglia"]])
)
ca_display_intended_dataset <- c(
  "Soma markers" = "neuron_soma",
  "Neuropil markers" = "neuron_neuropil",
  "Microglia/PVM markers" = "microglia"
)

ca_display_require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(label, " is missing required column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

ca_display_logical <- function(x, label) {
  if (is.logical(x)) return(x)
  key <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(key))
  out[key %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[key %in% c("false", "f", "0", "no", "n")] <- FALSE
  if (anyNA(out)) stop(label, " must contain only explicit TRUE/FALSE values.", call. = FALSE)
  out
}

ca_display_collapse <- function(x) {
  x <- trimws(as.character(x))
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x)) paste(x, collapse = " | ") else NA_character_
}

ca_validate_marker_display_config <- function(config) {
  required <- c(
    "marker_class", "official_gene_symbol", "label_order", "rank_label",
    "rank_label_dataset", "rank_label_order", "selection_timing",
    "selection_independence_from_observed_results",
    "observed_abundance_influenced_marker_selection", "selection_basis",
    "display_rationale", "rank_label_reason"
  )
  ca_display_require_columns(config, required, "Marker-display configuration")
  config <- as.data.frame(config, stringsAsFactors = FALSE, check.names = FALSE)
  config$marker_class <- trimws(as.character(config$marker_class))
  config$official_gene_symbol <- trimws(as.character(config$official_gene_symbol))
  config$marker_gene_key <- toupper(config$official_gene_symbol)
  config$rank_label <- ca_display_logical(config$rank_label, "rank_label")
  config$label_order <- suppressWarnings(as.integer(config$label_order))
  config$rank_label_order <- suppressWarnings(as.integer(config$rank_label_order))
  if (anyNA(config$label_order) || any(config$label_order < 1L) ||
      anyNA(config$rank_label_order) || any(config$rank_label_order < 1L)) {
    stop("Marker display and rank-label order values must be positive integers.", call. = FALSE)
  }
  if (anyDuplicated(config[c("marker_class", "marker_gene_key")])) {
    stop("Marker-display configuration must contain one row per marker class and gene.", call. = FALSE)
  }
  if (!setequal(unique(config$marker_class), ca_display_class_order)) {
    stop("Marker-display configuration must use exactly the three fidelity classes.", call. = FALSE)
  }
  intended <- unname(ca_display_intended_dataset[config$marker_class])
  if (any(config$rank_label & config$rank_label_dataset != intended)) {
    stop("Rank labels must be assigned only to the marker's intended compartment.", call. = FALSE)
  }
  label_counts <- table(config$rank_label_dataset[config$rank_label])
  if (length(label_counts) && any(label_counts > 4L)) {
    stop("Deterministic rank-label configuration permits at most four labels per facet.", call. = FALSE)
  }
  if (any(grepl("prespecif", config$display_rationale, ignore.case = TRUE))) {
    stop("Display rationale must not describe technically screened markers as prespecified.", call. = FALSE)
  }
  config
}

ca_prepare_configured_marker_matches <- function(config, marker_audit) {
  ca_display_require_columns(marker_audit, c(
    "marker_gene", "fidelity_marker_class", "ProteinGroupID", "matched",
    "primary_score_eligible", "gene_claim_eligible", "conflicting_marker_classes",
    "protein_group_ambiguity_class", "n_protein_groups_for_marker_gene"
  ), "Complete marker audit")
  audit <- as.data.frame(marker_audit, stringsAsFactors = FALSE, check.names = FALSE)
  audit$marker_gene_key <- toupper(trimws(as.character(audit$marker_gene)))
  eligible <- audit[
    audit$matched %in% TRUE & audit$primary_score_eligible %in% TRUE &
      audit$gene_claim_eligible %in% TRUE & !(audit$conflicting_marker_classes %in% TRUE),
    , drop = FALSE
  ]
  joined <- dplyr::inner_join(
    config,
    eligible,
    by = c("marker_class" = "fidelity_marker_class", "marker_gene_key"),
    relationship = "many-to-many"
  )
  mapping <- joined |>
    dplyr::transmute(
      marker_class, official_gene_symbol, label_order, rank_label,
      rank_label_dataset, rank_label_order, selection_timing,
      selection_independence_from_observed_results,
      observed_abundance_influenced_marker_selection, selection_basis,
      display_rationale, rank_label_reason, ProteinGroupID,
      protein_group_ambiguity_class, n_protein_groups_for_marker_gene,
      mapping_eligible = .data$matched,
      score_eligible = .data$primary_score_eligible,
      gene_claim_eligible = .data$gene_claim_eligible,
      conflicting_marker_classes = .data$conflicting_marker_classes
    ) |>
    dplyr::distinct()
  counts <- mapping |>
    dplyr::count(.data$marker_class, .data$official_gene_symbol, name = "n_mapped_protein_groups")
  if (nrow(counts) != nrow(config) || any(counts$n_mapped_protein_groups != 1L)) {
    stop("Every configured display gene must map to exactly one eligible canonical ProteinGroupID.", call. = FALSE)
  }
  mapping
}

ca_marker_selection_provenance <- function(config, marker_audit, mapping) {
  audit <- as.data.frame(marker_audit, stringsAsFactors = FALSE, check.names = FALSE)
  audit$marker_gene_key <- toupper(trimws(as.character(audit$marker_gene)))
  selected_ids <- unique(mapping$ProteinGroupID)
  eligible <- audit[
    audit$matched %in% TRUE & audit$primary_score_eligible %in% TRUE &
      audit$gene_claim_eligible %in% TRUE & !(audit$conflicting_marker_classes %in% TRUE),
    , drop = FALSE
  ] |>
    dplyr::distinct(.data$fidelity_marker_class, .data$marker_gene, .data$ProteinGroupID, .keep_all = TRUE)

  alternatives <- lapply(seq_len(nrow(mapping)), function(i) {
    selected <- mapping[i, , drop = FALSE]
    pool <- eligible[
      eligible$fidelity_marker_class == selected$marker_class &
        eligible$ProteinGroupID != selected$ProteinGroupID,
      , drop = FALSE
    ] |>
      dplyr::distinct(.data$marker_gene, .data$ProteinGroupID)
    if (!nrow(pool)) return(NULL)
    data.frame(
      selected_marker_class = selected$marker_class,
      selected_official_gene_symbol = selected$official_gene_symbol,
      selected_ProteinGroupID = selected$ProteinGroupID,
      alternative_marker_gene = pool$marker_gene,
      alternative_ProteinGroupID = pool$ProteinGroupID,
      alternative_is_also_configured = pool$ProteinGroupID %in% selected_ids,
      stringsAsFactors = FALSE
    )
  }) |>
    dplyr::bind_rows()

  source_rows <- dplyr::inner_join(
    config,
    audit,
    by = c("marker_class" = "fidelity_marker_class", "marker_gene_key"),
    relationship = "many-to-many"
  )
  provenance <- source_rows |>
    dplyr::group_by(.data$marker_class, .data$official_gene_symbol) |>
    dplyr::summarise(
      fidelity_subpanel = ca_display_collapse(.data$fidelity_subpanel),
      source_type = ca_display_collapse(.data$source_type),
      source_name = ca_display_collapse(.data$source_name),
      source_reference = ca_display_collapse(.data$source_reference),
      source_term_or_category = ca_display_collapse(.data$source_term_or_category),
      evidence_level = ca_display_collapse(.data$evidence_level),
      external_selection_rule = ca_display_collapse(.data$selection_rule),
      external_confidence = ca_display_collapse(.data$confidence),
      external_notes = ca_display_collapse(.data$notes),
      .groups = "drop"
    ) |>
    dplyr::left_join(mapping, by = c("marker_class", "official_gene_symbol"), relationship = "one-to-one")

  alt_summary <- alternatives |>
    dplyr::group_by(.data$selected_marker_class, .data$selected_official_gene_symbol) |>
    dplyr::summarise(
      n_total_technically_eligible_alternatives = dplyr::n(),
      n_other_configured_alternatives = sum(.data$alternative_is_also_configured),
      n_not_configured_eligible_alternatives = sum(!.data$alternative_is_also_configured),
      all_technically_eligible_alternative_genes = ca_display_collapse(.data$alternative_marker_gene),
      not_configured_eligible_alternative_genes = ca_display_collapse(
        .data$alternative_marker_gene[!.data$alternative_is_also_configured]
      ),
      .groups = "drop"
    )
  provenance <- provenance |>
    dplyr::left_join(
      alt_summary,
      by = c(
        "marker_class" = "selected_marker_class",
        "official_gene_symbol" = "selected_official_gene_symbol"
      ),
      relationship = "one-to-one"
    ) |>
    dplyr::arrange(match(.data$marker_class, ca_display_class_order), .data$label_order)
  list(summary = provenance, alternatives = alternatives)
}

ca_terminal_rank_drop_audit <- function(rank_table, final_fraction = 0.01) {
  ca_display_require_columns(rank_table, c(
    "dataset", "ProteinGroupID", "Rank", "N", "primary_region_balanced_mean_log2",
    "valid_animal_count", "observed_animal_count", "contributing_animal_count"
  ), "Rank-abundance source")
  if (!is.numeric(final_fraction) || length(final_fraction) != 1L ||
      !is.finite(final_fraction) || final_fraction <= 0 || final_fraction > 0.1) {
    stop("final_fraction must be a single number in (0, 0.1].", call. = FALSE)
  }
  primary <- rank_table
  if ("RankGroup" %in% names(primary)) {
    primary <- primary[primary$RankGroup == "CON_region_balanced_mean", , drop = FALSE]
  }
  parts <- lapply(split(primary, primary$dataset), function(x) {
    x <- x[order(x$Rank, x$ProteinGroupID, method = "radix"), , drop = FALSE]
    x$drop_from_previous_log2 <- c(NA_real_, head(x$primary_region_balanced_mean_log2, -1L) -
      tail(x$primary_region_balanced_mean_log2, -1L))
    window_n <- max(10L, ceiling(nrow(x) * final_fraction))
    candidate <- seq.int(max(2L, nrow(x) - window_n + 1L), nrow(x))
    split_row <- candidate[which.max(x$drop_from_previous_log2[candidate])]
    terminal <- x[seq.int(split_row, nrow(x)), , drop = FALSE]
    terminal$terminal_window_fraction <- final_fraction
    terminal$terminal_window_n <- window_n
    terminal$terminal_segment_start_rank <- x$Rank[[split_row]]
    terminal$terminal_segment_n <- nrow(terminal)
    terminal$terminal_segment_initiating_drop_log2 <- x$drop_from_previous_log2[[split_row]]
    terminal$terminal_drop_definition <- paste0(
      "proteins_at_or_after_largest_adjacent_abundance_drop_within_final_",
      format(final_fraction * 100, trim = TRUE), "_percent_of_ranks"
    )
    terminal
  })
  dplyr::bind_rows(parts) |>
    dplyr::arrange(match(.data$dataset, names(ca_display_dataset_labels)), .data$Rank)
}

ca_prepare_control_rendering_sources <- function(dot_source, rank_table, marker_audit, config) {
  config <- ca_validate_marker_display_config(config)
  mapping <- ca_prepare_configured_marker_matches(config, marker_audit)
  provenance <- ca_marker_selection_provenance(config, marker_audit, mapping)

  ca_display_require_columns(dot_source, c(
    "dataset", "ProteinGroupID", "median_robust_standardized_abundance",
    "zero_MAD_excluded", "raw_CON_valid_animal_detection_fraction"
  ), "Dot-heatmap source")
  primary <- rank_table
  if ("RankGroup" %in% names(primary)) {
    primary <- primary[primary$RankGroup == "CON_region_balanced_mean", , drop = FALSE]
  }
  rank_counts <- primary |>
    dplyr::select(
      "dataset", "ProteinGroupID", "valid_animal_count",
      "observed_animal_count", "contributing_animal_count"
    )
  dot_quant <- dot_source |>
    dplyr::select(
      "dataset", "ProteinGroupID", "median_robust_standardized_abundance",
      "zero_MAD_excluded", "raw_CON_valid_animal_detection_fraction"
    ) |>
    dplyr::distinct()
  dot_display <- dot_quant |>
    dplyr::inner_join(mapping, by = "ProteinGroupID", relationship = "many-to-one") |>
    dplyr::left_join(rank_counts, by = c("dataset", "ProteinGroupID"), relationship = "one-to-one") |>
    dplyr::mutate(
      displayed_robust_z = pmax(-3, pmin(3, .data$median_robust_standardized_abundance)),
      robust_z_display_cap = 3,
      dataset_label = unname(ca_display_dataset_labels[.data$dataset]),
      marker_class_heading = unname(ca_display_class_headings[.data$marker_class])
    )
  expected_rows <- nrow(config) * length(ca_display_dataset_labels)
  if (nrow(dot_display) != expected_rows || anyDuplicated(dot_display[c("dataset", "ProteinGroupID")])) {
    stop("Configured tile source must contain one row per dataset and configured ProteinGroupID.", call. = FALSE)
  }
  if (any(dot_display$valid_animal_count != 3L) ||
      any(abs(dot_display$raw_CON_valid_animal_detection_fraction - 1) > sqrt(.Machine$double.eps))) {
    stop(
      "Equal tile size is valid only when every configured protein has 3/3 valid-animal detection in every compartment.",
      call. = FALSE
    )
  }
  if (any(dot_display$zero_MAD_excluded %in% TRUE) ||
      any(!is.finite(dot_display$median_robust_standardized_abundance))) {
    stop("Configured tile source contains a zero-MAD exclusion or non-finite robust z-score.", call. = FALSE)
  }

  rank_representatives <- primary |>
    dplyr::inner_join(
      mapping |> dplyr::rename(configured_official_gene_symbol = "official_gene_symbol"),
      by = "ProteinGroupID", relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      dataset_label = unname(ca_display_dataset_labels[.data$dataset]),
      display_rank_label = .data$rank_label & .data$dataset == .data$rank_label_dataset
    ) |>
    dplyr::arrange(match(.data$dataset, names(ca_display_dataset_labels)), .data$rank_label_order)
  if (nrow(rank_representatives) != expected_rows ||
      anyDuplicated(rank_representatives[c("dataset", "ProteinGroupID")])) {
    stop("Configured rank display must contain one row per dataset and configured ProteinGroupID.", call. = FALSE)
  }
  label_counts <- table(rank_representatives$dataset[rank_representatives$display_rank_label])
  if (length(label_counts) && any(label_counts > 4L)) {
    stop("Prepared rank display exceeds four deterministic labels in a facet.", call. = FALSE)
  }

  primary$dataset_label <- factor(
    unname(ca_display_dataset_labels[primary$dataset]),
    levels = unname(ca_display_dataset_labels)
  )
  rank_representatives$dataset_label <- factor(
    unname(ca_display_dataset_labels[rank_representatives$dataset]),
    levels = unname(ca_display_dataset_labels)
  )
  dot_display$dataset_label <- factor(
    unname(ca_display_dataset_labels[dot_display$dataset]),
    levels = unname(ca_display_dataset_labels)
  )
  dot_display$marker_class <- factor(dot_display$marker_class, levels = ca_display_class_order)
  dot_display$marker_class_heading <- factor(
    dot_display$marker_class_heading,
    levels = unname(ca_display_class_headings[ca_display_class_order])
  )
  ordered_markers <- config$official_gene_symbol[
    order(match(config$marker_class, ca_display_class_order), config$label_order)
  ]
  dot_display$marker_label <- factor(dot_display$official_gene_symbol, levels = rev(ordered_markers))

  list(
    config = config,
    mapping = mapping,
    provenance = provenance$summary,
    alternatives = provenance$alternatives,
    dot = dot_display,
    rank_curve = primary,
    rank_representatives = rank_representatives,
    terminal_drop = ca_terminal_rank_drop_audit(primary)
  )
}

ca_build_control_marker_tile_heatmap <- function(rendering) {
  ggplot2::ggplot(
    rendering$dot,
    ggplot2::aes(.data$dataset_label, .data$marker_label, fill = .data$displayed_robust_z)
  ) +
    ggplot2::geom_tile(width = 0.9, height = 0.86, colour = "white", linewidth = 0.22) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(marker_class_heading),
      scales = "free_y", space = "free_y", switch = "y"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#3B6FB6", mid = "#F2F2F2", high = "#B24A4A", midpoint = 0,
      limits = c(-3, 3), breaks = c(-3, 0, 3),
      name = "Relative abundance (robust z-score)"
    ) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(
      direction = "horizontal", title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(31, "mm"), barheight = grid::unit(2.2, "mm"),
      ticks = TRUE
    )) +
    ggplot2::labs(x = NULL, y = NULL) +
    joint_pub_theme(base_size = 6.2) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 6.0),
      axis.text.y = ggplot2::element_text(size = 6.0, face = "italic"),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(
        angle = 0, hjust = 1, face = "bold", size = 5.8,
        margin = ggplot2::margin(r = 1.2)
      ),
      panel.spacing.y = grid::unit(0.8, "mm"),
      legend.position = "bottom", legend.justification = "center",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(size = 6.0),
      legend.text = ggplot2::element_text(size = 5.6),
      plot.margin = ggplot2::margin(2.2, 2.2, 1.5, 2.2)
    )
}

ca_rank_label_layout <- function(labels) {
  required <- c(
    "dataset", "configured_official_gene_symbol", "Rank",
    "primary_region_balanced_mean_log2"
  )
  ca_display_require_columns(labels, required, "Configured rank labels")
  labels <- as.data.frame(labels, stringsAsFactors = FALSE, check.names = FALSE)
  labels$label_x <- labels$Rank * 1.12
  labels$label_y <- labels$primary_region_balanced_mean_log2 + 0.14
  labels$label_hjust <- 0

  # Explicit positions make the published SVG deterministic and leave clear,
  # short leader lines in the dense soma and MG/PVM portions of the curve.
  manual <- data.frame(
    dataset = c(
      "neuron_soma", "neuron_soma", "neuron_soma", "neuron_soma",
      "neuron_neuropil", "neuron_neuropil", "neuron_neuropil", "neuron_neuropil",
      "microglia", "microglia", "microglia"
    ),
    configured_official_gene_symbol = c(
      "Npm1", "Anp32a", "Ptbp2", "Hdac5",
      "Syp", "Dlg4", "Snap25", "Camk2a",
      "P2ry12", "C1qa", "Ctss"
    ),
    label_x = c(
      290, 96, 760, 3250,
      42, 420, 118, 14,
      1020, 1830, 2420
    ),
    label_y = c(
      20.88, 21.66, 19.66, 15.13,
      22.92, 20.95, 22.33, 24.22,
      18.47, 17.66, 16.72
    ),
    label_hjust = c(
      0, 1, 0, 1,
      1, 0, 0, 1,
      1, 0, 1
    ),
    stringsAsFactors = FALSE
  )
  key <- paste(labels$dataset, labels$configured_official_gene_symbol, sep = "\r")
  manual_key <- paste(manual$dataset, manual$configured_official_gene_symbol, sep = "\r")
  idx <- match(key, manual_key)
  use_manual <- !is.na(idx)
  labels$label_x[use_manual] <- manual$label_x[idx[use_manual]]
  labels$label_y[use_manual] <- manual$label_y[idx[use_manual]]
  labels$label_hjust[use_manual] <- manual$label_hjust[idx[use_manual]]
  labels
}

ca_build_control_rank_abundance_plot <- function(rendering) {
  labels <- rendering$rank_representatives[
    rendering$rank_representatives$display_rank_label %in% TRUE,
    , drop = FALSE
  ]
  labels <- ca_rank_label_layout(labels)
  ggplot2::ggplot(
    rendering$rank_curve,
    ggplot2::aes(.data$Rank, .data$primary_region_balanced_mean_log2)
  ) +
    ggplot2::geom_line(colour = "#AFAFAF", linewidth = 0.19) +
    ggplot2::geom_point(
      data = rendering$rank_representatives,
      ggplot2::aes(colour = .data$marker_class),
      size = 0.95, alpha = 0.98
    ) +
    ggplot2::geom_segment(
      data = labels,
      ggplot2::aes(
        x = .data$Rank, y = .data$primary_region_balanced_mean_log2,
        xend = .data$label_x, yend = .data$label_y
      ),
      inherit.aes = FALSE, colour = "#8C8C8C", linewidth = 0.10,
      lineend = "round", show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(
        x = .data$label_x, y = .data$label_y,
        label = .data$configured_official_gene_symbol, hjust = .data$label_hjust
      ),
      inherit.aes = FALSE, colour = "#303030", size = 2.14,
      family = "sans", fontface = "plain", show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~dataset_label, nrow = 1, scales = "free_x") +
    ggplot2::scale_x_log10(
      labels = scales::label_number(),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_colour_manual(
      values = ca_display_class_palette,
      breaks = ca_display_class_order,
      labels = c(
        "Soma-associated", "Neuropil-associated", "Microglia/PVM-associated"
      ),
      name = "Representative markers"
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE, title.position = "left",
      keywidth = grid::unit(2.2, "mm"), keyheight = grid::unit(2.2, "mm")
    )) +
    ggplot2::labs(
      x = expression("Protein abundance rank (" * log[10] * " scale)"),
      y = expression("Mean log"[2] * " abundance")
    ) +
    joint_pub_theme(base_size = 6.0) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 7.5, face = "bold", colour = "#303030"),
      axis.text = ggplot2::element_text(size = 6.15),
      axis.title = ggplot2::element_text(size = 6.6),
      panel.spacing.x = grid::unit(1.4, "mm"),
      legend.position = "bottom", legend.justification = "center",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(size = 6.15, margin = ggplot2::margin(r = 0.55)),
      legend.text = ggplot2::element_text(size = 5.95),
      legend.spacing.x = grid::unit(0.35, "mm"),
      legend.box.spacing = grid::unit(0.25, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(1.2, 1.7, 0.2, 3.8)
    )
}

ca_control_rendering_legends <- function() {
  list(
    main = paste(
      "Configured marker proteins are shown across neuronal soma, neuronal neuropil and",
      "microglia/PVM-enriched ROI; these columns are compartments, not three cell types.",
      "Tile fill is the median robust z-score of region-balanced mean log2 abundance across",
      "CON biological units. Robust standardization was performed protein-wise using the median",
      "and 1.4826 x MAD. The colour display is capped at -3 and +3, while exported robust-z",
      "values remain uncapped. Tile size is constant because every configured protein had valid",
      "abundance in all three CON animals in every compartment. No significance testing is shown.",
      "Hexa is retained in the marker-selection audit despite its discordant observed pattern."
    ),
    extended = paste(
      "Control-only rank-abundance curves for neuronal soma, neuronal neuropil and",
      "microglia/PVM-enriched ROI. The grey curve contains every protein in the dataset-specific",
      "primary observed CON universe; configured representative proteins are coloured and only",
      "deterministically configured representatives are labelled in their intended compartment.",
      "Mean log2 abundance is the equal-weight arithmetic mean on the log2 scale across observed",
      "Left/Right ROI replicates within animal x region x layer, then contributing layers within",
      "region, then CA1, CA2, CA3 and DG within animal, and finally contributing CON animals.",
      "A valid animal-level value requires at least three regions and primary rank eligibility",
      "requires at least two of three CON animals. Terminal low-abundance proteins are shown without",
      "truncation. Values are relative jointly normalized abundance, not total hippocampal abundance",
      "or absolute protein quantity."
    )
  )
}

ca_prepare_control_rendering_sources_v2 <- function(dot_source, rank_table, display_selection) {
  ca_display_require_columns(display_selection, c(
    "marker_class", "ProteinGroupID", "marker_gene",
    "display_biological_subpanel", "selection_order_within_class",
    "rank_label", "rank_label_dataset", "intended_dataset"
  ), "V2 display selection")
  ca_display_require_columns(dot_source, c(
    "dataset", "ProteinGroupID", "median_centered_log2",
    "valid_animal_count", "valid_animal_fraction",
    "reliability_status", "intended_compartment"
  ), "V2 dot-heatmap source")
  primary <- rank_table
  if ("RankGroup" %in% names(primary)) {
    primary <- primary[primary$RankGroup == "CON_region_balanced_mean", , drop = FALSE]
  }
  mapping <- display_selection |>
    dplyr::distinct(.data$marker_class, .data$ProteinGroupID, .keep_all = TRUE)
  expected_rows <- nrow(mapping) * length(ca_display_dataset_labels)
  dot_display <- dot_source |>
    dplyr::inner_join(mapping, by = "ProteinGroupID", relationship = "many-to-one") |>
    dplyr::mutate(
      displayed_centered_log2 = ifelse(
        is.finite(.data$median_centered_log2),
        pmax(-3, pmin(3, .data$median_centered_log2)),
        NA_real_
      ),
      centered_log2_display_cap = 3,
      reliably_detected = .data$reliability_status == "reliably_detected",
      missing_symbol = ifelse(.data$reliably_detected, "", "\u00d7"),
      outline_colour = ifelse(
        .data$intended_compartment,
        unname(ca_dataset_palette[c(
          neuron_soma = "soma", neuron_neuropil = "neuropil", microglia = "microglia"
        )[.data$dataset]]),
        "#737373"
      ),
      dataset_label = unname(ca_display_dataset_labels[.data$dataset]),
      marker_class_heading = unname(ca_display_class_headings[.data$marker_class]),
      subpanel_order = ca_marker_subpanel_order(
        .data$marker_class, .data$display_biological_subpanel
      )
    )
  if (nrow(dot_display) != expected_rows ||
      anyDuplicated(dot_display[c("dataset", "ProteinGroupID")])) {
    stop("V2 dot source must contain one row per selected marker and dataset.", call. = FALSE)
  }
  if (any(dot_display$valid_animal_count < 0L | dot_display$valid_animal_count > 3L)) {
    stop("V2 dot source has invalid CON animal counts.", call. = FALSE)
  }

  ordered <- mapping |>
    dplyr::mutate(
      class_order = match(.data$marker_class, ca_display_class_order),
      subpanel_order = ca_marker_subpanel_order(
        .data$marker_class, .data$display_biological_subpanel
      )
    ) |>
    dplyr::arrange(
      .data$class_order, .data$subpanel_order,
      .data$selection_order_within_class, toupper(.data$marker_gene)
    )
  marker_levels <- rev(ordered$marker_gene)
  dot_display$marker_label <- factor(dot_display$marker_gene, levels = marker_levels)
  dot_display$marker_class_heading <- factor(
    dot_display$marker_class_heading,
    levels = unname(ca_display_class_headings[ca_display_class_order])
  )
  dot_display$dataset_label <- factor(
    dot_display$dataset_label,
    levels = unname(ca_display_dataset_labels)
  )

  rank_representatives <- primary |>
    dplyr::inner_join(mapping, by = "ProteinGroupID", relationship = "many-to-one") |>
    dplyr::mutate(
      configured_official_gene_symbol = .data$marker_gene,
      dataset_label = factor(
        unname(ca_display_dataset_labels[.data$dataset]),
        levels = unname(ca_display_dataset_labels)
      ),
      display_rank_label = .data$rank_label &
        .data$dataset == .data$rank_label_dataset
    )
  if (anyDuplicated(rank_representatives[c("dataset", "ProteinGroupID")])) {
    stop("V2 rank display is duplicated by dataset and selected ProteinGroupID.", call. = FALSE)
  }
  intended_keys <- paste(mapping$intended_dataset, mapping$ProteinGroupID, sep = "\r")
  observed_rank_keys <- paste(
    rank_representatives$dataset, rank_representatives$ProteinGroupID, sep = "\r"
  )
  if (!all(intended_keys %in% observed_rank_keys)) {
    stop("Every selected marker must be present in its intended rank-abundance facet.", call. = FALSE)
  }
  label_counts <- table(rank_representatives$dataset[rank_representatives$display_rank_label])
  if (length(label_counts) && any(label_counts > 4L)) {
    stop("V2 rank display exceeds four deterministic labels in a facet.", call. = FALSE)
  }
  primary$dataset_label <- factor(
    unname(ca_display_dataset_labels[primary$dataset]),
    levels = unname(ca_display_dataset_labels)
  )
  list(
    config = mapping,
    mapping = mapping,
    provenance = display_selection,
    alternatives = display_selection[0, , drop = FALSE],
    dot = dot_display,
    rank_curve = primary,
    rank_representatives = rank_representatives,
    terminal_drop = ca_terminal_rank_drop_audit(primary)
  )
}

ca_build_control_marker_dot_heatmap_v2 <- function(rendering) {
  missing <- rendering$dot[!(rendering$dot$reliably_detected %in% TRUE), , drop = FALSE]
  ggplot2::ggplot(
    rendering$dot,
    ggplot2::aes(.data$dataset_label, .data$marker_label, fill = .data$displayed_centered_log2)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(colour = .data$outline_colour),
      width = 0.94, height = 0.90, linewidth = 0.30
    ) +
    ggplot2::geom_text(
      data = missing,
      ggplot2::aes(label = .data$missing_symbol),
      colour = "#777777", size = 2.4, family = "sans"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(marker_class_heading),
      scales = "free_y", space = "free_y", switch = "y"
    ) +
    ggplot2::scale_fill_gradient2(
      low = nature_palette("signed")[["low"]],
      mid = nature_palette("signed")[["mid"]],
      high = nature_palette("signed")[["high"]],
      midpoint = 0, limits = c(-3, 3), breaks = c(-3, 0, 3),
      name = expression("Median centered log"[2]), na.value = "#F0F0F0"
    ) +
    ggplot2::scale_colour_identity(guide = "none") +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        direction = "horizontal", title.position = "top",
        barwidth = grid::unit(29, "mm"), barheight = grid::unit(2.2, "mm"),
        theme = ggplot2::theme(legend.text = ggplot2::element_text(size = nature_manuscript_text_sizes_pt()[["normal"]] / 0.7))
      )
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature_manuscript_panel(base_size = nature_manuscript_text_sizes_pt()[["normal"]], base_family = "Arial", axes = FALSE, publication_legible = TRUE) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = nature_manuscript_text_sizes_pt()[["normal"]]),
      axis.text.y = ggplot2::element_text(size = nature_manuscript_text_sizes_pt()[["dense"]], face = "italic"),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(
        angle = 0, hjust = 1, face = "bold", size = nature_manuscript_text_sizes_pt()[["normal"]],
        margin = ggplot2::margin(r = 1.2)
      ),
      panel.spacing.y = grid::unit(0.45, "mm"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing.x = grid::unit(2, "mm"),
      legend.title = ggplot2::element_text(
        size = nature_manuscript_text_sizes_pt()[["normal"]], margin = ggplot2::margin(b = 0.8, r = 1.2)
      ),
      legend.text = ggplot2::element_text(size = nature_manuscript_text_sizes_pt()[["normal"]]),
      plot.margin = ggplot2::margin(1.2, 1.2, 0.8, 1.2)
    )
}

ca_build_control_rank_abundance_plot_v2 <- function(rendering) {
  labels <- rendering$rank_representatives[
    rendering$rank_representatives$display_rank_label %in% TRUE,
    ,
    drop = FALSE
  ]
  ggplot2::ggplot(
    rendering$rank_curve,
    ggplot2::aes(.data$Rank, .data$primary_region_balanced_mean_log2)
  ) +
    ggplot2::geom_line(colour = "#AFAFAF", linewidth = 0.22) +
    ggplot2::geom_point(
      data = rendering$rank_representatives,
      ggplot2::aes(colour = .data$marker_class),
      size = 1.05, alpha = 0.98
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(label = .data$configured_official_gene_symbol),
      seed = 20260724L,
      direction = "both",
      box.padding = 0.28,
      point.padding = 0.12,
      min.segment.length = 0,
      max.overlaps = Inf,
      force = 1.6,
      segment.colour = "#8C8C8C",
      segment.size = 0.12,
      colour = "#303030",
      size = 2.35,
      family = "sans",
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~dataset_label, nrow = 1, scales = "free_x") +
    ggplot2::scale_x_log10(
      labels = scales::label_number(),
      expand = ggplot2::expansion(mult = c(0.03, 0.12))
    ) +
    ggplot2::scale_colour_manual(
      values = ca_display_class_palette,
      breaks = ca_display_class_order,
      labels = c(
        "Soma-associated", "Neuropil-associated", "Microglia/PVM-associated"
      ),
      name = "Representative markers"
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE, title.position = "left",
      keywidth = grid::unit(2.2, "mm"), keyheight = grid::unit(2.2, "mm")
    )) +
    ggplot2::labs(
      x = expression("Protein abundance rank (" * log[10] * " scale)"),
      y = expression("Mean log"[2] * " abundance")
    ) +
    theme_nature_manuscript_panel(base_size = nature_manuscript_text_sizes_pt()[["normal"]], base_family = "Arial", publication_legible = TRUE) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        size = nature_manuscript_text_sizes_pt()[["axis_title"]], face = "bold", colour = "#303030"
      ),
      axis.text = ggplot2::element_text(size = nature_manuscript_text_sizes_pt()[["normal"]]),
      axis.title = ggplot2::element_text(size = nature_manuscript_text_sizes_pt()[["axis_title"]]),
      panel.spacing.x = grid::unit(1.2, "mm"),
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.just = "center",
      legend.title = ggplot2::element_text(
        size = nature_manuscript_text_sizes_pt()[["normal"]], margin = ggplot2::margin(r = 0.55)
      ),
      legend.text = ggplot2::element_text(size = nature_manuscript_text_sizes_pt()[["normal"]]),
      plot.margin = ggplot2::margin(1.2, 1.5, 0.3, 2.2)
    )
}

ca_control_rendering_legends_v2 <- function() {
  list(
    main = paste(
      "Externally defined markers are shown across neuronal soma, neuronal neuropil and",
      "microglia/PVM-enriched ROI; these are spatial preparations, not three purified cell types.",
      "Fill is median within-protein centered log2 abundance across valid CON animals and is",
      "display-capped at -3 and +3 without changing exported values. Point size is the valid-CON-animal",
      "fraction; a cross denotes fewer than two valid animals, for which no quantitative abundance",
      "contrast is made. Dark outlines mark the intended compartment. Animal estimates use observed,",
      "non-imputed values with equal Left/Right hemisphere weight after layer and region aggregation.",
      "No inferential p-values are shown."
    ),
    extended = paste(
      "Control-only rank-abundance curves retain all canonical non-contaminant proteins with at least",
      "two of three valid CON animals in each dataset. The grey curve is the complete eligible universe;",
      "display markers are highlighted, with at most four deterministic labels per intended facet.",
      "Ranks use non-imputed values normalized by the joint shared-core sample offsets and equal-weight",
      "animal estimates. Marker selection used external provenance/rank, intended-compartment detection,",
      "subpanel balance and gene-symbol tie-breaking, never observed compartment effect magnitude."
    )
  )
}

ca_requested_intuitive_markers_v2 <- function() {
  data.frame(
    marker_gene = c(
      "Npm1", "Ptbp2", "Anp32a", "Hdac5", "Rbfox3",
      "Camk2a", "Dlg4", "Snap25", "Syp",
      "P2ry12", "C1qa", "Ctss", "Hexa"
    ),
    marker_class = c(
      rep("Soma markers", 5L), rep("Neuropil markers", 4L),
      rep("Microglia/PVM markers", 4L)
    ),
    intended_dataset = c(
      rep("neuron_soma", 5L), rep("neuron_neuropil", 4L),
      rep("microglia", 4L)
    ),
    intended_compartment = c(
      rep("Neuron soma", 5L), rep("Neuron neuropil", 4L),
      rep("Microglia/PVM-enriched ROI", 4L)
    ),
    stringsAsFactors = FALSE
  )
}

ca_compact_recognizable_marker_config_v2 <- function() {
  requested <- ca_requested_intuitive_markers_v2()
  keep <- requested$marker_gene %in% c(
    "Npm1", "Ptbp2", "Anp32a", "Hdac5",
    "Camk2a", "Snap25", "Syp",
    "P2ry12", "C1qa", "Ctss"
  )
  out <- requested[keep, , drop = FALSE]
  out$selection_order_within_class <- ave(
    seq_len(nrow(out)), out$marker_class, FUN = seq_along
  )
  out$candidate_selection_reason <- paste(
    "explicit_biologically_motivated_recognizable_marker_configuration;",
    "authoritative_external_registry_and_unique_canonical_mapping;",
    "joint_qc_gene_claim_and_intended_compartment_primary_and_strict_eligibility;",
    "selection_frozen_without_observed_cross_compartment_direction_or_effect_magnitude"
  )
  out
}

ca_requested_marker_eligibility_audit_v2 <- function(
    marker_mapping, marker_detection, marker_direction, rank_table) {
  ca_display_require_columns(marker_mapping, c(
    "marker_gene", "ProteinGroupID", "matched", "canonical_marker_eligible",
    "conflicting_marker_classes", "gene_claim_eligible", "joint_qc_eligible",
    "joint_qc_exclusion_reason", "feature_official_gene_symbol",
    "intended_primary_eligible",
    "intended_strict_eligible", "exclusion_reason_primary"
  ), "V2 marker mapping")
  ca_display_require_columns(marker_detection, c(
    "marker_class", "ProteinGroupID", "marker_gene", "dataset",
    "valid_animal_count", "valid_animal_fraction", "primary_eligible",
    "strict_eligible", "reliability_status"
  ), "V2 marker detection")
  ca_display_require_columns(marker_direction, c(
    "marker_class", "ProteinGroupID", "expected_direction_classification"
  ), "V2 marker direction")
  ca_display_require_columns(rank_table, c(
    "dataset", "ProteinGroupID", "official_gene_symbol", "RankGroup",
    "valid_animal_count", "valid_animal_fraction", "primary_eligible",
    "strict_eligible", "reliability_status", "joint_qc_eligible",
    "joint_qc_exclusion_reason", "gene_level_claim_allowed"
  ), "V2 rank-abundance table")

  requested <- ca_requested_intuitive_markers_v2()
  rows <- lapply(seq_len(nrow(requested)), function(i) {
    target <- requested[i, , drop = FALSE]
    mapping_rows <- marker_mapping[
      toupper(as.character(marker_mapping$marker_gene)) == toupper(target$marker_gene),
      , drop = FALSE
    ]
    registry_present <- nrow(mapping_rows) > 0L
    if (registry_present) {
      protein_ids <- unique(as.character(mapping_rows$ProteinGroupID[
        !is.na(mapping_rows$ProteinGroupID) & nzchar(mapping_rows$ProteinGroupID)
      ]))
      if (length(protein_ids) != 1L) {
        stop("Requested marker does not resolve to one ProteinGroupID: ", target$marker_gene, call. = FALSE)
      }
      pid <- protein_ids[[1]]
      map <- mapping_rows[mapping_rows$ProteinGroupID == pid, , drop = FALSE]
      detection <- marker_detection[
        marker_detection$ProteinGroupID == pid &
          marker_detection$dataset == target$intended_dataset,
        , drop = FALSE
      ]
      direction <- marker_direction[
        marker_direction$ProteinGroupID == pid &
          marker_direction$marker_class == target$marker_class,
        , drop = FALSE
      ]
      if (nrow(detection) != 1L || nrow(direction) != 1L) {
        stop("Requested marker audit is not one-to-one for: ", target$marker_gene, call. = FALSE)
      }
      logical_value <- function(name) {
        value <- unique(map[[name]])
        if (length(value) != 1L) stop("Inconsistent marker mapping field: ", name, call. = FALSE)
        value %in% TRUE
      }
      joint_reason <- unique(as.character(map$joint_qc_exclusion_reason))
      joint_reason <- joint_reason[!is.na(joint_reason) & nzchar(joint_reason)]
      data.frame(
        marker_gene = target$marker_gene,
        found_in_canonical_proteome = TRUE,
        found_in_authoritative_marker_registry = TRUE,
        ProteinGroupID = pid,
        canonical_gene_symbol = unique(as.character(map$feature_official_gene_symbol))[[1]],
        intended_compartment = target$intended_compartment,
        valid_CON_animal_count = detection$valid_animal_count,
        valid_CON_animal_fraction = detection$valid_animal_fraction,
        passes_primary_detection_threshold = detection$primary_eligible %in% TRUE,
        passes_strict_detection_threshold = detection$strict_eligible %in% TRUE,
        passes_primary_marker_eligibility = logical_value("intended_primary_eligible"),
        passes_strict_marker_eligibility = logical_value("intended_strict_eligible"),
        excluded_for_joint_qc_or_contamination = !logical_value("joint_qc_eligible"),
        excluded_for_mapping = !logical_value("matched") || !logical_value("gene_claim_eligible"),
        excluded_for_conflicting_marker_classes = logical_value("conflicting_marker_classes"),
        excluded_for_reliability = !(detection$primary_eligible %in% TRUE),
        reliability_status = as.character(detection$reliability_status),
        expected_direction_classification = as.character(direction$expected_direction_classification),
        exclusion_reason = unique(as.character(map$exclusion_reason_primary))[[1]],
        joint_qc_exclusion_reason = if (length(joint_reason)) paste(joint_reason, collapse = ";") else NA_character_,
        exact_source_tables = paste(
          "v2_04_all_marker_canonical_mapping_provenance.csv;",
          "v2_05_marker_by_dataset_detection_abundance.csv;",
          "v2_06_protein_intended_direction.csv"
        ),
        stringsAsFactors = FALSE
      )
    } else {
      rank_rows <- rank_table[
        rank_table$RankGroup %in% "CON_region_balanced_mean" &
          rank_table$dataset %in% target$intended_dataset &
          toupper(as.character(rank_table$official_gene_symbol)) %in%
            toupper(target$marker_gene),
        , drop = FALSE
      ]
      if (nrow(rank_rows) != 1L) {
        stop("Unregistered requested marker is not uniquely resolved in the intended rank universe: ", target$marker_gene, call. = FALSE)
      }
      data.frame(
        marker_gene = target$marker_gene,
        found_in_canonical_proteome = TRUE,
        found_in_authoritative_marker_registry = FALSE,
        ProteinGroupID = as.character(rank_rows$ProteinGroupID),
        canonical_gene_symbol = as.character(rank_rows$official_gene_symbol),
        intended_compartment = target$intended_compartment,
        valid_CON_animal_count = rank_rows$valid_animal_count,
        valid_CON_animal_fraction = rank_rows$valid_animal_fraction,
        passes_primary_detection_threshold = rank_rows$primary_eligible %in% TRUE,
        passes_strict_detection_threshold = rank_rows$strict_eligible %in% TRUE,
        passes_primary_marker_eligibility = FALSE,
        passes_strict_marker_eligibility = FALSE,
        excluded_for_joint_qc_or_contamination = !(rank_rows$joint_qc_eligible %in% TRUE),
        excluded_for_mapping = !(rank_rows$gene_level_claim_allowed %in% TRUE),
        excluded_for_conflicting_marker_classes = FALSE,
        excluded_for_reliability = !(rank_rows$primary_eligible %in% TRUE),
        reliability_status = as.character(rank_rows$reliability_status),
        expected_direction_classification = NA_character_,
        exclusion_reason = "not_in_authoritative_external_marker_registry",
        joint_qc_exclusion_reason = as.character(rank_rows$joint_qc_exclusion_reason),
        exact_source_tables = "v2_12_rank_abundance_data.csv; absent from v2_04 marker mapping",
        stringsAsFactors = FALSE
      )
    }
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ca_prepare_compact_recognizable_markers_v2 <- function(
    animal_abundance, marker_mapping, marker_detection, marker_direction,
    eligible_alternatives, rank_table) {
  audit <- ca_requested_marker_eligibility_audit_v2(
    marker_mapping, marker_detection, marker_direction, rank_table
  )
  config <- ca_compact_recognizable_marker_config_v2()
  selected_audit <- dplyr::inner_join(
    config, audit, by = c("marker_gene", "intended_compartment"),
    relationship = "one-to-one"
  )
  if (nrow(selected_audit) != nrow(config) ||
      !all(selected_audit$passes_primary_marker_eligibility) ||
      !all(selected_audit$passes_strict_marker_eligibility)) {
    stop(
      "Compact candidate markers must meet primary/strict eligibility in the intended compartment.",
      call. = FALSE
    )
  }

  selected_keys <- paste(selected_audit$marker_class, selected_audit$ProteinGroupID, sep = "\r")
  alternative_keys <- paste(
    eligible_alternatives$marker_class, eligible_alternatives$ProteinGroupID, sep = "\r"
  )
  provenance <- eligible_alternatives[alternative_keys %in% selected_keys, , drop = FALSE]
  if (nrow(provenance) != nrow(config) ||
      anyDuplicated(provenance[c("marker_class", "ProteinGroupID")])) {
    stop("Compact candidate markers are not unique authoritative eligible alternatives.", call. = FALSE)
  }
  provenance <- dplyr::left_join(
    provenance,
    selected_audit[, c(
      "marker_class", "ProteinGroupID", "intended_compartment",
      "selection_order_within_class", "candidate_selection_reason",
      "expected_direction_classification"
    )],
    by = c("marker_class", "ProteinGroupID"), relationship = "one-to-one"
  )
  provenance$candidate_metric <- "median_within_protein_centered_log2_across_valid_CON_animals"
  provenance$candidate_is_presentation_only <- TRUE
  provenance$marker_selection_used_observed_cross_compartment_direction <- FALSE
  provenance$marker_selection_used_observed_cross_compartment_effect_magnitude <- FALSE
  provenance$expected_direction_classification_role <-
    "post_selection_validation_only"
  provenance$candidate_source_tables <- paste(
    "v2_03_animal_level_abundance.csv;",
    "v2_04_all_marker_canonical_mapping_provenance.csv;",
    "v2_05_marker_by_dataset_detection_abundance.csv;",
    "v2_06_protein_intended_direction.csv;",
    "v2_11a_display_selection_eligible_alternatives.csv;",
    "v2_12_rank_abundance_data.csv"
  )

  primary_animal <- animal_abundance[
    animal_abundance$aggregation == "primary_median_hierarchy" &
      animal_abundance$ProteinGroupID %in% selected_audit$ProteinGroupID,
    , drop = FALSE
  ]
  centered <- ca_center_animal_abundance(primary_animal) |>
    dplyr::group_by(.data$dataset, .data$ProteinGroupID) |>
    dplyr::summarise(
      median_centered_log2 = stats::median(.data$centered_log2, na.rm = TRUE),
      centered_valid_animal_count = sum(is.finite(.data$centered_log2)),
      .groups = "drop"
    )
  detection <- marker_detection[
    marker_detection$ProteinGroupID %in% selected_audit$ProteinGroupID,
    , drop = FALSE
  ]
  detection <- detection[, c(
    "marker_class", "ProteinGroupID", "dataset", "valid_animal_count",
    "valid_animal_fraction", "reliability_status", "primary_eligible",
    "strict_eligible", "intended_compartment", "normalization_source",
    "feature_universe", "animal_detection_threshold", "region_threshold",
    "hemisphere_handling", "marker_source_policy", "quantitative_value_status",
    "selection_policy"
  ), drop = FALSE]
  dot <- dplyr::inner_join(
    detection,
    selected_audit[, c(
      "marker_class", "ProteinGroupID", "marker_gene", "intended_dataset",
      "selection_order_within_class", "candidate_selection_reason",
      "expected_direction_classification"
    )],
    by = c("marker_class", "ProteinGroupID"), relationship = "many-to-one"
  ) |>
    dplyr::left_join(
      centered, by = c("dataset", "ProteinGroupID"), relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      intended_compartment = .data$dataset == .data$intended_dataset,
      quantitative_source_contract = .data$quantitative_value_status,
      quantitatively_reliable =
        .data$primary_eligible %in% TRUE &
        .data$reliability_status == "reliably_detected" &
        is.finite(.data$median_centered_log2),
      quantitative_value_status = ifelse(
        .data$quantitatively_reliable,
        "quantitatively_estimable",
        "not_quantitatively_estimable"
      ),
      displayed_centered_log2 = ifelse(
        .data$quantitatively_reliable,
        pmax(-3, pmin(3, .data$median_centered_log2)), NA_real_
      ),
      centered_log2_display_cap = 3,
      reliably_detected = .data$quantitatively_reliable,
      missing_symbol = ifelse(.data$reliably_detected, "", "\u00d7"),
      outline_colour = ifelse(
        .data$intended_compartment,
        unname(ca_dataset_palette[c(
          neuron_soma = "soma", neuron_neuropil = "neuropil", microglia = "microglia"
        )[.data$dataset]]),
        "#737373"
      ),
      dataset_label = factor(
        unname(ca_manuscript_dataset_labels[.data$dataset]),
        levels = unname(ca_manuscript_dataset_labels)
      ),
      marker_class_heading = factor(
        unname(ca_manuscript_marker_class_headings[.data$marker_class]),
        levels = unname(ca_manuscript_marker_class_headings)
      ),
      analytical_recomputation = FALSE,
      marker_selection_used_observed_cross_compartment_direction = FALSE,
      marker_selection_used_observed_cross_compartment_effect_magnitude = FALSE,
      expected_direction_classification_role =
        "post_selection_validation_only",
      presentation_transform =
        "existing ca_center_animal_abundance; median across valid CON animals; display cap only"
    )
  ordered <- selected_audit |>
    dplyr::mutate(class_order = match(.data$marker_class, ca_display_class_order)) |>
    dplyr::arrange(.data$class_order, .data$selection_order_within_class)
  dot$marker_label <- factor(dot$marker_gene, levels = rev(ordered$marker_gene))
  expected_rows <- nrow(config) * 3L
  if (nrow(dot) != expected_rows || anyDuplicated(dot[c("dataset", "ProteinGroupID")])) {
    stop("Compact candidate source must contain one row per marker and compartment.", call. = FALSE)
  }
  intended_rows <- dot$intended_compartment %in% TRUE
  if (!all(dot$primary_eligible[intended_rows] %in% TRUE) ||
      !all(dot$strict_eligible[intended_rows] %in% TRUE) ||
      !all(dot$quantitatively_reliable[intended_rows] %in% TRUE)) {
    stop("Compact candidate contains an unreliable intended-compartment value.", call. = FALSE)
  }
  list(audit = audit, provenance = provenance, dot = dot, config = config)
}

ca_build_compact_recognizable_markers_v2 <- function(rendering) {
  ca_build_control_marker_dot_heatmap_v2(rendering) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(marker_class_heading),
      scales = "free_y", space = "free_y", switch = "y",
      labeller = ggplot2::labeller(
        marker_class_heading = ca_manuscript_marker_strip_labels
      )
    ) +
    ggplot2::scale_x_discrete(labels = c(
      "Neuron soma" = "Neuron soma",
      "Neuron neuropil" = "Neuron neuropil",
      "Microglia/PVM-enriched ROI" = "Microglia/PVM-\nenriched ROI"
    )) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = nature_manuscript_text_sizes_pt()[["normal"]], lineheight = 0.92
      ),
      strip.text.y.left = ggplot2::element_text(
        angle = 0, hjust = 1, face = "plain",
        size = nature_manuscript_text_sizes_pt()[["normal"]],
        margin = ggplot2::margin(r = 1.2)
      )
    )
}

ca_prepare_marker_enrichment_differences_v2 <- function(
    compact_source, source_path = NA_character_, source_sha256 = NA_character_) {
  required <- c(
    "marker_class", "ProteinGroupID", "marker_gene", "dataset",
    "intended_dataset", "median_centered_log2",
    "selection_order_within_class", "candidate_selection_reason",
    "expected_direction_classification", "normalization_source",
    "feature_universe", "quantitative_value_status", "primary_eligible",
    "strict_eligible", "reliability_status",
    "analytical_recomputation"
  )
  ca_display_require_columns(compact_source, required, "Compact Figure 2d source")
  if (nrow(compact_source) != 30L ||
      length(unique(compact_source$ProteinGroupID)) != 10L ||
      anyDuplicated(compact_source[c("dataset", "ProteinGroupID")]) ||
      any(compact_source$analytical_recomputation %in% TRUE)) {
    stop(
      "Compact Figure 2d source must retain the completed 10-marker x 3-preparation grid.",
      call. = FALSE
    )
  }

  marker_keys <- c(
    "marker_class", "ProteinGroupID", "marker_gene", "intended_dataset",
    "selection_order_within_class", "candidate_selection_reason",
    "expected_direction_classification", "normalization_source",
    "feature_universe"
  )
  intended <- compact_source |>
    dplyr::filter(.data$dataset == .data$intended_dataset) |>
    dplyr::select(dplyr::all_of(c(
      marker_keys, "dataset", "median_centered_log2", "primary_eligible",
      "strict_eligible", "reliability_status", "quantitative_value_status"
    ))) |>
    dplyr::rename(
      intended_preparation = "dataset",
      intended_median_centered_log2 = "median_centered_log2",
      intended_primary_eligible = "primary_eligible",
      intended_strict_eligible = "strict_eligible",
      intended_reliability_status = "reliability_status",
      intended_quantitative_value_status = "quantitative_value_status"
    )
  if (nrow(intended) != 10L || anyDuplicated(intended$ProteinGroupID) ||
      !all(intended$intended_primary_eligible %in% TRUE) ||
      !all(intended$intended_strict_eligible %in% TRUE) ||
      !all(intended$intended_reliability_status == "reliably_detected") ||
      any(!is.finite(intended$intended_median_centered_log2))) {
    stop(
      "Each compact marker must have one primary/strict-eligible intended preparation.",
      call. = FALSE
    )
  }
  comparators <- compact_source |>
    dplyr::filter(.data$dataset != .data$intended_dataset) |>
    dplyr::select(dplyr::all_of(c(
      marker_keys, "dataset", "median_centered_log2", "primary_eligible",
      "strict_eligible", "reliability_status", "quantitative_value_status"
    ))) |>
    dplyr::rename(
      comparator_preparation = "dataset",
      comparator_median_centered_log2 = "median_centered_log2",
      comparator_primary_eligible = "primary_eligible",
      comparator_strict_eligible = "strict_eligible",
      comparator_reliability_status = "reliability_status",
      comparator_quantitative_value_status = "quantitative_value_status"
    )
  differences <- dplyr::inner_join(
    intended, comparators, by = marker_keys, relationship = "one-to-many"
  ) |>
    dplyr::mutate(
      intended_quantitatively_reliable =
        .data$intended_primary_eligible %in% TRUE &
        .data$intended_reliability_status == "reliably_detected" &
        is.finite(.data$intended_median_centered_log2),
      comparator_quantitatively_reliable =
        .data$comparator_primary_eligible %in% TRUE &
        .data$comparator_reliability_status == "reliably_detected" &
        is.finite(.data$comparator_median_centered_log2),
      comparison_status = ifelse(
        .data$intended_quantitatively_reliable &
          .data$comparator_quantitatively_reliable,
        "quantitatively_estimable",
        "not_quantitatively_estimable"
      ),
      intended_minus_comparator_log2 = ifelse(
        .data$comparison_status == "quantitatively_estimable",
        .data$intended_median_centered_log2 -
          .data$comparator_median_centered_log2,
        NA_real_
      )
    )
  expected <- differences$intended_median_centered_log2 -
    differences$comparator_median_centered_log2
  estimable <- differences$comparison_status == "quantitatively_estimable"
  if (nrow(differences) != 20L ||
      any(table(differences$ProteinGroupID) != 2L) ||
      !isTRUE(all.equal(
        differences$intended_minus_comparator_log2[estimable],
        expected[estimable], tolerance = 0
      )) ||
      any(!is.na(differences$intended_minus_comparator_log2[!estimable]))) {
    stop(
      paste(
        "Figure 2d estimable differences must be exact intended-minus-comparator",
        "subtractions and non-estimable comparisons must remain missing."
      ),
      call. = FALSE
    )
  }

  order_table <- intended |>
    dplyr::mutate(class_order = match(.data$marker_class, ca_display_class_order)) |>
    dplyr::arrange(.data$class_order, .data$selection_order_within_class)
  differences <- differences |>
    dplyr::mutate(
      intended_preparation_label =
        unname(ca_manuscript_dataset_labels[.data$intended_preparation]),
      comparator_preparation_label =
        unname(ca_manuscript_dataset_labels[.data$comparator_preparation]),
      comparator_legend_label = paste("vs", .data$comparator_preparation_label),
      marker_class_heading = factor(
        unname(ca_manuscript_marker_class_headings[.data$marker_class]),
        levels = unname(ca_manuscript_marker_class_headings)
      ),
      marker_label = factor(.data$marker_gene, levels = rev(order_table$marker_gene)),
      difference_plot_value = ifelse(
        .data$comparison_status == "quantitatively_estimable",
        .data$intended_minus_comparator_log2,
        0
      ),
      difference_plot_value_is_quantitative =
        .data$comparison_status == "quantitatively_estimable",
      difference_definition =
        "intended preparation median centered log2 minus comparator preparation median centered log2",
      analytical_recomputation = FALSE,
      inferential_statistics = FALSE,
      marker_selection_used_difference = FALSE,
      marker_selection_used_observed_cross_compartment_direction = FALSE,
      marker_selection_used_observed_cross_compartment_effect_magnitude = FALSE,
      source_path = as.character(source_path),
      source_sha256 = as.character(source_sha256)
    ) |>
    dplyr::arrange(
      match(.data$marker_class, ca_display_class_order),
      .data$selection_order_within_class,
      match(.data$comparator_preparation, names(ca_manuscript_dataset_labels))
    )

  provenance <- intended |>
    dplyr::transmute(
      .data$marker_class, .data$ProteinGroupID, .data$marker_gene,
      .data$intended_dataset, .data$selection_order_within_class,
      .data$candidate_selection_reason,
      .data$expected_direction_classification,
      .data$normalization_source, .data$feature_universe,
      .data$intended_quantitative_value_status,
      source_path = as.character(source_path),
      source_sha256 = as.character(source_sha256),
      source_row_count = nrow(compact_source),
      difference_definition =
        "intended preparation median centered log2 minus comparator preparation median centered log2",
      analytical_recomputation = FALSE,
      inferential_statistics = FALSE,
      marker_selection_used_difference = FALSE,
      marker_selection_used_observed_cross_compartment_direction = FALSE,
      marker_selection_used_observed_cross_compartment_effect_magnitude = FALSE,
      expected_direction_classification_role =
        "post_selection_validation_only",
      presentation_only = TRUE
    )
  list(source = differences, provenance = provenance)
}

ca_build_marker_enrichment_differences_v2 <- function(rendering) {
  palette <- c(
    "vs Neuron soma" = ca_dataset_palette[["soma"]],
    "vs Neuron neuropil" = ca_dataset_palette[["neuropil"]],
    "vs Microglia/PVM-enriched ROI" = ca_dataset_palette[["microglia"]]
  )
  estimable <- rendering$source[
    rendering$source$comparison_status == "quantitatively_estimable", , drop = FALSE
  ]
  not_estimable <- rendering$source[
    rendering$source$comparison_status == "not_quantitatively_estimable", , drop = FALSE
  ]
  ggplot2::ggplot(
    rendering$source,
    ggplot2::aes(
      x = .data$difference_plot_value,
      y = .data$marker_label,
      colour = .data$comparator_legend_label
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.35) +
    ggplot2::geom_point(
      data = estimable,
      position = ggplot2::position_dodge(width = 0.48), size = 2.15,
      alpha = 0.95
    ) +
    ggplot2::geom_point(
      data = not_estimable,
      position = ggplot2::position_dodge(width = 0.48),
      shape = 4, size = 2.1, stroke = 0.65, alpha = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(marker_class_heading),
      scales = "free_y", space = "free_y", switch = "y",
      labeller = ggplot2::labeller(
        marker_class_heading = ca_manuscript_marker_strip_labels
      )
    ) +
    ggplot2::scale_colour_manual(values = palette, breaks = names(palette)) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE, title.position = "top"
    )) +
    ggplot2::labs(
      x = "Δ median centered log₂ abundance\n(intended − comparator)",
      y = NULL, colour = "Comparator preparation",
      caption = "Cross at zero = not estimable (status cue; not a value)."
    ) +
    theme_nature_manuscript_panel(
      base_size = nature_manuscript_text_sizes_pt()[["normal"]],
      base_family = "Arial", publication_legible = TRUE
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "left",
      strip.placement = "outside",
      axis.text.y = ggplot2::element_text(face = "italic"),
      strip.text.y.left = ggplot2::element_text(
        angle = 0, hjust = 1, face = "plain",
        size = nature_manuscript_text_sizes_pt()[["normal"]],
        margin = ggplot2::margin(r = 1.2)
      ),
      panel.spacing.y = grid::unit(1.4, "mm"),
      legend.key.height = grid::unit(3.5, "mm"),
      plot.caption = ggplot2::element_text(
        size = nature_manuscript_text_sizes_pt()[["dense"]],
        hjust = 0, colour = "#555555"
      ),
      plot.margin = ggplot2::margin(1.5, 1.5, 4, 1.5)
    )
}
