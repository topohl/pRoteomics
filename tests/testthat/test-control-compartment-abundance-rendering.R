testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "plotting_nature.R"))
source(testthat::test_path("..", "..", "R", "joint_compartment_qc_plotting.R"))
source(testthat::test_path("..", "..", "R", "compartment_abundance_utils.R"))
source(testthat::test_path("..", "..", "R", "control_compartment_abundance_rendering.R"))

ca_render_test_config <- function() {
  data.frame(
    marker_class = ca_display_class_order,
    official_gene_symbol = c("Soma1", "Neur1", "Mg1"),
    label_order = 1L,
    rank_label = TRUE,
    rank_label_dataset = unname(ca_display_intended_dataset[ca_display_class_order]),
    rank_label_order = 1L,
    selection_timing = "created_with_analysis_independence_not_documented",
    selection_independence_from_observed_results = "not_established",
    observed_abundance_influenced_marker_selection = "not_documented",
    selection_basis = "external_registry_plus_technical_screen",
    display_rationale = "Technically screened display representative",
    rank_label_reason = "Label in intended compartment",
    stringsAsFactors = FALSE
  )
}

ca_render_test_audit <- function() {
  selected <- data.frame(
    marker_gene = c("Soma1", "Neur1", "Mg1"),
    fidelity_marker_class = ca_display_class_order,
    ProteinGroupID = c("P1", "P2", "P3"),
    stringsAsFactors = FALSE
  )
  alternatives <- data.frame(
    marker_gene = c("Soma2", "Neur2", "Mg2"),
    fidelity_marker_class = ca_display_class_order,
    ProteinGroupID = c("Q1", "Q2", "Q3"),
    stringsAsFactors = FALSE
  )
  out <- rbind(selected, alternatives)
  out$matched <- TRUE
  out$primary_score_eligible <- TRUE
  out$gene_claim_eligible <- TRUE
  out$conflicting_marker_classes <- FALSE
  out$protein_group_ambiguity_class <- "single_accession_single_gene"
  out$n_protein_groups_for_marker_gene <- 1L
  out$fidelity_subpanel <- c("nucleus", "synapse", "allen", "nucleus", "synapse", "allen")
  out$source_type <- "external_cached"
  out$source_name <- "test source"
  out$source_reference <- "test reference"
  out$source_term_or_category <- "test term"
  out$evidence_level <- "test evidence"
  out$selection_rule <- "test external rule"
  out$confidence <- "test confidence"
  out$notes <- "test note"
  out
}

ca_render_test_sources <- function() {
  datasets <- names(ca_display_dataset_labels)
  selected <- data.frame(
    marker_class = ca_display_class_order,
    gene = c("Soma1", "Neur1", "Mg1"),
    ProteinGroupID = c("P1", "P2", "P3"),
    stringsAsFactors = FALSE
  )
  dot <- merge(data.frame(dataset = datasets), selected, by = NULL) |>
    dplyr::transmute(
      dataset, ProteinGroupID,
      median_robust_standardized_abundance = c(-5, 0, 5, -2, 0, 2, -1, 0, 1),
      zero_MAD_excluded = FALSE,
      raw_CON_valid_animal_detection_fraction = 1
    )
  rank <- merge(
    data.frame(dataset = datasets, stringsAsFactors = FALSE),
    data.frame(ProteinGroupID = c("P1", "P2", "P3", "B"), stringsAsFactors = FALSE),
    by = NULL
  ) |>
    dplyr::group_by(.data$dataset) |>
    dplyr::arrange(.data$ProteinGroupID, .by_group = TRUE) |>
    dplyr::mutate(
      RankGroup = "CON_region_balanced_mean",
      Rank = dplyr::row_number(), N = dplyr::n(),
      primary_region_balanced_mean_log2 = 20 - .data$Rank,
      valid_animal_count = 3L, observed_animal_count = 3L,
      contributing_animal_count = 3L
    ) |>
    dplyr::ungroup()
  list(dot = dot, rank = rank)
}

testthat::test_that("display config is deterministic, intended-compartment only, and not called prespecified", {
  config_path <- testthat::test_path(
    "..", "..", "config", "marker_panels", "control_compartment_publication_labels.csv"
  )
  config <- ca_validate_marker_display_config(utils::read.csv(config_path, stringsAsFactors = FALSE))
  testthat::expect_equal(nrow(config), 12L)
  testthat::expect_true("Hexa" %in% config$official_gene_symbol)
  testthat::expect_false(config$rank_label[config$official_gene_symbol == "Hexa"])
  testthat::expect_false(any(grepl("prespecif", config$display_rationale, ignore.case = TRUE)))
  testthat::expect_lte(max(table(config$rank_label_dataset[config$rank_label])), 4L)
  testthat::expect_true(all(
    config$rank_label_dataset[config$rank_label] ==
      unname(ca_display_intended_dataset[config$marker_class[config$rank_label]])
  ))
})

testthat::test_that("tile source retains uncapped values and uses a fixed symmetric display cap", {
  source <- ca_render_test_sources()
  rendering <- ca_prepare_control_rendering_sources(
    source$dot, source$rank, ca_render_test_audit(), ca_render_test_config()
  )
  testthat::expect_equal(nrow(rendering$dot), 9L)
  testthat::expect_equal(range(rendering$dot$displayed_robust_z), c(-3, 3))
  testthat::expect_equal(range(rendering$dot$median_robust_standardized_abundance), c(-5, 5))
  testthat::expect_true(all(rendering$dot$valid_animal_count == 3L))
  bad <- source$dot
  bad$raw_CON_valid_animal_detection_fraction[[1]] <- 2 / 3
  testthat::expect_error(
    ca_prepare_control_rendering_sources(bad, source$rank, ca_render_test_audit(), ca_render_test_config()),
    "Equal tile size"
  )
})

testthat::test_that("heatmap has equal tiles, one fill bar, and no size encoding", {
  source <- ca_render_test_sources()
  rendering <- ca_prepare_control_rendering_sources(
    source$dot, source$rank, ca_render_test_audit(), ca_render_test_config()
  )
  plot <- ca_build_control_marker_tile_heatmap(rendering)
  testthat::expect_s3_class(plot$layers[[1]]$geom, "GeomTile")
  testthat::expect_null(plot$scales$get_scales("size"))
  testthat::expect_identical(
    plot$scales$get_scales("fill")$name,
    "Relative abundance (robust z-score)"
  )
  testthat::expect_equal(plot$scales$get_scales("fill")$limits, c(-3, 3))
})

testthat::test_that("rank plot emphasizes only configured representatives and labels intended facets", {
  source <- ca_render_test_sources()
  rendering <- ca_prepare_control_rendering_sources(
    source$dot, source$rank, ca_render_test_audit(), ca_render_test_config()
  )
  plot <- ca_build_control_rank_abundance_plot(rendering)
  testthat::expect_s3_class(plot$layers[[1]]$geom, "GeomLine")
  testthat::expect_s3_class(plot$layers[[2]]$geom, "GeomPoint")
  testthat::expect_s3_class(plot$layers[[3]]$geom, "GeomSegment")
  testthat::expect_s3_class(plot$layers[[4]]$geom, "GeomText")
  testthat::expect_equal(nrow(plot$layers[[2]]$data), 9L)
  testthat::expect_equal(nrow(plot$layers[[4]]$data), 3L)
  testthat::expect_true(all(
    plot$layers[[4]]$data$dataset == plot$layers[[4]]$data$rank_label_dataset
  ))
  testthat::expect_equal(plot$layers[[1]]$aes_params$colour, "#AFAFAF")
  testthat::expect_equal(plot$layers[[1]]$aes_params$linewidth, 0.19)
  testthat::expect_equal(plot$layers[[2]]$aes_params$size, 0.95)
  testthat::expect_equal(plot$layers[[4]]$aes_params$colour, "#303030")
  testthat::expect_equal(plot$layers[[4]]$aes_params$size, 2.14)
  testthat::expect_true(is.expression(plot$labels$x))
  testthat::expect_true(is.expression(plot$labels$y))
  testthat::expect_identical(plot$scales$get_scales("colour")$name, "Representative markers")
})

testthat::test_that("manual rank-label layout deterministically separates the MG/PVM labels", {
  labels <- data.frame(
    dataset = rep("microglia", 3),
    configured_official_gene_symbol = c("P2ry12", "C1qa", "Ctss"),
    Rank = c(1345, 1536, 3119),
    primary_region_balanced_mean_log2 = c(18.174, 17.925, 16.312),
    stringsAsFactors = FALSE
  )
  laid_out <- ca_rank_label_layout(labels)
  testthat::expect_equal(laid_out$label_x, c(1020, 1830, 2420))
  testthat::expect_equal(laid_out$label_y, c(18.47, 17.66, 16.72))
  testthat::expect_identical(laid_out$label_hjust, c(1, 0, 1))
  testthat::expect_equal(length(unique(paste(laid_out$label_x, laid_out$label_y))), 3L)
})

testthat::test_that("terminal-drop audit is deterministic and preserves the complete terminal segment", {
  rank <- data.frame(
    dataset = "neuron_soma", ProteinGroupID = paste0("P", seq_len(100)),
    Rank = seq_len(100), N = 100L,
    primary_region_balanced_mean_log2 = 200 - seq_len(100),
    valid_animal_count = 2L, observed_animal_count = 2L,
    contributing_animal_count = 2L,
    stringsAsFactors = FALSE
  )
  rank$primary_region_balanced_mean_log2[98:100] <- c(90, 89, 88)
  audit_a <- ca_terminal_rank_drop_audit(rank)
  audit_b <- ca_terminal_rank_drop_audit(rank)
  testthat::expect_identical(audit_a, audit_b)
  testthat::expect_equal(nrow(audit_a), 3L)
  testthat::expect_equal(audit_a$Rank, 98:100)
  testthat::expect_true(all(audit_a$valid_animal_count == 2L))
})

testthat::test_that("04e exposes a rendering-only path based on completed source tables", {
  script <- readLines(testthat::test_path(
    "..", "..", "R", "control_compartment_abundance_workflow_v2.R"
  ), warn = FALSE)
  text <- paste(script, collapse = "\n")
  testthat::expect_match(text, "--render-only", fixed = TRUE)
  testthat::expect_match(text, "v2_17_marker_detection_dot_heatmap_source.csv", fixed = TRUE)
  testthat::expect_match(text, "v2_12_rank_abundance_data.csv", fixed = TRUE)
  testthat::expect_match(text, "v2_11_display_selection_provenance.csv", fixed = TRUE)
  testthat::expect_match(
    text,
    "PROTEOMICS_CONTROL_ABUNDANCE_V2_RENDER_ALLOW_OVERWRITE",
    fixed = TRUE
  )
  testthat::expect_match(text, "--compact-recognizable-candidate", fixed = TRUE)
  testthat::expect_match(text, "v2_03_animal_level_abundance.csv", fixed = TRUE)
  testthat::expect_match(
    text,
    "PROTEOMICS_CONTROL_ABUNDANCE_COMPACT_CANDIDATE_ALLOW_OVERWRITE",
    fixed = TRUE
  )
  testthat::expect_match(text, "--compact-enrichment-differences", fixed = TRUE)
  testthat::expect_match(
    text,
    "figure2d_compact_recognizable_markers_v2_source_data.csv",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "PROTEOMICS_CONTROL_ABUNDANCE_DIFFERENCE_CANDIDATE_ALLOW_OVERWRITE",
    fixed = TRUE
  )
})

testthat::test_that("compact recognizable marker candidate is explicit and bounded", {
  config <- ca_compact_recognizable_marker_config_v2()
  testthat::expect_identical(
    config$marker_gene,
    c(
      "Npm1", "Ptbp2", "Anp32a", "Hdac5",
      "Camk2a", "Snap25", "Syp",
      "P2ry12", "C1qa", "Ctss"
    )
  )
  testthat::expect_identical(
    as.integer(table(factor(config$marker_class, levels = ca_display_class_order))),
    c(4L, 3L, 3L)
  )
  testthat::expect_false(any(c("Rbfox3", "Dlg4", "Hexa") %in% config$marker_gene))
  testthat::expect_identical(
    unique(config$intended_compartment),
    c("Neuron soma", "Neuron neuropil", "Microglia/PVM-enriched ROI")
  )
})

ca_compact_candidate_fixture <- function() {
  requested <- ca_requested_intuitive_markers_v2()
  requested$ProteinGroupID <- paste0("PG", seq_len(nrow(requested)))
  marker_mapping <- transform(
    requested,
    matched = TRUE,
    canonical_marker_eligible = TRUE,
    conflicting_marker_classes = FALSE,
    gene_claim_eligible = TRUE,
    joint_qc_eligible = TRUE,
    joint_qc_exclusion_reason = NA_character_,
    feature_official_gene_symbol = marker_gene,
    intended_primary_eligible = TRUE,
    intended_strict_eligible = TRUE,
    exclusion_reason_primary = NA_character_
  )
  datasets <- names(ca_manuscript_dataset_labels)
  marker_detection <- merge(
    requested[c(
      "marker_class", "ProteinGroupID", "marker_gene", "intended_dataset"
    )],
    data.frame(dataset = datasets, stringsAsFactors = FALSE),
    by = NULL
  )
  marker_detection$valid_animal_count <- 3L
  marker_detection$valid_animal_fraction <- 1
  marker_detection$primary_eligible <- TRUE
  marker_detection$strict_eligible <- TRUE
  marker_detection$reliability_status <- "reliably_detected"
  marker_detection$intended_compartment <-
    marker_detection$dataset == marker_detection$intended_dataset
  marker_detection$normalization_source <- "completed centered source"
  marker_detection$feature_universe <- "completed authoritative v2"
  marker_detection$animal_detection_threshold <- ">=2_of_3_CON_animals"
  marker_detection$region_threshold <- ">=3_of_4_regions_per_valid_hemisphere"
  marker_detection$hemisphere_handling <- "equal_weight_valid_hemispheres"
  marker_detection$marker_source_policy <- "authoritative_external_registry"
  marker_detection$quantitative_value_status <- "observed_nonimputed"
  marker_detection$selection_policy <-
    "external_configuration_plus_intended_detection_no_observed_effect"

  compact <- ca_compact_recognizable_marker_config_v2()
  compact_ids <- marker_mapping$ProteinGroupID[
    match(compact$marker_gene, marker_mapping$marker_gene)
  ]
  unreliable_pid <- compact_ids[[2]]
  unreliable_dataset <- setdiff(
    datasets, compact$intended_dataset[[2]]
  )[[1]]
  unreliable <- marker_detection$ProteinGroupID == unreliable_pid &
    marker_detection$dataset == unreliable_dataset
  marker_detection$valid_animal_count[unreliable] <- 1L
  marker_detection$valid_animal_fraction[unreliable] <- 1 / 3
  marker_detection$primary_eligible[unreliable] <- FALSE
  marker_detection$strict_eligible[unreliable] <- FALSE
  marker_detection$reliability_status[unreliable] <- "not_reliably_detected"

  marker_direction <- requested[c("marker_class", "ProteinGroupID")]
  marker_direction$expected_direction_classification <- "intended_highest"
  discordant_pid <- compact_ids[[1]]
  marker_direction$expected_direction_classification[
    marker_direction$ProteinGroupID == discordant_pid
  ] <- "discordant_intended_not_highest"

  selected <- merge(
    data.frame(
      marker_gene = compact$marker_gene,
      ProteinGroupID = compact_ids,
      stringsAsFactors = FALSE
    ),
    compact, by = "marker_gene", sort = FALSE
  )
  eligible_alternatives <- selected[c(
    "marker_class", "ProteinGroupID", "marker_gene", "intended_dataset"
  )]

  animal_abundance <- merge(
    selected[c("ProteinGroupID")],
    merge(
      data.frame(dataset = datasets, stringsAsFactors = FALSE),
      data.frame(AnimalID = paste0("A", 1:3), stringsAsFactors = FALSE),
      by = NULL
    ),
    by = NULL
  )
  animal_abundance$animal_log2_abundance <-
    match(animal_abundance$dataset, datasets) +
    as.integer(sub("^PG", "", animal_abundance$ProteinGroupID)) / 100
  animal_abundance$aggregation <- "primary_median_hierarchy"
  animal_abundance <- animal_abundance[!(
    animal_abundance$ProteinGroupID == unreliable_pid &
      animal_abundance$dataset == unreliable_dataset &
      animal_abundance$AnimalID != "A1"
  ), , drop = FALSE]

  rank_table <- data.frame(
    dataset = "neuron_soma", ProteinGroupID = "unused",
    official_gene_symbol = "unused", RankGroup = "CON_region_balanced_mean",
    valid_animal_count = 3L, valid_animal_fraction = 1,
    primary_eligible = TRUE, strict_eligible = TRUE,
    reliability_status = "reliably_detected", joint_qc_eligible = TRUE,
    joint_qc_exclusion_reason = NA_character_, gene_level_claim_allowed = TRUE,
    stringsAsFactors = FALSE
  )
  list(
    animal_abundance = animal_abundance,
    marker_mapping = marker_mapping,
    marker_detection = marker_detection,
    marker_direction = marker_direction,
    eligible_alternatives = eligible_alternatives,
    rank_table = rank_table,
    discordant_pid = discordant_pid,
    unreliable_pid = unreliable_pid,
    unreliable_dataset = unreliable_dataset
  )
}

testthat::test_that("compact marker selection ignores observed direction and retains unreliable off-target cells", {
  x <- ca_compact_candidate_fixture()
  rendering <- ca_prepare_compact_recognizable_markers_v2(
    x$animal_abundance, x$marker_mapping, x$marker_detection,
    x$marker_direction, x$eligible_alternatives, x$rank_table
  )
  testthat::expect_true(x$discordant_pid %in% rendering$provenance$ProteinGroupID)
  discordant <- rendering$provenance$ProteinGroupID == x$discordant_pid
  testthat::expect_identical(
    rendering$provenance$expected_direction_classification[discordant],
    "discordant_intended_not_highest"
  )
  testthat::expect_false(any(
    rendering$provenance$marker_selection_used_observed_cross_compartment_direction
  ))
  unreliable <- rendering$dot$ProteinGroupID == x$unreliable_pid &
    rendering$dot$dataset == x$unreliable_dataset
  testthat::expect_equal(sum(unreliable), 1L)
  testthat::expect_identical(
    rendering$dot$quantitative_value_status[unreliable],
    "not_quantitatively_estimable"
  )
  testthat::expect_true(is.na(rendering$dot$displayed_centered_log2[unreliable]))
  testthat::expect_identical(rendering$dot$missing_symbol[unreliable], "\u00d7")

  difference_rendering <- ca_prepare_marker_enrichment_differences_v2(rendering$dot)
  differences <- difference_rendering$source
  affected <- differences$ProteinGroupID == x$unreliable_pid &
    differences$comparator_preparation == x$unreliable_dataset
  testthat::expect_equal(sum(affected), 1L)
  testthat::expect_identical(
    differences$comparison_status[affected], "not_quantitatively_estimable"
  )
  testthat::expect_true(is.na(differences$intended_minus_comparator_log2[affected]))
  testthat::expect_silent(
    ggplot2::ggplot_build(ca_build_marker_enrichment_differences_v2(
      difference_rendering
    ))
  )
})

ca_render_v2_fixture <- function() {
  classes <- ca_display_class_order
  selection <- data.frame(
    marker_class = classes,
    ProteinGroupID = c("P1", "P2", "P3"),
    marker_gene = c("Soma1", "Neur1", "Mg1"),
    display_biological_subpanel = c("chromatin", "postsynapse", "Allen microglia/PVM"),
    selection_order_within_class = 1L,
    rank_label = TRUE,
    rank_label_dataset = unname(ca_display_intended_dataset[classes]),
    intended_dataset = unname(ca_display_intended_dataset[classes]),
    stringsAsFactors = FALSE
  )
  datasets <- names(ca_display_dataset_labels)
  dot <- merge(
    data.frame(dataset = datasets, stringsAsFactors = FALSE),
    selection[c("ProteinGroupID")],
    by = NULL
  )
  dot$median_centered_log2 <- c(-5, 0, 5, -2, 0, 2, -1, 0, 1)
  dot$valid_animal_count <- c(3L, 1L, 3L, 3L, 3L, 3L, 3L, 3L, 3L)
  dot$valid_animal_fraction <- dot$valid_animal_count / 3
  dot$reliability_status <- ifelse(
    dot$valid_animal_count >= 2L,
    "reliably_detected",
    "not_reliably_detected"
  )
  intended_lookup <- selection$intended_dataset[match(dot$ProteinGroupID, selection$ProteinGroupID)]
  dot$intended_compartment <- dot$dataset == intended_lookup
  rank <- merge(
    data.frame(dataset = datasets, stringsAsFactors = FALSE),
    data.frame(ProteinGroupID = c("P1", "P2", "P3", "B"), stringsAsFactors = FALSE),
    by = NULL
  ) |>
    dplyr::group_by(.data$dataset) |>
    dplyr::arrange(.data$ProteinGroupID, .by_group = TRUE) |>
    dplyr::mutate(
      RankGroup = "CON_region_balanced_mean",
      Rank = dplyr::row_number(), N = dplyr::n(),
      primary_region_balanced_mean_log2 = 20 - .data$Rank,
      valid_animal_count = 3L, observed_animal_count = 3L,
      contributing_animal_count = 3L
    ) |>
    dplyr::ungroup()
  list(selection = selection, dot = dot, rank = rank)
}

testthat::test_that("v2 dot source has one row per selected marker and dataset", {
  x <- ca_render_v2_fixture()
  rendering <- ca_prepare_control_rendering_sources_v2(
    x$dot, x$rank, x$selection
  )
  testthat::expect_equal(nrow(rendering$dot), 9L)
  testthat::expect_false(anyDuplicated(rendering$dot[c("dataset", "ProteinGroupID")]) > 0)
  testthat::expect_equal(range(rendering$dot$displayed_centered_log2, na.rm = TRUE), c(-3, 3))
  testthat::expect_equal(sum(!rendering$dot$reliably_detected), 1L)
  bad <- rbind(x$dot, x$dot[1, ])
  testthat::expect_error(
    ca_prepare_control_rendering_sources_v2(bad, x$rank, x$selection),
    "one row per selected marker and dataset"
  )
})

testthat::test_that("v2 marker panel uses uniform tiles with non-area detection cues", {
  x <- ca_render_v2_fixture()
  rendering <- ca_prepare_control_rendering_sources_v2(
    x$dot, x$rank, x$selection
  )
  plot <- ca_build_control_marker_dot_heatmap_v2(rendering)
  testthat::expect_s3_class(plot$layers[[1]]$geom, "GeomTile")
  testthat::expect_s3_class(plot$layers[[2]]$geom, "GeomText")
  testthat::expect_null(plot$scales$get_scales("size"))
  testthat::expect_equal(plot$scales$get_scales("fill")$limits, c(-3, 3))
  testthat::expect_null(plot$scales$get_scales("x")$limits)
})

testthat::test_that("v2 rank plot uses deterministic repel labels", {
  x <- ca_render_v2_fixture()
  rendering <- ca_prepare_control_rendering_sources_v2(
    x$dot, x$rank, x$selection
  )
  plot <- ca_build_control_rank_abundance_plot_v2(rendering)
  testthat::expect_s3_class(plot$layers[[1]]$geom, "GeomLine")
  testthat::expect_s3_class(plot$layers[[2]]$geom, "GeomPoint")
  testthat::expect_s3_class(plot$layers[[3]]$geom, "GeomTextRepel")
  testthat::expect_lte(nrow(plot$layers[[3]]$data), 3L)
})

ca_compact_difference_fixture <- function() {
  config <- ca_compact_recognizable_marker_config_v2()
  config$ProteinGroupID <- paste0("P", seq_len(nrow(config)))
  source <- merge(
    config,
    data.frame(
      dataset = c("neuron_soma", "neuron_neuropil", "microglia"),
      stringsAsFactors = FALSE
    ),
    by = NULL
  )
  source$median_centered_log2 <-
    unname(c(neuron_soma = 1.25, neuron_neuropil = -0.75, microglia = 0.125)[source$dataset]) +
    source$selection_order_within_class / 100
  source$expected_direction_classification <- "intended_highest"
  source$normalization_source <- "completed centered source"
  source$feature_universe <- "completed authoritative v2"
  source$quantitative_value_status <- "quantitatively_estimable"
  source$primary_eligible <- TRUE
  source$strict_eligible <- TRUE
  source$reliability_status <- "reliably_detected"
  source$analytical_recomputation <- FALSE
  source
}

testthat::test_that("Figure 2d enrichment differences are exact completed-value subtractions", {
  compact <- ca_compact_difference_fixture()
  before <- compact
  rendering <- ca_prepare_marker_enrichment_differences_v2(
    compact, source_path = "completed.csv", source_sha256 = "abc123"
  )
  source <- rendering$source
  testthat::expect_identical(compact, before)
  testthat::expect_equal(nrow(source), 20L)
  testthat::expect_identical(
    as.numeric(source$intended_minus_comparator_log2),
    as.numeric(
      source$intended_median_centered_log2 -
        source$comparator_median_centered_log2
    )
  )
  testthat::expect_identical(
    as.integer(table(source$ProteinGroupID)), rep(2L, 10L)
  )
  testthat::expect_false(any(source$analytical_recomputation))
  testthat::expect_false(any(source$inferential_statistics))
  testthat::expect_false(any(source$marker_selection_used_difference))
  testthat::expect_identical(unique(source$source_sha256), "abc123")
})

testthat::test_that("Figure 2d difference candidate uses unscaled horizontal points without inference", {
  rendering <- ca_prepare_marker_enrichment_differences_v2(
    ca_compact_difference_fixture()
  )
  plot <- ca_build_marker_enrichment_differences_v2(rendering)
  testthat::expect_s3_class(plot$layers[[1]]$geom, "GeomVline")
  testthat::expect_s3_class(plot$layers[[2]]$geom, "GeomPoint")
  testthat::expect_s3_class(plot$layers[[3]]$geom, "GeomPoint")
  testthat::expect_equal(plot$layers[[1]]$data$xintercept, 0)
  testthat::expect_null(plot$scales$get_scales("x"))
  testthat::expect_null(plot$scales$get_scales("size"))
  testthat::expect_null(plot$scales$get_scales("fill"))
  testthat::expect_identical(
    plot$scales$get_scales("colour")$name, ggplot2::waiver()
  )
  testthat::expect_identical(plot$labels$colour, "Comparator preparation")
})
