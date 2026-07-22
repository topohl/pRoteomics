testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "joint_compartment_qc_plotting.R"))
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
  testthat::expect_equal(plot$layers[[1]]$aes_params$colour, "#B8B8B8")
  testthat::expect_equal(plot$layers[[1]]$aes_params$linewidth, 0.19)
  testthat::expect_equal(plot$layers[[2]]$aes_params$size, 0.95)
  testthat::expect_equal(plot$layers[[4]]$aes_params$colour, "#303030")
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
    "..", "..", "03_qc_exploration", "04e_control_compartment_abundance_publication_figures.r"
  ), warn = FALSE)
  text <- paste(script, collapse = "\n")
  testthat::expect_match(text, "--render-only", fixed = TRUE)
  testthat::expect_match(text, "07a_control_marker_dot_heatmap_source.csv", fixed = TRUE)
  testthat::expect_match(text, "04_rank_abundance_with_percentile.csv", fixed = TRUE)
  testthat::expect_match(text, "05_marker_matches_and_exclusions.csv", fixed = TRUE)
})
