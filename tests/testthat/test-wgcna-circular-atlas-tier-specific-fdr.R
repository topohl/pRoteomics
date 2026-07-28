test_root <- function(...) {
  testthat::test_path("..", "..", ...)
}

read_circular_csv <- function(...) {
  path <- test_root(...)
  testthat::expect_true(file.exists(path), info = path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

testthat::test_that("circular atlas uses the canonical Stage 05 consumer handoff", {
  source(test_root("R", "paths.R"))
  source(test_root("R", "wgcna_group_effect_consumer_utils.R"))

  datasets <- c("neuron_neuropil", "neuron_soma", "microglia")
  levels <- c("supermodule", "module")
  effects <- dplyr::bind_rows(lapply(datasets, function(dataset) {
    dplyr::bind_rows(lapply(levels, function(level) {
      path <- path_results(
        "tables", "06_modules_WGCNA", "group_effects", dataset,
        paste0(level, "_group_effects.csv")
      )
      wgcna_group_effect_consumer_adapt(
        readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
      )
    }))
  }))

  expected_fdr <- dplyr::case_when(
    effects$analysis_tier == "primary_wgcna_global" ~
      effects$FDR_primary_global,
    effects$analysis_tier == "secondary_contextual_global" ~
      effects$FDR_secondary_global,
    effects$analysis_tier == "secondary_spatial_heterogeneity" ~
      effects$FDR_interaction_omnibus,
    effects$analysis_tier == "exploratory_spatial_localization" ~
      effects$FDR_local_exploratory,
    TRUE ~ NA_real_
  )
  testthat::expect_equal(effects$tier_specific_fdr, expected_fdr)

  independent <- effects$independent_hypothesis %in% TRUE
  testthat::expect_true(all(is.finite(
    effects$tier_specific_fdr[independent]
  )))
  aliases <- effects$claim_entity_role == "compatibility_alias"
  testthat::expect_true(all(!effects$independent_hypothesis[aliases]))
  testthat::expect_true(all(is.na(effects$tier_specific_fdr[aliases])))
  testthat::expect_true(all(is.na(effects$tier_specific_family_id[aliases])))
  testthat::expect_true(all(is.na(effects$tier_specific_family_size[aliases])))
  testthat::expect_true(all(
    effects$statistical_support_status[aliases] ==
      "inherited_from_canonical_entity"
  ))

  primary <- effects$analysis_tier == "primary_wgcna_global"
  contextual <- effects$analysis_tier == "secondary_contextual_global"
  interaction <- effects$analysis_tier ==
    "secondary_spatial_heterogeneity"
  local <- effects$analysis_tier == "exploratory_spatial_localization"
  testthat::expect_true(all(effects$contrast[primary] == "SUS - RES"))
  testthat::expect_setequal(
    unique(effects$contrast[contextual]),
    c("RES - CON", "SUS - CON")
  )
  testthat::expect_true(all(
    effects$effect_scope[interaction] ==
      "stress_by_spatial_interaction"
  ))
  testthat::expect_true(all(
    effects$effect_scope[local] == "within_spatial_unit"
  ))
})

testthat::test_that("circular group-effect plotting sources retain complete effects and tier metadata", {
  source_dir <- file.path(
    "results", "source_data", "10_biological_integration",
    "wgcna_circular_atlas", "global"
  )
  required <- c(
    "analysis_tier", "contrast", "effect_scope", "spatial_unit",
    "independent_hypothesis", "estimate", "SE", "p_value",
    "tier_specific_fdr", "tier_specific_family_id",
    "tier_specific_family_size", "statistical_support_status",
    "model_valid_for_inference", "model_stability_status"
  )
  sources <- list(
    segments = read_circular_csv(
      source_dir, "wgcna_circular_atlas_segments.csv"
    ),
    plot_source = read_circular_csv(
      source_dir, "wgcna_circular_atlas_plot_source.csv"
    ),
    supermodule = read_circular_csv(
      source_dir, "wgcna_circular_heatmap_source_supermodule.csv"
    ),
    module = read_circular_csv(
      source_dir, "wgcna_circular_heatmap_source_module.csv"
    )
  )
  for (name in names(sources)) {
    testthat::expect_setequal(
      intersect(required, names(sources[[name]])),
      required
    )
  }
  standardized_required <- c(
    "raw_stage05_estimate", "CI_low", "CI_high", "response_SD",
    "standardized_effect_unclipped", "standardized_effect_display",
    "standardization_scope", "effect_clipped_for_display",
    "FDR_local_exploratory", "q_local_le_0_05",
    "q_local_le_0_10", "q_local_le_0_25",
    "figure_interpretation_guard"
  )
  for (name in c("supermodule", "module")) {
    testthat::expect_setequal(
      intersect(standardized_required, names(sources[[name]])),
      standardized_required
    )
  }

  testthat::expect_setequal(
    unique(sources$supermodule$contrast),
    c("RES - CON", "SUS - CON", "SUS - RES")
  )
  testthat::expect_setequal(
    unique(sources$module$contrast),
    c("RES - CON", "SUS - CON", "SUS - RES")
  )
  testthat::expect_true(all(
    sources$supermodule$effect_scope == "within_spatial_unit"
  ))
  testthat::expect_true(all(
    sources$module$effect_scope == "within_spatial_unit"
  ))
  testthat::expect_true(all(
    sources$supermodule$analysis_tier ==
      "exploratory_spatial_localization"
  ))
  testthat::expect_true(all(
    sources$module$analysis_tier ==
      "exploratory_spatial_localization"
  ))
  testthat::expect_equal(
    sources$supermodule$tier_specific_fdr,
    sources$supermodule$FDR_local_exploratory,
    tolerance = 0
  )
  testthat::expect_equal(
    sources$module$tier_specific_fdr,
    sources$module$FDR_local_exploratory,
    tolerance = 0
  )
  for (source_name in c("supermodule", "module")) {
    source <- sources[[source_name]]
    q <- source$FDR_local_exploratory
    testthat::expect_identical(
      source$q_local_le_0_05,
      dplyr::if_else(is.na(q), NA, q <= 0.05)
    )
    testthat::expect_identical(
      source$q_local_le_0_10,
      dplyr::if_else(is.na(q), NA, q <= 0.10)
    )
    testthat::expect_identical(
      source$q_local_le_0_25,
      dplyr::if_else(is.na(q), NA, q <= 0.25)
    )
  }

  expected_support <- function(x) {
    independent <- x$independent_hypothesis %in% TRUE
    stable <- x$model_valid_for_inference %in% TRUE &
      x$model_stability_status %in%
        c("stable_animal_level_lm", "stable_mixed_model")
    dplyr::case_when(
      !independent ~ "none",
      !stable ~ "invalid_or_unstable",
      !is.na(x$tier_specific_fdr) &
        x$tier_specific_fdr <= 0.05 ~ "FDR05",
      !is.na(x$tier_specific_fdr) &
        x$tier_specific_fdr <= 0.10 ~ "FDR10",
      TRUE ~ "none"
    )
  }
  testthat::expect_identical(
    sources$supermodule$support_class,
    expected_support(sources$supermodule)
  )
  testthat::expect_identical(
    sources$module$support_class,
    expected_support(sources$module)
  )

  aliases <- sources$supermodule$source_claim_entity_role ==
    "compatibility_alias"
  testthat::expect_true(all(
    sources$supermodule$support_class[aliases] == "none"
  ))
  testthat::expect_true(all(is.na(
    sources$supermodule$tier_specific_fdr[aliases]
  )))
  testthat::expect_true(all(
    sources$supermodule$statistical_support_status[aliases] ==
      "inherited_from_canonical_entity"
  ))
  testthat::expect_true(all(is.na(
    sources$supermodule$q_local_le_0_25[aliases]
  )))
})

testthat::test_that("circular heatmap estimates are copied exactly from Stage 05", {
  source(test_root("R", "paths.R"))
  datasets <- c("neuron_neuropil", "neuron_soma", "microglia")
  source_dir <- file.path(
    "results", "source_data", "10_biological_integration",
    "wgcna_circular_atlas", "global"
  )

  for (level in c("supermodule", "module")) {
    circular <- read_circular_csv(
      source_dir, paste0("wgcna_circular_heatmap_source_", level, ".csv")
    )
    source_rows <- dplyr::bind_rows(lapply(datasets, function(dataset) {
      readr::read_csv(
        path_results(
          "tables", "06_modules_WGCNA", "group_effects", dataset,
          paste0(level, "_group_effects.csv")
        ),
        show_col_types = FALSE, progress = FALSE
      ) |>
        dplyr::filter(.data$effect_scope == "within_spatial_unit")
    }))
    endpoint_column <- paste0(level, "_id")
    source_key <- paste(
      source_rows$dataset, source_rows[[endpoint_column]],
      source_rows$effect_scope, source_rows$spatial_unit,
      source_rows$contrast, sep = "||"
    )
    circular_key <- paste(
      circular$dataset, circular[[endpoint_column]],
      circular$effect_scope, circular$spatial_unit,
      circular$contrast, sep = "||"
    )
    testthat::expect_identical(anyDuplicated(source_key), 0L)
    testthat::expect_identical(anyDuplicated(circular_key), 0L)
    testthat::expect_setequal(circular_key, source_key)
    matched <- match(circular_key, source_key)
    testthat::expect_equal(
      circular$estimate, source_rows$estimate[matched], tolerance = 0
    )
    testthat::expect_equal(
      circular$raw_stage05_estimate,
      source_rows$estimate[matched],
      tolerance = 0
    )
    testthat::expect_equal(
      circular$SE, source_rows$SE[matched], tolerance = 0
    )
    testthat::expect_equal(
      circular$CI_low, source_rows$CI_low[matched], tolerance = 1e-15
    )
    testthat::expect_equal(
      circular$CI_high, source_rows$CI_high[matched], tolerance = 1e-15
    )
    testthat::expect_equal(
      circular$p_value, source_rows$p_value[matched], tolerance = 0
    )
    testthat::expect_equal(
      circular$FDR_local_exploratory,
      source_rows$FDR_local_exploratory[matched],
      tolerance = 0
    )
  }
})

testthat::test_that("outer supermodule support track uses exact Stage 05 global FDR fields", {
  source(test_root("R", "paths.R"))
  source(test_root("R", "wgcna_group_effect_consumer_utils.R"))
  datasets <- c("neuron_neuropil", "neuron_soma", "microglia")
  source_dir <- file.path(
    "results", "source_data", "10_biological_integration",
    "wgcna_circular_atlas", "global"
  )
  circular <- read_circular_csv(
    source_dir, "wgcna_circular_global_supermodule_support_source.csv"
  )
  required <- c(
    "dataset", "supermodule_id", "analysis_tier", "contrast",
    "effect_scope", "independent_hypothesis", "estimate", "SE",
    "CI_low", "CI_high", "p_value", "FDR_primary_global",
    "FDR_secondary_global", "global_support_fdr",
    "global_fdr_source", "tier_specific_fdr",
    "tier_specific_family_id", "tier_specific_family_size",
    "statistical_support_status", "model_valid_for_inference",
    "model_stability_status", "source_claim_entity_role",
    "global_support_class", "global_marker_visible"
  )
  testthat::expect_setequal(
    intersect(required, names(circular)),
    required
  )
  testthat::expect_setequal(unique(circular$dataset), datasets)
  testthat::expect_setequal(
    unique(circular$contrast),
    c("RES - CON", "SUS - CON", "SUS - RES")
  )
  testthat::expect_true(all(
    circular$effect_scope == "spatial_adjusted_global"
  ))
  rows_per_sector <- circular |>
    dplyr::count(.data$dataset, .data$supermodule_id)
  testthat::expect_true(all(rows_per_sector$n == 3L))

  source_rows <- dplyr::bind_rows(lapply(datasets, function(dataset) {
    path <- path_results(
      "tables", "06_modules_WGCNA", "group_effects", dataset,
      "supermodule_group_effects.csv"
    )
    wgcna_group_effect_consumer_adapt(readr::read_csv(
      path, show_col_types = FALSE, progress = FALSE
    )) |>
      dplyr::filter(
        .data$effect_scope == "spatial_adjusted_global",
        .data$contrast %in%
          c("RES - CON", "SUS - CON", "SUS - RES")
      )
  }))
  source_key <- paste(
    source_rows$dataset, source_rows$supermodule_id,
    source_rows$contrast, sep = "||"
  )
  circular_key <- paste(
    circular$dataset, circular$supermodule_id,
    circular$contrast, sep = "||"
  )
  testthat::expect_identical(anyDuplicated(source_key), 0L)
  testthat::expect_identical(anyDuplicated(circular_key), 0L)
  testthat::expect_setequal(circular_key, source_key)
  matched <- match(circular_key, source_key)
  for (column in c(
    "estimate", "SE", "CI_low", "CI_high", "p_value",
    "FDR_primary_global", "FDR_secondary_global",
    "tier_specific_fdr"
  )) {
    testthat::expect_equal(
      circular[[column]], source_rows[[column]][matched],
      tolerance = 0, info = column
    )
  }

  sus_res <- circular$contrast == "SUS - RES"
  contextual <- circular$contrast %in% c("RES - CON", "SUS - CON")
  testthat::expect_true(all(
    circular$analysis_tier[sus_res] == "primary_wgcna_global"
  ))
  testthat::expect_true(all(
    circular$analysis_tier[contextual] ==
      "secondary_contextual_global"
  ))
  testthat::expect_equal(
    circular$global_support_fdr[sus_res],
    circular$FDR_primary_global[sus_res],
    tolerance = 0
  )
  testthat::expect_equal(
    circular$global_support_fdr[contextual],
    circular$FDR_secondary_global[contextual],
    tolerance = 0
  )
  testthat::expect_true(all(
    circular$global_fdr_source[sus_res] == "FDR_primary_global"
  ))
  testthat::expect_true(all(
    circular$global_fdr_source[contextual] ==
      "FDR_secondary_global"
  ))

  invalid <- !(circular$model_valid_for_inference %in% TRUE) |
    circular$model_stability_status == "invalid"
  expected_support <- dplyr::case_when(
    !(circular$independent_hypothesis %in% TRUE) ~ "none",
    invalid ~ "invalid",
    !is.na(circular$global_support_fdr) &
      circular$global_support_fdr <= 0.05 ~ "FDR05",
    !is.na(circular$global_support_fdr) &
      circular$global_support_fdr <= 0.10 ~ "FDR10",
    TRUE ~ "none"
  )
  testthat::expect_identical(
    circular$global_support_class,
    expected_support
  )
  testthat::expect_identical(
    circular$global_marker_visible,
    expected_support != "none" &
      circular$independent_hypothesis %in% TRUE
  )
  aliases <- circular$source_claim_entity_role ==
    "compatibility_alias"
  testthat::expect_true(all(!circular$independent_hypothesis[aliases]))
  testthat::expect_true(all(is.na(
    circular$global_support_fdr[aliases]
  )))
  testthat::expect_true(all(
    circular$global_support_class[aliases] == "none"
  ))
  testthat::expect_true(all(
    !circular$global_marker_visible[aliases]
  ))
})

testthat::test_that("local standardized effects use exact Stage 05 animal-spatial response subsets", {
  source(test_root("R", "paths.R"))
  datasets <- c("neuron_neuropil", "neuron_soma", "microglia")
  source_dir <- file.path(
    "results", "source_data", "10_biological_integration",
    "wgcna_circular_atlas", "global"
  )
  circular <- dplyr::bind_rows(lapply(
    c("supermodule", "module"),
    function(level) {
      read_circular_csv(
        source_dir,
        paste0("wgcna_circular_heatmap_source_", level, ".csv")
      )
    }
  ))
  response_values <- dplyr::bind_rows(lapply(datasets, function(dataset) {
    readr::read_csv(
      path_results(
        "tables", "06_modules_WGCNA", "group_effects", dataset,
        "WGCNA_group_effect_animal_spatial_unit_values.csv"
      ),
      show_col_types = FALSE, progress = FALSE
    ) |>
      dplyr::filter(
        is.finite(.data$eigengene),
        !is.na(.data$StressGroup),
        nzchar(.data$StressGroup)
      )
  }))
  animal_key <- paste(
    response_values$dataset, response_values$level,
    response_values$endpoint_id,
    response_values$canonical_spatial_unit,
    response_values$AnimalID,
    sep = "||"
  )
  testthat::expect_identical(anyDuplicated(animal_key), 0L)

  response_lookup <- response_values |>
    dplyr::group_by(
      .data$dataset, .data$level, .data$endpoint_id,
      .data$canonical_spatial_unit
    ) |>
    dplyr::summarise(
      expected_response_SD = stats::sd(.data$eigengene),
      expected_n_rows = dplyr::n(),
      expected_n_animals = dplyr::n_distinct(.data$AnimalID),
      .groups = "drop"
    )
  response_key <- paste(
    response_lookup$dataset, response_lookup$level,
    response_lookup$endpoint_id,
    response_lookup$canonical_spatial_unit,
    sep = "||"
  )
  circular_key <- paste(
    circular$dataset, circular$standardization_level,
    circular$standardization_endpoint_id, circular$spatial_unit,
    sep = "||"
  )
  matched <- match(circular_key, response_key)
  testthat::expect_false(anyNA(matched))
  testthat::expect_equal(
    circular$response_SD,
    response_lookup$expected_response_SD[matched],
    tolerance = 1e-14
  )
  testthat::expect_equal(
    circular$response_n_rows,
    response_lookup$expected_n_rows[matched],
    tolerance = 0
  )
  testthat::expect_equal(
    circular$response_n_animals,
    response_lookup$expected_n_animals[matched],
    tolerance = 0
  )
  testthat::expect_true(all(
    circular$response_subset_unique_animal_rows
  ))
  testthat::expect_true(all(
    circular$response_n_rows == circular$response_n_animals
  ))
  expected_unclipped <- circular$estimate / circular$response_SD
  expected_display <- pmax(-2.5, pmin(2.5, expected_unclipped))
  testthat::expect_equal(
    circular$standardized_effect_unclipped,
    expected_unclipped,
    tolerance = 1e-14
  )
  testthat::expect_equal(
    circular$standardized_effect_display,
    expected_display,
    tolerance = 1e-14
  )
  testthat::expect_identical(
    circular$effect_clipped_for_display,
    abs(expected_unclipped) > 2.5
  )
  testthat::expect_true(all(circular$display_scale_min == -2.5))
  testthat::expect_true(all(circular$display_scale_max == 2.5))
  testthat::expect_true(all(
    circular$effect_metric ==
      "standardized model effect (SD units)"
  ))

  aliases <- circular$source_claim_entity_role ==
    "compatibility_alias"
  testthat::expect_true(all(
    circular$alias_standardization_status[aliases] ==
      "inherited_from_canonical_module_response"
  ))
  testthat::expect_true(all(grepl(
    "no independent alias denominator",
    circular$standardization_scope[aliases],
    fixed = TRUE
  )))
  testthat::expect_true(all(
    circular$support_class[aliases] == "none"
  ))
})

testthat::test_that("standardized circular SVG and PDF outputs exist", {
  figure_dir <- test_root(
    "results", "figures", "10_biological_integration",
    "wgcna_circular_atlas", "global"
  )
  stems <- c(
    "wgcna_circular_heatmap_all_datasets",
    "wgcna_circular_heatmap_neuron_neuropil",
    "wgcna_circular_heatmap_neuron_soma",
    "wgcna_circular_heatmap_microglia"
  )
  paths <- unlist(lapply(stems, function(stem) {
    file.path(figure_dir, paste0(stem, c(".svg", ".pdf")))
  }))
  testthat::expect_true(all(file.exists(paths)))
  testthat::expect_true(all(file.info(paths)$size > 0))
  canonical_svg <- paste(
    readLines(
      file.path(
        figure_dir,
        "wgcna_circular_heatmap_all_datasets.svg"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  testthat::expect_match(
    canonical_svg,
    "standardized model effect (SD units)",
    fixed = TRUE
  )
  testthat::expect_false(grepl(
    "Cohen's d", canonical_svg, fixed = TRUE
  ))
  testthat::expect_false(grepl(
    "not used", canonical_svg, fixed = TRUE
  ))
  testthat::expect_false(grepl(
    "q <= .25", canonical_svg, fixed = TRUE
  ))
  testthat::expect_match(
    canonical_svg,
    "Global SM effect support",
    fixed = TRUE
  )
  testthat::expect_match(
    canonical_svg,
    "Local spatial BH support",
    fixed = TRUE
  )
  testthat::expect_match(canonical_svg, "RES-CON", fixed = TRUE)
  testthat::expect_match(canonical_svg, "SUS-CON", fixed = TRUE)
  testthat::expect_match(canonical_svg, "SUS-RES", fixed = TRUE)
})

testthat::test_that("circular script does not use broad legacy FDR for support", {
  script <- paste(
    readLines(
      test_root("10_biological_integration", "04_wgcna_circular_atlas.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  testthat::expect_match(
    script,
    'source(repo_path("R", "wgcna_group_effect_consumer_utils.R"))',
    fixed = TRUE
  )
  testthat::expect_false(grepl("source_fdr <- function", script, fixed = TRUE))
  testthat::expect_false(grepl(
    "effect_support_class(p_value, fdr_within, fdr_global)",
    script,
    fixed = TRUE
  ))
  testthat::expect_match(
    script,
    "legacy broad FDR columns are provenance only",
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    'effect_column = "standardized_effect_display"',
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    "fixed_effect_limit = LOCAL_EFFECT_DISPLAY_LIMIT",
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    "Local circular support markers use only FDR_local_exploratory",
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    "global RES-CON and SUS-CON markers use FDR_secondary_global",
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    "global SUS-RES markers use FDR_primary_global",
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    "font = if (isTRUE(compact_local))",
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    "inner_radius = 0.52",
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    'border = "white"',
    fixed = TRUE
  )
  testthat::expect_match(
    script,
    '"alias: inherited / non-independent"',
    fixed = TRUE
  )
  testthat::expect_false(grepl(
    '"not used"', script, fixed = TRUE
  ))
})
