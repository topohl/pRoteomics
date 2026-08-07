test_root <- function(...) {
  testthat::test_path("..", "..", ...)
}

read_circular_csv <- function(...) {
  path <- test_root(...)
  testthat::skip_if_not(
    file.exists(path),
    paste("Generated circular-atlas output unavailable:", path)
  )
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

testthat::test_that("circular atlas uses the Stage 07 inferential handoff", {
  script <- paste(
    readLines(
      test_root("10_biological_integration", "04_wgcna_circular_atlas.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  testthat::expect_match(
    script, "WGCNA_inferential_handoff.csv", fixed = TRUE
  )
  testthat::expect_match(script, "tier_specific_family_id", fixed = TRUE)
  testthat::expect_match(script, "tier_specific_family_size", fixed = TRUE)
  testthat::expect_match(script, "claim_gate", fixed = TRUE)
  testthat::expect_match(script, "source_artifact", fixed = TRUE)
  testthat::expect_match(script, "source_key", fixed = TRUE)
  testthat::expect_match(
    script, "display_is_independent_endpoint", fixed = TRUE
  )
  testthat::expect_match(script, "display_support_origin", fixed = TRUE)
  prohibited <- c(
    "module_group_effects.csv", "supermodule_group_effects.csv",
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "FDR_global", "FDR_within_dataset_level"
  )
  for (token in prohibited) {
    testthat::expect_false(grepl(token, script, fixed = TRUE), info = token)
  }
})

testthat::test_that("live handoffs reconstruct explicit supermodule display semantics", {
  source(test_root("R", "paths.R"))
  source(test_root("R", "wgcna_group_effect_consumer_utils.R"))
  for (dataset in c("neuron_neuropil", "neuron_soma", "microglia")) {
    handoff <- wgcna_inferential_handoff_read(path_results(
      "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
      "WGCNA_inferential_handoff.csv"
    ))
    lookup <- readr::read_csv(
      path_results(
        "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
        "WGCNA_final_label_lookup.csv"
      ),
      show_col_types = FALSE, progress = FALSE
    )
    super <- lookup |>
      dplyr::filter(.data$level == "supermodule") |>
      dplyr::transmute(
        dataset, supermodule_id = .data$entity_id,
        n_member_modules, supermodule_label = .data$final_plot_label
      )
    display_lookup <- lookup |>
      dplyr::filter(.data$level == "module") |>
      dplyr::transmute(
        dataset, module_id = .data$entity_id,
        supermodule_id = .data$parent_entity_id
      ) |>
      dplyr::left_join(
        super, by = c("dataset", "supermodule_id"),
        relationship = "many-to-one"
      )
    display <- wgcna_inferential_handoff_supermodule_display(
      handoff, display_lookup
    )
    aliases <- display$display_is_compatibility_alias %in% TRUE
    testthat::expect_true(all(!display$display_is_independent_endpoint[aliases]))
    testthat::expect_true(all(
      display$display_support_origin[aliases] == "inherited_canonical_module"
    ))
    testthat::expect_true(all(
      display$display_is_independent_endpoint[!aliases]
    ))
    testthat::expect_true(all(
      display$display_support_origin[!aliases] ==
        "independent_supermodule_endpoint"
    ))
  }
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
  testthat::skip_if_not(
    all(vapply(
      sources[c("supermodule", "module")],
      function(x) all(c(
        "claim_gate", "source_artifact", "source_key",
        "display_is_independent_endpoint", "display_support_origin"
      ) %in% names(x)),
      logical(1)
    )),
    "Generated circular sources predate the current handoff/display provenance schema"
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
    "q_local_le_0_05", "q_local_le_0_10", "q_local_le_0_25",
    "claim_gate", "source_artifact", "source_key",
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
  for (source_name in c("supermodule", "module")) {
    source <- sources[[source_name]]
    q <- source$tier_specific_fdr
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
        c("stable_animal_level_lm", "stable_mixed_model") &
      x$claim_gate == "eligible_for_readiness_assessment"
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

  aliases <- sources$supermodule$display_is_compatibility_alias %in% TRUE
  testthat::expect_true(all(
    sources$supermodule$source_claim_entity_role[aliases] ==
      "canonical_module"
  ))
  testthat::expect_true(all(
    sources$supermodule$source_entity_level[aliases] == "module"
  ))
})

testthat::test_that("circular heatmap estimates retain exact handoff source rows", {
  source(test_root("R", "paths.R"))
  source(test_root("R", "wgcna_group_effect_consumer_utils.R"))
  datasets <- c("neuron_neuropil", "neuron_soma", "microglia")
  source_dir <- file.path(
    "results", "source_data", "10_biological_integration",
    "wgcna_circular_atlas", "global"
  )

  for (level in c("supermodule", "module")) {
    circular <- read_circular_csv(
      source_dir, paste0("wgcna_circular_heatmap_source_", level, ".csv")
    )
    testthat::skip_if_not(
      all(c(
        "source_key", "display_is_independent_endpoint",
        "display_support_origin"
      ) %in% names(circular)),
      "Generated circular heatmap source predates current display provenance"
    )
    source_rows <- dplyr::bind_rows(lapply(datasets, function(dataset) {
      handoff <- wgcna_inferential_handoff_read(
        path_results(
          "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
          "WGCNA_inferential_handoff.csv"
        )
      )
      if (level == "module") {
        handoff[handoff$entity_level == "module", , drop = FALSE]
      } else {
        lookup <- readr::read_csv(
          path_results(
            "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
            "WGCNA_final_label_lookup.csv"
          ),
          show_col_types = FALSE, progress = FALSE
        )
        super <- lookup |>
          dplyr::filter(.data$level == "supermodule") |>
          dplyr::transmute(
            dataset, supermodule_id = .data$entity_id,
            n_member_modules, supermodule_label = .data$final_plot_label
          )
        display_lookup <- lookup |>
          dplyr::filter(.data$level == "module") |>
          dplyr::transmute(
            dataset, module_id = .data$entity_id,
            supermodule_id = .data$parent_entity_id
          ) |>
          dplyr::left_join(
            super, by = c("dataset", "supermodule_id"),
            relationship = "many-to-one"
          )
        wgcna_inferential_handoff_supermodule_display(
          handoff, display_lookup
        )
      }
    })) |>
      dplyr::filter(.data$effect_scope == "within_spatial_unit")
    testthat::expect_identical(anyDuplicated(source_rows$source_key), 0L)
    testthat::expect_identical(anyDuplicated(circular$source_key), 0L)
    testthat::expect_setequal(circular$source_key, source_rows$source_key)
    matched <- match(circular$source_key, source_rows$source_key)
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
      circular$tier_specific_fdr,
      source_rows$tier_specific_fdr[matched],
      tolerance = 0
    )
    testthat::expect_identical(
      circular$tier_specific_family_id,
      source_rows$tier_specific_family_id[matched]
    )
    testthat::expect_identical(
      circular$source_artifact,
      source_rows$source_artifact[matched]
    )
  }
})

testthat::test_that("outer support track uses exact handoff FDR families", {
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
  testthat::skip_if_not(
    all(c(
      "claim_gate", "source_artifact", "source_key",
      "display_is_independent_endpoint", "display_support_origin"
    ) %in% names(circular)),
    "Generated circular global support source predates current display provenance"
  )
  required <- c(
    "dataset", "supermodule_id", "analysis_tier", "contrast",
    "effect_scope", "independent_hypothesis", "estimate", "SE",
    "CI_low", "CI_high", "p_value", "global_support_fdr",
    "global_fdr_source", "tier_specific_fdr",
    "tier_specific_family_id", "tier_specific_family_size",
    "statistical_support_status", "model_valid_for_inference",
    "model_stability_status", "source_claim_entity_role",
    "global_support_class", "global_marker_visible", "claim_gate",
    "source_artifact", "source_key"
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
    handoff <- wgcna_inferential_handoff_read(
      path_results(
        "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
        "WGCNA_inferential_handoff.csv"
      )
    )
    lookup <- readr::read_csv(
      path_results(
        "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
        "WGCNA_final_label_lookup.csv"
      ),
      show_col_types = FALSE, progress = FALSE
    )
    super <- lookup |>
      dplyr::filter(.data$level == "supermodule") |>
      dplyr::transmute(
        dataset, supermodule_id = .data$entity_id,
        n_member_modules, supermodule_label = .data$final_plot_label
      )
    display_lookup <- lookup |>
      dplyr::filter(.data$level == "module") |>
      dplyr::transmute(
        dataset, module_id = .data$entity_id,
        supermodule_id = .data$parent_entity_id
      ) |>
      dplyr::left_join(
        super, by = c("dataset", "supermodule_id"),
        relationship = "many-to-one"
      )
    wgcna_inferential_handoff_supermodule_display(
      handoff, display_lookup
    ) |>
      dplyr::filter(
        .data$effect_scope == "spatial_adjusted_global",
        .data$contrast %in%
          c("RES - CON", "SUS - CON", "SUS - RES")
      )
  }))
  testthat::expect_identical(anyDuplicated(source_rows$source_key), 0L)
  testthat::expect_identical(anyDuplicated(circular$source_key), 0L)
  testthat::expect_setequal(circular$source_key, source_rows$source_key)
  matched <- match(circular$source_key, source_rows$source_key)
  for (column in c(
    "estimate", "SE", "CI_low", "CI_high", "p_value",
    "tier_specific_fdr"
  )) {
    testthat::expect_equal(
      circular[[column]], source_rows[[column]][matched],
      tolerance = 0, info = column
    )
  }

  testthat::expect_equal(
    circular$global_support_fdr,
    circular$tier_specific_fdr,
    tolerance = 0
  )
  testthat::expect_true(all(
    circular$global_fdr_source == "tier_specific_fdr"
  ))

  invalid <- circular$claim_gate !=
    "eligible_for_readiness_assessment" |
    !(circular$model_valid_for_inference %in% TRUE) |
    circular$model_stability_status == "invalid"
  expected_support <- dplyr::case_when(
    !(circular$display_is_independent_endpoint %in% TRUE) ~ "none",
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
      circular$display_is_independent_endpoint %in% TRUE
  )
  aliases <- circular$display_is_compatibility_alias %in% TRUE
  testthat::expect_true(all(
    circular$source_claim_entity_role[aliases] == "canonical_module"
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
  testthat::skip_if_not(
    all(c(
      "display_is_compatibility_alias", "source_entity_level",
      "display_is_independent_endpoint", "display_support_origin"
    ) %in% names(circular)),
    "Generated circular standardized sources predate current display provenance"
  )
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

  aliases <- circular$display_is_compatibility_alias %in% TRUE
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
    circular$source_entity_level[aliases] == "module"
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
  testthat::skip_if_not(
    all(file.exists(paths)),
    "Generated circular-atlas figures are unavailable."
  )
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
    script, "WGCNA_inferential_handoff.csv", fixed = TRUE
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
    "Local circular support markers use only the Stage 07 tier-specific FDR",
    fixed = TRUE
  )
  for (token in c(
    "module_group_effects.csv", "supermodule_group_effects.csv",
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "FDR_global", "FDR_within_dataset_level"
  )) {
    testthat::expect_false(grepl(token, script, fixed = TRUE), info = token)
  }
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
