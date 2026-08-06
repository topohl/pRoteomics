testthat::local_edition(3)

source(testthat::test_path(
  "..", "..", "R", "wgcna_group_effect_consumer_utils.R"
))
source(testthat::test_path(
  "..", "..", "R", "wgcna_stage07_semantic_utils.R"
))

wgcna_consumer_fixture <- function() {
  data.frame(
    source_row_id = paste0("row_", seq_len(6L)),
    dataset = rep("microglia", 6L),
    level = c("module", "module", "module", "module", "module", "supermodule"),
    endpoint_id = c("m1", "m1", "m1", "m1", "m1", "sm_alias"),
    effect_scope = c(
      "spatial_adjusted_global",
      "spatial_adjusted_global",
      "stress_by_spatial_interaction",
      "within_spatial_unit",
      "stress_by_spatial_interaction",
      "spatial_adjusted_global"
    ),
    spatial_unit = c(
      "global_spatial_adjusted",
      "global_spatial_adjusted",
      "all_spatial_units",
      "PFC",
      "PFC",
      "global_spatial_adjusted"
    ),
    contrast = c(
      "SUS - RES", "RES - CON", "StressGroup x SpatialUnit",
      "SUS - RES", "SUS - RES", "SUS - RES"
    ),
    test_type = c(
      "named_contrast", "named_contrast", "interaction_omnibus",
      "named_contrast", "conditional_interaction_followup",
      "named_contrast"
    ),
    canonical_claim_entity_id = rep("m1", 6L),
    analysis_tier = c(
      "primary_wgcna_global",
      "secondary_contextual_global",
      "secondary_spatial_heterogeneity",
      "exploratory_spatial_localization",
      "exploratory_spatial_localization",
      "primary_wgcna_global"
    ),
    claim_entity_role = c(
      rep("canonical_module", 5L), "compatibility_alias"
    ),
    support_source_entity_id = rep("m1", 6L),
    independent_hypothesis = c(rep(TRUE, 5L), FALSE),
    statistical_support_status = c(
      "FDR_supported", "suggestive_FDR10", "not_supported",
      "nominal_exploratory", "nominal_exploratory",
      "inherited_from_canonical_entity"
    ),
    model_valid_for_inference = rep(TRUE, 6L),
    model_stability_status = rep("stable_mixed_model", 6L),
    claim_allowed_model = rep(TRUE, 6L),
    primary_model_stable = rep(TRUE, 6L),
    estimate = c(1, 2, NA, 4, 5, 1),
    SE = c(0.1, 0.2, NA, 0.4, 0.5, 0.1),
    p_value = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.01),
    FDR_primary_global = c(0.04, rep(NA_real_, 5L)),
    FDR_primary_family_id = c("primary_family", rep(NA_character_, 5L)),
    n_tests_FDR_primary = c(45L, rep(NA_integer_, 5L)),
    FDR_secondary_global = c(NA_real_, 0.08, rep(NA_real_, 4L)),
    FDR_secondary_family_id = c(
      NA_character_, "secondary_family", rep(NA_character_, 4L)
    ),
    n_tests_FDR_secondary_global = c(
      NA_integer_, 30L, rep(NA_integer_, 4L)
    ),
    FDR_interaction_omnibus = c(
      NA_real_, NA_real_, 0.12, rep(NA_real_, 3L)
    ),
    FDR_interaction_family_id = c(
      NA_character_, NA_character_, "interaction_family",
      rep(NA_character_, 3L)
    ),
    n_tests_FDR_interaction_omnibus = c(
      NA_integer_, NA_integer_, 45L, rep(NA_integer_, 3L)
    ),
    FDR_local_exploratory = c(
      rep(NA_real_, 3L), 0.20, NA_real_, NA_real_
    ),
    FDR_local_family_id = c(
      rep(NA_character_, 3L), "local_family", NA_character_, NA_character_
    ),
    n_tests_FDR_local_exploratory = c(
      rep(NA_integer_, 3L), 180L, NA_integer_, NA_integer_
    ),
    FDR_conservative_all_tests = c(
      0.90, 0.80, 0.70, NA_real_, NA_real_, NA_real_
    ),
    FDR_global = c(
      0.90, 0.80, 0.70, NA_real_, NA_real_, NA_real_
    ),
    FDR_within_dataset_level = rep(NA_real_, 6L),
    FDR_dataset_all_levels = c(
      0.60, 0.50, 0.40, NA_real_, NA_real_, NA_real_
    ),
    FDR_conservative_family_id = rep(NA_character_, 6L),
    n_tests_FDR_conservative_all_tests = rep(NA_integer_, 6L),
    FDR_family_within_level_id = rep(NA_character_, 6L),
    n_tests_FDR_within_dataset_level = rep(NA_integer_, 6L),
    FDR_family_dataset_id = rep(NA_character_, 6L),
    n_tests_FDR_dataset_all_levels = rep(NA_integer_, 6L),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("every authoritative family maps to typed adapter fields", {
  observed <- wgcna_group_effect_consumer_adapt(wgcna_consumer_fixture())

  testthat::expect_equal(
    observed$tier_specific_fdr[1:4],
    c(0.04, 0.08, 0.12, 0.20)
  )
  testthat::expect_identical(
    observed$tier_specific_family_id[1:4],
    c(
      "primary_family", "secondary_family",
      "interaction_family", "local_family"
    )
  )
  testthat::expect_identical(
    observed$tier_specific_family_size[1:4],
    c(45L, 30L, 45L, 180L)
  )
  testthat::expect_identical(
    observed$result_scope[1:4],
    c(
      "primary_global_vulnerability",
      "contextual_global_control",
      "spatial_heterogeneity_omnibus",
      "exploratory_localization"
    )
  )
  testthat::expect_silent(
    wgcna_group_effect_consumer_validate_adapted(observed)
  )
})

testthat::test_that("alias and conditional precedence suppress independent FDR", {
  observed <- wgcna_group_effect_consumer_adapt(wgcna_consumer_fixture())

  testthat::expect_identical(
    observed$result_scope[5:6],
    c("exploratory_conditional_followup", "compatibility_alias_display")
  )
  testthat::expect_true(all(is.na(observed$tier_specific_fdr[5:6])))
  testthat::expect_true(all(is.na(
    observed$tier_specific_family_id[5:6]
  )))
  testthat::expect_true(all(is.na(
    observed$tier_specific_family_size[5:6]
  )))
  testthat::expect_silent(
    wgcna_group_effect_consumer_validate_adapted(observed)
  )
})

testthat::test_that("adapted validation rejects every derived-field mutation", {
  observed <- wgcna_group_effect_consumer_adapt(wgcna_consumer_fixture())
  cases <- c(
    primary = 1L,
    interaction_omnibus = 3L,
    conditional_followup = 5L,
    compatibility_alias = 6L
  )
  columns <- wgcna_group_effect_consumer_added_columns()

  for (case in names(cases)) {
    row <- cases[[case]]
    for (column in columns) {
      mutated <- observed
      if (column == "tier_specific_fdr") {
        mutated[[column]][row] <- if (
          is.na(mutated[[column]][row])
        ) 0.01 else mutated[[column]][row] + 0.001
      } else if (column == "tier_specific_family_id") {
        mutated[[column]][row] <- "mutated_family"
      } else if (column == "tier_specific_family_size") {
        mutated[[column]][row] <- if (
          is.na(mutated[[column]][row])
        ) 1L else mutated[[column]][row] + 1L
      } else {
        mutated[[column]][row] <- "mutated_scope"
      }

      testthat::expect_error(
        wgcna_group_effect_consumer_validate_adapted(mutated),
        paste0("column ", column, " does not exactly match"),
        info = paste(case, column)
      )
    }
  }

  primary_mutation <- observed
  primary_mutation$tier_specific_fdr[1] <- 0.041
  testthat::expect_error(
    wgcna_group_effect_consumer_validate_primary_selection(
      primary_mutation[1, , drop = FALSE]
    ),
    "tier_specific_fdr does not exactly match"
  )
  testthat::expect_error(
    wgcna_group_effect_consumer_select_primary(primary_mutation),
    "tier_specific_fdr does not exactly match"
  )
})

testthat::test_that("adapter preserves source rows, order, values, and types", {
  input <- wgcna_consumer_fixture()
  source_names <- names(input)
  source_types <- vapply(input, typeof, character(1))
  source_rows <- row.names(input)

  observed <- wgcna_group_effect_consumer_adapt(input)

  testthat::expect_equal(nrow(observed), nrow(input))
  testthat::expect_identical(
    names(observed),
    c(source_names, wgcna_group_effect_consumer_added_columns())
  )
  testthat::expect_identical(observed[source_names], input)
  testthat::expect_identical(
    vapply(observed[source_names], typeof, character(1)),
    source_types
  )
  testthat::expect_identical(row.names(observed), source_rows)
  testthat::expect_identical(
    observed$statistical_support_status,
    input$statistical_support_status
  )
  testthat::expect_type(observed$tier_specific_fdr, "double")
  testthat::expect_type(observed$tier_specific_family_id, "character")
  testthat::expect_type(observed$tier_specific_family_size, "integer")
  testthat::expect_type(observed$result_scope, "character")
})

testthat::test_that("base adapter does not apply model or support filters", {
  input <- wgcna_consumer_fixture()
  input$model_valid_for_inference[1] <- FALSE
  input$claim_allowed_model[1] <- FALSE
  input$primary_model_stable[1] <- FALSE
  input$statistical_support_status[1] <- "not_supported"

  observed <- wgcna_group_effect_consumer_adapt(input)

  testthat::expect_equal(nrow(observed), nrow(input))
  testthat::expect_false(observed$model_valid_for_inference[1])
  testthat::expect_false(observed$claim_allowed_model[1])
  testthat::expect_false(observed$primary_model_stable[1])
  testthat::expect_identical(
    observed$statistical_support_status[1], "not_supported"
  )
})

testthat::test_that("legacy FDR columns are preserved but never selected", {
  input <- wgcna_consumer_fixture()
  input$FDR_global[1:3] <- c(0.001, 0.002, 0.003)
  input$FDR_dataset_all_levels[1:3] <- c(0.004, 0.005, 0.006)
  input$FDR_conservative_all_tests[1:3] <- c(0.007, 0.008, 0.009)
  testthat::expect_true(all(is.na(input$FDR_within_dataset_level)))

  observed <- wgcna_group_effect_consumer_adapt(input)

  testthat::expect_equal(
    observed$tier_specific_fdr[1:4],
    c(0.04, 0.08, 0.12, 0.20)
  )
  testthat::expect_identical(observed$FDR_global, input$FDR_global)
  testthat::expect_identical(
    observed$FDR_conservative_all_tests,
    input$FDR_conservative_all_tests
  )
})

testthat::test_that("blank CSV family metadata is treated as absent", {
  input <- wgcna_consumer_fixture()
  family_id_columns <- c(
    "FDR_primary_family_id", "FDR_secondary_family_id",
    "FDR_interaction_family_id", "FDR_local_family_id"
  )
  for (column in family_id_columns) {
    input[[column]][is.na(input[[column]])] <- ""
  }

  testthat::expect_silent(wgcna_group_effect_consumer_validate(input))
  observed <- wgcna_group_effect_consumer_adapt(input)
  testthat::expect_equal(
    observed$tier_specific_fdr[1:4],
    c(0.04, 0.08, 0.12, 0.20)
  )
  testthat::expect_true(all(is.na(
    observed$tier_specific_family_id[5:6]
  )))
})

testthat::test_that("all-NA CSV-inferred conditional fields are preserved", {
  input <- wgcna_consumer_fixture()[5, , drop = FALSE]
  inferred_logical <- c(
    "statistical_support_status",
    "FDR_primary_family_id", "FDR_secondary_family_id",
    "FDR_interaction_family_id", "FDR_local_family_id",
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "n_tests_FDR_primary", "n_tests_FDR_secondary_global",
    "n_tests_FDR_interaction_omnibus",
    "n_tests_FDR_local_exploratory"
  )
  for (column in inferred_logical) {
    input[[column]] <- NA
  }
  source_types <- vapply(input, typeof, character(1))

  observed <- wgcna_group_effect_consumer_adapt(input)

  testthat::expect_identical(
    vapply(observed[names(input)], typeof, character(1)),
    source_types
  )
  testthat::expect_identical(
    observed$result_scope, "exploratory_conditional_followup"
  )
  testthat::expect_true(is.na(observed$statistical_support_status))
})

wgcna_inferential_handoff_fixture <- function() {
  input <- wgcna_consumer_fixture()
  primary <- input[1, , drop = FALSE]
  module_rows <- do.call(rbind, lapply(seq_len(3L), function(i) {
    row <- primary
    row$source_row_id <- paste0("module_", i)
    row$endpoint_id <- paste0("m", i)
    row$canonical_claim_entity_id <- paste0("m", i)
    row$support_source_entity_id <- paste0("m", i)
    row
  }))
  alias <- input[6, , drop = FALSE]
  alias$source_row_id <- "singleton_alias"
  alias$endpoint_id <- "SM01"
  alias$canonical_claim_entity_id <- "m1"
  alias$support_source_entity_id <- "m1"
  block <- primary
  block$source_row_id <- "multi_supermodule"
  block$level <- "supermodule"
  block$endpoint_id <- "SM02"
  block$canonical_claim_entity_id <- "SM02"
  block$claim_entity_role <- "higher_order_block"
  block$support_source_entity_id <- "SM02"

  module_rows$module_id <- module_rows$endpoint_id
  module_rows$supermodule_id <- NA_character_
  super_rows <- rbind(alias, block)
  super_rows$module_id <- NA_character_
  super_rows$supermodule_id <- super_rows$endpoint_id
  module_effects <- wgcna_group_effect_consumer_adapt(module_rows)
  supermodule_effects <- wgcna_group_effect_consumer_adapt(super_rows)
  for (data_name in c("module_effects", "supermodule_effects")) {
    data <- get(data_name)
    data$CI_low <- data$estimate - 1.96 * data$SE
    data$CI_high <- data$estimate + 1.96 * data$SE
    data$n_animals_total <- 9L
    data$interpretation_sentence <- paste("Safe interpretation", data$endpoint_id)
    assign(data_name, data)
  }
  module_interpretable <- module_effects
  module_interpretable$ModulePlotLabel <- paste(
    "Module", module_interpretable$module_id
  )
  module_interpretable$microenvironment_class <- "contrast_blind_module_class"
  supermodule_interpretable <- supermodule_effects
  supermodule_interpretable$Supermodule_PlotLabel <- paste(
    "Supermodule", supermodule_interpretable$supermodule_id
  )
  supermodule_interpretable$dominant_microenvironment_class <-
    "contrast_blind_supermodule_class"
  membership <- data.frame(
    ModuleID = c("m1", "m2", "m3"),
    SupermoduleID = c("SM01", "SM02", "SM02"),
    stringsAsFactors = FALSE
  )
  list(
    module_effects = module_effects,
    supermodule_effects = supermodule_effects,
    module_interpretable = module_interpretable,
    supermodule_interpretable = supermodule_interpretable,
    membership = membership
  )
}

testthat::test_that(
  "inferential handoff keeps every module and only genuine supermodules",
  {
    fixture <- wgcna_inferential_handoff_fixture()
    observed <- do.call(
      wgcna_stage07_build_inferential_handoff, fixture
    )

    testthat::expect_setequal(
      observed$entity_id[observed$entity_level == "module"],
      c("m1", "m2", "m3")
    )
    testthat::expect_true(
      "m1" %in% observed$entity_id[observed$entity_level == "module"]
    )
    testthat::expect_false("SM01" %in% observed$entity_id)
    testthat::expect_true("SM02" %in% observed$entity_id)
    testthat::expect_equal(
      observed$n_member_modules[observed$entity_id == "SM02"], 2L
    )
    testthat::expect_false(any(
      observed$claim_entity_role == "compatibility_alias"
    ))
    testthat::expect_true(all(observed$independent_hypothesis))
    testthat::expect_false(any(
      observed$test_type == "conditional_interaction_followup"
    ))
    testthat::expect_true(all(
      observed$tier_specific_family_size == 45L
    ))
    testthat::expect_equal(
      observed$tier_specific_fdr,
      observed$FDR_primary_global,
      tolerance = 0
    )
    testthat::expect_identical(
      observed$source_key,
      wgcna_inferential_handoff_source_key(observed)
    )
    testthat::expect_identical(
      observed$source_artifact,
      wgcna_inferential_handoff_source_artifact(
        observed$dataset, observed$entity_level
      )
    )
    testthat::expect_silent(
      wgcna_stage07_validate_inferential_handoff(observed)
    )
  }
)

testthat::test_that(
  "annotation changes cannot alter inferential row selection",
  {
    fixture <- wgcna_inferential_handoff_fixture()
    first <- do.call(wgcna_stage07_build_inferential_handoff, fixture)
    fixture$module_interpretable$microenvironment_class <- rev(
      fixture$module_interpretable$microenvironment_class
    )
    fixture$module_interpretable$ModulePlotLabel <- paste0(
      "Changed ", seq_len(nrow(fixture$module_interpretable))
    )
    second <- do.call(wgcna_stage07_build_inferential_handoff, fixture)
    key <- c(
      "dataset", "entity_level", "entity_id", "analysis_tier", "contrast",
      "effect_scope", "spatial_unit", "test_type"
    )
    testthat::expect_identical(first[key], second[key])
    testthat::expect_identical(
      first$tier_specific_fdr, second$tier_specific_fdr
    )
  }
)

testthat::test_that(
  "singleton supermodule display reuses one canonical module endpoint",
  {
    fixture <- wgcna_inferential_handoff_fixture()
    handoff <- do.call(
      wgcna_stage07_build_inferential_handoff, fixture
    )
    lookup <- data.frame(
      dataset = "microglia",
      module_id = c("m1", "m2", "m3"),
      supermodule_id = c("SM01", "SM02", "SM02"),
      n_member_modules = c(1L, 2L, 2L),
      supermodule_label = c("Singleton", "Block", "Block"),
      stringsAsFactors = FALSE
    )
    display <- wgcna_inferential_handoff_supermodule_display(
      handoff, lookup
    )
    singleton <- display[display$supermodule_id == "SM01", , drop = FALSE]
    testthat::expect_equal(nrow(singleton), 1L)
    testthat::expect_identical(singleton$source_entity_level, "module")
    testthat::expect_identical(singleton$source_entity_id, "m1")
    testthat::expect_true(singleton$display_is_compatibility_alias)
    testthat::expect_identical(
      singleton$source_key,
      handoff$source_key[handoff$entity_id == "m1"]
    )
    testthat::expect_false(any(
      handoff$claim_entity_role == "compatibility_alias"
    ))
  }
)

testthat::test_that(
  "migrated manuscript consumers contain no legacy broad-FDR decisions",
  {
    root <- testthat::test_path("..", "..")
    consumers <- c(
      "06_modules_WGCNA/08_wgcna_publication_figures.R",
      "06_modules_WGCNA/08b_microglia_wgcna_readiness_publication_figures.R",
      "06_modules_WGCNA/09_microglia_neuropil_independence.R",
      "06_modules_WGCNA/11_module_robustness_sensitivity.r",
      "06_modules_WGCNA/13_wgcna_claim_readiness.R",
      "08_behavior_physio_coupling/03_module_behavior_coupling.r",
      "09_export_pride_journal/07_make_biological_claims_table.R",
      "10_biological_integration/01_cross_compartment_program_atlas.r",
      "10_biological_integration/04_wgcna_cross_compartment_overview.R",
      "10_biological_integration/04_wgcna_circular_atlas.R",
      "R/final_evidence_bundle_utils.R"
    )
    for (consumer in consumers) {
      text <- readLines(file.path(root, consumer), warn = FALSE)
      testthat::expect_true(
        any(grepl("WGCNA_inferential_handoff", text, fixed = TRUE)),
        info = consumer
      )
      testthat::expect_false(any(grepl(
        paste(
          "module_group_effects\\.csv|supermodule_group_effects\\.csv|",
          "FDR_primary_global|FDR_secondary_global|",
          "FDR_interaction_omnibus|FDR_local_exploratory|",
          "FDR_global|FDR_within_dataset_level",
          sep = ""
        ),
        text
      )), info = consumer)
    }
    semantic <- readLines(
      file.path(root, "R", "wgcna_stage07_semantic_utils.R"),
      warn = FALSE
    )
    testthat::expect_false(any(grepl("p\\.adjust\\s*\\(", semantic)))
  }
)

testthat::test_that("source and semantic keys are explicit and validated", {
  testthat::expect_identical(
    wgcna_group_effect_consumer_source_key(),
    c(
      "dataset", "level", "endpoint_id", "effect_scope",
      "spatial_unit", "contrast", "test_type"
    )
  )
  testthat::expect_identical(
    wgcna_group_effect_consumer_semantic_key(),
    c(
      "dataset", "level", "canonical_claim_entity_id", "analysis_tier",
      "contrast", "effect_scope", "spatial_unit", "test_type"
    )
  )

  input <- wgcna_consumer_fixture()
  duplicate_source <- rbind(input, input[1, , drop = FALSE])
  testthat::expect_error(
    wgcna_group_effect_consumer_validate(duplicate_source),
    "source keys are duplicated"
  )

  duplicate_semantic <- rbind(input, input[1, , drop = FALSE])
  duplicate_semantic$endpoint_id[nrow(duplicate_semantic)] <- "m2"
  testthat::expect_error(
    wgcna_group_effect_consumer_validate(duplicate_semantic),
    "semantic keys are duplicated"
  )
})

testthat::test_that("canonical claim IDs need not be globally unique", {
  input <- wgcna_consumer_fixture()
  testthat::expect_gt(sum(input$canonical_claim_entity_id == "m1"), 1L)
  testthat::expect_silent(wgcna_group_effect_consumer_validate(input))
  testthat::expect_silent(wgcna_group_effect_consumer_adapt(input))
})

testthat::test_that("primary claim selection is explicit and scoped", {
  observed <- wgcna_group_effect_consumer_adapt(wgcna_consumer_fixture())
  selected <- wgcna_group_effect_consumer_select_primary(observed)

  testthat::expect_identical(selected$source_row_id, "row_1")
  testthat::expect_true(all(selected$independent_hypothesis))
  testthat::expect_true(all(
    selected$result_scope == "primary_global_vulnerability"
  ))
  testthat::expect_silent(
    wgcna_group_effect_consumer_validate_primary_selection(selected)
  )

  alias <- observed[observed$claim_entity_role == "compatibility_alias", ]
  testthat::expect_error(
    wgcna_group_effect_consumer_validate_primary_selection(alias),
    "only independent SUS - RES primary"
  )

  duplicate_claim <- rbind(
    wgcna_consumer_fixture(),
    wgcna_consumer_fixture()[1, , drop = FALSE]
  )
  duplicate_claim$endpoint_id[nrow(duplicate_claim)] <- "m2"
  duplicate_claim$test_type[nrow(duplicate_claim)] <- "named_contrast_variant"
  duplicate_claim <- wgcna_group_effect_consumer_adapt(duplicate_claim)
  testthat::expect_error(
    wgcna_group_effect_consumer_select_primary(duplicate_claim),
    "duplicate scoped canonical keys"
  )
})

testthat::test_that("invalid cross-family and missing family rows fail closed", {
  input <- wgcna_consumer_fixture()
  input$FDR_secondary_global[1] <- 0.05
  input$FDR_secondary_family_id[1] <- "wrong_family"
  input$n_tests_FDR_secondary_global[1] <- 30L
  testthat::expect_error(
    wgcna_group_effect_consumer_validate(input),
    "cannot populate another tier-specific FDR family"
  )

  contracts <- list(
    c(1L, "FDR_primary_global"),
    c(2L, "FDR_secondary_global"),
    c(3L, "FDR_interaction_omnibus"),
    c(4L, "FDR_local_exploratory")
  )
  for (contract in contracts) {
    invalid <- wgcna_consumer_fixture()
    row <- as.integer(contract[[1]])
    column <- contract[[2]]
    invalid[[column]][row] <- NA_real_
    testthat::expect_error(
      wgcna_group_effect_consumer_validate(invalid),
      "require a finite",
      info = column
    )
  }

  missing_id <- wgcna_consumer_fixture()
  missing_id$FDR_primary_family_id[1] <- NA_character_
  testthat::expect_error(
    wgcna_group_effect_consumer_validate(missing_id),
    "require FDR_primary_family_id"
  )

  missing_size <- wgcna_consumer_fixture()
  missing_size$n_tests_FDR_primary[1] <- NA_integer_
  testthat::expect_error(
    wgcna_group_effect_consumer_validate(missing_size),
    "require a positive integer-valued n_tests_FDR_primary"
  )
})

testthat::test_that("special rows cannot carry independent q-values or metadata", {
  q_columns <- c(
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "FDR_conservative_all_tests", "FDR_global",
    "FDR_within_dataset_level", "FDR_dataset_all_levels"
  )
  metadata_values <- list(
    FDR_conservative_family_id = "conservative_family",
    n_tests_FDR_conservative_all_tests = 45L,
    FDR_family_within_level_id = "within_level_family",
    n_tests_FDR_within_dataset_level = 45L,
    FDR_family_dataset_id = "dataset_family",
    n_tests_FDR_dataset_all_levels = 90L
  )

  for (column in q_columns) {
    alias <- wgcna_consumer_fixture()
    alias[[column]][6] <- 0.04
    testthat::expect_error(
      wgcna_group_effect_consumer_validate(alias),
      "Compatibility aliases cannot carry",
      info = paste("alias", column)
    )

    conditional <- wgcna_consumer_fixture()
    conditional[[column]][5] <- 0.04
    testthat::expect_error(
      wgcna_group_effect_consumer_validate(conditional),
      "Conditional interaction follow-ups cannot carry",
      info = paste("conditional", column)
    )
  }

  for (column in names(metadata_values)) {
    alias <- wgcna_consumer_fixture()
    alias[[column]][6] <- metadata_values[[column]]
    testthat::expect_error(
      wgcna_group_effect_consumer_validate(alias),
      "Compatibility aliases cannot carry",
      info = paste("alias", column)
    )

    conditional <- wgcna_consumer_fixture()
    conditional[[column]][5] <- metadata_values[[column]]
    testthat::expect_error(
      wgcna_group_effect_consumer_validate(conditional),
      "Conditional interaction follow-ups cannot carry",
      info = paste("conditional", column)
    )
  }
})

testthat::test_that("adapted CSV round-trip validates and selects primary rows", {
  observed <- wgcna_group_effect_consumer_adapt(wgcna_consumer_fixture())
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(observed, path, row.names = FALSE, na = "")
  roundtrip <- utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = ""
  )

  testthat::expect_silent(
    wgcna_group_effect_consumer_validate_adapted(roundtrip)
  )
  selected <- wgcna_group_effect_consumer_select_primary(roundtrip)
  testthat::expect_identical(selected$source_row_id, "row_1")
  testthat::expect_identical(
    selected$tier_specific_fdr,
    selected$FDR_primary_global
  )
})

testthat::test_that("typed zero-row adapter and selections are valid", {
  empty <- wgcna_consumer_fixture()[0, , drop = FALSE]
  observed <- wgcna_group_effect_consumer_adapt(empty)

  testthat::expect_equal(nrow(observed), 0L)
  testthat::expect_type(observed$tier_specific_fdr, "double")
  testthat::expect_type(observed$tier_specific_family_id, "character")
  testthat::expect_type(observed$tier_specific_family_size, "integer")
  testthat::expect_type(observed$result_scope, "character")
  testthat::expect_silent(
    wgcna_group_effect_consumer_validate_adapted(observed)
  )
  testthat::expect_silent(
    wgcna_group_effect_consumer_validate_primary_selection(observed)
  )
  testthat::expect_equal(
    nrow(wgcna_group_effect_consumer_select_primary(observed)), 0L
  )
})

testthat::test_that("zero supported primary results are valid", {
  input <- wgcna_consumer_fixture()
  input$statistical_support_status[
    input$analysis_tier == "primary_wgcna_global" &
      input$independent_hypothesis
  ] <- "not_supported"
  testthat::expect_true(all(is.na(input$FDR_within_dataset_level)))

  observed <- wgcna_group_effect_consumer_adapt(input)

  testthat::expect_equal(
    nrow(wgcna_group_effect_consumer_select_primary(
      observed, "FDR_supported"
    )),
    0L
  )
  testthat::expect_equal(
    nrow(wgcna_group_effect_consumer_select_primary(
      observed, "suggestive_FDR10"
    )),
    0L
  )
  testthat::expect_equal(
    nrow(wgcna_group_effect_consumer_select_primary(
      observed, "optional_status_not_present"
    )),
    0L
  )
})

testthat::test_that("adapter source contains no statistical recomputation", {
  helper <- readLines(testthat::test_path(
    "..", "..", "R", "wgcna_group_effect_consumer_utils.R"
  ), warn = FALSE)
  testthat::expect_false(any(grepl(
    "p\\.adjust\\s*\\(",
    helper
  )))
  testthat::expect_false(any(grepl(
    "statistical_support_status\\s*<-",
    helper
  )))
  testthat::expect_false(any(grepl(
    "evidence_status\\s*<-",
    helper
  )))
})
