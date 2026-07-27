testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"))
source(testthat::test_path(
  "..", "..", "R", "wgcna_group_effect_consumer_utils.R"
))
source(testthat::test_path(
  "..", "..", "R", "wgcna_stage07_semantic_utils.R"
))

wave1a_effects <- function(dataset = "microglia", level = "module") {
  path <- testthat::test_path(
    "..", "..", "results", "tables", "06_modules_WGCNA",
    "group_effects", dataset, paste0(level, "_group_effects.csv")
  )
  testthat::skip_if_not(file.exists(path), "Committed Stage 05 v5 fixture absent.")
  wgcna_group_effect_consumer_adapt(readr::read_csv(
    path, show_col_types = FALSE
  ))
}

wave1a_conditional <- function(dataset = "microglia") {
  path <- testthat::test_path(
    "..", "..", "results", "tables", "06_modules_WGCNA",
    "group_effects", dataset,
    "WGCNA_group_effect_interaction_conditional_followup.csv"
  )
  testthat::skip_if_not(file.exists(path), "Committed conditional fixture absent.")
  wgcna_group_effect_consumer_adapt(readr::read_csv(
    path, show_col_types = FALSE
  ))
}

wave1a_entity_rows <- function(dataset = "microglia", level = "module") {
  effects <- wave1a_effects(dataset, level)
  primary <- wgcna_group_effect_consumer_select_primary(effects)
  testthat::expect_gt(nrow(primary), 0L)
  entity <- primary$canonical_claim_entity_id[[1]]
  effects <- effects[
    effects$canonical_claim_entity_id == entity,
    ,
    drop = FALSE
  ]
  effects$module_biological_label <- paste("Reviewed program", entity)
  conditional <- wave1a_conditional(dataset)
  conditional <- conditional[
    conditional$level == level &
      conditional$canonical_claim_entity_id == entity,
    ,
    drop = FALSE
  ]
  list(effects = effects, conditional = conditional, entity = entity)
}

wave1a_set_support <- function(
    effects, primary = "not_supported", interaction = "not_supported"
) {
  primary_row <- effects$analysis_tier == "primary_wgcna_global" &
    effects$contrast == "SUS - RES" &
    effects$effect_scope == "spatial_adjusted_global" &
    effects$spatial_unit == "global_spatial_adjusted" &
    effects$test_type == "named_contrast"
  interaction_row <-
    effects$analysis_tier == "secondary_spatial_heterogeneity" &
    effects$test_type == "interaction_omnibus"
  effects$statistical_support_status[primary_row] <- primary
  effects$statistical_support_status[interaction_row] <- interaction
  effects
}

wave1a_stage06_guard <- function() {
  expressions <- parse(testthat::test_path(
    "..", "..", "06_modules_WGCNA",
    "06_annotate_module_microenvironment.r"
  ))
  is_guard <- vapply(expressions, function(expression) {
    is.call(expression) &&
      identical(expression[[1]], as.name("<-")) &&
      identical(
        as.character(expression[[2]]),
        "wgcna_stage06_validate_contrast_blind_annotation"
      )
  }, logical(1))
  testthat::expect_equal(sum(is_guard), 1L)
  environment <- new.env(parent = baseenv())
  eval(expressions[[which(is_guard)]], envir = environment)
  environment$wgcna_stage06_validate_contrast_blind_annotation
}

testthat::test_that("Stage 06 is contrast-blind and rejects statistical fields", {
  script <- paste(readLines(testthat::test_path(
    "..", "..", "06_modules_WGCNA",
    "06_annotate_module_microenvironment.r"
  ), warn = FALSE), collapse = "\n")

  testthat::expect_false(grepl(
    paste(
      "module_group_effects\\.csv",
      "supermodule_group_effects\\.csv",
      "\\$FDR_global",
      "\\$FDR_within_dataset_level",
      "\\$tier_specific_fdr",
      "\\$statistical_support_status",
      "module_annot\\$Changed_in_group_contrasts",
      sep = "|"
    ),
    script
  ))

  annotation <- data.frame(
    dataset = "microglia",
    ModuleID = "WGCNA_#123456",
    biological_label = "Microglia-associated ROI program",
    microenvironment_class = "microglia_supported",
    annotation_confidence = "high",
    stringsAsFactors = FALSE
  )
  mutated_effects <- data.frame(
    contrast = sample(c("SUS - RES", "RES - CON")),
    p_value = c(1e-12, 0.99),
    FDR_global = c(1e-10, 1),
    stringsAsFactors = FALSE
  )
  shuffled_effects <- mutated_effects[rev(seq_len(nrow(
    mutated_effects
  ))), , drop = FALSE]

  annotation_after_support_mutation <- annotation
  testthat::expect_identical(annotation_after_support_mutation, annotation)
  testthat::expect_false(identical(mutated_effects, shuffled_effects))
  testthat::expect_match(
    script, "wgcna_stage06_validate_contrast_blind_annotation"
  )

  guard <- wave1a_stage06_guard()
  prohibited_examples <- list(
    family_specific_FDR = list(FDR_primary_global = 0.01),
    family_provenance = list(FDR_primary_family_id = "primary:test"),
    family_size = list(n_tests_FDR_primary = 10L),
    tier_scope = list(analysis_tier = "primary_wgcna_global"),
    result_scope = list(result_scope = "primary_global_vulnerability"),
    model_validity = list(model_valid_for_inference = TRUE),
    model_stability = list(model_stability_status = "stable_mixed_model")
  )
  for (category in names(prohibited_examples)) {
    candidate <- annotation
    field <- names(prohibited_examples[[category]])
    candidate[[field]] <- prohibited_examples[[category]][[field]]
    testthat::expect_error(
      guard(candidate),
      "contains group-effect-derived field",
      info = category
    )
  }

  targeted_signature_annotation <- transform(
    annotation,
    contrast_class = "stress_context",
    comparison = "curated signature versus module",
    signature_source = "targeted_panel",
    padj = 0.03,
    NES = 1.4
  )
  testthat::expect_silent(guard(targeted_signature_annotation))
})

testthat::test_that("Stage 07 uses the adapter without legacy inference", {
  script <- paste(readLines(testthat::test_path(
    "..", "..", "06_modules_WGCNA",
    "07_wgcna_interpretable_summary.r"
  ), warn = FALSE), collapse = "\n")

  testthat::expect_true(grepl(
    'source\\(repo_path\\("R", "wgcna_group_effect_consumer_utils\\.R"\\)\\)',
    script
  ))
  testthat::expect_gte(
    lengths(regmatches(
      script,
      gregexpr(
        "wgcna_group_effect_consumer_adapt\\(",
        script, fixed = FALSE
      )
    )),
    3L
  )
  testthat::expect_false(grepl(
    "FDR_global|FDR_within_dataset_level|p\\.adjust\\(",
    script
  ))
  testthat::expect_false(grepl(
    "slice_min\\([^)]*p_value|arrange\\([^)]*p_value|best_p",
    script
  ))
})

testthat::test_that("Stage 07 annotation joins preserve every Stage 05 value", {
  for (level in c("module", "supermodule")) {
    effects <- wave1a_effects("microglia", level)
    id_column <- if (level == "module") "module_id" else "supermodule_id"
    annotation_id <- if (level == "module") "ModuleID" else "SupermoduleID"
    annotation <- data.frame(
      dataset = unique(effects$dataset),
      entity = unique(effects[[id_column]]),
      biological_label = paste("Reviewed", unique(effects[[id_column]])),
      stringsAsFactors = FALSE
    )
    names(annotation)[names(annotation) == "entity"] <- annotation_id
    collision <- if (level == "module") "module_label" else "supermodule_label"
    annotation[[collision]] <- paste("Annotation", annotation[[annotation_id]])

    observed <- wgcna_stage07_join_annotations(
      effects,
      annotation,
      by = stats::setNames(c("dataset", annotation_id), c(
        "dataset", id_column
      )),
      artifact = paste(level, "test join")
    )
    testthat::expect_identical(observed[names(effects)], effects)
    testthat::expect_true(
      paste0(collision, "_annotation") %in% names(observed)
    )
  }
})

testthat::test_that("spatial class uses only primary and omnibus FDR support", {
  fixture <- wave1a_entity_rows()
  effects <- fixture$effects
  effects$p_value[] <- 1e-12
  effects$estimate <- rep(c(-2, 2), length.out = nrow(effects))

  cases <- list(
    no_support = c("not_supported", "not_supported"),
    primary = c("FDR_supported", "not_supported"),
    interaction = c("not_supported", "FDR_supported"),
    both = c("FDR_supported", "FDR_supported")
  )
  expected <- c(
    no_support = "no_FDR_supported_group_organization",
    primary = "primary_global_shift",
    interaction =
      "group_dependent_spatial_heterogeneity_without_primary_global_shift",
    both =
      "primary_global_shift_with_group_dependent_spatial_heterogeneity"
  )
  for (case in names(cases)) {
    selected <- wave1a_set_support(
      effects, cases[[case]][[1]], cases[[case]][[2]]
    )
    observed <- wgcna_stage07_build_spatial_organization(
      selected, fixture$conditional
    )
    testthat::expect_identical(
      observed$spatial_organization_class[[1]],
      expected[[case]],
      info = case
    )
  }
})

testthat::test_that("omnibus attribution and localization remain conservative", {
  fixture <- wave1a_entity_rows()
  effects <- wave1a_set_support(
    fixture$effects, "not_supported", "FDR_supported"
  )
  local <- effects$analysis_tier == "exploratory_spatial_localization"
  effects$estimate[local] <- rep(c(-1, 1), length.out = sum(local))
  effects$statistical_support_status[local] <- "nominal_exploratory"

  observed <- wgcna_stage07_build_spatial_organization(
    effects, fixture$conditional
  )
  testthat::expect_identical(
    observed$spatial_heterogeneity_contrast_attribution,
    "not_attributed_to_specific_contrast"
  )
  testthat::expect_true(observed$localization_evidence_available)
  testthat::expect_match(
    observed$classification_reason,
    "not classified as model instability"
  )
  testthat::expect_true(all(is.na(
    fixture$conditional$tier_specific_fdr
  )))
  testthat::expect_true(all(
    fixture$conditional$result_scope ==
      "exploratory_conditional_followup"
  ))
})

testthat::test_that("local support attributes contrast only with omnibus support", {
  fixture <- wave1a_entity_rows()
  local <- fixture$effects$analysis_tier ==
    "exploratory_spatial_localization" &
    fixture$effects$test_type != "conditional_interaction_followup"
  local_contrasts <- sort(unique(as.character(
    fixture$effects$contrast[local]
  )))
  testthat::expect_gte(length(local_contrasts), 2L)

  one_contrast <- fixture$effects
  one_contrast$statistical_support_status[local] <- "not_supported"
  one_contrast$statistical_support_status[
    local & one_contrast$contrast == local_contrasts[[1]]
  ] <- "FDR_supported"

  unsupported_omnibus <- wave1a_set_support(
    one_contrast, "not_supported", "not_supported"
  )
  observed_unsupported <- wgcna_stage07_build_spatial_organization(
    unsupported_omnibus, fixture$conditional
  )
  testthat::expect_identical(
    observed_unsupported$spatial_organization_class,
    "no_FDR_supported_group_organization"
  )
  testthat::expect_identical(
    observed_unsupported$spatial_heterogeneity_contrast_attribution,
    "not_attributed_to_specific_contrast"
  )
  testthat::expect_match(
    observed_unsupported$localization_support_summary,
    "local_FDR_supported=[1-9]"
  )
  invalid_without_omnibus <- observed_unsupported
  invalid_without_omnibus$spatial_heterogeneity_contrast_attribution <-
    "attributed_by_FDR_supported_localization:SUS-RES"
  testthat::expect_error(
    wgcna_stage07_validate_spatial_organization(invalid_without_omnibus),
    "unless the independent interaction omnibus is FDR_supported"
  )

  supported_omnibus <- wave1a_set_support(
    one_contrast, "not_supported", "FDR_supported"
  )
  observed_one <- wgcna_stage07_build_spatial_organization(
    supported_omnibus, fixture$conditional
  )
  testthat::expect_identical(
    observed_one$spatial_heterogeneity_contrast_attribution,
    paste0(
      "attributed_by_FDR_supported_localization:",
      gsub("\\s+", "", local_contrasts[[1]])
    )
  )

  multiple_contrasts <- supported_omnibus
  multiple_contrasts$statistical_support_status[
    local & multiple_contrasts$contrast == local_contrasts[[2]]
  ] <- "FDR_supported"
  observed_multiple <- wgcna_stage07_build_spatial_organization(
    multiple_contrasts, fixture$conditional
  )
  testthat::expect_identical(
    observed_multiple$spatial_heterogeneity_contrast_attribution,
    "multiple_contrasts_supported_by_localization"
  )

  no_supported_local <- wave1a_set_support(
    fixture$effects, "not_supported", "FDR_supported"
  )
  no_supported_local$statistical_support_status[local] <- "not_supported"
  conditional_only <- wgcna_stage07_build_spatial_organization(
    no_supported_local, fixture$conditional
  )
  testthat::expect_true(conditional_only$localization_evidence_available)
  testthat::expect_identical(
    conditional_only$spatial_heterogeneity_contrast_attribution,
    "not_attributed_to_specific_contrast"
  )
})

testthat::test_that("spatial validation controls attribution and omnibus type", {
  fixture <- wave1a_entity_rows()
  effects <- wave1a_set_support(
    fixture$effects, "not_supported", "FDR_supported"
  )
  observed <- wgcna_stage07_build_spatial_organization(
    effects, fixture$conditional
  )
  invalid_token <- observed
  invalid_token$spatial_heterogeneity_contrast_attribution <-
    "attributed_by_nominal_followup:SUS-RES"
  testthat::expect_error(
    wgcna_stage07_validate_spatial_organization(invalid_token),
    "invalid spatial-heterogeneity contrast attribution"
  )
  invalid_type <- observed
  invalid_type$interaction_test_type <- "conditional_interaction_followup"
  testthat::expect_error(
    wgcna_stage07_validate_spatial_organization(invalid_type),
    "only independent interaction omnibus rows"
  )
})

testthat::test_that("aliases cannot duplicate spatial findings", {
  effects <- wave1a_effects("microglia", "supermodule")
  primary <- wgcna_group_effect_consumer_select_primary(effects)
  aliases <- effects$claim_entity_role == "compatibility_alias"
  testthat::expect_true(any(aliases))
  testthat::expect_false(any(effects$independent_hypothesis[aliases]))
  testthat::expect_true(all(is.na(effects$tier_specific_fdr[aliases])))
  alias_sentence <- wgcna_stage07_effect_sentence(
    effects[which(aliases)[[1]], , drop = FALSE], "supermodule"
  )
  testthat::expect_match(alias_sentence, "compatibility display alias")
  testthat::expect_match(alias_sentence, "not a separate finding")

  observed <- wgcna_stage07_build_spatial_organization(
    effects, effects[0, , drop = FALSE]
  )
  testthat::expect_equal(nrow(observed), nrow(primary))
  testthat::expect_identical(anyDuplicated(observed[c(
    "dataset", "level", "canonical_claim_entity_id"
  )]), 0L)
})

testthat::test_that("boundary models retain validity but not primary stability", {
  fixture <- wave1a_entity_rows()
  effects <- wave1a_set_support(
    fixture$effects, "FDR_supported", "not_supported"
  )
  primary <- effects$analysis_tier == "primary_wgcna_global"
  effects$model_valid_for_inference[primary] <- TRUE
  effects$claim_allowed_model[primary] <- TRUE
  effects$model_stability_status[primary] <-
    "boundary_random_intercept_zero"
  effects$primary_model_stable[primary] <- FALSE

  observed <- wgcna_stage07_build_spatial_organization(
    effects, fixture$conditional
  )
  sentence <- wgcna_stage07_effect_sentence(
    effects[which(primary)[[1]], , drop = FALSE], "module"
  )
  testthat::expect_identical(
    observed$spatial_organization_class,
    "primary_global_shift"
  )
  testthat::expect_match(sentence, "valid and claim-eligible")
  testthat::expect_match(sentence, "boundary_random_intercept_zero")
})

testthat::test_that("interpretation wording preserves inferential scope", {
  fixture <- wave1a_entity_rows()
  effects <- fixture$effects
  primary <- which(effects$analysis_tier == "primary_wgcna_global")[[1]]
  contextual <- which(
    effects$analysis_tier == "secondary_contextual_global"
  )[[1]]
  interaction <- which(effects$test_type == "interaction_omnibus")[[1]]
  local <- which(
    effects$analysis_tier == "exploratory_spatial_localization" &
      effects$test_type != "conditional_interaction_followup"
  )[[1]]

  effects$statistical_support_status[primary] <- "not_supported"
  effects$statistical_support_status[interaction] <- "FDR_supported"
  sentences <- c(
    wgcna_stage07_effect_sentence(effects[primary, , drop = FALSE], "module"),
    wgcna_stage07_effect_sentence(
      effects[contextual, , drop = FALSE], "module"
    ),
    wgcna_stage07_effect_sentence(
      effects[interaction, , drop = FALSE], "module"
    ),
    wgcna_stage07_effect_sentence(effects[local, , drop = FALSE], "module"),
    wgcna_stage07_effect_sentence(
      fixture$conditional[1, , drop = FALSE], "module"
    )
  )
  testthat::expect_match(sentences[[1]], "not supported after correction")
  testthat::expect_match(
    sentences[[2]], "does not establish equivalence or normalization"
  )
  testthat::expect_match(
    sentences[[3]], "responsible contrast was not established"
  )
  testthat::expect_match(sentences[[4]], "exploratory localization")
  testthat::expect_match(sentences[[5]], "no independent q-value")
  testthat::expect_false(any(grepl(
    "restored|rescued|equivalent to controls|returned to baseline",
    sentences, ignore.case = TRUE
  )))
  testthat::expect_true(all(grepl(
    "microglia-enriched ROI/local microenvironment",
    sentences, fixed = TRUE
  )))
})

testthat::test_that("omnibus sentences preserve controlled support granularity", {
  fixture <- wave1a_entity_rows()
  interaction <- which(
    fixture$effects$test_type == "interaction_omnibus"
  )[[1]]
  expected <- c(
    FDR_supported =
      "FDR-supported in the prespecified interaction family",
    suggestive_FDR10 =
      "suggestive at 10% FDR, but not an FDR-supported result at 5%",
    nominal_exploratory = "nominal exploratory evidence only",
    not_supported =
      "not supported after correction in the prespecified interaction family"
  )
  for (status in names(expected)) {
    row <- fixture$effects[interaction, , drop = FALSE]
    row$statistical_support_status <- status
    sentence <- wgcna_stage07_effect_sentence(row, "module")
    testthat::expect_match(sentence, expected[[status]], fixed = TRUE)
    testthat::expect_match(
      sentence,
      "responsible contrast was not established by the omnibus test alone",
      fixed = TRUE
    )
  }
})

testthat::test_that("typed zero-row Wave 1A outputs are valid", {
  effects <- wave1a_effects()
  unsupported <- wgcna_group_effect_consumer_select_primary(
    effects, support_status = "FDR_supported"
  )
  empty <- effects[0, , drop = FALSE]
  interpreted <- wgcna_stage07_add_interpretation_sentences(empty)
  spatial <- wgcna_stage07_build_spatial_organization(empty, empty)

  testthat::expect_identical(interpreted$interpretation_sentence, character())
  testthat::expect_equal(nrow(unsupported), 0L)
  testthat::expect_silent(
    wgcna_group_effect_consumer_validate_adapted(interpreted)
  )
  testthat::expect_silent(
    wgcna_stage07_validate_interpretable(interpreted)
  )
  testthat::expect_silent(
    wgcna_stage07_validate_spatial_organization(spatial)
  )
  testthat::expect_identical(spatial$primary_tier_specific_fdr, double())
})
