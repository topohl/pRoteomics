# Stage 07-only semantic handoff helpers.
#
# These helpers consume the Wave 0 adapter output. They never calculate FDR,
# rebuild statistical support, decide readiness, or alter Stage 05 values.

.wgcna_stage07_required_interpretable_columns <- function() {
  c(
    "dataset", "level", "endpoint_id", "canonical_claim_entity_id",
    "claim_entity_role", "support_source_entity_id", "independent_hypothesis",
    "analysis_tier", "tier_specific_fdr", "tier_specific_family_id",
    "tier_specific_family_size", "result_scope", "statistical_support_status",
    "model_valid_for_inference", "model_stability_status",
    "claim_allowed_model", "primary_model_stable", "test_type",
    "effect_scope", "spatial_unit", "contrast", "estimate", "SE",
    "CI_low", "CI_high", "p_value", "interpretation_sentence"
  )
}

wgcna_stage07_validate_interpretable <- function(
    data, artifact = "Stage 07 interpretable output"
) {
  missing <- setdiff(
    .wgcna_stage07_required_interpretable_columns(), names(data)
  )
  if (length(missing)) {
    stop(
      artifact, " is missing required column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  wgcna_group_effect_consumer_validate_adapted(data)
  invisible(TRUE)
}

wgcna_stage07_join_annotations <- function(
    effects, annotation, by, artifact = "Stage 07 annotation join"
) {
  wgcna_group_effect_consumer_validate_adapted(effects)
  if (!is.data.frame(annotation)) {
    stop(artifact, " annotation input must be a data frame.", call. = FALSE)
  }
  if (is.null(names(by)) || any(!nzchar(names(by)))) {
    left_keys <- as.character(by)
    right_keys <- as.character(by)
  } else {
    left_keys <- names(by)
    right_keys <- unname(as.character(by))
  }
  missing <- c(
    setdiff(left_keys, names(effects)),
    setdiff(right_keys, names(annotation))
  )
  if (length(missing)) {
    stop(
      artifact, " is missing join key(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  if (nrow(annotation) && anyDuplicated(annotation[right_keys])) {
    stop(
      artifact,
      " annotation keys are duplicated; a many-to-one join is required.",
      call. = FALSE
    )
  }

  source_columns <- names(effects)
  collisions <- intersect(
    setdiff(names(annotation), right_keys), source_columns
  )
  if (length(collisions)) {
    replacement <- paste0(collisions, "_annotation")
    unavailable <- union(source_columns, names(annotation))
    for (i in seq_along(replacement)) {
      candidate <- replacement[[i]]
      suffix <- 2L
      while (candidate %in% unavailable) {
        candidate <- paste0(replacement[[i]], "_", suffix)
        suffix <- suffix + 1L
      }
      replacement[[i]] <- candidate
      unavailable <- c(unavailable, candidate)
    }
    names(annotation)[match(collisions, names(annotation))] <- replacement
  }

  out <- dplyr::left_join(effects, annotation, by = by)
  if (nrow(out) != nrow(effects)) {
    stop(artifact, " changed the Stage 05 source row count.", call. = FALSE)
  }
  wgcna_stage07_validate_source_preserved(effects, out, artifact)
  out
}

wgcna_stage07_validate_source_preserved <- function(
    source, output, artifact = "Stage 07 interpretable output"
) {
  wgcna_group_effect_consumer_validate_adapted(source)
  wgcna_group_effect_consumer_validate_adapted(output)
  if (nrow(source) != nrow(output)) {
    stop(artifact, " changed the Stage 05 source row count.", call. = FALSE)
  }
  missing <- setdiff(names(source), names(output))
  if (length(missing)) {
    stop(
      artifact, " removed Stage 05 source column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  changed <- names(source)[!vapply(
    names(source),
    function(column) identical(output[[column]], source[[column]]),
    logical(1)
  )]
  if (length(changed)) {
    stop(
      artifact, " changed Stage 05 source value(s) or type(s): ",
      paste(changed, collapse = ", "), ".", call. = FALSE
    )
  }
  invisible(TRUE)
}

.wgcna_stage07_source_key_value <- function(data) {
  key <- wgcna_group_effect_consumer_source_key()
  missing <- setdiff(key, names(data))
  if (length(missing)) {
    stop(
      "Stage 07 inferential source key is missing column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  do.call(paste, c(lapply(data[key], as.character), sep = "\r"))
}

wgcna_stage07_source_statistical_keys <- function(data) {
  key <- wgcna_group_effect_consumer_source_key()
  missing <- setdiff(key, names(data))
  if (length(missing)) {
    stop(
      "Stage 07 statistical source key is missing column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  values <- lapply(key, function(column) {
    paste0(column, "=", as.character(data[[column]]))
  })
  do.call(paste, c(values, sep = "|"))
}

.wgcna_stage07_first_character <- function(data, candidates) {
  out <- rep(NA_character_, nrow(data))
  for (column in intersect(candidates, names(data))) {
    value <- trimws(as.character(data[[column]]))
    value[is.na(value) | !nzchar(value)] <- NA_character_
    replace <- is.na(out) & !is.na(value)
    out[replace] <- value[replace]
  }
  out
}

.wgcna_stage07_first_integer <- function(data, candidates) {
  out <- rep(NA_integer_, nrow(data))
  for (column in intersect(candidates, names(data))) {
    value <- suppressWarnings(as.integer(data[[column]]))
    replace <- is.na(out) & !is.na(value)
    out[replace] <- value[replace]
  }
  out
}

.wgcna_stage07_member_counts <- function(membership) {
  if (!is.data.frame(membership)) {
    stop("Stage 07 inferential membership must be a data frame.", call. = FALSE)
  }
  module_column <- intersect(c("ModuleID", "module_id"), names(membership))
  supermodule_column <- intersect(
    c("SupermoduleID", "supermodule_id"), names(membership)
  )
  if (!length(module_column) || !length(supermodule_column)) {
    stop(
      "Stage 07 inferential membership requires module and supermodule IDs.",
      call. = FALSE
    )
  }
  map <- data.frame(
    module_id = as.character(membership[[module_column[[1]]]]),
    supermodule_id = as.character(membership[[supermodule_column[[1]]]]),
    stringsAsFactors = FALSE
  )
  map <- unique(map[
    !is.na(map$module_id) & nzchar(map$module_id) &
      !is.na(map$supermodule_id) & nzchar(map$supermodule_id),
    ,
    drop = FALSE
  ])
  if (!nrow(map) || anyDuplicated(map$module_id)) {
    stop(
      "Stage 07 inferential membership must map every module once.",
      call. = FALSE
    )
  }
  counts <- aggregate(
    map$module_id,
    by = list(supermodule_id = map$supermodule_id),
    FUN = function(x) length(unique(x))
  )
  names(counts)[[2]] <- "n_member_modules"
  counts$n_member_modules <- as.integer(counts$n_member_modules)
  counts
}

wgcna_stage07_select_inferential_source <- function(
    module_effects, supermodule_effects, membership
) {
  wgcna_group_effect_consumer_validate_adapted(module_effects)
  wgcna_group_effect_consumer_validate_adapted(supermodule_effects)
  counts <- .wgcna_stage07_member_counts(membership)

  module_keep <- module_effects$level == "module" &
    module_effects$independent_hypothesis %in% TRUE &
    module_effects$claim_entity_role == "canonical_module" &
    module_effects$test_type != "conditional_interaction_followup"
  selected_modules <- module_effects[module_keep, , drop = FALSE]
  selected_modules$.inferential_n_member_modules <- 1L

  super_ids <- .wgcna_stage07_first_character(
    supermodule_effects, c("supermodule_id", "endpoint_id")
  )
  super_count <- counts$n_member_modules[
    match(super_ids, counts$supermodule_id)
  ]
  super_keep <- supermodule_effects$level == "supermodule" &
    supermodule_effects$independent_hypothesis %in% TRUE &
    supermodule_effects$claim_entity_role == "higher_order_block" &
    supermodule_effects$test_type != "conditional_interaction_followup" &
    !is.na(super_count) & super_count >= 2L
  selected_supermodules <- supermodule_effects[super_keep, , drop = FALSE]
  selected_supermodules$.inferential_n_member_modules <-
    as.integer(super_count[super_keep])

  out <- dplyr::bind_rows(selected_modules, selected_supermodules)
  if (any(out$claim_entity_role == "compatibility_alias") ||
      any(!out$independent_hypothesis %in% TRUE) ||
      any(out$test_type == "conditional_interaction_followup")) {
    stop(
      "Stage 07 inferential selection retained a non-independent row.",
      call. = FALSE
    )
  }
  expected_module_keys <- .wgcna_stage07_source_key_value(
    module_effects[module_keep, , drop = FALSE]
  )
  observed_module_keys <- .wgcna_stage07_source_key_value(
    out[out$level == "module", , drop = FALSE]
  )
  if (!identical(expected_module_keys, observed_module_keys)) {
    stop(
      "Stage 07 inferential selection removed or reordered a valid module row.",
      call. = FALSE
    )
  }
  out
}

wgcna_stage07_validate_inferential_handoff <- function(
    data, artifact = "Stage 07 inferential handoff"
) {
  required <- c(
    "dataset", "entity_level", "entity_id", "display_label",
    "n_member_modules", "contrast", "analysis_tier", "spatial_scope",
    "estimate", "SE", "CI_low", "CI_high", "p_value",
    "tier_specific_fdr", "tier_specific_family_id",
    "tier_specific_family_size", "model_valid", "model_stability_status",
    "biological_n", "annotation_class", "support_class", "claim_gate",
    "safe_interpretation", "independent_hypothesis", "claim_entity_role",
    "test_type", "effect_scope", "spatial_unit", "result_scope",
    "source_artifact", "source_key", "source_key_contract"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      artifact, " is missing required column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  if (nrow(data)) {
    if (any(!data$independent_hypothesis %in% TRUE) ||
        any(data$claim_entity_role == "compatibility_alias") ||
        any(data$test_type == "conditional_interaction_followup")) {
      stop(artifact, " contains a non-independent endpoint.", call. = FALSE)
    }
    bad_supermodule <- data$entity_level == "supermodule" &
      (is.na(data$n_member_modules) | data$n_member_modules < 2L)
    if (any(bad_supermodule)) {
      stop(
        artifact,
        " contains a singleton-supermodule compatibility alias.",
        call. = FALSE
      )
    }
    key <- c(
      "dataset", "entity_level", "entity_id", "analysis_tier", "contrast",
      "effect_scope", "spatial_unit", "test_type"
    )
    if (anyDuplicated(data[key])) {
      stop(artifact, " contains duplicated inferential endpoints.", call. = FALSE)
    }
  }
  wgcna_group_effect_consumer_validate_adapted(data)
  wgcna_inferential_handoff_validate(data, artifact = artifact)
  invisible(TRUE)
}

wgcna_stage07_build_inferential_handoff <- function(
    module_effects, supermodule_effects,
    module_interpretable, supermodule_interpretable, membership
) {
  selected <- wgcna_stage07_select_inferential_source(
    module_effects, supermodule_effects, membership
  )
  wgcna_stage07_validate_interpretable(module_interpretable)
  wgcna_stage07_validate_interpretable(supermodule_interpretable)

  interpreted <- dplyr::bind_rows(
    module_interpretable, supermodule_interpretable
  )
  interpreted_keys <- .wgcna_stage07_source_key_value(interpreted)
  selected_keys <- .wgcna_stage07_source_key_value(selected)
  if (anyDuplicated(interpreted_keys)) {
    stop(
      "Stage 07 interpretable rows have duplicated source keys.",
      call. = FALSE
    )
  }
  positions <- match(selected_keys, interpreted_keys)
  if (anyNA(positions)) {
    stop(
      "Stage 07 inferential selection could not be matched to annotation rows.",
      call. = FALSE
    )
  }
  annotated <- interpreted[positions, , drop = FALSE]
  n_member_modules <- selected$.inferential_n_member_modules
  selected$.inferential_n_member_modules <- NULL
  for (column in names(selected)) {
    if (!identical(annotated[[column]], selected[[column]])) {
      stop(
        "Stage 07 annotation changed selected statistical column ",
        column, ".", call. = FALSE
      )
    }
  }

  entity_id <- ifelse(
    annotated$level == "module",
    .wgcna_stage07_first_character(
      annotated, c("module_id", "endpoint_id")
    ),
    .wgcna_stage07_first_character(
      annotated, c("supermodule_id", "endpoint_id")
    )
  )
  display_label <- ifelse(
    annotated$level == "module",
    .wgcna_stage07_first_character(
      annotated,
      c(
        "ModulePlotLabel", "interpretable_module_label",
        "module_label_display", "module_biological_label", "endpoint_label",
        "module_id"
      )
    ),
    .wgcna_stage07_first_character(
      annotated,
      c(
        "Supermodule_PlotLabel", "interpretable_supermodule_label",
        "Supermodule_DisplayLabel", "supermodule_label", "endpoint_label",
        "supermodule_id"
      )
    )
  )
  annotation_class <- ifelse(
    annotated$level == "module",
    .wgcna_stage07_first_character(
      annotated,
      c("microenvironment_class", "annotation_class", "annotation_confidence")
    ),
    .wgcna_stage07_first_character(
      annotated,
      c(
        "dominant_microenvironment_class", "annotation_class",
        "annotation_confidence"
      )
    )
  )
  model_valid <- as.logical(annotated$model_valid_for_inference)
  claim_allowed <- as.logical(annotated$claim_allowed_model)
  primary_stable <- as.logical(annotated$primary_model_stable)
  claim_gate <- ifelse(
    model_valid %in% TRUE & claim_allowed %in% TRUE &
      primary_stable %in% TRUE,
    "eligible_for_readiness_assessment",
    ifelse(
      !model_valid %in% TRUE | !primary_stable %in% TRUE,
      "diagnostic_only_model",
      "not_claim_allowed_model"
    )
  )
  canonical <- data.frame(
    dataset = as.character(annotated$dataset),
    entity_level = as.character(annotated$level),
    entity_id = as.character(entity_id),
    display_label = as.character(display_label),
    n_member_modules = as.integer(n_member_modules),
    contrast = as.character(annotated$contrast),
    analysis_tier = as.character(annotated$analysis_tier),
    spatial_scope = as.character(annotated$effect_scope),
    estimate = suppressWarnings(as.numeric(annotated$estimate)),
    SE = suppressWarnings(as.numeric(annotated$SE)),
    CI_low = suppressWarnings(as.numeric(
      if ("CI_low" %in% names(annotated)) annotated$CI_low else NA_real_
    )),
    CI_high = suppressWarnings(as.numeric(
      if ("CI_high" %in% names(annotated)) annotated$CI_high else NA_real_
    )),
    p_value = suppressWarnings(as.numeric(annotated$p_value)),
    tier_specific_fdr = suppressWarnings(as.numeric(
      annotated$tier_specific_fdr
    )),
    tier_specific_family_id = as.character(
      annotated$tier_specific_family_id
    ),
    tier_specific_family_size = as.integer(
      annotated$tier_specific_family_size
    ),
    model_valid = model_valid,
    model_stability_status = as.character(
      annotated$model_stability_status
    ),
    biological_n = .wgcna_stage07_first_integer(
      annotated, c("n_animals_total", "n_animals", "biological_n")
    ),
    annotation_class = as.character(annotation_class),
    support_class = as.character(annotated$statistical_support_status),
    claim_gate = as.character(claim_gate),
    safe_interpretation = as.character(annotated$interpretation_sentence),
    source_artifact = wgcna_inferential_handoff_source_artifact(
      annotated$dataset, annotated$level
    ),
    source_key = wgcna_inferential_handoff_source_key(annotated),
    source_key_contract = rep(
      wgcna_inferential_handoff_source_key_contract(), nrow(annotated)
    ),
    stringsAsFactors = FALSE
  )
  provenance <- annotated[
    , setdiff(names(annotated), names(canonical)),
    drop = FALSE
  ]
  out <- cbind(canonical, provenance, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  wgcna_stage07_validate_inferential_handoff(out)
  out
}

.wgcna_stage07_row_value <- function(row, name, default = NA) {
  if (name %in% names(row) && length(row[[name]])) row[[name]][[1]] else default
}

.wgcna_stage07_nonempty_value <- function(row, candidates) {
  for (name in candidates) {
    value <- as.character(.wgcna_stage07_row_value(
      row, name, NA_character_
    ))
    if (!is.na(value) && nzchar(trimws(value))) return(trimws(value))
  }
  NA_character_
}

wgcna_stage07_support_text <- function(
    status, family_label = "the prespecified family"
) {
  status <- as.character(status)
  out <- rep("Its statistical support status was unavailable.", length(status))
  out[status == "FDR_supported"] <-
    paste0("It was FDR-supported in ", family_label, ".")
  out[status == "suggestive_FDR10"] <-
    paste0(
      "It was suggestive at 10% FDR, but not an FDR-supported result at 5%, ",
      "in ", family_label, "."
    )
  out[status == "not_supported"] <-
    paste0("It was not supported after correction in ", family_label, ".")
  out[status == "inherited_from_canonical_entity"] <-
    "Its support status was inherited for compatibility display only."
  out
}

wgcna_stage07_effect_sentence <- function(row, level = NULL) {
  if (is.null(level)) {
    level <- as.character(.wgcna_stage07_row_value(row, "level", "entity"))
  }
  dataset <- as.character(.wgcna_stage07_row_value(
    row, "dataset", "unknown_dataset"
  ))
  role <- as.character(.wgcna_stage07_row_value(
    row, "claim_entity_role", NA_character_
  ))
  independent <- isTRUE(.wgcna_stage07_row_value(
    row, "independent_hypothesis", FALSE
  ))
  canonical_id <- as.character(.wgcna_stage07_row_value(
    row, "canonical_claim_entity_id", NA_character_
  ))
  entity_id <- .wgcna_stage07_nonempty_value(
    row,
    if (identical(level, "supermodule")) {
      c("supermodule_id", "endpoint_id")
    } else {
      c("module_id", "endpoint_id")
    }
  )

  if (identical(role, "compatibility_alias") || !independent) {
    return(paste0(
      "This singleton supermodule is a compatibility display alias of its ",
      "canonical member module ", canonical_id,
      "; it carries no independent q-value and is not a separate finding."
    ))
  }

  estimate <- suppressWarnings(as.numeric(
    .wgcna_stage07_row_value(row, "estimate", NA_real_)
  ))
  direction <- if (!is.finite(estimate)) {
    "not estimable"
  } else if (estimate > 0) {
    "positive"
  } else if (estimate < 0) {
    "negative"
  } else {
    "zero"
  }
  contrast <- as.character(.wgcna_stage07_row_value(
    row, "contrast", "unspecified contrast"
  ))
  tier <- as.character(.wgcna_stage07_row_value(
    row, "analysis_tier", NA_character_
  ))
  test_type <- as.character(.wgcna_stage07_row_value(
    row, "test_type", NA_character_
  ))
  support <- as.character(.wgcna_stage07_row_value(
    row, "statistical_support_status", NA_character_
  ))
  valid <- isTRUE(.wgcna_stage07_row_value(
    row, "model_valid_for_inference", FALSE
  ))
  claim_allowed <- isTRUE(.wgcna_stage07_row_value(
    row, "claim_allowed_model", FALSE
  ))
  stable <- isTRUE(.wgcna_stage07_row_value(
    row, "primary_model_stable", FALSE
  ))
  stability <- as.character(.wgcna_stage07_row_value(
    row, "model_stability_status", NA_character_
  ))

  model_text <- if (!valid || !claim_allowed) {
    "The source model was not valid or claim-eligible."
  } else if (!stable) {
    paste0(
      "The source model was valid and claim-eligible but not primarily ",
      "stable (", stability, ")."
    )
  } else {
    paste0("The source model was valid and stable (", stability, ").")
  }

  sentence <- if (identical(
    test_type, "conditional_interaction_followup"
  )) {
    paste0(
      "The estimated ", contrast, " conditional interaction difference was ",
      direction,
      ". This is exploratory localization evidence, carries no independent ",
      "q-value, and is not a separate discovery."
    )
  } else if (identical(test_type, "interaction_omnibus")) {
    paste0(
      "The all-group StressGroup by SpatialUnit omnibus interaction was ",
      "assessed. ",
      wgcna_stage07_support_text(
        support, family_label = "the prespecified interaction family"
      ),
      " The responsible contrast was not established by the omnibus test alone."
    )
  } else if (identical(tier, "primary_wgcna_global")) {
    paste0(
      "The estimated SUS-RES global difference for ", level, " ", entity_id,
      " was ", direction, ". ", wgcna_stage07_support_text(support)
    )
  } else if (identical(tier, "secondary_contextual_global")) {
    pattern <- if (identical(support, "FDR_supported") &&
                   identical(contrast, "SUS - CON")) {
      "susceptible deviation pattern"
    } else if (identical(support, "FDR_supported") &&
               identical(contrast, "RES - CON")) {
      "resilient-associated difference"
    } else {
      "inconclusive contextual pattern"
    }
    paste0(
      "The estimated ", contrast, " difference was ", direction,
      " and is contextual evidence only (", pattern, "). ",
      wgcna_stage07_support_text(support),
      " This comparison does not establish equivalence or normalization."
    )
  } else if (identical(tier, "exploratory_spatial_localization")) {
    paste0(
      "The estimated ", contrast, " regional difference was ", direction,
      ". This is exploratory localization evidence and was not an independent ",
      "global discovery. ", wgcna_stage07_support_text(support)
    )
  } else {
    paste0(
      "The estimated ", contrast, " difference was ", direction, ". ",
      wgcna_stage07_support_text(support)
    )
  }

  if (identical(dataset, "microglia")) {
    sentence <- paste0(
      sentence,
      " The biological annotation remains microglia-enriched ROI/local ",
      "microenvironment context and does not establish purified-cell or ",
      "cell-intrinsic regulation."
    )
  }
  paste(sentence, model_text)
}

wgcna_stage07_add_interpretation_sentences <- function(data, level = NULL) {
  if (!nrow(data)) {
    data$interpretation_sentence <- character()
    return(data)
  }
  data$interpretation_sentence <- vapply(
    seq_len(nrow(data)),
    function(i) wgcna_stage07_effect_sentence(
      data[i, , drop = FALSE], level = level
    ),
    character(1)
  )
  data
}

wgcna_stage07_empty_spatial_organization <- function() {
  data.frame(
    dataset = character(),
    level = character(),
    canonical_claim_entity_id = character(),
    entity_id = character(),
    biological_label = character(),
    primary_statistical_support_status = character(),
    primary_estimate = double(),
    primary_confidence_interval = character(),
    primary_tier_specific_fdr = double(),
    interaction_statistical_support_status = character(),
    interaction_tier_specific_fdr = double(),
    interaction_test_type = character(),
    spatial_organization_class = character(),
    spatial_heterogeneity_contrast_attribution = character(),
    localization_evidence_available = logical(),
    localization_support_summary = character(),
    classification_reason = character(),
    primary_source_statistical_key = character(),
    interaction_source_statistical_key = character(),
    stringsAsFactors = FALSE
  )
}

.wgcna_stage07_source_statistical_key <- function(row) {
  key <- wgcna_group_effect_consumer_source_key()
  values <- vapply(
    key,
    function(column) as.character(.wgcna_stage07_row_value(
      row, column, NA_character_
    )),
    character(1)
  )
  paste(paste0(key, "=", values), collapse = "|")
}

.wgcna_stage07_confidence_interval <- function(row) {
  low <- suppressWarnings(as.numeric(
    .wgcna_stage07_row_value(row, "CI_low", NA_real_)
  ))
  high <- suppressWarnings(as.numeric(
    .wgcna_stage07_row_value(row, "CI_high", NA_real_)
  ))
  if (!is.finite(low) || !is.finite(high)) return(NA_character_)
  paste0(
    "[", format(low, digits = 15, scientific = FALSE, trim = TRUE),
    ", ", format(high, digits = 15, scientific = FALSE, trim = TRUE), "]"
  )
}

.wgcna_stage07_localization_summary <- function(local, conditional) {
  count_status <- function(data, status) {
    if (!nrow(data)) return(0L)
    as.integer(sum(
      as.character(data$statistical_support_status) == status,
      na.rm = TRUE
    ))
  }
  paste0(
    "local_FDR_supported=", count_status(local, "FDR_supported"),
    "; local_suggestive_FDR10=", count_status(local, "suggestive_FDR10"),
    "; local_not_supported=",
    count_status(local, "not_supported"),
    "; conditional_followups=", nrow(conditional),
    " (exploratory_no_independent_q)"
  )
}

wgcna_stage07_build_spatial_organization <- function(
    effects, conditional_followups = NULL
) {
  wgcna_group_effect_consumer_validate_adapted(effects)
  if (is.null(conditional_followups)) {
    conditional_followups <- effects[0, , drop = FALSE]
  }
  wgcna_group_effect_consumer_validate_adapted(conditional_followups)

  primary <- wgcna_group_effect_consumer_select_primary(effects)
  if (!nrow(primary)) return(wgcna_stage07_empty_spatial_organization())

  rows <- lapply(seq_len(nrow(primary)), function(i) {
    primary_row <- primary[i, , drop = FALSE]
    same_entity <- as.character(effects$dataset) ==
      as.character(primary_row$dataset[[1]]) &
      as.character(effects$level) ==
      as.character(primary_row$level[[1]]) &
      as.character(effects$canonical_claim_entity_id) ==
      as.character(primary_row$canonical_claim_entity_id[[1]])
    interaction_idx <- which(
      same_entity &
        effects$independent_hypothesis %in% TRUE &
        effects$analysis_tier == "secondary_spatial_heterogeneity" &
        effects$test_type == "interaction_omnibus"
    )
    if (length(interaction_idx) > 1L) {
      stop(
        "Spatial organization found multiple independent interaction omnibus ",
        "rows for one canonical entity.", call. = FALSE
      )
    }
    interaction <- if (length(interaction_idx)) {
      effects[interaction_idx, , drop = FALSE]
    } else {
      effects[0, , drop = FALSE]
    }
    local <- effects[
      which(
        same_entity &
          effects$independent_hypothesis %in% TRUE &
          effects$analysis_tier == "exploratory_spatial_localization" &
          effects$test_type != "conditional_interaction_followup"
      ),
      ,
      drop = FALSE
    ]
    conditional_same <- as.character(conditional_followups$dataset) ==
      as.character(primary_row$dataset[[1]]) &
      as.character(conditional_followups$level) ==
      as.character(primary_row$level[[1]]) &
      as.character(conditional_followups$canonical_claim_entity_id) ==
      as.character(primary_row$canonical_claim_entity_id[[1]])
    conditional <- conditional_followups[
      which(
        conditional_same &
          conditional_followups$test_type ==
            "conditional_interaction_followup" &
          conditional_followups$claim_entity_role != "compatibility_alias"
      ),
      ,
      drop = FALSE
    ]

    primary_allowed <- isTRUE(primary_row$model_valid_for_inference[[1]]) &&
      isTRUE(primary_row$claim_allowed_model[[1]])
    interaction_allowed <- nrow(interaction) == 1L &&
      isTRUE(interaction$model_valid_for_inference[[1]]) &&
      isTRUE(interaction$claim_allowed_model[[1]])
    primary_supported <- identical(
      as.character(primary_row$statistical_support_status[[1]]),
      "FDR_supported"
    )
    interaction_supported <- nrow(interaction) == 1L && identical(
      as.character(interaction$statistical_support_status[[1]]),
      "FDR_supported"
    )
    spatial_class <- if (!primary_allowed || !interaction_allowed) {
      "model_or_evidence_not_claim_allowed"
    } else if (primary_supported && interaction_supported) {
      "primary_global_shift_with_group_dependent_spatial_heterogeneity"
    } else if (primary_supported) {
      "primary_global_shift"
    } else if (interaction_supported) {
      "group_dependent_spatial_heterogeneity_without_primary_global_shift"
    } else {
      "no_FDR_supported_group_organization"
    }

    supported_local <- local[
      local$statistical_support_status == "FDR_supported",
      ,
      drop = FALSE
    ]
    local_contrasts <- sort(unique(as.character(
      supported_local$contrast[
        !is.na(supported_local$contrast) &
          nzchar(as.character(supported_local$contrast))
      ]
    )))
    attribution <- if (!interaction_supported) {
      "not_attributed_to_specific_contrast"
    } else if (length(local_contrasts) == 1L) {
      paste0(
        "attributed_by_FDR_supported_localization:",
        gsub("\\s+", "", local_contrasts[[1]])
      )
    } else if (length(local_contrasts) > 1L) {
      "multiple_contrasts_supported_by_localization"
    } else {
      "not_attributed_to_specific_contrast"
    }

    localization_available <- nrow(local) > 0L || nrow(conditional) > 0L
    reason <- paste0(
      "Class uses only the prespecified primary-family status and the ",
      "independent all-group interaction-omnibus status. Conditional and ",
      "within-spatial-unit rows do not determine the class."
    )
    if (localization_available) {
      reason <- paste0(
        reason,
        " Opposite localization directions are biological context and are ",
        "not classified as model instability."
      )
    }
    if (!isTRUE(primary_row$primary_model_stable[[1]]) ||
        (nrow(interaction) &&
         !isTRUE(interaction$primary_model_stable[[1]]))) {
      reason <- paste0(
        reason,
        " At least one valid source model was not primarily stable; this ",
        "remains explicit for later readiness review."
      )
    }

    level <- as.character(primary_row$level[[1]])
    data.frame(
      dataset = as.character(primary_row$dataset[[1]]),
      level = level,
      canonical_claim_entity_id = as.character(
        primary_row$canonical_claim_entity_id[[1]]
      ),
      entity_id = .wgcna_stage07_nonempty_value(
        primary_row,
        if (identical(level, "module")) {
          c("module_id", "endpoint_id")
        } else {
          c("supermodule_id", "endpoint_id")
        }
      ),
      biological_label = .wgcna_stage07_nonempty_value(
        primary_row,
        if (identical(level, "module")) {
          c(
            "ModulePlotLabel", "module_biological_label",
            "module_label", "endpoint_label"
          )
        } else {
          c(
            "Supermodule_PlotLabel", "Supermodule_DisplayLabel",
            "supermodule_label", "endpoint_label"
          )
        }
      ),
      primary_statistical_support_status =
        as.character(primary_row$statistical_support_status[[1]]),
      primary_estimate = as.numeric(primary_row$estimate[[1]]),
      primary_confidence_interval =
        .wgcna_stage07_confidence_interval(primary_row),
      primary_tier_specific_fdr =
        as.numeric(primary_row$tier_specific_fdr[[1]]),
      interaction_statistical_support_status = if (nrow(interaction)) {
        as.character(interaction$statistical_support_status[[1]])
      } else {
        NA_character_
      },
      interaction_tier_specific_fdr = if (nrow(interaction)) {
        as.numeric(interaction$tier_specific_fdr[[1]])
      } else {
        NA_real_
      },
      interaction_test_type = if (nrow(interaction)) {
        as.character(interaction$test_type[[1]])
      } else {
        NA_character_
      },
      spatial_organization_class = spatial_class,
      spatial_heterogeneity_contrast_attribution = attribution,
      localization_evidence_available = localization_available,
      localization_support_summary =
        .wgcna_stage07_localization_summary(local, conditional),
      classification_reason = reason,
      primary_source_statistical_key =
        .wgcna_stage07_source_statistical_key(primary_row),
      interaction_source_statistical_key = if (nrow(interaction)) {
        .wgcna_stage07_source_statistical_key(interaction)
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(rows)
  out <- out[order(
    out$dataset, out$level, out$canonical_claim_entity_id
  ), , drop = FALSE]
  row.names(out) <- NULL
  out
}

wgcna_stage07_validate_spatial_organization <- function(
    data, artifact = "Stage 07 spatial organization summary"
) {
  required <- names(wgcna_stage07_empty_spatial_organization())
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      artifact, " is missing required column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  if (!nrow(data)) return(invisible(TRUE))
  key <- c("dataset", "level", "canonical_claim_entity_id")
  if (anyDuplicated(data[key])) {
    stop(
      artifact, " contains duplicate independent canonical entity rows.",
      call. = FALSE
    )
  }
  allowed <- c(
    "primary_global_shift",
    "primary_global_shift_with_group_dependent_spatial_heterogeneity",
    "group_dependent_spatial_heterogeneity_without_primary_global_shift",
    "no_FDR_supported_group_organization",
    "model_or_evidence_not_claim_allowed"
  )
  invalid <- setdiff(
    unique(as.character(data$spatial_organization_class)), allowed
  )
  if (length(invalid)) {
    stop(
      artifact, " contains invalid spatial_organization_class value(s): ",
      paste(invalid, collapse = ", "), ".", call. = FALSE
    )
  }
  if (any(
    data$interaction_test_type != "interaction_omnibus", na.rm = TRUE
  )) {
    stop(
      artifact, " must use only independent interaction omnibus rows.",
      call. = FALSE
    )
  }
  attribution <- as.character(
    data$spatial_heterogeneity_contrast_attribution
  )
  allowed_attribution <- attribution %in% c(
    "not_attributed_to_specific_contrast",
    "multiple_contrasts_supported_by_localization"
  ) | grepl(
    "^attributed_by_FDR_supported_localization:.+$", attribution
  )
  if (any(!allowed_attribution | is.na(allowed_attribution))) {
    invalid_attribution <- unique(attribution[
      !allowed_attribution | is.na(allowed_attribution)
    ])
    stop(
      artifact, " contains invalid spatial-heterogeneity contrast ",
      "attribution value(s): ",
      paste(invalid_attribution, collapse = ", "), ".", call. = FALSE
    )
  }
  contrast_attributed <- attribution !=
    "not_attributed_to_specific_contrast"
  omnibus_supported <- as.character(
    data$interaction_statistical_support_status
  ) == "FDR_supported"
  if (any(
    contrast_attributed & (!omnibus_supported | is.na(omnibus_supported))
  )) {
    stop(
      artifact, " cannot attribute spatial heterogeneity to a contrast ",
      "unless the independent interaction omnibus is FDR_supported.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------
# Presentation-only helpers.
#
# These shape what a figure renders and which stale files a zero-row
# result may leave behind. They never touch the source tables, the
# statistics, or any support classification.
# ---------------------------------------------------------------

# Keep only rows a contrast dot plot can actually draw: an estimable
# named contrast whose label is one of the plotted levels. Interaction
# omnibus rows are non-estimable by construction (estimate = NA, and a
# contrast label outside the plotted levels), so exclude them here
# rather than letting ggplot2 drop them silently and warn. This filters
# the plotting input only; the written source table is unaffected.
wgcna_stage07_filter_named_contrast_estimable <- function(data, contrast_levels) {
  if (!is.data.frame(data)) {
    stop(
      "Named-contrast plot filtering requires a data frame.", call. = FALSE
    )
  }
  if (!length(contrast_levels)) {
    stop(
      "Named-contrast plot filtering requires at least one contrast level.",
      call. = FALSE
    )
  }
  if (!nrow(data)) return(data)

  missing <- setdiff(c("test_type", "estimate", "contrast"), names(data))
  if (length(missing)) {
    stop(
      "Named-contrast plot filtering is missing required column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }

  keep <- as.character(data$test_type) == "named_contrast" &
    is.finite(suppressWarnings(as.numeric(data$estimate))) &
    as.character(data$contrast) %in% as.character(contrast_levels)

  data[which(keep), , drop = FALSE]
}

# A zero-row current result must not leave a historical figure on disk
# that a reader would mistake for the current canonical output. Remove
# the SVG and its PDF sibling deterministically.
wgcna_stage07_remove_stale_figure <- function(svg_path) {
  targets <- c(svg_path, sub("\\.svg$", ".pdf", svg_path))
  removed <- character(0)
  for (target in targets) {
    if (!file.exists(target)) next
    unlink(target)
    if (file.exists(target)) {
      stop(
        "Could not remove stale figure: ", target,
        ". A zero-row result must not leave a historical figure in place.",
        call. = FALSE
      )
    }
    removed <- c(removed, target)
  }
  if (length(removed)) {
    message(
      "Removed stale figure(s) for a zero-row current result: ",
      paste(basename(removed), collapse = ", ")
    )
  }
  invisible(removed)
}
