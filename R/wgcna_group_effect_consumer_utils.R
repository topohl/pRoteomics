# Stage 05 v5 WGCNA group-effect consumer contract.
#
# This helper is a semantic handoff only. It must not calculate FDR values,
# rebuild statistical support, decide manuscript readiness, or label biology.

if (!exists("wgcna_group_classify_statistical_support", mode = "function")) {
  support_utils <- c(
    if (exists("repo_path", mode = "function")) {
      repo_path("R", "wgcna_support_status_utils.R")
    } else character(),
    file.path("R", "wgcna_support_status_utils.R"),
    file.path("..", "R", "wgcna_support_status_utils.R"),
    file.path("..", "..", "R", "wgcna_support_status_utils.R")
  )
  support_utils <- support_utils[file.exists(support_utils)][[1]]
  source(support_utils)
  rm(support_utils)
}

wgcna_group_effect_consumer_added_columns <- function() {
  c(
    "tier_specific_fdr",
    "tier_specific_family_id",
    "tier_specific_family_size",
    "result_scope"
  )
}

wgcna_group_effect_consumer_source_key <- function() {
  c(
    "dataset", "level", "endpoint_id", "effect_scope",
    "spatial_unit", "contrast", "test_type"
  )
}

wgcna_group_effect_consumer_semantic_key <- function() {
  c(
    "dataset", "level", "canonical_claim_entity_id", "analysis_tier",
    "contrast", "effect_scope", "spatial_unit", "test_type"
  )
}

.wgcna_group_effect_consumer_family_contract <- function() {
  list(
    primary_wgcna_global = list(
      fdr = "FDR_primary_global",
      family_id = "FDR_primary_family_id",
      family_size = "n_tests_FDR_primary",
      result_scope = "primary_global_vulnerability"
    ),
    secondary_contextual_global = list(
      fdr = "FDR_secondary_global",
      family_id = "FDR_secondary_family_id",
      family_size = "n_tests_FDR_secondary_global",
      result_scope = "contextual_global_control"
    ),
    secondary_spatial_heterogeneity = list(
      fdr = "FDR_interaction_omnibus",
      family_id = "FDR_interaction_family_id",
      family_size = "n_tests_FDR_interaction_omnibus",
      result_scope = "spatial_heterogeneity_omnibus"
    ),
    exploratory_spatial_localization = list(
      fdr = "FDR_local_exploratory",
      family_id = "FDR_local_family_id",
      family_size = "n_tests_FDR_local_exploratory",
      result_scope = "exploratory_localization"
    )
  )
}

.wgcna_group_effect_consumer_family_columns <- function() {
  contract <- .wgcna_group_effect_consumer_family_contract()
  unique(unlist(lapply(
    contract,
    function(x) unname(unlist(x[c("fdr", "family_id", "family_size")]))
  ), use.names = FALSE))
}

.wgcna_group_effect_consumer_required_columns <- function() {
  unique(c(
    wgcna_group_effect_consumer_source_key(),
    wgcna_group_effect_consumer_semantic_key(),
    "claim_entity_role", "support_source_entity_id", "p_value",
    "independent_hypothesis", "statistical_support_status",
    "model_valid_for_inference", "model_stability_status",
    "claim_allowed_model", "primary_model_stable",
    .wgcna_group_effect_consumer_family_columns()
  ))
}

.wgcna_group_effect_consumer_independent_q_columns <- function(data) {
  intersect(
    c(
      "FDR_primary_global", "FDR_secondary_global",
      "FDR_interaction_omnibus", "FDR_local_exploratory",
      "FDR_conservative_all_tests", "FDR_global",
      "FDR_within_dataset_level", "FDR_dataset_all_levels"
    ),
    names(data)
  )
}

.wgcna_group_effect_consumer_independent_family_metadata_columns <- function(
    data
) {
  intersect(
    c(
      "FDR_conservative_family_id",
      "n_tests_FDR_conservative_all_tests",
      "FDR_family_within_level_id",
      "n_tests_FDR_within_dataset_level",
      "FDR_family_dataset_id",
      "n_tests_FDR_dataset_all_levels"
    ),
    names(data)
  )
}

.wgcna_group_effect_consumer_nonempty_character <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

.wgcna_group_effect_consumer_character_or_all_na <- function(x) {
  is.character(x) || all(is.na(x))
}

.wgcna_group_effect_consumer_numeric_or_all_na <- function(x) {
  is.numeric(x) || all(is.na(x))
}

.wgcna_group_effect_consumer_rows_text <- function(rows) {
  paste(utils::head(rows, 10L), collapse = ", ")
}

.wgcna_group_effect_consumer_stop_rows <- function(message, rows) {
  stop(
    message, " Offending row(s): ",
    .wgcna_group_effect_consumer_rows_text(rows), ".",
    call. = FALSE
  )
}

.wgcna_group_effect_consumer_require_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      "Stage 05 v5 group-effect input is missing required column(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.wgcna_group_effect_consumer_validate_types <- function(data) {
  character_columns <- unique(c(
    wgcna_group_effect_consumer_source_key(),
    wgcna_group_effect_consumer_semantic_key(),
    "claim_entity_role", "support_source_entity_id",
    "statistical_support_status",
    vapply(
      .wgcna_group_effect_consumer_family_contract(),
      `[[`, character(1), "family_id"
    )
  ))
  logical_columns <- c(
    "independent_hypothesis", "model_valid_for_inference",
    "claim_allowed_model", "primary_model_stable"
  )
  numeric_columns <- vapply(
    .wgcna_group_effect_consumer_family_contract(),
    `[[`, character(1), "fdr"
  )
  size_columns <- vapply(
    .wgcna_group_effect_consumer_family_contract(),
    `[[`, character(1), "family_size"
  )

  bad_character <- character_columns[
    !vapply(
      data[character_columns],
      .wgcna_group_effect_consumer_character_or_all_na,
      logical(1)
    )
  ]
  if (length(bad_character)) {
    stop(
      "Stage 05 v5 character column(s) have invalid types: ",
      paste(bad_character, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  bad_logical <- logical_columns[
    !vapply(data[logical_columns], is.logical, logical(1))
  ]
  if (length(bad_logical)) {
    stop(
      "Stage 05 v5 logical column(s) have invalid types: ",
      paste(bad_logical, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  bad_numeric <- numeric_columns[
    !vapply(
      data[numeric_columns],
      .wgcna_group_effect_consumer_numeric_or_all_na,
      logical(1)
    )
  ]
  if (length(bad_numeric)) {
    stop(
      "Stage 05 v5 FDR column(s) have invalid types: ",
      paste(bad_numeric, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  bad_size <- size_columns[
    !vapply(
      data[size_columns],
      .wgcna_group_effect_consumer_numeric_or_all_na,
      logical(1)
    )
  ]
  if (length(bad_size)) {
    stop(
      "Stage 05 v5 family-size column(s) have invalid types: ",
      paste(bad_size, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.character(data$model_stability_status)) {
    stop(
      "Stage 05 v5 model_stability_status must be character.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.wgcna_group_effect_consumer_any_populated <- function(data, rows, columns) {
  if (!length(rows) || !length(columns)) {
    return(logical(length(rows)))
  }
  Reduce(
    `|`,
    lapply(columns, function(column) {
      values <- data[[column]][rows]
      if (is.character(values)) {
        .wgcna_group_effect_consumer_nonempty_character(values)
      } else {
        !is.na(values)
      }
    }),
    init = rep(FALSE, length(rows))
  )
}

.wgcna_group_effect_consumer_validate_special_rows <- function(data) {
  alias <- data$claim_entity_role == "compatibility_alias"
  conditional <- data$test_type == "conditional_interaction_followup"

  inconsistent_alias <- which(alias & data$independent_hypothesis %in% TRUE)
  if (length(inconsistent_alias)) {
    .wgcna_group_effect_consumer_stop_rows(
      "Compatibility aliases must have independent_hypothesis = FALSE.",
      inconsistent_alias
    )
  }

  unidentified_nonindependent <- which(
    !data$independent_hypothesis %in% TRUE & !alias
  )
  if (length(unidentified_nonindependent)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(
        "Non-independent rows must be explicitly identified as ",
        "claim_entity_role = compatibility_alias."
      ),
      unidentified_nonindependent
    )
  }

  alias_conditional <- which(alias & conditional)
  if (length(alias_conditional)) {
    .wgcna_group_effect_consumer_stop_rows(
      "A row cannot be both a compatibility alias and a conditional follow-up.",
      alias_conditional
    )
  }

  independent_q <- .wgcna_group_effect_consumer_independent_q_columns(data)
  family_metadata <- c(
    .wgcna_group_effect_consumer_family_columns(),
    .wgcna_group_effect_consumer_independent_family_metadata_columns(data)
  )
  prohibited <- unique(c(independent_q, family_metadata))

  alias_rows <- which(alias)
  alias_has_q <- .wgcna_group_effect_consumer_any_populated(
    data, alias_rows, prohibited
  )
  if (any(alias_has_q)) {
    .wgcna_group_effect_consumer_stop_rows(
      "Compatibility aliases cannot carry an independent q-value or family.",
      alias_rows[alias_has_q]
    )
  }

  conditional_rows <- which(conditional)
  conditional_has_q <- .wgcna_group_effect_consumer_any_populated(
    data, conditional_rows, prohibited
  )
  if (any(conditional_has_q)) {
    .wgcna_group_effect_consumer_stop_rows(
      "Conditional interaction follow-ups cannot carry an independent q-value or family.",
      conditional_rows[conditional_has_q]
    )
  }
  invisible(TRUE)
}

.wgcna_group_effect_consumer_validate_family_rows <- function(data) {
  contract <- .wgcna_group_effect_consumer_family_contract()
  alias <- data$claim_entity_role == "compatibility_alias"
  conditional <- data$test_type == "conditional_interaction_followup"
  generic <- !alias & !conditional & data$independent_hypothesis %in% TRUE

  unsupported <- which(generic & !data$analysis_tier %in% names(contract))
  if (length(unsupported)) {
    .wgcna_group_effect_consumer_stop_rows(
      "Independent rows have unsupported analysis_tier values.",
      unsupported
    )
  }

  for (tier in names(contract)) {
    spec <- contract[[tier]]
    rows <- which(generic & data$analysis_tier == tier)
    if (!length(rows)) {
      next
    }

    fdr <- data[[spec$fdr]][rows]
    bad_fdr <- is.na(fdr) | !is.finite(fdr) | fdr < 0 | fdr > 1
    if (any(bad_fdr)) {
      .wgcna_group_effect_consumer_stop_rows(
        paste0(
          "Independent ", tier, " rows require a finite ", spec$fdr,
          " value in [0, 1]."
        ),
        rows[bad_fdr]
      )
    }

    family_id <- data[[spec$family_id]][rows]
    bad_family_id <- !.wgcna_group_effect_consumer_nonempty_character(
      family_id
    )
    if (any(bad_family_id)) {
      .wgcna_group_effect_consumer_stop_rows(
        paste0(
          "Independent ", tier, " rows require ", spec$family_id, "."
        ),
        rows[bad_family_id]
      )
    }

    family_size <- data[[spec$family_size]][rows]
    bad_family_size <- is.na(family_size) |
      !is.finite(family_size) |
      family_size < 1 |
      family_size != floor(family_size) |
      family_size > .Machine$integer.max
    if (any(bad_family_size)) {
      .wgcna_group_effect_consumer_stop_rows(
        paste0(
          "Independent ", tier, " rows require a positive integer-valued ",
          spec$family_size, "."
        ),
        rows[bad_family_size]
      )
    }

    other_tiers <- setdiff(names(contract), tier)
    other_columns <- unique(unlist(lapply(
      contract[other_tiers],
      function(x) unname(unlist(x[c("fdr", "family_id", "family_size")]))
    ), use.names = FALSE))
    cross_family <- .wgcna_group_effect_consumer_any_populated(
      data, rows, other_columns
    )
    if (any(cross_family)) {
      .wgcna_group_effect_consumer_stop_rows(
        paste0(
          "Independent ", tier,
          " rows cannot populate another tier-specific FDR family."
        ),
        rows[cross_family]
      )
    }
  }
  invisible(TRUE)
}

.wgcna_group_effect_consumer_validate_support_status <- function(data) {
  rows <- which(
    data$independent_hypothesis %in% TRUE &
      data$claim_entity_role != "compatibility_alias" &
      data$test_type != "conditional_interaction_followup"
  )
  if (!length(rows)) return(invisible(TRUE))

  applicable_fdr <- wgcna_group_applicable_fdr(data[rows, , drop = FALSE])
  expected <- wgcna_group_classify_statistical_support(
    data$p_value[rows], applicable_fdr
  )
  observed <- as.character(data$statistical_support_status[rows])
  mismatch <- which(is.na(observed) | observed != expected)
  if (length(mismatch)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(
        "statistical_support_status does not match the applicable ",
        "multiplicity-adjusted result; nominal p-values cannot restore ",
        "support after correction failure."
      ),
      rows[mismatch]
    )
  }
  invisible(TRUE)
}

wgcna_group_effect_consumer_validate <- function(data) {
  if (!is.data.frame(data)) {
    stop("Stage 05 v5 group-effect input must be a data frame.", call. = FALSE)
  }
  .wgcna_group_effect_consumer_require_columns(
    data, .wgcna_group_effect_consumer_required_columns()
  )
  .wgcna_group_effect_consumer_validate_types(data)

  source_key <- wgcna_group_effect_consumer_source_key()
  for (column in source_key) {
    bad <- which(
      !.wgcna_group_effect_consumer_nonempty_character(data[[column]])
    )
    if (length(bad)) {
      .wgcna_group_effect_consumer_stop_rows(
        paste0("Stage 05 source-key column ", column, " is missing."),
        bad
      )
    }
  }
  duplicated_source <- which(duplicated(data[source_key]))
  if (length(duplicated_source)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(
        "Stage 05 source keys are duplicated: ",
        paste(source_key, collapse = " + "), "."
      ),
      duplicated_source
    )
  }

  if (anyNA(data$independent_hypothesis)) {
    .wgcna_group_effect_consumer_stop_rows(
      "independent_hypothesis cannot be NA.",
      which(is.na(data$independent_hypothesis))
    )
  }
  if (any(!.wgcna_group_effect_consumer_nonempty_character(
    data$claim_entity_role
  ))) {
    .wgcna_group_effect_consumer_stop_rows(
      "claim_entity_role cannot be missing.",
      which(!.wgcna_group_effect_consumer_nonempty_character(
        data$claim_entity_role
      ))
    )
  }
  missing_support <- !.wgcna_group_effect_consumer_nonempty_character(
    data$statistical_support_status
  )
  missing_support[data$test_type == "conditional_interaction_followup"] <-
    FALSE
  if (any(missing_support)) {
    .wgcna_group_effect_consumer_stop_rows(
      "statistical_support_status cannot be missing.",
      which(missing_support)
    )
  }

  semantic_rows <- which(
    data$independent_hypothesis %in% TRUE &
      data$claim_entity_role != "compatibility_alias"
  )
  semantic_key <- wgcna_group_effect_consumer_semantic_key()
  for (column in semantic_key) {
    bad <- semantic_rows[
      !.wgcna_group_effect_consumer_nonempty_character(
        data[[column]][semantic_rows]
      )
    ]
    if (length(bad)) {
      .wgcna_group_effect_consumer_stop_rows(
        paste0(
          "Downstream semantic-key column ", column,
          " is missing for an independent row."
        ),
        bad
      )
    }
  }
  duplicated_semantic <- which(duplicated(data[semantic_rows, semantic_key]))
  if (length(duplicated_semantic)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(
        "Independent downstream semantic keys are duplicated: ",
        paste(semantic_key, collapse = " + "), "."
      ),
      semantic_rows[duplicated_semantic]
    )
  }

  .wgcna_group_effect_consumer_validate_special_rows(data)
  .wgcna_group_effect_consumer_validate_family_rows(data)
  .wgcna_group_effect_consumer_validate_support_status(data)
  invisible(TRUE)
}

.wgcna_group_effect_consumer_expected_adapter_columns <- function(data) {
  n <- nrow(data)
  expected <- data.frame(
    tier_specific_fdr = rep(NA_real_, n),
    tier_specific_family_id = rep(NA_character_, n),
    tier_specific_family_size = rep(NA_integer_, n),
    result_scope = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )

  alias <- data$claim_entity_role == "compatibility_alias"
  conditional <- data$test_type == "conditional_interaction_followup"
  expected$result_scope[alias] <- "compatibility_alias_display"
  expected$result_scope[conditional & !alias] <-
    "exploratory_conditional_followup"

  contract <- .wgcna_group_effect_consumer_family_contract()
  generic <- !alias & !conditional & data$independent_hypothesis %in% TRUE
  for (tier in names(contract)) {
    spec <- contract[[tier]]
    rows <- which(generic & data$analysis_tier == tier)
    if (!length(rows)) {
      next
    }
    expected$tier_specific_fdr[rows] <- data[[spec$fdr]][rows]
    expected$tier_specific_family_id[rows] <- data[[spec$family_id]][rows]
    expected$tier_specific_family_size[rows] <-
      as.integer(data[[spec$family_size]][rows])
    expected$result_scope[rows] <- spec$result_scope
  }

  if (anyNA(expected$result_scope)) {
    .wgcna_group_effect_consumer_stop_rows(
      "The adapter could not assign a controlled result_scope.",
      which(is.na(expected$result_scope))
    )
  }
  expected
}

.wgcna_group_effect_consumer_exact_matches <- function(observed, expected) {
  if (length(observed) != length(expected)) {
    return(rep(FALSE, max(length(observed), length(expected))))
  }
  vapply(
    seq_along(expected),
    function(i) identical(observed[i], expected[i]),
    logical(1)
  )
}

wgcna_group_effect_consumer_adapt <- function(data) {
  collisions <- intersect(
    wgcna_group_effect_consumer_added_columns(), names(data)
  )
  if (length(collisions)) {
    stop(
      "Stage 05 v5 group-effect input already contains adapter column(s): ",
      paste(collisions, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  wgcna_group_effect_consumer_validate(data)

  out <- data
  expected <- .wgcna_group_effect_consumer_expected_adapter_columns(out)
  for (column in wgcna_group_effect_consumer_added_columns()) {
    out[[column]] <- expected[[column]]
  }
  out
}

wgcna_group_effect_consumer_validate_adapted <- function(data) {
  wgcna_group_effect_consumer_validate(data)
  added <- wgcna_group_effect_consumer_added_columns()
  .wgcna_group_effect_consumer_require_columns(data, added)

  expected <- .wgcna_group_effect_consumer_expected_adapter_columns(data)
  for (column in added) {
    if (!identical(typeof(data[[column]]), typeof(expected[[column]]))) {
      stop(
        "Adapted Stage 05 column ", column,
        " has type ", typeof(data[[column]]),
        "; expected ", typeof(expected[[column]]), ".",
        call. = FALSE
      )
    }
    matches <- .wgcna_group_effect_consumer_exact_matches(
      data[[column]], expected[[column]]
    )
    if (any(!matches)) {
      .wgcna_group_effect_consumer_stop_rows(
        paste0(
          "Adapted Stage 05 column ", column,
          " does not exactly match its deterministic source-field mapping."
        ),
        which(!matches)
      )
    }
  }
  invisible(TRUE)
}

wgcna_group_effect_consumer_validate_primary_selection <- function(data) {
  wgcna_group_effect_consumer_validate_adapted(data)
  if (!nrow(data)) {
    return(invisible(TRUE))
  }

  valid_scope <- data$independent_hypothesis %in% TRUE &
    data$analysis_tier == "primary_wgcna_global" &
    data$contrast == "SUS - RES" &
    data$effect_scope == "spatial_adjusted_global" &
    data$spatial_unit == "global_spatial_adjusted" &
    data$result_scope == "primary_global_vulnerability"
  if (any(!valid_scope)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(
        "Primary claim selections may contain only independent SUS - RES ",
        "primary spatial-adjusted global rows."
      ),
      which(!valid_scope)
    )
  }

  claim_key <- c("dataset", "level", "canonical_claim_entity_id")
  duplicated_claim <- which(duplicated(data[claim_key]))
  if (length(duplicated_claim)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(
        "Primary claim selections contain duplicate scoped canonical keys: ",
        paste(claim_key, collapse = " + "), "."
      ),
      duplicated_claim
    )
  }
  invisible(TRUE)
}

wgcna_group_effect_consumer_select_primary <- function(
    data, support_status = NULL
) {
  wgcna_group_effect_consumer_validate_adapted(data)
  keep <- data$independent_hypothesis %in% TRUE &
    data$analysis_tier == "primary_wgcna_global" &
    data$contrast == "SUS - RES" &
    data$effect_scope == "spatial_adjusted_global" &
    data$spatial_unit == "global_spatial_adjusted"

  if (!is.null(support_status)) {
    if (!is.character(support_status) || anyNA(support_status)) {
      stop("support_status must be NULL or a non-missing character vector.",
           call. = FALSE)
    }
    keep <- keep & data$statistical_support_status %in% support_status
  }

  out <- data[keep, , drop = FALSE]
  wgcna_group_effect_consumer_validate_primary_selection(out)
  out
}

wgcna_inferential_handoff_source_key_contract <- function() {
  "dataset+level+endpoint_id+effect_scope+spatial_unit+contrast+test_type"
}

wgcna_inferential_handoff_source_artifact <- function(dataset, entity_level) {
  dataset <- as.character(dataset)
  entity_level <- as.character(entity_level)
  filename <- ifelse(
    entity_level == "module",
    "module_group_effects.csv",
    ifelse(
      entity_level == "supermodule",
      "supermodule_group_effects.csv",
      NA_character_
    )
  )
  file.path(
    "results", "tables", "06_modules_WGCNA", "group_effects",
    dataset, filename
  )
}

wgcna_inferential_handoff_source_key <- function(data) {
  key <- wgcna_group_effect_consumer_source_key()
  .wgcna_group_effect_consumer_require_columns(data, key)
  values <- lapply(key, function(column) {
    paste0(column, "=", as.character(data[[column]]))
  })
  do.call(paste, c(values, sep = "|"))
}

wgcna_inferential_handoff_required_columns <- function() {
  c(
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
}

wgcna_inferential_handoff_normalize_csv_types <- function(
    data, artifact = "WGCNA_inferential_handoff.csv"
) {
  if (!is.data.frame(data)) {
    stop(artifact, " must be a data frame.", call. = FALSE)
  }
  column <- "tier_specific_family_size"
  if (!column %in% names(data)) return(data)

  values <- data[[column]]
  if (!is.numeric(values)) {
    stop(
      artifact, " column ", column,
      " must be numeric before CSV type normalization.",
      call. = FALSE
    )
  }
  populated <- !is.na(values)
  invalid <- populated & (
    !is.finite(values) |
      values < 0 |
      values != trunc(values) |
      values > .Machine$integer.max
  )
  if (any(invalid)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(
        artifact, " column ", column,
        " must contain only finite, non-negative, whole values within ",
        "the R integer range"
      ),
      which(invalid)
    )
  }
  data[[column]] <- as.integer(values)
  data
}

wgcna_inferential_handoff_read <- function(
    path, artifact = as.character(path)
) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required to read the inferential handoff.",
         call. = FALSE)
  }
  data <- readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "unique"
  )
  data <- wgcna_inferential_handoff_normalize_csv_types(data, artifact)
  wgcna_inferential_handoff_validate(data, artifact)
  data
}

wgcna_inferential_handoff_validate <- function(
    data, artifact = "WGCNA_inferential_handoff.csv"
) {
  if (!is.data.frame(data)) {
    stop(artifact, " must be a data frame.", call. = FALSE)
  }
  missing <- setdiff(wgcna_inferential_handoff_required_columns(), names(data))
  if (length(missing)) {
    stop(
      artifact, " is missing required column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  wgcna_group_effect_consumer_validate_adapted(data)
  if (!nrow(data)) return(invisible(TRUE))

  invalid_independent <- !data$independent_hypothesis %in% TRUE |
    data$claim_entity_role == "compatibility_alias" |
    data$test_type == "conditional_interaction_followup"
  if (any(invalid_independent)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " contains a non-independent inferential endpoint."),
      which(invalid_independent)
    )
  }
  invalid_level <- !data$entity_level %in% c("module", "supermodule") |
    data$entity_level != data$level
  if (any(invalid_level)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " has inconsistent entity_level and Stage 05 level."),
      which(invalid_level)
    )
  }
  invalid_entity <- is.na(data$entity_id) | !nzchar(data$entity_id) |
    data$entity_id != data$endpoint_id
  if (any(invalid_entity)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " has an entity_id inconsistent with endpoint_id."),
      which(invalid_entity)
    )
  }
  invalid_count <- is.na(data$n_member_modules) |
    data$n_member_modules < 1L |
    (data$entity_level == "module" & data$n_member_modules != 1L) |
    (data$entity_level == "supermodule" & data$n_member_modules < 2L)
  if (any(invalid_count)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " has invalid module membership cardinality."),
      which(invalid_count)
    )
  }

  expected_artifact <- wgcna_inferential_handoff_source_artifact(
    data$dataset, data$entity_level
  )
  bad_artifact <- is.na(data$source_artifact) |
    data$source_artifact != expected_artifact
  if (any(bad_artifact)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " has incorrect Stage 05 source_artifact trace."),
      which(bad_artifact)
    )
  }
  expected_key <- wgcna_inferential_handoff_source_key(data)
  bad_key <- is.na(data$source_key) | data$source_key != expected_key
  if (any(bad_key)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " has incorrect Stage 05 source_key trace."),
      which(bad_key)
    )
  }
  expected_contract <- wgcna_inferential_handoff_source_key_contract()
  bad_contract <- is.na(data$source_key_contract) |
    data$source_key_contract != expected_contract
  if (any(bad_contract)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " has an unsupported source_key_contract."),
      which(bad_contract)
    )
  }

  source_model_valid <- as.logical(data$model_valid_for_inference)
  handoff_model_valid <- as.logical(data$model_valid)
  bad_model <- is.na(source_model_valid) | is.na(handoff_model_valid) |
    source_model_valid != handoff_model_valid
  if (any(bad_model)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " changed model validity."),
      which(bad_model)
    )
  }
  claim_allowed <- as.logical(data$claim_allowed_model)
  primary_stable <- as.logical(data$primary_model_stable)
  expected_gate <- ifelse(
    source_model_valid %in% TRUE & claim_allowed %in% TRUE &
      primary_stable %in% TRUE,
    "eligible_for_readiness_assessment",
    ifelse(
      !source_model_valid %in% TRUE | !primary_stable %in% TRUE,
      "diagnostic_only_model",
      "not_claim_allowed_model"
    )
  )
  bad_gate <- is.na(data$claim_gate) | data$claim_gate != expected_gate
  if (any(bad_gate)) {
    .wgcna_group_effect_consumer_stop_rows(
      paste0(artifact, " changed the deterministic claim gate."),
      which(bad_gate)
    )
  }

  semantic_key <- c(
    "dataset", "entity_level", "entity_id", "analysis_tier", "contrast",
    "effect_scope", "spatial_unit", "test_type"
  )
  if (anyDuplicated(data[semantic_key])) {
    stop(artifact, " contains duplicated inferential endpoints.", call. = FALSE)
  }
  invisible(TRUE)
}

wgcna_inferential_handoff_supermodule_display <- function(
    handoff, module_supermodule_lookup
) {
  wgcna_inferential_handoff_validate(handoff)
  required_lookup <- c(
    "dataset", "module_id", "supermodule_id", "n_member_modules",
    "supermodule_label"
  )
  missing <- setdiff(required_lookup, names(module_supermodule_lookup))
  if (length(missing)) {
    stop(
      "Supermodule display lookup is missing column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  lookup <- unique(module_supermodule_lookup[required_lookup])
  singleton_lookup <- lookup[
    !is.na(lookup$n_member_modules) & lookup$n_member_modules == 1L,
    ,
    drop = FALSE
  ]
  if (anyDuplicated(singleton_lookup[c("dataset", "supermodule_id")])) {
    stop(
      "Singleton supermodule display lookup is not one-to-one.",
      call. = FALSE
    )
  }

  genuine <- handoff[
    handoff$entity_level == "supermodule",
    ,
    drop = FALSE
  ]
  genuine$supermodule_id <- genuine$entity_id
  genuine$supermodule_label <- genuine$display_label
  genuine$display_n_member_modules <- genuine$n_member_modules
  genuine$source_entity_level <- genuine$entity_level
  genuine$source_entity_id <- genuine$entity_id
  genuine$display_is_compatibility_alias <- FALSE
  genuine$display_is_independent_endpoint <- TRUE
  genuine$display_support_origin <- "independent_supermodule_endpoint"
  genuine$display_entity_role <- "independent_supermodule_endpoint"

  module <- handoff[handoff$entity_level == "module", , drop = FALSE]
  singleton <- merge(
    module, singleton_lookup,
    by.x = c("dataset", "entity_id"),
    by.y = c("dataset", "module_id"),
    all = FALSE, sort = FALSE, suffixes = c("", "_lookup")
  )
  if (nrow(singleton)) {
    if ("supermodule_id_lookup" %in% names(singleton)) {
      singleton$supermodule_id <- singleton$supermodule_id_lookup
    }
    if ("supermodule_label_lookup" %in% names(singleton)) {
      singleton$supermodule_label <- singleton$supermodule_label_lookup
    }
    singleton$display_n_member_modules <- 1L
    singleton$source_entity_level <- singleton$entity_level
    singleton$source_entity_id <- singleton$entity_id
    singleton$display_is_compatibility_alias <- TRUE
    singleton$display_is_independent_endpoint <- FALSE
    singleton$display_support_origin <- "inherited_canonical_module"
    singleton$display_entity_role <- "singleton_compatibility_alias"
  }

  out <- dplyr::bind_rows(genuine, singleton)
  if (anyDuplicated(out[c(
    "dataset", "supermodule_id", "analysis_tier", "contrast",
    "effect_scope", "spatial_unit", "test_type"
  )])) {
    stop(
      "Supermodule display contains duplicated displayed endpoints.",
      call. = FALSE
    )
  }
  out
}
