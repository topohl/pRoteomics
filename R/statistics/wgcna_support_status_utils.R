# Canonical WGCNA multiplicity-support classification.

wgcna_group_classify_statistical_support <- function(
    p_value, applicable_fdr
) {
  if (!is.numeric(p_value) || !is.numeric(applicable_fdr)) {
    stop("p_value and applicable_fdr must be numeric.", call. = FALSE)
  }
  if (length(p_value) != length(applicable_fdr)) {
    stop("p_value and applicable_fdr must have equal lengths.", call. = FALSE)
  }
  invalid_fdr <- is.finite(applicable_fdr) &
    (applicable_fdr < 0 | applicable_fdr > 1)
  if (any(invalid_fdr)) {
    stop("Finite applicable_fdr values must lie in [0, 1].", call. = FALSE)
  }

  status <- rep("not_supported", length(applicable_fdr))
  finite_fdr <- is.finite(applicable_fdr)
  status[finite_fdr & applicable_fdr <= 0.10] <- "suggestive_FDR10"
  status[finite_fdr & applicable_fdr <= 0.05] <- "FDR_supported"
  status
}

wgcna_group_applicable_fdr <- function(data) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  required <- c(
    "analysis_tier", "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "WGCNA support input is missing required column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }

  tier <- as.character(data$analysis_tier)
  applicable_fdr <- rep(NA_real_, nrow(data))
  applicable_fdr[tier == "primary_wgcna_global"] <-
    data$FDR_primary_global[tier == "primary_wgcna_global"]
  applicable_fdr[tier == "secondary_contextual_global"] <-
    data$FDR_secondary_global[tier == "secondary_contextual_global"]
  applicable_fdr[tier == "secondary_spatial_heterogeneity"] <-
    data$FDR_interaction_omnibus[
      tier == "secondary_spatial_heterogeneity"
    ]
  applicable_fdr[tier == "exploratory_spatial_localization"] <-
    data$FDR_local_exploratory[
      tier == "exploratory_spatial_localization"
    ]
  applicable_fdr
}

wgcna_group_support_status_impact_audit <- function(
    data, source_artifact = NA_character_
) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  required <- c(
    "p_value", "statistical_support_status", "analysis_tier", "test_type",
    "claim_entity_role", "independent_hypothesis"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      "WGCNA impact-audit input is missing required column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }

  applicable_fdr <- wgcna_group_applicable_fdr(data)
  new_status <- wgcna_group_classify_statistical_support(
    data$p_value, applicable_fdr
  )
  auditable <- data$independent_hypothesis %in% TRUE &
    data$claim_entity_role != "compatibility_alias" &
    data$test_type != "conditional_interaction_followup"
  changed <- which(
    auditable &
      as.character(data$statistical_support_status) != new_status
  )
  identity_columns <- intersect(
    c(
      "dataset", "level", "endpoint_id", "canonical_claim_entity_id",
      "contrast", "analysis_tier", "test_type", "effect_scope",
      "spatial_unit", "p_value", "display_allowed",
      "manuscript_claim_ready"
    ),
    names(data)
  )
  out <- data[changed, identity_columns, drop = FALSE]
  out$applicable_fdr <- applicable_fdr[changed]
  out$old_status <- as.character(data$statistical_support_status[changed])
  out$new_status <- new_status[changed]
  out$source_artifact <- rep(as.character(source_artifact), nrow(out))
  rownames(out) <- NULL
  out
}
