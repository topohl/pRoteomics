# Additive diagnostics for the strict GSEA-program <-> WGCNA concordance.
#
# These helpers explain existing classifications and create explicitly
# descriptive direction summaries. They never modify source inference,
# concordance classes, p-values, FDRs, or claim gates.

gwwd_contrast_suffix <- function(x) {
  out <- c(
    "RES - CON" = "RES_CON",
    "SUS - CON" = "SUS_CON",
    "SUS - RES" = "SUS_RES"
  )[as.character(x)]
  unname(out)
}

gwwd_stable_model <- function(valid, stability) {
  valid %in% TRUE & !grepl(
    "invalid|unstable|fail|nonestimable",
    tolower(as.character(stability))
  )
}

gwwd_percentile <- function(x) {
  if (!length(x)) return(numeric())
  if (length(x) == 1L) return(1)
  (rank(x, ties.method = "average", na.last = "keep") - 1) /
    (sum(!is.na(x)) - 1)
}

gwwd_effect_distribution_lookup <- function(handoff) {
  gww_assert_columns(
    handoff,
    c(
      "dataset", "entity_level", "analysis_tier", "effect_scope",
      "contrast", "estimate", "CI_low", "CI_high", "model_valid",
      "source_key"
    ),
    "Canonical WGCNA handoff"
  )
  out <- handoff |>
    dplyr::filter(
      .data$contrast %in% unname(gww_formal_contrast_map()),
      .data$model_valid %in% TRUE,
      is.finite(.data$estimate), is.finite(.data$CI_low),
      is.finite(.data$CI_high)
    ) |>
    dplyr::mutate(
      absolute_estimate = abs(.data$estimate),
      CI_width = .data$CI_high - .data$CI_low
    ) |>
    dplyr::group_by(
      .data$dataset, .data$entity_level, .data$analysis_tier,
      .data$effect_scope
    ) |>
    dplyr::mutate(
      estimate_percentile = gwwd_percentile(.data$absolute_estimate),
      CI_width_percentile = gwwd_percentile(.data$CI_width),
      CI_width_extreme_threshold = as.numeric(stats::quantile(
        .data$CI_width, probs = 0.90, names = FALSE, type = 7,
        na.rm = TRUE
      )),
      effect_distribution_n = dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      wgcna_source_key = as.character(.data$source_key),
      absolute_estimate,
      estimate_percentile,
      CI_width,
      CI_width_percentile,
      CI_width_extreme_threshold,
      effect_distribution_n
    )
  gww_assert_unique(out, "wgcna_source_key", "WGCNA effect percentile lookup")
  out
}

gwwd_concordant_imprecise <- function(long, effect_lookup) {
  gww_assert_columns(
    long,
    c(
      "concordance_class", "gsea_evidence_id", "dataset", "contrast",
      "comparison_scope", "gsea_spatial_unit", "biological_program",
      "GO_ID", "GO_description", "NES", "GSEA_FDR",
      "wgcna_entity_level", "wgcna_entity_id", "wgcna_display_label",
      "wgcna_estimate", "wgcna_SE", "wgcna_CI_low", "wgcna_CI_high",
      "wgcna_p_value", "wgcna_tier_specific_fdr",
      "wgcna_biological_n", "wgcna_model_stability_status",
      "animal_instability_flag", "animal_instability_reason", "loao_status",
      "effect_scale_q25_abs", "n_overlap", "Jaccard",
      "overlap_fisher_p_value", "overlap_FDR", "wgcna_source_key"
    ),
    "Strict concordance long table"
  )
  counts_by_evidence <- long |>
    dplyr::group_by(
      .data$dataset, .data$contrast, .data$comparison_scope,
      .data$gsea_evidence_id
    ) |>
    dplyr::summarise(
      n_WGCNA_entities_matched_to_GSEA_evidence = dplyr::n_distinct(paste(
        .data$wgcna_entity_level, .data$wgcna_entity_id, sep = "|"
      )),
      .groups = "drop"
    )
  counts_by_entity <- long |>
    dplyr::group_by(
      .data$dataset, .data$contrast, .data$comparison_scope,
      .data$wgcna_entity_level, .data$wgcna_entity_id
    ) |>
    dplyr::summarise(
      n_GSEA_evidence_rows_matched_to_WGCNA_entity =
        dplyr::n_distinct(.data$gsea_evidence_id),
      .groups = "drop"
    )

  out <- long |>
    dplyr::filter(.data$concordance_class == "concordant_imprecise") |>
    dplyr::left_join(
      effect_lookup,
      by = "wgcna_source_key",
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      counts_by_evidence,
      by = c("dataset", "contrast", "comparison_scope", "gsea_evidence_id"),
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      counts_by_entity,
      by = c(
        "dataset", "contrast", "comparison_scope",
        "wgcna_entity_level", "wgcna_entity_id"
      ),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      CI_crosses_zero = .data$wgcna_CI_low <= 0 & .data$wgcna_CI_high >= 0,
      overlap_test_available = !is.na(.data$overlap_FDR),
      overlap_FDR_supported = .data$overlap_FDR <= 0.05,
      broad_CI_or_weak_precision_flag =
        .data$CI_width_percentile >= 0.90,
      stronger_imprecise_support_flag =
        .data$overlap_FDR_supported %in% TRUE,
      directional_only_imprecise_flag =
        !(.data$overlap_FDR_supported %in% TRUE) &
        !.data$broad_CI_or_weak_precision_flag,
      diagnostic_subgroup = dplyr::case_when(
        .data$stronger_imprecise_support_flag ~
          "stronger_imprecise_support",
        .data$broad_CI_or_weak_precision_flag ~
          "broad_CI_or_weak_precision",
        .data$directional_only_imprecise_flag ~
          "directional_only_imprecise",
        TRUE ~ "other"
      ),
      diagnostic_subgroup_rule = paste(
        "Official class remains concordant_imprecise. stronger=overlap",
        "FDR<=0.05; broad=CI-width percentile>=0.90 within the same",
        "dataset/entity-level/tier/effect-scope distribution; directional=",
        "no supported overlap and not broad; structural overlap is not",
        "animal-level group-effect evidence."
      )
    ) |>
    dplyr::transmute(
      dataset,
      contrast,
      local_global_scope = .data$comparison_scope,
      spatial_unit = .data$gsea_spatial_unit,
      biological_program,
      GSEA_GO_ID = .data$GO_ID,
      GSEA_GO_term = .data$GO_description,
      NES,
      GSEA_FDR,
      WGCNA_entity_level = .data$wgcna_entity_level,
      WGCNA_entity_ID = .data$wgcna_entity_id,
      WGCNA_entity_label = .data$wgcna_display_label,
      WGCNA_estimate = .data$wgcna_estimate,
      WGCNA_SE = .data$wgcna_SE,
      WGCNA_CI_low = .data$wgcna_CI_low,
      WGCNA_CI_high = .data$wgcna_CI_high,
      WGCNA_raw_p = .data$wgcna_p_value,
      WGCNA_tier_specific_FDR = .data$wgcna_tier_specific_fdr,
      biological_n = .data$wgcna_biological_n,
      model_stability = .data$wgcna_model_stability_status,
      animal_instability_flag,
      animal_stability_detail = .data$animal_instability_reason,
      LOAO_status = .data$loao_status,
      absolute_estimate,
      near_zero_threshold_used = .data$effect_scale_q25_abs,
      estimate_percentile,
      effect_distribution_n,
      CI_crosses_zero,
      CI_width,
      CI_width_percentile,
      CI_width_extreme_threshold,
      overlap_n = .data$n_overlap,
      overlap_Jaccard = .data$Jaccard,
      overlap_p = .data$overlap_fisher_p_value,
      overlap_FDR,
      overlap_test_available,
      n_WGCNA_entities_matched_to_GSEA_evidence,
      n_GSEA_evidence_rows_matched_to_WGCNA_entity,
      stronger_imprecise_support_flag,
      directional_only_imprecise_flag,
      broad_CI_or_weak_precision_flag,
      diagnostic_subgroup,
      official_concordance_class = .data$concordance_class,
      diagnostic_subgroup_rule,
      gsea_evidence_id,
      wgcna_source_key
    )
  gww_assert_unique(
    out,
    c("gsea_evidence_id", "WGCNA_entity_level", "WGCNA_entity_ID"),
    "Concordant-imprecise diagnostic"
  )
  out
}

gwwd_power_diagnostic <- function(long) {
  gww_assert_columns(
    long,
    c(
      "wgcna_source_key", "dataset", "contrast", "comparison_scope",
      "wgcna_spatial_unit", "wgcna_entity_level", "wgcna_entity_id",
      "wgcna_display_label", "wgcna_estimate", "wgcna_SE",
      "wgcna_CI_low", "wgcna_CI_high", "wgcna_p_value",
      "wgcna_tier_specific_fdr", "wgcna_tier_specific_family_id",
      "wgcna_tier_specific_family_size", "wgcna_biological_n",
      "wgcna_model_stability_status", "biological_program",
      "gsea_evidence_id", "overlap_FDR"
    ),
    "Strict concordance long table"
  )
  out <- long |>
    dplyr::group_by(
      .data$wgcna_source_key, .data$dataset, .data$contrast,
      .data$comparison_scope, .data$wgcna_spatial_unit,
      .data$wgcna_entity_level, .data$wgcna_entity_id,
      .data$wgcna_display_label, .data$wgcna_estimate, .data$wgcna_SE,
      .data$wgcna_CI_low, .data$wgcna_CI_high, .data$wgcna_p_value,
      .data$wgcna_tier_specific_fdr,
      .data$wgcna_tier_specific_family_id,
      .data$wgcna_tier_specific_family_size,
      .data$wgcna_biological_n, .data$wgcna_model_stability_status
    ) |>
    dplyr::summarise(
      n_matched_GSEA_evidence_rows = dplyr::n_distinct(.data$gsea_evidence_id),
      n_matched_programs = dplyr::n_distinct(.data$biological_program),
      matched_programs = paste(
        sort(unique(.data$biological_program)), collapse = ";"
      ),
      best_overlap_FDR = if (all(is.na(.data$overlap_FDR))) NA_real_ else
        min(.data$overlap_FDR, na.rm = TRUE),
      n_overlap_FDR_supported = sum(.data$overlap_FDR <= 0.05, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      absolute_estimate = abs(.data$wgcna_estimate),
      CI_excludes_zero = .data$wgcna_CI_low > 0 | .data$wgcna_CI_high < 0,
      raw_p_below_0_05 = .data$wgcna_p_value < 0.05,
      FDR_at_or_below_0_05 = .data$wgcna_tier_specific_fdr <= 0.05,
      FDR_at_or_below_0_10_descriptive =
        .data$wgcna_tier_specific_fdr <= 0.10,
      any_overlap_FDR_supported = .data$n_overlap_FDR_supported > 0,
      WGCNA_FDR_bin = dplyr::case_when(
        .data$wgcna_tier_specific_fdr <= 0.05 ~ "[0,0.05]",
        .data$wgcna_tier_specific_fdr <= 0.10 ~ "(0.05,0.10]",
        .data$wgcna_tier_specific_fdr <= 0.25 ~ "(0.10,0.25]",
        .data$wgcna_tier_specific_fdr <= 0.50 ~ "(0.25,0.50]",
        is.finite(.data$wgcna_tier_specific_fdr) ~ "(0.50,1.00]",
        TRUE ~ "missing"
      )
    ) |>
    dplyr::arrange(
      .data$wgcna_tier_specific_fdr,
      dplyr::desc(.data$absolute_estimate),
      .data$best_overlap_FDR,
      .data$dataset, .data$wgcna_entity_level, .data$wgcna_entity_id
    ) |>
    dplyr::mutate(
      strongest_matched_endpoint_rank = dplyr::row_number(),
      top_20_strongest_matched_endpoint =
        .data$strongest_matched_endpoint_rank <= 20L,
      ranking_rule = paste(
        "ascending unchanged tier-specific WGCNA FDR; descending absolute",
        "effect; ascending best program-specific overlap FDR"
      )
    )
  gww_assert_unique(out, "wgcna_source_key", "Matched WGCNA power diagnostic")
  out
}

gwwd_prepare_all_gsea <- function(term_table) {
  gww_assert_columns(
    term_table,
    c(
      "dataset", "phenotype_contrast", "spatial_unit", "GO_ID",
      "GO_description", "NES", "raw_p", "GSEA_FDR", "theme_id",
      "manuscript_theme", "theme_role", "theme_claim_eligible",
      "theme_assignment_id"
    ),
    "Ontology-aware all-contrast GSEA theme assignment table"
  )
  if (!"registry_version" %in% names(term_table)) {
    term_table$registry_version <- NA_character_
  }
  term_table |>
    dplyr::mutate(
      contrast = gww_formal_contrast(.data$phenotype_contrast),
      gsea_spatial_unit = gww_normalize_spatial_unit(.data$spatial_unit),
      biological_program = as.character(.data$theme_id),
      ID = as.character(.data$GO_ID),
      Description = as.character(.data$GO_description),
      NES = as.numeric(.data$NES),
      GSEA_FDR = as.numeric(.data$GSEA_FDR),
      raw_p = as.numeric(.data$raw_p)
    ) |>
    dplyr::filter(
      .data$contrast %in% unname(gww_formal_contrast_map()),
      .data$theme_claim_eligible %in% TRUE,
      .data$theme_role %in% c("primary", "supporting"),
      !is.na(.data$biological_program), nzchar(.data$biological_program),
      is.finite(.data$NES), is.finite(.data$GSEA_FDR)
    )
}

gwwd_local_program_effects <- function(term_table) {
  x <- gwwd_prepare_all_gsea(term_table)
  common_terms <- x |>
    dplyr::distinct(
      .data$dataset, .data$gsea_spatial_unit, .data$biological_program,
      .data$contrast, .data$ID
    ) |>
    dplyr::count(
      .data$dataset, .data$gsea_spatial_unit, .data$biological_program,
      .data$ID, name = "n_contrasts_available"
    ) |>
    dplyr::filter(.data$n_contrasts_available == 3L)

  anchor <- x |>
    dplyr::inner_join(
      common_terms,
      by = c("dataset", "gsea_spatial_unit", "biological_program", "ID"),
      relationship = "many-to-one"
    ) |>
    dplyr::group_by(
      .data$dataset, .data$gsea_spatial_unit,
      .data$biological_program, .data$ID, .data$Description
    ) |>
    dplyr::summarise(
      anchor_mean_FDR_across_contrasts = mean(.data$GSEA_FDR),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$dataset, .data$gsea_spatial_unit, .data$biological_program,
      .data$anchor_mean_FDR_across_contrasts, .data$ID, .data$Description
    ) |>
    dplyr::group_by(
      .data$dataset, .data$gsea_spatial_unit, .data$biological_program
    ) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::rename(
      anchor_GO_ID = "ID",
      anchor_GO_description = "Description"
    )

  effects <- x |>
    dplyr::inner_join(
      common_terms,
      by = c("dataset", "gsea_spatial_unit", "biological_program", "ID"),
      relationship = "many-to-one"
    ) |>
    dplyr::group_by(
      .data$dataset, .data$gsea_spatial_unit,
      .data$biological_program, .data$contrast
    ) |>
    dplyr::summarise(
      manuscript_theme = dplyr::first(.data$manuscript_theme),
      theme_role = dplyr::first(.data$theme_role),
      registry_version = dplyr::first(.data$registry_version),
      n_common_GO_terms = dplyr::n_distinct(.data$ID),
      gsea_program_NES = stats::median(.data$NES),
      gsea_min_FDR = min(.data$GSEA_FDR),
      gsea_median_FDR = stats::median(.data$GSEA_FDR),
      n_GSEA_FDR_supported_terms = sum(.data$GSEA_FDR < 0.05),
      fraction_common_terms_positive = mean(.data$NES > 0),
      fraction_common_terms_negative = mean(.data$NES < 0),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      anchor,
      by = c("dataset", "gsea_spatial_unit", "biological_program"),
      relationship = "many-to-one"
    )
  anchor_values <- x |>
    dplyr::inner_join(
      anchor |>
        dplyr::select(dplyr::all_of(c(
          "dataset", "gsea_spatial_unit", "biological_program",
          "anchor_GO_ID"
        ))),
      by = c(
        "dataset", "gsea_spatial_unit", "biological_program",
        "ID" = "anchor_GO_ID"
      ),
      relationship = "many-to-one"
    ) |>
    dplyr::transmute(
      dataset, gsea_spatial_unit, biological_program, contrast,
      anchor_NES = .data$NES,
      anchor_raw_p = .data$raw_p,
      anchor_GSEA_FDR = .data$GSEA_FDR
    )
  out <- effects |>
    dplyr::left_join(
      anchor_values,
      by = c(
        "dataset", "gsea_spatial_unit", "biological_program", "contrast"
      ),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      gsea_direction_sign = sign(.data$gsea_program_NES),
      gsea_representation_method = paste(
        "DESCRIPTIVE: median NES over GO IDs present in all three contrasts",
        "within the same dataset/spatial-unit/program. One anchor GO ID is",
        "chosen once across contrasts by lowest mean canonical FDR, then GO ID;",
        "the anchor does not determine the direction label."
      )
    )
  gww_assert_unique(
    out,
    c("dataset", "gsea_spatial_unit", "biological_program", "contrast"),
    "Complete local GSEA program effects"
  )
  out
}

gwwd_recurrent_cross_spatial_program_effects <- function(local_effects) {
  for (column in c("manuscript_theme", "theme_role", "registry_version")) {
    if (!column %in% names(local_effects)) local_effects[[column]] <- NA_character_
  }
  out <- local_effects |>
    dplyr::group_by(
      .data$dataset, .data$biological_program, .data$contrast
    ) |>
    dplyr::summarise(
      manuscript_theme = dplyr::first(.data$manuscript_theme),
      theme_role = dplyr::first(.data$theme_role),
      registry_version = dplyr::first(.data$registry_version),
      n_spatial_units = dplyr::n_distinct(.data$gsea_spatial_unit),
      n_units_positive = sum(.data$gsea_program_NES > 0),
      n_units_negative = sum(.data$gsea_program_NES < 0),
      n_units_with_any_FDR_supported_term = sum(
        .data$n_GSEA_FDR_supported_terms > 0
      ),
      common_GO_terms_min = min(.data$n_common_GO_terms),
      common_GO_terms_max = max(.data$n_common_GO_terms),
      gsea_program_NES = stats::median(.data$gsea_program_NES),
      gsea_min_FDR = min(.data$gsea_min_FDR),
      gsea_median_FDR = stats::median(.data$gsea_median_FDR),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      gsea_direction_sign = dplyr::case_when(
        .data$n_units_positive >= 2L &
          .data$n_units_positive > .data$n_units_negative ~ 1,
        .data$n_units_negative >= 2L &
          .data$n_units_negative > .data$n_units_positive ~ -1,
        TRUE ~ NA_real_
      ),
      recurrent_cross_spatial_direction_available = !is.na(
        .data$gsea_direction_sign
      ),
      gsea_representation_method = paste(
        "DESCRIPTIVE recurrent direction from the majority sign of local",
        "contrast-independent common-term program median NES values; at least",
        "two spatial units and a non-tied sign majority are required; FDR is",
        "reported but not used to assign direction."
      )
    )
  gww_assert_unique(
    out,
    c("dataset", "biological_program", "contrast"),
    "Complete recurrent-cross-spatial GSEA program effects"
  )
  out
}

gwwd_wgcna_program_effects <- function(handoff, matches, scope = c("global", "local")) {
  scope <- match.arg(scope)
  modules <- matches |>
    dplyr::filter(.data$entity_level == "module") |>
    dplyr::distinct(
      .data$dataset, .data$entity_id, .data$biological_program
    )
  x <- handoff |>
    dplyr::filter(
      .data$entity_level == "module",
      .data$contrast %in% unname(gww_formal_contrast_map())
    ) |>
    dplyr::inner_join(
      modules,
      by = c("dataset", "entity_id"),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      wgcna_spatial_unit = gww_normalize_spatial_unit(.data$spatial_unit),
      stable_model = gwwd_stable_model(
        .data$model_valid, .data$model_stability_status
      )
    )
  if (scope == "global") {
    x <- x |>
      dplyr::filter(
        .data$effect_scope == "spatial_adjusted_global",
        .data$analysis_tier %in% c(
          "primary_wgcna_global", "secondary_contextual_global"
        ),
        .data$wgcna_spatial_unit == "global_spatial_adjusted"
      )
    keys <- c("dataset", "biological_program", "contrast")
  } else {
    x <- x |>
      dplyr::filter(
        .data$effect_scope == "within_spatial_unit",
        .data$analysis_tier == "exploratory_spatial_localization"
      )
    keys <- c(
      "dataset", "wgcna_spatial_unit", "biological_program", "contrast"
    )
  }
  out <- x |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(
      n_matched_modules = dplyr::n_distinct(.data$entity_id),
      matched_module_IDs = paste(sort(unique(.data$entity_id)), collapse = ";"),
      wgcna_program_estimate = stats::median(.data$estimate),
      wgcna_CI_envelope_low = min(.data$CI_low),
      wgcna_CI_envelope_high = max(.data$CI_high),
      wgcna_min_raw_p = min(.data$p_value),
      wgcna_min_tier_FDR = min(.data$tier_specific_fdr),
      n_stable_models = sum(.data$stable_model),
      all_models_stable = all(.data$stable_model),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      wgcna_direction_sign = sign(.data$wgcna_program_estimate)
    )
  gww_assert_unique(out, keys, paste("WGCNA", scope, "program effects"))
  out
}

gwwd_pattern_label <- function(res, sus, direct, status_res, status_sus,
                               status_direct) {
  out <- rep("insufficient_data", length(res))
  incompatible <- status_res == "incompatible" |
    status_sus == "incompatible" | status_direct == "incompatible"
  out[incompatible %in% TRUE] <- "inconsistent"
  complete <- status_res == "compatible" & status_sus == "compatible" &
    status_direct == "compatible"
  out[complete & res == 1 & sus == -1 & direct == -1] <-
    "RES_gt_CON_gt_SUS"
  out[complete & res == -1 & sus == 1 & direct == 1] <-
    "RES_lt_CON_lt_SUS"
  out[complete & res == sus] <- "shared_same_direction"
  out[complete & res == -sus & !(
    (res == 1 & sus == -1 & direct == -1) |
      (res == -1 & sus == 1 & direct == 1)
  )] <- "opposing_RES_SUS"
  other_complete <- complete & out == "insufficient_data"
  out[other_complete] <- "inconsistent"
  out
}

gwwd_directional_patterns <- function(gsea_effects, wgcna_effects,
                                      scope = c("global", "local")) {
  scope <- match.arg(scope)
  if (scope == "global") {
    join_keys <- c("dataset", "biological_program", "contrast")
    key_columns <- c("dataset", "biological_program")
  } else {
    join_keys <- c(
      "dataset", "gsea_spatial_unit" = "wgcna_spatial_unit",
      "biological_program", "contrast"
    )
    key_columns <- c("dataset", "gsea_spatial_unit", "biological_program")
  }
  theme_metadata <- intersect(
    c("manuscript_theme", "theme_role", "registry_version"),
    names(gsea_effects)
  )
  key_columns <- c(key_columns, theme_metadata)
  joined <- gsea_effects |>
    dplyr::left_join(
      wgcna_effects,
      by = join_keys,
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      direction_compatibility = dplyr::case_when(
        is.na(.data$gsea_direction_sign) |
          is.na(.data$wgcna_direction_sign) ~ "unavailable",
        .data$gsea_direction_sign == .data$wgcna_direction_sign ~ "compatible",
        TRUE ~ "incompatible"
      ),
      integrated_direction_sign = dplyr::if_else(
        .data$direction_compatibility == "compatible",
        as.numeric(.data$gsea_direction_sign), NA_real_
      ),
      descriptive_evidence_label = paste(
        "DESCRIPTIVE / HYPOTHESIS-GENERATING; NOT AN INFERENTIAL CLAIM",
        "CLASSIFIER"
      )
    )
  value_columns <- setdiff(
    names(joined), c(key_columns, "contrast")
  )
  wide <- joined |>
    dplyr::mutate(contrast_suffix = gwwd_contrast_suffix(.data$contrast)) |>
    dplyr::select(-"contrast") |>
    tidyr::pivot_wider(
      names_from = "contrast_suffix",
      values_from = dplyr::all_of(value_columns),
      names_glue = "{.value}__{contrast_suffix}"
    )
  needed <- c("RES_CON", "SUS_CON", "SUS_RES")
  for (base in c("integrated_direction_sign", "direction_compatibility")) {
    for (suffix in needed) {
      column <- paste(base, suffix, sep = "__")
      if (!column %in% names(wide)) {
        wide[[column]] <- if (base == "direction_compatibility")
          "unavailable" else NA_real_
      }
    }
  }
  wide$directional_pattern <- gwwd_pattern_label(
    wide$integrated_direction_sign__RES_CON,
    wide$integrated_direction_sign__SUS_CON,
    wide$integrated_direction_sign__SUS_RES,
    wide$direction_compatibility__RES_CON,
    wide$direction_compatibility__SUS_CON,
    wide$direction_compatibility__SUS_RES
  )
  wide$pattern_method <- paste(
    "Directional labels require GSEA/WGCNA sign compatibility for all three",
    "contrasts. RES/SUS supported-shift labels are reserved and are not assigned",
    "from significance asymmetry or an invented near-zero GSEA threshold.",
    "FDRs and CIs remain descriptive columns and are not pattern gates."
  )
  gww_assert_unique(wide, key_columns, paste(scope, "directional patterns"))
  wide
}

gwwd_strict_gate_audit <- function(strict_long, strict_patterns) {
  recurrent_cross_spatial <- strict_long |>
    dplyr::filter(
      .data$comparison_scope == "recurrent_cross_spatial",
      .data$wgcna_entity_level == "module"
    ) |>
    dplyr::group_by(
      .data$dataset, .data$biological_program, .data$contrast
    ) |>
    dplyr::summarise(
      GSEA_available = dplyr::n() > 0,
      GSEA_direction_sign = dplyr::first(.data$gsea_direction_sign),
      GSEA_direction_status = dplyr::first(.data$gsea_direction_status),
      WGCNA_endpoint_available = dplyr::n_distinct(.data$wgcna_entity_id) > 0,
      n_WGCNA_endpoints = dplyr::n_distinct(.data$wgcna_entity_id),
      n_stable_concordant_endpoints = sum(.data$concordance_class %in% c(
        "FDR_supported_concordance", "concordant_imprecise"
      )),
      WGCNA_effects_stable = all(gwwd_stable_model(
        .data$wgcna_model_valid, .data$wgcna_model_stability_status
      )),
      animal_stability_pass = !any(.data$animal_instability_flag %in% TRUE),
      overlap_evidence_available = any(!is.na(.data$overlap_FDR)),
      overlap_FDR_supported = any(.data$overlap_FDR <= 0.05, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      suffix = gwwd_contrast_suffix(.data$contrast),
      stable_concordant_endpoint_available =
        .data$n_stable_concordant_endpoints > 0
    ) |>
    dplyr::select(-"contrast") |>
    tidyr::pivot_wider(
      names_from = "suffix",
      values_from = dplyr::all_of(c(
        "GSEA_available", "GSEA_direction_sign", "GSEA_direction_status",
        "WGCNA_endpoint_available", "n_WGCNA_endpoints",
        "n_stable_concordant_endpoints", "WGCNA_effects_stable",
        "animal_stability_pass", "overlap_evidence_available",
        "overlap_FDR_supported", "stable_concordant_endpoint_available"
      )),
      names_glue = "{.value}__{suffix}"
    )
  out <- strict_patterns |>
    dplyr::select(dplyr::all_of(c(
      "dataset", "biological_program", "adaptive_resilience_pattern"
    ))) |>
    dplyr::left_join(
      recurrent_cross_spatial,
      by = c("dataset", "biological_program"),
      relationship = "one-to-one"
    )
  bases_logical <- c(
    "GSEA_available", "WGCNA_endpoint_available", "WGCNA_effects_stable",
    "animal_stability_pass", "overlap_evidence_available",
    "overlap_FDR_supported", "stable_concordant_endpoint_available"
  )
  for (base in bases_logical) {
    for (suffix in c("RES_CON", "SUS_CON", "SUS_RES")) {
      column <- paste(base, suffix, sep = "__")
      if (!column %in% names(out)) out[[column]] <- FALSE
      out[[column]][is.na(out[[column]])] <- FALSE
    }
  }
  for (suffix in c("RES_CON", "SUS_CON", "SUS_RES")) {
    column <- paste("GSEA_direction_sign", suffix, sep = "__")
    if (!column %in% names(out)) out[[column]] <- NA_real_
  }
  res <- out$GSEA_direction_sign__RES_CON
  sus <- out$GSEA_direction_sign__SUS_CON
  direct <- out$GSEA_direction_sign__SUS_RES
  res_ok <- out$stable_concordant_endpoint_available__RES_CON
  sus_ok <- out$stable_concordant_endpoint_available__SUS_CON
  direct_ok <- out$stable_concordant_endpoint_available__SUS_RES
  out$recurrent_cross_spatial_GSEA_all_three_available <-
    out$GSEA_available__RES_CON & out$GSEA_available__SUS_CON &
    out$GSEA_available__SUS_RES
  out$program_module_mapping_available <- apply(
    cbind(
      out$WGCNA_endpoint_available__RES_CON,
      out$WGCNA_endpoint_available__SUS_CON,
      out$WGCNA_endpoint_available__SUS_RES
    ), 1L, any
  )
  out$direct_SUS_RES_direction_compatible <- !is.na(direct) & direct_ok
  out$RES_departure_from_CON_compatible <- !is.na(res) & res_ok
  out$SUS_departure_from_CON_compatible <- !is.na(sus) & sus_ok
  out$strict_RES_supported_shift_rule_pass <-
    res == 1 & direct == -1 & res_ok & direct_ok
  out$strict_SUS_supported_shift_rule_pass <-
    !is.na(sus) & sus != 0 & direct == sus & sus_ok & direct_ok
  out$strict_opposing_rule_pass <-
    !is.na(res) & !is.na(sus) & res == -sus & res_ok & sus_ok &
    direct == sus & direct_ok
  out$strict_graded_rule_pass <-
    !is.na(res) & !is.na(sus) & res == sus & !is.na(direct) &
    direct != 0 & res_ok & sus_ok & direct_ok
  out$strict_shared_rule_pass <-
    !is.na(res) & !is.na(sus) & res == sus & res_ok & sus_ok &
    (!direct_ok | is.na(direct))
  out$strict_classifier_significance_gated <- TRUE
  out$strict_GSEA_gate <- paste(
    "gww_prepare_supported_terms requires canonical GSEA p.adjust < 0.05",
    "before local/recurrent evidence enters the strict adaptive classifier"
  )
  out$overlap_is_strict_gate <- FALSE
  out$final_failed_gates <- vapply(seq_len(nrow(out)), function(i) {
    failures <- character()
    for (suffix in c("RES_CON", "SUS_CON", "SUS_RES")) {
      if (!isTRUE(out[[paste0("GSEA_available__", suffix)]][[i]])) {
        failures <- c(failures, paste0("missing_FDR_supported_recurrent_GSEA_", suffix))
      }
      if (!isTRUE(out[[paste0(
        "stable_concordant_endpoint_available__", suffix
      )]][[i]])) {
        failures <- c(failures, paste0("no_stable_concordant_WGCNA_", suffix))
      }
    }
    rules <- c(
      out$strict_RES_supported_shift_rule_pass[[i]],
      out$strict_SUS_supported_shift_rule_pass[[i]],
      out$strict_opposing_rule_pass[[i]],
      out$strict_graded_rule_pass[[i]],
      out$strict_shared_rule_pass[[i]]
    )
    if (!any(rules %in% TRUE)) failures <- c(failures, "no_strict_direction_rule_passed")
    paste(unique(failures), collapse = ";")
  }, character(1))
  gww_assert_unique(
    out, c("dataset", "biological_program"), "Strict adaptive gate audit"
  )
  out
}
