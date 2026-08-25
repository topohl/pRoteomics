TARGETED_EMPIRICAL_MARKER_CONTRACT <- "empirical_roi_marker_v2_animal_paired_limma"
TARGETED_CLAIM_READY_CONTRACT <- "microglia_targeted_claim_ready_v2"

targeted_effect_statistic_type <- function(enrichment_method) {
  dplyr::case_when(
    enrichment_method == "limma_ranked_geneSetTest" ~ "standardized_mean_rank_stat",
    enrichment_method %in% c("fgsea", "clusterProfiler_GSEA") ~ "NES",
    TRUE ~ "method_specific_effect_statistic"
  )
}

targeted_effect_statistic_display_label <- function(effect_statistic_type) {
  dplyr::case_when(
    effect_statistic_type == "standardized_mean_rank_stat" ~ "Standardized mean rank statistic",
    effect_statistic_type == "NES" ~ "Normalized enrichment score (NES)",
    TRUE ~ "Method-specific enrichment effect statistic"
  )
}

targeted_effect_statistic_panel_label <- function(effect_statistic_type) {
  types <- sort(unique(as.character(effect_statistic_type)))
  types <- types[!is.na(types) & nzchar(types)]
  if (length(types) == 1L) return(targeted_effect_statistic_display_label(types))
  "Method-specific enrichment effect statistic"
}

targeted_collapse_provenance <- function(x) {
  x <- sort(unique(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) paste(x, collapse = ";") else NA_character_
}

targeted_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  !is.na(x) & toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
}

targeted_nonblank <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))

build_canonical_empirical_roi_term2gene <- function(
    markers,
    expected_contract = TARGETED_EMPIRICAL_MARKER_CONTRACT) {
  required <- c(
    "marker_set", "marker_evidence_type", "inferential_test_valid",
    "claim_allowed", "FDR_pass", "GeneSymbol", "marker_contract_version"
  )
  missing <- setdiff(required, names(markers))
  if (length(missing)) {
    stop("Canonical empirical ROI marker table is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  contracts <- unique(as.character(markers$marker_contract_version))
  contracts <- contracts[targeted_nonblank(contracts)]
  if (!length(contracts) || !all(contracts == expected_contract)) {
    stop(
      "Canonical empirical ROI marker contract mismatch; expected only '",
      expected_contract, "', observed: ", paste(contracts, collapse = ", "),
      call. = FALSE
    )
  }

  enriched <- markers |>
    dplyr::filter(
      .data$marker_set == "empirical_microglia_roi_enriched",
      .data$marker_evidence_type == "inferential",
      targeted_flag(.data$inferential_test_valid),
      targeted_flag(.data$claim_allowed),
      targeted_flag(.data$FDR_pass),
      targeted_nonblank(.data$GeneSymbol)
    ) |>
    dplyr::transmute(
      term = "empirical_microglia_enriched",
      gene = toupper(trimws(as.character(.data$GeneSymbol))),
      signature_source = "canonical_empirical_roi_v2",
      signature_membership_contract = expected_contract,
      signature_membership_evidence = "animal_paired_compartment_marker_inferential_claim_safe",
      signature_claim_basis_eligible = TRUE,
      signature_identifier_type = "claim_safe_gene_symbol"
    ) |>
    dplyr::filter(targeted_nonblank(.data$gene)) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$gene)

  shared <- tibble::tibble()
  if ("mapped_gene_symbol" %in% names(markers)) {
    shared <- markers |>
      dplyr::filter(
        .data$marker_set == "empirical_microglia_neuropil_shared",
        .data$marker_evidence_type == "descriptive_only",
        !targeted_flag(.data$claim_allowed),
        targeted_nonblank(.data$mapped_gene_symbol)
      ) |>
      dplyr::transmute(
        term = "empirical_neuropil_shared",
        gene = toupper(trimws(as.character(.data$mapped_gene_symbol))),
        signature_source = "canonical_empirical_roi_v2_descriptive_shared_abundance",
        signature_membership_contract = expected_contract,
        signature_membership_evidence = "descriptive_shared_abundance_not_equivalence",
        signature_claim_basis_eligible = FALSE,
        signature_identifier_type = "descriptive_gene_symbol"
      ) |>
      dplyr::filter(targeted_nonblank(.data$gene)) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$gene)
  }

  term2gene <- dplyr::bind_rows(enriched, shared) |>
    dplyr::distinct(.data$term, .data$gene, .keep_all = TRUE) |>
    dplyr::arrange(.data$term, .data$gene)

  list(
    term2gene = term2gene,
    microglia = enriched,
    shared = shared,
    diagnostics = tibble::tibble(
      check = c(
        "canonical_empirical_roi_contract",
        "inferential_microglia_enriched_genes",
        "descriptive_shared_abundance_genes",
        "stress_rank_derived_inferential_genes"
      ),
      status = c("PASS", ifelse(nrow(enriched) > 0, "PASS", "WARN"),
                 ifelse(nrow(shared) > 0, "PASS", "WARN"), "PASS"),
      detail = c(expected_contract, as.character(nrow(enriched)),
                 as.character(nrow(shared)), "0")
    )
  )
}

read_canonical_empirical_roi_term2gene <- function(path, expected_contract = TARGETED_EMPIRICAL_MARKER_CONTRACT) {
  if (!file.exists(path)) {
    stop("Canonical empirical ROI marker table not found: ", path, call. = FALSE)
  }
  markers <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  out <- build_canonical_empirical_roi_term2gene(markers, expected_contract)
  out$input_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  out
}

targeted_finite_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

targeted_finite_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

targeted_finite_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

targeted_same_direction <- function(x, reference) {
  is.finite(x) & is.finite(reference) & sign(x) == sign(reference)
}

targeted_layer_results <- function(x) {
  if (!nrow(x)) return(NA_character_)
  x <- x |>
    dplyr::arrange(.data$comparison, .data$layer, .data$effect_statistic, .data$padj)
  layer <- ifelse(targeted_nonblank(x$layer), as.character(x$layer), "NA")
  effect <- ifelse(is.finite(x$effect_statistic), formatC(x$effect_statistic, digits = 8, format = "fg"), "NA")
  padj <- ifelse(is.finite(x$padj), formatC(x$padj, digits = 8, format = "fg"), "NA")
  paste0(
    x$comparison,
    "[layer=", layer,
    ";effect_statistic=", effect,
    ";effect_statistic_type=", x$effect_statistic_type,
    ";enrichment_method=", x$enrichment_method,
    ";padj=", padj,
    ";fdr_scope=", x$fdr_scope,
    "]",
    collapse = "|"
  )
}

targeted_reference_summary <- function(candidates, micro_effect_statistic) {
  sig <- is.finite(candidates$padj) & candidates$padj < 0.05
  same <- sig & targeted_same_direction(candidates$effect_statistic, micro_effect_statistic)
  opposite <- sig & is.finite(candidates$effect_statistic) & is.finite(micro_effect_statistic) &
    sign(candidates$effect_statistic) != sign(micro_effect_statistic)
  tibble::tibble(
    matched_neuropil_n_layers = dplyr::n_distinct(as.character(candidates$comparison)),
    matched_neuropil_n_significant_layers = sum(sig),
    matched_neuropil_any_significant = any(sig),
    matched_neuropil_any_sig_same_direction = any(same),
    matched_neuropil_any_sig_opposite_direction = any(opposite),
    matched_neuropil_median_effect_statistic = targeted_finite_median(candidates$effect_statistic),
    matched_neuropil_min_effect_statistic = targeted_finite_min(candidates$effect_statistic),
    matched_neuropil_max_effect_statistic = targeted_finite_max(candidates$effect_statistic),
    matched_neuropil_effect_statistic_types = targeted_collapse_provenance(candidates$effect_statistic_type),
    matched_neuropil_enrichment_methods = targeted_collapse_provenance(candidates$enrichment_method),
    matched_neuropil_fdr_scopes = targeted_collapse_provenance(candidates$fdr_scope),
    matched_neuropil_median_NES = targeted_finite_median(candidates$effect_statistic),
    matched_neuropil_min_NES = targeted_finite_min(candidates$effect_statistic),
    matched_neuropil_max_NES = targeted_finite_max(candidates$effect_statistic),
    matched_neuropil_min_padj = targeted_finite_min(candidates$padj),
    matched_neuropil_route_keys = paste(sort(unique(paste(
      candidates$left_region, candidates$right_region,
      candidates$left_condition, candidates$right_condition,
      sep = "|"
    ))), collapse = ";"),
    matched_neuropil_layer_results = targeted_layer_results(candidates),
    neuropil_reference_comparison = paste(sort(unique(as.character(candidates$comparison))), collapse = ";")
  )
}

targeted_empty_reference_summary <- function() {
  tibble::tibble(
    matched_neuropil_n_layers = 0L,
    matched_neuropil_n_significant_layers = 0L,
    matched_neuropil_any_significant = FALSE,
    matched_neuropil_any_sig_same_direction = FALSE,
    matched_neuropil_any_sig_opposite_direction = FALSE,
    matched_neuropil_median_effect_statistic = NA_real_,
    matched_neuropil_min_effect_statistic = NA_real_,
    matched_neuropil_max_effect_statistic = NA_real_,
    matched_neuropil_effect_statistic_types = NA_character_,
    matched_neuropil_enrichment_methods = NA_character_,
    matched_neuropil_fdr_scopes = NA_character_,
    matched_neuropil_median_NES = NA_real_,
    matched_neuropil_min_NES = NA_real_,
    matched_neuropil_max_NES = NA_real_,
    matched_neuropil_min_padj = NA_real_,
    matched_neuropil_route_keys = NA_character_,
    matched_neuropil_layer_results = NA_character_,
    neuropil_reference_comparison = NA_character_
  )
}

targeted_attach_neuropil_reference <- function(micro, neuropil, reference_dataset = "neuron_neuropil") {
  if (!nrow(micro)) return(micro)
  route_cols <- c("signature", "left_region", "right_region", "left_condition", "right_condition")
  effect_cols <- c("effect_statistic", "effect_statistic_type", "enrichment_method", "fdr_scope")
  required_micro <- c(route_cols, "NES", effect_cols)
  required_neuropil <- c(route_cols, "NES", effect_cols, "padj", "comparison", "layer")
  missing_micro <- setdiff(required_micro, names(micro))
  missing_neuropil <- setdiff(required_neuropil, names(neuropil))
  if (length(missing_micro)) stop("Microglia enrichment is missing route column(s): ", paste(missing_micro, collapse = ", "), call. = FALSE)
  if (nrow(neuropil) && length(missing_neuropil)) stop("Neuropil enrichment is missing route column(s): ", paste(missing_neuropil, collapse = ", "), call. = FALSE)

  new_cols <- c(
    "reference_dataset", "reference_match_type", "neuropil_reference_effect_statistic",
    "neuropil_reference_effect_statistic_type", "neuropil_reference_enrichment_method",
    "neuropil_reference_fdr_scope", "neuropil_reference_NES",
    "neuropil_reference_padj", "neuropil_reference_comparison",
    "matched_neuropil_n_layers", "matched_neuropil_n_significant_layers",
    "matched_neuropil_any_significant", "matched_neuropil_any_sig_same_direction",
    "matched_neuropil_any_sig_opposite_direction", "matched_neuropil_median_effect_statistic",
    "matched_neuropil_min_effect_statistic", "matched_neuropil_max_effect_statistic",
    "matched_neuropil_effect_statistic_types", "matched_neuropil_enrichment_methods",
    "matched_neuropil_fdr_scopes", "matched_neuropil_median_NES",
    "matched_neuropil_min_NES", "matched_neuropil_max_NES", "matched_neuropil_min_padj",
    "matched_neuropil_route_keys", "matched_neuropil_layer_results", "global_neuropil_effect_statistic",
    "global_neuropil_effect_statistic_type", "global_neuropil_enrichment_method",
    "global_neuropil_fdr_scope", "global_neuropil_NES", "global_neuropil_padj",
    "global_neuropil_min_padj", "global_neuropil_comparison", "global_reference_match_type",
    "old_region_partial_match_would_have_been_used", "old_global_fallback_would_have_been_used",
    "matched_reference_uses_global_fallback",
    "neuropil_reference_summary_contract"
  )
  micro <- dplyr::select(micro, -dplyr::any_of(new_cols))

  summaries <- lapply(seq_len(nrow(micro)), function(i) {
    row <- micro[i, , drop = FALSE]
    exact <- if (nrow(neuropil)) neuropil else neuropil[FALSE, , drop = FALSE]
    for (col in route_cols) {
      value <- as.character(row[[col]][[1]])
      exact <- exact[!is.na(exact[[col]]) & !is.na(value) & as.character(exact[[col]]) == value, , drop = FALSE]
    }

    exact_summary <- if (nrow(exact)) targeted_reference_summary(exact, row$effect_statistic[[1]]) else targeted_empty_reference_summary()

    global <- if (nrow(neuropil)) {
      neuropil |>
        dplyr::filter(
          .data$signature == as.character(row$signature[[1]]),
          .data$left_condition == as.character(row$left_condition[[1]]),
          .data$right_condition == as.character(row$right_condition[[1]])
        )
    } else neuropil
    global_type <- if (nrow(global)) "same_formal_contrast_all_regions_descriptive" else "missing_global_context"
    if (!nrow(global) && nrow(neuropil)) {
      global <- neuropil |> dplyr::filter(.data$signature == as.character(row$signature[[1]]))
      if (nrow(global)) global_type <- "signature_global_descriptive"
    }

    old_region <- neuropil[FALSE, , drop = FALSE]
    old_condition <- neuropil[FALSE, , drop = FALSE]
    old_signature <- neuropil[FALSE, , drop = FALSE]
    if (nrow(neuropil)) {
      old_signature <- neuropil |>
        dplyr::filter(.data$signature == as.character(row$signature[[1]]))
      if (all(c("condition_pair", "region_level_unit") %in% names(neuropil)) &&
          all(c("condition_pair", "region_level_unit") %in% names(row))) {
        old_region <- old_signature |>
          dplyr::filter(
            .data$condition_pair == as.character(row$condition_pair[[1]]),
            .data$region_level_unit == as.character(row$region_level_unit[[1]])
          )
        old_condition <- old_signature |>
          dplyr::filter(.data$condition_pair == as.character(row$condition_pair[[1]]))
      } else {
        old_condition <- global
      }
    }
    old_global_used <- !nrow(old_region) && (nrow(old_condition) > 0 || nrow(old_signature) > 0)
    old_partial_region_used <- nrow(old_region) > 0 && !nrow(exact)

    dplyr::bind_cols(
      exact_summary,
      tibble::tibble(
        reference_dataset = reference_dataset,
        reference_match_type = if (nrow(exact)) "exact_same_region_formal_contrast" else "missing_exact_neuropil_reference",
        neuropil_reference_effect_statistic = if (nrow(exact)) exact_summary$matched_neuropil_median_effect_statistic else NA_real_,
        neuropil_reference_effect_statistic_type = if (nrow(exact)) exact_summary$matched_neuropil_effect_statistic_types else NA_character_,
        neuropil_reference_enrichment_method = if (nrow(exact)) exact_summary$matched_neuropil_enrichment_methods else NA_character_,
        neuropil_reference_fdr_scope = if (nrow(exact)) exact_summary$matched_neuropil_fdr_scopes else NA_character_,
        neuropil_reference_NES = if (nrow(exact)) exact_summary$matched_neuropil_median_NES else NA_real_,
        neuropil_reference_padj = NA_real_,
        global_neuropil_effect_statistic = if (nrow(global)) targeted_finite_median(global$effect_statistic) else NA_real_,
        global_neuropil_effect_statistic_type = if (nrow(global)) targeted_collapse_provenance(global$effect_statistic_type) else NA_character_,
        global_neuropil_enrichment_method = if (nrow(global)) targeted_collapse_provenance(global$enrichment_method) else NA_character_,
        global_neuropil_fdr_scope = if (nrow(global)) targeted_collapse_provenance(global$fdr_scope) else NA_character_,
        global_neuropil_NES = if (nrow(global)) targeted_finite_median(global$effect_statistic) else NA_real_,
        global_neuropil_padj = NA_real_,
        global_neuropil_min_padj = if (nrow(global)) targeted_finite_min(global$padj) else NA_real_,
        global_neuropil_comparison = if (nrow(global)) paste(sort(unique(as.character(global$comparison))), collapse = ";") else NA_character_,
        global_reference_match_type = global_type,
        old_region_partial_match_would_have_been_used = old_partial_region_used,
        old_global_fallback_would_have_been_used = old_global_used,
        matched_reference_uses_global_fallback = FALSE,
        neuropil_reference_summary_contract = "exact_route_all_layers_median_effect_statistic_v3"
      )
    )
  })

  dplyr::bind_cols(micro, dplyr::bind_rows(summaries))
}

targeted_claim_eligible_classes <- function() {
  c("microglia_enriched_empirical", "microglia_enriched_reference_supported")
}

targeted_classify_signature_evidence <- function(
    microglia_significant,
    matched_sig_same_direction,
    matched_sig_opposite_direction,
    empirical_support,
    reference_support,
    curated_support) {
  dplyr::case_when(
    targeted_flag(microglia_significant) & targeted_flag(matched_sig_same_direction) ~ "neuropil_shared",
    targeted_flag(microglia_significant) & targeted_flag(matched_sig_opposite_direction) ~ "mixed_microenvironment",
    targeted_flag(microglia_significant) & targeted_flag(empirical_support) ~ "microglia_enriched_empirical",
    targeted_flag(microglia_significant) & targeted_flag(reference_support) ~ "microglia_enriched_reference_supported",
    targeted_flag(microglia_significant) & targeted_flag(curated_support) ~ "curated_microglia_program",
    TRUE ~ "ambiguous"
  )
}

targeted_add_claim_ready <- function(df) {
  if (!nrow(df)) return(df)
  required <- c(
    "padj", "contrast_class", "reference_celltype_support",
    "microglia_signature_class", "signature_source",
    "signature_membership_contract", "signature_claim_basis_eligible",
    "reference_match_type", "matched_reference_uses_global_fallback"
  )
  missing <- setdiff(required, names(df))
  if (length(missing)) stop("Claim-ready input is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)

  significant <- is.finite(df$padj) & df$padj < 0.05
  interpretable <- df$contrast_class %in% c("within_region_condition", "cross_region_same_condition")
  reference_supported <- !is.na(df$reference_celltype_support) & df$reference_celltype_support == "microglia_supported"
  class_eligible <- df$microglia_signature_class %in% targeted_claim_eligible_classes()
  membership_eligible <- targeted_flag(df$signature_claim_basis_eligible)
  empirical_rows <- !is.na(df$microglia_signature_class) &
    df$microglia_signature_class == "microglia_enriched_empirical"
  empirical_independent <- !empirical_rows |
    (!is.na(df$signature_source) & df$signature_source == "canonical_empirical_roi_v2" &
       !is.na(df$signature_membership_contract) &
       df$signature_membership_contract == TARGETED_EMPIRICAL_MARKER_CONTRACT)
  no_global_fallback <- !targeted_flag(df$matched_reference_uses_global_fallback) &
    !grepl("global", ifelse(is.na(df$reference_match_type), "", df$reference_match_type), ignore.case = TRUE)
  ready <- significant & interpretable & reference_supported & class_eligible &
    membership_eligible & empirical_independent & no_global_fallback

  failure <- vapply(seq_len(nrow(df)), function(i) {
    reasons <- character()
    if (!significant[[i]]) reasons <- c(reasons, "microglia_not_fdr_significant")
    if (!interpretable[[i]]) reasons <- c(reasons, "contrast_not_interpretable")
    if (!reference_supported[[i]]) reasons <- c(reasons, "reference_celltype_not_microglia_supported")
    if (!class_eligible[[i]]) reasons <- c(reasons, "signature_class_not_claim_eligible")
    if (!membership_eligible[[i]]) reasons <- c(reasons, "signature_membership_not_claim_eligible")
    if (!empirical_independent[[i]]) reasons <- c(reasons, "empirical_membership_not_canonical_roi_v2")
    if (!no_global_fallback[[i]]) reasons <- c(reasons, "global_neuropil_fallback_used")
    if (length(reasons)) paste(reasons, collapse = ";") else "claim_ready"
  }, character(1))

  dplyr::mutate(
    df,
    claim_ready = ready,
    claim_ready_failure_reason = failure,
    claim_ready_contract_version = TARGETED_CLAIM_READY_CONTRACT
  )
}
