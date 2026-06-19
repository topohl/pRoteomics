# Pure helpers for deterministic WGCNA label candidate scoring and selection.
# This file performs no I/O and does not alter WGCNA statistics or assignments.

wgcna_label_or <- function(x, y) if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x

wgcna_label_col <- function(df, name, default = NA_character_) {
  if (name %in% names(df)) return(df[[name]])
  rep(default, nrow(df))
}

wgcna_label_clean <- function(x) {
  x <- trimws(gsub("\\s+", " ", as.character(x)))
  x[x %in% c("", "NA", "NaN", "Unlabelled", "Unassigned")] <- NA_character_
  x
}

wgcna_biological_label <- function(x) {
  x <- wgcna_label_clean(x)
  x <- sub("^(microglia[- ]enriched )?shared ROI:\\s*", "", x, ignore.case = TRUE)
  x <- sub("^microglia ROI:\\s*", "", x, ignore.case = TRUE)
  x <- sub("^SM[0-9]+\\s*(?:·|:|-)\\s*", "", x, ignore.case = TRUE, perl = TRUE)
  x <- sub("^(?:singleton|dominant|mostly):?\\s*", "", x, ignore.case = TRUE, perl = TRUE)
  x <- sub("^mixed:\\s*", "", x, ignore.case = TRUE)
  wgcna_label_clean(x)
}

wgcna_candidate_rows <- function(df, level, id_col, source_columns, dataset = NULL) {
  if (is.null(df) || !nrow(df)) return(tibble::tibble())
  ds <- wgcna_label_clean(wgcna_label_col(df, "dataset", dataset))
  ds[is.na(ds)] <- as.character(wgcna_label_or(dataset, NA_character_))
  ids <- wgcna_label_clean(wgcna_label_col(df, id_col))
  hub_text <- dplyr::coalesce(
    wgcna_label_clean(wgcna_label_col(df, "raw_hub_proteins")),
    wgcna_label_clean(wgcna_label_col(df, "evidence_hubs")),
    wgcna_label_clean(wgcna_label_col(df, "MemberEvidenceHubs"))
  )
  marker_text <- dplyr::coalesce(
    wgcna_label_clean(wgcna_label_col(df, "raw_marker_or_signature_label")),
    wgcna_label_clean(wgcna_label_col(df, "evidence_marker_panels")),
    wgcna_label_clean(wgcna_label_col(df, "MemberEvidenceMarkerPanels")),
    wgcna_label_clean(wgcna_label_col(df, "marker_panels_supporting"))
  )
  mismatch_text <- tolower(paste(
    wgcna_label_col(df, "GO_label_relevance_flag"),
    wgcna_label_col(df, "label_warning"),
    wgcna_label_col(df, "GO_label_relevance_rationale")
  ))
  mismatch <- grepl("mismatch|artefact|non[- ]neural|incongruent|unsafe", mismatch_text)
  conflict_text <- tolower(paste(
    wgcna_label_col(df, "label_warning"),
    wgcna_label_col(df, "annotation_downgrade_reason"),
    wgcna_label_col(df, "label_downgrade_reason")
  ))
  conflict <- grepl("conflict|unstable|mixed|unresolved|mismatch", conflict_text)

  out <- lapply(names(source_columns), function(source) {
    spec <- source_columns[[source]]
    label <- wgcna_biological_label(wgcna_label_col(df, spec$column))
    tibble::tibble(
      dataset = ds,
      level = level,
      entity_id = ids,
      candidate_label = label,
      candidate_source = source,
      evidence_strength = as.numeric(spec$strength),
      hub_support = as.numeric(!is.na(hub_text) & nzchar(hub_text)),
      marker_context_support = as.numeric(!is.na(marker_text) & nzchar(marker_text)),
      genericity_penalty = as.numeric(is.na(label) | grepl("^(mixed|unresolved|unknown|unlabelled)|low-specificity|^SM[0-9]+$|^WGCNA_", label, ignore.case = TRUE)),
      ontology_mismatch_flag = mismatch,
      conflict_penalty = as.numeric(conflict)
    )
  })
  dplyr::bind_rows(out) |>
    dplyr::filter(!is.na(.data$dataset), !is.na(.data$entity_id), !is.na(.data$candidate_label), nzchar(.data$candidate_label)) |>
    dplyr::distinct(.data$dataset, .data$level, .data$entity_id, .data$candidate_source, .data$candidate_label, .keep_all = TRUE)
}

wgcna_make_label_candidates <- function(module_rows, supermodule_rows, dataset = NULL) {
  module_sources <- list(
    cleaned_biological_label = list(column = "cleaned_biological_label", strength = 5),
    module_specific_label = list(column = "module_specific_label", strength = 4.8),
    module_biological_label = list(column = "module_biological_label", strength = 4.6),
    safe_display_label = list(column = "safe_display_label", strength = 4.0),
    raw_top_GO_label = list(column = "raw_top_GO_label", strength = 3.4),
    ModuleLabel_Final = list(column = "ModuleLabel_Final", strength = 2.6),
    raw_module_label = list(column = "raw_module_label", strength = 1.5)
  )
  super_sources <- list(
    DominantMemberTheme = list(column = "DominantMemberTheme", strength = 5.2),
    cleaned_biological_label = list(column = "cleaned_biological_label", strength = 5.0),
    Supermodule_DataDrivenLabel = list(column = "Supermodule_DataDrivenLabel", strength = 4.6),
    Supermodule_CompositionLabel = list(column = "Supermodule_CompositionLabel", strength = 4.2),
    Supermodule_ConservativeLabel = list(column = "Supermodule_ConservativeLabel", strength = 3.8),
    raw_top_GO_label = list(column = "raw_top_GO_label", strength = 3.2),
    Supermodule_FinalLabel = list(column = "Supermodule_FinalLabel", strength = 2.4)
  )
  dplyr::bind_rows(
    wgcna_candidate_rows(module_rows, "module", "module_id", module_sources, dataset),
    wgcna_candidate_rows(supermodule_rows, "supermodule", "supermodule_id", super_sources, dataset)
  )
}

wgcna_score_label_candidates <- function(candidates) {
  if (is.null(candidates) || !nrow(candidates)) return(candidates)
  candidates |>
    dplyr::mutate(
      final_label_score = .data$evidence_strength + 0.75 * .data$hub_support +
        0.35 * .data$marker_context_support - 1.5 * .data$genericity_penalty -
        1.25 * as.numeric(.data$ontology_mismatch_flag) - 0.75 * .data$conflict_penalty
    )
}

wgcna_select_final_labels <- function(candidates) {
  if (is.null(candidates) || !nrow(candidates)) return(candidates)
  candidates |>
    dplyr::group_by(.data$dataset, .data$level, .data$entity_id) |>
    dplyr::arrange(dplyr::desc(.data$final_label_score), dplyr::desc(.data$evidence_strength), .data$candidate_source, .by_group = TRUE) |>
    dplyr::mutate(
      selected_label = dplyr::row_number() == 1L,
      rejection_reason = dplyr::case_when(
        .data$selected_label ~ "selected_highest_score",
        .data$ontology_mismatch_flag ~ "ontology_mismatch_penalty",
        .data$genericity_penalty > 0 ~ "generic_or_unresolved_label",
        .data$conflict_penalty > 0 ~ "conflicting_or_unstable_evidence",
        TRUE ~ "lower_final_label_score"
      )
    ) |>
    dplyr::ungroup()
}

wgcna_parse_named_fractions <- function(x) {
  parts <- trimws(unlist(strsplit(as.character(wgcna_label_or(x, "")), "\\s*;\\s*", perl = TRUE), use.names = FALSE))
  parts <- parts[nzchar(parts) & grepl("=", parts, fixed = TRUE)]
  vals <- suppressWarnings(as.numeric(sub("^.*=", "", parts)))
  names(vals) <- trimws(sub("=.*$", "", parts))
  vals <- vals[is.finite(vals) & !is.na(names(vals)) & nzchar(names(vals))]
  sort(vals, decreasing = TRUE)
}

wgcna_supermodule_plot_label <- function(id, n_members, dominant_label, dominant_fraction, fractions) {
  id <- as.character(id)
  dominant_label <- wgcna_biological_label(dominant_label)
  n_members <- suppressWarnings(as.integer(n_members))
  dominant_fraction <- suppressWarnings(as.numeric(dominant_fraction))
  vals <- wgcna_parse_named_fractions(fractions)
  eligible <- vals[vals >= 0.25 & !grepl("mixed|unresolved|low-specificity", names(vals), ignore.case = TRUE)]
  labels <- unique(names(eligible))
  if (is.finite(n_members) && n_members == 1L) {
    return(paste0(id, " · singleton: ", wgcna_label_or(dominant_label, "mixed / unresolved")))
  }
  if (is.finite(dominant_fraction) && dominant_fraction >= 0.67 && !is.na(dominant_label)) {
    return(paste0(id, " · dominant: ", dominant_label))
  }
  if (length(labels) >= 2L) {
    return(paste0(id, " · mixed: ", paste(utils::head(labels, 3L), collapse = "; ")))
  }
  paste0(id, " · mixed / unresolved")
}

wgcna_build_final_label_lookup <- function(selected_candidates, module_rows, supermodule_rows, dataset = NULL) {
  selected <- selected_candidates |>
    dplyr::filter(.data$selected_label) |>
    dplyr::select(
      "dataset", "level", "entity_id", best_data_driven_label = "candidate_label",
      label_score = "final_label_score", selected_source = "candidate_source",
      "hub_support", "marker_context_support", "ontology_mismatch_flag"
    )

  module_meta <- module_rows |>
    dplyr::transmute(
      dataset = dplyr::coalesce(wgcna_label_clean(wgcna_label_col(module_rows, "dataset")), dataset),
      level = "module",
      entity_id = wgcna_label_clean(wgcna_label_col(module_rows, "module_id")),
      parent_entity_id = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(module_rows, "supermodule_id")),
        wgcna_label_clean(wgcna_label_col(module_rows, "supermodule_id_for_module")),
        wgcna_label_clean(wgcna_label_col(module_rows, "module_supermodule_id"))
      ),
      n_member_modules = 1L,
      raw_top_GO_label = wgcna_label_clean(wgcna_label_col(module_rows, "raw_top_GO_label")),
      parent_program = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(module_rows, "module_program_primary")),
        wgcna_label_clean(wgcna_label_col(module_rows, "Macroprogram_Display"))
      ),
      microenvironment_context = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(module_rows, "microenvironment_caution_label")),
        wgcna_label_clean(wgcna_label_col(module_rows, "microenvironment_label")),
        wgcna_label_clean(wgcna_label_col(module_rows, "microenvironment_class"))
      ),
      annotation_confidence = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(module_rows, "label_confidence")),
        wgcna_label_clean(wgcna_label_col(module_rows, "annotation_confidence")),
        wgcna_label_clean(wgcna_label_col(module_rows, "cleaned_biological_label_confidence"))
      ),
      label_basis_upstream = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(module_rows, "label_basis")),
        wgcna_label_clean(wgcna_label_col(module_rows, "cleaned_biological_label_source"))
      ),
      unsafe_interpretation = wgcna_label_clean(wgcna_label_col(module_rows, "unsafe_interpretation")),
      independent_go = !is.na(wgcna_label_clean(wgcna_label_col(module_rows, "raw_top_GO_label")))
    ) |>
    dplyr::filter(!is.na(.data$entity_id)) |>
    dplyr::distinct(.data$dataset, .data$level, .data$entity_id, .keep_all = TRUE)

  super_meta <- supermodule_rows |>
    dplyr::transmute(
      dataset = dplyr::coalesce(wgcna_label_clean(wgcna_label_col(supermodule_rows, "dataset")), dataset),
      level = "supermodule",
      entity_id = wgcna_label_clean(wgcna_label_col(supermodule_rows, "supermodule_id")),
      parent_entity_id = NA_character_,
      n_member_modules = suppressWarnings(as.integer(dplyr::coalesce(
        wgcna_label_col(supermodule_rows, "n_member_modules", NA_integer_),
        wgcna_label_col(supermodule_rows, "DataDrivenClusterSize", NA_integer_)
      ))),
      raw_top_GO_label = wgcna_label_clean(wgcna_label_col(supermodule_rows, "raw_top_GO_label")),
      parent_program = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "module_program_primary")),
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "Macroprogram_Display"))
      ),
      microenvironment_context = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "microenvironment_caution_label")),
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "dominant_microenvironment_class"))
      ),
      annotation_confidence = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "label_confidence")),
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "annotation_confidence")),
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "Supermodule_LabelConfidence"))
      ),
      label_basis_upstream = dplyr::coalesce(
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "label_basis")),
        wgcna_label_clean(wgcna_label_col(supermodule_rows, "Supermodule_LabelSource"))
      ),
      unsafe_interpretation = wgcna_label_clean(wgcna_label_col(supermodule_rows, "unsafe_interpretation")),
      independent_go = !is.na(wgcna_label_clean(wgcna_label_col(supermodule_rows, "raw_top_GO_label"))),
      dominant_label = wgcna_label_clean(wgcna_label_col(supermodule_rows, "DominantMemberTheme")),
      dominant_fraction = suppressWarnings(as.numeric(wgcna_label_col(supermodule_rows, "DominantMemberThemeFraction", NA_real_))),
      member_fractions = wgcna_label_clean(wgcna_label_col(supermodule_rows, "MemberThemeFractions"))
    ) |>
    dplyr::filter(!is.na(.data$entity_id)) |>
    dplyr::distinct(.data$dataset, .data$level, .data$entity_id, .keep_all = TRUE)

  meta <- dplyr::bind_rows(module_meta, super_meta) |>
    dplyr::left_join(selected, by = c("dataset", "level", "entity_id")) |>
    dplyr::mutate(
      best_data_driven_label = dplyr::coalesce(.data$best_data_driven_label, "mixed / unresolved"),
      final_plot_label = dplyr::if_else(
        .data$level == "supermodule",
        mapply(wgcna_supermodule_plot_label, .data$entity_id, .data$n_member_modules,
               dplyr::coalesce(.data$dominant_label, .data$best_data_driven_label),
               .data$dominant_fraction, .data$member_fractions, USE.NAMES = FALSE),
        paste0(.data$entity_id, " · ", .data$best_data_driven_label)
      ),
      label_confidence = dplyr::case_when(
        .data$ontology_mismatch_flag & !(.data$hub_support > 0 | .data$independent_go) ~ "low",
        .data$label_score >= 5.5 ~ "high",
        .data$label_score >= 3.5 ~ "moderate",
        TRUE ~ dplyr::coalesce(.data$annotation_confidence, "low")
      ),
      label_basis = paste0("selected_source=", .data$selected_source, "; upstream=", dplyr::coalesce(.data$label_basis_upstream, "not_available")),
      manual_review_required = .data$ontology_mismatch_flag & !(.data$hub_support > 0 | .data$independent_go) |
        .data$label_confidence == "low" | grepl("unresolved", .data$final_plot_label, ignore.case = TRUE),
      independent_evidence_recorded = .data$hub_support > 0 | .data$independent_go,
      label_rationale = paste0(
        "Candidate selected by deterministic evidence score; score=", sprintf("%.2f", .data$label_score),
        "; biological label kept separate from microenvironment context."
      ),
      unsafe_interpretation = dplyr::coalesce(
        .data$unsafe_interpretation,
        "Do not interpret a WGCNA label as causal, cell-pure, or cell-intrinsic evidence."
      )
    ) |>
    dplyr::select(
      "dataset", "level", "entity_id", "parent_entity_id", "n_member_modules",
      "raw_top_GO_label", "best_data_driven_label", "parent_program",
      "microenvironment_context", "final_plot_label", "label_confidence",
      "label_basis", "label_score", "manual_review_required",
      "ontology_mismatch_flag", "independent_evidence_recorded", "label_rationale", "unsafe_interpretation"
    ) |>
    dplyr::arrange(.data$dataset, .data$level, .data$entity_id)
  meta
}

wgcna_validate_label_lookup <- function(lookup) {
  required <- c(
    "dataset", "level", "entity_id", "parent_entity_id", "n_member_modules",
    "raw_top_GO_label", "best_data_driven_label", "parent_program",
    "microenvironment_context", "final_plot_label", "label_confidence",
    "label_basis", "label_score", "manual_review_required",
    "ontology_mismatch_flag", "label_rationale", "unsafe_interpretation"
  )
  missing <- setdiff(required, names(lookup))
  if (length(missing)) stop("WGCNA final label lookup missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  dup <- duplicated(lookup[c("dataset", "level", "entity_id")])
  if (any(dup)) stop("WGCNA final label lookup has duplicate dataset x level x entity_id rows.", call. = FALSE)
  if (any(is.na(lookup$final_plot_label) | !nzchar(trimws(lookup$final_plot_label)))) stop("WGCNA final label lookup contains empty final_plot_label values.", call. = FALSE)
  singleton_mostly <- lookup$level == "supermodule" & lookup$n_member_modules == 1L & grepl("mostly", lookup$final_plot_label, ignore.case = TRUE)
  if (any(singleton_mostly, na.rm = TRUE)) stop("Singleton supermodule labels must not contain 'mostly'.", call. = FALSE)
  independent <- if ("independent_evidence_recorded" %in% names(lookup)) lookup$independent_evidence_recorded %in% TRUE else FALSE
  mismatch_unreviewed <- lookup$ontology_mismatch_flag %in% TRUE & !independent & !(lookup$manual_review_required %in% TRUE)
  if (any(mismatch_unreviewed, na.rm = TRUE)) stop("Ontology-context mismatch labels require manual review unless independent evidence cleared the mismatch before lookup construction.", call. = FALSE)
  invisible(TRUE)
}
