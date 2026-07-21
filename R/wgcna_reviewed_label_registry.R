# Reviewed WGCNA label-registry and canonical lookup helpers.
# These helpers are deterministic and do not alter network or statistical state.

wgcna_reviewed_registry_columns <- function() {
  c(
    "dataset", "level", "entity_id", "reviewed_biological_label",
    "reviewed_short_label", "subcellular_context", "roi_context",
    "label_confidence", "manual_review_required",
    "aggregation_evidence_class", "structural_status", "rationale",
    "reviewer", "review_date", "adjudication_status"
  )
}

wgcna_canonical_lookup_columns <- function() {
  c(
    "dataset", "level", "entity_id", "parent_entity_id",
    "canonical_biological_label", "canonical_short_label",
    "canonical_plot_label", "biological_label_confidence",
    "manual_review_required", "structural_status",
    "aggregation_evidence_class", "structural_coherence_class",
    "subcellular_context", "roi_context", "n_member_modules",
    "member_ModuleIDs", "member_modules_fingerprint",
    "annotation_evidence_fingerprint", "label_contract_version",
    "label_source", "label_rationale", "unsafe_interpretation",
    "final_plot_label"
  )
}

wgcna_sha256_text <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for WGCNA label fingerprints.", call. = FALSE)
  }
  x <- as.character(x)
  vapply(x, digest::digest, character(1), algo = "sha256", serialize = FALSE)
}

wgcna_clean_registry_text <- function(x) {
  x <- trimws(gsub("\\s+", " ", as.character(x)))
  x[is.na(x) | !nzchar(x) | toupper(x) %in% c("NA", "NAN", "NULL")] <- NA_character_
  x
}

wgcna_normalize_current_member_map <- function(member_map, dataset) {
  if (is.null(member_map) || !is.data.frame(member_map) || !nrow(member_map)) {
    stop("Current WGCNA module-to-supermodule member map is missing or empty.", call. = FALSE)
  }
  module_hits <- intersect(c("ModuleID", "module_id", "module_key"), names(member_map))
  super_hits <- intersect(c("SupermoduleID", "supermodule_id", "module_supermodule_id"), names(member_map))
  if (!length(module_hits) || !length(super_hits)) {
    stop("Current member map requires stable ModuleID and SupermoduleID columns.", call. = FALSE)
  }
  module_col <- module_hits[[1]]
  super_col <- super_hits[[1]]
  out <- tibble::tibble(
    dataset = if ("dataset" %in% names(member_map)) as.character(member_map$dataset) else dataset,
    ModuleID = wgcna_clean_registry_text(member_map[[module_col]]),
    SupermoduleID = wgcna_clean_registry_text(member_map[[super_col]])
  ) |>
    dplyr::filter(.data$dataset == .env$dataset, !is.na(.data$ModuleID), !is.na(.data$SupermoduleID)) |>
    dplyr::distinct()
  if (any(!grepl("^WGCNA_m[0-9]{2}$", out$ModuleID))) {
    stop("Current member map contains stale or invalid ModuleID values.", call. = FALSE)
  }
  if (any(!grepl("^SM[0-9]{2}$", out$SupermoduleID))) {
    stop("Current member map contains stale or invalid SupermoduleID values.", call. = FALSE)
  }
  conflict <- out |>
    dplyr::count(.data$dataset, .data$ModuleID, name = "n_memberships") |>
    dplyr::filter(.data$n_memberships != 1L)
  if (nrow(conflict)) {
    stop("Each stable ModuleID must map to exactly one current SupermoduleID.", call. = FALSE)
  }
  if (anyDuplicated(out[c("dataset", "ModuleID")])) {
    stop("Current member map contains duplicated dataset + ModuleID keys.", call. = FALSE)
  }
  out |>
    dplyr::arrange(.data$dataset, .data$SupermoduleID, .data$ModuleID)
}

wgcna_validate_reviewed_registry <- function(registry, dataset, member_map,
                                             expected_n_modules = 13L,
                                             expected_n_supermodules = 9L) {
  if (!is.data.frame(registry)) stop("Reviewed WGCNA label registry must be a data frame.", call. = FALSE)
  expected_columns <- wgcna_reviewed_registry_columns()
  if (!identical(names(registry), expected_columns)) {
    stop(
      "Reviewed WGCNA label registry schema mismatch; missing: ",
      paste(setdiff(expected_columns, names(registry)), collapse = ", "),
      "; unexpected: ", paste(setdiff(names(registry), expected_columns), collapse = ", "), ".",
      call. = FALSE
    )
  }
  current_map <- wgcna_normalize_current_member_map(member_map, dataset)
  registry <- registry |>
    dplyr::mutate(
      dataset = as.character(.data$dataset),
      level = as.character(.data$level),
      entity_id = as.character(.data$entity_id),
      manual_review_required = as.logical(.data$manual_review_required)
    )
  if (any(registry$dataset != dataset)) stop("Reviewed registry contains rows for the wrong dataset.", call. = FALSE)
  if (anyDuplicated(registry[c("dataset", "level", "entity_id")])) {
    stop("Reviewed registry has duplicate dataset + level + entity_id keys.", call. = FALSE)
  }
  modules <- registry |> dplyr::filter(.data$level == "module")
  supers <- registry |> dplyr::filter(.data$level == "supermodule")
  current_modules <- sort(unique(current_map$ModuleID))
  current_supers <- sort(unique(current_map$SupermoduleID))
  if (nrow(modules) != expected_n_modules || length(current_modules) != expected_n_modules) {
    stop("Reviewed registry/current map must contain exactly ", expected_n_modules, " modules.", call. = FALSE)
  }
  if (nrow(supers) != expected_n_supermodules || length(current_supers) != expected_n_supermodules) {
    stop("Reviewed registry/current map must contain exactly ", expected_n_supermodules, " supermodules.", call. = FALSE)
  }
  if (!identical(sort(modules$entity_id), current_modules) || !identical(sort(supers$entity_id), current_supers)) {
    stop("Reviewed registry contains stale identities or omits current stable identities.", call. = FALSE)
  }
  if (any(registry$adjudication_status != "reviewed")) {
    stop("All reviewed registry rows must have adjudication_status = reviewed.", call. = FALSE)
  }
  required_text <- c("reviewed_biological_label", "reviewed_short_label", "subcellular_context", "roi_context", "label_confidence", "rationale")
  for (nm in required_text) {
    if (any(is.na(wgcna_clean_registry_text(registry[[nm]])))) stop("Reviewed registry contains empty ", nm, " values.", call. = FALSE)
  }
  directional <- grepl(
    "\\b(SUS|RES|CON|higher|lower|increased?|decreased?|up-?regulated|down-?regulated)\\b",
    registry$reviewed_biological_label,
    ignore.case = TRUE,
    perl = TRUE
  )
  if (any(directional)) stop("Biological labels must not include group-effect direction.", call. = FALSE)

  singleton_map <- current_map |>
    dplyr::group_by(.data$SupermoduleID) |>
    dplyr::summarise(n = dplyr::n(), member = dplyr::first(.data$ModuleID), .groups = "drop") |>
    dplyr::filter(.data$n == 1L)
  if (nrow(singleton_map) != 6L) stop("Current microglia map must contain exactly six singleton supermodules.", call. = FALSE)
  singleton_check <- singleton_map |>
    dplyr::left_join(
      modules |>
        dplyr::select(
          member = "entity_id", member_bio = "reviewed_biological_label",
          member_short = "reviewed_short_label", member_conf = "label_confidence",
          member_manual = "manual_review_required"
        ),
      by = "member", relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      supers |>
        dplyr::select(
          SupermoduleID = "entity_id", super_bio = "reviewed_biological_label",
          super_short = "reviewed_short_label", super_conf = "label_confidence",
          super_manual = "manual_review_required", "structural_status",
          "aggregation_evidence_class"
        ),
      by = "SupermoduleID", relationship = "many-to-one"
    )
  inherited_ok <- with(singleton_check,
    member_bio == super_bio & member_short == super_short & member_conf == super_conf &
      member_manual == super_manual & structural_status == "singleton" &
      aggregation_evidence_class == "not_applicable_singleton"
  )
  if (!isTRUE(all(inherited_ok))) stop("Singleton supermodules must inherit the exact current member label, confidence, and review state.", call. = FALSE)
  invisible(TRUE)
}

wgcna_read_reviewed_registry <- function(dataset, member_map,
                                         path = file.path("config", "wgcna_labels", paste0(dataset, ".csv"))) {
  if (!file.exists(path)) stop("Missing reviewed WGCNA label registry: ", path, call. = FALSE)
  registry <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  wgcna_validate_reviewed_registry(registry, dataset, member_map)
  registry
}

wgcna_annotation_fingerprints <- function(registry, label_candidates, module_rows, supermodule_rows) {
  evidence_columns <- c(
    "raw_GO_BP_terms", "raw_GO_CC_terms", "raw_GO_MF_terms", "raw_top_GO_label",
    "raw_hub_proteins", "evidence_BP", "evidence_CC", "evidence_MF", "evidence_hubs",
    "microenvironment_class", "microenvironment_label", "dominant_microenvironment_class",
    "MemberThemeCounts", "MemberThemeFractions", "Supermodule_CompositionLabel"
  )
  collapse_rows <- function(df, level, id_candidates) {
    if (is.null(df) || !nrow(df)) return(tibble::tibble(level = character(), entity_id = character(), evidence = character()))
    id_hits <- intersect(id_candidates, names(df))
    cols <- intersect(evidence_columns, names(df))
    if (!length(id_hits)) return(tibble::tibble(level = character(), entity_id = character(), evidence = character()))
    id_col <- id_hits[[1]]
    if (!length(cols)) cols <- id_col
    tibble::tibble(
      level = level,
      entity_id = as.character(df[[id_col]]),
      evidence = apply(df[cols], 1L, function(z) paste(names(z), wgcna_clean_registry_text(z), sep = "=", collapse = "|"))
    ) |>
      dplyr::group_by(.data$level, .data$entity_id) |>
      dplyr::summarise(evidence = paste(sort(unique(.data$evidence)), collapse = "||"), .groups = "drop")
  }
  annotation <- dplyr::bind_rows(
    collapse_rows(module_rows, "module", c("module_id", "ModuleID")),
    collapse_rows(supermodule_rows, "supermodule", c("supermodule_id", "SupermoduleID"))
  )
  candidates <- if (!is.null(label_candidates) && nrow(label_candidates)) {
    label_candidates |>
      dplyr::transmute(
        level = as.character(.data$level), entity_id = as.character(.data$entity_id),
        evidence = paste(.data$candidate_source, .data$candidate_label, .data$final_label_score, sep = "=")
      ) |>
      dplyr::group_by(.data$level, .data$entity_id) |>
      dplyr::summarise(candidate_evidence = paste(sort(unique(.data$evidence)), collapse = "||"), .groups = "drop")
  } else {
    tibble::tibble(level = character(), entity_id = character(), candidate_evidence = character())
  }
  registry |>
    dplyr::select("level", "entity_id") |>
    dplyr::left_join(annotation, by = c("level", "entity_id"), relationship = "many-to-one") |>
    dplyr::left_join(candidates, by = c("level", "entity_id"), relationship = "many-to-one") |>
    dplyr::mutate(
      fingerprint_text = paste(
        dplyr::coalesce(.data$evidence, "no_annotation_evidence"),
        dplyr::coalesce(.data$candidate_evidence, "no_automatic_candidates"),
        sep = "||"
      ),
      annotation_evidence_fingerprint = wgcna_sha256_text(.data$fingerprint_text)
    ) |>
    dplyr::select("level", "entity_id", "annotation_evidence_fingerprint")
}

wgcna_structural_coherence <- function(supermodule_rows, ids) {
  if (is.null(supermodule_rows) || !nrow(supermodule_rows)) {
    return(tibble::tibble(entity_id = ids, structural_coherence_class = "not_available"))
  }
  id_hits <- intersect(c("supermodule_id", "SupermoduleID"), names(supermodule_rows))
  if (!length(id_hits)) return(tibble::tibble(entity_id = ids, structural_coherence_class = "not_available"))
  id_col <- id_hits[[1]]
  n_hits <- intersect(c("n_member_modules", "DataDrivenClusterSize"), names(supermodule_rows))
  n_col <- if (length(n_hits)) n_hits[[1]] else NULL
  get_num <- function(candidates) {
    hit <- intersect(candidates, names(supermodule_rows))
    if (length(hit)) suppressWarnings(as.numeric(supermodule_rows[[hit[[1]]]])) else rep(NA_real_, nrow(supermodule_rows))
  }
  out <- tibble::tibble(
    entity_id = as.character(supermodule_rows[[id_col]]),
    n = if (!is.null(n_col)) suppressWarnings(as.integer(supermodule_rows[[n_col]])) else NA_integer_,
    min_cor = get_num(c("adjusted_signed_min_pairwise_eigengene_correlation", "signed_min_pairwise_eigengene_correlation")),
    pc1 = get_num("pc1_variance_explained"),
    stability = get_num("cut_height_stability_fraction_stable")
  ) |>
    dplyr::distinct(.data$entity_id, .keep_all = TRUE) |>
    dplyr::mutate(
      structural_coherence_class = dplyr::case_when(
        .data$n == 1L ~ "not_applicable_singleton",
        is.finite(.data$min_cor) & .data$min_cor >= 0.50 & is.finite(.data$pc1) & .data$pc1 >= 0.70 & is.finite(.data$stability) & .data$stability >= 0.50 ~ "coherent_stable",
        TRUE ~ "mixed_or_limited_coherence"
      )
    ) |>
    dplyr::select("entity_id", "structural_coherence_class")
  tibble::tibble(entity_id = ids) |>
    dplyr::left_join(out, by = "entity_id", relationship = "many-to-one") |>
    dplyr::mutate(structural_coherence_class = dplyr::coalesce(.data$structural_coherence_class, "not_available"))
}

wgcna_build_reviewed_canonical_lookup <- function(registry, member_map, label_candidates,
                                                   module_rows, supermodule_rows,
                                                   dataset,
                                                   label_source = file.path("config", "wgcna_labels", paste0(dataset, ".csv"))) {
  wgcna_validate_reviewed_registry(registry, dataset, member_map)
  current_map <- wgcna_normalize_current_member_map(member_map, dataset)
  module_meta <- current_map |>
    dplyr::transmute(
      level = "module", entity_id = .data$ModuleID, parent_entity_id = .data$SupermoduleID,
      n_member_modules = 1L, member_ModuleIDs = .data$ModuleID,
      member_modules_fingerprint = wgcna_sha256_text(paste(.data$dataset, .data$ModuleID, .data$SupermoduleID, sep = "|")),
      structural_coherence_class = "not_applicable_module"
    )
  super_meta <- current_map |>
    dplyr::group_by(.data$dataset, .data$SupermoduleID) |>
    dplyr::summarise(
      level = "supermodule", entity_id = dplyr::first(.data$SupermoduleID),
      parent_entity_id = NA_character_, n_member_modules = dplyr::n(),
      member_ModuleIDs = paste(sort(.data$ModuleID), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      member_modules_fingerprint = wgcna_sha256_text(paste(.data$dataset, .data$entity_id, .data$member_ModuleIDs, sep = "|"))
    ) |>
    dplyr::left_join(wgcna_structural_coherence(supermodule_rows, unique(current_map$SupermoduleID)), by = "entity_id", relationship = "many-to-one") |>
    dplyr::mutate(
      structural_coherence_class = dplyr::if_else(
        .data$n_member_modules == 1L,
        "not_applicable_singleton",
        .data$structural_coherence_class
      )
    ) |>
    dplyr::select(-"SupermoduleID")
  fingerprints <- wgcna_annotation_fingerprints(registry, label_candidates, module_rows, supermodule_rows)
  meta <- dplyr::bind_rows(module_meta, super_meta) |>
    dplyr::select(-dplyr::any_of("dataset"))
  out <- registry |>
    dplyr::left_join(meta, by = c("level", "entity_id"), relationship = "one-to-one") |>
    dplyr::left_join(fingerprints, by = c("level", "entity_id"), relationship = "one-to-one") |>
    dplyr::mutate(
      parent_entity_id = as.character(.data$parent_entity_id),
      canonical_biological_label = as.character(.data$reviewed_biological_label),
      canonical_short_label = as.character(.data$reviewed_short_label),
      canonical_plot_label = paste0(.data$entity_id, "\n", .data$canonical_short_label),
      biological_label_confidence = as.character(.data$label_confidence),
      structural_coherence_class = as.character(.data$structural_coherence_class),
      label_contract_version = "microglia_wgcna_reviewed_labels_v1",
      label_source = label_source,
      label_rationale = as.character(.data$rationale),
      unsafe_interpretation = dplyr::case_when(
        .data$level == "module" ~ "Do not interpret an ROI-derived WGCNA module label as causal, cell-pure, or microglia-intrinsic evidence.",
        TRUE ~ "Do not interpret an ROI-derived WGCNA supermodule label as causal, cell-pure, or microglia-intrinsic evidence."
      ),
      final_plot_label = .data$canonical_plot_label
    ) |>
    dplyr::select(dplyr::all_of(wgcna_canonical_lookup_columns())) |>
    dplyr::arrange(factor(.data$level, levels = c("module", "supermodule")), .data$entity_id)
  wgcna_validate_canonical_lookup(out, dataset, current_map)
  out
}

wgcna_validate_canonical_lookup <- function(lookup, dataset = NULL, member_map = NULL) {
  expected <- wgcna_canonical_lookup_columns()
  if (!is.data.frame(lookup) || !identical(names(lookup), expected)) {
    stop("Canonical WGCNA label lookup does not match the reviewed-label contract.", call. = FALSE)
  }
  if (anyDuplicated(lookup[c("dataset", "level", "entity_id")])) stop("Canonical lookup has duplicate entity keys.", call. = FALSE)
  if (nrow(lookup) != 22L || sum(lookup$level == "module") != 13L || sum(lookup$level == "supermodule") != 9L) {
    stop("Canonical microglia lookup must contain exactly 13 module and 9 supermodule rows.", call. = FALSE)
  }
  if (!is.null(dataset) && any(lookup$dataset != dataset)) stop("Canonical lookup contains the wrong dataset.", call. = FALSE)
  if (!identical(lookup$canonical_plot_label, lookup$final_plot_label)) stop("final_plot_label must be an exact compatibility alias.", call. = FALSE)
  if (any(is.na(wgcna_clean_registry_text(lookup$canonical_biological_label))) || any(is.na(wgcna_clean_registry_text(lookup$canonical_short_label)))) {
    stop("Canonical lookup contains empty reviewed labels.", call. = FALSE)
  }
  if (anyDuplicated(lookup[lookup$level == "supermodule", c("dataset", "entity_id")])) stop("Duplicate supermodule rows are not allowed.", call. = FALSE)
  unstable_claim <- lookup$level == "supermodule" & lookup$structural_coherence_class %in% c("model_singular", "model_unstable") & lookup$aggregation_evidence_class == "claim_supported"
  if (any(unstable_claim, na.rm = TRUE)) stop("Structurally unstable supermodules cannot be claim-supported.", call. = FALSE)
  if (!is.null(member_map)) {
    current_dataset <- if (!is.null(dataset) && length(dataset) && !is.na(dataset)) dataset else unique(lookup$dataset)[[1]]
    current_map <- wgcna_normalize_current_member_map(member_map, current_dataset)
    expected_fingerprints <- current_map |>
      dplyr::group_by(.data$dataset, .data$SupermoduleID) |>
      dplyr::summarise(member_ModuleIDs = paste(sort(.data$ModuleID), collapse = ";"), .groups = "drop") |>
      dplyr::mutate(expected = wgcna_sha256_text(paste(.data$dataset, .data$SupermoduleID, .data$member_ModuleIDs, sep = "|")))
    observed <- lookup |>
      dplyr::filter(.data$level == "supermodule") |>
      dplyr::select(SupermoduleID = "entity_id", observed = "member_modules_fingerprint") |>
      dplyr::left_join(expected_fingerprints, by = "SupermoduleID", relationship = "one-to-one")
    if (any(observed$observed != observed$expected)) stop("Canonical lookup membership fingerprints do not match current membership.", call. = FALSE)
  }
  invisible(TRUE)
}
