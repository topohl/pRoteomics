# Deterministic downstream GSEA-program <-> WGCNA concordance helpers.
#
# These functions consume existing inferential outputs. They do not fit models,
# recompute WGCNA modules, alter source p-values/FDRs, or create a combined test.

gww_formal_contrast_map <- function() {
  c(
    RES_vs_CON = "RES - CON",
    SUS_vs_CON = "SUS - CON",
    SUS_vs_RES = "SUS - RES"
  )
}

gww_formal_contrast <- function(x) {
  x <- as.character(x)
  map <- gww_formal_contrast_map()
  out <- unname(map[x])
  already <- x %in% unname(map)
  out[already] <- x[already]
  out
}

gww_normalize_spatial_unit <- function(x) {
  tolower(trimws(gsub("\\s+", "_", as.character(x))))
}

gww_assert_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(label, " is missing required column(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

gww_assert_unique <- function(x, key, label) {
  gww_assert_columns(x, key, label)
  duplicated_rows <- which(duplicated(x[key]))
  if (length(duplicated_rows)) {
    stop(label, " contains duplicated key(s): ", paste(key, collapse = " + "),
         "; first duplicated row ", duplicated_rows[[1]], ".", call. = FALSE)
  }
  invisible(TRUE)
}

gww_validate_gsea_terms <- function(x) {
  gww_assert_columns(
    x,
    c(
      "dataset", "phenotype_contrast", "spatial_unit", "program_class",
      "Comparison", "ID", "Description", "NES", "pvalue", "p.adjust",
      "core_enrichment", "evidence_source_family"
    ),
    "Canonical spatial ranked-GSEA term table"
  )
  allowed <- names(gww_formal_contrast_map())
  bad <- !is.na(x$phenotype_contrast) & !x$phenotype_contrast %in% allowed
  if (any(bad)) {
    stop("Canonical spatial ranked-GSEA term table contains unsupported contrast(s): ",
         paste(unique(x$phenotype_contrast[bad]), collapse = ", "), ".",
         call. = FALSE)
  }
  if (any(as.character(x$evidence_source_family) != "ranked_GSEA", na.rm = TRUE)) {
    stop("All GSEA rows must retain evidence_source_family = ranked_GSEA.",
         call. = FALSE)
  }
  invisible(TRUE)
}

gww_prepare_supported_terms <- function(x) {
  gww_validate_gsea_terms(x)
  x |>
    dplyr::mutate(
      phenotype_contrast = as.character(.data$phenotype_contrast),
      contrast = gww_formal_contrast(.data$phenotype_contrast),
      gsea_spatial_unit_source = as.character(.data$spatial_unit),
      gsea_spatial_unit = gww_normalize_spatial_unit(.data$spatial_unit),
      biological_program = as.character(.data$program_class),
      term_direction_sign = dplyr::case_when(
        is.finite(.data$NES) & .data$NES > 0 ~ 1L,
        is.finite(.data$NES) & .data$NES < 0 ~ -1L,
        TRUE ~ 0L
      ),
      term_occurrence_key = paste(
        .data$dataset, .data$Comparison, .data$ID, sep = "|"
      )
    ) |>
    dplyr::filter(
      .data$phenotype_contrast %in% names(gww_formal_contrast_map()),
      !is.na(.data$p.adjust), is.finite(.data$p.adjust), .data$p.adjust < 0.05,
      .data$term_direction_sign != 0L,
      !is.na(.data$biological_program),
      nzchar(.data$biological_program),
      .data$biological_program != "Other"
    )
}

gww_representative_rows <- function(x, group_columns, direction_sign = NULL,
                                    prefix = "") {
  if (!is.null(direction_sign)) {
    x <- dplyr::filter(x, .data$term_direction_sign == direction_sign)
  }
  out <- x |>
    dplyr::arrange(
      dplyr::across(dplyr::all_of(group_columns)),
      .data$p.adjust, dplyr::desc(abs(.data$NES)),
      .data$ID, .data$Description, .data$Comparison
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(group_columns)),
      GO_ID = as.character(.data$ID),
      GO_description = as.character(.data$Description),
      NES = as.numeric(.data$NES),
      raw_p = as.numeric(.data$pvalue),
      GSEA_FDR = as.numeric(.data$p.adjust),
      leading_edge_accessions = as.character(.data$core_enrichment),
      source_comparison = as.character(.data$Comparison),
      gsea_source_key = as.character(.data$term_occurrence_key),
      source_supplementary_file = as.character(.data$source_supplementary_file)
    )
  if (nzchar(prefix)) {
    names(out)[!names(out) %in% group_columns] <- paste0(
      prefix, names(out)[!names(out) %in% group_columns]
    )
  }
  out
}

gww_build_local_gsea_evidence <- function(term_table) {
  supported <- gww_prepare_supported_terms(term_table)
  group_columns <- c(
    "dataset", "phenotype_contrast", "contrast",
    "gsea_spatial_unit_source", "gsea_spatial_unit", "biological_program"
  )
  total_units <- term_table |>
    dplyr::filter(.data$phenotype_contrast %in% names(gww_formal_contrast_map())) |>
    dplyr::mutate(
      contrast = gww_formal_contrast(.data$phenotype_contrast),
      gsea_spatial_unit = gww_normalize_spatial_unit(.data$spatial_unit),
      biological_program = as.character(.data$program_class)
    ) |>
    dplyr::filter(.data$biological_program != "Other") |>
    dplyr::distinct(
      .data$dataset, .data$contrast, .data$biological_program,
      .data$gsea_spatial_unit
    ) |>
    dplyr::count(
      .data$dataset, .data$contrast, .data$biological_program,
      name = "n_total_spatial_units"
    )

  summary <- supported |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) |>
    dplyr::summarise(
      n_supported_term_occurrences = dplyr::n_distinct(.data$term_occurrence_key),
      n_supported_terms = dplyr::n_distinct(.data$ID),
      n_positive_supported_terms = dplyr::n_distinct(
        .data$ID[.data$term_direction_sign == 1L]
      ),
      n_negative_supported_terms = dplyr::n_distinct(
        .data$ID[.data$term_direction_sign == -1L]
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      gsea_direction_status = dplyr::case_when(
        .data$n_positive_supported_terms > 0L &
          .data$n_negative_supported_terms > 0L ~ "mixed_direction",
        .data$n_positive_supported_terms > 0L ~ "positive_only",
        .data$n_negative_supported_terms > 0L ~ "negative_only",
        TRUE ~ "unresolved"
      ),
      gsea_direction_sign = dplyr::case_when(
        .data$gsea_direction_status == "positive_only" ~ 1L,
        .data$gsea_direction_status == "negative_only" ~ -1L,
        TRUE ~ NA_integer_
      )
    )

  overall <- gww_representative_rows(supported, group_columns)
  positive <- gww_representative_rows(
    supported, group_columns, direction_sign = 1L, prefix = "positive_"
  )
  negative <- gww_representative_rows(
    supported, group_columns, direction_sign = -1L, prefix = "negative_"
  )

  out <- summary |>
    dplyr::left_join(overall, by = group_columns, relationship = "one-to-one") |>
    dplyr::left_join(positive, by = group_columns, relationship = "one-to-one") |>
    dplyr::left_join(negative, by = group_columns, relationship = "one-to-one") |>
    dplyr::left_join(
      total_units,
      by = c("dataset", "contrast", "biological_program"),
      relationship = "many-to-one"
    )

  recurrence <- out |>
    dplyr::group_by(.data$dataset, .data$contrast, .data$biological_program) |>
    dplyr::summarise(
      n_supported_spatial_units = dplyr::n_distinct(.data$gsea_spatial_unit),
      n_positive_direction_units = dplyr::n_distinct(
        .data$gsea_spatial_unit[.data$gsea_direction_status == "positive_only"]
      ),
      n_negative_direction_units = dplyr::n_distinct(
        .data$gsea_spatial_unit[.data$gsea_direction_status == "negative_only"]
      ),
      n_mixed_direction_units = dplyr::n_distinct(
        .data$gsea_spatial_unit[.data$gsea_direction_status == "mixed_direction"]
      ),
      .groups = "drop"
    )

  out <- out |>
    dplyr::left_join(
      recurrence,
      by = c("dataset", "contrast", "biological_program"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      n_units_supporting_same_direction = dplyr::case_when(
        .data$gsea_direction_sign == 1L ~ .data$n_positive_direction_units,
        .data$gsea_direction_sign == -1L ~ .data$n_negative_direction_units,
        TRUE ~ NA_integer_
      ),
      fraction_spatial_units_same_direction =
        .data$n_units_supporting_same_direction / .data$n_total_spatial_units,
      comparison_scope = "local_local",
      gsea_evidence_id = paste(
        .data$dataset, .data$contrast, .data$gsea_spatial_unit,
        .data$biological_program, "local", sep = "|"
      ),
      representative_selection_rule = paste(
        "original GSEA p.adjust < 0.05; lowest p.adjust;",
        "largest abs(NES); GO ID; description; comparison"
      ),
      evidence_source_family = "ranked_GSEA",
      gsea_representation_role = "canonical_spatial_term_summary"
    )
  gww_assert_unique(out, "gsea_evidence_id", "Local GSEA program evidence")
  out
}

gww_build_global_gsea_evidence <- function(term_table, local_evidence,
                                           min_recurrent_units = 2L) {
  supported <- gww_prepare_supported_terms(term_table)
  recurrence <- local_evidence |>
    dplyr::group_by(.data$dataset, .data$phenotype_contrast, .data$contrast,
                    .data$biological_program) |>
    dplyr::summarise(
      n_total_spatial_units = dplyr::first(.data$n_total_spatial_units),
      n_supported_spatial_units = dplyr::n_distinct(.data$gsea_spatial_unit),
      n_positive_direction_units = dplyr::n_distinct(
        .data$gsea_spatial_unit[.data$gsea_direction_status == "positive_only"]
      ),
      n_negative_direction_units = dplyr::n_distinct(
        .data$gsea_spatial_unit[.data$gsea_direction_status == "negative_only"]
      ),
      n_mixed_direction_units = dplyr::n_distinct(
        .data$gsea_spatial_unit[.data$gsea_direction_status == "mixed_direction"]
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      gsea_direction_status = dplyr::case_when(
        .data$n_positive_direction_units >= min_recurrent_units &
          .data$n_positive_direction_units > .data$n_negative_direction_units ~
          "positive_recurrent",
        .data$n_negative_direction_units >= min_recurrent_units &
          .data$n_negative_direction_units > .data$n_positive_direction_units ~
          "negative_recurrent",
        .data$n_positive_direction_units >= min_recurrent_units &
          .data$n_negative_direction_units >= min_recurrent_units &
          .data$n_positive_direction_units == .data$n_negative_direction_units ~
          "mixed_recurrent",
        TRUE ~ "insufficient_recurrence"
      ),
      gsea_direction_sign = dplyr::case_when(
        .data$gsea_direction_status == "positive_recurrent" ~ 1L,
        .data$gsea_direction_status == "negative_recurrent" ~ -1L,
        TRUE ~ NA_integer_
      ),
      n_units_supporting_same_direction = pmax(
        .data$n_positive_direction_units,
        .data$n_negative_direction_units
      ),
      fraction_spatial_units_same_direction =
        .data$n_units_supporting_same_direction / .data$n_total_spatial_units
    ) |>
    dplyr::filter(.data$gsea_direction_status != "insufficient_recurrence")

  group_columns <- c(
    "dataset", "phenotype_contrast", "contrast", "biological_program"
  )
  term_summary <- supported |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) |>
    dplyr::summarise(
      n_supported_term_occurrences = dplyr::n_distinct(.data$term_occurrence_key),
      n_supported_terms = dplyr::n_distinct(.data$ID),
      n_positive_supported_terms = dplyr::n_distinct(
        .data$ID[.data$term_direction_sign == 1L]
      ),
      n_negative_supported_terms = dplyr::n_distinct(
        .data$ID[.data$term_direction_sign == -1L]
      ),
      .groups = "drop"
    )

  candidate <- supported |>
    dplyr::inner_join(
      recurrence[group_columns |> c("gsea_direction_sign")],
      by = group_columns,
      relationship = "many-to-one"
    ) |>
    dplyr::filter(
      is.na(.data$gsea_direction_sign) |
        .data$term_direction_sign == .data$gsea_direction_sign
    )
  overall <- gww_representative_rows(candidate, group_columns)
  positive <- gww_representative_rows(
    supported, group_columns, direction_sign = 1L, prefix = "positive_"
  )
  negative <- gww_representative_rows(
    supported, group_columns, direction_sign = -1L, prefix = "negative_"
  )

  out <- recurrence |>
    dplyr::inner_join(term_summary, by = group_columns, relationship = "one-to-one") |>
    dplyr::left_join(overall, by = group_columns, relationship = "one-to-one") |>
    dplyr::left_join(positive, by = group_columns, relationship = "one-to-one") |>
    dplyr::left_join(negative, by = group_columns, relationship = "one-to-one") |>
    dplyr::mutate(
      gsea_spatial_unit_source = "recurrent_across_spatial_units",
      gsea_spatial_unit = "global_spatial_adjusted",
      comparison_scope = "global_cross_spatial",
      gsea_evidence_id = paste(
        .data$dataset, .data$contrast, .data$gsea_spatial_unit,
        .data$biological_program, "global", sep = "|"
      ),
      representative_selection_rule = paste(
        "recurrent direction in >=", min_recurrent_units,
        "spatial units; original GSEA p.adjust < 0.05; lowest p.adjust;",
        "largest abs(NES); GO ID; description; comparison"
      ),
      evidence_source_family = "ranked_GSEA",
      gsea_representation_role = "recurrent_cross_spatial_summary"
    )
  gww_assert_unique(out, "gsea_evidence_id", "Global GSEA program evidence")
  out
}

gww_read_program_mapping <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  gww_assert_columns(
    x, c("annotation_token", "biological_program", "mapping_rationale"),
    "Program-module mapping registry"
  )
  x |>
    dplyr::mutate(
      registry_annotation_token = as.character(.data$annotation_token),
      annotation_token_normalized = tolower(trimws(.data$annotation_token)),
      biological_program = as.character(.data$biological_program)
    ) |>
    dplyr::select(-"annotation_token") |>
    dplyr::distinct(
      .data$annotation_token_normalized, .data$biological_program,
      .keep_all = TRUE
    )
}

gww_build_entity_program_matches <- function(handoff, registry) {
  gww_assert_columns(
    handoff,
    c(
      "dataset", "entity_level", "entity_id", "display_label",
      "module_program_primary", "module_program_secondary",
      "independent_hypothesis", "claim_entity_role"
    ),
    "WGCNA inferential handoff"
  )
  if (any(!handoff$independent_hypothesis %in% TRUE) ||
      any(handoff$claim_entity_role == "compatibility_alias")) {
    stop("Program matching accepts only independent non-alias WGCNA entities.",
         call. = FALSE)
  }
  entity_columns <- c(
    "dataset", "entity_level", "entity_id", "display_label",
    "module_program_primary", "module_program_secondary",
    "cleaned_biological_label"
  )
  for (column in setdiff(entity_columns, names(handoff))) {
    handoff[[column]] <- NA_character_
  }
  entities <- handoff |>
    dplyr::distinct(dplyr::across(dplyr::all_of(entity_columns)))
  counts <- entities |>
    dplyr::count(.data$dataset, .data$entity_level, .data$entity_id)
  if (any(counts$n != 1L)) {
    stop("WGCNA biological annotations are not constant within entity.",
         call. = FALSE)
  }

  annotations <- dplyr::bind_rows(
    entities |>
      dplyr::transmute(
        dplyr::across(dplyr::all_of(c(
          "dataset", "entity_level", "entity_id", "display_label",
          "cleaned_biological_label"
        ))),
        annotation_field = "module_program_primary",
        annotation_priority = 1L,
        annotation_token = as.character(.data$module_program_primary)
      ),
    entities |>
      dplyr::transmute(
        dplyr::across(dplyr::all_of(c(
          "dataset", "entity_level", "entity_id", "display_label",
          "cleaned_biological_label"
        ))),
        annotation_field = "module_program_secondary",
        annotation_priority = 2L,
        annotation_token = as.character(.data$module_program_secondary)
      )
  ) |>
    dplyr::mutate(
      annotation_token_normalized = tolower(trimws(.data$annotation_token))
    ) |>
    dplyr::filter(
      !is.na(.data$annotation_token_normalized),
      nzchar(.data$annotation_token_normalized)
    )

  out <- annotations |>
    dplyr::inner_join(
      registry,
      by = "annotation_token_normalized",
      relationship = "many-to-many"
    ) |>
    dplyr::arrange(
      .data$dataset, .data$entity_level, .data$entity_id,
      .data$biological_program, .data$annotation_priority,
      .data$annotation_token_normalized
    ) |>
    dplyr::group_by(
      .data$dataset, .data$entity_level, .data$entity_id,
      .data$biological_program
    ) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      program_module_match_rule = paste(
        "exact normalized Stage 07 curated primary/secondary annotation token",
        "to prespecified registry"
      ),
      module_annotation_evidence_role =
        "biological_identity_only_not_group_effect_evidence"
    )
  gww_assert_unique(
    out,
    c("dataset", "entity_level", "entity_id", "biological_program"),
    "Program-module matches"
  )
  out
}

gww_normalize_accession <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("_MOUSE$", "", x)
}

gww_build_universe_bundle <- function(universe, dataset) {
  gww_assert_columns(
    universe,
    c("ProteinGroupID", "included_in_wgcna", "ModuleID", "member_accessions"),
    paste(dataset, "canonical WGCNA feature universe")
  )
  x <- universe |>
    dplyr::filter(.data$included_in_wgcna %in% TRUE) |>
    dplyr::transmute(
      dataset = dataset,
      ProteinGroupID = as.character(.data$ProteinGroupID),
      ModuleID = as.character(.data$ModuleID),
      member_accessions = as.character(.data$member_accessions),
      RepresentativeUniProt = if ("RepresentativeUniProt" %in% names(universe))
        as.character(.data$RepresentativeUniProt) else NA_character_,
      MemberUniProts = if ("MemberUniProts" %in% names(universe))
        as.character(.data$MemberUniProts) else NA_character_,
      representative_accession = if ("representative_accession" %in% names(universe))
        as.character(.data$representative_accession) else NA_character_
    )
  if (any(is.na(x$ProteinGroupID) | !nzchar(x$ProteinGroupID)) ||
      any(is.na(x$ModuleID) | !nzchar(x$ModuleID))) {
    stop(dataset, " canonical WGCNA universe has missing identity/module values.",
         call. = FALSE)
  }
  gww_assert_unique(x, "ProteinGroupID", paste(dataset, "WGCNA universe"))

  token_rows <- lapply(seq_len(nrow(x)), function(i) {
    values <- paste(
      x$member_accessions[[i]], x$RepresentativeUniProt[[i]],
      x$MemberUniProts[[i]], x$representative_accession[[i]], sep = ";"
    )
    tokens <- unique(gww_normalize_accession(
      unlist(strsplit(values, "[/;,|[:space:]]+"), use.names = FALSE)
    ))
    tokens <- tokens[!is.na(tokens) & nzchar(tokens) & tokens != "NA"]
    data.frame(
      dataset = dataset,
      ProteinGroupID = x$ProteinGroupID[[i]],
      token = tokens,
      stringsAsFactors = FALSE
    )
  })
  token_map <- dplyr::bind_rows(token_rows) |>
    dplyr::distinct(.data$dataset, .data$ProteinGroupID, .data$token)
  ambiguous <- token_map |>
    dplyr::count(.data$dataset, .data$token) |>
    dplyr::filter(.data$n > 1L)
  if (nrow(ambiguous)) {
    stop(dataset, " WGCNA universe maps accession token(s) to multiple ",
         "ProteinGroupIDs; first token: ", ambiguous$token[[1]], ".",
         call. = FALSE)
  }
  list(
    universe = x,
    token_map = token_map,
    module_sets = split(x$ProteinGroupID, x$ModuleID)
  )
}

gww_map_leading_edges <- function(evidence, token_map) {
  gww_assert_unique(evidence, "gsea_evidence_id", "GSEA evidence")
  expanded <- lapply(seq_len(nrow(evidence)), function(i) {
    tokens <- unique(gww_normalize_accession(
      unlist(strsplit(
        as.character(evidence$leading_edge_accessions[[i]]),
        "[/;,|[:space:]]+"
      ), use.names = FALSE)
    ))
    tokens <- tokens[!is.na(tokens) & nzchar(tokens) & tokens != "NA"]
    data.frame(
      dataset = evidence$dataset[[i]],
      gsea_evidence_id = evidence$gsea_evidence_id[[i]],
      leading_edge_token = tokens,
      stringsAsFactors = FALSE
    )
  })
  expanded <- dplyr::bind_rows(expanded)
  mapped <- expanded |>
    dplyr::left_join(
      token_map,
      by = c("dataset", "leading_edge_token" = "token"),
      relationship = "many-to-one"
    )
  summary <- mapped |>
    dplyr::group_by(.data$dataset, .data$gsea_evidence_id) |>
    dplyr::summarise(
      n_leading_edge_source = dplyr::n_distinct(.data$leading_edge_token),
      n_leading_edge_mapped =
        dplyr::n_distinct(.data$ProteinGroupID, na.rm = TRUE),
      leading_edge_mapping_fraction =
        dplyr::n_distinct(.data$leading_edge_token[!is.na(.data$ProteinGroupID)]) /
        .data$n_leading_edge_source,
      leading_edge_proteins = paste(
        sort(unique(stats::na.omit(.data$ProteinGroupID))), collapse = ";"
      ),
      unmapped_leading_edge_tokens = paste(
        sort(unique(.data$leading_edge_token[is.na(.data$ProteinGroupID)])),
        collapse = ";"
      ),
      .groups = "drop"
    )
  out <- evidence |>
    dplyr::left_join(
      summary,
      by = c("dataset", "gsea_evidence_id"),
      relationship = "one-to-one"
    )
  if (any(is.na(out$n_leading_edge_mapped) |
          out$n_leading_edge_mapped < 1L)) {
    stop("At least one representative leading edge has no proteins in the ",
         "canonical WGCNA universe.", call. = FALSE)
  }
  out
}

gww_split_proteins <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(x)) return(character())
  unique(trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)))
}

gww_overlap_stat <- function(module_set, leading_edge_set, universe_set) {
  universe_set <- unique(as.character(universe_set))
  module_set <- intersect(unique(as.character(module_set)), universe_set)
  leading_edge_set <- intersect(unique(as.character(leading_edge_set)), universe_set)
  overlap <- intersect(module_set, leading_edge_set)
  a <- length(overlap)
  b <- length(leading_edge_set) - a
  c <- length(module_set) - a
  d <- length(universe_set) - a - b - c
  if (d < 0L) stop("Invalid overlap contingency table.", call. = FALSE)
  p <- stats::fisher.test(
    matrix(c(a, b, c, d), nrow = 2L), alternative = "greater"
  )$p.value
  data.frame(
    n_universe = length(universe_set),
    n_module = length(module_set),
    n_leading_edge = length(leading_edge_set),
    n_overlap = a,
    overlap_proteins = paste(sort(overlap), collapse = ";"),
    Jaccard = if (length(union(module_set, leading_edge_set)))
      a / length(union(module_set, leading_edge_set)) else NA_real_,
    fisher_p_value = as.numeric(p),
    stringsAsFactors = FALSE
  )
}

gww_program_specific_overlap <- function(evidence, matches, bundles) {
  modules <- matches |>
    dplyr::filter(.data$entity_level == "module") |>
    dplyr::select(dplyr::all_of(c(
      "dataset", "entity_id", "biological_program", "annotation_field",
      "annotation_token", "program_module_match_rule"
    )))
  hypotheses <- evidence |>
    dplyr::inner_join(
      modules,
      by = c("dataset", "biological_program"),
      relationship = "many-to-many"
    ) |>
    dplyr::arrange(
      .data$dataset, .data$contrast, .data$gsea_spatial_unit,
      .data$biological_program, .data$entity_id
    )
  rows <- lapply(seq_len(nrow(hypotheses)), function(i) {
    row <- hypotheses[i, , drop = FALSE]
    bundle <- bundles[[row$dataset[[1]]]]
    module_set <- bundle$module_sets[[row$entity_id[[1]]]]
    if (is.null(module_set)) {
      stop("Matched module is absent from canonical universe: ",
           row$dataset[[1]], " / ", row$entity_id[[1]], ".", call. = FALSE)
    }
    stat <- gww_overlap_stat(
      module_set,
      gww_split_proteins(row$leading_edge_proteins[[1]]),
      bundle$universe$ProteinGroupID
    )
    cbind(
      row[c(
        "gsea_evidence_id", "dataset", "phenotype_contrast", "contrast",
        "comparison_scope", "gsea_spatial_unit", "biological_program",
        "GO_ID", "GO_description", "NES", "GSEA_FDR",
        "gsea_source_key", "entity_id", "annotation_field",
        "annotation_token", "program_module_match_rule"
      )],
      stat,
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows) |>
    dplyr::group_by(
      .data$dataset, .data$contrast, .data$gsea_spatial_unit
    ) |>
    dplyr::mutate(
      overlap_family_id = paste(
        "program_module_overlap", .data$dataset,
        gsub(" ", "", gsub("-", "_", .data$contrast)),
        .data$gsea_spatial_unit, sep = "::"
      ),
      overlap_family_size = dplyr::n(),
      overlap_FDR = stats::p.adjust(.data$fisher_p_value, method = "BH")
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      overlap_universe = "canonical_dataset_WGCNA_feature_universe",
      overlap_test_scope = "program_specific_representative_term_x_module"
    )
  gww_assert_unique(
    out,
    c("gsea_evidence_id", "dataset", "entity_id"),
    "Program-specific overlap hypotheses"
  )
  family_check <- out |>
    dplyr::count(.data$overlap_family_id, name = "observed_family_size") |>
    dplyr::left_join(
      out |>
        dplyr::distinct(.data$overlap_family_id, .data$overlap_family_size),
      by = "overlap_family_id",
      relationship = "one-to-one"
    )
  if (any(family_check$observed_family_size != family_check$overlap_family_size)) {
    stop("Overlap family sizes do not reconcile.", call. = FALSE)
  }
  out
}

gww_effect_scale_thresholds <- function(effects) {
  gww_assert_columns(
    effects,
    c(
      "dataset", "entity_level", "analysis_tier", "effect_scope",
      "estimate", "model_valid", "contrast"
    ),
    "WGCNA effects"
  )
  effects |>
    dplyr::filter(
      .data$contrast %in% unname(gww_formal_contrast_map()),
      .data$model_valid %in% TRUE,
      is.finite(.data$estimate)
    ) |>
    dplyr::group_by(
      .data$dataset, .data$entity_level, .data$analysis_tier,
      .data$effect_scope
    ) |>
    dplyr::summarise(
      effect_scale_n = dplyr::n(),
      effect_scale_q25_abs = as.numeric(stats::quantile(
        abs(.data$estimate), probs = 0.25, names = FALSE, type = 7
      )),
      effect_scale_median_abs = as.numeric(stats::median(abs(.data$estimate))),
      effect_scale_rule = paste(
        "within dataset x entity_level x analysis_tier x effect_scope;",
        "q25(abs estimate)=near-zero boundary; median(abs estimate)=",
        "meaningful-effect compatibility reference"
      ),
      .groups = "drop"
    )
}

gww_classify_concordance <- function(x) {
  gww_assert_columns(
    x,
    c(
      "gsea_direction_sign", "wgcna_estimate", "wgcna_CI_low",
      "wgcna_CI_high", "wgcna_tier_specific_fdr", "wgcna_model_valid",
      "wgcna_model_stability_status", "animal_instability_flag",
      "effect_scale_q25_abs", "effect_scale_median_abs"
    ),
    "Concordance classification input"
  )
  sign_w <- sign(x$wgcna_estimate)
  sign_g <- x$gsea_direction_sign
  concordant <- !is.na(sign_g) & sign_w == sign_g & sign_w != 0
  opposite <- !is.na(sign_g) & sign_w == -sign_g & sign_w != 0
  near_zero <- is.finite(x$wgcna_estimate) &
    is.finite(x$effect_scale_q25_abs) &
    abs(x$wgcna_estimate) <= x$effect_scale_q25_abs
  ci_reaches_median <- dplyr::case_when(
    sign_g == 1L ~ x$wgcna_CI_high >= x$effect_scale_median_abs,
    sign_g == -1L ~ x$wgcna_CI_low <= -x$effect_scale_median_abs,
    TRUE ~ FALSE
  )
  clearly_opposite <- dplyr::case_when(
    sign_g == 1L ~ x$wgcna_CI_high < 0,
    sign_g == -1L ~ x$wgcna_CI_low > 0,
    TRUE ~ FALSE
  )
  stable <- x$wgcna_model_valid %in% TRUE &
    !grepl(
      "invalid|unstable|fail|nonestimable",
      tolower(as.character(x$wgcna_model_stability_status))
    )
  fdr_supported <- stable & is.finite(x$wgcna_tier_specific_fdr) &
    x$wgcna_tier_specific_fdr <= 0.05
  animal_sensitive <- x$animal_instability_flag %in% TRUE

  x$effect_direction_concordant <- concordant
  x$wgcna_effect_near_zero <- near_zero
  x$wgcna_ci_reaches_median_effect_in_gsea_direction <- ci_reaches_median
  x$wgcna_ci_clearly_opposes_gsea_direction <- clearly_opposite
  x$wgcna_endpoint_FDR_supported <- fdr_supported
  x$concordance_class <- dplyr::case_when(
    is.na(sign_g) ~ "unresolved",
    concordant & animal_sensitive ~ "animal_sensitive",
    concordant & fdr_supported ~ "FDR_supported_concordance",
    opposite & clearly_opposite & stable ~ "discordant",
    concordant & stable & !fdr_supported & !near_zero &
      ci_reaches_median & !animal_sensitive ~ "concordant_imprecise",
    stable & near_zero ~ "weak_or_near_zero_module_support",
    TRUE ~ "unresolved"
  )
  x$concordance_class_reason <- dplyr::case_when(
    x$concordance_class == "FDR_supported_concordance" ~
      "GSEA direction concordant with a valid tier-specific FDR-supported WGCNA endpoint.",
    x$concordance_class == "concordant_imprecise" ~
      "Direction concordant; estimate above the within-tier near-zero boundary and CI reaches the median observed effect, but WGCNA is not FDR-supported.",
    x$concordance_class == "weak_or_near_zero_module_support" ~
      "WGCNA absolute estimate is at or below the within-tier lower-quartile effect boundary.",
    x$concordance_class == "discordant" ~
      "WGCNA CI lies entirely in the direction opposite to the supported GSEA program.",
    x$concordance_class == "animal_sensitive" ~
      "Apparent directional concordance has an existing leave-one-animal-out or equivalent animal-instability flag.",
    TRUE ~
      "Direction is mixed, model/robustness evidence is insufficient, or the CI/effect geometry does not meet another deterministic class."
  )
  x
}

gww_validate_local_global_semantics <- function(x) {
  gww_assert_columns(
    x,
    c(
      "comparison_scope", "gsea_spatial_unit", "wgcna_spatial_unit",
      "wgcna_effect_scope", "wgcna_analysis_tier"
    ),
    "Concordance table"
  )
  local <- x$comparison_scope == "local_local"
  bad_local <- local & (
    x$gsea_spatial_unit != x$wgcna_spatial_unit |
      x$wgcna_analysis_tier != "exploratory_spatial_localization" |
      x$wgcna_effect_scope != "within_spatial_unit"
  )
  global <- x$comparison_scope == "global_cross_spatial"
  bad_global <- global & (
    x$gsea_spatial_unit != "global_spatial_adjusted" |
      x$wgcna_spatial_unit != "global_spatial_adjusted" |
      x$wgcna_effect_scope != "spatial_adjusted_global" |
      !x$wgcna_analysis_tier %in% c(
        "primary_wgcna_global", "secondary_contextual_global"
      )
  )
  if (any(bad_local | bad_global)) {
    stop("Local/global concordance semantics are inconsistent.", call. = FALSE)
  }
  invisible(TRUE)
}

gww_class_rank <- function(x) {
  match(
    x,
    c(
      "FDR_supported_concordance", "concordant_imprecise",
      "weak_or_near_zero_module_support", "animal_sensitive",
      "discordant", "unresolved"
    )
  )
}

gww_program_summary <- function(long) {
  long |>
    dplyr::group_by(
      .data$dataset, .data$contrast, .data$comparison_scope,
      .data$biological_program
    ) |>
    dplyr::summarise(
      gsea_direction_status = dplyr::first(.data$gsea_direction_status),
      gsea_direction_sign = dplyr::first(.data$gsea_direction_sign),
      n_gsea_evidence_rows = dplyr::n_distinct(.data$gsea_evidence_id),
      n_wgcna_entities = dplyr::n_distinct(paste(
        .data$wgcna_entity_level, .data$wgcna_entity_id, sep = "|"
      )),
      n_wgcna_modules = dplyr::n_distinct(
        .data$wgcna_entity_id[.data$wgcna_entity_level == "module"]
      ),
      n_wgcna_supermodules = dplyr::n_distinct(
        .data$wgcna_entity_id[.data$wgcna_entity_level == "supermodule"]
      ),
      median_wgcna_estimate = stats::median(.data$wgcna_estimate, na.rm = TRUE),
      n_FDR_supported_concordance = sum(
        .data$concordance_class == "FDR_supported_concordance"
      ),
      n_concordant_imprecise = sum(
        .data$concordance_class == "concordant_imprecise"
      ),
      n_weak_or_near_zero = sum(
        .data$concordance_class == "weak_or_near_zero_module_support"
      ),
      n_discordant = sum(.data$concordance_class == "discordant"),
      n_animal_sensitive = sum(.data$concordance_class == "animal_sensitive"),
      n_unresolved = sum(.data$concordance_class == "unresolved"),
      best_overlap_FDR = if (all(is.na(.data$overlap_FDR))) NA_real_ else
        min(.data$overlap_FDR, na.rm = TRUE),
      best_concordance_class = {
        classes <- unique(.data$concordance_class)
        classes[[which.min(gww_class_rank(classes))]]
      },
      evidence_source_family = "ranked_GSEA;WGCNA_group_effect;WGCNA_overlap",
      independence_note = paste(
        "GSEA/compareGO/program summaries are one ranked_GSEA lineage;",
        "module annotation is identity only; shared-group contrasts are not",
        "independent replications"
      ),
      .groups = "drop"
    )
}

gww_adaptive_pattern_summary <- function(long) {
  global_modules <- long |>
    dplyr::filter(
      .data$comparison_scope == "global_cross_spatial",
      .data$wgcna_entity_level == "module"
    ) |>
    dplyr::group_by(
      .data$dataset, .data$biological_program, .data$contrast
    ) |>
    dplyr::summarise(
      gsea_direction_sign = dplyr::first(.data$gsea_direction_sign),
      gsea_direction_status = dplyr::first(.data$gsea_direction_status),
      n_matched_modules = dplyr::n_distinct(.data$wgcna_entity_id),
      median_wgcna_estimate = stats::median(.data$wgcna_estimate, na.rm = TRUE),
      wgcna_effect_min = min(.data$wgcna_estimate, na.rm = TRUE),
      wgcna_effect_max = max(.data$wgcna_estimate, na.rm = TRUE),
      wgcna_CI_envelope_low = min(.data$wgcna_CI_low, na.rm = TRUE),
      wgcna_CI_envelope_high = max(.data$wgcna_CI_high, na.rm = TRUE),
      n_stable_concordant = sum(.data$concordance_class %in% c(
        "FDR_supported_concordance", "concordant_imprecise"
      )),
      n_animal_sensitive = sum(.data$concordance_class == "animal_sensitive"),
      n_discordant = sum(.data$concordance_class == "discordant"),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = "contrast",
      values_from = c(
        "gsea_direction_sign", "gsea_direction_status", "n_matched_modules",
        "median_wgcna_estimate", "wgcna_effect_min", "wgcna_effect_max",
        "wgcna_CI_envelope_low", "wgcna_CI_envelope_high",
        "n_stable_concordant", "n_animal_sensitive", "n_discordant"
      ),
      names_glue = "{.value}__{contrast}"
    )

  needed_contrasts <- unname(gww_formal_contrast_map())
  bases <- c(
    "gsea_direction_sign", "gsea_direction_status", "n_matched_modules",
    "median_wgcna_estimate", "wgcna_effect_min", "wgcna_effect_max",
    "wgcna_CI_envelope_low", "wgcna_CI_envelope_high",
    "n_stable_concordant", "n_animal_sensitive", "n_discordant"
  )
  for (base in bases) {
    for (contrast in needed_contrasts) {
      column <- paste(base, contrast, sep = "__")
      if (!column %in% names(global_modules)) global_modules[[column]] <- NA
    }
  }
  get_col <- function(base, contrast) global_modules[[paste(base, contrast, sep = "__")]]
  res <- get_col("gsea_direction_sign", "RES - CON")
  sus <- get_col("gsea_direction_sign", "SUS - CON")
  direct <- get_col("gsea_direction_sign", "SUS - RES")
  res_ok <- get_col("n_stable_concordant", "RES - CON") > 0
  sus_ok <- get_col("n_stable_concordant", "SUS - CON") > 0
  direct_ok <- get_col("n_stable_concordant", "SUS - RES") > 0
  res_ok[is.na(res_ok)] <- FALSE
  sus_ok[is.na(sus_ok)] <- FALSE
  direct_ok[is.na(direct_ok)] <- FALSE

  res_gt_con_gt_sus <- res == 1L & sus == -1L & direct == -1L &
    res_ok & sus_ok & direct_ok
  global_modules$RES_gt_CON_gt_SUS_supported <-
    dplyr::coalesce(res_gt_con_gt_sus, FALSE)
  global_modules$adaptive_resilience_pattern <- dplyr::case_when(
    res == 1L & direct == -1L & res_ok & direct_ok ~
      "RES_specific_adaptive_candidate",
    !is.na(sus) & sus != 0L & direct == sus & sus_ok & direct_ok ~
      "SUS_specific_candidate",
    !is.na(res) & !is.na(sus) & res == -sus & res_ok & sus_ok &
      direct == sus & direct_ok ~ "opposing_RES_vs_SUS",
    !is.na(res) & !is.na(sus) & res == sus & !is.na(direct) &
      direct != 0L & res_ok & sus_ok & direct_ok ~ "graded_stress_response",
    !is.na(res) & !is.na(sus) & res == sus & res_ok & sus_ok &
      (!direct_ok | is.na(direct)) ~ "shared_stress_response",
    TRUE ~ "unresolved"
  )
  global_modules$pattern_rule <- paste(
    "Requires recurrent ranked-GSEA direction plus at least one stable",
    "direction-concordant global module endpoint for each contrast used;",
    "SUS-RES is the direct group-difference endpoint; significance/non-",
    "significance asymmetry is never used. CI envelope is descriptive, not pooled."
  )
  global_modules
}
