# Shared module contracts and identifier helpers.

normalize_module_identifier <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(";.*$", "", x)
  x <- gsub("\\|.*$", "", x)
  x[is.na(x)] <- ""
  x
}

normalize_gene_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x)] <- ""
  x
}

wgcna_publication_module_palette <- function() {
  c(
    deep_steel_blue = "#486A8A",
    muted_sky_blue = "#7FA6C1",
    slate_blue = "#5D7894",
    blue_grey = "#8AA0AF",
    deep_teal = "#2F6F73",
    sea_teal = "#5B9992",
    pale_teal = "#9CC9BF",
    dark_sea_green = "#4E7D66",
    sage_green = "#7FA37A",
    olive_green = "#8A9364",
    moss_green = "#637A52",
    lichen_green = "#A8B58A",
    muted_ochre = "#B08A4A",
    sand_ochre = "#C5A46D",
    clay = "#A66E5A",
    terracotta = "#B67861",
    umber = "#7F6654",
    taupe = "#9A8977",
    warm_grey = "#8E8A82",
    cool_grey = "#7E8A8F",
    charcoal_grey = "#5D6468",
    light_grey = "#C8CDD0",
    mauve = "#9B7895",
    dusty_lilac = "#B39AB2",
    heather = "#7D6F8F",
    plum_grey = "#6F5F73",
    muted_rose = "#B9878A",
    dusty_raspberry = "#9B6672",
    spruce = "#3F6B5B",
    petrol_blue = "#3F6E82",
    celadon = "#86B49D",
    graphite = "#4F5457"
  )
}

wgcna_module_color_metadata <- function(module_colors, palette = wgcna_publication_module_palette()) {
  module_colors <- sort(unique(as.character(module_colors)))
  module_colors <- module_colors[!is.na(module_colors) & nzchar(module_colors)]
  if (length(module_colors) > length(palette)) {
    stop(
      "WGCNA publication palette has ", length(palette), " colours but ",
      length(module_colors), " modules were detected; add unique colours instead of recycling.",
      call. = FALSE
    )
  }
  idx <- seq_along(module_colors)
  module_ids <- sprintf("WGCNA_m%02d", idx)
  color_names <- names(palette)[idx]
  color_labels <- gsub("_", " ", color_names)
  data.frame(
    ModuleID = module_ids,
    ModuleLegacyID = paste0("WGCNA_", module_colors),
    ModuleColor = unname(palette[idx]),
    ModuleColorName = color_names,
    ModuleColorLabel = color_labels,
    WGCNAInternalColor = module_colors,
    stringsAsFactors = FALSE
  )
}

validate_wgcna_module_color_metadata <- function(module_meta, observed_internal_colors = NULL,
                                                 artifact = "WGCNA module colour metadata") {
  required <- c(
    "WGCNAInternalColor", "ModuleID", "ModuleLegacyID", "ModuleColor",
    "ModuleColorName", "ModuleColorLabel"
  )
  require_module_contract_columns(module_meta, required, artifact)
  for (nm in c("WGCNAInternalColor", "ModuleID", "ModuleColor")) {
    values <- trimws(as.character(module_meta[[nm]]))
    if (any(is.na(values) | !nzchar(values))) {
      stop(artifact, " contains missing ", nm, " values.", call. = FALSE)
    }
  }
  if (any(!grepl("^WGCNA_m[0-9]{2,}$", as.character(module_meta$ModuleID)))) {
    stop(artifact, " must use stable ModuleID values such as WGCNA_m01.", call. = FALSE)
  }
  if (anyDuplicated(as.character(module_meta$WGCNAInternalColor))) {
    stop(artifact, " must map WGCNAInternalColor one-to-one to ModuleID.", call. = FALSE)
  }
  if (anyDuplicated(as.character(module_meta$ModuleID))) {
    stop(artifact, " contains duplicated ModuleID values.", call. = FALSE)
  }
  if (anyDuplicated(as.character(module_meta$ModuleColor))) {
    stop(artifact, " must map ModuleID one-to-one to publication ModuleColor.", call. = FALSE)
  }
  if (!is.null(observed_internal_colors)) {
    observed <- unique(as.character(observed_internal_colors))
    observed <- observed[!is.na(observed) & nzchar(observed)]
    missing <- setdiff(observed, as.character(module_meta$WGCNAInternalColor))
    if (length(missing)) {
      stop(
        artifact, " is missing internal colour(s) observed in mergedColors: ",
        paste(missing, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

wgcna_integrate_feature_module_metadata <- function(protein_group_ids, merged_colors,
                                                    module_meta, feature_table = NULL) {
  protein_group_ids <- as.character(protein_group_ids)
  if (anyNA(protein_group_ids) || any(!nzchar(protein_group_ids)) || anyDuplicated(protein_group_ids)) {
    stop("WGCNA feature integration requires unique, nonmissing ProteinGroupID values.", call. = FALSE)
  }
  internal_colors <- if (!is.null(names(merged_colors))) {
    as.character(merged_colors[protein_group_ids])
  } else {
    if (length(merged_colors) != length(protein_group_ids)) {
      stop("mergedColors and ProteinGroupID must have the same length when mergedColors is unnamed.", call. = FALSE)
    }
    as.character(merged_colors)
  }
  if (anyNA(internal_colors) || any(!nzchar(internal_colors))) {
    stop("mergedColors contains missing WGCNA internal colours for one or more ProteinGroupID values.", call. = FALSE)
  }
  validate_wgcna_module_color_metadata(module_meta, internal_colors)
  out <- tibble::tibble(
    ProteinGroupID = protein_group_ids,
    WGCNAInternalColor = internal_colors
  ) |>
    dplyr::left_join(
      module_meta,
      by = "WGCNAInternalColor",
      relationship = "many-to-one"
    )
  if (nrow(out) != length(protein_group_ids)) {
    stop("WGCNA module metadata integration changed the feature row count.", call. = FALSE)
  }
  if (anyNA(out$ModuleID) || any(!grepl("^WGCNA_m[0-9]{2,}$", out$ModuleID))) {
    stop("WGCNA feature module metadata contains missing or unstable ModuleID values.", call. = FALSE)
  }
  if (!is.null(feature_table)) {
    require_module_contract_columns(feature_table, "ProteinGroupID", "WGCNA network feature table")
    if (anyDuplicated(as.character(feature_table$ProteinGroupID))) {
      stop("WGCNA network feature table contains duplicated ProteinGroupID values.", call. = FALSE)
    }
    out <- out |>
      dplyr::left_join(feature_table, by = "ProteinGroupID", relationship = "many-to-one")
    if (nrow(out) != length(protein_group_ids)) {
      stop("WGCNA feature annotation join changed the feature row count.", call. = FALSE)
    }
  }
  out
}

validate_wgcna_module_label_table <- function(module_label_table,
                                              artifact = "WGCNA module label table") {
  validate_wgcna_module_color_metadata(module_label_table, artifact = artifact)
  if (anyDuplicated(as.character(module_label_table$ModuleID))) {
    stop(artifact, " must contain exactly one row per ModuleID.", call. = FALSE)
  }
  invisible(TRUE)
}

wgcna_join_module_labels <- function(feature_module_tbl, module_label_table) {
  require_module_contract_columns(feature_module_tbl, "ModuleID", "WGCNA feature module table")
  validate_wgcna_module_label_table(module_label_table)
  label_cols <- setdiff(
    names(module_label_table),
    c("WGCNAInternalColor", "ModuleLegacyID", "ModuleColor", "ModuleColorName", "ModuleColorLabel")
  )
  out <- feature_module_tbl |>
    dplyr::left_join(
      module_label_table[, label_cols, drop = FALSE],
      by = "ModuleID",
      relationship = "many-to-one"
    )
  if (nrow(out) != nrow(feature_module_tbl)) {
    stop("Joining WGCNA module labels changed the feature row count.", call. = FALSE)
  }
  out
}

wgcna_apply_module_metadata <- function(df, module_meta,
                                        internal_color_col = "ModuleColor",
                                        dataset_col = NULL) {
  if (!internal_color_col %in% names(df)) return(df)
  meta <- module_meta
  if (!is.null(dataset_col) && dataset_col %in% names(df) && dataset_col %in% names(meta)) {
    by_cols <- stats::setNames(c(dataset_col, "WGCNAInternalColor"), c(dataset_col, internal_color_col))
  } else {
    by_cols <- stats::setNames("WGCNAInternalColor", internal_color_col)
  }
  drop_cols <- intersect(c("ModuleID", "ModuleLegacyID", "ModuleColor", "ModuleColorName", "ModuleColorLabel"), names(df))
  df |>
    dplyr::select(-dplyr::any_of(drop_cols)) |>
    dplyr::left_join(meta, by = by_cols) |>
    dplyr::select(-dplyr::any_of("WGCNAInternalColor"))
}

wgcna_normalize_module_ids <- function(df, module_lookup = NULL, id_col = "ModuleID",
                                       legacy_col = "ModuleLegacyID") {
  if (!id_col %in% names(df) && legacy_col %in% names(df)) df[[id_col]] <- df[[legacy_col]]
  if (!id_col %in% names(df)) return(df)
  df[[id_col]] <- as.character(df[[id_col]])
  if (!legacy_col %in% names(df)) df[[legacy_col]] <- NA_character_
  legacy_like <- grepl("^WGCNA_#", df[[id_col]])
  df[[legacy_col]][legacy_like & (is.na(df[[legacy_col]]) | !nzchar(df[[legacy_col]]))] <- df[[id_col]][legacy_like]
  if (!is.null(module_lookup) && all(c("ModuleID", "ModuleLegacyID") %in% names(module_lookup))) {
    hit <- match(df[[id_col]], module_lookup$ModuleLegacyID)
    replace <- legacy_like & !is.na(hit)
    df[[id_col]][replace] <- module_lookup$ModuleID[hit[replace]]
    miss_legacy <- is.na(df[[legacy_col]]) | !nzchar(df[[legacy_col]])
    hit2 <- match(df[[id_col]], module_lookup$ModuleID)
    fill <- miss_legacy & !is.na(hit2)
    df[[legacy_col]][fill] <- module_lookup$ModuleLegacyID[hit2[fill]]
  }
  df
}

require_module_contract_columns <- function(df, cols, artifact = "artifact") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing)) {
    stop(
      artifact, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

wgcna_join_supermodule_hub_handoff <- function(wgcna_modules_long,
                                               supermodule_annotation,
                                               merged_me_names) {
  module_artifact <- "WGCNA hub-overlap module rows"
  lookup_artifact <- "WGCNA hub-overlap supermodule lookup"
  require_module_contract_columns(
    wgcna_modules_long,
    c(
      "ProteinGroupID", "ModuleID", "WGCNAInternalColor", "ModuleColor",
      "abs_kME"
    ),
    module_artifact
  )
  require_module_contract_columns(
    supermodule_annotation,
    c(
      "dataset", "ModuleID", "module_eigengene", "SupermoduleID",
      "Supermodule_DisplayLabel", "present_in_dataset"
    ),
    lookup_artifact
  )
  if (anyDuplicated(as.character(wgcna_modules_long$ProteinGroupID))) {
    stop(module_artifact, " contains duplicated ProteinGroupID rows.", call. = FALSE)
  }
  if (anyNA(wgcna_modules_long$ModuleID) || any(!nzchar(as.character(wgcna_modules_long$ModuleID)))) {
    stop(module_artifact, " contains missing ModuleID values.", call. = FALSE)
  }

  lookup <- supermodule_annotation |>
    dplyr::filter(.data$present_in_dataset %in% TRUE) |>
    dplyr::select(
      "dataset", "ModuleID", "module_eigengene",
      "SupermoduleID", "Supermodule_DisplayLabel"
    )
  if (!nrow(lookup)) {
    stop(lookup_artifact, " contains no present modules.", call. = FALSE)
  }
  for (nm in c("dataset", "ModuleID", "module_eigengene", "SupermoduleID")) {
    values <- trimws(as.character(lookup[[nm]]))
    if (any(is.na(values) | !nzchar(values))) {
      stop(lookup_artifact, " contains missing ", nm, " values.", call. = FALSE)
    }
  }
  lookup_key <- paste(as.character(lookup$dataset), as.character(lookup$ModuleID), sep = "::")
  if (anyDuplicated(lookup_key)) {
    stop(lookup_artifact, " contains duplicated dataset + ModuleID keys.", call. = FALSE)
  }
  eigengenes_per_module <- tapply(
    as.character(lookup$module_eigengene),
    as.character(lookup$ModuleID),
    function(x) length(unique(x))
  )
  if (any(eigengenes_per_module != 1L)) {
    stop(lookup_artifact, " must map every ModuleID to exactly one module_eigengene.", call. = FALSE)
  }
  supermodules_per_module <- tapply(
    as.character(lookup$SupermoduleID),
    as.character(lookup$ModuleID),
    function(x) length(unique(x))
  )
  if (any(supermodules_per_module != 1L)) {
    stop(lookup_artifact, " must map every ModuleID to exactly one SupermoduleID.", call. = FALSE)
  }

  merged_me_names <- unique(as.character(merged_me_names))
  if (!length(merged_me_names) || anyNA(merged_me_names) || any(!nzchar(merged_me_names))) {
    stop("WGCNA hub-overlap handoff requires valid merged eigengene column names.", call. = FALSE)
  }
  unexpected_eigengenes <- setdiff(unique(as.character(lookup$module_eigengene)), merged_me_names)
  if (length(unexpected_eigengenes)) {
    stop(
      lookup_artifact, " contains module_eigengene values absent from mergedMEs: ",
      paste(unexpected_eigengenes, collapse = ", "), ".",
      call. = FALSE
    )
  }

  join_cols <- if ("dataset" %in% names(wgcna_modules_long)) {
    dataset_values <- trimws(as.character(wgcna_modules_long$dataset))
    if (any(is.na(dataset_values) | !nzchar(dataset_values))) {
      stop(module_artifact, " contains missing dataset values.", call. = FALSE)
    }
    c("dataset", "ModuleID")
  } else {
    lookup_datasets <- unique(as.character(lookup$dataset))
    if (length(lookup_datasets) != 1L) {
      stop(
        "ModuleID-only hub-overlap joins require exactly one dataset in the authoritative lookup.",
        call. = FALSE
      )
    }
    "ModuleID"
  }

  before_n <- nrow(wgcna_modules_long)
  joined <- wgcna_modules_long |>
    dplyr::left_join(lookup, by = join_cols, relationship = "many-to-one")
  if (nrow(joined) != before_n) {
    stop("WGCNA hub-overlap lookup join changed the protein row count.", call. = FALSE)
  }
  missing_mapping <- is.na(joined$module_eigengene) |
    !nzchar(as.character(joined$module_eigengene)) |
    is.na(joined$SupermoduleID) |
    !nzchar(as.character(joined$SupermoduleID))
  if (any(missing_mapping)) {
    missing_ids <- sort(unique(as.character(joined$ModuleID[missing_mapping])))
    stop(
      "WGCNA hub-overlap lookup is missing module mappings for ModuleID: ",
      paste(missing_ids, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (any(!as.character(joined$module_eigengene) %in% merged_me_names)) {
    stop("WGCNA hub-overlap join produced module_eigengene values absent from mergedMEs.", call. = FALSE)
  }
  joined
}

wgcna_build_hub_module_sets <- function(wgcna_modules_long,
                                        supermodule_annotation,
                                        merged_me_names,
                                        top_n = 25L) {
  top_n <- suppressWarnings(as.integer(top_n))
  if (length(top_n) != 1L || is.na(top_n) || top_n < 1L) {
    stop("top_n must be one positive integer.", call. = FALSE)
  }
  joined <- wgcna_join_supermodule_hub_handoff(
    wgcna_modules_long,
    supermodule_annotation,
    merged_me_names
  )
  out <- joined |>
    dplyr::group_by(
      .data$dataset, .data$ModuleID, .data$module_eigengene,
      .data$WGCNAInternalColor, .data$ModuleColor, .data$SupermoduleID
    ) |>
    dplyr::arrange(dplyr::desc(.data$abs_kME), .by_group = TRUE) |>
    dplyr::mutate(
      .hub_keep = if ("is_top_hub_25" %in% names(joined)) {
        dplyr::coalesce(.data$is_top_hub_25, FALSE) | dplyr::row_number() <= top_n
      } else {
        dplyr::row_number() <= top_n
      }
    ) |>
    dplyr::filter(.data$.hub_keep) |>
    dplyr::summarise(
      Supermodule_DisplayLabel = dplyr::first(.data$Supermodule_DisplayLabel),
      top_hub_proteins = list(unique(as.character(.data$ProteinGroupID))),
      .groups = "drop"
    )
  expected_groups <- length(unique(paste(joined$dataset, joined$ModuleID, sep = "::")))
  if (nrow(out) != expected_groups) {
    stop("WGCNA hub-overlap grouping did not produce exactly one hub set per module.", call. = FALSE)
  }
  out
}

wgcna_supermodule_sensitivity_grid <- function() {
  c(0.25, 0.35, 0.45, 0.50, 0.55, 0.65)
}

supermodule_contract_key <- function(dataset, supermodule_id) {
  paste(as.character(dataset), as.character(supermodule_id), sep = "::")
}

validate_supermodule_member_map <- function(df, expected_modules = NULL,
                                            artifact = "authoritative supermodule member map",
                                            dataset_col = "dataset",
                                            module_col = "module_eigengene",
                                            id_col = "SupermoduleID",
                                            display_col = NULL) {
  require_module_contract_columns(df, c(dataset_col, module_col, id_col), artifact)
  dataset <- trimws(as.character(df[[dataset_col]]))
  modules <- trimws(as.character(df[[module_col]]))
  ids <- trimws(as.character(df[[id_col]]))
  if (any(is.na(dataset) | !nzchar(dataset))) stop(artifact, " contains missing dataset keys.", call. = FALSE)
  if (any(is.na(modules) | !nzchar(modules))) stop(artifact, " contains missing member modules.", call. = FALSE)
  if (any(is.na(ids) | !nzchar(ids))) stop(artifact, " contains missing SupermoduleID values.", call. = FALSE)
  if (any(!grepl("^SM[0-9]{2,}$", ids))) {
    stop(artifact, " contains unexpected SupermoduleID values; expected stable IDs such as SM01.", call. = FALSE)
  }
  member_key <- paste(dataset, modules, sep = "::")
  if (anyDuplicated(member_key)) {
    stop(artifact, " contains duplicate dataset + member-module keys.", call. = FALSE)
  }
  for (ds in unique(dataset)) {
    observed_ids <- sort(unique(ids[dataset == ds]))
    expected_ids <- sprintf("SM%02d", seq_along(observed_ids))
    if (!identical(observed_ids, expected_ids)) {
      stop(
        artifact, " contains unexpected IDs for ", ds, ": ", paste(observed_ids, collapse = ", "),
        "; expected exactly ", paste(expected_ids, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  if (!is.null(expected_modules)) {
    expected_modules <- sort(unique(as.character(expected_modules)))
    observed_modules <- sort(unique(modules))
    missing <- setdiff(expected_modules, observed_modules)
    unexpected <- setdiff(observed_modules, expected_modules)
    if (length(missing) || length(unexpected)) {
      stop(
        artifact, " does not exactly match the authoritative module universe; missing: ",
        paste(missing, collapse = ", "), "; unexpected: ", paste(unexpected, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  if (!is.null(display_col)) {
    require_module_contract_columns(df, display_col, artifact)
    labels <- trimws(as.character(df[[display_col]]))
    keep <- !is.na(labels) & nzchar(labels)
    label_map <- unique(data.frame(dataset = dataset[keep], SupermoduleID = ids[keep], label = labels[keep], stringsAsFactors = FALSE))
    if (nrow(label_map)) {
      id_key <- supermodule_contract_key(label_map$dataset, label_map$SupermoduleID)
      labels_per_id <- tapply(label_map$label, id_key, function(x) length(unique(x)))
      if (any(labels_per_id > 1L)) stop(artifact, " contains multiple display labels for one dataset + SupermoduleID.", call. = FALSE)
      label_key <- paste(label_map$dataset, label_map$label, sep = "::")
      ids_per_label <- tapply(label_map$SupermoduleID, label_key, function(x) length(unique(x)))
      if (any(ids_per_label > 1L)) stop(artifact, " contains a display-label collision across distinct SupermoduleID values.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

validate_supermodule_summary_ids <- function(summary_df, authoritative_map,
                                             artifact = "supermodule summary") {
  require_module_contract_columns(summary_df, c("dataset", "SupermoduleID"), artifact)
  require_module_contract_columns(authoritative_map, c("dataset", "SupermoduleID"), "authoritative supermodule map")
  summary_keys <- supermodule_contract_key(summary_df$dataset, summary_df$SupermoduleID)
  if (anyDuplicated(summary_keys)) stop(artifact, " contains duplicate dataset + SupermoduleID rows.", call. = FALSE)
  expected <- sort(unique(supermodule_contract_key(authoritative_map$dataset, authoritative_map$SupermoduleID)))
  observed <- sort(unique(summary_keys))
  if (!identical(observed, expected)) {
    stop(
      artifact, " SupermoduleID set does not exactly match the authoritative cluster map; missing: ",
      paste(setdiff(expected, observed), collapse = ", "), "; unexpected: ",
      paste(setdiff(observed, expected), collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

require_supermodule_composition_columns <- function(df,
                                                    artifact = "supermodule biological annotation") {
  required <- c(
    "dataset", "SupermoduleID", "Supermodule_CompositionLabel",
    "Supermodule_CompositionDisplayLabel", "Supermodule_CompositionConfidence",
    "Supermodule_CompositionRationale"
  )
  require_module_contract_columns(df, required, artifact)
  for (nm in c("Supermodule_CompositionLabel", "Supermodule_CompositionDisplayLabel")) {
    values <- trimws(as.character(df[[nm]]))
    if (nrow(df) && all(is.na(values) | !nzchar(values))) {
      stop(artifact, " has a stale or empty ", nm, " column; regenerate module annotation before interpretation.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

wgcna_supermodule_go_support <- function(go_df, member_map, fdr_threshold = 0.05) {
  go_names <- if (is.null(go_df)) character() else names(go_df)
  module_key <- if ("ModuleID" %in% names(member_map) && "ModuleID" %in% go_names) {
    "ModuleID"
  } else if ("WGCNAInternalColor" %in% names(member_map) && "WGCNAInternalColor" %in% go_names) {
    "WGCNAInternalColor"
  } else {
    "ModuleColor"
  }
  required_map <- c("dataset", "module_eigengene", module_key, "SupermoduleID")
  require_module_contract_columns(member_map, required_map, "supermodule GO member map")
  validate_supermodule_member_map(
    member_map, artifact = "supermodule GO member map",
    display_col = if ("Supermodule_DisplayLabel" %in% names(member_map)) "Supermodule_DisplayLabel" else NULL
  )
  empty <- tibble::tibble(
    dataset = character(), SupermoduleID = character(), Ontology = character(),
    Description = character(), ModuleProteinSetType = character(),
    n_member_modules = integer(), n_modules_supporting_term = integer(),
    fraction_member_modules_supporting_term = numeric(),
    best_member_module_fdr = numeric(), worst_member_module_fdr = numeric(),
    supporting_member_modules = character(), go_support_confidence = character(),
    recurring_significant_term = logical(), support_fdr_threshold = numeric()
  )
  if (is.null(go_df) || !nrow(go_df)) return(empty)
  require_module_contract_columns(
    go_df,
    c("ModuleProteinSetType", module_key, "Ontology", "Description", "p.adjust"),
    "module GO enrichment"
  )
  map <- member_map |>
    dplyr::select(dplyr::all_of(required_map)) |>
    dplyr::distinct()
  map_key <- c("dataset", module_key)
  if (anyDuplicated(map[map_key])) {
    stop(
      "supermodule GO member map has duplicated ",
      paste(map_key, collapse = " + "),
      " keys; GO evidence requires a many-to-one module lookup.",
      call. = FALSE
    )
  }
  module_counts <- map |>
    dplyr::count(.data$dataset, .data$SupermoduleID, name = "n_member_modules")
  if (!"dataset" %in% names(go_df) && dplyr::n_distinct(map$dataset) != 1L) {
    stop(
      "module GO enrichment lacks dataset, but the supermodule member map contains multiple datasets; ",
      "dataset + ModuleID is required.",
      call. = FALSE
    )
  }
  go_join_cols <- if ("dataset" %in% names(go_df)) c("dataset", module_key) else module_key
  lookup <- map
  go_filtered <- go_df |>
    dplyr::filter(
      .data$ModuleProteinSetType == "all",
      !is.na(.data$Description), nzchar(as.character(.data$Description))
    ) |>
    dplyr::mutate(p.adjust = suppressWarnings(as.numeric(.data$p.adjust)))
  n_go_filtered <- nrow(go_filtered)
  module_terms <- go_filtered |>
    dplyr::left_join(lookup, by = go_join_cols, relationship = "many-to-one")
  if (nrow(module_terms) != n_go_filtered) {
    stop("supermodule GO lookup changed the number of module GO rows.", call. = FALSE)
  }
  if (anyNA(module_terms$SupermoduleID) || anyNA(module_terms$module_eigengene)) {
    missing_keys <- unique(module_terms[[module_key]][
      is.na(module_terms$SupermoduleID) | is.na(module_terms$module_eigengene)
    ])
    stop(
      "module GO enrichment contains modules absent from the supermodule member map: ",
      paste(missing_keys, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  module_terms <- module_terms |>
    dplyr::group_by(
      .data$dataset, .data$SupermoduleID, .data$module_eigengene,
      .data$Ontology, .data$Description
    ) |>
    dplyr::summarise(
      member_module_fdr = {
        z <- .data$p.adjust[is.finite(.data$p.adjust)]
        if (length(z)) min(z) else NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(significant_module_support = is.finite(.data$member_module_fdr) & .data$member_module_fdr <= fdr_threshold)
  if (!nrow(module_terms)) return(empty)
  module_terms |>
    dplyr::group_by(.data$dataset, .data$SupermoduleID, .data$Ontology, .data$Description) |>
    dplyr::summarise(
      n_modules_supporting_term = sum(.data$significant_module_support),
      best_member_module_fdr = {
        z <- .data$member_module_fdr[.data$significant_module_support]
        if (length(z)) min(z) else NA_real_
      },
      worst_member_module_fdr = {
        z <- .data$member_module_fdr[.data$significant_module_support]
        if (length(z)) max(z) else NA_real_
      },
      supporting_member_modules = paste(sort(unique(.data$module_eigengene[.data$significant_module_support])), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      module_counts,
      by = c("dataset", "SupermoduleID"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      ModuleProteinSetType = "all",
      fraction_member_modules_supporting_term = .data$n_modules_supporting_term / .data$n_member_modules,
      go_support_confidence = dplyr::case_when(
        .data$n_member_modules <= 1L ~ "singleton",
        .data$n_modules_supporting_term == .data$n_member_modules ~ "high",
        .data$n_modules_supporting_term >= 2L & .data$fraction_member_modules_supporting_term >= 0.5 ~ "medium",
        TRUE ~ "not_recurring"
      ),
      recurring_significant_term = .data$go_support_confidence %in% c("high", "medium"),
      support_fdr_threshold = fdr_threshold
    ) |>
    dplyr::select(dplyr::all_of(c(
      "dataset", "SupermoduleID", "Ontology", "Description",
      "ModuleProteinSetType", "n_member_modules",
      "n_modules_supporting_term", "fraction_member_modules_supporting_term",
      "best_member_module_fdr", "worst_member_module_fdr",
      "supporting_member_modules", "go_support_confidence",
      "recurring_significant_term", "support_fdr_threshold"
    ))) |>
    dplyr::arrange(
      .data$dataset, .data$SupermoduleID,
      factor(.data$go_support_confidence, levels = c("high", "medium", "not_recurring", "singleton")),
      dplyr::desc(.data$fraction_member_modules_supporting_term),
      .data$worst_member_module_fdr, .data$Ontology, .data$Description
    )
}

wgcna_supermodule_go_evidence_summary <- function(go_df, member_map, fdr_threshold = 0.05) {
  require_module_contract_columns(
    member_map,
    c("dataset", "ModuleID", "module_eigengene", "SupermoduleID"),
    "supermodule GO evidence member map"
  )
  members <- member_map |>
    dplyr::group_by(.data$dataset, .data$SupermoduleID) |>
    dplyr::summarise(
      member_ModuleIDs = paste(sort(unique(as.character(.data$ModuleID))), collapse = ";"),
      n_member_modules = dplyr::n_distinct(.data$ModuleID),
      .groups = "drop"
    )
  support <- wgcna_supermodule_go_support(
    go_df = go_df,
    member_map = member_map,
    fdr_threshold = fdr_threshold
  )
  recurring <- support |>
    dplyr::filter(.data$recurring_significant_term) |>
    dplyr::group_by(.data$dataset, .data$SupermoduleID) |>
    dplyr::summarise(
      n_recurring_significant_GO_terms = dplyr::n(),
      n_modules_with_GO_support = max(.data$n_modules_supporting_term),
      n_high_confidence_GO_terms = sum(.data$go_support_confidence == "high"),
      n_medium_confidence_GO_terms = sum(.data$go_support_confidence == "medium"),
      .groups = "drop"
    )
  members |>
    dplyr::left_join(
      recurring,
      by = c("dataset", "SupermoduleID"),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      n_recurring_significant_GO_terms = dplyr::coalesce(.data$n_recurring_significant_GO_terms, 0L),
      n_modules_with_GO_support = dplyr::coalesce(.data$n_modules_with_GO_support, 0L),
      n_high_confidence_GO_terms = dplyr::coalesce(.data$n_high_confidence_GO_terms, 0L),
      n_medium_confidence_GO_terms = dplyr::coalesce(.data$n_medium_confidence_GO_terms, 0L),
      GO_label_confidence_class = dplyr::case_when(
        .data$n_member_modules <= 1L ~ "singleton",
        .data$n_high_confidence_GO_terms > 0L ~ "high",
        .data$n_medium_confidence_GO_terms > 0L ~ "medium",
        TRUE ~ "mixed_or_unresolved"
      ),
      Supermodule_NameSource = dplyr::case_when(
        .data$n_member_modules <= 1L ~ "singleton",
        .data$n_recurring_significant_GO_terms > 0L ~ "recurring_significant_GO",
        TRUE ~ "no_recurring_significant_GO"
      ),
      ManualReviewRequired = .data$Supermodule_NameSource != "recurring_significant_GO"
    ) |>
    dplyr::arrange(.data$dataset, .data$SupermoduleID)
}

wgcna_supermodule_coherence <- function(mergedMEs, member_map,
                                        module_members = NULL,
                                        sample_info = NULL,
                                        sensitivity = NULL) {
  require_module_contract_columns(member_map, c("dataset", "module_eigengene", "SupermoduleID"), "supermodule coherence member map")
  module_names <- colnames(mergedMEs)
  validate_supermodule_member_map(
    member_map,
    expected_modules = module_names,
    artifact = "supermodule coherence member map",
    display_col = if ("Supermodule_DisplayLabel" %in% names(member_map)) "Supermodule_DisplayLabel" else NULL
  )
  me <- as.data.frame(mergedMEs[, member_map$module_eigengene, drop = FALSE], check.names = FALSE)
  signed_corr <- stats::cor(me, use = "pairwise.complete.obs", method = "pearson")

  adjusted <- matrix(NA_real_, nrow = nrow(me), ncol = ncol(me), dimnames = dimnames(as.matrix(me)))
  adjustment_variables <- character()
  adjustment_status <- "not_available"
  if (!is.null(sample_info) && nrow(sample_info) == nrow(me)) {
    candidates <- names(sample_info)[tolower(names(sample_info)) %in% c("region", "layer")]
    candidates <- candidates[vapply(sample_info[candidates], function(x) length(unique(stats::na.omit(as.character(x)))) > 1L, logical(1))]
    if (length(candidates)) {
      covars <- sample_info[, candidates, drop = FALSE]
      keep <- stats::complete.cases(covars)
      if (sum(keep) >= max(5L, length(candidates) + 2L)) {
        design <- stats::model.matrix(~ ., data = covars[keep, , drop = FALSE])
        for (nm in names(me)) {
          y <- suppressWarnings(as.numeric(me[[nm]][keep]))
          ok <- is.finite(y) & apply(design, 1L, function(z) all(is.finite(z)))
          if (sum(ok) > ncol(design)) {
            adjusted[which(keep)[ok], nm] <- stats::lm.fit(design[ok, , drop = FALSE], y[ok])$residuals
          }
        }
        adjustment_variables <- candidates
        adjustment_status <- "adjusted"
      }
    }
  }
  adjusted_corr <- if (identical(adjustment_status, "adjusted")) {
    stats::cor(adjusted, use = "pairwise.complete.obs", method = "pearson")
  } else {
    matrix(NA_real_, ncol = ncol(me), nrow = ncol(me), dimnames = list(names(me), names(me)))
  }

  protein_counts <- NULL
  if (!is.null(module_members) && nrow(module_members)) {
    require_module_contract_columns(module_members, "ProteinGroupID", "WGCNA module membership")
    join_col <- if ("module_eigengene" %in% names(module_members)) {
      "module_eigengene"
    } else if ("ModuleID" %in% names(module_members) && "ModuleID" %in% names(member_map)) {
      "ModuleID"
    } else if ("WGCNAInternalColor" %in% names(module_members) && "WGCNAInternalColor" %in% names(member_map)) {
      "WGCNAInternalColor"
    } else if ("ModuleColor" %in% names(module_members) && "ModuleColor" %in% names(member_map)) {
      "ModuleColor"
    } else {
      NA_character_
    }
    if (is.na(join_col)) stop("WGCNA module membership requires module_eigengene, ModuleID, or WGCNAInternalColor for supermodule protein counts.", call. = FALSE)
    map_for_join <- member_map |>
      dplyr::select(dplyr::all_of(c("dataset", join_col, "SupermoduleID")))
    join_cols <- if ("dataset" %in% names(module_members)) c("dataset", join_col) else join_col
    protein_joined <- dplyr::inner_join(module_members, map_for_join, by = join_cols)
    if (!"dataset" %in% names(protein_joined)) protein_joined$dataset <- unique(member_map$dataset)[[1]]
    protein_counts <- protein_joined |>
      dplyr::group_by(.data$dataset, .data$SupermoduleID) |>
      dplyr::summarise(n_member_proteins = dplyr::n_distinct(.data$ProteinGroupID), .groups = "drop")
  }

  groups <- split(seq_len(nrow(member_map)), supermodule_contract_key(member_map$dataset, member_map$SupermoduleID))
  rows <- lapply(groups, function(idx) {
    ds <- as.character(member_map$dataset[idx[[1]]])
    sid <- as.character(member_map$SupermoduleID[idx[[1]]])
    members <- as.character(member_map$module_eigengene[idx])
    pair_values <- if (length(members) >= 2L) signed_corr[members, members][upper.tri(signed_corr[members, members])] else numeric()
    adjusted_values <- if (length(members) >= 2L) adjusted_corr[members, members][upper.tri(adjusted_corr[members, members])] else numeric()
    pair_values <- pair_values[is.finite(pair_values)]
    adjusted_values <- adjusted_values[is.finite(adjusted_values)]
    pc1 <- if (length(members) == 1L) 1 else tryCatch({
      fit <- stats::prcomp(me[, members, drop = FALSE], center = TRUE, scale. = TRUE)
      fit$sdev[[1]]^2 / sum(fit$sdev^2)
    }, error = function(e) NA_real_)
    stability <- NULL
    if (!is.null(sensitivity) && nrow(sensitivity) && all(c("dataset", "cut_height", "matched_primary_supermodule_id", "jaccard_to_primary_supermodule", "stable_primary_match") %in% names(sensitivity))) {
      stability <- sensitivity |>
        dplyr::filter(.data$dataset == ds, .data$matched_primary_supermodule_id == sid) |>
        dplyr::group_by(.data$cut_height) |>
        dplyr::summarise(
          jaccard = max(.data$jaccard_to_primary_supermodule, na.rm = TRUE),
          stable = any(.data$stable_primary_match %in% TRUE),
          .groups = "drop"
        )
    }
    tibble::tibble(
      dataset = ds,
      SupermoduleID = sid,
      n_member_modules = length(members),
      n_member_module_pairs = choose(length(members), 2L),
      member_modules = paste(members, collapse = ";"),
      signed_min_pairwise_eigengene_correlation = if (length(pair_values)) min(pair_values, na.rm = TRUE) else NA_real_,
      signed_mean_pairwise_eigengene_correlation = if (length(pair_values)) mean(pair_values, na.rm = TRUE) else NA_real_,
      signed_median_pairwise_eigengene_correlation = if (length(pair_values)) stats::median(pair_values, na.rm = TRUE) else NA_real_,
      adjusted_signed_min_pairwise_eigengene_correlation = if (length(adjusted_values)) min(adjusted_values) else NA_real_,
      adjusted_signed_mean_pairwise_eigengene_correlation = if (length(adjusted_values)) mean(adjusted_values) else NA_real_,
      adjusted_signed_median_pairwise_eigengene_correlation = if (length(adjusted_values)) stats::median(adjusted_values) else NA_real_,
      correlation_adjustment_variables = paste(adjustment_variables, collapse = ";"),
      correlation_adjustment_status = adjustment_status,
      pc1_variance_explained = pc1,
      pc1_interpretation_note = "PC1 variance is a data-reduction coherence metric, not proof of a shared biological pathway.",
      cut_height_grid_n = if (!is.null(stability)) nrow(stability) else 0L,
      cut_height_stability_min_jaccard = if (!is.null(stability) && any(is.finite(stability$jaccard))) min(stability$jaccard[is.finite(stability$jaccard)]) else NA_real_,
      cut_height_stability_mean_jaccard = if (!is.null(stability) && any(is.finite(stability$jaccard))) mean(stability$jaccard[is.finite(stability$jaccard)]) else NA_real_,
      cut_height_stability_fraction_stable = if (!is.null(stability) && nrow(stability)) mean(stability$stable) else NA_real_,
      correlation_metric = "signed_pearson",
      construction_basis = "average_linkage_on_1_minus_signed_module_eigengene_correlation"
    )
  })
  out <- dplyr::bind_rows(rows)
  if (!is.null(protein_counts)) out <- dplyr::left_join(out, protein_counts, by = c("dataset", "SupermoduleID"))
  if (!"n_member_proteins" %in% names(out)) out$n_member_proteins <- NA_integer_
  out
}

validate_wgcna_module_definitions <- function(df, artifact = "WGCNA module definitions") {
  require_module_contract_columns(
    df,
    c(
      "ModuleSet", "WGCNAInternalColor", "ModuleID", "ModuleLegacyID", "ModuleColor", "ModuleColorName", "ModuleColorLabel",
      "ProteinGroupID", "ProteinID", "member_accessions", "member_gene_symbols",
      "original_identifier", "representative_accession", "representative_gene_symbol",
      "protein_group_ambiguity_class", "n_mapped_accessions", "n_gene_symbols",
      "gene_level_claim_allowed", "protein_level_claim_allowed", "mapping_status", "FeatureDisplayLabel"
    ),
    artifact
  )
  bad_id <- !grepl("^WGCNA_m[0-9]{2,}$", as.character(df$ModuleID))
  if (any(bad_id, na.rm = TRUE)) {
    stop(artifact, " must use stable ModuleID values such as WGCNA_m01, not colour labels.", call. = FALSE)
  }
  if (any(duplicated(unique(df[, c("ModuleID", "ModuleColor"), drop = FALSE])$ModuleID))) {
    stop(artifact, " has non-unique ModuleID to ModuleColor mappings.", call. = FALSE)
  }
  module_map <- unique(df[, c(
    "WGCNAInternalColor", "ModuleID", "ModuleLegacyID", "ModuleColor",
    "ModuleColorName", "ModuleColorLabel"
  ), drop = FALSE])
  validate_wgcna_module_color_metadata(
    module_map,
    observed_internal_colors = df$WGCNAInternalColor,
    artifact = paste0(artifact, " module metadata")
  )
  if (!any(c("kME", "Weight") %in% colnames(df))) {
    stop(artifact, " must contain kME or Weight.", call. = FALSE)
  }
  if (anyNA(df$ProteinGroupID) || any(!nzchar(as.character(df$ProteinGroupID)))) {
    stop(artifact, " contains missing ProteinGroupID values.", call. = FALSE)
  }
  if (any(as.character(df$ProteinID) != as.character(df$ProteinGroupID), na.rm = TRUE)) {
    stop(artifact, " ProteinID must be a deprecated alias of ProteinGroupID.", call. = FALSE)
  }
  invisible(TRUE)
}

wgcna_feature_key_contract_version <- function() "protein_group_id_v1"

wgcna_feature_key_fingerprint <- function(protein_group_ids) {
  ids <- as.character(protein_group_ids)
  if (anyNA(ids) || any(!nzchar(ids))) stop("Cannot fingerprint missing ProteinGroupID values.", call. = FALSE)
  if (anyDuplicated(ids)) stop("Cannot fingerprint duplicate ProteinGroupID values.", call. = FALSE)
  if (!exists("stable_pg_hash", mode = "function")) stop("stable_pg_hash() must be loaded before fingerprinting WGCNA features.", call. = FALSE)
  stable_pg_hash(paste(ids, collapse = "\n"))
}

validate_wgcna_feature_alignment <- function(expression_data, feature_table) {
  require_module_contract_columns(feature_table, "ProteinGroupID", "WGCNA feature table")
  ids <- as.character(feature_table$ProteinGroupID[feature_table$included_in_wgcna %in% TRUE])
  if (anyNA(ids) || any(!nzchar(ids))) stop("WGCNA feature table contains missing ProteinGroupID values.", call. = FALSE)
  if (anyDuplicated(ids)) stop("WGCNA feature table contains duplicate ProteinGroupID values; make.unique() repair is forbidden.", call. = FALSE)
  if (!identical(colnames(expression_data), ids)) stop("Expression features and canonical WGCNA annotation rows are misaligned.", call. = FALSE)
  invisible(TRUE)
}

validate_wgcna_cached_state <- function(state, expected_feature_ids = NULL) {
  rerun <- "Full WGCNA rerun required under the Phase 1B ProteinGroupID feature-key contract."
  if (!is.list(state) || !identical(state$feature_key_contract_version, wgcna_feature_key_contract_version())) {
    stop("Legacy WGCNA cached state lacks the Phase 1B feature-key contract. ", rerun, call. = FALSE)
  }
  expr <- state$expression.data
  ids <- if (!is.null(expr)) colnames(expr) else character()
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Cached WGCNA state lacks valid canonical ProteinGroupID columns. ", rerun, call. = FALSE)
  }
  if (any(!grepl("^PG:", ids))) {
    stop("Cached WGCNA state contains legacy accession or repaired feature identifiers. ", rerun, call. = FALSE)
  }
  fingerprint <- wgcna_feature_key_fingerprint(ids)
  if (!identical(state$feature_key_fingerprint, fingerprint)) {
    stop("Cached WGCNA feature-key fingerprint is invalid. ", rerun, call. = FALSE)
  }
  if (!is.null(expected_feature_ids) && !identical(ids, as.character(expected_feature_ids))) {
    stop("Cached WGCNA feature keys differ from current canonical features. ", rerun, call. = FALSE)
  }
  invisible(TRUE)
}

validate_wgcna_hydrated_go_enrichment <- function(go_enrichment_long,
                                                  module_label_table,
                                                  saved_go_evidence = FALSE) {
  if (is.null(go_enrichment_long)) {
    if (isTRUE(saved_go_evidence)) {
      stop(
        "Cached WGCNA state contains GO enrichment evidence but hydrated GO_enrichment_long is NULL.",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }
  require_module_contract_columns(
    go_enrichment_long,
    c(
      "WGCNAInternalColor", "ModuleID", "ModuleColor", "ModuleProteinSetType",
      "Ontology", "Description", "p.adjust"
    ),
    "hydrated WGCNA GO enrichment"
  )
  validate_wgcna_module_label_table(module_label_table, artifact = "hydrated WGCNA module label table")
  if (nrow(go_enrichment_long)) {
    go_ids <- as.character(go_enrichment_long$ModuleID)
    if (anyNA(go_ids) || any(!grepl("^WGCNA_m[0-9]{2,}$", go_ids))) {
      stop(
        "hydrated WGCNA GO enrichment must use stable ModuleID values, not legacy colour labels.",
        call. = FALSE
      )
    }
    if (anyNA(go_enrichment_long$ModuleColor) || anyNA(go_enrichment_long$WGCNAInternalColor)) {
      stop("hydrated WGCNA GO enrichment contains missing module colour metadata.", call. = FALSE)
    }
  }
  label_lookup <- module_label_table |>
    dplyr::select(dplyr::all_of(c("ModuleID", "WGCNAInternalColor", "ModuleColor"))) |>
    dplyr::rename(
      expected_WGCNAInternalColor = "WGCNAInternalColor",
      expected_ModuleColor = "ModuleColor"
    ) |>
    dplyr::distinct()
  if (anyDuplicated(label_lookup$ModuleID)) {
    stop("hydrated WGCNA module label table has duplicated ModuleID values.", call. = FALSE)
  }
  n_before <- nrow(go_enrichment_long)
  checked <- go_enrichment_long |>
    dplyr::left_join(label_lookup, by = "ModuleID", relationship = "many-to-one")
  if (nrow(checked) != n_before) {
    stop("hydrated WGCNA GO validation changed the number of GO rows.", call. = FALSE)
  }
  if (nrow(checked) && anyNA(checked$expected_ModuleColor)) {
    missing_ids <- unique(checked$ModuleID[is.na(checked$expected_ModuleColor)])
    stop(
      "hydrated WGCNA GO enrichment contains ModuleID values absent from module_label_table: ",
      paste(missing_ids, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  color_mismatch <- as.character(checked$ModuleColor) != as.character(checked$expected_ModuleColor)
  internal_mismatch <- as.character(checked$WGCNAInternalColor) != as.character(checked$expected_WGCNAInternalColor)
  if (any(color_mismatch, na.rm = TRUE)) {
    stop("hydrated WGCNA GO publication ModuleColor disagrees with module_label_table by ModuleID.", call. = FALSE)
  }
  if (any(internal_mismatch, na.rm = TRUE)) {
    stop("hydrated WGCNA GO WGCNAInternalColor disagrees with module_label_table by ModuleID.", call. = FALSE)
  }
  if (isTRUE(saved_go_evidence) && !nrow(go_enrichment_long)) {
    stop(
      "Cached WGCNA state reports GO enrichment evidence but hydrated GO_enrichment_long has no rows.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_wgcna_hydrated_go_qc <- function(go_enrichment_qc, module_label_table) {
  require_module_contract_columns(
    go_enrichment_qc,
    c(
      "WGCNAInternalColor", "ModuleID", "ModuleColor", "ModuleProteinSetType",
      "Ontology", "status"
    ),
    "hydrated WGCNA GO enrichment QC"
  )
  if (!nrow(go_enrichment_qc)) return(invisible(TRUE))
  qc_ids <- as.character(go_enrichment_qc$ModuleID)
  if (anyNA(qc_ids) || any(!grepl("^WGCNA_m[0-9]{2,}$", qc_ids))) {
    stop("hydrated WGCNA GO enrichment QC must use stable ModuleID values.", call. = FALSE)
  }
  label_index <- match(qc_ids, as.character(module_label_table$ModuleID))
  if (anyNA(label_index)) {
    stop("hydrated WGCNA GO enrichment QC contains ModuleID values absent from module_label_table.", call. = FALSE)
  }
  expected_color <- as.character(module_label_table$ModuleColor[label_index])
  expected_internal <- as.character(module_label_table$WGCNAInternalColor[label_index])
  if (anyNA(go_enrichment_qc$ModuleColor) ||
      any(as.character(go_enrichment_qc$ModuleColor) != expected_color)) {
    stop("hydrated WGCNA GO enrichment QC ModuleColor disagrees with module_label_table by ModuleID.", call. = FALSE)
  }
  if (anyNA(go_enrichment_qc$WGCNAInternalColor) ||
      any(as.character(go_enrichment_qc$WGCNAInternalColor) != expected_internal)) {
    stop("hydrated WGCNA GO enrichment QC WGCNAInternalColor disagrees with module_label_table by ModuleID.", call. = FALSE)
  }
  invisible(TRUE)
}

wgcna_hydrate_cached_state <- function(state, expected_feature_ids = NULL) {
  validate_wgcna_cached_state(state, expected_feature_ids = expected_feature_ids)
  saved_to_runtime <- list(
    expression.data = "expression.data",
    wgcna_feature_table = "wgcna_feature_table",
    wgcna_member_bridge = "wgcna_member_bridge",
    WGCNA_feature_universe = "WGCNA_feature_universe",
    sample_info = "sample_info",
    mergedColors = "mergedColors",
    mergedMEs = "mergedMEs",
    kME = "kME",
    WGCNA_modules_long = "WGCNA_modules_long",
    WGCNA_module_summary = c("WGCNA_module_summary", "module_summary"),
    GO_enrichment_long = c("GO_enrichment_long", "GO_enrichment"),
    GO_enrichment_QC = "GO_enrichment_QC",
    module_name_map = "module_name_map",
    module_label_table = "module_label_table",
    color_to_MEcol = "color_to_MEcol",
    ME_names_stable = "ME_names_stable",
    module_preservation_long = c("module_preservation_long", "module_preservation"),
    geneTree = "geneTree",
    softPower = "softPower",
    parameters = "parameters"
  )
  resolve_saved <- function(saved_names) {
    for (saved_name in saved_names) {
      if (saved_name %in% names(state) && !is.null(state[[saved_name]])) return(state[[saved_name]])
    }
    NULL
  }
  hydrated <- lapply(saved_to_runtime, resolve_saved)
  missing_runtime <- names(hydrated)[vapply(hydrated, is.null, logical(1))]

  qc <- hydrated$GO_enrichment_QC
  qc_reports_evidence <- !is.null(qc) && nrow(qc) && "status" %in% names(qc) &&
    any(tolower(as.character(qc$status)) == "ok", na.rm = TRUE)
  saved_go <- resolve_saved(c("GO_enrichment_long", "GO_enrichment"))
  saved_go_evidence <- (!is.null(saved_go) && nrow(saved_go) > 0L) || qc_reports_evidence
  if ("GO_enrichment_long" %in% missing_runtime && saved_go_evidence) {
    stop(
      "Cached WGCNA state contains GO enrichment evidence but hydrated GO_enrichment_long is NULL.",
      call. = FALSE
    )
  }

  required_runtime <- setdiff(names(saved_to_runtime), "GO_enrichment_long")
  missing_required <- intersect(missing_runtime, required_runtime)
  if (length(missing_required)) {
    stop(
      "Cached WGCNA state cannot reconstruct required downstream runtime objects: ",
      paste(missing_required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  module_color_metadata <- resolve_saved("module_color_metadata")
  if (is.null(module_color_metadata)) {
    metadata_columns <- c(
      "WGCNAInternalColor", "ModuleID", "ModuleLegacyID", "ModuleColor",
      "ModuleColorName", "ModuleColorLabel"
    )
    require_module_contract_columns(
      hydrated$module_label_table,
      metadata_columns,
      "cached WGCNA module label table"
    )
    module_color_metadata <- hydrated$module_label_table |>
      dplyr::select(dplyr::all_of(metadata_columns)) |>
      dplyr::distinct()
  }
  validate_wgcna_module_color_metadata(
    module_color_metadata,
    observed_internal_colors = hydrated$mergedColors,
    artifact = "hydrated WGCNA module colour metadata"
  )
  hydrated$module_color_metadata <- module_color_metadata
  validate_wgcna_hydrated_go_enrichment(
    hydrated$GO_enrichment_long,
    hydrated$module_label_table,
    saved_go_evidence = saved_go_evidence
  )
  validate_wgcna_hydrated_go_qc(
    hydrated$GO_enrichment_QC,
    hydrated$module_label_table
  )
  hydrated
}

wgcna_feature_universe_audit <- function(feature_table, dataset) {
  cls <- as.character(feature_table$protein_group_ambiguity_class)
  included <- feature_table$included_in_wgcna %in% TRUE
  tibble::tibble(
    dataset = dataset,
    total_expression_features = nrow(feature_table),
    features_with_ProteinGroupID = sum(!is.na(feature_table$ProteinGroupID) & nzchar(feature_table$ProteinGroupID)),
    ProteinGroupID_collisions = sum(duplicated(feature_table$ProteinGroupID)),
    single_accession_groups = sum(cls == "single_accession_single_gene"),
    same_gene_multi_accession_groups = sum(cls == "multi_accession_same_gene"),
    multi_gene_groups = sum(cls == "multi_gene_indistinguishable"),
    partially_mapped_groups = sum(cls == "partially_mapped_group"),
    unresolved_groups = sum(cls == "unresolved_group"),
    mixed_species_or_contaminant_groups = sum(cls == "mixed_species_or_contaminant"),
    features_included_in_WGCNA = sum(included),
    features_excluded_from_WGCNA = sum(!included),
    features_eligible_for_gene_level_annotation = sum(included & feature_table$gene_level_claim_allowed),
    features_eligible_for_protein_level_claims = sum(included & feature_table$protein_level_claim_allowed),
    features_excluded_only_from_annotation = sum(included & !feature_table$gene_level_claim_allowed)
  )
}

validate_curated_overlap_programs <- function(df, artifact = "curated overlap programs") {
  require_module_contract_columns(
    df,
    c("ModuleSet", "ModuleID", "UniProt", "GeneSymbol", "Source"),
    artifact
  )
  invisible(TRUE)
}

validate_module_score_output <- function(df, artifact = "module score output") {
  require_module_contract_columns(
    df,
    c(
      "dataset", "module_definition_source", "ModuleID",
      "ModuleScore", "ScoreType", "n_found_in_matrix", "coverage_fraction"
    ),
    artifact
  )
  if (!any(c("Sample", "AnimalID") %in% colnames(df))) {
    stop(
      artifact, " is missing required key column(s): expected Sample or AnimalID.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_wgcna_group_effects <- function(df, artifact = "WGCNA group effects") {
  require_module_contract_columns(
    df,
    c(
      "dataset", "level", "spatial_unit", "contrast", "estimate", "SE",
      "p_value", "FDR_within_dataset_level", "FDR_global", "direction",
      "effect_scope", "SpatialUnitType", "model_type", "has_repeated_animals",
      "n_animals", "n_animals_total", "n_animals_per_group", "min_animals_per_group",
      "n_samples", "n_samples_total", "n_samples_per_group",
      "animal_level_status", "pseudoreplication_guard",
      "formula_requested", "formula_used", "dropped_covariates",
      "model_family", "model_formula", "primary_model_stable", "claim_allowed_model",
      "model_downgrade_reason", "fallback_used", "fallback_type",
      "rank_deficient_model", "singular_model", "emmeans_success",
      "animal_random_effect_used", "biological_replicate_unit", "model_warning"
    ),
    artifact
  )
  if (!any(c("module_id", "supermodule_id") %in% colnames(df))) {
    stop(artifact, " must contain module_id or supermodule_id.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_wgcna_module_annotation <- function(df, artifact = "WGCNA module biological annotation") {
  require_module_contract_columns(
    df,
    c(
      "dataset", "ModuleID", "ModuleLegacyID", "ModuleColor", "ModuleColorName", "ModuleColorLabel",
      "n_proteins", "microenvironment_class",
      "microglia_evidence", "neuropil_evidence", "other_cellular_evidence",
      "canonical_microglia_evidence", "empirical_microglia_roi_evidence",
      "canonical_neuropil_evidence", "empirical_neuropil_evidence",
      "empirical_shared_microenvironment_evidence",
      "microglia_state_or_activation_evidence",
      "peripheral_myeloid_caution_evidence",
      "classification_threshold", "classification_rationale",
      "marker_registry_version", "empirical_marker_set_version",
      "interpretation_note"
    ),
    artifact
  )
  invisible(TRUE)
}

validate_wgcna_interpretable_summary <- function(df, artifact = "WGCNA interpretable summary") {
  require_module_contract_columns(
    df,
    c("dataset", "level", "contrast", "estimate", "p_value", "FDR_global", "interpretation_sentence"),
    artifact
  )
  invisible(TRUE)
}

validate_cross_compartment_program_atlas <- function(df, artifact = "cross-compartment program atlas") {
  require_module_contract_columns(
    df,
    c(
      "dataset", "evidence_domain", "evidence_id", "program_label",
      "entity_type", "entity_id", "source_file", "evidence_status",
      "interpretation_note", "qc_flag", "canonical_claim_entity_id",
      "claim_entity_role", "separate_manuscript_claim_allowed",
      "wgcna_architecture_status", "wgcna_group_effect_status",
      "wgcna_allowed_claim_scope", "wgcna_prohibited_claim_scope",
      "readiness_contract_version", "counts_toward_convergence",
      "evidence_semantic_class"
    ),
    artifact
  )
  invisible(TRUE)
}

validate_manuscript_program_summary <- function(df, artifact = "manuscript program summary") {
  require_module_contract_columns(
    df,
    c(
      "program_key", "manuscript_claim_scope", "datasets_supported",
      "evidence_domains", "strongest_evidence", "safe_manuscript_sentence",
      "main_limitation", "qc_flag", "n_evidence_rows_total",
      "n_evidence_rows_counting_toward_convergence", "n_wgcna_architecture_rows",
      "n_wgcna_stress_effect_rows", "n_wgcna_alias_rows_excluded",
      "claim_semantic_scope"
    ),
    artifact
  )
  invisible(TRUE)
}

validate_evidence_priority_matrix <- function(df, artifact = "evidence priority matrix") {
  require_module_contract_columns(
    df,
    c(
      "priority_id", "program_key", "dataset", "priority_tier",
      "evidence_domain_count", "evidence_rows_counting_toward_convergence",
      "strongest_fdr", "robustness_flag",
      "behavior_flag", "qc_flag", "recommended_use"
    ),
    artifact
  )
  invisible(TRUE)
}

validate_final_evidence_bundle <- function(bundle_dir, artifact = "final biological evidence bundle") {
  required_sheets <- c(
    "README", "input_status", "manuscript_program_summary",
    "evidence_priority_matrix", "cross_compartment_program_atlas",
    "wgcna_key_modules", "wgcna_key_supermodules",
    "microglia_wgcna_claim_readiness",
    "microglia_roi_signature_drivers", "qc_flags", "biological_claims"
  )
  if (!dir.exists(bundle_dir)) {
    stop(artifact, " directory does not exist: ", bundle_dir, call. = FALSE)
  }
  missing <- required_sheets[!file.exists(file.path(bundle_dir, paste0(required_sheets, ".csv")))]
  if (length(missing)) {
    stop(artifact, " is missing CSV sheet mirror(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  read_sheet <- function(sheet) {
    utils::read.csv(file.path(bundle_dir, paste0(sheet, ".csv")), check.names = FALSE, stringsAsFactors = FALSE)
  }

  readme <- read_sheet("README")
  require_module_contract_columns(readme, c("sheet", "produced_from", "meaning", "manuscript_safe_columns"), paste(artifact, "README sheet"))
  missing_readme_rows <- setdiff(required_sheets, as.character(readme$sheet))
  if (length(missing_readme_rows)) {
    stop(artifact, " README sheet does not document: ", paste(missing_readme_rows, collapse = ", "), call. = FALSE)
  }

  input_status <- read_sheet("input_status")
  require_module_contract_columns(input_status, c("input_name", "path", "status", "n_rows"), paste(artifact, "input_status sheet"))

  claims <- read_sheet("biological_claims")
  require_module_contract_columns(
    claims,
    c(
      "claim_id", "dataset", "biological_program", "evidence_type", "claim_type",
      "claim_use_class", "raw_top_GO_term", "representative_GO_terms",
      "semantic_parent_label", "safe_program_label", "term_label_risk",
      "label_confidence", "label_basis", "label_downgrade_reason",
      "claim_grade", "primary_evidence", "orthogonal_support",
      "major_limitation", "safe_interpretation", "unsafe_overinterpretation",
      "claim_allowed", "claim_gate_status", "claim_downgrade_reason",
      "model_fit_status", "statistical_evidence_status", "claim_gate_model_status",
      "primary_model_status", "animal_level_gate", "qc_gate",
      "missingness_gate", "batch_confound_gate", "marker_contamination_gate",
      "microglia_roi_gate", "neuropil_independence_gate", "robustness_gate",
      "evidence_independence_gate",
      "missingness_confounded", "batch_or_plate_confounded",
      "region_layer_imbalance_risk", "animal_pseudoreplication_risk",
      "marker_contamination_or_roi_mixture_flag", "qc_interpretation_flag"
    ),
    paste(artifact, "biological_claims sheet")
  )

  modules <- read_sheet("wgcna_key_modules")
  require_module_contract_columns(
    modules,
    c(
      "dataset", "ModuleID", "ModuleColor", "targeted_signature_primary_driver",
      "targeted_signature_driver_class", "targeted_signature_driver_signature",
      "targeted_signature_driver_padj", "targeted_signature_driver_NES",
      "targeted_signature_driver_overlap_proteins", "canonical_claim_entity_id",
      "claim_entity_role", "separate_manuscript_claim_allowed",
      "primary_architecture_status", "group_effect_status",
      "readiness_contract_version", "wgcna_claim_semantic_status"
    ),
    paste(artifact, "wgcna_key_modules sheet")
  )

  readiness <- read_sheet("microglia_wgcna_claim_readiness")
  require_module_contract_columns(
    readiness,
    c(
      "dataset", "level", "entity_id", "canonical_claim_entity_id",
      "claim_entity_role", "separate_manuscript_claim_allowed",
      "primary_architecture_status", "group_effect_status",
      "allowed_claim_scope", "prohibited_claim_scope",
      "readiness_contract_version", "stage13_source_file"
    ),
    paste(artifact, "microglia_wgcna_claim_readiness sheet")
  )
  if (nrow(readiness) != 22L) stop(artifact, " Stage 13 readiness sheet must contain exactly 22 rows.", call. = FALSE)

  drivers <- read_sheet("microglia_roi_signature_drivers")
  drivers_unavailable <- all(c("status", "message") %in% names(drivers)) &&
    nrow(drivers) == 1L && identical(as.character(drivers$status[[1]]), "missing_optional")
  if (!drivers_unavailable) {
    require_module_contract_columns(
      drivers,
      c(
        "ModuleID", "ModuleColor", "microenvironment_label",
        "targeted_signature_primary_driver", "targeted_signature_driver_class",
        "targeted_signature_driver_signature", "targeted_signature_driver_padj",
        "targeted_signature_driver_NES", "targeted_signature_driver_overlap_proteins"
      ),
      paste(artifact, "microglia_roi_signature_drivers sheet")
    )
  }

  invisible(TRUE)
}

write_contract_validation_status <- function(path, artifact, ok, message = "") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      artifact = artifact,
      ok = isTRUE(ok),
      message = as.character(message),
      stringsAsFactors = FALSE
    ),
    path,
    row.names = FALSE
  )
  invisible(path)
}
