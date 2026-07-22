# Shared calculations for control-only compartment abundance validation.
# All abundance summaries operate on observed, non-imputed log2 values.

ca_group_levels <- function() c("CON", "RES", "SUS")

ca_normalize_group <- function(x, fail_unknown = TRUE) {
  token <- toupper(trimws(as.character(x)))
  token[is.na(x)] <- NA_character_
  out <- rep(NA_character_, length(token))
  out[token %in% c("1", "CON", "CTRL", "CONTROL")] <- "CON"
  out[token %in% c("2", "RES", "RESILIENT")] <- "RES"
  out[token %in% c("3", "SUS", "SUSCEPTIBLE")] <- "SUS"
  bad <- sort(unique(token[!is.na(token) & nzchar(token) & is.na(out)]))
  if (isTRUE(fail_unknown) && length(bad)) {
    stop("Unsupported experimental-group value(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  out
}

ca_resolve_experimental_group_column <- function(metadata) {
  nms <- names(metadata)
  norm <- gsub("[^a-z0-9]", "", tolower(nms))
  preferred <- c("expgroup", "experimentalgroup", "stressgroup", "group")
  for (key in preferred) {
    hits <- which(norm == key)
    if (!length(hits)) next
    for (hit in hits) {
      values <- as.character(metadata[[hit]])
      nonempty <- !is.na(values) & nzchar(trimws(values))
      mapped <- ca_normalize_group(values, fail_unknown = FALSE)
      if (any(nonempty) && all(!nonempty | !is.na(mapped))) return(nms[[hit]])
      if (key == "expgroup") {
        stop(
          "ExpGroup exists but contains values outside the supported 1/2/3 or CON/RES/SUS mapping.",
          call. = FALSE
        )
      }
    }
  }
  stop(
    "No valid experimental-group column was found. ExpGroup is required when an anatomical group column is present.",
    call. = FALSE
  )
}

ca_resolve_metadata_columns <- function(metadata, require_dataset = FALSE, require_group = TRUE,
                                        require_animal = TRUE) {
  choose <- function(lower, upper = NULL, required = TRUE) {
    if (lower %in% names(metadata)) return(lower)
    if (!is.null(upper) && upper %in% names(metadata)) return(upper)
    if (required) stop("Metadata is missing required column: ", lower, call. = FALSE)
    NA_character_
  }
  dataset_col <- if ("dataset" %in% names(metadata)) {
    "dataset"
  } else if ("celltype_layer" %in% names(metadata)) {
    "celltype_layer"
  } else if (isTRUE(require_dataset)) {
    stop("Metadata is missing dataset/celltype_layer.", call. = FALSE)
  } else {
    NA_character_
  }
  list(
    sample = choose("Sample", "sample"),
    group = if (isTRUE(require_group)) ca_resolve_experimental_group_column(metadata) else NA_character_,
    animal = choose("AnimalID", "animal_id", required = require_animal),
    region = choose("region", "Region"),
    layer = choose("layer", "Layer"),
    replicate = choose("ReplicateGroup", "replicate", required = FALSE),
    plate = choose("plate", "Plate", required = FALSE),
    dataset = dataset_col
  )
}

ca_filter_expression_group <- function(mat, metadata, requested_group = "") {
  mat <- as.matrix(mat)
  requested_group <- if (length(requested_group) && !is.na(requested_group[[1]])) {
    trimws(as.character(requested_group[[1]]))
  } else {
    ""
  }
  cols <- ca_resolve_metadata_columns(
    metadata, require_group = nzchar(requested_group), require_animal = nzchar(requested_group)
  )
  sample_ids <- as.character(metadata[[cols$sample]])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("Metadata Sample values must be unique and nonblank.", call. = FALSE)
  }
  idx <- match(colnames(mat), sample_ids)
  if (anyNA(idx) || ncol(mat) != nrow(metadata)) {
    stop("Expression matrix and metadata are not exactly aligned by Sample.", call. = FALSE)
  }
  metadata <- metadata[idx, , drop = FALSE]
  if (!identical(colnames(mat), as.character(metadata[[cols$sample]]))) {
    stop("Expression matrix and metadata order differs after alignment.", call. = FALSE)
  }
  if (!nzchar(requested_group)) {
    return(list(
      mat = mat, meta = metadata, requested_group = "", canonical_group = "",
      resolved_group_column = NA_character_, included_samples = ncol(mat),
      included_animals = if (is.na(cols$animal)) NA_integer_ else length(unique(as.character(metadata[[cols$animal]]))),
      columns = cols
    ))
  }
  canonical <- ca_normalize_group(requested_group)
  observed <- ca_normalize_group(metadata[[cols$group]])
  keep <- !is.na(observed) & observed == canonical
  if (!any(keep)) {
    stop(
      "Requested group '", requested_group, "' (", canonical,
      ") is absent from resolved column '", cols$group, "'.",
      call. = FALSE
    )
  }
  list(
    mat = mat[, keep, drop = FALSE], meta = metadata[keep, , drop = FALSE],
    requested_group = requested_group, canonical_group = canonical,
    resolved_group_column = cols$group, included_samples = sum(keep),
    included_animals = length(unique(as.character(metadata[[cols$animal]][keep]))), columns = cols
  )
}

ca_namespace_paths <- function(paths, canonical_group = "") {
  if (!nzchar(canonical_group)) return(paths)
  lapply(paths, function(path) file.path(path, paste0("group_", canonical_group)))
}

ca_validate_bundle_reconstruction <- function(bundle, tolerance = 1e-12) {
  if (!identical(bundle$contract_version, "joint_compartment_qc_v1")) {
    stop("Unsupported joint-QC bundle contract.", call. = FALSE)
  }
  required_top <- c("metadata", "primary", "broad_union", "feature_table", "feature_filter_audit")
  missing_top <- setdiff(required_top, names(bundle))
  if (length(missing_top)) stop("Joint-QC bundle is missing: ", paste(missing_top, collapse = ", "), call. = FALSE)
  meta <- bundle$metadata
  cols <- ca_resolve_metadata_columns(meta, require_dataset = TRUE)
  sample_ids <- as.character(meta[[cols$sample]])
  feature_table <- bundle$feature_table
  if (!"ProteinGroupID" %in% names(feature_table)) stop("Bundle feature_table lacks ProteinGroupID.", call. = FALSE)
  protein_ids <- as.character(feature_table$ProteinGroupID)
  if (anyNA(protein_ids) || any(!nzchar(protein_ids)) || anyDuplicated(protein_ids)) {
    stop("Bundle feature_table has invalid ProteinGroupID values.", call. = FALSE)
  }
  if (!all(sample_ids %in% names(feature_table))) {
    stop("Bundle feature_table does not retain every raw quantitative sample column.", call. = FALSE)
  }
  raw <- as.matrix(data.frame(lapply(feature_table[sample_ids], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  }), check.names = FALSE))
  rownames(raw) <- protein_ids
  colnames(raw) <- sample_ids
  if (any(is.finite(raw) & raw < 0)) stop("Raw bundle values contain negative intensities.", call. = FALSE)
  raw[is.finite(raw) & raw == 0] <- NA_real_
  raw_log2 <- log2(raw)

  primary_ids <- as.character(bundle$primary$feature_ids)
  primary_raw <- raw_log2[primary_ids, , drop = FALSE]
  if (!identical(dim(primary_raw), dim(bundle$primary$unnormalized_log2)) ||
      !isTRUE(all.equal(primary_raw, bundle$primary$unnormalized_log2, tolerance = 0, check.attributes = TRUE))) {
    stop("Raw-positive log2 reconstruction is not exact for the primary shared core.", call. = FALSE)
  }
  broad_ids <- as.character(bundle$broad_union$feature_ids)
  broad_mask <- (!is.na(raw_log2[broad_ids, , drop = FALSE])) * 1L
  if (!identical(broad_mask, bundle$broad_union$detected_binary)) {
    stop("Raw-observation reconstruction does not reproduce the broad-union detection mask.", call. = FALSE)
  }
  sample_medians <- apply(primary_raw, 2, stats::median, na.rm = TRUE)
  if (any(!is.finite(sample_medians))) stop("A sample has no observed shared-core abundance.", call. = FALSE)
  target_median <- stats::median(sample_medians)
  offsets <- sample_medians - target_median
  normalized_nonimputed <- sweep(raw_log2, 2, offsets, FUN = "-")
  observed_primary <- !is.na(primary_raw)
  reconstructed_primary <- normalized_nonimputed[primary_ids, , drop = FALSE]
  if (!isTRUE(all.equal(
    reconstructed_primary[observed_primary], bundle$primary$matrix[observed_primary],
    tolerance = tolerance, check.attributes = FALSE
  ))) {
    stop("Normalized observed values do not reproduce the primary bundle matrix.", call. = FALSE)
  }
  list(
    metadata = meta,
    columns = cols,
    feature_table = feature_table,
    raw_log2 = raw_log2,
    observed = !is.na(raw_log2),
    normalized_nonimputed = normalized_nonimputed,
    primary_imputation_mask = is.na(primary_raw),
    primary_feature_ids = primary_ids,
    broad_feature_ids = broad_ids,
    normalization_offsets = offsets,
    normalization_target_median = target_median
  )
}

ca_row_summary <- function(x, method) {
  if (method == "mean") {
    out <- rowMeans(x, na.rm = TRUE)
  } else {
    out <- apply(x, 1, stats::median, na.rm = TRUE)
  }
  out[!is.finite(out)] <- NA_real_
  out
}

ca_group_matrix <- function(mat, keys, method = c("mean", "median")) {
  method <- match.arg(method)
  mat <- as.matrix(mat)
  keys <- as.data.frame(keys, stringsAsFactors = FALSE)
  if (nrow(keys) != ncol(mat)) stop("Grouping keys and matrix columns are not aligned.", call. = FALSE)
  key_text <- do.call(paste, c(lapply(keys, as.character), sep = "\r"))
  levels <- unique(key_text)
  group_rows <- match(levels, key_text)
  group_meta <- keys[group_rows, , drop = FALSE]
  values <- matrix(NA_real_, nrow(mat), length(levels), dimnames = list(rownames(mat), levels))
  counts <- matrix(0L, nrow(mat), length(levels), dimnames = list(rownames(mat), levels))
  for (i in seq_along(levels)) {
    keep <- key_text == levels[[i]]
    block <- mat[, keep, drop = FALSE]
    counts[, i] <- rowSums(is.finite(block))
    values[, i] <- ca_row_summary(block, method)
  }
  list(values = values, observed_count = counts, group_meta = group_meta)
}

ca_aggregate_abundance <- function(mat, metadata, method = c("mean", "median"), min_regions = 3L) {
  method <- match.arg(method)
  min_regions <- as.integer(min_regions)
  if (!is.finite(min_regions) || min_regions < 1L || min_regions > 4L) stop("min_regions must be an integer from 1 to 4.", call. = FALSE)
  cols <- ca_resolve_metadata_columns(metadata, require_dataset = TRUE, require_group = FALSE)
  idx <- match(colnames(mat), as.character(metadata[[cols$sample]]))
  if (anyNA(idx) || ncol(mat) != nrow(metadata)) stop("Aggregation matrix and metadata are not aligned.", call. = FALSE)
  metadata <- metadata[idx, , drop = FALSE]
  canonical <- data.frame(
    dataset = as.character(metadata[[cols$dataset]]),
    AnimalID = as.character(metadata[[cols$animal]]),
    region = toupper(as.character(metadata[[cols$region]])),
    layer = tolower(as.character(metadata[[cols$layer]])),
    stringsAsFactors = FALSE
  )
  allowed_regions <- c("CA1", "CA2", "CA3", "DG")
  if (anyNA(canonical) || any(!nzchar(as.matrix(canonical)))) stop("Aggregation metadata contains missing dataset/animal/region/layer values.", call. = FALSE)
  if (length(setdiff(unique(canonical$region), allowed_regions))) stop("Aggregation metadata contains unsupported regions.", call. = FALSE)

  stratum <- ca_group_matrix(mat, canonical, method)
  region_keys <- stratum$group_meta[, c("dataset", "AnimalID", "region"), drop = FALSE]
  region <- ca_group_matrix(stratum$values, region_keys, method)
  region$contributing_layer_count <- region$observed_count
  animal_keys <- region$group_meta[, c("dataset", "AnimalID"), drop = FALSE]
  animal <- ca_group_matrix(region$values, animal_keys, method)
  animal$contributing_region_count <- animal$observed_count
  animal$values[animal$contributing_region_count < min_regions] <- NA_real_
  dataset_keys <- animal$group_meta[, "dataset", drop = FALSE]
  dataset <- ca_group_matrix(animal$values, dataset_keys, method)
  dataset$contributing_animal_count <- dataset$observed_count
  list(
    method = method, min_regions = min_regions, metadata = canonical,
    stratum = stratum, region = region, animal = animal, dataset = dataset
  )
}

ca_matrix_long <- function(values, group_meta, value_name = "abundance", extra_matrices = list()) {
  nr <- nrow(values)
  nc <- ncol(values)
  out <- data.frame(
    ProteinGroupID = rep(rownames(values), times = nc),
    group_meta[rep(seq_len(nc), each = nr), , drop = FALSE],
    value = as.vector(values),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out)[names(out) == "value"] <- value_name
  for (nm in names(extra_matrices)) out[[nm]] <- as.vector(extra_matrices[[nm]])
  rownames(out) <- NULL
  out
}

ca_animal_abundance_table <- function(aggregation, value_name = "animal_log2_abundance") {
  ca_matrix_long(
    aggregation$animal$values, aggregation$animal$group_meta, value_name,
    list(contributing_region_count = aggregation$animal$contributing_region_count)
  )
}

ca_detection_audit <- function(aggregation, min_animals = 2L, strict_animals = 3L) {
  min_animals <- as.integer(min_animals)
  strict_animals <- as.integer(strict_animals)
  datasets <- as.character(aggregation$dataset$group_meta$dataset)
  parts <- lapply(datasets, function(ds) {
    s <- which(aggregation$stratum$group_meta$dataset == ds)
    r <- which(aggregation$region$group_meta$dataset == ds)
    a <- which(aggregation$animal$group_meta$dataset == ds)
    valid_animals <- rowSums(is.finite(aggregation$animal$values[, a, drop = FALSE]))
    observed_animals <- rowSums(aggregation$animal$contributing_region_count[, a, drop = FALSE] > 0L)
    total_samples <- sum(aggregation$metadata$dataset == ds)
    data.frame(
      dataset = ds,
      ProteinGroupID = rownames(aggregation$animal$values),
      observed_sample_count = rowSums(aggregation$stratum$observed_count[, s, drop = FALSE]),
      total_sample_count = total_samples,
      observed_stratum_count = rowSums(is.finite(aggregation$stratum$values[, s, drop = FALSE])),
      total_stratum_count = length(s),
      observed_region_count = rowSums(is.finite(aggregation$region$values[, r, drop = FALSE])),
      total_region_count = length(r),
      observed_animal_count = observed_animals,
      valid_animal_count = valid_animals,
      total_animal_count = length(a),
      observed_animal_fraction = observed_animals / length(a),
      valid_animal_fraction = valid_animals / length(a),
      detection_fraction = valid_animals / length(a),
      primary_eligible = valid_animals >= min_animals,
      strict_eligible = valid_animals >= strict_animals,
      eligibility_reason = ifelse(valid_animals >= min_animals, "passes_primary_valid_animal_count", "fewer_than_primary_valid_animals"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parts)
}

ca_dataset_estimates <- function(aggregation, estimate_name = "region_balanced_log2") {
  ca_matrix_long(
    aggregation$dataset$values, aggregation$dataset$group_meta, estimate_name,
    list(contributing_animal_count = aggregation$dataset$contributing_animal_count)
  )
}

ca_rank_table <- function(estimates, detection_audit, eligibility_col = "primary_eligible",
                          estimate_col = "region_balanced_log2", universe = "primary_2_of_3",
                          rank_group = "region_balanced_global") {
  required <- c("dataset", "ProteinGroupID", estimate_col)
  if (length(setdiff(required, names(estimates)))) stop("Estimate table lacks required columns.", call. = FALSE)
  if (!eligibility_col %in% names(detection_audit)) stop("Detection audit lacks eligibility column.", call. = FALSE)
  merged <- merge(
    estimates, detection_audit,
    by = c("dataset", "ProteinGroupID"), all.x = TRUE, sort = FALSE
  )
  merged <- merged[merged[[eligibility_col]] %in% TRUE & is.finite(merged[[estimate_col]]), , drop = FALSE]
  parts <- lapply(split(merged, merged$dataset), function(x) {
    ord <- order(-x[[estimate_col]], x$ProteinGroupID, method = "radix")
    x <- x[ord, , drop = FALSE]
    n <- nrow(x)
    x$RankGroup <- rank_group
    x$Rank <- seq_len(n)
    x$N <- n
    x$RankPercentile <- if (n == 1L) 1 else 1 - (x$Rank - 1) / (n - 1)
    x$universe <- universe
    x$rank_tie_rule <- "ordinal_desc_abundance_then_ProteinGroupID"
    x
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  if (anyDuplicated(out[, c("dataset", "RankGroup", "ProteinGroupID")])) {
    stop("Base rank table is duplicated by dataset/RankGroup/ProteinGroupID.", call. = FALSE)
  }
  out
}

ca_collapse <- function(x) paste(sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))), collapse = ";")

ca_as_logical <- function(x) {
  token <- toupper(trimws(as.character(x)))
  token %in% c("TRUE", "T", "1", "YES", "Y")
}

ca_match_fidelity_markers <- function(registry, feature_table, shared_feature_ids) {
  required_registry <- c(
    "marker_set", "fidelity_marker_class", "fidelity_subpanel", "gene_symbol",
    "include_in_fidelity_score"
  )
  required_features <- c(
    "ProteinGroupID", "member_gene_symbols", "official_gene_symbol",
    "protein_group_ambiguity_class", "gene_level_claim_allowed", "original_identifier"
  )
  missing <- c(setdiff(required_registry, names(registry)), setdiff(required_features, names(feature_table)))
  if (length(missing)) stop("Marker matching input is missing: ", paste(unique(missing), collapse = ", "), call. = FALSE)
  classes <- c("Soma markers", "Neuropil markers", "Microglia/PVM markers")
  registry <- registry[
    registry$fidelity_marker_class %in% classes & ca_as_logical(registry$include_in_fidelity_score),
    , drop = FALSE
  ]
  registry$marker_gene <- as.character(registry$gene_symbol)
  registry$marker_panel <- as.character(registry$marker_set)
  registry$marker_gene_key <- toupper(trimws(registry$marker_gene))
  registry <- registry[!is.na(registry$marker_gene_key) & nzchar(registry$marker_gene_key), , drop = FALSE]
  registry <- unique(registry)

  member_rows <- lapply(seq_len(nrow(feature_table)), function(i) {
    genes <- trimws(unlist(strsplit(as.character(feature_table$member_gene_symbols[[i]]), ";", fixed = TRUE)))
    genes <- unique(genes[!is.na(genes) & nzchar(genes) & genes != "NA"])
    if (!length(genes)) return(NULL)
    data.frame(
      ProteinGroupID = as.character(feature_table$ProteinGroupID[[i]]),
      matched_official_gene = genes,
      marker_gene_key = toupper(genes),
      stringsAsFactors = FALSE
    )
  })
  member_rows <- Filter(Negate(is.null), member_rows)
  members <- if (length(member_rows)) do.call(rbind, member_rows) else data.frame(
    ProteinGroupID = character(), matched_official_gene = character(), marker_gene_key = character()
  )
  matches <- merge(registry, members, by = "marker_gene_key", all.x = TRUE, sort = FALSE)
  feature_meta <- feature_table[, required_features, drop = FALSE]
  names(feature_meta)[names(feature_meta) == "official_gene_symbol"] <- "feature_official_gene_symbol"
  matches <- merge(matches, feature_meta, by = "ProteinGroupID", all.x = TRUE, sort = FALSE)
  matches$matched <- !is.na(matches$ProteinGroupID) & nzchar(matches$ProteinGroupID)
  matches$in_shared_core <- matches$matched & matches$ProteinGroupID %in% shared_feature_ids
  matches$canonical_marker_eligible <- matches$matched &
    matches$protein_group_ambiguity_class != "mixed_species_or_contaminant" &
    !is.na(matches$original_identifier) & nzchar(trimws(matches$original_identifier))
  class_count <- aggregate(
    fidelity_marker_class ~ ProteinGroupID,
    data = matches[matches$matched, , drop = FALSE],
    FUN = function(x) length(unique(x[!is.na(x) & nzchar(x)]))
  )
  names(class_count)[[2]] <- "n_fidelity_marker_classes"
  matches <- merge(matches, class_count, by = "ProteinGroupID", all.x = TRUE, sort = FALSE)
  matches$n_fidelity_marker_classes[is.na(matches$n_fidelity_marker_classes)] <- 0L
  matches$conflicting_marker_classes <- matches$n_fidelity_marker_classes > 1L
  gene_match_count <- aggregate(
    ProteinGroupID ~ marker_gene_key,
    data = unique(matches[matches$matched, c("marker_gene_key", "ProteinGroupID")]),
    FUN = length
  )
  names(gene_match_count)[[2]] <- "n_protein_groups_for_marker_gene"
  matches <- merge(matches, gene_match_count, by = "marker_gene_key", all.x = TRUE, sort = FALSE)
  matches$n_protein_groups_for_marker_gene[is.na(matches$n_protein_groups_for_marker_gene)] <- 0L
  matches$gene_claim_eligible <- matches$matched & matches$gene_level_claim_allowed %in% TRUE &
    matches$n_protein_groups_for_marker_gene == 1L &
    toupper(as.character(matches$feature_official_gene_symbol)) == matches$marker_gene_key
  matches$primary_score_eligible <- matches$matched & matches$in_shared_core &
    matches$canonical_marker_eligible & !matches$conflicting_marker_classes
  matches$exclusion_reason <- ifelse(
    !matches$matched, "no_canonical_protein_match",
    ifelse(!matches$in_shared_core, "not_in_joint_shared_core",
      ifelse(!matches$canonical_marker_eligible, "canonical_marker_ineligible",
        ifelse(matches$conflicting_marker_classes, "conflicting_fidelity_marker_classes", "included_primary_score")))
  )
  rownames(matches) <- NULL
  matches
}

ca_join_marker_annotations <- function(rank_table, marker_table) {
  if (!nrow(marker_table)) return(rank_table)
  collapsed <- marker_table |>
    dplyr::group_by(.data$ProteinGroupID) |>
    dplyr::summarise(
      marker_classes = ca_collapse(.data$marker_class),
      marker_panels = ca_collapse(.data$marker_panel),
      marker_genes = ca_collapse(.data$marker_gene),
      .groups = "drop"
    )
  out <- dplyr::left_join(rank_table, collapsed, by = "ProteinGroupID", relationship = "many-to-one")
  if (nrow(out) != nrow(rank_table) || anyDuplicated(out[, c("dataset", "RankGroup", "ProteinGroupID")])) {
    stop("Marker annotation duplicated the base rank curve.", call. = FALSE)
  }
  out
}

ca_robust_standardize <- function(animal_abundance) {
  required <- c("ProteinGroupID", "dataset", "AnimalID", "animal_log2_abundance")
  if (length(setdiff(required, names(animal_abundance)))) stop("Animal abundance table lacks robust-z columns.", call. = FALSE)
  parts <- lapply(split(animal_abundance, animal_abundance$ProteinGroupID), function(x) {
    values <- x$animal_log2_abundance
    center <- stats::median(values, na.rm = TRUE)
    raw_mad <- stats::median(abs(values - center), na.rm = TRUE)
    scale <- 1.4826 * raw_mad
    zero <- !is.finite(scale) || scale <= .Machine$double.eps
    x$protein_median_log2 <- center
    x$raw_MAD <- raw_mad
    x$robust_scale <- scale
    x$zero_MAD_excluded <- zero
    x$robust_z <- if (zero) NA_real_ else (values - center) / scale
    x
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

ca_prepare_marker_map <- function(marker_matches) {
  marker_matches |>
    dplyr::filter(.data$primary_score_eligible %in% TRUE, !is.na(.data$marker_class), nzchar(.data$marker_class)) |>
    dplyr::group_by(.data$marker_class, .data$ProteinGroupID) |>
    dplyr::summarise(
      marker_gene = ca_collapse(.data$marker_gene),
      marker_panel = ca_collapse(.data$marker_panel),
      all_fidelity_subpanels = ca_collapse(.data$fidelity_subpanel),
      fidelity_subpanel = {
        z <- sort(unique(as.character(.data$fidelity_subpanel[!is.na(.data$fidelity_subpanel) & nzchar(.data$fidelity_subpanel)])))
        if (length(z)) z[[1]] else "unspecified"
      },
      n_subpanel_memberships = dplyr::n_distinct(.data$fidelity_subpanel[!is.na(.data$fidelity_subpanel) & nzchar(.data$fidelity_subpanel)]),
      .groups = "drop"
    )
}

ca_score_markers <- function(robust_table, marker_map, min_proteins = 3L, min_fraction = 0.5) {
  min_proteins <- as.integer(min_proteins)
  min_fraction <- as.numeric(min_fraction)
  if (min_proteins < 1L || !is.finite(min_fraction) || min_fraction < 0 || min_fraction > 1) {
    stop("Marker-score thresholds are invalid.", call. = FALSE)
  }
  marker_map <- ca_prepare_marker_map(marker_map)
  eligible_counts <- marker_map |>
    dplyr::group_by(.data$marker_class) |>
    dplyr::summarise(n_eligible_proteins = dplyr::n_distinct(.data$ProteinGroupID), .groups = "drop")
  score_grid <- merge(
    unique(robust_table[, c("dataset", "AnimalID"), drop = FALSE]),
    unique(marker_map["marker_class"]),
    by = NULL
  )
  joined <- dplyr::inner_join(robust_table, marker_map, by = "ProteinGroupID", relationship = "many-to-many")
  protein_scores <- joined |>
    dplyr::filter(is.finite(.data$robust_z)) |>
    dplyr::distinct(.data$dataset, .data$AnimalID, .data$marker_class, .data$ProteinGroupID, .keep_all = TRUE)
  observed_scores <- protein_scores |>
    dplyr::group_by(.data$dataset, .data$AnimalID, .data$marker_class) |>
    dplyr::summarise(
      n_contributing_proteins = dplyr::n_distinct(.data$ProteinGroupID),
      mean_robust_z_score = mean(.data$robust_z),
      median_robust_z_score = stats::median(.data$robust_z),
      .groups = "drop"
    )
  class_scores <- score_grid |>
    dplyr::left_join(observed_scores, by = c("dataset", "AnimalID", "marker_class")) |>
    dplyr::left_join(eligible_counts, by = "marker_class") |>
    dplyr::mutate(
      n_contributing_proteins = dplyr::coalesce(.data$n_contributing_proteins, 0L),
      fraction_contributing_proteins = .data$n_contributing_proteins / .data$n_eligible_proteins,
      score_eligible = .data$n_contributing_proteins >= min_proteins & .data$fraction_contributing_proteins >= min_fraction,
      mean_robust_z_score = ifelse(.data$score_eligible, .data$mean_robust_z_score, NA_real_),
      median_robust_z_score = ifelse(.data$score_eligible, .data$median_robust_z_score, NA_real_),
      precision_interpretation = "descriptive_only_marker_classes_have_unequal_protein_counts"
    )
  subpanel <- protein_scores |>
    dplyr::group_by(.data$dataset, .data$AnimalID, .data$marker_class, .data$fidelity_subpanel) |>
    dplyr::summarise(
      subpanel_mean_robust_z = mean(.data$robust_z),
      n_subpanel_proteins = dplyr::n_distinct(.data$ProteinGroupID),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$dataset, .data$AnimalID, .data$marker_class) |>
    dplyr::summarise(
      subpanel_balanced_mean_robust_z = mean(.data$subpanel_mean_robust_z),
      n_contributing_subpanels = dplyr::n_distinct(.data$fidelity_subpanel),
      .groups = "drop"
    )
  class_scores <- dplyr::left_join(class_scores, subpanel, by = c("dataset", "AnimalID", "marker_class"))
  list(scores = class_scores, protein_contributions = protein_scores, marker_map = marker_map)
}

ca_exact_paired_bootstrap_draws <- function(animals) {
  animals <- sort(unique(as.character(animals)))
  if (length(animals) != 3L) stop("Exact paired bootstrap requires exactly three animals.", call. = FALSE)
  grid <- expand.grid(draw_1 = animals, draw_2 = animals, draw_3 = animals, stringsAsFactors = FALSE)
  grid$bootstrap_draw_id <- sprintf("B%02d", seq_len(nrow(grid)))
  grid <- grid[, c("bootstrap_draw_id", "draw_1", "draw_2", "draw_3")]
  if (nrow(grid) != 27L || anyDuplicated(grid[, c("draw_1", "draw_2", "draw_3")])) {
    stop("Exact paired bootstrap did not produce 27 unique ordered draws.", call. = FALSE)
  }
  grid
}

ca_leave_one_animal_out_cases <- function(animals) {
  animals <- sort(unique(as.character(animals)))
  if (length(animals) != 3L) stop("Leave-one-animal-out requires exactly three animals.", call. = FALSE)
  data.frame(omission_id = sprintf("LOAO_%s", animals), omitted_AnimalID = animals, stringsAsFactors = FALSE)
}

ca_rank_concordance <- function(a, b, label_a, label_b) {
  joined <- merge(
    a[, c("dataset", "ProteinGroupID", "Rank", "RankPercentile")],
    b[, c("dataset", "ProteinGroupID", "Rank", "RankPercentile")],
    by = c("dataset", "ProteinGroupID"), suffixes = c("_a", "_b")
  )
  summary <- do.call(rbind, lapply(split(joined, joined$dataset), function(x) data.frame(
    dataset = x$dataset[[1]], comparison = paste(label_a, "vs", label_b),
    n_common_proteins = nrow(x),
    spearman_rank = stats::cor(x$Rank_a, x$Rank_b, method = "spearman"),
    spearman_rank_percentile = stats::cor(x$RankPercentile_a, x$RankPercentile_b, method = "spearman"),
    stringsAsFactors = FALSE
  )))
  list(by_protein = joined, summary = summary)
}

ca_save_vector_figure <- function(plot, root, stem, width_mm, height_mm) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  svg <- file.path(root, paste0(stem, ".svg"))
  pdf <- file.path(root, paste0(stem, ".pdf"))
  ggplot2::ggsave(svg, plot, width = width_mm, height = height_mm, units = "mm",
                  device = function(...) svglite::svglite(..., pointsize = 7, fix_text_size = FALSE),
                  limitsize = FALSE, bg = "white")
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  args <- list(filename = pdf, plot = plot, width = width_mm, height = height_mm,
               units = "mm", device = pdf_device, limitsize = FALSE, bg = "white")
  if (identical(pdf_device, grDevices::pdf)) args$useDingbats <- FALSE
  do.call(ggplot2::ggsave, args)
  c(svg = svg, pdf = pdf)
}
