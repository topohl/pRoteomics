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
    hemisphere = ca_resolve_hemisphere_column(metadata),
    plate = choose("plate", "Plate", required = FALSE),
    dataset = dataset_col
  )
}

ca_normalize_hemisphere <- function(x, fail_unknown = TRUE) {
  token <- toupper(trimws(as.character(x)))
  token[is.na(x)] <- NA_character_
  out <- rep(NA_character_, length(token))
  out[token %in% c("L", "LEFT")] <- "Left"
  out[token %in% c("R", "RIGHT")] <- "Right"
  bad <- sort(unique(token[!is.na(token) & nzchar(token) & is.na(out)]))
  if (isTRUE(fail_unknown) && length(bad)) {
    stop("Unsupported hemisphere value(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  out
}

ca_resolve_hemisphere_column <- function(metadata) {
  candidates <- intersect(
    c("hemisphere", "Hemisphere", "side", "Side", "ReplicateGroup", "replicate"),
    names(metadata)
  )
  for (candidate in candidates) {
    values <- as.character(metadata[[candidate]])
    nonempty <- !is.na(values) & nzchar(trimws(values))
    normalized <- ca_normalize_hemisphere(values, fail_unknown = FALSE)
    if (any(nonempty) && all(!nonempty | !is.na(normalized))) return(candidate)
  }
  NA_character_
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
  optional_features <- intersect(
    c("joint_qc_eligible", "joint_qc_exclusion_reason"),
    names(feature_table)
  )
  feature_meta <- feature_table[, unique(c(required_features, optional_features)), drop = FALSE]
  names(feature_meta)[names(feature_meta) == "official_gene_symbol"] <- "feature_official_gene_symbol"
  matches <- merge(matches, feature_meta, by = "ProteinGroupID", all.x = TRUE, sort = FALSE)
  matches$matched <- !is.na(matches$ProteinGroupID) & nzchar(matches$ProteinGroupID)
  matches$in_shared_core <- matches$matched & matches$ProteinGroupID %in% shared_feature_ids
  if (!"joint_qc_eligible" %in% names(matches)) {
    matches$joint_qc_eligible <- matches$matched &
      matches$protein_group_ambiguity_class != "mixed_species_or_contaminant"
  }
  if (!"joint_qc_exclusion_reason" %in% names(matches)) {
    matches$joint_qc_exclusion_reason <- ifelse(
      matches$joint_qc_eligible %in% TRUE, NA_character_, "mixed_species_or_contaminant"
    )
  }
  matches$canonical_marker_eligible <- matches$matched &
    matches$joint_qc_eligible %in% TRUE &
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
  matches$primary_score_eligible <- matches$matched &
    matches$canonical_marker_eligible & !matches$conflicting_marker_classes
  matches$exclusion_reason <- ifelse(
    !matches$matched, "no_canonical_protein_match",
    ifelse(!matches$canonical_marker_eligible, "canonical_marker_ineligible",
      ifelse(matches$conflicting_marker_classes, "conflicting_fidelity_marker_classes",
        "canonical_marker_eligible_detection_pending"))
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

# -----------------------------------------------------------------------------
# Control compartment-abundance v2 contracts
# -----------------------------------------------------------------------------

ca_aggregate_abundance_v2 <- function(mat, metadata, method = c("median", "mean"),
                                      min_regions = 3L, require_both_hemispheres = FALSE) {
  method <- match.arg(method)
  min_regions <- as.integer(min_regions)
  if (!is.finite(min_regions) || min_regions < 1L || min_regions > 4L) {
    stop("min_regions must be an integer from 1 to 4.", call. = FALSE)
  }
  cols <- ca_resolve_metadata_columns(metadata, require_dataset = TRUE, require_group = FALSE)
  idx <- match(colnames(mat), as.character(metadata[[cols$sample]]))
  if (anyNA(idx) || ncol(mat) != nrow(metadata)) {
    stop("Aggregation matrix and metadata are not aligned.", call. = FALSE)
  }
  metadata <- metadata[idx, , drop = FALSE]
  has_hemisphere <- !is.na(cols$hemisphere)
  hemisphere <- if (has_hemisphere) {
    ca_normalize_hemisphere(metadata[[cols$hemisphere]])
  } else {
    rep("not_available", nrow(metadata))
  }
  canonical <- data.frame(
    dataset = as.character(metadata[[cols$dataset]]),
    AnimalID = as.character(metadata[[cols$animal]]),
    hemisphere = hemisphere,
    region = toupper(as.character(metadata[[cols$region]])),
    layer = tolower(as.character(metadata[[cols$layer]])),
    stringsAsFactors = FALSE
  )
  allowed_regions <- c("CA1", "CA2", "CA3", "DG")
  if (anyNA(canonical) || any(!nzchar(as.matrix(canonical)))) {
    stop("Aggregation metadata contains missing dataset/animal/hemisphere/region/layer values.", call. = FALSE)
  }
  if (length(setdiff(unique(canonical$region), allowed_regions))) {
    stop("Aggregation metadata contains unsupported regions.", call. = FALSE)
  }

  technical <- ca_group_matrix(mat, canonical, method)
  region_keys <- technical$group_meta[, c("dataset", "AnimalID", "hemisphere", "region"), drop = FALSE]
  region <- ca_group_matrix(technical$values, region_keys, method)
  region$contributing_layer_count <- region$observed_count
  hemisphere_keys <- region$group_meta[, c("dataset", "AnimalID", "hemisphere"), drop = FALSE]
  hemisphere_level <- ca_group_matrix(region$values, hemisphere_keys, method)
  hemisphere_level$contributing_region_count <- hemisphere_level$observed_count
  hemisphere_level$values[hemisphere_level$contributing_region_count < min_regions] <- NA_real_

  animal_keys <- hemisphere_level$group_meta[, c("dataset", "AnimalID"), drop = FALSE]
  animal <- ca_group_matrix(hemisphere_level$values, animal_keys, method)
  animal$contributing_hemisphere_count <- animal$observed_count
  if (isTRUE(require_both_hemispheres) && has_hemisphere) {
    animal$values[animal$contributing_hemisphere_count < 2L] <- NA_real_
  }
  dataset_keys <- animal$group_meta[, "dataset", drop = FALSE]
  dataset <- ca_group_matrix(animal$values, dataset_keys, method)
  dataset$contributing_animal_count <- dataset$observed_count

  list(
    method = method,
    min_regions = min_regions,
    require_both_hemispheres = isTRUE(require_both_hemispheres),
    hemisphere_column = if (has_hemisphere) cols$hemisphere else NA_character_,
    hemisphere_handling = if (has_hemisphere) {
      paste0(
        "technical_replicates_then_layers_then_regions_within_hemisphere;",
        "hemisphere_valid_if_at_least_", min_regions,
        "_regions;valid_hemispheres_equal_weight_within_animal;",
        if (isTRUE(require_both_hemispheres)) "two_valid_hemispheres_required" else "one_or_two_valid_hemispheres_allowed"
      )
    } else {
      "hemisphere_field_not_available;legacy_region_then_animal_hierarchy"
    },
    metadata = canonical,
    technical = technical,
    stratum = technical,
    region = region,
    hemisphere = hemisphere_level,
    animal = animal,
    dataset = dataset
  )
}

ca_hemisphere_abundance_table <- function(aggregation, value_name = "hemisphere_log2_abundance") {
  ca_matrix_long(
    aggregation$hemisphere$values,
    aggregation$hemisphere$group_meta,
    value_name,
    list(contributing_region_count = aggregation$hemisphere$contributing_region_count)
  )
}

ca_animal_abundance_table_v2 <- function(aggregation, value_name = "animal_log2_abundance") {
  ca_matrix_long(
    aggregation$animal$values,
    aggregation$animal$group_meta,
    value_name,
    list(contributing_hemisphere_count = aggregation$animal$contributing_hemisphere_count)
  )
}

ca_detection_audit_v2 <- function(aggregation, min_animals = 2L, strict_animals = 3L) {
  min_animals <- as.integer(min_animals)
  strict_animals <- as.integer(strict_animals)
  datasets <- as.character(aggregation$dataset$group_meta$dataset)
  parts <- lapply(datasets, function(ds) {
    s <- which(aggregation$technical$group_meta$dataset == ds)
    r <- which(aggregation$region$group_meta$dataset == ds)
    h <- which(aggregation$hemisphere$group_meta$dataset == ds)
    a <- which(aggregation$animal$group_meta$dataset == ds)
    animal_values <- aggregation$animal$values[, a, drop = FALSE]
    valid_animals <- rowSums(is.finite(animal_values))
    observed_animals <- rowSums(
      aggregation$animal$contributing_hemisphere_count[, a, drop = FALSE] > 0L
    )
    animal_median <- apply(animal_values, 1, stats::median, na.rm = TRUE)
    animal_mean <- rowMeans(animal_values, na.rm = TRUE)
    animal_median[!is.finite(animal_median)] <- NA_real_
    animal_mean[!is.finite(animal_mean)] <- NA_real_
    total_samples <- sum(aggregation$metadata$dataset == ds)
    data.frame(
      dataset = ds,
      ProteinGroupID = rownames(aggregation$animal$values),
      observed_sample_count = rowSums(aggregation$technical$observed_count[, s, drop = FALSE]),
      total_sample_count = total_samples,
      observed_stratum_count = rowSums(is.finite(aggregation$technical$values[, s, drop = FALSE])),
      total_stratum_count = length(s),
      observed_region_count = rowSums(is.finite(aggregation$region$values[, r, drop = FALSE])),
      total_region_count = length(r),
      observed_hemisphere_count = rowSums(
        aggregation$hemisphere$contributing_region_count[, h, drop = FALSE] > 0L
      ),
      valid_hemisphere_count = rowSums(
        is.finite(aggregation$hemisphere$values[, h, drop = FALSE])
      ),
      total_hemisphere_count = length(h),
      observed_animal_count = observed_animals,
      valid_animal_count = valid_animals,
      total_animal_count = length(a),
      observed_animal_fraction = observed_animals / length(a),
      valid_animal_fraction = valid_animals / length(a),
      median_animal_log2_abundance = animal_median,
      mean_animal_log2_abundance = animal_mean,
      primary_eligible = valid_animals >= min_animals,
      strict_eligible = valid_animals >= strict_animals,
      reliability_status = ifelse(
        valid_animals >= min_animals,
        "reliably_detected",
        "not_reliably_detected"
      ),
      animal_detection_threshold = paste0(">=", min_animals, "_valid_animals"),
      region_threshold = paste0(">=", aggregation$min_regions, "_regions_per_valid_hemisphere"),
      hemisphere_handling = aggregation$hemisphere_handling,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

ca_standardize_marker_registry_ranks <- function(registry) {
  registry <- as.data.frame(registry, stringsAsFactors = FALSE, check.names = FALSE)
  registry$.registry_row <- seq_len(nrow(registry))
  if (!"source_rank" %in% names(registry)) registry$source_rank <- NA_integer_
  if (!"source_rank_scope" %in% names(registry)) registry$source_rank_scope <- NA_character_
  if (!"source_ranking_method" %in% names(registry)) registry$source_ranking_method <- NA_character_
  registry$source_rank <- suppressWarnings(as.integer(registry$source_rank))
  eligible <- ca_as_logical(registry$include_in_fidelity_score)

  assign_rank <- function(rows, order_rows, scope, method) {
    rows <- rows[is.na(registry$source_rank[rows])]
    if (!length(rows)) return(invisible(NULL))
    ordered <- intersect(order_rows, rows)
    registry$source_rank[ordered] <<- seq_along(ordered)
    registry$source_rank_scope[ordered] <<- scope
    registry$source_ranking_method[ordered] <<- method
  }

  syngo_rows <- which(eligible & registry$source_name == "SynGO")
  syngo_groups <- split(
    syngo_rows,
    paste(
      registry$fidelity_marker_class[syngo_rows],
      registry$fidelity_subpanel[syngo_rows],
      sep = "\r"
    )
  )
  for (rows in syngo_groups) {
    evidence_priority <- ifelse(registry$evidence_level[rows] == "traceable_experimental", 1L, 2L)
    ord <- rows[order(
      evidence_priority,
      toupper(registry$gene_symbol[rows]),
      registry$.registry_row[rows],
      method = "radix"
    )]
    assign_rank(
      rows, ord,
      paste0("fidelity_subpanel:", registry$fidelity_subpanel[rows[[1]]]),
      "rank_within_subpanel_by_SynGO_evidence_priority_then_gene_symbol"
    )
  }

  go_priority <- c("EXP", "IDA", "IMP", "IPI", "IGI")
  go_rows <- which(eligible & registry$source_name == "GO_MGI")
  go_groups <- split(
    go_rows,
    paste(
      registry$fidelity_marker_class[go_rows],
      registry$fidelity_subpanel[go_rows],
      sep = "\r"
    )
  )
  for (rows in go_groups) {
    evidence_priority <- match(registry$evidence_level[rows], go_priority)
    evidence_priority[is.na(evidence_priority)] <- length(go_priority) + 1L
    ord <- rows[order(
      evidence_priority,
      registry$source_term_or_category[rows],
      toupper(registry$gene_symbol[rows]),
      registry$.registry_row[rows],
      method = "radix"
    )]
    assign_rank(
      rows, ord,
      paste0("fidelity_subpanel:", registry$fidelity_subpanel[rows[[1]]]),
      "rank_after_documented_GO_MGI_evidence_priority_then_term_then_gene"
    )
  }

  allen_rows <- which(
    eligible & grepl("allen", registry$source_name, ignore.case = TRUE)
  )
  allen_groups <- split(
    allen_rows,
    paste(
      registry$source_name[allen_rows],
      registry$fidelity_subpanel[allen_rows],
      sep = "\r"
    )
  )
  for (rows in allen_groups) {
    ord <- rows[order(registry$.registry_row[rows], method = "radix")]
    assign_rank(
      rows, ord,
      paste0("source:", registry$source_name[rows[[1]]], ";subpanel:", registry$fidelity_subpanel[rows[[1]]]),
      "preserved_04b_registry_order_from_external_rank_within_source"
    )
  }
  registry$.registry_row <- NULL
  registry
}

ca_intended_dataset <- function(marker_class) {
  out <- c(
    "Soma markers" = "neuron_soma",
    "Neuropil markers" = "neuron_neuropil",
    "Microglia/PVM markers" = "microglia"
  )[as.character(marker_class)]
  unname(out)
}

ca_display_biological_subpanel <- function(marker_class, fidelity_subpanel) {
  x <- tolower(as.character(fidelity_subpanel))
  cls <- as.character(marker_class)
  out <- x
  soma <- cls == "Soma markers"
  neuropil <- cls == "Neuropil markers"
  microglia <- cls == "Microglia/PVM markers"
  out[soma & grepl("chromatin|nuclear_organization", x)] <- "chromatin/nuclear organization"
  out[soma & grepl("nuclear_matrix", x)] <- "nuclear matrix"
  out[soma & grepl("nuclear_speck", x)] <- "nuclear speck"
  out[soma & grepl("spliceosome", x)] <- "spliceosome"
  out[soma & grepl("rnp|rna", x)] <- "RNP/RNA processing"
  out[neuropil & grepl("presynaptic_vesicle|synaptic_vesicle", x)] <- "presynaptic vesicle"
  out[neuropil & grepl("active_zone|presynapse", x)] <- "active zone/presynapse"
  out[neuropil & grepl("postsynaptic_density", x)] <- "postsynaptic density"
  out[neuropil & grepl("postsynapse", x) & !grepl("density", x)] <- "postsynapse"
  out[neuropil & grepl("dendritic_spine|projection", x)] <- "dendritic spine/projection"
  out[microglia] <- "Allen microglia/PVM"
  out[is.na(out) | !nzchar(out)] <- "unspecified"
  out
}

ca_marker_subpanel_order <- function(marker_class, subpanel) {
  priorities <- list(
    "Soma markers" = c(
      "chromatin/nuclear organization", "nuclear matrix", "nuclear speck",
      "spliceosome", "RNP/RNA processing"
    ),
    "Neuropil markers" = c(
      "presynaptic vesicle", "active zone/presynapse", "postsynapse",
      "postsynaptic density", "dendritic spine/projection"
    ),
    "Microglia/PVM markers" = "Allen microglia/PVM"
  )
  vapply(seq_along(marker_class), function(i) {
    p <- priorities[[as.character(marker_class[[i]])]]
    hit <- match(as.character(subpanel[[i]]), p)
    ifelse(is.na(hit), length(p) + 1L, hit)
  }, integer(1))
}

ca_apply_marker_eligibility_v2 <- function(marker_matches, detection_audit,
                                           shared_feature_ids, broad_feature_ids) {
  matches <- as.data.frame(marker_matches, stringsAsFactors = FALSE, check.names = FALSE)
  matches$intended_dataset <- ca_intended_dataset(matches$fidelity_marker_class)
  matches$in_joint_shared_core <- matches$matched &
    matches$ProteinGroupID %in% shared_feature_ids
  matches$in_broad_union <- matches$matched &
    matches$ProteinGroupID %in% broad_feature_ids
  intended_detection <- detection_audit[, c(
    "dataset", "ProteinGroupID", "valid_animal_count", "valid_animal_fraction",
    "primary_eligible", "strict_eligible"
  )]
  names(intended_detection) <- c(
    "intended_dataset", "ProteinGroupID", "intended_valid_animal_count",
    "intended_valid_animal_fraction", "intended_detection_primary",
    "intended_detection_strict"
  )
  matches <- dplyr::left_join(
    matches,
    intended_detection,
    by = c("intended_dataset", "ProteinGroupID"),
    relationship = "many-to-one"
  )
  base_eligible <- matches$matched %in% TRUE &
    matches$canonical_marker_eligible %in% TRUE &
    !(matches$conflicting_marker_classes %in% TRUE)
  matches$intended_primary_eligible <- base_eligible &
    matches$intended_detection_primary %in% TRUE
  matches$intended_strict_eligible <- base_eligible &
    matches$intended_detection_strict %in% TRUE
  matches$shared_core_sensitivity_eligible <- matches$intended_primary_eligible &
    matches$in_joint_shared_core
  matches$primary_score_eligible <- matches$intended_primary_eligible
  matches$exclusion_reason_primary <- ifelse(
    !matches$matched, "no_canonical_protein_match",
    ifelse(
      !(matches$canonical_marker_eligible %in% TRUE),
      "canonical_marker_ineligible",
      ifelse(
        matches$conflicting_marker_classes %in% TRUE,
        "conflicting_fidelity_marker_classes",
        ifelse(
          !(matches$intended_detection_primary %in% TRUE),
          "intended_compartment_fewer_than_2_valid_CON_animals",
          "included_primary_intended_detection"
        )
      )
    )
  )
  matches$exclusion_reason_shared_core_sensitivity <- ifelse(
    !matches$intended_primary_eligible,
    matches$exclusion_reason_primary,
    ifelse(
      !matches$in_joint_shared_core,
      "not_in_joint_shared_core",
      "included_shared_core_sensitivity"
    )
  )
  matches
}

ca_prepare_marker_map_v2 <- function(marker_matches, eligibility_col = NULL) {
  x <- as.data.frame(marker_matches, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- x$matched %in% TRUE & x$canonical_marker_eligible %in% TRUE &
    !(x$conflicting_marker_classes %in% TRUE)
  if (!is.null(eligibility_col)) {
    if (!eligibility_col %in% names(x)) stop("Unknown marker eligibility column.", call. = FALSE)
    keep <- keep & x[[eligibility_col]] %in% TRUE
  }
  x <- x[keep, , drop = FALSE]
  if (!nrow(x)) {
    return(data.frame(
      marker_class = character(), ProteinGroupID = character(),
      marker_gene = character(), marker_panel = character(),
      fidelity_subpanel = character(), display_biological_subpanel = character(),
      intended_dataset = character(), source_rank = integer(),
      source_rank_scope = character(), source_ranking_method = character(),
      intended_primary_eligible = logical(), intended_strict_eligible = logical(),
      shared_core_sensitivity_eligible = logical(), gene_claim_eligible = logical(),
      stringsAsFactors = FALSE
    ))
  }
  x$display_biological_subpanel <- ca_display_biological_subpanel(
    x$fidelity_marker_class, x$fidelity_subpanel
  )
  x$subpanel_priority <- ca_marker_subpanel_order(
    x$fidelity_marker_class, x$display_biological_subpanel
  )
  rank_value <- suppressWarnings(as.numeric(x$source_rank))
  rank_value[!is.finite(rank_value)] <- Inf
  x <- x[order(
    match(x$fidelity_marker_class, c(
      "Soma markers", "Neuropil markers", "Microglia/PVM markers"
    )),
    x$ProteinGroupID,
    x$subpanel_priority,
    rank_value,
    toupper(x$marker_gene),
    method = "radix"
  ), , drop = FALSE]
  first_rows <- !duplicated(x[, c("fidelity_marker_class", "ProteinGroupID")])
  representative <- x[first_rows, c(
    "fidelity_marker_class", "ProteinGroupID", "marker_gene", "marker_panel",
    "fidelity_subpanel", "display_biological_subpanel", "intended_dataset",
    "source_rank", "source_rank_scope", "source_ranking_method", "source_name",
    "source_reference", "source_term_or_category", "evidence_level",
    "selection_rule", "confidence", "gene_claim_eligible",
    "intended_primary_eligible", "intended_strict_eligible",
    "shared_core_sensitivity_eligible", "in_joint_shared_core", "in_broad_union"
  ), drop = FALSE]
  names(representative)[names(representative) == "fidelity_marker_class"] <- "marker_class"
  collapsed <- x |>
    dplyr::group_by(.data$fidelity_marker_class, .data$ProteinGroupID) |>
    dplyr::summarise(
      all_marker_genes = ca_collapse(.data$marker_gene),
      all_marker_panels = ca_collapse(.data$marker_panel),
      all_fidelity_subpanels = ca_collapse(.data$fidelity_subpanel),
      all_source_names = ca_collapse(.data$source_name),
      .groups = "drop"
    )
  names(collapsed)[names(collapsed) == "fidelity_marker_class"] <- "marker_class"
  out <- dplyr::left_join(
    representative,
    collapsed,
    by = c("marker_class", "ProteinGroupID"),
    relationship = "one-to-one"
  )
  if (anyDuplicated(out[c("marker_class", "ProteinGroupID")])) {
    stop("Marker map is duplicated by marker_class and ProteinGroupID.", call. = FALSE)
  }
  out
}

ca_center_animal_abundance <- function(animal_abundance) {
  required <- c("ProteinGroupID", "dataset", "AnimalID", "animal_log2_abundance")
  if (length(setdiff(required, names(animal_abundance)))) {
    stop("Animal abundance table lacks centered-abundance columns.", call. = FALSE)
  }
  animal_abundance |>
    dplyr::group_by(.data$ProteinGroupID) |>
    dplyr::mutate(
      protein_center_log2 = {
        value <- stats::median(.data$animal_log2_abundance, na.rm = TRUE)
        ifelse(is.finite(value), value, NA_real_)
      },
      centered_log2 = .data$animal_log2_abundance - .data$protein_center_log2
    ) |>
    dplyr::ungroup()
}

ca_score_markers_v2 <- function(centered_abundance, marker_map,
                                min_proteins = 3L, min_fraction = 0.5) {
  min_proteins <- as.integer(min_proteins)
  min_fraction <- as.numeric(min_fraction)
  if (min_proteins < 1L || !is.finite(min_fraction) ||
      min_fraction < 0 || min_fraction > 1) {
    stop("Marker-score thresholds are invalid.", call. = FALSE)
  }
  if (anyDuplicated(marker_map[c("marker_class", "ProteinGroupID")])) {
    stop("Subpanel balancing requires one marker-map row per class and ProteinGroupID.", call. = FALSE)
  }
  units <- unique(centered_abundance[c("dataset", "AnimalID")])
  classes <- c("Soma markers", "Neuropil markers", "Microglia/PVM markers")
  grid <- merge(units, data.frame(marker_class = classes), by = NULL)
  eligible <- marker_map |>
    dplyr::group_by(.data$marker_class) |>
    dplyr::summarise(
      n_eligible_proteins = dplyr::n_distinct(.data$ProteinGroupID),
      n_eligible_subpanels = dplyr::n_distinct(.data$display_biological_subpanel),
      .groups = "drop"
    )
  contributions <- dplyr::inner_join(
    centered_abundance,
    marker_map,
    by = "ProteinGroupID",
    relationship = "many-to-many"
  ) |>
    dplyr::filter(is.finite(.data$centered_log2)) |>
    dplyr::distinct(
      .data$dataset, .data$AnimalID, .data$marker_class,
      .data$ProteinGroupID, .keep_all = TRUE
    )
  subpanel_scores <- contributions |>
    dplyr::group_by(
      .data$dataset, .data$AnimalID, .data$marker_class,
      .data$display_biological_subpanel
    ) |>
    dplyr::summarise(
      subpanel_median_centered_log2 = stats::median(.data$centered_log2),
      n_subpanel_contributing_proteins = dplyr::n_distinct(.data$ProteinGroupID),
      .groups = "drop"
    )
  observed <- contributions |>
    dplyr::group_by(.data$dataset, .data$AnimalID, .data$marker_class) |>
    dplyr::summarise(
      median_centered_log2 = stats::median(.data$centered_log2),
      n_contributing_proteins = dplyr::n_distinct(.data$ProteinGroupID),
      .groups = "drop"
    )
  balanced <- subpanel_scores |>
    dplyr::group_by(.data$dataset, .data$AnimalID, .data$marker_class) |>
    dplyr::summarise(
      subpanel_balanced_median_centered_log2 =
        stats::median(.data$subpanel_median_centered_log2),
      n_contributing_subpanels =
        dplyr::n_distinct(.data$display_biological_subpanel),
      .groups = "drop"
    )
  scores <- grid |>
    dplyr::left_join(observed, by = c("dataset", "AnimalID", "marker_class")) |>
    dplyr::left_join(balanced, by = c("dataset", "AnimalID", "marker_class")) |>
    dplyr::left_join(eligible, by = "marker_class") |>
    dplyr::mutate(
      n_eligible_proteins = dplyr::coalesce(.data$n_eligible_proteins, 0L),
      n_eligible_subpanels = dplyr::coalesce(.data$n_eligible_subpanels, 0L),
      n_contributing_proteins = dplyr::coalesce(.data$n_contributing_proteins, 0L),
      n_contributing_subpanels = dplyr::coalesce(.data$n_contributing_subpanels, 0L),
      fraction_contributing = ifelse(
        .data$n_eligible_proteins > 0L,
        .data$n_contributing_proteins / .data$n_eligible_proteins,
        NA_real_
      ),
      score_eligible = .data$n_contributing_proteins >= min_proteins &
        .data$fraction_contributing >= min_fraction,
      median_centered_log2 = ifelse(
        .data$score_eligible, .data$median_centered_log2, NA_real_
      ),
      subpanel_balanced_median_centered_log2 = ifelse(
        .data$score_eligible,
        .data$subpanel_balanced_median_centered_log2,
        NA_real_
      ),
      primary_statistic = "median_within_protein_centered_log2",
      class_balance_statistic = "median_of_subpanel_medians"
    )
  list(
    scores = scores,
    subpanel_scores = subpanel_scores,
    protein_contributions = contributions
  )
}

ca_marker_dataset_audit_v2 <- function(marker_map, detection_audit) {
  datasets <- c("neuron_soma", "neuron_neuropil", "microglia")
  if (!nrow(marker_map)) {
    return(data.frame(
      marker_class = character(), ProteinGroupID = character(),
      dataset = character(), stringsAsFactors = FALSE
    ))
  }
  grid <- merge(marker_map, data.frame(dataset = datasets), by = NULL)
  detection_columns <- c(
    "dataset", "ProteinGroupID", "valid_animal_count", "valid_animal_fraction",
    "observed_animal_count", "observed_sample_count", "observed_region_count",
    "observed_hemisphere_count", "valid_hemisphere_count",
    "median_animal_log2_abundance", "mean_animal_log2_abundance",
    "primary_eligible", "strict_eligible", "reliability_status",
    "animal_detection_threshold", "region_threshold", "hemisphere_handling"
  )
  out <- dplyr::left_join(
    grid,
    detection_audit[, detection_columns],
    by = c("dataset", "ProteinGroupID"),
    relationship = "many-to-one"
  ) |>
    dplyr::mutate(
      intended_compartment = .data$dataset == .data$intended_dataset,
      off_target_reliability_status = dplyr::case_when(
        .data$intended_compartment ~ "intended_compartment",
        .data$primary_eligible %in% TRUE ~ "off_target_reliably_detected",
        TRUE ~ "off_target_not_reliably_detected"
      ),
      abundance_observation_status = ifelse(
        .data$primary_eligible %in% TRUE,
        "observed_nonimputed_reliable_animal_abundance",
        "detection_only_no_quantitative_contrast"
      ),
      normalization_source =
        "joint_shared_core_sample_median_offsets_applied_to_raw_positive_log2",
      feature_universe = "external_fidelity_registry_canonical_noncontaminant",
      marker_source_policy = "external_SynGO_GO_MGI_Allen_with_documented_fallbacks"
    )
  if (anyDuplicated(out[c("marker_class", "ProteinGroupID", "dataset")])) {
    stop("Marker-by-dataset audit is duplicated at its intended key.", call. = FALSE)
  }
  out
}

ca_protein_direction_v2 <- function(marker_dataset_audit) {
  eligible <- marker_dataset_audit[
    marker_dataset_audit$intended_primary_eligible %in% TRUE,
    ,
    drop = FALSE
  ]
  keys <- unique(eligible[c("marker_class", "ProteinGroupID")])
  parts <- lapply(seq_len(nrow(keys)), function(i) {
    cls <- keys$marker_class[[i]]
    pid <- keys$ProteinGroupID[[i]]
    x <- eligible[
      eligible$marker_class == cls & eligible$ProteinGroupID == pid,
      ,
      drop = FALSE
    ]
    intended <- x[x$intended_compartment %in% TRUE, , drop = FALSE]
    off <- x[!(x$intended_compartment %in% TRUE), , drop = FALSE]
    comparable <- off[
      off$primary_eligible %in% TRUE &
        is.finite(off$median_animal_log2_abundance),
      ,
      drop = FALSE
    ]
    strongest_abundance <- if (nrow(comparable)) {
      max(comparable$median_animal_log2_abundance)
    } else {
      NA_real_
    }
    strongest_detection <- if (nrow(off)) {
      max(off$valid_animal_fraction, na.rm = TRUE)
    } else {
      NA_real_
    }
    if (!is.finite(strongest_detection)) strongest_detection <- NA_real_
    intended_abundance <- intended$median_animal_log2_abundance[[1]]
    margin <- if (is.finite(intended_abundance) && is.finite(strongest_abundance)) {
      intended_abundance - strongest_abundance
    } else {
      NA_real_
    }
    intended_highest <- if (is.finite(margin)) margin > 0 else NA
    classification <- if (!nrow(comparable)) {
      "reliably_detected_only_in_intended_compartment"
    } else if (isTRUE(intended_highest)) {
      "intended_highest"
    } else {
      "discordant_intended_not_highest"
    }
    data.frame(
      marker_class = cls,
      ProteinGroupID = pid,
      marker_gene = intended$marker_gene[[1]],
      display_biological_subpanel = intended$display_biological_subpanel[[1]],
      intended_dataset = intended$intended_dataset[[1]],
      intended_median_log2 = intended_abundance,
      strongest_offtarget_median_log2 = strongest_abundance,
      intended_margin_log2 = margin,
      intended_is_highest = intended_highest,
      n_reliably_comparable_offtarget_datasets = nrow(comparable),
      intended_valid_animal_fraction = intended$valid_animal_fraction[[1]],
      strongest_offtarget_valid_animal_fraction = strongest_detection,
      detection_advantage = intended$valid_animal_fraction[[1]] - strongest_detection,
      all_offtargets_not_reliably_detected = !nrow(comparable),
      expected_direction_classification = classification,
      stringsAsFactors = FALSE
    )
  })
  if (!length(parts)) {
    return(data.frame(
      marker_class = character(), ProteinGroupID = character(),
      marker_gene = character(), display_biological_subpanel = character(),
      intended_dataset = character(), intended_median_log2 = numeric(),
      strongest_offtarget_median_log2 = numeric(), intended_margin_log2 = numeric(),
      intended_is_highest = logical(),
      n_reliably_comparable_offtarget_datasets = integer(),
      intended_valid_animal_fraction = numeric(),
      strongest_offtarget_valid_animal_fraction = numeric(),
      detection_advantage = numeric(),
      all_offtargets_not_reliably_detected = logical(),
      expected_direction_classification = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

ca_round_robin_select <- function(candidates, n_select) {
  if (!nrow(candidates) || n_select < 1L) return(candidates[0, , drop = FALSE])
  panels <- unique(candidates$display_biological_subpanel[
    order(candidates$subpanel_order, candidates$display_biological_subpanel, method = "radix")
  ])
  selected <- integer()
  depth <- 1L
  while (length(selected) < n_select) {
    added <- FALSE
    for (panel in panels) {
      rows <- which(candidates$display_biological_subpanel == panel)
      if (length(rows) >= depth) {
        row <- rows[[depth]]
        if (row %in% selected) next
        selected <- c(selected, row)
        added <- TRUE
        if (length(selected) >= n_select) break
      }
    }
    if (!added) break
    depth <- depth + 1L
  }
  if (length(selected) < n_select) {
    selected <- c(selected, head(setdiff(seq_len(nrow(candidates)), selected), n_select - length(selected)))
  }
  candidates[selected, , drop = FALSE]
}

ca_select_display_markers_v2 <- function(marker_map, n_per_class = 6L) {
  n_per_class <- as.integer(n_per_class)
  classes <- c("Soma markers", "Neuropil markers", "Microglia/PVM markers")
  candidates <- marker_map |>
    dplyr::filter(
      .data$intended_primary_eligible %in% TRUE,
      .data$gene_claim_eligible %in% TRUE,
      is.finite(.data$source_rank)
    ) |>
    dplyr::mutate(
      technical_detection_gate = ifelse(
        .data$intended_strict_eligible %in% TRUE,
        "3_of_3_intended_CON_animals",
        "2_of_3_intended_CON_animals_fallback"
      ),
      detection_priority = ifelse(.data$intended_strict_eligible %in% TRUE, 1L, 2L),
      subpanel_order = ca_marker_subpanel_order(
        .data$marker_class, .data$display_biological_subpanel
      )
    ) |>
    dplyr::arrange(
      match(.data$marker_class, classes), .data$detection_priority,
      .data$source_rank, .data$subpanel_order,
      toupper(.data$marker_gene), .data$ProteinGroupID
    )
  selected_parts <- lapply(classes, function(cls) {
    pool <- candidates[candidates$marker_class == cls, , drop = FALSE]
    strict <- pool[pool$intended_strict_eligible %in% TRUE, , drop = FALSE]
    use <- if (nrow(strict) >= n_per_class) {
      strict
    } else {
      rbind(strict, pool[!(pool$intended_strict_eligible %in% TRUE), , drop = FALSE])
    }
    use <- use[order(
      use$source_rank, use$subpanel_order, toupper(use$marker_gene),
      use$ProteinGroupID, method = "radix"
    ), , drop = FALSE]
    chosen <- ca_round_robin_select(use, n_per_class)
    if (nrow(chosen) < n_per_class) {
      stop("Fewer than ", n_per_class, " display markers are available for ", cls, ".", call. = FALSE)
    }
    chosen$selection_order_within_class <- seq_len(nrow(chosen))
    chosen$used_2_of_3_fallback <- !(chosen$intended_strict_eligible %in% TRUE)
    chosen
  })
  selected <- dplyr::bind_rows(selected_parts) |>
    dplyr::mutate(
      display_selection_policy = paste(
        "external_registry;unique_canonical_mapping;intended_detection;",
        "external_source_rank;round_robin_subpanel_balance;gene_tiebreak"
      ),
      observed_effect_magnitude_used = FALSE,
      rank_label = .data$selection_order_within_class <= 4L,
      rank_label_dataset = .data$intended_dataset
    )
  selected_keys <- paste(selected$marker_class, selected$ProteinGroupID, sep = "\r")
  alternatives <- candidates |>
    dplyr::mutate(
      selected_for_display =
        paste(.data$marker_class, .data$ProteinGroupID, sep = "\r") %in% selected_keys,
      observed_effect_magnitude_used = FALSE,
      display_selection_policy = paste(
        "external_registry;unique_canonical_mapping;intended_detection;",
        "external_source_rank;round_robin_subpanel_balance;gene_tiebreak"
      )
    )
  if (anyDuplicated(selected[c("marker_class", "ProteinGroupID")])) {
    stop("Display selection duplicated a ProteinGroupID within a marker class.", call. = FALSE)
  }
  list(selected = selected, alternatives = alternatives)
}
