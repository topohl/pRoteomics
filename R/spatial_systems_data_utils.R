# Hemisphere-resolved data contract for the spatial systems layer.
#
# WHY THIS FILE EXISTS
#   The biological replicate in this study is AnimalID. Left and right
#   hemispheres are repeated tissue samples WITHIN an animal, never independent
#   replicates. Several existing objects average the two sides before any
#   consumer sees them - including one historically returned as `hemisphere_mat`,
#   which is keyed AnimalID x SpatialUnit and has no side component at all. A
#   bilateral analysis built on such an object is degenerate by construction: it
#   returns an exactly-zero left-right difference and looks like a real result.
#
#   This file is the single place that builds a genuinely side-resolved layer,
#   and the aggregation hierarchy is explicit:
#
#     LEVEL 0  canonical source samples
#     LEVEL 1  ProteinGroupID x AnimalID x Hemisphere x SpatialUnit
#              (technical rows collapsed WITHIN a side, sides kept apart)
#     LEVEL 2  ProteinGroupID x AnimalID x SpatialUnit
#              (equal-weight mean of the available sides)
#     LEVEL 3  explicit left-only / right-only / bilateral matrices
#
# THE ONE-SIDED POLICY IS NOT INVENTED HERE.
#   Both canonical animal-level paths already agree, and this file reuses their
#   vocabulary rather than adding a third rule:
#     R/wgcna_group_effects_utils.R  "equal_weight_mean_available_LR_after_within_hemisphere_mean"
#                                    / "one_sided_observed_no_imputation"
#     R/protigy_input_utils.R        "single_observed_hemisphere_no_imputation"
#   A one-sided cell contributes its observed side unchanged. Nothing is imputed.
#
# NAMING
#   The historical metadata column is `ReplicateGroup` with values Left/Right.
#   That name invites exactly the error this layer exists to prevent, so the
#   biological variable is called `Hemisphere` with canonical values L / R, and
#   `source_hemisphere_field` records where it came from.

# `%||%` comes from the canonical R/null_coalescing.R, loaded via R/paths.R.

sps_contract_version <- function() "spatial_systems_data_v1"

# The metadata column the hemisphere is derived FROM. Recorded as provenance on
# every level so the lineage back to the raw metadata stays visible.
sps_hemisphere_source_field <- function() "ReplicateGroup"

sps_hemisphere_levels <- function() c("L", "R")

# Normalize the historical Left/Right vocabulary exactly once.
sps_normalize_hemisphere <- function(x) {
  raw <- toupper(trimws(as.character(x)))
  out <- ifelse(raw %in% c("L", "LEFT"), "L",
                ifelse(raw %in% c("R", "RIGHT"), "R", NA_character_))
  bad <- !is.na(raw) & nzchar(raw) & is.na(out)
  if (any(bad)) {
    stop("Unrecognised ", sps_hemisphere_source_field(), " value(s): ",
         paste(unique(raw[bad]), collapse = ", "),
         ". Expected Left/Right (or L/R).", call. = FALSE)
  }
  if (any(is.na(out))) {
    stop("Hemisphere is missing for ", sum(is.na(out)),
         " sample(s). A side-resolved layer cannot be built with an unknown side.",
         call. = FALSE)
  }
  factor(out, levels = sps_hemisphere_levels())
}

# The anatomical unit that defines a spatial node, per dataset.
# neuropil resolves Region x Layer; soma and microglia resolve Region.
sps_spatial_unit_column <- function(dataset) {
  if (identical(dataset, "neuron_neuropil")) "RegionLayer" else "Region"
}

sps_spatial_unit_type <- function(dataset) {
  if (identical(dataset, "neuron_neuropil")) "region_layer" else "region"
}

.sps_key <- function(...) paste(..., sep = "\037")

.sps_require_cols <- function(meta, cols, what) {
  missing <- setdiff(cols, names(meta))
  if (length(missing)) {
    stop(what, " is missing required column(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

# Collapse matrix columns by a grouping key using an arithmetic mean.
# Returns the aggregated matrix plus the per-key source counts, so a caller can
# always see how many rows went into a value.
sps_aggregate_columns <- function(mat, keys) {
  if (ncol(mat) != length(keys)) {
    stop("Aggregation matrix/key lengths disagree.", call. = FALSE)
  }
  keys <- as.character(keys)
  idx <- split(seq_along(keys), keys)
  # deterministic column order, independent of the input sample order
  idx <- idx[order(names(idx), method = "radix")]
  out <- vapply(idx, function(i) rowMeans(mat[, i, drop = FALSE], na.rm = TRUE),
                numeric(nrow(mat)))
  if (is.null(dim(out))) out <- matrix(out, nrow = nrow(mat))
  dimnames(out) <- list(rownames(mat), names(idx))
  list(mat = out, n_source_rows = vapply(idx, length, integer(1)))
}

# ---------------------------------------------------------------- LEVEL 1..3

# Build the full hemisphere-resolved hierarchy for one dataset.
#
# `canonical` is the object returned by qc_load_canonical_expression(); passing
# it in lets callers (and tests) avoid re-reading the matrix.
sps_build_spatial_levels <- function(dataset, canonical = NULL,
                                     metadata = NULL, mat = NULL) {
  dataset <- as.character(dataset)[[1]]
  if (is.null(mat) || is.null(metadata)) {
    if (is.null(canonical)) {
      stop("Provide either `canonical`, or both `mat` and `metadata`.", call. = FALSE)
    }
    mat <- canonical$mat
    metadata <- canonical$meta
  }
  unit_col <- sps_spatial_unit_column(dataset)
  .sps_require_cols(metadata,
                    c("Sample", "AnimalID", "StressGroup", "Region",
                      sps_hemisphere_source_field(), unit_col),
                    "Spatial systems metadata")

  # align metadata to the matrix columns
  keep <- colnames(mat) %in% metadata$Sample
  mat <- mat[, keep, drop = FALSE]
  metadata <- metadata[match(colnames(mat), metadata$Sample), , drop = FALSE]
  if (anyNA(metadata$Sample)) {
    stop("Canonical matrix and metadata failed to align on Sample.", call. = FALSE)
  }
  if (!nrow(metadata)) stop("No samples remain for ", dataset, ".", call. = FALSE)

  hemi <- sps_normalize_hemisphere(metadata[[sps_hemisphere_source_field()]])
  animal <- as.character(metadata$AnimalID)
  unit <- as.character(metadata[[unit_col]])
  if (any(is.na(animal) | !nzchar(animal))) {
    stop("AnimalID is missing for some samples; it is the biological replicate.",
         call. = FALSE)
  }
  if (any(is.na(unit) | !nzchar(unit))) {
    stop("SpatialUnit (", unit_col, ") is missing for some samples.", call. = FALSE)
  }

  level0 <- data.frame(
    dataset = dataset,
    Sample = as.character(metadata$Sample),
    AnimalID = animal,
    StressGroup = as.character(metadata$StressGroup),
    Hemisphere = as.character(hemi),
    source_hemisphere_field = sps_hemisphere_source_field(),
    source_hemisphere_value = as.character(metadata[[sps_hemisphere_source_field()]]),
    Region = as.character(metadata$Region),
    Layer = if ("Layer" %in% names(metadata)) as.character(metadata$Layer) else NA_character_,
    SpatialUnit = unit,
    SpatialUnitType = sps_spatial_unit_type(dataset),
    stringsAsFactors = FALSE
  )

  # ---- LEVEL 1: collapse technical rows WITHIN a side; sides stay apart -----
  k1 <- .sps_key(animal, as.character(hemi), unit)
  agg1 <- sps_aggregate_columns(mat, k1)
  parts <- do.call(rbind, strsplit(colnames(agg1$mat), "\037", fixed = TRUE))
  level1_meta <- data.frame(
    dataset = dataset,
    key = colnames(agg1$mat),
    AnimalID = parts[, 1],
    Hemisphere = parts[, 2],
    SpatialUnit = parts[, 3],
    source_hemisphere_field = sps_hemisphere_source_field(),
    n_source_rows = as.integer(agg1$n_source_rows),
    SpatialUnitType = sps_spatial_unit_type(dataset),
    aggregation_method = "arithmetic_mean_within_hemisphere",
    stringsAsFactors = FALSE
  )
  lut <- level0[!duplicated(.sps_key(level0$AnimalID, level0$Hemisphere, level0$SpatialUnit)), ]
  m <- match(level1_meta$key,
             .sps_key(lut$AnimalID, lut$Hemisphere, lut$SpatialUnit))
  level1_meta$StressGroup <- lut$StressGroup[m]
  level1_meta$Region <- lut$Region[m]
  level1_meta$Layer <- lut$Layer[m]
  level1_meta <- level1_meta[, c("dataset", "key", "AnimalID", "StressGroup",
                                 "Hemisphere", "source_hemisphere_field",
                                 "Region", "Layer", "SpatialUnit",
                                 "SpatialUnitType", "n_source_rows",
                                 "aggregation_method")]

  # ---- LEVEL 2: equal-weight mean over the AVAILABLE sides -----------------
  cell <- .sps_key(level1_meta$AnimalID, level1_meta$SpatialUnit)
  agg2 <- sps_aggregate_columns(agg1$mat, cell)
  p2 <- do.call(rbind, strsplit(colnames(agg2$mat), "\037", fixed = TRUE))

  side_of <- split(level1_meta$Hemisphere, cell)
  n_left <- vapply(side_of, function(s) sum(s == "L"), integer(1))
  n_right <- vapply(side_of, function(s) sum(s == "R"), integer(1))
  ord <- match(colnames(agg2$mat), names(side_of))
  n_left <- n_left[ord]; n_right <- n_right[ord]

  left_key <- .sps_key(p2[, 1], "L", p2[, 2])
  right_key <- .sps_key(p2[, 1], "R", p2[, 2])
  has_l <- left_key %in% colnames(agg1$mat)
  has_r <- right_key %in% colnames(agg1$mat)

  # Reuse the canonical hemisphere-status vocabulary rather than adding a third.
  status <- if (exists("protigy_unit_status", mode = "function")) {
    vapply(seq_along(n_left),
           function(i) protigy_unit_status(as.integer(n_left[[i]]),
                                           as.integer(n_right[[i]])),
           character(1))
  } else {
    ifelse(n_left == 1L & n_right == 1L, "bilateral_complete",
           ifelse(n_left == 1L & n_right == 0L, "left_only_observed",
                  ifelse(n_left == 0L & n_right == 1L, "right_only_observed",
                         "invalid_not_output")))
  }
  n_hemispheres <- as.integer(has_l) + as.integer(has_r)
  level2_meta <- data.frame(
    dataset = dataset,
    key = colnames(agg2$mat),
    AnimalID = p2[, 1],
    SpatialUnit = p2[, 2],
    SpatialUnitType = sps_spatial_unit_type(dataset),
    n_hemispheres = n_hemispheres,
    hemisphere_status = status,
    source_hemisphere_field = sps_hemisphere_source_field(),
    aggregation_method = ifelse(
      n_hemispheres == 2L,
      "equal_weight_mean_available_LR_after_within_hemisphere_mean",
      "single_observed_hemisphere_no_imputation"),
    stringsAsFactors = FALSE
  )
  lut2 <- level1_meta[!duplicated(.sps_key(level1_meta$AnimalID, level1_meta$SpatialUnit)), ]
  m2 <- match(level2_meta$key, .sps_key(lut2$AnimalID, lut2$SpatialUnit))
  level2_meta$StressGroup <- lut2$StressGroup[m2]
  level2_meta$Region <- lut2$Region[m2]
  level2_meta$Layer <- lut2$Layer[m2]

  # ---- LEVEL 3: explicit side-specific matrices ----------------------------
  # All three share ProteinGroupID rows and AnimalID::SpatialUnit columns, so a
  # consumer can subtract them directly without re-aligning.
  cell_names <- colnames(agg2$mat)
  take_side <- function(side) {
    want <- .sps_key(p2[, 1], side, p2[, 2])
    out <- matrix(NA_real_, nrow = nrow(agg1$mat), ncol = length(want),
                  dimnames = list(rownames(agg1$mat), cell_names))
    hit <- want %in% colnames(agg1$mat)
    if (any(hit)) out[, hit] <- agg1$mat[, want[hit], drop = FALSE]
    out
  }
  left_mat <- take_side("L")
  right_mat <- take_side("R")

  structure(list(
    dataset = dataset,
    contract_version = sps_contract_version(),
    spatial_unit_column = unit_col,
    spatial_unit_type = sps_spatial_unit_type(dataset),
    source_hemisphere_field = sps_hemisphere_source_field(),
    level0 = level0,
    level1 = list(mat = agg1$mat, meta = level1_meta),
    level2 = list(mat = agg2$mat, meta = level2_meta),
    level3 = list(left = left_mat, right = right_mat, bilateral = agg2$mat)
  ), class = "sps_spatial_levels")
}

# Long form of LEVEL 1, materialised on demand.
#
# The eager long form is ~1e6 rows per dataset (5k proteins x 9 animals x 2
# sides x up to 10 units), which no consumer needs in full, so it is built only
# for the requested proteins.
sps_level1_long <- function(levels, proteins = NULL) {
  if (!inherits(levels, "sps_spatial_levels")) {
    stop("`levels` must come from sps_build_spatial_levels().", call. = FALSE)
  }
  m <- levels$level1$mat
  meta <- levels$level1$meta
  if (!is.null(proteins)) {
    proteins <- intersect(as.character(proteins), rownames(m))
    m <- m[proteins, , drop = FALSE]
  }
  if (!nrow(m)) {
    return(data.frame(dataset = character(), ProteinGroupID = character(),
                      AnimalID = character(), StressGroup = character(),
                      Hemisphere = character(), source_hemisphere_field = character(),
                      Region = character(), Layer = character(),
                      SpatialUnit = character(), n_source_rows = integer(),
                      value = numeric(), stringsAsFactors = FALSE))
  }
  idx <- rep(seq_len(ncol(m)), each = nrow(m))
  data.frame(
    dataset = levels$dataset,
    ProteinGroupID = rep(rownames(m), times = ncol(m)),
    AnimalID = meta$AnimalID[idx],
    StressGroup = meta$StressGroup[idx],
    Hemisphere = meta$Hemisphere[idx],
    source_hemisphere_field = meta$source_hemisphere_field[idx],
    Region = meta$Region[idx],
    Layer = meta$Layer[idx],
    SpatialUnit = meta$SpatialUnit[idx],
    n_source_rows = meta$n_source_rows[idx],
    value = as.vector(m),
    stringsAsFactors = FALSE
  )
}

# Long form of LEVEL 2, carrying both sides beside the bilateral value so a
# reviewer can see what was averaged.
sps_level2_long <- function(levels, proteins = NULL) {
  if (!inherits(levels, "sps_spatial_levels")) {
    stop("`levels` must come from sps_build_spatial_levels().", call. = FALSE)
  }
  keep <- if (is.null(proteins)) rownames(levels$level2$mat) else
    intersect(as.character(proteins), rownames(levels$level2$mat))
  b <- levels$level2$mat[keep, , drop = FALSE]
  l <- levels$level3$left[keep, , drop = FALSE]
  r <- levels$level3$right[keep, , drop = FALSE]
  meta <- levels$level2$meta
  if (!nrow(b)) {
    return(data.frame(dataset = character(), ProteinGroupID = character(),
                      AnimalID = character(), SpatialUnit = character(),
                      left_value = numeric(), right_value = numeric(),
                      n_hemispheres = integer(), hemisphere_status = character(),
                      bilateral_value = numeric(), stringsAsFactors = FALSE))
  }
  idx <- rep(seq_len(ncol(b)), each = nrow(b))
  data.frame(
    dataset = levels$dataset,
    ProteinGroupID = rep(rownames(b), times = ncol(b)),
    AnimalID = meta$AnimalID[idx],
    StressGroup = meta$StressGroup[idx],
    Region = meta$Region[idx],
    Layer = meta$Layer[idx],
    SpatialUnit = meta$SpatialUnit[idx],
    left_value = as.vector(l),
    right_value = as.vector(r),
    n_hemispheres = meta$n_hemispheres[idx],
    hemisphere_status = meta$hemisphere_status[idx],
    aggregation_method = meta$aggregation_method[idx],
    bilateral_value = as.vector(b),
    stringsAsFactors = FALSE
  )
}

# Convenience: the side-specific sample metadata for a one-sided model fit.
# Returns ONLY the requested side, so a left model can never see a right sample.
sps_side_samples <- function(levels, side) {
  side <- match.arg(toupper(side), sps_hemisphere_levels())
  levels$level0[levels$level0$Hemisphere == side, , drop = FALSE]
}

# Guard used by consumers: assert a frame never treats a side as a replicate.
sps_assert_animal_is_replicate <- function(df, label = "table") {
  if (!"AnimalID" %in% names(df)) {
    stop(label, " has no AnimalID column; AnimalID is the biological replicate.",
         call. = FALSE)
  }
  if (sps_hemisphere_source_field() %in% names(df)) {
    stop(label, " carries a raw ", sps_hemisphere_source_field(),
         " column. Use the normalized Hemisphere variable so a side is never ",
         "mistaken for a replicate.", call. = FALSE)
  }
  invisible(TRUE)
}
