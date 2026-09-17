# Shared helpers for historical dataset-specific ProTigy inputs and the
# animal-level ProTigy handoff.

protigy_prepare_legacy_expression <- function(df) {
  description_candidates <- intersect(
    c("Description", "First.Protein.Description"),
    colnames(df)
  )
  if (length(description_candidates) > 1L) {
    description_values <- lapply(
      df[description_candidates],
      function(x) as.character(x)
    )
    if (!all(vapply(
      description_values[-1L],
      identical,
      logical(1),
      description_values[[1L]]
    ))) {
      stop("Conflicting protein description columns were supplied.", call. = FALSE)
    }
  }
  description_source <- if (length(description_candidates)) {
    description_candidates[[1L]]
  } else {
    "id_compatibility_fallback"
  }
  source_descriptions <- if (length(description_candidates)) {
    as.character(df[[description_candidates[[1L]]]])
  } else {
    NULL
  }

  name_col <- intersect(c("T: Protein.Names", "Protein.Names"), colnames(df))
  if (length(name_col) == 1L) {
    colnames(df)[colnames(df) == name_col] <- "id"
    df <- df[, c("id", setdiff(colnames(df), "id")), drop = FALSE]
  }

  group_col <- intersect(c("T: Protein.Group", "Protein.Group"), colnames(df))
  if (length(group_col) == 1L) {
    idx <- which(colnames(df) == group_col)
    if (idx > 1L) {
      df <- df[, seq_len(idx - 1L), drop = FALSE]
    } else {
      df <- df[, FALSE, drop = FALSE]
    }
  }
  if ("id" %in% colnames(df)) {
    descriptions <- if (is.null(source_descriptions)) {
      as.character(df$id)
    } else {
      source_descriptions
    }
    fallback <- is.na(descriptions) | !nzchar(trimws(descriptions))
    descriptions[fallback] <- as.character(df$id[fallback])
    df$Description <- descriptions
    df <- df[, c("id", "Description", setdiff(colnames(df), c("id", "Description"))), drop = FALSE]
    attr(df, "protigy_description_source") <- description_source
    attr(df, "protigy_description_fallback_rows") <- fallback
  }
  df
}

gct_v1.3_fields <- function(line) {
  tabs <- gregexpr("\t", line, fixed = TRUE)[[1]]
  if (identical(tabs, -1L)) return(line)
  starts <- c(1L, tabs + 1L)
  ends <- c(tabs - 1L, nchar(line))
  substring(line, starts, ends)
}

validate_gct_v1.3 <- function(file,
                              expected_rows = NULL,
                              expected_sample_cols = NULL,
                              expected_row_descriptor_cols = NULL,
                              expected_col_metadata_rows = NULL) {
  lines <- readLines(file, warn = FALSE)
  if (length(lines) < 3L || !identical(lines[[1]], "#1.3")) {
    stop("Invalid GCT marker in written file: ", file, call. = FALSE)
  }
  dimension_fields <- gct_v1.3_fields(lines[[2]])
  if (
    length(dimension_fields) != 4L ||
      any(!grepl("^[0-9]+$", dimension_fields))
  ) {
    stop("Invalid GCT v1.3 dimension line in written file: ", file, call. = FALSE)
  }
  dims <- suppressWarnings(as.integer(dimension_fields))
  if (anyNA(dims)) {
    stop("Invalid GCT v1.3 dimension line in written file: ", file, call. = FALSE)
  }
  names(dims) <- c(
    "matrix_rows",
    "sample_columns",
    "row_descriptor_columns",
    "column_metadata_rows"
  )
  expected <- c(
    matrix_rows = expected_rows,
    sample_columns = expected_sample_cols,
    row_descriptor_columns = expected_row_descriptor_cols,
    column_metadata_rows = expected_col_metadata_rows
  )
  expected <- expected[!vapply(expected, is.null, logical(1))]
  if (length(expected) && any(dims[names(expected)] != as.integer(expected))) {
    stop(
      "GCT v1.3 dimension mismatch in written file: ", file,
      ". Expected ", paste(unname(expected), collapse = "\t"),
      " but found ", paste(unname(dims), collapse = "\t"),
      call. = FALSE
    )
  }

  expected_line_count <- 3L + dims[["column_metadata_rows"]] +
    dims[["matrix_rows"]]
  if (length(lines) != expected_line_count) {
    stop(
      "GCT v1.3 line-count mismatch in written file: ", file,
      ". Expected ", expected_line_count, " lines but found ", length(lines), ".",
      call. = FALSE
    )
  }

  expected_field_count <- 1L + dims[["row_descriptor_columns"]] +
    dims[["sample_columns"]]
  structural_line_numbers <- c(
    3L,
    if (dims[["column_metadata_rows"]] > 0L) {
      seq.int(4L, 3L + dims[["column_metadata_rows"]])
    },
    if (dims[["matrix_rows"]] > 0L) {
      seq.int(
        4L + dims[["column_metadata_rows"]],
        expected_line_count
      )
    }
  )
  structural_fields <- lapply(lines[structural_line_numbers], gct_v1.3_fields)
  observed_field_counts <- lengths(structural_fields)
  bad_width <- observed_field_counts != expected_field_count
  if (any(bad_width)) {
    stop(
      "GCT v1.3 field-count mismatch in written file: ", file,
      ". Expected ", expected_field_count, " fields on every structural row; line(s) ",
      paste(structural_line_numbers[bad_width], collapse = ", "), " contained ",
      paste(observed_field_counts[bad_width], collapse = ", "), ".",
      call. = FALSE
    )
  }

  header_fields <- structural_fields[[1]]
  if (!identical(header_fields[[1]], "id")) {
    stop("GCT v1.3 header must begin with reserved field 'id': ", file, call. = FALSE)
  }
  if (
    dims[["row_descriptor_columns"]] == 1L &&
      !identical(header_fields[[2]], "Description")
  ) {
    stop(
      "ProTigy-compatible GCT v1.3 with one row descriptor must name it 'Description': ",
      file,
      call. = FALSE
    )
  }
  sample_start <- 2L + dims[["row_descriptor_columns"]]
  sample_end <- expected_field_count
  sample_ids <- header_fields[seq.int(sample_start, sample_end)]
  if (any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("GCT v1.3 sample IDs must be non-empty and unique: ", file, call. = FALSE)
  }

  data_offset <- 1L + dims[["column_metadata_rows"]]
  data_fields <- structural_fields[
    seq.int(data_offset + 1L, data_offset + dims[["matrix_rows"]])
  ]
  row_ids <- vapply(data_fields, `[[`, character(1), 1L)
  if (any(!nzchar(row_ids)) || anyDuplicated(row_ids)) {
    stop("GCT v1.3 row IDs must be non-empty and unique: ", file, call. = FALSE)
  }
  if (dims[["row_descriptor_columns"]] > 0L) {
    descriptor_fields <- unlist(lapply(
      data_fields,
      function(fields) fields[seq.int(2L, 1L + dims[["row_descriptor_columns"]])]
    ), use.names = FALSE)
    if (any(is.na(descriptor_fields)) || any(!nzchar(descriptor_fields))) {
      stop(
        "GCT v1.3 row-descriptor values must align one-to-one with non-empty protein rows: ",
        file,
        call. = FALSE
      )
    }
  }
  matrix_fields <- unlist(lapply(
    data_fields,
    function(fields) fields[seq.int(sample_start, sample_end)]
  ), use.names = FALSE)
  numeric_matrix <- suppressWarnings(as.numeric(matrix_fields))
  if (any(is.na(numeric_matrix)) || any(!is.finite(numeric_matrix))) {
    stop("GCT v1.3 matrix cells must all be finite numeric values: ", file, call. = FALSE)
  }
  invisible(TRUE)
}

write_gct_v1.3 <- function(df, file, metadata) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"id" %in% colnames(df)) {
    stop("ProTigy-targeted GCT input must contain reserved column 'id'.", call. = FALSE)
  }
  if (!"Description" %in% colnames(df)) {
    df$Description <- as.character(df$id)
  }
  missing_description <- is.na(df$Description) | !nzchar(trimws(as.character(df$Description)))
  df$Description[missing_description] <- as.character(df$id[missing_description])
  if (any(grepl("[\t\r\n]", as.character(df$Description)))) {
    stop("Description values cannot contain tab or newline characters.", call. = FALSE)
  }
  meta_row_idx <- which(df$id %in% c(
    names(metadata),
    "region_layer_ExpGroup",
    "phenotypeWithinUnit"
  ))
  data_rows <- df[-meta_row_idx, , drop = FALSE]
  data_rows <- data_rows[!is.na(data_rows$id) & data_rows$id != "", , drop = FALSE]
  sample_cols <- intersect(setdiff(colnames(df), c("id", "Description")), metadata$sample_id)
  sample_cols <- intersect(sample_cols, colnames(data_rows))
  unexpected_descriptor_cols <- setdiff(
    colnames(data_rows),
    c("id", "Description", sample_cols)
  )
  if (length(unexpected_descriptor_cols)) {
    stop(
      "Unexpected additional row descriptor column(s) in ProTigy-targeted GCT input: ",
      paste(unexpected_descriptor_cols, collapse = ", "),
      call. = FALSE
    )
  }
  row_descriptor_cols <- "Description"
  if (length(sample_cols) == 0L || nrow(data_rows) == 0L) {
    warning(sprintf("No data to write for file: %s", file))
    return()
  }
  is_numeric_row <- function(row) {
    values <- suppressWarnings(as.numeric(row[sample_cols]))
    all(!is.na(values) & is.finite(values))
  }
  numeric_rows_idx <- apply(data_rows, 1, is_numeric_row)
  data_rows <- data_rows[numeric_rows_idx, , drop = FALSE]
  for (col in sample_cols) {
    data_rows[[col]] <- as.numeric(data_rows[[col]])
  }
  if (nrow(data_rows) == 0L) {
    warning(sprintf("No numeric data to write for file: %s", file))
    return()
  }
  physical_cols <- c("id", row_descriptor_cols, sample_cols)
  meta_rows <- df[meta_row_idx, physical_cols, drop = FALSE]
  meta_rows <- meta_rows[!is.na(meta_rows$id) & meta_rows$id != "", , drop = FALSE]
  if (length(row_descriptor_cols) && nrow(meta_rows)) {
    meta_rows[, row_descriptor_cols] <- "na"
  }
  con <- file(file, "wt")
  on.exit(close(con))
  writeLines("#1.3", con)
  writeLines(
    sprintf(
      "%d\t%d\t%d\t%d",
      nrow(data_rows),
      length(sample_cols),
      length(row_descriptor_cols),
      nrow(meta_rows)
    ),
    con
  )
  writeLines(paste(physical_cols, collapse = "\t"), con)
  if (nrow(meta_rows) > 0L) {
    for (i in seq_len(nrow(meta_rows))) {
      writeLines(
        paste(unlist(meta_rows[i, physical_cols]), collapse = "\t"),
        con
      )
    }
  }
  for (i in seq_len(nrow(data_rows))) {
    writeLines(
      paste(unlist(data_rows[i, physical_cols]), collapse = "\t"),
      con
    )
  }
  close(con)
  on.exit()
  validate_gct_v1.3(
    file,
    nrow(data_rows),
    length(sample_cols),
    length(row_descriptor_cols),
    nrow(meta_rows)
  )
}

self_test_write_gct_v1.3 <- function() {
  tmp <- tempfile(fileext = ".gct")
  test_metadata <- data.frame(
    sample_id = c("sample1", "sample2", "sample3"),
    group = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )
  test_df <- data.frame(
    id = c("group", "protein_a", "protein_b"),
    sample1 = c("A", 1, 4),
    sample2 = c("B", 2, 5),
    sample3 = c("A", 3, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  write_gct_v1.3(test_df, tmp, test_metadata)
  lines <- readLines(tmp, warn = FALSE)
  stopifnot(identical(lines[[1]], "#1.3"))
  stopifnot(identical(lines[[2]], "2\t3\t1\t1"))
  invisible(TRUE)
}

protigy_excluded <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

protigy_resolve_hemisphere <- function(replicate_group, sample_id) {
  if (length(replicate_group) != length(sample_id)) {
    stop("ReplicateGroup and sample_id must have the same length.", call. = FALSE)
  }
  raw <- trimws(as.character(replicate_group))
  hemisphere <- rep(NA_character_, length(raw))
  hemisphere[tolower(raw) %in% c("left", "l")] <- "Left"
  hemisphere[tolower(raw) %in% c("right", "r")] <- "Right"
  unresolved <- is.na(hemisphere)
  if (any(unresolved)) {
    stop(
      "Unresolved or unsupported ReplicateGroup value(s): ",
      paste(unique(raw[unresolved]), collapse = ", "),
      ". Only Left/L and Right/R are accepted.",
      call. = FALSE
    )
  }

  sample_id <- as.character(sample_id)
  has_left <- grepl("_L_", sample_id, fixed = TRUE)
  has_right <- grepl("_R_", sample_id, fixed = TRUE)
  unresolved_name <- has_left == has_right
  if (any(unresolved_name)) {
    stop(
      "Hemisphere cannot be resolved uniquely from sample_id(s): ",
      paste(sample_id[unresolved_name], collapse = "; "),
      call. = FALSE
    )
  }
  from_name <- ifelse(has_left, "Left", "Right")
  conflict <- hemisphere != from_name
  if (any(conflict)) {
    stop(
      "ReplicateGroup conflicts with sample_id hemisphere for: ",
      paste(sample_id[conflict], collapse = "; "),
      call. = FALSE
    )
  }
  hemisphere
}

protigy_required_metadata_columns <- function() {
  c(
    "sample_id", "AnimalID", "ExpGroup", "region", "layer", "celltype",
    "celltype_layer", "ReplicateGroup"
  )
}

protigy_single_value <- function(x, label, key) {
  values <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(values) != 1L) {
    stop(
      "Expected exactly one ", label, " for ", key,
      "; found: ", if (length(values)) paste(values, collapse = ", ") else "none",
      call. = FALSE
    )
  }
  values[[1]]
}

protigy_unit_status <- function(n_left, n_right) {
  if (n_left > 1L || n_right > 1L) {
    issues <- character()
    if (n_left > 1L) issues <- c(issues, "duplicate_left")
    if (n_right > 1L) issues <- c(issues, "duplicate_right")
    if (n_left == 0L) issues <- c(issues, "missing_left")
    if (n_right == 0L) issues <- c(issues, "missing_right")
    return(paste(issues, collapse = "_"))
  }
  if (n_left == 1L && n_right == 1L) return("bilateral_complete")
  if (n_left == 1L && n_right == 0L) return("left_only_observed")
  if (n_left == 0L && n_right == 1L) return("right_only_observed")
  if (n_left == 0L && n_right == 0L) return("missing_both")
  issues <- character()
  if (n_left < 0L || n_right < 0L) issues <- c(issues, "negative_count")
  paste(c("invalid_hemisphere_cardinality", issues), collapse = "_")
}

protigy_primary_hemisphere_statuses <- function() {
  c("bilateral_complete", "left_only_observed", "right_only_observed")
}

protigy_aggregate_expression_columns <- function(expression_matrix, aggregation_audit) {
  expression_matrix <- as.matrix(expression_matrix)
  if (is.null(rownames(expression_matrix)) || anyNA(rownames(expression_matrix)) ||
      any(!nzchar(rownames(expression_matrix)))) {
    stop("Expression matrix requires complete, non-empty row names before animal aggregation.", call. = FALSE)
  }
  if (anyDuplicated(rownames(expression_matrix))) {
    stop("Expression matrix row names must be unique before animal aggregation.", call. = FALSE)
  }
  if (is.null(colnames(expression_matrix)) || anyDuplicated(colnames(expression_matrix))) {
    stop("Expression matrix column names must be present and unique before animal aggregation.", call. = FALSE)
  }

  required <- c(
    "left_sample", "right_sample", "n_left_source_samples", "n_right_source_samples",
    "output_column_name", "hemisphere_status"
  )
  missing <- setdiff(required, names(aggregation_audit))
  if (length(missing)) {
    stop("Aggregation audit is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  audit <- aggregation_audit
  if ("inclusion_status" %in% names(audit)) {
    audit <- audit[audit$inclusion_status == "included_primary", , drop = FALSE]
  }
  if (!nrow(audit)) stop("Aggregation audit contains no primary animal-level units.", call. = FALSE)
  if (anyDuplicated(audit$output_column_name)) {
    stop("Aggregation audit output column names must be unique.", call. = FALSE)
  }

  audit_samples <- function(x) {
    if (!length(x) || is.na(x) || !nzchar(trimws(as.character(x)))) return(character())
    samples <- trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1]])
    samples[nzchar(samples)]
  }

  aggregated <- matrix(
    NA_real_,
    nrow = nrow(expression_matrix),
    ncol = nrow(audit),
    dimnames = list(rownames(expression_matrix), as.character(audit$output_column_name))
  )
  for (i in seq_len(nrow(audit))) {
    left <- audit_samples(audit$left_sample[[i]])
    right <- audit_samples(audit$right_sample[[i]])
    if (length(left) != as.integer(audit$n_left_source_samples[[i]]) ||
        length(right) != as.integer(audit$n_right_source_samples[[i]])) {
      stop(
        "Aggregation audit source-sample counts disagree for ", audit$output_column_name[[i]], ".",
        call. = FALSE
      )
    }
    source_samples <- c(left, right)
    missing_samples <- setdiff(source_samples, colnames(expression_matrix))
    if (length(missing_samples)) {
      stop(
        "Aggregation audit references sample column(s) absent from the gene matrix: ",
        paste(missing_samples, collapse = "; "),
        call. = FALSE
      )
    }
    status <- as.character(audit$hemisphere_status[[i]])
    if (identical(status, "bilateral_complete") && length(left) == 1L && length(right) == 1L) {
      aggregated[, i] <- rowMeans(expression_matrix[, source_samples, drop = FALSE])
    } else if (status %in% c("left_only_observed", "right_only_observed") && length(source_samples) == 1L) {
      aggregated[, i] <- expression_matrix[, source_samples[[1]]]
    } else {
      stop(
        "Unsupported hemisphere status/cardinality for ", audit$output_column_name[[i]],
        ": ", status, " (Left=", length(left), ", Right=", length(right), ").",
        call. = FALSE
      )
    }
  }

  if (!identical(rownames(aggregated), rownames(expression_matrix))) {
    stop("Animal aggregation changed expression-matrix row identity.", call. = FALSE)
  }
  aggregated
}

protigy_build_annotated_table <- function(expression, sample_metadata) {
  sample_cols <- as.character(sample_metadata$sample_id)
  if (!"Description" %in% names(expression)) {
    stop("Animal-level expression must contain the ProTigy Description row descriptor.", call. = FALSE)
  }
  if (!identical(setdiff(names(expression), c("id", "Description")), sample_cols)) {
    stop("Animal-level expression and metadata columns are not aligned.", call. = FALSE)
  }
  meta_rows <- lapply(names(sample_metadata), function(meta_col) {
    row <- rep(NA_character_, ncol(expression))
    names(row) <- names(expression)
    row[["id"]] <- meta_col
    row[sample_cols] <- as.character(sample_metadata[[meta_col]])
    row
  })
  meta_df <- as.data.frame(do.call(rbind, meta_rows), stringsAsFactors = FALSE)
  names(meta_df) <- names(expression)
  rbind(meta_df, expression)
}

protigy_prepare_animal_level <- function(expression,
                                         metadata,
                                         dataset,
                                         validate_e9_design = TRUE,
                                         fail_on_invalid = TRUE) {
  dataset <- validate_dataset(dataset)
  spatial_unit <- dataset_spatial_unit(dataset)
  expression <- protigy_prepare_legacy_expression(as.data.frame(
    expression,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  if (!"id" %in% names(expression)) {
    stop(
      "Expression must contain exactly one historical ProTigy ID column: ",
      "T: Protein.Names or Protein.Names.",
      call. = FALSE
    )
  }
  description_source <- attr(expression, "protigy_description_source", exact = TRUE)
  description_fallback_rows <- attr(
    expression,
    "protigy_description_fallback_rows",
    exact = TRUE
  )
  sample_cols <- setdiff(names(expression), c("id", "Description"))
  if (!length(sample_cols)) stop("Expression contains no sample columns.", call. = FALSE)
  if (anyDuplicated(sample_cols)) stop("Expression sample columns are duplicated.", call. = FALSE)

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  missing_metadata_cols <- setdiff(protigy_required_metadata_columns(), names(metadata))
  if (length(missing_metadata_cols)) {
    stop(
      "Metadata is missing required column(s): ",
      paste(missing_metadata_cols, collapse = ", "),
      call. = FALSE
    )
  }
  missing_meta <- setdiff(sample_cols, as.character(metadata$sample_id))
  if (length(missing_meta)) {
    stop(
      "Expression sample column(s) have no metadata row: ",
      paste(missing_meta, collapse = "; "),
      call. = FALSE
    )
  }
  used_meta_all <- metadata[metadata$sample_id %in% sample_cols, , drop = FALSE]
  duplicate_used <- unique(used_meta_all$sample_id[duplicated(used_meta_all$sample_id)])
  if (length(duplicate_used)) {
    stop(
      "Metadata sample_id is duplicated for used sample(s): ",
      paste(duplicate_used, collapse = "; "),
      call. = FALSE
    )
  }
  used_meta_all <- used_meta_all[match(sample_cols, used_meta_all$sample_id), , drop = FALSE]
  excluded <- if ("exclude" %in% names(used_meta_all)) {
    protigy_excluded(used_meta_all$exclude)
  } else if ("Exclude" %in% names(used_meta_all)) {
    protigy_excluded(used_meta_all$Exclude)
  } else {
    rep(FALSE, nrow(used_meta_all))
  }
  included_meta <- used_meta_all[!excluded, , drop = FALSE]
  included_samples <- as.character(included_meta$sample_id)
  if (!nrow(included_meta)) stop("No samples remain after applying exclusions.", call. = FALSE)

  dataset_match <- metadata_matches_dataset(included_meta, dataset)
  if (any(!dataset_match)) {
    wrong <- included_meta$sample_id[!dataset_match]
    stop(
      "Matched samples belong to the wrong dataset family for '", dataset, "': ",
      paste(wrong, collapse = "; "),
      call. = FALSE
    )
  }
  included_meta$hemisphere <- protigy_resolve_hemisphere(
    included_meta$ReplicateGroup,
    included_meta$sample_id
  )
  included_meta$AnimalID <- trimws(as.character(included_meta$AnimalID))
  included_meta$ExpGroup <- trimws(as.character(included_meta$ExpGroup))
  included_meta$region <- trimws(as.character(included_meta$region))
  included_meta$layer <- trimws(as.character(included_meta$layer))
  included_meta$celltype <- trimws(as.character(included_meta$celltype))
  included_meta$celltype_layer <- trimws(as.character(included_meta$celltype_layer))
  required_values <- c("AnimalID", "ExpGroup", "region")
  if (!identical(dataset, "microglia")) required_values <- c(required_values, "layer")
  for (column in required_values) {
    bad <- is.na(included_meta[[column]]) | !nzchar(included_meta[[column]])
    if (any(bad)) {
      stop(
        "Missing ", column, " for included sample(s): ",
        paste(included_meta$sample_id[bad], collapse = "; "),
        call. = FALSE
      )
    }
  }
  if (identical(dataset, "microglia")) {
    included_meta$layer <- "microglia"
    included_meta$celltype_layer <- "microglia"
  }

  if (identical(spatial_unit, "region_layer")) {
    included_meta$canonical_spatial_unit <- paste(
      included_meta$region,
      included_meta$layer,
      sep = "_"
    )
  } else if (identical(spatial_unit, "region")) {
    included_meta$canonical_spatial_unit <- included_meta$region
  } else {
    stop("Unsupported dataset spatial unit: ", spatial_unit, call. = FALSE)
  }
  included_meta$unit_key <- paste(
    included_meta$AnimalID,
    included_meta$canonical_spatial_unit,
    sep = "\r"
  )

  validation_issues <- character()
  animal_groups <- split(included_meta$ExpGroup, included_meta$AnimalID)
  group_counts_by_animal <- vapply(animal_groups, function(x) length(unique(x)), integer(1))
  if (any(group_counts_by_animal != 1L)) {
    validation_issues <- c(
      validation_issues,
      paste0(
        "AnimalID maps to more than one ExpGroup: ",
        paste(names(group_counts_by_animal)[group_counts_by_animal != 1L], collapse = ", ")
      )
    )
  }

  unit_keys <- unique(included_meta$unit_key)
  unit_parts <- lapply(unit_keys, function(key) included_meta[included_meta$unit_key == key, , drop = FALSE])
  unit_audit <- lapply(unit_parts, function(part) {
    animal <- unique(part$AnimalID)
    region <- unique(part$region)
    canonical <- unique(part$canonical_spatial_unit)
    layer_by_hemisphere <- split(part$layer, part$hemisphere)
    bad_hemisphere_layers <- names(layer_by_hemisphere)[vapply(
      layer_by_hemisphere,
      function(x) length(unique(x)) != 1L,
      logical(1)
    )]
    if (length(bad_hemisphere_layers)) {
      stop(
        "Soma/descriptive layer is not unique within AnimalID x region x hemisphere for ",
        animal[[1]], " / ", region[[1]], ": ",
        paste(bad_hemisphere_layers, collapse = ", "),
        call. = FALSE
      )
    }
    layer <- protigy_single_value(part$layer, "descriptive layer", paste(animal[[1]], canonical[[1]], sep = " / "))
    exp_group <- protigy_single_value(part$ExpGroup, "ExpGroup", paste(animal[[1]], canonical[[1]], sep = " / "))
    celltype <- protigy_single_value(part$celltype, "celltype", paste(animal[[1]], canonical[[1]], sep = " / "))
    celltype_layer <- protigy_single_value(part$celltype_layer, "celltype_layer", paste(animal[[1]], canonical[[1]], sep = " / "))
    left <- as.character(part$sample_id[part$hemisphere == "Left"])
    right <- as.character(part$sample_id[part$hemisphere == "Right"])
    output_name <- if (identical(dataset, "microglia")) {
      paste(animal[[1]], region[[1]], "microglia", sep = "_")
    } else {
      paste(animal[[1]], region[[1]], layer, sep = "_")
    }
    data.frame(
      dataset = dataset,
      AnimalID = animal[[1]],
      ExpGroup = exp_group,
      region = region[[1]],
      layer = layer,
      canonical_spatial_unit = canonical[[1]],
      left_sample = paste(left, collapse = ";"),
      right_sample = paste(right, collapse = ";"),
      n_left_source_samples = length(left),
      n_right_source_samples = length(right),
      output_column_name = output_name,
      hemisphere_status = protigy_unit_status(length(left), length(right)),
      aggregation_method = if (length(left) == 1L && length(right) == 1L) {
        "equal_weight_mean_LR_on_existing_imputed_log2_values"
      } else if (length(left) + length(right) == 1L) {
        "single_observed_hemisphere_no_imputation"
      } else {
        NA_character_
      },
      inclusion_status = if (
        protigy_unit_status(length(left), length(right)) %in%
          protigy_primary_hemisphere_statuses()
      ) "included_primary" else "invalid_not_output",
      celltype = celltype,
      celltype_layer = celltype_layer,
      stringsAsFactors = FALSE
    )
  })
  aggregation_audit <- do.call(rbind, unit_audit)
  aggregation_audit <- aggregation_audit[order(
    aggregation_audit$AnimalID,
    aggregation_audit$region,
    aggregation_audit$layer
  ), , drop = FALSE]
  rownames(aggregation_audit) <- NULL

  if (anyDuplicated(aggregation_audit$output_column_name)) {
    validation_issues <- c(validation_issues, "Output column names are not unique.")
  }
  invalid_cardinality <- !aggregation_audit$hemisphere_status %in%
    protigy_primary_hemisphere_statuses()
  if (any(invalid_cardinality)) {
    details <- paste0(
      aggregation_audit$AnimalID[invalid_cardinality], "/",
      aggregation_audit$canonical_spatial_unit[invalid_cardinality],
      " (Left=", aggregation_audit$n_left_source_samples[invalid_cardinality],
      ", Right=", aggregation_audit$n_right_source_samples[invalid_cardinality], ")"
    )
    validation_issues <- c(
      validation_issues,
      paste0("Invalid hemisphere cardinality unit(s): ", paste(details, collapse = "; "))
    )
  }

  if (isTRUE(validate_e9_design)) {
    animals <- sort(unique(included_meta$AnimalID))
    if (length(animals) != 9L) {
      validation_issues <- c(
        validation_issues,
        paste0("Expected exactly 9 unique AnimalIDs; found ", length(animals), ".")
      )
    }
    observed_groups <- sort(unique(included_meta$ExpGroup))
    if (!identical(observed_groups, c("1", "2", "3"))) {
      validation_issues <- c(
        validation_issues,
        paste0("Expected ExpGroup codes 1/2/3; found: ", paste(observed_groups, collapse = ", "))
      )
    }
    animal_group <- unique(included_meta[, c("AnimalID", "ExpGroup"), drop = FALSE])
    group_sizes <- table(animal_group$ExpGroup)
    for (group in c("1", "2", "3")) {
      n_group <- if (group %in% names(group_sizes)) as.integer(group_sizes[[group]]) else 0L
      if (n_group != 3L) {
        validation_issues <- c(
          validation_issues,
          paste0("Expected 3 animals in ExpGroup ", group, "; found ", n_group, ".")
        )
      }
    }
  }

  assignment_included <- included_meta[, c(
    "sample_id", "AnimalID", "ExpGroup", "region", "layer", "hemisphere",
    "canonical_spatial_unit", "unit_key"
  ), drop = FALSE]
  assignment_match <- match(assignment_included$unit_key, paste(
      aggregation_audit$AnimalID,
      aggregation_audit$canonical_spatial_unit,
      sep = "\r"
    ))
  assignment_unit_status <- aggregation_audit$hemisphere_status[assignment_match]
  assignment_primary <- assignment_unit_status %in% protigy_primary_hemisphere_statuses()
  assignment_included$output_column_name <- ifelse(
    assignment_primary,
    aggregation_audit$output_column_name[assignment_match],
    NA_character_
  )
  assignment_included$inclusion_status <- ifelse(
    assignment_primary,
    "included",
    "invalid_unit_not_output"
  )
  assignment_included$assigned_output_count <- as.integer(assignment_primary)
  assignment_included$dataset <- dataset
  assignment_included <- assignment_included[, c(
    "dataset", "sample_id", "AnimalID", "ExpGroup", "region", "layer",
    "hemisphere", "canonical_spatial_unit", "output_column_name",
    "inclusion_status", "assigned_output_count"
  ), drop = FALSE]
  excluded_meta <- used_meta_all[excluded, , drop = FALSE]
  if (nrow(excluded_meta)) {
    assignment_excluded <- data.frame(
      dataset = dataset,
      sample_id = as.character(excluded_meta$sample_id),
      AnimalID = as.character(excluded_meta$AnimalID),
      ExpGroup = as.character(excluded_meta$ExpGroup),
      region = as.character(excluded_meta$region),
      layer = as.character(excluded_meta$layer),
      hemisphere = NA_character_,
      canonical_spatial_unit = NA_character_,
      output_column_name = NA_character_,
      inclusion_status = "excluded_by_metadata",
      assigned_output_count = 0L,
      stringsAsFactors = FALSE
    )
    source_assignment <- rbind(assignment_included, assignment_excluded)
  } else {
    source_assignment <- assignment_included
  }
  source_assignment <- source_assignment[order(
    source_assignment$inclusion_status,
    source_assignment$sample_id
  ), , drop = FALSE]
  rownames(source_assignment) <- NULL
  if (any(source_assignment$assigned_output_count[source_assignment$inclusion_status == "included"] != 1L)) {
    validation_issues <- c(
      validation_issues,
      "At least one included input sample does not map to exactly one output unit."
    )
  }
  if (anyDuplicated(source_assignment$sample_id[source_assignment$inclusion_status == "included"])) {
    validation_issues <- c(validation_issues, "An included input sample is reused across output units.")
  }

  primary_audit <- aggregation_audit[
    aggregation_audit$hemisphere_status %in% protigy_primary_hemisphere_statuses(),
    ,
    drop = FALSE
  ]
  feature_numeric <- lapply(expression[included_samples], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  feature_numeric <- as.data.frame(feature_numeric, check.names = FALSE)
  numeric_rows <- if (nrow(expression)) {
    apply(feature_numeric, 1, function(x) all(!is.na(x)))
  } else {
    logical()
  }
  feature_keep <- !is.na(expression$id) & nzchar(as.character(expression$id)) & numeric_rows
  expression_numeric <- feature_numeric[feature_keep, , drop = FALSE]
  feature_ids <- as.character(expression$id[feature_keep])
  feature_descriptions <- as.character(expression$Description[feature_keep])
  description_fallback_count <- if (length(description_fallback_rows)) {
    sum(description_fallback_rows[feature_keep])
  } else {
    length(feature_ids)
  }

  animal_expression <- data.frame(
    id = feature_ids,
    Description = feature_descriptions,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (nrow(primary_audit)) {
    for (i in seq_len(nrow(primary_audit))) {
      source_samples <- c(
        primary_audit$left_sample[[i]],
        primary_audit$right_sample[[i]]
      )
      source_samples <- source_samples[nzchar(source_samples)]
      animal_expression[[primary_audit$output_column_name[[i]]]] <- if (
        length(source_samples) == 2L
      ) {
        rowMeans(expression_numeric[, source_samples, drop = FALSE])
      } else {
        expression_numeric[[source_samples[[1]]]]
      }
    }
  }
  output_metadata <- data.frame(
    sample_id = primary_audit$output_column_name,
    AnimalID = primary_audit$AnimalID,
    ExpGroup = primary_audit$ExpGroup,
    region = primary_audit$region,
    layer = primary_audit$layer,
    celltype = primary_audit$celltype,
    celltype_layer = primary_audit$celltype_layer,
    phenotypeWithinUnit = paste(
      primary_audit$region,
      primary_audit$layer,
      primary_audit$ExpGroup,
      sep = "_"
    ),
    region_layer_ExpGroup = paste(
      primary_audit$region,
      primary_audit$layer,
      primary_audit$ExpGroup,
      sep = "_"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  annotated_table <- protigy_build_annotated_table(animal_expression, output_metadata)
  strict_audit <- primary_audit[
    primary_audit$hemisphere_status == "bilateral_complete",
    ,
    drop = FALSE
  ]
  strict_columns <- c("id", "Description", strict_audit$output_column_name)
  strict_bilateral_expression <- animal_expression[, strict_columns, drop = FALSE]
  strict_bilateral_metadata <- output_metadata[
    match(strict_audit$output_column_name, output_metadata$sample_id),
    ,
    drop = FALSE
  ]
  rownames(strict_bilateral_metadata) <- NULL
  strict_bilateral_annotated_table <- protigy_build_annotated_table(
    strict_bilateral_expression,
    strict_bilateral_metadata
  )

  status <- if (length(validation_issues)) {
    paste0("FAIL: ", paste(unique(validation_issues), collapse = " | "))
  } else {
    "PASS"
  }
  result <- list(
    dataset = dataset,
    spatial_unit = spatial_unit,
    n_source_expression_rows = nrow(expression),
    n_proteins = length(feature_ids),
    n_input_samples = length(included_samples),
    feature_ids = feature_ids,
    feature_descriptions = feature_descriptions,
    description_source = if (length(description_source)) {
      description_source
    } else {
      "id_compatibility_fallback"
    },
    description_fallback_count = as.integer(description_fallback_count),
    animal_expression = animal_expression,
    output_metadata = output_metadata,
    annotated_table = annotated_table,
    strict_bilateral_expression = strict_bilateral_expression,
    strict_bilateral_metadata = strict_bilateral_metadata,
    strict_bilateral_annotated_table = strict_bilateral_annotated_table,
    aggregation_audit = aggregation_audit,
    source_sample_assignment = source_assignment,
    validation_status = status,
    validation_issues = unique(validation_issues),
    valid = !length(validation_issues)
  )
  if (isTRUE(fail_on_invalid) && !isTRUE(result$valid)) {
    stop(result$validation_status, call. = FALSE)
  }
  result
}
