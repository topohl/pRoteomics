# Shared, raw-derived helpers for global joint compartment QC.
# These helpers deliberately keep biological dataset labels out of normalization
# and imputation; labels are used only for documented feature eligibility/audits.

if (!exists("repo_path", mode = "function")) source(file.path("R", "paths.R"))
if (!exists("valid_datasets", mode = "function")) source(repo_path("R", "dataset_config.R"))
if (!exists("build_canonical_protein_group_tables", mode = "function")) source(repo_path("R", "protein_mapping_utils.R"))
if (!exists("qc_build_mapping_context", mode = "function")) source(repo_path("R", "qc_exploration_utils.R"))

joint_qc_datasets <- function() valid_datasets()

joint_qc_env_number <- function(name, default, lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value) || value < lower || value > upper) {
    stop(name, " must be a finite number in [", lower, ", ", upper, "].", call. = FALSE)
  }
  value
}

joint_qc_first_col <- function(df, candidates, required = FALSE, context = "table") {
  norm <- function(x) tolower(gsub("[^a-z0-9]", "", x))
  hit <- match(norm(candidates), norm(names(df)))
  hit <- hit[!is.na(hit)]
  if (length(hit)) return(names(df)[hit[[1]]])
  if (required) stop("Missing required column in ", context, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  NA_character_
}

joint_qc_truthy <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x %in% c("1", "true", "t", "yes", "y", "excluded", "exclude")
}

joint_qc_normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\\\", "/", x)
  tolower(x)
}

joint_qc_metadata_dataset <- function(meta) {
  celltype_layer <- joint_qc_first_col(meta, c("dataset", "celltype_layer", "CellTypeLayer"))
  celltype <- joint_qc_first_col(meta, c("celltype", "CellType"))
  layer <- joint_qc_first_col(meta, c("layer", "Layer"))
  out <- rep(NA_character_, nrow(meta))
  token <- function(x) tolower(gsub("[[:space:]-]+", "_", trimws(as.character(x))))
  if (!is.na(celltype_layer)) out <- normalize_dataset(meta[[celltype_layer]])
  if (!is.na(celltype)) {
    unresolved <- is.na(out) | !out %in% joint_qc_datasets()
    candidate <- normalize_dataset(meta[[celltype]])
    out[unresolved & candidate %in% joint_qc_datasets()] <- candidate[unresolved & candidate %in% joint_qc_datasets()]
  }
  if (!is.na(layer)) {
    lay <- token(meta[[layer]])
    unresolved <- is.na(out) | !out %in% joint_qc_datasets()
    out[unresolved & lay %in% c("so", "sr", "slm", "mo", "po")] <- "neuron_neuropil"
    out[unresolved & lay %in% c("sp", "sg")] <- "neuron_soma"
  }
  out
}

joint_qc_prepare_metadata <- function(metadata, matrix_sample_ids) {
  sample_col <- joint_qc_first_col(metadata, c("sample_id", "Sample", "sample", "SampleID"), TRUE, "sample metadata")
  metadata$.joint_sample_id <- as.character(metadata[[sample_col]])
  metadata$.joint_sample_key <- joint_qc_normalize_sample_id(metadata$.joint_sample_id)
  matrix_key <- joint_qc_normalize_sample_id(matrix_sample_ids)
  if (anyDuplicated(matrix_key)) stop("Duplicated matrix sample IDs after exact canonical normalization.", call. = FALSE)
  if (anyDuplicated(metadata$.joint_sample_key)) stop("Duplicated metadata sample IDs after exact canonical normalization.", call. = FALSE)
  metadata$dataset <- joint_qc_metadata_dataset(metadata)
  exclude_col <- joint_qc_first_col(metadata, c("exclude", "Exclude", "excluded", "QC_status"))
  metadata$exclusion_status <- if (is.na(exclude_col)) "included" else ifelse(joint_qc_truthy(metadata[[exclude_col]]), "excluded", "included")
  background_cols <- intersect(c("sample_type", "SampleType", "celltype", "CellType", "celltype_layer", "CellTypeLayer", "group", "Group"), names(metadata))
  background_text <- if (length(background_cols)) apply(metadata[background_cols], 1, function(z) paste(z, collapse = " ")) else ""
  metadata$background_or_blank <- grepl("blank|background", background_text, ignore.case = TRUE)
  idx <- match(matrix_key, metadata$.joint_sample_key)
  unmatched_matrix <- matrix_sample_ids[is.na(idx)]
  unmatched_metadata <- metadata$.joint_sample_id[!metadata$.joint_sample_key %in% matrix_key]
  if (length(unmatched_matrix) || length(unmatched_metadata)) {
    stop("Matrix/metadata matching is not one-to-one. Unmatched matrix samples: ", length(unmatched_matrix), "; unmatched metadata rows: ", length(unmatched_metadata), ".", call. = FALSE)
  }
  metadata <- metadata[idx, , drop = FALSE]
  rownames(metadata) <- matrix_sample_ids
  metadata$Sample <- matrix_sample_ids
  metadata
}

joint_qc_select_technical_block <- function(metadata) {
  candidates <- c("batch", "acquisition_batch", "plate", "run", "acquisition_date", "run_date")
  col <- joint_qc_first_col(metadata, candidates)
  if (!is.na(col)) {
    value <- trimws(as.character(metadata[[col]]))
    value[is.na(value) | !nzchar(value)] <- NA_character_
    if (sum(!is.na(value)) == nrow(metadata) && length(unique(value)) >= 2L) {
      return(list(variable = col, block = value, mode = "metadata_technical_block", warning = NA_character_))
    }
  }
  list(variable = "dataset", block = as.character(metadata$dataset), mode = "dataset_only_fallback",
       warning = "No complete technical batch/plate/run/date field with at least two levels; used dataset-only blocks.")
}

joint_qc_detection_long <- function(observed, metadata, protein_ids, technical_block) {
  groups <- unique(data.frame(dataset = metadata$dataset, technical_block = technical_block, stringsAsFactors = FALSE))
  out <- lapply(seq_len(nrow(groups)), function(i) {
    keep <- metadata$dataset == groups$dataset[i] & technical_block == groups$technical_block[i]
    data.frame(ProteinGroupID = protein_ids, dataset = groups$dataset[i], technical_block = groups$technical_block[i],
               n_samples = sum(keep), observed_detection_rate_raw = rowMeans(observed[, keep, drop = FALSE]), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

joint_qc_filter_universes <- function(observed, metadata, protein_ids, technical_block, min_block = .70, min_union = .30) {
  block_long <- joint_qc_detection_long(observed, metadata, protein_ids, technical_block)
  dataset_long <- do.call(rbind, lapply(joint_qc_datasets(), function(ds) {
    keep <- metadata$dataset == ds
    data.frame(ProteinGroupID = protein_ids, dataset = ds, n_samples = sum(keep), observed_detection_rate_raw = rowMeans(observed[, keep, drop = FALSE]), stringsAsFactors = FALSE)
  }))
  by_id <- split(block_long, block_long$ProteinGroupID)
  primary <- vapply(protein_ids, function(id) all(by_id[[id]]$observed_detection_rate_raw >= min_block), logical(1))
  complete <- rowSums(observed) == ncol(observed)
  by_dataset <- split(dataset_long, dataset_long$ProteinGroupID)
  broad <- vapply(protein_ids, function(id) any(by_dataset[[id]]$observed_detection_rate_raw >= min_union), logical(1))
  audit <- data.frame(ProteinGroupID = protein_ids, primary_shared_core_inclusion = primary,
                      complete_case_inclusion = complete, broad_union_inclusion = broad,
                      n_missing_before_imputation = rowSums(!observed), fraction_missing_before_imputation = rowMeans(!observed),
                      exclusion_reason = ifelse(primary, NA_character_, "below_min_detection_in_at_least_one_dataset_x_technical_block"), stringsAsFactors = FALSE)
  for (ds in joint_qc_datasets()) {
    x <- dataset_long[dataset_long$dataset == ds, c("ProteinGroupID", "observed_detection_rate_raw")]
    audit[[paste0("observed_detection_rate_raw_", ds)]] <- x$observed_detection_rate_raw[match(protein_ids, x$ProteinGroupID)]
  }
  list(primary = primary, complete = complete, broad = broad, audit = audit, block_long = block_long, dataset_long = dataset_long)
}

joint_qc_report_feature_subset_diagnostics <- function(ids, mat, label) {
  matrix_ids <- rownames(mat)
  fmt <- function(x) if (is.null(x)) "<NULL>" else paste(head(x, 20L), collapse = " | ")
  message(label, " length(primary_ids): ", length(ids))
  message(label, " nrow(matrix): ", nrow(mat))
  message(label, " head(requested IDs): ", fmt(ids))
  message(label, " head(matrix row names): ", fmt(matrix_ids))
  message(label, " requested IDs in matrix: ", sum(ids %in% matrix_ids))
  message(label, " requested IDs absent from matrix: ", fmt(setdiff(ids, matrix_ids)))
  message(label, " matrix IDs absent from requested IDs: ", fmt(setdiff(matrix_ids, ids)))
  message(label, " anyDuplicated(requested IDs): ", anyDuplicated(ids))
  message(label, " anyDuplicated(matrix row names): ", if (is.null(matrix_ids)) "<NULL>" else anyDuplicated(matrix_ids))
  invisible(NULL)
}

joint_validate_feature_subset <- function(ids, mat, label) {
  ids <- as.character(ids)
  matrix_ids <- rownames(mat)
  if (is.null(matrix_ids)) stop(label, ": matrix has no row names.", call. = FALSE)
  if (anyNA(ids) || any(!nzchar(ids))) stop(label, ": requested feature IDs contain missing or blank values.", call. = FALSE)
  if (anyDuplicated(ids)) stop(label, ": requested feature IDs are duplicated.", call. = FALSE)
  if (anyNA(matrix_ids) || any(!nzchar(matrix_ids)) || anyDuplicated(matrix_ids)) stop(label, ": matrix row identifiers are invalid or duplicated.", call. = FALSE)
  missing_ids <- setdiff(ids, matrix_ids)
  if (length(missing_ids)) stop(label, ": ", length(missing_ids), " requested ProteinGroupID values are absent from the matrix. Examples: ", paste(head(missing_ids, 10L), collapse = ", "), call. = FALSE)
  idx <- match(ids, matrix_ids)
  if (anyNA(idx)) stop(label, ": internal feature matching failed.", call. = FALSE)
  out <- mat[idx, , drop = FALSE]
  if (!identical(rownames(out), ids)) stop(label, ": selected matrix row names do not exactly match requested ProteinGroupID order.", call. = FALSE)
  out
}

joint_validate_canonical_feature_matrix <- function(feature_table, mat, label) {
  if (!"ProteinGroupID" %in% names(feature_table)) stop(label, ": feature table lacks ProteinGroupID.", call. = FALSE)
  feature_ids <- as.character(feature_table$ProteinGroupID)
  if (nrow(mat) != length(feature_ids)) stop(label, ": matrix row count (", nrow(mat), ") differs from canonical feature-table row count (", length(feature_ids), ").", call. = FALSE)
  if (anyNA(feature_ids) || any(!nzchar(feature_ids)) || anyDuplicated(feature_ids)) stop(label, ": canonical feature-table ProteinGroupID values are invalid or duplicated.", call. = FALSE)
  matrix_ids <- rownames(mat)
  if (is.null(matrix_ids) || anyNA(matrix_ids) || any(!nzchar(matrix_ids)) || anyDuplicated(matrix_ids)) stop(label, ": matrix row identifiers are missing, invalid, or duplicated.", call. = FALSE)
  if (!identical(matrix_ids, feature_ids)) stop(label, ": matrix row names are not identical to canonical ProteinGroupID source-row order.", call. = FALSE)
  invisible(TRUE)
}

joint_qc_log2_positive <- function(raw) {
  if (is.null(dim(raw))) stop("Raw intensity input must be a matrix.", call. = FALSE)
  raw <- matrix(
    suppressWarnings(as.numeric(raw)), nrow = nrow(raw), ncol = ncol(raw),
    dimnames = dimnames(raw)
  )
  if (is.null(rownames(raw)) || anyNA(rownames(raw)) || any(!nzchar(rownames(raw))) || anyDuplicated(rownames(raw))) {
    stop("Raw intensity matrix must carry unique non-blank canonical row names before positive-value filtering.", call. = FALSE)
  }
  zero <- is.finite(raw) & raw == 0
  negative <- is.finite(raw) & raw < 0
  if (any(negative)) stop("Raw intensity matrix contains negative values; log2 transformation is invalid.", call. = FALSE)
  raw[zero] <- NA_real_
  log2_matrix <- log2(raw)
  rownames(log2_matrix) <- rownames(raw)
  colnames(log2_matrix) <- colnames(raw)
  list(raw = raw, log2 = log2_matrix, zero = zero, nonfinite = !is.finite(raw) & !is.na(raw))
}

joint_qc_joint_median_normalize <- function(log2_matrix) {
  med <- apply(log2_matrix, 2, stats::median, na.rm = TRUE)
  if (any(!is.finite(med))) stop("At least one retained sample has no positive observed intensity.", call. = FALSE)
  target <- stats::median(med)
  normalized <- sweep(log2_matrix, 2, med - target, FUN = "-")
  list(matrix = normalized, pre_median = med, post_median = apply(normalized, 2, stats::median, na.rm = TRUE), target_median = target)
}

joint_qc_median_impute <- function(matrix, protein_ids = rownames(matrix)) {
  med <- apply(matrix, 1, stats::median, na.rm = TRUE)
  if (any(!is.finite(med))) stop("Protein-wise median imputation encountered a protein with no observed values.", call. = FALSE)
  miss <- is.na(matrix)
  out <- matrix
  out[miss] <- med[row(out)[miss]]
  list(matrix = out, protein_median = med, missing = miss,
       audit = data.frame(ProteinGroupID = protein_ids, n_imputed = rowSums(miss), fraction_imputed = rowMeans(miss), stringsAsFactors = FALSE))
}

joint_qc_write_matrix_tsv <- function(matrix, path) {
  dir_create(dirname(path))
  utils::write.table(data.frame(ProteinGroupID = rownames(matrix), matrix, check.names = FALSE), path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

joint_qc_write_gct_v13 <- function(matrix, metadata, feature_table, path) {
  if (any(!is.finite(matrix))) stop("GCT export requires a complete finite matrix.", call. = FALSE)
  if (!identical(colnames(matrix), as.character(metadata$Sample))) stop("GCT matrix and metadata columns are not aligned.", call. = FALSE)
  desc_col <- if ("FeatureDisplayLabel" %in% names(feature_table)) "FeatureDisplayLabel" else "original_identifier"
  row_meta <- data.frame(ProteinGroupID = rownames(matrix), Description = as.character(feature_table[[desc_col]][match(rownames(matrix), feature_table$ProteinGroupID)]), stringsAsFactors = FALSE)
  row_meta$Description[is.na(row_meta$Description)] <- ""
  meta_cols <- intersect(c("dataset", "compartment", "cell_type", "celltype", "region", "layer", "region_layer", "experimental_group", "ExpGroup", "group", "animal_id", "AnimalID", "hemisphere", "batch", "plate", "run", "technical_block", "exclusion_status"), names(metadata))
  dir_create(dirname(path)); con <- file(path, "wt", encoding = "UTF-8"); on.exit(close(con), add = TRUE)
  writeLines("#1.3", con)
  writeLines(paste(nrow(matrix), ncol(matrix), ncol(row_meta), length(meta_cols), sep = "\t"), con)
  writeLines(paste(c(names(row_meta), colnames(matrix)), collapse = "\t"), con)
  for (nm in meta_cols) writeLines(paste(c(nm, rep("", ncol(row_meta) - 1L), as.character(metadata[[nm]])), collapse = "\t"), con)
  for (i in seq_len(nrow(matrix))) writeLines(paste(c(as.character(row_meta[i, ]), format(matrix[i, ], scientific = FALSE, trim = TRUE, digits = 15)), collapse = "\t"), con)
  close(con); on.exit(NULL, add = FALSE)
  joint_qc_validate_gct_v13(path, nrow(matrix), ncol(matrix), ncol(row_meta), length(meta_cols))
  invisible(path)
}

joint_qc_validate_gct_v13 <- function(path, n_rows = NULL, n_cols = NULL, n_row_meta = NULL, n_col_meta = NULL) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 3L || lines[[1]] != "#1.3") stop("Invalid GCT v1.3 header.", call. = FALSE)
  dims <- suppressWarnings(as.integer(strsplit(lines[[2]], "\t", fixed = TRUE)[[1]]))
  if (length(dims) != 4L || anyNA(dims)) stop("Invalid GCT v1.3 dimension line.", call. = FALSE)
  expected <- c(n_rows, n_cols, n_row_meta, n_col_meta)
  if (all(!is.null(expected)) && any(!is.na(expected)) && !identical(dims, as.integer(expected))) stop("GCT v1.3 dimensions do not match written data.", call. = FALSE)
  data_lines <- lines[(4L + dims[[4]]):(3L + dims[[4]] + dims[[1]])]
  numeric_start <- dims[[3]] + 1L
  values <- suppressWarnings(as.numeric(unlist(lapply(strsplit(data_lines, "\t", fixed = TRUE), function(x) x[numeric_start:(numeric_start + dims[[2]] - 1L)]))))
  if (length(values) != dims[[1]] * dims[[2]] || any(!is.finite(values))) stop("GCT v1.3 numeric contents are invalid.", call. = FALSE)
  invisible(TRUE)
}

joint_qc_observed_detection_provenance <- function(protein_ids, dataset, root = path_processed("01_preprocessing", "joint_compartment_qc", "global")) {
  file <- file.path(root, "observed_detection_by_dataset.csv")
  out <- data.frame(ProteinGroupID = protein_ids, observed_detection_rate_raw = NA_real_, detectability_source = "raw_unified_output_not_available", stringsAsFactors = FALSE)
  if (!file.exists(file)) return(out)
  tab <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  tab <- tab[tab$dataset == dataset, c("ProteinGroupID", "observed_detection_rate_raw"), drop = FALSE]
  out$observed_detection_rate_raw <- tab$observed_detection_rate_raw[match(out$ProteinGroupID, tab$ProteinGroupID)]
  out$detectability_source <- ifelse(is.finite(out$observed_detection_rate_raw), "raw_unified_pre_imputation", "raw_unified_feature_not_present")
  out
}
