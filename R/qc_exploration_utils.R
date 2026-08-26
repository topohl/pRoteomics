# Shared helpers for dataset-aware QC and exploration scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("validate_dataset", mode = "function")) {
  source(repo_path("R", "dataset_config.R"))
}
if (!exists("resolve_dataset_inputs", mode = "function")) {
  source(repo_path("R", "dataset_inputs.R"))
}
if (!exists("build_canonical_protein_group_tables", mode = "function")) {
  source(repo_path("R", "protein_mapping_utils.R"))
}
if (!exists("wgcna_feature_key_fingerprint", mode = "function")) {
  source(repo_path("R", "module_contracts.R"))
}

qc_args <- function(default_dataset = "neuron_neuropil") {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag, default = "") {
    hit <- which(args == flag)
    if (!length(hit) || hit[[1]] == length(args)) return(default)
    args[[hit[[1]] + 1L]]
  }
  dataset_cli <- value_after("--dataset", "")
  if (nzchar(dataset_cli)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset_cli, source = "--dataset"))
  list(
    args = args,
    dataset = current_dataset(default_dataset),
    dry_run = is_dry_run(),
    run_embeddings = any(args %in% c("--run-embeddings", "--run-umap", "--run-tsne")) ||
      tolower(Sys.getenv("PROTEOMICS_QC_RUN_EMBEDDINGS", unset = "false")) %in% c("1", "true", "yes", "y"),
    run_clustering = "--run-clustering" %in% args ||
      tolower(Sys.getenv("PROTEOMICS_QC_RUN_CLUSTERING", unset = "false")) %in% c("1", "true", "yes", "y")
  )
}

qc_paths <- function(substep, dataset) {
  create_module_dirs("03_qc_exploration", file.path(substep, dataset))
}

qc_latest <- function(root, pattern) {
  latest_matching_file(root, pattern, recursive = TRUE)
}

qc_resolve_matrix <- function(dataset, env = "PROTEOMICS_QC_MATRIX_FILE") {
  override <- Sys.getenv(env, unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))

  dataset_inputs <- resolve_dataset_inputs(dataset, purpose = "wgcna")
  candidates <- c(
    dataset_inputs$expression_file,
    qc_latest(path_processed("01_preprocessing", "impute"), paste0("pgmatrix.*", dataset, ".*\\.(xlsx|csv|tsv)$")),
    qc_latest(path_processed("01_preprocessing"), paste0("pgmatrix.*", dataset, ".*\\.(xlsx|csv|tsv|gct)$")),
    qc_latest(path_processed("01_preprocessing", "excel_convert"), paste0(".*", dataset, ".*\\.(gct|xlsx|csv|tsv)$")),
    path_processed("01_preprocessing", "impute", paste0("pgmatrix_imputed_", dataset, ".xlsx"))
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) && !is.na(hit)) normalizePath(hit, winslash = "/", mustWork = FALSE) else candidates[1]
}

qc_resolve_metadata <- function(dataset, env = "PROTEOMICS_QC_METADATA_FILE") {
  override <- Sys.getenv(env, unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))

  dataset_inputs <- resolve_dataset_inputs(dataset, purpose = "wgcna")
  candidates <- c(
    dataset_inputs$metadata_file,
    path_metadata("sample_metadata_clean.tsv"),
    path_metadata("sample_metadata_clean.csv"),
    path_metadata("TPE9_sample_metadata_males.xlsx")
  )
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) && !is.na(hit)) normalizePath(hit, winslash = "/", mustWork = FALSE) else candidates[1]
}

qc_read_table <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required to read: ", path, call. = FALSE)
    return(as.data.frame(readxl::read_excel(path, sheet = sheet %||% 1), check.names = FALSE))
  }
  if (ext == "tsv") {
    return(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (ext == "csv") {
    return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  stop("Unsupported tabular file extension for: ", path, call. = FALSE)
}

qc_first_col <- function(df, candidates) {
  norm <- function(x) tolower(gsub("[^a-z0-9]", "", x))
  hit <- match(norm(candidates), norm(names(df)))
  hit <- hit[!is.na(hit)]
  if (!length(hit)) return(NA_character_)
  names(df)[hit[[1]]]
}

qc_sample_metadata_from_names <- function(samples) {
  data.frame(
    Sample = samples,
    sample_id = samples,
    Region = ifelse(grepl("CA1", samples, ignore.case = TRUE), "CA1",
      ifelse(grepl("CA2", samples, ignore.case = TRUE), "CA2",
        ifelse(grepl("CA3", samples, ignore.case = TRUE), "CA3",
          ifelse(grepl("DG", samples, ignore.case = TRUE), "DG", NA_character_)))),
    Layer = ifelse(grepl("SLM", samples, ignore.case = TRUE), "SLM",
      ifelse(grepl("\\bSO\\b|_SO|SO_", samples, ignore.case = TRUE), "SO",
        ifelse(grepl("\\bSR\\b|_SR|SR_", samples, ignore.case = TRUE), "SR",
          ifelse(grepl("\\bSP\\b|_SP|SP_", samples, ignore.case = TRUE), "SP",
            ifelse(grepl("\\bSG\\b|_SG|SG_", samples, ignore.case = TRUE), "SG", NA_character_))))),
    stringsAsFactors = FALSE
  )
}

qc_make_feature_ids <- function(ids, prefix = "unannotated_protein") {
  # Legacy display-only helper for non-gene-aware diagnostics. These repaired
  # labels are never a quantitative identity and must not be joined to genes.
  ids <- trimws(as.character(ids))
  missing <- is.na(ids) | !nzchar(ids) | tolower(ids) %in% c("na", "nan")
  ids[missing] <- paste0(prefix, "_", which(missing))
  make.unique(ids)
}

qc_read_gct_v13 <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 4L || trimws(lines[[1]]) != "#1.3") stop("Expected strict GCT v1.3 file: ", path, call. = FALSE)
  dims <- suppressWarnings(as.integer(strsplit(lines[[2]], "\t")[[1]]))
  if (length(dims) < 4L || anyNA(dims[1:4])) stop("Invalid GCT v1.3 dimension line in: ", path, call. = FALSE)
  n_rows <- dims[[1]]
  n_samples <- dims[[2]]
  n_row_meta <- dims[[3]]
  n_col_meta <- dims[[4]]
  header <- strsplit(lines[[3]], "\t", fixed = TRUE)[[1]]
  sample_ids <- header[(n_row_meta + 1L):(n_row_meta + n_samples)]
  meta <- qc_sample_metadata_from_names(sample_ids)
  if (n_col_meta > 0L) {
    meta_lines <- lines[4L:(3L + n_col_meta)]
    for (ln in meta_lines) {
      parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      key <- parts[[1]]
      vals <- parts[(n_row_meta + 1L):(n_row_meta + n_samples)]
      meta[[key]] <- vals
    }
  }
  expr_lines <- lines[(4L + n_col_meta):(3L + n_col_meta + n_rows)]
  split <- strsplit(expr_lines, "\t", fixed = TRUE)
  ids <- vapply(split, `[`, character(1), 1L)
  mat <- do.call(rbind, lapply(split, function(x) suppressWarnings(as.numeric(x[(n_row_meta + 1L):(n_row_meta + n_samples)]))))
  rownames(mat) <- qc_make_feature_ids(ids)
  colnames(mat) <- sample_ids
  list(mat = mat, meta = meta, source = "gct")
}

qc_read_expression <- function(matrix_file, metadata_file = NA_character_, dataset = current_dataset()) {
  ext <- tolower(tools::file_ext(matrix_file))
  if (ext == "gct") {
    parsed <- qc_read_gct_v13(matrix_file)
    mat <- parsed$mat
    meta <- parsed$meta
  } else {
    df <- qc_read_table(matrix_file)
    id_col <- qc_first_col(df, c("Genes", "gene_symbol", "Protein.Group", "ProteinID", "Protein", "id", "T: Protein.Names"))
    if (is.na(id_col)) id_col <- names(df)[[1]]
    numeric_cols <- names(df)[vapply(df, function(x) {
      suppressWarnings(mean(!is.na(as.numeric(as.character(x)))) > 0.5)
    }, logical(1))]
    sample_cols <- setdiff(numeric_cols, id_col)
    if (length(sample_cols) < 2L) stop("Could not detect at least two numeric sample columns in: ", matrix_file, call. = FALSE)
    mat <- as.matrix(data.frame(lapply(df[sample_cols], function(x) as.numeric(as.character(x))), check.names = FALSE))
    rownames(mat) <- qc_make_feature_ids(df[[id_col]])
    colnames(mat) <- sample_cols
    meta <- qc_sample_metadata_from_names(sample_cols)
  }

  if (!is.na(metadata_file) && file.exists(metadata_file)) {
    meta0 <- qc_read_table(metadata_file)
    sample_col <- qc_first_col(meta0, c("Sample", "sample", "sample_id", "SampleID", "SampleColumn", "row.names", "shortname", "sampleNumber"))
    if (!is.na(sample_col)) {
      meta0[[sample_col]] <- as.character(meta0[[sample_col]])
      keep <- match(colnames(mat), meta0[[sample_col]])
      if (sum(!is.na(keep)) > 0L) {
        merged <- meta
        meta_match <- meta0[keep, , drop = FALSE]
        for (nm in names(meta_match)) merged[[nm]] <- meta_match[[nm]]
        meta <- merged
      }
    }
  }
  keep_meta <- metadata_matches_dataset(meta, dataset)
  if (any(keep_meta) && sum(keep_meta) < nrow(meta)) {
    mat <- mat[, keep_meta, drop = FALSE]
    meta <- meta[keep_meta, , drop = FALSE]
  }
  rownames(meta) <- colnames(mat)
  meta$Sample <- colnames(mat)
  list(
    mat = mat,
    meta = meta,
    feature_identity_contract = "legacy_display_labels_only",
    source_row_provenance = data.frame(source_row_id = seq_len(nrow(mat)), display_label = rownames(mat))
  )
}

qc_detect_sample_columns <- function(df) {
  annotation_columns <- unique(c(
    detect_source_provenance_columns(df),
    "First.Protein.Description", "Description", ".row_id", "source_row_id",
    "member_identifiers_original", "member_identifiers_canonical",
    "member_accessions", "member_gene_symbols", "representative_accession",
    "representative_gene_symbol", "protein_group_ambiguity_class",
    "gene_level_claim_allowed", "protein_level_claim_allowed", "mapping_status"
  ))
  candidates <- setdiff(names(df), annotation_columns)
  candidates[vapply(df[candidates], function(x) {
    vals <- suppressWarnings(as.numeric(as.character(x)))
    mean(!is.na(vals)) >= 0.7
  }, logical(1))]
}

qc_build_mapping_context <- function(idmap_file = path_external("MOUSE_10090_idmapping.dat"),
                                     manual_mapping_file = Sys.getenv(
                                       "PROTEOMICS_MANUAL_MAPPING_FILE",
                                       unset = path_metadata("manual_mapping.xlsx")
                                     ),
                                     manual_gene_annotation_file = Sys.getenv(
                                       "PROTEOMICS_MANUAL_GENE_ANNOTATION_FILE",
                                       unset = path_metadata("manual_gene_annotation_overrides.csv")
                                     )) {
  required <- c("dplyr", "AnnotationDbi", "org.Mm.eg.db")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Canonical QC mapping requires package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  idmap <- load_mouse_idmapping(idmap_file)
  maps <- build_mouse_maps(idmap)
  orgdb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
  annotation_maps <- build_mouse_gene_annotation_maps(
    orgdb,
    accessions = unique(maps$accession_gene_map$UNIPROT)
  )
  list(
    entry_map = maps$entry_map,
    gene_map = maps$gene_map,
    accession_gene_map = maps$accession_gene_map,
    reviewed_map = maps$reviewed_map,
    manual_mapping = read_manual_mapping_table(manual_mapping_file),
    gene_annotation_maps = annotation_maps,
    manual_gene_annotation_overrides = read_manual_gene_annotation_overrides(manual_gene_annotation_file),
    idmap_file = normalizePath(idmap_file, winslash = "/", mustWork = FALSE),
    idmap_sha256 = file_hash_sha256(idmap_file),
    manual_mapping_file = normalizePath(manual_mapping_file, winslash = "/", mustWork = FALSE),
    manual_gene_annotation_file = normalizePath(manual_gene_annotation_file, winslash = "/", mustWork = FALSE)
  )
}

qc_align_sample_metadata <- function(mat, metadata_file = NA_character_, dataset = current_dataset()) {
  meta <- qc_sample_metadata_from_names(colnames(mat))
  if (!is.na(metadata_file) && file.exists(metadata_file)) {
    meta0 <- qc_read_table(metadata_file)
    sample_col <- qc_first_col(meta0, c(
      "Sample", "sample", "sample_id", "SampleID", "SampleColumn", "row.names",
      "shortname", "sampleNumber"
    ))
    if (!is.na(sample_col)) {
      meta0[[sample_col]] <- as.character(meta0[[sample_col]])
      keep <- match(colnames(mat), meta0[[sample_col]])
      if (any(!is.na(keep))) {
        matched <- meta0[keep, , drop = FALSE]
        for (nm in names(matched)) meta[[nm]] <- matched[[nm]]
      }
    }
  }
  keep <- metadata_matches_dataset(meta, dataset)
  if (any(keep) && sum(keep) < nrow(meta)) {
    mat <- mat[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
  }
  rownames(meta) <- colnames(mat)
  meta$Sample <- colnames(mat)
  list(mat = mat, meta = meta)
}

qc_validate_canonical_expression <- function(mat, feature_table) {
  if (!"ProteinGroupID" %in% names(feature_table)) {
    stop("Canonical QC feature table lacks ProteinGroupID.", call. = FALSE)
  }
  ids <- as.character(feature_table$ProteinGroupID)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Canonical QC feature table has a missing or duplicate ProteinGroupID; label repair is forbidden.", call. = FALSE)
  }
  if (!identical(rownames(mat), ids)) {
    stop("Canonical QC matrix and feature table are not aligned by ordered ProteinGroupID.", call. = FALSE)
  }
  invisible(TRUE)
}

qc_validate_source_row_preservation <- function(source_rows, feature_table) {
  if (nrow(source_rows) != nrow(feature_table)) {
    stop("Canonical QC construction lost quantitative source rows before eligibility filtering.", call. = FALSE)
  }
  if (!"source_row_id" %in% names(feature_table) ||
      !identical(as.integer(feature_table$source_row_id), seq_len(nrow(source_rows)))) {
    stop("Canonical QC source-row provenance is incomplete or reordered.", call. = FALSE)
  }
  invisible(TRUE)
}

qc_feature_universe_reconciliation <- function(feature_table, dataset) {
  data.frame(
    dataset = dataset,
    quantitative_features = nrow(feature_table),
    unique_ProteinGroupID = length(unique(feature_table$ProteinGroupID)),
    wgcna_eligible_features = sum(feature_table$wgcna_eligible %in% TRUE),
    marker_eligible_features = sum(feature_table$marker_eligible %in% TRUE),
    enrichment_eligible_features = sum(feature_table$enrichment_eligible %in% TRUE),
    multi_member_groups = sum(feature_table$n_members_canonical > 1L),
    missing_legacy_display_identifiers = sum(is.na(feature_table$original_identifier) | !nzchar(feature_table$original_identifier)),
    stringsAsFactors = FALSE
  )
}

qc_load_canonical_expression <- function(matrix_file, metadata_file = NA_character_,
                                         dataset = current_dataset(), mapping_context = NULL,
                                         sample_columns = NULL, strict = TRUE) {
  if (!file.exists(matrix_file)) stop("Canonical QC expression source not found: ", matrix_file, call. = FALSE)
  dataset <- validate_dataset(dataset)
  raw <- qc_read_table(matrix_file)
  raw$.source_row_id <- seq_len(nrow(raw))
  if (is.null(sample_columns)) sample_columns <- qc_detect_sample_columns(raw)
  if (length(sample_columns) < 2L) {
    stop("Could not detect at least two quantitative sample columns in canonical QC source: ", matrix_file, call. = FALSE)
  }
  if (is.null(mapping_context)) mapping_context <- qc_build_mapping_context()

  identifier_col <- if ("original_identifier" %in% names(raw)) {
    "original_identifier"
  } else if ("T: Protein.Names" %in% names(raw)) {
    "T: Protein.Names"
  } else if ("gene_symbol" %in% names(raw)) {
    "gene_symbol"
  } else if ("Protein.Group" %in% names(raw)) {
    "Protein.Group"
  } else names(raw)[[1]]
  feature_col <- if ("ProteinGroupID" %in% names(raw)) {
    "ProteinGroupID"
  } else if ("source_feature_id" %in% names(raw)) {
    "source_feature_id"
  } else if ("Protein.Group" %in% names(raw)) {
    "Protein.Group"
  } else NA_character_

  canonical <- build_canonical_protein_group_tables(
    df_raw = raw,
    dataset = dataset,
    source_file = matrix_file,
    entry_map = mapping_context$entry_map,
    gene_map = mapping_context$gene_map,
    accession_gene_map = mapping_context$accession_gene_map,
    reviewed_map = mapping_context$reviewed_map,
    manual_mapping = mapping_context$manual_mapping,
    gene_annotation_maps = mapping_context$gene_annotation_maps,
    manual_gene_annotation_overrides = mapping_context$manual_gene_annotation_overrides,
    uniprot_mapping_file_hash = mapping_context$idmap_sha256 %||% NA_character_,
    strict = strict,
    identifier_col = identifier_col,
    feature_col = feature_col
  )

  canonical_columns <- setdiff(names(canonical$wide), c(sample_columns, ".source_row_id"))
  feature_table <- canonical$wide[, canonical_columns, drop = FALSE]
  qc_validate_source_row_preservation(raw, feature_table)
  feature_table$FeatureDisplayLabel <- wgcna_feature_display_label(feature_table)
  feature_table$wgcna_eligible <- feature_table$protein_group_ambiguity_class != "mixed_species_or_contaminant"
  feature_table$wgcna_exclusion_reason <- dplyr::case_when(
    feature_table$wgcna_eligible ~ NA_character_,
    nzchar(feature_table$source_provenance_exclusion_reason) ~ feature_table$source_provenance_exclusion_reason,
    TRUE ~ "mixed_species_or_contaminant_member_evidence"
  )

  eligible_member <- canonical$bridge |>
    dplyr::mutate(
      eligible_official_mouse_member = .data$gene_annotation_status == "resolved" &
        !is.na(.data$member_gene_symbol) & nzchar(.data$member_gene_symbol) &
        .data$contaminant_status != "contaminant" &
        (is.na(.data$member_species) | .data$member_species == "mouse")
    ) |>
    dplyr::group_by(.data$ProteinGroupID) |>
    dplyr::summarise(has_eligible_official_mouse_member = any(.data$eligible_official_mouse_member), .groups = "drop")
  feature_table <- feature_table |>
    dplyr::left_join(eligible_member, by = "ProteinGroupID") |>
    dplyr::mutate(
      has_eligible_official_mouse_member = dplyr::coalesce(.data$has_eligible_official_mouse_member, FALSE),
      annotation_source_identifier_eligible = !is.na(.data$original_identifier) & nzchar(trimws(.data$original_identifier)),
      marker_eligible = .data$wgcna_eligible & .data$annotation_source_identifier_eligible &
        .data$has_eligible_official_mouse_member,
      marker_exclusion_reason = dplyr::case_when(
        .data$marker_eligible ~ NA_character_,
        !.data$wgcna_eligible ~ .data$wgcna_exclusion_reason,
        !.data$annotation_source_identifier_eligible ~ "blank_legacy_display_identifier_not_gene_claim_eligible",
        TRUE ~ "no_eligible_official_mouse_member"
      ),
      enrichment_eligible = .data$wgcna_eligible & .data$annotation_source_identifier_eligible &
        .data$gene_level_claim_allowed,
      enrichment_exclusion_reason = dplyr::case_when(
        .data$enrichment_eligible ~ NA_character_,
        !.data$wgcna_eligible ~ .data$wgcna_exclusion_reason,
        !.data$annotation_source_identifier_eligible ~ "blank_legacy_display_identifier_not_gene_claim_eligible",
        TRUE ~ "gene_level_claim_not_allowed"
      )
    )

  mat <- as.matrix(data.frame(lapply(raw[sample_columns], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  }), check.names = FALSE))
  rownames(mat) <- feature_table$ProteinGroupID
  colnames(mat) <- sample_columns
  aligned <- qc_align_sample_metadata(mat, metadata_file, dataset)
  mat <- aligned$mat
  meta <- aligned$meta
  qc_validate_canonical_expression(mat, feature_table)

  contract_version <- wgcna_feature_key_contract_version()
  feature_fingerprint <- wgcna_feature_key_fingerprint(feature_table$ProteinGroupID)
  wgcna_fingerprint <- wgcna_feature_key_fingerprint(feature_table$ProteinGroupID[feature_table$wgcna_eligible])
  source_path <- normalizePath(matrix_file, winslash = "/", mustWork = TRUE)
  source_sha256 <- file_hash_sha256(source_path)
  input_manifest <- data.frame(
    input_role = c("post_filter_imputed_quantitative_matrix", "sample_metadata", "mouse_uniprot_idmapping", "manual_protein_mapping", "manual_gene_annotation"),
    resolved_path = c(
      source_path,
      normalizePath(metadata_file, winslash = "/", mustWork = FALSE),
      mapping_context$idmap_file %||% NA_character_,
      mapping_context$manual_mapping_file %||% NA_character_,
      mapping_context$manual_gene_annotation_file %||% NA_character_
    ),
    sha256 = vapply(c(
      source_path, metadata_file,
      mapping_context$idmap_file %||% NA_character_,
      mapping_context$manual_mapping_file %||% NA_character_,
      mapping_context$manual_gene_annotation_file %||% NA_character_
    ), file_hash_sha256, character(1)),
    exists = file.exists(c(
      source_path, metadata_file,
      mapping_context$idmap_file %||% NA_character_,
      mapping_context$manual_mapping_file %||% NA_character_,
      mapping_context$manual_gene_annotation_file %||% NA_character_
    )),
    stringsAsFactors = FALSE
  )
  exclusion_audit <- dplyr::bind_rows(
    feature_table |> dplyr::transmute(ProteinGroupID, source_row_id, eligibility_layer = "wgcna", included = wgcna_eligible, exclusion_reason = wgcna_exclusion_reason),
    feature_table |> dplyr::transmute(ProteinGroupID, source_row_id, eligibility_layer = "marker", included = marker_eligible, exclusion_reason = marker_exclusion_reason),
    feature_table |> dplyr::transmute(ProteinGroupID, source_row_id, eligibility_layer = "enrichment", included = enrichment_eligible, exclusion_reason = enrichment_exclusion_reason)
  )
  source_row_provenance <- feature_table |>
    dplyr::select(dplyr::all_of(c(
      "ProteinGroupID", "source_file", "source_feature_id", "source_row_id",
      "original_identifier", "member_identifier_source", "member_identifiers_original",
      "source_provenance_columns", "source_provenance_values",
      "contaminant_source_evidence", "non_mouse_source_evidence"
    )))

  list(
    mat = mat,
    meta = meta,
    feature_table = feature_table,
    member_bridge = canonical$bridge,
    collision_audit = canonical$collision_audit,
    exclusion_audit = exclusion_audit,
    source_path = source_path,
    source_sha256 = source_sha256,
    contract_version = contract_version,
    feature_fingerprint = feature_fingerprint,
    wgcna_feature_fingerprint = wgcna_fingerprint,
    source_row_provenance = source_row_provenance,
    input_manifest = input_manifest,
    feature_universe_reconciliation = qc_feature_universe_reconciliation(feature_table, dataset),
    sample_columns = sample_columns
  )
}

qc_official_gene_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x)] <- ""
  x
}

qc_match_markers_to_protein_groups <- function(marker_registry, member_bridge, feature_table,
                                               panel_col = "marker_panel", gene_col = "marker_gene",
                                               class_col = NULL) {
  required_marker <- c(panel_col, gene_col)
  if (length(setdiff(required_marker, names(marker_registry)))) {
    stop("Marker registry lacks panel/gene columns required for canonical matching.", call. = FALSE)
  }
  qc_validate_canonical_expression(
    matrix(nrow = nrow(feature_table), ncol = 0L, dimnames = list(feature_table$ProteinGroupID, NULL)),
    feature_table
  )
  registry <- marker_registry |>
    dplyr::mutate(
      marker_panel = as.character(.data[[panel_col]]),
      requested_marker = as.character(.data[[gene_col]]),
      marker_gene_key = qc_official_gene_key(.data[[gene_col]])
    ) |>
    dplyr::filter(nzchar(.data$marker_panel), nzchar(.data$marker_gene_key)) |>
    dplyr::distinct(.data$marker_panel, .data$marker_gene_key, .keep_all = TRUE)
  eligible_ids <- feature_table$ProteinGroupID[feature_table$marker_eligible %in% TRUE]
  members <- member_bridge |>
    dplyr::filter(
      .data$ProteinGroupID %in% eligible_ids,
      .data$gene_annotation_status == "resolved",
      !is.na(.data$member_gene_symbol), nzchar(.data$member_gene_symbol),
      .data$contaminant_status != "contaminant",
      is.na(.data$member_species) | .data$member_species == "mouse"
    ) |>
    dplyr::mutate(marker_gene_key = qc_official_gene_key(.data$member_gene_symbol))
  expanded <- registry |>
    dplyr::inner_join(members, by = "marker_gene_key", relationship = "many-to-many") |>
    dplyr::distinct(.data$marker_panel, .data$ProteinGroupID, .data$marker_gene_key,
                    .data$member_identifier_original, .data$member_accession, .keep_all = TRUE)

  collapse_values <- function(x) paste(sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))), collapse = ";")
  matches <- expanded |>
    dplyr::group_by(.data$marker_panel, .data$ProteinGroupID) |>
    dplyr::summarise(
      requested_markers = collapse_values(.data$requested_marker),
      matched_official_genes = collapse_values(.data$member_gene_symbol),
      matched_member_accessions = collapse_values(.data$member_accession),
      matched_member_identifiers = collapse_values(.data$member_identifier_original),
      n_matched_official_genes = dplyr::n_distinct(.data$member_gene_symbol),
      n_matched_members = dplyr::n_distinct(.data$member_identifier_original),
      marker_classes = if (!is.null(class_col) && class_col %in% names(expanded)) collapse_values(.data[[class_col]]) else "",
      .groups = "drop"
    )
  conflicts <- matches |>
    dplyr::group_by(.data$ProteinGroupID) |>
    dplyr::summarise(
      n_marker_panels = dplyr::n_distinct(.data$marker_panel),
      matched_marker_panels = collapse_values(.data$marker_panel),
      n_marker_classes = length(unique(unlist(strsplit(.data$marker_classes[nzchar(.data$marker_classes)], ";", fixed = TRUE)))),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      conflicting_marker_panels = .data$n_marker_panels > 1L,
      conflicting_marker_classes = .data$n_marker_classes > 1L
    )
  matches <- matches |>
    dplyr::left_join(conflicts, by = "ProteinGroupID") |>
    dplyr::mutate(primary_score_eligible = !.data$conflicting_marker_classes)
  unmatched <- registry |>
    dplyr::anti_join(expanded |> dplyr::distinct(.data$marker_panel, .data$marker_gene_key),
                     by = c("marker_panel", "marker_gene_key"))
  list(matches = matches, expanded = expanded, conflicts = conflicts, unmatched = unmatched)
}

qc_validate_optional_wgcna_bridge <- function(state_file, canonical, strict = FALSE) {
  result <- data.frame(
    bridge_path = as.character(state_file),
    validation_status = "not_available",
    contract_version = canonical$contract_version,
    expected_feature_fingerprint = canonical$wgcna_feature_fingerprint,
    observed_feature_fingerprint = NA_character_,
    reason = "WGCNA state not found",
    stringsAsFactors = FALSE
  )
  if (is.na(state_file) || !file.exists(state_file)) return(result)
  state <- tryCatch(readRDS(state_file), error = function(e) e)
  if (inherits(state, "error")) {
    result$validation_status <- "rejected"
    result$reason <- conditionMessage(state)
  } else {
    error <- tryCatch({
      validate_wgcna_cached_state(
        state,
        expected_feature_ids = canonical$feature_table$ProteinGroupID[canonical$feature_table$wgcna_eligible]
      )
      NULL
    }, error = function(e) e)
    if (is.null(error)) {
      result$validation_status <- "validated"
      result$observed_feature_fingerprint <- state$feature_key_fingerprint
      result$reason <- NA_character_
    } else {
      result$validation_status <- "rejected"
      result$reason <- conditionMessage(error)
    }
  }
  if (isTRUE(strict) && result$validation_status != "validated") stop(result$reason, call. = FALSE)
  result
}

qc_add_input_manifest_paths <- function(canonical, paths) {
  if (!length(paths)) return(canonical)
  roles <- names(paths)
  if (is.null(roles) || any(!nzchar(roles))) stop("Additional canonical QC inputs require named roles.", call. = FALSE)
  paths <- as.character(paths)
  extra <- data.frame(
    input_role = roles,
    resolved_path = vapply(paths, normalizePath, character(1), winslash = "/", mustWork = FALSE),
    sha256 = vapply(paths, file_hash_sha256, character(1)),
    exists = file.exists(paths),
    stringsAsFactors = FALSE
  )
  canonical$input_manifest <- dplyr::bind_rows(canonical$input_manifest, extra) |>
    dplyr::distinct(.data$input_role, .data$resolved_path, .keep_all = TRUE)
  canonical
}

qc_write_canonical_feature_artifacts <- function(canonical, table_dir, marker_matches = NULL,
                                                 wgcna_validation = NULL) {
  qc_write_csv(canonical$feature_universe_reconciliation, file.path(table_dir, "feature_universe_reconciliation.csv"))
  qc_write_csv(canonical$feature_table, file.path(table_dir, "canonical_feature_table.csv"))
  qc_write_csv(canonical$member_bridge, file.path(table_dir, "protein_group_member_bridge.csv"))
  qc_write_csv(canonical$source_row_provenance, file.path(table_dir, "source_row_provenance.csv"))
  qc_write_csv(canonical$exclusion_audit, file.path(table_dir, "protein_group_exclusion_audit.csv"))
  qc_write_csv(
    canonical$feature_table |> dplyr::filter(.data$n_members_canonical > 1L),
    file.path(table_dir, "multi_member_group_audit.csv")
  )
  qc_write_csv(canonical$input_manifest, file.path(table_dir, "input_path_hash_manifest.csv"))
  qc_write_csv(data.frame(
    contract_version = canonical$contract_version,
    ordered_quantitative_feature_fingerprint = canonical$feature_fingerprint,
    ordered_wgcna_eligible_feature_fingerprint = canonical$wgcna_feature_fingerprint,
    n_quantitative_features = nrow(canonical$feature_table),
    n_wgcna_eligible_features = sum(canonical$feature_table$wgcna_eligible),
    stringsAsFactors = FALSE
  ), file.path(table_dir, "feature_contract.csv"))
  if (!is.null(marker_matches)) {
    qc_write_csv(marker_matches$matches, file.path(table_dir, "marker_to_ProteinGroupID_provenance.csv"))
    qc_write_csv(
      marker_matches$expanded,
      file.path(table_dir, "marker_to_ProteinGroupID_member_provenance.csv")
    )
    qc_write_csv(
      marker_matches$conflicts |>
        dplyr::filter(.data$conflicting_marker_panels) |>
        dplyr::left_join(
          canonical$feature_table |>
          dplyr::select(dplyr::all_of(c(
            "ProteinGroupID",
            "source_feature_id",
            "member_accessions",
            "member_gene_symbols",
            "official_gene_symbol"
          ))),
          by = "ProteinGroupID"
        ),
      file.path(table_dir, "conflicting_marker_panel_audit.csv")
    )
  }
  if (!is.null(wgcna_validation)) {
    qc_write_csv(wgcna_validation, file.path(table_dir, "wgcna_feature_contract_validation.csv"))
  }
  invisible(table_dir)
}

qc_write_csv <- function(x, path) {
  dir_create(dirname(path))
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

qc_write_xlsx <- function(x, path) {
  dir_create(dirname(path))
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(x, path)
  } else if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(x, path, overwrite = TRUE)
  } else {
    warning("No XLSX writer package installed; skipped: ", path, call. = FALSE)
  }
  invisible(path)
}

qc_embedding_palette <- function(n) {
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::hcl.colors(n, palette = "Dark 3")
}

qc_legend_rows <- function(x, max_per_row = 4L) {
  n <- length(unique(stats::na.omit(as.character(x))))
  max(1L, ceiling(n / max_per_row))
}

qc_embedding_theme <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1),
      legend.position = "bottom",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = grid::unit(3, "mm"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.spacing.x = grid::unit(1.5, "mm"),
      legend.margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0),
      plot.title = ggplot2::element_text(size = base_size, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, hjust = 0),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      aspect.ratio = 1
    )
}

qc_save_square_svg <- function(filename, plot, size_mm = 90) {
  ggplot2::ggsave(filename, plot, width = size_mm, height = size_mm, units = "mm", device = svglite::svglite)
  invisible(filename)
}

qc_safe_lm_p <- function(y, x) {
  ok <- is.finite(y) & !is.na(x) & nzchar(as.character(x))
  if (sum(ok) < 4L || length(unique(x[ok])) < 2L) return(NA_real_)
  fit <- try(stats::lm(y[ok] ~ factor(x[ok])), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  a <- stats::anova(fit)
  if (!"Pr(>F)" %in% names(a)) return(NA_real_)
  a[1, "Pr(>F)"]
}

qc_eta2 <- function(y, x) {
  ok <- is.finite(y) & !is.na(x) & nzchar(as.character(x))
  if (sum(ok) < 4L || length(unique(x[ok])) < 2L) return(NA_real_)
  fit <- try(stats::lm(y[ok] ~ factor(x[ok])), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  a <- stats::anova(fit)
  ss <- a[["Sum Sq"]]
  if (length(ss) < 2L || !is.finite(sum(ss))) return(NA_real_)
  ss[[1]] / sum(ss)
}

qc_metadata_terms <- function(meta) {
  candidates <- c("Group", "group", "ExpGroup", "group2", "Region", "region", "Layer", "layer",
    "AnimalID", "plate", "batch", "ReplicateGroup", "Sex", "sex", "run", "order", "sample_prep")
  intersect(candidates, names(meta))
}

qc_impute_for_pca <- function(mat) {
  mat[!is.finite(mat)] <- NA_real_
  row_missing <- rowMeans(is.na(mat))
  mat <- mat[row_missing < 0.8, , drop = FALSE]
  row_var <- apply(mat, 1L, stats::var, na.rm = TRUE)
  mat <- mat[is.finite(row_var) & row_var > 0, , drop = FALSE]
  med <- apply(mat, 1L, stats::median, na.rm = TRUE)
  idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(idx)) mat[idx] <- med[idx[, 1L]]
  mat
}

qc_dry_run_contract <- function(script, dataset, matrix_file = NULL, metadata_file = NULL, paths = NULL, extra = character()) {
  dry_run_line("Script", script)
  dry_run_line("Dataset", dataset)
  if (!is.null(matrix_file)) dry_run_line("Expression matrix", matrix_file, if (file.exists(matrix_file)) "PASS" else "FAIL")
  if (!is.null(metadata_file)) dry_run_line("Metadata", metadata_file, if (file.exists(metadata_file)) "PASS" else "WARN")
  if (!is.null(paths)) {
    invisible(lapply(unlist(paths), dir_create))
    dry_run_line("Output directories", paste(unlist(paths), collapse = "; "), "PASS")
  }
  if (length(extra)) for (line in extra) dry_run_line("Check", line)
  invisible(if (is.null(matrix_file) || file.exists(matrix_file)) 0L else 1L)
}
