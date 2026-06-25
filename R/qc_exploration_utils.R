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

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

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
  list(mat = mat, meta = meta)
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
