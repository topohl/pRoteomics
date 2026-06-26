#!/usr/bin/env Rscript

# Directionality audit: compare ProTigy forward logFC signs against raw/imputed abundance means.
#
# This is a diagnostic handoff check. It does not change downstream data.
# It answers whether a ProTigy forward contrast behaves like left-minus-right
# in the pre-ProTigy abundance matrix, or whether the exported logFC appears
# reversed relative to the biological comparison label.

suppressPackageStartupMessages({
  required_pkgs <- c("readxl", "readr", "dplyr", "tidyr", "stringr", "tibble", "purrr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
  }
})

source(if (file.exists(file.path("R", "script_runtime.R"))) file.path("R", "script_runtime.R") else file.path("..", "R", "script_runtime.R"))
source(repo_path("R", "dataset_inputs.R"))

runtime <- init_script_runtime(
  script = "01_preprocessing/03b_directionality_audit_raw_vs_protigy.r",
  stage = "core",
  default_dataset = "neuron_neuropil"
)
dataset <- runtime$dataset
Sys.setenv(PROTEOMICS_DATASET = dataset)
Sys.setenv(PROTEOMICS_SCRIPT_ID = runtime$script)

module_id <- "01_preprocessing/03b_directionality_audit_raw_vs_protigy"
out_dir <- path_results("tables", module_id, dataset)
log_dir <- path_results("logs", module_id, dataset)
dir_create(out_dir)
dir_create(log_dir)

message("Directionality audit dataset: ", dataset)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

norm_token <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:].-]+", "_", x)
  x <- gsub("_+", "_", x)
  x[x %in% c("", "na", "nan", "none", "null")] <- NA_character_
  x
}

normalize_group <- function(x) {
  x0 <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x0 %in% c("1", "CON", "CTRL", "CONTROL") ~ "con",
    x0 %in% c("2", "RES", "RESILIENT") ~ "res",
    x0 %in% c("3", "SUS", "SUSCEPTIBLE") ~ "sus",
    TRUE ~ norm_token(x0)
  )
}

first_existing_col <- function(df, candidates) {
  nms_clean <- tolower(gsub("[^a-z0-9]", "", names(df)))
  cand_clean <- tolower(gsub("[^a-z0-9]", "", candidates))
  hit <- match(cand_clean, nms_clean)
  hit <- hit[!is.na(hit)]
  if (!length(hit)) return(NA_character_)
  names(df)[hit[[1]]]
}

file_hash_or_na <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  file_hash_sha256(path) %||% file_hash(path)
}

parse_forward_fallback <- function(label) {
  label <- norm_token(label)
  m <- stringr::str_match(label, "^([a-z0-9]+?)(con|res|sus)_([a-z0-9]+?)(con|res|sus)$")
  if (is.na(m[1, 1])) return(NULL)
  split_prefix <- function(prefix) {
    region <- stringr::str_match(prefix, "^(ca[123]|dg)(.*)$")
    if (!is.na(region[1, 1])) {
      list(region = region[1, 2], token = ifelse(nzchar(region[1, 3]), region[1, 3], NA_character_))
    } else {
      list(region = NA_character_, token = prefix)
    }
  }
  left <- split_prefix(m[1, 2])
  right <- split_prefix(m[1, 4])
  list(
    parsed = TRUE,
    parse_mode = "parsed_forward_fallback",
    left_region = left$region,
    left_token = left$token,
    left_group = m[1, 3],
    right_region = right$region,
    right_token = right$token,
    right_group = m[1, 5]
  )
}

parse_comparison <- function(comparison, parsed_forward_comparison = NA_character_) {
  label_map <- c("1" = "con", "2" = "res", "3" = "sus", "con" = "con", "res" = "res", "sus" = "sus")
  comp <- as.character(comparison %||% NA_character_)
  m <- stringr::str_match(
    comp,
    "^([A-Za-z0-9]+)_([A-Za-z0-9]+)_([123]|con|res|sus)\\.over\\.([A-Za-z0-9]+)_([A-Za-z0-9]+)_([123]|con|res|sus)$"
  )
  if (!is.na(m[1, 1])) {
    return(list(
      parsed = TRUE,
      parse_mode = "index_comparison",
      left_region = norm_token(m[1, 2]),
      left_token = norm_token(m[1, 3]),
      left_group = unname(label_map[[tolower(m[1, 4])]]),
      right_region = norm_token(m[1, 5]),
      right_token = norm_token(m[1, 6]),
      right_group = unname(label_map[[tolower(m[1, 7])]])
    ))
  }
  fallback <- parse_forward_fallback(parsed_forward_comparison)
  if (!is.null(fallback)) return(fallback)
  list(
    parsed = FALSE,
    parse_mode = "failed",
    left_region = NA_character_,
    left_token = NA_character_,
    left_group = NA_character_,
    right_region = NA_character_,
    right_token = NA_character_,
    right_group = NA_character_
  )
}

prepare_metadata <- function(metadata, sample_candidates) {
  sample_col <- first_existing_col(metadata, sample_candidates)
  group_col <- first_existing_col(metadata, c("ExpGroup", "StressGroup", "group", "Group", "condition", "Condition"))
  region_col <- first_existing_col(metadata, c("region", "Region"))
  layer_col <- first_existing_col(metadata, c("layer", "Layer"))
  celltype_col <- first_existing_col(metadata, c("celltype", "CellType"))
  celltype_layer_col <- first_existing_col(metadata, c("celltype_layer", "CellTypeLayer", "cell_type_layer", "Cell.Type.Layer"))

  if (is.na(sample_col)) stop("Could not find sample identifier column in metadata.", call. = FALSE)
  if (is.na(group_col)) stop("Could not find group/ExpGroup column in metadata.", call. = FALSE)

  out <- tibble::tibble(
    sample_id = as.character(metadata[[sample_col]]),
    group = normalize_group(metadata[[group_col]]),
    region = if (!is.na(region_col)) norm_token(metadata[[region_col]]) else NA_character_,
    layer = if (!is.na(layer_col)) norm_token(metadata[[layer_col]]) else NA_character_,
    celltype = if (!is.na(celltype_col)) norm_token(metadata[[celltype_col]]) else NA_character_,
    celltype_layer = if (!is.na(celltype_layer_col)) norm_token(metadata[[celltype_layer_col]]) else NA_character_
  ) |>
    dplyr::filter(!is.na(.data$sample_id), nzchar(.data$sample_id)) |>
    dplyr::distinct(.data$sample_id, .keep_all = TRUE)

  if ("exclude" %in% names(metadata)) {
    keep <- is.na(metadata$exclude) | metadata$exclude != TRUE
    out <- out[keep[match(out$sample_id, as.character(metadata[[sample_col]]))] %||% TRUE, , drop = FALSE]
  }
  out
}

match_side_samples <- function(meta, side_region, side_token, side_group) {
  matched_fields <- c("group")
  warnings <- character()
  keep <- !is.na(meta$sample_id) & meta$group == side_group

  if (!is.na(side_region) && nzchar(side_region)) {
    if (any(meta$region == side_region, na.rm = TRUE)) {
      keep <- keep & meta$region == side_region
      matched_fields <- c(matched_fields, "region")
    } else {
      warnings <- c(warnings, paste0("region_not_found_in_metadata:", side_region))
    }
  }

  if (!is.na(side_token) && nzchar(side_token)) {
    token_matched <- FALSE
    if (any(meta$layer == side_token, na.rm = TRUE)) {
      keep <- keep & meta$layer == side_token
      matched_fields <- c(matched_fields, "layer")
      token_matched <- TRUE
    } else if (any(meta$celltype_layer == side_token, na.rm = TRUE)) {
      keep <- keep & meta$celltype_layer == side_token
      matched_fields <- c(matched_fields, "celltype_layer")
      token_matched <- TRUE
    } else if (any(meta$celltype == side_token, na.rm = TRUE)) {
      keep <- keep & meta$celltype == side_token
      matched_fields <- c(matched_fields, "celltype")
      token_matched <- TRUE
    } else if (side_token %in% c("neuron", "neuropil", "neuron_neuropil", "soma", "neuron_soma", "microglia", "microglial")) {
      token_norm <- dplyr::case_when(
        side_token == "neuropil" ~ "neuron_neuropil",
        side_token == "soma" ~ "neuron_soma",
        side_token == "microglial" ~ "microglia",
        TRUE ~ side_token
      )
      if (any(meta$celltype == token_norm, na.rm = TRUE)) {
        keep <- keep & meta$celltype == token_norm
        matched_fields <- c(matched_fields, "celltype")
        token_matched <- TRUE
      }
    }
    if (!token_matched) {
      warnings <- c(warnings, paste0("stratum_token_not_matched_to_metadata:", side_token))
    }
  }

  list(
    sample_ids = meta$sample_id[keep],
    n_samples = sum(keep, na.rm = TRUE),
    matched_fields = paste(unique(matched_fields), collapse = ";"),
    warning = paste(unique(warnings), collapse = ";")
  )
}

summarise_audit <- function(audit, comparison_label, parsed, left, right) {
  sig_col <- if ("padj" %in% names(audit) && any(is.finite(audit$padj))) "padj" else if ("pval" %in% names(audit) && any(is.finite(audit$pval))) "pval" else NA_character_
  audit <- audit |>
    dplyr::mutate(
      significant = if (!is.na(sig_col)) .data[[sig_col]] <= 0.05 else FALSE,
      valid_sign = is.finite(.data$log2fc) & is.finite(.data$raw_delta_left_minus_right) & .data$log2fc != 0 & .data$raw_delta_left_minus_right != 0,
      sign_agrees = dplyr::if_else(.data$valid_sign, sign(.data$log2fc) == sign(.data$raw_delta_left_minus_right), NA)
    )

  n_valid <- sum(audit$valid_sign, na.rm = TRUE)
  frac_agree <- if (n_valid > 0L) mean(audit$sign_agrees[audit$valid_sign], na.rm = TRUE) else NA_real_
  cor_pearson <- if (n_valid >= 3L) suppressWarnings(stats::cor(audit$log2fc, audit$raw_delta_left_minus_right, use = "pairwise.complete.obs", method = "pearson")) else NA_real_
  cor_spearman <- if (n_valid >= 3L) suppressWarnings(stats::cor(audit$log2fc, audit$raw_delta_left_minus_right, use = "pairwise.complete.obs", method = "spearman")) else NA_real_

  orientation <- dplyr::case_when(
    n_valid < 10L ~ "insufficient_overlap",
    is.finite(cor_pearson) & cor_pearson >= 0.80 & is.finite(frac_agree) & frac_agree >= 0.80 ~ "forward_matches_raw_left_minus_right",
    is.finite(cor_pearson) & cor_pearson <= -0.80 & is.finite(frac_agree) & frac_agree <= 0.20 ~ "forward_appears_reversed_vs_raw_left_minus_right",
    TRUE ~ "mixed_or_model_dependent"
  )

  tibble::tibble(
    dataset = dataset,
    comparison = comparison_label,
    parsed = parsed$parsed,
    parse_mode = parsed$parse_mode,
    left_label = paste(stats::na.omit(c(parsed$left_region, parsed$left_token, parsed$left_group)), collapse = "_"),
    right_label = paste(stats::na.omit(c(parsed$right_region, parsed$right_token, parsed$right_group)), collapse = "_"),
    left_group = parsed$left_group,
    right_group = parsed$right_group,
    left_region = parsed$left_region,
    right_region = parsed$right_region,
    left_token = parsed$left_token,
    right_token = parsed$right_token,
    left_n_samples = left$n_samples,
    right_n_samples = right$n_samples,
    left_matched_fields = left$matched_fields,
    right_matched_fields = right$matched_fields,
    matching_warning = paste(unique(c(left$warning, right$warning)[nzchar(c(left$warning, right$warning))]), collapse = ";"),
    n_overlap_proteins = nrow(audit),
    n_valid_signs = n_valid,
    frac_sign_agree_forward_vs_raw_delta = frac_agree,
    cor_log2fc_raw_delta_pearson = cor_pearson,
    cor_log2fc_raw_delta_spearman = cor_spearman,
    inferred_forward_orientation = orientation,
    n_significant = sum(audit$significant, na.rm = TRUE),
    significant_metric = sig_col %||% NA_character_,
    n_significant_protigy_positive = sum(audit$significant & audit$log2fc > 0, na.rm = TRUE),
    n_significant_protigy_negative = sum(audit$significant & audit$log2fc < 0, na.rm = TRUE),
    n_significant_raw_left_higher = sum(audit$significant & audit$raw_delta_left_minus_right > 0, na.rm = TRUE),
    n_significant_raw_right_higher = sum(audit$significant & audit$raw_delta_left_minus_right < 0, na.rm = TRUE),
    note = "Positive raw_delta_left_minus_right means the left comparison group has higher raw/imputed abundance than the right comparison group."
  )
}

dataset_inputs <- resolve_dataset_inputs(dataset, purpose = "wgcna", script = runtime$script, stage = runtime$stage)
expr_xlsx <- Sys.getenv("PROTEOMICS_DIRECTIONALITY_RAW_XLSX", unset = dataset_inputs$expression_file)
meta_xlsx <- Sys.getenv("PROTEOMICS_DIRECTIONALITY_METADATA_XLSX", unset = dataset_inputs$metadata_file)
gct_extract_dir <- Sys.getenv(
  "PROTEOMICS_DIRECTIONALITY_GCT_EXTRACT_DIR",
  unset = path_processed("01_preprocessing", "gct_extractR", dataset)
)
index_path <- file.path(gct_extract_dir, "indexComparisons.csv")
forward_dir <- file.path(gct_extract_dir, "forward")

if (isTRUE(runtime$dry_run)) {
  message("[dry-run] expression workbook: ", expr_xlsx, " exists=", file.exists(expr_xlsx))
  message("[dry-run] metadata workbook: ", meta_xlsx, " exists=", file.exists(meta_xlsx))
  message("[dry-run] GCT extract index: ", index_path, " exists=", file.exists(index_path))
  message("[dry-run] forward contrast directory: ", forward_dir, " exists=", dir.exists(forward_dir))
  message("[dry-run] outputs: ", out_dir)
  quit(status = if (file.exists(expr_xlsx) && file.exists(meta_xlsx) && file.exists(index_path) && dir.exists(forward_dir)) 0 else 1, save = "no")
}

if (!file.exists(expr_xlsx)) stop("Expression workbook not found: ", expr_xlsx, call. = FALSE)
if (!file.exists(meta_xlsx)) stop("Metadata workbook not found: ", meta_xlsx, call. = FALSE)
if (!file.exists(index_path)) stop("GCT extraction index not found: ", index_path, call. = FALSE)
if (!dir.exists(forward_dir)) stop("Forward GCT extraction directory not found: ", forward_dir, call. = FALSE)

expr <- as.data.frame(readxl::read_excel(expr_xlsx))
metadata <- as.data.frame(readxl::read_excel(meta_xlsx))
index <- readr::read_csv(index_path, show_col_types = FALSE)

protein_col <- first_existing_col(expr, dataset_inputs$protein_id_col_candidates)
if (is.na(protein_col)) {
  stop("Could not find a protein identifier column in expression workbook: ", expr_xlsx, call. = FALSE)
}

non_sample_cols <- unique(c(dataset_inputs$protein_id_col_candidates, "Protein.Group", "First.Protein.Description", protein_col))
candidate_sample_cols <- setdiff(names(expr), non_sample_cols)
numeric_fraction <- vapply(candidate_sample_cols, function(cc) {
  vals <- suppressWarnings(as.numeric(as.character(expr[[cc]])))
  mean(!is.na(vals))
}, numeric(1))
sample_cols <- candidate_sample_cols[numeric_fraction >= 0.70]
if (length(sample_cols) < 4L) {
  stop("Could not detect expression sample columns in: ", expr_xlsx, call. = FALSE)
}

expr_num <- expr[, sample_cols, drop = FALSE]
expr_num[] <- lapply(expr_num, function(x) suppressWarnings(as.numeric(as.character(x))))
expr_num <- as.data.frame(expr_num, check.names = FALSE)
expr_num$gene_symbol <- as.character(expr[[protein_col]])
expr_num <- expr_num[!is.na(expr_num$gene_symbol) & nzchar(expr_num$gene_symbol), , drop = FALSE]

meta_prepared <- prepare_metadata(metadata, dataset_inputs$sample_id_col_candidates)
matched_samples <- intersect(sample_cols, meta_prepared$sample_id)
if (!length(matched_samples)) {
  stop("No expression sample columns matched metadata sample IDs.", call. = FALSE)
}
meta_prepared <- meta_prepared |>
  dplyr::filter(.data$sample_id %in% sample_cols)

if (!all(c("comparison", "parsed_forward_comparison", "forward_file") %in% names(index))) {
  stop("indexComparisons.csv must contain comparison, parsed_forward_comparison, and forward_file columns.", call. = FALSE)
}

comparison_arg <- script_arg_value("--comparison", default = "", args = runtime$args)
if (nzchar(comparison_arg)) {
  index <- index |>
    dplyr::filter(
      stringr::str_detect(.data$comparison, stringr::fixed(comparison_arg)) |
        stringr::str_detect(.data$parsed_forward_comparison, stringr::fixed(comparison_arg)) |
        stringr::str_detect(basename(.data$forward_file), stringr::fixed(comparison_arg))
    )
}
if (!nrow(index)) stop("No GCT forward comparisons matched the requested audit selection.", call. = FALSE)

all_audits <- list()
all_summaries <- list()

for (i in seq_len(nrow(index))) {
  row <- index[i, , drop = FALSE]
  forward_file <- as.character(row$forward_file)
  if (!file.exists(forward_file)) {
    forward_file <- file.path(forward_dir, basename(forward_file))
  }
  if (!file.exists(forward_file)) {
    warning("Skipping missing forward contrast file: ", as.character(row$forward_file), call. = FALSE)
    next
  }

  parsed <- parse_comparison(row$comparison, row$parsed_forward_comparison)
  if (!isTRUE(parsed$parsed)) {
    warning("Could not parse comparison strata for: ", row$comparison, call. = FALSE)
    next
  }

  left <- match_side_samples(meta_prepared, parsed$left_region, parsed$left_token, parsed$left_group)
  right <- match_side_samples(meta_prepared, parsed$right_region, parsed$right_token, parsed$right_group)
  left_samples <- intersect(left$sample_ids, sample_cols)
  right_samples <- intersect(right$sample_ids, sample_cols)

  if (length(left_samples) < 1L || length(right_samples) < 1L) {
    warning("Skipping comparison with missing left/right raw samples: ", row$parsed_forward_comparison, call. = FALSE)
    next
  }

  left_mean <- rowMeans(as.matrix(expr_num[, left_samples, drop = FALSE]), na.rm = TRUE)
  right_mean <- rowMeans(as.matrix(expr_num[, right_samples, drop = FALSE]), na.rm = TRUE)
  left_nonmiss <- rowSums(!is.na(as.matrix(expr_num[, left_samples, drop = FALSE])))
  right_nonmiss <- rowSums(!is.na(as.matrix(expr_num[, right_samples, drop = FALSE])))
  raw_delta <- tibble::tibble(
    gene_symbol = expr_num$gene_symbol,
    raw_left_mean = left_mean,
    raw_right_mean = right_mean,
    raw_delta_left_minus_right = left_mean - right_mean,
    raw_left_n_nonmissing = left_nonmiss,
    raw_right_n_nonmissing = right_nonmiss
  ) |>
    dplyr::filter(is.finite(.data$raw_delta_left_minus_right)) |>
    dplyr::group_by(.data$gene_symbol) |>
    dplyr::summarise(
      raw_left_mean = mean(.data$raw_left_mean, na.rm = TRUE),
      raw_right_mean = mean(.data$raw_right_mean, na.rm = TRUE),
      raw_delta_left_minus_right = mean(.data$raw_delta_left_minus_right, na.rm = TRUE),
      raw_left_n_nonmissing = max(.data$raw_left_n_nonmissing, na.rm = TRUE),
      raw_right_n_nonmissing = max(.data$raw_right_n_nonmissing, na.rm = TRUE),
      .groups = "drop"
    )

  de <- readr::read_csv(forward_file, show_col_types = FALSE)
  log_col <- first_existing_col(de, c("log2fc", "logFC", "RawlogFC", "rawlog2fc"))
  if (is.na(log_col)) {
    warning("Skipping contrast without logFC/log2fc column: ", forward_file, call. = FALSE)
    next
  }
  de_gene_col <- first_existing_col(de, c("gene_symbol", "T: Protein.Names", "Genes", "Protein.Group", "ProteinID", "UniProt"))
  if (is.na(de_gene_col)) {
    warning("Skipping contrast without protein ID column: ", forward_file, call. = FALSE)
    next
  }
  pval_col <- first_existing_col(de, c("pval", "P.Value", "p.value", "PValue"))
  padj_col <- first_existing_col(de, c("padj", "adj.P.Val", "FDR", "qvalue", "q.value"))

  de_small <- tibble::tibble(
    gene_symbol = as.character(de[[de_gene_col]]),
    log2fc = suppressWarnings(as.numeric(as.character(de[[log_col]]))),
    pval = if (!is.na(pval_col)) suppressWarnings(as.numeric(as.character(de[[pval_col]]))) else NA_real_,
    padj = if (!is.na(padj_col)) suppressWarnings(as.numeric(as.character(de[[padj_col]]))) else NA_real_
  ) |>
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol))

  audit <- de_small |>
    dplyr::inner_join(raw_delta, by = "gene_symbol") |>
    dplyr::mutate(
      dataset = dataset,
      comparison = as.character(row$comparison),
      parsed_forward_comparison = as.character(row$parsed_forward_comparison),
      forward_file = normalizePath(forward_file, winslash = "/", mustWork = FALSE),
      left_label = paste(stats::na.omit(c(parsed$left_region, parsed$left_token, parsed$left_group)), collapse = "_"),
      right_label = paste(stats::na.omit(c(parsed$right_region, parsed$right_token, parsed$right_group)), collapse = "_"),
      protigy_sign = sign(.data$log2fc),
      raw_sign_left_minus_right = sign(.data$raw_delta_left_minus_right),
      valid_sign = is.finite(.data$log2fc) & is.finite(.data$raw_delta_left_minus_right) & .data$log2fc != 0 & .data$raw_delta_left_minus_right != 0,
      sign_agrees_forward_vs_raw_delta = dplyr::if_else(.data$valid_sign, .data$protigy_sign == .data$raw_sign_left_minus_right, NA)
    ) |>
    dplyr::select(
      .data$dataset, .data$comparison, .data$parsed_forward_comparison, .data$forward_file,
      .data$left_label, .data$right_label, .data$gene_symbol, .data$log2fc, .data$pval, .data$padj,
      .data$raw_left_mean, .data$raw_right_mean, .data$raw_delta_left_minus_right,
      .data$raw_left_n_nonmissing, .data$raw_right_n_nonmissing,
      .data$protigy_sign, .data$raw_sign_left_minus_right,
      .data$valid_sign, .data$sign_agrees_forward_vs_raw_delta
    )

  if (!nrow(audit)) {
    warning("No overlapping protein IDs between ProTigy and raw matrix for: ", row$parsed_forward_comparison, call. = FALSE)
    next
  }

  all_audits[[length(all_audits) + 1L]] <- audit
  all_summaries[[length(all_summaries) + 1L]] <- summarise_audit(audit, as.character(row$comparison), parsed, left, right)
}

if (!length(all_audits)) {
  stop("Directionality audit produced no comparable contrast rows. Check protein IDs, comparison parsing, and metadata sample IDs.", call. = FALSE)
}

audit_long <- dplyr::bind_rows(all_audits)
summary_tbl <- dplyr::bind_rows(all_summaries) |>
  dplyr::arrange(.data$inferred_forward_orientation, .data$comparison)

decision_tbl <- summary_tbl |>
  dplyr::count(.data$inferred_forward_orientation, name = "n_comparisons") |>
  dplyr::mutate(
    dataset = dataset,
    interpretation = dplyr::case_when(
      .data$inferred_forward_orientation == "forward_matches_raw_left_minus_right" ~ "ProTigy forward logFC agrees with raw left-minus-right abundance means for these comparisons.",
      .data$inferred_forward_orientation == "forward_appears_reversed_vs_raw_left_minus_right" ~ "ProTigy forward logFC appears reversed relative to raw left-minus-right abundance means for these comparisons.",
      .data$inferred_forward_orientation == "mixed_or_model_dependent" ~ "Direction cannot be resolved globally; inspect matching warnings, spatial strata, model design, and protein-level rows.",
      .data$inferred_forward_orientation == "insufficient_overlap" ~ "Insufficient protein overlap to infer direction reliably.",
      TRUE ~ "Unclassified orientation."
    )
  ) |>
  dplyr::select(.data$dataset, .data$inferred_forward_orientation, .data$n_comparisons, .data$interpretation)

audit_path <- file.path(out_dir, "directionality_audit_raw_vs_protigy_by_protein.csv")
summary_path <- file.path(out_dir, "directionality_audit_summary.csv")
decision_path <- file.path(out_dir, "directionality_audit_decision.csv")
input_status_path <- file.path(log_dir, "input_status.csv")
manifest_path <- file.path(log_dir, "run_manifest.yml")

readr::write_csv(audit_long, audit_path, na = "")
readr::write_csv(summary_tbl, summary_path, na = "")
readr::write_csv(decision_tbl, decision_path, na = "")

input_rows <- dplyr::bind_rows(
  input_status_row("pre_protigy_expression_matrix", expr_xlsx, dataset, required = TRUE, n_rows = nrow(expr)),
  input_status_row("sample_metadata", meta_xlsx, dataset, required = TRUE, n_rows = nrow(metadata)),
  input_status_row("gct_extract_index", index_path, dataset, required = TRUE, n_rows = nrow(index)),
  input_status_row("gct_extract_forward_dir", forward_dir, dataset, required = TRUE, n_rows = length(list.files(forward_dir, pattern = "\\.csv$", full.names = TRUE)))
)
write_input_status(input_rows, input_status_path, dry_run = FALSE)

finish_script_runtime(
  runtime,
  manifest_path = manifest_path,
  inputs = c(
    pre_protigy_expression_matrix = expr_xlsx,
    sample_metadata = meta_xlsx,
    gct_extract_index = index_path,
    gct_extract_forward_dir = forward_dir
  ),
  outputs = c(
    protein_level_audit = audit_path,
    comparison_summary = summary_path,
    decision_summary = decision_path,
    input_status = input_status_path
  ),
  notes = c(
    "Diagnostic-only directionality audit; no downstream data are modified.",
    "Positive raw_delta_left_minus_right means the left comparison group has higher raw/imputed abundance than the right comparison group.",
    "Use this to verify whether ProTigy forward logFC signs agree with raw/imputed abundance means before interpreting GSEA direction."
  )
)

message("Wrote protein-level audit: ", audit_path)
message("Wrote comparison summary: ", summary_path)
message("Wrote decision summary: ", decision_path)
message("Finished directionality audit for dataset: ", dataset)
