# ================================================================
# Script: 04_differential_expression_enrichment/01b_gsea_protein_direction_audit.r
# Stage: enrichment
# Scope: dataset_specific
# Consumes: clusterProfiler manifest/output tables and mapped forward protein files.
# Produces: diagnostic GSEA/protein direction audit CSVs under results/tables and run_manifest.yml under results/logs.
# Dataset behavior: runs for neuron_neuropil, neuron_soma, microglia via --dataset.
# Notes: Read-only diagnostic. Does not modify ProTigy, mapped, clusterProfiler, or compareGO outputs.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "script_runtime.R"))

MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- "01b_gsea_protein_direction_audit"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

runtime <- init_script_runtime(
  script = "04_differential_expression_enrichment/01b_gsea_protein_direction_audit.r",
  stage = "enrichment",
  default_dataset = "neuron_neuropil"
)
args <- runtime$args
DATASET <- runtime$dataset
ONTOLOGY <- toupper(script_arg_value("--ontology", "BP", args = args))
COMPARISON_FILTER <- script_arg_value("--comparison", "", args = args)

if (!ONTOLOGY %in% c("BP", "MF", "CC")) {
  stop("--ontology must be one of BP, MF, or CC.", call. = FALSE)
}

tables_dir <- path_results("tables", MODULE_ID, SUBSTEP_ID, DATASET, ONTOLOGY)
logs_dir <- path_results("logs", MODULE_ID, SUBSTEP_ID, DATASET, ONTOLOGY)
output_paths <- list(
  term_audit = file.path(tables_dir, "gsea_term_direction_audit.csv"),
  contrast_summary = file.path(tables_dir, "gsea_contrast_direction_summary.csv"),
  ora_warning = file.path(tables_dir, "ora_pooled_direction_warning.csv"),
  input_status = file.path(tables_dir, "input_status.csv"),
  run_manifest = file.path(logs_dir, "run_manifest.yml")
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

read_csv_base <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

as_logical_manifest <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

canonical_manifest_path <- function(dataset) {
  path_processed(MODULE_ID, "clusterProfiler", dataset, "clusterProfiler_manifest.csv")
}

discover_manifest <- function(dataset) {
  canonical <- canonical_manifest_path(dataset)
  if (file.exists(canonical)) return(normalizePath(canonical, winslash = "/", mustWork = FALSE))

  candidates <- c(
    list.files(
      path_processed(MODULE_ID, "clusterProfiler", dataset),
      pattern = "^clusterProfiler_manifest.*\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    ),
    list.files(
      path_results("tables", MODULE_ID, "01_clusterProfiler", dataset),
      pattern = "^clusterProfiler_manifest.*\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    ),
    list.files(
      path_results("tables", MODULE_ID, "clusterProfiler", dataset),
      pattern = "^clusterProfiler_manifest.*\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(normalizePath(canonical, winslash = "/", mustWork = FALSE))
  info <- file.info(candidates)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

fill_missing_manifest_dataset <- function(manifest, dataset, manifest_path) {
  if (!"dataset" %in% names(manifest)) return(manifest)
  dataset_missing <- is.na(manifest$dataset) |
    !nzchar(trimws(as.character(manifest$dataset))) |
    toupper(trimws(as.character(manifest$dataset))) == "NA"
  manifest_parent_dataset <- basename(dirname(normalizePath(manifest_path, winslash = "/", mustWork = FALSE)))
  if (any(dataset_missing) && identical(manifest_parent_dataset, dataset)) {
    manifest$dataset[dataset_missing] <- dataset
  }
  manifest
}

normalize_col_key <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(trimws(as.character(x))))
}

rename_first_alias <- function(df, target, aliases) {
  keys <- normalize_col_key(names(df))
  alias_keys <- normalize_col_key(aliases)
  hit <- which(keys %in% alias_keys)
  if (length(hit)) names(df)[hit[[1]]] <- target
  df
}

normalize_protein_table <- function(df, path) {
  df <- rename_first_alias(df, "gene_symbol", c("gene_symbol", "gene", "genes", "symbol", "UNIPROT", "id"))
  df <- rename_first_alias(df, "log2fc", c("log2fc", "logFC", "log2foldchange", "avg_log2FC"))
  df <- rename_first_alias(df, "pvalue", c("pvalue", "pval", "p.value"))
  df <- rename_first_alias(df, "padj", c("padj", "adj.P.Val", "p.adjust", "FDR", "qvalue"))

  required <- c("gene_symbol", "log2fc")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("Mapped protein table is missing required normalized column(s) ",
         paste(missing, collapse = ", "), ": ", path, call. = FALSE)
  }

  df$gene_symbol <- trimws(as.character(df$gene_symbol))
  df$log2fc <- suppressWarnings(as.numeric(df$log2fc))
  df <- df[!is.na(df$gene_symbol) & nzchar(df$gene_symbol), , drop = FALSE]
  df <- df[!is.na(df$log2fc), , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol), , drop = FALSE]
  df
}

validate_gsea_table <- function(df, path) {
  required <- c("ID", "Description", "NES", "p.adjust", "core_enrichment")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("GSEA table is missing required column(s) ",
         paste(missing, collapse = ", "), ": ", path, call. = FALSE)
  }
  df$NES <- suppressWarnings(as.numeric(df$NES))
  df$`p.adjust` <- suppressWarnings(as.numeric(df$`p.adjust`))
  df$core_enrichment <- as.character(df$core_enrichment)
  df
}

sign_label <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x > 0] <- "positive"
  out[!is.na(x) & x < 0] <- "negative"
  out[!is.na(x) & x == 0] <- "zero"
  out
}

parse_comparison_sides <- function(comparison) {
  comp <- tolower(as.character(comparison))
  compact <- gsub("[^a-z0-9.]+", "_", comp)
  if (grepl("2\\.over\\.1|res[_-]?con|con[_-]?res", compact)) {
    return(list(formal_contrast = comparison, positive_side_label = "RES", negative_side_label = "CON"))
  }
  if (grepl("3\\.over\\.2|sus[_-]?res|res[_-]?sus", compact)) {
    return(list(formal_contrast = comparison, positive_side_label = "SUS", negative_side_label = "RES"))
  }
  if (grepl("3\\.over\\.1|sus[_-]?con|con[_-]?sus", compact)) {
    return(list(formal_contrast = comparison, positive_side_label = "SUS", negative_side_label = "CON"))
  }
  list(formal_contrast = comparison, positive_side_label = "positive_log2fc_side", negative_side_label = "negative_log2fc_side")
}

side_from_sign <- function(sign_value, positive_side, negative_side) {
  if (is.na(sign_value) || identical(sign_value, "zero")) return(NA_character_)
  if (identical(sign_value, "positive")) positive_side else negative_side
}

load_threshold <- function() {
  cfg_paths <- c(
    repo_path("config", "clusterProfiler_config.local.yml"),
    repo_path("config", "clusterProfiler_config.yml")
  )
  cfg_path <- cfg_paths[file.exists(cfg_paths)][1]
  if (length(cfg_path) && !is.na(cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
    value <- cfg$analysis$top_gene_abs_log2fc %||% NA_real_
    value <- suppressWarnings(as.numeric(value))
    if (!is.na(value)) return(value)
  }
  1
}

input_status_rows <- function(manifest_path, manifest_filtered, output_paths) {
  rows <- list(
    input_status_row("clusterProfiler_manifest", manifest_path, dataset = DATASET, required = TRUE,
                     n_rows = if (file.exists(manifest_path)) nrow(read_csv_base(manifest_path)) else NA_integer_)
  )
  if (!is.null(manifest_filtered) && nrow(manifest_filtered)) {
    rows <- c(rows, lapply(seq_len(nrow(manifest_filtered)), function(i) {
      input_status_row(
        paste0("gsea_output_table:", manifest_filtered$comparison[[i]]),
        manifest_filtered$output_table[[i]],
        dataset = DATASET,
        required = TRUE
      )
    }))
    rows <- c(rows, lapply(seq_len(nrow(manifest_filtered)), function(i) {
      input_status_row(
        paste0("mapped_input_file:", manifest_filtered$comparison[[i]]),
        manifest_filtered$input_gene_file[[i]],
        dataset = DATASET,
        required = TRUE
      )
    }))
  }
  do.call(rbind, rows)
}

manifest_path <- discover_manifest(DATASET)
manifest_exists <- file.exists(manifest_path)
manifest <- if (manifest_exists) read_csv_base(manifest_path) else data.frame()
if (manifest_exists) manifest <- fill_missing_manifest_dataset(manifest, DATASET, manifest_path)

required_manifest_cols <- c("dataset", "ontology", "result_type", "used_for_plot", "comparison", "input_gene_file", "output_table")
missing_manifest_cols <- if (manifest_exists) setdiff(required_manifest_cols, names(manifest)) else required_manifest_cols

manifest_filtered <- data.frame()
if (manifest_exists && !length(missing_manifest_cols)) {
  manifest$used_for_plot <- as_logical_manifest(manifest$used_for_plot)
  manifest_filtered <- manifest[
    manifest$dataset == DATASET &
      toupper(as.character(manifest$ontology)) == ONTOLOGY &
      manifest$result_type == "GSEA_GO" &
      manifest$used_for_plot %in% TRUE,
    ,
    drop = FALSE
  ]
  if (nzchar(COMPARISON_FILTER)) {
    manifest_filtered <- manifest_filtered[
      grepl(COMPARISON_FILTER, manifest_filtered$comparison, fixed = TRUE),
      ,
      drop = FALSE
    ]
  }
}

status_df <- input_status_rows(manifest_path, manifest_filtered, output_paths)

if (isTRUE(runtime$dry_run)) {
  dry_run_line("Dataset", DATASET)
  dry_run_line("Ontology", ONTOLOGY)
  dry_run_line("Comparison filter", if (nzchar(COMPARISON_FILTER)) COMPARISON_FILTER else "<none>")
  dry_run_line("Manifest", manifest_path, if (manifest_exists) "PASS" else "FAIL")
  dry_run_line("Manifest rows", if (manifest_exists) nrow(manifest) else 0L)
  dry_run_line("Filtered GSEA_GO plot rows", nrow(manifest_filtered))
  if (length(missing_manifest_cols)) dry_run_line("Missing manifest columns", paste(missing_manifest_cols, collapse = ", "), "FAIL")
  dry_run_line("Output term audit", output_paths$term_audit)
  dry_run_line("Output contrast summary", output_paths$contrast_summary)
  dry_run_line("Output ORA warning", output_paths$ora_warning)
  dry_run_line("Output input status", output_paths$input_status)
  dry_run_line("Output run manifest", output_paths$run_manifest)
  if (nrow(status_df)) {
    for (i in seq_len(nrow(status_df))) {
      dry_run_line(status_df$input_name[[i]], status_df$path[[i]], if (status_df$status[[i]] == "present") "PASS" else "FAIL")
    }
  }
  quit(save = "no", status = 0L)
}

missing_pkgs <- c("dplyr", "tidyr")[!vapply(c("dplyr", "tidyr"), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       ". Install them explicitly before running this script.", call. = FALSE)
}

if (!manifest_exists) {
  stop("clusterProfiler manifest not found: ", manifest_path,
       "\nRun 04_differential_expression_enrichment/01_clusterProfiler.r first.", call. = FALSE)
}
if (length(missing_manifest_cols)) {
  stop("clusterProfiler manifest is missing required column(s): ",
       paste(missing_manifest_cols, collapse = ", "), call. = FALSE)
}
if (!nrow(manifest_filtered)) {
  stop("No manifest rows matched dataset=", DATASET, ", ontology=", ONTOLOGY,
       ", result_type=GSEA_GO, used_for_plot=TRUE",
       if (nzchar(COMPARISON_FILTER)) paste0(", comparison containing '", COMPARISON_FILTER, "'") else "",
       ".", call. = FALSE)
}

missing_tables <- manifest_filtered$output_table[!file.exists(manifest_filtered$output_table)]
missing_inputs <- manifest_filtered$input_gene_file[!file.exists(manifest_filtered$input_gene_file)]
if (length(missing_tables) || length(missing_inputs)) {
  write_input_status(status_df, output_paths$input_status, dry_run = FALSE)
  stop("Required audit input(s) are missing.\nMissing GSEA table(s): ",
       paste(missing_tables, collapse = "; "),
       "\nMissing mapped protein file(s): ", paste(missing_inputs, collapse = "; "),
       call. = FALSE)
}

top_abs_threshold <- load_threshold()
term_rows <- list()
ora_rows <- list()

for (i in seq_len(nrow(manifest_filtered))) {
  row <- manifest_filtered[i, , drop = FALSE]
  comparison <- as.character(row$comparison[[1]])
  sides <- parse_comparison_sides(comparison)

  gsea <- validate_gsea_table(read_csv_base(row$output_table[[1]]), row$output_table[[1]])
  protein <- normalize_protein_table(read_csv_base(row$input_gene_file[[1]]), row$input_gene_file[[1]])

  if (!nrow(gsea)) next

  core_long <- gsea |>
    dplyr::mutate(
      contrast = comparison,
      core_gene = strsplit(as.character(.data$core_enrichment), "/", fixed = TRUE)
    ) |>
    tidyr::unnest(core_gene) |>
    dplyr::mutate(core_gene = trimws(as.character(.data$core_gene))) |>
    dplyr::filter(!is.na(.data$core_gene), nzchar(.data$core_gene)) |>
    dplyr::left_join(
      protein[, c("gene_symbol", "log2fc"), drop = FALSE],
      by = c("core_gene" = "gene_symbol")
    )

  audit_one <- core_long |>
    dplyr::group_by(.data$contrast, .data$ID, .data$Description) |>
    dplyr::summarise(
      NES = dplyr::first(.data$NES),
      p.adjust = dplyr::first(.data$`p.adjust`),
      n_core = dplyr::n_distinct(.data$core_gene),
      n_core_mapped = sum(!is.na(.data$log2fc)),
      median_core_log2fc = ifelse(all(is.na(.data$log2fc)), NA_real_, stats::median(.data$log2fc, na.rm = TRUE)),
      mean_core_log2fc = ifelse(all(is.na(.data$log2fc)), NA_real_, mean(.data$log2fc, na.rm = TRUE)),
      frac_core_positive = ifelse(sum(!is.na(.data$log2fc)) > 0, mean(.data$log2fc > 0, na.rm = TRUE), NA_real_),
      frac_core_negative = ifelse(sum(!is.na(.data$log2fc)) > 0, mean(.data$log2fc < 0, na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dataset = DATASET,
      ontology = ONTOLOGY,
      formal_contrast = sides$formal_contrast,
      positive_side_label = sides$positive_side_label,
      negative_side_label = sides$negative_side_label,
      NES_sign = sign_label(.data$NES),
      median_core_sign = sign_label(.data$median_core_log2fc),
      direction_consistent = !is.na(.data$NES_sign) & !is.na(.data$median_core_sign) & .data$NES_sign == .data$median_core_sign,
      direction_strength = pmax(.data$frac_core_positive, .data$frac_core_negative, na.rm = TRUE),
      direction_strength = ifelse(is.infinite(.data$direction_strength), NA_real_, .data$direction_strength),
      classification = dplyr::case_when(
        .data$n_core_mapped < 3 | is.na(.data$direction_strength) | .data$direction_strength < 0.60 ~ "weak_or_mixed_core",
        .data$direction_consistent ~ "consistent",
        TRUE ~ "inconsistent_check_mapping_or_contrast"
      ),
      NES_biological_side = vapply(.data$NES_sign, side_from_sign, character(1),
                                   positive_side = sides$positive_side_label,
                                   negative_side = sides$negative_side_label),
      core_median_biological_side = vapply(.data$median_core_sign, side_from_sign, character(1),
                                           positive_side = sides$positive_side_label,
                                           negative_side = sides$negative_side_label),
      input_gene_file = row$input_gene_file[[1]],
      output_table = row$output_table[[1]]
    ) |>
    dplyr::select(dplyr::all_of(c(
      "dataset", "ontology", "contrast", "formal_contrast",
      "positive_side_label", "negative_side_label",
      "NES_biological_side", "core_median_biological_side",
      "ID", "Description", "NES", "p.adjust",
      "n_core", "n_core_mapped", "median_core_log2fc",
      "mean_core_log2fc", "frac_core_positive", "frac_core_negative",
      "NES_sign", "median_core_sign", "direction_consistent",
      "direction_strength", "classification",
      "input_gene_file", "output_table"
    )))

  term_rows[[length(term_rows) + 1L]] <- audit_one

  ora_rows[[length(ora_rows) + 1L]] <- data.frame(
    dataset = DATASET,
    ontology = ONTOLOGY,
    contrast = comparison,
    formal_contrast = sides$formal_contrast,
    positive_side_label = sides$positive_side_label,
    negative_side_label = sides$negative_side_label,
    n_positive_log2fc = sum(protein$log2fc > 0, na.rm = TRUE),
    n_negative_log2fc = sum(protein$log2fc < 0, na.rm = TRUE),
    top_abs_log2fc_threshold = top_abs_threshold,
    n_top_abs_log2fc = sum(abs(protein$log2fc) > top_abs_threshold, na.rm = TRUE),
    n_top_positive = sum(protein$log2fc > top_abs_threshold, na.rm = TRUE),
    n_top_negative = sum(protein$log2fc < -top_abs_threshold, na.rm = TRUE),
    note = "Current ORA/top-regulated GO is based on abs(log2fc) and is not direction-specific unless split into positive and negative log2fc sets.",
    input_gene_file = row$input_gene_file[[1]],
    stringsAsFactors = FALSE
  )
}

term_audit <- dplyr::bind_rows(term_rows)
ora_warning <- dplyr::bind_rows(ora_rows)

if (!nrow(term_audit)) {
  stop("No GSEA terms were available for audit after reading matched manifest rows.", call. = FALSE)
}

summary <- term_audit |>
  dplyr::group_by(.data$dataset, .data$ontology, .data$contrast, .data$formal_contrast,
                  .data$positive_side_label, .data$negative_side_label) |>
  dplyr::summarise(
    n_terms = dplyr::n(),
    n_consistent = sum(.data$classification == "consistent", na.rm = TRUE),
    n_weak_or_mixed_core = sum(.data$classification == "weak_or_mixed_core", na.rm = TRUE),
    n_inconsistent = sum(.data$classification == "inconsistent_check_mapping_or_contrast", na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    frac_consistent = ifelse(.data$n_terms > 0, .data$n_consistent / .data$n_terms, NA_real_),
    dominant_problem = dplyr::case_when(
      .data$n_consistent >= .data$n_weak_or_mixed_core & .data$n_consistent >= .data$n_inconsistent ~ "mostly_consistent",
      .data$n_weak_or_mixed_core >= .data$n_inconsistent ~ "weak_or_mixed_core",
      TRUE ~ "inconsistent_check_mapping_or_contrast"
    ),
    interpretation = dplyr::case_when(
      .data$dominant_problem == "mostly_consistent" ~ "GSEA NES direction agrees with the direction of its leading/core proteins.",
      .data$dominant_problem == "weak_or_mixed_core" ~ "Many GO terms have mixed core proteins; inspect leading-edge genes rather than term label only.",
      TRUE ~ "Potential mapping, contrast, or sign mismatch; inspect affected terms."
    )
  )

dir_create(tables_dir)
dir_create(logs_dir)
utils::write.csv(term_audit, output_paths$term_audit, row.names = FALSE, na = "")
utils::write.csv(summary, output_paths$contrast_summary, row.names = FALSE, na = "")
utils::write.csv(ora_warning, output_paths$ora_warning, row.names = FALSE, na = "")
write_input_status(status_df, output_paths$input_status, dry_run = FALSE)

finish_script_runtime(
  runtime,
  manifest_path = output_paths$run_manifest,
  outputs = unlist(output_paths[c("term_audit", "contrast_summary", "ora_warning", "input_status")], use.names = TRUE),
  inputs = c(
    clusterProfiler_manifest = manifest_path,
    stats::setNames(unique(manifest_filtered$output_table), paste0("gsea_table_", seq_along(unique(manifest_filtered$output_table)))),
    stats::setNames(unique(manifest_filtered$input_gene_file), paste0("mapped_input_", seq_along(unique(manifest_filtered$input_gene_file))))
  ),
  status = "completed",
  notes = c(
    "Read-only diagnostic audit. No upstream enrichment, mapped, or ProTigy files are modified.",
    "NES > 0 corresponds to the positive log2fc side of the formal contrast; NES < 0 corresponds to the negative side."
  )
)

message("[INFO] GSEA/protein direction audit complete.")
message("[INFO] Audited contrasts: ", length(unique(term_audit$contrast)))
message("[INFO] Audited GSEA terms: ", nrow(term_audit))
message("[INFO] Outputs: ", tables_dir)
