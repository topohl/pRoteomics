# ================================================================
# Script: 04_differential_expression_enrichment/01b_gsea_protein_direction_audit.r
# Stage: enrichment
# Scope: dataset_specific
# Consumes: canonical clusterProfiler manifest, collapsed-gene inputs, and term-gene provenance.
# Produces: diagnostic GSEA/protein direction audit CSVs under results/tables and run_manifest.yml under results/logs.
# Dataset behavior: runs for neuron_neuropil, neuron_soma, microglia via --dataset.
# Notes: Read-only diagnostic. Does not modify ProTigy, mapped, clusterProfiler, or compareGO outputs.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "enrichment_io.R"))

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

write_csv_in_chunks <- function(x, path, chunk_size = 100000L) {
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be one positive integer.", call. = FALSE)
  }
  write_chunk <- function(chunk, append, column_names) {
    utils::write.table(
      chunk, file = path, sep = ",", dec = ".", quote = TRUE,
      qmethod = "double", row.names = FALSE, col.names = column_names,
      append = append, na = ""
    )
  }

  if (!nrow(x)) {
    write_chunk(x, append = FALSE, column_names = TRUE)
    return(invisible(path))
  }

  starts <- seq.int(1L, nrow(x), by = chunk_size)
  for (i in seq_along(starts)) {
    first <- starts[[i]]
    last <- min(first + chunk_size - 1L, nrow(x))
    write_chunk(x[first:last, , drop = FALSE], append = i > 1L, column_names = i == 1L)
  }
  invisible(path)
}

canonical_manifest_path <- function(dataset) {
  canonical_clusterprofiler_manifest_path(dataset)
}

sign_label <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x > 0] <- "positive"
  out[!is.na(x) & x < 0] <- "negative"
  out[!is.na(x) & x == 0] <- "zero"
  out
}

parse_comparison_sides <- function(comparison) {
  comparison <- as.character(comparison)
  out <- data.frame(
    formal_contrast = comparison,
    positive_side_label = rep("positive_log2fc_side", length(comparison)),
    negative_side_label = rep("negative_log2fc_side", length(comparison)),
    stringsAsFactors = FALSE
  )
  if (!length(comparison)) return(out)

  comp <- tolower(comparison)
  compact <- gsub("[^a-z0-9.]+", "_", comp)
  has_label_pair <- function(x, left, right) {
    left_right <- paste0("(^|_)[a-z0-9.]*", left, "(_|$).*[a-z0-9.]*", right, "(_|$)")
    right_left <- paste0("(^|_)[a-z0-9.]*", right, "(_|$).*[a-z0-9.]*", left, "(_|$)")
    grepl(left_right, x) | grepl(right_left, x)
  }
  res_con <- !is.na(compact) & (grepl("2\\.over\\.1", compact) | has_label_pair(compact, "res", "con"))
  sus_res <- !is.na(compact) & (grepl("3\\.over\\.2", compact) | has_label_pair(compact, "sus", "res"))
  sus_con <- !is.na(compact) & (grepl("3\\.over\\.1", compact) | has_label_pair(compact, "sus", "con"))

  out$positive_side_label[res_con] <- "RES"
  out$negative_side_label[res_con] <- "CON"
  out$positive_side_label[sus_res] <- "SUS"
  out$negative_side_label[sus_res] <- "RES"
  out$positive_side_label[sus_con] <- "SUS"
  out$negative_side_label[sus_con] <- "CON"
  out
}

make_status_summary_rows <- function(status) {
  status_only <- status |>
    dplyr::filter(.data$analysis_status != "success_with_terms") |>
    dplyr::select(dataset, ontology, comparison, analysis_status)
  status_only <- dplyr::bind_cols(
    status_only,
    parse_comparison_sides(status_only$comparison)
  )
  n_status <- nrow(status_only)
  interpretation <- rep(NA_character_, n_status)
  failed <- !is.na(status_only$analysis_status) & status_only$analysis_status == "failed"
  completed <- !is.na(status_only$analysis_status) & !failed
  interpretation[failed] <- "Upstream GSEA analysis failed; excluded from direction audit."
  interpretation[completed] <- "Upstream GSEA completed successfully with zero terms."

  status_only |>
    dplyr::transmute(
      dataset, ontology, contrast = .data$comparison,
      formal_contrast, positive_side_label, negative_side_label,
      analysis_status,
      n_terms = rep(0L, n_status),
      n_consistent = rep(0L, n_status),
      n_weak_or_mixed_core = rep(0L, n_status),
      n_inconsistent = rep(0L, n_status),
      frac_consistent = rep(NA_real_, n_status),
      dominant_problem = rep(NA_character_, n_status),
      interpretation = interpretation
    )
}

assert_bind_rows_compatible <- function(left, right, context) {
  if (!identical(names(left), names(right))) {
    stop(context, " tables have different column names or order.", call. = FALSE)
  }
  left_types <- vapply(left, typeof, character(1))
  right_types <- vapply(right, typeof, character(1))
  if (!identical(left_types, right_types)) {
    incompatible <- names(left_types)[left_types != right_types]
    detail <- paste0(incompatible, " (", left_types[incompatible], " vs ", right_types[incompatible], ")")
    stop(context, " tables have incompatible column types: ", paste(detail, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
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
    successful <- manifest_filtered$analysis_status %in% c("success_with_terms", "success_zero_terms")
    required_paths <- c("output_table", "collapsed_gene_input_file", "term_gene_provenance_file")
    for (column in required_paths) {
      rows <- c(rows, lapply(which(successful), function(i) {
        input_status_row(
          paste0(column, ":", manifest_filtered$comparison[[i]]),
          manifest_filtered[[column]][[i]], dataset = DATASET, required = TRUE
        )
      }))
    }
  }
  do.call(rbind, rows)
}

manifest_path <- canonical_manifest_path(DATASET)
manifest_exists <- file.exists(manifest_path)
manifest <- if (manifest_exists) {
  read_canonical_clusterprofiler_manifest(
    manifest_path, DATASET, strict = TRUE,
    require_files = !isTRUE(runtime$dry_run)
  )
} else {
  data.frame()
}
manifest_filtered <- if (nrow(manifest)) {
  manifest[
    toupper(as.character(manifest$ontology)) == ONTOLOGY & manifest$result_type == "GSEA_GO",
    , drop = FALSE
  ]
} else {
  data.frame()
}
if (nrow(manifest_filtered)) {
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
  dry_run_line("Filtered canonical GSEA_GO rows", nrow(manifest_filtered))
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
if (!nrow(manifest_filtered)) {
  stop("No manifest rows matched dataset=", DATASET, ", ontology=", ONTOLOGY,
       ", result_type=GSEA_GO",
       if (nzchar(COMPARISON_FILTER)) paste0(", comparison containing '", COMPARISON_FILTER, "'") else "",
       ".", call. = FALSE)
}

bundle <- read_canonical_clusterprofiler_bundle(
  manifest_path, DATASET, result_types = "GSEA_GO", ontology = ONTOLOGY,
  comparisons = manifest_filtered$comparison, strict = TRUE
)

top_abs_threshold <- load_threshold()
identity_columns <- c("dataset", "comparison", "result_type", "ontology")
term_identity_columns <- c(identity_columns, "term_id")
manifest_identity <- bundle$manifest |>
  dplyr::select(dplyr::all_of(c(identity_columns, "analysis_status", "n_terms", "route_category",
    "route_unit", "output_table", "collapsed_gene_input_file", "term_gene_provenance_file")))

term_meta <- bundle$terms |>
  dplyr::transmute(
    dataset, comparison, result_type, ontology, term_id,
    NES = suppressWarnings(as.numeric(.data$NES)),
    p.adjust = suppressWarnings(as.numeric(.data$`p.adjust`))
  ) |>
  dplyr::left_join(manifest_identity, by = identity_columns)

if (nrow(term_meta)) {
  missing_provenance <- dplyr::anti_join(
    term_meta |> dplyr::distinct(dplyr::across(dplyr::all_of(term_identity_columns))),
    bundle$provenance |> dplyr::distinct(dplyr::across(dplyr::all_of(term_identity_columns))),
    by = term_identity_columns
  )
  if (nrow(missing_provenance)) {
    stop("Canonical GSEA term is missing required term-gene provenance.", call. = FALSE)
  }
}

contributors <- join_term_provenance_to_collapsed_genes(
  bundle$provenance |>
    dplyr::filter(.data$core_enrichment_member, .data$gene_level_claim_allowed),
  bundle$collapsed,
  strict = TRUE
) |>
  dplyr::left_join(term_meta, by = term_identity_columns)

term_gene_values <- contributors |>
  dplyr::distinct(dplyr::across(dplyr::all_of(c(term_identity_columns, "official_gene_symbol"))),
    .keep_all = TRUE)
term_metrics <- term_gene_values |>
  dplyr::group_by(dplyr::across(dplyr::all_of(term_identity_columns))) |>
  dplyr::summarise(
    NES = dplyr::first(.data$NES), p.adjust = dplyr::first(.data$p.adjust),
    n_core = dplyr::n_distinct(.data$official_gene_symbol),
    n_core_mapped = sum(is.finite(.data$collapsed_logfc)),
    median_core_log2fc = ifelse(all(!is.finite(.data$collapsed_logfc)), NA_real_, stats::median(.data$collapsed_logfc, na.rm = TRUE)),
    mean_core_log2fc = ifelse(all(!is.finite(.data$collapsed_logfc)), NA_real_, mean(.data$collapsed_logfc, na.rm = TRUE)),
    frac_core_positive = ifelse(any(is.finite(.data$collapsed_logfc)), mean(.data$collapsed_logfc > 0, na.rm = TRUE), NA_real_),
    frac_core_negative = ifelse(any(is.finite(.data$collapsed_logfc)), mean(.data$collapsed_logfc < 0, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

term_audit <- contributors |>
  dplyr::select(-dplyr::any_of(c("NES", "p.adjust"))) |>
  dplyr::left_join(term_metrics, by = term_identity_columns)
term_audit <- dplyr::bind_cols(
  term_audit,
  parse_comparison_sides(term_audit$comparison)
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    contrast = .data$comparison,
    NES_sign = sign_label(.data$NES),
    log2fc_sign = sign_label(.data$collapsed_logfc),
    median_core_sign = sign_label(.data$median_core_log2fc),
    direction_consistent = !is.na(.data$NES_sign) && !is.na(.data$log2fc_sign) && .data$NES_sign == .data$log2fc_sign,
    direction_strength = max(.data$frac_core_positive, .data$frac_core_negative, na.rm = TRUE),
    direction_strength = dplyr::if_else(is.infinite(.data$direction_strength), NA_real_, as.numeric(.data$direction_strength), ptype = double()),
    classification = dplyr::case_when(
      .data$n_core_mapped < 3 || is.na(.data$direction_strength) || .data$direction_strength < 0.60 ~ "weak_or_mixed_core",
      sign_label(.data$NES) == sign_label(.data$median_core_log2fc) ~ "consistent",
      TRUE ~ "inconsistent_check_mapping_or_contrast"
    ),
    concordance_status = dplyr::case_when(
      is.na(.data$NES_sign) || is.na(.data$log2fc_sign) ~ "direction_unavailable",
      .data$direction_consistent ~ "concordant",
      TRUE ~ "discordant"
    ),
    NES_biological_side = side_from_sign(.data$NES_sign, .data$positive_side_label, .data$negative_side_label),
    core_median_biological_side = side_from_sign(.data$median_core_sign, .data$positive_side_label, .data$negative_side_label),
    source_manifest = bundle$manifest_source,
    collapsed_gene_input_file = .data$source_file,
    term_provenance_file = .data$term_gene_provenance_file
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(.data$dataset, .data$comparison, .data$ontology, .data$term_id,
    .data$official_gene_symbol, .data$ProteinGroupID)

term_summary <- term_metrics |>
  dplyr::left_join(manifest_identity, by = identity_columns)
term_summary <- dplyr::bind_cols(
  term_summary,
  parse_comparison_sides(term_summary$comparison)
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    contrast = .data$comparison,
    direction_strength = max(.data$frac_core_positive, .data$frac_core_negative, na.rm = TRUE),
    direction_strength = dplyr::if_else(is.infinite(.data$direction_strength), NA_real_, as.numeric(.data$direction_strength), ptype = double()),
    classification = dplyr::case_when(
      .data$n_core_mapped < 3 || is.na(.data$direction_strength) || .data$direction_strength < 0.60 ~ "weak_or_mixed_core",
      sign_label(.data$NES) == sign_label(.data$median_core_log2fc) ~ "consistent",
      TRUE ~ "inconsistent_check_mapping_or_contrast"
    )
  ) |>
  dplyr::ungroup()

summary_with_terms <- term_summary |>
  dplyr::group_by(.data$dataset, .data$ontology, .data$contrast, .data$formal_contrast,
                  .data$positive_side_label, .data$negative_side_label, .data$analysis_status) |>
  dplyr::summarise(
    n_terms = dplyr::n(),
    n_consistent = sum(.data$classification == "consistent", na.rm = TRUE),
    n_weak_or_mixed_core = sum(.data$classification == "weak_or_mixed_core", na.rm = TRUE),
    n_inconsistent = sum(.data$classification == "inconsistent_check_mapping_or_contrast", na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    frac_consistent = dplyr::if_else(.data$n_terms > 0L, .data$n_consistent / .data$n_terms, NA_real_, ptype = double()),
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

status_only <- make_status_summary_rows(bundle$status)
assert_bind_rows_compatible(summary_with_terms, status_only, "GSEA direction summary")
summary <- dplyr::bind_rows(summary_with_terms, status_only) |>
  dplyr::arrange(.data$dataset, .data$ontology, .data$contrast)

ora_warning <- bundle$manifest |>
  dplyr::filter(.data$analysis_status != "failed")
ora_warning <- dplyr::bind_cols(
  ora_warning,
  parse_comparison_sides(ora_warning$comparison)
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    values = list(bundle$collapsed$collapsed_logfc[bundle$collapsed$comparison == .data$comparison & bundle$collapsed$ontology == .data$ontology])
  ) |>
  dplyr::transmute(
    dataset, ontology, contrast = .data$comparison, formal_contrast,
    positive_side_label, negative_side_label,
    n_positive_log2fc = sum(unlist(.data$values) > 0, na.rm = TRUE),
    n_negative_log2fc = sum(unlist(.data$values) < 0, na.rm = TRUE),
    top_abs_log2fc_threshold = top_abs_threshold,
    n_top_abs_log2fc = sum(abs(unlist(.data$values)) > top_abs_threshold, na.rm = TRUE),
    n_top_positive = sum(unlist(.data$values) > top_abs_threshold, na.rm = TRUE),
    n_top_negative = sum(unlist(.data$values) < -top_abs_threshold, na.rm = TRUE),
    note = "Directional diagnostic uses canonical collapsed log2fc values; enrichment results are unchanged.",
    collapsed_gene_input_file
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(.data$dataset, .data$ontology, .data$contrast)

dir_create(tables_dir)
dir_create(logs_dir)
write_csv_in_chunks(term_audit, output_paths$term_audit)
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
    stats::setNames(unique(manifest_filtered$collapsed_gene_input_file), paste0("collapsed_gene_input_", seq_along(unique(manifest_filtered$collapsed_gene_input_file)))),
    stats::setNames(unique(manifest_filtered$term_gene_provenance_file), paste0("term_gene_provenance_", seq_along(unique(manifest_filtered$term_gene_provenance_file))))
  ),
  status = "completed",
  notes = c(
    "Read-only diagnostic audit. No upstream enrichment, mapped, or ProTigy files are modified.",
    "NES > 0 corresponds to the positive log2fc side of the formal contrast; NES < 0 corresponds to the negative side."
  )
)

message("[INFO] GSEA/protein direction audit complete.")
message("[INFO] Audited contrasts: ", length(unique(summary$contrast)))
message("[INFO] Audited GSEA terms: ", nrow(term_summary))
message("[INFO] Outputs: ", tables_dir)
