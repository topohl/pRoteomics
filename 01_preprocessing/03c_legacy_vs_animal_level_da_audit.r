#!/usr/bin/env Rscript

# Direct audit of historical hemisphere-level versus corrected animal-level ProTigy DA statistics.
# This script reads the six ProTigy statistical GCTs directly and does not run mapping or enrichment.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "protigy_stat_gct_utils.R"))

SCRIPT_ID <- "01_preprocessing/03c_legacy_vs_animal_level_da_audit.r"
MODULE_ID <- "01_preprocessing/03c_legacy_vs_animal_level_da_audit"
runtime <- init_script_runtime(
  script = SCRIPT_ID,
  stage = "diagnostic",
  default_dataset = "all",
  allow_all = TRUE
)
datasets <- if (identical(runtime$dataset, "all")) valid_datasets() else runtime$dataset

legacy_root <- resolve_protigy_root(
  Sys.getenv("PROTEOMICS_LEGACY_GCT_ROOT", unset = ""),
  path_processed("01_preprocessing", "protigy_output")
)
animal_root <- resolve_protigy_root(
  Sys.getenv("PROTEOMICS_ANIMAL_LEVEL_GCT_ROOT", unset = ""),
  path_processed("01_preprocessing", "protigy_output_animal_level")
)

table_root <- path_results("tables", MODULE_ID)
report_root <- path_results("reports", MODULE_ID)
log_root <- path_results("logs", MODULE_ID)

safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L) return(NA_real_)
  if (stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

format_num <- function(x, digits = 3L) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

top_rank_overlap <- function(legacy, animal, n) {
  legacy_ids <- legacy$protein_id[is.finite(legacy$t) & legacy$rank <= n]
  animal_ids <- animal$protein_id[is.finite(animal$t) & animal$rank <= n]
  length(intersect(legacy_ids, animal_ids))
}

build_unmatched_table <- function(dataset, legacy, animal, matches) {
  legacy_desc <- stats::setNames(legacy$description, legacy$ids)
  animal_desc <- stats::setNames(animal$description, animal$ids)
  rbind(
    data.frame(
      dataset = rep(dataset, length(matches$legacy_only)),
      source_only = rep("legacy_only", length(matches$legacy_only)),
      protein_id = matches$legacy_only,
      description = unname(legacy_desc[matches$legacy_only]),
      stringsAsFactors = FALSE
    ),
    data.frame(
      dataset = rep(dataset, length(matches$animal_only)),
      source_only = rep("animal_level_only", length(matches$animal_only)),
      protein_id = matches$animal_only,
      description = unname(animal_desc[matches$animal_only]),
      stringsAsFactors = FALSE
    )
  )
}

build_comparison_validation <- function(gct, source_label) {
  diag <- protigy_gct_diagnostic_row(gct)
  out <- gct$comparison_validation
  out$source <- source_label
  out$gct_path <- gct$path
  out$gct_sha256 <- gct$sha256
  out$nrmat <- diag$nrmat
  out$ncmat <- diag$ncmat
  out$nrhd <- diag$nrhd
  out$nchd <- diag$nchd
  out$n_protein_ids <- diag$n_protein_ids
  out$protein_ids_unique <- diag$protein_ids_unique
  out$description_present <- diag$description_present
  out$description_aligned <- diag$description_aligned
  out$metrics_found <- diag$metrics_found
  out$metrics_missing <- diag$metrics_missing
  out$n_missing_numeric_values <- diag$n_missing_numeric_values
  out$n_non_numeric_values <- diag$n_non_numeric_values
  out$physical_dimension_match <- diag$physical_dimension_match
  out
}

audit_dataset <- function(dataset) {
  legacy_resolved <- resolve_single_protigy_gct(legacy_root, dataset)
  animal_resolved <- resolve_single_protigy_gct(animal_root, dataset)
  legacy <- read_protigy_stat_gct(legacy_resolved$path, dataset, strict_primary = FALSE)
  animal <- read_protigy_stat_gct(animal_resolved$path, dataset, strict_primary = TRUE)
  animal_contract <- validate_corrected_protigy_gct_contract(animal)

  expected_proteins <- c(neuron_neuropil = 5045L, neuron_soma = 5529L)
  if (dataset %in% names(expected_proteins) && length(animal$ids) != expected_proteins[[dataset]]) {
    stop(
      "Corrected ", dataset, " protein count mismatch. Expected ", expected_proteins[[dataset]],
      "; observed ", length(animal$ids), ".",
      call. = FALSE
    )
  }

  legacy_valid <- legacy$comparison_validation$comparison[legacy$comparison_validation$valid_stress_comparison]
  animal_valid <- animal$comparison_validation$comparison[animal$comparison_validation$valid_stress_comparison]
  if (!setequal(legacy_valid, animal_valid)) {
    stop(
      "Legacy and animal-level within-unit stress comparison sets differ for ", dataset,
      ". Legacy-only: ", paste(setdiff(legacy_valid, animal_valid), collapse = ", "),
      "; animal-only: ", paste(setdiff(animal_valid, legacy_valid), collapse = ", "),
      call. = FALSE
    )
  }
  comparisons <- sort(animal_valid, method = "radix")
  metrics <- c("logFC", "t", "P.Value", "adj.P.Val", "B", "AveExpr")
  legacy_long <- extract_protigy_metric_long(legacy, comparisons, metrics)
  animal_long <- extract_protigy_metric_long(animal, comparisons, metrics)
  legacy_long$rank <- ave(
    seq_len(nrow(legacy_long)), legacy_long$comparison,
    FUN = function(i) deterministic_signed_rank(legacy_long$t[i], legacy_long$protein_id[i])
  )
  animal_long$rank <- ave(
    seq_len(nrow(animal_long)), animal_long$comparison,
    FUN = function(i) deterministic_signed_rank(animal_long$t[i], animal_long$protein_id[i])
  )

  matches <- match_protigy_proteins(legacy$ids, animal$ids)
  unmatched <- build_unmatched_table(dataset, legacy, animal, matches)

  merged <- merge(
    legacy_long, animal_long,
    by = c("dataset", "spatial_unit", "comparison", "biological_contrast", "protein_id"),
    suffixes = c("_legacy", "_animal"),
    all = FALSE, sort = FALSE
  )
  merged <- merged[order(merged$comparison, merged$protein_id, method = "radix"), , drop = FALSE]
  protein <- data.frame(
    dataset = merged$dataset,
    spatial_unit = merged$spatial_unit,
    comparison = merged$comparison,
    biological_contrast = merged$biological_contrast,
    protein_id = merged$protein_id,
    description_legacy = merged$description_legacy,
    description_animal_level = merged$description_animal,
    legacy_logFC = merged$logFC_legacy,
    animal_logFC = merged$logFC_animal,
    delta_logFC = merged$logFC_animal - merged$logFC_legacy,
    legacy_t = merged$t_legacy,
    animal_t = merged$t_animal,
    delta_t = merged$t_animal - merged$t_legacy,
    legacy_P.Value = merged$P.Value_legacy,
    animal_P.Value = merged$P.Value_animal,
    legacy_adj.P.Val = merged$adj.P.Val_legacy,
    animal_adj.P.Val = merged$adj.P.Val_animal,
    legacy_B = merged$B_legacy,
    animal_B = merged$B_animal,
    legacy_AveExpr = merged$AveExpr_legacy,
    animal_AveExpr = merged$AveExpr_animal,
    legacy_significant_FDR05 = is.finite(merged$adj.P.Val_legacy) & merged$adj.P.Val_legacy <= 0.05,
    animal_significant_FDR05 = is.finite(merged$adj.P.Val_animal) & merged$adj.P.Val_animal <= 0.05,
    effect_sign_same = ifelse(
      is.finite(merged$logFC_legacy) & is.finite(merged$logFC_animal),
      sign(merged$logFC_legacy) == sign(merged$logFC_animal), NA
    ),
    rank_legacy = as.integer(merged$rank_legacy),
    rank_animal = as.integer(merged$rank_animal),
    rank_delta = as.integer(merged$rank_animal - merged$rank_legacy),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  summary_rows <- lapply(comparisons, function(comp) {
    legacy_comp <- legacy_long[legacy_long$comparison == comp, , drop = FALSE]
    animal_comp <- animal_long[animal_long$comparison == comp, , drop = FALSE]
    matched_comp <- protein[protein$comparison == comp, , drop = FALSE]
    legacy_sig <- legacy_comp$protein_id[is.finite(legacy_comp$adj.P.Val) & legacy_comp$adj.P.Val <= 0.05]
    animal_sig <- animal_comp$protein_id[is.finite(animal_comp$adj.P.Val) & animal_comp$adj.P.Val <= 0.05]
    union_sig <- union(legacy_sig, animal_sig)
    intersection_sig <- intersect(legacy_sig, animal_sig)
    jaccard_status <- if (!length(union_sig)) "both_sets_empty_defined_as_1" else "observed_union"
    data.frame(
      dataset = dataset,
      spatial_unit = matched_comp$spatial_unit[[1]],
      comparison = comp,
      biological_contrast = matched_comp$biological_contrast[[1]],
      n_legacy_proteins = nrow(legacy_comp),
      n_animal_level_proteins = nrow(animal_comp),
      n_matched_proteins = nrow(matched_comp),
      n_legacy_only = length(matches$legacy_only),
      n_animal_only = length(matches$animal_only),
      pearson_logFC = safe_cor(matched_comp$legacy_logFC, matched_comp$animal_logFC, "pearson"),
      spearman_logFC = safe_cor(matched_comp$legacy_logFC, matched_comp$animal_logFC, "spearman"),
      spearman_t = safe_cor(matched_comp$legacy_t, matched_comp$animal_t, "spearman"),
      effect_direction_concordance = mean(matched_comp$effect_sign_same, na.rm = TRUE),
      legacy_FDR05_count = length(legacy_sig),
      animal_level_FDR05_count = length(animal_sig),
      retained_significant = length(intersection_sig),
      lost_significant = length(setdiff(legacy_sig, animal_sig)),
      gained_significant = length(setdiff(animal_sig, legacy_sig)),
      significant_set_union_size = length(union_sig),
      significant_set_intersection_size = length(intersection_sig),
      significant_set_Jaccard = if (!length(union_sig)) 1 else length(intersection_sig) / length(union_sig),
      significant_set_Jaccard_status = jaccard_status,
      top20_t_rank_overlap = top_rank_overlap(legacy_comp, animal_comp, 20L),
      top50_t_rank_overlap = top_rank_overlap(legacy_comp, animal_comp, 50L),
      top100_t_rank_overlap = top_rank_overlap(legacy_comp, animal_comp, 100L),
      median_abs_delta_logFC = safe_median(abs(matched_comp$delta_logFC)),
      median_abs_delta_t = safe_median(abs(matched_comp$delta_t)),
      max_abs_delta_logFC = safe_max(abs(matched_comp$delta_logFC)),
      max_abs_delta_t = safe_max(abs(matched_comp$delta_t)),
      stringsAsFactors = FALSE
    )
  })
  contrast <- do.call(rbind, summary_rows)
  rownames(contrast) <- NULL

  validation <- rbind(
    build_comparison_validation(legacy, "legacy"),
    build_comparison_validation(animal, "animal_level")
  )
  validation$expected_corrected_comparison_count <- animal_contract$expected_comparison_count
  validation$corrected_contract_status <- animal_contract$status

  dataset_dir <- file.path(table_root, dataset)
  if (!runtime$dry_run) {
    dir_create(dataset_dir)
    utils::write.csv(protein, file.path(dataset_dir, "protein_level_comparison.csv"), row.names = FALSE, na = "")
    utils::write.csv(contrast, file.path(dataset_dir, "contrast_level_summary.csv"), row.names = FALSE, na = "")
    utils::write.csv(unmatched, file.path(dataset_dir, "unmatched_proteins.csv"), row.names = FALSE, na = "")
    utils::write.csv(validation, file.path(dataset_dir, "comparison_validation.csv"), row.names = FALSE, na = "")
  }

  list(
    dataset = dataset,
    legacy = legacy,
    animal = animal,
    contract = animal_contract,
    legacy_long = legacy_long,
    animal_long = animal_long,
    protein = protein,
    contrast = contrast,
    unmatched = unmatched,
    validation = validation
  )
}

dataset_summary_row <- function(result) {
  protein <- result$protein
  contrast <- result$contrast
  legacy_unique_sig <- unique(result$legacy_long$protein_id[
    is.finite(result$legacy_long$adj.P.Val) & result$legacy_long$adj.P.Val <= 0.05
  ])
  animal_unique_sig <- unique(result$animal_long$protein_id[
    is.finite(result$animal_long$adj.P.Val) & result$animal_long$adj.P.Val <= 0.05
  ])
  spatial_t <- aggregate(abs(protein$delta_t), list(spatial_unit = protein$spatial_unit), safe_median)
  names(spatial_t)[[2]] <- "median_abs_delta_t"
  spatial_fc <- aggregate(abs(protein$delta_logFC), list(spatial_unit = protein$spatial_unit), safe_median)
  names(spatial_fc)[[2]] <- "median_abs_delta_logFC"
  spatial_t <- spatial_t[order(-spatial_t$median_abs_delta_t, spatial_t$spatial_unit, method = "radix"), , drop = FALSE]
  spatial_fc <- spatial_fc[order(-spatial_fc$median_abs_delta_logFC, spatial_fc$spatial_unit, method = "radix"), , drop = FALSE]
  sus <- contrast[contrast$biological_contrast == "SUS_vs_RES", , drop = FALSE]
  context <- contrast[contrast$biological_contrast != "SUS_vs_RES", , drop = FALSE]
  unmatched_legacy <- sum(result$unmatched$source_only == "legacy_only")
  unmatched_animal <- sum(result$unmatched$source_only == "animal_level_only")

  data.frame(
    dataset = result$dataset,
    legacy_gct = result$legacy$path,
    legacy_sha256 = result$legacy$sha256,
    animal_level_gct = result$animal$path,
    animal_level_sha256 = result$animal$sha256,
    legacy_proteins = length(result$legacy$ids),
    animal_level_proteins = length(result$animal$ids),
    matched_proteins = length(intersect(result$legacy$ids, result$animal$ids)),
    legacy_only_proteins = unmatched_legacy,
    animal_level_only_proteins = unmatched_animal,
    legacy_detected_comparisons = nrow(result$legacy$comparison_validation),
    legacy_valid_stress_comparisons = sum(result$legacy$comparison_validation$valid_stress_comparison),
    animal_level_detected_comparisons = nrow(result$animal$comparison_validation),
    animal_level_valid_stress_comparisons = sum(result$animal$comparison_validation$valid_stress_comparison),
    expected_animal_level_comparisons = result$contract$expected_comparison_count,
    n_spatial_units = length(result$contract$expected_spatial_units),
    median_abs_delta_logFC = safe_median(abs(protein$delta_logFC)),
    max_abs_delta_logFC = safe_max(abs(protein$delta_logFC)),
    median_abs_delta_t = safe_median(abs(protein$delta_t)),
    max_abs_delta_t = safe_max(abs(protein$delta_t)),
    median_contrast_pearson_logFC = safe_median(contrast$pearson_logFC),
    median_contrast_spearman_logFC = safe_median(contrast$spearman_logFC),
    median_contrast_spearman_t = safe_median(contrast$spearman_t),
    median_abs_t_rank_delta = safe_median(abs(protein$rank_delta)),
    max_abs_t_rank_delta = safe_max(abs(protein$rank_delta)),
    median_top20_t_rank_overlap = safe_median(contrast$top20_t_rank_overlap),
    median_top50_t_rank_overlap = safe_median(contrast$top50_t_rank_overlap),
    median_top100_t_rank_overlap = safe_median(contrast$top100_t_rank_overlap),
    overall_effect_direction_concordance = mean(protein$effect_sign_same, na.rm = TRUE),
    legacy_FDR05_protein_contrast_count = sum(contrast$legacy_FDR05_count),
    animal_level_FDR05_protein_contrast_count = sum(contrast$animal_level_FDR05_count),
    legacy_unique_FDR05_proteins = length(legacy_unique_sig),
    animal_level_unique_FDR05_proteins = length(animal_unique_sig),
    retained_significant_count = sum(contrast$retained_significant),
    lost_significant_count = sum(contrast$lost_significant),
    gained_significant_count = sum(contrast$gained_significant),
    most_changed_spatial_unit_by_delta_t = spatial_t$spatial_unit[[1]],
    most_changed_spatial_unit_median_abs_delta_t = spatial_t$median_abs_delta_t[[1]],
    most_changed_spatial_unit_by_delta_logFC = spatial_fc$spatial_unit[[1]],
    most_changed_spatial_unit_median_abs_delta_logFC = spatial_fc$median_abs_delta_logFC[[1]],
    sus_res_legacy_FDR05_count = sum(sus$legacy_FDR05_count),
    sus_res_animal_level_FDR05_count = sum(sus$animal_level_FDR05_count),
    sus_res_retained_significant = sum(sus$retained_significant),
    sus_res_lost_significant = sum(sus$lost_significant),
    sus_res_gained_significant = sum(sus$gained_significant),
    sus_res_median_abs_delta_logFC = safe_median(sus$median_abs_delta_logFC),
    contextual_median_abs_delta_logFC = safe_median(context$median_abs_delta_logFC),
    sus_res_minus_contextual_median_abs_delta_logFC = safe_median(sus$median_abs_delta_logFC) - safe_median(context$median_abs_delta_logFC),
    sus_res_median_abs_delta_t = safe_median(sus$median_abs_delta_t),
    contextual_median_abs_delta_t = safe_median(context$median_abs_delta_t),
    sus_res_minus_contextual_median_abs_delta_t = safe_median(sus$median_abs_delta_t) - safe_median(context$median_abs_delta_t),
    sus_res_median_spearman_t = safe_median(sus$spearman_t),
    contextual_median_spearman_t = safe_median(context$spearman_t),
    sus_res_median_abs_t_rank_delta = safe_median(abs(protein$rank_delta[protein$biological_contrast == "SUS_vs_RES"])),
    contextual_median_abs_t_rank_delta = safe_median(abs(protein$rank_delta[protein$biological_contrast != "SUS_vs_RES"])),
    stringsAsFactors = FALSE
  )
}

if (runtime$dry_run) {
  dry_run_line("Script", SCRIPT_ID)
  dry_run_line("Datasets", paste(datasets, collapse = ", "))
  dry_run_line("Legacy root", legacy_root, if (dir.exists(legacy_root)) "PASS" else "FAIL")
  dry_run_line("Animal-level root", animal_root, if (dir.exists(animal_root)) "PASS" else "FAIL")
  for (dataset in datasets) {
    legacy <- resolve_single_protigy_gct(legacy_root, dataset)
    animal <- resolve_single_protigy_gct(animal_root, dataset)
    dry_run_line(paste(dataset, "legacy GCT"), legacy$path, "PASS")
    dry_run_line(paste(dataset, "animal-level GCT"), animal$path, "PASS")
  }
  dry_run_line("Table root", table_root)
  dry_run_line("Report", file.path(report_root, "legacy_vs_animal_level_summary.md"))
  dry_run_line("Log root", log_root)
  quit(status = 0, save = "no")
}

results <- lapply(datasets, audit_dataset)
global_summary <- do.call(rbind, lapply(results, dataset_summary_row))
rownames(global_summary) <- NULL

dir_create(table_root)
dir_create(report_root)
dir_create(log_root)
global_summary_path <- file.path(table_root, "legacy_vs_animal_level_global_summary.csv")
utils::write.csv(global_summary, global_summary_path, row.names = FALSE, na = "")

source_manifest <- do.call(rbind, lapply(results, function(x) {
  rbind(
    transform(protigy_gct_diagnostic_row(x$legacy), source = "legacy"),
    transform(protigy_gct_diagnostic_row(x$animal), source = "animal_level")
  )
}))
source_manifest <- source_manifest[, c("dataset", "source", setdiff(names(source_manifest), c("dataset", "source"))), drop = FALSE]
source_manifest_path <- file.path(log_root, "source_gct_manifest.csv")
utils::write.csv(source_manifest, source_manifest_path, row.names = FALSE, na = "")
write_session_info(file.path(log_root, "sessionInfo.txt"))

report_path <- file.path(report_root, "legacy_vs_animal_level_summary.md")
report <- c(
  "# Legacy versus animal-level ProTigy DA audit",
  "",
  "## Technical summary",
  "",
  "This audit compares the historical hemisphere-level and corrected animal-level ProTigy statistical-result GCTs directly. Proteins are matched only by the shared ProTigy `id`; comparisons are canonicalized from historical `.over.` and corrected `_over_` syntax and restricted to within-unit 2/1, 3/2, and 3/1 stress contrasts. No mapping, enrichment, WGCNA, or EWCE output was read or changed.",
  "",
  "The principal rank statistic is moderated `t`, matching the current GSEA statistic-selection contract, which prefers `t` over log fold-change.",
  "",
  "## Dataset-level findings",
  "",
  "| Dataset | Proteins legacy / animal / matched | Contrasts | Median abs delta logFC | Median abs delta t | Median Spearman t | Median abs rank shift | Direction concordance | FDR<=0.05 rows legacy / animal | Unique DAPs legacy / animal | Most changed unit by delta t |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
  vapply(seq_len(nrow(global_summary)), function(i) {
    x <- global_summary[i, ]
    paste0(
      "| ", x$dataset, " | ", x$legacy_proteins, " / ", x$animal_level_proteins, " / ", x$matched_proteins,
      " | ", x$animal_level_valid_stress_comparisons,
      " | ", format_num(x$median_abs_delta_logFC),
      " | ", format_num(x$median_abs_delta_t),
      " | ", format_num(x$median_contrast_spearman_t),
      " | ", format_num(x$median_abs_t_rank_delta, 0L),
      " | ", format_num(x$overall_effect_direction_concordance),
      " | ", x$legacy_FDR05_protein_contrast_count, " / ", x$animal_level_FDR05_protein_contrast_count,
      " | ", x$legacy_unique_FDR05_proteins, " / ", x$animal_level_unique_FDR05_proteins,
      " | ", x$most_changed_spatial_unit_by_delta_t, " (", format_num(x$most_changed_spatial_unit_median_abs_delta_t), ") |"
    )
  }, character(1)),
  "",
  "Counts of significant proteins are protein-by-contrast counts summed over the validated spatial contrasts; they are not unique-protein counts across the dataset.",
  "",
  "## SUS versus RES",
  "",
  "| Dataset | FDR<=0.05 legacy / animal | Retained / lost / gained | Median abs delta logFC | Contextual median | Median abs delta t | Contextual median | Median Spearman t |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(global_summary)), function(i) {
    x <- global_summary[i, ]
    paste0(
      "| ", x$dataset,
      " | ", x$sus_res_legacy_FDR05_count, " / ", x$sus_res_animal_level_FDR05_count,
      " | ", x$sus_res_retained_significant, " / ", x$sus_res_lost_significant, " / ", x$sus_res_gained_significant,
      " | ", format_num(x$sus_res_median_abs_delta_logFC),
      " | ", format_num(x$contextual_median_abs_delta_logFC),
      " | ", format_num(x$sus_res_median_abs_delta_t),
      " | ", format_num(x$contextual_median_abs_delta_t),
      " | ", format_num(x$sus_res_median_spearman_t), " |"
    )
  }, character(1)),
  "",
  "Contextual values pool RES-vs-CON and SUS-vs-CON contrast summaries. The table reports numerical differences only; it does not test whether the rerun effect differs by biological contrast.",
  "",
  "## Scope and validation",
  "",
  paste0("- Corrected comparison counts: ", paste(vapply(results, function(x) paste0(x$dataset, "=", x$contract$expected_comparison_count), character(1)), collapse = "; "), "."),
  paste0("- Corrected protein counts: ", paste(vapply(results, function(x) paste0(x$dataset, "=", length(x$animal$ids)), character(1)), collapse = "; "), "."),
  paste0("- Unmatched IDs (legacy-only / animal-only): ", paste(vapply(results, function(x) paste0(x$dataset, "=", sum(x$unmatched$source_only == "legacy_only"), "/", sum(x$unmatched$source_only == "animal_level_only")), character(1)), collapse = "; "), "."),
  "- Every corrected spatial unit contains exactly 2/1, 3/2, and 3/1, with no cross-unit comparisons, duplicate metric/comparison fields, duplicate protein IDs, missing required DA fields, or non-numeric required statistics.",
  "- Historical GCTs are audited as found. Historical `B` and `significant` fields are absent, so `legacy_B` is left missing and significance is calculated directly from historical `adj.P.Val <= 0.05`.",
  "- A significant-set Jaccard is defined as 1 when both sets are empty; `significant_set_Jaccard_status` marks those rows explicitly.",
  "",
  "## Methods and limitations",
  "",
  "Ranks are deterministic descending signed moderated-t ranks, with protein ID as the tie-breaker. Correlations use matched proteins with finite values. DAP retention, loss, gain, union, and intersection are calculated from each source's complete protein-ID set, including unmatched IDs. This is a descriptive rerun audit and does not establish a new inferential test of the difference between ProTigy analyses.",
  "",
  "The historical neuropil GCT declares fewer physical fields than are present and contains repeated physical metric/comparison fields. The reader records this discrepancy and, consistently with the historical extractor, uses the first physical occurrence of each metric/comparison field.",
  "",
  "## Output inventory",
  "",
  paste0("- Global dataset summary: `", normalizePath(global_summary_path, winslash = "/", mustWork = FALSE), "`"),
  paste0("- Source GCT manifest and hashes: `", normalizePath(source_manifest_path, winslash = "/", mustWork = FALSE), "`"),
  paste0("- Dataset audit tables: `", normalizePath(table_root, winslash = "/", mustWork = FALSE), "/<dataset>/`"),
  "",
  "## Next step boundary",
  "",
  "This audit stops before ID mapping and enrichment. Promotion should use the separate `gct_extractR_animal_level` namespace first; canonical pipeline contracts remain unchanged until a separately authorized downstream migration."
)
writeLines(report, report_path, useBytes = TRUE)

finish_script_runtime(
  runtime,
  manifest_path = file.path(log_root, "run_manifest.yml"),
  inputs = stats::setNames(source_manifest$path, paste(source_manifest$dataset, source_manifest$source, sep = "_")),
  outputs = c(
    global_summary = global_summary_path,
    report = report_path,
    source_gct_manifest = source_manifest_path,
    dataset_tables = table_root
  ),
  status = "completed",
  notes = c(
    "Direct ProTigy statistical-GCT comparison only.",
    "Protein identity is exact GCT id; no Description or gene-symbol fallback.",
    "Principal ranking statistic is moderated t, matching the current GSEA preference.",
    "No mapping, enrichment, WGCNA, or EWCE stage was run."
  )
)

message("Legacy versus animal-level DA audit completed: ", global_summary_path)
