#!/usr/bin/env Rscript

# Additive before/after audit for the finite-df CI and ontology-theme
# correction. This script never fits a model or recalculates a source p/FDR.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source(file.path("R", "paths.R"))

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
baseline_arg <- grep("^--baseline-dir=", args, value = TRUE)
baseline_dir <- if (length(baseline_arg)) {
  sub("^--baseline-dir=", "", baseline_arg[[1L]])
} else {
  path_results(
    "reviewer_audit", "gsea_wgcna_concordance_correction_20260813",
    "baseline"
  )
}
baseline_dir <- normalizePath(baseline_dir, winslash = "/", mustWork = FALSE)

table_dir <- path_results(
  "tables", "10_biological_integration",
  "gsea_wgcna_concordance_diagnostics", "global"
)
source_dir <- path_results(
  "source_data", "10_biological_integration",
  "gsea_wgcna_concordance_diagnostics", "global"
)
datasets <- c("microglia", "neuron_neuropil", "neuron_soma")
stage05_files <- c(
  "module_group_effects.csv", "supermodule_group_effects.csv",
  "WGCNA_group_effect_interaction_conditional_followup.csv"
)

if (dry_run) {
  cat("[DRY-RUN] Baseline directory:", baseline_dir, "\n")
  cat("[DRY-RUN] Reads nine Stage-05 baseline/current table pairs.\n")
  cat("[DRY-RUN] Writes finite-df CI and protected-value audits to:",
      table_dir, "and", source_dir, "\n")
  quit(save = "no", status = 0L)
}

if (!dir.exists(baseline_dir)) {
  stop("Correction baseline directory is missing: ", baseline_dir, call. = FALSE)
}

key_columns <- c(
  "dataset", "level", "endpoint_id", "effect_scope", "spatial_unit",
  "contrast", "test_type"
)
protected_base <- c("estimate", "SE", "statistic", "df_num", "df_den", "p_value")

read_pair <- function(dataset, filename) {
  before_file <- file.path(baseline_dir, paste0(dataset, "__", filename))
  after_file <- path_results(
    "tables", "06_modules_WGCNA", "group_effects", dataset, filename
  )
  if (!file.exists(before_file) || !file.exists(after_file)) {
    stop("Missing Stage-05 audit pair: ", before_file, " / ", after_file,
         call. = FALSE)
  }
  before <- read_csv(before_file, show_col_types = FALSE, progress = FALSE,
                     guess_max = Inf)
  after <- read_csv(after_file, show_col_types = FALSE, progress = FALSE,
                    guess_max = Inf)
  missing_keys <- setdiff(key_columns, intersect(names(before), names(after)))
  if (length(missing_keys)) {
    stop("Stage-05 audit key is incomplete in ", filename, ": ",
         paste(missing_keys, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(before[key_columns]) || anyDuplicated(after[key_columns])) {
    stop("Duplicated canonical Stage-05 inferential key in ", filename,
         call. = FALSE)
  }
  before_keys <- before |> select(all_of(key_columns)) |> arrange(across(everything()))
  after_keys <- after |> select(all_of(key_columns)) |> arrange(across(everything()))
  if (!identical(before_keys, after_keys)) {
    stop("Stage-05 canonical inferential keys changed in ", filename,
         call. = FALSE)
  }
  list(before = before, after = after, filename = filename)
}

pairs <- unlist(lapply(datasets, function(dataset) {
  lapply(stage05_files, function(filename) read_pair(dataset, filename))
}), recursive = FALSE)

same_value <- function(before, after, tolerance = 1e-12) {
  if (is.numeric(before) || is.numeric(after)) {
    before <- suppressWarnings(as.numeric(before))
    after <- suppressWarnings(as.numeric(after))
    return((is.na(before) & is.na(after)) |
      (is.finite(before) & is.finite(after) & abs(before - after) <= tolerance))
  }
  (is.na(before) & is.na(after)) | as.character(before) == as.character(after)
}

protected_audit <- bind_rows(lapply(pairs, function(pair) {
  before <- pair$before
  after <- pair$after
  protected <- unique(c(
    protected_base,
    grep("^(FDR_|n_tests_FDR)", intersect(names(before), names(after)), value = TRUE)
  ))
  protected <- intersect(protected, intersect(names(before), names(after)))
  joined <- before |>
    select(all_of(c(key_columns, protected))) |>
    rename_with(~ paste0(.x, "__before"), all_of(protected)) |>
    inner_join(
      after |>
        select(all_of(c(key_columns, protected))) |>
        rename_with(~ paste0(.x, "__after"), all_of(protected)),
      by = key_columns, relationship = "one-to-one"
    )
  bind_rows(lapply(protected, function(field) {
    before_value <- joined[[paste0(field, "__before")]]
    after_value <- joined[[paste0(field, "__after")]]
    tibble(
      joined[key_columns], source_table = pair$filename,
      protected_field = field,
      before_value = as.character(before_value),
      after_value = as.character(after_value),
      unchanged = same_value(before_value, after_value)
    )
  }))
}))

if (any(!protected_audit$unchanged)) {
  first <- protected_audit |> filter(!.data$unchanged) |> slice_head(n = 1L)
  stop("Protected Stage-05 value changed: ", first$source_table, " / ",
       first$protected_field, call. = FALSE)
}

ci_audit <- bind_rows(lapply(pairs, function(pair) {
  before <- pair$before
  after <- pair$after
  joined <- before |>
    select(all_of(key_columns), CI_low_before = "CI_low",
           CI_high_before = "CI_high", p_value_before = "p_value") |>
    inner_join(
      after |>
        select(
          all_of(key_columns), CI_low_after = "CI_low",
          CI_high_after = "CI_high", p_value_after = "p_value",
          any_of(c("CI_method", "CI_level", "CI_df_method"))
        ),
      by = key_columns, relationship = "one-to-one"
    ) |>
    mutate(source_table = pair$filename)
  for (column in c("CI_method", "CI_level", "CI_df_method")) {
    if (!column %in% names(joined)) joined[[column]] <- NA
  }
  joined |>
    mutate(
      CI_excludes_zero_before = is.finite(.data$CI_low_before) &
        is.finite(.data$CI_high_before) &
        (.data$CI_low_before > 0 | .data$CI_high_before < 0),
      CI_excludes_zero_after = is.finite(.data$CI_low_after) &
        is.finite(.data$CI_high_after) &
        (.data$CI_low_after > 0 | .data$CI_high_after < 0),
      raw_p_below_0_05 = is.finite(.data$p_value_after) & .data$p_value_after < 0.05,
      ordinary_two_sided_contrast = .data$test_type %in% c(
        "named_contrast", "conditional_interaction_followup"
      ),
      CI_p_reconcile = !.data$ordinary_two_sided_contrast |
        is.na(.data$p_value_after) |
        .data$CI_excludes_zero_after == .data$raw_p_below_0_05
    )
}))

if (any(!ci_audit$CI_p_reconcile, na.rm = TRUE)) {
  stop("Finite-df Stage-05 CI exclusion does not reconcile with raw p<0.05.",
       call. = FALSE)
}

write_mirror <- function(x, filename) {
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(x, file.path(table_dir, filename), na = "")
  write_csv(x, file.path(source_dir, filename), na = "")
}

write_mirror(ci_audit, "wgcna_ci_correction_audit.csv")
write_mirror(
  protected_audit,
  "wgcna_protected_inferential_before_after.csv"
)

summary <- bind_rows(
  tibble(
    audit_section = "stage05_all_named_and_conditional_contrasts",
    metric = c(
      "CI_excludes_zero_before", "CI_excludes_zero_after",
      "raw_p_below_0_05", "protected_fields_changed"
    ),
    value = c(
      sum(ci_audit$CI_excludes_zero_before & ci_audit$ordinary_two_sided_contrast,
          na.rm = TRUE),
      sum(ci_audit$CI_excludes_zero_after & ci_audit$ordinary_two_sided_contrast,
          na.rm = TRUE),
      sum(ci_audit$raw_p_below_0_05 & ci_audit$ordinary_two_sided_contrast,
          na.rm = TRUE),
      sum(!protected_audit$unchanged)
    ),
    interpretation = c(
      "Legacy normal-multiplier intervals.",
      "Finite-df emmeans intervals from the same contrast inference as p-values.",
      "Unadjusted raw p; no threshold or FDR change.",
      "Must equal zero."
    )
  )
)

integration_dir <- path_results(
  "tables", "10_biological_integration", "gsea_wgcna_concordance", "global"
)
diagnostic_dir <- table_dir
old_long_file <- file.path(
  baseline_dir, "integration__gsea_wgcna_concordance_long.csv"
)
new_long_file <- file.path(integration_dir, "gsea_wgcna_concordance_long.csv")
old_overlap_file <- file.path(
  baseline_dir, "integration__program_specific_leading_edge_module_overlap.csv"
)
new_overlap_file <- file.path(
  integration_dir, "program_specific_leading_edge_module_overlap.csv"
)
if (all(file.exists(c(
  old_long_file, new_long_file, old_overlap_file, new_overlap_file
)))) {
  old_long <- read_csv(old_long_file, show_col_types = FALSE, progress = FALSE)
  new_long <- read_csv(new_long_file, show_col_types = FALSE, progress = FALSE)
  old_overlap <- read_csv(
    old_overlap_file, show_col_types = FALSE, progress = FALSE
  )
  new_overlap <- read_csv(
    new_overlap_file, show_col_types = FALSE, progress = FALSE
  )
  classes <- c(
    "FDR_supported_concordance", "concordant_imprecise",
    "weak_or_near_zero_module_support", "discordant", "animal_sensitive",
    "unresolved"
  )
  class_summary <- bind_rows(lapply(classes, function(class_name) {
    tibble(
      audit_section = "strict_concordance_class",
      metric = paste0(class_name, c("__before", "__after")),
      value = c(
        sum(old_long$concordance_class == class_name, na.rm = TRUE),
        sum(new_long$concordance_class == class_name, na.rm = TRUE)
      ),
      interpretation = c(
        "Legacy broad-program identity and legacy CI reporting.",
        "Ontology-aware theme identity and finite-df CI reporting."
      )
    )
  }))
  overlap_summary <- tibble(
    audit_section = "structural_theme_module_overlap",
    metric = c(
      "tests_before", "FDR_le_0_05_before", "tests_after",
      "FDR_le_0_05_after", "BH_families_after"
    ),
    value = c(
      nrow(old_overlap), sum(old_overlap$overlap_FDR <= 0.05, na.rm = TRUE),
      nrow(new_overlap), sum(new_overlap$overlap_FDR <= 0.05, na.rm = TRUE),
      n_distinct(new_overlap$overlap_family_id)
    ),
    interpretation = paste(
      "Structural convergence only; BH families remain dataset x formal",
      "contrast x spatial unit and do not alter animal-level inference."
    )
  )
  summary <- bind_rows(summary, class_summary, overlap_summary)
}

theme_file <- file.path(
  integration_dir, "ontology_aware_gsea_theme_assignments_all_contrasts.csv"
)
if (file.exists(theme_file)) {
  themes <- read_csv(theme_file, show_col_types = FALSE, progress = FALSE) |>
    mutate(source_gsea_row_id = sub("\\|.*$", "", .data$theme_assignment_id))
  legacy_change <- themes |>
    group_by(.data$source_gsea_row_id, .data$legacy_program_class) |>
    summarise(
      ontology_themes = paste(
        sort(unique(stats::na.omit(.data$theme_id))), collapse = ";"
      ),
      ontology_mapped = any(!is.na(.data$theme_id)),
      .groups = "drop"
    )
  expected <- c(
    Mitochondria_OXPHOS_Metabolism =
      "mitochondrial_respiration_oxphos",
    Synapse_Vesicle_Organization = "synaptic_signaling_vesicle",
    RNA_RNP_processing = "rna_processing_splicing_rnp",
    Ribosome_Translation = "ribosome_translation",
    Autophagy_Lysosome = "autophagy_lysosome_endosome"
  )
  legacy_change$expected_authoritative_theme <- unname(
    expected[legacy_change$legacy_program_class]
  )
  retained <- !is.na(legacy_change$expected_authoritative_theme) &
    mapply(
      grepl, legacy_change$expected_authoritative_theme,
      legacy_change$ontology_themes, MoreArgs = list(fixed = TRUE)
    )
  legacy_change$assignment_change_status <- ifelse(
    legacy_change$legacy_program_class == "Other" &
      !legacy_change$ontology_mapped,
    "remains_unmapped",
    ifelse(
      legacy_change$legacy_program_class == "Other" &
        legacy_change$ontology_mapped,
      "legacy_Other_now_registry_mapped",
      ifelse(
        retained, "broad_identity_retained_by_ontology",
        "legacy_broad_assignment_not_retained"
      )
    )
  )
  change_audit <- legacy_change |>
    count(
      .data$legacy_program_class, .data$expected_authoritative_theme,
      .data$ontology_themes, .data$ontology_mapped,
      .data$assignment_change_status,
      name = "n_source_GSEA_occurrences"
    ) |>
    arrange(
      .data$legacy_program_class, .data$assignment_change_status,
      .data$ontology_themes
    ) |>
    mutate(
      change_rule = paste(
        "A legacy non-Other identity is retained only when its prespecified",
        "ontology theme is present. Legacy Other becomes changed only when",
        "the registry supplies an authoritative theme. No text regex is used."
      )
    )
  write_mirror(
    change_audit, "legacy_program_theme_assignment_change_audit.csv"
  )
  status_counts <- legacy_change |>
    count(.data$assignment_change_status)
  summary <- bind_rows(
    summary,
    status_counts |>
      transmute(
        audit_section = "ontology_identity_change",
        metric = .data$assignment_change_status,
        value = .data$n,
        interpretation = paste(
          "Source GSEA occurrence grain; identity audit only; NES/p/FDR",
          "remain unchanged."
        )
      )
  )
}
write_mirror(summary, "gsea_wgcna_correction_summary.csv")

cat("Correction audit PASS:", nrow(ci_audit), "Stage-05 rows and",
    nrow(protected_audit), "protected field comparisons.\n")
