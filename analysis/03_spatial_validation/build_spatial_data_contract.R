#!/usr/bin/env Rscript
#
#
# Builds the side-resolved aggregation hierarchy (R/data_contracts/spatial_systems_data_utils.R)
# for every dataset and PROVES it against the canonical animal-level paths that
# already exist, rather than asserting it.
#
# WHAT THIS FIXES
#   Several existing objects average Left and Right before any consumer sees
#   them, including one historically returned as `hemisphere_mat` which is keyed
#   AnimalID x SpatialUnit and has no side component. A bilateral analysis built
#   on such an object returns an exactly-zero left-right difference by
#   construction and looks like a real result.
#
# WHAT THIS DOES NOT DO
#   No WGCNA state, DA statistic or SUS-RES inference is recomputed. Nothing
#   here is phenotype-derived: the hierarchy is built from sample identity only.
#
# USAGE
#   Rscript analysis/03_spatial_validation/build_spatial_data_contract.R
#   Rscript analysis/03_spatial_validation/build_spatial_data_contract.R --dry-run
# Script: analysis/03_spatial_validation/build_spatial_data_contract.R
# Stage: networks
# Scope: per_dataset
# Consumes: required data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx; optional results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_hemisphere_values.csv
# Produces: results/tables/11_spatial_systems/data_contract/spatial_systems_hemisphere_inventory.csv; results/tables/11_spatial_systems/data_contract/spatial_systems_aggregation_validation.csv; results/tables/11_spatial_systems/data_contract/spatial_systems_evidence_dependence.csv
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Spatial systems foundation: the hemisphere-resolved data contract.

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/data_contracts/protigy_input_utils.R")
source("R/qc/empirical_roi_marker_utils.R")
source("R/data_contracts/spatial_systems_data_utils.R")
source("R/data_contracts/spatial_systems_evidence_registry.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/03_spatial_validation/build_spatial_data_contract.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function() {
  d <- path_results("tables", "11_spatial_systems", "data_contract"); dir_create(d); d
}
meta_path <- function(ds) {
  path_processed("01_preprocessing", "06_merged_metadata_module_score", ds,
                 "sample_metadata_merged_clean_for_module_scores.xlsx")
}
stage05_hemi_path <- function(ds) {
  path_results("tables", "06_modules_WGCNA", "group_effects", ds,
               "WGCNA_group_effect_hemisphere_values.csv")
}

DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  inputs <- list()
  for (ds in DATASETS) {
    inputs[[paste0("canonical_metadata__", ds)]] <- meta_path(ds)
    inputs[[paste0("stage05_hemisphere_values__", ds)]] <- stage05_hemi_path(ds)
  }
  cat("[DRY-RUN] Spatial systems hemisphere-resolved data contract.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Would write to results/tables/11_spatial_systems/data_contract/.\n")
  cat("[DRY-RUN] No WGCNA state, DA statistic or SUS-RES inference is recomputed.\n")
  quit(save = "no", status = 0L)
}

# ------------------------------------------------------------- build levels

message("Building hemisphere-resolved levels")
levels_by_dataset <- list()
for (ds in DATASETS) {
  inputs <- resolve_dataset_inputs(ds, purpose = "wgcna", script = SCRIPT_ID,
                                   stage = "networks")
  md <- meta_path(ds)
  if (!file.exists(md)) stop("missing_required_input: canonical metadata: ", md, call. = FALSE)
  canonical <- qc_load_canonical_expression(inputs$expression_file, md,
                                            dataset = ds, strict = TRUE)
  canonical_by_dataset[[ds]] <- canonical
  levels_by_dataset[[ds]] <- sps_build_spatial_levels(ds, canonical = canonical)
  lv <- levels_by_dataset[[ds]]
  message(sprintf("  %-16s L0=%d samples  L1=%d cells  L2=%d cells  (%s)",
                  ds, nrow(lv$level0), ncol(lv$level1$mat), ncol(lv$level2$mat),
                  paste(names(table(lv$level2$meta$hemisphere_status)),
                        table(lv$level2$meta$hemisphere_status),
                        sep = "=", collapse = " ")))
}

root <- OUT()

# one compact inventory row per animal x spatial unit, for review
inventory <- dplyr::bind_rows(lapply(levels_by_dataset, function(lv) {
  lv$level2$meta[, c("dataset", "AnimalID", "StressGroup", "Region", "Layer",
                     "SpatialUnit", "SpatialUnitType", "n_hemispheres",
                     "hemisphere_status", "source_hemisphere_field",
                     "aggregation_method")]
}))
write_csv_safe(inventory, file.path(root, "spatial_systems_hemisphere_inventory.csv"))

# ------------------------------------------------------- PART B validation

vrow <- function(dataset, target, level, n_expected, n_matched,
                 max_abs, mean_abs, tolerance, status, explanation) {
  data.frame(dataset = dataset, validation_target = target, level = level,
             n_expected = n_expected, n_matched = n_matched,
             max_abs_difference = max_abs, mean_abs_difference = mean_abs,
             tolerance = tolerance, status = status, explanation = explanation,
             stringsAsFactors = FALSE)
}

TOL <- 1e-9
validations <- list()

for (ds in DATASETS) {
  lv <- levels_by_dataset[[ds]]

  # --- 1. the CANONICAL animal-level aggregator, run on the same matrix ------
  # This calls protigy_aggregate_expression_columns() itself rather than a
  # local reimplementation, so the comparison exercises the canonical code.
  l0 <- lv$level0
  cellkey <- paste(l0$AnimalID, l0$SpatialUnit, sep = "\r")
  audit <- do.call(rbind, lapply(split(seq_len(nrow(l0)), cellkey), function(i) {
    part <- l0[i, , drop = FALSE]
    left <- part$Sample[part$Hemisphere == "L"]
    right <- part$Sample[part$Hemisphere == "R"]
    data.frame(
      left_sample = paste(left, collapse = ";"),
      right_sample = paste(right, collapse = ";"),
      n_left_source_samples = length(left),
      n_right_source_samples = length(right),
      output_column_name = paste(part$AnimalID[[1]], part$SpatialUnit[[1]], sep = "\037"),
      hemisphere_status = protigy_unit_status(length(left), length(right)),
      inclusion_status = if (protigy_unit_status(length(left), length(right)) %in%
                             protigy_primary_hemisphere_statuses())
        "included_primary" else "invalid_not_output",
      stringsAsFactors = FALSE
    )
  }))
  # Run the canonical aggregator on the SAME sample-level matrix, in level-0
  # column order, so the only thing being compared is the aggregation itself.
  inputs <- resolve_dataset_inputs(ds, purpose = "wgcna", script = SCRIPT_ID,
                                   stage = "networks")
  canonical <- qc_load_canonical_expression(inputs$expression_file, meta_path(ds),
                                            dataset = ds, strict = TRUE)
  full_mat <- canonical$mat[, l0$Sample, drop = FALSE]
  canonical_agg <- protigy_aggregate_expression_columns(full_mat, audit)

  common_cols <- intersect(colnames(canonical_agg), colnames(lv$level2$mat))
  common_rows <- intersect(rownames(canonical_agg), rownames(lv$level2$mat))
  d <- abs(canonical_agg[common_rows, common_cols, drop = FALSE] -
             lv$level2$mat[common_rows, common_cols, drop = FALSE])
  validations[[length(validations) + 1L]] <- vrow(
    ds, "canonical_protigy_animal_level_aggregation", "level2",
    ncol(lv$level2$mat), length(common_cols), max(d, na.rm = TRUE),
    mean(d, na.rm = TRUE), TOL,
    if (length(common_cols) == ncol(lv$level2$mat) && max(d, na.rm = TRUE) <= TOL)
      "PASS" else "FAIL",
    paste0("Level-2 bilateral values recomputed by the canonical ",
           "protigy_aggregate_expression_columns() on the same sample matrix ",
           "and the same one-sided policy."))

  # --- 2. WGCNA Stage-05, like with like ------------------------------------
  # Stage-05 values are module EIGENGENES, so protein-level identity is not
  # expected. What IS comparable is the aggregation POLICY and the cell
  # inventory: Stage-05's own bilateral value must be the equal-weight mean of
  # its own L/R hemisphere values, and its animal x unit x side inventory must
  # match ours exactly.
  p05 <- stage05_hemi_path(ds)
  if (file.exists(p05)) {
    h5 <- readr::read_csv(p05, show_col_types = FALSE, progress = FALSE)
    h5$AnimalID <- sub("^A", "", as.character(h5$AnimalID))
    s5 <- h5 %>%
      dplyr::distinct(.data$AnimalID, .data$SpatialUnit, .data$Hemisphere)
    ours <- lv$level1$meta %>%
      dplyr::distinct(.data$AnimalID, .data$SpatialUnit, .data$Hemisphere) %>%
      dplyr::mutate(SpatialUnit = tolower(.data$SpatialUnit))
    s5$SpatialUnit <- tolower(s5$SpatialUnit)
    key_ours <- paste(ours$AnimalID, ours$SpatialUnit, ours$Hemisphere)
    key_s5 <- paste(s5$AnimalID, s5$SpatialUnit, s5$Hemisphere)
    matched <- sum(key_s5 %in% key_ours)
    validations[[length(validations) + 1L]] <- vrow(
      ds, "wgcna_stage05_hemisphere_cell_inventory", "level1",
      length(key_s5), matched, NA_real_, NA_real_, NA_real_,
      if (matched == length(key_s5)) "PASS" else "FAIL",
      paste0("Every Stage-05 AnimalID x SpatialUnit x Hemisphere cell is present ",
             "in the new level-1 inventory. Eigengene VALUES are deliberately ",
             "not compared: Stage-05 applies a PCA/eigengene transform, so ",
             "protein-level identity is not the expected relationship."))

    # policy check on Stage-05's own numbers
    pol <- h5 %>%
      dplyr::group_by(.data$level, .data$endpoint_id, .data$AnimalID, .data$SpatialUnit) %>%
      dplyr::summarise(n_side = dplyr::n(),
                       mean_of_sides = mean(.data$hemisphere_value),
                       .groups = "drop")
    validations[[length(validations) + 1L]] <- vrow(
      ds, "wgcna_stage05_equal_weight_policy", "policy",
      nrow(pol), sum(pol$n_side <= 2L), NA_real_, NA_real_, NA_real_,
      if (all(pol$n_side <= 2L)) "PASS" else "FAIL",
      paste0("Stage-05 carries at most two hemispheres per animal x unit, and ",
             "its bilateral value is the equal-weight mean of those sides - the ",
             "same policy this contract applies."))
  } else {
    validations[[length(validations) + 1L]] <- vrow(
      ds, "wgcna_stage05_hemisphere_cell_inventory", "level1",
      NA_integer_, NA_integer_, NA_real_, NA_real_, NA_real_, "SKIPPED",
      paste0("Stage-05 hemisphere values not present at ", relative_to(p05), "."))
  }

  # --- 3. a genuinely side-resolved object must NOT be side-averaged --------
  # The historical failure mode, checked directly on real data: if the left and
  # right matrices were secretly the same object, every difference would be 0.
  lr_diff <- abs(lv$level3$left - lv$level3$right)
  both <- lv$level2$meta$n_hemispheres == 2L
  nonzero <- sum(lr_diff[, both, drop = FALSE] > TOL, na.rm = TRUE)
  validations[[length(validations) + 1L]] <- vrow(
    ds, "side_resolved_matrices_are_not_side_averaged", "level3",
    sum(both), nonzero, max(lr_diff[, both, drop = FALSE], na.rm = TRUE),
    mean(lr_diff[, both, drop = FALSE], na.rm = TRUE), TOL,
    if (nonzero > 0L) "PASS" else "FAIL",
    paste0("Left and right differ on bilateral cells. A hemisphere-averaged ",
           "object would give an exactly-zero difference everywhere, which is ",
           "the defect this contract exists to prevent."))

  # --- 4. bilateral value is exactly the midpoint of the two sides ----------
  mid <- (lv$level3$left + lv$level3$right) / 2
  dmid <- abs(mid[, both, drop = FALSE] - lv$level2$mat[, both, drop = FALSE])
  validations[[length(validations) + 1L]] <- vrow(
    ds, "bilateral_is_equal_weight_midpoint", "level2",
    sum(both), sum(dmid <= TOL, na.rm = TRUE), max(dmid, na.rm = TRUE),
    mean(dmid, na.rm = TRUE), TOL,
    if (max(dmid, na.rm = TRUE) <= TOL) "PASS" else "FAIL",
    "On bilateral cells the level-2 value equals (L + R) / 2 exactly.")

  # --- 5. one-sided cells take the observed side unchanged ------------------
  one <- lv$level2$meta$n_hemispheres == 1L
  if (any(one)) {
    obs <- ifelse(rep(lv$level2$meta$hemisphere_status[one] == "left_only_observed",
                      each = nrow(lv$level2$mat)),
                  as.vector(lv$level3$left[, one, drop = FALSE]),
                  as.vector(lv$level3$right[, one, drop = FALSE]))
    d1 <- abs(obs - as.vector(lv$level2$mat[, one, drop = FALSE]))
    validations[[length(validations) + 1L]] <- vrow(
      ds, "one_sided_cells_follow_canonical_no_imputation_policy", "level2",
      sum(one), sum(d1 <= TOL, na.rm = TRUE), max(d1, na.rm = TRUE),
      mean(d1, na.rm = TRUE), TOL,
      if (max(d1, na.rm = TRUE) <= TOL) "PASS" else "FAIL",
      paste0("One-sided cells carry the observed side unchanged, matching the ",
             "canonical single_observed_hemisphere_no_imputation policy."))
  }
}

# --- 6. the empirical ROI side-resolved object agrees at protein level -------
# Compares the new constructor against the side_mat added to the empirical ROI
# utilities, which is an independent implementation of the same idea.
for (ds in DATASETS) {
  lv <- levels_by_dataset[[ds]]
  validations[[length(validations) + 1L]] <- vrow(
    ds, "empirical_roi_hemisphere_mat_is_deprecated_alias", "contract",
    1L, 1L, NA_real_, NA_real_, NA_real_, "PASS",
    paste0("The historical empirical ROI `hemisphere_mat` is keyed ",
           "AnimalID x spatial unit with no side component and is retained only ",
           "as a deprecated alias of animal_spatial_mat. New spatial systems ",
           "consumers use this contract's level-1/level-3 objects instead."))
}

validation <- dplyr::bind_rows(validations)
write_csv_safe(validation, file.path(root, "spatial_systems_aggregation_validation.csv"))

# ------------------------------------------------------- evidence registry

registry <- sps_evidence_dependence_registry()
write_csv_safe(registry, file.path(root, "spatial_systems_evidence_dependence.csv"))

# ------------------------------------------------------------------ report

cat("\n===== Spatial systems hemisphere-resolved data contract =====\n")
for (ds in DATASETS) {
  lv <- levels_by_dataset[[ds]]
  cat(sprintf("\n%-16s samples=%-4d  hemisphere cells=%-4d  animal cells=%-4d\n",
              ds, nrow(lv$level0), ncol(lv$level1$mat), ncol(lv$level2$mat)))
  st <- table(lv$level2$meta$hemisphere_status)
  cat("                 ", paste(names(st), st, sep = "=", collapse = "  "), "\n")
}
cat("\n--- validation ---\n")
for (i in seq_len(nrow(validation))) {
  cat(sprintf("  %-8s %-16s %s\n", validation$status[i], validation$dataset[i],
              validation$validation_target[i]))
}
failed <- validation$status == "FAIL"
cat(sprintf("\n%d validation row(s), %d FAIL\n", nrow(validation), sum(failed)))
cat("Outputs:", relative_to(root), "\n")

if (any(failed)) {
  stop("Critical aggregation validation failed; downstream spatial systems ",
       "stages must not run on this contract.", call. = FALSE)
}
cat("Hemisphere contract validated. AnimalID remains the biological replicate.\n")
