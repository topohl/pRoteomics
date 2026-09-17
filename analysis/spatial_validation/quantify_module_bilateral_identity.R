#!/usr/bin/env Rscript
#
#
# Consumes the EXISTING Stage-05 hemisphere-level module values. No WGCNA state,
# membership or eigengene is recomputed here.
#
# TWO DIFFERENT QUESTIONS, KEPT APART
#   1. ABSOLUTE VALUE reproducibility - across matched AnimalID x SpatialUnit
#      cells, does the left value equal the right value?
#   2. SPATIAL PATTERN reproducibility - within one animal, does the module's
#      profile ACROSS spatial units have the same shape on both sides?
#
#   These can disagree, and the distinction matters: a module with a constant
#   side offset scores poorly on (1) while being perfectly preserved on (2).
#   Collapsing them into one number would hide exactly that case.
#
# A LOW CORRELATION IS NOT AUTOMATICALLY TECHNICAL FAILURE. It may be real
# hemispheric asymmetry. Nothing here is labelled as a quality verdict.
#
# USAGE
#   Rscript analysis/spatial_validation/quantify_module_bilateral_identity.R
#   Rscript analysis/spatial_validation/quantify_module_bilateral_identity.R --dry-run
# Script: analysis/spatial_validation/quantify_module_bilateral_identity.R
# Stage: networks
# Scope: per_dataset
# Consumes: required results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_hemisphere_values.csv; optional none declared in pipeline.yml
# Produces: results/tables/11_spatial_systems/bilateral/WGCNA_module_bilateral_reproducibility.csv; results/tables/11_spatial_systems/bilateral/WGCNA_module_bilateral_spatial_profile_reproducibility.csv; results/tables/11_spatial_systems/bilateral/WGCNA_module_bilateral_spatial_profile_summary.csv
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: WGCNA module bilateral reproducibility.

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/data_contracts/spatial_systems_data_utils.R")
source("R/spatial/spatial_systems_bilateral_utils.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/spatial_validation/quantify_module_bilateral_identity.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function() {
  d <- path_results("tables", "11_spatial_systems", "bilateral"); dir_create(d); d
}
stage05_hemi_path <- function(ds) {
  path_results("tables", "06_modules_WGCNA", "group_effects", ds,
               "WGCNA_group_effect_hemisphere_values.csv")
}
DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  inputs <- stats::setNames(lapply(DATASETS, stage05_hemi_path),
                            paste0("stage05_hemisphere_values__", DATASETS))
  cat("[DRY-RUN] WGCNA module bilateral reproducibility.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Consumes Stage-05 hemisphere values; recomputes no WGCNA state.\n")
  quit(save = "no", status = 0L)
}

read_stage05 <- function(ds) {
  p <- stage05_hemi_path(ds)
  if (!file.exists(p)) stop("missing_required_input: Stage-05 hemisphere values: ", p,
                            call. = FALSE)
  x <- readr::read_csv(p, show_col_types = FALSE, progress = FALSE, guess_max = Inf)
  need <- c("dataset", "level", "endpoint_id", "AnimalID", "StressGroup",
            "SpatialUnit", "Hemisphere", "hemisphere_value")
  missing <- setdiff(need, names(x))
  if (length(missing)) {
    stop("Stage-05 hemisphere values missing column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  x$Hemisphere <- as.character(sps_normalize_hemisphere(x$Hemisphere))
  as.data.frame(x)
}

absolute_rows <- list()
profile_rows <- list()
profile_summary_rows <- list()

for (ds in DATASETS) {
  h <- read_stage05(ds)
  message(sprintf("%-16s %d rows, %d endpoints, %d spatial units",
                  ds, nrow(h), length(unique(h$endpoint_id)),
                  length(unique(h$SpatialUnit))))

  for (ep in sort(unique(h$endpoint_id))) {
    e <- h[h$endpoint_id == ep, , drop = FALSE]
    lvl <- e$level[[1]]

    # ---- 1. ABSOLUTE VALUE reproducibility -------------------------------
    # Pair strictly on AnimalID + SpatialUnit: a left value may only ever be
    # compared with the right value of the SAME animal and the SAME unit.
    w <- sps_pair_sides(e, id_cols = c("AnimalID", "SpatialUnit"),
                        side_col = "Hemisphere", value_col = "hemisphere_value")
    absolute_rows[[length(absolute_rows) + 1L]] <- cbind(
      data.frame(dataset = ds, level = lvl, endpoint_id = ep,
                 stringsAsFactors = FALSE),
      sps_paired_agreement(w$left, w$right),
      data.frame(n_animals = length(unique(w$AnimalID)),
                 n_spatial_units = length(unique(w$SpatialUnit)),
                 stringsAsFactors = FALSE))

    # ---- 2. SPATIAL PATTERN reproducibility ------------------------------
    # Within one animal, correlate the across-unit profile on L with the
    # across-unit profile on R. Units are matched exactly and in the same order.
    for (an in sort(unique(w$AnimalID))) {
      wa <- w[w$AnimalID == an, , drop = FALSE]
      wa <- wa[order(wa$SpatialUnit), , drop = FALSE]
      ok <- is.finite(wa$left) & is.finite(wa$right)
      n_units <- sum(ok)
      pr <- sps_safe_cor(wa$left[ok], wa$right[ok], "pearson")
      sr <- sps_safe_cor(wa$left[ok], wa$right[ok], "spearman")
      profile_rows[[length(profile_rows) + 1L]] <- data.frame(
        dataset = ds, level = lvl, endpoint_id = ep, AnimalID = an,
        StressGroup = e$StressGroup[match(an, e$AnimalID)],
        n_spatial_units = n_units,
        spatial_profile_Pearson = pr,
        spatial_profile_Spearman = sr,
        matched_units = paste(wa$SpatialUnit[ok], collapse = ";"),
        stringsAsFactors = FALSE)
    }
  }
}

absolute <- dplyr::bind_rows(absolute_rows)
profiles <- dplyr::bind_rows(profile_rows)

# summarise the per-animal profile correlations by module
profile_summary <- profiles %>%
  dplyr::group_by(.data$dataset, .data$level, .data$endpoint_id) %>%
  dplyr::summarise(
    n_animals = dplyr::n(),
    n_spatial_units = stats::median(.data$n_spatial_units),
    median_profile_Pearson = stats::median(.data$spatial_profile_Pearson, na.rm = TRUE),
    IQR_profile_Pearson = stats::IQR(.data$spatial_profile_Pearson, na.rm = TRUE),
    min_profile_Pearson = suppressWarnings(min(.data$spatial_profile_Pearson, na.rm = TRUE)),
    max_profile_Pearson = suppressWarnings(max(.data$spatial_profile_Pearson, na.rm = TRUE)),
    median_profile_Spearman = stats::median(.data$spatial_profile_Spearman, na.rm = TRUE),
    IQR_profile_Spearman = stats::IQR(.data$spatial_profile_Spearman, na.rm = TRUE),
    .groups = "drop")

# A descriptive class that names WHICH of the two reproducibilities holds.
# Deliberately not a quality score: poor bilateral agreement can be real
# hemispheric asymmetry rather than a measurement problem.
absolute <- absolute %>%
  dplyr::left_join(profile_summary[, c("dataset", "level", "endpoint_id",
                                       "median_profile_Pearson")],
                   by = c("dataset", "level", "endpoint_id")) %>%
  dplyr::mutate(
    bilateral_reproducibility_class = sps_reproducibility_class(
      absolute_r = .data$pearson_r,
      profile_r = .data$median_profile_Pearson,
      mean_signed_difference = .data$mean_signed_L_minus_R,
      median_abs_difference = .data$median_abs_difference))

root <- OUT()
write_csv_safe(absolute, file.path(root, "WGCNA_module_bilateral_reproducibility.csv"))
write_csv_safe(profiles, file.path(root, "WGCNA_module_bilateral_spatial_profile_reproducibility.csv"))
write_csv_safe(profile_summary, file.path(root, "WGCNA_module_bilateral_spatial_profile_summary.csv"))

cat("\n===== WGCNA module bilateral reproducibility =====\n")
for (ds in DATASETS) {
  a <- absolute[absolute$dataset == ds & absolute$level == "module", , drop = FALSE]
  if (!nrow(a)) next
  cat(sprintf("\n--- %s (%d modules) ---\n", ds, nrow(a)))
  cat(sprintf("  %-12s %8s %8s %8s %9s  %s\n",
              "module", "abs_r", "prof_r", "MAE", "signedLR", "class"))
  for (i in seq_len(nrow(a))) {
    cat(sprintf("  %-12s %8.3f %8.3f %8.4f %9.4f  %s\n",
                a$endpoint_id[i], a$pearson_r[i], a$median_profile_Pearson[i],
                a$MAE[i], a$mean_signed_L_minus_R[i],
                a$bilateral_reproducibility_class[i]))
  }
}
cat("\nClass counts:\n"); print(table(absolute$bilateral_reproducibility_class))
cat("\nOutputs:", relative_to(root), "\n")
cat("Stage-05 hemisphere values consumed; no WGCNA state recomputed.\n")
