#!/usr/bin/env Rscript
#
#
# The point is NOT to show left and right correlate. It is to estimate how much
# of the variance in an animal-level estimate is between-animal signal and how
# much is within-animal hemispheric variation, and therefore whether averaging
# the two sides measurably improves precision.
#
# THE MODEL
#   value ~ Hemisphere + (1 | AnimalID) + (1 | SpatialUnit) + (1 | AnimalID:SpatialUnit)
#
#   Hemisphere is FIXED, so a systematic side offset is estimated rather than
#   absorbed into noise. SpatialUnit is a random effect so the anatomical
#   gradient is accounted for - without it the residual would be inflated by
#   spatial structure and the reliability gain would be overstated.
#   AnimalID:SpatialUnit is the cell each L/R pair sits in, so the RESIDUAL is
#   the within-animal hemispheric variance rather than a mixture.
#
# WHY THE VARIANCE IS NOT CALLED MEASUREMENT ERROR
#   Left-right variation is partly real anatomical asymmetry. It is reported as
#   `within_animal_hemispheric_variance` throughout.
#
# THE RELIABILITY FORMULA IS GUARDED
#   ICC_single   = s2_animal / (s2_animal + s2_hemi)
#   ICC_bilateral= s2_animal / (s2_animal + s2_hemi / 2)
#   These assume independent side noise. They are reported ONLY when the model
#   converged, was non-singular, and spatial structure was in the model. Every
#   row records the formula, the assumption status and the fit diagnostics.
#
# USAGE
#   Rscript analysis/spatial_validation/decompose_bilateral_variance.R
#   Rscript analysis/spatial_validation/decompose_bilateral_variance.R --dry-run
# Script: analysis/spatial_validation/decompose_bilateral_variance.R
# Stage: networks
# Scope: per_dataset
# Consumes: required results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_hemisphere_values.csv; optional config/marker_panels/wgcna_reference_marker_sets.csv; results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv
# Produces: results/spatial_validation/decompose_bilateral_variance/global/tables/bilateral_variance_decomposition.csv; results/spatial_validation/decompose_bilateral_variance/global/tables/bilateral_precision_gain.csv
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Bilateral variance decomposition: what does measuring both hemispheres buy?
#  
#  

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/data_contracts/spatial_systems_data_utils.R")
source("R/spatial/spatial_systems_bilateral_utils.R")
source("R/data_contracts/spatial_systems_endpoint_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))

# Phase 6G.3: destinations resolve through the normalized output contract,
# addressed by this analysis's own identity rather than by the historical
# 11_spatial_systems stage directory. Outputs already written there stay
# exactly where they are and are read, never rewritten.
ANALYSIS_ID <- "decompose_bilateral_variance"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/spatial_validation/decompose_bilateral_variance.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function() {
  d <- CANONICAL_PATHS$tables; dir_create(d); d
}
DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  inputs <- list(
    reference_marker_panels = repo_path("config", "marker_panels",
                                        "wgcna_reference_marker_sets.csv"),
    empirical_marker_sets = path_results("tables", "03_qc_exploration",
                                         "05_empirical_roi_marker_discovery",
                                         "empirical_roi_marker_sets.csv"))
  for (ds in DATASETS) {
    inputs[[paste0("stage05_hemisphere_values__", ds)]] <-
      path_results("tables", "06_modules_WGCNA", "group_effects", ds,
                   "WGCNA_group_effect_hemisphere_values.csv")
  }
  cat("[DRY-RUN] Bilateral variance decomposition and precision gain.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Endpoint classes: WGCNA modules, reference marker scores, empirical compartment scores.\n")
  quit(save = "no", status = 0L)
}

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("missing_required_input: lme4 is required for the variance decomposition.",
       call. = FALSE)
}

# --------------------------------------------------------------- endpoints

endpoints <- list()
for (ds in DATASETS) {
  message("Assembling endpoints for ", ds)
  # 1. WGCNA modules, from the accepted Stage-05 hemisphere values
  endpoints[[length(endpoints) + 1L]] <- sps_endpoints_wgcna_modules(ds)
  # 2 + 3. marker and compartment scores, built on the hemisphere-resolved
  #        level-1 matrix from the data contract
  lv <- sps_levels_for_dataset(ds)
  endpoints[[length(endpoints) + 1L]] <- sps_endpoints_reference_marker_scores(ds, lv)
  endpoints[[length(endpoints) + 1L]] <- sps_endpoints_empirical_compartment_scores(ds, lv)
}
endpoints <- dplyr::bind_rows(endpoints)
endpoints <- endpoints[is.finite(endpoints$value), , drop = FALSE]
sps_assert_animal_is_replicate(endpoints, "endpoint table")

message(sprintf("Endpoints assembled: %d rows, %d endpoint(s) across %d class(es)",
                nrow(endpoints), length(unique(paste(endpoints$dataset, endpoints$endpoint_id))),
                length(unique(endpoints$endpoint_class))))

# ------------------------------------------------------------ decomposition

vc_rows <- list()
for (key in unique(paste(endpoints$dataset, endpoints$endpoint_class,
                         endpoints$endpoint_id, sep = "\037"))) {
  parts <- strsplit(key, "\037", fixed = TRUE)[[1]]
  e <- endpoints[endpoints$dataset == parts[1] &
                   endpoints$endpoint_class == parts[2] &
                   endpoints$endpoint_id == parts[3], , drop = FALSE]
  vc_rows[[length(vc_rows) + 1L]] <- sps_variance_components(e)
}
vc <- dplyr::bind_rows(vc_rows)

precision <- sps_precision_gain(vc)

root <- OUT()
write_csv_safe(vc, file.path(root, "bilateral_variance_decomposition.csv"))
write_csv_safe(precision, file.path(root, "bilateral_precision_gain.csv"))

# ------------------------------------------------------------------ report

cat("\n===== Bilateral variance decomposition =====\n")
usable <- precision[precision$assumption_status == "assumptions_met", , drop = FALSE]
cat(sprintf("\n%d endpoint(s) fitted; %d usable for a reliability statement\n",
            nrow(vc), nrow(usable)))
cat("  model convergence:  ", paste(names(table(vc$model_convergence)),
                                    table(vc$model_convergence), sep = "=",
                                    collapse = "  "), "\n")
cat("  singular fits:      ", sum(vc$is_singular %in% TRUE), "\n")

for (cl in unique(vc$endpoint_class)) {
  z <- vc[vc$endpoint_class == cl, , drop = FALSE]
  u <- usable[usable$endpoint_class == cl, , drop = FALSE]
  cat(sprintf("\n--- %s (%d endpoints, %d usable) ---\n", cl, nrow(z), nrow(u)))
  if (!nrow(u)) { cat("  no endpoint met the assumption guards\n"); next }
  cat(sprintf("  median between-animal variance      : %.5f\n",
              stats::median(u$between_animal_variance, na.rm = TRUE)))
  cat(sprintf("  median within-animal hemispheric var: %.5f\n",
              stats::median(u$within_animal_hemispheric_variance, na.rm = TRUE)))
  cat(sprintf("  median ICC single side              : %.3f\n",
              stats::median(u$ICC_single_side, na.rm = TRUE)))
  cat(sprintf("  median ICC bilateral mean           : %.3f\n",
              stats::median(u$ICC_bilateral_mean, na.rm = TRUE)))
  cat(sprintf("  median measurement-variance reduction: %.1f%%\n",
              100 * stats::median(u$relative_measurement_variance_reduction, na.rm = TRUE)))
  cat(sprintf("  endpoints with a systematic side effect (p<0.05): %d of %d\n",
              sum(u$hemisphere_fixed_effect_p < 0.05, na.rm = TRUE), nrow(u)))
}

cat("\nOutputs:", relative_to(root), "\n")
cat("Left-right variation is reported as within-animal hemispheric variance,\n")
cat("not as measurement error, and n remains the number of animals.\n")
