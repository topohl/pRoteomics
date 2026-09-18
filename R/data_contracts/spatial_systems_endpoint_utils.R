# Endpoint assembly and repeated-measures variance decomposition for the
# spatial systems layer.
#
# An "endpoint" here is one scalar per AnimalID x SpatialUnit x Hemisphere:
# a module eigengene, a reference marker score or an empirical compartment
# score. Keeping them in one tidy shape lets a single variance model serve all
# three classes without special-casing.

# ------------------------------------------------------------ level cache

.sps_level_cache <- new.env(parent = emptyenv())

# Build (and memoise) the hemisphere-resolved levels for a dataset.

## qc_find() resolves QC outputs normalized-first with a historical
## fallback; it lives in its own small file so this library does not have to
## load the whole QC utility stack.
if (!exists("qc_find", mode = "function")) {
  if (!exists("repo_path", mode = "function")) {
    paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
    source(paths_file)
  }
  source(repo_path("R", "qc_result_paths.R"))
}
if (!exists("preprocessing_find", mode = "function")) source(repo_path("R", "preprocessing_paths.R"))

sps_levels_for_dataset <- function(dataset) {
  key <- as.character(dataset)
  if (!is.null(.sps_level_cache[[key]])) return(.sps_level_cache[[key]])
  inputs <- resolve_dataset_inputs(dataset, purpose = "wgcna",
                                   script = Sys.getenv("PROTEOMICS_SCRIPT_ID",
                                                       "spatial_systems"),
                                   stage = "networks")
  ## Normalized first, historical second (Phase 6G.7).
  md <- preprocessing_find("sample_metadata_merged_clean_for_module_scores.xlsx",
                           owner = "build_module_score_metadata",
                           legacy_stage = "01_preprocessing",
                           legacy_substep = "06_merged_metadata_module_score",
                           scope = dataset, child = "tables",
                           legacy_family = "processed")
  canonical <- qc_load_canonical_expression(inputs$expression_file, md,
                                            dataset = dataset, strict = TRUE)
  lv <- sps_build_spatial_levels(dataset, canonical = canonical)
  lv$feature_table <- canonical$feature_table
  .sps_level_cache[[key]] <- lv
  lv
}

.sps_endpoint_frame <- function(dataset, endpoint_class, endpoint_id, meta, value) {
  data.frame(
    dataset = dataset,
    endpoint_class = endpoint_class,
    endpoint_id = endpoint_id,
    AnimalID = meta$AnimalID,
    StressGroup = meta$StressGroup,
    SpatialUnit = meta$SpatialUnit,
    Hemisphere = meta$Hemisphere,
    value = as.numeric(value),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------- endpoint builders

# 1. WGCNA modules: the ACCEPTED Stage-05 hemisphere values, not a recomputation.
sps_endpoints_wgcna_modules <- function(dataset, level_filter = "module") {
  p <- path_results("tables", "06_modules_WGCNA", "group_effects", dataset,
                    "WGCNA_group_effect_hemisphere_values.csv")
  if (!file.exists(p)) return(NULL)
  x <- as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                     guess_max = Inf))
  if (!is.null(level_filter)) x <- x[x$level %in% level_filter, , drop = FALSE]
  if (!nrow(x)) return(NULL)
  data.frame(
    dataset = dataset,
    endpoint_class = "wgcna_module_eigengene",
    endpoint_id = as.character(x$endpoint_id),
    AnimalID = sub("^A", "", as.character(x$AnimalID)),
    StressGroup = as.character(x$StressGroup),
    SpatialUnit = toupper(as.character(x$SpatialUnit)),
    Hemisphere = as.character(sps_normalize_hemisphere(x$Hemisphere)),
    value = as.numeric(x$hemisphere_value),
    stringsAsFactors = FALSE
  )
}

# Score a protein set on the hemisphere-resolved level-1 matrix.
#
# Each protein is z-scored ACROSS the level-1 columns first, so a panel score is
# not dominated by whichever member happens to be most abundant; the score is
# then the mean z across available members.
sps_score_protein_set <- function(levels, protein_ids, min_members = 3L) {
  m <- levels$level1$mat
  ids <- intersect(unique(as.character(protein_ids)), rownames(m))
  if (length(ids) < min_members) return(NULL)
  sub <- m[ids, , drop = FALSE]
  mu <- rowMeans(sub, na.rm = TRUE)
  sdv <- apply(sub, 1, stats::sd, na.rm = TRUE)
  keep <- is.finite(sdv) & sdv > 0
  if (sum(keep) < min_members) return(NULL)
  z <- sweep(sweep(sub[keep, , drop = FALSE], 1, mu[keep], "-"), 1, sdv[keep], "/")
  list(score = colMeans(z, na.rm = TRUE), n_members = sum(keep))
}

# 2. Canonical reference marker panels (external/curated), mapped to this
#    dataset's protein groups through the canonical feature table.
sps_endpoints_reference_marker_scores <- function(dataset, levels,
                                                  min_members = 3L) {
  p <- repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")
  if (!file.exists(p)) return(NULL)
  panels <- as.data.frame(readr::read_csv(p, show_col_types = FALSE,
                                          progress = FALSE))
  if (!all(c("marker_set", "gene_symbol") %in% names(panels))) return(NULL)
  ft <- levels$feature_table
  sym_col <- intersect(c("official_gene_symbol", "gene_symbol",
                         "representative_gene_symbol"), names(ft))
  if (!length(sym_col)) return(NULL)
  lut <- data.frame(
    sym = toupper(trimws(as.character(ft[[sym_col[[1]]]]))),
    pg = as.character(ft$ProteinGroupID), stringsAsFactors = FALSE)
  lut <- lut[nzchar(lut$sym) & !is.na(lut$sym), , drop = FALSE]

  out <- list()
  for (set in sort(unique(panels$marker_set))) {
    syms <- toupper(trimws(panels$gene_symbol[panels$marker_set == set]))
    pg <- unique(lut$pg[lut$sym %in% syms])
    sc <- sps_score_protein_set(levels, pg, min_members = min_members)
    if (is.null(sc)) next
    out[[length(out) + 1L]] <- .sps_endpoint_frame(
      dataset, "reference_marker_score", set, levels$level1$meta, sc$score)
  }
  if (!length(out)) return(NULL)
  dplyr::bind_rows(out)
}

# 3. Empirical (experiment-derived) compartment marker sets.
sps_endpoints_empirical_compartment_scores <- function(dataset, levels,
                                                       min_members = 3L) {
  p <- qc_find("empirical_roi_marker_sets.csv",
               owner = "discover_empirical_roi_markers",
               legacy_substep = "05_empirical_roi_marker_discovery")
  if (!file.exists(p)) return(NULL)
  sets <- as.data.frame(readr::read_csv(p, show_col_types = FALSE,
                                        progress = FALSE, guess_max = Inf))
  if (!"marker_set" %in% names(sets)) return(NULL)
  # prefer the dataset-specific protein group id, fall back to the shared key
  id_col <- intersect(c(paste0("ProteinGroupID_", dataset), "ProteinGroupID"),
                      names(sets))
  if (!length(id_col)) return(NULL)
  out <- list()
  for (set in sort(unique(sets$marker_set))) {
    pg <- unique(unlist(lapply(id_col, function(cc)
      as.character(sets[[cc]][sets$marker_set == set]))))
    pg <- pg[!is.na(pg) & nzchar(pg)]
    sc <- sps_score_protein_set(levels, pg, min_members = min_members)
    if (is.null(sc)) next
    out[[length(out) + 1L]] <- .sps_endpoint_frame(
      dataset, "empirical_compartment_score", set, levels$level1$meta, sc$score)
  }
  if (!length(out)) return(NULL)
  dplyr::bind_rows(out)
}

# ------------------------------------------------- variance decomposition

# Fit the repeated-measures model for ONE endpoint and return its components.
#
# Every guard the reliability formula depends on is recorded on the row, so a
# reader can see exactly why a number is or is not reported.
sps_variance_components <- function(e) {
  base <- data.frame(
    dataset = e$dataset[[1]], endpoint_class = e$endpoint_class[[1]],
    endpoint_id = e$endpoint_id[[1]],
    n_observations = nrow(e),
    n_animals = length(unique(e$AnimalID)),
    n_spatial_units = length(unique(e$SpatialUnit)),
    n_hemispheres = length(unique(e$Hemisphere)),
    formula_used = "value ~ Hemisphere + (1|AnimalID) + (1|SpatialUnit) + (1|AnimalID:SpatialUnit)",
    stringsAsFactors = FALSE)

  insufficient <- function(reason) {
    cbind(base, data.frame(
      between_animal_variance = NA_real_, spatial_unit_variance = NA_real_,
      animal_spatial_cell_variance = NA_real_,
      within_animal_hemispheric_variance = NA_real_, total_variance = NA_real_,
      hemisphere_fixed_effect = NA_real_, hemisphere_fixed_effect_se = NA_real_,
      hemisphere_fixed_effect_p = NA_real_,
      model_convergence = "not_fitted", is_singular = NA,
      assumption_status = reason, stringsAsFactors = FALSE))
  }
  if (base$n_animals < 3L) return(insufficient("too_few_animals"))
  if (base$n_hemispheres < 2L) return(insufficient("single_hemisphere_only"))
  if (stats::sd(e$value, na.rm = TRUE) == 0) return(insufficient("zero_variance_endpoint"))

  e$AnimalID <- factor(e$AnimalID)
  e$SpatialUnit <- factor(e$SpatialUnit)
  e$Hemisphere <- factor(e$Hemisphere, levels = c("L", "R"))

  msgs <- character()
  fit <- withCallingHandlers(
    tryCatch(
      lme4::lmer(value ~ Hemisphere + (1 | AnimalID) + (1 | SpatialUnit) +
                   (1 | AnimalID:SpatialUnit), data = e, REML = TRUE),
      error = function(err) NULL),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(fit)) return(insufficient("model_failed"))

  vc <- as.data.frame(lme4::VarCorr(fit))
  getv <- function(grp) {
    v <- vc$vcov[vc$grp == grp]
    if (!length(v)) NA_real_ else v[[1]]
  }
  s2_animal <- getv("AnimalID")
  s2_unit <- getv("SpatialUnit")
  s2_cell <- getv("AnimalID:SpatialUnit")
  s2_res <- getv("Residual")

  fe <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  hemi_row <- if (!is.null(fe) && "HemisphereR" %in% rownames(fe)) fe["HemisphereR", ] else NULL
  hemi_est <- if (!is.null(hemi_row)) unname(hemi_row[["Estimate"]]) else NA_real_
  hemi_se <- if (!is.null(hemi_row)) unname(hemi_row[["Std. Error"]]) else NA_real_
  # Wald p: adequate here because Hemisphere is a single balanced fixed contrast
  hemi_p <- if (is.finite(hemi_est) && is.finite(hemi_se) && hemi_se > 0) {
    2 * stats::pnorm(-abs(hemi_est / hemi_se))
  } else NA_real_

  singular <- isTRUE(lme4::isSingular(fit, tol = 1e-4))
  converged <- !length(grep("converge", msgs, ignore.case = TRUE))
  status <- if (singular) "singular_fit_reliability_not_reported"
  else if (!converged) "did_not_converge"
  else if (!is.finite(s2_animal) || !is.finite(s2_res)) "variance_not_estimable"
  else if (s2_animal + s2_res <= 0) "degenerate_variance"
  else "assumptions_met"

  cbind(base, data.frame(
    between_animal_variance = s2_animal,
    spatial_unit_variance = s2_unit,
    animal_spatial_cell_variance = s2_cell,
    within_animal_hemispheric_variance = s2_res,
    total_variance = sum(c(s2_animal, s2_unit, s2_cell, s2_res), na.rm = TRUE),
    hemisphere_fixed_effect = hemi_est,
    hemisphere_fixed_effect_se = hemi_se,
    hemisphere_fixed_effect_p = hemi_p,
    model_convergence = if (converged) "converged" else "convergence_warning",
    is_singular = singular,
    assumption_status = status,
    stringsAsFactors = FALSE))
}

# Reliability and expected precision gain, reported ONLY where the guards pass.
sps_precision_gain <- function(vc) {
  ok <- vc$assumption_status == "assumptions_met"
  a <- vc$between_animal_variance
  h <- vc$within_animal_hemispheric_variance
  icc1 <- ifelse(ok, a / (a + h), NA_real_)
  icc2 <- ifelse(ok, a / (a + h / 2), NA_real_)
  out <- vc[, c("dataset", "endpoint_class", "endpoint_id", "n_animals",
                "n_spatial_units", "between_animal_variance",
                "within_animal_hemispheric_variance",
                "hemisphere_fixed_effect", "hemisphere_fixed_effect_p",
                "model_convergence", "is_singular", "assumption_status",
                "formula_used")]
  out$ICC_single_side <- icc1
  out$ICC_bilateral_mean <- icc2
  # averaging two independent sides halves the hemispheric variance component
  out$relative_measurement_variance_reduction <- ifelse(ok, 0.5, NA_real_)
  out$expected_reliability_gain <- icc2 - icc1
  out$reliability_formula <- ifelse(
    ok,
    "ICC1 = s2_animal/(s2_animal + s2_hemi); ICC2 = s2_animal/(s2_animal + s2_hemi/2)",
    NA_character_)
  out$reliability_assumption <- ifelse(
    ok,
    paste0("independent side noise, spatial structure modelled as a random ",
           "effect, non-singular converged fit"),
    NA_character_)
  out
}
