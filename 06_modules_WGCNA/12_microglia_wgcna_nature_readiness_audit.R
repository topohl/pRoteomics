#!/usr/bin/env Rscript

# Additive, fixed-membership robustness audit for the current microglia WGCNA
# network. This script never writes to the primary WGCNA, Stage 05, registry,
# or publication-figure paths.

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c(
  "DBI", "digest", "ggplot2", "jsonlite", "lme4", "RSQLite", "svglite",
  "variancePartition", "WGCNA"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}
repo_root <- normalizePath(arg_value("--repo-root", getwd()), winslash = "/", mustWork = TRUE)
output_dir <- file.path(
  repo_root,
  arg_value("--output-dir", "results/reviewer_audit/microglia_wgcna_nature_readiness")
)
n_permutations <- as.integer(arg_value("--permutations", "100"))
if (is.na(n_permutations) || n_permutations < 0L) stop("--permutations must be a non-negative integer")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

atomic_write <- function(path, writer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writer(tmp)
  if (!file.exists(tmp) || file.info(tmp)$size <= 0) stop("Incomplete write: ", path)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not publish completed output: ", path)
  invisible(path)
}

write_csv_atomic <- function(x, path) {
  atomic_write(path, function(tmp) utils::write.csv(x, tmp, row.names = FALSE, na = ""))
}

write_gz_csv_atomic <- function(x, path) {
  atomic_write(path, function(tmp) {
    con <- gzfile(tmp, open = "wt", compression = 6)
    on.exit(close(con), add = TRUE)
    utils::write.csv(x, con, row.names = FALSE, na = "")
  })
}

write_lines_atomic <- function(x, path) {
  atomic_write(path, function(tmp) writeLines(enc2utf8(x), tmp, useBytes = TRUE))
}

save_plot_atomic <- function(plot, path, width, height) {
  atomic_write(path, function(tmp) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "svg") {
      ggplot2::ggsave(tmp, plot = plot, device = svglite::svglite,
                      width = width, height = height, units = "in")
    } else if (ext == "pdf") {
      ggplot2::ggsave(tmp, plot = plot, device = grDevices::cairo_pdf,
                      width = width, height = height, units = "in")
    } else stop("Unsupported figure extension: ", ext)
  })
}

rel <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(repo_root, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}

state_path <- file.path(repo_root, "data/processed/06_modules_WGCNA/01_WGCNA/microglia/wgcna_final_model_state.rds")
definitions_path <- file.path(repo_root, "results/tables/06_modules_WGCNA/01_WGCNA/microglia/modules/WGCNA_module_definitions_for_downstream.csv")
super_values_path <- file.path(repo_root, "results/tables/06_modules_WGCNA/group_effects/microglia/all_supermodule_eigengene_group_values.csv")
annotation_path <- file.path(repo_root, "results/tables/06_modules_WGCNA/module_annotation/microglia/WGCNA_module_biological_annotation.csv")
registry_path <- file.path(repo_root, "config/wgcna_labels/microglia.csv")
joint_path <- file.path(repo_root, "data/processed/01_preprocessing/joint_compartment_qc/global/joint_compartment_qc_matrices.rds")
reference_markers_path <- file.path(repo_root, "config/marker_panels/wgcna_reference_marker_sets.csv")
microenvironment_markers_path <- file.path(repo_root, "config/marker_panels/microenvironment_marker_panels.csv")
empirical_markers_path <- file.path(repo_root, "results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv")

required_inputs <- c(
  state_path, definitions_path, super_values_path, annotation_path,
  registry_path, joint_path, reference_markers_path,
  microenvironment_markers_path, empirical_markers_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Required inputs missing: ", paste(rel(missing_inputs), collapse = ", "))
reference_state_md5 <- unname(tools::md5sum(state_path))
preservation_cache_contract <- "fixed_membership_sensitivity_v1"

protected_roots <- c(
  state_path,
  file.path(repo_root, "results/tables/06_modules_WGCNA/group_effects/microglia"),
  registry_path,
  file.path(repo_root, "results/figures/06_modules_WGCNA/wgcna_publication_figures")
)
protected_files <- unique(unlist(lapply(protected_roots, function(x) {
  if (!file.exists(x)) return(character())
  if (dir.exists(x)) list.files(x, recursive = TRUE, full.names = TRUE, all.files = TRUE) else x
}), use.names = FALSE))
protected_files <- protected_files[file.exists(protected_files) & !dir.exists(protected_files)]
protected_before <- data.frame(
  protected_path = vapply(protected_files, rel, character(1)),
  bytes = unname(file.info(protected_files)$size),
  md5_before = unname(tools::md5sum(protected_files)),
  stringsAsFactors = FALSE
)

state <- readRDS(state_path)
expression_primary <- as.matrix(state$expression.data)
storage.mode(expression_primary) <- "double"
if (nrow(expression_primary) != 72L || ncol(expression_primary) != 5201L) {
  stop("Current primary expression matrix is not the verified 72 x 5201 state")
}
if (anyDuplicated(colnames(expression_primary))) stop("Primary ProteinGroupID keys are duplicated")
if (anyNA(expression_primary)) stop("Primary expression matrix contains missing values")

module_long <- state$WGCNA_modules_long
module_long <- module_long[match(colnames(expression_primary), module_long$ProteinGroupID), , drop = FALSE]
if (anyNA(module_long$ModuleID) || !identical(module_long$ProteinGroupID, colnames(expression_primary))) {
  stop("Module membership is not exactly aligned to the active ProteinGroupID universe")
}
module_ids <- sprintf("WGCNA_m%02d", seq_len(13))
if (!setequal(unique(module_long$ModuleID), module_ids)) stop("Current module identities are not WGCNA_m01..WGCNA_m13")
module_long$ModuleID <- factor(module_long$ModuleID, levels = module_ids)

definitions <- utils::read.csv(definitions_path, check.names = FALSE)
module_map <- unique(definitions[c("ModuleID", "module_eigengene", "ModuleColor", "SupermoduleID")])
module_map <- module_map[match(module_ids, module_map$ModuleID), , drop = FALSE]
if (anyNA(module_map$module_eigengene) || anyDuplicated(module_map$ModuleID)) stop("Module definition identity mismatch")

expected_blocks <- list(
  SM01 = c("WGCNA_m05", "WGCNA_m10"),
  SM03 = c("WGCNA_m04", "WGCNA_m08", "WGCNA_m12"),
  SM09 = c("WGCNA_m07", "WGCNA_m11")
)
observed_members <- split(module_map$ModuleID, module_map$SupermoduleID)
for (sm in names(expected_blocks)) {
  if (!setequal(observed_members[[sm]], expected_blocks[[sm]])) stop("Higher-order membership mismatch for ", sm)
}
singleton_ids <- setdiff(sprintf("SM%02d", seq_len(9)), names(expected_blocks))
if (any(lengths(observed_members[singleton_ids]) != 1L)) stop("The six standalone supermodule identities are not singletons")

sample_names <- rownames(expression_primary)
sample_info <- state$sample_info[match(sample_names, rownames(state$sample_info)), , drop = FALSE]
extract_one <- function(pattern, x, replacement = "\\1") {
  out <- sub(paste0(".*", pattern, ".*"), replacement, x, perl = TRUE)
  bad <- identical(out, x) | !grepl(pattern, x, perl = TRUE)
  out[bad] <- NA_character_
  out
}
metadata <- data.frame(
  Sample = sample_names,
  AnimalID = paste0("A", as.integer(extract_one("_A0*([0-9]+)_", sample_names))),
  Region = toupper(sample_info$region),
  StressGroup = toupper(sample_info$ExpGroup),
  Hemisphere = extract_one("_(L|R)_(?:CA1|CA2|CA3|DG)_", sample_names),
  AcquisitionBatch = extract_one("Bluto_([0-9]{8})_", sample_names),
  stringsAsFactors = FALSE
)
if (anyNA(metadata)) stop("Sample metadata could not be reconstructed deterministically from the current state")
rownames(metadata) <- sample_names
metadata[] <- lapply(metadata, function(x) if (is.character(x)) x else x)
for (nm in setdiff(names(metadata), "Sample")) metadata[[nm]] <- factor(metadata[[nm]])
if (length(unique(metadata$AnimalID)) != 9L || any(table(metadata$AnimalID) != 8L)) {
  stop("Repeated-measures contract is not 9 animals x 8 ROIs")
}

current_me <- state$mergedMEs[, module_map$module_eigengene, drop = FALSE]
colnames(current_me) <- module_map$ModuleID
if (!identical(rownames(current_me), sample_names)) stop("Primary module eigengenes are not sample-aligned")

super_values <- utils::read.csv(super_values_path, check.names = FALSE)
super_values <- super_values[super_values$dataset == "microglia", , drop = FALSE]
if (nrow(super_values) != 72L * 9L || anyDuplicated(super_values[c("Sample", "supermodule_id")])) {
  stop("Stage 05 supermodule eigengene table is not one row per sample and SM identity")
}
super_ids <- sprintf("SM%02d", seq_len(9))
super_me <- sapply(super_ids, function(sm) {
  z <- super_values[super_values$supermodule_id == sm, ]
  z$eigengene[match(sample_names, z$Sample)]
})
colnames(super_me) <- super_ids
rownames(super_me) <- sample_names
if (anyNA(super_me)) stop("Stage 05 supermodule eigengenes are incompletely aligned")

metadata_export <- metadata
metadata_export$AnimalID <- as.character(metadata_export$AnimalID)
metadata_export$Region <- as.character(metadata_export$Region)
metadata_export$StressGroup <- as.character(metadata_export$StressGroup)
metadata_export$Hemisphere <- as.character(metadata_export$Hemisphere)
metadata_export$AcquisitionBatch <- as.character(metadata_export$AcquisitionBatch)
metadata_export$repeated_measure_unit <- "AnimalID"
metadata_export$n_ROI_per_animal <- as.integer(table(metadata$AnimalID)[metadata$AnimalID])
write_csv_atomic(metadata_export, file.path(output_dir, "sample_metadata_audit.csv"))

# -------------------------------------------------------------------------
# Part 1: mixed-model variance decomposition.
# AcquisitionBatch is fully nested in AnimalID and therefore excluded from
# attribution rather than receiving an arbitrary component.
# -------------------------------------------------------------------------

batch_per_animal <- tapply(as.character(metadata$AcquisitionBatch), metadata$AnimalID, function(x) length(unique(x)))
animal_per_batch <- table(metadata$AcquisitionBatch, metadata$AnimalID)
batch_nested <- all(batch_per_animal == 1L)
vp_formula_text <- "~ (1|AnimalID) + (1|Region) + (1|StressGroup) + (1|Hemisphere)"
vp_formula <- stats::as.formula(vp_formula_text)

estimability_base <- data.frame(
  term = c("AnimalID", "Region", "StressGroup", "Hemisphere", "AcquisitionBatch", "Residuals"),
  requested = TRUE,
  included_in_formula = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  n_levels = c(nlevels(metadata$AnimalID), nlevels(metadata$Region), nlevels(metadata$StressGroup),
               nlevels(metadata$Hemisphere), nlevels(metadata$AcquisitionBatch), NA_integer_),
  estimability = c(
    "estimable_random_intercept",
    "estimable_random_intercept",
    "estimable_but_low_level_and_nested_within_animal",
    "estimable_within_animal",
    if (batch_nested) "not_estimable_separately_from_animal" else "not_included",
    "estimable"
  ),
  handling = c(
    "variance_component_reported",
    "variance_component_reported",
    "variance_component_reported_with_low_level_caution",
    "variance_component_reported",
    "excluded_and_variance_fraction_left_NA",
    "variance_component_reported"
  ),
  limitation = c(
    "nine biological units; ROI observations are repeated measures",
    "four hippocampal regions observed within animals",
    "three group levels; group is constant within animal and estimates are descriptive",
    "two levels observed within animal",
    "two acquisition dates; each animal occurs in exactly one date, so separate attribution is aliased",
    "unexplained ROI-level variance"
  ),
  stringsAsFactors = FALSE
)

fit_variance_partition <- function(me_matrix, level_name, entity_ids) {
  vp <- variancePartition::fitExtractVarPartModel(t(me_matrix), vp_formula, metadata)
  vp <- as.data.frame(vp)
  component_order <- c("AnimalID", "Region", "StressGroup", "Hemisphere", "Residuals")
  missing_components <- setdiff(component_order, names(vp))
  if (length(missing_components)) stop("Variance partition omitted components: ", paste(missing_components, collapse = ", "))
  rows <- do.call(rbind, lapply(seq_along(entity_ids), function(i) {
    values <- as.numeric(vp[i, component_order])
    negative <- values < -sqrt(.Machine$double.eps)
    values[values < 0] <- 0
    total <- sum(values)
    if (!is.finite(total) || total <= 0) stop("Invalid variance decomposition for ", entity_ids[[i]])
    values <- values / total
    data.frame(
      dataset = "microglia", level = level_name, entity_id = entity_ids[[i]],
      component = c(component_order, "AcquisitionBatch"),
      variance_fraction = c(values, NA_real_),
      estimable = c(rep(TRUE, length(values)), FALSE),
      boundary_estimate = c(values <= 1e-8, NA),
      negative_estimate_detected = c(negative, NA),
      negative_estimate_handling = c(rep("truncated_to_zero_then_renormalized_if_present", length(values)),
                                     "not_applicable_non_estimable"),
      model_formula = vp_formula_text,
      method = "variancePartition::fitExtractVarPartModel; categorical random-intercept variance fractions",
      biological_unit = "AnimalID (n=9)",
      stringsAsFactors = FALSE
    )
  }))
  rows
}

module_vp <- fit_variance_partition(current_me, "module", module_ids)
super_vp <- fit_variance_partition(super_me, "supermodule", super_ids)

fit_model_diagnostics <- function(me_matrix, level_name, entity_ids) {
  do.call(rbind, lapply(seq_along(entity_ids), function(i) {
    d <- metadata
    d$y <- me_matrix[, i]
    warnings <- character()
    fit <- withCallingHandlers(
      lme4::lmer(stats::as.formula(paste("y", vp_formula_text)), data = d, REML = TRUE),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
    vc <- as.data.frame(lme4::VarCorr(fit))
    conv <- fit@optinfo$conv$lme4$messages
    data.frame(
      dataset = "microglia", level = level_name, entity_id = entity_ids[[i]],
      model_formula = paste("eigengene", vp_formula_text),
      n_observations = nrow(d), n_animals = nlevels(d$AnimalID),
      singular_fit = lme4::isSingular(fit, tol = 1e-5),
      convergence_message = if (is.null(conv)) "none" else paste(conv, collapse = " | "),
      warning = if (length(warnings)) paste(unique(warnings), collapse = " | ") else "none",
      boundary_components = paste(vc$grp[vc$vcov <= 1e-8], collapse = ";"),
      negative_variance_components = 0L,
      acquisition_batch_estimable = FALSE,
      acquisition_batch_reason = "fully nested within AnimalID; excluded from attribution",
      stress_group_caution = "three levels and constant within AnimalID; descriptive variance component only",
      stringsAsFactors = FALSE
    )
  }))
}

rbind_fill <- function(...) {
  parts <- list(...)
  all_names <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[all_names]
  })
  do.call(rbind, parts)
}
model_audit <- rbind_fill(
  transform(estimability_base, dataset = "microglia", level = "global", entity_id = "all"),
  fit_model_diagnostics(current_me, "module", module_ids),
  fit_model_diagnostics(super_me, "supermodule", super_ids)
)
write_csv_atomic(module_vp, file.path(output_dir, "module_eigengene_variance_partition.csv"))
write_csv_atomic(super_vp, file.path(output_dir, "supermodule_eigengene_variance_partition.csv"))
write_csv_atomic(model_audit, file.path(output_dir, "model_estimability_audit.csv"))

vp_plot_data <- rbind(module_vp, super_vp)
vp_plot_data <- vp_plot_data[vp_plot_data$estimable, ]
vp_plot_data$entity_id <- factor(vp_plot_data$entity_id, levels = rev(c(module_ids, super_ids)))
vp_plot_data$component <- factor(vp_plot_data$component,
                                 levels = c("AnimalID", "StressGroup", "Region", "Hemisphere", "Residuals"))
vp_plot <- ggplot2::ggplot(vp_plot_data, ggplot2::aes(component, entity_id, fill = variance_fraction)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
  ggplot2::facet_grid(level ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "C", name = "Variance\nfraction") +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_family = "Arial", base_size = 8) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
                 axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
save_plot_atomic(vp_plot, file.path(output_dir, "variance_partition_heatmap.svg"), 6.4, 7.4)
save_plot_atomic(vp_plot, file.path(output_dir, "variance_partition_heatmap.pdf"), 6.4, 7.4)

# -------------------------------------------------------------------------
# Part 2: fixed-universe sensitivity matrices.
# -------------------------------------------------------------------------

matrix_rank <- function(x, tol = 1e-8) qr(x, tol = tol)$rank
projection <- function(x) {
  q <- qr(x)
  qr.Q(q)[, seq_len(q$rank), drop = FALSE]
}

collapse_key <- interaction(metadata$AnimalID, metadata$Region, drop = TRUE, lex.order = TRUE)
collapse_levels <- levels(collapse_key)
expression_collapsed <- do.call(rbind, lapply(collapse_levels, function(k) {
  colMeans(expression_primary[collapse_key == k, , drop = FALSE])
}))
collapsed_meta <- do.call(rbind, lapply(collapse_levels, function(k) {
  ix <- which(collapse_key == k)
  data.frame(
    Sample = paste(as.character(metadata$AnimalID[ix[1]]), as.character(metadata$Region[ix[1]]), sep = "__"),
    AnimalID = as.character(metadata$AnimalID[ix[1]]), Region = as.character(metadata$Region[ix[1]]),
    StressGroup = as.character(metadata$StressGroup[ix[1]]), Hemisphere = "collapsed_LR",
    AcquisitionBatch = as.character(metadata$AcquisitionBatch[ix[1]]), n_primary_ROI = length(ix),
    stringsAsFactors = FALSE
  )
}))
rownames(expression_collapsed) <- collapsed_meta$Sample
colnames(expression_collapsed) <- colnames(expression_primary)

x_protected <- stats::model.matrix(~ 0 + AnimalID, metadata)
x_nuisance_raw <- stats::model.matrix(~ Region + Hemisphere + AcquisitionBatch, metadata)
q_protected <- projection(x_protected)
x_nuisance_orth <- x_nuisance_raw - q_protected %*% crossprod(q_protected, x_nuisance_raw)
keep_nuisance <- colSums(x_nuisance_orth^2) > 1e-10
x_nuisance_orth <- x_nuisance_orth[, keep_nuisance, drop = FALSE]
q_nuisance <- projection(x_nuisance_orth)
expression_spatial <- expression_primary - q_nuisance %*% crossprod(q_nuisance, expression_primary)
rownames(expression_spatial) <- rownames(expression_primary)
colnames(expression_spatial) <- colnames(expression_primary)

x_strict <- stats::model.matrix(~ 0 + AnimalID + Region, metadata)
q_strict <- projection(x_strict)
expression_strict <- expression_primary - q_strict %*% crossprod(q_strict, expression_primary)
rownames(expression_strict) <- rownames(expression_primary)
colnames(expression_strict) <- colnames(expression_primary)

protected_projection_delta <- max(abs(crossprod(q_protected, expression_spatial) - crossprod(q_protected, expression_primary)))
group_design <- stats::model.matrix(~ 0 + StressGroup, metadata)
q_group <- projection(group_design)
group_projection_delta <- max(abs(crossprod(q_group, expression_spatial) - crossprod(q_group, expression_primary)))

sensitivity_audit <- data.frame(
  sensitivity_id = c("primary", "hemisphere_collapsed", "spatial_adjusted", "within_animal_spatial_adjusted"),
  purpose = c("immutable reference", "bilateral sensitivity", "spatial nuisance sensitivity preserving animal/group", "strict within-animal spatial sensitivity"),
  n_samples = c(nrow(expression_primary), nrow(expression_collapsed), nrow(expression_spatial), nrow(expression_strict)),
  n_animals = c(9L, 9L, 9L, 9L),
  n_features = c(ncol(expression_primary), ncol(expression_collapsed), ncol(expression_spatial), ncol(expression_strict)),
  feature_order_exact = c(TRUE, identical(colnames(expression_primary), colnames(expression_collapsed)),
                          identical(colnames(expression_primary), colnames(expression_spatial)),
                          identical(colnames(expression_primary), colnames(expression_strict))),
  missing_values = c(sum(is.na(expression_primary)), sum(is.na(expression_collapsed)), sum(is.na(expression_spatial)), sum(is.na(expression_strict))),
  residualization_formula = c("none", "mean within AnimalID + Region", "orthogonalized nuisance: Region + Hemisphere + AcquisitionBatch | protected AnimalID", "residuals of ~ 0 + AnimalID + Region"),
  design_rank = c(NA_integer_, NA_integer_, matrix_rank(x_nuisance_orth), matrix_rank(x_strict)),
  design_columns = c(NA_integer_, NA_integer_, ncol(x_nuisance_orth), ncol(x_strict)),
  full_rank = c(NA, NA, matrix_rank(x_nuisance_orth) == ncol(x_nuisance_orth), matrix_rank(x_strict) == ncol(x_strict)),
  animal_projection_max_abs_delta = c(NA, NA, protected_projection_delta, NA),
  stress_group_projection_max_abs_delta = c(NA, NA, group_projection_delta, NA),
  group_effect_leakage = c("not_applicable", "group preserved through animal-level averaging", if (group_projection_delta < 1e-8) "none_detected_linear_projection" else "detected", "removed_by_design_strict_sensitivity"),
  permitted_use = c("reference_only", "sensitivity_only", "sensitivity_only", "diagnostic_only_never_primary"),
  stringsAsFactors = FALSE
)
write_csv_atomic(sensitivity_audit, file.path(output_dir, "sensitivity_matrix_validation.csv"))

residualization_estimability <- data.frame(
  sensitivity_id = c(
    "hemisphere_collapsed", "hemisphere_collapsed",
    "spatial_adjusted", "spatial_adjusted", "spatial_adjusted", "spatial_adjusted",
    "within_animal_spatial_adjusted", "within_animal_spatial_adjusted"
  ),
  term = c(
    "AnimalID", "Region", "Region", "Hemisphere", "AcquisitionBatch", "AnimalID_protected_space",
    "AnimalID", "Region"
  ),
  requested_role = c(
    "averaging_key", "averaging_key", "nuisance", "nuisance", "nuisance", "protected",
    "strict_nuisance", "strict_nuisance"
  ),
  estimability = c(
    "estimable", "estimable", "estimable_after_orthogonalization", "estimable_after_orthogonalization",
    "not_estimable_after_orthogonalization_to_AnimalID", "estimable_full_rank",
    "estimable_as_fixed_effect_design", "estimable_conditional_on_AnimalID"
  ),
  degrees_of_freedom = c(8L, 3L, 3L, 1L, 0L, 9L, 9L, 3L),
  included_in_residualization = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
  handling = c(
    "left/right averaged within AnimalID+Region", "left/right averaged within AnimalID+Region",
    "orthogonal nuisance component subtracted", "orthogonal nuisance component subtracted",
    "zero after projection because batch is nested in AnimalID; no variance removed",
    "projection preserved exactly; group is contained in this space",
    "removed only in strict diagnostic matrix", "removed only in strict diagnostic matrix"
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(residualization_estimability, file.path(output_dir, "sensitivity_residualization_estimability.csv"))

matrix_to_df <- function(mat, meta) {
  cbind(meta, as.data.frame(mat, check.names = FALSE))
}
write_gz_csv_atomic(matrix_to_df(expression_collapsed, collapsed_meta), file.path(output_dir, "sensitivity_hemisphere_collapsed.csv.gz"))
write_gz_csv_atomic(matrix_to_df(expression_spatial, metadata_export), file.path(output_dir, "sensitivity_spatial_adjusted.csv.gz"))
write_gz_csv_atomic(matrix_to_df(expression_strict, metadata_export), file.path(output_dir, "sensitivity_within_animal_spatial_adjusted.csv.gz"))

feature_audit <- data.frame(
  feature_index = seq_along(colnames(expression_primary)),
  ProteinGroupID = colnames(expression_primary),
  ModuleID = as.character(module_long$ModuleID),
  present_primary = TRUE, present_hemisphere_collapsed = TRUE,
  present_spatial_adjusted = TRUE, present_within_animal_spatial_adjusted = TRUE,
  order_identical = TRUE,
  stringsAsFactors = FALSE
)
write_csv_atomic(feature_audit, file.path(output_dir, "sensitivity_feature_identity_audit.csv"))

# -------------------------------------------------------------------------
# Part 3: immutable-set module robustness and leave-one-animal-out checks.
# -------------------------------------------------------------------------

bicor_safe <- function(x, y = NULL) {
  if (is.null(y)) WGCNA::bicor(x, use = "pairwise.complete.obs", maxPOutliers = 0.05)
  else as.numeric(WGCNA::bicor(x, y, use = "pairwise.complete.obs", maxPOutliers = 0.05))
}

pc1_scores <- function(x) {
  fit <- stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = 1)
  as.numeric(fit$x[, 1])
}

fixed_module_mes <- function(expr) {
  out <- sapply(module_ids, function(mid) pc1_scores(expr[, as.character(module_long$ModuleID) == mid, drop = FALSE]))
  colnames(out) <- module_ids
  rownames(out) <- rownames(expr)
  out
}

align_score <- function(score, reference) {
  r <- suppressWarnings(stats::cor(score, reference, use = "pairwise.complete.obs"))
  if (is.finite(r) && r < 0) -score else score
}

connectivity <- function(expr_module, power = state$softPower) {
  cor_mat <- bicor_safe(expr_module)
  adj <- ((1 + cor_mat) / 2)^power
  diag(adj) <- 0
  rowSums(adj)
}

primary_connectivity <- setNames(vector("list", length(module_ids)), module_ids)
primary_abs_kme <- setNames(vector("list", length(module_ids)), module_ids)
for (mid in module_ids) {
  ix <- as.character(module_long$ModuleID) == mid
  primary_connectivity[[mid]] <- connectivity(expression_primary[, ix, drop = FALSE])
  primary_abs_kme[[mid]] <- abs(apply(expression_primary[, ix, drop = FALSE], 2, bicor_safe, y = current_me[, mid]))
  names(primary_connectivity[[mid]]) <- colnames(expression_primary)[ix]
  names(primary_abs_kme[[mid]]) <- colnames(expression_primary)[ix]
}

primary_collapsed_me <- do.call(rbind, lapply(collapse_levels, function(k) colMeans(current_me[collapse_key == k, , drop = FALSE])))
rownames(primary_collapsed_me) <- rownames(expression_collapsed)

sensitivity_expr <- list(
  hemisphere_collapsed = expression_collapsed,
  spatial_adjusted = expression_spatial,
  within_animal_spatial_adjusted = expression_strict
)
sensitivity_ref_me <- list(
  hemisphere_collapsed = primary_collapsed_me,
  spatial_adjusted = current_me,
  within_animal_spatial_adjusted = current_me
)
sensitivity_mes <- lapply(sensitivity_expr, fixed_module_mes)

metric_rows <- list()
for (sid in names(sensitivity_expr)) {
  expr <- sensitivity_expr[[sid]]
  mes <- sensitivity_mes[[sid]]
  ref_mes <- sensitivity_ref_me[[sid]]
  for (mid in module_ids) {
    ix <- as.character(module_long$ModuleID) == mid
    aligned <- align_score(mes[, mid], ref_mes[, mid])
    mes[, mid] <- aligned
    ktest <- connectivity(expr[, ix, drop = FALSE])
    names(ktest) <- colnames(expr)[ix]
    kref <- primary_connectivity[[mid]][names(ktest)]
    top_ref <- names(sort(kref, decreasing = TRUE))[seq_len(min(25L, length(kref)))]
    top_test <- names(sort(ktest, decreasing = TRUE))[seq_len(min(25L, length(ktest)))]
    rank_r <- suppressWarnings(stats::cor(rank(-kref[top_ref]), rank(-ktest[top_ref]), method = "spearman"))
    metric_rows[[length(metric_rows) + 1L]] <- data.frame(
      dataset = "microglia", sensitivity_id = sid, ModuleID = mid,
      feature_count = sum(ix),
      module_eigengene_correlation = bicor_safe(aligned, ref_mes[, mid]),
      intramodular_connectivity_spearman = suppressWarnings(stats::cor(kref, ktest, method = "spearman")),
      primary_top25_rank_spearman = rank_r,
      primary_top25_retained_fraction = length(intersect(top_ref, top_test)) / length(top_ref),
      preservation_Zsummary = NA_real_, preservation_median_rank = NA_real_,
      preservation_metric_reason = "pending_modulePreservation",
      stringsAsFactors = FALSE
    )
  }
  sensitivity_mes[[sid]] <- mes
}
module_robustness <- do.call(rbind, metric_rows)

run_preservation <- function(test_expr, test_name) {
  cache_path <- file.path(output_dir, paste0("preservation_", test_name, "_", n_permutations, "perm.csv"))
  if (file.exists(cache_path) && file.info(cache_path)$size > 0) {
    cached <- utils::read.csv(cache_path, check.names = FALSE)
    if (nrow(cached) == 13L && setequal(cached$ModuleID, module_ids) &&
        all(cached$preservation_n_permutations == n_permutations) &&
        "reference_state_md5" %in% names(cached) &&
        all(cached$reference_state_md5 == reference_state_md5) &&
        "preservation_cache_contract" %in% names(cached) &&
        all(cached$preservation_cache_contract == preservation_cache_contract)) return(cached)
  }
  common <- intersect(colnames(expression_primary), colnames(test_expr))
  colors <- as.character(module_long$ModuleID[match(common, module_long$ProteinGroupID)])
  multi_expr <- list(
    primary = list(data = expression_primary[, common, drop = FALSE]),
    sensitivity = list(data = test_expr[, common, drop = FALSE])
  )
  multi_color <- list(primary = colors, sensitivity = colors)
  set.seed(12345)
  mp <- WGCNA::modulePreservation(
    multiData = multi_expr, multiColor = multi_color,
    referenceNetworks = 1, nPermutations = n_permutations,
    networkType = "signed", corFnc = "bicor",
    corOptions = "use='p', maxPOutliers=0.05",
    randomSeed = 12345, parallelCalculation = FALSE,
    maxModuleSize = ncol(expression_primary),
    greyName = "__grey_unused", goldName = "__gold_unused",
    savePermutedStatistics = FALSE, useInterpolation = FALSE,
    verbose = 1
  )
  z <- mp$preservation$Z[[1]][[2]]
  observed <- mp$preservation$observed[[1]][[2]]
  keep_modules <- intersect(module_ids, rownames(z))
  out <- data.frame(
    sensitivity_id = test_name,
    ModuleID = keep_modules,
    preservation_Zsummary = z[keep_modules, "Zsummary.pres"],
    preservation_median_rank = observed[keep_modules, "medianRank.pres"],
    preservation_metric_reason = if (n_permutations > 0) "available" else "permutations_disabled",
    preservation_n_permutations = n_permutations,
    reference_state_md5 = reference_state_md5,
    preservation_cache_contract = preservation_cache_contract,
    stringsAsFactors = FALSE
  )
  write_csv_atomic(out, cache_path)
  out
}

preservation_rows <- lapply(names(sensitivity_expr), function(sid) run_preservation(sensitivity_expr[[sid]], sid))
preservation <- do.call(rbind, preservation_rows)
key_rob <- paste(module_robustness$sensitivity_id, module_robustness$ModuleID)
key_pre <- paste(preservation$sensitivity_id, preservation$ModuleID)
match_pre <- match(key_rob, key_pre)
module_robustness$preservation_Zsummary <- preservation$preservation_Zsummary[match_pre]
module_robustness$preservation_median_rank <- preservation$preservation_median_rank[match_pre]
module_robustness$preservation_metric_reason <- preservation$preservation_metric_reason[match_pre]

status_one <- function(me, conn, hubs, z) {
  if (is.finite(me) && is.finite(conn) && is.finite(hubs) && is.finite(z) &&
      me >= 0.80 && conn >= 0.60 && hubs >= 0.60 && z >= 10) return("pass")
  if (is.finite(me) && is.finite(conn) && is.finite(hubs) && is.finite(z) &&
      me >= 0.60 && conn >= 0.30 && hubs >= 0.40 && z >= 2) return("suggestive")
  "fail"
}
module_robustness$status <- mapply(
  status_one,
  module_robustness$module_eigengene_correlation,
  module_robustness$intramodular_connectivity_spearman,
  module_robustness$primary_top25_retained_fraction,
  module_robustness$preservation_Zsummary
)
module_robustness$status_rule <- "pass: ME>=0.80, connectivity>=0.60, top25>=0.60, Z>=10; suggestive: ME>=0.60, connectivity>=0.30, top25>=0.40, Z>=2; otherwise fail"
write_csv_atomic(module_robustness, file.path(output_dir, "module_sensitivity_robustness.csv"))

loao_rows <- list()
for (animal in levels(metadata$AnimalID)) {
  keep <- metadata$AnimalID != animal
  expr <- expression_primary[keep, , drop = FALSE]
  mes <- fixed_module_mes(expr)
  for (mid in module_ids) {
    ix <- as.character(module_long$ModuleID) == mid
    aligned <- align_score(mes[, mid], current_me[keep, mid])
    test_kme <- abs(apply(expr[, ix, drop = FALSE], 2, bicor_safe, y = aligned))
    names(test_kme) <- colnames(expr)[ix]
    ref_kme <- primary_abs_kme[[mid]][names(test_kme)]
    top_ref <- names(sort(ref_kme, decreasing = TRUE))[seq_len(min(25L, length(ref_kme)))]
    top_test <- names(sort(test_kme, decreasing = TRUE))[seq_len(min(25L, length(test_kme)))]
    loao_rows[[length(loao_rows) + 1L]] <- data.frame(
      dataset = "microglia", omitted_AnimalID = animal, ModuleID = mid,
      n_remaining_animals = 8L, n_remaining_ROI = sum(keep), feature_count = sum(ix),
      eigengene_correlation = bicor_safe(aligned, current_me[keep, mid]),
      hub_rank_spearman = suppressWarnings(stats::cor(ref_kme, test_kme, method = "spearman")),
      primary_top25_retained_fraction = length(intersect(top_ref, top_test)) / length(top_ref),
      membership_status = "immutable_reference_set_no_reassignment",
      stringsAsFactors = FALSE
    )
  }
}
loao <- do.call(rbind, loao_rows)
loao$influential_animal_flag <- with(loao, eigengene_correlation < 0.80 | hub_rank_spearman < 0.60 | primary_top25_retained_fraction < 0.60)
loao_summary <- do.call(rbind, lapply(split(loao, loao$ModuleID), function(z) {
  weakest <- z$omitted_AnimalID[which.min(pmin(z$eigengene_correlation, z$hub_rank_spearman, z$primary_top25_retained_fraction))]
  data.frame(
    dataset = "microglia", ModuleID = z$ModuleID[1], n_leave_one_animal_out = nrow(z),
    min_eigengene_correlation = min(z$eigengene_correlation),
    median_eigengene_correlation = median(z$eigengene_correlation),
    max_eigengene_correlation = max(z$eigengene_correlation),
    min_hub_rank_spearman = min(z$hub_rank_spearman),
    median_hub_rank_spearman = median(z$hub_rank_spearman),
    min_primary_top25_retained_fraction = min(z$primary_top25_retained_fraction),
    n_influential_animal_flags = sum(z$influential_animal_flag),
    weakest_omitted_AnimalID = weakest,
    member_stability_interpretation = "fixed membership; hub/kME stability only, no reassignment claim",
    stringsAsFactors = FALSE
  )
}))
write_csv_atomic(loao, file.path(output_dir, "module_leave_one_animal_out.csv"))
write_csv_atomic(loao_summary, file.path(output_dir, "module_leave_one_animal_out_summary.csv"))

module_consensus <- do.call(rbind, lapply(split(module_robustness, module_robustness$ModuleID), function(z) {
  rank_status <- c(pass = 3L, suggestive = 2L, fail = 1L)
  worst <- names(which.min(rank_status[z$status]))
  lz <- loao_summary[loao_summary$ModuleID == z$ModuleID[1], ]
  data.frame(
    dataset = "microglia", ModuleID = z$ModuleID[1], feature_count = z$feature_count[1],
    worst_sensitivity_status = worst,
    min_sensitivity_ME_correlation = min(z$module_eigengene_correlation, na.rm = TRUE),
    min_sensitivity_connectivity_spearman = min(z$intramodular_connectivity_spearman, na.rm = TRUE),
    min_sensitivity_top25_retained = min(z$primary_top25_retained_fraction, na.rm = TRUE),
    min_preservation_Zsummary = min(z$preservation_Zsummary, na.rm = TRUE),
    max_preservation_median_rank = max(z$preservation_median_rank, na.rm = TRUE),
    loao_min_ME_correlation = lz$min_eigengene_correlation,
    loao_min_hub_rank_spearman = lz$min_hub_rank_spearman,
    loao_min_top25_retained = lz$min_primary_top25_retained_fraction,
    influential_animal = lz$n_influential_animal_flags > 0,
    stringsAsFactors = FALSE
  )
}))
write_csv_atomic(module_consensus, file.path(output_dir, "module_robustness_consensus.csv"))

# -------------------------------------------------------------------------
# Part 4: higher-order blocks; all other identities remain standalone.
# -------------------------------------------------------------------------

matrix_metadata <- list(
  primary = metadata,
  hemisphere_collapsed = transform(collapsed_meta,
                                   AnimalID = factor(AnimalID), Region = factor(Region),
                                   StressGroup = factor(StressGroup), Hemisphere = factor(Hemisphere),
                                   AcquisitionBatch = factor(AcquisitionBatch)),
  spatial_adjusted = metadata,
  within_animal_spatial_adjusted = metadata
)
all_mes <- c(list(primary = current_me), sensitivity_mes)

partial_residual <- function(y, d) {
  usable <- intersect(c("Region", "Hemisphere"), names(d))
  usable <- usable[vapply(d[usable], function(x) length(unique(x)) > 1L, logical(1))]
  if (!length(usable)) return(y)
  stats::residuals(stats::lm(stats::reformulate(usable, response = "y"), data = cbind(d, y = y)))
}

cluster_support <- function(me, members, cut_height = 0.40) {
  cor_mat <- bicor_safe(me)
  hc <- stats::hclust(stats::as.dist(1 - cor_mat), method = "average")
  clusters <- stats::cutree(hc, h = cut_height)
  same <- length(unique(clusters[members])) == 1L
  exclusive <- same && all(names(clusters)[clusters == clusters[members[1]]] %in% members)
  within <- cor_mat[members, members, drop = FALSE]
  min_within <- if (length(members) == 2L) within[1, 2] else min(within[upper.tri(within)])
  outsiders <- setdiff(colnames(me), members)
  max_out <- max(cor_mat[members, outsiders, drop = FALSE])
  c(same_group = same, exclusive_group = exclusive, min_within = min_within,
    max_member_outsider = max_out, separation_margin = min_within - max_out)
}

block_rows <- list()
for (sid in names(all_mes)) {
  me <- all_mes[[sid]]
  d <- matrix_metadata[[sid]]
  for (sm in names(expected_blocks)) {
    members <- expected_blocks[[sm]]
    cor_mat <- bicor_safe(me[, members, drop = FALSE])
    pair_values <- cor_mat[upper.tri(cor_mat)]
    adj_me <- apply(me[, members, drop = FALSE], 2, partial_residual, d = d)
    adj_cor <- bicor_safe(adj_me)
    adj_values <- adj_cor[upper.tri(adj_cor)]
    pc <- stats::prcomp(me[, members, drop = FALSE], center = TRUE, scale. = TRUE)
    support <- cluster_support(me, members)
    block_rows[[length(block_rows) + 1L]] <- data.frame(
      dataset = "microglia", sensitivity_id = sid, SupermoduleID = sm,
      member_ModuleIDs = paste(members, collapse = ";"), n_member_modules = length(members),
      min_member_eigengene_correlation = min(pair_values),
      median_member_eigengene_correlation = median(pair_values),
      min_spatial_adjusted_member_correlation = min(adj_values),
      median_spatial_adjusted_member_correlation = median(adj_values),
      PC1_variance_explained = summary(pc)$importance[2, 1],
      same_group_at_cut_0_40 = as.logical(support["same_group"]),
      exclusive_group_at_cut_0_40 = as.logical(support["exclusive_group"]),
      min_within_correlation = as.numeric(support["min_within"]),
      max_member_to_outsider_correlation = as.numeric(support["max_member_outsider"]),
      cut_height_independent_separation_margin = as.numeric(support["separation_margin"]),
      cut_height_independent_support = as.numeric(support["separation_margin"]) > 0 && as.numeric(support["min_within"]) > 0,
      cut_height_use = "descriptive_only_no_optimization",
      stringsAsFactors = FALSE
    )
  }
}
block_robustness <- do.call(rbind, block_rows)

block_loao_rows <- list()
for (animal in levels(metadata$AnimalID)) {
  keep <- metadata$AnimalID != animal
  me <- fixed_module_mes(expression_primary[keep, , drop = FALSE])
  for (mid in module_ids) me[, mid] <- align_score(me[, mid], current_me[keep, mid])
  for (sm in names(expected_blocks)) {
    members <- expected_blocks[[sm]]
    cor_mat <- bicor_safe(me[, members, drop = FALSE])
    pair_values <- cor_mat[upper.tri(cor_mat)]
    pc <- stats::prcomp(me[, members, drop = FALSE], center = TRUE, scale. = TRUE)
    support <- cluster_support(me, members)
    block_loao_rows[[length(block_loao_rows) + 1L]] <- data.frame(
      omitted_AnimalID = animal, SupermoduleID = sm,
      min_member_eigengene_correlation = min(pair_values),
      PC1_variance_explained = summary(pc)$importance[2, 1],
      same_group_at_cut_0_40 = as.logical(support["same_group"]),
      separation_margin = as.numeric(support["separation_margin"]),
      stringsAsFactors = FALSE
    )
  }
}
block_loao <- do.call(rbind, block_loao_rows)
block_loao_summary <- do.call(rbind, lapply(split(block_loao, block_loao$SupermoduleID), function(z) {
  data.frame(
    SupermoduleID = z$SupermoduleID[1], n_leave_one_animal_out = nrow(z),
    min_member_correlation = min(z$min_member_eigengene_correlation),
    median_member_correlation = median(z$min_member_eigengene_correlation),
    min_PC1_variance = min(z$PC1_variance_explained),
    grouping_persistence_fraction = mean(z$same_group_at_cut_0_40),
    min_separation_margin = min(z$separation_margin),
    weakest_omitted_AnimalID = z$omitted_AnimalID[which.min(z$separation_margin)],
    stringsAsFactors = FALSE
  )
}))
standalone_audit <- data.frame(
  SupermoduleID = singleton_ids,
  member_ModuleID = vapply(observed_members[singleton_ids], `[`, character(1), 1L),
  structural_status = "standalone_singleton_module",
  higher_order_block_metrics = "not_applicable",
  reason = "one member module; not treated as a multi-module higher-order block",
  stringsAsFactors = FALSE
)
write_csv_atomic(block_robustness, file.path(output_dir, "higher_order_block_robustness.csv"))
write_csv_atomic(block_loao, file.path(output_dir, "higher_order_block_leave_one_animal_out.csv"))
write_csv_atomic(block_loao_summary, file.path(output_dir, "higher_order_block_leave_one_animal_out_summary.csv"))
write_csv_atomic(standalone_audit, file.path(output_dir, "standalone_supermodule_identity_audit.csv"))

# -------------------------------------------------------------------------
# Part 5: same-study compartment and marker-panel context.
# -------------------------------------------------------------------------

annotation <- utils::read.csv(annotation_path, check.names = FALSE)
annotation <- annotation[match(module_ids, annotation$ModuleID), , drop = FALSE]
registry <- utils::read.csv(registry_path, check.names = FALSE)
registry_modules <- registry[registry$dataset == "microglia" & registry$level == "module", ]
registry_modules <- registry_modules[match(module_ids, registry_modules$entity_id), , drop = FALSE]
if (nrow(registry_modules) != 13L || anyNA(registry_modules$reviewed_biological_label)) stop("Reviewed module registry contract mismatch")

joint <- readRDS(joint_path)
joint_matrix <- as.matrix(joint$primary$matrix)
joint_features <- joint$feature_table
active_map <- module_long[c("ProteinGroupID", "ModuleID", "member_accessions", "GeneSymbol")]
joint_key_counts <- table(joint_features$member_accessions)
joint_unique <- joint_features[!is.na(joint_features$member_accessions) & joint_key_counts[joint_features$member_accessions] == 1L, ]
joint_match <- match(active_map$member_accessions, joint_unique$member_accessions)
active_map$joint_ProteinGroupID <- joint_unique$ProteinGroupID[joint_match]
active_map$joint_match_status <- ifelse(is.na(joint_match), "no_unique_exact_member_accessions_match", "unique_exact_member_accessions_match")
joint_row <- match(active_map$joint_ProteinGroupID, rownames(joint_matrix))
active_map$joint_primary_matrix_present <- !is.na(joint_row)

compartments <- c("microglia", "neuron_neuropil", "neuron_soma")
comp_means <- matrix(NA_real_, nrow(active_map), length(compartments), dimnames = list(NULL, compartments))
for (comp in compartments) {
  cols <- joint$metadata$dataset == comp
  ok <- !is.na(joint_row)
  comp_means[ok, comp] <- rowMeans(joint_matrix[joint_row[ok], cols, drop = FALSE], na.rm = TRUE)
}
active_map <- cbind(active_map, as.data.frame(comp_means))
active_map$microglia_vs_neuropil <- active_map$microglia - active_map$neuron_neuropil
active_map$microglia_vs_soma <- active_map$microglia - active_map$neuron_soma
active_map$dominant_compartment <- ifelse(
  apply(comp_means, 1, function(x) all(is.na(x))), NA_character_,
  compartments[max.col(comp_means, ties.method = "first")]
)
write_csv_atomic(active_map, file.path(output_dir, "protein_compartment_specificity_source.csv"))

split_genes <- function(x) unique(toupper(trimws(unlist(strsplit(paste(x, collapse = ";"), "[;,|]")))))
universe_genes <- split_genes(module_long$GeneSymbol)
reference_markers <- utils::read.csv(reference_markers_path, check.names = FALSE)
microenv_markers <- utils::read.csv(microenvironment_markers_path, check.names = FALSE)
empirical_markers <- utils::read.csv(empirical_markers_path, check.names = FALSE)

panel_sets <- list(
  canonical_microglia = unique(toupper(reference_markers$gene_symbol[grepl("microglia|pvm|myeloid", paste(reference_markers$marker_set, reference_markers$cell_class), ignore.case = TRUE)])),
  neuronal_neuropil = unique(toupper(reference_markers$gene_symbol[grepl("neuron|synap|neuropil", paste(reference_markers$marker_set, reference_markers$cell_class), ignore.case = TRUE)])),
  vascular_ecm = unique(c(
    toupper(reference_markers$gene_symbol[grepl("vascular|endothelial|pericyte|ecm", paste(reference_markers$marker_set, reference_markers$cell_class), ignore.case = TRUE)]),
    toupper(microenv_markers$marker_symbol[grepl("vascular|perivascular|basement|ecm|endothelial|pericyte", microenv_markers$panel_id, ignore.case = TRUE)])
  )),
  astrocyte_endfoot = unique(c(
    toupper(reference_markers$gene_symbol[grepl("astro", paste(reference_markers$marker_set, reference_markers$cell_class), ignore.case = TRUE)]),
    toupper(microenv_markers$marker_symbol[grepl("astro|endfoot|gliovascular", microenv_markers$panel_id, ignore.case = TRUE)])
  )),
  oligodendrocyte_myelin = unique(c(
    toupper(reference_markers$gene_symbol[grepl("oligo|myelin", paste(reference_markers$marker_set, reference_markers$cell_class), ignore.case = TRUE)]),
    toupper(microenv_markers$marker_symbol[grepl("oligo|myelin", microenv_markers$panel_id, ignore.case = TRUE)])
  )),
  empirical_microglia_ROI = unique(toupper(empirical_markers$GeneSymbol[grepl("microglia", empirical_markers$marker_set, ignore.case = TRUE)])),
  empirical_neuropil = unique(toupper(empirical_markers$GeneSymbol[grepl("neuropil", empirical_markers$marker_set, ignore.case = TRUE)])),
  empirical_soma = unique(toupper(empirical_markers$GeneSymbol[grepl("soma", empirical_markers$marker_set, ignore.case = TRUE)]))
)
panel_sets <- lapply(panel_sets, function(x) intersect(x[!is.na(x) & nzchar(x)], universe_genes))

marker_rows <- list()
for (mid in module_ids) {
  genes <- intersect(split_genes(module_long$GeneSymbol[as.character(module_long$ModuleID) == mid]), universe_genes)
  for (panel in names(panel_sets)) {
    markers <- panel_sets[[panel]]
    overlap <- intersect(genes, markers)
    p <- stats::phyper(length(overlap) - 1L, length(markers), length(universe_genes) - length(markers), length(genes), lower.tail = FALSE)
    marker_rows[[length(marker_rows) + 1L]] <- data.frame(
      ModuleID = mid, marker_panel = panel, module_gene_count = length(genes),
      panel_gene_count_in_background = length(markers), overlap_count = length(overlap),
      overlap_fraction_module = if (length(genes)) length(overlap) / length(genes) else NA_real_,
      hypergeometric_p = p, overlap_genes = paste(sort(overlap), collapse = ";"),
      interpretation_scope = "ROI/cellular context only; not evidence of cell-intrinsic process",
      stringsAsFactors = FALSE
    )
  }
}
marker_enrichment <- do.call(rbind, marker_rows)
marker_enrichment$FDR <- ave(marker_enrichment$hypergeometric_p, marker_enrichment$marker_panel, FUN = p.adjust, method = "BH")
write_csv_atomic(marker_enrichment, file.path(output_dir, "module_marker_panel_enrichment.csv"))

compartment_summary <- do.call(rbind, lapply(split(active_map, active_map$ModuleID), function(z) {
  usable <- z$joint_primary_matrix_present
  data.frame(
    ModuleID = as.character(z$ModuleID[1]), feature_count = nrow(z),
    n_joint_exact_matches = sum(usable), joint_exact_match_fraction = mean(usable),
    median_microglia_vs_neuropil = median(z$microglia_vs_neuropil[usable], na.rm = TRUE),
    median_microglia_vs_soma = median(z$microglia_vs_soma[usable], na.rm = TRUE),
    fraction_dominant_microglia = mean(z$dominant_compartment[usable] == "microglia", na.rm = TRUE),
    fraction_dominant_neuropil = mean(z$dominant_compartment[usable] == "neuron_neuropil", na.rm = TRUE),
    fraction_dominant_soma = mean(z$dominant_compartment[usable] == "neuron_soma", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
compartment_summary <- compartment_summary[match(module_ids, compartment_summary$ModuleID), ]

context_rows <- lapply(seq_along(module_ids), function(i) {
  mid <- module_ids[i]
  ce <- compartment_summary[compartment_summary$ModuleID == mid, ]
  me <- marker_enrichment[marker_enrichment$ModuleID == mid, ]
  sig <- me$marker_panel[me$FDR <= 0.10 & me$overlap_count >= 2]
  ann <- annotation[i, ]
  empirical_both <- all(c("empirical_microglia_ROI", "empirical_neuropil") %in% sig)
  dominant <- if ("oligodendrocyte_myelin" %in% sig) {
    "oligodendrocyte_myelin-associated_ROI"
  } else if ("vascular_ecm" %in% sig) {
    "vascular_ECM-associated_ROI"
  } else if (all(c("astrocyte_endfoot", "canonical_microglia") %in% sig)) {
    "mixed_astrocyte_microglia_ROI"
  } else if (empirical_both) {
    "shared_microglia_neuropil_ROI"
  } else if ("astrocyte_endfoot" %in% sig) {
    "astrocyte_endfoot-associated_ROI"
  } else if ("canonical_microglia" %in% sig ||
             ("empirical_microglia_ROI" %in% sig && ce$median_microglia_vs_neuropil > 0.25)) {
    "microglia-associated_ROI"
  } else if ("neuronal_neuropil" %in% sig || "empirical_neuropil" %in% sig ||
             ce$median_microglia_vs_neuropil < -0.25) {
    "neuropil-associated_ROI"
  } else {
    "shared_or_unresolved_ROI"
  }
  confidence <- if (ce$joint_exact_match_fraction >= 0.75 && length(sig)) "moderate" else if (ce$joint_exact_match_fraction >= 0.5) "low" else "low"
  limitation <- switch(dominant,
    "microglia-associated_ROI" = "Microglia-associated ROI evidence; not purified-cell or microglia-intrinsic validation.",
    "shared_microglia_neuropil_ROI" = "Shared ROI signal; cannot assign the program exclusively to microglia.",
    "vascular_ECM-associated_ROI" = "Vascular/ECM context may reflect local microenvironment and must not be interpreted as microglia-intrinsic biology.",
    "astrocyte_endfoot-associated_ROI" = "Astrocyte/endfoot context must not be interpreted as microglia-intrinsic.",
    "mixed_astrocyte_microglia_ROI" = "Mixed astrocyte/endfoot and microglia-associated ROI evidence must not be assigned exclusively to microglia.",
    "oligodendrocyte_myelin-associated_ROI" = "Myelin/oligodendrocyte context must not be interpreted as microglia-intrinsic.",
    "neuropil-associated_ROI" = "Neuropil association may reflect adjacent neuronal material rather than microglia-intrinsic biology.",
    "Context remains unresolved; biological-process claims require orthogonal cell-resolved evidence."
  )
  data.frame(
    dataset = "microglia", ModuleID = mid,
    molecular_process = registry_modules$reviewed_biological_label[i],
    subcellular_localization = registry_modules$subcellular_context[i],
    reviewed_roi_context = registry_modules$roi_context[i],
    dominant_context_class = dominant, context_confidence = confidence,
    marker_panels_FDR_le_0_10 = paste(sig, collapse = ";"),
    empirical_microglia_vs_neuropil_support = ce$median_microglia_vs_neuropil,
    joint_exact_match_fraction = ce$joint_exact_match_fraction,
    network_structure_status = module_consensus$worst_sensitivity_status[match(mid, module_consensus$ModuleID)],
    claim_limitation = limitation,
    stringsAsFactors = FALSE
  )
})
context_validation <- do.call(rbind, context_rows)
write_csv_atomic(compartment_summary, file.path(output_dir, "module_compartment_specificity.csv"))
write_csv_atomic(context_validation, file.path(output_dir, "module_biological_context_validation.csv"))

# -------------------------------------------------------------------------
# Part 6: transcriptomic validation readiness only. No concordance is run.
# -------------------------------------------------------------------------

bridge <- module_long[c(
  "ProteinGroupID", "ModuleID", "member_accessions", "GeneSymbol",
  "protein_group_ambiguity_class", "gene_level_claim_allowed", "mapping_status"
)]
bridge$GeneSymbol <- toupper(bridge$GeneSymbol)
gene_counts <- table(bridge$GeneSymbol[nzchar(bridge$GeneSymbol)])
bridge$duplicate_gene_across_protein_groups <- !is.na(bridge$GeneSymbol) & nzchar(bridge$GeneSymbol) & gene_counts[bridge$GeneSymbol] > 1L
bridge$eligible_for_future_gene_set <- bridge$gene_level_claim_allowed %in% TRUE &
  bridge$mapping_status == "mapped" & !is.na(bridge$GeneSymbol) & nzchar(bridge$GeneSymbol) &
  !bridge$duplicate_gene_across_protein_groups
bridge$gene_claim_exclusion_reason <- ifelse(
  bridge$eligible_for_future_gene_set, "eligible_unique_gene_mapping",
  ifelse(!(bridge$gene_level_claim_allowed %in% TRUE), paste0("ambiguous_protein_group:", bridge$protein_group_ambiguity_class),
         ifelse(bridge$mapping_status != "mapped", paste0("mapping_status:", bridge$mapping_status),
                ifelse(bridge$duplicate_gene_across_protein_groups, "duplicate_gene_across_protein_groups", "missing_gene_symbol")))
)
bridge$transcriptomic_background_measured <- FALSE
bridge$transcriptomic_concordance_run <- FALSE
write_csv_atomic(bridge, file.path(output_dir, "module_to_gene_bridge.csv"))

gene_set_audit <- do.call(rbind, lapply(split(bridge, bridge$ModuleID), function(z) {
  data.frame(
    ModuleID = as.character(z$ModuleID[1]), protein_group_count = nrow(z),
    mapped_protein_group_count = sum(z$mapping_status == "mapped", na.rm = TRUE),
    ambiguous_protein_group_count = sum(!(z$gene_level_claim_allowed %in% TRUE)),
    duplicate_gene_protein_group_count = sum(z$duplicate_gene_across_protein_groups, na.rm = TRUE),
    eligible_unique_gene_count = length(unique(z$GeneSymbol[z$eligible_for_future_gene_set])),
    measured_transcriptomic_background_count = NA_integer_,
    region_matching_status = "not_validated_missing_transcriptomic_contract",
    contrast_matching_status = "not_validated_missing_transcriptomic_contract",
    concordance_status = "not_run_invalid_matching_contract",
    stringsAsFactors = FALSE
  )
}))
gene_set_audit <- gene_set_audit[match(module_ids, gene_set_audit$ModuleID), ]
write_csv_atomic(gene_set_audit, file.path(output_dir, "module_transcriptomic_gene_set_readiness.csv"))

transcript_contract <- data.frame(
  contract_element = c("stable_module_identity", "ProteinGroupID_to_gene_mapping", "measured_transcriptomic_background", "region_matching", "contrast_matching", "duplicate_gene_policy", "ambiguous_protein_group_policy", "concordance_execution"),
  status = c("valid", "partially_valid", "missing", "missing", "missing", "defined", "defined", "not_run"),
  evidence = c(
    "13 current ModuleIDs and fixed memberships",
    paste0(sum(bridge$eligible_for_future_gene_set), " of ", nrow(bridge), " ProteinGroupIDs have eligible unique gene mappings"),
    "no measured hippocampal transcriptomic background resolved in current repository inputs",
    "no protein-to-transcript region crosswalk resolved",
    "no exact proteomic-to-transcriptomic contrast crosswalk resolved",
    "duplicate genes excluded until an explicit collapse rule is approved",
    "gene_level_claim_allowed=FALSE protein groups excluded",
    "contract incomplete; no transcriptomic concordance calculated or claimed"
  ),
  required_for_concordance = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  claim_safe = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
write_csv_atomic(transcript_contract, file.path(output_dir, "transcriptomic_matching_contract_audit.csv"))

# -------------------------------------------------------------------------
# Part 7: answer-first recommendations and validation.
# -------------------------------------------------------------------------

status_counts <- table(factor(module_consensus$worst_sensitivity_status, levels = c("pass", "suggestive", "fail")))
robust_n <- unname(status_counts["pass"])
suggestive_n <- unname(status_counts["suggestive"])
fail_n <- unname(status_counts["fail"])
block_primary <- block_robustness[block_robustness$sensitivity_id == "primary", ]
block_all_support <- tapply(block_robustness$cut_height_independent_support, block_robustness$SupermoduleID, all)

recommendations <- data.frame(
  output = c(
    "fixed module architecture and robust co-variation programs",
    "eigengene variance decomposition",
    "sensitivity preservation and leave-one-animal-out tables",
    "higher-order blocks SM01/SM03/SM09",
    "cellular/ROI context validation",
    "group-effect panels",
    "within-animal spatial-adjusted network",
    "transcriptomic bridge readiness"
  ),
  placement = c(
    if (fail_n == 0L) "main-text candidate" else "Extended Data",
    "Extended Data", "Supplementary Data",
    if (all(block_all_support)) "main-text candidate" else "Extended Data",
    "Extended Data", "Extended Data", "diagnostic only", "Supplementary Data"
  ),
  rationale = c(
    paste0(robust_n, " pass, ", suggestive_n, " suggestive, ", fail_n, " fail under the conservative all-sensitivity rule"),
    "descriptive mixed-model decomposition respects nine animal units and exposes spatial dominance",
    "full numerical robustness evidence is necessary for reproducibility but too detailed for the main narrative",
    "only three identities are multi-module blocks; support must persist beyond the descriptive 0.40 cut",
    "distinguishes process from ROI/cellular context and limits intrinsic-cell claims",
    "group contrasts are not independent validation of network construction and retain the existing Stage 05 FDR policy",
    "AnimalID and region are removed by design; strict diagnostic sensitivity cannot replace the primary network",
    "stable bridge is prepared, but the transcriptomic background/region/contrast matching contract is incomplete"
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(recommendations, file.path(output_dir, "output_placement_recommendations.csv"))

context_counts <- table(context_validation$dominant_context_class)
context_text <- paste(paste(names(context_counts), as.integer(context_counts), sep = "="), collapse = "; ")
answers <- data.frame(
  question_id = 1:5,
  question = c(
    "Are the modules robust to repeated-measure and spatial structure?",
    "Which programs are microglia-associated versus neuropil, vascular or shared?",
    "Is WGCNA suitable for a main-text biological claim?",
    "Should the group-effect panels remain Extended Data?",
    "Which single orthogonal validation would add the most value?"
  ),
  answer = c(
    paste0("Conservative fixed-membership result: ", robust_n, " pass, ", suggestive_n, " suggestive, and ", fail_n, " fail. Interpret per-module results; 72 ROIs are not treated as independent animals."),
    paste0("Context classes: ", context_text, ". These are ROI associations, not intrinsic-cell assignments."),
    if (fail_n == 0L) "Yes, for fixed co-variation programs with explicit spatial/context limitations; not as proof of cell-intrinsic mechanism or group causality." else "Only selectively: modules passing robustness and context gates may support a main-text co-variation claim; the network as a whole should not be generalized.",
    "Yes. Group-effect panels should remain Extended Data because they do not validate network construction and retain the existing Stage 05 multiplicity policy.",
    "Matched hippocampal single-nucleus or spatial transcriptomic gene-set validation with an explicit measured background, region crosswalk and identical contrasts would add the most orthogonal value."
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(answers, file.path(output_dir, "primary_questions_answers.csv"))

protected_after <- protected_before
if (nrow(protected_after)) {
  after_paths <- file.path(repo_root, protected_after$protected_path)
  protected_after$md5_after <- unname(tools::md5sum(after_paths))
  protected_after$unchanged <- protected_after$md5_before == protected_after$md5_after
}
write_csv_atomic(protected_after, file.path(output_dir, "protected_output_hash_audit.csv"))

validation <- data.frame(
  check_id = c(
    "current_module_ids", "current_supermodule_ids", "no_membership_duplicates",
    "repeated_measure_unit", "sensitivity_dimensions", "feature_alignment",
    "variance_fractions", "preservation_availability", "loao_coverage",
    "higher_order_membership", "context_claim_limits", "transcriptomic_contract",
    "protected_outputs_unchanged", "no_row_multiplication"
  ),
  passed = c(
    identical(sort(unique(as.character(module_long$ModuleID))), module_ids),
    identical(sort(unique(module_map$SupermoduleID)), super_ids),
    !anyDuplicated(module_long$ProteinGroupID),
    length(unique(metadata$AnimalID)) == 9L && all(table(metadata$AnimalID) == 8L),
    identical(c(nrow(expression_collapsed), nrow(expression_spatial), nrow(expression_strict)), c(36L, 72L, 72L)),
    all(sensitivity_audit$feature_order_exact),
    all(abs(aggregate(variance_fraction ~ level + entity_id, rbind(module_vp, super_vp), sum, na.rm = TRUE)$variance_fraction - 1) < 1e-8),
    nrow(module_robustness) == 39L && all(module_robustness$preservation_metric_reason == "available"),
    nrow(loao) == 9L * 13L && all(table(loao$ModuleID) == 9L),
    all(vapply(names(expected_blocks), function(sm) setequal(observed_members[[sm]], expected_blocks[[sm]]), logical(1))),
    all(grepl("not|cannot|unresolved|require", context_validation$claim_limitation, ignore.case = TRUE)),
    all(transcript_contract$status[transcript_contract$contract_element == "concordance_execution"] == "not_run"),
    all(protected_after$unchanged),
    nrow(context_validation) == 13L && !anyDuplicated(context_validation$ModuleID)
  ),
  detail = c(
    "exact WGCNA_m01..WGCNA_m13", "exact SM01..SM09",
    "one immutable module assignment per active ProteinGroupID",
    "nine animals with eight repeated ROI observations each",
    "36 hemisphere-collapsed; 72 spatial-adjusted; 72 strict residualized",
    "all matrices retain the exact ordered 5201-feature universe",
    "estimable components sum to one per eigengene; non-estimable batch is NA",
    paste0("39 module-matrix rows; ", n_permutations, " preservation permutations"),
    "nine omissions for every module", "SM01/SM03/SM09 match current membership",
    "all module context rows contain explicit claim limitations",
    "concordance not run because matching contract is incomplete",
    "pre/post MD5 hashes for protected paths", "one row per current module"
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(validation, file.path(output_dir, "validation_table.csv"))
if (!all(validation$passed)) stop("One or more audit validation checks failed: ", paste(validation$check_id[!validation$passed], collapse = ", "))

input_hashes <- data.frame(
  input_path = vapply(required_inputs, rel, character(1)),
  bytes = unname(file.info(required_inputs)$size),
  md5 = unname(tools::md5sum(required_inputs)),
  stringsAsFactors = FALSE
)
write_csv_atomic(input_hashes, file.path(output_dir, "input_hashes.csv"))

report_lines <- c(
  "# Microglia WGCNA Nature-readiness audit",
  "",
  "## Answer first",
  "",
  paste0("The current fixed network contains 13 modules from 72 ROI observations nested within nine animals. Under a conservative rule spanning bilateral collapse, spatial adjustment and strict within-animal residualization, **", robust_n, " modules pass, ", suggestive_n, " are suggestive and ", fail_n, " fail**. The detailed module table, not the aggregate count, controls interpretation."),
  "",
  paste0("Only SM01, SM03 and SM09 are evaluated as multi-module higher-order blocks. Cut height 0.40 is descriptive and was not optimized. Cut-height-independent support across all matrices was: ", paste(names(block_all_support), ifelse(block_all_support, "supported", "not consistently supported"), sep = "=", collapse = "; "), "."),
  "",
  "The network can support a main-text **co-variation** claim only for modules and blocks passing the fixed-membership robustness gates. It does not establish microglia-intrinsic mechanism, causal stress effects or independent replication.",
  "",
  "## Repeated-measures and estimability",
  "",
  paste0("Variance fractions were estimated with `", vp_formula_text, "` using `variancePartition::fitExtractVarPartModel`. AnimalID is the biological unit. StressGroup has only three levels and is constant within animals, so its component is descriptive. AcquisitionBatch has two levels and is fully nested in AnimalID; it is non-estimable separately, excluded from the formula and reported as NA rather than assigned variance. Boundary components are flagged. Negative variance estimates, if returned numerically, are truncated to zero and the estimable fractions renormalized."),
  "",
  "## Sensitivity design",
  "",
  "The three sensitivity matrices preserve the exact ordered 5,201 ProteinGroupID universe. Hemisphere collapse produces 36 AnimalID-region observations. Spatial adjustment removes the component of region, hemisphere and acquisition design orthogonal to AnimalID, preserving the complete animal/group linear projection. The strict matrix removes AnimalID and region and is diagnostic only.",
  "",
  "## Biological context",
  "",
  paste0("Same-study jointly normalized proteomics and marker panels classify module context as: ", context_text, ". Molecular-process labels, subcellular localization, ROI/cellular context and network robustness remain separate columns. ROI association is never converted into a microglia-intrinsic claim."),
  "",
  "## Transcriptomic readiness",
  "",
  "A stable ProteinGroupID-to-gene bridge is exported. Transcriptomic concordance was not run: the measured transcriptomic background, region crosswalk and exact contrast crosswalk were not resolved. Duplicate genes and ambiguous protein groups are explicitly excluded pending an approved collapse/matching contract.",
  "",
  "## Publication placement",
  "",
  paste0("- Fixed architecture/programs: ", recommendations$placement[1], "."),
  "- Variance decomposition and cellular-context evidence: Extended Data.",
  "- Preservation, leave-one-animal-out and bridge tables: Supplementary Data.",
  "- Strict within-animal residualization: diagnostic only.",
  "- Existing group-effect panels: remain Extended Data.",
  "",
  "## Highest-value next validation",
  "",
  "Matched hippocampal single-nucleus or spatial transcriptomic gene-set validation, using a measured background and explicit region/contrast crosswalk, would add the most orthogonal value.",
  "",
  "## Reproducibility",
  "",
  paste0("Audit script: `06_modules_WGCNA/12_microglia_wgcna_nature_readiness_audit.R`; preservation permutations: ", n_permutations, "; generated: ", format(Sys.time(), tz = "Europe/Berlin", usetz = TRUE), ". See `input_hashes.csv`, `protected_output_hash_audit.csv`, `validation_table.csv` and `session_info.txt`." )
)
write_lines_atomic(report_lines, file.path(output_dir, "WGCNA_nature_readiness_report.md"))

readme_lines <- c(
  "# Microglia WGCNA Nature-readiness audit outputs",
  "",
  "This directory is additive. It contains fixed-membership sensitivity analyses and does not replace the primary network.",
  "",
  "Key entry points:",
  "",
  "- `WGCNA_nature_readiness_report.html`: portable answer-first report.",
  "- `WGCNA_nature_readiness_report.md`: report source.",
  "- `module_robustness_consensus.csv`: one row per current module.",
  "- `higher_order_block_robustness.csv`: SM01/SM03/SM09 only.",
  "- `module_biological_context_validation.csv`: process/context/claim-limit separation.",
  "- `transcriptomic_matching_contract_audit.csv`: readiness gate; no concordance was run.",
  "- `validation_table.csv`: machine-readable checks.",
  "- `audit_manifest.csv`: output inventory and hashes.",
  "",
  "Rerun from the repository root:",
  "",
  "```powershell",
  "& 'C:\\Users\\topohl\\AppData\\Local\\Programs\\R\\R-4.5.1\\bin\\Rscript.exe' '06_modules_WGCNA/12_microglia_wgcna_nature_readiness_audit.R' --permutations 100",
  "& 'C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\MSBuild\\Microsoft\\VisualStudio\\NodeJs\\node.exe' 'C:\\Users\\topohl\\.codex\\plugins\\cache\\openai-curated-remote\\data-analytics\\0.2.8-13ceeea1f599\\skills\\build-report\\scripts\\deliver_portable_artifact.mjs' --input 'results/reviewer_audit/microglia_wgcna_nature_readiness/WGCNA_nature_readiness_report_artifact.json' --output 'results/reviewer_audit/microglia_wgcna_nature_readiness/WGCNA_nature_readiness_report.html'",
  "& 'C:\\Users\\topohl\\AppData\\Local\\Programs\\R\\R-4.5.1\\bin\\Rscript.exe' '06_modules_WGCNA/12b_finalize_microglia_wgcna_nature_readiness_audit.R'",
  "& 'C:\\Users\\topohl\\AppData\\Local\\Programs\\R\\R-4.5.1\\bin\\Rscript.exe' -e \"testthat::test_file('tests/testthat/test-microglia-wgcna-nature-readiness-audit.R', reporter='summary')\"",
  "```"
)
write_lines_atomic(readme_lines, file.path(output_dir, "README.md"))

session_lines <- capture.output(utils::sessionInfo())
write_lines_atomic(session_lines, file.path(output_dir, "session_info.txt"))

# The canonical report artifact is deliberately compact; detailed evidence is
# linked by safe repository-relative paths and remains in the source-data CSVs.
variance_chart_sql <- paste(
  "SELECT entity_id, component, variance_fraction",
  "FROM module_variance_partition",
  "WHERE estimable = 1 AND variance_fraction IS NOT NULL",
  "ORDER BY component, entity_id",
  sep = "\n"
)
sql_connection <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
on.exit(DBI::dbDisconnect(sql_connection), add = TRUE)
DBI::dbWriteTable(sql_connection, "module_variance_partition", module_vp, overwrite = TRUE)
vp_chart_summary <- DBI::dbGetQuery(sql_connection, variance_chart_sql)
rows_as_lists <- function(x) lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
artifact <- list(
  surface = "report",
  manifest = list(
    version = 1,
    title = "Microglia WGCNA Nature-readiness audit",
    description = "Fixed-membership robustness, repeated-measures, spatial and biological-context audit",
    blocks = list(
      list(id = "title", type = "markdown", body = paste(report_lines[1:3], collapse = "\n")),
      list(id = "answer", type = "markdown", body = paste(report_lines[4:10], collapse = "\n"), sourceId = "audit-summary"),
      list(id = "variance", type = "chart", chartId = "variance_chart"),
      list(id = "caveats", type = "markdown", body = paste(report_lines[11:length(report_lines)], collapse = "\n"), sourceId = "audit-summary")
    ),
    charts = list(
      list(
        id = "variance_chart", title = "Mean module eigengene variance fraction",
        subtitle = "Region dominates most module eigengenes; acquisition batch is excluded as non-estimable",
        type = "bar", dataset = "variance_summary", sourceId = "variance-query",
        encodings = list(
          x = list(field = "component", type = "nominal", label = "Variance component"),
          y = list(field = "variance_fraction", type = "quantitative", aggregate = "avg", format = "percent", label = "Mean variance fraction")
        ),
        xAxisTitle = "Variance component", yAxisTitle = "Mean variance fraction",
        maxRows = 100
      )
    ),
    sources = list(
      list(
        id = "audit-summary", name = "Microglia WGCNA fixed-membership audit",
        type = "file", path = "results/reviewer_audit/microglia_wgcna_nature_readiness",
        description = "Generated by the additive audit script from hashed current-state inputs"
      ),
      list(
        id = "variance-query", name = "Module variance partition chart source",
        type = "query", path = "results/reviewer_audit/microglia_wgcna_nature_readiness/module_eigengene_variance_partition.csv",
        description = "SQLite query executed in-memory against the completed module variance-partition table",
        query = list(
          engine = "SQLite", language = "sql", sql = variance_chart_sql,
          description = "Select all estimable module-level variance fractions for chart aggregation",
          tables_used = list("module_variance_partition"),
          filters = list("estimable = 1", "variance_fraction IS NOT NULL")
        )
      )
    )
  ),
  snapshot = list(
    version = 1,
    status = "ready",
    datasets = list(
      variance_summary = rows_as_lists(vp_chart_summary)
    ),
    access_issues = list()
  ),
  sources = list()
)
artifact_path <- file.path(output_dir, "WGCNA_nature_readiness_report_artifact.json")
atomic_write(artifact_path, function(tmp) jsonlite::write_json(artifact, tmp, auto_unbox = TRUE, pretty = TRUE, na = "null"))

outputs <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
outputs <- outputs[file.exists(outputs) & !dir.exists(outputs)]
outputs <- outputs[basename(outputs) != "audit_manifest.csv"]
manifest <- data.frame(
  relative_path = vapply(outputs, rel, character(1)),
  bytes = unname(file.info(outputs)$size),
  md5 = unname(tools::md5sum(outputs)),
  generated_at = format(Sys.time(), tz = "Europe/Berlin", usetz = TRUE),
  audit_scope = "additive_fixed_membership_microglia_WGCNA",
  stringsAsFactors = FALSE
)
write_csv_atomic(manifest, file.path(output_dir, "audit_manifest.csv"))

message("Audit computation complete: ", rel(output_dir))
