#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/05_module_supermodule_group_effects.r
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/; optional results/tables/06_modules_WGCNA/module_score/<dataset>/.
# Produces: results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv; results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Uses existing WGCNA state/supermodule annotation; safe to run without recomputing WGCNA.
# ================================================================

#
# Test WGCNA module and supermodule eigengene differences between CON, RES, and SUS.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

run <- wgcna_cli()
DATASET <- run$dataset
LEVEL <- match.arg(tolower(run$level), c("module", "supermodule", "both"))
PATHS <- wgcna_downstream_paths("group_effects", DATASET)
FILES <- resolve_wgcna_files(DATASET)
SpatialUnitType <- if (DATASET == "neuron_neuropil") "region_layer" else "region"

if (run$dry_run) {
  invisible(lapply(unlist(PATHS), dir_create))
  dry_run_line("Script", "06_modules_WGCNA/05_module_supermodule_group_effects.r")
  dry_run_line("Dataset", DATASET)
  dry_run_line("Level", LEVEL)
  dry_run_line("Effect scopes", "within_spatial_unit; spatial_adjusted_global; stress_by_spatial_interaction")
  dry_run_line("Mixed models", "Use lmerTest::lmer when AnimalID repeats and lmerTest is installed")
  dry_run_line("WGCNA final state", FILES$state, if (file.exists(FILES$state)) "PASS" else "WARN")
  dry_run_line("Module definitions", FILES$definitions, if (file.exists(FILES$definitions)) "PASS" else "WARN")
  dry_run_line("Supermodule annotation", FILES$supermodule_annotation, if (file.exists(FILES$supermodule_annotation)) "PASS" else "WARN")
  dry_run_line("Optional marker traits", FILES$marker_traits, if (file.exists(FILES$marker_traits)) "PASS" else "WARN")
  dry_run_line("Output tables", PATHS$tables)
  dry_run_line("Supermodule PCA eigenvalue tables", file.path(PATHS$tables, "supermodule_pca_eigenvalues.csv"))
  dry_run_line("Supermodule PCA eigenvalue source data", file.path(PATHS$source_data, "supermodule_pca_eigenvalues.csv"))
  dry_run_line("Supermodule PCA eigenvalue scree SVG", file.path(PATHS$figures, "supermodule_pca_eigenvalue_scree.svg"))
  dry_run_line("Supermodule PCA input audit", file.path(PATHS$tables, "supermodule_pca_input_audit.csv"))
  dry_run_line("Selected SUS-RES PCA eigenvalues", file.path(PATHS$tables, "selected_sus_res_supermodule_pca_eigenvalues.csv"))
  dry_run_line("Selected SUS-RES PCA scree SVG", file.path(PATHS$figures, "selected_sus_res_supermodule_pca_eigenvalue_scree.svg"))
  dry_run_line("Selected SUS-RES PCA selection audit", file.path(PATHS$tables, "selected_sus_res_supermodule_pca_selection_audit.csv"))
  quit(status = 0, save = "no")
}

if (!length(missing_pkgs)) theme_set(theme_classic(base_size = 8))

has_repeats <- function(dat) {
  "AnimalID" %in% names(dat) &&
    any(!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))) &&
    any(table(dat$AnimalID[!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))]) > 1L)
}

collapse_counts <- function(x, group) {
  if (is.null(x) || is.null(group) || !length(x)) return(NA_character_)
  ok <- !is.na(group) & nzchar(as.character(group))
  if (!any(ok)) return(NA_character_)
  counts <- tapply(x[ok], group[ok], function(z) length(unique(stats::na.omit(as.character(z)))))
  if (!length(counts)) return(NA_character_)
  paste(paste0(names(counts), "=", as.integer(counts)), collapse = ";")
}

animal_support_summary <- function(dat, model_type = NA_character_, min_animals_threshold = 3L) {
  if (is.null(dat) || !nrow(dat) || !"AnimalID" %in% names(dat)) {
    return(list(
      n_animals_total = NA_integer_, n_animals_per_group = NA_character_,
      min_animals_per_group = NA_integer_, n_samples_total = if (is.null(dat)) 0L else nrow(dat),
      n_samples_per_group = NA_character_, animal_level_status = "missing_animal_id",
      pseudoreplication_guard = "missing_animal_id",
      biological_replicate_unit = "sample_or_unknown"
    ))
  }
  animal <- as.character(dat$AnimalID)
  usable <- !is.na(animal) & nzchar(animal)
  if (!any(usable)) {
    return(list(
      n_animals_total = NA_integer_, n_animals_per_group = NA_character_,
      min_animals_per_group = NA_integer_, n_samples_total = nrow(dat),
      n_samples_per_group = collapse_counts(seq_len(nrow(dat)), dat$StressGroup),
      animal_level_status = "missing_animal_id",
      pseudoreplication_guard = "missing_animal_id",
      biological_replicate_unit = "sample_or_unknown"
    ))
  }
  dat2 <- dat[usable & !is.na(dat$StressGroup), , drop = FALSE]
  n_total <- dplyr::n_distinct(dat2$AnimalID)
  group_counts <- if (nrow(dat2)) tapply(dat2$AnimalID, dat2$StressGroup, function(x) length(unique(x))) else integer()
  min_group <- if (length(group_counts)) min(as.integer(group_counts), na.rm = TRUE) else NA_integer_
  repeated <- has_repeats(dat2)
  used_mixed <- identical(as.character(model_type), "lmerTest_lmer")
  status <- dplyr::case_when(
    is.na(n_total) || n_total == 0L ~ "missing_animal_id",
    !is.na(min_group) && min_group < min_animals_threshold ~ "insufficient_animals",
    repeated && used_mixed ~ "repeated_sample_mixed_model",
    repeated && !used_mixed ~ "sample_level_or_unclear",
    TRUE ~ "animal_level"
  )
  guard <- dplyr::case_when(
    status %in% c("animal_level", "repeated_sample_mixed_model") ~ "pass",
    status == "insufficient_animals" ~ "insufficient_animals",
    status == "missing_animal_id" ~ "missing_animal_id",
    TRUE ~ "sample_level_or_unclear"
  )
  list(
    n_animals_total = as.integer(n_total),
    n_animals_per_group = collapse_counts(dat2$AnimalID, dat2$StressGroup),
    min_animals_per_group = as.integer(min_group),
    n_samples_total = nrow(dat),
    n_samples_per_group = collapse_counts(seq_len(nrow(dat)), dat$StressGroup),
    animal_level_status = status,
    pseudoreplication_guard = guard,
    biological_replicate_unit = if (status %in% c("animal_level", "repeated_sample_mixed_model")) "animal" else "sample_or_unknown"
  )
}

add_claim_model_fields <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  warn <- tolower(as.character(df$model_warning %||% ""))
  warn[is.na(warn)] <- ""
  model_type <- tolower(as.character(df$model_type %||% ""))
  formula_used <- as.character(df$formula_used %||% NA_character_)
  rank_def <- suppressWarnings(as.logical(df$rank_deficient_model))
  rank_def[is.na(rank_def)] <- FALSE
  fallback <- grepl("fallback|t_test|t-test", model_type) | grepl("fallback|t-test", warn)
  emmeans_ok <- !(grepl("emmeans.*failed|emmeans unavailable|emmeans.*unavailable", warn) | fallback)
  singular <- grepl("singular", warn)
  primary_stable <- !fallback & !rank_def & !singular & emmeans_ok & !is.na(df$p_value)
  downgrade <- character(nrow(df))
  for (i in seq_len(nrow(df))) {
    reasons <- c(
      if (isTRUE(fallback[[i]])) "diagnostic_only_model_fallback" else NULL,
      if (isTRUE(rank_def[[i]])) "rank_deficient_model" else NULL,
      if (isTRUE(singular[[i]])) "singular_model" else NULL,
      if (!isTRUE(emmeans_ok[[i]])) "emmeans_failed_or_unavailable" else NULL,
      if (is.na(df$p_value[[i]])) "missing_model_p_value" else NULL
    )
    downgrade[[i]] <- if (length(reasons)) paste(unique(reasons), collapse = ";") else "none"
  }
  df$model_family <- dplyr::case_when(
    grepl("lmertest_lmer", model_type) ~ "linear_mixed_model",
    grepl("lm", model_type) ~ "linear_model",
    TRUE ~ as.character(df$model_type)
  )
  df$model_formula <- dplyr::coalesce(formula_used, as.character(df$formula_requested))
  df$primary_model_stable <- primary_stable
  df$claim_allowed_model <- primary_stable
  df$model_downgrade_reason <- downgrade
  df$fallback_used <- fallback
  df$fallback_type <- ifelse(fallback, "two_group_t_test_diagnostic", "none")
  df$singular_model <- singular
  df$emmeans_success <- emmeans_ok
  df$animal_random_effect_used <- grepl("lmertest_lmer", model_type)
  df
}

status_row <- function(level, endpoint_id, endpoint_label, spatial_unit, effect_scope, reason,
                       formula_requested = NA_character_, formula_used = NA_character_,
                       model_type = NA_character_, dat = NULL) {
  animal_support <- animal_support_summary(dat, model_type = model_type)
  tibble::tibble(
    dataset = DATASET,
    level = level,
    endpoint_id = endpoint_id,
    endpoint_label = endpoint_label,
    module_id = if (level == "module") endpoint_id else NA_character_,
    supermodule_id = if (level == "supermodule") endpoint_id else NA_character_,
    module_label = if (level == "module") endpoint_label else NA_character_,
    supermodule_label = if (level == "supermodule") endpoint_label else NA_character_,
    spatial_unit = spatial_unit,
    effect_scope = effect_scope,
    SpatialUnitType = SpatialUnitType,
    model_type = model_type,
    has_repeated_animals = if (is.null(dat)) NA else has_repeats(dat),
    n_animals = if (is.null(dat) || !"AnimalID" %in% names(dat)) NA_integer_ else dplyr::n_distinct(dat$AnimalID[!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))]),
    n_animals_total = animal_support$n_animals_total,
    n_animals_per_group = animal_support$n_animals_per_group,
    min_animals_per_group = animal_support$min_animals_per_group,
    n_samples_total = animal_support$n_samples_total,
    n_samples_per_group = animal_support$n_samples_per_group,
    animal_level_status = animal_support$animal_level_status,
    pseudoreplication_guard = animal_support$pseudoreplication_guard,
    contrast = NA_character_,
    estimate = NA_real_,
    SE = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    FDR_within_dataset_level = NA_real_,
    FDR_global = NA_real_,
    evidence_status = "not_supported",
    direction = NA_character_,
    n_samples = if (is.null(dat)) 0L else nrow(dat),
    formula_requested = formula_requested,
    formula_used = formula_used,
    dropped_covariates = NA_character_,
    rank_deficient_model = NA,
    singular_model = grepl("singular", tolower(reason)),
    emmeans_success = FALSE,
    animal_random_effect_used = identical(model_type, "lmerTest_lmer"),
    biological_replicate_unit = animal_support$biological_replicate_unit,
    model_warning = reason
  ) |> add_claim_model_fields()
}

valid_covariates <- function(dat, covars, protect = character()) {
  missing <- setdiff(covars, names(dat))
  present <- intersect(covars, names(dat))
  nonvarying <- present[vapply(present, function(v) dplyr::n_distinct(dat[[v]][!is.na(dat[[v]])]) < 2L, logical(1))]
  dropped <- unique(c(missing, setdiff(nonvarying, protect)))
  list(keep = setdiff(present, dropped), dropped = dropped)
}

fit_model <- function(dat, rhs_terms, random_animal = TRUE) {
  repeated <- has_repeats(dat)
  warning_text <- character()
  use_lmer <- isTRUE(random_animal) && repeated
  if (use_lmer && !requireNamespace("lmerTest", quietly = TRUE)) {
    warning_text <- c(warning_text, "AnimalID repeats detected but lmerTest is unavailable; used lm")
    use_lmer <- FALSE
  }
  rhs <- paste(rhs_terms, collapse = " + ")
  formula_requested <- paste0("eigengene ~ ", rhs, if (use_lmer) " + (1 | AnimalID)" else "")
  fit_formula <- stats::as.formula(formula_requested)
  fit <- tryCatch(
    if (use_lmer) lmerTest::lmer(fit_formula, data = dat, REML = FALSE) else stats::lm(fit_formula, data = dat),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(fit = NULL, warning = c(warning_text, conditionMessage(fit)), model_type = if (use_lmer) "lmerTest_lmer" else "lm", formula_requested = formula_requested))
  rank_deficient <- tryCatch({
    X <- if (use_lmer && requireNamespace("lme4", quietly = TRUE)) lme4::getME(fit, "X") else stats::model.matrix(fit)
    qr(X)$rank < ncol(X)
  }, error = function(e) NA)
  singular <- tryCatch(use_lmer && requireNamespace("lme4", quietly = TRUE) && lme4::isSingular(fit), error = function(e) NA)
  list(
    fit = fit,
    warning = warning_text,
    model_type = if (use_lmer) "lmerTest_lmer" else "lm",
    formula_requested = formula_requested,
    formula_used = paste(deparse(stats::formula(fit)), collapse = ""),
    rank_deficient = rank_deficient,
    singular = singular
  )
}

contrast_rows <- function(fit_info, dat, level, endpoint_id, endpoint_label, spatial_unit, effect_scope, by_spatial = FALSE, dropped = character()) {
  animal_support <- animal_support_summary(dat, model_type = fit_info$model_type)
  if (is.null(fit_info$fit)) {
    return(status_row(level, endpoint_id, endpoint_label, spatial_unit, effect_scope, paste(fit_info$warning, collapse = "; "), fit_info$formula_requested, NA_character_, fit_info$model_type, dat))
  }
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emm <- tryCatch({
      if (by_spatial) emmeans::emmeans(fit_info$fit, specs = stats::as.formula("~ StressGroup | SpatialLabel"))
      else emmeans::emmeans(fit_info$fit, specs = "StressGroup")
    }, error = function(e) e)
    if (!inherits(emm, "error")) {
      contr <- as.data.frame(emmeans::contrast(emm, method = list(
        "RES - CON" = c(-1, 1, 0),
        "SUS - CON" = c(-1, 0, 1),
        "SUS - RES" = c(0, -1, 1)
      )))
      spatial_vals <- if (by_spatial && "SpatialLabel" %in% names(contr)) as.character(contr$SpatialLabel) else spatial_unit
      stat_col <- intersect(c("t.ratio", "z.ratio"), names(contr))[1]
      return(tibble::tibble(
        dataset = DATASET,
        level = level,
        endpoint_id = endpoint_id,
        endpoint_label = endpoint_label,
        module_id = if (level == "module") endpoint_id else NA_character_,
        supermodule_id = if (level == "supermodule") endpoint_id else NA_character_,
        module_label = if (level == "module") endpoint_label else NA_character_,
        supermodule_label = if (level == "supermodule") endpoint_label else NA_character_,
        spatial_unit = spatial_vals,
        effect_scope = effect_scope,
        SpatialUnitType = SpatialUnitType,
        model_type = fit_info$model_type,
        has_repeated_animals = has_repeats(dat),
        n_animals = if ("AnimalID" %in% names(dat)) dplyr::n_distinct(dat$AnimalID[!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))]) else NA_integer_,
        n_animals_total = animal_support$n_animals_total,
        n_animals_per_group = animal_support$n_animals_per_group,
        min_animals_per_group = animal_support$min_animals_per_group,
        n_samples_total = animal_support$n_samples_total,
        n_samples_per_group = animal_support$n_samples_per_group,
        animal_level_status = animal_support$animal_level_status,
        pseudoreplication_guard = animal_support$pseudoreplication_guard,
        contrast = as.character(contr$contrast),
        estimate = as.numeric(contr$estimate),
        SE = as.numeric(contr$SE),
        statistic = if (!is.na(stat_col)) as.numeric(contr[[stat_col]]) else NA_real_,
        p_value = as.numeric(contr$p.value),
        FDR_within_dataset_level = NA_real_,
        FDR_global = NA_real_,
        evidence_status = NA_character_,
        direction = dplyr::case_when(estimate > 0 ~ "higher", estimate < 0 ~ "lower", TRUE ~ "zero"),
        n_samples = nrow(dat),
        formula_requested = fit_info$formula_requested,
        formula_used = fit_info$formula_used,
        dropped_covariates = paste(dropped, collapse = ";"),
        rank_deficient_model = fit_info$rank_deficient,
        singular_model = fit_info$singular,
        emmeans_success = TRUE,
        animal_random_effect_used = identical(fit_info$model_type, "lmerTest_lmer"),
        biological_replicate_unit = animal_support$biological_replicate_unit,
        model_warning = paste(fit_info$warning, collapse = "; ")
      ) |> add_claim_model_fields())
    }
  }

  warning_note <- c(fit_info$warning, "emmeans unavailable or failed; used two-group t-test contrasts")
  contrasts <- c("RES - CON", "SUS - CON", "SUS - RES")
  units <- if (by_spatial) sort(unique(stats::na.omit(dat$SpatialLabel))) else spatial_unit
  dplyr::bind_rows(lapply(units, function(unit) {
    subdat <- if (by_spatial) dat[dat$SpatialLabel == unit, , drop = FALSE] else dat
    fallback_animal_support <- animal_support_summary(subdat, model_type = paste0(fit_info$model_type, "_fallback_t_test"))
    means <- stats::aggregate(subdat$eigengene, list(StressGroup = subdat$StressGroup), mean, na.rm = TRUE)
    names(means)[2] <- "mean"
    dplyr::bind_rows(lapply(contrasts, function(contrast) {
      parts <- strsplit(contrast, " - ", fixed = TRUE)[[1]]
      dd <- subdat[subdat$StressGroup %in% parts, , drop = FALSE]
      tt <- tryCatch(stats::t.test(eigengene ~ StressGroup, data = dd), error = function(e) NULL)
      est <- means$mean[match(parts[[1]], means$StressGroup)] - means$mean[match(parts[[2]], means$StressGroup)]
      tibble::tibble(
        dataset = DATASET, level = level,
        endpoint_id = endpoint_id,
        endpoint_label = endpoint_label,
        module_id = if (level == "module") endpoint_id else NA_character_,
        supermodule_id = if (level == "supermodule") endpoint_id else NA_character_,
        module_label = if (level == "module") endpoint_label else NA_character_,
        supermodule_label = if (level == "supermodule") endpoint_label else NA_character_,
        spatial_unit = unit, effect_scope = effect_scope, SpatialUnitType = SpatialUnitType,
        model_type = paste0(fit_info$model_type, "_fallback_t_test"), has_repeated_animals = has_repeats(subdat),
        n_animals = if ("AnimalID" %in% names(subdat)) dplyr::n_distinct(subdat$AnimalID[!is.na(subdat$AnimalID) & nzchar(as.character(subdat$AnimalID))]) else NA_integer_,
        n_animals_total = fallback_animal_support$n_animals_total,
        n_animals_per_group = fallback_animal_support$n_animals_per_group,
        min_animals_per_group = fallback_animal_support$min_animals_per_group,
        n_samples_total = fallback_animal_support$n_samples_total,
        n_samples_per_group = fallback_animal_support$n_samples_per_group,
        animal_level_status = fallback_animal_support$animal_level_status,
        pseudoreplication_guard = fallback_animal_support$pseudoreplication_guard,
        contrast = contrast, estimate = est,
        SE = if (!is.null(tt)) unname(diff(tt$conf.int)) / (2 * 1.96) else NA_real_,
        statistic = if (!is.null(tt)) unname(tt$statistic) else NA_real_,
        p_value = if (!is.null(tt)) tt$p.value else NA_real_,
        FDR_within_dataset_level = NA_real_, FDR_global = NA_real_,
        evidence_status = NA_character_,
        direction = dplyr::case_when(est > 0 ~ "higher", est < 0 ~ "lower", TRUE ~ "zero"),
        n_samples = nrow(dat), formula_requested = fit_info$formula_requested,
        formula_used = fit_info$formula_used, dropped_covariates = paste(dropped, collapse = ";"),
        rank_deficient_model = fit_info$rank_deficient,
        singular_model = fit_info$singular,
        emmeans_success = FALSE,
        animal_random_effect_used = identical(fit_info$model_type, "lmerTest_lmer"),
        biological_replicate_unit = fallback_animal_support$biological_replicate_unit,
        model_warning = paste(warning_note, collapse = "; ")
      ) |> add_claim_model_fields()
    }))
  }))
}

fit_endpoint_scope <- function(dat_all, endpoint_col, endpoint_id, endpoint_label, level, effect_scope, spatial_unit = NA_character_) {
  dat <- dat_all |>
    dplyr::mutate(eigengene = as.numeric(.data[[endpoint_col]])) |>
    dplyr::filter(is.finite(.data$eigengene), !is.na(.data$StressGroup), .data$StressGroup %in% c("CON", "RES", "SUS"))
  if (!is.na(spatial_unit)) dat <- dat |> dplyr::filter(.data$SpatialLabel == spatial_unit)
  if (nrow(dat) < 4L || dplyr::n_distinct(dat$StressGroup) < 2L) {
    return(status_row(level, endpoint_id, endpoint_label, spatial_unit %||% "global", effect_scope, "too few samples/groups", dat = dat))
  }
  dat$StressGroup <- factor(dat$StressGroup, levels = c("CON", "RES", "SUS"))
  dat$SpatialLabel <- factor(dat$SpatialLabel)

  if (effect_scope == "within_spatial_unit") {
    cov <- valid_covariates(dat, c("StressGroup", "Sex", "Batch"), protect = "StressGroup")
    fit_info <- fit_model(dat, cov$keep, random_animal = TRUE)
    return(contrast_rows(fit_info, dat, level, endpoint_id, endpoint_label, spatial_unit, effect_scope, FALSE, cov$dropped))
  }

  if (effect_scope == "spatial_adjusted_global") {
    cov <- valid_covariates(dat, c("StressGroup", "SpatialLabel", "Sex", "Batch"), protect = c("StressGroup", "SpatialLabel"))
    fit_info <- fit_model(dat, cov$keep, random_animal = TRUE)
    return(contrast_rows(fit_info, dat, level, endpoint_id, endpoint_label, "global_spatial_adjusted", effect_scope, FALSE, cov$dropped))
  }

  if (effect_scope == "stress_by_spatial_interaction") {
    combo <- table(dat$StressGroup, dat$SpatialLabel)
    if (nrow(combo) < 2L || ncol(combo) < 2L || any(combo == 0L)) {
      return(status_row(level, endpoint_id, endpoint_label, "all_spatial_units", effect_scope, "StressGroup x SpatialLabel interaction is not estimable because at least one group/spatial cell is empty.", "eigengene ~ StressGroup * SpatialLabel + Sex + Batch", dat = dat))
    }
    cov <- valid_covariates(dat, c("Sex", "Batch"))
    rhs <- c("StressGroup * SpatialLabel", cov$keep)
    fit_info <- fit_model(dat, rhs, random_animal = TRUE)
    return(contrast_rows(fit_info, dat, level, endpoint_id, endpoint_label, "all_spatial_units", effect_scope, TRUE, cov$dropped))
  }

  status_row(level, endpoint_id, endpoint_label, "unknown", effect_scope, "unknown effect scope", dat = dat)
}

run_effects <- function(eigengenes, endpoint_map, level, state) {
  meta <- standardize_wgcna_metadata(state$sample_info, DATASET)
  dat <- dplyr::inner_join(meta, eigengenes, by = "Sample")
  endpoint_cols <- intersect(endpoint_map$endpoint_col, names(dat))
  spatial_units <- sort(unique(stats::na.omit(dat$SpatialLabel)))
  dplyr::bind_rows(lapply(endpoint_cols, function(col) {
    map <- endpoint_map[match(col, endpoint_map$endpoint_col), , drop = FALSE]
    within <- dplyr::bind_rows(lapply(spatial_units, function(unit) {
      fit_endpoint_scope(dat, col, map$endpoint_id, map$endpoint_label, level, "within_spatial_unit", unit)
    }))
    global <- fit_endpoint_scope(dat, col, map$endpoint_id, map$endpoint_label, level, "spatial_adjusted_global")
    interaction <- fit_endpoint_scope(dat, col, map$endpoint_id, map$endpoint_label, level, "stress_by_spatial_interaction")
    dplyr::bind_rows(within, global, interaction)
  }))
}

make_endpoint_maps <- function(module_eig, definitions, super_ann) {
  module_map <- data.frame(endpoint_col = setdiff(names(module_eig), "Sample"), stringsAsFactors = FALSE) |>
    dplyr::mutate(module_eigengene = .data$endpoint_col, ModuleColor = module_col_to_id(.data$endpoint_col))
  if (nrow(definitions)) {
    defs <- definitions[, intersect(c("ModuleColor", "ModuleID", "module_eigengene", "ModuleLabel_Final"), names(definitions)), drop = FALSE] |> dplyr::distinct()
    module_map <- module_map |> dplyr::left_join(defs, by = "module_eigengene")
  }
  for (nm in c("ModuleID", "ModuleColor.x", "ModuleColor.y", "ModuleLabel_Final")) if (!nm %in% names(module_map)) module_map[[nm]] <- NA_character_
  module_map <- module_map |> dplyr::mutate(endpoint_id = dplyr::coalesce(.data$ModuleID, .data$ModuleColor.x, .data$ModuleColor.y, .data$module_eigengene), endpoint_label = dplyr::coalesce(.data$ModuleLabel_Final, .data$endpoint_id))

  if (nrow(super_ann)) {
    super_map0 <- super_ann
    if ("present_in_dataset" %in% names(super_map0)) super_map0 <- super_map0 |> dplyr::filter(.data$present_in_dataset %in% c(TRUE, "TRUE", "true", 1))
    for (nm in c("Supermodule_DataDrivenID", "Supermodule_DataDriven", "SupermoduleID", "Supermodule_DisplayLabel", "Macroprogram_Display", "Supermodule_LongLabel", "Supermodule_FinalLabel", "Supermodule", "Supermodule_DataDrivenLabel", "supermodule_merge_rule")) if (!nm %in% names(super_map0)) super_map0[[nm]] <- NA_character_
    if (!"supermodule_merge_cut_height" %in% names(super_map0)) super_map0$supermodule_merge_cut_height <- NA_real_
    super_map <- super_map0 |>
      dplyr::mutate(
        module_eigengene = as.character(.data$module_eigengene),
        SupermoduleID = dplyr::coalesce(as.character(.data$Supermodule_DataDrivenID), as.character(.data$Supermodule_DataDriven), as.character(.data$SupermoduleID), as.character(.data$Supermodule)),
        SupermoduleLabel = dplyr::coalesce(as.character(.data$Supermodule_DisplayLabel), as.character(.data$Supermodule_FinalLabel), as.character(.data$Macroprogram_Display), as.character(.data$Supermodule_DataDrivenLabel), as.character(.data$Supermodule), .data$SupermoduleID),
        Supermodule_LongLabel = dplyr::coalesce(as.character(.data$Supermodule_LongLabel), as.character(.data$Supermodule_FinalLabel), as.character(.data$Supermodule)),
        Macroprogram_Display = as.character(.data$Macroprogram_Display)
      )
  } else {
    super_map <- data.frame(module_eigengene = character(), SupermoduleID = character(), SupermoduleLabel = character())
  }
  list(module_map = module_map, super_map = super_map)
}

empty_supermodule_composition <- function() {
  data.frame(
    supermodule_id = character(),
    supermodule_eigengene = character(),
    n_member_modules = integer(),
    member_modules = character(),
    stringsAsFactors = FALSE
  )
}

supermodule_annotation_meta <- function(super_map) {
  for (nm in c("SupermoduleID", "SupermoduleLabel", "supermodule_merge_cut_height", "supermodule_merge_rule")) {
    if (!nm %in% names(super_map)) super_map[[nm]] <- NA
  }
  super_map |>
    dplyr::distinct(
      .data$SupermoduleID,
      .data$SupermoduleLabel,
      .data$supermodule_merge_cut_height,
      .data$supermodule_merge_rule
    )
}

normalize_module_eigengene_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^WGCNA_", "", x, ignore.case = TRUE)
  x <- sub("^ME", "", x, ignore.case = TRUE)
  toupper(x)
}

collapse_examples <- function(x, n = 8L) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  paste(utils::head(x, n), collapse = ";")
}

reconcile_supermodule_map_to_module_eigengenes <- function(super_map, module_eigengenes, definitions = data.frame()) {
  if (is.null(super_map) || !nrow(super_map)) return(super_map)
  me_cols <- setdiff(names(module_eigengenes), "Sample")
  if (!length(me_cols) || !"module_eigengene" %in% names(super_map)) return(super_map)
  super_map$module_eigengene_original <- as.character(super_map$module_eigengene)
  exact <- match(super_map$module_eigengene_original, me_cols)
  norm_cols <- normalize_module_eigengene_key(me_cols)
  norm_map <- match(normalize_module_eigengene_key(super_map$module_eigengene_original), norm_cols)
  label_resolved <- rep(NA_character_, nrow(super_map))
  if (nrow(definitions) && all(c("module_eigengene", "ModuleLabel_Final") %in% names(definitions)) && "ModuleLabel_Final" %in% names(super_map)) {
    label_lookup <- definitions |>
      dplyr::filter(.data$module_eigengene %in% me_cols) |>
      dplyr::transmute(
        module_label_key = trimws(as.character(.data$ModuleLabel_Final)),
        module_eigengene_lookup = as.character(.data$module_eigengene)
      ) |>
      dplyr::filter(nzchar(.data$module_label_key), !is.na(.data$module_label_key)) |>
      dplyr::distinct(.data$module_label_key, .data$module_eigengene_lookup) |>
      dplyr::group_by(.data$module_label_key) |>
      dplyr::filter(dplyr::n_distinct(.data$module_eigengene_lookup) == 1L) |>
      dplyr::ungroup() |>
      dplyr::distinct(.data$module_label_key, .keep_all = TRUE)
    label_hit <- match(trimws(as.character(super_map$ModuleLabel_Final)), label_lookup$module_label_key)
    label_resolved <- label_lookup$module_eigengene_lookup[label_hit]
  }
  resolved <- dplyr::coalesce(
    ifelse(!is.na(exact), me_cols[exact], NA_character_),
    ifelse(!is.na(norm_map), me_cols[norm_map], NA_character_),
    label_resolved
  )
  super_map$module_eigengene_match_method <- dplyr::case_when(
    !is.na(exact) ~ "exact",
    is.na(exact) & !is.na(norm_map) ~ "normalized_ME_prefix",
    is.na(exact) & is.na(norm_map) & !is.na(label_resolved) ~ "module_label_lookup",
    TRUE ~ "unmatched"
  )
  super_map$module_eigengene <- dplyr::coalesce(resolved, super_map$module_eigengene_original)
  super_map
}

supermodule_pca_input_audit <- function(module_eigengenes, super_map, dataset) {
  me_cols <- setdiff(names(module_eigengenes), "Sample")
  if (is.null(super_map)) super_map <- data.frame()
  map_vals <- if ("module_eigengene" %in% names(super_map)) super_map$module_eigengene else character()
  original_vals <- if ("module_eigengene_original" %in% names(super_map)) super_map$module_eigengene_original else map_vals
  tibble::tibble(
    dataset = dataset,
    n_module_eigengene_columns = length(me_cols),
    n_supermodule_map_rows = nrow(super_map),
    n_supermodule_map_rows_matched_to_module_eigengenes = sum(map_vals %in% me_cols, na.rm = TRUE),
    example_module_eigengene_columns = collapse_examples(me_cols),
    example_supermodule_map_module_eigengenes = collapse_examples(original_vals)
  )
}

supermodule_pca_eigenvalues <- function(module_eigengenes, super_map, dataset) {
  me_cols <- setdiff(names(module_eigengenes), "Sample")
  super_map <- super_map[super_map$module_eigengene %in% me_cols & !is.na(super_map$SupermoduleID), , drop = FALSE]
  if (!nrow(super_map)) {
    return(tibble::tibble(
      dataset = character(), supermodule_id = character(), supermodule_label = character(),
      pc = integer(), eigenvalue = numeric(), variance_explained = numeric(),
      cumulative_variance_explained = numeric(), n_member_modules = integer(),
      n_variable_member_modules = integer(), n_samples_used = integer(),
      member_modules = character(), pca_status = character()
    ))
  }

  dplyr::bind_rows(lapply(unique(super_map$SupermoduleID), function(sid) {
    sm <- super_map[super_map$SupermoduleID == sid, , drop = FALSE]
    members <- unique(sm$module_eigengene)
    label <- dplyr::coalesce(na.omit(as.character(sm$SupermoduleLabel))[1], sid)
    vals <- module_eigengenes[, members, drop = FALSE]
    variable <- vapply(vals, function(z) stats::var(as.numeric(z), na.rm = TRUE) > 0, logical(1))
    vals <- vals[, variable, drop = FALSE]
    n_complete <- if (ncol(vals)) sum(stats::complete.cases(vals)) else 0L
    base <- tibble::tibble(
      dataset = dataset,
      supermodule_id = sid,
      supermodule_label = label,
      n_member_modules = length(members),
      n_variable_member_modules = ncol(vals),
      n_samples_used = n_complete,
      member_modules = paste(members, collapse = ";")
    )
    if (!ncol(vals)) {
      return(dplyr::mutate(
        base,
        pc = NA_integer_,
        eigenvalue = NA_real_,
        variance_explained = NA_real_,
        cumulative_variance_explained = NA_real_,
        pca_status = "no_variable_member_modules"
      ))
    }
    if (ncol(vals) == 1L) {
      return(dplyr::mutate(
        base,
        pc = 1L,
        eigenvalue = 1,
        variance_explained = 1,
        cumulative_variance_explained = 1,
        pca_status = "singleton_or_single_variable_module"
      ))
    }
    vals_complete <- vals[stats::complete.cases(vals), , drop = FALSE]
    if (nrow(vals_complete) < 2L) {
      return(dplyr::mutate(
        base,
        pc = NA_integer_,
        eigenvalue = NA_real_,
        variance_explained = NA_real_,
        cumulative_variance_explained = NA_real_,
        pca_status = "too_few_complete_samples"
      ))
    }
    pca <- tryCatch(stats::prcomp(vals_complete, center = TRUE, scale. = TRUE), error = function(e) e)
    if (inherits(pca, "error")) {
      return(dplyr::mutate(
        base,
        pc = NA_integer_,
        eigenvalue = NA_real_,
        variance_explained = NA_real_,
        cumulative_variance_explained = NA_real_,
        pca_status = paste0("pca_failed: ", conditionMessage(pca))
      ))
    }
    eigenvalues <- pca$sdev^2
    variance <- eigenvalues / sum(eigenvalues)
    dplyr::bind_cols(
      base[rep(1L, length(eigenvalues)), , drop = FALSE],
      tibble::tibble(
        pc = seq_along(eigenvalues),
        eigenvalue = eigenvalues,
        variance_explained = variance,
        cumulative_variance_explained = cumsum(variance),
        pca_status = "ok"
      )
    )
  }))
}

plot_supermodule_pca_scree <- function(pca_eigenvalues, path) {
  plot_df <- pca_eigenvalues |>
    dplyr::filter(!is.na(.data$pc), is.finite(.data$variance_explained)) |>
    dplyr::mutate(
      pc_label = paste0("PC", .data$pc),
      facet_label = paste0(.data$supermodule_id, ": ", .data$supermodule_label)
    )
  if (!nrow(plot_df)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No variable supermodule PCA components available\nSee supermodule_pca_input_audit.csv", size = 3) +
      ggplot2::labs(title = "Supermodule PCA eigenvalue scree") +
      ggplot2::theme_void()
  } else {
    line_df <- plot_df |>
      dplyr::group_by(.data$facet_label) |>
      dplyr::filter(dplyr::n() > 1L) |>
      dplyr::ungroup()
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$pc_label)) +
      ggplot2::geom_col(ggplot2::aes(y = .data$variance_explained), fill = "#4E79A7", width = 0.72)
    if (nrow(line_df)) {
      p <- p +
        ggplot2::geom_line(
          data = line_df,
          ggplot2::aes(y = .data$cumulative_variance_explained, group = 1),
          color = "#D95F02",
          linewidth = 0.35
        )
    }
    p <- p +
      ggplot2::geom_point(ggplot2::aes(y = .data$cumulative_variance_explained), color = "#D95F02", size = 0.75) +
      ggplot2::facet_wrap(~facet_label, scales = "free_x") +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
      ggplot2::labs(x = NULL, y = "Variance explained", title = "Supermodule PCA eigenvalue scree") +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5),
        strip.text = ggplot2::element_text(size = 6)
      )
  }
  ggplot2::ggsave(path, p, width = 210, height = 180, units = "mm", device = svglite::svglite)
}

standardize_contrast_key <- function(x) {
  toupper(gsub("[[:space:]]+", "", gsub("–", "-", as.character(x))))
}

select_sus_res_supermodules <- function(super_effects, dataset, max_n = 2L) {
  no_selection_row <- function(message) tibble::tibble(
    dataset = dataset,
    supermodule_id = NA_character_,
    supermodule_label = NA_character_,
    contrast = "SUS-RES",
    estimate = NA_real_,
    p_value = NA_real_,
    FDR_within_dataset_level = NA_real_,
    FDR_global = NA_real_,
    selection_support = "none",
    selection_rank = NA_integer_,
    selection_message = message
  )
  if (is.null(super_effects) || !nrow(super_effects)) {
    return(no_selection_row("no supermodule group-effect rows available"))
  }
  sus_res <- super_effects |>
    dplyr::filter(standardize_contrast_key(.data$contrast) == "SUS-RES") |>
    dplyr::mutate(
      p_value = suppressWarnings(as.numeric(.data$p_value)),
      FDR_within_dataset_level = suppressWarnings(as.numeric(.data$FDR_within_dataset_level)),
      FDR_global = suppressWarnings(as.numeric(.data$FDR_global)),
      estimate = suppressWarnings(as.numeric(.data$estimate)),
      support_rank = dplyr::case_when(
        !is.na(.data$FDR_within_dataset_level) & .data$FDR_within_dataset_level <= 0.05 ~ 1L,
        !is.na(.data$FDR_within_dataset_level) & .data$FDR_within_dataset_level <= 0.10 ~ 2L,
        !is.na(.data$p_value) & .data$p_value < 0.05 ~ 3L,
        TRUE ~ 99L
      ),
      selection_support = dplyr::case_when(
        .data$support_rank == 1L ~ "FDR_within_dataset_level <= 0.05",
        .data$support_rank == 2L ~ "FDR_within_dataset_level <= 0.10",
        .data$support_rank == 3L ~ "nominal p_value < 0.05",
        TRUE ~ "not_selected"
      ),
      sort_fdr = dplyr::coalesce(.data$FDR_within_dataset_level, .data$FDR_global, Inf),
      sort_p = dplyr::coalesce(.data$p_value, Inf),
      supermodule_label = dplyr::coalesce(.data$supermodule_label, .data$endpoint_label, .data$supermodule_id)
    ) |>
    dplyr::filter(.data$support_rank < 99L) |>
    dplyr::arrange(.data$support_rank, .data$sort_fdr, .data$sort_p, dplyr::desc(abs(.data$estimate))) |>
    dplyr::group_by(.data$supermodule_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::slice_head(n = max_n) |>
    dplyr::mutate(selection_rank = dplyr::row_number(), selection_message = "selected") |>
    dplyr::select(
      "dataset", "supermodule_id", "supermodule_label", "contrast", "estimate", "p_value",
      "FDR_within_dataset_level", "FDR_global", "selection_support",
      "selection_rank", "selection_message"
    )
  if (!nrow(sus_res)) {
    return(no_selection_row("No SUS-RES significant supermodule eigengenes selected"))
  }
  sus_res
}

selected_sus_res_pca_eigenvalues <- function(pca_eigenvalues, selected) {
  if (is.null(selected) || !nrow(selected) || !"supermodule_id" %in% names(selected)) {
    return(pca_eigenvalues[0, , drop = FALSE])
  }
  selected <- selected |> dplyr::filter(.data$selection_message == "selected")
  if (!nrow(selected)) {
    return(pca_eigenvalues[0, , drop = FALSE])
  }
  pca_eigenvalues |>
    dplyr::inner_join(
      selected |>
        dplyr::select(
          "supermodule_id",
          selected_supermodule_label = "supermodule_label",
          sus_res_estimate = "estimate",
          sus_res_p_value = "p_value",
          sus_res_FDR_within_dataset_level = "FDR_within_dataset_level",
          sus_res_FDR_global = "FDR_global",
          "selection_support",
          "selection_rank"
        ),
      by = "supermodule_id"
    ) |>
    dplyr::arrange(.data$selection_rank, .data$pc)
}

plot_selected_sus_res_pca_scree <- function(selected_pca, selected, path) {
  if (!nrow(selected) || all(selected$selection_message != "selected")) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No SUS-RES significant supermodule eigengenes selected", size = 3) +
      ggplot2::labs(title = "Selected SUS-RES supermodule PCA eigenvalue scree") +
      ggplot2::theme_void()
  } else {
    plot_df <- selected_pca |>
      dplyr::filter(!is.na(.data$pc), is.finite(.data$variance_explained)) |>
      dplyr::mutate(
        pc_label = paste0("PC", .data$pc),
        facet_label = paste0(
          .data$supermodule_id, ": ", .data$supermodule_label,
          "\nSUS-RES est=", signif(.data$sus_res_estimate, 3),
          ", p=", signif(.data$sus_res_p_value, 3),
          ", FDR=", signif(.data$sus_res_FDR_within_dataset_level, 3)
        )
      )
    if (!nrow(plot_df)) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No PCA eigenvalues available for selected SUS-RES supermodules", size = 3) +
        ggplot2::labs(title = "Selected SUS-RES supermodule PCA eigenvalue scree") +
        ggplot2::theme_void()
    } else {
      line_df <- plot_df |>
        dplyr::group_by(.data$facet_label) |>
        dplyr::filter(dplyr::n() > 1L) |>
        dplyr::ungroup()
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$pc_label)) +
        ggplot2::geom_col(ggplot2::aes(y = .data$variance_explained), fill = "#4E79A7", width = 0.72)
      if (nrow(line_df)) {
        p <- p +
          ggplot2::geom_line(
            data = line_df,
            ggplot2::aes(y = .data$cumulative_variance_explained, group = 1),
            color = "#D95F02",
            linewidth = 0.35
          )
      }
      p <- p +
        ggplot2::geom_point(ggplot2::aes(y = .data$cumulative_variance_explained), color = "#D95F02", size = 0.75) +
        ggplot2::facet_wrap(~facet_label, scales = "free_x") +
        ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
        ggplot2::labs(
          x = NULL,
          y = "Variance explained",
          title = "Selected SUS-RES supermodule PCA eigenvalue scree",
          subtitle = "Selection: up to two SUS-RES supermodule effects prioritized by within-dataset FDR, nominal p-value, and absolute estimate"
        ) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6),
          strip.text = ggplot2::element_text(size = 6)
        )
    }
  }
  ggplot2::ggsave(path, p, width = 210, height = 180, units = "mm", device = svglite::svglite)
}

correlate_marker_traits <- function(eigengenes, endpoint_map, marker_traits, level) {
  if (is.null(marker_traits) || !nrow(marker_traits)) return(tibble::tibble())
  trait_cols <- grep("^(z|raw)_.+_(score|ratio)$|^z_microglia_minus_neuropil_score$|^z_microglia_to_neuropil_ratio$|^raw_microglia_minus_neuropil_score$|^raw_microglia_to_neuropil_ratio$", names(marker_traits), value = TRUE)
  trait_cols <- trait_cols[vapply(marker_traits[trait_cols], is.numeric, logical(1))]
  dat <- dplyr::inner_join(eigengenes, marker_traits[, c("Sample", trait_cols), drop = FALSE], by = "Sample")
  endpoint_cols <- intersect(endpoint_map$endpoint_col, names(dat))
  if (!length(endpoint_cols) || !length(trait_cols)) return(tibble::tibble(
    dataset = character(), level = character(), module_id = character(), supermodule_id = character(),
    marker_trait = character(), correlation = numeric(), p_value = numeric(), FDR = numeric(),
    n_samples = integer(), interpretation_note = character()
  ))
  rows <- lapply(endpoint_cols, function(ec) {
    map <- endpoint_map[match(ec, endpoint_map$endpoint_col), , drop = FALSE]
    dplyr::bind_rows(lapply(trait_cols, function(tc) {
      ok <- is.finite(dat[[ec]]) & is.finite(dat[[tc]])
      ct <- if (sum(ok) >= 4L) tryCatch(stats::cor.test(dat[[ec]][ok], dat[[tc]][ok], method = "spearman"), error = function(e) NULL) else NULL
      tibble::tibble(
        dataset = DATASET, level = level,
        module_id = if (level == "module") map$endpoint_id else NA_character_,
        supermodule_id = if (level == "supermodule") map$endpoint_id else NA_character_,
        marker_trait = tc,
        correlation = if (!is.null(ct)) unname(ct$estimate) else NA_real_,
        p_value = if (!is.null(ct)) ct$p.value else NA_real_,
        FDR = NA_real_,
        n_samples = sum(ok),
        interpretation_note = "Marker traits are annotation/sensitivity layers only; not default covariates."
      )
    }))
  })
  out <- dplyr::bind_rows(rows)
  out$FDR <- stats::p.adjust(out$p_value, method = "BH")
  out
}

classify_group_effect_evidence <- function(df) {
  warn <- if ("model_warning" %in% names(df)) tolower(as.character(df$model_warning)) else rep("", nrow(df))
  warn[is.na(warn)] <- ""
  rank_def <- if ("rank_deficient_model" %in% names(df)) suppressWarnings(as.logical(df$rank_deficient_model)) else rep(FALSE, nrow(df))
  rank_def[is.na(rank_def)] <- FALSE
  unstable <- rank_def |
    grepl("rank|singular|not estimable|failed|unavailable|t-test|too few|empty", warn) |
    grepl("fallback_t_test", tolower(as.character(df$model_type %||% "")))
  dplyr::case_when(
    is.na(df$p_value) ~ "not_supported",
    unstable ~ "model_unstable",
    !is.na(df$FDR_global) & df$FDR_global <= 0.05 ~ "robust_FDR",
    !is.na(df$FDR_global) & df$FDR_global <= 0.10 ~ "suggestive_FDR10",
    !is.na(df$FDR_within_dataset_level) & df$FDR_within_dataset_level <= 0.05 ~ "robust_FDR",
    !is.na(df$FDR_within_dataset_level) & df$FDR_within_dataset_level <= 0.10 ~ "suggestive_FDR10",
    !is.na(df$p_value) & df$p_value < 0.05 ~ "nominal_only",
    TRUE ~ "not_supported"
  )
}

rank_group_effects <- function(df) {
  for (nm in required_group_effect_columns) if (!nm %in% names(df)) df[[nm]] <- NA
  out <- df |>
    dplyr::mutate(
      endpoint_id = dplyr::coalesce(.data$endpoint_id, ifelse(.data$level == "module", .data$module_id, .data$supermodule_id)),
      endpoint_label = dplyr::coalesce(.data$endpoint_label, ifelse(.data$level == "module", .data$module_label, .data$supermodule_label))
    )
  out$evidence_status <- classify_group_effect_evidence(out)
  out |>
    dplyr::mutate(
      evidence_rank = dplyr::case_when(
        .data$evidence_status == "robust_FDR" ~ 1L,
        .data$evidence_status == "suggestive_FDR10" ~ 2L,
        .data$evidence_status == "nominal_only" ~ 3L,
        .data$evidence_status == "model_unstable" ~ 4L,
        TRUE ~ 5L
      )
    ) |>
    dplyr::arrange(.data$evidence_rank, .data$FDR_global, .data$FDR_within_dataset_level, .data$p_value, dplyr::desc(abs(.data$estimate))) |>
    dplyr::select(-"evidence_rank")
}

state <- tryCatch(load_wgcna_state(FILES$state), error = function(e) e)
if (inherits(state, "error")) {
  msg <- conditionMessage(state)
  module_out <- empty_group_effects(DATASET, "module", msg)
  super_out <- empty_group_effects(DATASET, "supermodule", msg)
  comp <- data.frame(dataset = DATASET, supermodule_id = NA_character_, n_member_modules = 0L, member_modules = NA_character_, stringsAsFactors = FALSE)
  module_marker <- tibble::tibble()
  super_marker <- tibble::tibble()
  super_pca_input_audit <- supermodule_pca_input_audit(data.frame(Sample = character()), data.frame(), DATASET)
  super_pca_eigenvalues <- supermodule_pca_eigenvalues(data.frame(Sample = character()), data.frame(), DATASET)
} else {
  module_eig <- extract_module_eigengenes(state)
  definitions <- safe_read_csv(FILES$definitions) %||% data.frame()
  super_ann <- safe_read_csv(FILES$supermodule_annotation) %||% data.frame()
  maps <- make_endpoint_maps(module_eig, definitions, super_ann)
  maps$super_map <- reconcile_supermodule_map_to_module_eigengenes(maps$super_map, module_eig, definitions)
  super_pca_input_audit <- supermodule_pca_input_audit(module_eig, maps$super_map, DATASET)
  super_pca_eigenvalues <- supermodule_pca_eigenvalues(module_eig, maps$super_map, DATASET)
  super <- make_supermodule_eigengenes(module_eig, maps$super_map)
  comp0 <- super$composition
  if (is.null(comp0) || !"supermodule_id" %in% names(comp0)) comp0 <- empty_supermodule_composition()
  comp <- comp0 |>
    dplyr::mutate(dataset = DATASET, .before = "supermodule_id") |>
    dplyr::left_join(
      supermodule_annotation_meta(maps$super_map),
      by = c("supermodule_id" = "SupermoduleID")
    ) |>
    dplyr::rename(supermodule_label = "SupermoduleLabel")

  module_out <- if (LEVEL %in% c("module", "both")) run_effects(module_eig, maps$module_map, "module", state) else empty_group_effects(DATASET, "module", "not requested")
  super_endpoint_map <- comp |> dplyr::transmute(endpoint_col = .data$supermodule_eigengene, endpoint_id = .data$supermodule_id, endpoint_label = .data$supermodule_label)
  super_out <- if (LEVEL %in% c("supermodule", "both") && nrow(super_endpoint_map)) run_effects(super$eigengenes, super_endpoint_map, "supermodule", state) else empty_group_effects(DATASET, "supermodule", "not requested or no supermodules")

  marker_traits <- safe_read_csv(FILES$marker_traits)
  module_marker <- correlate_marker_traits(module_eig, maps$module_map, marker_traits, "module")
  super_marker <- correlate_marker_traits(super$eigengenes, super_endpoint_map, marker_traits, "supermodule")
}

all_p <- c(module_out$p_value, super_out$p_value)
all_fdr <- stats::p.adjust(all_p, method = "BH")
module_out$FDR_global <- all_fdr[seq_len(nrow(module_out))]
super_out$FDR_global <- all_fdr[seq_len(nrow(super_out)) + nrow(module_out)]
module_out$FDR_within_dataset_level <- stats::p.adjust(module_out$p_value, method = "BH")
super_out$FDR_within_dataset_level <- stats::p.adjust(super_out$p_value, method = "BH")
module_out <- rank_group_effects(module_out)
super_out <- rank_group_effects(super_out)
module_out <- module_out[, required_group_effect_columns]
super_out <- super_out[, required_group_effect_columns]
selected_sus_res_audit <- select_sus_res_supermodules(super_out, DATASET, max_n = 2L)
selected_sus_res_pca <- selected_sus_res_pca_eigenvalues(super_pca_eigenvalues, selected_sus_res_audit)

write_table_and_source(module_out, PATHS$tables, PATHS$source_data, "module_group_effects.csv")
write_table_and_source(super_out, PATHS$tables, PATHS$source_data, "supermodule_group_effects.csv")
write_table_and_source(comp, PATHS$tables, PATHS$source_data, "supermodule_composition.csv")
write_table_and_source(super_pca_input_audit, PATHS$tables, PATHS$source_data, "supermodule_pca_input_audit.csv")
write_table_and_source(super_pca_eigenvalues, PATHS$tables, PATHS$source_data, "supermodule_pca_eigenvalues.csv")
write_table_and_source(selected_sus_res_audit, PATHS$tables, PATHS$source_data, "selected_sus_res_supermodule_pca_selection_audit.csv")
write_table_and_source(selected_sus_res_pca, PATHS$tables, PATHS$source_data, "selected_sus_res_supermodule_pca_eigenvalues.csv")
write_table_and_source(module_marker, PATHS$tables, PATHS$source_data, "module_marker_trait_correlations.csv")
write_table_and_source(super_marker, PATHS$tables, PATHS$source_data, "supermodule_marker_trait_correlations.csv")

map_out <- safe_read_csv(FILES$supermodule_annotation)
if (is.null(map_out)) map_out <- data.frame(dataset = DATASET, status = "missing_supermodule_annotation")
write_table_and_source(map_out, PATHS$tables, PATHS$source_data, "module_to_supermodule_map_with_annotations.csv")

plot_effects <- function(df, path, title) {
  plot_df <- df |> dplyr::filter(!is.na(.data$p_value))
  if (!nrow(plot_df)) return(invisible(NULL))
  label_col <- if (unique(plot_df$level)[[1]] == "module") "module_id" else "supermodule_id"
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data[[label_col]], color = .data$estimate, size = -log10(.data$p_value))) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_grid(effect_scope ~ spatial_unit, scales = "free_y", space = "free_y") +
    ggplot2::scale_color_gradient2(low = "#3B6FB6", mid = "grey92", high = "#C84C5A") +
    ggplot2::labs(x = NULL, y = NULL, color = "Estimate", size = "-log10 P", title = title) +
    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(path, p, width = 190, height = 140, units = "mm", device = svglite::svglite)
}
plot_effects(module_out, file.path(PATHS$figures, "module_group_effect_dotplot.svg"), "Module group effects")
plot_effects(super_out, file.path(PATHS$figures, "supermodule_effects_by_spatial_unit.svg"), "Supermodule group effects")
plot_effects(super_out, file.path(PATHS$figures, "supermodule_group_effect_dotplot.svg"), "Supermodule group effects")
plot_supermodule_pca_scree(super_pca_eigenvalues, file.path(PATHS$figures, "supermodule_pca_eigenvalue_scree.svg"))
plot_selected_sus_res_pca_scree(selected_sus_res_pca, selected_sus_res_audit, file.path(PATHS$figures, "selected_sus_res_supermodule_pca_eigenvalue_scree.svg"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = FILES,
  outputs = list(tables = PATHS$tables, source_data = PATHS$source_data, figures = PATHS$figures),
  parameters = list(dataset = DATASET, level = LEVEL, effect_scopes = c("within_spatial_unit", "spatial_adjusted_global", "stress_by_spatial_interaction")),
  notes = "Primary results are not filtered by marker or microenvironment class. Marker traits, if present, are annotation only."
)

message("WGCNA module/supermodule group effects complete for dataset: ", DATASET)
