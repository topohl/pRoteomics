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
  dry_run_line("Supermodule PCA member loadings", file.path(PATHS$tables, "supermodule_pca_member_loadings.csv"))
  dry_run_line("Supermodule PCA input audit", file.path(PATHS$tables, "supermodule_pca_input_audit.csv"))
  dry_run_line("Selected SUS-RES PCA eigenvalues", file.path(PATHS$tables, "selected_sus_res_supermodule_pca_eigenvalues.csv"))
  dry_run_line("Selected SUS-RES PCA scree SVG", file.path(PATHS$figures, "selected_sus_res_supermodule_pca_eigenvalue_scree.svg"))
  dry_run_line("Selected SUS-RES PCA selection audit", file.path(PATHS$tables, "selected_sus_res_supermodule_pca_selection_audit.csv"))
  dry_run_line("Selected SUS-RES supermodule contents", file.path(PATHS$tables, "selected_sus_res_supermodule_contents.csv"))
  dry_run_line("Selected SUS-RES supermodule interpretation", file.path(PATHS$tables, "selected_sus_res_supermodule_interpretation.csv"))
  dry_run_line("Selected SUS-RES interpretation summary", file.path(PATHS$tables, "selected_sus_res_supermodule_interpretation_summary.csv"))
  dry_run_line("Selected SUS-RES member loading SVG", file.path(PATHS$figures, "selected_sus_res_supermodule_member_loading_plot.svg"))
  dry_run_line("All supermodule eigengene group values", file.path(PATHS$tables, "all_supermodule_eigengene_group_values.csv"))
  dry_run_line("All supermodule eigengene group SVG", file.path(PATHS$figures, "all_supermodule_eigengene_group_plot.svg"))
  dry_run_line("All supermodule eigengene spatial group values", file.path(PATHS$tables, "all_supermodule_eigengene_spatial_group_values.csv"))
  dry_run_line("All supermodule eigengene spatial group SVG", file.path(PATHS$figures, "all_supermodule_eigengene_spatial_group_plot.svg"))
  dry_run_line("All supermodule member loading SVG", file.path(PATHS$figures, "all_supermodule_member_loading_plot.svg"))
  dry_run_line("All supermodule contents summary", file.path(PATHS$tables, "all_supermodule_contents_summary.csv"))
  dry_run_line("All supermodule region/layer effects", file.path(PATHS$tables, "all_supermodule_region_layer_effects.csv"))
  dry_run_line("All supermodule region/layer Cohen's d effects", file.path(PATHS$tables, "all_supermodule_region_layer_cohend_effects.csv"))
  dry_run_line("All supermodule region/layer heatmap SVG", file.path(PATHS$figures, "all_supermodule_region_layer_effect_heatmap.svg"))
  dry_run_line("All SUS-RES region/layer effects SVG", file.path(PATHS$figures, "all_sus_res_supermodule_region_layer_effects.svg"))
  dry_run_line("All SUS-RES region/layer effects source", file.path(PATHS$source_data, "all_sus_res_supermodule_region_layer_effects.csv"))
  dry_run_line("Selected SUS-RES effect summary SVG", file.path(PATHS$figures, "selected_sus_res_supermodule_effect_summary.svg"))
  dry_run_line("Selected SUS-RES contents overview SVG", file.path(PATHS$figures, "selected_sus_res_supermodule_contents_overview.svg"))
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

collapse_unique_values <- function(x, max_n = 12L) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  if (length(x) > max_n) return(paste(c(utils::head(x, max_n), "..."), collapse = "; "))
  paste(x, collapse = "; ")
}

col_if_present <- function(df, nm, default = NA_character_) {
  if (!is.null(df) && nm %in% names(df)) return(df[[nm]])
  rep(default, if (is.null(df)) 0L else nrow(df))
}

first_value <- function(x, default = NA) {
  if (is.null(x) || !length(x)) return(default)
  x[[1]]
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

empty_supermodule_eigengene_group_values <- function() {
  tibble::tibble(
    dataset = character(),
    Sample = character(),
    supermodule_id = character(),
    supermodule_label = character(),
    supermodule_eigengene = character(),
    eigengene = numeric(),
    eigengene_z = numeric(),
    StressGroup = character(),
    AnimalID = character(),
    Sex = character(),
    Batch = character(),
    SpatialLabel = character(),
    SpatialUnitType = character()
  )
}

all_supermodule_eigengene_group_values <- function(super_eigengenes, super_endpoint_map, state, dataset, spatial_unit_type) {
  if (is.null(super_eigengenes) || !nrow(super_eigengenes) || is.null(super_endpoint_map) || !nrow(super_endpoint_map)) {
    return(empty_supermodule_eigengene_group_values())
  }
  meta <- standardize_wgcna_metadata(state$sample_info, dataset)
  for (nm in c("AnimalID", "Sex", "Batch", "SpatialLabel", "StressGroup")) {
    if (!nm %in% names(meta)) meta[[nm]] <- NA_character_
  }
  endpoint_map <- super_endpoint_map |>
    dplyr::filter(.data$endpoint_col %in% names(super_eigengenes)) |>
    dplyr::transmute(
      supermodule_eigengene = as.character(.data$endpoint_col),
      supermodule_id = as.character(.data$endpoint_id),
      supermodule_label = as.character(dplyr::coalesce(.data$endpoint_label, .data$endpoint_id))
    ) |>
    dplyr::distinct(.data$supermodule_eigengene, .keep_all = TRUE)
  if (!nrow(endpoint_map)) return(empty_supermodule_eigengene_group_values())

  joined <- dplyr::inner_join(
    meta |>
      dplyr::select("Sample", "StressGroup", "AnimalID", "Sex", "Batch", "SpatialLabel"),
    super_eigengenes |>
      dplyr::select("Sample", dplyr::all_of(endpoint_map$supermodule_eigengene)),
    by = "Sample"
  )
  if (!nrow(joined)) return(empty_supermodule_eigengene_group_values())

  joined |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(endpoint_map$supermodule_eigengene),
      names_to = "supermodule_eigengene",
      values_to = "eigengene"
    ) |>
    dplyr::left_join(endpoint_map, by = "supermodule_eigengene") |>
    dplyr::group_by(.data$supermodule_id) |>
    dplyr::mutate(
      eigengene_z = {
        mu <- mean(.data$eigengene, na.rm = TRUE)
        sig <- stats::sd(.data$eigengene, na.rm = TRUE)
        if (is.finite(sig) && sig > 0) (.data$eigengene - mu) / sig else NA_real_
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      dataset = dataset,
      StressGroup = factor(as.character(.data$StressGroup), levels = c("CON", "RES", "SUS")),
      SpatialUnitType = spatial_unit_type
    ) |>
    dplyr::arrange(.data$supermodule_id, .data$Sample) |>
    dplyr::select(
      "dataset", "Sample", "supermodule_id", "supermodule_label", "supermodule_eigengene",
      "eigengene", "eigengene_z", "StressGroup", "AnimalID", "Sex", "Batch", "SpatialLabel", "SpatialUnitType"
    )
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

supermodule_pca_member_loadings <- function(module_eigengenes, super_map, dataset) {
  me_cols <- setdiff(names(module_eigengenes), "Sample")
  super_map <- super_map[super_map$module_eigengene %in% me_cols & !is.na(super_map$SupermoduleID), , drop = FALSE]
  empty <- tibble::tibble(
    dataset = character(), supermodule_id = character(), supermodule_label = character(),
    module_eigengene = character(), module_eigengene_original = character(),
    module_id = character(), module_label = character(), pc = integer(),
    loading = numeric(), abs_loading = numeric(), loading_rank = integer(),
    n_member_modules = integer(), n_variable_member_modules = integer(),
    n_samples_used = integer(), pca_status = character()
  )
  if (!nrow(super_map)) return(empty)

  dplyr::bind_rows(lapply(unique(super_map$SupermoduleID), function(sid) {
    sm <- super_map[super_map$SupermoduleID == sid, , drop = FALSE]
    members <- unique(sm$module_eigengene)
    label <- dplyr::coalesce(na.omit(as.character(sm$SupermoduleLabel))[1], sid)
    vals <- module_eigengenes[, members, drop = FALSE]
    variable <- vapply(vals, function(z) stats::var(as.numeric(z), na.rm = TRUE) > 0, logical(1))
    variable_members <- names(variable)[variable]
    n_complete <- if (length(variable_members)) sum(stats::complete.cases(vals[, variable_members, drop = FALSE])) else 0L
    module_meta_src <- sm |>
      dplyr::distinct(.data$module_eigengene, .keep_all = TRUE)
    for (nm in c("module_eigengene_original", "ModuleID", "ModuleLabel_Final")) {
      if (!nm %in% names(module_meta_src)) module_meta_src[[nm]] <- NA_character_
    }
    module_meta <- module_meta_src |>
      dplyr::transmute(
        module_eigengene = .data$module_eigengene,
        module_eigengene_original = dplyr::coalesce(.data$module_eigengene_original, .data$module_eigengene),
        module_id = dplyr::coalesce(.data$ModuleID, paste0("WGCNA_", sub("^ME", "", .data$module_eigengene))),
        module_label = dplyr::coalesce(.data$ModuleLabel_Final, .data$module_id)
      )
    base <- module_meta |>
      dplyr::mutate(
        dataset = dataset,
        supermodule_id = sid,
        supermodule_label = label,
        pc = 1L,
        n_member_modules = length(members),
        n_variable_member_modules = length(variable_members),
        n_samples_used = n_complete,
        .before = "module_eigengene"
      )
    if (!length(variable_members)) {
      return(base |>
        dplyr::mutate(loading = NA_real_, abs_loading = NA_real_, loading_rank = NA_integer_, pca_status = "no_variable_member_modules") |>
        dplyr::select(dplyr::all_of(names(empty))))
    }
    if (length(variable_members) == 1L) {
      out <- base |>
        dplyr::filter(.data$module_eigengene == variable_members[[1]]) |>
        dplyr::mutate(loading = 1, abs_loading = 1, loading_rank = 1L, pca_status = "singleton_or_single_variable_module")
      return(out |> dplyr::select(dplyr::all_of(names(empty))))
    }
    vals_complete <- vals[stats::complete.cases(vals[, variable_members, drop = FALSE]), variable_members, drop = FALSE]
    if (nrow(vals_complete) < 2L) {
      return(base |>
        dplyr::filter(.data$module_eigengene %in% variable_members) |>
        dplyr::mutate(loading = NA_real_, abs_loading = NA_real_, loading_rank = NA_integer_, pca_status = "too_few_complete_samples") |>
        dplyr::select(dplyr::all_of(names(empty))))
    }
    pca <- tryCatch(stats::prcomp(vals_complete, center = TRUE, scale. = TRUE), error = function(e) e)
    if (inherits(pca, "error")) {
      return(base |>
        dplyr::filter(.data$module_eigengene %in% variable_members) |>
        dplyr::mutate(loading = NA_real_, abs_loading = NA_real_, loading_rank = NA_integer_, pca_status = paste0("pca_failed: ", conditionMessage(pca))) |>
        dplyr::select(dplyr::all_of(names(empty))))
    }
    sign_flip <- tryCatch({
      mean_vec <- rowMeans(vals_complete, na.rm = TRUE)
      stats::cor(pca$x[, 1L], mean_vec, use = "pairwise.complete.obs") < 0
    }, error = function(e) FALSE)
    loading <- pca$rotation[, 1L]
    if (isTRUE(sign_flip)) loading <- -loading
    load_df <- tibble::tibble(module_eigengene = names(loading), loading = as.numeric(loading))
    base |>
      dplyr::filter(.data$module_eigengene %in% variable_members) |>
      dplyr::left_join(load_df, by = "module_eigengene") |>
      dplyr::mutate(abs_loading = abs(.data$loading), pca_status = "ok") |>
      dplyr::arrange(dplyr::desc(.data$abs_loading), .data$module_eigengene) |>
      dplyr::mutate(loading_rank = dplyr::row_number()) |>
      dplyr::select(dplyr::all_of(names(empty)))
  }))
}

compact_supermodule_label <- function(supermodule_id, supermodule_label, width = 26L) {
  supermodule_id <- as.character(supermodule_id)
  label <- mapply(
    function(sid, lbl) gsub(paste0("^", sid, "\\s*(/|:|-|\\|)\\s*"), "", as.character(lbl), ignore.case = TRUE),
    supermodule_id,
    as.character(supermodule_label),
    USE.NAMES = FALSE
  )
  missing_label <- is.na(label) | !nzchar(label)
  label[missing_label] <- supermodule_id[missing_label]
  paste0(supermodule_id, "\n", wrap_plot_label(label, width = width))
}

plot_supermodule_pca_scree <- function(pca_eigenvalues, path, max_pc = 5L, title = "Supermodule PCA eigenvalue scree") {
  plot_df <- pca_eigenvalues |>
    dplyr::filter(!is.na(.data$pc), .data$pc <= max_pc, is.finite(.data$variance_explained)) |>
    dplyr::mutate(
      pc_label = paste0("PC", .data$pc),
      facet_label = compact_supermodule_label(.data$supermodule_id, .data$supermodule_label, width = 24L)
    )
  if (!nrow(plot_df)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No variable supermodule PCA components available\nSee supermodule_pca_input_audit.csv", size = 3) +
      ggplot2::labs(title = title) +
      ggplot2::theme_void()
  } else {
    line_df <- plot_df |>
      dplyr::group_by(.data$facet_label) |>
      dplyr::filter(dplyr::n() > 1L) |>
      dplyr::ungroup()
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$pc_label)) +
      ggplot2::geom_col(ggplot2::aes(y = .data$variance_explained), fill = "#7FA6C9", width = 0.64)
    if (nrow(line_df)) {
      p <- p +
        ggplot2::geom_line(
          data = line_df,
          ggplot2::aes(y = .data$cumulative_variance_explained, group = 1),
          color = "#C66A2B",
          linewidth = 0.3
        )
    }
    p <- p +
      ggplot2::geom_point(ggplot2::aes(y = .data$cumulative_variance_explained), color = "#C66A2B", size = 1.05) +
      ggplot2::facet_wrap(~facet_label, ncol = 3) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::labs(x = NULL, y = "Variance explained", title = title, subtitle = paste0("Bars show per-PC variance; points/line show cumulative variance. Display limited to first ", max_pc, " PCs.")) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.spacing = grid::unit(5, "mm"),
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 7),
        axis.text.y = ggplot2::element_text(size = 7),
        strip.text = ggplot2::element_text(size = 8, face = "bold", lineheight = 0.95),
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey35")
      )
  }
  ggplot2::ggsave(path, p, width = 210, height = 170, units = "mm", device = svglite::svglite)
}

standardize_contrast_key <- function(x) {
  x <- gsub("[\u2010-\u2015]", "-", as.character(x))
  toupper(gsub("[[:space:]]+", "", x))
}

is_sus_res_contrast <- function(x) {
  standardize_contrast_key(x) %in% c("SUS-RES", "RES-SUS")
}

orient_sus_minus_res_estimate <- function(estimate, contrast) {
  key <- standardize_contrast_key(contrast)
  out <- suppressWarnings(as.numeric(estimate))
  out[key == "RES-SUS"] <- -out[key == "RES-SUS"]
  out
}

contrast_block <- function(x) {
  key <- standardize_contrast_key(x)
  dplyr::case_when(
    key %in% c("RES-CON", "CON-RES") ~ "RES-CON",
    key %in% c("SUS-CON", "CON-SUS") ~ "SUS-CON",
    key %in% c("SUS-RES", "RES-SUS") ~ "SUS-RES",
    TRUE ~ as.character(x)
  )
}

contrast_block_order <- function(x) {
  block <- contrast_block(x)
  dplyr::case_when(
    block == "RES-CON" ~ 1L,
    block == "SUS-CON" ~ 2L,
    block == "SUS-RES" ~ 3L,
    TRUE ~ 99L
  )
}

parse_spatial_unit_for_effects <- function(x) {
  raw <- gsub("[^a-z0-9]+", "_", tolower(as.character(x)))
  raw <- gsub("^_+|_+$", "", raw)
  raw[!nzchar(raw) | is.na(raw) | raw %in% c("global", "global_spatial_adjusted", "all_spatial_units", "na")] <- "no_local_support"
  display <- toupper(gsub("_", " ", raw))
  display[raw == "no_local_support"] <- "No local support"
  region <- dplyr::case_when(
    grepl("^ca1($|_)", raw) ~ "CA1",
    grepl("^ca2($|_)", raw) ~ "CA2",
    grepl("^ca3($|_)", raw) ~ "CA3",
    grepl("^dg($|_)", raw) ~ "DG",
    raw == "no_local_support" ~ "Global/no local support",
    TRUE ~ "Other"
  )
  layer <- dplyr::case_when(
    raw == "no_local_support" ~ "Layer not available",
    grepl("_slm$", raw) ~ "SLM",
    grepl("_so$", raw) ~ "SO",
    grepl("_sr$", raw) ~ "SR",
    grepl("_sp$", raw) ~ "SP",
    grepl("_mo$", raw) ~ "MO",
    grepl("_po$", raw) ~ "PO",
    grepl("_ml$", raw) ~ "ML",
    grepl("_gcl$", raw) ~ "GCL",
    grepl("_hilus$", raw) ~ "Hilus",
    TRUE ~ "Layer not available"
  )
  tibble::tibble(
    parsed_region = region,
    parsed_layer_or_unit = layer,
    parsed_spatial_unit_display = display
  )
}

effect_support_class <- function(p_value, fdr_within, fdr_global) {
  fdr <- dplyr::coalesce(suppressWarnings(as.numeric(fdr_within)), suppressWarnings(as.numeric(fdr_global)))
  p <- suppressWarnings(as.numeric(p_value))
  dplyr::case_when(
    !is.na(fdr) & fdr <= 0.05 ~ "FDR05",
    !is.na(fdr) & fdr <= 0.10 ~ "FDR10",
    !is.na(p) & p < 0.05 ~ "nominal",
    TRUE ~ "none"
  )
}

spatial_order_value <- function(region, layer_or_unit, unit) {
  unit <- tolower(as.character(unit))
  layer <- toupper(as.character(layer_or_unit))
  dplyr::case_when(
    region == "CA1" & layer == "SO" ~ 101L,
    region == "CA1" & layer == "SP" ~ 102L,
    region == "CA1" & layer == "SR" ~ 103L,
    region == "CA1" & layer == "SLM" ~ 104L,
    region == "CA1" ~ 109L,
    region == "CA2" ~ 201L,
    region == "CA3" ~ 221L,
    region == "DG" & layer == "MO" ~ 301L,
    region == "DG" & layer == "ML" ~ 302L,
    region == "DG" & layer == "GCL" ~ 303L,
    region == "DG" & layer == "PO" ~ 304L,
    region == "DG" ~ 309L,
    region == "Global/no local support" ~ 998L,
    TRUE ~ 999L
  )
}

order_spatial_labels_for_plot <- function(labels) {
  labels <- unique(as.character(labels))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (!length(labels)) return(labels)
  parsed <- parse_spatial_unit_for_effects(labels)
  order_df <- tibble::tibble(SpatialLabel = labels) |>
    dplyr::bind_cols(parsed) |>
    dplyr::mutate(spatial_order = spatial_order_value(.data$parsed_region, .data$parsed_layer_or_unit, .data$SpatialLabel)) |>
    dplyr::arrange(.data$spatial_order, .data$SpatialLabel)
  order_df$SpatialLabel
}

plot_all_supermodule_eigengene_group <- function(values, path) {
  plot_df <- values |>
    dplyr::filter(
      is.finite(.data$eigengene_z),
      !is.na(.data$StressGroup),
      as.character(.data$StressGroup) %in% c("CON", "RES", "SUS")
    ) |>
    dplyr::mutate(
      StressGroup = factor(as.character(.data$StressGroup), levels = c("CON", "RES", "SUS")),
      supermodule_plot_label = compact_supermodule_label(.data$supermodule_id, .data$supermodule_label, width = 22L)
    )
  if (!nrow(plot_df)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No all-supermodule eigengene group values available", size = 3) +
      ggplot2::labs(title = "All supermodule eigengenes by stress group") +
      ggplot2::theme_void()
    ggplot2::ggsave(path, p, width = 180, height = 105, units = "mm", device = svglite::svglite)
    return(invisible(path))
  }

  n_supermodules <- dplyr::n_distinct(plot_df$supermodule_id)
  ncol <- dplyr::case_when(
    n_supermodules <= 4L ~ 2L,
    n_supermodules <= 9L ~ 3L,
    TRUE ~ 4L
  )
  height_mm <- max(115, 48 * ceiling(n_supermodules / ncol))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$StressGroup, y = .data$eigengene_z)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = .data$StressGroup),
      width = 0.58,
      outlier.shape = NA,
      alpha = 0.42,
      linewidth = 0.25,
      color = "grey30"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$StressGroup),
      position = ggplot2::position_jitter(width = 0.12, height = 0),
      size = 0.75,
      alpha = 0.72,
      stroke = 0
    ) +
    ggplot2::facet_wrap(~supermodule_plot_label, ncol = ncol) +
    ggplot2::scale_fill_manual(values = c("CON" = "#8B8F96", "RES" = "#4F7CAC", "SUS" = "#B8664B"), drop = FALSE) +
    ggplot2::scale_color_manual(values = c("CON" = "#62666D", "RES" = "#315E8A", "SUS" = "#8F4330"), drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = "Supermodule eigengene (z)",
      title = "All supermodule eigengenes by stress group"
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.2),
      axis.line = ggplot2::element_line(color = "grey35", linewidth = 0.2),
      axis.ticks = ggplot2::element_line(color = "grey35", linewidth = 0.2),
      axis.text.x = ggplot2::element_text(size = 6.5),
      axis.text.y = ggplot2::element_text(size = 6.5),
      strip.text = ggplot2::element_text(size = 7, face = "bold", lineheight = 0.95),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    )
  ggplot2::ggsave(path, p, width = 215, height = height_mm, units = "mm", device = svglite::svglite, limitsize = FALSE)
}

plot_all_supermodule_eigengene_spatial_group <- function(values, path) {
  plot_df <- values |>
    dplyr::filter(
      is.finite(.data$eigengene_z),
      !is.na(.data$StressGroup),
      as.character(.data$StressGroup) %in% c("CON", "RES", "SUS"),
      !is.na(.data$SpatialLabel),
      nzchar(as.character(.data$SpatialLabel))
    ) |>
    dplyr::mutate(
      StressGroup = factor(as.character(.data$StressGroup), levels = c("CON", "RES", "SUS")),
      supermodule_plot_label = compact_supermodule_label(.data$supermodule_id, .data$supermodule_label, width = 22L)
    )
  if (!nrow(plot_df)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No all-supermodule spatial eigengene group values available", size = 3) +
      ggplot2::labs(title = "All supermodule eigengenes by spatial unit and stress group") +
      ggplot2::theme_void()
    ggplot2::ggsave(path, p, width = 190, height = 105, units = "mm", device = svglite::svglite)
    return(invisible(path))
  }

  spatial_levels <- order_spatial_labels_for_plot(plot_df$SpatialLabel)
  plot_df$SpatialLabel <- factor(as.character(plot_df$SpatialLabel), levels = spatial_levels)
  n_supermodules <- dplyr::n_distinct(plot_df$supermodule_id)
  n_spatial <- length(spatial_levels)
  ncol <- if (n_spatial >= 6L || n_supermodules <= 4L) 2L else 3L
  height_mm <- max(125, 52 * ceiling(n_supermodules / ncol))
  width_mm <- max(205, min(290, 140 + 12 * n_spatial))
  dodge <- ggplot2::position_dodge(width = 0.72)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$SpatialLabel, y = .data$eigengene_z, fill = .data$StressGroup, color = .data$StressGroup)) +
    ggplot2::geom_boxplot(
      position = dodge,
      width = 0.62,
      outlier.shape = NA,
      alpha = 0.36,
      linewidth = 0.23,
      color = "grey35"
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, jitter.height = 0, dodge.width = 0.72),
      size = 0.55,
      alpha = 0.62,
      stroke = 0
    ) +
    ggplot2::facet_wrap(~supermodule_plot_label, ncol = ncol) +
    ggplot2::scale_fill_manual(values = c("CON" = "#8B8F96", "RES" = "#4F7CAC", "SUS" = "#B8664B"), drop = FALSE) +
    ggplot2::scale_color_manual(values = c("CON" = "#62666D", "RES" = "#315E8A", "SUS" = "#8F4330"), drop = FALSE) +
    ggplot2::labs(
      x = "Spatial unit",
      y = "Supermodule eigengene (z)",
      fill = NULL,
      color = NULL,
      title = "All supermodule eigengenes by spatial unit and stress group"
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.2),
      axis.line = ggplot2::element_line(color = "grey35", linewidth = 0.2),
      axis.ticks = ggplot2::element_line(color = "grey35", linewidth = 0.2),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5.6),
      axis.text.y = ggplot2::element_text(size = 6),
      strip.text = ggplot2::element_text(size = 7, face = "bold", lineheight = 0.95),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    )
  ggplot2::ggsave(path, p, width = width_mm, height = height_mm, units = "mm", device = svglite::svglite, limitsize = FALSE)
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
    evidence_status = NA_character_,
    selection_support = "none",
    selection_rank = NA_integer_,
    selection_message = message
  )
  if (is.null(super_effects) || !nrow(super_effects)) {
    return(no_selection_row("no supermodule group-effect rows available"))
  }
  sus_res <- super_effects |>
    dplyr::filter(is_sus_res_contrast(.data$contrast)) |>
    dplyr::mutate(
      p_value = suppressWarnings(as.numeric(.data$p_value)),
      FDR_within_dataset_level = suppressWarnings(as.numeric(.data$FDR_within_dataset_level)),
      FDR_global = suppressWarnings(as.numeric(.data$FDR_global)),
      estimate = orient_sus_minus_res_estimate(.data$estimate, .data$contrast),
      contrast = "SUS - RES",
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
      "FDR_within_dataset_level", "FDR_global", "evidence_status", "selection_support",
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
      dplyr::filter(!is.na(.data$pc), .data$pc <= 5L, is.finite(.data$variance_explained)) |>
      dplyr::mutate(
        pc_label = paste0("PC", .data$pc),
        facet_label = paste0(
          compact_supermodule_label(.data$supermodule_id, .data$supermodule_label, width = 24L),
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
        ggplot2::geom_col(ggplot2::aes(y = .data$variance_explained), fill = "#7FA6C9", width = 0.64)
      if (nrow(line_df)) {
        p <- p +
          ggplot2::geom_line(
            data = line_df,
            ggplot2::aes(y = .data$cumulative_variance_explained, group = 1),
            color = "#C66A2B",
            linewidth = 0.3
          )
      }
      p <- p +
        ggplot2::geom_point(ggplot2::aes(y = .data$cumulative_variance_explained), color = "#C66A2B", size = 1.05) +
        ggplot2::facet_wrap(~facet_label, ncol = 2) +
        ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::labs(
          x = NULL,
          y = "Variance explained",
          title = "Selected SUS-RES supermodule PCA eigenvalue scree",
          subtitle = "Bars show per-PC variance; points/line show cumulative variance. Display limited to first 5 PCs."
        ) +
        ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.spacing = grid::unit(5, "mm"),
          legend.position = "none",
          axis.text.x = ggplot2::element_text(size = 7),
          axis.text.y = ggplot2::element_text(size = 7),
          strip.text = ggplot2::element_text(size = 8, face = "bold", lineheight = 0.95),
          plot.title = ggplot2::element_text(face = "bold", size = 11),
          plot.subtitle = ggplot2::element_text(size = 8, color = "grey35")
        )
    }
  }
  ggplot2::ggsave(path, p, width = 185, height = 115, units = "mm", device = svglite::svglite)
}

selected_sus_res_member_loadings <- function(member_loadings, selected, contents) {
  if (is.null(selected) || !nrow(selected) || all(selected$selection_message != "selected")) {
    return(member_loadings[0, , drop = FALSE])
  }
  selected <- selected |> dplyr::filter(.data$selection_message == "selected")
  content_meta <- if (!is.null(contents) && nrow(contents)) {
    contents |>
      dplyr::select("supermodule_id", "contrast", "estimate", "p_value", "FDR_within_dataset_level", "FDR_global", "evidence_status", "selection_support")
  } else {
    tibble::tibble(
      supermodule_id = character(), contrast = character(), estimate = numeric(), p_value = numeric(),
      FDR_within_dataset_level = numeric(), FDR_global = numeric(), evidence_status = character(), selection_support = character()
    )
  }
  member_loadings |>
    dplyr::inner_join(selected |> dplyr::select("supermodule_id", "selection_rank"), by = "supermodule_id") |>
    dplyr::left_join(content_meta, by = "supermodule_id") |>
    dplyr::arrange(.data$selection_rank, .data$loading_rank)
}

plot_selected_sus_res_member_loadings <- function(selected_loadings, path) {
  plot_df <- selected_loadings |>
    dplyr::filter(is.finite(.data$abs_loading)) |>
    dplyr::group_by(.data$supermodule_id) |>
    dplyr::slice_max(.data$abs_loading, n = 8, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      module_plot_label = dplyr::coalesce(.data$module_label, .data$module_id, .data$module_eigengene),
      module_plot_label = paste0(.data$module_plot_label, "\n", .data$module_id),
      facet_label = paste0(
        compact_supermodule_label(.data$supermodule_id, .data$supermodule_label, width = 24L),
        "\nSUS-RES est=", signif(.data$estimate, 3), ", p=", signif(.data$p_value, 3),
        ", FDR=", signif(.data$FDR_within_dataset_level, 3)
      ),
      signed_direction = ifelse(.data$loading >= 0, "positive PC1 loading", "negative PC1 loading")
    )
  if (!nrow(plot_df)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No PCA member loadings available for selected SUS-RES supermodules", size = 3) +
      ggplot2::labs(title = "Selected SUS-RES supermodule PC1 member loadings") +
      ggplot2::theme_void()
  } else {
    plot_df <- plot_df |>
      dplyr::group_by(.data$facet_label) |>
      dplyr::mutate(module_plot_label = stats::reorder(.data$module_plot_label, .data$abs_loading)) |>
      dplyr::ungroup()
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$module_plot_label, y = .data$abs_loading, fill = .data$signed_direction)) +
      ggplot2::geom_col(width = 0.66, color = "grey30", linewidth = 0.15) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~facet_label, scales = "free_y", ncol = 1) +
      ggplot2::scale_fill_manual(values = c("positive PC1 loading" = "#7FA6C9", "negative PC1 loading" = "#D8A070")) +
      ggplot2::labs(
        x = NULL,
        y = "|PC1 loading|",
        fill = NULL,
        title = "Selected SUS-RES supermodule PC1 member loadings",
        subtitle = "Higher absolute loading indicates stronger contribution of a member module to the supermodule eigengene."
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 8, face = "bold", lineheight = 0.95),
        axis.text.y = ggplot2::element_text(size = 7),
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey35")
      )
  }
  ggplot2::ggsave(path, p, width = 190, height = 125, units = "mm", device = svglite::svglite)
}

plot_all_supermodule_member_loadings <- function(member_loadings, path, max_modules_per_supermodule = 5L) {
  plot_df <- member_loadings |>
    dplyr::filter(is.finite(.data$abs_loading)) |>
    dplyr::group_by(.data$supermodule_id) |>
    dplyr::slice_max(.data$abs_loading, n = max_modules_per_supermodule, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      module_plot_label = dplyr::coalesce(.data$module_label, .data$module_id, .data$module_eigengene),
      facet_label = compact_supermodule_label(.data$supermodule_id, .data$supermodule_label, width = 20L),
      signed_direction = ifelse(.data$loading >= 0, "positive PC1 loading", "negative PC1 loading")
    )
  if (!nrow(plot_df)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No PCA member loadings available", size = 3) +
      ggplot2::labs(title = "All supermodule PC1 member loadings") +
      ggplot2::theme_void()
  } else {
    plot_df <- plot_df |>
      dplyr::group_by(.data$facet_label) |>
      dplyr::mutate(module_plot_label = stats::reorder(.data$module_plot_label, .data$abs_loading)) |>
      dplyr::ungroup()
    n_supermodules <- dplyr::n_distinct(plot_df$supermodule_id)
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$module_plot_label, y = .data$abs_loading, fill = .data$signed_direction)) +
      ggplot2::geom_col(width = 0.64, color = "grey35", linewidth = 0.12) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~facet_label, scales = "free_y", ncol = 3) +
      ggplot2::scale_fill_manual(values = c("positive PC1 loading" = "#7FA6C9", "negative PC1 loading" = "#D8A070")) +
      ggplot2::labs(
        x = NULL,
        y = "|PC1 loading|",
        fill = NULL,
        title = "All supermodule PC1 member loadings",
        subtitle = paste0("Top ", max_modules_per_supermodule, " member modules per supermodule by absolute PC1 loading.")
      ) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 7, face = "bold", lineheight = 0.95),
        axis.text.y = ggplot2::element_text(size = 5.5),
        axis.text.x = ggplot2::element_text(size = 6),
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey35")
      )
  }
  height_mm <- if (exists("n_supermodules", inherits = FALSE)) max(150, 34 * ceiling(n_supermodules / 3)) else 120
  ggplot2::ggsave(path, p, width = 230, height = height_mm, units = "mm", device = svglite::svglite, limitsize = FALSE)
}

all_supermodule_region_layer_effect_source <- function(super_effects) {
  if (is.null(super_effects) || !nrow(super_effects)) {
    return(tibble::tibble(
      dataset = character(), supermodule_id = character(), supermodule_label = character(),
      spatial_unit = character(), spatial_unit_label = character(), effect_scope = character(),
      SpatialUnitType = character(), parsed_region = character(), parsed_layer_or_unit = character(),
      contrast = character(), contrast_block = character(), contrast_block_order = integer(),
      estimate = numeric(), estimate_sus_minus_res = numeric(), p_value = numeric(),
      FDR_within_dataset_level = numeric(), FDR_global = numeric(),
      evidence_status = character(), support_class = character(), spatial_order = integer()
    ))
  }
  parsed <- parse_spatial_unit_for_effects(super_effects$spatial_unit)
  super_effects |>
    dplyr::bind_cols(parsed) |>
    dplyr::mutate(
      estimate = suppressWarnings(as.numeric(.data$estimate)),
      estimate_sus_minus_res = dplyr::if_else(is_sus_res_contrast(.data$contrast), orient_sus_minus_res_estimate(.data$estimate, .data$contrast), NA_real_),
      p_value = suppressWarnings(as.numeric(.data$p_value)),
      FDR_within_dataset_level = suppressWarnings(as.numeric(.data$FDR_within_dataset_level)),
      FDR_global = suppressWarnings(as.numeric(.data$FDR_global)),
      support_class = effect_support_class(.data$p_value, .data$FDR_within_dataset_level, .data$FDR_global),
      spatial_unit = dplyr::coalesce(.data$spatial_unit, "global"),
      spatial_unit_label = .data$parsed_spatial_unit_display,
      spatial_order = spatial_order_value(.data$parsed_region, .data$parsed_layer_or_unit, .data$spatial_unit),
      supermodule_label = dplyr::coalesce(.data$supermodule_label, .data$endpoint_label, .data$supermodule_id),
      contrast_block = contrast_block(.data$contrast),
      contrast_block_order = contrast_block_order(.data$contrast)
    ) |>
    dplyr::select(
      "dataset", "supermodule_id", "supermodule_label", "spatial_unit", "spatial_unit_label",
      "effect_scope", "SpatialUnitType", "parsed_region", "parsed_layer_or_unit",
      "contrast", "contrast_block", "contrast_block_order", "estimate", "estimate_sus_minus_res",
      "p_value", "FDR_within_dataset_level", "FDR_global", "evidence_status", "support_class", "spatial_order"
    )
}

supermodule_score_effects_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "module_score", dataset, "wgcna", "supermodule_directional_effects.csv")
}

publication_supermodule_region_layer_effect_source <- function(dataset, supermodule_contents) {
  path <- supermodule_score_effects_path(dataset)
  score <- safe_read_csv(path)
  empty <- tibble::tibble(
    dataset = character(), supermodule_id = character(), supermodule_label = character(),
    spatial_unit = character(), spatial_unit_label = character(), effect_scope = character(),
    SpatialUnitType = character(), parsed_region = character(), parsed_layer_or_unit = character(),
    contrast = character(), contrast_block = character(), contrast_block_order = integer(),
    estimate = numeric(), estimate_sus_minus_res = numeric(), Cohen_d = numeric(),
    p_value = numeric(), FDR_within_dataset_level = numeric(), FDR_global = numeric(),
    p_adj_within_model_BH = numeric(), p_adj_global_BH = numeric(),
    evidence_status = character(), support_class = character(), spatial_order = integer(),
    metric_used = character(), source_name = character(), source_path = character(),
    analysis = character(), significance_marker_column = character(), significance_rule = character()
  )
  if (is.null(score) || !nrow(score)) return(empty)
  required <- c("Analysis", "RegionLayer", "Module", "Cohen_d", "p_adj_within_model_BH")
  if (length(setdiff(required, names(score)))) return(empty)

  label_lookup <- if (!is.null(supermodule_contents) && nrow(supermodule_contents)) {
    supermodule_contents |>
      dplyr::transmute(
        supermodule_id = as.character(.data$supermodule_id),
        supermodule_label_lookup = as.character(.data$supermodule_label)
      ) |>
      dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)
  } else {
    tibble::tibble(supermodule_id = character(), supermodule_label_lookup = character())
  }

  within_flag <- if ("within_BH_significant" %in% names(score)) {
    trimws(as.character(score$within_BH_significant)) %in% c("TRUE", "true", "True", "1")
  } else {
    rep(FALSE, nrow(score))
  }
  contrast_source <- dplyr::coalesce(
    col_if_present(score, "Contrast"),
    col_if_present(score, "contrast_label"),
    col_if_present(score, "contrast")
  )
  p_nominal <- dplyr::coalesce(
    suppressWarnings(as.numeric(col_if_present(score, "p.value"))),
    suppressWarnings(as.numeric(col_if_present(score, "p_nominal")))
  )
  parsed <- parse_spatial_unit_for_effects(score$RegionLayer)

  score |>
    dplyr::bind_cols(parsed) |>
    dplyr::mutate(
      dataset = dataset,
      supermodule_id = as.character(.data$Module),
      supermodule_label = .data$supermodule_id,
      spatial_unit = as.character(.data$RegionLayer),
      spatial_unit_label = .data$parsed_spatial_unit_display,
      effect_scope = "publication_score_region_layer",
      SpatialUnitType = SpatialUnitType,
      contrast = as.character(.env$contrast_source),
      contrast_block = contrast_block(.data$contrast),
      contrast_block_order = contrast_block_order(.data$contrast),
      estimate = suppressWarnings(as.numeric(.data$Cohen_d)),
      Cohen_d = .data$estimate,
      estimate_sus_minus_res = dplyr::if_else(.data$contrast_block == "SUS-RES", orient_sus_minus_res_estimate(.data$estimate, .data$contrast), NA_real_),
      p_value = .env$p_nominal,
      p_adj_within_model_BH = suppressWarnings(as.numeric(.data$p_adj_within_model_BH)),
      p_adj_global_BH = suppressWarnings(as.numeric(col_if_present(score, "p_adj_global_BH"))),
      FDR_within_dataset_level = .data$p_adj_within_model_BH,
      FDR_global = .data$p_adj_global_BH,
      support_class = dplyr::if_else(.env$within_flag | (!is.na(.data$p_adj_within_model_BH) & .data$p_adj_within_model_BH <= 0.05), "FDR05", "none"),
      evidence_status = dplyr::if_else(.data$support_class == "FDR05", "robust_FDR", "not_supported"),
      spatial_order = spatial_order_value(.data$parsed_region, .data$parsed_layer_or_unit, .data$spatial_unit),
      metric_used = "Cohen_d",
      source_name = "publication_score_supermodule_directional_effects",
      source_path = path,
      analysis = as.character(.data$Analysis),
      significance_marker_column = "within_BH_significant/p_adj_within_model_BH",
      significance_rule = "within_BH_significant == TRUE or p_adj_within_model_BH <= 0.05"
    ) |>
    dplyr::filter(.data$analysis == "primary_all_replicates") |>
    dplyr::left_join(label_lookup, by = "supermodule_id") |>
    dplyr::mutate(supermodule_label = dplyr::coalesce(.data$supermodule_label_lookup, .data$supermodule_label)) |>
    dplyr::select(
      "dataset", "supermodule_id", "supermodule_label", "spatial_unit", "spatial_unit_label",
      "effect_scope", "SpatialUnitType", "parsed_region", "parsed_layer_or_unit",
      "contrast", "contrast_block", "contrast_block_order", "estimate", "estimate_sus_minus_res",
      "Cohen_d", "p_value", "FDR_within_dataset_level", "FDR_global",
      "p_adj_within_model_BH", "p_adj_global_BH", "evidence_status", "support_class",
      "spatial_order", "metric_used", "source_name", "source_path", "analysis",
      "significance_marker_column", "significance_rule"
    ) |>
    dplyr::filter(is.finite(.data$estimate), nzchar(.data$supermodule_id), nzchar(.data$spatial_unit))
}

sus_res_region_layer_effect_source <- function(super_effects) {
  all_supermodule_region_layer_effect_source(super_effects) |>
    dplyr::filter(.data$contrast_block == "SUS-RES") |>
    dplyr::mutate(contrast = "SUS - RES")
}

plot_all_sus_res_region_layer_effects <- function(effect_source, path) {
  plot_df <- effect_source |>
    dplyr::filter(is.finite(.data$estimate_sus_minus_res), !is.na(.data$p_value)) |>
    dplyr::mutate(
      supermodule_plot_label = compact_supermodule_label(.data$supermodule_id, .data$supermodule_label, width = 22L),
      spatial_unit_label = dplyr::coalesce(.data$spatial_unit_label, "GLOBAL"),
      support_class = factor(.data$support_class, levels = c("FDR05", "FDR10", "nominal", "none")),
      neg_log10_p = pmin(8, -log10(pmax(.data$p_value, .Machine$double.xmin)))
    )
  if (!nrow(plot_df)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No SUS-RES region/layer supermodule effects available", size = 3) +
      ggplot2::labs(title = "All SUS-RES supermodule effects by region/layer") +
      ggplot2::theme_void()
  } else {
    uses_cohend <- "metric_used" %in% names(plot_df) && any(plot_df$metric_used == "Cohen_d", na.rm = TRUE)
    fill_label <- if (uses_cohend) "SUS-RES\nCohen's d" else "SUS-RES\nestimate"
    subtitle <- if (uses_cohend) {
      "Values are Cohen's d oriented as SUS minus RES; markers use within-model BH adjusted p <= 0.05."
    } else {
      "Estimates are oriented as SUS minus RES; positive values are higher in SUS."
    }
    spatial_levels <- plot_df |>
      dplyr::distinct(.data$spatial_unit_label) |>
      dplyr::pull(.data$spatial_unit_label)
    supermodule_levels <- plot_df |>
      dplyr::group_by(.data$supermodule_plot_label) |>
      dplyr::summarise(best = max(abs(.data$estimate_sus_minus_res), na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.data$best, .data$supermodule_plot_label) |>
      dplyr::pull(.data$supermodule_plot_label)
    plot_df$spatial_unit_label <- factor(plot_df$spatial_unit_label, levels = spatial_levels)
    plot_df$supermodule_plot_label <- factor(plot_df$supermodule_plot_label, levels = supermodule_levels)
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$spatial_unit_label, y = .data$supermodule_plot_label)) +
      ggplot2::geom_point(ggplot2::aes(fill = .data$estimate_sus_minus_res, size = .data$neg_log10_p, shape = .data$support_class), color = "grey25", stroke = 0.22, alpha = 0.92) +
      ggplot2::facet_wrap(~effect_scope, scales = "free_x") +
      ggplot2::scale_fill_gradient2(low = "#3B6EA8", mid = "grey96", high = "#B84A4A", midpoint = 0, name = fill_label) +
      ggplot2::scale_size_continuous(range = c(1.4, 4.6), name = "-log10 p") +
      ggplot2::scale_shape_manual(values = c("FDR05" = 21, "FDR10" = 22, "nominal" = 24, "none" = 21), labels = c("FDR <= 0.05", "FDR <= 0.10", "nominal p < 0.05", "not supported"), drop = FALSE, name = "Support") +
      ggplot2::labs(
        x = "Region / layer",
        y = NULL,
        title = "All SUS-RES supermodule effects by region/layer",
        subtitle = subtitle
      ) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        panel.grid = ggplot2::element_line(color = "grey90", linewidth = 0.2),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = ggplot2::element_text(size = 5.5),
        strip.text = ggplot2::element_text(size = 8, face = "bold"),
        legend.position = "right",
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey35")
      )
  }
  ggplot2::ggsave(path, p, width = 235, height = 165, units = "mm", device = svglite::svglite)
}

plot_all_supermodule_region_layer_heatmap <- function(effect_source, path) {
  plot_df <- effect_source |>
    dplyr::filter(is.finite(.data$estimate), !is.na(.data$p_value)) |>
    dplyr::mutate(
      supermodule_plot_label = compact_supermodule_label(.data$supermodule_id, .data$supermodule_label, width = 22L),
      spatial_plot_label = paste(.data$parsed_region, .data$parsed_layer_or_unit, sep = " / "),
      spatial_plot_label = gsub(" / Layer not available$", "", .data$spatial_plot_label),
      support_marker = dplyr::case_when(
        .data$support_class == "FDR05" ~ "FDR <= 0.05",
        .data$support_class == "FDR10" ~ "FDR <= 0.10",
        .data$support_class == "nominal" ~ "nominal p < 0.05",
        TRUE ~ NA_character_
      )
    )
  if (!nrow(plot_df)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No supermodule region/layer effects available", size = 3) +
      ggplot2::labs(title = "All supermodule effects by region/layer") +
      ggplot2::theme_void()
  } else {
    uses_cohend <- "metric_used" %in% names(plot_df) && any(plot_df$metric_used == "Cohen_d", na.rm = TRUE)
    fill_label <- if (uses_cohend) "Cohen's d" else "Estimate"
    subtitle <- if (uses_cohend) {
      "Tile fill is Cohen's d from the publication score source; markers use within-model BH adjusted p <= 0.05."
    } else {
      "Tile fill is the source model estimate; markers use the same FDR/nominal support classes as the circular heatmap source."
    }
    spatial_levels <- plot_df |>
      dplyr::distinct(.data$spatial_plot_label, .data$spatial_order) |>
      dplyr::arrange(.data$spatial_order, .data$spatial_plot_label) |>
      dplyr::pull(.data$spatial_plot_label)
    supermodule_levels <- plot_df |>
      dplyr::group_by(.data$supermodule_plot_label) |>
      dplyr::summarise(best = max(abs(.data$estimate), na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(.data$best, .data$supermodule_plot_label) |>
      dplyr::pull(.data$supermodule_plot_label)
    plot_df$spatial_plot_label <- factor(plot_df$spatial_plot_label, levels = unique(spatial_levels))
    plot_df$supermodule_plot_label <- factor(plot_df$supermodule_plot_label, levels = unique(supermodule_levels))
    lim <- stats::quantile(abs(plot_df$estimate[is.finite(plot_df$estimate)]), 0.95, na.rm = TRUE, names = FALSE)
    lim <- max(lim, 1e-6)
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$spatial_plot_label, y = .data$supermodule_plot_label)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$estimate), color = "white", linewidth = 0.12) +
      ggplot2::geom_point(
        data = plot_df |> dplyr::filter(!is.na(.data$support_marker)),
        ggplot2::aes(shape = .data$support_marker),
        size = 1.15,
        color = "black",
        stroke = 0.25
      ) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$contrast_block), cols = ggplot2::vars(.data$effect_scope), scales = "free_x", space = "free_x") +
      ggplot2::scale_fill_gradient2(low = "#3B6EA8", mid = "grey96", high = "#B84A4A", midpoint = 0, limits = c(-lim, lim), oob = scales::squish, name = fill_label) +
      ggplot2::scale_shape_manual(values = c("FDR <= 0.05" = 16, "FDR <= 0.10" = 1, "nominal p < 0.05" = 4), name = "Support") +
      ggplot2::labs(
        x = "Region / layer",
        y = NULL,
        title = "All supermodule effects by region/layer",
        subtitle = subtitle
      ) +
      ggplot2::theme_minimal(base_size = 7) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 5.5),
        axis.text.y = ggplot2::element_text(size = 5.2),
        strip.text = ggplot2::element_text(size = 7, face = "bold"),
        legend.position = "right",
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey35")
      )
  }
  ggplot2::ggsave(path, p, width = 250, height = 180, units = "mm", device = svglite::svglite, limitsize = FALSE)
}

split_member_modules <- function(x) {
  z <- unlist(strsplit(as.character(x %||% ""), ";", fixed = TRUE))
  z <- trimws(z)
  z[nzchar(z)]
}

selected_sus_res_supermodule_contents <- function(selected, pca_eigenvalues, comp, definitions, super_ann, super_summary, modules_long, module_summary, go_enrichment, module_name_map, dataset) {
  empty <- tibble::tibble(
    dataset = character(), supermodule_id = character(), supermodule_label = character(),
    contrast = character(), estimate = numeric(), p_value = numeric(),
    FDR_within_dataset_level = numeric(), FDR_global = numeric(), evidence_status = character(),
    pca_PC1_variance_explained = numeric(), n_member_modules = integer(),
    member_modules = character(), member_module_labels = character(),
    top_GO_BP_terms = character(), top_GO_MF_terms = character(), top_GO_CC_terms = character(),
    top_hub_symbols = character(), top_hub_proteins = character(), top_core_kME_proteins = character(),
    selection_support = character(), selection_message = character()
  )
  if (is.null(selected) || !nrow(selected) || all(!selected$selection_message %in% c("selected", "all_supermodules"))) return(empty)
  selected <- selected |>
    dplyr::filter(.data$selection_message %in% c("selected", "all_supermodules")) |>
    dplyr::mutate(
      effect_support_label = dplyr::case_when(
        !is.na(.data$FDR_within_dataset_level) & .data$FDR_within_dataset_level <= 0.05 ~ "FDR <= 0.05",
        !is.na(.data$FDR_within_dataset_level) & .data$FDR_within_dataset_level <= 0.10 ~ "FDR <= 0.10",
        !is.na(.data$p_value) & .data$p_value < 0.05 ~ "nominal-only p < 0.05",
        TRUE ~ "not selected"
      )
    )

  pc1 <- pca_eigenvalues |>
    dplyr::filter(.data$pc == 1L) |>
    dplyr::select("supermodule_id", pca_PC1_variance_explained = "variance_explained")

  summary_meta <- if (!is.null(super_summary) && nrow(super_summary) && "SupermoduleID" %in% names(super_summary)) {
    super_summary |>
      dplyr::transmute(
        supermodule_id = as.character(.data$SupermoduleID),
        summary_label = dplyr::coalesce(col_if_present(super_summary, "Supermodule_DisplayLabel"), col_if_present(super_summary, "Supermodule_FinalLabel"), as.character(.data$SupermoduleID)),
        summary_member_modules = col_if_present(super_summary, "member_modules"),
        summary_n_member_modules = suppressWarnings(as.integer(col_if_present(super_summary, "n_modules", NA_integer_))),
        top_GO_BP_terms = col_if_present(super_summary, "top_GO_BP_terms"),
        top_GO_MF_terms = col_if_present(super_summary, "top_GO_MF_terms"),
        top_GO_CC_terms = col_if_present(super_summary, "top_GO_CC_terms"),
        top_hub_symbols = col_if_present(super_summary, "top_hub_symbols"),
        top_hub_proteins_summary = dplyr::coalesce(col_if_present(super_summary, "top_hub_proteins"), col_if_present(super_summary, "top_hub_uniprot"))
      ) |>
      dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)
  } else {
    tibble::tibble(
      supermodule_id = character(), summary_label = character(), summary_member_modules = character(),
      summary_n_member_modules = integer(), top_GO_BP_terms = character(), top_GO_MF_terms = character(),
      top_GO_CC_terms = character(), top_hub_symbols = character(), top_hub_proteins_summary = character()
    )
  }

  definition_lookup <- if (!is.null(definitions) && nrow(definitions)) {
    definitions |>
      dplyr::transmute(
        module_eigengene = as.character(col_if_present(definitions, "module_eigengene")),
        ModuleID = as.character(col_if_present(definitions, "ModuleID")),
        ModuleLabel_Final = as.character(col_if_present(definitions, "ModuleLabel_Final"))
      ) |>
      dplyr::filter(nzchar(.data$module_eigengene) | nzchar(.data$ModuleID)) |>
      dplyr::distinct()
  } else {
    tibble::tibble(module_eigengene = character(), ModuleID = character(), ModuleLabel_Final = character())
  }
  name_map_lookup <- if (!is.null(module_name_map) && nrow(module_name_map) && "ModuleID" %in% names(module_name_map)) {
    module_name_map |>
      dplyr::transmute(
        ModuleID = as.character(.data$ModuleID),
        ModuleLabel_Final = dplyr::coalesce(col_if_present(module_name_map, "ModuleLabel_Final"), col_if_present(module_name_map, "final_label"), col_if_present(module_name_map, "primary_label"))
      ) |>
      dplyr::filter(nzchar(.data$ModuleID)) |>
      dplyr::distinct(.data$ModuleID, .keep_all = TRUE)
  } else {
    tibble::tibble(ModuleID = character(), ModuleLabel_Final = character())
  }

  ann_module_meta <- if (!is.null(super_ann) && nrow(super_ann) && "SupermoduleID" %in% names(super_ann)) {
    super_ann |>
      dplyr::transmute(
        supermodule_id = as.character(.data$SupermoduleID),
        module_eigengene = as.character(col_if_present(super_ann, "module_eigengene")),
        ModuleID = as.character(col_if_present(super_ann, "ModuleID")),
        ModuleLabel_Final = as.character(col_if_present(super_ann, "ModuleLabel_Final"))
      ) |>
      dplyr::filter(nzchar(.data$supermodule_id)) |>
      dplyr::distinct()
  } else {
    tibble::tibble(supermodule_id = character(), module_eigengene = character(), ModuleID = character(), ModuleLabel_Final = character())
  }

  comp_long <- if (!is.null(comp) && nrow(comp) && all(c("supermodule_id", "member_modules") %in% names(comp))) {
    dplyr::bind_rows(lapply(seq_len(nrow(comp)), function(i) {
      tibble::tibble(
        supermodule_id = as.character(comp$supermodule_id[[i]]),
        module_eigengene = split_member_modules(comp$member_modules[[i]])
      )
    })) |>
      dplyr::left_join(definition_lookup, by = "module_eigengene") |>
      dplyr::left_join(name_map_lookup, by = "ModuleID", suffix = c("", "_name_map")) |>
      dplyr::filter(nzchar(.data$supermodule_id)) |>
      dplyr::mutate(
        module_key = dplyr::coalesce(.data$ModuleID, .data$module_eigengene),
        module_label = dplyr::coalesce(.data$ModuleLabel_Final, .data$ModuleLabel_Final_name_map, .data$module_key)
      ) |>
      dplyr::distinct(.data$supermodule_id, .data$module_key, .keep_all = TRUE)
  } else {
    ann_module_meta |>
      dplyr::mutate(
        module_key = dplyr::coalesce(.data$ModuleID, .data$module_eigengene),
        module_label = dplyr::coalesce(.data$ModuleLabel_Final, .data$module_key)
      )
  }

  go_terms_for <- function(module_ids, ontology, max_n = 5L) {
    if (is.null(go_enrichment) || !nrow(go_enrichment) || !all(c("ModuleID", "Ontology", "Description") %in% names(go_enrichment))) {
      return(NA_character_)
    }
    out <- go_enrichment |>
      dplyr::filter(.data$ModuleID %in% module_ids, .data$Ontology == ontology) |>
      dplyr::mutate(p_adjust_num = if ("p.adjust" %in% names(go_enrichment)) suppressWarnings(as.numeric(.data[["p.adjust"]])) else NA_real_) |>
      dplyr::arrange(.data$p_adjust_num, .data$Description) |>
      dplyr::distinct(.data$Description, .keep_all = TRUE) |>
      dplyr::slice_head(n = max_n) |>
      dplyr::pull(.data$Description)
    collapse_unique_values(out, max_n = max_n)
  }

  module_summary_hubs <- if (!is.null(module_summary) && nrow(module_summary) && "ModuleID" %in% names(module_summary)) {
    module_summary |>
      dplyr::transmute(ModuleID = as.character(.data$ModuleID), top_hub_proteins_module = col_if_present(module_summary, "top_hub_proteins"))
  } else {
    tibble::tibble(ModuleID = character(), top_hub_proteins_module = character())
  }

  module_long_ranked <- if (!is.null(modules_long) && nrow(modules_long) && "ModuleID" %in% names(modules_long)) {
    modules_long |>
      dplyr::mutate(
        ModuleID = as.character(.data$ModuleID),
        abs_kME_num = suppressWarnings(as.numeric(col_if_present(modules_long, "abs_kME", NA_real_))),
        is_core = as.character(col_if_present(modules_long, "is_core_kME_0.6", FALSE)) %in% c("TRUE", "true", "1"),
        is_hub = as.character(col_if_present(modules_long, "is_top_hub_25", FALSE)) %in% c("TRUE", "true", "1"),
        protein_label = dplyr::coalesce(
          col_if_present(modules_long, "GeneSymbol"),
          col_if_present(modules_long, "ProteinID"),
          col_if_present(modules_long, "UniProt"),
          col_if_present(modules_long, "resolved_uniprot")
        )
      )
  } else {
    tibble::tibble(ModuleID = character(), abs_kME_num = numeric(), is_core = logical(), is_hub = logical(), protein_label = character())
  }

  rows <- lapply(selected$supermodule_id, function(sid) {
    sel <- selected[selected$supermodule_id == sid, , drop = FALSE][1, , drop = FALSE]
    members <- comp_long |> dplyr::filter(.data$supermodule_id == sid)
    member_ids <- unique(members$module_key)
    member_ids <- member_ids[!is.na(member_ids) & nzchar(member_ids)]
    member_labels <- unique(members$module_label)
    member_labels <- member_labels[!is.na(member_labels) & nzchar(member_labels)]
    module_hubs <- module_summary_hubs |> dplyr::filter(.data$ModuleID %in% member_ids)
    long_hits <- module_long_ranked |> dplyr::filter(.data$ModuleID %in% member_ids)
    core <- long_hits |>
      dplyr::filter(.data$is_core %in% TRUE) |>
      dplyr::arrange(dplyr::desc(.data$abs_kME_num)) |>
      dplyr::pull(.data$protein_label)
    if (!length(core)) {
      core <- long_hits |>
        dplyr::arrange(dplyr::desc(.data$abs_kME_num)) |>
        dplyr::slice_head(n = 12) |>
        dplyr::pull(.data$protein_label)
    }
    sm <- summary_meta |> dplyr::filter(.data$supermodule_id == sid)
    pc1_value <- pc1 |>
      dplyr::filter(.data$supermodule_id == sid) |>
      dplyr::pull(.data$pca_PC1_variance_explained) |>
      first_value(NA_real_)
    tibble::tibble(
      dataset = dataset,
      supermodule_id = sid,
      supermodule_label = dplyr::coalesce(first_value(sel$supermodule_label, NA_character_), first_value(sm$summary_label, NA_character_), sid),
      contrast = sel$contrast[[1]],
      estimate = sel$estimate[[1]],
      p_value = sel$p_value[[1]],
      FDR_within_dataset_level = sel$FDR_within_dataset_level[[1]],
      FDR_global = sel$FDR_global[[1]],
      evidence_status = sel$evidence_status[[1]],
      pca_PC1_variance_explained = pc1_value,
      n_member_modules = if (length(member_ids)) length(member_ids) else first_value(sm$summary_n_member_modules, NA_integer_),
      member_modules = if (length(member_ids)) paste(member_ids, collapse = ";") else first_value(sm$summary_member_modules, NA_character_),
      member_module_labels = collapse_unique_values(member_labels, max_n = 16L),
      top_GO_BP_terms = dplyr::coalesce(first_value(sm$top_GO_BP_terms, NA_character_), go_terms_for(member_ids, "BP")),
      top_GO_MF_terms = dplyr::coalesce(first_value(sm$top_GO_MF_terms, NA_character_), go_terms_for(member_ids, "MF")),
      top_GO_CC_terms = dplyr::coalesce(first_value(sm$top_GO_CC_terms, NA_character_), go_terms_for(member_ids, "CC")),
      top_hub_symbols = first_value(sm$top_hub_symbols, NA_character_),
      top_hub_proteins = dplyr::coalesce(first_value(sm$top_hub_proteins_summary, NA_character_), collapse_unique_values(module_hubs$top_hub_proteins_module, max_n = 20L)),
      top_core_kME_proteins = collapse_unique_values(core, max_n = 20L),
      selection_support = sel$selection_support[[1]],
      selection_message = sel$selection_message[[1]]
    )
  })
  dplyr::bind_rows(rows) |>
    dplyr::select(dplyr::all_of(names(empty)))
}

all_supermodule_contents_summary <- function(pca_eigenvalues, comp, definitions, super_ann, super_summary, modules_long, module_summary, go_enrichment, module_name_map, dataset) {
  pc1 <- pca_eigenvalues |>
    dplyr::filter(.data$pc == 1L) |>
    dplyr::select("supermodule_id", pca_PC1_variance_explained = "variance_explained")
  selected_like <- if (!is.null(comp) && nrow(comp) && "supermodule_id" %in% names(comp)) {
    comp |>
      dplyr::transmute(
        dataset = dataset,
        supermodule_id = as.character(.data$supermodule_id),
        supermodule_label = dplyr::coalesce(col_if_present(comp, "supermodule_label"), as.character(.data$supermodule_id)),
        contrast = NA_character_,
        estimate = NA_real_,
        p_value = NA_real_,
        FDR_within_dataset_level = NA_real_,
        FDR_global = NA_real_,
        evidence_status = "not_applicable_composition_summary",
        selection_support = "all_supermodules",
        selection_rank = dplyr::row_number(),
        selection_message = "all_supermodules"
      )
  } else if (!is.null(super_summary) && nrow(super_summary) && "SupermoduleID" %in% names(super_summary)) {
    super_summary |>
      dplyr::transmute(
        dataset = dataset,
        supermodule_id = as.character(.data$SupermoduleID),
        supermodule_label = dplyr::coalesce(col_if_present(super_summary, "Supermodule_DisplayLabel"), col_if_present(super_summary, "Supermodule_FinalLabel"), as.character(.data$SupermoduleID)),
        contrast = NA_character_,
        estimate = NA_real_,
        p_value = NA_real_,
        FDR_within_dataset_level = NA_real_,
        FDR_global = NA_real_,
        evidence_status = "not_applicable_composition_summary",
        selection_support = "all_supermodules",
        selection_rank = dplyr::row_number(),
        selection_message = "all_supermodules"
      )
  } else {
    tibble::tibble(
      dataset = character(), supermodule_id = character(), supermodule_label = character(),
      contrast = character(), estimate = numeric(), p_value = numeric(),
      FDR_within_dataset_level = numeric(), FDR_global = numeric(),
      evidence_status = character(), selection_support = character(),
      selection_rank = integer(), selection_message = character()
    )
  }
  out <- selected_sus_res_supermodule_contents(
    selected_like,
    pca_eigenvalues,
    comp,
    definitions,
    super_ann,
    super_summary,
    modules_long,
    module_summary,
    go_enrichment,
    module_name_map,
    dataset
  )
  out |>
    dplyr::select(-dplyr::any_of(c("selection_support", "selection_message"))) |>
    dplyr::left_join(pc1, by = "supermodule_id", suffix = c("", "_from_pc1")) |>
    dplyr::mutate(pca_PC1_variance_explained = dplyr::coalesce(.data$pca_PC1_variance_explained, .data$pca_PC1_variance_explained_from_pc1)) |>
    dplyr::select(-dplyr::any_of("pca_PC1_variance_explained_from_pc1"))
}

selected_sus_res_supermodule_interpretation_summary <- function(contents, selected_loadings = NULL) {
  if (is.null(contents) || !nrow(contents)) {
    return(tibble::tibble(
      dataset = character(), supermodule_id = character(), supermodule_label = character(),
      selection_support = character(), interpretation_sentence = character()
    ))
  }
  loading_summary <- if (!is.null(selected_loadings) && nrow(selected_loadings)) {
    selected_loadings |>
      dplyr::filter(is.finite(.data$abs_loading)) |>
      dplyr::arrange(.data$supermodule_id, dplyr::desc(.data$abs_loading), .data$module_label) |>
      dplyr::group_by(.data$supermodule_id) |>
      dplyr::summarise(
        strongest_contributing_modules = collapse_unique_values(dplyr::coalesce(.data$module_label, .data$module_id, .data$module_eigengene), max_n = 3L),
        .groups = "drop"
      )
  } else {
    tibble::tibble(supermodule_id = character(), strongest_contributing_modules = character())
  }
  contents |>
    dplyr::left_join(loading_summary, by = "supermodule_id") |>
    dplyr::mutate(
      support_phrase = dplyr::case_when(
        !is.na(.data$FDR_within_dataset_level) & .data$FDR_within_dataset_level <= 0.05 ~ "FDR-significant",
        !is.na(.data$FDR_within_dataset_level) & .data$FDR_within_dataset_level <= 0.10 ~ "FDR-suggestive",
        !is.na(.data$p_value) & .data$p_value < 0.05 ~ "nominal-only",
        TRUE ~ "not significant"
      ),
      direction_phrase = dplyr::case_when(
        !is.na(.data$estimate) & .data$estimate > 0 ~ "higher in SUS than RES",
        !is.na(.data$estimate) & .data$estimate < 0 ~ "lower in SUS than RES",
        TRUE ~ "different between SUS and RES"
      ),
      pc1_pct = ifelse(is.na(.data$pca_PC1_variance_explained), "NA", paste0(round(100 * .data$pca_PC1_variance_explained), "%")),
      interpretation_sentence = paste0(
        .data$supermodule_id, " (", .data$supermodule_label, ") is a ", .data$support_phrase,
        " SUS-RES supermodule eigengene effect (estimate=", signif(.data$estimate, 3),
        ", p=", signif(.data$p_value, 3), ", FDR=", signif(.data$FDR_within_dataset_level, 3),
        "), with the eigengene ", .data$direction_phrase, ". PC1 explains ", .data$pc1_pct,
        " of member-module variance. The strongest contributing modules are ",
        dplyr::coalesce(.data$strongest_contributing_modules, .data$member_module_labels, "not available"),
        ". Top biology: ", dplyr::coalesce(.data$top_GO_BP_terms, "not available"),
        ". Top hubs/core proteins: ",
        dplyr::coalesce(.data$top_hub_symbols, .data$top_core_kME_proteins, "not available"), "."
      )
    ) |>
    dplyr::select(
      "dataset", "supermodule_id", "supermodule_label", "contrast",
      "estimate", "p_value", "FDR_within_dataset_level", "FDR_global",
      "evidence_status", "selection_support", "pca_PC1_variance_explained",
      "n_member_modules", "member_module_labels", "strongest_contributing_modules", "top_GO_BP_terms",
      "top_hub_symbols", "top_core_kME_proteins", "interpretation_sentence"
    )
}

shorten_terms <- function(x, max_chars = 95L) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "not available"
  vapply(x, function(z) {
    first_terms <- strsplit(z, ";", fixed = TRUE)[[1]]
    first_terms <- trimws(first_terms)
    out <- paste(utils::head(first_terms[nzchar(first_terms)], 2), collapse = "; ")
    if (!nzchar(out)) out <- "not available"
    if (nchar(out) > max_chars) paste0(substr(out, 1, max_chars - 3L), "...") else out
  }, character(1))
}

wrap_plot_label <- function(x, width = 36L) {
  vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

plot_selected_sus_res_supermodule_effect_summary <- function(contents, path) {
  if (is.null(contents) || !nrow(contents)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No SUS-RES significant supermodule eigengenes selected", size = 3) +
      ggplot2::labs(title = "Selected SUS-RES supermodule effects") +
      ggplot2::theme_void()
  } else {
    plot_df <- contents |>
      dplyr::mutate(
        support_class = dplyr::case_when(
          !is.na(.data$FDR_within_dataset_level) & .data$FDR_within_dataset_level <= 0.05 ~ "FDR <= 0.05",
          !is.na(.data$FDR_within_dataset_level) & .data$FDR_within_dataset_level <= 0.10 ~ "FDR <= 0.10",
          !is.na(.data$p_value) & .data$p_value < 0.05 ~ "nominal-only p < 0.05",
          TRUE ~ "not significant"
        ),
        label = paste0(.data$supermodule_id, "\n", wrap_plot_label(.data$supermodule_label, width = 28L)),
        stat_label = paste0("p=", signif(.data$p_value, 2), "\nFDR=", signif(.data$FDR_within_dataset_level, 2))
      ) |>
      dplyr::arrange(.data$estimate)
    plot_df$label <- factor(plot_df$label, levels = plot_df$label)
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$label, y = .data$estimate, fill = .data$support_class)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.25) +
      ggplot2::geom_col(width = 0.64, color = "grey25", linewidth = 0.2) +
      ggplot2::geom_text(ggplot2::aes(label = .data$stat_label), hjust = ifelse(plot_df$estimate < 0, 1.05, -0.05), size = 2.4) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::scale_fill_manual(values = c("FDR <= 0.05" = "#0B6E4F", "FDR <= 0.10" = "#65A30D", "nominal-only p < 0.05" = "#F59E0B", "not significant" = "#BDBDBD"), drop = FALSE) +
      ggplot2::labs(x = NULL, y = "SUS-RES eigengene estimate", fill = "Selection support", title = "Selected SUS-RES supermodule eigengene effects") +
      ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 6), plot.margin = ggplot2::margin(5.5, 36, 5.5, 5.5))
  }
  ggplot2::ggsave(path, p, width = 185, height = 105, units = "mm", device = svglite::svglite)
}

plot_selected_sus_res_supermodule_contents_overview <- function(contents, path) {
  if (is.null(contents) || !nrow(contents)) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No SUS-RES significant supermodule eigengenes selected", size = 3) +
      ggplot2::labs(title = "Selected SUS-RES supermodule contents") +
      ggplot2::theme_void()
  } else {
    plot_df <- contents |>
      dplyr::mutate(
        label = paste0(.data$supermodule_id, "\n", wrap_plot_label(.data$supermodule_label, width = 28L)),
        pc1_pct = pmax(0, pmin(1, .data$pca_PC1_variance_explained)),
        pc1_label = paste0(round(100 * .data$pc1_pct), "% PC1"),
        biology_label = wrap_plot_label(paste0("Modules: ", .data$member_module_labels, "\nGO: ", shorten_terms(.data$top_GO_BP_terms)), width = 58L)
      ) |>
      dplyr::arrange(dplyr::desc(.data$pc1_pct))
    plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$label, y = .data$pc1_pct)) +
      ggplot2::geom_col(fill = "#4E79A7", width = 0.62) +
      ggplot2::geom_text(ggplot2::aes(label = .data$pc1_label), hjust = -0.08, size = 2.5) +
      ggplot2::geom_text(ggplot2::aes(y = 0.02, label = .data$biology_label), hjust = 0, size = 2.35, color = "grey15") +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1.22)) +
      ggplot2::labs(x = NULL, y = "PC1 variance explained", title = "Selected SUS-RES supermodule contents") +
      ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 6), plot.margin = ggplot2::margin(5.5, 34, 5.5, 5.5))
  }
  ggplot2::ggsave(path, p, width = 205, height = 115, units = "mm", device = svglite::svglite)
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

definitions <- safe_read_csv(FILES$definitions) %||% data.frame()
super_ann <- safe_read_csv(FILES$supermodule_annotation) %||% data.frame()
module_summary <- safe_read_csv(FILES$module_summary) %||% data.frame()
super_summary <- safe_read_csv(FILES$supermodule_summary) %||% data.frame()
modules_long <- safe_read_csv(file.path(dirname(FILES$module_summary), "WGCNA_modules_long.csv")) %||% data.frame()
go_enrichment <- safe_read_csv(FILES$go) %||% data.frame()
module_name_map <- safe_read_csv(file.path(dirname(FILES$module_summary), "module_name_map.csv")) %||% data.frame()

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
  super_pca_member_loadings <- supermodule_pca_member_loadings(data.frame(Sample = character()), data.frame(), DATASET)
  all_supermodule_eigengene_group_values_tbl <- empty_supermodule_eigengene_group_values()
  all_supermodule_eigengene_spatial_group_values_tbl <- empty_supermodule_eigengene_group_values()
} else {
  module_eig <- extract_module_eigengenes(state)
  maps <- make_endpoint_maps(module_eig, definitions, super_ann)
  maps$super_map <- reconcile_supermodule_map_to_module_eigengenes(maps$super_map, module_eig, definitions)
  super_pca_input_audit <- supermodule_pca_input_audit(module_eig, maps$super_map, DATASET)
  super_pca_eigenvalues <- supermodule_pca_eigenvalues(module_eig, maps$super_map, DATASET)
  super_pca_member_loadings <- supermodule_pca_member_loadings(module_eig, maps$super_map, DATASET)
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
  all_supermodule_eigengene_group_values_tbl <- all_supermodule_eigengene_group_values(super$eigengenes, super_endpoint_map, state, DATASET, SpatialUnitType)
  all_supermodule_eigengene_spatial_group_values_tbl <- all_supermodule_eigengene_group_values_tbl

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
selected_sus_res_contents <- selected_sus_res_supermodule_contents(
  selected_sus_res_audit,
  super_pca_eigenvalues,
  comp,
  definitions,
  super_ann,
  super_summary,
  modules_long,
  module_summary,
  go_enrichment,
  module_name_map,
  DATASET
)
all_supermodule_contents <- all_supermodule_contents_summary(
  super_pca_eigenvalues,
  comp,
  definitions,
  super_ann,
  super_summary,
  modules_long,
  module_summary,
  go_enrichment,
  module_name_map,
  DATASET
)
selected_sus_res_loadings <- selected_sus_res_member_loadings(super_pca_member_loadings, selected_sus_res_audit, selected_sus_res_contents)
selected_sus_res_interpretation_summary <- selected_sus_res_supermodule_interpretation_summary(selected_sus_res_contents, selected_sus_res_loadings)
all_supermodule_region_layer_effects <- all_supermodule_region_layer_effect_source(super_out)
all_supermodule_region_layer_cohend_effects <- publication_supermodule_region_layer_effect_source(DATASET, all_supermodule_contents)
all_supermodule_region_layer_heatmap_source <- if (nrow(all_supermodule_region_layer_cohend_effects)) {
  all_supermodule_region_layer_cohend_effects
} else {
  all_supermodule_region_layer_effects
}
all_sus_res_region_layer_effects <- all_supermodule_region_layer_heatmap_source |>
  dplyr::filter(.data$contrast_block == "SUS-RES") |>
  dplyr::mutate(contrast = "SUS - RES")

write_table_and_source(module_out, PATHS$tables, PATHS$source_data, "module_group_effects.csv")
write_table_and_source(super_out, PATHS$tables, PATHS$source_data, "supermodule_group_effects.csv")
write_table_and_source(comp, PATHS$tables, PATHS$source_data, "supermodule_composition.csv")
write_table_and_source(super_pca_input_audit, PATHS$tables, PATHS$source_data, "supermodule_pca_input_audit.csv")
write_table_and_source(super_pca_eigenvalues, PATHS$tables, PATHS$source_data, "supermodule_pca_eigenvalues.csv")
write_table_and_source(super_pca_member_loadings, PATHS$tables, PATHS$source_data, "supermodule_pca_member_loadings.csv")
write_table_and_source(selected_sus_res_audit, PATHS$tables, PATHS$source_data, "selected_sus_res_supermodule_pca_selection_audit.csv")
write_table_and_source(selected_sus_res_pca, PATHS$tables, PATHS$source_data, "selected_sus_res_supermodule_pca_eigenvalues.csv")
write_table_and_source(selected_sus_res_contents, PATHS$tables, PATHS$source_data, "selected_sus_res_supermodule_contents.csv")
write_table_and_source(all_supermodule_contents, PATHS$tables, PATHS$source_data, "all_supermodule_contents_summary.csv")
write_table_and_source(selected_sus_res_interpretation_summary, PATHS$tables, PATHS$source_data, "selected_sus_res_supermodule_interpretation_summary.csv")
write_table_and_source(selected_sus_res_interpretation_summary, PATHS$tables, PATHS$source_data, "selected_sus_res_supermodule_interpretation.csv")
write_table_and_source(all_supermodule_eigengene_group_values_tbl, PATHS$tables, PATHS$source_data, "all_supermodule_eigengene_group_values.csv")
write_table_and_source(all_supermodule_eigengene_spatial_group_values_tbl, PATHS$tables, PATHS$source_data, "all_supermodule_eigengene_spatial_group_values.csv")
write_table_and_source(all_supermodule_region_layer_effects, PATHS$tables, PATHS$source_data, "all_supermodule_region_layer_effects.csv")
write_table_and_source(all_supermodule_region_layer_cohend_effects, PATHS$tables, PATHS$source_data, "all_supermodule_region_layer_cohend_effects.csv")
write_table_and_source(all_sus_res_region_layer_effects, PATHS$tables, PATHS$source_data, "all_sus_res_supermodule_region_layer_effects.csv")
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
plot_all_supermodule_eigengene_group(all_supermodule_eigengene_group_values_tbl, file.path(PATHS$figures, "all_supermodule_eigengene_group_plot.svg"))
plot_all_supermodule_eigengene_spatial_group(all_supermodule_eigengene_spatial_group_values_tbl, file.path(PATHS$figures, "all_supermodule_eigengene_spatial_group_plot.svg"))
plot_selected_sus_res_pca_scree(selected_sus_res_pca, selected_sus_res_audit, file.path(PATHS$figures, "selected_sus_res_supermodule_pca_eigenvalue_scree.svg"))
plot_all_supermodule_member_loadings(super_pca_member_loadings, file.path(PATHS$figures, "all_supermodule_member_loading_plot.svg"))
plot_selected_sus_res_member_loadings(selected_sus_res_loadings, file.path(PATHS$figures, "selected_sus_res_supermodule_member_loading_plot.svg"))
plot_all_supermodule_region_layer_heatmap(all_supermodule_region_layer_heatmap_source, file.path(PATHS$figures, "all_supermodule_region_layer_effect_heatmap.svg"))
plot_all_sus_res_region_layer_effects(all_sus_res_region_layer_effects, file.path(PATHS$figures, "all_sus_res_supermodule_region_layer_effects.svg"))
plot_selected_sus_res_supermodule_effect_summary(selected_sus_res_contents, file.path(PATHS$figures, "selected_sus_res_supermodule_effect_summary.svg"))
plot_selected_sus_res_supermodule_contents_overview(selected_sus_res_contents, file.path(PATHS$figures, "selected_sus_res_supermodule_contents_overview.svg"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = FILES,
  outputs = list(tables = PATHS$tables, source_data = PATHS$source_data, figures = PATHS$figures),
  parameters = list(dataset = DATASET, level = LEVEL, effect_scopes = c("within_spatial_unit", "spatial_adjusted_global", "stress_by_spatial_interaction")),
  notes = "Primary results are not filtered by marker or microenvironment class. Marker traits, if present, are annotation only."
)

message("WGCNA module/supermodule group effects complete for dataset: ", DATASET)
