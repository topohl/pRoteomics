#!/usr/bin/env Rscript
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
  quit(status = 0, save = "no")
}

if (!length(missing_pkgs)) theme_set(theme_classic(base_size = 8))

has_repeats <- function(dat) {
  "AnimalID" %in% names(dat) &&
    any(!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))) &&
    any(table(dat$AnimalID[!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))]) > 1L)
}

status_row <- function(level, endpoint_id, endpoint_label, spatial_unit, effect_scope, reason,
                       formula_requested = NA_character_, formula_used = NA_character_,
                       model_type = NA_character_, dat = NULL) {
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
    model_warning = reason
  )
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
  list(
    fit = fit,
    warning = warning_text,
    model_type = if (use_lmer) "lmerTest_lmer" else "lm",
    formula_requested = formula_requested,
    formula_used = paste(deparse(stats::formula(fit)), collapse = ""),
    rank_deficient = rank_deficient
  )
}

contrast_rows <- function(fit_info, dat, level, endpoint_id, endpoint_label, spatial_unit, effect_scope, by_spatial = FALSE, dropped = character()) {
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
        model_warning = paste(fit_info$warning, collapse = "; ")
      ))
    }
  }

  warning_note <- c(fit_info$warning, "emmeans unavailable or failed; used two-group t-test contrasts")
  contrasts <- c("RES - CON", "SUS - CON", "SUS - RES")
  units <- if (by_spatial) sort(unique(stats::na.omit(dat$SpatialLabel))) else spatial_unit
  dplyr::bind_rows(lapply(units, function(unit) {
    subdat <- if (by_spatial) dat[dat$SpatialLabel == unit, , drop = FALSE] else dat
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
        model_type = fit_info$model_type, has_repeated_animals = has_repeats(dat),
        n_animals = if ("AnimalID" %in% names(dat)) dplyr::n_distinct(dat$AnimalID[!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))]) else NA_integer_,
        contrast = contrast, estimate = est,
        SE = if (!is.null(tt)) unname(diff(tt$conf.int)) / (2 * 1.96) else NA_real_,
        statistic = if (!is.null(tt)) unname(tt$statistic) else NA_real_,
        p_value = if (!is.null(tt)) tt$p.value else NA_real_,
        FDR_within_dataset_level = NA_real_, FDR_global = NA_real_,
        evidence_status = NA_character_,
        direction = dplyr::case_when(est > 0 ~ "higher", est < 0 ~ "lower", TRUE ~ "zero"),
        n_samples = nrow(dat), formula_requested = fit_info$formula_requested,
        formula_used = fit_info$formula_used, dropped_covariates = paste(dropped, collapse = ";"),
        rank_deficient_model = fit_info$rank_deficient, model_warning = paste(warning_note, collapse = "; ")
      )
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
    for (nm in c("Supermodule_DataDrivenID", "Supermodule_DataDriven", "SupermoduleID", "Supermodule_DisplayLabel", "Macroprogram_Display", "Supermodule_LongLabel", "Supermodule_FinalLabel", "Supermodule", "Supermodule_DataDrivenLabel")) if (!nm %in% names(super_map0)) super_map0[[nm]] <- NA_character_
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

correlate_marker_traits <- function(eigengenes, endpoint_map, marker_traits, level) {
  if (is.null(marker_traits) || !nrow(marker_traits)) return(tibble::tibble())
  trait_cols <- grep("^(z|raw)_.+_(score|ratio)$|^z_microglia_minus_neuropil_score$|^z_microglia_to_neuropil_ratio$|^raw_microglia_minus_neuropil_score$|^raw_microglia_to_neuropil_ratio$", names(marker_traits), value = TRUE)
  trait_cols <- trait_cols[vapply(marker_traits[trait_cols], is.numeric, logical(1))]
  dat <- dplyr::inner_join(eigengenes, marker_traits[, c("Sample", trait_cols), drop = FALSE], by = "Sample")
  endpoint_cols <- intersect(endpoint_map$endpoint_col, names(dat))
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
    grepl("rank|singular|not estimable|failed|unavailable|t-test|too few|empty", warn)
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
} else {
  module_eig <- extract_module_eigengenes(state)
  definitions <- safe_read_csv(FILES$definitions) %||% data.frame()
  super_ann <- safe_read_csv(FILES$supermodule_annotation) %||% data.frame()
  maps <- make_endpoint_maps(module_eig, definitions, super_ann)
  super <- make_supermodule_eigengenes(module_eig, maps$super_map)
  comp <- super$composition |>
    dplyr::mutate(dataset = DATASET, .before = "supermodule_id") |>
    dplyr::left_join(maps$super_map |> dplyr::distinct(SupermoduleID, SupermoduleLabel), by = c("supermodule_id" = "SupermoduleID")) |>
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

write_table_and_source(module_out, PATHS$tables, PATHS$source_data, "module_group_effects.csv")
write_table_and_source(super_out, PATHS$tables, PATHS$source_data, "supermodule_group_effects.csv")
write_table_and_source(comp, PATHS$tables, PATHS$source_data, "supermodule_composition.csv")
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
plot_effects(super_out, file.path(PATHS$figures, "supermodule_group_effect_heatmap.svg"), "Supermodule group effects")

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = FILES,
  outputs = list(tables = PATHS$tables, source_data = PATHS$source_data, figures = PATHS$figures),
  parameters = list(dataset = DATASET, level = LEVEL, effect_scopes = c("within_spatial_unit", "spatial_adjusted_global", "stress_by_spatial_interaction")),
  notes = "Primary results are not filtered by marker or microenvironment class. Marker traits, if present, are annotation only."
)

message("WGCNA module/supermodule group effects complete for dataset: ", DATASET)
