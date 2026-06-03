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

default_formula <- switch(
  DATASET,
  neuron_neuropil = "eigengene ~ StressGroup * RegionLayer + Sex + Batch",
  neuron_soma = "eigengene ~ StressGroup * Region + Sex + Batch",
  microglia = "eigengene ~ StressGroup * Region + Sex + Batch"
)
formula_override <- Sys.getenv("PROTEOMICS_WGCNA_GROUP_FORMULA", unset = "")
global_formula_text <- if (nzchar(formula_override)) formula_override else default_formula

if (run$dry_run) {
  invisible(lapply(unlist(PATHS), dir_create))
  dry_run_line("Script", "06_modules_WGCNA/05_module_supermodule_group_effects.r")
  dry_run_line("Dataset", DATASET)
  dry_run_line("Level", LEVEL)
  dry_run_line("WGCNA final state", FILES$state, if (file.exists(FILES$state)) "PASS" else "WARN")
  dry_run_line("Module definitions", FILES$definitions, if (file.exists(FILES$definitions)) "PASS" else "WARN")
  dry_run_line("Supermodule annotation", FILES$supermodule_annotation, if (file.exists(FILES$supermodule_annotation)) "PASS" else "WARN")
  dry_run_line("Optional marker traits", FILES$marker_traits, if (file.exists(FILES$marker_traits)) "PASS" else "WARN")
  dry_run_line("Formula", global_formula_text)
  dry_run_line("Output tables", PATHS$tables)
  quit(status = 0, save = "no")
}

if (!length(missing_pkgs)) {
  theme_set(theme_classic(base_size = 8))
}

fit_one_endpoint <- function(df, endpoint_col, endpoint_id, endpoint_label, level, spatial_label = "global", formula_text) {
  dat <- df |>
    dplyr::mutate(eigengene = as.numeric(.data[[endpoint_col]])) |>
    dplyr::filter(is.finite(.data$eigengene), !is.na(.data$StressGroup), .data$StressGroup %in% c("CON", "RES", "SUS"))
  if (nrow(dat) < 4L || dplyr::n_distinct(dat$StressGroup) < 2L) {
    out <- empty_group_effects(DATASET, level, "too few samples/groups")
    out$spatial_unit <- spatial_label
    if (level == "module") out$module_id <- endpoint_id else out$supermodule_id <- endpoint_id
    return(out)
  }
  dat$StressGroup <- factor(dat$StressGroup, levels = c("CON", "RES", "SUS"))
  vars_in_formula <- all.vars(stats::as.formula(formula_text))
  vars_in_formula <- setdiff(vars_in_formula, "eigengene")
  dropped <- vars_in_formula[!vars_in_formula %in% names(dat)]
  keep_vars <- vars_in_formula[vars_in_formula %in% names(dat)]
  invalid <- keep_vars[vapply(keep_vars, function(v) dplyr::n_distinct(dat[[v]][!is.na(dat[[v]])]) < 2L, logical(1))]
  dropped <- unique(c(dropped, invalid))
  rhs_vars <- setdiff(keep_vars, dropped)
  if (!"StressGroup" %in% rhs_vars) rhs_vars <- c("StressGroup", rhs_vars)
  rhs_vars <- unique(rhs_vars)
  simple_formula <- stats::as.formula(paste("eigengene ~", paste(rhs_vars, collapse = " + ")))
  warning_text <- NA_character_
  fit <- tryCatch(
    stats::lm(simple_formula, data = dat),
    warning = function(w) {
      warning_text <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    },
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    out <- empty_group_effects(DATASET, level, conditionMessage(fit))
    out$spatial_unit <- spatial_label
    if (level == "module") out$module_id <- endpoint_id else out$supermodule_id <- endpoint_id
    return(out)
  }
  rank_deficient <- fit$rank < ncol(stats::model.matrix(fit))
  formula_used <- paste(deparse(stats::formula(fit)), collapse = "")

  contrasts <- c("RES - CON", "SUS - CON", "RES - SUS")
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emm <- tryCatch(emmeans::emmeans(fit, specs = "StressGroup"), error = function(e) NULL)
    if (!is.null(emm)) {
      contr <- as.data.frame(emmeans::contrast(emm, method = list(
        "RES - CON" = c(-1, 1, 0),
        "SUS - CON" = c(-1, 0, 1),
        "RES - SUS" = c(0, -1, 1)
      )))
      names(contr)[names(contr) == "emmean"] <- "estimate"
      return(tibble::tibble(
        dataset = DATASET,
        level = level,
        module_id = if (level == "module") endpoint_id else NA_character_,
        supermodule_id = if (level == "supermodule") endpoint_id else NA_character_,
        module_label = if (level == "module") endpoint_label else NA_character_,
        supermodule_label = if (level == "supermodule") endpoint_label else NA_character_,
        spatial_unit = spatial_label,
        contrast = as.character(contr$contrast),
        estimate = as.numeric(contr$estimate),
        SE = as.numeric(contr$SE),
        statistic = as.numeric(contr$t.ratio %||% contr$z.ratio %||% NA_real_),
        p_value = as.numeric(contr$p.value),
        FDR_within_dataset_level = NA_real_,
        FDR_global = NA_real_,
        direction = dplyr::case_when(estimate > 0 ~ "higher", estimate < 0 ~ "lower", TRUE ~ "zero"),
        n_samples = nrow(dat),
        formula_used = formula_used,
        dropped_covariates = paste(dropped, collapse = ";"),
        rank_deficient_model = rank_deficient,
        model_warning = warning_text
      ))
    }
  }

  means <- stats::aggregate(dat$eigengene, list(StressGroup = dat$StressGroup), mean, na.rm = TRUE)
  names(means)[2] <- "mean"
  rows <- lapply(contrasts, function(contrast) {
    parts <- strsplit(contrast, " - ", fixed = TRUE)[[1]]
    sub <- dat[dat$StressGroup %in% parts, , drop = FALSE]
    tt <- tryCatch(stats::t.test(eigengene ~ StressGroup, data = sub), error = function(e) NULL)
    est <- means$mean[match(parts[[1]], means$StressGroup)] - means$mean[match(parts[[2]], means$StressGroup)]
    tibble::tibble(
      dataset = DATASET, level = level,
      module_id = if (level == "module") endpoint_id else NA_character_,
      supermodule_id = if (level == "supermodule") endpoint_id else NA_character_,
      module_label = if (level == "module") endpoint_label else NA_character_,
      supermodule_label = if (level == "supermodule") endpoint_label else NA_character_,
      spatial_unit = spatial_label, contrast = contrast, estimate = est,
      SE = if (!is.null(tt)) unname(diff(tt$conf.int)) / (2 * 1.96) else NA_real_,
      statistic = if (!is.null(tt)) unname(tt$statistic) else NA_real_,
      p_value = if (!is.null(tt)) tt$p.value else NA_real_,
      FDR_within_dataset_level = NA_real_, FDR_global = NA_real_,
      direction = dplyr::case_when(est > 0 ~ "higher", est < 0 ~ "lower", TRUE ~ "zero"),
      n_samples = nrow(dat), formula_used = formula_used,
      dropped_covariates = paste(dropped, collapse = ";"),
      rank_deficient_model = rank_deficient,
      model_warning = paste(stats::na.omit(c(warning_text, "emmeans unavailable; used two-group t-test contrasts")), collapse = "; ")
    )
  })
  dplyr::bind_rows(rows)
}

run_effects <- function(eigengenes, endpoint_map, level, formula_text) {
  meta <- standardize_wgcna_metadata(state$sample_info, DATASET)
  dat <- dplyr::inner_join(meta, eigengenes, by = "Sample")
  endpoint_cols <- intersect(endpoint_map$endpoint_col, names(dat))
  spatial_col <- if (DATASET == "neuron_neuropil" && "RegionLayer" %in% names(dat) && any(!is.na(dat$RegionLayer))) "RegionLayer" else "Region"
  within <- lapply(endpoint_cols, function(col) {
    map <- endpoint_map[match(col, endpoint_map$endpoint_col), , drop = FALSE]
    dplyr::bind_rows(lapply(sort(unique(stats::na.omit(dat[[spatial_col]]))), function(unit) {
      fit_one_endpoint(dat[dat[[spatial_col]] == unit, , drop = FALSE], col, map$endpoint_id, map$endpoint_label, level, unit, "eigengene ~ StressGroup + Sex + Batch")
    }))
  }) |> dplyr::bind_rows()
  global <- lapply(endpoint_cols, function(col) {
    map <- endpoint_map[match(col, endpoint_map$endpoint_col), , drop = FALSE]
    fit_one_endpoint(dat, col, map$endpoint_id, map$endpoint_label, level, "global_spatial_adjusted", formula_text)
  }) |> dplyr::bind_rows()
  out <- dplyr::bind_rows(within, global)
  out$FDR_within_dataset_level <- stats::p.adjust(out$p_value, method = "BH")
  out
}

state <- tryCatch(load_wgcna_state(FILES$state), error = function(e) e)
if (inherits(state, "error")) {
  msg <- conditionMessage(state)
  module_out <- empty_group_effects(DATASET, "module", msg)
  super_out <- empty_group_effects(DATASET, "supermodule", msg)
  comp <- data.frame(dataset = DATASET, supermodule_id = NA_character_, n_member_modules = 0L, member_modules = NA_character_, stringsAsFactors = FALSE)
} else {
  module_eig <- extract_module_eigengenes(state)
  definitions <- safe_read_csv(FILES$definitions)
  super_ann <- safe_read_csv(FILES$supermodule_annotation)
  if (is.null(definitions)) definitions <- data.frame()
  if (is.null(super_ann)) super_ann <- data.frame()

  module_map <- data.frame(endpoint_col = setdiff(names(module_eig), "Sample"), stringsAsFactors = FALSE) |>
    dplyr::mutate(module_eigengene = .data$endpoint_col, ModuleColor = module_col_to_id(.data$endpoint_col))
  if (nrow(definitions)) {
    label_cols <- c("ModuleColor", "ModuleID", "module_eigengene", "ModuleLabel_Final")
    defs <- definitions[, intersect(label_cols, names(definitions)), drop = FALSE] |> dplyr::distinct()
    module_map <- module_map |>
      dplyr::left_join(defs, by = "module_eigengene")
    for (nm in c("ModuleID", "ModuleColor.x", "ModuleColor.y", "ModuleLabel_Final")) {
      if (!nm %in% names(module_map)) module_map[[nm]] <- NA_character_
    }
    module_map <- module_map |>
      dplyr::mutate(
        endpoint_id = dplyr::coalesce(.data$ModuleID, .data$ModuleColor.x, .data$ModuleColor.y, .data$module_eigengene),
        endpoint_label = dplyr::coalesce(.data$ModuleLabel_Final, .data$endpoint_id)
      )
  } else {
    module_map <- module_map |> dplyr::mutate(endpoint_id = .data$ModuleColor, endpoint_label = .data$ModuleColor)
  }

  if (nrow(super_ann)) {
    super_map0 <- super_ann
    if ("present_in_dataset" %in% names(super_map0)) {
      super_map0 <- super_map0 |>
        dplyr::filter(.data$present_in_dataset %in% c(TRUE, "TRUE", "true", 1))
    }
    for (nm in c("Supermodule_DataDriven", "Supermodule", "Supermodule_Manual")) {
      if (!nm %in% names(super_map0)) super_map0[[nm]] <- NA_character_
    }
    super_map <- super_map0 |>
      dplyr::mutate(
        module_eigengene = as.character(.data$module_eigengene),
        SupermoduleID = dplyr::coalesce(as.character(.data$Supermodule_DataDriven), as.character(.data$Supermodule), as.character(.data$Supermodule_Manual)),
        SupermoduleLabel = dplyr::coalesce(as.character(.data$Supermodule), .data$SupermoduleID)
      )
  } else {
    super_map <- data.frame(module_eigengene = character(), SupermoduleID = character(), SupermoduleLabel = character())
  }
  super <- make_supermodule_eigengenes(module_eig, super_map)
  comp <- super$composition |>
    dplyr::mutate(dataset = DATASET, .before = .data$supermodule_id) |>
    dplyr::left_join(
      super_map |> dplyr::distinct(SupermoduleID, SupermoduleLabel),
      by = c("supermodule_id" = "SupermoduleID")
    ) |>
    dplyr::rename(supermodule_label = .data$SupermoduleLabel)

  module_out <- if (LEVEL %in% c("module", "both")) run_effects(module_eig, module_map, "module", global_formula_text) else empty_group_effects(DATASET, "module", "not requested")
  super_endpoint_map <- comp |>
    dplyr::transmute(endpoint_col = .data$supermodule_eigengene, endpoint_id = .data$supermodule_id, endpoint_label = .data$supermodule_label)
  super_out <- if (LEVEL %in% c("supermodule", "both") && nrow(super_endpoint_map)) run_effects(super$eigengenes, super_endpoint_map, "supermodule", global_formula_text) else empty_group_effects(DATASET, "supermodule", "not requested or no supermodules")
}

all_p <- c(module_out$p_value, super_out$p_value)
all_fdr <- stats::p.adjust(all_p, method = "BH")
module_out$FDR_global <- all_fdr[seq_len(nrow(module_out))]
super_out$FDR_global <- all_fdr[seq_len(nrow(super_out)) + nrow(module_out)]
module_out <- module_out[, required_group_effect_columns]
super_out <- super_out[, required_group_effect_columns]

write_table_and_source(module_out, PATHS$tables, PATHS$source_data, "module_group_effects.csv")
write_table_and_source(super_out, PATHS$tables, PATHS$source_data, "supermodule_group_effects.csv")
write_table_and_source(comp, PATHS$tables, PATHS$source_data, "supermodule_composition.csv")

map_out <- safe_read_csv(FILES$supermodule_annotation)
if (is.null(map_out)) map_out <- data.frame(dataset = DATASET, status = "missing_supermodule_annotation")
write_table_and_source(map_out, PATHS$tables, PATHS$source_data, "module_to_supermodule_map_with_annotations.csv")

plot_effects <- function(df, path, title) {
  plot_df <- df |> dplyr::filter(!is.na(.data$p_value), .data$spatial_unit != "global_spatial_adjusted")
  if (!nrow(plot_df)) return(invisible(NULL))
  label_col <- if (unique(plot_df$level)[[1]] == "module") "module_id" else "supermodule_id"
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$contrast, y = .data[[label_col]], color = .data$estimate, size = -log10(.data$p_value))) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_wrap(~ spatial_unit, scales = "free_y") +
    ggplot2::scale_color_gradient2(low = "#3B6FB6", mid = "grey92", high = "#C84C5A") +
    ggplot2::labs(x = NULL, y = NULL, color = "Estimate", size = "-log10 P", title = title) +
    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(path, p, width = 170, height = 120, units = "mm", device = svglite::svglite)
}
plot_effects(module_out, file.path(PATHS$figures, "module_group_effect_dotplot.svg"), "Module group effects")
plot_effects(super_out, file.path(PATHS$figures, "supermodule_effects_by_spatial_unit.svg"), "Supermodule group effects")
plot_effects(super_out |> dplyr::filter(.data$spatial_unit == "global_spatial_adjusted"), file.path(PATHS$figures, "supermodule_group_effect_heatmap.svg"), "Spatial-adjusted supermodule effects")

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = FILES,
  outputs = list(tables = PATHS$tables, source_data = PATHS$source_data, figures = PATHS$figures),
  parameters = list(dataset = DATASET, level = LEVEL, formula = global_formula_text),
  notes = "Primary results are not filtered by marker or microenvironment class. Marker traits, if present, are annotation only."
)

message("WGCNA module/supermodule group effects complete for dataset: ", DATASET)
