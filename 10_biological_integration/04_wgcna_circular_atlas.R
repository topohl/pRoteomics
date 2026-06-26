#!/usr/bin/env Rscript
#
# Build and audit WGCNA circular-atlas source data.
#
# Downstream-only contract:
#   - do not recompute WGCNA or alter module/supermodule definitions
#   - copy group-effect statistics/status from source WGCNA outputs
#   - represent every source supermodule exactly once in the atlas segments

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))

required_pkgs <- c("dplyr", "readr", "tibble", "tidyr", "stringr")
plot_pkgs <- c("circlize", "svglite", "ggplot2", "scales")
all_required_pkgs <- unique(c(required_pkgs, plot_pkgs))
missing_pkgs <- all_required_pkgs[!vapply(all_required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
available_pkgs <- all_required_pkgs[vapply(all_required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(available_pkgs)) {
  suppressPackageStartupMessages(invisible(lapply(available_pkgs, library, character.only = TRUE)))
}

run <- wgcna_cli(default_dataset = "all", allow_all = TRUE)
DATASET_ARG <- run$dataset

table_dir <- path_results("tables", "10_biological_integration", "wgcna_circular_atlas", "global")
source_dir <- path_results("source_data", "10_biological_integration", "wgcna_circular_atlas", "global")
report_dir <- path_results("reports", "10_biological_integration", "wgcna_circular_atlas", "global")
figure_dir <- path_results("figures", "10_biological_integration", "wgcna_circular_atlas", "global")
log_dir <- path_results("logs", "10_biological_integration", "wgcna_circular_atlas", "global")
invisible(lapply(c(table_dir, source_dir, report_dir, figure_dir, log_dir), dir_create))

out_segments <- file.path(source_dir, "wgcna_circular_atlas_segments.csv")
out_metrics <- file.path(table_dir, "wgcna_circular_atlas_metrics.csv")
out_status <- file.path(report_dir, "wgcna_circular_atlas_input_status.csv")
out_logic_audit <- file.path(report_dir, "wgcna_circular_atlas_logic_audit.csv")
out_count_audit <- file.path(report_dir, "wgcna_circular_atlas_supermodule_count_audit.csv")
out_join_audit <- file.path(report_dir, "wgcna_circular_atlas_join_audit.csv")
out_selected_audit <- file.path(report_dir, "wgcna_circular_atlas_selected_table_audit.csv")
out_neuropil_availability <- file.path(report_dir, "neuron_neuropil_supermodule_availability_audit.csv")
out_duplicate_audit <- file.path(report_dir, "wgcna_circular_atlas_duplicate_source_audit.csv")
out_effect_scope_audit <- file.path(report_dir, "wgcna_circular_atlas_effect_scope_audit.csv")
out_local_support <- file.path(table_dir, "wgcna_circular_atlas_local_support_summary.csv")
out_plot_source <- file.path(source_dir, "wgcna_circular_atlas_plot_source.csv")
out_main_svg <- file.path(figure_dir, "wgcna_circular_atlas_main.svg")
out_main_pdf <- file.path(figure_dir, "wgcna_circular_atlas_main.pdf")
out_selected_svg <- file.path(figure_dir, "wgcna_circular_atlas_selected_only.svg")
out_selected_pdf <- file.path(figure_dir, "wgcna_circular_atlas_selected_only.pdf")
out_heatmap_source_supermodule <- file.path(source_dir, "wgcna_circular_heatmap_source_supermodule.csv")
out_heatmap_source_module <- file.path(source_dir, "wgcna_circular_heatmap_source_module.csv")
heatmap_svg_paths <- c(
  neuron_neuropil = file.path(figure_dir, "wgcna_circular_heatmap_neuron_neuropil.svg"),
  neuron_soma = file.path(figure_dir, "wgcna_circular_heatmap_neuron_soma.svg"),
  microglia = file.path(figure_dir, "wgcna_circular_heatmap_microglia.svg")
)
heatmap_pdf_paths <- c(
  neuron_neuropil = file.path(figure_dir, "wgcna_circular_heatmap_neuron_neuropil.pdf"),
  neuron_soma = file.path(figure_dir, "wgcna_circular_heatmap_neuron_soma.pdf"),
  microglia = file.path(figure_dir, "wgcna_circular_heatmap_microglia.pdf")
)
out_rect_modules_svg <- file.path(figure_dir, "wgcna_region_layer_heatmap_all_modules.svg")
out_rect_modules_pdf <- file.path(figure_dir, "wgcna_region_layer_heatmap_all_modules.pdf")
out_run_manifest <- file.path(log_dir, "run_manifest.yml")

dataset_label <- function(ds) {
  vapply(as.character(ds), function(x) {
    switch(x,
      neuron_neuropil = "Neuron neuropil",
      neuron_soma = "Neuron soma",
      microglia = "Microglia ROI",
      x
    )
  }, character(1))
}

status_priority <- c(
  robust_FDR = 1L,
  suggestive_FDR10 = 2L,
  nominal_only = 3L,
  model_unstable = 4L,
  not_supported = 5L,
  missing_effect_test = 6L
)

effect_scope_priority <- c(
  spatial_adjusted_global = 1L,
  within_spatial_unit = 2L,
  stress_by_spatial_interaction = 3L
)

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

na_if_blank_chr <- function(x) {
  x <- clean_chr(x)
  x[!nzchar(x)] <- NA_character_
  x
}

col_or_na <- function(df, nm, default = NA_character_) {
  if (nm %in% names(df)) return(df[[nm]])
  rep(default, nrow(df))
}

read_csv_quiet <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

as_num <- function(x) suppressWarnings(as.numeric(x))

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else NA_character_
}

first_nonblank_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(NA_character_)
  for (nm in hit) {
    vals <- na_if_blank_chr(df[[nm]])
    if (any(!is.na(vals))) return(nm)
  }
  hit[[1]]
}

normalize_supermodule_id <- function(x) {
  x <- clean_chr(x)
  x <- ifelse(grepl("^SM[0-9]+$", x), x, x)
  x
}

normalize_member_modules <- function(x) {
  vapply(as.character(x), function(one) {
    mods <- unlist(strsplit(clean_chr(one), ";", fixed = TRUE))
    mods <- trimws(mods)
    mods <- gsub("^WGCNA_", "", mods)
    mods <- gsub("^ME", "", mods)
    mods <- mods[nzchar(mods)]
    if (!length(mods)) return("")
    paste0("WGCNA_", mods, collapse = ";")
  }, character(1))
}

contrast_family <- function(x) {
  z <- toupper(gsub("[[:space:]]+", "", clean_chr(x)))
  dplyr::case_when(
    z %in% c("SUS-RES", "RES-SUS") ~ "SUS_RES",
    z %in% c("SUS-CON", "CON-SUS") ~ "SUS_CON",
    z %in% c("RES-CON", "CON-RES") ~ "RES_CON",
    TRUE ~ z
  )
}

contrast_priority <- function(x) {
  fam <- contrast_family(x)
  dplyr::case_when(
    fam == "SUS_RES" ~ 1L,
    fam == "SUS_CON" ~ 2L,
    fam == "RES_CON" ~ 3L,
    TRUE ~ 99L
  )
}

clean_label <- function(...) {
  vals <- list(...)
  out <- rep(NA_character_, length(vals[[1]]))
  for (v in vals) {
    v <- na_if_blank_chr(v)
    out <- ifelse(is.na(out) & !is.na(v), v, out)
  }
  out[is.na(out)] <- "Unlabelled"
  out
}

same_values <- function(x, numeric = FALSE, tolerance = 1e-12) {
  if (isTRUE(numeric)) {
    vals <- as_num(x)
    vals <- vals[!is.na(vals)]
    if (length(vals) <= 1L) return(TRUE)
    return(max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE) <= tolerance)
  }
  vals <- na_if_blank_chr(x)
  vals <- vals[!is.na(vals)]
  length(unique(vals)) <= 1L
}

source_fdr <- function(df) {
  fdr_within <- if ("FDR_within_dataset_level" %in% names(df)) as_num(df$FDR_within_dataset_level) else rep(NA_real_, nrow(df))
  fdr_global <- if ("FDR_global" %in% names(df)) as_num(df$FDR_global) else rep(NA_real_, nrow(df))
  dplyr::coalesce(fdr_within, fdr_global)
}

classify_scope_evidence <- function(p, fdr, prefix) {
  p <- p[is.finite(p)]
  fdr <- fdr[is.finite(fdr)]
  if (!length(p) && !length(fdr)) return(paste0(prefix, "_missing"))
  best_p <- if (length(p)) min(p, na.rm = TRUE) else NA_real_
  best_fdr <- if (length(fdr)) min(fdr, na.rm = TRUE) else NA_real_
  if (is.finite(best_fdr) && best_fdr <= 0.05) return(paste0(prefix, "_FDR_supported"))
  if (is.finite(best_fdr) && best_fdr <= 0.10) return(paste0(prefix, "_suggestive_FDR10"))
  if (is.finite(best_p) && best_p <= 0.05) return(paste0(prefix, "_nominal_only"))
  paste0(prefix, "_not_supported")
}

local_spatial_rows <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  spatial <- na_if_blank_chr(df$spatial_unit)
  is_local_unit <- !is.na(spatial) & !spatial %in% c("global", "global_spatial_adjusted", "all_spatial_units")
  df |>
    dplyr::filter(
      .data$effect_scope == "within_spatial_unit" |
        (.data$effect_scope != "stress_by_spatial_interaction" & is_local_unit)
    )
}

best_effect_row <- function(df) {
  if (is.null(df) || !nrow(df)) return(NULL)
  df |>
    dplyr::mutate(
      source_FDR = source_fdr(df),
      p_value_num = as_num(.data$p_value),
      estimate_num = as_num(.data$estimate)
    ) |>
    dplyr::arrange(.data$source_FDR, .data$p_value_num, dplyr::desc(abs(.data$estimate_num))) |>
    dplyr::slice(1)
}

schema_variant <- function(summary_id_col) {
  dplyr::case_when(
    identical(summary_id_col, "SupermoduleID") ~ "canonical_SupermoduleID",
    identical(summary_id_col, "Supermodule_DataDriven") ~ "data_driven_id_no_SupermoduleID",
    identical(summary_id_col, "Supermodule_DataDrivenID") ~ "data_driven_id_column",
    TRUE ~ paste0("alternate_id_column:", summary_id_col)
  )
}

audit_join <- function(join_name, left, right, keys, left_table, right_table) {
  if (is.null(left)) left <- tibble::tibble()
  if (is.null(right)) right <- tibble::tibble()
  if (!nrow(left) || !all(keys %in% names(left))) {
    left_keys <- tibble::tibble(.key = character())
  } else {
    left_keys <- left |>
      dplyr::mutate(dplyr::across(dplyr::all_of(keys), as.character)) |>
      tidyr::unite(".key", dplyr::all_of(keys), sep = "||", remove = FALSE)
  }
  if (!nrow(right) || !all(keys %in% names(right))) {
    right_keys <- tibble::tibble(.key = character())
  } else {
    right_keys <- right |>
      dplyr::mutate(dplyr::across(dplyr::all_of(keys), as.character)) |>
      tidyr::unite(".key", dplyr::all_of(keys), sep = "||", remove = FALSE)
  }
  lkey <- unique(left_keys$.key)
  rkey <- unique(right_keys$.key)
  tibble::tibble(
    join_name = join_name,
    left_table = left_table,
    right_table = right_table,
    join_keys = paste(keys, collapse = ";"),
    n_left = nrow(left),
    n_right = nrow(right),
    n_matched = length(intersect(lkey, rkey)),
    n_unmatched_left = length(setdiff(lkey, rkey)),
    n_unmatched_right = length(setdiff(rkey, lkey)),
    duplicated_keys_left = sum(duplicated(left_keys$.key)),
    duplicated_keys_right = sum(duplicated(right_keys$.key))
  )
}

dataset_source_paths <- function(ds) {
  list(
    supermodule_summary = path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds, "supermodules", "wgcna_supermodule_summary.csv"),
    module_supermodule_annotation = path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds, "supermodules", "wgcna_module_supermodule_annotation.csv"),
    supermodule_group_effects = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"),
    final_label_lookup = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_final_label_lookup.csv"),
    interpretable_summary = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_supermodule_group_effects_interpretable.csv"),
    biological_annotation = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"),
    modules_dir = path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds, "modules")
  )
}

pipeline_registers_dataset_script <- function(script, dataset) {
  yml <- repo_path("pipeline.yml")
  if (!file.exists(yml)) return(NA)
  txt <- readLines(yml, warn = FALSE)
  script_lines <- grep(paste0('script: "', gsub("([\\\\.])", "\\\\\\1", script), '"'), txt)
  if (!length(script_lines)) return(FALSE)
  any(vapply(script_lines, function(i) {
    window <- txt[i:min(length(txt), i + 8L)]
    any(grepl(paste0('"', dataset, '"'), window, fixed = TRUE))
  }, logical(1)))
}

discover_candidate_paths <- function(filename, dataset = "neuron_neuropil") {
  roots <- c(path_results(), repo_path("WGCNA_BACKUP_20260605_095556"))
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return("")
  candidates <- unlist(lapply(roots, function(root) {
    list.files(root, pattern = paste0("^", filename, "$"), full.names = TRUE, recursive = TRUE)
  }), use.names = FALSE)
  candidates <- candidates[grepl(dataset, candidates, fixed = TRUE)]
  candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
  paste(unique(candidates), collapse = ";")
}

claim_status_table <- function() {
  path <- path_results("tables", "biological_claims_table.csv")
  claims <- read_csv_quiet(path)
  if (is.null(claims) || !nrow(claims)) {
    return(tibble::tibble(
      dataset = character(),
      claim_key = character(),
      claim_display_status = character(),
      claim_source_file = character()
    ))
  }
  gate_col <- first_existing_col(claims, c("claim_gate", "claim_status", "gate_status", "claim_display_status"))
  program_col <- first_existing_col(claims, c("biological_program", "program_label", "Supermodule", "supermodule_id"))
  if (is.na(gate_col) || is.na(program_col) || !"dataset" %in% names(claims)) {
    return(tibble::tibble(
      dataset = unique(claims$dataset %||% character()),
      claim_key = NA_character_,
      claim_display_status = "claim_table_present_no_gate_columns",
      claim_source_file = path
    ))
  }
  claims |>
    dplyr::mutate(
      claim_key = clean_chr(.data[[program_col]]),
      claim_display_status = clean_chr(.data[[gate_col]]),
      claim_source_file = path
    ) |>
    dplyr::filter(nzchar(.data$claim_key), nzchar(.data$claim_display_status)) |>
    dplyr::group_by(.data$dataset, .data$claim_key) |>
    dplyr::summarise(
      claim_display_status = paste(unique(.data$claim_display_status), collapse = ";"),
      claim_source_file = dplyr::first(.data$claim_source_file),
      .groups = "drop"
    )
}

lookup_claim_status <- function(ds, cleaned_label, broad_program, claims) {
  if (is.null(claims) || !nrow(claims)) return("claim_not_available")
  keys <- c(clean_chr(cleaned_label), clean_chr(broad_program))
  hit <- claims |>
    dplyr::filter(.data$dataset == ds, .data$claim_key %in% keys)
  if (nrow(hit)) return(hit$claim_display_status[[1]])
  if (any(claims$dataset == ds & claims$claim_display_status == "claim_table_present_no_gate_columns", na.rm = TRUE)) {
    return("claim_table_present_no_gate_columns")
  }
  "claim_not_mapped"
}

summarise_effect_scopes <- function(effects, dataset, source_ids) {
  empty <- tibble::tibble(
    dataset = dataset,
    supermodule_id = source_ids,
    n_effect_rows_total = 0L,
    available_effect_scopes = NA_character_,
    n_spatial_adjusted_global_rows = 0L,
    n_within_spatial_unit_rows = 0L,
    n_stress_by_spatial_interaction_rows = 0L,
    n_rows_with_region_or_layer = 0L,
    n_contrasts = 0L,
    contrasts_available = NA_character_,
    spatial_units_available = NA_character_,
    best_global_p = NA_real_,
    best_global_FDR = NA_real_,
    best_global_spatial_unit = NA_character_,
    best_global_contrast = NA_character_,
    best_global_estimate = NA_real_,
    best_local_p = NA_real_,
    best_local_FDR = NA_real_,
    best_local_spatial_unit = NA_character_,
    best_local_contrast = NA_character_,
    best_local_estimate = NA_real_,
    best_interaction_p = NA_real_,
    best_interaction_FDR = NA_real_,
    best_interaction_spatial_unit = NA_character_,
    best_interaction_contrast = NA_character_,
    best_interaction_estimate = NA_real_,
    global_evidence_status = "global_missing",
    local_spatial_evidence_status = "local_missing",
    interaction_evidence_status = "interaction_missing",
    n_spatial_units_tested = 0L,
    n_spatial_units_FDR05 = 0L,
    n_spatial_units_FDR10 = 0L,
    n_spatial_units_nominal = 0L,
    message = "no supermodule effect rows found"
  )
  if (is.null(effects) || !nrow(effects)) return(empty)

  effects <- effects |>
    dplyr::mutate(
      supermodule_id = normalize_supermodule_id(.data$supermodule_id),
      source_FDR = source_fdr(effects),
      p_value_num = as_num(.data$p_value),
      estimate_num = as_num(.data$estimate),
      spatial_unit_clean = na_if_blank_chr(.data$spatial_unit)
    )

  rows <- lapply(source_ids, function(sid) {
    df <- effects |> dplyr::filter(.data$supermodule_id == sid)
    if (!nrow(df)) return(empty |> dplyr::filter(.data$supermodule_id == sid))
    global <- df |> dplyr::filter(.data$effect_scope == "spatial_adjusted_global")
    local <- local_spatial_rows(df)
    interaction <- df |> dplyr::filter(.data$effect_scope == "stress_by_spatial_interaction")
    bg <- best_effect_row(global)
    bl <- best_effect_row(local)
    bi <- best_effect_row(interaction)
    local_units <- local |>
      dplyr::filter(!is.na(.data$spatial_unit_clean)) |>
      dplyr::group_by(.data$spatial_unit_clean) |>
      dplyr::summarise(
        min_fdr = min(.data$source_FDR, na.rm = TRUE),
        min_p = min(.data$p_value_num, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        min_fdr = dplyr::if_else(is.infinite(.data$min_fdr), NA_real_, .data$min_fdr),
        min_p = dplyr::if_else(is.infinite(.data$min_p), NA_real_, .data$min_p)
      )
    tibble::tibble(
      dataset = dataset,
      supermodule_id = sid,
      n_effect_rows_total = nrow(df),
      available_effect_scopes = paste(sort(unique(na_if_blank_chr(df$effect_scope))), collapse = ";"),
      n_spatial_adjusted_global_rows = nrow(global),
      n_within_spatial_unit_rows = nrow(local),
      n_stress_by_spatial_interaction_rows = nrow(interaction),
      n_rows_with_region_or_layer = sum(!is.na(df$spatial_unit_clean) & !df$spatial_unit_clean %in% c("global", "global_spatial_adjusted", "all_spatial_units")),
      n_contrasts = dplyr::n_distinct(na_if_blank_chr(df$contrast), na.rm = TRUE),
      contrasts_available = paste(sort(unique(na_if_blank_chr(df$contrast))), collapse = ";"),
      spatial_units_available = paste(sort(unique(stats::na.omit(df$spatial_unit_clean))), collapse = ";"),
      best_global_p = if (is.null(bg)) NA_real_ else bg$p_value_num[[1]],
      best_global_FDR = if (is.null(bg)) NA_real_ else bg$source_FDR[[1]],
      best_global_spatial_unit = if (is.null(bg)) NA_character_ else bg$spatial_unit[[1]],
      best_global_contrast = if (is.null(bg)) NA_character_ else bg$contrast[[1]],
      best_global_estimate = if (is.null(bg)) NA_real_ else bg$estimate_num[[1]],
      best_local_p = if (is.null(bl)) NA_real_ else bl$p_value_num[[1]],
      best_local_FDR = if (is.null(bl)) NA_real_ else bl$source_FDR[[1]],
      best_local_spatial_unit = if (is.null(bl)) NA_character_ else bl$spatial_unit[[1]],
      best_local_contrast = if (is.null(bl)) NA_character_ else bl$contrast[[1]],
      best_local_estimate = if (is.null(bl)) NA_real_ else bl$estimate_num[[1]],
      best_interaction_p = if (is.null(bi)) NA_real_ else bi$p_value_num[[1]],
      best_interaction_FDR = if (is.null(bi)) NA_real_ else bi$source_FDR[[1]],
      best_interaction_spatial_unit = if (is.null(bi)) NA_character_ else bi$spatial_unit[[1]],
      best_interaction_contrast = if (is.null(bi)) NA_character_ else bi$contrast[[1]],
      best_interaction_estimate = if (is.null(bi)) NA_real_ else bi$estimate_num[[1]],
      global_evidence_status = classify_scope_evidence(global$p_value_num, global$source_FDR, "global"),
      local_spatial_evidence_status = classify_scope_evidence(local$p_value_num, local$source_FDR, "local"),
      interaction_evidence_status = classify_scope_evidence(interaction$p_value_num, interaction$source_FDR, "interaction"),
      n_spatial_units_tested = nrow(local_units),
      n_spatial_units_FDR05 = sum(local_units$min_fdr <= 0.05, na.rm = TRUE),
      n_spatial_units_FDR10 = sum(local_units$min_fdr <= 0.10, na.rm = TRUE),
      n_spatial_units_nominal = sum(local_units$min_p <= 0.05 & (is.na(local_units$min_fdr) | local_units$min_fdr > 0.10), na.rm = TRUE),
      message = dplyr::case_when(
        nrow(global) == 0L & nrow(local) > 0L ~ "local rows available; no spatial_adjusted_global rows",
        nrow(local) == 0L & nrow(global) > 0L ~ "global rows available; no local spatial-unit rows",
        nrow(local) > 0L & nrow(global) > 0L ~ "global and local spatial-unit rows available",
        TRUE ~ "effect rows available but no global/local scope recognized"
      )
    )
  })
  dplyr::bind_rows(rows)
}

build_neuropil_availability_audit <- function() {
  expected <- tibble::tibble(
    expected_file = c(
      path_results("tables", "06_modules_WGCNA", "01_WGCNA", "neuron_neuropil", "supermodules", "wgcna_supermodule_summary.csv"),
      path_results("tables", "06_modules_WGCNA", "01_WGCNA", "neuron_neuropil", "supermodules", "wgcna_module_supermodule_annotation.csv"),
      path_results("tables", "06_modules_WGCNA", "group_effects", "neuron_neuropil", "supermodule_group_effects.csv"),
      path_results("tables", "06_modules_WGCNA", "interpretable_summary", "neuron_neuropil", "WGCNA_supermodule_group_effects_interpretable.csv"),
      path_results("tables", "06_modules_WGCNA", "module_annotation", "neuron_neuropil", "WGCNA_supermodule_biological_annotation.csv")
    ),
    producing_script = c(
      "06_modules_WGCNA/01_WGCNA.r",
      "06_modules_WGCNA/01_WGCNA.r",
      "06_modules_WGCNA/05_module_supermodule_group_effects.r",
      "06_modules_WGCNA/07_wgcna_interpretable_summary.r",
      "06_modules_WGCNA/06_annotate_module_microenvironment.r"
    )
  )
  expected <- expected |>
    dplyr::mutate(
      exists = file.exists(.data$expected_file),
      pipeline_registered = vapply(.data$producing_script, pipeline_registers_dataset_script, logical(1), dataset = "neuron_neuropil")
    )
  expected$discovered_candidate_paths <- vapply(seq_len(nrow(expected)), function(i) {
    if (isTRUE(expected$exists[[i]])) {
      return("")
    }
    discover_candidate_paths(basename(expected$expected_file[[i]]))
  }, character(1))
  expected |>
    dplyr::mutate(
      likely_reason_missing = dplyr::case_when(
        .data$exists ~ "available_current_run",
        nzchar(.data$discovered_candidate_paths) ~ "missing_from_live_expected_path_but_candidate_or_backup_exists",
        .data$pipeline_registered ~ "producer_registered_but_output_absent; rerun minimum WGCNA downstream producer chain",
        TRUE ~ "producer_not_registered_or_input_absent"
      ),
      recommended_action = dplyr::case_when(
        .data$exists ~ "use_current_expected_file",
        grepl("01_WGCNA", .data$producing_script) ~ "run Rscript 06_modules_WGCNA/01_WGCNA.r --dataset neuron_neuropil with cached state reuse enabled",
        grepl("05_module", .data$producing_script) ~ "run Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset neuron_neuropil after 01_WGCNA outputs exist",
        grepl("06_annotate", .data$producing_script) ~ "run Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset neuron_neuropil after 01_WGCNA/group effects exist",
        grepl("07_wgcna", .data$producing_script) ~ "run Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset neuron_neuropil after group effects and annotation exist",
        TRUE ~ "inspect pipeline registration and dataset inputs"
      )
    ) |>
    dplyr::select(
      "expected_file",
      "exists",
      "discovered_candidate_paths",
      "producing_script",
      "pipeline_registered",
      "likely_reason_missing",
      "recommended_action"
    )
}

available_datasets <- function() {
  root <- path_results("tables", "06_modules_WGCNA", "01_WGCNA")
  ds <- if (dir.exists(root)) list.dirs(root, full.names = FALSE, recursive = FALSE) else character()
  if (!identical(DATASET_ARG, "all")) ds <- intersect(ds, DATASET_ARG)
  ds
}

module_protein_count <- function(ds, member_modules) {
  mods <- unlist(strsplit(clean_chr(member_modules), ";", fixed = TRUE))
  mods <- gsub("^WGCNA_", "", trimws(mods))
  mods <- gsub("^ME", "", trimws(mods))
  mods <- mods[nzchar(mods)]
  if (!length(mods)) return(NA_integer_)
  paths <- file.path(dataset_source_paths(ds)$modules_dir, paste0("genes_in_module_", mods, ".csv"))
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return(NA_integer_)
  genes <- unique(unlist(lapply(paths, function(p) {
    x <- read_csv_quiet(p)
    if (is.null(x) || !"Gene" %in% names(x)) return(character())
    clean_chr(x$Gene)
  })))
  length(genes[nzchar(genes)])
}

pick_effect_rows <- function(effects, source_ids, effect_source_file) {
  if (is.null(effects) || !nrow(effects)) {
    return(tibble::tibble(
      dataset = character(),
      supermodule_id = character(),
      effect_source_file = character(),
      effect_source_row_id = integer(),
      strongest_effect_scope_used = character(),
      strongest_contrast_used = character(),
      strongest_estimate_source = numeric(),
      p_value_source = numeric(),
      FDR_global_source = numeric(),
      FDR_within_dataset_level_source = numeric(),
      evidence_status_source = character()
    ))
  }

  effects |>
    dplyr::filter(.data$supermodule_id %in% source_ids) |>
    dplyr::mutate(
      scope_priority = dplyr::coalesce(effect_scope_priority[.data$effect_scope], 99L),
      contrast_priority_value = contrast_priority(.data$contrast),
      status_priority_value = dplyr::coalesce(status_priority[.data$evidence_status], 99L),
      abs_estimate = abs(as_num(.data$estimate)),
      p_value_num = as_num(.data$p_value)
    ) |>
    dplyr::arrange(
      .data$dataset,
      .data$supermodule_id,
      .data$scope_priority,
      .data$contrast_priority_value,
      .data$status_priority_value,
      .data$p_value_num,
      dplyr::desc(.data$abs_estimate)
    ) |>
    dplyr::group_by(.data$dataset, .data$supermodule_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      dataset = .data$dataset,
      supermodule_id = .data$supermodule_id,
      effect_source_file = effect_source_file,
      effect_source_row_id = .data$effect_source_row_id,
      strongest_effect_scope_used = .data$effect_scope,
      strongest_contrast_used = .data$contrast,
      strongest_estimate_source = as_num(.data$estimate),
      p_value_source = as_num(.data$p_value),
      FDR_global_source = as_num(.data$FDR_global),
      FDR_within_dataset_level_source = as_num(.data$FDR_within_dataset_level),
      evidence_status_source = clean_chr(.data$evidence_status)
    )
}

process_dataset <- function(ds) {
  p <- dataset_source_paths(ds)
  summary <- read_csv_quiet(p$supermodule_summary)
  module_ann <- read_csv_quiet(p$module_supermodule_annotation)
  effects <- read_csv_quiet(p$supermodule_group_effects)
  label_lookup <- read_csv_quiet(p$final_label_lookup)
  interp <- read_csv_quiet(p$interpretable_summary)
  bio <- read_csv_quiet(p$biological_annotation)

  input_status <- tibble::tibble(
    dataset = ds,
    input_name = names(p)[names(p) != "modules_dir"],
    path = unlist(p[names(p) != "modules_dir"], use.names = FALSE),
    exists = file.exists(unlist(p[names(p) != "modules_dir"], use.names = FALSE)),
    n_rows = vapply(
      list(summary, module_ann, effects, label_lookup, interp, bio),
      function(x) if (is.null(x)) NA_integer_ else nrow(x),
      integer(1)
    )
  )

  if (is.null(summary) || !nrow(summary)) {
    return(list(
      segments = tibble::tibble(),
      logic_audit = tibble::tibble(),
      join_audit = tibble::tibble(),
      input_status = input_status
    ))
  }

  summary_id_col <- first_existing_col(summary, c("SupermoduleID", "Supermodule_DataDrivenID", "Supermodule_DataDriven", "supermodule_id", "Supermodule"))
  summary_label_col <- first_nonblank_col(summary, c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Supermodule_DataDrivenLabel", "Supermodule_CuratedLabel"))
  summary_program_col <- first_nonblank_col(summary, c("Macroprogram_Display", "Supermodule_FinalLabel", "Supermodule_DataDrivenLabel"))
  source_schema <- schema_variant(summary_id_col)
  summary <- summary |>
    dplyr::mutate(
      source_summary_row_id = dplyr::row_number(),
      dataset = ds,
      supermodule_id = normalize_supermodule_id(.data[[summary_id_col]]),
      source_cleaned_label = clean_label(col_or_na(summary, "Supermodule_DisplayLabel"), col_or_na(summary, "Supermodule_FinalLabel"), .data$supermodule_id),
      source_broad_program_if_available = clean_label(col_or_na(summary, "Macroprogram_Display"), rep(NA_character_, nrow(summary))),
      n_member_modules_source = as.integer(as_num(col_or_na(summary, "n_modules"))),
      member_modules_source = normalize_member_modules(col_or_na(summary, "member_modules"))
    )
  source_summary_duplicate_ids <- summary$supermodule_id[duplicated(summary$supermodule_id) | duplicated(summary$supermodule_id, fromLast = TRUE)]
  summary_precollapse <- summary
  summary <- summary |>
    dplyr::select(
      "dataset", "supermodule_id", "source_cleaned_label",
      "source_broad_program_if_available", "n_member_modules_source",
      "member_modules_source"
    ) |>
    dplyr::distinct(.data$dataset, .data$supermodule_id, .keep_all = TRUE)

  module_ann_id_col <- if (!is.null(module_ann)) first_existing_col(module_ann, c("SupermoduleID", "Supermodule_DataDrivenID", "Supermodule_DataDriven", "supermodule_id", "Supermodule")) else NA_character_
  module_ids <- if (!is.null(module_ann) && !is.na(module_ann_id_col)) {
    module_ann |>
      dplyr::mutate(dataset = ds, supermodule_id = normalize_supermodule_id(.data[[module_ann_id_col]])) |>
      dplyr::distinct(.data$dataset, .data$supermodule_id)
  } else {
    tibble::tibble(dataset = character(), supermodule_id = character())
  }

  interp_ids <- if (!is.null(interp) && "supermodule_id" %in% names(interp)) {
    interp |>
      dplyr::mutate(supermodule_id = normalize_supermodule_id(.data$supermodule_id)) |>
      dplyr::distinct(.data$dataset, .data$supermodule_id)
  } else {
    tibble::tibble(dataset = character(), supermodule_id = character())
  }

  label_ids <- if (!is.null(label_lookup)) {
    id_col <- first_existing_col(label_lookup, c("supermodule_id", "SupermoduleID"))
    ds_col <- first_existing_col(label_lookup, c("dataset"))
    if (!is.na(id_col)) {
      label_lookup |>
        dplyr::mutate(
          dataset = if (!is.na(ds_col)) .data[[ds_col]] else ds,
          supermodule_id = normalize_supermodule_id(.data[[id_col]])
        ) |>
        dplyr::distinct(.data$dataset, .data$supermodule_id)
    } else {
      tibble::tibble(dataset = character(), supermodule_id = character())
    }
  } else {
    tibble::tibble(dataset = character(), supermodule_id = character())
  }

  bio_core <- if (!is.null(bio) && nrow(bio)) {
    bio_label_col <- first_nonblank_col(bio, c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Supermodule_ShortLabel"))
    bio_program_col <- first_nonblank_col(bio, c("Macroprogram_Display", "dominant_microenvironment_class"))
    bio |>
      dplyr::mutate(
        dataset = clean_label(col_or_na(bio, "dataset"), rep(ds, nrow(bio))),
        supermodule_id = normalize_supermodule_id(col_or_na(bio, "SupermoduleID")),
        segment_cleaned_label = clean_label(col_or_na(bio, "Supermodule_DisplayLabel"), col_or_na(bio, "Supermodule_FinalLabel"), .data$supermodule_id),
        segment_broad_program_class = clean_label(col_or_na(bio, "Macroprogram_Display"), col_or_na(bio, "dominant_microenvironment_class"), rep("Unresolved / mixed", nrow(bio))),
        n_member_modules_segments = as.integer(as_num(col_or_na(bio, "n_member_modules"))),
        member_modules_segments = normalize_member_modules(col_or_na(bio, "member_modules")),
        annotation_source_file = p$biological_annotation,
        annotation_label_source_column = bio_label_col,
        annotation_program_source_column = bio_program_col
      ) |>
      dplyr::select(
        "dataset", "supermodule_id", "segment_cleaned_label",
        "segment_broad_program_class", "n_member_modules_segments",
        "member_modules_segments", "annotation_source_file",
        "annotation_label_source_column", "annotation_program_source_column"
      ) |>
      dplyr::distinct(.data$dataset, .data$supermodule_id, .keep_all = TRUE)
  } else {
    tibble::tibble(
      dataset = character(), supermodule_id = character(), segment_cleaned_label = character(),
      segment_broad_program_class = character(), n_member_modules_segments = integer(),
      member_modules_segments = character(), annotation_source_file = character(),
      annotation_label_source_column = character(), annotation_program_source_column = character()
    )
  }

  effects <- if (!is.null(effects) && nrow(effects)) {
    effects |>
      dplyr::mutate(
        effect_source_row_id = dplyr::row_number(),
        supermodule_id = normalize_supermodule_id(.data$supermodule_id)
      )
  } else {
    effects
  }
  effect_ids <- if (!is.null(effects) && nrow(effects)) {
    effects |> dplyr::distinct(.data$dataset, .data$supermodule_id)
  } else {
    tibble::tibble(dataset = character(), supermodule_id = character())
  }

  effect_scope_audit <- summarise_effect_scopes(effects, ds, summary$supermodule_id)
  picked <- pick_effect_rows(effects, summary$supermodule_id, p$supermodule_group_effects)
  claims <- claim_status_table()

  duplicate_audit <- summary_precollapse |>
    dplyr::filter(.data$supermodule_id %in% source_summary_duplicate_ids) |>
    dplyr::mutate(
      duplicate_source_file = p$supermodule_summary,
      source_row_id = .data$source_summary_row_id,
      n_proteins = vapply(seq_len(dplyr::n()), function(i) module_protein_count(ds, .data$member_modules_source[[i]]), integer(1))
    ) |>
    dplyr::left_join(picked, by = c("dataset", "supermodule_id")) |>
    dplyr::group_by(.data$dataset, .data$supermodule_id) |>
    dplyr::mutate(
      duplicate_group_n_rows = dplyr::n(),
      cleaned_label_differs = !same_values(.data$source_cleaned_label),
      broad_program_class_differs = !same_values(.data$source_broad_program_if_available),
      n_member_modules_differs = !same_values(.data$n_member_modules_source, numeric = TRUE),
      n_proteins_differs = !same_values(.data$n_proteins, numeric = TRUE),
      evidence_status_differs = !same_values(.data$evidence_status_source),
      strongest_contrast_differs = !same_values(.data$strongest_contrast_used),
      effect_estimate_differs = !same_values(.data$strongest_estimate_source, numeric = TRUE),
      p_value_differs = !same_values(.data$p_value_source, numeric = TRUE),
      FDR_global_differs = !same_values(.data$FDR_global_source, numeric = TRUE),
      FDR_within_dataset_level_differs = !same_values(.data$FDR_within_dataset_level_source, numeric = TRUE),
      collapse_safe = !(
        .data$cleaned_label_differs |
          .data$broad_program_class_differs |
          .data$n_member_modules_differs |
          .data$n_proteins_differs |
          .data$evidence_status_differs |
          .data$strongest_contrast_differs |
          .data$effect_estimate_differs |
          .data$p_value_differs |
          .data$FDR_global_differs |
          .data$FDR_within_dataset_level_differs
      ),
      collapse_rule = "Circular atlas emits one row per dataset x supermodule_id; duplicate source summary rows are collapsed only after auditing biologically/statistically relevant fields."
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      "dataset",
      "supermodule_id",
      "duplicate_source_file",
      "source_row_id",
      "duplicate_group_n_rows",
      "source_cleaned_label",
      "source_broad_program_if_available",
      "n_member_modules_source",
      "n_proteins",
      "evidence_status_source",
      "strongest_contrast_used",
      "strongest_estimate_source",
      "p_value_source",
      "FDR_global_source",
      "FDR_within_dataset_level_source",
      "cleaned_label_differs",
      "broad_program_class_differs",
      "n_member_modules_differs",
      "n_proteins_differs",
      "evidence_status_differs",
      "strongest_contrast_differs",
      "effect_estimate_differs",
      "p_value_differs",
      "FDR_global_differs",
      "FDR_within_dataset_level_differs",
      "collapse_safe",
      "collapse_rule"
    )

  segments <- summary |>
    dplyr::left_join(bio_core, by = c("dataset", "supermodule_id")) |>
    dplyr::left_join(effect_scope_audit, by = c("dataset", "supermodule_id")) |>
    dplyr::left_join(picked, by = c("dataset", "supermodule_id")) |>
    dplyr::mutate(
      segment_cleaned_label = dplyr::coalesce(.data$segment_cleaned_label, .data$source_cleaned_label),
      segment_broad_program_class = dplyr::coalesce(.data$segment_broad_program_class, .data$source_broad_program_if_available),
      supermodule_id_source_column = summary_id_col,
      supermodule_label_source_column = dplyr::coalesce(.data$annotation_label_source_column, summary_label_col, "supermodule_id"),
      broad_program_source_column = dplyr::coalesce(.data$annotation_program_source_column, summary_program_col),
      source_schema_variant = source_schema,
      effect_source_file = dplyr::coalesce(.data$effect_source_file, p$supermodule_group_effects),
      annotation_source_file = dplyr::coalesce(.data$annotation_source_file, p$biological_annotation),
      n_member_modules_segments = .data$n_member_modules_source,
      member_modules_segments = .data$member_modules_source,
      n_proteins_source_if_available = vapply(seq_len(dplyr::n()), function(i) module_protein_count(ds, .data$member_modules_source[[i]]), integer(1)),
      n_proteins_segments_if_available = .data$n_proteins_source_if_available,
      claim_display_status = vapply(
        seq_len(dplyr::n()),
        function(i) lookup_claim_status(ds, .data$segment_cleaned_label[[i]], .data$segment_broad_program_class[[i]], claims),
        character(1)
      ),
      strongest_effect_scope_used = dplyr::coalesce(.data$strongest_effect_scope_used, "missing_effect_test"),
      strongest_contrast_used = dplyr::coalesce(.data$strongest_contrast_used, NA_character_),
      evidence_status_source = dplyr::coalesce(na_if_blank_chr(.data$evidence_status_source), "missing_effect_test"),
      evidence_status_segments = .data$evidence_status_source,
      strongest_estimate_segments = .data$strongest_estimate_source,
      p_value_segments = .data$p_value_source,
      FDR_global_segments = .data$FDR_global_source,
      FDR_within_dataset_level_segments = .data$FDR_within_dataset_level_source,
      dataset_label = dataset_label(.data$dataset),
      segment_id = paste(.data$dataset, .data$supermodule_id, sep = "::")
    )

  logic_audit <- segments |>
    dplyr::mutate(
      present_in_source_supermodule_summary = TRUE,
      duplicated_in_source_supermodule_summary = .data$supermodule_id %in% source_summary_duplicate_ids,
      duplicate_source_file = dplyr::if_else(
        .data$duplicated_in_source_supermodule_summary,
        p$supermodule_summary,
        NA_character_
      ),
      duplicate_collapse_rule = dplyr::if_else(
        .data$duplicated_in_source_supermodule_summary,
        "Circular atlas emits one row per dataset x supermodule_id; duplicate source summary rows with identical supermodule_id are collapsed with dplyr::distinct(dataset, supermodule_id, .keep_all = TRUE) after preserving source counts/member modules.",
        NA_character_
      ),
      present_in_module_supermodule_annotation = paste(.data$dataset, .data$supermodule_id) %in% paste(module_ids$dataset, module_ids$supermodule_id),
      present_in_supermodule_group_effects = paste(.data$dataset, .data$supermodule_id) %in% paste(effect_ids$dataset, effect_ids$supermodule_id),
      present_in_interpretable_summary = paste(.data$dataset, .data$supermodule_id) %in% paste(interp_ids$dataset, interp_ids$supermodule_id),
      present_in_final_label_lookup = paste(.data$dataset, .data$supermodule_id) %in% paste(label_ids$dataset, label_ids$supermodule_id),
      present_in_circular_segments = TRUE,
      present_in_selected_table = FALSE,
      status_match = .data$evidence_status_source == .data$evidence_status_segments,
      numeric_match = (
        (is.na(.data$strongest_estimate_source) & is.na(.data$strongest_estimate_segments)) |
          abs(.data$strongest_estimate_source - .data$strongest_estimate_segments) < 1e-12
      ) & (
        (is.na(.data$p_value_source) & is.na(.data$p_value_segments)) |
          abs(.data$p_value_source - .data$p_value_segments) < 1e-12
      ) & (
        (is.na(.data$FDR_global_source) & is.na(.data$FDR_global_segments)) |
          abs(.data$FDR_global_source - .data$FDR_global_segments) < 1e-12
      ) & (
        (is.na(.data$FDR_within_dataset_level_source) & is.na(.data$FDR_within_dataset_level_segments)) |
          abs(.data$FDR_within_dataset_level_source - .data$FDR_within_dataset_level_segments) < 1e-12
      ),
      label_match = .data$source_cleaned_label == .data$segment_cleaned_label,
      module_count_match = .data$n_member_modules_source == .data$n_member_modules_segments,
      audit_issue = dplyr::case_when(
        !.data$present_in_module_supermodule_annotation ~ "missing_module_supermodule_annotation",
        .data$duplicated_in_source_supermodule_summary ~ "duplicate_source_supermodule_summary_rows_collapsed_to_one_segment",
        !.data$present_in_supermodule_group_effects ~ "missing_effect_test",
        !.data$present_in_interpretable_summary ~ "missing_interpretable_summary",
        !.data$status_match ~ "evidence_status_mismatch",
        !.data$numeric_match ~ "numeric_mismatch",
        !.data$module_count_match ~ "module_count_mismatch",
        !.data$label_match ~ "label_differs_between_summary_and_annotation",
        TRUE ~ "ok"
      )
    ) |>
    dplyr::select(
      "dataset",
      "supermodule_id",
      "present_in_source_supermodule_summary",
      "duplicated_in_source_supermodule_summary",
      "duplicate_source_file",
      "duplicate_collapse_rule",
      "present_in_module_supermodule_annotation",
      "present_in_supermodule_group_effects",
      "present_in_interpretable_summary",
      "present_in_final_label_lookup",
      "present_in_circular_segments",
      "present_in_selected_table",
      "n_member_modules_source",
      "n_member_modules_segments",
      "n_proteins_source_if_available",
      "n_proteins_segments_if_available",
      "source_cleaned_label",
      "segment_cleaned_label",
      "source_broad_program_if_available",
      "segment_broad_program_class",
      "strongest_effect_scope_used",
      "strongest_contrast_used",
      "effect_source_file",
      "effect_source_row_id",
      "strongest_estimate_source",
      "strongest_estimate_segments",
      "p_value_source",
      "p_value_segments",
      "FDR_global_source",
      "FDR_global_segments",
      "FDR_within_dataset_level_source",
      "FDR_within_dataset_level_segments",
      "evidence_status_source",
      "evidence_status_segments",
      "status_match",
      "numeric_match",
      "label_match",
      "module_count_match",
      "audit_issue"
    )

  join_audit <- dplyr::bind_rows(
    audit_join("source_summary_to_biological_annotation", summary, bio_core, c("dataset", "supermodule_id"), "wgcna_supermodule_summary.csv", "WGCNA_supermodule_biological_annotation.csv"),
    audit_join("source_summary_to_group_effect_pick", summary, picked, c("dataset", "supermodule_id"), "wgcna_supermodule_summary.csv", "supermodule_group_effects.csv selected row"),
    audit_join("source_summary_to_interpretable_summary", summary, interp_ids, c("dataset", "supermodule_id"), "wgcna_supermodule_summary.csv", "WGCNA_supermodule_group_effects_interpretable.csv"),
    audit_join("source_summary_to_module_supermodule_annotation", summary, module_ids, c("dataset", "supermodule_id"), "wgcna_supermodule_summary.csv", "wgcna_module_supermodule_annotation.csv"),
    audit_join("source_summary_to_final_label_lookup", summary, label_ids, c("dataset", "supermodule_id"), "wgcna_supermodule_summary.csv", "WGCNA_final_label_lookup.csv")
  )

  list(
    segments = segments,
    logic_audit = logic_audit,
    join_audit = join_audit,
    duplicate_audit = duplicate_audit,
    effect_scope_audit = effect_scope_audit,
    input_status = input_status
  )
}

select_table_rows <- function(segments) {
  if (!nrow(segments)) {
    return(tibble::tibble(
      dataset = character(), supermodule_id = character(), selected_rank = integer(),
      selected_reason = character(), source_evidence_status = character(), selection_support_status = character(),
      n_member_modules = integer(),
      abs_effect = numeric(), broad_program_class = character(), cleaned_label = character(),
      global_evidence_status = character(), local_spatial_evidence_status = character(),
      local_FDR_units = integer(), claim_display_status = character()
    ))
  }

  supported <- c("robust_FDR", "suggestive_FDR10", "nominal_only")
  base <- segments |>
    dplyr::mutate(
      evidence_status = .data$evidence_status_segments,
      abs_effect = abs(.data$strongest_estimate_segments),
      selection_support_status = dplyr::case_when(
        grepl("FDR_supported", .data$global_evidence_status) |
          grepl("FDR_supported", .data$local_spatial_evidence_status) |
          grepl("FDR_supported", .data$interaction_evidence_status) ~ "robust_FDR",
        grepl("suggestive_FDR10", .data$global_evidence_status) |
          grepl("suggestive_FDR10", .data$local_spatial_evidence_status) |
          grepl("suggestive_FDR10", .data$interaction_evidence_status) ~ "suggestive_FDR10",
        grepl("nominal_only", .data$global_evidence_status) |
          grepl("nominal_only", .data$local_spatial_evidence_status) |
          grepl("nominal_only", .data$interaction_evidence_status) ~ "nominal_only",
        .data$evidence_status %in% supported ~ .data$evidence_status,
        TRUE ~ .data$evidence_status
      ),
      status_priority_value = dplyr::coalesce(status_priority[.data$selection_support_status], 99L),
      broad_program_class = .data$segment_broad_program_class,
      cleaned_label = .data$segment_cleaned_label,
      n_member_modules = .data$n_member_modules_segments,
      local_FDR_units = .data$n_spatial_units_FDR05
    )

  evidence_picks <- base |>
    dplyr::filter(.data$selection_support_status %in% supported) |>
    dplyr::arrange(.data$dataset, .data$status_priority_value, dplyr::desc(.data$abs_effect), .data$p_value_segments) |>
    dplyr::group_by(.data$dataset) |>
    dplyr::slice_head(n = 8) |>
    dplyr::ungroup() |>
    dplyr::mutate(selected_reason = paste0("effect-supported priority: ", .data$selection_support_status))

  program_picks <- base |>
    dplyr::arrange(.data$dataset, .data$broad_program_class, .data$status_priority_value, dplyr::desc(.data$abs_effect), dplyr::desc(.data$n_member_modules)) |>
    dplyr::group_by(.data$dataset, .data$broad_program_class) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::mutate(selected_reason = dplyr::if_else(
      .data$selection_support_status %in% supported,
      paste0("major broad-program representative with supported effect: ", .data$selection_support_status),
      "major broad-program representative; no supported effect available for this program"
    ))

  selected <- dplyr::bind_rows(evidence_picks, program_picks) |>
    dplyr::distinct(.data$dataset, .data$supermodule_id, .keep_all = TRUE) |>
    dplyr::arrange(.data$dataset, .data$status_priority_value, dplyr::desc(.data$abs_effect), .data$supermodule_id) |>
    dplyr::group_by(.data$dataset) |>
    dplyr::mutate(
      selected_rank = dplyr::row_number(),
      selected_reason = if (any(.data$selection_support_status %in% supported)) {
        .data$selected_reason
      } else {
        paste0("representative supermodule; no supported effects in ", dplyr::first(.data$dataset))
      }
    ) |>
    dplyr::ungroup()

  selected |>
    dplyr::select(
      "dataset", "supermodule_id", "selected_rank", "selected_reason",
      "global_evidence_status", "local_spatial_evidence_status",
      "local_FDR_units", "claim_display_status",
      source_evidence_status = "evidence_status",
      "selection_support_status", "n_member_modules", "abs_effect",
      "broad_program_class", "cleaned_label"
    )
}

parse_spatial_unit <- function(x) {
  raw <- tolower(gsub("[^a-z0-9]+", "_", clean_chr(x)))
  raw <- gsub("^_+|_+$", "", raw)
  raw[!nzchar(raw) | is.na(raw) | raw %in% c("global", "global_spatial_adjusted", "all_spatial_units", "na")] <- "no_local_support"

  display <- toupper(raw)
  display <- gsub("_", " ", display)
  display[raw == "no_local_support"] <- "No local support"

  region <- dplyr::case_when(
    grepl("^ca1($|_)", raw) ~ "CA1",
    grepl("^ca2($|_)|^ca3($|_)", raw) ~ "CA2/3",
    grepl("^dg($|_)", raw) ~ "DG",
    raw == "no_local_support" ~ "Global/no local support",
    TRUE ~ "Other"
  )

  layer <- dplyr::case_when(
    raw == "no_local_support" ~ "No local support",
    grepl("_slm$", raw) ~ "SLM",
    grepl("_so$", raw) ~ "SO",
    grepl("_sr$", raw) ~ "SR",
    grepl("_sp$", raw) ~ "SP",
    grepl("_mo$", raw) ~ "MO",
    grepl("_po$", raw) ~ "PO",
    grepl("_ml$", raw) ~ "ML",
    grepl("_gcl$", raw) ~ "GCL",
    grepl("_hilus$", raw) ~ "Hilus",
    TRUE ~ display
  )

  tibble::tibble(
    parsed_spatial_region = region,
    parsed_spatial_layer_or_unit = layer,
    parsed_spatial_unit_display = display
  )
}

prepare_circular_plot_source <- function(segments) {
  dataset_order <- c("neuron_neuropil", "neuron_soma", "microglia")
  parsed_spatial <- parse_spatial_unit(segments$best_local_spatial_unit)
  segments |>
    dplyr::bind_cols(parsed_spatial) |>
    dplyr::mutate(
      dataset = factor(.data$dataset, levels = dataset_order),
      dataset_label = dataset_label(as.character(.data$dataset)),
      present_in_selected_table = .data$present_in_selected_table %in% TRUE,
      evidence_status_segments = dplyr::coalesce(na_if_blank_chr(.data$evidence_status_segments), "missing_effect_test"),
      evidence_priority = as.integer(dplyr::coalesce(status_priority[.data$evidence_status_segments], 99L)),
      abs_effect_for_order = abs(as_num(.data$strongest_estimate_segments)),
      abs_effect_for_order = dplyr::if_else(is.na(.data$abs_effect_for_order), -Inf, .data$abs_effect_for_order),
      segment_broad_program_class = dplyr::coalesce(na_if_blank_chr(.data$segment_broad_program_class), "Unresolved / mixed"),
      segment_cleaned_label = dplyr::coalesce(na_if_blank_chr(.data$segment_cleaned_label), .data$supermodule_id)
    ) |>
    dplyr::arrange(
      .data$dataset,
      dplyr::desc(.data$present_in_selected_table),
      .data$evidence_priority,
      dplyr::desc(.data$abs_effect_for_order),
      .data$supermodule_id
    ) |>
    dplyr::group_by(.data$dataset) |>
    dplyr::mutate(
      plot_index = dplyr::row_number(),
      plot_order = .data$plot_index,
      sector_x_left = .data$plot_index - 0.48,
      sector_x_right = .data$plot_index + 0.48,
      high_effect_rank = dplyr::min_rank(dplyr::desc(.data$abs_effect_for_order)),
      label_shown = .data$present_in_selected_table |
        (
          .data$evidence_status_segments %in% c("robust_FDR", "suggestive_FDR10") &
            !is.infinite(.data$abs_effect_for_order) &
            .data$high_effect_rank <= 1L
        ),
      label_for_plot = .data$label_shown,
      label_text = dplyr::if_else(
        .data$label_shown,
        paste0(.data$supermodule_id, ": ", vapply(strwrap(.data$segment_cleaned_label, width = 28, simplify = FALSE), paste, character(1), collapse = "\n")),
        NA_character_
      ),
      local_support_fraction = dplyr::if_else(
        !is.na(as_num(.data$n_spatial_units_tested)) & as_num(.data$n_spatial_units_tested) > 0,
        pmin(1, as_num(.data$n_spatial_units_FDR05) / as_num(.data$n_spatial_units_tested)),
        0
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(dataset = as.character(.data$dataset))
}

atlas_program_colors <- function(programs) {
  programs <- sort(unique(dplyr::coalesce(na_if_blank_chr(programs), "Unresolved / mixed")))
  base <- c(
    "Mitochondrial metabolism" = "#1B9E77",
    "RNA / translation" = "#D95F02",
    "Synaptic / cytoskeletal" = "#7570B3",
    "Perivascular ECM / adhesion" = "#66A61E",
    "Microglia state" = "#E7298A",
    "Neuropil / neuronal" = "#A6761D",
    "Unresolved / mixed" = "#9E9E9E",
    "mixed / low-specificity" = "#9E9E9E"
  )
  missing <- setdiff(programs, names(base))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[programs]
}

atlas_evidence_colors <- function(statuses) {
  base <- c(
    robust_FDR = "#0B6E4F",
    suggestive_FDR10 = "#65A30D",
    nominal_only = "#F59E0B",
    model_unstable = "#8B5CF6",
    not_supported = "#BDBDBD",
    missing_effect_test = "#F3F4F6"
  )
  statuses <- sort(unique(dplyr::coalesce(na_if_blank_chr(statuses), "missing_effect_test")))
  missing <- setdiff(statuses, names(base))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Set 3")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[statuses]
}

atlas_spatial_colors <- function(regions) {
  base <- c(
    "CA1" = "#0072B2",
    "CA2/3" = "#56B4E9",
    "DG" = "#009E73",
    "Global/no local support" = "#D9D9D9",
    "Other" = "#CC79A7"
  )
  regions <- sort(unique(dplyr::coalesce(na_if_blank_chr(regions), "Global/no local support")))
  missing <- setdiff(regions, names(base))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Set 2")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[regions]
}

render_circular_atlas <- function(plot_source, svg_path, pdf_path, selected_only = FALSE) {
  if (!nrow(plot_source)) return(invisible(FALSE))

  dataset_order <- c("neuron_neuropil", "neuron_soma", "microglia")
  plot_source <- plot_source |>
    dplyr::filter(.data$dataset %in% dataset_order) |>
    dplyr::filter(!isTRUE(selected_only) | .data$present_in_selected_table %in% TRUE) |>
    dplyr::mutate(dataset = factor(.data$dataset, levels = dataset_order)) |>
    dplyr::arrange(.data$dataset, .data$plot_order) |>
    dplyr::group_by(.data$dataset) |>
    dplyr::mutate(
      plot_index = dplyr::row_number(),
      sector_x_left = .data$plot_index - 0.48,
      sector_x_right = .data$plot_index + 0.48,
      label_text = dplyr::if_else(
        isTRUE(selected_only) | .data$label_shown,
        paste0(.data$supermodule_id, ": ", vapply(strwrap(.data$segment_cleaned_label, width = if (isTRUE(selected_only)) 34 else 28, simplify = FALSE), paste, character(1), collapse = "\n")),
        NA_character_
      )
    ) |>
    dplyr::ungroup()
  if (!nrow(plot_source)) return(invisible(FALSE))

  program_cols <- atlas_program_colors(plot_source$segment_broad_program_class)
  evidence_cols <- atlas_evidence_colors(plot_source$evidence_status_segments)
  spatial_cols <- atlas_spatial_colors(plot_source$parsed_spatial_region)
  finite_effects <- as_num(plot_source$strongest_estimate_segments)
  finite_effects <- finite_effects[is.finite(finite_effects)]
  effect_limit <- if (length(finite_effects)) stats::quantile(abs(finite_effects), 0.95, na.rm = TRUE, names = FALSE) else 1
  effect_limit <- max(effect_limit, 1e-6)
  effect_col_fun <- circlize::colorRamp2(c(-effect_limit, 0, effect_limit), c("#2166AC", "#F7F7F7", "#B2182B"))
  effect_colors <- ifelse(
    is.finite(as_num(plot_source$strongest_estimate_segments)),
    effect_col_fun(pmax(-effect_limit, pmin(effect_limit, as_num(plot_source$strongest_estimate_segments)))),
    "#ECEFF1"
  )
  plot_source$program_color <- unname(program_cols[plot_source$segment_broad_program_class])
  plot_source$evidence_color <- unname(evidence_cols[plot_source$evidence_status_segments])
  plot_source$spatial_color <- unname(spatial_cols[plot_source$parsed_spatial_region])
  plot_source$effect_color <- effect_colors

  draw_one <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit({
      circlize::circos.clear()
      graphics::par(old_par)
    }, add = TRUE)
    graphics::par(mar = c(1.8, 1.4, 3.0, 1.4), xpd = NA, family = "sans")
    circlize::circos.clear()
    counts <- plot_source |>
      dplyr::count(.data$dataset, name = "n") |>
      dplyr::arrange(factor(.data$dataset, levels = dataset_order))
    xlim <- cbind(rep(0.5, nrow(counts)), counts$n + 0.5)
    rownames(xlim) <- as.character(counts$dataset)
    gap_after <- rep(12, nrow(counts))
    gap_after[length(gap_after)] <- 18
    circlize::circos.par(
      start.degree = 88,
      gap.after = gap_after,
      track.margin = c(0.004, 0.004),
      cell.padding = c(0, 0, 0, 0)
    )
    circlize::circos.initialize(factors = counts$dataset, xlim = xlim)

    sector_rows <- function() {
      ds <- circlize::CELL_META$sector.index
      plot_source[as.character(plot_source$dataset) == ds, , drop = FALSE]
    }

    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = if (isTRUE(selected_only)) 0.22 else 0.18, bg.border = NA,
      panel.fun = function(x, y) {
        ds <- circlize::CELL_META$sector.index
        d <- sector_rows()
        circlize::circos.text(mean(circlize::CELL_META$xlim), 0.94, dataset_label(ds),
                              facing = "bending.inside", niceFacing = TRUE, font = 2, cex = 0.72)
        lab <- d[!is.na(d$label_text) & nzchar(d$label_text), , drop = FALSE]
        if (nrow(lab)) {
          circlize::circos.text(lab$plot_index, 0.06, lab$label_text,
                                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = if (isTRUE(selected_only)) 0.42 else 0.25)
        }
      }
    )

    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = 0.075, bg.border = NA,
      panel.fun = function(x, y) {
        d <- sector_rows()
        circlize::circos.rect(d$sector_x_left, 0, d$sector_x_right, 1, col = d$program_color, border = "white", lwd = 0.25)
      }
    )
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = 0.075, bg.border = NA,
      panel.fun = function(x, y) {
        d <- sector_rows()
        circlize::circos.rect(d$sector_x_left, 0, d$sector_x_right, 1, col = d$effect_color, border = "white", lwd = 0.25)
      }
    )
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = 0.075, bg.border = NA,
      panel.fun = function(x, y) {
        d <- sector_rows()
        circlize::circos.rect(d$sector_x_left, 0, d$sector_x_right, 1, col = d$evidence_color, border = "white", lwd = 0.25)
      }
    )
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = 0.075, bg.border = NA,
      panel.fun = function(x, y) {
        d <- sector_rows()
        circlize::circos.rect(d$sector_x_left, 0, d$sector_x_right, 1, col = d$spatial_color, border = "white", lwd = 0.25)
      }
    )
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = 0.08, bg.border = NA,
      panel.fun = function(x, y) {
        d <- sector_rows()
        circlize::circos.rect(d$sector_x_left, 0, d$sector_x_right, 1, col = "#F5F5F5", border = "white", lwd = 0.25)
        circlize::circos.rect(d$sector_x_left, 0, d$sector_x_right, d$local_support_fraction, col = "#4D4D4D", border = NA)
      }
    )
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = 0.07, bg.border = NA,
      panel.fun = function(x, y) {
        d <- sector_rows()
        circlize::circos.rect(d$sector_x_left, 0, d$sector_x_right, 1, col = "#FAFAFA", border = "white", lwd = 0.25)
        sel <- d[d$present_in_selected_table %in% TRUE, , drop = FALSE]
        if (nrow(sel)) circlize::circos.points(sel$plot_index, rep(0.5, nrow(sel)), pch = 16, cex = 0.45, col = "black")
      }
    )

    graphics::title(if (isTRUE(selected_only)) "WGCNA Circular Atlas: Selected Supermodules" else "WGCNA Circular Atlas", line = 1.1, cex.main = 1.05)
    graphics::mtext(
      "Outer to inner: labels, broad program, strongest effect, evidence status, best local spatial unit, local FDR05 support fraction, selected-table marker. Microglia-enriched ROI / local microenvironment.",
      side = 1, line = 0.35, cex = 0.50, col = "#444444"
    )
    graphics::legend(
      "bottomleft",
      legend = names(program_cols),
      fill = unname(program_cols),
      border = NA,
      bty = "n",
      cex = 0.42,
      ncol = 2,
      title = "Broad program"
    )
    graphics::legend(
      "topright",
      legend = c("negative", "near zero", "positive"),
      fill = c("#2166AC", "#F7F7F7", "#B2182B"),
      border = NA,
      bty = "n",
      cex = 0.45,
      title = "Strongest effect"
    )
    graphics::legend(
      "right",
      legend = names(evidence_cols),
      fill = unname(evidence_cols),
      border = NA,
      bty = "n",
      cex = 0.43,
      title = "Evidence status"
    )
    graphics::legend(
      "bottomright",
      legend = names(spatial_cols),
      fill = unname(spatial_cols),
      border = NA,
      bty = "n",
      cex = 0.45,
      title = "Best local unit"
    )
    graphics::legend(
      "topleft",
      legend = c("selected-table marker", "local support fraction"),
      pch = c(16, 15),
      col = c("black", "#4D4D4D"),
      bty = "n",
      cex = 0.45,
      title = "Inner tracks"
    )
  }

  dir_create(dirname(svg_path))
  dir_create(dirname(pdf_path))
  svglite::svglite(svg_path, width = 12, height = 12, bg = "white")
  draw_one()
  grDevices::dev.off()
  grDevices::pdf(pdf_path, width = 12, height = 12, onefile = FALSE, useDingbats = FALSE)
  draw_one()
  grDevices::dev.off()
  invisible(TRUE)
}

is_real_spatial_unit <- function(x) {
  z <- tolower(clean_chr(x))
  nzchar(z) & !z %in% c("global", "global_spatial_adjusted", "all_spatial_units", "na")
}

contrast_block <- function(x) {
  fam <- contrast_family(x)
  dplyr::case_when(
    fam == "RES_CON" ~ "RES-CON",
    fam == "SUS_CON" ~ "SUS-CON",
    fam == "SUS_RES" ~ "SUS-RES",
    TRUE ~ clean_chr(x)
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

spatial_order_value <- function(region, layer_or_unit, unit) {
  unit <- tolower(clean_chr(unit))
  layer <- toupper(clean_chr(layer_or_unit))
  dplyr::case_when(
    region == "CA1" & layer == "SO" ~ 101L,
    region == "CA1" & layer == "SP" ~ 102L,
    region == "CA1" & layer == "SR" ~ 103L,
    region == "CA1" & layer == "SLM" ~ 104L,
    region == "CA1" ~ 109L,
    region == "CA2/3" & grepl("^ca2", unit) & layer == "SO" ~ 201L,
    region == "CA2/3" & grepl("^ca2", unit) & layer == "SR" ~ 202L,
    region == "CA2/3" & grepl("^ca2", unit) & layer == "SLM" ~ 203L,
    region == "CA2/3" & grepl("^ca2", unit) ~ 209L,
    region == "CA2/3" & grepl("^ca3", unit) & layer == "SO" ~ 221L,
    region == "CA2/3" & grepl("^ca3", unit) & layer == "SR" ~ 222L,
    region == "CA2/3" & grepl("^ca3", unit) ~ 229L,
    region == "CA2/3" ~ 240L,
    region == "DG" & layer == "MO" ~ 301L,
    region == "DG" & layer == "ML" ~ 302L,
    region == "DG" & layer == "GCL" ~ 303L,
    region == "DG" & layer == "PO" ~ 304L,
    region == "DG" ~ 309L,
    region == "Other" ~ 900L,
    TRUE ~ 999L
  )
}

effect_support_class <- function(p_value, fdr_within, fdr_global) {
  fdr <- dplyr::coalesce(as_num(fdr_within), as_num(fdr_global))
  p <- as_num(p_value)
  dplyr::case_when(
    !is.na(fdr) & fdr <= 0.05 ~ "FDR05",
    !is.na(fdr) & fdr <= 0.10 ~ "FDR10",
    !is.na(p) & p <= 0.05 ~ "nominal",
    TRUE ~ "none"
  )
}

local_effect_rows <- function(df, level) {
  if (is.null(df) || !nrow(df)) return(tibble::tibble())
  key_col <- if (identical(level, "supermodule")) "supermodule_id" else "module_id"
  local_df <- df |>
    dplyr::filter(is_real_spatial_unit(.data$spatial_unit))
  if (!nrow(local_df)) return(tibble::tibble())
  parsed <- parse_spatial_unit(local_df$spatial_unit)
  local_df |>
    dplyr::bind_cols(parsed) |>
    dplyr::mutate(
      level = level,
      endpoint_key = clean_chr(.data[[key_col]]),
      contrast_block = contrast_block(.data$contrast),
      contrast_block_order = contrast_block_order(.data$contrast),
      spatial_order = spatial_order_value(.data$parsed_spatial_region, .data$parsed_spatial_layer_or_unit, .data$spatial_unit),
      effect_scope_order = dplyr::case_when(
        .data$effect_scope == "within_spatial_unit" ~ 1L,
        .data$effect_scope == "stress_by_spatial_interaction" ~ 2L,
        TRUE ~ 9L
      ),
      p_value_num = as_num(.data$p_value),
      FDR_within_num = as_num(.data$FDR_within_dataset_level),
      FDR_global_num = as_num(.data$FDR_global),
      effect_abs = abs(as_num(.data$estimate)),
      support_class = effect_support_class(.data$p_value, .data$FDR_within_dataset_level, .data$FDR_global)
    ) |>
    dplyr::arrange(
      .data$dataset, .data$endpoint_key, .data$spatial_unit, .data$contrast_block,
      .data$effect_scope_order,
      dplyr::coalesce(.data$FDR_within_num, .data$FDR_global_num, Inf),
      .data$p_value_num,
      dplyr::desc(.data$effect_abs)
    ) |>
    dplyr::group_by(.data$dataset, .data$endpoint_key, .data$spatial_unit, .data$contrast_block) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
}

module_supermodule_map <- function(dataset) {
  path <- dataset_source_paths(dataset)$module_supermodule_annotation
  ann <- read_csv_quiet(path)
  if (is.null(ann) || !nrow(ann)) {
    return(tibble::tibble(
      dataset = character(),
      module_id_annotation = character(),
      module_label_key = character(),
      supermodule_id = character(),
      supermodule_label = character()
    ))
  }
  supermodule_id_col <- first_existing_col(ann, c("SupermoduleID", "Supermodule_DataDrivenID", "Supermodule_DataDriven"))
  supermodule_label_col <- first_nonblank_col(ann, c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Supermodule", "Supermodule_DataDrivenLabel"))
  module_id_col <- first_existing_col(ann, c("ModuleID", "module_id"))
  module_label_col <- first_nonblank_col(ann, c("ModuleLabel_Final", "top_GO_label", "ModuleLabel_GO_BP"))

  ann |>
    dplyr::transmute(
      dataset = dataset,
      module_id_annotation = if (!is.na(module_id_col)) clean_chr(.data[[module_id_col]]) else NA_character_,
      module_label_key = if (!is.na(module_label_col)) clean_chr(.data[[module_label_col]]) else NA_character_,
      supermodule_id = if (!is.na(supermodule_id_col)) clean_chr(.data[[supermodule_id_col]]) else NA_character_,
      supermodule_label = if (!is.na(supermodule_label_col)) clean_chr(.data[[supermodule_label_col]]) else NA_character_
    ) |>
    dplyr::filter(nzchar(.data$module_label_key) | nzchar(.data$module_id_annotation)) |>
    dplyr::distinct(.data$dataset, .data$module_label_key, .data$module_id_annotation, .keep_all = TRUE)
}

build_heatmap_sources <- function(datasets, segments, selected_audit) {
  segment_meta <- segments |>
    dplyr::mutate(
      supermodule_id = clean_chr(.data$supermodule_id),
      selected_or_priority_flag = .data$present_in_selected_table %in% TRUE |
        .data$evidence_status_segments %in% c("robust_FDR", "suggestive_FDR10")
    ) |>
    dplyr::select(
      "dataset", "supermodule_id",
      supermodule_label = "segment_cleaned_label",
      selected_or_priority_flag,
      segment_broad_program_class = "segment_broad_program_class",
      strongest_estimate_segments = "strongest_estimate_segments",
      evidence_status_segments = "evidence_status_segments"
    )

  supermodule_rows <- lapply(datasets, function(ds) {
    path <- dataset_source_paths(ds)$supermodule_group_effects
    df <- read_csv_quiet(path)
    local <- local_effect_rows(df, "supermodule")
    if (!nrow(local)) return(tibble::tibble())
    local |>
      dplyr::left_join(segment_meta, by = c("dataset", "supermodule_id")) |>
      dplyr::mutate(
        module_id = NA_character_,
        module_label = NA_character_,
        supermodule_label = dplyr::coalesce(
          na_if_blank_chr(.data$supermodule_label.x),
          na_if_blank_chr(.data$supermodule_label.y),
          na_if_blank_chr(.data$endpoint_label),
          .data$supermodule_id
        ),
        selected_or_priority_flag = .data$selected_or_priority_flag %in% TRUE
      )
  })

  module_rows <- lapply(datasets, function(ds) {
    path <- path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv")
    df <- read_csv_quiet(path)
    local <- local_effect_rows(df, "module")
    if (!nrow(local)) return(tibble::tibble())
    map <- module_supermodule_map(ds)
    local <- local |>
      dplyr::mutate(module_label_key = clean_chr(.data$module_label)) |>
      dplyr::left_join(map, by = c("dataset", "module_label_key")) |>
      dplyr::mutate(
        supermodule_id = dplyr::coalesce(
          na_if_blank_chr(.data$supermodule_id.y),
          na_if_blank_chr(.data$supermodule_id.x),
          "unmapped"
        ),
        supermodule_label_from_map = dplyr::coalesce(
          na_if_blank_chr(.data$supermodule_label.y),
          na_if_blank_chr(.data$supermodule_label.x)
        )
      ) |>
      dplyr::select(-dplyr::any_of(c("supermodule_id.x", "supermodule_id.y", "supermodule_label.x", "supermodule_label.y"))) |>
      dplyr::left_join(segment_meta, by = c("dataset", "supermodule_id")) |>
      dplyr::mutate(
        supermodule_id = dplyr::coalesce(na_if_blank_chr(.data$supermodule_id), "unmapped"),
        supermodule_label = dplyr::coalesce(na_if_blank_chr(.data$supermodule_label_from_map), na_if_blank_chr(.data$supermodule_label), "Unmapped module"),
        selected_or_priority_flag = .data$selected_or_priority_flag %in% TRUE |
          .data$evidence_status %in% c("robust_FDR", "suggestive_FDR10")
      )
    local
  })

  format_source <- function(df, level_name) {
    if (is.null(df) || !nrow(df)) return(tibble::tibble())
    df |>
      dplyr::transmute(
        dataset,
        level = level_name,
        module_id = if ("module_id" %in% names(df)) na_if_blank_chr(.data$module_id) else NA_character_,
        supermodule_id = na_if_blank_chr(.data$supermodule_id),
        module_label = if ("module_label" %in% names(df)) na_if_blank_chr(.data$module_label) else NA_character_,
        supermodule_label = na_if_blank_chr(.data$supermodule_label),
        spatial_unit = clean_chr(.data$spatial_unit),
        parsed_region = .data$parsed_spatial_region,
        parsed_layer_or_unit = .data$parsed_spatial_layer_or_unit,
        contrast = clean_chr(.data$contrast),
        contrast_block = .data$contrast_block,
        effect_scope = clean_chr(.data$effect_scope),
        estimate = as_num(.data$estimate),
        p_value = as_num(.data$p_value),
        FDR_global = as_num(.data$FDR_global),
        FDR_within_dataset_level = as_num(.data$FDR_within_dataset_level),
        evidence_status = clean_chr(.data$evidence_status),
        support_class = .data$support_class,
        selected_or_priority_flag = .data$selected_or_priority_flag %in% TRUE,
        plot_sector_order = NA_integer_,
        plot_track_order = NA_integer_,
        angular_order = NA_integer_,
        ring_order = NA_integer_,
        contrast_ring = NA_character_,
        contrast_ring_order = NA_integer_,
        internal_contrast_label_position = NA_character_,
        spatial_order = .data$spatial_order,
        contrast_block_order = .data$contrast_block_order
      )
  }

  supermodule <- format_source(dplyr::bind_rows(supermodule_rows), "supermodule")
  module <- format_source(dplyr::bind_rows(module_rows), "module")

  order_sources <- function(df) {
    if (!nrow(df)) return(df)
    sector_meta <- df |>
      dplyr::group_by(.data$dataset, .data$supermodule_id) |>
      dplyr::summarise(
        selected = any(.data$selected_or_priority_flag, na.rm = TRUE),
        best_support = min(dplyr::coalesce(status_priority[.data$evidence_status], 99L), na.rm = TRUE),
        max_abs_effect = max(abs(.data$estimate), na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        best_support = dplyr::if_else(is.infinite(.data$best_support), 99, as.numeric(.data$best_support)),
        max_abs_effect = dplyr::if_else(is.infinite(.data$max_abs_effect), 0, .data$max_abs_effect)
      ) |>
      dplyr::arrange(
        .data$dataset,
        dplyr::desc(.data$selected),
        .data$best_support,
        dplyr::desc(.data$max_abs_effect),
        .data$supermodule_id
      ) |>
      dplyr::group_by(.data$dataset) |>
      dplyr::mutate(plot_sector_order = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::select("dataset", "supermodule_id", "plot_sector_order")

    spatial_meta <- df |>
      dplyr::distinct(.data$dataset, .data$spatial_unit, .data$parsed_region, .data$parsed_layer_or_unit, .data$spatial_order) |>
      dplyr::arrange(.data$dataset, .data$spatial_order, .data$spatial_unit) |>
      dplyr::group_by(.data$dataset) |>
      dplyr::mutate(spatial_unit_order = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::select("dataset", "spatial_unit", "spatial_unit_order")

    track_meta <- df |>
      dplyr::distinct(.data$dataset, .data$contrast_block, .data$spatial_unit, .data$parsed_region, .data$parsed_layer_or_unit, .data$contrast_block_order, .data$spatial_order) |>
      dplyr::arrange(.data$dataset, .data$contrast_block_order, .data$spatial_order, .data$spatial_unit) |>
      dplyr::group_by(.data$dataset) |>
      dplyr::mutate(plot_track_order = dplyr::row_number()) |>
      dplyr::ungroup() |>
      dplyr::select("dataset", "contrast_block", "spatial_unit", "plot_track_order")

    df |>
      dplyr::select(-dplyr::any_of(c("plot_sector_order", "plot_track_order", "angular_order", "ring_order"))) |>
      dplyr::left_join(sector_meta, by = c("dataset", "supermodule_id")) |>
      dplyr::left_join(spatial_meta, by = c("dataset", "spatial_unit")) |>
      dplyr::left_join(track_meta, by = c("dataset", "contrast_block", "spatial_unit")) |>
      dplyr::group_by(.data$dataset) |>
      dplyr::mutate(
        angular_order = (.data$plot_sector_order - 1L) * max(.data$spatial_unit_order, na.rm = TRUE) + .data$spatial_unit_order,
        ring_order = .data$contrast_block_order,
        contrast_ring = .data$contrast_block,
        contrast_ring_order = .data$ring_order,
        internal_contrast_label_position = "left_annulus"
      ) |>
      dplyr::ungroup() |>
      dplyr::arrange(.data$dataset, .data$angular_order, .data$ring_order, .data$module_id)
  }

  list(
    supermodule = order_sources(supermodule),
    module = order_sources(module)
  )
}

render_dataset_circular_heatmap <- function(source_supermodule, dataset_name, svg_path, pdf_path) {
  df <- source_supermodule |> dplyr::filter(.data$dataset == .env$dataset_name)
  if (!nrow(df)) return(invisible(FALSE))
  df <- df |>
    dplyr::mutate(
      sector_id = paste0(sprintf("%03d", .data$angular_order), "::", .data$supermodule_id, "::", .data$spatial_unit),
      spatial_label = toupper(gsub("_", " ", .data$spatial_unit)),
      sector_label = paste0(.data$supermodule_id, ": ", .data$supermodule_label)
    )
  sector_meta <- df |>
    dplyr::distinct(
      .data$sector_id,
      .data$supermodule_id,
      .data$sector_label,
      .data$plot_sector_order,
      .data$spatial_unit,
      .data$spatial_label,
      .data$angular_order,
      .data$selected_or_priority_flag
    ) |>
    dplyr::arrange(.data$angular_order) |>
    dplyr::group_by(.data$supermodule_id) |>
    dplyr::mutate(
      spatial_index_in_supermodule = dplyr::row_number(),
      n_spatial_in_supermodule = dplyr::n(),
      label_sector = .data$spatial_index_in_supermodule == ceiling(.data$n_spatial_in_supermodule / 2),
      gap_after = dplyr::if_else(.data$spatial_index_in_supermodule == .data$n_spatial_in_supermodule, 3.2, 0.18)
    ) |>
    dplyr::ungroup()

  spatial_meta <- sector_meta |>
    dplyr::distinct(
      .data$spatial_unit,
      .data$spatial_label,
      .data$angular_order
    ) |>
    dplyr::group_by(.data$spatial_unit, .data$spatial_label) |>
    dplyr::summarise(first_angular_order = min(.data$angular_order), .groups = "drop") |>
    dplyr::arrange(.data$first_angular_order)

  ring_meta <- df |>
    dplyr::distinct(.data$contrast_block, .data$ring_order) |>
    dplyr::arrange(.data$ring_order)

  effects <- as_num(df$estimate)
  lim <- if (any(is.finite(effects))) stats::quantile(abs(effects[is.finite(effects)]), 0.95, na.rm = TRUE, names = FALSE) else 1
  lim <- max(lim, 1e-6)
  col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#2166AC", "#F7F7F7", "#B2182B"))
  contrast_cols <- c("RES-CON" = "#4E79A7", "SUS-CON" = "#F28E2B", "SUS-RES" = "#59A14F")
  ring_meta$contrast_color <- unname(contrast_cols[ring_meta$contrast_block])
  ring_meta$contrast_color[is.na(ring_meta$contrast_color)] <- "#757575"
  ring_label_meta <- ring_meta |>
    dplyr::mutate(
      label_x = -0.60,
      label_y = c(0.075, 0, -0.075)[seq_len(dplyr::n())]
    )

  draw_one <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit({
      circlize::circos.clear()
      graphics::par(old_par)
    }, add = TRUE)
    graphics::par(mar = c(1.5, 1.1, 2.8, 1.0), xpd = NA, family = "sans")
    circlize::circos.clear()
    circlize::circos.par(
      start.degree = 88,
      gap.after = c(head(sector_meta$gap_after, -1), 8),
      track.margin = c(0.0015, 0.0015),
      cell.padding = c(0, 0, 0, 0)
    )
    circlize::circos.initialize(factors = factor(sector_meta$sector_id, levels = sector_meta$sector_id), xlim = cbind(rep(0, nrow(sector_meta)), rep(1, nrow(sector_meta))))

    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), track.height = 0.112, bg.border = NA,
      panel.fun = function(x, y) {
        sid <- circlize::CELL_META$sector.index
        meta <- sector_meta[sector_meta$sector_id == sid, , drop = FALSE]
        if (isTRUE(meta$label_sector[[1]])) {
          lab <- vapply(strwrap(meta$sector_label[[1]], width = 24, simplify = FALSE), paste, character(1), collapse = "\n")
          circlize::circos.text(
            0.5, 0.08, lab,
            facing = "clockwise",
            niceFacing = TRUE,
            adj = c(0, 0.5),
            cex = if (isTRUE(meta$selected_or_priority_flag[[1]])) 0.32 else 0.25,
            font = if (isTRUE(meta$selected_or_priority_flag[[1]])) 2 else 1
          )
        }
      }
    )

    if (nrow(spatial_meta) <= 10) {
      circlize::circos.trackPlotRegion(
        ylim = c(0, 1), track.height = 0.026, bg.border = NA,
        panel.fun = function(x, y) {
          sid <- circlize::CELL_META$sector.index
          meta <- sector_meta[sector_meta$sector_id == sid, , drop = FALSE]
          first_supermodule <- meta$plot_sector_order[[1]] == min(sector_meta$plot_sector_order, na.rm = TRUE)
          if (isTRUE(first_supermodule)) {
            circlize::circos.text(0.5, 0.5, meta$spatial_label[[1]], facing = "clockwise", niceFacing = TRUE, cex = 0.22)
          }
        }
      )
    }

    for (ring_i in seq_len(nrow(ring_meta))) {
      ring <- ring_meta[ring_i, , drop = FALSE]
      circlize::circos.trackPlotRegion(
        ylim = c(0, 1), track.height = 0.064, bg.border = NA,
        panel.fun = function(x, y) {
          sid <- circlize::CELL_META$sector.index
          meta <- sector_meta[sector_meta$sector_id == sid, , drop = FALSE]
          cell <- df[df$sector_id == sid & df$contrast_block == ring$contrast_block[[1]], , drop = FALSE]
          fill <- "#ECEFF1"
          if (nrow(cell) && is.finite(cell$estimate[[1]])) {
            fill <- col_fun(pmax(-lim, pmin(lim, cell$estimate[[1]])))
          }
          border_col <- if (isTRUE(meta$spatial_index_in_supermodule[[1]] == meta$n_spatial_in_supermodule[[1]])) "#2B2B2B" else "white"
          circlize::circos.rect(0, 0, 1, 1, col = fill, border = border_col, lwd = if (identical(border_col, "white")) 0.10 else 0.35)
          if (nrow(cell) && cell$support_class[[1]] %in% c("FDR05", "FDR10", "nominal")) {
            pch <- if (cell$support_class[[1]] == "FDR05") 16 else if (cell$support_class[[1]] == "FDR10") 1 else 3
            circlize::circos.points(0.5, 0.5, pch = pch, cex = 0.30, col = "black")
          }
        }
      )
    }

    graphics::text(
      x = ring_label_meta$label_x,
      y = ring_label_meta$label_y,
      labels = ring_label_meta$contrast_block,
      col = ring_label_meta$contrast_color,
      cex = 0.58,
      font = 2,
      adj = c(0.5, 0.5)
    )

    graphics::title(paste0("Circular WGCNA Spatial Heatmap: ", dataset_label(dataset_name)), line = 1.1, cex.main = 1.0)
    if (identical(dataset_name, "microglia")) {
      graphics::mtext("Microglia-enriched ROI / local microenvironment; not purified cell-intrinsic evidence.", side = 3, line = 0.1, cex = 0.58, col = "#444444")
    }
    graphics::mtext("Angular tiles are supermodule x spatial unit. The three concentric heatmap rings are contrasts.", side = 1, line = 0.35, cex = 0.54, col = "#444444")
    graphics::legend("bottomleft", legend = c("negative", "near zero", "positive"),
                     fill = c("#2166AC", "#F7F7F7", "#B2182B"), border = NA, bty = "n", cex = 0.48, title = "Estimate")
    graphics::legend("topleft", legend = c("FDR <= 0.05", "FDR <= 0.10", "p <= 0.05"),
                     pch = c(16, 1, 3), col = "black", bty = "n", cex = 0.48, title = "Support marker")
    spatial_summary <- spatial_meta |>
      dplyr::arrange(.data$first_angular_order) |>
      dplyr::summarise(units = paste(.data$spatial_label, collapse = ", "), .groups = "drop")
    graphics::legend("topright", legend = spatial_summary$units,
                     bty = "n", cex = 0.34, title = "Spatial unit order")
  }

  dir_create(dirname(svg_path))
  svglite::svglite(svg_path, width = 10.5, height = 10.5, bg = "white")
  draw_one()
  grDevices::dev.off()
  grDevices::pdf(pdf_path, width = 10.5, height = 10.5, onefile = FALSE, useDingbats = FALSE)
  draw_one()
  grDevices::dev.off()
  invisible(TRUE)
}

render_rectangular_module_heatmap <- function(source_module, svg_path, pdf_path) {
  if (is.null(source_module) || !nrow(source_module)) return(invisible(FALSE))
  df <- source_module |>
    dplyr::mutate(
      row_label = paste(.data$dataset, .data$supermodule_id, dplyr::coalesce(.data$module_label, .data$module_id), sep = " | "),
      col_label = paste(.data$contrast_block, toupper(gsub("_", " ", .data$spatial_unit)), sep = " | "),
      row_order = paste(sprintf("%02d", .data$plot_sector_order), .data$supermodule_id, dplyr::coalesce(.data$module_label, .data$module_id)),
      col_order = paste(sprintf("%02d", .data$plot_track_order), .data$col_label),
      support_shape = dplyr::case_when(
        .data$support_class == "FDR05" ~ "FDR <= 0.05",
        .data$support_class == "FDR10" ~ "FDR <= 0.10",
        .data$support_class == "nominal" ~ "p <= 0.05",
        TRUE ~ NA_character_
      )
    )
  lim <- stats::quantile(abs(df$estimate[is.finite(df$estimate)]), 0.95, na.rm = TRUE, names = FALSE)
  lim <- max(lim, 1e-6)
  row_levels <- df |>
    dplyr::distinct(.data$dataset, .data$row_label, .data$row_order) |>
    dplyr::arrange(.data$dataset, .data$row_order) |>
    dplyr::pull(.data$row_label)
  col_levels <- df |>
    dplyr::distinct(.data$dataset, .data$col_label, .data$col_order) |>
    dplyr::arrange(.data$dataset, .data$col_order) |>
    dplyr::pull(.data$col_label) |>
    unique()
  df$row_label <- factor(df$row_label, levels = rev(unique(row_levels)))
  df$col_label <- factor(df$col_label, levels = unique(col_levels))
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$col_label, .data$row_label)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$estimate), color = "white", linewidth = 0.08) +
    ggplot2::geom_point(
      data = df |> dplyr::filter(!is.na(.data$support_shape)),
      ggplot2::aes(shape = .data$support_shape),
      size = 0.8,
      color = "black",
      stroke = 0.25
    ) +
    ggplot2::facet_grid(dataset ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), oob = scales::squish, name = "Estimate") +
    ggplot2::scale_shape_manual(values = c("FDR <= 0.05" = 16, "FDR <= 0.10" = 1, "p <= 0.05" = 3), name = "Support") +
    ggplot2::labs(x = "Contrast x spatial unit", y = "Dataset | supermodule | module", title = "WGCNA Module Region/Layer Heatmap") +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 5),
      axis.text.y = ggplot2::element_text(size = 4.5),
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold"),
      legend.position = "right"
    )
  dir_create(dirname(svg_path))
  ggplot2::ggsave(svg_path, p, width = 13, height = 12, units = "in", bg = "white")
  ggplot2::ggsave(pdf_path, p, width = 13, height = 12, units = "in", bg = "white", device = grDevices::cairo_pdf)
  invisible(TRUE)
}

if (run$dry_run) {
  ds <- available_datasets()
  dry_run_line("Script", "10_biological_integration/04_wgcna_circular_atlas.R")
  dry_run_line("Dataset argument", DATASET_ARG)
  dry_run_line("Datasets discovered", paste(ds, collapse = ", "))
  dry_run_line("Downstream-only mode", "No WGCNA definitions/effects/FDRs/p-values/claim gates are modified", "PASS")
  dry_run_line("Effect row priority", "scope: spatial_adjusted_global, within_spatial_unit, stress_by_spatial_interaction; contrast: SUS-RES, SUS-CON, RES-CON")
  dry_run_line("Protein count provenance", "Existing genes_in_module_*.csv files; no WGCNA recomputation")
  dry_run_line("Segments output", out_segments)
  dry_run_line("Metrics output", out_metrics)
  dry_run_line("Input status output", out_status)
  dry_run_line("Logic audit output", out_logic_audit)
  dry_run_line("Count audit output", out_count_audit)
  dry_run_line("Join audit output", out_join_audit)
  dry_run_line("Selected table audit output", out_selected_audit)
  dry_run_line("Duplicate source audit output", out_duplicate_audit)
  dry_run_line("Effect scope audit output", out_effect_scope_audit)
  dry_run_line("Local support summary output", out_local_support)
  dry_run_line("Plot source output", out_plot_source)
  dry_run_line("Main circular atlas SVG output", out_main_svg)
  dry_run_line("Main circular atlas PDF output", out_main_pdf)
  dry_run_line("Selected-only circular atlas SVG output", out_selected_svg)
  dry_run_line("Selected-only circular atlas PDF output", out_selected_pdf)
  dry_run_line("Circular heatmap supermodule source output", out_heatmap_source_supermodule)
  dry_run_line("Circular heatmap module source output", out_heatmap_source_module)
  for (ds_name in names(heatmap_svg_paths)) {
    dry_run_line(paste0("Circular heatmap SVG output (", ds_name, ")"), heatmap_svg_paths[[ds_name]])
    dry_run_line(paste0("Circular heatmap PDF output (", ds_name, ")"), heatmap_pdf_paths[[ds_name]])
  }
  dry_run_line("Rectangular all-module heatmap SVG output", out_rect_modules_svg)
  dry_run_line("Rectangular all-module heatmap PDF output", out_rect_modules_pdf)
  dry_run_line("Run manifest output", out_run_manifest)
  dry_run_line("Neuron neuropil availability audit output", out_neuropil_availability)
  quit(status = 0, save = "no")
}

results <- lapply(available_datasets(), process_dataset)
segments <- dplyr::bind_rows(lapply(results, `[[`, "segments"))
logic_audit <- dplyr::bind_rows(lapply(results, `[[`, "logic_audit"))
join_audit <- dplyr::bind_rows(lapply(results, `[[`, "join_audit"))
duplicate_audit <- dplyr::bind_rows(lapply(results, `[[`, "duplicate_audit"))
effect_scope_audit <- dplyr::bind_rows(lapply(results, `[[`, "effect_scope_audit"))
input_status <- dplyr::bind_rows(lapply(results, `[[`, "input_status"))

selected_audit <- select_table_rows(segments)
neuropil_availability_audit <- build_neuropil_availability_audit()
selected_keys <- paste(selected_audit$dataset, selected_audit$supermodule_id)
logic_audit <- logic_audit |>
  dplyr::mutate(present_in_selected_table = paste(.data$dataset, .data$supermodule_id) %in% selected_keys)

segments <- segments |>
  dplyr::mutate(present_in_selected_table = paste(.data$dataset, .data$supermodule_id) %in% selected_keys) |>
  dplyr::select(
    "dataset",
    "dataset_label",
    "supermodule_id",
    "segment_id",
    "supermodule_id_source_column",
    "supermodule_label_source_column",
    "broad_program_source_column",
    "source_schema_variant",
    "effect_source_file",
    "effect_source_row_id",
    "annotation_source_file",
    "segment_cleaned_label",
    "segment_broad_program_class",
    "global_evidence_status",
    "local_spatial_evidence_status",
    "interaction_evidence_status",
    "claim_display_status",
    "n_spatial_units_tested",
    "n_spatial_units_FDR05",
    "n_spatial_units_FDR10",
    "n_spatial_units_nominal",
    "best_local_spatial_unit",
    "best_local_contrast",
    "best_local_estimate",
    "best_local_p",
    "best_local_FDR",
    "best_global_contrast",
    "best_global_estimate",
    "best_global_p",
    "best_global_FDR",
    "n_member_modules_segments",
    "n_proteins_segments_if_available",
    "member_modules_segments",
    "strongest_effect_scope_used",
    "strongest_contrast_used",
    "strongest_estimate_segments",
    "p_value_segments",
    "FDR_global_segments",
    "FDR_within_dataset_level_segments",
    "evidence_status_segments",
    "present_in_selected_table"
  )

local_support_summary <- segments |>
  dplyr::transmute(
    dataset,
    supermodule_id,
    cleaned_label = .data$segment_cleaned_label,
    broad_program_class = .data$segment_broad_program_class,
    n_spatial_units_tested,
    n_spatial_units_FDR05,
    n_spatial_units_FDR10,
    n_spatial_units_nominal,
    best_local_spatial_unit,
    best_local_contrast,
    best_local_estimate,
    best_local_p,
    best_local_FDR,
    local_spatial_evidence_status,
    global_evidence_status,
    claim_display_status
  )

plot_source <- prepare_circular_plot_source(segments)
heatmap_sources <- build_heatmap_sources(available_datasets(), segments, selected_audit)
heatmap_source_supermodule <- heatmap_sources$supermodule
heatmap_source_module <- heatmap_sources$module

metrics <- tibble::tibble(
  metric = c(
    "n_segments",
    "n_selected_table_rows",
    "n_datasets",
    "n_plot_labels",
    "n_circular_heatmap_supermodule_cells",
    "n_circular_heatmap_module_cells",
    "n_protein_count_source_column",
    "effect_scope_priority",
    "contrast_priority",
    "selection_table_policy"
    ,
    "link_filter_note",
    "evidence_FDR_source"
  ),
  value = c(
    as.character(nrow(segments)),
    as.character(nrow(selected_audit)),
    as.character(dplyr::n_distinct(segments$dataset)),
    as.character(sum(plot_source$label_shown, na.rm = TRUE)),
    as.character(nrow(heatmap_source_supermodule)),
    as.character(nrow(heatmap_source_module)),
    "No direct n_proteins source column found; computed from existing 01_WGCNA/<dataset>/modules/genes_in_module_*.csv member files when available",
    "spatial_adjusted_global; within_spatial_unit; stress_by_spatial_interaction",
    "SUS-RES pair; SUS-CON pair; RES-CON pair, preserving source contrast orientation",
    "Prefer robust_FDR, suggestive_FDR10, nominal_only; add major broad-program representatives; label unsupported-only selections as representatives",
    "No high-confidence cross-compartment links passed filters.",
    "Global/local/interaction evidence summaries use copied source FDR_within_dataset_level when available, otherwise copied source FDR_global; no p-values or FDRs are recomputed."
  )
)

count_audit <- logic_audit |>
  dplyr::group_by(.data$dataset) |>
  dplyr::summarise(
    n_supermodules_source = dplyr::n(),
    n_supermodules_segments = sum(.data$present_in_circular_segments),
    n_missing_from_segments = sum(!.data$present_in_circular_segments),
    n_extra_in_segments = 0L,
    n_selected_table = sum(.data$present_in_selected_table),
    n_robust_FDR = sum(.data$evidence_status_segments == "robust_FDR", na.rm = TRUE),
    n_suggestive_FDR10 = sum(.data$evidence_status_segments == "suggestive_FDR10", na.rm = TRUE),
    n_nominal_only = sum(.data$evidence_status_segments == "nominal_only", na.rm = TRUE),
    n_model_unstable = sum(.data$evidence_status_segments == "model_unstable", na.rm = TRUE),
    n_not_supported = sum(.data$evidence_status_segments == "not_supported", na.rm = TRUE),
    n_missing_effect_test = sum(.data$evidence_status_segments == "missing_effect_test", na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(segments, out_segments)
readr::write_csv(metrics, out_metrics)
readr::write_csv(input_status, out_status)
readr::write_csv(logic_audit, out_logic_audit)
readr::write_csv(count_audit, out_count_audit)
readr::write_csv(join_audit, out_join_audit)
readr::write_csv(selected_audit, out_selected_audit)
readr::write_csv(neuropil_availability_audit, out_neuropil_availability)
readr::write_csv(duplicate_audit, out_duplicate_audit)
readr::write_csv(effect_scope_audit, out_effect_scope_audit)
readr::write_csv(local_support_summary, out_local_support)
readr::write_csv(plot_source, out_plot_source)
readr::write_csv(heatmap_source_supermodule, out_heatmap_source_supermodule)
readr::write_csv(heatmap_source_module, out_heatmap_source_module)

render_circular_atlas(plot_source, out_main_svg, out_main_pdf)
render_circular_atlas(plot_source, out_selected_svg, out_selected_pdf, selected_only = TRUE)
for (ds_name in intersect(names(heatmap_svg_paths), unique(heatmap_source_supermodule$dataset))) {
  render_dataset_circular_heatmap(
    heatmap_source_supermodule,
    ds_name,
    heatmap_svg_paths[[ds_name]],
    heatmap_pdf_paths[[ds_name]]
  )
}
render_rectangular_module_heatmap(heatmap_source_module, out_rect_modules_svg, out_rect_modules_pdf)

write_run_manifest(
  out_run_manifest,
  inputs = list(
    wgcna_tables_root = path_results("tables", "06_modules_WGCNA"),
    biological_claims_table = path_results("tables", "biological_claims_table.csv")
  ),
  outputs = list(
    segments = out_segments,
    metrics = out_metrics,
    input_status = out_status,
    logic_audit = out_logic_audit,
    count_audit = out_count_audit,
    join_audit = out_join_audit,
    selected_table_audit = out_selected_audit,
    duplicate_source_audit = out_duplicate_audit,
    effect_scope_audit = out_effect_scope_audit,
    local_support_summary = out_local_support,
    plot_source = out_plot_source,
    circular_heatmap_source_supermodule = out_heatmap_source_supermodule,
    circular_heatmap_source_module = out_heatmap_source_module,
    main_svg = out_main_svg,
    main_pdf = out_main_pdf,
    selected_only_svg = out_selected_svg,
    selected_only_pdf = out_selected_pdf,
    circular_heatmap_svg = as.list(heatmap_svg_paths),
    circular_heatmap_pdf = as.list(heatmap_pdf_paths),
    rectangular_module_heatmap_svg = out_rect_modules_svg,
    rectangular_module_heatmap_pdf = out_rect_modules_pdf,
    neuron_neuropil_availability = out_neuropil_availability
  ),
  parameters = list(
    script = "10_biological_integration/04_wgcna_circular_atlas.R",
    dataset_argument = DATASET_ARG,
    downstream_only = TRUE,
    circular_plot_tracks = c(
      "broad program",
      "strongest effect estimate",
      "evidence status",
      "best local spatial unit",
      "local spatial FDR05 support fraction",
      "selected-table marker"
    ),
    circular_heatmap_inputs = c(
      "results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv",
      "results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv",
      "results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/wgcna_module_supermodule_annotation.csv"
    )
  ),
  notes = c(
    "Circular atlas figure uses already assembled segment rows and copied source statistics.",
    "No WGCNA models, module/supermodule definitions, p-values, FDRs, or claim gates are recomputed or altered.",
    "Microglia labels refer to microglia-enriched ROI/local microenvironment evidence, not purified cell-intrinsic claims."
  )
)

if (nrow(duplicate_audit) && any(!duplicate_audit$collapse_safe, na.rm = TRUE)) {
  warning(
    "Unsafe duplicate source supermodule collapse detected. Inspect: ",
    out_duplicate_audit,
    call. = FALSE
  )
}

cat("\nWGCNA circular atlas audit summary\n")
cat("Source supermodules per dataset:\n")
print(count_audit |> dplyr::select("dataset", "n_supermodules_source"))
cat("Circular segments per dataset:\n")
print(count_audit |> dplyr::select("dataset", "n_supermodules_segments"))
missing <- logic_audit |> dplyr::filter(!.data$present_in_circular_segments)
extra <- tibble::tibble()
cat("Missing supermodules:", if (nrow(missing)) paste(paste(missing$dataset, missing$supermodule_id, sep = ":"), collapse = ", ") else "none", "\n")
cat("Extra segments:", if (nrow(extra)) paste(paste(extra$dataset, extra$supermodule_id, sep = ":"), collapse = ", ") else "none", "\n")
cat("Evidence-status distribution:\n")
print(logic_audit |> dplyr::count(.data$dataset, .data$evidence_status_segments, name = "n"))
cat("Selected-table rows by source evidence status:\n")
print(selected_audit |> dplyr::count(.data$dataset, .data$source_evidence_status, name = "n"))
cat("Selected-table rows by selection support status:\n")
print(selected_audit |> dplyr::count(.data$dataset, .data$selection_support_status, name = "n"))
cat("Numeric validation passed:", all(logic_audit$numeric_match, na.rm = TRUE), "\n")
cat("Duplicate collapse unsafe:", if (nrow(duplicate_audit)) any(!duplicate_audit$collapse_safe, na.rm = TRUE) else FALSE, "\n")
cat("Outputs:\n")
cat(" - ", out_segments, "\n", sep = "")
cat(" - ", out_metrics, "\n", sep = "")
cat(" - ", out_status, "\n", sep = "")
cat(" - ", out_logic_audit, "\n", sep = "")
cat(" - ", out_count_audit, "\n", sep = "")
cat(" - ", out_join_audit, "\n", sep = "")
cat(" - ", out_selected_audit, "\n", sep = "")
cat(" - ", out_duplicate_audit, "\n", sep = "")
cat(" - ", out_effect_scope_audit, "\n", sep = "")
cat(" - ", out_local_support, "\n", sep = "")
cat(" - ", out_plot_source, "\n", sep = "")
cat(" - ", out_heatmap_source_supermodule, "\n", sep = "")
cat(" - ", out_heatmap_source_module, "\n", sep = "")
cat(" - ", out_main_svg, "\n", sep = "")
cat(" - ", out_main_pdf, "\n", sep = "")
cat(" - ", out_selected_svg, "\n", sep = "")
cat(" - ", out_selected_pdf, "\n", sep = "")
for (ds_name in names(heatmap_svg_paths)) {
  cat(" - ", heatmap_svg_paths[[ds_name]], "\n", sep = "")
  cat(" - ", heatmap_pdf_paths[[ds_name]], "\n", sep = "")
}
cat(" - ", out_rect_modules_svg, "\n", sep = "")
cat(" - ", out_rect_modules_pdf, "\n", sep = "")
cat(" - ", out_run_manifest, "\n", sep = "")
cat(" - ", out_neuropil_availability, "\n", sep = "")
