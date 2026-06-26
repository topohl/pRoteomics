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
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
if (!length(missing_pkgs)) {
  suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))
}

run <- wgcna_cli(allow_all = TRUE)
DATASET_ARG <- run$dataset

table_dir <- path_results("tables", "10_biological_integration", "wgcna_circular_atlas", "global")
source_dir <- path_results("source_data", "10_biological_integration", "wgcna_circular_atlas", "global")
report_dir <- path_results("reports", "10_biological_integration", "wgcna_circular_atlas", "global")
invisible(lapply(c(table_dir, source_dir, report_dir), dir_create))

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

metrics <- tibble::tibble(
  metric = c(
    "n_segments",
    "n_selected_table_rows",
    "n_datasets",
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
cat(" - ", out_neuropil_availability, "\n", sep = "")
