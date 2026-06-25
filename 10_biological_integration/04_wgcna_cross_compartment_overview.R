#!/usr/bin/env Rscript
# ================================================================
# Script: 10_biological_integration/04_wgcna_cross_compartment_overview.R
# Stage: integration
# Scope: global
# Consumes: existing WGCNA module/supermodule summaries, downstream group effects,
#           cleaned labels, biological annotations, and optional claims table.
# Produces: conservative cross-compartment WGCNA/supermodule overview.
# Notes: Read-only reporting layer. Does not recompute WGCNA or alter module,
#        supermodule, model, FDR, effect-size, or claim-gate state.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "wgcna_labeling_utils.R"))

SCRIPT_ID <- "10_biological_integration/04_wgcna_cross_compartment_overview.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "readr", "stringr", "scales", "svglite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

run <- integration_cli(default_dataset = "all", allow_all = TRUE)
DATASETS <- integration_datasets(run$dataset)
paths <- integration_paths("wgcna_cross_compartment_overview", "global")

input_spec <- function(ds) {
  base01 <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds)
  list(
    module_summary = list(path = file.path(base01, "modules", "WGCNA_module_summary.csv"), required = TRUE),
    module_definitions = list(path = file.path(base01, "modules", "WGCNA_module_definitions_for_downstream.csv"), required = TRUE),
    supermodule_summary = list(path = file.path(base01, "supermodules", "wgcna_supermodule_summary.csv"), required = TRUE),
    supermodule_clustering_sensitivity = list(path = file.path(base01, "supermodules", "supermodule_clustering_sensitivity.csv"), required = FALSE),
    module_supermodule_annotation = list(path = file.path(base01, "supermodules", "wgcna_module_supermodule_annotation.csv"), required = TRUE),
    module_group_effects = list(path = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"), required = TRUE),
    supermodule_group_effects = list(path = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"), required = TRUE),
    final_label_lookup = list(path = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_final_label_lookup.csv"), required = TRUE),
    supermodule_group_effects_interpretable = list(path = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_supermodule_group_effects_interpretable.csv"), required = TRUE),
    module_biological_annotation = list(path = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"), required = FALSE),
    supermodule_biological_annotation = list(path = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"), required = FALSE)
  )
}

global_inputs <- list(
  biological_claims_table = list(path = path_results("tables", "biological_claims_table.csv"), required = FALSE),
  final_biological_evidence_bundle = list(path = path_results("tables", "10_biological_integration", "final_evidence_bundle", "global", "final_biological_evidence_bundle.xlsx"), required = FALSE)
)

all_input_paths <- c(
  unlist(lapply(DATASETS, function(ds) vapply(input_spec(ds), `[[`, character(1), "path")), use.names = TRUE),
  vapply(global_inputs, `[[`, character(1), "path")
)
if (run$dry_run) {
  invisible(lapply(unlist(paths), dir_create))
  dry_run_inputs(SCRIPT_ID, as.list(all_input_paths))
  quit(status = 0, save = "no")
}

read_input <- function(ds, input_name, spec) {
  path <- spec$path
  exists <- file.exists(path)
  status <- tibble::tibble(
    dataset = ds,
    input_name = input_name,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    required = isTRUE(spec$required),
    status = if (exists) "present" else if (isTRUE(spec$required)) "missing_required" else "missing_optional",
    n_rows = NA_integer_,
    message = if (exists) "loaded" else if (isTRUE(spec$required)) "required input missing" else "optional input missing"
  )
  record_input_resolution(
    script = SCRIPT_ID,
    dataset = ds,
    stage = "integration",
    input_name = input_name,
    expected_path = path,
    resolved_path = path,
    resolution_mode = status$status[[1]],
    strict_mode = strict_inputs_enabled(),
    allowed_in_strict_mode = TRUE,
    producer_script_or_artifact_id = "wgcna_cross_compartment_overview",
    warning = if (!exists && isTRUE(spec$required)) status$message[[1]] else NA_character_
  )
  if (!exists) return(list(data = NULL, status = status))
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("csv", "tsv", "txt")) {
    status$message <- "present; not parsed by this read-only CSV overview"
    return(list(data = NULL, status = status))
  }
  data <- tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), error = function(e) e)
  if (inherits(data, "error")) {
    status$status <- if (isTRUE(spec$required)) "missing_required" else "missing_optional"
    status$message <- paste("read error:", conditionMessage(data))
    return(list(data = NULL, status = status))
  }
  status$n_rows <- nrow(data)
  list(data = data, status = status)
}

require_cols <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) stop(label, " missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

label_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

clean_chr <- function(x, fallback = NA_character_) {
  x <- stringr::str_squish(as.character(x))
  x[is.na(x) | !nzchar(x) | toupper(x) %in% c("NA", "NAN", "NULL")] <- fallback
  x
}

dataset_label <- function(ds) vapply(ds, function(x) dataset_contracts()[[x]]$label %||% x, character(1))

effect_class_counts <- function(df) {
  statuses <- c("robust_FDR", "suggestive_FDR10", "nominal_only", "model_unstable", "not_supported")
  out <- stats::setNames(rep(0L, length(statuses)), statuses)
  if (!is.null(df) && nrow(df) && "evidence_status" %in% names(df)) {
    tab <- table(as.character(df$evidence_status))
    out[names(tab)] <- as.integer(tab[names(tab)])
  }
  out
}

selected_cut_height <- function(sens, ann) {
  if (!is.null(sens) && nrow(sens) && "selected_supermodule_cut_height" %in% names(sens)) {
    val <- suppressWarnings(as.numeric(stats::na.omit(sens$selected_supermodule_cut_height)[[1]]))
    if (is.finite(val)) return(val)
  }
  if (!is.null(ann) && nrow(ann) && "supermodule_merge_cut_height" %in% names(ann)) {
    val <- suppressWarnings(as.numeric(stats::na.omit(ann$supermodule_merge_cut_height)[[1]]))
    if (is.finite(val)) return(val)
  }
  NA_real_
}

supermodule_id_col <- function(df) label_col(df, c("supermodule_id", "SupermoduleID", "Supermodule_DataDriven", "Supermodule_DataDrivenID", "endpoint_id"))

broad_program_bin <- function(label, dataset = NA_character_) {
  z <- tolower(as.character(label))
  out <- dplyr::case_when(
    grepl("ecm|extracellular matrix|adhesion|collagen|laminin|basement|perivascular", z) ~ "ECM / adhesion",
    grepl("synap|vesicle|cytoskeleton|actin|microtubule|postsynap|presynap", z) ~ "synaptic / vesicle / cytoskeletal",
    grepl("mitochond|respirat|oxidative|\\batp\\b|tca|acetyl|electron transport", z) ~ "mitochondrial / energy metabolism",
    grepl("\\brna\\b|translation|ribosom|splice|proteostasis|proteasom|rnp|chaperon", z) ~ "RNA / translation / proteostasis",
    grepl("microglia|immune|phago|lysosom|complement|inflamm|antigen|interferon", z) ~ "immune / phagolysosomal / microglia-associated ROI",
    grepl("vascular|bbb|endothelial|pericyte|blood vessel|mural", z) ~ "vascular / BBB / microenvironment",
    TRUE ~ "mixed / unresolved"
  )
  out[dataset == "microglia" & out == "immune / phagolysosomal / microglia-associated ROI"] <- "immune / phagolysosomal / microglia-associated ROI"
  out
}

source_label_vector <- function(summary, annotation, lookup) {
  sid_sum <- supermodule_id_col(summary)
  ids <- clean_chr(if (!is.na(sid_sum)) summary[[sid_sum]] else character())
  out <- tibble::tibble(supermodule_id = ids, source_label = NA_character_, source_label_field = NA_character_, label_confidence = NA_character_)

  add_source <- function(current, df, id_candidates, label_candidates, field_prefix) {
    if (is.null(df) || !nrow(df)) return(current)
    idc <- supermodule_id_col(df)
    if (is.na(idc)) idc <- label_col(df, id_candidates)
    labc <- label_col(df, label_candidates)
    if (is.na(idc) || is.na(labc)) return(current)
    tmp <- tibble::tibble(
      supermodule_id = clean_chr(df[[idc]]),
      candidate = clean_chr(df[[labc]]),
      field = labc,
      conf = clean_chr(if ("label_confidence" %in% names(df)) df$label_confidence else if ("annotation_confidence" %in% names(df)) df$annotation_confidence else if ("Supermodule_CompositionConfidence" %in% names(df)) df$Supermodule_CompositionConfidence else NA_character_)
    ) |>
      dplyr::filter(!is.na(.data$supermodule_id), !is.na(.data$candidate)) |>
      dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)
    idx <- match(current$supermodule_id, tmp$supermodule_id)
    fill <- is.na(current$source_label) & !is.na(idx)
    current$source_label[fill] <- tmp$candidate[idx[fill]]
    current$source_label_field[fill] <- paste0(field_prefix, ":", tmp$field[idx[fill]])
    current$label_confidence[fill] <- tmp$conf[idx[fill]]
    current
  }

  out <- add_source(out, lookup, "entity_id", c("final_plot_label", "best_data_driven_label", "parent_program"), "WGCNA_final_label_lookup.csv")
  out <- add_source(out, annotation, "SupermoduleID", c("Supermodule_CompositionDisplayLabel", "Supermodule_CompositionLabel", "Macroprogram_Display", "DominantMemberTheme", "cleaned_biological_label", "module_program_primary", "microenvironment_caution_label"), "WGCNA_supermodule_biological_annotation.csv")
  out <- add_source(out, summary, "SupermoduleID", c("Macroprogram_Display", "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Supermodule_DataDrivenLabel"), "wgcna_supermodule_summary.csv")
  out
}

load_one_dataset <- function(ds) {
  specs <- input_spec(ds)
  loaded <- lapply(names(specs), function(nm) read_input(ds, nm, specs[[nm]]))
  names(loaded) <- names(specs)
  status <- dplyr::bind_rows(lapply(loaded, `[[`, "status"))
  if (any(status$status == "missing_required")) return(list(status = status, missing_required = TRUE))
  data <- lapply(loaded, `[[`, "data")

  require_cols(data$module_summary, c("ModuleID"), paste(ds, "module_summary"))
  require_cols(data$module_definitions, c("ModuleID"), paste(ds, "module_definitions"))
  require_cols(data$supermodule_summary, c("SupermoduleID"), paste(ds, "supermodule_summary"))
  require_cols(data$module_group_effects, c("evidence_status"), paste(ds, "module_group_effects"))
  require_cols(data$supermodule_group_effects, c("supermodule_id", "contrast", "estimate", "p_value", "FDR_global", "FDR_within_dataset_level", "evidence_status"), paste(ds, "supermodule_group_effects"))
  require_cols(data$final_label_lookup, c("dataset", "level", "entity_id", "final_plot_label"), paste(ds, "WGCNA_final_label_lookup"))
  wgcna_validate_label_lookup(data$final_label_lookup)

  list(status = status, missing_required = FALSE, data = data)
}

loaded_global <- lapply(names(global_inputs), function(nm) read_input("global", nm, global_inputs[[nm]]))
names(loaded_global) <- names(global_inputs)
global_status <- dplyr::bind_rows(lapply(loaded_global, `[[`, "status"))
claims <- loaded_global$biological_claims_table$data

loaded <- lapply(DATASETS, load_one_dataset)
names(loaded) <- DATASETS
input_status <- dplyr::bind_rows(global_status, dplyr::bind_rows(lapply(loaded, `[[`, "status")))
input_status <- input_status |>
  dplyr::transmute(
    dataset = as.character(.data$dataset),
    input_name = as.character(.data$input_name),
    path = as.character(.data$path),
    required = as.logical(.data$required),
    status = as.character(.data$status),
    n_rows = suppressWarnings(as.integer(.data$n_rows)),
    message = as.character(.data$message)
  )
invisible(write_integration_table(input_status, paths, "wgcna_cross_compartment_input_status.csv"))
if (any(vapply(loaded, `[[`, logical(1), "missing_required"))) {
  stop("Missing required WGCNA overview inputs. See ", file.path(paths$tables, "wgcna_cross_compartment_input_status.csv"), call. = FALSE)
}

dataset_bundle <- lapply(loaded, `[[`, "data")

metric_rows <- lapply(names(dataset_bundle), function(ds) {
  x <- dataset_bundle[[ds]]
  mod_size <- suppressWarnings(as.numeric(x$module_summary[[label_col(x$module_summary, c("n_features", "n_proteins", "n_mapped_features"))]]))
  sm_size <- suppressWarnings(as.numeric(x$supermodule_summary[[label_col(x$supermodule_summary, c("n_modules", "DataDrivenClusterSize", "n_member_modules"))]]))
  super_eff <- x$supermodule_group_effects
  module_eff <- x$module_group_effects
  s_counts <- effect_class_counts(super_eff)
  m_counts <- effect_class_counts(module_eff)
  claims_ds <- if (!is.null(claims) && nrow(claims) && "dataset" %in% names(claims)) claims[claims$dataset == ds, , drop = FALSE] else NULL
  claims_text <- if (!is.null(claims_ds)) tolower(paste(claims_ds$evidence_type %||% "", claims_ds$source_file %||% "", claims_ds$claim_type %||% "", sep = " ")) else character()
  wgcna_claim <- if (!is.null(claims_ds)) grepl("wgcna|module|supermodule", claims_text) else logical()
  allowed <- if (!is.null(claims_ds) && "claim_allowed" %in% names(claims_ds)) as.logical(claims_ds$claim_allowed) else rep(FALSE, length(wgcna_claim))
  n_animals <- suppressWarnings(as.integer(stats::na.omit(super_eff$n_animals_total)[[1]] %||% NA_integer_))
  n_samples <- suppressWarnings(as.integer(stats::na.omit(super_eff$n_samples_total)[[1]] %||% NA_integer_))
  tibble::tibble(
    dataset = ds,
    dataset_label = dataset_label(ds),
    n_samples_total = n_samples,
    n_animals_total = n_animals,
    n_animals_per_group = stats::na.omit(super_eff$n_animals_per_group)[[1]] %||% NA_character_,
    min_animals_per_group = suppressWarnings(as.integer(stats::na.omit(super_eff$min_animals_per_group)[[1]] %||% NA_integer_)),
    n_proteins_in_module_definitions = dplyr::n_distinct(x$module_definitions$ProteinID %||% seq_len(nrow(x$module_definitions))),
    n_wgcna_modules = dplyr::n_distinct(x$module_summary$ModuleID),
    median_module_size = stats::median(mod_size, na.rm = TRUE),
    min_module_size = suppressWarnings(min(mod_size, na.rm = TRUE)),
    max_module_size = suppressWarnings(max(mod_size, na.rm = TRUE)),
    n_supermodules = dplyr::n_distinct(x$supermodule_summary$SupermoduleID),
    n_singleton_supermodules = sum(sm_size == 1, na.rm = TRUE),
    fraction_singleton_supermodules = mean(sm_size == 1, na.rm = TRUE),
    median_supermodule_size = stats::median(sm_size, na.rm = TRUE),
    max_supermodule_size = suppressWarnings(max(sm_size, na.rm = TRUE)),
    selected_supermodule_cut_height = selected_cut_height(x$supermodule_clustering_sensitivity, x$module_supermodule_annotation),
    n_supermodule_effect_tests = nrow(super_eff),
    n_supermodule_effects_robust_FDR = unname(s_counts["robust_FDR"]),
    n_supermodule_effects_suggestive_FDR10 = unname(s_counts["suggestive_FDR10"]),
    n_supermodule_effects_nominal_only = unname(s_counts["nominal_only"]),
    n_supermodule_effects_model_unstable = unname(s_counts["model_unstable"]),
    n_supermodule_effects_not_supported = unname(s_counts["not_supported"]),
    n_module_effects_robust_FDR = unname(m_counts["robust_FDR"]),
    n_claim_allowed_wgcna_rows = if (!is.null(claims_ds)) sum(wgcna_claim & allowed, na.rm = TRUE) else NA_integer_,
    n_claim_blocked_or_diagnostic_wgcna_rows = if (!is.null(claims_ds)) sum(wgcna_claim & !allowed, na.rm = TRUE) else NA_integer_,
    main_limitation_note = if (ds == "microglia") "microglia-enriched ROI/local microenvironment; not purified or cell-intrinsic microglia evidence" else "supermodules are data-reduction/interpretation objects, not independent discoveries"
  )
})
metrics <- dplyr::bind_rows(metric_rows)

effect_overview <- dplyr::bind_rows(lapply(names(dataset_bundle), function(ds) {
  x <- dataset_bundle[[ds]]
  labels <- x$final_label_lookup |>
    dplyr::filter(.data$level == "supermodule") |>
    dplyr::transmute(supermodule_id = .data$entity_id, final_plot_label = .data$final_plot_label, label_source = "WGCNA_final_label_lookup.csv")
  x$supermodule_group_effects |>
    dplyr::filter(.data$effect_scope == "spatial_adjusted_global") |>
    dplyr::left_join(labels, by = "supermodule_id") |>
    dplyr::transmute(
      dataset = ds,
      supermodule_id,
      final_plot_label,
      contrast,
      effect_scope,
      estimate = as.numeric(.data$estimate),
      p_value = as.numeric(.data$p_value),
      FDR_global = as.numeric(.data$FDR_global),
      FDR_within_dataset_level = as.numeric(.data$FDR_within_dataset_level),
      evidence_status,
      direction,
      n_animals_total,
      n_samples_total,
      model_family,
      primary_model_stable,
      claim_allowed_model,
      microglia_roi_note = ifelse(ds == "microglia", "microglia-enriched ROI/local microenvironment; not purified or cell-intrinsic microglia evidence", NA_character_)
    )
}))

if (any(is.na(effect_overview$final_plot_label) | !nzchar(effect_overview$final_plot_label))) {
  stop("Displayed supermodule labels must come from WGCNA_final_label_lookup.csv; missing label(s) detected.", call. = FALSE)
}
if (any(effect_overview$dataset == "microglia" & !grepl("microglia-enriched ROI/local microenvironment", effect_overview$microglia_roi_note))) {
  stop("Microglia overview rows must use conservative ROI/local microenvironment wording.", call. = FALSE)
}

program_composition <- dplyr::bind_rows(lapply(names(dataset_bundle), function(ds) {
  x <- dataset_bundle[[ds]]
  src <- source_label_vector(x$supermodule_summary, x$supermodule_biological_annotation, x$final_label_lookup)
  size_col <- label_col(x$supermodule_summary, c("n_modules", "DataDrivenClusterSize", "n_member_modules"))
  size <- if (!is.na(size_col)) suppressWarnings(as.integer(x$supermodule_summary[[size_col]])) else rep(NA_integer_, nrow(x$supermodule_summary))
  id_col <- supermodule_id_col(x$supermodule_summary)
  tibble::tibble(supermodule_id = clean_chr(x$supermodule_summary[[id_col]]), n_member_modules = size) |>
    dplyr::left_join(src, by = "supermodule_id") |>
    dplyr::mutate(
      dataset = ds,
      broad_program_bin = broad_program_bin(.data$source_label, ds),
      confidence_status = dplyr::case_when(
        is.na(.data$source_label) ~ "unresolved_no_label",
        grepl("mixed|unresolved", .data$source_label, ignore.case = TRUE) ~ "mixed_or_unresolved",
        grepl("WGCNA_final_label_lookup", .data$source_label_field) ~ "cleaned_label_lookup",
        TRUE ~ "annotation_label"
      ),
      source_label = dplyr::coalesce(.data$source_label, "mixed / unresolved"),
      source_label_field_used = dplyr::coalesce(.data$source_label_field, "none_available")
    ) |>
    dplyr::select(dataset, supermodule_id, n_member_modules, broad_program_bin, source_label, source_label_field_used, label_confidence, confidence_status)
}))

numeric_qc <- dplyr::bind_rows(lapply(names(dataset_bundle), function(ds) {
  raw <- dataset_bundle[[ds]]$supermodule_group_effects |>
    dplyr::filter(.data$effect_scope == "spatial_adjusted_global") |>
    dplyr::arrange(.data$dataset, .data$supermodule_id, .data$contrast, .data$effect_scope)
  out <- effect_overview |>
    dplyr::filter(.data$dataset == ds) |>
    dplyr::arrange(.data$dataset, .data$supermodule_id, .data$contrast, .data$effect_scope)
  checks <- c("estimate", "p_value", "FDR_global", "FDR_within_dataset_level")
  dplyr::bind_rows(lapply(checks, function(nm) {
    a <- suppressWarnings(as.numeric(raw[[nm]]))
    b <- suppressWarnings(as.numeric(out[[nm]]))
    tibble::tibble(
      dataset = ds,
      numeric_field = nm,
      max_abs_diff = max(abs(a - b), na.rm = TRUE),
      status = if (identical(length(a), length(b)) && isTRUE(all.equal(a, b, tolerance = 0, check.attributes = FALSE))) "unchanged" else "changed"
    )
  }))
}))
if (any(numeric_qc$status != "unchanged")) stop("Numeric validation failed: overview changed copied group-effect values.", call. = FALSE)

invisible(write_integration_table(metrics, paths, "wgcna_cross_compartment_metrics.csv"))
invisible(write_integration_table(effect_overview, paths, "wgcna_supermodule_effect_overview.csv"))
invisible(write_integration_table(program_composition, paths, "wgcna_supermodule_program_composition.csv"))
write_csv_safe(numeric_qc, file.path(paths$reports, "numeric_validation.csv"))
write_csv_safe(numeric_qc, file.path(paths$source_data, "wgcna_cross_compartment_numeric_validation.csv"))

theme_overview <- function(base_size = 7) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(color = "#333333"),
      strip.text = ggplot2::element_text(face = "plain", color = "#333333"),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = base_size + 1, face = "plain"),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}
dataset_levels <- dataset_label(c("microglia", "neuron_soma", "neuron_neuropil"))
save_plot_both <- function(plot, stem, width = 180, height = 95) {
  ggplot2::ggsave(file.path(paths$figures, paste0(stem, ".svg")), plot = plot, width = width, height = height, units = "mm", device = svglite::svglite)
  ggplot2::ggsave(file.path(paths$figures, paste0(stem, ".pdf")), plot = plot, width = width, height = height, units = "mm", device = grDevices::cairo_pdf)
}
wrap_lab <- function(x, width = 38) vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))

arch_source <- metrics |>
  dplyr::mutate(dataset_label = factor(.data$dataset_label, levels = dataset_levels)) |>
  tidyr::pivot_longer(c("n_wgcna_modules", "n_supermodules", "n_singleton_supermodules"), names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(.data$metric, levels = c("n_wgcna_modules", "n_supermodules", "n_singleton_supermodules"), labels = c("WGCNA modules", "Supermodules", "Singleton supermodules")))
write_csv_safe(arch_source, file.path(paths$source_data, "wgcna_cross_compartment_overview_architecture_source.csv"))
p1 <- ggplot2::ggplot(arch_source, ggplot2::aes(x = .data$dataset_label, y = .data$value, fill = .data$metric)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.62, color = "white", linewidth = 0.15) +
  ggplot2::geom_text(
    data = metrics |> dplyr::mutate(dataset_label = factor(.data$dataset_label, levels = dataset_levels)),
    ggplot2::aes(x = .data$dataset_label, y = pmax(.data$n_wgcna_modules, .data$n_supermodules, na.rm = TRUE) * 1.08, label = paste0("proteins=", .data$n_proteins_in_module_definitions, "\nmedian module=", round(.data$median_module_size))),
    inherit.aes = FALSE, size = 2.1, color = "#333333"
  ) +
  ggplot2::scale_fill_manual(values = c("WGCNA modules" = "#3B6FB6", "Supermodules" = "#2A9D8F", "Singleton supermodules" = "#8A8A84")) +
  ggplot2::labs(title = "WGCNA module-to-supermodule overview", x = NULL, y = "Count", fill = NULL) +
  theme_overview()
save_plot_both(p1, "wgcna_cross_compartment_overview_architecture", 178, 92)

heat_source <- effect_overview |>
  dplyr::mutate(
    dataset_label = factor(dataset_label(.data$dataset), levels = dataset_levels),
    contrast = factor(.data$contrast, levels = c("RES - CON", "SUS - CON", "SUS - RES")),
    final_plot_label_wrapped = wrap_lab(.data$final_plot_label, 34),
    evidence_marker = dplyr::case_when(
      .data$evidence_status == "robust_FDR" ~ "FDR <= 0.05",
      .data$evidence_status == "suggestive_FDR10" ~ "FDR <= 0.10",
      .data$evidence_status == "nominal_only" ~ "nominal",
      .data$evidence_status == "model_unstable" ~ "model unstable",
      TRUE ~ "not supported"
    )
  )
write_csv_safe(heat_source, file.path(paths$source_data, "wgcna_supermodule_effect_heatmap_global_source.csv"))
max_est <- max(abs(heat_source$estimate), na.rm = TRUE)
max_est <- if (is.finite(max_est)) max(0.1, stats::quantile(abs(heat_source$estimate), 0.95, na.rm = TRUE)) else 1
p2 <- ggplot2::ggplot(heat_source, ggplot2::aes(x = .data$contrast, y = .data$final_plot_label_wrapped, fill = .data$estimate)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.18) +
  ggplot2::geom_point(ggplot2::aes(shape = .data$evidence_marker), size = 1.45, stroke = 0.25, color = "#222222", fill = "white") +
  ggplot2::facet_grid(dataset_label ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "#F8FAFC", high = "#C84C5A", midpoint = 0, limits = c(-max_est, max_est), oob = scales::squish) +
  ggplot2::scale_shape_manual(values = c("FDR <= 0.05" = 21, "FDR <= 0.10" = 22, "nominal" = 4, "model unstable" = 1, "not supported" = 46), drop = FALSE) +
  ggplot2::labs(title = "Global spatial-adjusted supermodule effects", x = NULL, y = NULL, fill = "Estimate", shape = NULL) +
  theme_overview(6.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), axis.text.y = ggplot2::element_text(size = 5.5))
save_plot_both(p2, "wgcna_supermodule_effect_heatmap_global", 183, 150)

program_source <- program_composition |>
  dplyr::mutate(dataset_label = factor(dataset_label(.data$dataset), levels = dataset_levels)) |>
  dplyr::group_by(.data$dataset, .data$dataset_label, .data$broad_program_bin) |>
  dplyr::summarise(n_supermodules = dplyr::n(), n_member_modules = sum(.data$n_member_modules, na.rm = TRUE), .groups = "drop")
write_csv_safe(program_source, file.path(paths$source_data, "wgcna_supermodule_program_composition_source.csv"))
program_cols <- c(
  "ECM / adhesion" = "#8C510A",
  "synaptic / vesicle / cytoskeletal" = "#238B45",
  "mitochondrial / energy metabolism" = "#6A51A3",
  "RNA / translation / proteostasis" = "#3182BD",
  "immune / phagolysosomal / microglia-associated ROI" = "#C51B7D",
  "vascular / BBB / microenvironment" = "#E6550D",
  "mixed / unresolved" = "#969696"
)
p3 <- ggplot2::ggplot(program_source, ggplot2::aes(x = .data$dataset_label, y = .data$n_member_modules, fill = .data$broad_program_bin)) +
  ggplot2::geom_col(width = 0.66, color = "white", linewidth = 0.16) +
  ggplot2::scale_fill_manual(values = program_cols, drop = FALSE) +
  ggplot2::labs(title = "Supermodule broad program composition", x = NULL, y = "Member modules", fill = NULL) +
  theme_overview(6.8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
save_plot_both(p3, "wgcna_supermodule_program_composition", 178, 100)

quality_source <- metrics |>
  dplyr::select(dataset, dataset_label, n_supermodule_effects_robust_FDR, n_supermodule_effects_suggestive_FDR10, n_supermodule_effects_nominal_only, n_supermodule_effects_model_unstable, n_supermodule_effects_not_supported) |>
  tidyr::pivot_longer(-c("dataset", "dataset_label"), names_to = "evidence_status", values_to = "n_tests") |>
  dplyr::mutate(
    dataset_label = factor(.data$dataset_label, levels = dataset_levels),
    evidence_status = factor(.data$evidence_status, levels = c("n_supermodule_effects_robust_FDR", "n_supermodule_effects_suggestive_FDR10", "n_supermodule_effects_nominal_only", "n_supermodule_effects_model_unstable", "n_supermodule_effects_not_supported"), labels = c("robust_FDR", "suggestive_FDR10", "nominal_only", "model_unstable", "not_supported"))
  )
write_csv_safe(quality_source, file.path(paths$source_data, "wgcna_evidence_quality_summary_source.csv"))
p4 <- ggplot2::ggplot(quality_source, ggplot2::aes(x = .data$dataset_label, y = .data$n_tests, fill = .data$evidence_status)) +
  ggplot2::geom_col(width = 0.66, color = "white", linewidth = 0.16) +
  ggplot2::scale_fill_manual(values = c(robust_FDR = "#2A9D8F", suggestive_FDR10 = "#8AB17D", nominal_only = "#E9C46A", model_unstable = "#F4A261", not_supported = "#B8B8B8"), drop = FALSE) +
  ggplot2::labs(title = "Supermodule evidence quality summary", x = NULL, y = "Effect tests", fill = NULL) +
  theme_overview(6.8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
save_plot_both(p4, "wgcna_evidence_quality_summary", 178, 95)

manifest <- list(
  script = SCRIPT_ID,
  script_version = "1.0.0",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  datasets = DATASETS,
  read_only_note = "Read-only downstream reporting layer; does not recompute WGCNA or alter module assignments, supermodule assignments, labels, estimates, p-values, FDRs, or claim gates.",
  inputs = c(lapply(DATASETS, function(ds) lapply(input_spec(ds), `[[`, "path")), list(global = lapply(global_inputs, `[[`, "path"))),
  outputs = list(tables = paths$tables, source_data = paths$source_data, figures = paths$figures, logs = paths$logs, reports = paths$reports)
)
manifest_path <- file.path(paths$tables, "wgcna_cross_compartment_overview_manifest.yml")
if (requireNamespace("yaml", quietly = TRUE)) writeLines(yaml::as.yaml(manifest), manifest_path) else capture.output(str(manifest, max.level = 4), file = manifest_path)
invisible(file.copy(manifest_path, file.path(paths$logs, "wgcna_cross_compartment_overview_manifest.yml"), overwrite = TRUE))

write_run_manifest(
  file.path(paths$logs, "run_manifest.yml"),
  inputs = as.list(all_input_paths),
  outputs = list(tables = paths$tables, source_data = paths$source_data, figures = paths$figures, reports = paths$reports, manifest = manifest_path),
  parameters = list(dataset = run$dataset, effect_scope = "spatial_adjusted_global"),
  notes = "Conservative read-only WGCNA/supermodule cross-compartment overview. Supermodules are data-reduction/interpretation objects, not independent discoveries."
)

message("WGCNA cross-compartment overview complete: ", paths$tables)
