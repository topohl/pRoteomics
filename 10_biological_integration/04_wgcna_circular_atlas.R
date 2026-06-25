#!/usr/bin/env Rscript
# ================================================================
# Script: 10_biological_integration/04_wgcna_circular_atlas.R
# Stage: integration
# Scope: global
# Consumes: existing WGCNA module/supermodule summaries, labels,
#           biological annotations, group effects, and optional claims.
# Produces: read-only circular WGCNA atlas figures and source data.
# Notes: Manuscript visualization layer only. Does not recompute WGCNA,
#        alter module/supermodule definitions, group-effect statistics,
#        eigengenes, labels, or claim gates.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "wgcna_labeling_utils.R"))

SCRIPT_ID <- "10_biological_integration/04_wgcna_circular_atlas.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

required_pkgs <- c("circlize", "ComplexHeatmap", "grid", "dplyr", "tidyr", "tibble", "readr", "stringr", "scales", "svglite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "),
       ". Install circlize/ComplexHeatmap to render the circular atlas; no fallback plotting engine is used.",
       call. = FALSE)
}
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

run <- integration_cli(default_dataset = "all", allow_all = TRUE)
DATASETS <- integration_datasets(run$dataset)
paths <- integration_paths("wgcna_circular_atlas", "global")

input_spec <- function(ds) {
  base01 <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds)
  list(
    module_summary = list(path = file.path(base01, "modules", "WGCNA_module_summary.csv"), required = TRUE),
    module_definitions = list(path = file.path(base01, "modules", "WGCNA_module_definitions_for_downstream.csv"), required = TRUE),
    supermodule_summary = list(path = file.path(base01, "supermodules", "wgcna_supermodule_summary.csv"), required = TRUE),
    module_supermodule_annotation = list(path = file.path(base01, "supermodules", "wgcna_module_supermodule_annotation.csv"), required = TRUE),
    supermodule_clustering_sensitivity = list(path = file.path(base01, "supermodules", "supermodule_clustering_sensitivity.csv"), required = TRUE),
    supermodule_group_effects = list(path = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "supermodule_group_effects.csv"), required = TRUE),
    module_group_effects = list(path = path_results("tables", "06_modules_WGCNA", "group_effects", ds, "module_group_effects.csv"), required = TRUE),
    final_label_lookup = list(path = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_final_label_lookup.csv"), required = TRUE),
    supermodule_group_effects_interpretable = list(path = path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds, "WGCNA_supermodule_group_effects_interpretable.csv"), required = TRUE),
    supermodule_biological_annotation = list(path = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv"), required = TRUE),
    module_biological_annotation = list(path = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"), required = TRUE)
  )
}

global_inputs <- list(
  biological_claims_table = list(path = path_results("tables", "biological_claims_table.csv"), required = FALSE),
  evidence_priority_matrix = list(path = path_results("tables", "10_biological_integration", "evidence_priority_matrix", "global", "evidence_priority_matrix.csv"), required = FALSE)
)

all_input_paths <- c(
  unlist(lapply(DATASETS, function(ds) vapply(input_spec(ds), `[[`, character(1), "path")), use.names = TRUE),
  vapply(global_inputs, `[[`, character(1), "path")
)

if (run$dry_run) {
  invisible(lapply(unlist(paths), dir_create))
  dry_run_inputs(SCRIPT_ID, as.list(all_input_paths))
  dry_run_line("Output figures", paths$figures)
  dry_run_line("Output tables", paths$tables)
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
    message = if (exists) "loaded" else if (isTRUE(spec$required)) "required input missing" else "optional input missing; recorded as status"
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
    producer_script_or_artifact_id = "wgcna_circular_atlas",
    warning = if (!exists && isTRUE(spec$required)) status$message[[1]] else NA_character_
  )
  if (!exists) return(list(data = NULL, status = status))
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

clean_chr <- function(x, fallback = NA_character_) {
  x <- stringr::str_squish(as.character(x))
  x[is.na(x) | !nzchar(x) | toupper(x) %in% c("NA", "NAN", "NULL")] <- fallback
  x
}

label_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

dataset_label <- function(ds) {
  dplyr::case_when(
    ds == "microglia" ~ "Microglia ROI",
    ds == "neuron_neuropil" ~ "Neuron neuropil",
    ds == "neuron_soma" ~ "Neuron soma",
    TRUE ~ as.character(ds)
  )
}

program_levels <- c(
  "immune / microglia-associated ROI",
  "ECM / adhesion / vascular",
  "mitochondrial / energy metabolism",
  "RNA / translation / proteostasis",
  "synaptic / vesicle / cytoskeletal",
  "mixed / unresolved"
)

evidence_levels <- c("robust_FDR", "suggestive_FDR10", "nominal_only", "model_unstable", "not_supported", "missing_effect_test")

program_palette <- c(
  "immune / microglia-associated ROI" = "#C51B7D",
  "ECM / adhesion / vascular" = "#8C510A",
  "mitochondrial / energy metabolism" = "#6A51A3",
  "RNA / translation / proteostasis" = "#3182BD",
  "synaptic / vesicle / cytoskeletal" = "#238B45",
  "mixed / unresolved" = "#969696"
)

evidence_palette <- c(
  robust_FDR = "#2A9D8F",
  suggestive_FDR10 = "#8AB17D",
  nominal_only = "#E9C46A",
  model_unstable = "#F4A261",
  not_supported = "#B8B8B8",
  missing_effect_test = "#E5E5E5"
)

contrast_priority <- c("SUS - RES", "SUS - CON", "RES - CON")
evidence_rank <- setNames(seq_along(evidence_levels), evidence_levels)

broad_program_class <- function(label, dataset) {
  z <- tolower(as.character(label))
  out <- dplyr::case_when(
    grepl("microglia|immune|phago|lysosom|complement|inflamm|antigen|interferon|myeloid", z) ~ "immune / microglia-associated ROI",
    grepl("ecm|extracellular matrix|adhesion|collagen|laminin|basement|perivascular|vascular|bbb|endothelial|pericyte|blood vessel|mural", z) ~ "ECM / adhesion / vascular",
    grepl("mitochond|respirat|oxidative|\\batp\\b|tca|acetyl|electron transport|precursor metabolites", z) ~ "mitochondrial / energy metabolism",
    grepl("\\brna\\b|translation|ribosom|splice|proteostasis|proteasom|rnp|chaperon", z) ~ "RNA / translation / proteostasis",
    grepl("synap|vesicle|cytoskeleton|actin|microtubule|postsynap|presynap|axon|dendrit", z) ~ "synaptic / vesicle / cytoskeletal",
    TRUE ~ "mixed / unresolved"
  )
  out[is.na(out) | !nzchar(out)] <- "mixed / unresolved"
  factor(out, levels = program_levels)
}

source_label_rows <- function(summary, annotation, lookup) {
  ids <- clean_chr(summary$SupermoduleID)
  out <- tibble::tibble(supermodule_id = ids, cleaned_label = NA_character_, source_label_field = NA_character_, label_confidence = NA_character_)
  add_source <- function(current, df, id_col, label_candidates, prefix) {
    if (is.null(df) || !nrow(df) || !id_col %in% names(df)) return(current)
    labc <- label_col(df, label_candidates)
    if (is.na(labc)) return(current)
    confc <- label_col(df, c("label_confidence", "annotation_confidence", "Supermodule_LabelConfidence", "cleaned_biological_label_confidence", "Supermodule_CompositionConfidence"))
    tmp <- tibble::tibble(
      supermodule_id = clean_chr(df[[id_col]]),
      label = clean_chr(df[[labc]]),
      source = paste0(prefix, ":", labc),
      confidence = if (!is.na(confc)) clean_chr(df[[confc]]) else NA_character_
    ) |>
      dplyr::filter(!is.na(.data$supermodule_id), !is.na(.data$label)) |>
      dplyr::distinct(.data$supermodule_id, .keep_all = TRUE)
    idx <- match(current$supermodule_id, tmp$supermodule_id)
    fill <- is.na(current$cleaned_label) & !is.na(idx)
    current$cleaned_label[fill] <- tmp$label[idx[fill]]
    current$source_label_field[fill] <- tmp$source[idx[fill]]
    current$label_confidence[fill] <- tmp$confidence[idx[fill]]
    current
  }
  out <- add_source(out, lookup |> dplyr::filter(.data$level == "supermodule"), "entity_id",
                    c("final_plot_label", "best_data_driven_label", "parent_program"), "WGCNA_final_label_lookup.csv")
  out <- add_source(out, annotation, "SupermoduleID",
                    c("Supermodule_CompositionDisplayLabel", "Supermodule_CleanPlotLabel", "Supermodule_PlotLabel",
                      "cleaned_biological_label", "DominantMemberTheme", "Macroprogram_Display", "Supermodule_DisplayLabel"),
                    "WGCNA_supermodule_biological_annotation.csv")
  out <- add_source(out, summary, "SupermoduleID",
                    c("Supermodule_DisplayLabel", "Macroprogram_Display", "Supermodule_FinalLabel", "Supermodule_DataDrivenLabel"),
                    "wgcna_supermodule_summary.csv")
  out |>
    dplyr::mutate(
      cleaned_label = dplyr::coalesce(.data$cleaned_label, "mixed / unresolved"),
      source_label_field = dplyr::coalesce(.data$source_label_field, "fallback:mixed_unresolved")
    )
}

select_strongest_effect <- function(effects, super_ids) {
  if (is.null(effects) || !nrow(effects)) {
    return(tibble::tibble(
      supermodule_id = super_ids,
      strongest_contrast = NA_character_,
      strongest_effect_estimate = NA_real_,
      strongest_effect_direction = NA_character_,
      p_value = NA_real_,
      FDR_global = NA_real_,
      FDR_within_dataset_level = NA_real_,
      evidence_status = "missing_effect_test",
      effect_selection_status = "supermodule_group_effects_missing_or_empty"
    ))
  }
  scope <- if ("effect_scope" %in% names(effects) && any(effects$effect_scope == "spatial_adjusted_global", na.rm = TRUE)) {
    "spatial_adjusted_global"
  } else {
    NA_character_
  }
  eff <- effects
  if (!is.na(scope)) {
    eff <- eff |> dplyr::filter(.data$effect_scope == "spatial_adjusted_global")
    selection_status <- "spatial_adjusted_global"
  } else {
    eff <- eff |> dplyr::filter(!.data$effect_scope %in% c("within_spatial_unit", "stress_by_spatial_interaction"))
    if (!nrow(eff)) eff <- effects
    selection_status <- "fallback_non_spatial_adjusted_or_all_rows"
  }
  for (nm in c("supermodule_id", "contrast", "estimate", "p_value", "FDR_global", "FDR_within_dataset_level", "evidence_status", "direction")) {
    if (!nm %in% names(eff)) eff[[nm]] <- NA
  }
  eff |>
    dplyr::mutate(
      supermodule_id = clean_chr(.data$supermodule_id),
      contrast_priority = match(.data$contrast, contrast_priority),
      contrast_priority = ifelse(is.na(.data$contrast_priority), length(contrast_priority) + 1L, .data$contrast_priority),
      evidence_rank = evidence_rank[as.character(.data$evidence_status)],
      evidence_rank = ifelse(is.na(.data$evidence_rank), length(evidence_levels), .data$evidence_rank),
      estimate = suppressWarnings(as.numeric(.data$estimate)),
      p_value = suppressWarnings(as.numeric(.data$p_value)),
      FDR_global = suppressWarnings(as.numeric(.data$FDR_global)),
      FDR_within_dataset_level = suppressWarnings(as.numeric(.data$FDR_within_dataset_level))
    ) |>
    dplyr::filter(.data$supermodule_id %in% super_ids) |>
    dplyr::group_by(.data$supermodule_id) |>
    dplyr::arrange(.data$evidence_rank, .data$contrast_priority, .data$FDR_global, .data$FDR_within_dataset_level, .data$p_value, dplyr::desc(abs(.data$estimate)), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      supermodule_id,
      strongest_contrast = .data$contrast,
      strongest_effect_estimate = .data$estimate,
      strongest_effect_direction = dplyr::coalesce(clean_chr(.data$direction), dplyr::case_when(.data$estimate > 0 ~ "higher", .data$estimate < 0 ~ "lower", TRUE ~ "zero")),
      p_value,
      FDR_global,
      FDR_within_dataset_level,
      evidence_status = dplyr::coalesce(clean_chr(.data$evidence_status), "missing_effect_test"),
      effect_selection_status = selection_status
    ) |>
    dplyr::right_join(tibble::tibble(supermodule_id = super_ids), by = "supermodule_id") |>
    dplyr::mutate(
      evidence_status = dplyr::coalesce(.data$evidence_status, "missing_effect_test"),
      effect_selection_status = dplyr::coalesce(.data$effect_selection_status, "no_matching_effect_row")
    )
}

protein_counts_by_supermodule <- function(defs, ann) {
  if (is.null(defs) || !nrow(defs) || is.null(ann) || !nrow(ann)) return(tibble::tibble(supermodule_id = character(), n_proteins_if_available = integer()))
  mod_col <- label_col(defs, c("ModuleID", "module_id", "ModuleColor", "module_eigengene"))
  prot_col <- label_col(defs, c("ProteinID", "protein_id", "UniProt", "uniprot", "GeneSymbol", "gene_symbol"))
  ann_mod_col <- label_col(ann, c("ModuleID", "module_id", "module_eigengene", "ModuleColor"))
  ann_sm_col <- label_col(ann, c("SupermoduleID", "Supermodule_DataDrivenID", "Supermodule_DataDriven", "Supermodule"))
  if (any(is.na(c(mod_col, prot_col, ann_mod_col, ann_sm_col)))) {
    return(tibble::tibble(supermodule_id = character(), n_proteins_if_available = integer()))
  }
  defs |>
    dplyr::transmute(module_key = clean_chr(.data[[mod_col]]), protein_key = clean_chr(.data[[prot_col]])) |>
    dplyr::left_join(ann |> dplyr::transmute(module_key = clean_chr(.data[[ann_mod_col]]), supermodule_id = clean_chr(.data[[ann_sm_col]])) |> dplyr::distinct(), by = "module_key") |>
    dplyr::filter(!is.na(.data$supermodule_id), !is.na(.data$protein_key)) |>
    dplyr::group_by(.data$supermodule_id) |>
    dplyr::summarise(n_proteins_if_available = dplyr::n_distinct(.data$protein_key), .groups = "drop")
}

claims_for_supermodules <- function(claims, ds, ids) {
  if (is.null(claims) || !nrow(claims)) return(tibble::tibble(supermodule_id = ids, claim_allowed_if_available = NA_character_))
  cdf <- claims
  if ("dataset" %in% names(cdf)) cdf <- cdf |> dplyr::filter(.data$dataset == ds)
  if (!nrow(cdf)) return(tibble::tibble(supermodule_id = ids, claim_allowed_if_available = "no_dataset_claim_rows"))
  text <- apply(as.data.frame(cdf, stringsAsFactors = FALSE), 1, paste, collapse = " ")
  allowed <- if ("claim_allowed" %in% names(cdf)) as.character(cdf$claim_allowed) else NA_character_
  tibble::tibble(supermodule_id = ids) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      claim_allowed_if_available = {
        hit <- grepl(supermodule_id, text, fixed = TRUE)
        if (!any(hit)) "no_supermodule_specific_claim_row" else paste(unique(allowed[hit]), collapse = ";")
      }
    ) |>
    dplyr::ungroup()
}

load_one_dataset <- function(ds) {
  specs <- input_spec(ds)
  loaded <- lapply(names(specs), function(nm) read_input(ds, nm, specs[[nm]]))
  names(loaded) <- names(specs)
  status <- dplyr::bind_rows(lapply(loaded, `[[`, "status"))
  if (any(status$status == "missing_required")) return(list(status = status, missing_required = TRUE))
  data <- lapply(loaded, `[[`, "data")
  require_cols(data$supermodule_summary, c("SupermoduleID"), paste(ds, "supermodule_summary"))
  require_cols(data$module_summary, c("ModuleID"), paste(ds, "module_summary"))
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
input_status <- dplyr::bind_rows(global_status, dplyr::bind_rows(lapply(loaded, `[[`, "status"))) |>
  dplyr::transmute(dataset, input_name, path, required, status, n_rows = suppressWarnings(as.integer(.data$n_rows)), message)
write_csv_safe(input_status, file.path(paths$tables, "wgcna_circular_atlas_input_status.csv"))
write_csv_safe(input_status, file.path(paths$reports, "wgcna_circular_atlas_input_status.csv"))
write_csv_safe(input_status, file.path(paths$source_data, "wgcna_circular_atlas_input_status.csv"))
if (any(vapply(loaded, `[[`, logical(1), "missing_required"))) {
  stop("Missing required WGCNA circular atlas inputs. See ", file.path(paths$tables, "wgcna_circular_atlas_input_status.csv"), call. = FALSE)
}

dataset_bundle <- lapply(loaded, `[[`, "data")

segments <- dplyr::bind_rows(lapply(names(dataset_bundle), function(ds) {
  x <- dataset_bundle[[ds]]
  summary <- x$supermodule_summary
  size_col <- label_col(summary, c("n_modules", "DataDrivenClusterSize", "n_member_modules"))
  super_ids <- clean_chr(summary$SupermoduleID)
  label_rows <- source_label_rows(summary, x$supermodule_biological_annotation, x$final_label_lookup)
  effects <- select_strongest_effect(x$supermodule_group_effects, super_ids)
  prot <- protein_counts_by_supermodule(x$module_definitions, x$module_supermodule_annotation)
  claim_rows <- claims_for_supermodules(claims, ds, super_ids)
  tibble::tibble(
    dataset = ds,
    compartment_label = dataset_label(ds),
    supermodule_id = super_ids,
    supermodule_display_id = super_ids,
    n_member_modules = if (!is.na(size_col)) suppressWarnings(as.integer(summary[[size_col]])) else NA_integer_
  ) |>
    dplyr::left_join(label_rows, by = "supermodule_id") |>
    dplyr::left_join(prot, by = "supermodule_id") |>
    dplyr::left_join(effects, by = "supermodule_id") |>
    dplyr::left_join(claim_rows, by = "supermodule_id") |>
    dplyr::mutate(
      broad_program_class = as.character(broad_program_class(.data$cleaned_label, ds)),
      singleton_supermodule = .data$n_member_modules <= 1L,
      label_confidence_if_available = .data$label_confidence,
      microglia_roi_wording_status = dplyr::case_when(
        .data$dataset != "microglia" ~ "not_applicable",
        grepl("purified|cell-intrinsic", .data$cleaned_label, ignore.case = TRUE) ~ "unsafe_wording_detected",
        TRUE ~ "conservative_microglia_roi_local_microenvironment_wording"
      )
    ) |>
    dplyr::select(
      dataset, compartment_label, supermodule_id, supermodule_display_id, cleaned_label,
      broad_program_class, n_member_modules, n_proteins_if_available, singleton_supermodule,
      strongest_contrast, strongest_effect_estimate, strongest_effect_direction, p_value,
      FDR_global, FDR_within_dataset_level, evidence_status, claim_allowed_if_available,
      label_confidence_if_available, source_label_field, microglia_roi_wording_status,
      effect_selection_status
    )
}))

segments <- segments |>
  dplyr::mutate(
    evidence_status = ifelse(.data$evidence_status %in% evidence_levels, .data$evidence_status, "missing_effect_test"),
    broad_program_class = ifelse(.data$broad_program_class %in% program_levels, .data$broad_program_class, "mixed / unresolved"),
    evidence_order = evidence_rank[.data$evidence_status],
    abs_effect = abs(suppressWarnings(as.numeric(.data$strongest_effect_estimate))),
    n_member_modules = suppressWarnings(as.integer(.data$n_member_modules)),
    segment_index = dplyr::row_number()
  ) |>
  dplyr::arrange(factor(.data$dataset, levels = c("microglia", "neuron_neuropil", "neuron_soma")),
                 factor(.data$broad_program_class, levels = program_levels),
                 .data$evidence_order, dplyr::desc(.data$abs_effect), dplyr::desc(.data$n_member_modules), .data$supermodule_id) |>
  dplyr::group_by(.data$dataset) |>
  dplyr::mutate(seg_start = dplyr::row_number() - 1, seg_end = dplyr::row_number(), seg_mid = (.data$seg_start + .data$seg_end) / 2) |>
  dplyr::ungroup()

metrics <- dplyr::bind_rows(lapply(names(dataset_bundle), function(ds) {
  x <- dataset_bundle[[ds]]
  counts <- table(factor(segments$evidence_status[segments$dataset == ds], levels = evidence_levels))
  tibble::tibble(
    dataset = ds,
    n_proteins = dplyr::n_distinct(x$module_definitions[[label_col(x$module_definitions, c("ProteinID", "protein_id", "UniProt", "GeneSymbol"))]]),
    n_modules = dplyr::n_distinct(x$module_summary$ModuleID),
    n_supermodules = dplyr::n_distinct(x$supermodule_summary$SupermoduleID),
    n_singleton_supermodules = sum(segments$dataset == ds & segments$singleton_supermodule, na.rm = TRUE),
    fraction_singleton_supermodules = mean(segments$singleton_supermodule[segments$dataset == ds], na.rm = TRUE),
    n_robust_supermodule_effects = unname(counts["robust_FDR"]),
    n_suggestive_supermodule_effects = unname(counts["suggestive_FDR10"]),
    n_nominal_supermodule_effects = unname(counts["nominal_only"]),
    n_not_supported_supermodule_effects = unname(counts["not_supported"])
  )
}))

make_links <- function(seg) {
  eligible <- seg |>
    dplyr::filter(.data$broad_program_class != "mixed / unresolved", .data$evidence_status %in% c("robust_FDR", "suggestive_FDR10")) |>
    dplyr::group_by(.data$broad_program_class) |>
    dplyr::filter(dplyr::n_distinct(.data$dataset) >= 2L) |>
    dplyr::arrange(.data$evidence_order, dplyr::desc(.data$abs_effect), .by_group = TRUE) |>
    dplyr::slice_head(n = 3) |>
    dplyr::ungroup()
  if (!nrow(eligible)) {
    return(tibble::tibble(
      source_dataset = character(), source_supermodule_id = character(),
      target_dataset = character(), target_supermodule_id = character(),
      shared_program_class = character(), link_rule = character(),
      correlation_or_similarity_if_available = numeric(), evidence_status = character(),
      include_in_sparse_plot = logical()
    ))
  }
  out <- list()
  k <- 1L
  for (pc in unique(eligible$broad_program_class)) {
    d <- eligible |> dplyr::filter(.data$broad_program_class == pc)
    pairs <- utils::combn(seq_len(nrow(d)), 2)
    for (i in seq_len(ncol(pairs))) {
      a <- d[pairs[1, i], ]
      b <- d[pairs[2, i], ]
      if (a$dataset == b$dataset) next
      out[[k]] <- tibble::tibble(
        source_dataset = a$dataset,
        source_supermodule_id = a$supermodule_id,
        target_dataset = b$dataset,
        target_supermodule_id = b$supermodule_id,
        shared_program_class = pc,
        link_rule = "same_broad_program_class_and_both_robust_or_suggestive; max_sparse_links_per_class",
        correlation_or_similarity_if_available = NA_real_,
        evidence_status = paste(a$evidence_status, b$evidence_status, sep = ";"),
        include_in_sparse_plot = TRUE
      )
      k <- k + 1L
    }
  }
  dplyr::bind_rows(out) |>
    dplyr::group_by(.data$shared_program_class) |>
    dplyr::slice_head(n = 2) |>
    dplyr::ungroup()
}

links <- make_links(segments)
if (!nrow(links)) {
  links <- tibble::tibble(
    source_dataset = NA_character_, source_supermodule_id = NA_character_,
    target_dataset = NA_character_, target_supermodule_id = NA_character_,
    shared_program_class = NA_character_,
    link_rule = "no sparse high-confidence cross-compartment shared-program links passed filters",
    correlation_or_similarity_if_available = NA_real_,
    evidence_status = NA_character_,
    include_in_sparse_plot = FALSE
  )
}

validate_atlas <- function(seg) {
  numeric_qc <- dplyr::bind_rows(lapply(names(dataset_bundle), function(ds) {
    raw <- dataset_bundle[[ds]]$supermodule_group_effects |>
      dplyr::filter(.data$effect_scope == "spatial_adjusted_global") |>
      dplyr::mutate(supermodule_id = clean_chr(.data$supermodule_id)) |>
      dplyr::select(
        supermodule_id,
        strongest_contrast = "contrast",
        source_estimate = "estimate",
        source_p_value = "p_value",
        source_FDR_global = "FDR_global",
        source_FDR_within_dataset_level = "FDR_within_dataset_level",
        source_evidence_status = "evidence_status"
      )
    selected <- seg |> dplyr::filter(.data$dataset == ds, .data$effect_selection_status == "spatial_adjusted_global")
    keyed <- selected |>
      dplyr::left_join(raw, by = c("supermodule_id", "strongest_contrast"))
    validation_pairs <- c(
      strongest_effect_estimate = "source_estimate",
      p_value = "source_p_value",
      FDR_global = "source_FDR_global",
      FDR_within_dataset_level = "source_FDR_within_dataset_level"
    )
    dplyr::bind_rows(lapply(names(validation_pairs), function(out_nm) {
      source_nm <- unname(validation_pairs[[out_nm]])
      diff <- abs(suppressWarnings(as.numeric(keyed[[out_nm]])) - suppressWarnings(as.numeric(keyed[[source_nm]])))
      max_diff <- suppressWarnings(max(diff, na.rm = TRUE))
      if (!is.finite(max_diff)) max_diff <- 0
      missing_key <- sum(is.na(keyed[[source_nm]]) & !is.na(keyed[[out_nm]]))
      tibble::tibble(
        dataset = ds,
        numeric_field = out_nm,
        max_abs_diff = max_diff,
        missing_source_key_rows = missing_key,
        status = if (missing_key == 0L && max_diff <= 1e-12) "copied_from_existing_group_effects" else "mismatch"
      )
    }))
  }))
  if (any(numeric_qc$status == "mismatch")) stop("Validation failed: numeric group-effect fields were not copied exactly.", call. = FALSE)
  if (any(!seg$evidence_status %in% evidence_levels)) stop("Validation failed: invalid evidence_status value.", call. = FALSE)
  if (any(is.na(seg$source_label_field) | !nzchar(seg$source_label_field))) stop("Validation failed: every segment must have a source label field.", call. = FALSE)
  if (any(seg$dataset == "microglia" & seg$microglia_roi_wording_status == "unsafe_wording_detected")) {
    stop("Validation failed: unsafe microglia wording detected.", call. = FALSE)
  }
  write_csv_safe(numeric_qc, file.path(paths$reports, "wgcna_circular_atlas_numeric_validation.csv"))
  invisible(numeric_qc)
}
numeric_qc <- validate_atlas(segments)

segment_out <- segments |> dplyr::select(-"evidence_order", -"abs_effect", -"seg_start", -"seg_end", -"seg_mid")
write_csv_safe(segment_out, file.path(paths$tables, "wgcna_circular_atlas_segments.csv"))
write_csv_safe(segment_out, file.path(paths$source_data, "wgcna_circular_atlas_segments.csv"))
write_csv_safe(links, file.path(paths$tables, "wgcna_circular_atlas_links.csv"))
write_csv_safe(links, file.path(paths$source_data, "wgcna_circular_atlas_links.csv"))
write_csv_safe(metrics, file.path(paths$tables, "wgcna_circular_atlas_metrics.csv"))
write_csv_safe(metrics, file.path(paths$source_data, "wgcna_circular_atlas_metrics.csv"))

circos_setup <- function(seg) {
  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 90,
    gap.after = c(3, 3, 6),
    track.margin = c(0.002, 0.002),
    cell.padding = c(0, 0, 0, 0)
  )
  sectors <- c("microglia", "neuron_neuropil", "neuron_soma")
  xlim <- cbind(rep(0, length(sectors)), vapply(sectors, function(ds) max(seg$seg_end[seg$dataset == ds], na.rm = TRUE), numeric(1)))
  circlize::circos.initialize(factors = sectors, xlim = xlim)
}

draw_atlas <- function(seg, links_df = NULL, compact = FALSE, show_id_track = TRUE) {
  circos_setup(seg)
  max_modules <- max(seg$n_member_modules, na.rm = TRUE)
  if (!is.finite(max_modules) || max_modules <= 0) max_modules <- 1
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = if (compact) 0.06 else 0.075, bg.border = NA,
    panel.fun = function(x, y) {
      ds <- circlize::get.cell.meta.data("sector.index")
      lim <- circlize::get.cell.meta.data("xlim")
      circlize::circos.text(
        mean(lim), 0.42, dataset_label(ds),
        facing = "bending.outside", cex = if (compact) 0.78 else 0.72,
        niceFacing = TRUE, col = "#2F2F2F"
      )
    })
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = if (compact) 0.19 else 0.22, bg.border = NA,
    panel.fun = function(x, y) {
      ds <- circlize::get.cell.meta.data("sector.index")
      d <- seg |> dplyr::filter(.data$dataset == ds)
      for (i in seq_len(nrow(d))) {
        circlize::circos.rect(
          d$seg_start[i], 0, d$seg_end[i], 1,
          col = program_palette[d$broad_program_class[i]],
          border = "white", lwd = 0.18
        )
      }
    })
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = if (compact) 0.055 else 0.065, bg.border = NA,
    panel.fun = function(x, y) {
      ds <- circlize::get.cell.meta.data("sector.index")
      d <- seg |> dplyr::filter(.data$dataset == ds)
      for (i in seq_len(nrow(d))) {
        col <- scales::alpha(evidence_palette[d$evidence_status[i]], ifelse(d$evidence_status[i] %in% c("not_supported", "missing_effect_test"), 0.42, 0.82))
        circlize::circos.rect(d$seg_start[i], 0.18, d$seg_end[i], 0.82, col = col, border = "white", lwd = 0.16)
      }
    })
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = if (compact) 0.05 else 0.06, bg.border = NA,
    panel.fun = function(x, y) {
      ds <- circlize::get.cell.meta.data("sector.index")
      d <- seg |> dplyr::filter(.data$dataset == ds)
      for (i in seq_len(nrow(d))) {
        h <- 0.15 + 0.7 * pmin(d$n_member_modules[i] / max_modules, 1)
        circlize::circos.lines(c(d$seg_mid[i], d$seg_mid[i]), c(0.15, h), col = "#4A4A4A", lwd = 0.32)
        pch <- c(robust_FDR = 16, suggestive_FDR10 = 15, nominal_only = 4, model_unstable = 1, not_supported = 46, missing_effect_test = 3)[d$evidence_status[i]]
        circlize::circos.points(d$seg_mid[i], 0.92, pch = pch, cex = 0.23, col = "#404040")
      }
    })
  if (isTRUE(show_id_track)) {
    circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = if (compact) 0.07 else 0.09, bg.border = NA,
      panel.fun = function(x, y) {
        ds <- circlize::get.cell.meta.data("sector.index")
        d <- seg |> dplyr::filter(.data$dataset == ds)
        keep <- d |>
          dplyr::filter(.data$evidence_status %in% c("robust_FDR", "suggestive_FDR10") |
                          .data$n_member_modules >= stats::quantile(seg$n_member_modules, 0.8, na.rm = TRUE)) |>
          dplyr::slice_head(n = if (compact) 5 else 7)
        for (i in seq_len(nrow(keep))) {
          circlize::circos.text(
            keep$seg_mid[i], 0.45, keep$supermodule_display_id[i],
            cex = if (compact) 0.36 else 0.40,
            facing = "clockwise", niceFacing = TRUE, col = "#333333"
          )
        }
      })
  }
  if (!is.null(links_df) && nrow(links_df) && any(links_df$include_in_sparse_plot %in% TRUE)) {
    show_links <- links_df |> dplyr::filter(.data$include_in_sparse_plot %in% TRUE)
    for (i in seq_len(nrow(show_links))) {
      a <- seg |> dplyr::filter(.data$dataset == show_links$source_dataset[i], .data$supermodule_id == show_links$source_supermodule_id[i])
      b <- seg |> dplyr::filter(.data$dataset == show_links$target_dataset[i], .data$supermodule_id == show_links$target_supermodule_id[i])
      if (nrow(a) != 1L || nrow(b) != 1L) next
      circlize::circos.link(a$dataset, c(a$seg_mid, a$seg_mid), b$dataset, c(b$seg_mid, b$seg_mid),
                            col = scales::alpha(program_palette[a$broad_program_class], 0.22), border = NA, h.ratio = 0.58)
    }
  }
}

draw_label_table <- function(seg) {
  top <- seg |>
    dplyr::arrange(.data$evidence_order, dplyr::desc(.data$abs_effect), dplyr::desc(.data$n_member_modules)) |>
    dplyr::slice_head(n = 12) |>
    dplyr::mutate(
      dataset_table = dplyr::case_when(
        .data$dataset == "microglia" ~ "Microglia ROI",
        .data$dataset == "neuron_neuropil" ~ "Neuropil",
        .data$dataset == "neuron_soma" ~ "Soma",
        TRUE ~ .data$compartment_label
      ),
      evidence_table = dplyr::case_when(
        .data$evidence_status == "robust_FDR" ~ "robust",
        .data$evidence_status == "suggestive_FDR10" ~ "suggestive",
        .data$evidence_status == "nominal_only" ~ "nominal",
        .data$evidence_status == "model_unstable" ~ "unstable",
        .data$evidence_status == "not_supported" ~ "not supported",
        TRUE ~ "missing"
      ),
      short_label = stringr::str_trunc(.data$cleaned_label, width = 40, side = "right", ellipsis = "...")
    )
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::text(0.02, 0.965, "Selected supermodules", adj = c(0, 1), cex = 0.88, col = "#2F2F2F")
  headers <- c("", "Dataset", "ID", "Short cleaned label", "n", "Evidence")
  xpos <- c(0.03, 0.10, 0.27, 0.37, 0.78, 0.86)
  y_header <- 0.895
  for (j in seq_along(headers)) {
    graphics::text(xpos[j], y_header, headers[j], adj = c(0, 0.5), cex = 0.56, font = 2, col = "#333333")
  }
  graphics::segments(0.02, 0.865, 0.98, 0.865, col = "#D8D8D8", lwd = 0.7)
  row_h <- 0.062
  y <- 0.825
  for (i in seq_len(nrow(top))) {
    if (i %% 2 == 0) graphics::rect(0.015, y - row_h / 2, 0.985, y + row_h / 2, col = "#F7F7F7", border = NA)
    graphics::rect(0.027, y - 0.014, 0.052, y + 0.014, col = program_palette[top$broad_program_class[i]], border = NA)
    graphics::text(0.10, y, top$dataset_table[i], adj = c(0, 0.5), cex = 0.52, col = "#333333")
    graphics::text(0.27, y, top$supermodule_id[i], adj = c(0, 0.5), cex = 0.54, col = "#333333")
    graphics::text(0.37, y, top$short_label[i], adj = c(0, 0.5), cex = 0.45, col = "#333333")
    graphics::text(0.78, y, as.character(top$n_member_modules[i]), adj = c(0, 0.5), cex = 0.52, col = "#333333")
    graphics::text(0.86, y, top$evidence_table[i], adj = c(0, 0.5), cex = 0.45, col = "#333333")
    y <- y - row_h
  }
}

legend_grob <- function() {
  lgd1 <- ComplexHeatmap::Legend(
    title = "Broad program", at = names(program_palette),
    legend_gp = grid::gpar(fill = unname(program_palette)),
    ncol = 2, title_gp = grid::gpar(fontsize = 7), labels_gp = grid::gpar(fontsize = 5.8),
    grid_height = grid::unit(3.2, "mm"), grid_width = grid::unit(3.2, "mm")
  )
  lgd2 <- ComplexHeatmap::Legend(
    title = "Evidence/status", at = names(evidence_palette),
    legend_gp = grid::gpar(fill = scales::alpha(unname(evidence_palette), c(0.82, 0.82, 0.82, 0.82, 0.42, 0.42))),
    ncol = 2, title_gp = grid::gpar(fontsize = 7), labels_gp = grid::gpar(fontsize = 5.8),
    grid_height = grid::unit(3.2, "mm"), grid_width = grid::unit(3.2, "mm")
  )
  ComplexHeatmap::packLegend(lgd1, lgd2, direction = "horizontal", gap = grid::unit(5, "mm"))
}

save_base_figure <- function(stem, plot_fun, width = 7.2, height = 7.2) {
  svg_path <- file.path(paths$figures, paste0(stem, ".svg"))
  pdf_path <- file.path(paths$figures, paste0(stem, ".pdf"))
  png_path <- file.path(paths$figures, paste0(stem, ".png"))
  svglite::svglite(svg_path, width = width, height = height); plot_fun(); grDevices::dev.off()
  grDevices::cairo_pdf(pdf_path, width = width, height = height); plot_fun(); grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = 300, type = "cairo"); plot_fun(); grDevices::dev.off()
  c(svg = svg_path, pdf = pdf_path, png = png_path)
}

rendered <- c(
  save_base_figure("wgcna_circular_atlas_no_ribbons", function() draw_atlas(segments, NULL, compact = FALSE), 7.2, 7.2),
  save_base_figure("wgcna_circular_atlas_sparse_links", function() {
    draw_atlas(segments, links, compact = FALSE)
    if (!any(links$include_in_sparse_plot %in% TRUE)) {
      grid::grid.text(
        "No high-confidence cross-compartment links passed filters",
        x = grid::unit(0.5, "npc"), y = grid::unit(0.035, "npc"),
        gp = grid::gpar(fontsize = 7, col = "#4A4A4A")
      )
    }
  }, 7.2, 7.2),
  save_base_figure("wgcna_circular_atlas_compact_with_label_table", function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(fig = c(0.00, 0.64, 0.00, 1.00), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(0, 0, 0, 0))
    draw_atlas(segments, NULL, compact = TRUE, show_id_track = TRUE)
    graphics::par(fig = c(0.64, 1.00, 0.00, 1.00), mar = c(0.3, 0.1, 0.3, 0.2), new = TRUE)
    draw_label_table(segments)
  }, 14.0, 8.4),
  save_base_figure("wgcna_circular_atlas_compact_with_label_table_clean", function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(fig = c(0.00, 0.64, 0.00, 1.00), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(0, 0, 0, 0))
    draw_atlas(segments, NULL, compact = TRUE, show_id_track = TRUE)
    graphics::par(fig = c(0.64, 1.00, 0.00, 1.00), mar = c(0.3, 0.1, 0.3, 0.2), new = TRUE)
    draw_label_table(segments)
  }, 14.0, 8.4),
  save_base_figure("wgcna_circular_atlas_track_legend", function() {
    grid::grid.newpage()
    ComplexHeatmap::draw(legend_grob(), x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"), just = c("center", "center"))
  }, 6.2, 1.65)
)
circlize::circos.clear()

manifest_path <- file.path(paths$logs, "run_manifest.yml")
write_run_manifest(
  manifest_path,
  inputs = as.list(all_input_paths),
  outputs = list(figures = paths$figures, tables = paths$tables, source_data = paths$source_data, reports = paths$reports),
  parameters = list(dataset = run$dataset, effect_scope_priority = "spatial_adjusted_global", contrast_priority = paste(contrast_priority, collapse = ";")),
  notes = "Read-only global circular WGCNA atlas; estimates, p-values, FDRs, and evidence_status are copied from existing group-effect outputs."
)
invisible(file.copy(manifest_path, file.path(paths$tables, "run_manifest.yml"), overwrite = TRUE))

message("WGCNA circular atlas complete: ", paths$figures)
