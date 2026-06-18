#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/06_annotate_module_microenvironment.r
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/; optional config/marker_panels/wgcna_reference_marker_sets.csv; results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv; +1 more.
# Produces: results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv; results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv; results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_targeted_signature_overlap_details.csv.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Consumes marker registry, empirical markers, WGCNA outputs, and neuropil annotations when present.
# ================================================================

#
# Annotate WGCNA modules and supermodules with biological and microenvironment evidence.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "readr", "stringr", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

run <- wgcna_cli()
DATASET <- run$dataset
PATHS <- wgcna_downstream_paths("module_annotation", DATASET)
FILES <- resolve_wgcna_files(DATASET)
force_microglia <- tolower(Sys.getenv("PROTEOMICS_FORCE_MICROGLIA_MODULE_ANNOTATION", unset = "false")) %in% c("1", "true", "yes")
classification_threshold <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_WGCNA_MARKER_FRACTION_THRESHOLD", unset = "0.10")))
if (!is.finite(classification_threshold)) classification_threshold <- 0.10
marker_threshold_sensitivity_values <- c(0.05, 0.10, 0.20)
supplemental_marker_panel_file <- Sys.getenv(
  "PROTEOMICS_WGCNA_MICROENV_MARKER_PANEL_FILE",
  unset = repo_path("config", "marker_panels", "microenvironment_marker_panels.csv")
)
supplemental_marker_panel_file <- normalizePath(supplemental_marker_panel_file, winslash = "/", mustWork = FALSE)
supplemental_marker_panel_hash <- if (file.exists(supplemental_marker_panel_file)) unname(tools::md5sum(supplemental_marker_panel_file)) else NA_character_

as_count0 <- function(x) {
  y <- suppressWarnings(as.integer(x %||% 0L))
  y[is.na(y)] <- 0L
  y
}

read_microenvironment_marker_panels <- function(path = supplemental_marker_panel_file) {
  required <- c(
    "panel_id", "marker_symbol", "source_type", "source_reference",
    "allowed_use", "claim_role", "caution_note"
  )
  panels <- safe_read_csv(path)
  if (is.null(panels) || !nrow(panels)) {
    stop("Missing supplemental microenvironment marker panel config: ", path, call. = FALSE)
  }
  missing <- setdiff(required, names(panels))
  if (length(missing)) {
    stop("Supplemental microenvironment marker panel config is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  panels <- panels |>
    dplyr::mutate(
      panel_id = as.character(.data$panel_id),
      marker_symbol = as.character(.data$marker_symbol),
      source_type = as.character(.data$source_type),
      allowed_use = as.character(.data$allowed_use),
      claim_role = as.character(.data$claim_role),
      gene_token = normalize_gene_token(.data$marker_symbol)
    ) |>
    dplyr::filter(nzchar(.data$panel_id), nzchar(.data$gene_token))
  valid_source <- c("external_reference", "empirical_this_study", "curated_local", "diagnostic")
  valid_use <- c("QC_only", "annotation", "claim_support")
  valid_role <- c("primary_support", "orthogonal_support", "context_only", "caution_only")
  if (any(!panels$source_type %in% valid_source)) stop("Invalid source_type in supplemental marker panel config.", call. = FALSE)
  if (any(!panels$allowed_use %in% valid_use)) stop("Invalid allowed_use in supplemental marker panel config.", call. = FALSE)
  if (any(!panels$claim_role %in% valid_role)) stop("Invalid claim_role in supplemental marker panel config.", call. = FALSE)
  panels
}

if (run$dry_run) {
  invisible(lapply(unlist(PATHS), dir_create))
  dry_run_line("Script", "06_modules_WGCNA/06_annotate_module_microenvironment.r")
  dry_run_line("Dataset", DATASET)
  dry_run_line("Module definitions", FILES$definitions, if (file.exists(FILES$definitions)) "PASS" else "WARN")
  dry_run_line("GO enrichment", FILES$go, if (file.exists(FILES$go)) "PASS" else "WARN")
  dry_run_line("Supermodule annotation", FILES$supermodule_annotation, if (file.exists(FILES$supermodule_annotation)) "PASS" else "WARN")
  dry_run_line("Group effects", path_results("tables", "06_modules_WGCNA", "group_effects", DATASET), "INFO")
  if (DATASET == "microglia" || force_microglia) dry_run_line("Neuropil reference annotation", FILES$neuropil_annotation, if (file.exists(FILES$neuropil_annotation)) "PASS" else "WARN")
  if (DATASET == "microglia" || force_microglia) dry_run_line("Targeted microglia signature enrichment", path_results("tables", "04_differential_expression_enrichment", "microglia_targeted_signature_enrichment", "microglia", "microglia_signature_enrichment_with_contrast_class.csv"), if (file.exists(path_results("tables", "04_differential_expression_enrichment", "microglia_targeted_signature_enrichment", "microglia", "microglia_signature_enrichment_with_contrast_class.csv"))) "PASS" else "WARN")
  dry_run_line("Marker registry", Sys.getenv("PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE", unset = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")), "INFO")
  dry_run_line("Empirical ROI marker sets", Sys.getenv("PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE", unset = path_results("tables", "03_qc_exploration", "05_empirical_roi_marker_discovery", "empirical_roi_marker_sets.csv")), "INFO")
  dry_run_line("Supplemental microenvironment marker panels", supplemental_marker_panel_file, if (file.exists(supplemental_marker_panel_file)) "PASS" else "FAIL")
  dry_run_line("Microenvironment threshold sensitivity audit", path_results("reviewer_audit", "wgcna_microenvironment_threshold_sensitivity.csv"), "INFO")
  dry_run_line("Output tables", PATHS$tables)
  quit(status = 0, save = "no")
}

if (DATASET != "microglia" && force_microglia) {
  warning("Forcing microglia-specific module annotation outside dataset == microglia.", call. = FALSE)
}

definitions <- safe_read_csv(FILES$definitions)
if (is.null(definitions) || !nrow(definitions)) {
  module_annot <- data.frame(dataset = DATASET, ModuleID = NA_character_, ModuleColor = NA_character_, n_proteins = 0L, microenvironment_class = "missing_module_definitions", interpretation_note = WGCNA_ROI_NOTE)
  super_annot <- data.frame(dataset = DATASET, SupermoduleID = NA_character_, supermodule_id = NA_character_, n_member_modules = 0L, dominant_microenvironment_class = "missing_module_definitions", interpretation_note = WGCNA_ROI_NOTE)
  write_table_and_source(module_annot, PATHS$tables, PATHS$source_data, "WGCNA_module_biological_annotation.csv")
  write_table_and_source(super_annot, PATHS$tables, PATHS$source_data, "WGCNA_supermodule_biological_annotation.csv")
  write_run_manifest(file.path(PATHS$logs, "run_manifest.yml"), inputs = FILES, outputs = list(tables = PATHS$tables), parameters = list(dataset = DATASET), notes = paste(WGCNA_ROI_NOTE, "Improved annotation separates vascular basement membrane/BBB/mural, astrocyte/endfoot, oligodendrocyte/myelin, neuropil, and microglia-supported ROI evidence; labels are annotation only."))
  quit(status = 0, save = "no")
}

validate_wgcna_module_definitions(definitions, "WGCNA downstream definitions")
marker_sets <- load_wgcna_marker_sets()
marker_registry_version <- attr(marker_sets, "marker_registry_version") %||% NA_character_
empirical_marker_set_version <- attr(marker_sets, "empirical_marker_set_version") %||% NA_character_
supplemental_marker_panels <- read_microenvironment_marker_panels(supplemental_marker_panel_file)
supplemental_marker_sets <- split(as.character(supplemental_marker_panels$marker_symbol), as.character(supplemental_marker_panels$panel_id)) |>
  lapply(function(x) unique(x[nzchar(normalize_gene_token(x))]))
module_summary <- safe_read_csv(FILES$module_summary)
go <- safe_read_csv(FILES$go)
super_ann <- safe_read_csv(FILES$supermodule_annotation)
super_summary <- safe_read_csv(FILES$supermodule_summary)
module_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", DATASET, "module_group_effects.csv"))
super_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", DATASET, "supermodule_group_effects.csv"))
neuropil_ref <- if (DATASET == "microglia" || force_microglia) safe_read_csv(FILES$neuropil_annotation) else NULL
targeted_signature_file <- path_results("tables", "04_differential_expression_enrichment", "microglia_targeted_signature_enrichment", "microglia", "microglia_signature_enrichment_with_contrast_class.csv")
targeted_signature_ref <- if (DATASET == "microglia" || force_microglia) safe_read_csv(targeted_signature_file) else NULL

for (nm in names(supplemental_marker_sets)) {
  if (!nm %in% names(marker_sets)) marker_sets[[nm]] <- supplemental_marker_sets[[nm]]
}
panel_names <- names(marker_sets)

compact_label <- function(x, n = 3L) {
  x <- unique(trimws(as.character(x)))
  x <- x[nzchar(x) & !is.na(x)]
  paste(utils::head(x, n), collapse = "; ")
}

clean_term_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- gsub("\\bGO\\b", "", x, ignore.case = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

generic_supermodule_label <- function(x) {
  x <- tolower(as.character(x))
  is.na(x) | !nzchar(x) |
    grepl("hub-supported|module cluster|ambiguous|unlabelled|unknown|mixed$", x)
}

informative_member_label <- function(x) {
  x <- stringr::str_squish(as.character(x))
  bad <- is.na(x) | !nzchar(x) |
    grepl("mixed|low-specificity|unresolved|unlabelled|unassigned|unknown|^module\\b|^ME[A-Za-z0-9]+$", x, ignore.case = TRUE)
  x[bad] <- NA_character_
  x
}

first_informative_value <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- informative_member_label(vals)
  vals <- vals[!is.na(vals)]
  if (length(vals)) vals[[1]] else NA_character_
}

split_label_terms <- function(x) {
  wgcna_split_label_terms(x)
}

member_theme_from_label <- function(x) {
  z <- tolower(clean_term_label(x))
  cleaned <- wgcna_clean_semantic_label(raw_label = x)
  dplyr::case_when(
    is.na(z) | !nzchar(z) ~ NA_character_,
    !is.na(cleaned$cleaned_biological_label) & nzchar(cleaned$cleaned_biological_label) ~ cleaned$cleaned_biological_label,
    TRUE ~ shorten_supermodule_label(clean_term_label(x), max_chars = 42)
  )
}

semantic_marker_panel_values <- function(row) {
  candidate_cols <- grep("(_marker_fraction|_fraction)$", names(row), value = TRUE)
  candidate_cols <- setdiff(candidate_cols, c(
    "fraction_modules_neuropil_sensitive", "fraction_modules_vascular_basement_membrane_ecm",
    "fraction_modules_vascular_bbb_mural", "fraction_modules_astrocyte_or_endfoot_sensitive",
    "fraction_modules_oligodendrocyte_or_myelin_sensitive", "fraction_modules_peripheral_myeloid_caution",
    "fraction_modules_ambiguous_or_mixed", "fraction_member_modules_with_informative_labels"
  ))
  if (!length(candidate_cols)) return(stats::setNames(numeric(), character()))
  vals <- suppressWarnings(as.numeric(unlist(row[candidate_cols], use.names = FALSE)))
  stats::setNames(vals, candidate_cols)
}

composition_short_label <- function(x) {
  x <- as.character(x)
  x <- sub("^mostly\\s+", "", x, ignore.case = TRUE)
  x <- sub("^mixed:\\s*", "mixed: ", x, ignore.case = TRUE)
  shorten_supermodule_label(x, max_chars = 38)
}

microglia_guarded_composition_label <- function(label, dominant_class) {
  label
}

label_from_evidence <- function(row) {
  frac <- function(candidates) {
    hit <- candidates[candidates %in% names(row)][1]
    if (length(hit) && !is.na(hit)) suppressWarnings(as.numeric(row[[hit]])) else NA_real_
  }
  max_frac <- function(...) {
    vals <- c(...)
    vals <- vals[is.finite(vals)]
    if (length(vals)) max(vals) else NA_real_
  }

  bm <- max_frac(frac(c("basement_membrane_perivascular_ecm_fraction", "basement_membrane_perivascular_ecm_marker_fraction")))
  bbb <- max_frac(frac(c("bbb_endothelial_transport_fraction", "bbb_endothelial_transport_marker_fraction")), frac(c("canonical_endothelial_vascular_fraction", "endothelial_pericyte_vascular_marker_fraction")))
  mural <- max_frac(frac(c("pericyte_mural_ng2_fraction", "pericyte_mural_ng2_marker_fraction")), frac(c("canonical_pericyte_vascular_fraction")))
  adhesion <- max_frac(frac(c("integrin_ecm_adhesion_fraction", "integrin_ecm_adhesion_marker_fraction")))
  micro <- max_frac(frac(c("empirical_microglia_roi_high_confidence_fraction", "empirical_microglia_roi_enriched_fraction")), frac(c("canonical_microglia_homeostatic_fraction", "microglia_marker_fraction")))
  micro_state <- frac(c("canonical_microglia_phagolysosomal_state_fraction"))
  targeted_driver <- as.character(row$targeted_signature_primary_driver %||% NA_character_)
  targeted_micro <- as_count0(row$n_targeted_microglia_signature_overlaps)
  targeted_claim <- as_count0(row$n_targeted_claim_ready_signature_overlaps)
  targeted_shared <- as_count0(row$n_targeted_neuropil_shared_overlaps %||% row$n_targeted_neuropil_shared_signature_overlaps)
  targeted_curated <- as_count0(row$n_targeted_curated_microglia_program_overlaps)
  curated_warning <- as.character(row$curated_program_overlap_warning %||% NA_character_)
  unique_curated_proteins <- as_count0(row$n_unique_curated_program_overlap_proteins)
  driver_is_caution <- grepl("caution|single-protein|mitochondrial", targeted_driver, ignore.case = TRUE) ||
    (!is.na(curated_warning) && nzchar(curated_warning))
  neuro <- max_frac(frac(c("empirical_neuropil_sensitive_high_confidence_fraction", "empirical_neuropil_enriched_fraction")), frac(c("canonical_neuronal_synaptic_neuropil_fraction", "neuropil_synaptic_neuronal_marker_fraction")))
  astro <- max_frac(frac(c("astrocyte_endfoot_gliovascular_fraction", "astrocyte_endfoot_gliovascular_marker_fraction")), frac(c("canonical_astrocyte_fraction", "astrocyte_marker_fraction")))
  oligo <- max_frac(frac(c("canonical_oligodendrocyte_myelin_fraction", "oligodendrocyte_myelin_marker_fraction")), frac(c("canonical_opc_fraction")))
  mito <- frac(c("canonical_mitochondrial_oxphos_fraction", "mitochondrial_oxphos_marker_fraction"))
  ribo <- frac(c("canonical_ribosomal_translation_fraction", "ribosomal_translation_marker_fraction"))
  rnp <- frac(c("canonical_rnp_rna_processing_fraction", "rnp_rna_processing_marker_fraction"))

  vascular_ecm <- max_frac(bm, bbb, mural, adhesion)
  if (is.finite(vascular_ecm) && vascular_ecm >= classification_threshold && max_frac(bm, adhesion) >= max_frac(bbb, mural, 0)) return("perivascular basement membrane / ECM")
  if (is.finite(max_frac(bbb, mural)) && max_frac(bbb, mural) >= classification_threshold) return("BBB / vascular mural")
  if (!is.na(targeted_driver) && nzchar(targeted_driver) && !driver_is_caution) return(targeted_driver)
  if (targeted_shared > 0 && (targeted_micro > 0 || targeted_curated > 0)) return("shared microglia-neuropil signature overlap")
  if (targeted_claim > 0 || targeted_micro > 1) return("microglia-enriched signature overlap")
  if (!is.na(curated_warning) && grepl("non_specific_mitochondrial_overlap", curated_warning)) return("curated mitochondrial/oxidative-stress overlap")
  if (!is.na(curated_warning) && grepl("single_protein_curated_overlap", curated_warning)) return("curated microglia-relevant program overlap (single protein)")
  if (targeted_curated > 0 && unique_curated_proteins <= 1L) return("curated microglia-relevant program overlap (single protein)")
  if (targeted_curated > 0) return("curated microglia-relevant program overlap")
  if (targeted_shared > 0) return("shared microglia-neuropil signature overlap")
  if (is.finite(micro_state) && micro_state >= classification_threshold) return("microglia phagolysosomal / activation")
  if (is.finite(micro) && micro >= classification_threshold && !is.finite(neuro)) return("microglia-enriched ROI")
  if (is.finite(micro) && micro >= classification_threshold && is.finite(neuro) && neuro >= classification_threshold) return("shared microglia-neuropil microenvironment")
  if (is.finite(neuro) && neuro >= classification_threshold) return("neuronal synaptic / neuropil")
  if (is.finite(astro) && astro >= classification_threshold) return("astrocyte / endfoot")
  if (is.finite(oligo) && oligo >= classification_threshold) return("oligodendrocyte / myelin")
  if (is.finite(mito) && mito >= classification_threshold) return("mitochondrial / OXPHOS")
  if (is.finite(ribo) && ribo >= classification_threshold) return("ribosomal translation")
  if (is.finite(rnp) && rnp >= classification_threshold) return("RNA processing / RNP")
  "mixed / low-specificity"
}

confidence_from_evidence <- function(row) {
  vals <- suppressWarnings(as.numeric(unlist(row[grep("(_marker_fraction|_fraction)$", names(row), value = TRUE)])))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return("low")
  m <- max(vals)
  if (m >= 0.30) "high" else if (m >= classification_threshold) "moderate" else "low"
}

targeted_signature_caution_confidence <- function(row, confidence) {
  warning <- as.character(row$curated_program_overlap_warning %||% NA_character_)
  driver <- as.character(row$targeted_signature_primary_driver %||% NA_character_)
  targeted_micro <- as_count0(row$n_targeted_microglia_signature_overlaps)
  targeted_claim <- as_count0(row$n_targeted_claim_ready_signature_overlaps)
  if ((!is.na(warning) && nzchar(warning)) &&
      targeted_claim == 0L && targeted_micro == 0L &&
      !grepl("shared microglia-neuropil", driver, ignore.case = TRUE)) {
    return("low")
  }
  confidence
}


marker_stats_for <- function(genes) {
  genes_key <- normalize_gene_token(genes)
  out <- lapply(panel_names, function(panel) {
    markers <- normalize_gene_token(marker_sets[[panel]] %||% character())
    hits <- genes[genes_key %in% markers]
    data.frame(
      panel = panel,
      fraction = if (length(genes)) length(unique(normalize_gene_token(hits))) / length(unique(genes_key)) else NA_real_,
      hits = paste(unique(hits), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }) |> dplyr::bind_rows()
  frac <- stats::setNames(out$fraction, paste0(out$panel, "_fraction"))
  hits <- stats::setNames(out$hits, paste0(out$panel, "_hits"))
  legacy_frac <- stats::setNames(out$fraction, paste0(out$panel, "_marker_fraction"))
  legacy_hits <- stats::setNames(out$hits, paste0(out$panel, "_marker_hits"))
  as.list(c(frac, hits, legacy_frac, legacy_hits))
}

top_go_for_module <- function(module_color) {
  wgcna_go_terms_for_module(go, module_color)
}

classify_module <- function(row) {
  frac <- function(candidates) {
    hit <- candidates[candidates %in% names(row)][1]
    if (length(hit) && !is.na(hit)) as.numeric(row[[hit]]) else NA_real_
  }
  max_frac <- function(...) {
    vals <- c(...)
    vals <- vals[is.finite(vals)]
    if (length(vals)) max(vals) else NA_real_
  }

  micro <- max_frac(frac(c("empirical_microglia_roi_high_confidence_fraction", "empirical_microglia_roi_enriched_fraction")), frac(c("canonical_microglia_homeostatic_fraction", "microglia_marker_fraction")))
  micro_state <- frac(c("canonical_microglia_phagolysosomal_state_fraction"))
  neuro <- max_frac(frac(c("empirical_neuropil_sensitive_high_confidence_fraction", "empirical_neuropil_enriched_fraction")), frac(c("canonical_neuronal_synaptic_neuropil_fraction", "neuropil_synaptic_neuronal_marker_fraction")))
  shared <- frac(c("empirical_microglia_neuropil_shared_fraction"))

  astro <- max_frac(frac(c("astrocyte_endfoot_gliovascular_fraction", "astrocyte_endfoot_gliovascular_marker_fraction")), frac(c("canonical_astrocyte_fraction", "astrocyte_marker_fraction")))
  oligo <- max_frac(frac(c("canonical_oligodendrocyte_myelin_fraction", "oligodendrocyte_myelin_marker_fraction")), frac(c("canonical_opc_fraction")))
  vascular <- max_frac(frac(c("canonical_endothelial_vascular_fraction", "endothelial_pericyte_vascular_marker_fraction")), frac(c("canonical_pericyte_vascular_fraction")))
  bm_ecm <- max_frac(frac(c("basement_membrane_perivascular_ecm_fraction", "basement_membrane_perivascular_ecm_marker_fraction")), frac(c("integrin_ecm_adhesion_fraction", "integrin_ecm_adhesion_marker_fraction")))
  bbb_mural <- max_frac(frac(c("bbb_endothelial_transport_fraction", "bbb_endothelial_transport_marker_fraction")), frac(c("pericyte_mural_ng2_fraction", "pericyte_mural_ng2_marker_fraction")))
  peripheral <- frac(c("canonical_peripheral_myeloid_caution_fraction"))

  robust <- as_count0(row$n_microglia_robust_term_overlaps)
  sensitive <- as_count0(row$n_neuropil_sensitive_term_overlaps)
  mixed <- as_count0(row$n_mixed_microenvironment_term_overlaps)
  targeted_micro <- as_count0(row$n_targeted_microglia_signature_overlaps)
  targeted_claim <- as_count0(row$n_targeted_claim_ready_signature_overlaps)
  targeted_shared <- as_count0(row$n_targeted_neuropil_shared_signature_overlaps)
  unique_curated_proteins <- as_count0(row$n_unique_curated_program_overlap_proteins)
  curated_warning <- as.character(row$curated_program_overlap_warning %||% NA_character_)

  microglia_evidence <- (is.finite(micro) && micro >= classification_threshold) || robust > 0 || targeted_micro > 0 || targeted_claim > 0
  microglia_state_evidence <- is.finite(micro_state) && micro_state >= classification_threshold
  neuropil_evidence <- (is.finite(neuro) && neuro >= classification_threshold) || sensitive > 0 || targeted_shared > 0
  shared_evidence <- (is.finite(shared) && shared >= classification_threshold) || mixed > 0 || (targeted_micro > 0 && targeted_shared > 0)
  vascular_ecm_evidence <- is.finite(max_frac(vascular, bm_ecm, bbb_mural)) && max_frac(vascular, bm_ecm, bbb_mural) >= classification_threshold
  astro_evidence <- is.finite(astro) && astro >= classification_threshold
  oligo_evidence <- is.finite(oligo) && oligo >= classification_threshold
  peripheral_evidence <- is.finite(peripheral) && peripheral >= classification_threshold

  if ((microglia_evidence && neuropil_evidence) || shared_evidence) return("shared_microenvironment")
  if (is.finite(bm_ecm) && bm_ecm >= classification_threshold) return("vascular_basement_membrane_ecm")
  if (vascular_ecm_evidence) return("vascular_bbb_mural")
  if (microglia_state_evidence && !neuropil_evidence) return("microglia_state_or_activation_supported")
  if (microglia_evidence && !neuropil_evidence && !vascular_ecm_evidence && !astro_evidence && !oligo_evidence) return("microglia_supported")
  if (neuropil_evidence && !microglia_evidence) return("neuropil_sensitive")
  if (astro_evidence) return("astrocyte_or_endfoot_sensitive")
  if (oligo_evidence) return("oligodendrocyte_or_myelin_sensitive")
  if (peripheral_evidence && !microglia_evidence) return("peripheral_myeloid_caution")
  "ambiguous_or_mixed"
}

classify_module_at_threshold <- function(row, threshold) {
  old_threshold <- classification_threshold
  assign("classification_threshold", threshold, envir = environment(classify_module))
  on.exit(assign("classification_threshold", old_threshold, envir = environment(classify_module)), add = TRUE)
  classify_module(row)
}

supporting_marker_panels <- function(row, threshold = classification_threshold) {
  vals <- semantic_marker_panel_values(row)
  vals <- vals[is.finite(vals) & vals >= threshold]
  if (!length(vals)) return(NA_character_)
  paste(names(sort(vals, decreasing = TRUE)), collapse = ";")
}

primary_marker_fraction <- function(row) {
  vals <- semantic_marker_panel_values(row)
  vals <- vals[is.finite(vals)]
  if (length(vals)) max(vals) else NA_real_
}

annotation_basis_for <- function(row) {
  marker_frac <- primary_marker_fraction(row)
  targeted <- as_count0(row$n_targeted_microglia_signature_overlaps) +
    as_count0(row$n_targeted_claim_ready_signature_overlaps) +
    as_count0(row$n_targeted_curated_microglia_program_overlaps)
  go_hits <- as_count0(row$n_microglia_robust_term_overlaps) +
    as_count0(row$n_neuropil_sensitive_term_overlaps) +
    as_count0(row$n_mixed_microenvironment_term_overlaps)
  dplyr::case_when(
    is.finite(marker_frac) & marker_frac >= classification_threshold & targeted > 0 ~ "mixed",
    targeted > 0 ~ "targeted_signature",
    is.finite(marker_frac) & marker_frac >= classification_threshold ~ "marker_panel",
    go_hits > 0 ~ "GO",
    TRUE ~ "insufficient"
  )
}

annotation_confidence_for <- function(confidence, stable, basis) {
  confidence <- as.character(confidence)
  basis <- as.character(basis)
  dplyr::case_when(
    !isTRUE(stable) ~ "low",
    basis == "insufficient" ~ "unresolved",
    confidence == "high" ~ "high",
    confidence %in% c("moderate", "medium") ~ "moderate",
    TRUE ~ "low"
  )
}

annotation_downgrade_for <- function(class_at_0.05, class_at_0.10, class_at_0.20, basis, confidence) {
  reasons <- character()
  if (length(unique(c(class_at_0.05, class_at_0.10, class_at_0.20))) > 1L) reasons <- c(reasons, "threshold_unstable")
  if (identical(as.character(basis), "insufficient")) reasons <- c(reasons, "insufficient_annotation_evidence")
  if (as.character(confidence) %in% c("low", "unresolved")) reasons <- c(reasons, "low_annotation_confidence")
  if (length(reasons)) paste(unique(reasons), collapse = ";") else "not_applicable"
}

classify_rationale <- function(row) {
  frac <- function(candidates) {
    hit <- candidates[candidates %in% names(row)][1]
    if (length(hit) && !is.na(hit)) as.numeric(row[[hit]]) else NA_real_
  }
  max_frac <- function(...) {
    vals <- c(...)
    vals <- vals[is.finite(vals)]
    if (length(vals)) max(vals) else NA_real_
  }
  canonical_micro <- frac(c("canonical_microglia_homeostatic_fraction", "microglia_marker_fraction"))
  empirical_micro <- max_frac(frac(c("empirical_microglia_roi_high_confidence_fraction")), frac(c("empirical_microglia_roi_enriched_fraction")))
  canonical_neuro <- frac(c("canonical_neuronal_synaptic_neuropil_fraction", "neuropil_synaptic_neuronal_marker_fraction"))
  empirical_neuro <- max_frac(frac(c("empirical_neuropil_sensitive_high_confidence_fraction")), frac(c("empirical_neuropil_enriched_fraction")))
  shared <- frac(c("empirical_microglia_neuropil_shared_fraction"))
  micro_state <- frac(c("canonical_microglia_phagolysosomal_state_fraction"))
  other_cellular <- max_frac(
    frac(c("canonical_astrocyte_fraction", "astrocyte_marker_fraction")),
    frac(c("canonical_oligodendrocyte_myelin_fraction", "oligodendrocyte_myelin_marker_fraction")),
    frac(c("canonical_opc_fraction")),
    frac(c("canonical_endothelial_vascular_fraction", "endothelial_pericyte_vascular_marker_fraction")),
    frac(c("canonical_pericyte_vascular_fraction"))
  )
  peripheral <- frac(c("canonical_peripheral_myeloid_caution_fraction"))
  robust <- as_count0(row$n_microglia_robust_term_overlaps)
  sensitive <- as_count0(row$n_neuropil_sensitive_term_overlaps)
  mixed <- as_count0(row$n_mixed_microenvironment_term_overlaps)
  targeted_micro <- as_count0(row$n_targeted_microglia_signature_overlaps)
  targeted_claim <- as_count0(row$n_targeted_claim_ready_signature_overlaps)
  targeted_shared <- as_count0(row$n_targeted_neuropil_shared_signature_overlaps)
  unique_curated_proteins <- as_count0(row$n_unique_curated_program_overlap_proteins)
  curated_warning <- as.character(row$curated_program_overlap_warning %||% NA_character_)
  canonical_microglia_evidence <- (is.finite(canonical_micro) && canonical_micro >= classification_threshold) || robust > 0
  empirical_microglia_roi_evidence <- (is.finite(empirical_micro) && empirical_micro >= classification_threshold) || targeted_micro > 0 || targeted_claim > 0
  canonical_neuropil_evidence <- (is.finite(canonical_neuro) && canonical_neuro >= classification_threshold) || sensitive > 0 || targeted_shared > 0
  empirical_neuropil_evidence <- is.finite(empirical_neuro) && empirical_neuro >= classification_threshold
  empirical_shared_microenvironment_evidence <- (is.finite(shared) && shared >= classification_threshold) || mixed > 0 || (targeted_micro > 0 && targeted_shared > 0)
  other_cellular_or_vascular_evidence <- is.finite(other_cellular) && other_cellular >= classification_threshold
  microglia_state_or_activation_evidence <- is.finite(micro_state) && micro_state >= classification_threshold
  peripheral_myeloid_caution_evidence <- is.finite(peripheral) && peripheral >= classification_threshold
  list(
    canonical_microglia_evidence = canonical_microglia_evidence,
    empirical_microglia_roi_evidence = empirical_microglia_roi_evidence,
    canonical_neuropil_evidence = canonical_neuropil_evidence,
    empirical_neuropil_evidence = empirical_neuropil_evidence,
    empirical_shared_microenvironment_evidence = empirical_shared_microenvironment_evidence,
    other_cellular_or_vascular_evidence = other_cellular_or_vascular_evidence,
    microglia_state_or_activation_evidence = microglia_state_or_activation_evidence,
    peripheral_myeloid_caution_evidence = peripheral_myeloid_caution_evidence,
    microglia_evidence = canonical_microglia_evidence || empirical_microglia_roi_evidence,
    neuropil_evidence = canonical_neuropil_evidence || empirical_neuropil_evidence,
    other_cellular_evidence = other_cellular_or_vascular_evidence,
    classification_rationale = paste0(
      "threshold=", classification_threshold,
      "; canonical_microglia=", signif(canonical_micro, 3), "; empirical_microglia=", signif(empirical_micro, 3), "; microglia_terms=", robust,
      "; targeted_microglia_signatures=", targeted_micro, "; targeted_claim_ready_signatures=", targeted_claim,
      "; unique_curated_program_overlap_proteins=", unique_curated_proteins,
      "; curated_program_overlap_warning=", dplyr::coalesce(curated_warning, "none"),
      "; canonical_neuropil=", signif(canonical_neuro, 3), "; empirical_neuropil=", signif(empirical_neuro, 3), "; neuropil_terms=", sensitive,
      "; shared=", signif(shared, 3), "; mixed_terms=", mixed, "; targeted_neuropil_shared_signatures=", targeted_shared,
      "; microglia_state=", signif(micro_state, 3),
      "; other_cellular_or_vascular=", signif(other_cellular, 3),
      "; peripheral_myeloid_caution=", signif(peripheral, 3)
    )
  )
}

module_neuropil_reference_counts <- function(module_color) {
  empty <- data.frame(
    n_microglia_robust_term_overlaps = 0L,
    n_neuropil_sensitive_term_overlaps = 0L,
    n_mixed_microenvironment_term_overlaps = 0L,
    n_ambiguous_term_overlaps = 0L,
    best_overlapping_microglia_terms = NA_character_,
    best_overlapping_neuropil_terms = NA_character_,
    stringsAsFactors = FALSE
  )
  if (is.null(go) || !nrow(go) || is.null(neuropil_ref) || !nrow(neuropil_ref)) return(empty)
  color_col <- first_present_col(go, c("ModuleColor", "module", "ModuleID"))
  desc_col <- first_present_col(go, c("Description", "description", "term_description", "ID"))
  ref_desc_col <- first_present_col(neuropil_ref, c("term_description", "Description", "description", "term_id", "ID"))
  class_col <- first_present_col(neuropil_ref, c("interpretation_class", "microenvironment_class"))
  if (is.na(color_col) || is.na(desc_col) || is.na(ref_desc_col) || is.na(class_col)) return(empty)
  module_terms <- go |>
    dplyr::filter(as.character(.data[[color_col]]) %in% c(as.character(module_color), paste0("ME", module_color))) |>
    dplyr::mutate(term_key = toupper(trimws(as.character(.data[[desc_col]])))) |>
    dplyr::filter(nzchar(.data$term_key))
  ref_terms <- neuropil_ref |>
    dplyr::mutate(term_key = toupper(trimws(as.character(.data[[ref_desc_col]])))) |>
    dplyr::filter(nzchar(.data$term_key))
  hit <- dplyr::inner_join(module_terms, ref_terms, by = "term_key", relationship = "many-to-many") |>
    dplyr::mutate(.ref_class = .data[[class_col]]) |>
    dplyr::distinct(term_key, .ref_class, .keep_all = TRUE)
  if (!nrow(hit)) return(empty)
  cls <- as.character(hit[[class_col]])
  data.frame(
    n_microglia_robust_term_overlaps = length(unique(hit$term_key[cls == "microglia_robust"])),
    n_neuropil_sensitive_term_overlaps = length(unique(hit$term_key[cls %in% c("neuropil_sensitive", "neuropil_marker_enriched")])),
    n_mixed_microenvironment_term_overlaps = length(unique(hit$term_key[cls == "mixed_microenvironment"])),
    n_ambiguous_term_overlaps = length(unique(hit$term_key[cls == "ambiguous"])),
    best_overlapping_microglia_terms = paste(utils::head(unique(hit$term_key[cls == "microglia_robust"]), 5), collapse = ";"),
    best_overlapping_neuropil_terms = paste(utils::head(unique(hit$term_key[cls %in% c("neuropil_sensitive", "neuropil_marker_enriched")]), 5), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

signature_tokens <- function(...) {
  fields <- c(...)
  fields <- fields[!is.na(fields) & nzchar(fields)]
  if (!length(fields)) return(character())
  normalize_gene_token(unlist(strsplit(paste(fields, collapse = ";"), "[/;,|[:space:]]+"), use.names = FALSE))
}

is_mito_oxidative_signature <- function(signature) {
  grepl("mitochond|oxidative|oxphos|respirat|electron_transport", as.character(signature), ignore.case = TRUE)
}

generic_mitochondrial_tokens <- normalize_gene_token(c(
  "VDAC1", "VDAC2", "VDAC3", "Q60932",
  "ATP5F1A", "ATP5F1B", "ATP5A1", "ATP5B",
  "NDUFA9", "NDUFS1", "NDUFV1", "UQCRC1", "UQCRC2",
  "COX4I1", "CYCS", "SOD1", "SOD2", "PRDX1", "PRDX3",
  "GPX1", "TXN", "TXNRD1", "PARK7"
))

curated_program_overlap_warning_for <- function(curated_rows) {
  if (is.null(curated_rows) || !nrow(curated_rows)) return(NA_character_)
  proteins <- unique(signature_tokens(curated_rows$overlap_proteins))
  flags <- character()
  if (length(proteins) == 1L) flags <- c(flags, "single_protein_curated_overlap")
  if (any(is_mito_oxidative_signature(curated_rows$signature)) &&
      length(intersect(proteins, generic_mitochondrial_tokens)) > 0L) {
    flags <- c(flags, "non_specific_mitochondrial_overlap")
  }
  if (!length(flags)) NA_character_ else paste(unique(flags), collapse = ";")
}

targeted_signature_interpretation_note <- function(curated_warning, shared_n, microglia_enriched_n, claim_ready_n) {
  if (!is.na(curated_warning) && nzchar(curated_warning)) {
    return(paste(
      "Curated microglia signatures are microglia-relevant gene sets, not necessarily microglia-specific markers.",
      "Single ubiquitous proteins such as VDAC1/Q60932 support mitochondrial/oxidative-stress interpretation, not microglia specificity."
    ))
  }
  if (shared_n > 0L) {
    return("Targeted signature evidence overlaps neuropil-shared signal; interpret as local microenvironment/shared biology unless independent microglia-enriched evidence is present.")
  }
  if (claim_ready_n > 0L || microglia_enriched_n > 0L) {
    return("Microglia-targeted evidence is driven by empirical/reference-supported signatures; claim-ready status remains restricted to these classes.")
  }
  NA_character_
}

targeted_signature_detail_columns <- c(
  "dataset", "ModuleID", "ModuleColor", "signature", "signature_source",
  "microglia_signature_class", "contrast_class", "comparison", "NES", "padj",
  "neuropil_reference_NES", "neuropil_reference_padj", "claim_ready",
  "overlap_proteins", "n_overlap_proteins", "leading_edge", "matched_genes"
)

empty_targeted_signature_details <- function() {
  out <- data.frame(
    dataset = character(),
    ModuleID = character(),
    ModuleColor = character(),
    signature = character(),
    signature_source = character(),
    microglia_signature_class = character(),
    contrast_class = character(),
    comparison = character(),
    NES = numeric(),
    padj = numeric(),
    neuropil_reference_NES = numeric(),
    neuropil_reference_padj = numeric(),
    claim_ready = logical(),
    overlap_proteins = character(),
    n_overlap_proteins = integer(),
    leading_edge = character(),
    matched_genes = character(),
    stringsAsFactors = FALSE
  )
  out[targeted_signature_detail_columns]
}

empty_targeted_signature_summary <- function() {
  data.frame(
    n_targeted_microglia_signature_overlaps = 0L,
    n_targeted_claim_ready_signature_overlaps = 0L,
    n_targeted_neuropil_shared_signature_overlaps = 0L,
    n_targeted_signature_overlap_proteins = 0L,
    best_targeted_microglia_signatures = NA_character_,
    n_targeted_microglia_enriched_empirical_overlaps = 0L,
    n_targeted_microglia_enriched_reference_supported_overlaps = 0L,
    n_targeted_curated_microglia_program_overlaps = 0L,
    n_targeted_neuropil_shared_overlaps = 0L,
    n_targeted_ambiguous_overlaps = 0L,
    n_unique_targeted_signatures = 0L,
    n_unique_targeted_overlap_proteins = 0L,
    n_unique_curated_program_signatures = 0L,
    n_unique_curated_program_overlap_proteins = 0L,
    curated_program_overlap_proteins = NA_character_,
    curated_program_overlap_warning = NA_character_,
    targeted_signature_overlap_interpretation_note = NA_character_,
    targeted_signature_driver_evidence_tier = NA_character_,
    targeted_signature_primary_driver = NA_character_,
    targeted_signature_driver_class = NA_character_,
    targeted_signature_driver_signature = NA_character_,
    targeted_signature_driver_padj = NA_real_,
    targeted_signature_driver_NES = NA_real_,
    targeted_signature_driver_overlap_proteins = NA_character_,
    best_targeted_microglia_enriched_signatures = NA_character_,
    best_targeted_curated_microglia_programs = NA_character_,
    best_targeted_neuropil_shared_signatures = NA_character_,
    best_targeted_signature_overlap_proteins = NA_character_,
    stringsAsFactors = FALSE
  )
}

targeted_signature_claim_ready <- function(padj, contrast_class, microglia_signature_class) {
  is.finite(padj) & padj < 0.05 &
    contrast_class %in% c("within_region_condition", "cross_region_same_condition") &
    microglia_signature_class %in% c("microglia_enriched_empirical", "microglia_enriched_reference_supported")
}

targeted_signature_driver_label <- function(empirical_n, reference_n, curated_n, shared_n, ambiguous_n) {
  microglia_enriched_n <- empirical_n + reference_n
  if (microglia_enriched_n > 0 && microglia_enriched_n >= max(curated_n, shared_n, ambiguous_n)) {
    return("microglia-enriched signature overlap")
  }
  if (shared_n > 0 && (microglia_enriched_n > 0 || curated_n > 0 || shared_n >= max(microglia_enriched_n, curated_n, ambiguous_n))) {
    return("shared microglia-neuropil signature overlap (caution)")
  }
  if (curated_n > 0 && curated_n >= max(microglia_enriched_n, shared_n, ambiguous_n)) {
    return("curated microglia-relevant program overlap (caution)")
  }
  if (shared_n > 0) return("shared microglia-neuropil signature overlap (caution)")
  if (curated_n > 0) return("curated microglia-relevant program overlap (caution)")
  if (microglia_enriched_n > 0) return("microglia-enriched signature overlap")
  if (ambiguous_n > 0) return("ambiguous targeted signature overlap")
  NA_character_
}

build_module_targeted_signature_details <- function(module_rows_tbl) {
  if (is.null(targeted_signature_ref) || !nrow(targeted_signature_ref)) return(empty_targeted_signature_details())
  needed <- c("signature", "microglia_signature_class", "contrast_class", "matched_genes", "leading_edge", "padj")
  if (!all(needed %in% names(targeted_signature_ref))) return(empty_targeted_signature_details())

  sig <- targeted_signature_ref
  optional_cols <- c("signature_source", "comparison", "NES", "neuropil_reference_NES", "neuropil_reference_padj")
  for (nm in optional_cols) if (!nm %in% names(sig)) sig[[nm]] <- NA
  sig <- sig |>
    dplyr::mutate(
      padj = suppressWarnings(as.numeric(.data$padj)),
      NES = suppressWarnings(as.numeric(.data$NES)),
      neuropil_reference_NES = suppressWarnings(as.numeric(.data$neuropil_reference_NES)),
      neuropil_reference_padj = suppressWarnings(as.numeric(.data$neuropil_reference_padj)),
      significant = !is.na(.data$padj) & .data$padj < 0.05,
      claim_ready = targeted_signature_claim_ready(.data$padj, .data$contrast_class, .data$microglia_signature_class)
    ) |>
    dplyr::filter(.data$significant)
  if (!nrow(sig)) return(empty_targeted_signature_details())

  detail_rows <- lapply(seq_len(nrow(module_rows_tbl)), function(m) {
    genes <- normalize_gene_token(module_rows_tbl$proteins[[m]])
    genes <- unique(genes[nzchar(genes) & !is.na(genes)])
    if (!length(genes)) return(NULL)
    rows <- lapply(seq_len(nrow(sig)), function(i) {
      tokens <- signature_tokens(sig$leading_edge[[i]], sig$matched_genes[[i]])
      hits <- intersect(genes, tokens)
      if (!length(hits)) return(NULL)
      data.frame(
        dataset = DATASET,
        ModuleID = as.character(module_rows_tbl$ModuleID[[m]]),
        ModuleColor = as.character(module_rows_tbl$ModuleColor[[m]]),
        signature = as.character(sig$signature[[i]]),
        signature_source = as.character(sig$signature_source[[i]]),
        microglia_signature_class = as.character(sig$microglia_signature_class[[i]]),
        contrast_class = as.character(sig$contrast_class[[i]]),
        comparison = as.character(sig$comparison[[i]]),
        NES = sig$NES[[i]],
        padj = sig$padj[[i]],
        neuropil_reference_NES = sig$neuropil_reference_NES[[i]],
        neuropil_reference_padj = sig$neuropil_reference_padj[[i]],
        claim_ready = isTRUE(sig$claim_ready[[i]]),
        overlap_proteins = paste(hits, collapse = ";"),
        n_overlap_proteins = length(unique(hits)),
        leading_edge = as.character(sig$leading_edge[[i]]),
        matched_genes = as.character(sig$matched_genes[[i]]),
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(rows)
  })
  out <- dplyr::bind_rows(detail_rows)
  if (!nrow(out)) return(empty_targeted_signature_details())
  out[targeted_signature_detail_columns]
}

summarise_module_targeted_signature_details <- function(details, module_rows_tbl) {
  empty <- empty_targeted_signature_summary()
  if (is.null(details) || !nrow(details)) {
    return(dplyr::bind_cols(module_rows_tbl |> dplyr::select("ModuleID", "ModuleColor"), empty[rep(1, nrow(module_rows_tbl)), , drop = FALSE]))
  }
  detail_summary <- details |>
    dplyr::mutate(
      overlap_key = paste(.data$signature, .data$comparison, sep = "||"),
      class_group = dplyr::case_when(
        .data$microglia_signature_class == "microglia_enriched_empirical" ~ "microglia_enriched_empirical",
        .data$microglia_signature_class == "microglia_enriched_reference_supported" ~ "microglia_enriched_reference_supported",
        .data$microglia_signature_class == "curated_microglia_program" ~ "curated_microglia_program",
        .data$microglia_signature_class == "neuropil_shared" ~ "neuropil_shared",
        TRUE ~ "ambiguous"
      )
    ) |>
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) |>
    dplyr::group_modify(function(.x, .y) {
      empirical_n <- dplyr::n_distinct(.x$overlap_key[.x$class_group == "microglia_enriched_empirical"])
      reference_n <- dplyr::n_distinct(.x$overlap_key[.x$class_group == "microglia_enriched_reference_supported"])
      curated_n <- dplyr::n_distinct(.x$overlap_key[.x$class_group == "curated_microglia_program"])
      shared_n <- dplyr::n_distinct(.x$overlap_key[.x$class_group == "neuropil_shared"])
      ambiguous_n <- dplyr::n_distinct(.x$overlap_key[.x$class_group == "ambiguous"])
      microglia_enriched <- .x |> dplyr::filter(.data$class_group %in% c("microglia_enriched_empirical", "microglia_enriched_reference_supported"))
      curated <- .x |> dplyr::filter(.data$class_group == "curated_microglia_program")
      shared <- .x |> dplyr::filter(.data$class_group == "neuropil_shared")
      proteins <- signature_tokens(.x$overlap_proteins)
      microglia_enriched_proteins <- signature_tokens(microglia_enriched$overlap_proteins)
      curated_proteins <- signature_tokens(curated$overlap_proteins)
      curated_warning <- curated_program_overlap_warning_for(curated)
      claim_ready_n <- dplyr::n_distinct(.x$overlap_key[.x$claim_ready])
      claim_ready_proteins <- signature_tokens(.x$overlap_proteins[.x$claim_ready])
      microglia_enriched_n <- empirical_n + reference_n
      has_claim_ready_multi <- claim_ready_n > 0L && length(unique(claim_ready_proteins)) >= 2L
      has_microglia_enriched_multi <- microglia_enriched_n > 0L && length(unique(microglia_enriched_proteins)) >= 2L
      primary_driver <- targeted_signature_driver_label(empirical_n, reference_n, curated_n, shared_n, ambiguous_n)
      if (has_claim_ready_multi || has_microglia_enriched_multi) {
        primary_driver <- "microglia-enriched signature overlap"
      } else if (shared_n > 0L) {
        primary_driver <- "shared microglia-neuropil signature overlap (caution)"
      } else if (curated_n > 0L) {
        primary_driver <- dplyr::case_when(
          !is.na(curated_warning) & grepl("non_specific_mitochondrial_overlap", curated_warning) ~ "curated mitochondrial/oxidative-stress overlap (caution)",
          !is.na(curated_warning) & grepl("single_protein_curated_overlap", curated_warning) ~ "curated microglia-relevant program overlap (single-protein caution)",
          TRUE ~ "curated microglia-relevant program overlap (caution)"
        )
      } else if (microglia_enriched_n > 0L) {
        primary_driver <- "microglia-enriched signature overlap (single-protein caution)"
      }
      driver_pool <- .x |>
        dplyr::mutate(
          driver_priority = dplyr::case_when(
            .data$claim_ready & .data$n_overlap_proteins >= 2L ~ 1L,
            .data$class_group %in% c("microglia_enriched_empirical", "microglia_enriched_reference_supported") & .data$n_overlap_proteins >= 2L ~ 2L,
            .data$class_group == "neuropil_shared" ~ 3L,
            .data$class_group %in% c("microglia_enriched_empirical", "microglia_enriched_reference_supported") ~ 4L,
            .data$class_group == "curated_microglia_program" & (is.na(curated_warning) | !nzchar(curated_warning)) ~ 5L,
            .data$class_group == "curated_microglia_program" ~ 6L,
            primary_driver == "ambiguous targeted signature overlap" & .data$class_group == "ambiguous" ~ 7L,
            TRUE ~ 8L
          )
        ) |>
        dplyr::arrange(.data$driver_priority, .data$padj, dplyr::desc(abs(.data$NES)), dplyr::desc(.data$n_overlap_proteins))
      driver <- driver_pool[1, , drop = FALSE]
      evidence_tier <- dplyr::case_when(
        has_claim_ready_multi ~ "claim_ready_empirical_or_reference_supported",
        has_microglia_enriched_multi ~ "empirical_or_reference_supported_multi_protein",
        shared_n > 0L ~ "caution_neuropil_shared",
        curated_n > 0L & !is.na(curated_warning) & nzchar(curated_warning) ~ "caution_curated_single_or_low_specificity",
        curated_n > 0L ~ "caution_curated_microglia_relevant_program",
        microglia_enriched_n > 0L ~ "caution_single_protein_empirical_or_reference_supported",
        TRUE ~ "ambiguous"
      )
      note <- targeted_signature_interpretation_note(curated_warning, shared_n, microglia_enriched_n, claim_ready_n)
      # Curated microglia programs are microglia-relevant gene sets, not purified
      # microglia-specific markers; single ubiquitous proteins such as VDAC1/Q60932
      # support mitochondrial/oxidative-stress interpretation, not microglia specificity.
      data.frame(
        n_targeted_microglia_signature_overlaps = empirical_n + reference_n,
        n_targeted_claim_ready_signature_overlaps = claim_ready_n,
        n_targeted_neuropil_shared_signature_overlaps = shared_n,
        n_targeted_signature_overlap_proteins = length(unique(proteins)),
        best_targeted_microglia_signatures = paste(utils::head(unique(microglia_enriched$signature[order(microglia_enriched$padj)]), 5), collapse = ";"),
        n_targeted_microglia_enriched_empirical_overlaps = empirical_n,
        n_targeted_microglia_enriched_reference_supported_overlaps = reference_n,
        n_targeted_curated_microglia_program_overlaps = curated_n,
        n_targeted_neuropil_shared_overlaps = shared_n,
        n_targeted_ambiguous_overlaps = ambiguous_n,
        n_unique_targeted_signatures = dplyr::n_distinct(.x$signature),
        n_unique_targeted_overlap_proteins = length(unique(proteins)),
        n_unique_curated_program_signatures = dplyr::n_distinct(curated$signature),
        n_unique_curated_program_overlap_proteins = length(unique(curated_proteins)),
        curated_program_overlap_proteins = paste(utils::head(unique(curated_proteins), 20), collapse = ";"),
        curated_program_overlap_warning = curated_warning,
        targeted_signature_overlap_interpretation_note = note,
        targeted_signature_driver_evidence_tier = evidence_tier,
        targeted_signature_primary_driver = primary_driver,
        targeted_signature_driver_class = as.character(driver$microglia_signature_class[[1]]),
        targeted_signature_driver_signature = as.character(driver$signature[[1]]),
        targeted_signature_driver_padj = driver$padj[[1]],
        targeted_signature_driver_NES = driver$NES[[1]],
        targeted_signature_driver_overlap_proteins = as.character(driver$overlap_proteins[[1]]),
        best_targeted_microglia_enriched_signatures = paste(utils::head(unique(microglia_enriched$signature[order(microglia_enriched$padj)]), 5), collapse = ";"),
        best_targeted_curated_microglia_programs = paste(utils::head(unique(curated$signature[order(curated$padj)]), 5), collapse = ";"),
        best_targeted_neuropil_shared_signatures = paste(utils::head(unique(shared$signature[order(shared$padj)]), 5), collapse = ";"),
        best_targeted_signature_overlap_proteins = paste(utils::head(unique(proteins), 20), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }) |>
    dplyr::ungroup()
  out <- module_rows_tbl |>
    dplyr::select("ModuleID", "ModuleColor") |>
    dplyr::left_join(detail_summary, by = c("ModuleID", "ModuleColor"))
  summary_cols <- names(empty)
  for (nm in summary_cols) {
    if (!nm %in% names(out)) out[[nm]] <- empty[[nm]][[1]]
    if (is.integer(empty[[nm]])) out[[nm]][is.na(out[[nm]])] <- 0L
    if (is.numeric(empty[[nm]])) out[[nm]][is.na(out[[nm]])] <- NA_real_
  }
  out[summary_cols] <- Map(function(x, template) {
    if (is.integer(template)) {
      x <- suppressWarnings(as.integer(x))
      x[is.na(x)] <- 0L
    } else if (is.numeric(template)) {
      x <- suppressWarnings(as.numeric(x))
    } else {
      x <- as.character(x)
      x[is.na(x) | !nzchar(x)] <- NA_character_
    }
    x
  }, out[summary_cols], empty[summary_cols])
  out
}

module_targeted_signature_counts <- function(module_genes) {
  empty <- empty_targeted_signature_summary()
  if (is.null(targeted_signature_ref) || !nrow(targeted_signature_ref)) return(empty)
  needed <- c("signature", "microglia_signature_class", "contrast_class", "matched_genes", "leading_edge", "padj")
  if (!all(needed %in% names(targeted_signature_ref))) return(empty)

  genes <- normalize_gene_token(module_genes)
  genes <- unique(genes[nzchar(genes) & !is.na(genes)])
  if (!length(genes)) return(empty)

  sig <- targeted_signature_ref |>
    dplyr::mutate(
      padj = suppressWarnings(as.numeric(.data$padj)),
      significant = !is.na(.data$padj) & .data$padj < 0.05,
      claim_ready = .data$significant &
        .data$contrast_class %in% c("within_region_condition", "cross_region_same_condition") &
        .data$microglia_signature_class %in% c("microglia_enriched_empirical", "microglia_enriched_reference_supported")
    ) |>
    dplyr::filter(.data$significant)
  if (!nrow(sig)) return(empty)

  overlap_rows <- lapply(seq_len(nrow(sig)), function(i) {
    tokens <- signature_tokens(sig$leading_edge[[i]], sig$matched_genes[[i]])
    hits <- intersect(genes, tokens)
    if (!length(hits)) return(NULL)
    data.frame(
      signature = as.character(sig$signature[[i]]),
      microglia_signature_class = as.character(sig$microglia_signature_class[[i]]),
      claim_ready = isTRUE(sig$claim_ready[[i]]),
      padj = sig$padj[[i]],
      overlap_proteins = paste(hits, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  hit <- dplyr::bind_rows(overlap_rows)
  if (!nrow(hit)) return(empty)

  robust <- hit |>
    dplyr::filter(.data$microglia_signature_class %in% c("microglia_enriched_empirical", "microglia_enriched_reference_supported"))
  shared <- hit |> dplyr::filter(.data$microglia_signature_class == "neuropil_shared")
  claim <- hit |> dplyr::filter(.data$claim_ready)
  proteins <- signature_tokens(hit$overlap_proteins)

  data.frame(
    n_targeted_microglia_signature_overlaps = length(unique(robust$signature)),
    n_targeted_claim_ready_signature_overlaps = length(unique(claim$signature)),
    n_targeted_neuropil_shared_signature_overlaps = length(unique(shared$signature)),
    n_targeted_signature_overlap_proteins = length(unique(proteins)),
    best_targeted_microglia_signatures = paste(utils::head(unique(robust$signature[order(robust$padj)]), 5), collapse = ";"),
    best_targeted_signature_overlap_proteins = paste(utils::head(unique(proteins), 20), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

module_rows <- definitions |>
  dplyr::mutate(
    ProteinToken = dplyr::coalesce(as.character(.data$GeneSymbol), as.character(.data$ProteinID), as.character(.data$UniProt))
  ) |>
  dplyr::group_by(.data$ModuleID, .data$ModuleColor) |>
  dplyr::summarise(
    n_proteins = dplyr::n(),
    proteins = list(unique(.data$ProteinToken)),
    top_hub_proteins = paste(utils::head(.data$ProteinToken[order(abs(as.numeric(.data$kME %||% .data$Weight)), decreasing = TRUE)], 25), collapse = ";"),
    n_top_hub_25 = min(25L, dplyr::n()),
    n_core_kME_0_6 = sum(abs(as.numeric(.data$kME %||% .data$Weight)) >= 0.6, na.rm = TRUE),
    module_eigengene = dplyr::first(as.character(.data$module_eigengene %||% paste0("ME", .data$ModuleColor))),
    module_label = dplyr::first(as.character(.data$ModuleLabel_Final %||% .data$ModuleID)),
    .groups = "drop"
  )

stats_list <- lapply(seq_len(nrow(module_rows)), function(i) {
  c(marker_stats_for(module_rows$proteins[[i]]), top_go_for_module(module_rows$ModuleColor[[i]]))
})
stats_df <- dplyr::bind_rows(lapply(stats_list, as.data.frame, stringsAsFactors = FALSE))

module_annot <- dplyr::bind_cols(module_rows |> dplyr::select(-"proteins"), stats_df) |>
  dplyr::mutate(dataset = DATASET, .before = "ModuleID")
for (nm in c("ModuleLabel_Final", "ModuleLabel_GO_BP", "best_GO_BP", "best_GO_MF", "best_GO_CC")) {
  if (!nm %in% names(module_annot)) module_annot[[nm]] <- NA_character_
}
fraction_cols <- grep("(_marker_fraction|_fraction)$", names(module_annot), value = TRUE)
module_annot[fraction_cols] <- lapply(module_annot[fraction_cols], function(x) suppressWarnings(as.numeric(x)))

if (DATASET == "microglia" || force_microglia) {
  ref_counts <- dplyr::bind_rows(lapply(module_annot$ModuleColor, module_neuropil_reference_counts))
  targeted_signature_details <- build_module_targeted_signature_details(module_rows)
  targeted_counts <- summarise_module_targeted_signature_details(targeted_signature_details, module_rows) |>
    dplyr::select(-"ModuleID", -"ModuleColor")
  module_annot <- module_annot |>
    dplyr::bind_cols(ref_counts, targeted_counts)
} else {
  targeted_signature_details <- empty_targeted_signature_details()
  module_annot <- module_annot |>
    dplyr::mutate(
      n_microglia_robust_term_overlaps = NA_integer_,
      n_neuropil_sensitive_term_overlaps = NA_integer_,
      n_mixed_microenvironment_term_overlaps = NA_integer_,
      n_ambiguous_term_overlaps = NA_integer_,
      best_overlapping_microglia_terms = NA_character_,
      best_overlapping_neuropil_terms = NA_character_,
      n_targeted_microglia_signature_overlaps = NA_integer_,
      n_targeted_claim_ready_signature_overlaps = NA_integer_,
      n_targeted_neuropil_shared_signature_overlaps = NA_integer_,
      n_targeted_signature_overlap_proteins = NA_integer_,
      best_targeted_microglia_signatures = NA_character_,
      n_targeted_microglia_enriched_empirical_overlaps = NA_integer_,
      n_targeted_microglia_enriched_reference_supported_overlaps = NA_integer_,
      n_targeted_curated_microglia_program_overlaps = NA_integer_,
      n_targeted_neuropil_shared_overlaps = NA_integer_,
      n_targeted_ambiguous_overlaps = NA_integer_,
      n_unique_targeted_signatures = NA_integer_,
      n_unique_targeted_overlap_proteins = NA_integer_,
      n_unique_curated_program_signatures = NA_integer_,
      n_unique_curated_program_overlap_proteins = NA_integer_,
      curated_program_overlap_proteins = NA_character_,
      curated_program_overlap_warning = NA_character_,
      targeted_signature_overlap_interpretation_note = NA_character_,
      targeted_signature_driver_evidence_tier = NA_character_,
      targeted_signature_primary_driver = NA_character_,
      targeted_signature_driver_class = NA_character_,
      targeted_signature_driver_signature = NA_character_,
      targeted_signature_driver_padj = NA_real_,
      targeted_signature_driver_NES = NA_real_,
      targeted_signature_driver_overlap_proteins = NA_character_,
      best_targeted_microglia_enriched_signatures = NA_character_,
      best_targeted_curated_microglia_programs = NA_character_,
      best_targeted_neuropil_shared_signatures = NA_character_,
      best_targeted_signature_overlap_proteins = NA_character_
    )
}

evidence_df <- dplyr::bind_rows(lapply(seq_len(nrow(module_annot)), function(i) as.data.frame(classify_rationale(module_annot[i, , drop = FALSE]), stringsAsFactors = FALSE)))
module_annot <- dplyr::bind_cols(module_annot, evidence_df)
module_annot$classification_threshold <- classification_threshold
module_annot$microenvironment_class <- vapply(seq_len(nrow(module_annot)), function(i) classify_module(module_annot[i, , drop = FALSE]), character(1))
module_annot$microenvironment_label <- vapply(seq_len(nrow(module_annot)), function(i) label_from_evidence(module_annot[i, , drop = FALSE]), character(1))
module_annot$microenvironment_confidence <- vapply(seq_len(nrow(module_annot)), function(i) confidence_from_evidence(module_annot[i, , drop = FALSE]), character(1))
module_annot$microenvironment_confidence <- vapply(seq_len(nrow(module_annot)), function(i) targeted_signature_caution_confidence(module_annot[i, , drop = FALSE], module_annot$microenvironment_confidence[[i]]), character(1))
module_annot$module_display_label <- paste0(module_annot$ModuleID, " — ", module_annot$microenvironment_label)
module_threshold_classes <- dplyr::bind_cols(
  module_annot |> dplyr::select("dataset", module_or_supermodule_id = "ModuleID"),
  dplyr::bind_cols(lapply(marker_threshold_sensitivity_values, function(thr) {
    out <- vapply(seq_len(nrow(module_annot)), function(i) classify_module_at_threshold(module_annot[i, , drop = FALSE], thr), character(1))
    tibble::tibble(!!paste0("class_at_", sprintf("%.2f", thr)) := out)
  }))
)
module_annot$annotation_basis <- vapply(seq_len(nrow(module_annot)), function(i) annotation_basis_for(module_annot[i, , drop = FALSE]), character(1))
module_annot$marker_fraction_primary <- vapply(seq_len(nrow(module_annot)), function(i) primary_marker_fraction(module_annot[i, , drop = FALSE]), numeric(1))
module_annot$marker_panels_supporting <- vapply(seq_len(nrow(module_annot)), function(i) supporting_marker_panels(module_annot[i, , drop = FALSE], classification_threshold), character(1))
module_annot$annotation_stable_across_thresholds <- module_threshold_classes$`class_at_0.05` == module_threshold_classes$`class_at_0.10` &
  module_threshold_classes$`class_at_0.10` == module_threshold_classes$`class_at_0.20`
module_annot$annotation_confidence <- vapply(seq_len(nrow(module_annot)), function(i) annotation_confidence_for(module_annot$microenvironment_confidence[[i]], module_annot$annotation_stable_across_thresholds[[i]], module_annot$annotation_basis[[i]]), character(1))
module_annot$annotation_downgrade_reason <- vapply(seq_len(nrow(module_annot)), function(i) annotation_downgrade_for(module_threshold_classes$`class_at_0.05`[[i]], module_threshold_classes$`class_at_0.10`[[i]], module_threshold_classes$`class_at_0.20`[[i]], module_annot$annotation_basis[[i]], module_annot$annotation_confidence[[i]]), character(1))
module_annot$unsafe_interpretation <- "Do not treat marker-panel microenvironment annotation as a purified cell-type or causal claim without claim-gate support."
module_annot$raw_GO_BP_terms <- dplyr::coalesce(as.character(module_annot$raw_GO_BP_terms), as.character(module_annot$top_GO_BP_labels))
module_annot$raw_GO_MF_terms <- dplyr::coalesce(as.character(module_annot$raw_GO_MF_terms), as.character(module_annot$top_GO_MF_labels))
module_annot$raw_GO_CC_terms <- dplyr::coalesce(as.character(module_annot$raw_GO_CC_terms), as.character(module_annot$top_GO_CC_labels))
module_annot$raw_top_GO_label <- dplyr::coalesce(
  as.character(module_annot$raw_top_GO_label),
  vapply(seq_len(nrow(module_annot)), function(i) {
    wgcna_collapse_terms(c(
      split_label_terms(module_annot$raw_GO_BP_terms[[i]]),
      split_label_terms(module_annot$raw_GO_MF_terms[[i]]),
      split_label_terms(module_annot$raw_GO_CC_terms[[i]])
    ), n = 1L)
  }, character(1))
)
module_annot$raw_module_label <- dplyr::coalesce(
  informative_member_label(module_annot$ModuleLabel_Final),
  informative_member_label(module_annot$module_label),
  informative_member_label(module_annot$raw_top_GO_label),
  as.character(module_annot$ModuleID)
)
module_annot$raw_hub_proteins <- as.character(module_annot$top_hub_proteins)
module_annot$raw_marker_or_signature_label <- as.character(module_annot$microenvironment_label)
semantic_df <- dplyr::bind_rows(lapply(seq_len(nrow(module_annot)), function(i) {
  row <- module_annot[i, , drop = FALSE]
  bp <- split_label_terms(row$raw_GO_BP_terms)
  mf <- split_label_terms(row$raw_GO_MF_terms)
  cc <- split_label_terms(row$raw_GO_CC_terms)
  as.data.frame(wgcna_clean_semantic_label(
    raw_label = dplyr::coalesce(row$raw_module_label, row$raw_top_GO_label),
    go_terms = c(bp, mf, cc),
    hubs = split_tokens(row$raw_hub_proteins),
    marker_label = row$raw_marker_or_signature_label,
    bp_terms = bp,
    mf_terms = mf,
    cc_terms = cc,
    marker_panel_fractions = semantic_marker_panel_values(row)
  ), stringsAsFactors = FALSE)
}))
module_annot <- dplyr::bind_cols(module_annot, semantic_df)
caution_df <- dplyr::bind_rows(lapply(seq_len(nrow(module_annot)), function(i) {
  row <- module_annot[i, , drop = FALSE]
  as.data.frame(wgcna_microenvironment_caution(
    cls = row$microenvironment_class,
    label = row$microenvironment_label,
    rationale = row$classification_rationale,
    dataset = DATASET
  ), stringsAsFactors = FALSE)
}))
module_annot <- dplyr::bind_cols(module_annot, caution_df)
module_annot$Module_CleanPlotLabel <- paste0(module_annot$ModuleID, " | ", module_annot$cleaned_biological_label_short)
module_annot$module_biological_label <- dplyr::coalesce(
  informative_member_label(module_annot$cleaned_biological_label),
  informative_member_label(module_annot$ModuleLabel_Final),
  informative_member_label(module_annot$module_label),
  informative_member_label(module_annot$raw_top_GO_label),
  as.character(module_annot$ModuleID)
)
module_annot$module_biological_label_short <- dplyr::coalesce(
  informative_member_label(module_annot$cleaned_biological_label_short),
  informative_member_label(module_annot$module_biological_label)
)
module_annot$module_label_display <- dplyr::coalesce(
  informative_member_label(module_annot$Module_CleanPlotLabel),
  paste(module_annot$ModuleID, module_annot$module_biological_label_short, sep = " | ")
)
module_annot$raw_annotation_label <- module_annot$raw_marker_or_signature_label
module_annot$cleaned_annotation_label <- module_annot$cleaned_biological_label
module_annot$safe_display_label <- dplyr::coalesce(
  informative_member_label(module_annot$Module_CleanPlotLabel),
  informative_member_label(module_annot$module_label_display),
  informative_member_label(module_annot$module_biological_label),
  as.character(module_annot$ModuleID)
)
module_annot$label_confidence <- module_annot$annotation_confidence
module_annot$label_basis <- module_annot$annotation_basis
module_annot$label_downgrade_reason <- module_annot$annotation_downgrade_reason

module_annot$member_theme_label <- vapply(seq_len(nrow(module_annot)), function(i) {
  row <- module_annot[i, , drop = FALSE]
  first_informative_value(
    row$module_specific_label,
    row$cleaned_biological_label,
    row$ModuleLabel_Final,
    row$module_label,
    row$ModuleLabel_GO_BP,
    row$best_GO_BP,
    split_label_terms(row$top_GO_BP_labels),
    split_label_terms(row$top_GO_MF_labels),
    split_label_terms(row$top_GO_CC_labels)
  )
}, character(1))
module_annot$member_theme_source <- vapply(seq_len(nrow(module_annot)), function(i) {
  row <- module_annot[i, , drop = FALSE]
  candidates <- list(
    module_specific_label = row$module_specific_label,
    cleaned_biological_label = row$cleaned_biological_label,
    ModuleLabel_Final = row$ModuleLabel_Final,
    module_label = row$module_label,
    ModuleLabel_GO_BP = row$ModuleLabel_GO_BP,
    best_GO_BP = row$best_GO_BP,
    top_GO_BP_labels = split_label_terms(row$top_GO_BP_labels),
    top_GO_MF_labels = split_label_terms(row$top_GO_MF_labels),
    top_GO_CC_labels = split_label_terms(row$top_GO_CC_labels)
  )
  hit <- names(candidates)[vapply(candidates, function(x) any(!is.na(informative_member_label(x))), logical(1))]
  if (length(hit)) hit[[1]] else if (length(split_tokens(row$top_hub_proteins))) "top_hub_proteins_audit_only" else "unresolved"
}, character(1))
module_annot$member_theme <- vapply(seq_len(nrow(module_annot)), function(i) {
  program <- module_annot$module_program_primary[[i]]
  if (!is.na(program) && nzchar(program)) return(program)
  x <- module_annot$member_theme_label[[i]]
  theme <- member_theme_from_label(x)
  if (is.na(theme) || !nzchar(theme)) "mixed / low-specificity" else theme
}, character(1))
module_annot$has_informative_member_theme <- !is.na(informative_member_label(module_annot$member_theme_label))
module_annot$SemanticManualReviewRequired <- suppressWarnings(as.logical(module_annot$ManualReviewFromSemanticLabel))
module_annot$marker_registry_version <- marker_registry_version
module_annot$empirical_marker_set_version <- empirical_marker_set_version
module_annot$interpretation_note <- WGCNA_ROI_NOTE

if (!is.null(module_effects) && nrow(module_effects)) {
  changed <- module_effects |>
    dplyr::filter(!is.na(.data$FDR_global), .data$FDR_global < 0.10) |>
    dplyr::group_by(.data$module_id) |>
    dplyr::summarise(Changed_in_group_contrasts = paste(unique(paste(.data$spatial_unit, .data$contrast, sep = ":")), collapse = ";"), .groups = "drop")
  module_annot <- module_annot |> dplyr::left_join(changed, by = c("ModuleID" = "module_id"))
}

if (is.null(super_ann) || !nrow(super_ann)) {
  super_annot <- data.frame(dataset = DATASET, SupermoduleID = NA_character_, supermodule_id = NA_character_, Supermodule_FinalLabel = NA_character_, n_member_modules = 0L, dominant_microenvironment_class = "missing_supermodule_annotation", interpretation_note = WGCNA_ROI_NOTE)
} else {
  super_ann2 <- super_ann
  for (nm in c("Supermodule_DataDrivenID", "Supermodule_DisplayLabel", "Supermodule_LongLabel", "Macroprogram_Display", "Supermodule_DataDrivenLabel", "Supermodule_CuratedLabel", "Supermodule_FinalLabel", "Supermodule_LabelSource", "Supermodule_LabelConfidence", "Supermodule_LabelRationale", "GO_label_confidence_class", "annotation_scope", "manual_label_status", "ManualReviewRequired", "Supermodule_DataDriven", "Supermodule", "SupermoduleConfidence", "SupermoduleRationale")) {
    if (!nm %in% names(super_ann2)) super_ann2[[nm]] <- NA_character_
  }
  smap <- super_ann2 |>
    dplyr::mutate(
      SupermoduleID = dplyr::coalesce(as.character(.data$Supermodule_DataDrivenID), as.character(.data$Supermodule_DataDriven), as.character(.data$Supermodule)),
      Supermodule_FinalLabel = dplyr::coalesce(as.character(.data$Supermodule_FinalLabel), as.character(.data$Supermodule), .data$SupermoduleID),
      Supermodule_LongLabel = dplyr::coalesce(as.character(.data$Supermodule_LongLabel), .data$Supermodule_FinalLabel),
      Macroprogram_Display = dplyr::coalesce(as.character(.data$Macroprogram_Display), macroprogram_display(.data$Supermodule_LongLabel)),
      Supermodule_DisplayLabel = dplyr::coalesce(as.character(.data$Supermodule_DisplayLabel), compose_supermodule_display_label(.data$SupermoduleID, dplyr::coalesce(.data$Supermodule_FinalLabel, .data$Macroprogram_Display))),
      Supermodule_ShortLabel = .data$Supermodule_DisplayLabel
    ) |>
    dplyr::select(dplyr::any_of(c("ModuleColor", "module_eigengene", "SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_LongLabel", "Macroprogram_Display", "Supermodule_DataDrivenLabel", "Supermodule_CuratedLabel", "Supermodule_FinalLabel", "Supermodule_ShortLabel", "Supermodule_LabelSource", "Supermodule_LabelConfidence", "Supermodule_LabelRationale", "GO_label_confidence_class", "annotation_scope", "manual_label_status", "ManualReviewRequired", "SupermoduleConfidence", "SupermoduleRationale")))
  super_annot <- module_annot |>
    dplyr::left_join(smap, by = c("ModuleColor" = "ModuleColor")) |>
    dplyr::filter(!is.na(.data$SupermoduleID)) |>
    dplyr::group_by(.data$dataset, .data$SupermoduleID, .data$Supermodule_DisplayLabel, .data$Supermodule_LongLabel, .data$Macroprogram_Display, .data$Supermodule_FinalLabel, .data$Supermodule_ShortLabel) |>
    dplyr::summarise(
      n_member_modules = dplyr::n_distinct(.data$ModuleID),
      member_modules = paste(unique(.data$ModuleID), collapse = ";"),
      fraction_modules_microglia_supported = mean(as.character(.data$microenvironment_class) == "microglia_supported", na.rm = TRUE),
      fraction_modules_microglia_state_or_activation_supported = mean(as.character(.data$microenvironment_class) == "microglia_state_or_activation_supported", na.rm = TRUE),
      fraction_modules_shared_microenvironment = mean(as.character(.data$microenvironment_class) == "shared_microenvironment", na.rm = TRUE),
      fraction_modules_neuropil_sensitive = mean(as.character(.data$microenvironment_class) == "neuropil_sensitive", na.rm = TRUE),
      fraction_modules_vascular_basement_membrane_ecm = mean(as.character(.data$microenvironment_class) == "vascular_basement_membrane_ecm", na.rm = TRUE),
      fraction_modules_vascular_bbb_mural = mean(as.character(.data$microenvironment_class) == "vascular_bbb_mural", na.rm = TRUE),
      fraction_modules_astrocyte_or_endfoot_sensitive = mean(as.character(.data$microenvironment_class) == "astrocyte_or_endfoot_sensitive", na.rm = TRUE),
      fraction_modules_oligodendrocyte_or_myelin_sensitive = mean(as.character(.data$microenvironment_class) == "oligodendrocyte_or_myelin_sensitive", na.rm = TRUE),
      fraction_modules_peripheral_myeloid_caution = mean(as.character(.data$microenvironment_class) == "peripheral_myeloid_caution", na.rm = TRUE),
      fraction_modules_ambiguous_or_mixed = mean(as.character(.data$microenvironment_class) == "ambiguous_or_mixed", na.rm = TRUE),
      dominant_microenvironment_class = names(sort(table(.data$microenvironment_class), decreasing = TRUE))[1],
      dominant_module_labels = compact_label(.data$microenvironment_label, 3L),
      dominant_GO_terms = paste(utils::head(unique(c(split_label_terms(.data$raw_GO_BP_terms), split_label_terms(.data$raw_GO_MF_terms), split_label_terms(.data$raw_GO_CC_terms))), 10), collapse = ";"),
      top_hub_proteins = paste(utils::head(unique(split_tokens(.data$top_hub_proteins)), 30), collapse = ";"),
      raw_GO_BP_terms = wgcna_collapse_terms(split_label_terms(.data$raw_GO_BP_terms), n = 20L),
      raw_GO_MF_terms = wgcna_collapse_terms(split_label_terms(.data$raw_GO_MF_terms), n = 20L),
      raw_GO_CC_terms = wgcna_collapse_terms(split_label_terms(.data$raw_GO_CC_terms), n = 20L),
      raw_top_GO_label = wgcna_collapse_terms(split_label_terms(.data$raw_top_GO_label), n = 6L),
      raw_module_label = compact_label(.data$raw_module_label, 8L),
      raw_hub_proteins = paste(utils::head(unique(split_tokens(.data$raw_hub_proteins)), 30), collapse = ";"),
      raw_marker_or_signature_label = compact_label(.data$raw_marker_or_signature_label, 6L),
      MemberModuleSpecificLabels = compact_label(.data$module_specific_label, 8L),
      MemberProgramPrimaryCounts = wgcna_member_theme_summary(.data$module_program_primary)$counts_label,
      MemberProgramSecondary = compact_label(.data$module_program_secondary, 8L),
      MemberLabelDecisionRules = compact_label(.data$label_decision_rule, 8L),
      MemberLabelWarnings = compact_label(.data$label_warning, 8L),
      MemberAdhesionInterpretations = compact_label(.data$adhesion_interpretation, 8L),
      MemberEvidenceBP = wgcna_collapse_terms(split_label_terms(.data$evidence_BP), n = 20L),
      MemberEvidenceMF = wgcna_collapse_terms(split_label_terms(.data$evidence_MF), n = 20L),
      MemberEvidenceCC = wgcna_collapse_terms(split_label_terms(.data$evidence_CC), n = 20L),
      MemberEvidenceHubs = paste(utils::head(unique(split_tokens(.data$evidence_hubs)), 30), collapse = ";"),
      MemberEvidenceMarkerPanels = compact_label(.data$evidence_marker_panels, 10L),
      microenvironment_caution_label = names(sort(table(.data$microenvironment_caution_label), decreasing = TRUE))[1],
      microenvironment_caution_class = names(sort(table(.data$microenvironment_caution_class), decreasing = TRUE))[1],
      microenvironment_caution_rationale = compact_label(.data$microenvironment_caution_rationale, 4L),
      MemberThemeCounts = wgcna_member_theme_summary(.data$member_theme)$counts_label,
      MemberThemeFractions = wgcna_member_theme_summary(.data$member_theme)$fractions_label,
      n_distinct_member_themes = wgcna_member_theme_summary(.data$member_theme)$n_distinct,
      is_multi_theme_supermodule = wgcna_member_theme_summary(.data$member_theme)$is_multi,
      themes_above_display_threshold = wgcna_member_theme_summary(.data$member_theme)$themes_above_display_threshold,
      Supermodule_CompositionLabel_FromAllMemberThemes = wgcna_member_theme_summary(.data$member_theme)$display_label,
      themes_omitted_from_display_label = wgcna_member_theme_summary(.data$member_theme)$themes_omitted_from_display_label,
      DominantMemberTheme = {
        theme_summary <- wgcna_member_theme_summary(.data$member_theme)
        if (length(theme_summary$fractions)) names(theme_summary$fractions)[[1]] else NA_character_
      },
      DominantMemberThemeFraction = {
        theme_summary <- wgcna_member_theme_summary(.data$member_theme)
        if (length(theme_summary$fractions)) as.numeric(theme_summary$fractions[[1]]) else NA_real_
      },
      SecondMemberTheme = {
        theme_summary <- wgcna_member_theme_summary(.data$member_theme)
        if (length(theme_summary$fractions) >= 2L) names(theme_summary$fractions)[[2]] else NA_character_
      },
      SecondMemberThemeFraction = {
        theme_summary <- wgcna_member_theme_summary(.data$member_theme)
        if (length(theme_summary$fractions) >= 2L) as.numeric(theme_summary$fractions[[2]]) else NA_real_
      },
      TopMemberModuleLabels = compact_label(informative_member_label(.data$cleaned_biological_label), 8L),
      TopMemberGOTerms = paste(utils::head(unique(c(split_label_terms(.data$raw_GO_BP_terms), split_label_terms(.data$raw_GO_MF_terms), split_label_terms(.data$raw_GO_CC_terms))), 12), collapse = ";"),
      n_member_modules_with_informative_labels = sum(.data$has_informative_member_theme, na.rm = TRUE),
      fraction_member_modules_with_informative_labels = mean(as.logical(.data$has_informative_member_theme), na.rm = TRUE),
      any_member_semantic_manual_review = any(.data$SemanticManualReviewRequired, na.rm = TRUE),
      member_theme_fraction_summary = wgcna_member_theme_summary(.data$member_theme)$fractions_label,
      member_theme_label_sources = compact_label(.data$member_theme_source, 5L),
      Supermodule_DataDrivenLabel = dplyr::first(.data$Supermodule_DataDrivenLabel),
      Supermodule_CuratedLabel = dplyr::first(.data$Supermodule_CuratedLabel),
      Supermodule_LabelSource = dplyr::first(.data$Supermodule_LabelSource),
      Supermodule_LabelConfidence = dplyr::first(.data$Supermodule_LabelConfidence),
      Supermodule_LabelRationale = dplyr::first(.data$Supermodule_LabelRationale),
      GO_label_confidence_class = dplyr::first(.data$GO_label_confidence_class),
      annotation_scope = dplyr::first(.data$annotation_scope),
      manual_label_status = dplyr::first(.data$manual_label_status),
      ManualReviewRequired = dplyr::first(.data$ManualReviewRequired),
      label_confidence = dplyr::coalesce(dplyr::first(.data$Supermodule_LabelConfidence), dplyr::first(.data$SupermoduleConfidence %||% NA_character_)),
      label_decision_rule = "member_module_semantic_composition",
      label_warning = compact_label(c(.data$label_warning, .data$themes_omitted_from_display_label), 10L),
      adhesion_interpretation = compact_label(.data$adhesion_interpretation, 8L),
      evidence_BP = wgcna_collapse_terms(split_label_terms(.data$evidence_BP), n = 20L),
      evidence_MF = wgcna_collapse_terms(split_label_terms(.data$evidence_MF), n = 20L),
      evidence_CC = wgcna_collapse_terms(split_label_terms(.data$evidence_CC), n = 20L),
      evidence_hubs = paste(utils::head(unique(split_tokens(.data$evidence_hubs)), 30), collapse = ";"),
      evidence_marker_panels = compact_label(.data$evidence_marker_panels, 10L),
      marker_registry_version = dplyr::first(.data$marker_registry_version),
      empirical_marker_set_version = dplyr::first(.data$empirical_marker_set_version),
      interpretation_note = WGCNA_ROI_NOTE,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Supermodule_FinalLabel_Original = .data$Supermodule_FinalLabel,
      Supermodule_CompositionLabel = .data$Supermodule_CompositionLabel_FromAllMemberThemes,
      Supermodule_CompositionLabel = if (DATASET == "microglia" || force_microglia) microglia_guarded_composition_label(.data$Supermodule_CompositionLabel, .data$dominant_microenvironment_class) else .data$Supermodule_CompositionLabel,
      Supermodule_CompositionDisplayLabel = compose_supermodule_composition_display(.data$Supermodule_CompositionLabel, .data$microenvironment_caution_label, dataset = DATASET),
      Supermodule_CompositionShortLabel = composition_short_label(.data$Supermodule_CompositionLabel),
      Supermodule_CleanPlotLabel = compose_supermodule_display_label(.data$SupermoduleID, .data$Supermodule_CompositionDisplayLabel),
      Supermodule_CompositionLabelSource = dplyr::case_when(
        .data$member_theme_label_sources == "top_hub_proteins_audit_only" ~ "hub_audit_only",
        .data$n_member_modules_with_informative_labels > 0L ~ "member_module_biological_label_or_GO",
        TRUE ~ "unresolved"
      ),
      Supermodule_CompositionConfidence = dplyr::case_when(
        .data$n_member_modules <= 1L ~ "low",
        .data$Supermodule_CompositionLabel == "mixed / low-specificity" ~ "low",
        grepl("^mixed:", .data$Supermodule_CompositionLabel) ~ "low",
        .data$DominantMemberThemeFraction >= 0.60 & .data$fraction_member_modules_with_informative_labels >= 0.60 ~ "medium",
        TRUE ~ "low"
      ),
      Supermodule_CompositionRationale = paste0(
        "member_theme_fractions: ", .data$member_theme_fraction_summary,
        "; member_theme_counts: ", .data$MemberThemeCounts,
        "; themes_above_display_threshold=", dplyr::coalesce(.data$themes_above_display_threshold, ""),
        "; themes_omitted_from_display_label=", dplyr::coalesce(.data$themes_omitted_from_display_label, ""),
        "; labels=", dplyr::coalesce(.data$TopMemberModuleLabels, ""),
        "; GO=", dplyr::coalesce(.data$TopMemberGOTerms, ""),
        "; hubs_audit=", dplyr::coalesce(.data$top_hub_proteins, "")
      ),
      Supermodule_ConservativeLabel = dplyr::case_when(
        .data$n_member_modules <= 1L ~ "mixed / unresolved (singleton supermodule)",
        grepl("^mostly ", .data$Supermodule_CompositionLabel) & .data$Supermodule_CompositionConfidence %in% c("medium", "high") ~ .data$Supermodule_CompositionLabel,
        TRUE ~ "mixed / unresolved"
      ),
      cleaned_biological_label = .data$Supermodule_CompositionLabel,
      cleaned_biological_label_short = .data$Supermodule_CompositionShortLabel,
      cleaned_biological_label_source = .data$Supermodule_CompositionLabelSource,
      cleaned_biological_label_confidence = .data$Supermodule_CompositionConfidence,
      cleaned_biological_label_rationale = .data$Supermodule_CompositionRationale,
      module_specific_label = .data$Supermodule_CompositionLabel,
      module_program_primary = .data$DominantMemberTheme,
      module_program_secondary = .data$SecondMemberTheme,
      GO_label_relevance_flag = dplyr::case_when(
        grepl("skin development|dorsal ventral pattern|binding sperm|zona pellucida|fertilization", .data$raw_module_label, ignore.case = TRUE) ~ "cleaned_or_context_checked",
        grepl("barrier / cell-junction|ontology artefact|regulatory module", .data$Supermodule_CompositionLabel, ignore.case = TRUE) ~ "cleaned_or_context_checked",
        TRUE ~ "direct_or_suggestive"
      ),
      GO_label_relevance_rationale = .data$Supermodule_CompositionRationale,
      Supermodule_FinalLabel = dplyr::case_when(
        generic_supermodule_label(.data$Supermodule_FinalLabel) & !is.na(.data$dominant_module_labels) & nzchar(.data$dominant_module_labels) & !grepl("mixed|low-specificity|unresolved", .data$dominant_module_labels, ignore.case = TRUE) ~ .data$dominant_module_labels,
        TRUE ~ .data$Supermodule_FinalLabel
      ),
      Supermodule_LabelRationale = paste0(
        dplyr::coalesce(.data$Supermodule_LabelRationale, ""),
        ifelse(generic_supermodule_label(.data$Supermodule_FinalLabel_Original), "; generic_supermodule_label_retained_or_microenvironment_fallback; composition_label_from_member_module_biology_available", "")
      ),
      ManualReviewRequired = dplyr::case_when(
        generic_supermodule_label(.data$Supermodule_FinalLabel_Original) ~ TRUE,
        .data$dominant_microenvironment_class %in% c("ambiguous_or_mixed", "shared_microenvironment") ~ TRUE,
        .data$Supermodule_CompositionConfidence == "low" ~ TRUE,
        .data$any_member_semantic_manual_review ~ TRUE,
        TRUE ~ suppressWarnings(as.logical(.data$ManualReviewRequired))
      )
    )
  if (DATASET == "microglia" || force_microglia) {
    super_annot <- super_annot |>
      dplyr::mutate(
        Supermodule_FinalLabel = dplyr::case_when(
          .data$dominant_microenvironment_class == "vascular_basement_membrane_ecm" & !grepl("vascular|ECM|basement|BBB", .data$Supermodule_FinalLabel, ignore.case = TRUE) ~ paste0("perivascular ECM / ", .data$Supermodule_FinalLabel),
          .data$dominant_microenvironment_class == "vascular_bbb_mural" & !grepl("vascular|BBB|mural|pericyte", .data$Supermodule_FinalLabel, ignore.case = TRUE) ~ paste0("BBB / vascular mural / ", .data$Supermodule_FinalLabel),
          .data$dominant_microenvironment_class == "shared_microenvironment" & grepl("synap|vesicle|neur", .data$Supermodule_FinalLabel, ignore.case = TRUE) ~ paste0("shared local microenvironment / ", .data$Supermodule_FinalLabel),
          .data$dominant_microenvironment_class == "neuropil_sensitive" & grepl("synap|vesicle|neur", .data$Supermodule_FinalLabel, ignore.case = TRUE) ~ paste0("neuropil-sensitive ", .data$Supermodule_FinalLabel),
          .data$dominant_microenvironment_class == "microglia_supported" & grepl("phago|lyso|complement|immune", .data$Supermodule_FinalLabel, ignore.case = TRUE) ~ paste0("microglia-supported ", .data$Supermodule_FinalLabel),
          .data$dominant_microenvironment_class == "ambiguous_or_mixed" ~ paste0(.data$Supermodule_FinalLabel, " / mixed local ROI program"),
          TRUE ~ .data$Supermodule_FinalLabel
        )
      )
  }
  if (!"supermodule_id" %in% names(super_annot)) super_annot$supermodule_id <- NA_character_
  super_annot <- super_annot |>
    dplyr::mutate(
      SupermoduleID = dplyr::coalesce(as.character(.data$SupermoduleID), as.character(.data$supermodule_id)),
      supermodule_id = dplyr::coalesce(as.character(.data$supermodule_id), as.character(.data$SupermoduleID)),
      microenvironment_label_component = supermodule_microenvironment_label(.data$dominant_microenvironment_class, dataset = DATASET),
      has_active_manual_label = .data$manual_label_status == "manual_label_present_for_dataset_module" |
        (!is.na(.data$Supermodule_CuratedLabel) & nzchar(as.character(.data$Supermodule_CuratedLabel))),
      display_label_component = dplyr::case_when(
        .data$has_active_manual_label ~ shorten_supermodule_label(.data$Supermodule_FinalLabel, max_chars = 30),
        !is.na(.data$Supermodule_CompositionDisplayLabel) & nzchar(.data$Supermodule_CompositionDisplayLabel) & !grepl("mixed / low-specificity|mixed / unresolved", .data$Supermodule_CompositionDisplayLabel, ignore.case = TRUE) ~ .data$Supermodule_CompositionDisplayLabel,
        !is.na(.data$microenvironment_label_component) & nzchar(.data$microenvironment_label_component) ~ .data$microenvironment_label_component,
        !is.na(.data$Macroprogram_Display) & nzchar(.data$Macroprogram_Display) & .data$Macroprogram_Display != "Unresolved / mixed" ~ .data$Macroprogram_Display,
        TRUE ~ shorten_supermodule_label(.data$Supermodule_FinalLabel, max_chars = 30)
      ),
      Supermodule_DisplayLabel = compose_supermodule_display_label(.data$SupermoduleID, .data$display_label_component),
      Supermodule_ShortLabel = .data$Supermodule_DisplayLabel,
      Supermodule_CleanPlotLabel = dplyr::coalesce(.data$Supermodule_CleanPlotLabel, compose_supermodule_display_label(.data$SupermoduleID, .data$display_label_component)),
      Supermodule_PlotLabel = .data$Supermodule_CleanPlotLabel,
      Supermodule_LabelSource = dplyr::case_when(
        .data$has_active_manual_label ~ .data$Supermodule_LabelSource,
        !is.na(.data$microenvironment_label_component) & nzchar(.data$microenvironment_label_component) ~ paste0(dplyr::coalesce(.data$Supermodule_LabelSource, "data_driven"), "+microenvironment_annotation"),
        TRUE ~ .data$Supermodule_LabelSource
      ),
      has_coherent_hubs = lengths(strsplit(dplyr::coalesce(.data$top_hub_proteins, ""), "[;,]")) >= 3L,
      Supermodule_LabelConfidence = vapply(seq_len(dplyr::n()), function(i) {
        classify_supermodule_label_confidence(
          n_modules = n_member_modules[[i]],
          go_class = GO_label_confidence_class[[i]],
          has_coherent_hubs = has_coherent_hubs[[i]],
          microenvironment_class = dominant_microenvironment_class[[i]]
        )
      }, character(1)),
      Supermodule_LabelConfidence = dplyr::if_else(.data$has_active_manual_label & !is.na(.data$label_confidence), as.character(.data$label_confidence), .data$Supermodule_LabelConfidence),
      Supermodule_LabelConfidence = dplyr::if_else(.data$n_member_modules <= 1L & .data$Supermodule_LabelConfidence == "high", "low", .data$Supermodule_LabelConfidence),
      ManualReviewRequired = dplyr::case_when(
        .data$dominant_microenvironment_class %in% c("ambiguous_or_mixed", "shared_microenvironment") ~ TRUE,
        .data$n_member_modules <= 1L ~ TRUE,
        TRUE ~ suppressWarnings(as.logical(.data$ManualReviewRequired))
      ),
      Supermodule_LabelRationale = paste0(
        dplyr::coalesce(.data$Supermodule_LabelRationale, ""),
        ifelse(!is.na(.data$microenvironment_label_component) & nzchar(.data$microenvironment_label_component),
               paste0("; display_label_microenvironment_component=", .data$microenvironment_label_component),
               "")
      )
    ) |>
    dplyr::select(-dplyr::any_of(c("has_active_manual_label", "display_label_component", "has_coherent_hubs", "member_theme_fraction_summary", "member_theme_label_sources", "any_member_semantic_manual_review", "Supermodule_CompositionLabel_FromAllMemberThemes")))
}

if (nrow(super_annot) && exists("smap")) {
  super_threshold_classes <- module_threshold_classes |>
    dplyr::left_join(module_annot |> dplyr::select(module_or_supermodule_id = "ModuleID", ModuleColor), by = "module_or_supermodule_id") |>
    dplyr::left_join(smap |> dplyr::select("ModuleColor", "SupermoduleID"), by = "ModuleColor") |>
    dplyr::filter(!is.na(.data$SupermoduleID)) |>
    dplyr::group_by(.data$dataset, module_or_supermodule_id = .data$SupermoduleID) |>
    dplyr::summarise(
      class_at_0.05 = names(sort(table(.data$`class_at_0.05`), decreasing = TRUE))[1],
      class_at_0.10 = names(sort(table(.data$`class_at_0.10`), decreasing = TRUE))[1],
      class_at_0.20 = names(sort(table(.data$`class_at_0.20`), decreasing = TRUE))[1],
      .groups = "drop"
    )
  super_annot <- super_annot |>
    dplyr::left_join(super_threshold_classes, by = c("dataset", "SupermoduleID" = "module_or_supermodule_id")) |>
    dplyr::mutate(
      annotation_stable_across_thresholds = .data$`class_at_0.05` == .data$`class_at_0.10` & .data$`class_at_0.10` == .data$`class_at_0.20`,
      marker_fraction_primary = dplyr::coalesce(suppressWarnings(as.numeric(.data$DominantMemberThemeFraction)), suppressWarnings(as.numeric(.data$fraction_member_modules_with_informative_labels))),
      marker_panels_supporting = .data$MemberEvidenceMarkerPanels,
      annotation_basis = dplyr::case_when(
        !is.na(.data$marker_panels_supporting) & nzchar(as.character(.data$marker_panels_supporting)) ~ "marker_panel",
        !is.na(.data$TopMemberGOTerms) & nzchar(as.character(.data$TopMemberGOTerms)) ~ "GO",
        .data$n_member_modules_with_informative_labels > 0L ~ "hub/member_composition",
        TRUE ~ "insufficient"
      ),
      annotation_confidence = dplyr::case_when(
        !.data$annotation_stable_across_thresholds ~ "low",
        .data$Supermodule_CompositionConfidence %in% c("high") ~ "high",
        .data$Supermodule_CompositionConfidence %in% c("medium", "moderate") ~ "moderate",
        .data$annotation_basis == "insufficient" ~ "unresolved",
        TRUE ~ "low"
      ),
      annotation_downgrade_reason = vapply(seq_len(dplyr::n()), function(i) annotation_downgrade_for(`class_at_0.05`[[i]], `class_at_0.10`[[i]], `class_at_0.20`[[i]], annotation_basis[[i]], annotation_confidence[[i]]), character(1)),
      unsafe_interpretation = "Do not treat supermodule composition annotation as a purified cell-type or causal claim without claim-gate support.",
      raw_annotation_label = .data$raw_marker_or_signature_label,
      cleaned_annotation_label = .data$cleaned_biological_label,
      safe_display_label = dplyr::coalesce(.data$Supermodule_CleanPlotLabel, .data$Supermodule_CompositionDisplayLabel, .data$Supermodule_DisplayLabel, .data$SupermoduleID),
      label_confidence = .data$annotation_confidence,
      label_basis = .data$annotation_basis,
      label_downgrade_reason = .data$annotation_downgrade_reason
    )
}

supermodule_composition_columns <- c(
  "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms",
  "raw_top_GO_label", "raw_module_label", "raw_hub_proteins",
  "raw_marker_or_signature_label",
  "cleaned_biological_label", "cleaned_biological_label_short",
  "cleaned_biological_label_source", "cleaned_biological_label_confidence",
  "cleaned_biological_label_rationale", "GO_label_relevance_flag",
  "GO_label_relevance_rationale",
  "module_specific_label", "module_program_primary", "module_program_secondary",
  "label_decision_rule", "label_warning", "adhesion_interpretation",
  "evidence_BP", "evidence_MF", "evidence_CC", "evidence_hubs",
  "evidence_marker_panels",
  "microenvironment_caution_label", "microenvironment_caution_class",
  "microenvironment_caution_rationale",
  "Supermodule_ConservativeLabel", "Supermodule_CompositionLabel",
  "Supermodule_CompositionShortLabel", "Supermodule_CompositionDisplayLabel",
  "Supermodule_CompositionLabelSource",
  "Supermodule_CompositionConfidence", "Supermodule_CompositionRationale",
  "Supermodule_CleanPlotLabel", "Supermodule_PlotLabel",
  "DominantMemberTheme", "DominantMemberThemeFraction", "SecondMemberTheme",
  "SecondMemberThemeFraction", "TopMemberModuleLabels", "TopMemberGOTerms",
  "MemberThemeCounts", "MemberThemeFractions", "n_distinct_member_themes",
  "is_multi_theme_supermodule", "themes_above_display_threshold",
  "themes_omitted_from_display_label",
  "MemberModuleSpecificLabels", "MemberProgramPrimaryCounts",
  "MemberProgramSecondary", "MemberLabelDecisionRules", "MemberLabelWarnings",
  "MemberAdhesionInterpretations", "MemberEvidenceBP", "MemberEvidenceMF",
  "MemberEvidenceCC", "MemberEvidenceHubs", "MemberEvidenceMarkerPanels",
  "n_member_modules_with_informative_labels",
  "fraction_member_modules_with_informative_labels",
  "annotation_confidence", "annotation_basis", "annotation_downgrade_reason",
  "annotation_stable_across_thresholds", "unsafe_interpretation",
  "raw_annotation_label", "cleaned_annotation_label", "safe_display_label",
  "label_confidence", "label_basis", "label_downgrade_reason",
  "marker_fraction_primary", "marker_panels_supporting",
  "class_at_0.05", "class_at_0.10", "class_at_0.20"
)
for (nm in supermodule_composition_columns) {
  if (!nm %in% names(super_annot)) super_annot[[nm]] <- NA
}

audit_dir <- path_results("reviewer_audit")
dir_create(audit_dir)
replace_dataset_audit <- function(path, rows) {
  old <- safe_read_csv(path)
  if (!is.null(old) && nrow(old) && "dataset" %in% names(old)) {
    old <- old |> dplyr::filter(.data$dataset != DATASET)
  }
  if (!is.null(old)) old[] <- lapply(old, as.character)
  rows[] <- lapply(rows, as.character)
  out <- dplyr::bind_rows(old, rows)
  readr::write_csv(out, path, na = "")
  invisible(out)
}

module_threshold_audit <- module_threshold_classes |>
  dplyr::mutate(
    level = "module",
    primary_class = module_annot$microenvironment_class,
    annotation_stable_across_thresholds = module_annot$annotation_stable_across_thresholds,
    classification_flip_reason = dplyr::case_when(
      .data$annotation_stable_across_thresholds ~ "not_applicable",
      TRUE ~ paste0("threshold_sensitive: ", .data$`class_at_0.05`, " -> ", .data$`class_at_0.10`, " -> ", .data$`class_at_0.20`)
    ),
    marker_fraction_primary = module_annot$marker_fraction_primary,
    marker_panels_supporting = module_annot$marker_panels_supporting,
    claim_relevance = dplyr::case_when(
      module_annot$annotation_confidence %in% c("high", "moderate") & module_annot$annotation_stable_across_thresholds ~ "annotation_context",
      TRUE ~ "annotation_only_or_caution"
    )
  ) |>
  dplyr::select(
    "dataset", "module_or_supermodule_id", "level",
    "class_at_0.05", "class_at_0.10", "class_at_0.20",
    "primary_class", "annotation_stable_across_thresholds",
    "classification_flip_reason", "marker_fraction_primary",
    "marker_panels_supporting", "claim_relevance"
  )

super_threshold_audit <- if (nrow(super_annot)) {
  super_annot |>
    dplyr::transmute(
      dataset = .data$dataset,
      module_or_supermodule_id = .data$SupermoduleID,
      level = "supermodule",
      `class_at_0.05` = as.character(.data$`class_at_0.05`),
      `class_at_0.10` = as.character(.data$`class_at_0.10`),
      `class_at_0.20` = as.character(.data$`class_at_0.20`),
      primary_class = as.character(.data$dominant_microenvironment_class),
      annotation_stable_across_thresholds = suppressWarnings(as.logical(.data$annotation_stable_across_thresholds)),
      classification_flip_reason = dplyr::case_when(
        .data$annotation_stable_across_thresholds ~ "not_applicable",
        TRUE ~ paste0("threshold_sensitive: ", .data$`class_at_0.05`, " -> ", .data$`class_at_0.10`, " -> ", .data$`class_at_0.20`)
      ),
      marker_fraction_primary = suppressWarnings(as.numeric(.data$marker_fraction_primary)),
      marker_panels_supporting = as.character(.data$marker_panels_supporting),
      claim_relevance = dplyr::case_when(
        .data$annotation_confidence %in% c("high", "moderate") & .data$annotation_stable_across_thresholds ~ "annotation_context",
        TRUE ~ "annotation_only_or_caution"
      )
    )
} else {
  module_threshold_audit[0, ]
}
threshold_audit <- dplyr::bind_rows(module_threshold_audit, super_threshold_audit)
replace_dataset_audit(file.path(audit_dir, "wgcna_microenvironment_threshold_sensitivity.csv"), threshold_audit)

label_confidence_audit <- dplyr::bind_rows(
  module_annot |>
    dplyr::transmute(
      dataset = .data$dataset, module_or_supermodule_id = .data$ModuleID, level = "module",
      raw_annotation_label = .data$raw_annotation_label, cleaned_annotation_label = .data$cleaned_annotation_label,
      safe_display_label = .data$safe_display_label, label_confidence = .data$label_confidence,
      label_basis = .data$label_basis, label_downgrade_reason = .data$label_downgrade_reason,
      annotation_stable_across_thresholds = .data$annotation_stable_across_thresholds,
      unsafe_interpretation = .data$unsafe_interpretation
    ),
  super_annot |>
    dplyr::transmute(
      dataset = .data$dataset, module_or_supermodule_id = .data$SupermoduleID, level = "supermodule",
      raw_annotation_label = .data$raw_annotation_label, cleaned_annotation_label = .data$cleaned_annotation_label,
      safe_display_label = .data$safe_display_label, label_confidence = .data$label_confidence,
      label_basis = .data$label_basis, label_downgrade_reason = .data$label_downgrade_reason,
      annotation_stable_across_thresholds = suppressWarnings(as.logical(.data$annotation_stable_across_thresholds)),
      unsafe_interpretation = .data$unsafe_interpretation
    )
)
replace_dataset_audit(file.path(audit_dir, "wgcna_label_confidence_audit.csv"), label_confidence_audit)

annotation_source_audit <- supplemental_marker_panels |>
  dplyr::count(
    panel_version = dplyr::coalesce(as.character(.data$panel_version), "unversioned"),
    panel_id = .data$panel_id,
    source_type = .data$source_type,
    source_reference = .data$source_reference,
    allowed_use = .data$allowed_use,
    claim_role = .data$claim_role,
    caution_note = .data$caution_note,
    name = "n_markers"
  ) |>
  dplyr::mutate(
    dataset = DATASET,
    config_file = supplemental_marker_panel_file,
    config_hash = supplemental_marker_panel_hash,
    .before = "panel_version"
  )
replace_dataset_audit(file.path(audit_dir, "wgcna_annotation_source_audit.csv"), annotation_source_audit)

write_table_and_source(module_annot, PATHS$tables, PATHS$source_data, "WGCNA_module_biological_annotation.csv")
write_table_and_source(super_annot, PATHS$tables, PATHS$source_data, "WGCNA_supermodule_biological_annotation.csv")
write_table_and_source(targeted_signature_details, PATHS$tables, PATHS$source_data, "WGCNA_module_targeted_signature_overlap_details.csv")
if (nrow(super_annot)) {
  display_audit <- super_annot |>
    dplyr::mutate(
      is_singleton_supermodule = .data$n_member_modules <= 1L,
      manual_review_required = suppressWarnings(as.logical(.data$ManualReviewRequired)),
      label_rationale = .data$Supermodule_LabelRationale,
      top_GO_BP_terms = dplyr::coalesce(.data$dominant_GO_terms, NA_character_),
      top_hub_symbols = dplyr::coalesce(.data$top_hub_proteins, NA_character_)
    ) |>
    dplyr::select(dplyr::any_of(c(
      "dataset", "SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_FinalLabel",
      "Supermodule_CleanPlotLabel", "Supermodule_PlotLabel",
      "raw_GO_BP_terms", "raw_GO_MF_terms", "raw_GO_CC_terms",
      "raw_top_GO_label", "raw_module_label", "raw_hub_proteins",
      "raw_marker_or_signature_label",
      "cleaned_biological_label", "cleaned_biological_label_short",
      "cleaned_biological_label_source", "cleaned_biological_label_confidence",
      "cleaned_biological_label_rationale", "GO_label_relevance_flag",
      "GO_label_relevance_rationale",
      "module_specific_label", "module_program_primary", "module_program_secondary",
      "label_decision_rule", "label_warning", "adhesion_interpretation",
      "evidence_BP", "evidence_MF", "evidence_CC", "evidence_hubs",
      "evidence_marker_panels",
      "microenvironment_caution_label", "microenvironment_caution_class",
      "microenvironment_caution_rationale",
      "Supermodule_ConservativeLabel", "Supermodule_CompositionLabel",
      "Supermodule_CompositionShortLabel", "Supermodule_CompositionDisplayLabel",
      "Supermodule_CompositionLabelSource",
      "Supermodule_CompositionConfidence", "Supermodule_CompositionRationale",
      "DominantMemberTheme", "DominantMemberThemeFraction", "SecondMemberTheme",
      "SecondMemberThemeFraction", "TopMemberModuleLabels", "TopMemberGOTerms",
      "MemberThemeCounts", "MemberThemeFractions", "n_distinct_member_themes",
      "is_multi_theme_supermodule", "themes_above_display_threshold",
      "themes_omitted_from_display_label",
      "MemberModuleSpecificLabels", "MemberProgramPrimaryCounts",
      "MemberProgramSecondary", "MemberLabelDecisionRules", "MemberLabelWarnings",
      "MemberAdhesionInterpretations", "MemberEvidenceBP", "MemberEvidenceMF",
      "MemberEvidenceCC", "MemberEvidenceHubs", "MemberEvidenceMarkerPanels",
      "Supermodule_LabelSource", "Supermodule_LabelConfidence", "is_singleton_supermodule",
      "n_member_modules", "GO_label_confidence_class", "dominant_microenvironment_class",
      "n_member_modules_with_informative_labels", "fraction_member_modules_with_informative_labels",
      "top_GO_BP_terms", "top_hub_symbols", "label_rationale", "manual_review_required"
    ))) |>
    dplyr::rename(Supermodule_DataDriven = "SupermoduleID")
  write_table_and_source(display_audit, PATHS$tables, PATHS$source_data, "WGCNA_supermodule_display_label_audit.csv")
}
if (requireNamespace("writexl", quietly = TRUE)) {
  tryCatch(
    writexl::write_xlsx(list(modules = module_annot, supermodules = super_annot), file.path(PATHS$tables, "WGCNA_module_microenvironment_annotation.xlsx")),
    error = function(e) warning("Could not write WGCNA_module_microenvironment_annotation.xlsx: ", conditionMessage(e), call. = FALSE)
  )
}

marker_fraction_plot_cols <- paste0(names(marker_sets), "_fraction")
marker_fraction_plot_cols <- marker_fraction_plot_cols[marker_fraction_plot_cols %in% names(module_annot)]
marker_long <- module_annot |>
  dplyr::select("module_display_label", dplyr::all_of(marker_fraction_plot_cols)) |>
  tidyr::pivot_longer(cols = -dplyr::all_of("module_display_label"), names_to = "marker_panel", values_to = "fraction")
if (nrow(marker_long)) {
  p <- ggplot2::ggplot(marker_long, ggplot2::aes(x = .data$marker_panel, y = .data$module_display_label, fill = .data$fraction)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "grey97", high = "#2F6F73", na.value = "grey90") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Marker\nfraction", title = "Module marker evidence") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(file.path(PATHS$figures, "module_marker_fraction_heatmap.svg"), p, width = 170, height = 120, units = "mm", device = svglite::svglite)
}
if (nrow(super_annot) && "dominant_microenvironment_class" %in% names(super_annot)) {
  comp <- super_annot |>
    dplyr::mutate(
      Supermodule_DisplayLabel = dplyr::coalesce(.data$Supermodule_DisplayLabel, compose_supermodule_display_label(.data$SupermoduleID, dplyr::coalesce(.data$Supermodule_FinalLabel, .data$Macroprogram_Display)), .data$SupermoduleID),
      plot_label = vapply(as.character(.data$Supermodule_DisplayLabel), function(z) paste(strwrap(z, width = 34), collapse = "\n"), character(1))
    ) |>
    dplyr::select("plot_label", "SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_LongLabel", "Macroprogram_Display", dplyr::starts_with("fraction_modules_")) |>
    tidyr::pivot_longer(cols = dplyr::starts_with("fraction_modules_"), names_to = "class", values_to = "fraction") |>
    dplyr::mutate(
      class = gsub("^fraction_modules_", "", .data$class),
      class = gsub("_", " ", .data$class),
      fraction = suppressWarnings(as.numeric(.data$fraction)),
      supermodule_order = match(.data$SupermoduleID, sort(unique(.data$SupermoduleID)))
    ) |>
    dplyr::filter(is.finite(.data$fraction), .data$fraction > 0)
  write_table_and_source(comp, PATHS$tables, PATHS$source_data, "WGCNA_supermodule_microenvironment_composition_source.csv")
  if (nrow(comp)) {
    fig_h <- max(90, 14 * length(unique(comp$plot_label)) + 35)
    p2 <- ggplot2::ggplot(comp, ggplot2::aes(x = .data$fraction, y = stats::reorder(.data$plot_label, .data$supermodule_order), fill = .data$class)) +
      ggplot2::geom_col(width = 0.70) +
      ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = "Fraction of member modules", y = NULL, fill = NULL) +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 7.4))
    ggplot2::ggsave(file.path(PATHS$figures, "supermodule_microenvironment_composition.svg"), p2, width = 165, height = fig_h, units = "mm", device = svglite::svglite)
  }
}
if (DATASET == "microglia" || force_microglia) {
  x_col <- if ("empirical_neuropil_enriched_fraction" %in% names(module_annot)) "empirical_neuropil_enriched_fraction" else if ("canonical_neuronal_synaptic_neuropil_fraction" %in% names(module_annot)) "canonical_neuronal_synaptic_neuropil_fraction" else "neuropil_synaptic_neuronal_marker_fraction"
  y_col <- if ("empirical_microglia_roi_enriched_fraction" %in% names(module_annot)) "empirical_microglia_roi_enriched_fraction" else if ("canonical_microglia_homeostatic_fraction" %in% names(module_annot)) "canonical_microglia_homeostatic_fraction" else "microglia_marker_fraction"
  p3 <- ggplot2::ggplot(module_annot, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = .data$microenvironment_class)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "Neuropil marker fraction", y = "Microglia marker fraction", color = "Class") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(PATHS$figures, "microglia_vs_neuropil_module_evidence.svg"), p3, width = 120, height = 95, units = "mm", device = svglite::svglite)

  if ("n_targeted_microglia_signature_overlaps" %in% names(module_annot) &&
      any(module_annot$n_targeted_microglia_signature_overlaps > 0 |
            module_annot$n_targeted_curated_microglia_program_overlaps > 0 |
            module_annot$n_targeted_neuropil_shared_overlaps > 0, na.rm = TRUE)) {
    targeted_plot <- module_annot |>
      dplyr::mutate(
        microglia_enriched = suppressWarnings(as.numeric(.data$n_targeted_microglia_signature_overlaps)),
        curated_program = suppressWarnings(as.numeric(.data$n_targeted_curated_microglia_program_overlaps)),
        claim_ready = suppressWarnings(as.numeric(.data$n_targeted_claim_ready_signature_overlaps)),
        shared_neuropil = suppressWarnings(as.numeric(.data$n_targeted_neuropil_shared_overlaps)),
        module_display_label = dplyr::coalesce(.data$module_display_label, .data$ModuleID)
      ) |>
      dplyr::filter(.data$microglia_enriched > 0 | .data$curated_program > 0 | .data$claim_ready > 0 | .data$shared_neuropil > 0) |>
      dplyr::arrange(dplyr::desc(.data$claim_ready), dplyr::desc(.data$microglia_enriched), dplyr::desc(.data$shared_neuropil)) |>
      dplyr::mutate(module_display_label = stats::reorder(.data$module_display_label, .data$microglia_enriched + .data$curated_program + .data$shared_neuropil))
    targeted_long <- targeted_plot |>
      dplyr::select("module_display_label", microglia_enriched, curated_program, claim_ready, shared_neuropil) |>
      tidyr::pivot_longer(cols = -dplyr::all_of("module_display_label"), names_to = "evidence", values_to = "n")
    write_table_and_source(targeted_long, PATHS$tables, PATHS$source_data, "WGCNA_module_targeted_microglia_signature_overlap_source.csv")
    targeted_summary_source <- targeted_plot |>
      dplyr::transmute(
        dataset = .data$dataset,
        ModuleID = .data$ModuleID,
        ModuleColor = .data$ModuleColor,
        module_display_label = .data$module_display_label,
        n_signature_comparison_overlaps = .data$n_targeted_microglia_signature_overlaps + .data$n_targeted_curated_microglia_program_overlaps + .data$n_targeted_neuropil_shared_overlaps,
        n_claim_ready_signature_comparison_overlaps = .data$n_targeted_claim_ready_signature_overlaps,
        n_unique_targeted_signatures = .data$n_unique_targeted_signatures,
        n_unique_targeted_overlap_proteins = .data$n_unique_targeted_overlap_proteins,
        n_unique_curated_program_signatures = .data$n_unique_curated_program_signatures,
        n_unique_curated_program_overlap_proteins = .data$n_unique_curated_program_overlap_proteins,
        top_overlap_proteins = .data$best_targeted_signature_overlap_proteins,
        curated_program_overlap_proteins = .data$curated_program_overlap_proteins,
        caution_flag = .data$curated_program_overlap_warning,
        interpretation_note = .data$targeted_signature_overlap_interpretation_note,
        driver_evidence_tier = .data$targeted_signature_driver_evidence_tier,
        targeted_signature_primary_driver = .data$targeted_signature_primary_driver
      )
    write_table_and_source(targeted_summary_source, PATHS$tables, PATHS$source_data, "WGCNA_module_targeted_microglia_signature_overlap_summary_source.csv")
    p4 <- ggplot2::ggplot(targeted_long, ggplot2::aes(x = .data$n, y = .data$module_display_label, fill = .data$evidence)) +
      ggplot2::geom_col(position = "dodge", width = 0.68) +
      ggplot2::scale_fill_manual(
        values = c(microglia_enriched = "#3C5488", curated_program = "#7E6148", claim_ready = "#00A087", shared_neuropil = "#E64B35"),
        labels = c(microglia_enriched = "microglia-enriched", curated_program = "curated microglia-relevant", claim_ready = "claim-ready", shared_neuropil = "shared neuropil"),
        name = NULL
      ) +
      ggplot2::labs(x = "Overlapping targeted signatures", y = NULL) +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 7.2))
    ggplot2::ggsave(file.path(PATHS$figures, "microglia_module_targeted_signature_overlap.svg"), p4, width = 135, height = max(85, 7.5 * length(unique(targeted_long$module_display_label)) + 24), units = "mm", device = svglite::svglite, limitsize = FALSE)
  }
}

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = c(FILES, targeted_microglia_signature_enrichment = targeted_signature_file, supplemental_microenvironment_marker_panels = supplemental_marker_panel_file),
  outputs = list(
    tables = PATHS$tables,
    source_data = PATHS$source_data,
    figures = PATHS$figures,
    reviewer_audit = c(
      wgcna_microenvironment_threshold_sensitivity = file.path(audit_dir, "wgcna_microenvironment_threshold_sensitivity.csv"),
      wgcna_label_confidence_audit = file.path(audit_dir, "wgcna_label_confidence_audit.csv"),
      wgcna_annotation_source_audit = file.path(audit_dir, "wgcna_annotation_source_audit.csv")
    )
  ),
  parameters = list(
    dataset = DATASET,
    force_microglia_annotation = force_microglia,
    classification_threshold = classification_threshold,
    threshold_sensitivity_values = paste(marker_threshold_sensitivity_values, collapse = ","),
    marker_registry_version = marker_registry_version,
    empirical_marker_set_version = empirical_marker_set_version,
    supplemental_marker_panel_file = supplemental_marker_panel_file,
    supplemental_marker_panel_hash = supplemental_marker_panel_hash,
    supplemental_marker_panel_version = paste(unique(supplemental_marker_panels$panel_version %||% "unversioned"), collapse = ";")
  ),
  notes = paste(WGCNA_ROI_NOTE, "Improved annotation separates vascular basement membrane/BBB/mural, astrocyte/endfoot, oligodendrocyte/myelin, neuropil, and microglia-supported ROI evidence; labels are annotation only.")
)

message("WGCNA module/supermodule biological annotation complete for dataset: ", DATASET)
