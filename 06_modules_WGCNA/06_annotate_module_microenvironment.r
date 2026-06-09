#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/06_annotate_module_microenvironment.r
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/; optional config/marker_panels/wgcna_reference_marker_sets.csv; results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv; +1 more.
# Produces: results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv; results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv.
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
  dry_run_line("Output tables", PATHS$tables)
  quit(status = 0, save = "no")
}

if (DATASET != "microglia" && force_microglia) {
  warning("Forcing microglia-specific module annotation outside dataset == microglia.", call. = FALSE)
}

definitions <- safe_read_csv(FILES$definitions)
if (is.null(definitions) || !nrow(definitions)) {
  module_annot <- data.frame(dataset = DATASET, ModuleID = NA_character_, ModuleColor = NA_character_, n_proteins = 0L, microenvironment_class = "missing_module_definitions", interpretation_note = WGCNA_ROI_NOTE)
  super_annot <- data.frame(dataset = DATASET, SupermoduleID = NA_character_, n_member_modules = 0L, dominant_microenvironment_class = "missing_module_definitions", interpretation_note = WGCNA_ROI_NOTE)
  write_table_and_source(module_annot, PATHS$tables, PATHS$source_data, "WGCNA_module_biological_annotation.csv")
  write_table_and_source(super_annot, PATHS$tables, PATHS$source_data, "WGCNA_supermodule_biological_annotation.csv")
  write_run_manifest(file.path(PATHS$logs, "run_manifest.yml"), inputs = FILES, outputs = list(tables = PATHS$tables), parameters = list(dataset = DATASET), notes = paste(WGCNA_ROI_NOTE, "Improved annotation separates vascular basement membrane/BBB/mural, astrocyte/endfoot, oligodendrocyte/myelin, neuropil, and microglia-supported ROI evidence; labels are annotation only."))
  quit(status = 0, save = "no")
}

validate_wgcna_module_definitions(definitions, "WGCNA downstream definitions")
marker_sets <- load_wgcna_marker_sets()
marker_registry_version <- attr(marker_sets, "marker_registry_version") %||% NA_character_
empirical_marker_set_version <- attr(marker_sets, "empirical_marker_set_version") %||% NA_character_
module_summary <- safe_read_csv(FILES$module_summary)
go <- safe_read_csv(FILES$go)
super_ann <- safe_read_csv(FILES$supermodule_annotation)
super_summary <- safe_read_csv(FILES$supermodule_summary)
module_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", DATASET, "module_group_effects.csv"))
super_effects <- safe_read_csv(path_results("tables", "06_modules_WGCNA", "group_effects", DATASET, "supermodule_group_effects.csv"))
neuropil_ref <- if (DATASET == "microglia" || force_microglia) safe_read_csv(FILES$neuropil_annotation) else NULL
targeted_signature_file <- path_results("tables", "04_differential_expression_enrichment", "microglia_targeted_signature_enrichment", "microglia", "microglia_signature_enrichment_with_contrast_class.csv")
targeted_signature_ref <- if (DATASET == "microglia" || force_microglia) safe_read_csv(targeted_signature_file) else NULL


# Add a small, local compartment marker layer for microenvironment interpretation.
# These sets supplement the registry; they do not replace curated external marker sets.
extra_microenvironment_marker_sets <- list(
  basement_membrane_perivascular_ecm = c(
    "Agrn", "Hspg2", "Lamb1", "Lamb2", "Lamc1", "Lama2", "Lama4", "Lama5",
    "Col4a1", "Col4a2", "Nid1", "Nid2", "Bcam", "Tinagl1", "Serpinh1"
  ),
  bbb_endothelial_transport = c(
    "Slc2a1", "Abcb1a", "Abcb1b", "Abcg2", "Cldn5", "Pecam1", "Kdr", "Flt1",
    "Tek", "Vwf", "Esam", "Mfsd2a", "Ocln", "Cav1"
  ),
  pericyte_mural_ng2 = c(
    "Cspg4", "Pdgfrb", "Rgs5", "Kcnj8", "Abcc9", "Des", "Acta2", "Tagln",
    "Mcam", "Notch3", "Anpep"
  ),
  astrocyte_endfoot_gliovascular = c(
    "Aqp4", "Gfap", "Aldh1l1", "Slc1a2", "Slc1a3", "Gja1", "Glul", "Aldoc",
    "Dtna", "Dag1"
  ),
  integrin_ecm_adhesion = c(
    "Itga1", "Itga2", "Itga5", "Itga6", "Itgb1", "Tln1", "Tln2", "Vcl",
    "Parva", "Parvb", "Pxn"
  )
)

for (nm in names(extra_microenvironment_marker_sets)) {
  if (!nm %in% names(marker_sets)) marker_sets[[nm]] <- extra_microenvironment_marker_sets[[nm]]
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
  targeted_micro <- suppressWarnings(as.integer(row$n_targeted_microglia_signature_overlaps %||% 0))
  targeted_claim <- suppressWarnings(as.integer(row$n_targeted_claim_ready_signature_overlaps %||% 0))
  neuro <- max_frac(frac(c("empirical_neuropil_sensitive_high_confidence_fraction", "empirical_neuropil_enriched_fraction")), frac(c("canonical_neuronal_synaptic_neuropil_fraction", "neuropil_synaptic_neuronal_marker_fraction")))
  astro <- max_frac(frac(c("astrocyte_endfoot_gliovascular_fraction", "astrocyte_endfoot_gliovascular_marker_fraction")), frac(c("canonical_astrocyte_fraction", "astrocyte_marker_fraction")))
  oligo <- max_frac(frac(c("canonical_oligodendrocyte_myelin_fraction", "oligodendrocyte_myelin_marker_fraction")), frac(c("canonical_opc_fraction")))
  mito <- frac(c("canonical_mitochondrial_oxphos_fraction", "mitochondrial_oxphos_marker_fraction"))
  ribo <- frac(c("canonical_ribosomal_translation_fraction", "ribosomal_translation_marker_fraction"))
  rnp <- frac(c("canonical_rnp_rna_processing_fraction", "rnp_rna_processing_marker_fraction"))

  vascular_ecm <- max_frac(bm, bbb, mural, adhesion)
  if (is.finite(vascular_ecm) && vascular_ecm >= classification_threshold && max_frac(bm, adhesion) >= max_frac(bbb, mural, 0)) return("perivascular basement membrane / ECM")
  if (is.finite(max_frac(bbb, mural)) && max_frac(bbb, mural) >= classification_threshold) return("BBB / vascular mural")
  if (targeted_claim > 0 || targeted_micro > 1) return("targeted microglia signature overlap")
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
  if (is.null(go) || !nrow(go)) return(list(top_GO_BP_labels = NA_character_, top_GO_MF_labels = NA_character_, top_GO_CC_labels = NA_character_))
  color_col <- first_present_col(go, c("ModuleColor", "module", "ModuleID"))
  desc_col <- first_present_col(go, c("Description", "description", "term_description", "ID"))
  ont_col <- first_present_col(go, c("ONTOLOGY", "ontology"))
  padj_col <- first_present_col(go, c("p.adjust", "p_adj", "padj", "qvalue"))
  if (is.na(color_col) || is.na(desc_col)) return(list(top_GO_BP_labels = NA_character_, top_GO_MF_labels = NA_character_, top_GO_CC_labels = NA_character_))
  tab <- go |> dplyr::filter(as.character(.data[[color_col]]) %in% c(as.character(module_color), paste0("ME", module_color)))
  if (nrow(tab) && !is.na(padj_col)) tab <- tab |> dplyr::arrange(.data[[padj_col]])
  one <- function(ont) {
    x <- if (!is.na(ont_col)) tab |> dplyr::filter(toupper(as.character(.data[[ont_col]])) == ont) else tab
    paste(utils::head(unique(as.character(x[[desc_col]])), 5), collapse = ";")
  }
  list(top_GO_BP_labels = one("BP"), top_GO_MF_labels = one("MF"), top_GO_CC_labels = one("CC"))
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

  robust <- as.integer(row$n_microglia_robust_term_overlaps %||% 0)
  sensitive <- as.integer(row$n_neuropil_sensitive_term_overlaps %||% 0)
  mixed <- as.integer(row$n_mixed_microenvironment_term_overlaps %||% 0)
  targeted_micro <- as.integer(row$n_targeted_microglia_signature_overlaps %||% 0)
  targeted_claim <- as.integer(row$n_targeted_claim_ready_signature_overlaps %||% 0)
  targeted_shared <- as.integer(row$n_targeted_neuropil_shared_signature_overlaps %||% 0)

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
  robust <- as.integer(row$n_microglia_robust_term_overlaps %||% 0)
  sensitive <- as.integer(row$n_neuropil_sensitive_term_overlaps %||% 0)
  mixed <- as.integer(row$n_mixed_microenvironment_term_overlaps %||% 0)
  targeted_micro <- as.integer(row$n_targeted_microglia_signature_overlaps %||% 0)
  targeted_claim <- as.integer(row$n_targeted_claim_ready_signature_overlaps %||% 0)
  targeted_shared <- as.integer(row$n_targeted_neuropil_shared_signature_overlaps %||% 0)
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

module_targeted_signature_counts <- function(module_genes) {
  empty <- data.frame(
    n_targeted_microglia_signature_overlaps = 0L,
    n_targeted_claim_ready_signature_overlaps = 0L,
    n_targeted_neuropil_shared_signature_overlaps = 0L,
    n_targeted_signature_overlap_proteins = 0L,
    best_targeted_microglia_signatures = NA_character_,
    best_targeted_signature_overlap_proteins = NA_character_,
    stringsAsFactors = FALSE
  )
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
fraction_cols <- grep("(_marker_fraction|_fraction)$", names(module_annot), value = TRUE)
module_annot[fraction_cols] <- lapply(module_annot[fraction_cols], function(x) suppressWarnings(as.numeric(x)))

if (DATASET == "microglia" || force_microglia) {
  ref_counts <- dplyr::bind_rows(lapply(module_annot$ModuleColor, module_neuropil_reference_counts))
  targeted_counts <- dplyr::bind_rows(lapply(module_rows$proteins, module_targeted_signature_counts))
  module_annot <- module_annot |>
    dplyr::bind_cols(ref_counts, targeted_counts)
} else {
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
      best_targeted_signature_overlap_proteins = NA_character_
    )
}

evidence_df <- dplyr::bind_rows(lapply(seq_len(nrow(module_annot)), function(i) as.data.frame(classify_rationale(module_annot[i, , drop = FALSE]), stringsAsFactors = FALSE)))
module_annot <- dplyr::bind_cols(module_annot, evidence_df)
module_annot$classification_threshold <- classification_threshold
module_annot$microenvironment_class <- vapply(seq_len(nrow(module_annot)), function(i) classify_module(module_annot[i, , drop = FALSE]), character(1))
module_annot$microenvironment_label <- vapply(seq_len(nrow(module_annot)), function(i) label_from_evidence(module_annot[i, , drop = FALSE]), character(1))
module_annot$microenvironment_confidence <- vapply(seq_len(nrow(module_annot)), function(i) confidence_from_evidence(module_annot[i, , drop = FALSE]), character(1))
module_annot$module_display_label <- paste0(module_annot$ModuleID, " — ", module_annot$microenvironment_label)
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
  super_annot <- data.frame(dataset = DATASET, SupermoduleID = NA_character_, Supermodule_FinalLabel = NA_character_, n_member_modules = 0L, dominant_microenvironment_class = "missing_supermodule_annotation", interpretation_note = WGCNA_ROI_NOTE)
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
      dominant_GO_terms = paste(utils::head(unique(c(split_tokens(.data$top_GO_BP_labels), split_tokens(.data$top_GO_MF_labels), split_tokens(.data$top_GO_CC_labels))), 10), collapse = ";"),
      top_hub_proteins = paste(utils::head(unique(split_tokens(.data$top_hub_proteins)), 30), collapse = ";"),
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
      marker_registry_version = dplyr::first(.data$marker_registry_version),
      empirical_marker_set_version = dplyr::first(.data$empirical_marker_set_version),
      interpretation_note = WGCNA_ROI_NOTE,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Supermodule_FinalLabel_Original = .data$Supermodule_FinalLabel,
      Supermodule_FinalLabel = dplyr::case_when(
        generic_supermodule_label(.data$Supermodule_FinalLabel) & !is.na(.data$dominant_module_labels) & nzchar(.data$dominant_module_labels) ~ .data$dominant_module_labels,
        TRUE ~ .data$Supermodule_FinalLabel
      ),
      Supermodule_LabelRationale = paste0(
        dplyr::coalesce(.data$Supermodule_LabelRationale, ""),
        ifelse(generic_supermodule_label(.data$Supermodule_FinalLabel_Original), "; fallback_label_from_member_module_microenvironment_labels", "")
      ),
      ManualReviewRequired = dplyr::case_when(
        generic_supermodule_label(.data$Supermodule_FinalLabel_Original) ~ TRUE,
        .data$dominant_microenvironment_class %in% c("ambiguous_or_mixed", "shared_microenvironment") ~ TRUE,
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
  super_annot <- super_annot |>
    dplyr::mutate(
      microenvironment_label_component = supermodule_microenvironment_label(.data$dominant_microenvironment_class, dataset = DATASET),
      has_active_manual_label = .data$manual_label_status == "manual_label_present_for_dataset_module" |
        (!is.na(.data$Supermodule_CuratedLabel) & nzchar(as.character(.data$Supermodule_CuratedLabel))),
      display_label_component = dplyr::case_when(
        .data$has_active_manual_label ~ shorten_supermodule_label(.data$Supermodule_FinalLabel, max_chars = 30),
        !is.na(.data$microenvironment_label_component) & nzchar(.data$microenvironment_label_component) ~ .data$microenvironment_label_component,
        !is.na(.data$Macroprogram_Display) & nzchar(.data$Macroprogram_Display) & .data$Macroprogram_Display != "Unresolved / mixed" ~ .data$Macroprogram_Display,
        TRUE ~ shorten_supermodule_label(.data$Supermodule_FinalLabel, max_chars = 30)
      ),
      Supermodule_DisplayLabel = compose_supermodule_display_label(.data$SupermoduleID, .data$display_label_component),
      Supermodule_ShortLabel = .data$Supermodule_DisplayLabel,
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
    dplyr::select(-dplyr::any_of(c("has_active_manual_label", "display_label_component", "has_coherent_hubs")))
}

write_table_and_source(module_annot, PATHS$tables, PATHS$source_data, "WGCNA_module_biological_annotation.csv")
write_table_and_source(super_annot, PATHS$tables, PATHS$source_data, "WGCNA_supermodule_biological_annotation.csv")
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
      "Supermodule_LabelSource", "Supermodule_LabelConfidence", "is_singleton_supermodule",
      "n_member_modules", "GO_label_confidence_class", "dominant_microenvironment_class",
      "top_GO_BP_terms", "top_hub_symbols", "label_rationale", "manual_review_required"
    ))) |>
    dplyr::rename(Supermodule_DataDriven = "SupermoduleID")
  write_table_and_source(display_audit, PATHS$tables, PATHS$source_data, "WGCNA_supermodule_display_label_audit.csv")
}
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(modules = module_annot, supermodules = super_annot), file.path(PATHS$tables, "WGCNA_module_microenvironment_annotation.xlsx"))
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
      fraction = suppressWarnings(as.numeric(.data$fraction))
    ) |>
    dplyr::filter(is.finite(.data$fraction), .data$fraction > 0)
  write_table_and_source(comp, PATHS$tables, PATHS$source_data, "WGCNA_supermodule_microenvironment_composition_source.csv")
  if (nrow(comp)) {
    fig_h <- max(90, 14 * length(unique(comp$plot_label)) + 35)
    p2 <- ggplot2::ggplot(comp, ggplot2::aes(x = .data$fraction, y = stats::reorder(.data$plot_label, .data$SupermoduleID), fill = .data$class)) +
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
      any(module_annot$n_targeted_microglia_signature_overlaps > 0, na.rm = TRUE)) {
    targeted_plot <- module_annot |>
      dplyr::mutate(
        total_targeted = suppressWarnings(as.numeric(.data$n_targeted_microglia_signature_overlaps)),
        claim_ready = suppressWarnings(as.numeric(.data$n_targeted_claim_ready_signature_overlaps)),
        shared = suppressWarnings(as.numeric(.data$n_targeted_neuropil_shared_signature_overlaps)),
        module_display_label = dplyr::coalesce(.data$module_display_label, .data$ModuleID)
      ) |>
      dplyr::filter(.data$total_targeted > 0 | .data$claim_ready > 0 | .data$shared > 0) |>
      dplyr::arrange(dplyr::desc(.data$claim_ready), dplyr::desc(.data$total_targeted)) |>
      dplyr::mutate(module_display_label = stats::reorder(.data$module_display_label, .data$total_targeted))
    targeted_long <- targeted_plot |>
      dplyr::select("module_display_label", total_targeted, claim_ready, shared) |>
      tidyr::pivot_longer(cols = -dplyr::all_of("module_display_label"), names_to = "evidence", values_to = "n")
    write_table_and_source(targeted_long, PATHS$tables, PATHS$source_data, "WGCNA_module_targeted_microglia_signature_overlap_source.csv")
    p4 <- ggplot2::ggplot(targeted_long, ggplot2::aes(x = .data$n, y = .data$module_display_label, fill = .data$evidence)) +
      ggplot2::geom_col(position = "dodge", width = 0.68) +
      ggplot2::scale_fill_manual(
        values = c(total_targeted = "#3C5488", claim_ready = "#00A087", shared = "#E64B35"),
        labels = c(total_targeted = "microglia signature", claim_ready = "claim-ready", shared = "neuropil-shared"),
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
  inputs = c(FILES, targeted_microglia_signature_enrichment = targeted_signature_file),
  outputs = list(tables = PATHS$tables, source_data = PATHS$source_data, figures = PATHS$figures),
  parameters = list(
    dataset = DATASET,
    force_microglia_annotation = force_microglia,
    classification_threshold = classification_threshold,
    marker_registry_version = marker_registry_version,
    empirical_marker_set_version = empirical_marker_set_version
  ),
  notes = paste(WGCNA_ROI_NOTE, "Improved annotation separates vascular basement membrane/BBB/mural, astrocyte/endfoot, oligodendrocyte/myelin, neuropil, and microglia-supported ROI evidence; labels are annotation only.")
)

message("WGCNA module/supermodule biological annotation complete for dataset: ", DATASET)
