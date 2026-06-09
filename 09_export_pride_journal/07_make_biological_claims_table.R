#!/usr/bin/env Rscript
# Manuscript/journal biological claims index — not a PRIDE-required deposition artifact.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "schema_validation.R"))

required_pkgs <- c("dplyr", "readr", "tibble", "stringr", "tidyr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

claim_columns <- c(
  "claim_id", "dataset", "region", "layer_cell_compartment", "contrast",
  "biological_program", "direction", "key_proteins_genes", "evidence_type",
  "effect_size_NES", "raw_p", "FDR", "robustness_stability_metric",
  "source_file", "figure_table_target", "interpretation_note",
  "claim_grade", "primary_evidence", "orthogonal_support", "major_limitation",
  "safe_interpretation", "unsafe_overinterpretation"
)

numeric_claim_columns <- c("effect_size_NES", "raw_p", "FDR")
character_claim_columns <- setdiff(claim_columns, numeric_claim_columns)

empty_claims <- function() {
  out <- as.data.frame(setNames(rep(list(character()), length(claim_columns)), claim_columns), stringsAsFactors = FALSE)
  for (col in numeric_claim_columns) out[[col]] <- numeric()
  out
}

standardize_claims <- function(df) {
  if (is.null(df) || !nrow(df)) return(empty_claims())
  for (col in claim_columns) if (!col %in% names(df)) df[[col]] <- NA
  df <- df[, claim_columns, drop = FALSE]
  for (col in numeric_claim_columns) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  for (col in character_claim_columns) df[[col]] <- as.character(df[[col]])
  df
}

grade_claim <- function(evidence_type, fdr, interpretation_note) {
  evidence_type <- as.character(evidence_type)
  note <- tolower(as.character(interpretation_note))
  has_da <- grepl("de|differential|overlap", tolower(evidence_type))
  has_enrichment <- grepl("gsea|enrichment|program", tolower(evidence_type))
  has_module <- grepl("wgcna|module|network", tolower(evidence_type))
  fdr <- suppressWarnings(as.numeric(fdr))
  if (grepl("confound|unsafe|ambiguous", note)) return("X")
  if (isTRUE(has_da && has_enrichment && has_module)) return("A")
  if (isTRUE(has_da && has_enrichment)) return("B")
  if (isTRUE(has_enrichment || has_module)) return("C")
  if (!is.na(fdr) && fdr <= 0.05) return("C")
  "D"
}

primary_evidence_label <- function(evidence_type) {
  dplyr::case_when(
    grepl("WGCNA_DE_GSEA_overlap", evidence_type, ignore.case = TRUE) ~ "DA/GSEA overlap with WGCNA module evidence",
    grepl("GSEA|program|enrichment", evidence_type, ignore.case = TRUE) ~ "Enrichment/program summary",
    grepl("WGCNA|module", evidence_type, ignore.case = TRUE) ~ "Module/network evidence",
    grepl("behavior", evidence_type, ignore.case = TRUE) ~ "Behavior/network association",
    TRUE ~ evidence_type
  )
}

supermodule_annotation_for_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "module_annotation", dataset, "WGCNA_supermodule_biological_annotation.csv")
  ann <- read_csv_if_exists(f)
  if (is.null(ann) || !nrow(ann)) return(NULL)
  for (col in c("dataset", "SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "dominant_microenvironment_class")) {
    if (!col %in% names(ann)) ann[[col]] <- NA
  }
  ann %>%
    dplyr::mutate(
      dataset = dataset,
      Supermodule_DisplayLabel = dplyr::coalesce(
        as.character(.data$Supermodule_DisplayLabel),
        as.character(.data$Supermodule_FinalLabel),
        as.character(.data$Macroprogram_Display),
        as.character(.data$SupermoduleID)
      ),
      Supermodule_DisplayLabel_annotation = .data$Supermodule_DisplayLabel,
      Supermodule_FinalLabel_annotation = .data$Supermodule_FinalLabel,
      Macroprogram_Display_annotation = .data$Macroprogram_Display
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(c("SupermoduleID", "Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display")),
      names_to = "supermodule_claim_key_source",
      values_to = "supermodule_claim_key"
    ) %>%
    dplyr::filter(!is.na(.data$supermodule_claim_key), nzchar(.data$supermodule_claim_key)) %>%
    dplyr::distinct(.data$dataset, .data$supermodule_claim_key, .keep_all = TRUE) %>%
    dplyr::transmute(
      dataset = .data$dataset,
      supermodule_claim_key = .data$supermodule_claim_key,
      Supermodule_DisplayLabel = .data$Supermodule_DisplayLabel_annotation,
      Supermodule_FinalLabel = .data$Supermodule_FinalLabel_annotation,
      Macroprogram_Display = .data$Macroprogram_Display_annotation,
      dominant_microenvironment_class = .data$dominant_microenvironment_class
    )
}

latest_csv <- function(root, pattern) latest_file(root, pattern)

collect_program_claims <- function(dataset) {
  f <- path_results("tables", "04_differential_expression_enrichment", "biological_program_summary", dataset, "program_summary.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"biological_program" %in% names(df)) return(empty_claims())
  for (col in c("route_category", "route_unit", "min_raw_p", "min_fdr", "representative_NES", "key_genes", "core_genes", "n_terms", "top_term")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df %>%
    dplyr::transmute(
      dataset = dataset,
      region = .data$route_unit,
      layer_cell_compartment = .data$route_category,
      contrast = .data$comparison,
      direction = .data$direction,
      biological_program = .data$biological_program,
      key_proteins_genes = dplyr::coalesce(.data$key_genes, .data$core_genes),
      evidence_type = "GSEA_program_summary",
      effect_size_NES = .data$representative_NES,
      raw_p = .data$min_raw_p,
      FDR = .data$min_fdr,
      robustness_stability_metric = paste0("n_terms=", .data$n_terms),
      source_file = f,
      figure_table_target = "program_atlas_heatmap; manuscript program summary table",
      interpretation_note = paste0(
        vapply(seq_along(.data$min_fdr), function(i) interpretation_strength(fdr = .data$min_fdr[[i]], effect_size = .data$representative_NES[[i]], n = NA), character(1)),
        "; top term: ", dplyr::coalesce(.data$top_term, "NA"),
        "; regex-based thematic synthesis"
      )
    ) %>%
    standardize_claims()
}

collect_wgcna_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_evidence_rank.csv")
  if (!file.exists(f)) f <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_priority_summary.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"ModuleID" %in% names(df)) return(empty_claims())
  for (col in c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "SupermoduleID", "Supermodule", "ModuleLabel_Final", "strongest_condition_contrast", "condition_model_fdr",
                "condition_model_p", "strongest_condition_adjusted_delta", "evidence_warning")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  super_ann <- supermodule_annotation_for_claims(dataset)
  df <- df %>%
    dplyr::mutate(
      dataset = dataset,
      supermodule_claim_key = dplyr::coalesce(as.character(.data$SupermoduleID), as.character(.data$Supermodule), as.character(.data$Supermodule_FinalLabel), as.character(.data$Macroprogram_Display))
    )
  if (!is.null(super_ann)) {
    df <- df %>%
      dplyr::left_join(super_ann, by = c("dataset", "supermodule_claim_key"), suffix = c("", ".annotation")) %>%
      dplyr::mutate(
        Supermodule_DisplayLabel = dplyr::coalesce(.data$Supermodule_DisplayLabel.annotation, .data$Supermodule_DisplayLabel),
        Supermodule_FinalLabel = dplyr::coalesce(.data$Supermodule_FinalLabel.annotation, .data$Supermodule_FinalLabel),
        Macroprogram_Display = dplyr::coalesce(.data$Macroprogram_Display.annotation, .data$Macroprogram_Display)
      )
  }
  df %>%
    dplyr::transmute(
      dataset = dataset,
      contrast = .data$strongest_condition_contrast,
      biological_program = dplyr::coalesce(.data$Supermodule_DisplayLabel, .data$Supermodule_FinalLabel, .data$Macroprogram_Display, .data$Supermodule, .data$ModuleLabel_Final),
      direction = dplyr::case_when(
        suppressWarnings(as.numeric(.data$strongest_condition_adjusted_delta)) > 0 ~ "positive_effect",
        suppressWarnings(as.numeric(.data$strongest_condition_adjusted_delta)) < 0 ~ "negative_effect",
        TRUE ~ NA_character_
      ),
      evidence_type = "WGCNA_module",
      effect_size_NES = .data$strongest_condition_adjusted_delta,
      raw_p = .data$condition_model_p,
      FDR = .data$condition_model_fdr,
      robustness_stability_metric = if ("preservation_Zsummary_median" %in% names(df)) .data$preservation_Zsummary_median else NA,
      key_proteins_genes = paste(.data$ModuleID, .data$ModuleLabel_Final, sep = ": "),
      source_file = f,
      figure_table_target = "WGCNA_module_priority; WGCNA module evidence table",
      interpretation_note = paste0(
        vapply(seq_along(.data$condition_model_fdr), function(i) interpretation_strength(fdr = .data$condition_model_fdr[[i]], effect_size = .data$strongest_condition_adjusted_delta[[i]], n = NA), character(1)),
        "; ", dplyr::coalesce(.data$evidence_warning, "Low-n module evidence; prioritize replicated or convergent modules.")
      )
    ) %>%
    standardize_claims()
}

collect_overlap_claims <- function(dataset) {
  f <- path_results("tables", "06_modules_WGCNA", "05_wgcna_de_gsea_overlap", dataset, "WGCNA_vs_DE_GSEA_overlap.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"ModuleID" %in% names(df)) return(empty_claims())
  if ("status" %in% names(df)) return(empty_claims())
  for (col in c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Macroprogram_Display", "SupermoduleID", "Supermodule", "jaccard_DE", "n_DE_overlap", "top_overlap_proteins")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  super_ann <- supermodule_annotation_for_claims(dataset)
  df <- df %>%
    dplyr::mutate(
      dataset = dataset,
      supermodule_claim_key = dplyr::coalesce(as.character(.data$SupermoduleID), as.character(.data$Supermodule), as.character(.data$Supermodule_FinalLabel), as.character(.data$Macroprogram_Display))
    )
  if (!is.null(super_ann)) {
    df <- df %>%
      dplyr::left_join(super_ann, by = c("dataset", "supermodule_claim_key"), suffix = c("", ".annotation")) %>%
      dplyr::mutate(
        Supermodule_DisplayLabel = dplyr::coalesce(.data$Supermodule_DisplayLabel.annotation, .data$Supermodule_DisplayLabel),
        Supermodule_FinalLabel = dplyr::coalesce(.data$Supermodule_FinalLabel.annotation, .data$Supermodule_FinalLabel),
        Macroprogram_Display = dplyr::coalesce(.data$Macroprogram_Display.annotation, .data$Macroprogram_Display)
      )
  }
  df %>%
    dplyr::arrange(.data$fisher_fdr) %>%
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      dataset = dataset,
      contrast = .data$contrast,
      biological_program = dplyr::coalesce(.data$Supermodule_DisplayLabel, .data$Supermodule_FinalLabel, .data$Macroprogram_Display, .data$Supermodule),
      direction = NA_character_,
      evidence_type = "WGCNA_DE_GSEA_overlap",
      effect_size_NES = if ("jaccard_DE" %in% names(df)) .data$jaccard_DE else NA,
      raw_p = .data$fisher_p,
      FDR = .data$fisher_fdr,
      robustness_stability_metric = paste0("n_DE_overlap=", .data$n_DE_overlap),
      key_proteins_genes = .data$top_overlap_proteins,
      source_file = f,
      figure_table_target = "WGCNA_vs_DE_GSEA_overlap",
      interpretation_note = paste0(
        vapply(seq_along(.data$fisher_fdr), function(i) interpretation_strength(fdr = .data$fisher_fdr[[i]], effect_size = .data$jaccard_DE[[i]], n = .data$n_DE_overlap[[i]]), character(1)),
        "; overlap supports convergence, not causality."
      )
    ) %>%
    standardize_claims()
}

collect_behavior_claims <- function() {
  f <- path_results("tables", "08_behavior_physio_coupling", "network_behavior_coupling", "edge_behavior_figure_ready_table.csv")
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df)) return(empty_claims())
  df %>%
    dplyr::transmute(
      dataset = current_dataset(),
      contrast = .data$Edge,
      biological_program = .data$Outcome,
      direction = dplyr::case_when(.data$estimate > 0 ~ "positive_correlation", .data$estimate < 0 ~ "negative_correlation", TRUE ~ "neutral"),
      evidence_type = "network_behavior_coupling",
      effect_size_NES = .data$estimate,
      raw_p = .data$p.value,
      FDR = .data$fdr,
      robustness_stability_metric = paste0("n=", .data$n),
      source_file = f,
      figure_table_target = "edge_behavior_correlation_forest",
      interpretation_note = paste(.data$interpretation_strength, .data$limitations, sep = "; ")
    ) %>%
    standardize_claims()
}

collect_microglia_signature_claims <- function(dataset) {
  f <- path_results(
    "tables",
    "04_differential_expression_enrichment",
    "microglia_targeted_signature_enrichment",
    dataset,
    "microglia_signature_claims_ready.csv"
  )
  df <- read_csv_if_exists(f)
  if (is.null(df) || !nrow(df) || !"signature" %in% names(df)) return(empty_claims())
  for (col in c("comparison", "left_region", "right_region", "left_unit", "left_condition", "right_condition", "claim_type", "NES", "pvalue", "padj", "matched_genes", "contrast_class")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df %>%
    dplyr::transmute(
      dataset = dataset,
      region = ifelse(!is.na(.data$left_region) & !is.na(.data$right_region), paste(.data$left_region, .data$right_region, sep = "_vs_"), .data$left_region),
      layer_cell_compartment = .data$left_unit,
      contrast = .data$comparison,
      direction = dplyr::case_when(.data$NES > 0 ~ "positive_NES", .data$NES < 0 ~ "negative_NES", TRUE ~ "neutral"),
      biological_program = .data$signature,
      key_proteins_genes = .data$matched_genes,
      evidence_type = "microglia_signature_enrichment",
      effect_size_NES = .data$NES,
      raw_p = .data$pvalue,
      FDR = .data$padj,
      robustness_stability_metric = paste0("contrast_class=", .data$contrast_class, "; claim_type=", .data$claim_type),
      source_file = f,
      figure_table_target = "microglia_signature_claims_ready",
      interpretation_note = dplyr::case_when(
        .data$claim_type == "within_region_condition_microglia_program" ~ "Within-region condition contrast; conservative microglia-supported claim.",
        .data$claim_type == "regional_microglia_program" ~ "Cross-region same-condition contrast; conservative regional microglia program claim.",
        TRUE ~ "Microglia signature enrichment claim-ready row."
      )
    ) %>%
    standardize_claims()
}

if (is_dry_run()) {
  dry_run_line("Script", "09_export_pride_journal/07_make_biological_claims_table.R")
  dry_run_line("Datasets", paste(valid_datasets(), collapse = ", "))
  dry_run_line("Output CSV", path_results("tables", "biological_claims_table.csv"))
  dry_run_line("Output XLSX", path_results("tables", "biological_claims_table.xlsx"))
  quit(status = 0, save = "no")
}

claims <- dplyr::bind_rows(
  lapply(valid_datasets(), collect_program_claims),
  lapply(valid_datasets(), collect_wgcna_claims),
  lapply(valid_datasets(), collect_overlap_claims),
  lapply(valid_datasets(), collect_microglia_signature_claims),
  list(collect_behavior_claims())
) %>%
  standardize_claims() %>%
  dplyr::mutate(
    claim_id = sprintf("CLAIM_%04d", dplyr::row_number()),
    interpretation_note = dplyr::coalesce(.data$interpretation_note, "No interpretation note available; review source file before manuscript use."),
    claim_grade = vapply(seq_along(.data$evidence_type), function(i) grade_claim(.data$evidence_type[[i]], .data$FDR[[i]], .data$interpretation_note[[i]]), character(1)),
    primary_evidence = primary_evidence_label(.data$evidence_type),
    orthogonal_support = dplyr::case_when(
      grepl("overlap", .data$evidence_type, ignore.case = TRUE) ~ "WGCNA module and DA/GSEA convergence",
      grepl("microglia_signature", .data$evidence_type, ignore.case = TRUE) ~ "Microglia-targeted signature evidence",
      grepl("WGCNA_module", .data$evidence_type, ignore.case = TRUE) ~ "Module evidence; seek DA/enrichment support before strong biological claims",
      TRUE ~ "Review companion DA, enrichment, module, and QC outputs"
    ),
    major_limitation = dplyr::case_when(
      .data$dataset == "microglia" ~ "Microglia dataset is microglia-enriched ROI/local microenvironment, not purified microglia; region-only interpretation.",
      .data$claim_grade %in% c("C", "D") ~ "Single evidence stream or exploratory evidence; avoid causal or cell-intrinsic claims.",
      TRUE ~ "Observational proteomics; causal direction is not established."
    ),
    safe_interpretation = dplyr::case_when(
      .data$dataset == "microglia" ~ paste0("In microglia-enriched ROIs, evidence supports a local microenvironment-associated ", .data$biological_program, " signal."),
      TRUE ~ paste0("Evidence supports an association between ", .data$biological_program, " and the specified dataset/contrast.")
    ),
    unsafe_overinterpretation = dplyr::case_when(
      .data$dataset == "microglia" ~ "Do not claim purified microglial cell-intrinsic regulation or subtractive neuropil correction from these data alone.",
      TRUE ~ "Do not claim causality, mechanism, or cell-type specificity without independent validation."
    )
  )

validate_table_schema(claims, "biological_claims_table", strict = TRUE)

dir_create(path_results("tables"))
csv_out <- path_results("tables", "biological_claims_table.csv")
xlsx_out <- path_results("tables", "biological_claims_table.xlsx")
readr::write_csv(claims, csv_out, na = "")
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(biological_claims = claims), xlsx_out)
}

message("Biological claims table written: ", csv_out)

write_run_manifest(
  path_results("logs", "09_export_pride_journal", "biological_claims_table", "run_manifest.yml"),
  inputs = list(source_files = unique(claims$source_file)),
  outputs = list(csv = csv_out, xlsx = if (file.exists(xlsx_out)) xlsx_out else NA_character_),
  parameters = list(datasets = valid_datasets(), schema = "biological_claims_table"),
  notes = "Evidence-graded manuscript claim table. Missing statistics remain NA."
)
