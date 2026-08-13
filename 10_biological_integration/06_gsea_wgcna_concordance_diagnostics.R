#!/usr/bin/env Rscript
# Script: 10_biological_integration/06_gsea_wgcna_concordance_diagnostics.R
# Stage: additive diagnostic review of 05_gsea_wgcna_concordance.R
#
# The strict concordance classifications and adaptive pattern summary are read
# unchanged. This script explains their gates and writes a separate descriptive,
# hypothesis-generating direction layer from complete canonical GSEA effects.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "wgcna_group_effect_consumer_utils.R"))
source(repo_path("R", "manuscript_go_theme_utils.R"))
source(repo_path("R", "gsea_wgcna_concordance_utils.R"))
source(repo_path("R", "gsea_wgcna_concordance_diagnostic_utils.R"))

SCRIPT_ID <- "10_biological_integration/06_gsea_wgcna_concordance_diagnostics.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
run <- integration_cli(default_dataset = "all", allow_all = TRUE)
datasets <- integration_datasets(run$dataset)
if (!setequal(datasets, valid_datasets())) {
  stop("The diagnostic is a global three-dataset integration and requires --dataset all.",
       call. = FALSE)
}
paths <- integration_paths("gsea_wgcna_concordance_diagnostics", "global")

required_packages <- c(
  "dplyr", "tidyr", "tidyselect", "readr", "stringr", "tibble"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "),
       ".", call. = FALSE)
}

strict_dir <- path_results(
  "tables", "10_biological_integration", "gsea_wgcna_concordance", "global"
)
strict_long_file <- file.path(strict_dir, "gsea_wgcna_concordance_long.csv")
strict_pattern_file <- file.path(
  strict_dir, "adaptive_resilience_pattern_summary.csv"
)
strict_overlap_file <- file.path(
  strict_dir, "program_specific_leading_edge_module_overlap.csv"
)
strict_matching_file <- file.path(strict_dir, "program_module_matching_audit.csv")
gsea_terms_file <- file.path(
  strict_dir, "ontology_aware_gsea_theme_assignments_all_contrasts.csv"
)
theme_mapping_file <- repo_path(
  "config", "gsea_wgcna_theme_module_mapping.csv"
)
theme_registry_file <- repo_path("config", "manuscript_go_theme_registry.tsv")

dataset_inputs <- lapply(datasets, function(dataset) {
  list(
    wgcna_handoff = path_results(
      "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
      "WGCNA_inferential_handoff.csv"
    ),
    feature_universe = path_results(
      "tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules",
      "WGCNA_feature_universe.csv"
    ),
    module_annotation = path_results(
      "tables", "06_modules_WGCNA", "module_annotation", dataset,
      "WGCNA_module_biological_annotation.csv"
    ),
    supermodule_annotation = path_results(
      "tables", "06_modules_WGCNA", "module_annotation", dataset,
      "WGCNA_supermodule_biological_annotation.csv"
    ),
    module_go = path_results(
      "tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules",
      "WGCNA_module_GO_enrichment_long.csv"
    ),
    final_labels = path_results(
      "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
      "WGCNA_final_label_lookup.csv"
    ),
    group_sensitivity = path_results(
      "tables", "06_modules_WGCNA", "group_effects", dataset,
      "WGCNA_group_effect_sensitivity.csv"
    )
  )
})
names(dataset_inputs) <- datasets

required_inputs <- c(
  strict_long = strict_long_file,
  strict_patterns = strict_pattern_file,
  strict_overlap = strict_overlap_file,
  strict_matching = strict_matching_file,
  ontology_gsea_terms = gsea_terms_file,
  theme_entity_mapping = theme_mapping_file,
  manuscript_theme_registry = theme_registry_file,
  unlist(lapply(names(dataset_inputs), function(dataset) {
    setNames(
      unlist(dataset_inputs[[dataset]]),
      paste0(dataset, "_", names(dataset_inputs[[dataset]]))
    )
  }))
)

output_files <- c(
  "concordant_imprecise_diagnostic.csv",
  "concordance_wgcna_power_diagnostic.csv",
  "unmatched_gsea_program_mapping_audit.csv",
  "unmatched_gsea_theme_mapping_audit.csv",
  "adaptive_resilience_gate_audit.csv",
  "adaptive_resilience_directional_patterns.csv",
  "adaptive_resilience_local_spatial_patterns.csv",
  "neuropil_RES_CON_SUS_pattern_audit.csv",
  "diagnostic_summary.csv"
)

if (run$dry_run) {
  dry_run_inputs(SCRIPT_ID, required_inputs)
  for (filename in output_files) {
    dry_run_line("Output table", file.path(paths$tables, filename))
    dry_run_line("Source-data copy", file.path(paths$source_data, filename))
  }
  dry_run_line("Diagnostic report", file.path(paths$reports, "README.md"))
  quit(save = "no", status = 0L)
}

missing_required <- required_inputs[!file.exists(required_inputs)]
if (length(missing_required)) {
  stop("Required diagnostic input(s) missing: ", paste(
    names(missing_required), missing_required, sep = "=", collapse = "; "
  ), ".", call. = FALSE)
}

source_hashes <- function(paths_named) {
  data.frame(
    input_name = names(paths_named),
    source_file = normalizePath(
      unname(paths_named), winslash = "/", mustWork = TRUE
    ),
    source_md5 = as.character(tools::md5sum(unname(paths_named))),
    stringsAsFactors = FALSE
  )
}
hash_before <- source_hashes(required_inputs)

strict_long <- readr::read_csv(
  strict_long_file, show_col_types = FALSE, progress = FALSE, guess_max = Inf
)
strict_patterns <- readr::read_csv(
  strict_pattern_file, show_col_types = FALSE, progress = FALSE, guess_max = Inf
)
strict_overlap <- readr::read_csv(
  strict_overlap_file, show_col_types = FALSE, progress = FALSE, guess_max = Inf
)
strict_matching <- readr::read_csv(
  strict_matching_file, show_col_types = FALSE, progress = FALSE, guess_max = Inf
)
gsea_terms <- readr::read_csv(
  gsea_terms_file, show_col_types = FALSE, progress = FALSE, guess_max = Inf
)
gww_assert_columns(
  gsea_terms,
  c(
    "dataset", "phenotype_contrast", "spatial_unit", "GO_ID",
    "GO_description", "NES", "raw_p", "GSEA_FDR", "theme_id",
    "manuscript_theme", "theme_role", "theme_claim_eligible",
    "theme_assignment_id"
  ),
  "Ontology-aware all-contrast GSEA theme assignments"
)

handoff_by_dataset <- lapply(datasets, function(dataset) {
  wgcna_inferential_handoff_read(
    dataset_inputs[[dataset]]$wgcna_handoff,
    paste(dataset, "WGCNA_inferential_handoff.csv")
  )
})
names(handoff_by_dataset) <- datasets
handoff <- dplyr::bind_rows(handoff_by_dataset)
wgcna_inferential_handoff_validate(handoff, "Combined canonical WGCNA handoff")
if (any(!handoff$independent_hypothesis %in% TRUE) ||
    any(handoff$claim_entity_role == "compatibility_alias")) {
  stop("Canonical WGCNA handoff contains an alias/non-independent row.",
       call. = FALSE)
}

theme_registry <- read_manuscript_go_theme_registry(theme_registry_file)
theme_mapping <- gww_read_theme_entity_mapping(theme_mapping_file)
validated_theme_mapping <- gww_validate_theme_entity_mapping(
  handoff, theme_mapping, theme_registry
)
program_matches <- gww_build_entity_theme_matches(
  handoff, theme_mapping, theme_registry
)
effect_lookup <- gwwd_effect_distribution_lookup(handoff)
imprecise <- gwwd_concordant_imprecise(strict_long, effect_lookup)
power <- gwwd_power_diagnostic(strict_long)

# Reconstruct ontology-theme exclusions before making manual-review suggestions.
strict_theme_terms <- gsea_terms |>
  dplyr::transmute(
    dataset, phenotype_contrast, spatial_unit,
    Comparison = .data$source_comparison,
    ID = .data$GO_ID, Description = .data$GO_description,
    GO_ID, GO_description, NES,
    pvalue = .data$raw_p, p.adjust = .data$GSEA_FDR,
    core_enrichment = .data$leading_edge_proteins,
    core_enrichment_gene = .data$leading_edge_genes,
    evidence_source_family,
    theme_id, manuscript_theme, theme_role, theme_claim_eligible,
    anchor_GO_IDs = .data$anchor_GO_ID,
    mapping_method = .data$mapping_type,
    registry_version, theme_assignment_id,
    source_supplementary_file = NA_character_
  )
gww_validate_theme_terms(strict_theme_terms)
strict_local_gsea <- gww_build_local_gsea_evidence(strict_theme_terms)
strict_recurrent_cross_spatial_gsea <- gww_build_recurrent_cross_spatial_gsea_evidence(
  strict_theme_terms, strict_local_gsea, min_recurrent_units = 2L
)
strict_gsea_evidence <- dplyr::bind_rows(
  strict_local_gsea, strict_recurrent_cross_spatial_gsea
)
# QC-review themes remain visible for audit and manual mapping review, but are
# deliberately excluded from the claim-eligible concordance evidence above.
qc_context_evidence <- strict_theme_terms |>
  dplyr::filter(
    .data$theme_role == "qc_review",
    is.finite(.data$p.adjust), .data$p.adjust <= 0.05,
    is.finite(.data$NES), .data$NES != 0
  ) |>
  dplyr::mutate(
    direction = ifelse(.data$NES > 0, "positive", "negative")
  ) |>
  dplyr::group_by(
    .data$dataset, .data$phenotype_contrast, .data$spatial_unit,
    .data$theme_id, .data$direction
  ) |>
  dplyr::arrange(
    .data$p.adjust, .data$pvalue, dplyr::desc(abs(.data$NES)), .data$GO_ID,
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 1L) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    gsea_evidence_id = paste(
      .data$dataset, .data$phenotype_contrast, .data$spatial_unit,
      .data$theme_id, .data$direction, "qc_context", sep = "|"
    ),
    dataset, phenotype_contrast, contrast = .data$phenotype_contrast,
    comparison_scope = "local_spatial",
    gsea_spatial_unit = .data$spatial_unit,
    biological_program = .data$theme_id,
    theme_id, manuscript_theme, theme_role, theme_claim_eligible,
    registry_version,
    GO_ID, GO_description, NES,
    raw_p = .data$pvalue, GSEA_FDR = .data$p.adjust,
    leading_edge_accessions = .data$core_enrichment,
    direction
  )
unmatched <- strict_gsea_evidence |>
  dplyr::anti_join(
    program_matches |>
      dplyr::distinct(.data$dataset, .data$biological_program),
    by = c("dataset", "biological_program")
  ) |>
  dplyr::bind_rows(qc_context_evidence)
gww_assert_unique(unmatched, "gsea_evidence_id", "Unmatched strict GSEA evidence")

universe_by_dataset <- lapply(datasets, function(dataset) {
  readr::read_csv(
    dataset_inputs[[dataset]]$feature_universe,
    col_select = tidyselect::any_of(c(
      "ProteinGroupID", "included_in_wgcna", "ModuleID", "member_accessions",
      "RepresentativeUniProt", "MemberUniProts", "representative_accession"
    )),
    show_col_types = FALSE, progress = FALSE, guess_max = Inf
  )
})
names(universe_by_dataset) <- datasets
universe_bundles <- lapply(datasets, function(dataset) {
  gww_build_universe_bundle(universe_by_dataset[[dataset]], dataset)
})
names(universe_bundles) <- datasets
unmatched <- gww_map_leading_edges(
  unmatched,
  dplyr::bind_rows(lapply(universe_bundles, `[[`, "token_map"))
)

module_annotation <- dplyr::bind_rows(lapply(datasets, function(dataset) {
  readr::read_csv(
    dataset_inputs[[dataset]]$module_annotation,
    col_select = tidyselect::any_of(c(
      "dataset", "ModuleID", "cleaned_biological_label",
      "module_program_primary", "module_program_secondary", "label_confidence",
      "cleaned_biological_label_confidence", "evidence_BP", "evidence_CC",
      "evidence_MF", "SupermoduleID", "top_hub_proteins"
    )),
    show_col_types = FALSE, progress = FALSE, guess_max = Inf
  )
})) |>
  dplyr::mutate(
    candidate_annotation_token = dplyr::coalesce(
      as.character(.data$module_program_primary),
      as.character(.data$module_program_secondary)
    ),
    candidate_label_confidence = dplyr::coalesce(
      as.character(.data$cleaned_biological_label_confidence),
      as.character(.data$label_confidence)
    )
  )
gww_assert_unique(
  module_annotation, c("dataset", "ModuleID"), "Module biological annotation"
)

current_labels <- dplyr::bind_rows(lapply(datasets, function(dataset) {
  x <- readr::read_csv(
    dataset_inputs[[dataset]]$final_labels,
    col_select = tidyselect::any_of(c(
      "dataset", "level", "entity_id", "canonical_biological_label",
      "final_plot_label", "best_data_driven_label", "label_confidence",
      "biological_label_confidence", "label_basis", "label_source"
    )),
    show_col_types = FALSE, progress = FALSE, guess_max = Inf
  )
  for (column in c(
    "canonical_biological_label", "final_plot_label",
    "best_data_driven_label", "label_confidence",
    "biological_label_confidence", "label_basis", "label_source"
  )) {
    if (!column %in% names(x)) x[[column]] <- NA_character_
  }
  x |>
    dplyr::filter(.data$level == "module") |>
    dplyr::transmute(
      dataset,
      ModuleID = as.character(.data$entity_id),
      current_reviewed_WGCNA_label = dplyr::coalesce(
        as.character(.data$canonical_biological_label),
        as.character(.data$final_plot_label),
        as.character(.data$best_data_driven_label)
      ),
      current_reviewed_label_confidence = dplyr::coalesce(
        as.character(.data$biological_label_confidence),
        as.character(.data$label_confidence)
      ),
      current_reviewed_label_source = dplyr::coalesce(
        as.character(.data$label_source), as.character(.data$label_basis),
        "WGCNA_final_label_lookup.csv"
      )
    )
}))
gww_assert_unique(
  current_labels, c("dataset", "ModuleID"), "Current final WGCNA labels"
)
module_annotation <- module_annotation |>
  dplyr::left_join(
    current_labels,
    by = c("dataset", "ModuleID"),
    relationship = "one-to-one"
  )

supermodule_annotation <- dplyr::bind_rows(lapply(datasets, function(dataset) {
  readr::read_csv(
    dataset_inputs[[dataset]]$supermodule_annotation,
    col_select = tidyselect::any_of(c(
      "dataset", "SupermoduleID", "cleaned_biological_label",
      "Supermodule_CompositionLabel", "member_modules", "member_ModuleIDs",
      "module_program_primary", "module_program_secondary", "evidence_BP"
    )),
    show_col_types = FALSE, progress = FALSE, guess_max = Inf
  )
})) |>
  dplyr::transmute(
    dataset,
    SupermoduleID = as.character(.data$SupermoduleID),
    candidate_supermodule_label = dplyr::coalesce(
      as.character(.data$cleaned_biological_label),
      as.character(.data$Supermodule_CompositionLabel)
    ),
    candidate_supermodule_annotation = paste(
      dplyr::coalesce(as.character(.data$module_program_primary), ""),
      dplyr::coalesce(as.character(.data$module_program_secondary), ""),
      dplyr::coalesce(as.character(.data$evidence_BP), ""),
      sep = " | "
    )
  ) |>
  dplyr::distinct(.data$dataset, .data$SupermoduleID, .keep_all = TRUE)

module_go <- dplyr::bind_rows(lapply(datasets, function(dataset) {
  readr::read_csv(
    dataset_inputs[[dataset]]$module_go,
    col_select = tidyselect::any_of(c(
      "ModuleID", "ModuleProteinSetType", "Ontology", "ID", "Description",
      "pvalue", "p.adjust", "geneID", "Count"
    )),
    show_col_types = FALSE, progress = FALSE, guess_max = Inf
  ) |>
    dplyr::mutate(dataset = dataset)
})) |>
  dplyr::filter(
    .data$ModuleProteinSetType == "all", .data$Ontology == "BP"
  )

mapping_candidate_rows <- lapply(seq_len(nrow(unmatched)), function(i) {
  row <- unmatched[i, , drop = FALSE]
  dataset <- row$dataset[[1]]
  bundle <- universe_bundles[[dataset]]
  leading <- gww_split_proteins(row$leading_edge_proteins[[1]])
  candidates <- module_annotation |>
    dplyr::filter(.data$dataset == .env$dataset) |>
    dplyr::left_join(
      module_go |>
        dplyr::filter(
          .data$dataset == .env$dataset,
          .data$ID == row$GO_ID[[1]]
        ) |>
        dplyr::transmute(
          dataset, ModuleID,
          exact_GO_ID_match = TRUE,
          module_GO_description = as.character(.data$Description),
          module_GO_raw_p = as.numeric(.data$pvalue),
          module_GO_FDR = as.numeric(.data$p.adjust),
          module_GO_geneID = as.character(.data$geneID),
          module_GO_count = as.integer(.data$Count)
        ),
      by = c("dataset", "ModuleID"),
      relationship = "one-to-one"
    )
  composition <- lapply(seq_len(nrow(candidates)), function(j) {
    module_id <- candidates$ModuleID[[j]]
    module_set <- bundle$module_sets[[module_id]]
    if (is.null(module_set)) {
      return(data.frame(
        dataset = dataset, ModuleID = module_id,
        candidate_overlap_n = 0L,
        candidate_overlap_proteins = "",
        candidate_overlap_Jaccard = NA_real_, stringsAsFactors = FALSE
      ))
    }
    overlap <- intersect(module_set, leading)
    data.frame(
      dataset = dataset, ModuleID = module_id,
      candidate_overlap_n = length(overlap),
      candidate_overlap_proteins = paste(sort(overlap), collapse = ";"),
      candidate_overlap_Jaccard = if (length(union(module_set, leading)))
        length(overlap) / length(union(module_set, leading)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::distinct(.data$dataset, .data$ModuleID, .keep_all = TRUE)
  candidates <- candidates |>
    dplyr::left_join(
      composition, by = c("dataset", "ModuleID"),
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      supermodule_annotation,
      by = c("dataset", "SupermoduleID"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      exact_GO_ID_match = dplyr::coalesce(.data$exact_GO_ID_match, FALSE),
      exact_GO_FDR_supported =
        .data$exact_GO_ID_match & .data$module_GO_FDR <= 0.05,
      candidate_evidence_rank = dplyr::case_when(
        .data$exact_GO_FDR_supported & .data$candidate_overlap_n >= 3L &
          .data$candidate_overlap_Jaccard >= 0.01 &
          !is.na(.data$candidate_annotation_token) ~ 1L,
        .data$exact_GO_FDR_supported ~ 2L,
        .data$exact_GO_ID_match ~ 3L,
        .data$candidate_overlap_n >= 5L &
          .data$candidate_overlap_Jaccard >= 0.02 ~ 4L,
        TRUE ~ 5L
      )
    ) |>
    dplyr::arrange(
      .data$candidate_evidence_rank, .data$module_GO_FDR,
      dplyr::desc(.data$candidate_overlap_Jaccard), .data$ModuleID
    )
  best <- candidates[1L, , drop = FALSE]
  confidence <- dplyr::case_when(
    best$candidate_evidence_rank == 1L ~ "high",
    best$candidate_evidence_rank %in% c(2L, 3L, 4L) ~ "moderate",
    TRUE ~ "none"
  )
  action <- dplyr::case_when(
    row$theme_role[[1]] == "qc_review" ~
      "keep_unmatched_for_biological_claims",
    confidence == "high" ~ "add_mapping",
    confidence == "moderate" ~ "review_manually",
    TRUE ~ "keep_unmatched"
  )
  data.frame(
    dataset = dataset,
    contrast = row$contrast[[1]],
    comparison_scope = row$comparison_scope[[1]],
    spatial_unit = row$gsea_spatial_unit[[1]],
    biological_program = row$biological_program[[1]],
    theme_id = row$theme_id[[1]],
    manuscript_theme = row$manuscript_theme[[1]],
    theme_role = row$theme_role[[1]],
    GO_ID = row$GO_ID[[1]],
    GO_term = row$GO_description[[1]],
    NES = row$NES[[1]],
    GSEA_FDR = row$GSEA_FDR[[1]],
    leading_edge_proteins = row$leading_edge_proteins[[1]],
    reason_no_mapping_available = paste(
      "No approved checked-in mapping exists for this dataset and ontology-aware",
      "theme. QC-review themes are retained as context and are ineligible for",
      "primary biological concordance claims."
    ),
    candidate_WGCNA_entity = if (confidence == "none") NA_character_ else
      best$ModuleID[[1]],
    candidate_WGCNA_label = if (confidence == "none") NA_character_ else
      dplyr::coalesce(
        best$current_reviewed_WGCNA_label[[1]],
        best$cleaned_biological_label[[1]]
      ),
    candidate_reviewed_label_confidence = if (confidence == "none")
      NA_character_ else best$current_reviewed_label_confidence[[1]],
    candidate_reviewed_label_source = if (confidence == "none")
      NA_character_ else best$current_reviewed_label_source[[1]],
    candidate_annotation_token = if (confidence == "none") NA_character_ else
      best$candidate_annotation_token[[1]],
    candidate_parent_supermodule = if (confidence == "none") NA_character_ else
      best$SupermoduleID[[1]],
    candidate_supermodule_label = if (confidence == "none") NA_character_ else
      best$candidate_supermodule_label[[1]],
    exact_GO_ID_match = best$exact_GO_ID_match[[1]],
    module_GO_FDR = best$module_GO_FDR[[1]],
    candidate_overlap_n = best$candidate_overlap_n[[1]],
    candidate_overlap_proteins = best$candidate_overlap_proteins[[1]],
    candidate_overlap_Jaccard = best$candidate_overlap_Jaccard[[1]],
    evidence_for_mapping = if (confidence == "none")
      "No exact module GO-ID support or prespecified substantial leading-edge composition overlap."
    else paste0(
      "exact_GO_ID_match=", best$exact_GO_ID_match[[1]],
      "; module_GO_FDR=", signif(best$module_GO_FDR[[1]], 4),
      "; leading_edge_overlap_n=", best$candidate_overlap_n[[1]],
      "; Jaccard=", signif(best$candidate_overlap_Jaccard[[1]], 4),
      "; parent_supermodule=", dplyr::coalesce(
        best$candidate_supermodule_label[[1]], "none"
      )
    ),
    annotation_source = paste(
      "WGCNA_module_GO_enrichment_long.csv exact GO ID;",
      "WGCNA_module_biological_annotation.csv reviewed label/token;",
      "WGCNA_supermodule_biological_annotation.csv parent context;",
      "canonical WGCNA feature-universe composition"
    ),
    mapping_confidence = confidence,
    recommended_action = action,
    candidate_rule = paste(
      "No fuzzy label matching. high requires exact module GO FDR<=0.05,",
      ">=3 leading-edge proteins, Jaccard>=0.01 and an existing annotation",
      "token; moderate requires exact GO support or >=5 proteins with",
      "Jaccard>=0.02; otherwise keep unmatched. Manual review is required."
    ),
    gsea_evidence_id = row$gsea_evidence_id[[1]],
    stringsAsFactors = FALSE
  )
})
unmatched_audit <- dplyr::bind_rows(mapping_candidate_rows) |>
  dplyr::arrange(
    .data$dataset, .data$contrast, .data$comparison_scope,
    .data$spatial_unit, .data$biological_program, .data$GO_ID
  )
gww_assert_unique(
  unmatched_audit, "gsea_evidence_id", "Unmatched GSEA mapping audit"
)
candidate_keys <- unmatched_audit |>
  dplyr::filter(!is.na(.data$candidate_WGCNA_entity)) |>
  dplyr::distinct(
    .data$dataset, ModuleID = .data$candidate_WGCNA_entity
  )
invalid_candidate_keys <- candidate_keys |>
  dplyr::anti_join(
    module_annotation |>
      dplyr::distinct(.data$dataset, .data$ModuleID),
    by = c("dataset", "ModuleID")
  )
if (nrow(invalid_candidate_keys)) {
  stop(
    "Unmatched-theme audit proposed a WGCNA module outside the same dataset: ",
    paste(
      paste(invalid_candidate_keys$dataset, invalid_candidate_keys$ModuleID,
            sep = "|"),
      collapse = "; "
    ), ".",
    call. = FALSE
  )
}

# Preserve the current strict result at its generated grain and expose each gate.
gate_audit <- gwwd_strict_gate_audit(strict_long, strict_patterns)

# Second layer: complete effects, explicitly descriptive and independent of FDR
# for direction assignment.
local_gsea_all <- gwwd_local_program_effects(gsea_terms)
recurrent_cross_spatial_gsea_all <-
  gwwd_recurrent_cross_spatial_program_effects(local_gsea_all)
global_wgcna_program <- gwwd_wgcna_program_effects(
  handoff, program_matches, scope = "global"
)
local_wgcna_program <- gwwd_wgcna_program_effects(
  handoff, program_matches, scope = "local"
)
recurrent_cross_spatial_patterns <- gwwd_directional_patterns(
  recurrent_cross_spatial_gsea_all |>
    dplyr::semi_join(
      program_matches |>
        dplyr::distinct(.data$dataset, .data$biological_program),
      by = c("dataset", "biological_program")
    ),
  global_wgcna_program,
  scope = "global"
) |>
  dplyr::mutate(
    pattern_scope = "recurrent_cross_spatial_with_global_wgcna"
  ) |>
  dplyr::arrange(.data$dataset, .data$biological_program)
local_patterns <- gwwd_directional_patterns(
  local_gsea_all |>
    dplyr::semi_join(
      program_matches |>
        dplyr::distinct(.data$dataset, .data$biological_program),
      by = c("dataset", "biological_program")
    ),
  local_wgcna_program,
  scope = "local"
) |>
  dplyr::mutate(pattern_scope = "same_spatial_unit_local") |>
  dplyr::arrange(
    .data$dataset, .data$gsea_spatial_unit, .data$biological_program
  )

# Neuropil RES > CON > SUS audit. The available Stage 05 sensitivity table
# covers SUS-RES only; no full three-contrast LOAO claim is made.
neuropil_global <- handoff |>
  dplyr::filter(
    .data$dataset == "neuron_neuropil",
    .data$entity_level == "module",
    .data$effect_scope == "spatial_adjusted_global",
    .data$contrast %in% unname(gww_formal_contrast_map())
  ) |>
  dplyr::transmute(
    dataset = as.character(.data$dataset),
    ModuleID = .data$entity_id,
    WGCNA_label = .data$display_label,
    contrast_suffix = gwwd_contrast_suffix(.data$contrast),
    estimate = .data$estimate,
    SE = .data$SE,
    CI_low = .data$CI_low,
    CI_high = .data$CI_high,
    raw_p = .data$p_value,
    tier_FDR = .data$tier_specific_fdr,
    model_stability = .data$model_stability_status
  ) |>
  tidyr::pivot_wider(
    names_from = "contrast_suffix",
    values_from = dplyr::all_of(c(
      "estimate", "SE", "CI_low", "CI_high", "raw_p", "tier_FDR",
      "model_stability"
    )),
    names_glue = "{.value}__{contrast_suffix}"
  ) |>
  dplyr::mutate(
    WGCNA_RES_gt_CON_gt_SUS_point_pattern =
      .data$estimate__RES_CON > 0 & .data$estimate__SUS_CON < 0,
    direct_SUS_RES_point_direction_compatible = .data$estimate__SUS_RES < 0
  )

neuropil_sensitivity <- readr::read_csv(
  dataset_inputs$neuron_neuropil$group_sensitivity,
  show_col_types = FALSE, progress = FALSE, guess_max = Inf
) |>
  dplyr::filter(
    .data$level == "module", .data$contrast == "SUS - RES"
  ) |>
  dplyr::mutate(dataset = "neuron_neuropil") |>
  dplyr::group_by(.data$dataset, .data$endpoint_id) |>
  dplyr::summarise(
    n_available_SUS_RES_sensitivity_models = dplyr::n(),
    sensitivity_tiers = paste(sort(.data$analysis_tier), collapse = ";"),
    n_negative_SUS_RES_sensitivity_estimates = sum(.data$estimate < 0),
    available_SUS_RES_sensitivity_sign_stable = all(.data$estimate < 0),
    all_available_sensitivity_models_stable = all(
      .data$model_valid_for_inference %in% TRUE &
        .data$primary_model_stable %in% TRUE
    ),
    sensitivity_estimate_min = min(.data$estimate),
    sensitivity_estimate_max = max(.data$estimate),
    .groups = "drop"
  ) |>
  dplyr::rename(ModuleID = "endpoint_id")

neuropil_programs <- program_matches |>
  dplyr::filter(
    .data$dataset == "neuron_neuropil",
    .data$entity_level == "module"
  ) |>
  dplyr::transmute(
    dataset = as.character(.data$dataset),
    ModuleID = .data$entity_id,
    biological_program,
    program_mapping_token = .data$annotation_token,
    program_mapping_field = .data$annotation_field
  )
all_neuropil_modules <- neuropil_global |>
  dplyr::select("dataset", "ModuleID") |>
  dplyr::left_join(
    neuropil_programs, by = c("dataset", "ModuleID"),
    relationship = "one-to-many"
  )

neuropil_overlap <- strict_overlap |>
  dplyr::filter(.data$dataset == "neuron_neuropil") |>
  dplyr::group_by(
    .data$dataset, .data$entity_id, .data$biological_program
  ) |>
  dplyr::summarise(
    best_program_specific_overlap_FDR = min(.data$overlap_FDR),
    best_program_specific_overlap_Jaccard =
      .data$Jaccard[[which.min(.data$overlap_FDR)]],
    any_program_specific_overlap_FDR_supported =
      any(.data$overlap_FDR <= 0.05),
    overlap_scopes = paste(sort(unique(.data$comparison_scope)), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::rename(ModuleID = "entity_id")

neuropil_recurrent_cross_spatial_gsea <- recurrent_cross_spatial_patterns |>
  dplyr::filter(.data$dataset == "neuron_neuropil") |>
  dplyr::select(dplyr::all_of(c(
    "dataset", "biological_program", "gsea_direction_sign__RES_CON",
    "gsea_direction_sign__SUS_CON", "gsea_direction_sign__SUS_RES",
    "gsea_program_NES__RES_CON", "gsea_program_NES__SUS_CON",
    "gsea_program_NES__SUS_RES", "gsea_min_FDR__RES_CON",
    "gsea_min_FDR__SUS_CON", "gsea_min_FDR__SUS_RES",
    "recurrent_cross_spatial_direction_available__RES_CON",
    "recurrent_cross_spatial_direction_available__SUS_CON",
    "recurrent_cross_spatial_direction_available__SUS_RES", "directional_pattern"
  )))
neuropil_local_support <- local_patterns |>
  dplyr::filter(.data$dataset == "neuron_neuropil") |>
  dplyr::group_by(.data$dataset, .data$biological_program) |>
  dplyr::summarise(
    n_local_RES_gt_CON_gt_SUS_units = sum(
      .data$directional_pattern == "RES_gt_CON_gt_SUS"
    ),
    local_RES_gt_CON_gt_SUS_units = paste(
      sort(.data$gsea_spatial_unit[
        .data$directional_pattern == "RES_gt_CON_gt_SUS"
      ]), collapse = ";"
    ),
    .groups = "drop"
  )

neuropil_audit <- all_neuropil_modules |>
  dplyr::left_join(
    neuropil_global, by = c("dataset", "ModuleID"),
    relationship = "many-to-one"
  ) |>
  dplyr::left_join(
    neuropil_sensitivity, by = c("dataset", "ModuleID"),
    relationship = "many-to-one"
  ) |>
  dplyr::left_join(
    neuropil_recurrent_cross_spatial_gsea,
    by = c("dataset", "biological_program"),
    relationship = "many-to-one"
  ) |>
  dplyr::left_join(
    neuropil_local_support, by = c("dataset", "biological_program"),
    relationship = "many-to-one"
  ) |>
  dplyr::left_join(
    neuropil_overlap,
    by = c("dataset", "ModuleID", "biological_program"),
    relationship = "many-to-one"
  ) |>
  dplyr::mutate(
    full_three_contrast_LOAO_status =
      "not_available; only four existing SUS-RES sensitivity models assessed",
    recurrent_cross_spatial_GSEA_contradictory = dplyr::coalesce(
      (!is.na(.data$gsea_direction_sign__RES_CON) &
         .data$gsea_direction_sign__RES_CON != sign(.data$estimate__RES_CON)) |
        (!is.na(.data$gsea_direction_sign__SUS_CON) &
           .data$gsea_direction_sign__SUS_CON != sign(.data$estimate__SUS_CON)) |
        (!is.na(.data$gsea_direction_sign__SUS_RES) &
           .data$gsea_direction_sign__SUS_RES != sign(.data$estimate__SUS_RES)),
      FALSE
    ),
    all_three_recurrent_cross_spatial_GSEA_compatible = dplyr::coalesce(
      .data$recurrent_cross_spatial_direction_available__RES_CON %in% TRUE &
        .data$recurrent_cross_spatial_direction_available__SUS_CON %in% TRUE &
        .data$recurrent_cross_spatial_direction_available__SUS_RES %in% TRUE &
        .data$gsea_direction_sign__RES_CON == sign(.data$estimate__RES_CON) &
        .data$gsea_direction_sign__SUS_CON == sign(.data$estimate__SUS_CON) &
        .data$gsea_direction_sign__SUS_RES == sign(.data$estimate__SUS_RES),
      FALSE
    ),
    GSEA_support_scope = dplyr::case_when(
      is.na(.data$biological_program) ~ "unsupported_by_GSEA",
      .data$recurrent_cross_spatial_GSEA_contradictory ~ "contradictory",
      .data$all_three_recurrent_cross_spatial_GSEA_compatible ~
        "recurrent_cross_spatial",
      .data$n_local_RES_gt_CON_gt_SUS_units > 0 ~ "spatially_localized_only",
      TRUE ~ "unsupported_by_GSEA"
    ),
    interpretation_boundary = paste(
      "Point-estimate pattern only. RES-CON/SUS-CON significance asymmetry is",
      "not used. Sensitivity sign stability is partial because available",
      "sensitivity models cover SUS-RES only."
    )
  ) |>
  dplyr::arrange(.data$ModuleID, .data$biological_program)
gww_assert_unique(
  neuropil_audit,
  c("dataset", "ModuleID", "biological_program"),
  "Neuropil RES-CON-SUS pattern audit"
)

# Summary records are deliberately long and machine-auditable.
summary_row <- function(section, metric, value, detail = NA_character_) {
  data.frame(
    section = section, metric = metric, value = as.character(value),
    detail = detail, stringsAsFactors = FALSE
  )
}
strict_counts <- strict_long |>
  dplyr::count(.data$concordance_class)
imprecise_counts <- imprecise |>
  dplyr::count(.data$diagnostic_subgroup) |>
  tidyr::complete(
    diagnostic_subgroup = c(
      "stronger_imprecise_support", "directional_only_imprecise",
      "broad_CI_or_weak_precision", "other"
    ),
    fill = list(n = 0L)
  )
imprecise_by_dataset <- imprecise |>
  dplyr::count(.data$dataset, .data$diagnostic_subgroup)
imprecise_by_contrast <- imprecise |>
  dplyr::count(.data$contrast, .data$diagnostic_subgroup)
imprecise_by_program <- imprecise |>
  dplyr::count(.data$biological_program, .data$diagnostic_subgroup)
imprecise_by_scope <- imprecise |>
  dplyr::count(.data$local_global_scope, .data$diagnostic_subgroup)
fdr_bins <- power |>
  dplyr::count(.data$WGCNA_FDR_bin)
recurrent_cross_spatial_pattern_counts <- recurrent_cross_spatial_patterns |>
  dplyr::count(.data$directional_pattern)
local_pattern_counts <- local_patterns |>
  dplyr::count(.data$directional_pattern)
diagnostic_summary <- dplyr::bind_rows(
  lapply(seq_len(nrow(strict_counts)), function(i) summary_row(
    "strict_class_count", strict_counts$concordance_class[[i]],
    strict_counts$n[[i]], "Read unchanged from the corrected strict concordance output."
  )),
  lapply(seq_len(nrow(imprecise_counts)), function(i) summary_row(
    "concordant_imprecise_subgroup", imprecise_counts$diagnostic_subgroup[[i]],
    imprecise_counts$n[[i]], "Diagnostic descriptor; not an inferential class."
  )),
  lapply(seq_len(nrow(imprecise_by_dataset)), function(i) summary_row(
    "concordant_imprecise_by_dataset",
    paste(
      imprecise_by_dataset$dataset[[i]],
      imprecise_by_dataset$diagnostic_subgroup[[i]], sep = "|"
    ), imprecise_by_dataset$n[[i]]
  )),
  lapply(seq_len(nrow(imprecise_by_contrast)), function(i) summary_row(
    "concordant_imprecise_by_contrast",
    paste(
      imprecise_by_contrast$contrast[[i]],
      imprecise_by_contrast$diagnostic_subgroup[[i]], sep = "|"
    ), imprecise_by_contrast$n[[i]]
  )),
  lapply(seq_len(nrow(imprecise_by_program)), function(i) summary_row(
    "concordant_imprecise_by_program",
    paste(
      imprecise_by_program$biological_program[[i]],
      imprecise_by_program$diagnostic_subgroup[[i]], sep = "|"
    ), imprecise_by_program$n[[i]]
  )),
  lapply(seq_len(nrow(imprecise_by_scope)), function(i) summary_row(
    "concordant_imprecise_by_scope",
    paste(
      imprecise_by_scope$local_global_scope[[i]],
      imprecise_by_scope$diagnostic_subgroup[[i]], sep = "|"
    ), imprecise_by_scope$n[[i]]
  )),
  lapply(seq_len(nrow(fdr_bins)), function(i) summary_row(
    "WGCNA_FDR_distribution", fdr_bins$WGCNA_FDR_bin[[i]],
    fdr_bins$n[[i]], "Unchanged tier-specific WGCNA FDR bins."
  )),
  list(
    summary_row(
      "WGCNA_power", "minimum_tier_specific_FDR",
      min(power$wgcna_tier_specific_fdr, na.rm = TRUE)
    ),
    summary_row("WGCNA_power", "unique_matched_endpoints", nrow(power)),
    summary_row(
      "WGCNA_power", "raw_p_below_0_05",
      sum(power$raw_p_below_0_05, na.rm = TRUE)
    ),
    summary_row(
      "WGCNA_power", "FDR_at_or_below_0_05",
      sum(power$FDR_at_or_below_0_05, na.rm = TRUE)
    ),
    summary_row(
      "WGCNA_power", "FDR_at_or_below_0_10_descriptive",
      sum(power$FDR_at_or_below_0_10_descriptive, na.rm = TRUE),
      "Descriptive only; no threshold change."
    ),
    summary_row(
      "WGCNA_power", "CI_excludes_zero",
      sum(power$CI_excludes_zero, na.rm = TRUE)
    ),
    summary_row(
      "structural_overlap", "matched_endpoints_with_supported_overlap",
      sum(power$any_overlap_FDR_supported, na.rm = TRUE),
      "Structural convergence only; not animal-level confirmation."
    ),
    summary_row(
      "structural_overlap", "matched_comparisons_with_supported_overlap",
      sum(strict_long$overlap_FDR <= 0.05, na.rm = TRUE),
      "Comparison-level count; structural convergence only."
    ),
    summary_row(
      "concordant_imprecise", "overlap_FDR_at_or_below_0_05",
      sum(imprecise$overlap_FDR <= 0.05, na.rm = TRUE)
    ),
    summary_row("mapping", "unmatched_GSEA_evidence_rows", nrow(unmatched_audit)),
    summary_row(
      "mapping", "candidate_add_mapping",
      sum(unmatched_audit$recommended_action == "add_mapping")
    ),
    summary_row(
      "mapping", "candidate_review_manually",
      sum(unmatched_audit$recommended_action == "review_manually")
    ),
    summary_row(
      "mapping", "candidate_keep_unmatched",
      sum(unmatched_audit$recommended_action == "keep_unmatched")
    ),
    summary_row(
      "mapping", "qc_keep_unmatched_for_biological_claims",
      sum(
        unmatched_audit$recommended_action ==
          "keep_unmatched_for_biological_claims"
      )
    ),
    summary_row(
      "strict_adaptive", "significance_gated", TRUE,
      "Strict GSEA evidence is filtered at canonical FDR<0.05 before adaptive classification."
    ),
    summary_row("strict_adaptive", "unresolved_rows", nrow(gate_audit))
  ),
  lapply(seq_len(nrow(recurrent_cross_spatial_pattern_counts)), function(i) summary_row(
    "recurrent_cross_spatial_descriptive_pattern",
    recurrent_cross_spatial_pattern_counts$directional_pattern[[i]],
    recurrent_cross_spatial_pattern_counts$n[[i]],
    "DESCRIPTIVE / HYPOTHESIS-GENERATING; NOT AN INFERENTIAL CLAIM."
  )),
  lapply(seq_len(nrow(local_pattern_counts)), function(i) summary_row(
    "local_descriptive_pattern", local_pattern_counts$directional_pattern[[i]],
    local_pattern_counts$n[[i]],
    "Same-spatial-unit only; descriptive, not inferential."
  )),
  list(
    summary_row(
      "neuropil_pattern", "modules_RES_gt_CON_gt_SUS_point_estimates",
      dplyr::n_distinct(neuropil_audit$ModuleID[
        neuropil_audit$WGCNA_RES_gt_CON_gt_SUS_point_pattern %in% TRUE
      ])
    ),
    summary_row(
      "neuropil_pattern", "modules_sign_stable_available_SUS_RES_sensitivity",
      dplyr::n_distinct(neuropil_audit$ModuleID[
        neuropil_audit$WGCNA_RES_gt_CON_gt_SUS_point_pattern %in% TRUE &
          neuropil_audit$available_SUS_RES_sensitivity_sign_stable %in% TRUE
      ]),
      "Partial sensitivity only; full three-contrast LOAO is unavailable."
    )
  )
)

# Assertions: source classifications, exact expected grains, and no row
# multiplication in the principal diagnostics.
if (nrow(imprecise) != sum(
  strict_long$concordance_class == "concordant_imprecise"
)) {
  stop("Concordant-imprecise diagnostic lost or multiplied rows.", call. = FALSE)
}
if (nrow(gate_audit) != nrow(strict_patterns)) {
  stop("Strict adaptive gate audit does not preserve its input row grain.",
       call. = FALSE)
}
if (nrow(unmatched_audit) != nrow(unmatched)) {
  stop("Unmatched mapping audit does not preserve one row per exclusion.",
       call. = FALSE)
}
if (!identical(
  sort(unique(imprecise$official_concordance_class)),
  "concordant_imprecise"
)) {
  stop("An official concordance class was altered in diagnostics.", call. = FALSE)
}

hash_after <- source_hashes(required_inputs)
hash_audit <- hash_before |>
  dplyr::rename(source_md5_before = "source_md5") |>
  dplyr::left_join(
    hash_after |>
      dplyr::rename(source_md5_after = "source_md5"),
    by = c("input_name", "source_file"),
    relationship = "one-to-one"
  ) |>
  dplyr::mutate(unchanged = .data$source_md5_before == .data$source_md5_after)
if (any(!hash_audit$unchanged)) {
  stop("A protected input changed during diagnostic generation.", call. = FALSE)
}

tables <- list(
  concordant_imprecise_diagnostic.csv = imprecise,
  concordance_wgcna_power_diagnostic.csv = power,
  unmatched_gsea_program_mapping_audit.csv = unmatched_audit,
  unmatched_gsea_theme_mapping_audit.csv = unmatched_audit,
  adaptive_resilience_gate_audit.csv = gate_audit,
  adaptive_resilience_directional_patterns.csv = recurrent_cross_spatial_patterns,
  adaptive_resilience_local_spatial_patterns.csv = local_patterns,
  neuropil_RES_CON_SUS_pattern_audit.csv = neuropil_audit,
  diagnostic_summary.csv = diagnostic_summary
)
for (filename in names(tables)) {
  invisible(write_integration_table(tables[[filename]], paths, filename))
}
write_csv_safe(hash_audit, file.path(paths$reports, "protected_source_hash_audit.csv"))

imprecise_subgroups <- imprecise |>
  dplyr::count(.data$diagnostic_subgroup)
recurrent_cross_spatial_counts <- recurrent_cross_spatial_patterns |>
  dplyr::count(.data$directional_pattern)
local_counts <- local_patterns |>
  dplyr::count(.data$directional_pattern)
readme_lines <- c(
  "# GSEA/GO <-> WGCNA concordance diagnostics",
  "",
  paste0("Generated: ", format(Sys.time(), tz = "Europe/Berlin", usetz = TRUE)),
  "",
  "## Boundary",
  "",
  "The strict concordance classes and adaptive-resilience result are read unchanged from the corrected ontology-aware concordance output. This diagnostic does not refit enrichment or WGCNA, alter p-values/FDRs, relax a threshold, change a claim gate, combine p-values, or create a composite concordance statistic.",
  "",
  "GSEA and WGCNA use the same proteomics data and are not independent replications. Leading-edge/module overlap is structural convergence, not animal-level group-effect significance. WGCNA GO annotation supplies biological identity only.",
  "",
  "## Why strict FDR-supported concordance is zero",
  "",
  paste0(
    "Among ", nrow(power), " unique biologically matched Stage 07 endpoints, the minimum unchanged tier-specific WGCNA FDR is ",
    signif(min(power$wgcna_tier_specific_fdr, na.rm = TRUE), 4), ". Raw p<0.05 occurs for ",
    sum(power$raw_p_below_0_05, na.rm = TRUE), " endpoints; FDR<=0.05 occurs for ",
    sum(power$FDR_at_or_below_0_05, na.rm = TRUE), "; FDR<=0.10 occurs for ",
    sum(power$FDR_at_or_below_0_10_descriptive, na.rm = TRUE), " (descriptive only)."
  ),
  "",
  "The zero is therefore expected from the current WGCNA multiplicity/power rather than absence of structural biological matching; supported program-specific overlaps remain visible but do not change WGCNA inference.",
  "",
  "## concordant_imprecise diagnostic descriptors",
  "",
  paste0("- ", imprecise_subgroups$diagnostic_subgroup, ": ", imprecise_subgroups$n),
  "",
  "These are descriptive flags only. Extremely broad CI means the top decile of CI width within the same dataset/entity-level/tier/effect-scope distribution. Supermodules are not Fisher-tested as module unions, so missing overlap is retained explicitly.",
  "",
  "## Strict versus descriptive adaptive patterns",
  "",
  "The strict classifier is significance-gated because its input builder requires canonical GSEA FDR<0.05 before recurrent-cross-spatial evidence is constructed. That rule remains unchanged.",
  "",
  "The second layer is labelled DESCRIPTIVE / HYPOTHESIS-GENERATING; NOT AN INFERENTIAL CLAIM CLASSIFIER. It uses the median NES across a fixed set of GO IDs present in all three contrasts within each dataset/spatial-unit/program. The same GO set is used across contrasts. A single anchor GO term is selected once across contrasts for auditability, but it does not determine direction. Recurrent-cross-spatial direction uses at least two spatial units and a non-tied majority sign; it is not a formal across-all-spatial-units GSEA inferential test. FDR is reported but is not a direction gate.",
  "",
  "Recurrent-cross-spatial descriptive counts paired with global WGCNA:",
  paste0("- ", recurrent_cross_spatial_counts$directional_pattern, ": ", recurrent_cross_spatial_counts$n),
  "",
  "Local same-spatial-unit descriptive counts:",
  paste0("- ", local_counts$directional_pattern, ": ", local_counts$n),
  "",
  "RES-supported-shift and SUS-supported-shift labels are deliberately not assigned from significant/non-significant asymmetry or an invented GSEA near-zero threshold. Rows lacking three compatible directions remain insufficient_data; sign conflicts are inconsistent.",
  "",
  "## Neuropil sensitivity limitation",
  "",
  "Current Stage 05 sensitivity outputs provide four SUS-RES models (left-only, right-only, complete bilateral pairs, animal-wide mean). They do not provide full RES-CON/SUS-CON LOAO estimates. The Neuropil audit therefore reports available SUS-RES sign stability and labels the full three-contrast LOAO assessment unavailable.",
  ""
)
dir_create(paths$reports)
writeLines(readme_lines, file.path(paths$reports, "README.md"), useBytes = TRUE)

write_integration_manifest(
  paths,
  required_inputs,
  list(
    tables = paths$tables,
    source_data = paths$source_data,
    report = file.path(paths$reports, "README.md"),
    protected_hash_audit = file.path(
      paths$reports, "protected_source_hash_audit.csv"
    )
  ),
  list(
    analysis_role = "additive_diagnostic_and_descriptive_only",
    strict_classifications_modified = FALSE,
    descriptive_GSEA_method =
      "median NES over contrast-independent common GO term sets",
    p_values_combined = FALSE,
    thresholds_relaxed = FALSE,
    manuscript_figure_generated = FALSE
  ),
  notes = paste(
    "Additive diagnostic review only. Strict concordance/adaptive classes",
    "remain unchanged; complete-effect directional patterns are explicitly",
    "descriptive and hypothesis-generating."
  )
)

message("GSEA/WGCNA concordance diagnostics complete: ", paths$tables)
