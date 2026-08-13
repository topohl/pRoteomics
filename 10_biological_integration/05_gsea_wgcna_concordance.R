#!/usr/bin/env Rscript
# Script: 10_biological_integration/05_gsea_wgcna_concordance.R
# Stage: additive downstream biological integration
# Consumes: canonical spatial ranked-GSEA terms, Stage 07 inferential handoffs,
#   canonical WGCNA feature universes, and existing robustness/readiness outputs.
# Produces: transparent GSEA-program <-> WGCNA concordance and module-overlap
#   tables. No upstream inference, WGCNA state, p-value, FDR, or claim gate is
#   recomputed or modified.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "wgcna_group_effect_consumer_utils.R"))
source(repo_path("R", "gsea_wgcna_concordance_utils.R"))

SCRIPT_ID <- "10_biological_integration/05_gsea_wgcna_concordance.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
run <- integration_cli(default_dataset = "all", allow_all = TRUE)
datasets <- integration_datasets(run$dataset)
paths <- integration_paths("gsea_wgcna_concordance", "global")

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

gsea_terms_file <- path_results(
  "source_data", "04_differential_expression_enrichment",
  "compareGO_spatial_atlas", "spatial_atlas_enrichment_long.csv"
)
program_mapping_file <- repo_path(
  "config", "gsea_wgcna_program_module_mapping.csv"
)

dataset_inputs <- lapply(datasets, function(dataset) {
  list(
    wgcna_handoff = path_results(
      "tables", "06_modules_WGCNA", "interpretable_summary", dataset,
      "WGCNA_inferential_handoff.csv"
    ),
    wgcna_universe = path_results(
      "tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules",
      "WGCNA_feature_universe.csv"
    ),
    module_preservation = path_results(
      "tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules",
      "WGCNA_module_preservation_summary.csv"
    )
  )
})
names(dataset_inputs) <- datasets

microglia_inputs <- list(
  claim_readiness = path_results(
    "tables", "06_modules_WGCNA", "claim_readiness", "microglia",
    "WGCNA_entity_claim_readiness.csv"
  ),
  module_robustness = path_results(
    "reviewer_audit", "microglia_wgcna_nature_readiness",
    "module_robustness_consensus.csv"
  ),
  block_readiness = path_results(
    "reviewer_audit", "microglia_wgcna_nature_readiness",
    "higher_order_block_readiness_summary.csv"
  ),
  legacy_overlap = path_results(
    "tables", "06_modules_WGCNA", "04_wgcna_de_gsea_overlap", "microglia",
    "WGCNA_vs_DE_GSEA_overlap.csv"
  )
)

required_inputs <- c(
  gsea_terms = gsea_terms_file,
  program_mapping = program_mapping_file,
  unlist(lapply(names(dataset_inputs), function(dataset) {
    c(
      setNames(dataset_inputs[[dataset]]$wgcna_handoff,
               paste0(dataset, "_wgcna_handoff")),
      setNames(dataset_inputs[[dataset]]$wgcna_universe,
               paste0(dataset, "_wgcna_universe"))
    )
  }))
)
optional_inputs <- c(
  unlist(lapply(names(dataset_inputs), function(dataset) {
    setNames(
      dataset_inputs[[dataset]]$module_preservation,
      paste0(dataset, "_module_preservation")
    )
  })),
  setNames(unlist(microglia_inputs), paste0("microglia_", names(microglia_inputs)))
)

output_files <- c(
  "gsea_wgcna_concordance_long.csv",
  "gsea_wgcna_concordance_summary.csv",
  "adaptive_resilience_pattern_summary.csv",
  "program_specific_leading_edge_module_overlap.csv",
  "gsea_wgcna_concordance_audit.csv",
  "program_module_matching_audit.csv"
)

if (run$dry_run) {
  dry_run_inputs(SCRIPT_ID, c(required_inputs, optional_inputs))
  for (filename in output_files) {
    dry_run_line("Output table", file.path(paths$tables, filename))
    dry_run_line("Source-data copy", file.path(paths$source_data, filename))
  }
  dry_run_line(
    "Interpretation README",
    file.path(paths$reports, "README.md")
  )
  quit(save = "no", status = 0L)
}

missing_required <- required_inputs[!file.exists(required_inputs)]
if (length(missing_required)) {
  stop("Required canonical input(s) missing: ",
       paste(names(missing_required), missing_required, sep = "=", collapse = "; "),
       ".", call. = FALSE)
}

source_hashes <- function(paths_named) {
  values <- unname(paths_named)
  hashes <- as.character(tools::md5sum(values))
  data.frame(
    input_name = names(paths_named),
    source_file = normalizePath(values, winslash = "/", mustWork = FALSE),
    source_md5 = hashes,
    stringsAsFactors = FALSE
  )
}

protected_sources <- c(
  required_inputs,
  optional_inputs[file.exists(optional_inputs)]
)
hash_before <- source_hashes(protected_sources)

input_status <- dplyr::bind_rows(lapply(seq_along(c(required_inputs, optional_inputs)), function(i) {
  all_inputs <- c(required_inputs, optional_inputs)
  required <- names(all_inputs)[[i]] %in% names(required_inputs)
  path <- all_inputs[[i]]
  data.frame(
    input_name = names(all_inputs)[[i]],
    source_file = normalizePath(path, winslash = "/", mustWork = FALSE),
    required = required,
    status = if (file.exists(path)) "present" else if (required)
      "missing_required" else "missing_optional",
    stringsAsFactors = FALSE
  )
}))

gsea_terms <- readr::read_csv(
  gsea_terms_file,
  show_col_types = FALSE,
  progress = FALSE,
  guess_max = Inf
)
gww_validate_gsea_terms(gsea_terms)

wgcna_by_dataset <- lapply(datasets, function(dataset) {
  wgcna_inferential_handoff_read(
    dataset_inputs[[dataset]]$wgcna_handoff,
    paste(dataset, "WGCNA_inferential_handoff.csv")
  )
})
names(wgcna_by_dataset) <- datasets
wgcna <- dplyr::bind_rows(wgcna_by_dataset)
wgcna_inferential_handoff_validate(wgcna, "Combined WGCNA inferential handoff")
if (any(!wgcna$independent_hypothesis %in% TRUE) ||
    any(wgcna$claim_entity_role == "compatibility_alias")) {
  stop("Combined WGCNA handoff contains non-independent or alias rows.",
       call. = FALSE)
}

universe_by_dataset <- lapply(datasets, function(dataset) {
  readr::read_csv(
    dataset_inputs[[dataset]]$wgcna_universe,
    col_select = tidyselect::any_of(c(
      "ProteinGroupID", "included_in_wgcna", "ModuleID", "member_accessions",
      "RepresentativeUniProt", "MemberUniProts", "representative_accession"
    )),
    show_col_types = FALSE,
    progress = FALSE,
    guess_max = Inf
  )
})
names(universe_by_dataset) <- datasets
universe_bundles <- lapply(datasets, function(dataset) {
  gww_build_universe_bundle(universe_by_dataset[[dataset]], dataset)
})
names(universe_bundles) <- datasets
token_map <- dplyr::bind_rows(lapply(universe_bundles, `[[`, "token_map"))

registry <- gww_read_program_mapping(program_mapping_file)
program_matches <- gww_build_entity_program_matches(wgcna, registry)

local_gsea <- gww_build_local_gsea_evidence(gsea_terms)
global_gsea <- gww_build_global_gsea_evidence(
  gsea_terms, local_gsea, min_recurrent_units = 2L
)
gsea_evidence <- dplyr::bind_rows(local_gsea, global_gsea)
gww_assert_unique(gsea_evidence, "gsea_evidence_id", "Combined GSEA evidence")
gsea_evidence <- gww_map_leading_edges(gsea_evidence, token_map)

overlap <- gww_program_specific_overlap(
  gsea_evidence, program_matches, universe_bundles
)

for (column in c(
  "cleaned_biological_label", "module_program_primary",
  "module_program_secondary", "annotation_class"
)) {
  if (!column %in% names(wgcna)) wgcna[[column]] <- NA_character_
}

effects <- wgcna |>
  dplyr::filter(.data$contrast %in% unname(gww_formal_contrast_map())) |>
  dplyr::transmute(
    dataset = as.character(.data$dataset),
    wgcna_entity_level = as.character(.data$entity_level),
    wgcna_entity_id = as.character(.data$entity_id),
    wgcna_display_label = as.character(.data$display_label),
    wgcna_canonical_biological_label = dplyr::coalesce(
      as.character(.data$cleaned_biological_label),
      as.character(.data$display_label)
    ),
    module_program_primary = as.character(.data$module_program_primary),
    module_program_secondary = as.character(.data$module_program_secondary),
    contrast = as.character(.data$contrast),
    wgcna_analysis_tier = as.character(.data$analysis_tier),
    wgcna_effect_scope = as.character(.data$effect_scope),
    wgcna_spatial_unit = gww_normalize_spatial_unit(.data$spatial_unit),
    wgcna_estimate = as.numeric(.data$estimate),
    wgcna_SE = as.numeric(.data$SE),
    wgcna_CI_low = as.numeric(.data$CI_low),
    wgcna_CI_high = as.numeric(.data$CI_high),
    wgcna_p_value = as.numeric(.data$p_value),
    wgcna_tier_specific_fdr = as.numeric(.data$tier_specific_fdr),
    wgcna_tier_specific_family_id = as.character(.data$tier_specific_family_id),
    wgcna_tier_specific_family_size = as.integer(.data$tier_specific_family_size),
    wgcna_biological_n = as.numeric(.data$biological_n),
    wgcna_model_valid = as.logical(.data$model_valid),
    wgcna_model_stability_status = as.character(.data$model_stability_status),
    wgcna_support_class = as.character(.data$support_class),
    wgcna_claim_gate = as.character(.data$claim_gate),
    wgcna_independent_hypothesis = as.logical(.data$independent_hypothesis),
    wgcna_claim_entity_role = as.character(.data$claim_entity_role),
    wgcna_result_scope = as.character(.data$result_scope),
    wgcna_test_type = as.character(.data$test_type),
    wgcna_source_artifact = as.character(.data$source_artifact),
    wgcna_source_key = as.character(.data$source_key),
    wgcna_source_key_contract = as.character(.data$source_key_contract),
    wgcna_annotation_class = as.character(.data$annotation_class)
  )
gww_assert_unique(
  effects,
  c(
    "dataset", "wgcna_entity_level", "wgcna_entity_id", "contrast",
    "wgcna_analysis_tier", "wgcna_effect_scope", "wgcna_spatial_unit",
    "wgcna_test_type"
  ),
  "WGCNA named-contrast endpoints"
)

interaction <- wgcna |>
  dplyr::filter(.data$analysis_tier == "secondary_spatial_heterogeneity") |>
  dplyr::transmute(
    dataset = as.character(.data$dataset),
    wgcna_entity_level = as.character(.data$entity_level),
    wgcna_entity_id = as.character(.data$entity_id),
    interaction_p_value = as.numeric(.data$p_value),
    interaction_FDR = as.numeric(.data$tier_specific_fdr),
    interaction_family_id = as.character(.data$tier_specific_family_id),
    interaction_family_size = as.integer(.data$tier_specific_family_size),
    interaction_model_valid = as.logical(.data$model_valid),
    interaction_source_key = as.character(.data$source_key)
  )
gww_assert_unique(
  interaction,
  c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
  "WGCNA interaction endpoints"
)

preservation <- dplyr::bind_rows(lapply(datasets, function(dataset) {
  path <- dataset_inputs[[dataset]]$module_preservation
  if (!file.exists(path)) return(tibble::tibble())
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE) |>
    dplyr::transmute(
      dataset = dataset,
      wgcna_entity_level = "module",
      wgcna_entity_id = as.character(.data$ModuleID),
      preservation_min_Zsummary = as.numeric(.data$preservation_min_Zsummary),
      preservation_status = as.character(.data$preservation_interpretation),
      preservation_note = as.character(.data$preservation_warning)
    )
}))
if (nrow(preservation)) {
  gww_assert_unique(
    preservation,
    c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
    "Module preservation summary"
  )
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) return(tibble::tibble())
  readr::read_csv(
    path, show_col_types = FALSE, progress = FALSE, guess_max = Inf,
    na = c("", "NA")
  )
}

stage13 <- read_optional_csv(microglia_inputs$claim_readiness)
if (nrow(stage13)) {
  stage13 <- stage13 |>
    dplyr::filter(
      .data$level %in% c("module", "supermodule"),
      .data$claim_entity_role != "compatibility_alias"
    ) |>
    dplyr::transmute(
      dataset = "microglia",
      wgcna_entity_level = as.character(.data$level),
      wgcna_entity_id = as.character(.data$entity_id),
      primary_architecture_status = as.character(.data$primary_architecture_status),
      spatial_dependence_class = as.character(.data$spatial_dependence_class),
      animal_stability_status = as.character(.data$animal_stability_status),
      strict_nonspatial_sensitivity_status = as.character(
        .data$strict_nonspatial_sensitivity_status
      ),
      conventional_preservation_status = as.character(
        .data$conventional_preservation_status
      ),
      biological_process_confidence = as.character(.data$biological_process_confidence),
      context_assignment_class = as.character(.data$context_assignment_class),
      neuropil_independence_status = as.character(.data$neuropil_independence_status),
      stage13_group_effect_status = as.character(.data$group_effect_status),
      stage13_allowed_claim_scope = as.character(.data$allowed_claim_scope),
      stage13_readiness_contract_version = as.character(.data$readiness_contract_version)
    )
  gww_assert_unique(
    stage13,
    c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
    "Stage 13 readiness rows"
  )
}

module_robustness <- read_optional_csv(microglia_inputs$module_robustness)
if (nrow(module_robustness)) {
  module_robustness <- module_robustness |>
    dplyr::transmute(
      dataset = "microglia",
      wgcna_entity_level = "module",
      wgcna_entity_id = as.character(.data$ModuleID),
      loao_status = as.character(.data$animal_stability_status),
      loao_min_eigengene_correlation = as.numeric(.data$loao_min_ME_correlation),
      loao_min_hub_rank_spearman = as.numeric(.data$loao_min_hub_rank_spearman),
      loao_min_top25_retained = as.numeric(.data$loao_min_top25_retained),
      animal_bootstrap_severe_instability = as.logical(
        .data$bootstrap_severe_instability
      ),
      influential_animal = as.logical(.data$influential_animal),
      strict_nonspatial_detail = as.character(
        .data$strict_nonspatial_sensitivity_status
      ),
      sensitivity_architecture_status = as.character(
        .data$primary_architecture_status
      )
    )
}

block_robustness <- read_optional_csv(microglia_inputs$block_readiness)
if (nrow(block_robustness)) {
  block_robustness <- block_robustness |>
    dplyr::transmute(
      dataset = "microglia",
      wgcna_entity_level = "supermodule",
      wgcna_entity_id = as.character(.data$SupermoduleID),
      loao_status = as.character(.data$leave_one_animal_out_support),
      animal_cluster_bootstrap_support = as.character(
        .data$animal_cluster_bootstrap_support
      ),
      strict_nonspatial_detail = as.character(.data$strict_nonspatial_support),
      sensitivity_architecture_status = as.character(
        .data$higher_order_block_classification
      ),
      cut_height_independent_support = as.logical(
        .data$cut_height_independent_support
      )
    )
}
robustness_detail <- dplyr::bind_rows(module_robustness, block_robustness)
if (nrow(robustness_detail)) {
  gww_assert_unique(
    robustness_detail,
    c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
    "Microglia robustness detail"
  )
}

robustness <- effects |>
  dplyr::distinct(
    .data$dataset, .data$wgcna_entity_level, .data$wgcna_entity_id
  ) |>
  dplyr::left_join(
    preservation,
    by = c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
    relationship = "one-to-one"
  ) |>
  dplyr::left_join(
    stage13,
    by = c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
    relationship = "one-to-one"
  ) |>
  dplyr::left_join(
    robustness_detail,
    by = c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
    relationship = "one-to-one"
  ) |>
  dplyr::mutate(
    animal_instability_flag = dplyr::case_when(
      .data$dataset != "microglia" ~ NA,
      .data$wgcna_entity_level == "module" ~
        .data$loao_status %in% c("animal_sensitive", "fail") |
        .data$animal_bootstrap_severe_instability %in% TRUE |
        .data$influential_animal %in% TRUE,
      .data$wgcna_entity_level == "supermodule" ~
        .data$loao_status %in% c("animal_sensitive", "fail") |
        grepl("severe|catastrophic", tolower(dplyr::coalesce(
          .data$animal_cluster_bootstrap_support, ""
        ))),
      TRUE ~ NA
    ),
    animal_instability_reason = dplyr::case_when(
      .data$dataset != "microglia" ~ "leave-one-animal-out summary unavailable",
      .data$animal_instability_flag %in% TRUE ~ paste(
        "loao_status=", .data$loao_status,
        "; bootstrap=", .data$animal_cluster_bootstrap_support,
        "; influential_animal=", .data$influential_animal,
        sep = ""
      ),
      TRUE ~ "no strong animal-instability flag in available readiness summaries"
    )
  )

matches_for_join <- program_matches |>
  dplyr::transmute(
    dataset = as.character(.data$dataset),
    wgcna_entity_level = as.character(.data$entity_level),
    wgcna_entity_id = as.character(.data$entity_id),
    biological_program = as.character(.data$biological_program),
    matching_annotation_field = as.character(.data$annotation_field),
    matching_annotation_token = as.character(.data$annotation_token),
    matching_rationale = as.character(.data$mapping_rationale),
    program_module_match_rule = as.character(.data$program_module_match_rule),
    module_annotation_evidence_role = as.character(
      .data$module_annotation_evidence_role
    )
  )

matched_evidence <- gsea_evidence |>
  dplyr::inner_join(
    matches_for_join,
    by = c("dataset", "biological_program"),
    relationship = "many-to-many"
  )
unmatched_program_evidence <- gsea_evidence |>
  dplyr::anti_join(
    matches_for_join |>
      dplyr::distinct(.data$dataset, .data$biological_program),
    by = c("dataset", "biological_program")
  )

local_long <- matched_evidence |>
  dplyr::filter(.data$comparison_scope == "local_local") |>
  dplyr::inner_join(
    effects |>
      dplyr::filter(
        .data$wgcna_analysis_tier == "exploratory_spatial_localization",
        .data$wgcna_effect_scope == "within_spatial_unit"
      ),
    by = c(
      "dataset", "contrast", "wgcna_entity_level", "wgcna_entity_id",
      "gsea_spatial_unit" = "wgcna_spatial_unit"
    ),
    relationship = "many-to-one"
  ) |>
  dplyr::mutate(wgcna_spatial_unit = .data$gsea_spatial_unit)

global_long <- matched_evidence |>
  dplyr::filter(.data$comparison_scope == "global_cross_spatial") |>
  dplyr::inner_join(
    effects |>
      dplyr::filter(
        .data$wgcna_analysis_tier %in% c(
          "primary_wgcna_global", "secondary_contextual_global"
        ),
        .data$wgcna_effect_scope == "spatial_adjusted_global",
        .data$wgcna_spatial_unit == "global_spatial_adjusted"
      ) |>
      dplyr::select(-"wgcna_spatial_unit"),
    by = c(
      "dataset", "contrast", "wgcna_entity_level", "wgcna_entity_id"
    ),
    relationship = "many-to-one"
  ) |>
  dplyr::mutate(wgcna_spatial_unit = "global_spatial_adjusted")

long <- dplyr::bind_rows(local_long, global_long)
gww_assert_unique(
  long,
  c("gsea_evidence_id", "wgcna_entity_level", "wgcna_entity_id"),
  "GSEA-WGCNA comparison rows before contextual joins"
)
pre_context_rows <- nrow(long)
unmatched_scope_pairs <- matched_evidence |>
  dplyr::select(dplyr::all_of(c(
    "dataset", "comparison_scope", "gsea_evidence_id",
    "wgcna_entity_level", "wgcna_entity_id"
  ))) |>
  dplyr::anti_join(
    long |>
      dplyr::select(dplyr::all_of(c(
        "dataset", "comparison_scope", "gsea_evidence_id",
        "wgcna_entity_level", "wgcna_entity_id"
      ))),
    by = c(
      "dataset", "comparison_scope", "gsea_evidence_id",
      "wgcna_entity_level", "wgcna_entity_id"
    )
  )

effect_thresholds <- gww_effect_scale_thresholds(
  effects |>
    dplyr::transmute(
      dataset,
      entity_level = .data$wgcna_entity_level,
      analysis_tier = .data$wgcna_analysis_tier,
      effect_scope = .data$wgcna_effect_scope,
      estimate = .data$wgcna_estimate,
      model_valid = .data$wgcna_model_valid,
      contrast
    )
)

overlap_for_join <- overlap |>
  dplyr::transmute(
    gsea_evidence_id,
    dataset,
    wgcna_entity_id = as.character(.data$entity_id),
    n_universe,
    n_module,
    n_leading_edge,
    n_overlap,
    overlap_proteins,
    Jaccard,
    overlap_fisher_p_value = .data$fisher_p_value,
    overlap_FDR,
    overlap_family_id,
    overlap_family_size,
    overlap_universe,
    overlap_test_scope
  )

long <- long |>
  dplyr::left_join(
    interaction,
    by = c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
    relationship = "many-to-one"
  ) |>
  dplyr::left_join(
    robustness,
    by = c("dataset", "wgcna_entity_level", "wgcna_entity_id"),
    relationship = "many-to-one"
  ) |>
  dplyr::left_join(
    effect_thresholds,
    by = c(
      "dataset", "wgcna_entity_level" = "entity_level",
      "wgcna_analysis_tier" = "analysis_tier",
      "wgcna_effect_scope" = "effect_scope"
    ),
    relationship = "many-to-one"
  ) |>
  dplyr::left_join(
    overlap_for_join,
    by = c("gsea_evidence_id", "dataset", "wgcna_entity_id"),
    relationship = "many-to-one"
  ) |>
  dplyr::mutate(
    overlap_supermodule_exclusion = dplyr::if_else(
      .data$wgcna_entity_level == "supermodule",
      "not_tested_to_avoid_correlated_member_module_union_hypothesis",
      NA_character_
    ),
    evidence_source_family =
      "ranked_GSEA;WGCNA_group_effect;WGCNA_module_annotation;WGCNA_overlap",
    evidence_independence_note = paste(
      "ranked GSEA, compareGO and program summaries are one inferential lineage;",
      "WGCNA module annotation is identity only; overlaps are dependent",
      "convergence evidence; RES-CON, SUS-CON and SUS-RES share groups"
    ),
    legacy_overlap_role =
      "04_wgcna_de_gsea_overlap retained as context only; not used for matching, classification or FDR"
  )
if (nrow(long) != pre_context_rows) {
  stop("Contextual joins multiplied comparison rows.", call. = FALSE)
}

long_first <- gww_classify_concordance(long)
long_second <- gww_classify_concordance(long)
if (!identical(long_first$concordance_class, long_second$concordance_class) ||
    !identical(long_first$concordance_class_reason,
               long_second$concordance_class_reason)) {
  stop("Concordance classification is not deterministic.", call. = FALSE)
}
long <- long_first
gww_validate_local_global_semantics(long)
gww_assert_unique(
  long,
  c("gsea_evidence_id", "wgcna_entity_level", "wgcna_entity_id"),
  "Final GSEA-WGCNA concordance rows"
)

# Exact source-value traces: no inferential value may be silently rewritten.
gsea_source <- gsea_terms |>
  dplyr::transmute(
    gsea_source_key = paste(.data$dataset, .data$Comparison, .data$ID, sep = "|"),
    source_NES = as.numeric(.data$NES),
    source_raw_p = as.numeric(.data$pvalue),
    source_GSEA_FDR = as.numeric(.data$p.adjust)
  ) |>
  dplyr::distinct()
gww_assert_unique(gsea_source, "gsea_source_key", "GSEA source-value trace")
gsea_index <- match(long$gsea_source_key, gsea_source$gsea_source_key)
if (anyNA(gsea_index) ||
    !identical(long$NES, gsea_source$source_NES[gsea_index]) ||
    !identical(long$raw_p, gsea_source$source_raw_p[gsea_index]) ||
    !identical(long$GSEA_FDR, gsea_source$source_GSEA_FDR[gsea_index])) {
  stop("A selected GSEA NES/p/FDR value differs from its canonical source row.",
       call. = FALSE)
}

wgcna_source <- wgcna |>
  dplyr::transmute(
    wgcna_source_key = as.character(.data$source_key),
    source_estimate = as.numeric(.data$estimate),
    source_SE = as.numeric(.data$SE),
    source_CI_low = as.numeric(.data$CI_low),
    source_CI_high = as.numeric(.data$CI_high),
    source_p_value = as.numeric(.data$p_value),
    source_tier_specific_fdr = as.numeric(.data$tier_specific_fdr),
    source_family_id = as.character(.data$tier_specific_family_id),
    source_family_size = as.integer(.data$tier_specific_family_size)
  )
gww_assert_unique(wgcna_source, "wgcna_source_key", "WGCNA source-value trace")
wgcna_index <- match(long$wgcna_source_key, wgcna_source$wgcna_source_key)
if (anyNA(wgcna_index) ||
    !identical(long$wgcna_estimate, wgcna_source$source_estimate[wgcna_index]) ||
    !identical(long$wgcna_SE, wgcna_source$source_SE[wgcna_index]) ||
    !identical(long$wgcna_CI_low, wgcna_source$source_CI_low[wgcna_index]) ||
    !identical(long$wgcna_CI_high, wgcna_source$source_CI_high[wgcna_index]) ||
    !identical(long$wgcna_p_value, wgcna_source$source_p_value[wgcna_index]) ||
    !identical(
      long$wgcna_tier_specific_fdr,
      wgcna_source$source_tier_specific_fdr[wgcna_index]
    ) ||
    !identical(
      long$wgcna_tier_specific_family_id,
      wgcna_source$source_family_id[wgcna_index]
    ) ||
    !identical(
      long$wgcna_tier_specific_family_size,
      wgcna_source$source_family_size[wgcna_index]
    )) {
  stop("A WGCNA estimate/SE/CI/p/FDR/family value differs from Stage 07.",
       call. = FALSE)
}

if (any(long$contrast == "SUS - RES" & !is.na(long$gsea_direction_sign) &
        sign(long$NES) != long$gsea_direction_sign)) {
  stop("SUS - RES GSEA direction semantics are inconsistent.", call. = FALSE)
}
if (any(!long$wgcna_independent_hypothesis %in% TRUE) ||
    any(long$wgcna_claim_entity_role == "compatibility_alias")) {
  stop("A non-independent or singleton compatibility alias entered concordance.",
       call. = FALSE)
}

summary <- gww_program_summary(long)
patterns <- gww_adaptive_pattern_summary(long)

long <- long |>
  dplyr::arrange(
    .data$dataset, .data$contrast, .data$comparison_scope,
    .data$gsea_spatial_unit, .data$biological_program,
    .data$wgcna_entity_level, .data$wgcna_entity_id
  )
summary <- summary |>
  dplyr::arrange(
    .data$dataset, .data$contrast, .data$comparison_scope,
    .data$biological_program
  )
patterns <- patterns |>
  dplyr::arrange(.data$dataset, .data$biological_program)
overlap <- overlap |>
  dplyr::arrange(
    .data$dataset, .data$contrast, .data$gsea_spatial_unit,
    .data$overlap_FDR, .data$biological_program, .data$entity_id
  )

audit_row <- function(audit_type, item, status = "PASS", n = NA_integer_,
                      dataset = NA_character_, contrast = NA_character_,
                      spatial_unit = NA_character_, detail = NA_character_,
                      source_file = NA_character_) {
  data.frame(
    audit_type = audit_type,
    dataset = dataset,
    contrast = contrast,
    spatial_unit = spatial_unit,
    item = item,
    status = status,
    n = as.integer(n),
    detail = detail,
    source_file = source_file,
    stringsAsFactors = FALSE
  )
}

hash_after <- source_hashes(protected_sources)
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
  stop("A protected source changed during downstream integration.",
       call. = FALSE)
}

concordance_classes <- c(
  "FDR_supported_concordance", "concordant_imprecise",
  "weak_or_near_zero_module_support", "discordant", "animal_sensitive",
  "unresolved"
)
adaptive_pattern_classes <- c(
  "RES_specific_adaptive_candidate", "SUS_specific_candidate",
  "opposing_RES_vs_SUS", "graded_stress_response",
  "shared_stress_response", "unresolved"
)

audit <- dplyr::bind_rows(
  lapply(seq_len(nrow(input_status)), function(i) {
    audit_row(
      "source_file", input_status$input_name[[i]],
      status = if (input_status$status[[i]] == "present") "PASS" else "INFO",
      detail = paste(
        "required=", input_status$required[[i]],
        "; status=", input_status$status[[i]], sep = ""
      ),
      source_file = input_status$source_file[[i]]
    )
  }),
  lapply(seq_len(nrow(hash_audit)), function(i) {
    audit_row(
      "protected_source_hash", hash_audit$input_name[[i]],
      status = if (hash_audit$unchanged[[i]]) "PASS" else "FAIL",
      detail = paste(
        "before=", hash_audit$source_md5_before[[i]],
        "; after=", hash_audit$source_md5_after[[i]], sep = ""
      ),
      source_file = hash_audit$source_file[[i]]
    )
  }),
  list(
    audit_row(
      "matching_rule", "program_module_matching", detail = paste(
        "Exact normalized Stage 07 module_program_primary/secondary token",
        "to config/gsea_wgcna_program_module_mapping.csv; labels are identity",
        "annotations, not evidence for a group effect."
      )
    ),
    audit_row(
      "matching_rule", "local_local", detail =
        "Local GSEA is paired only to within_spatial_unit WGCNA endpoints with the same normalized spatial unit."
    ),
    audit_row(
      "matching_rule", "global_cross_spatial", detail =
        "Global WGCNA is paired only to GSEA programs recurring in one direction across at least two spatial units."
    ),
    audit_row(
      "lineage", "ranked_GSEA", detail = paste(
        "GSEA, compareGO and program-summary representations are one inferential",
        "lineage and are never counted as independent confirmation."
      )
    ),
    audit_row(
      "exclusion", "singleton_compatibility_aliases", n = 0L,
      detail = "Canonical Stage 07 handoff excludes aliases; assertion passed."
    ),
    audit_row(
      "exclusion", "supermodule_overlap_tests", n = sum(
        program_matches$entity_level == "supermodule"
      ), detail = paste(
        "Genuine supermodules remain independent WGCNA concordance endpoints,",
        "but no Fisher test is run on correlated unions of member modules."
      )
    ),
    audit_row(
      "exclusion", "legacy_saturated_overlap", detail = paste(
        "04_wgcna_de_gsea_overlap is context only and is not used for matching,",
        "classification or overlap FDR."
      ), source_file = normalizePath(
        microglia_inputs$legacy_overlap, winslash = "/", mustWork = FALSE
      )
    ),
    audit_row(
      "exclusion", "gsea_evidence_without_curated_module_match",
      n = nrow(unmatched_program_evidence), detail = paste(
        "GSEA evidence rows without an exact checked-in program-to-module",
        "annotation match are retained in the source lineage but excluded",
        "from concordance comparisons rather than matched by label similarity."
      )
    ),
    audit_row(
      "exclusion", "matched_pairs_without_compatible_wgcna_endpoint",
      n = nrow(unmatched_scope_pairs), detail = paste(
        "Curated program-entity pairs lacking the required exact local spatial",
        "endpoint or the prespecified spatial-adjusted global endpoint are",
        "excluded. This includes soma route units such as ca1_sp when the",
        "WGCNA local endpoint is region-only ca1."
      )
    ),
    audit_row(
      "join_cardinality", "contextual_joins", n = nrow(long),
      detail = paste0("rows_before=", pre_context_rows, "; rows_after=", nrow(long))
    ),
    audit_row(
      "classification_rule", "low_n_effect_geometry", detail = paste(
        "Near-zero boundary=q25(abs estimate); concordant_imprecise requires",
        "estimate above q25, CI reaching median(abs estimate) in the GSEA",
        "direction, valid/stable model, no animal-instability flag, and q>0.05."
      )
    )
  ),
  overlap |>
    dplyr::distinct(
      .data$dataset, .data$contrast, .data$gsea_spatial_unit,
      .data$overlap_family_id, .data$overlap_family_size
    ) |>
    dplyr::transmute(
      audit_type = "overlap_family",
      dataset = .data$dataset,
      contrast = .data$contrast,
      spatial_unit = .data$gsea_spatial_unit,
      item = .data$overlap_family_id,
      status = "PASS",
      n = as.integer(.data$overlap_family_size),
      detail = paste(
        "BH within dataset x formal contrast x GSEA spatial unit over all",
        "prespecified biologically matched program x module hypotheses."
      ),
      source_file = normalizePath(
        gsea_terms_file, winslash = "/", mustWork = TRUE
      )
    ),
  long |>
    dplyr::count(.data$dataset, .data$contrast, .data$concordance_class) |>
    tidyr::complete(
      dataset, contrast,
      concordance_class = concordance_classes,
      fill = list(n = 0L)
    ) |>
    dplyr::transmute(
      audit_type = "concordance_class_count",
      dataset = .data$dataset,
      contrast = .data$contrast,
      spatial_unit = NA_character_,
      item = .data$concordance_class,
      status = "PASS",
      n = as.integer(.data$n),
      detail = "Deterministic class count.",
      source_file = NA_character_
    ),
  patterns |>
    dplyr::count(.data$adaptive_resilience_pattern) |>
    tidyr::complete(
      adaptive_resilience_pattern = adaptive_pattern_classes,
      fill = list(n = 0L)
    ) |>
    dplyr::transmute(
      audit_type = "adaptive_pattern_count",
      dataset = NA_character_,
      contrast = NA_character_,
      spatial_unit = NA_character_,
      item = .data$adaptive_resilience_pattern,
      status = "PASS",
      n = as.integer(.data$n),
      detail = "Joint three-contrast pattern count; contrasts share groups.",
      source_file = NA_character_
    )
)

matching_audit <- program_matches |>
  dplyr::select(dplyr::all_of(c(
    "dataset", "entity_level", "entity_id", "display_label",
    "cleaned_biological_label", "biological_program", "annotation_field",
    "annotation_token", "mapping_rationale", "program_module_match_rule",
    "module_annotation_evidence_role"
  ))) |>
  dplyr::arrange(
    .data$dataset, .data$entity_level, .data$entity_id,
    .data$biological_program
  )

for (filename in output_files) {
  # All expected files are written explicitly below. This loop documents the
  # complete contract and guards accidental duplicate names.
  if (anyDuplicated(output_files)) stop("Duplicated output filename contract.")
}
invisible(write_integration_table(long, paths, "gsea_wgcna_concordance_long.csv"))
invisible(write_integration_table(summary, paths, "gsea_wgcna_concordance_summary.csv"))
invisible(write_integration_table(patterns, paths, "adaptive_resilience_pattern_summary.csv"))
invisible(write_integration_table(
  overlap, paths, "program_specific_leading_edge_module_overlap.csv"
))
invisible(write_integration_table(audit, paths, "gsea_wgcna_concordance_audit.csv"))
invisible(write_integration_table(
  matching_audit, paths, "program_module_matching_audit.csv"
))
write_csv_safe(input_status, file.path(paths$reports, "input_status.csv"))
write_csv_safe(hash_audit, file.path(paths$reports, "protected_source_hash_audit.csv"))
write_csv_safe(audit, file.path(paths$reports, "gsea_wgcna_concordance_audit.csv"))

class_counts <- tibble::tibble(concordance_class = concordance_classes) |>
  dplyr::left_join(
    long |> dplyr::count(.data$concordance_class),
    by = "concordance_class",
    relationship = "one-to-one"
  ) |>
  dplyr::mutate(n = dplyr::coalesce(.data$n, 0L))
pattern_counts <- tibble::tibble(
  adaptive_resilience_pattern = adaptive_pattern_classes
) |>
  dplyr::left_join(
    patterns |> dplyr::count(.data$adaptive_resilience_pattern),
    by = "adaptive_resilience_pattern",
    relationship = "one-to-one"
  ) |>
  dplyr::mutate(n = dplyr::coalesce(.data$n, 0L))
readme_lines <- c(
  "# GSEA/GO <-> WGCNA concordance under low biological n",
  "",
  paste0("Generated: ", format(Sys.time(), tz = "Europe/Berlin", usetz = TRUE)),
  "",
  "## Interpretation boundary",
  "",
  "This is a purely downstream integration of existing ranked-GSEA and WGCNA outputs. It does not refit WGCNA, recalculate any source p-value/FDR, loosen any threshold, alter a claim gate, or create a combined significance statistic.",
  "",
  "GSEA, compareGO and program-summary views derived from the same ranked analysis are one inferential lineage. WGCNA module annotations establish biological matching only; they are not evidence for a group effect. RES-CON, SUS-CON and SUS-RES share animals/groups and are not independent replications.",
  "",
  "## Matching and multiple testing",
  "",
  "Local GSEA is compared only with the same local WGCNA spatial unit. Spatial-adjusted global WGCNA is compared only with a program recurring in one direction in at least two spatial units. Program-module identity uses the checked-in exact annotation registry. Fisher tests use the canonical dataset-specific WGCNA ProteinGroupID universe and are BH-adjusted within dataset x formal contrast x GSEA spatial unit over all prespecified matched program x module hypotheses.",
  "",
  "Genuine multimodule supermodules remain independent WGCNA endpoints in the concordance table, but no Fisher test is run on unions of their member modules because that would duplicate correlated module hypotheses. The legacy 04_wgcna_de_gsea_overlap output is retained as context only.",
  "",
  "## Low-n classes",
  "",
  "The near-zero boundary is the lower quartile of absolute WGCNA estimates within dataset x entity level x analysis tier x effect scope. Concordant-imprecise additionally requires a point estimate above that boundary, a CI reaching the median absolute effect in the GSEA-supported direction, a valid/stable model, and no existing strong animal-instability flag. Mixed GSEA direction is unresolved. WGCNA FDR support always uses the unchanged tier-specific FDR <= 0.05.",
  "",
  "Class counts:",
  paste0("- ", class_counts$concordance_class, ": ", class_counts$n),
  "",
  "Adaptive-pattern counts:",
  paste0("- ", pattern_counts$adaptive_resilience_pattern, ": ", pattern_counts$n),
  "",
  "## Limitations",
  "",
  "Biological n is nine animals total (three per group). Wide CIs are expected and remain visible. Pattern summaries require recurrent GSEA direction plus a stable direction-concordant global module endpoint for every contrast used. A significant/non-significant asymmetry is never interpreted as a difference. Microglia refers to a microglia-enriched ROI and not purified microglia. Neuronal datasets lack the dedicated Stage 12 leave-one-animal-out readiness summary, so that field is explicitly unavailable rather than imputed.",
  ""
)
dir_create(paths$reports)
writeLines(readme_lines, file.path(paths$reports, "README.md"), useBytes = TRUE)

write_integration_manifest(
  paths,
  c(required_inputs, optional_inputs),
  list(
    tables = paths$tables,
    source_data = paths$source_data,
    report = file.path(paths$reports, "README.md"),
    audit = file.path(paths$reports, "gsea_wgcna_concordance_audit.csv")
  ),
  list(
    datasets = datasets,
    gsea_fdr_threshold = 0.05,
    wgcna_fdr_threshold = 0.05,
    global_gsea_min_recurrent_units = 2L,
    overlap_family = "dataset x formal contrast x GSEA spatial unit",
    overlap_entities = "canonical modules only",
    near_zero_rule = "within-tier q25(abs estimate)",
    meaningful_ci_reference = "within-tier median(abs estimate)"
  ),
  paste(
    "Additive downstream integration only. Source ranked-GSEA and Stage 07",
    "inferential values are copied exactly; source hashes are protected.",
    "No combined p-value and no supermodule-union overlap test."
  )
)

message("GSEA/WGCNA concordance complete: ", paths$tables)
