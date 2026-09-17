#!/usr/bin/env Rscript
#
#
# WHAT THIS SCRIPT DOES
#   Joins frozen Stage 01 WGCNA module membership to frozen manifest-selected
#   protein-level differential abundance, to the canonical Stage 07 module-level
#   inferential handoff, and to the existing GSEA/WGCNA leading-edge overlap.
#   It then derives interpretable ranks and boolean flags and exports tables for
#   manual inspection.
#
# WHAT THIS SCRIPT DOES NOT DO
#   * No WGCNA is recomputed. kME, is_core_kME_0.6, is_top_hub_25,
#     GeneSignificanceP and GeneSignificanceFDR are copied unchanged.
#   * No differential abundance is recomputed. log2FC, p and BH FDR are copied
#     unchanged from the manifest-selected Stage 02 mapped contrast files.
#   * No module-level model is refitted. Estimates and tier-specific FDRs are
#     copied unchanged from the Stage 07 handoff.
#   * No frozen state, module identity, module label, Stage 01-13 scientific
#     output or manuscript result is modified.
#   * No weighted composite "importance score" is produced. Every criterion
#     stays a separately inspectable boolean.
#
# INTERPRETATION GUARDRAILS
#   * GeneSignificanceP/FDR is an OMNIBUS CON/RES/SUS ANOVA. It is unsigned and
#     is not a SUS-vs-RES phenotype correlation.
#   * GSEA leading-edge membership is NOT independent of the DA columns: the
#     ranked GSEA that produced it ranks on the moderated t statistic of the
#     same limma fit.
#   * Hippocampal spatial contexts are repeated measurements from the same
#     animals. Agreement across them is spatial consistency, not replication.
#
# USAGE
#   Rscript analysis/08_integration/build_candidate_protein_shortlist.R
#   Rscript analysis/08_integration/build_candidate_protein_shortlist.R --dataset microglia
#   Rscript analysis/08_integration/build_candidate_protein_shortlist.R --dataset neuron_neuropil --module WGCNA_m02
#   Rscript analysis/08_integration/build_candidate_protein_shortlist.R --dry-run
# Script: analysis/08_integration/build_candidate_protein_shortlist.R
# Stage: integration
# Scope: per_dataset_and_global
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_modules_long.csv; results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_inferential_handoff.csv; data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv; optional results/tables/10_biological_integration/gsea_wgcna_concordance/global/program_specific_leading_edge_module_overlap.csv
# Produces: results/tables/10_biological_integration/wgcna_candidate_protein_shortlist/<dataset>/wgcna_candidate_proteins_all.csv; results/tables/10_biological_integration/wgcna_candidate_protein_shortlist/<dataset>/wgcna_candidate_proteins_long.csv; results/tables/10_biological_integration/wgcna_candidate_protein_shortlist/<dataset>/wgcna_candidate_proteins_shortlist.csv; +5 more
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Additive downstream candidate-protein shortlist for frozen WGCNA modules.

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/enrichment/enrichment_io.R")
source("R/statistics/sus_res_spatial_dap_atlas_utils.R")
source("R/statistics/wgcna_candidate_protein_utils.R")
source("R/utilities/xlsx_package_utils.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

SCRIPT_ID <- "analysis/08_integration/build_candidate_protein_shortlist.R"
SUBSTEP <- "wgcna_candidate_protein_shortlist"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

cli <- integration_cli(default_dataset = "all")
cli_value <- function(flag, default = "") {
  hit <- which(cli$args == flag)
  if (!length(hit) || hit[[1]] == length(cli$args)) return(default)
  cli$args[[hit[[1]] + 1L]]
}
module_focus <- trimws(cli_value("--module", ""))

# ------------------------------------------------------------------- inputs

membership_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules",
               "WGCNA_modules_long.csv")
}
handoff_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset,
               "WGCNA_inferential_handoff.csv")
}
leading_edge_path <- function() {
  path_results("tables", "10_biological_integration", "gsea_wgcna_concordance",
               "global", "program_specific_leading_edge_module_overlap.csv")
}

datasets <- integration_datasets(cli$dataset)

if (isTRUE(cli$dry_run)) {
  inputs <- list()
  for (dataset in datasets) {
    inputs[[paste0("module_membership__", dataset)]] <- membership_path(dataset)
    inputs[[paste0("stage07_inferential_handoff__", dataset)]] <- handoff_path(dataset)
    inputs[[paste0("clusterprofiler_manifest__", dataset)]] <-
      canonical_clusterprofiler_manifest_path(dataset)
  }
  inputs[["gsea_wgcna_leading_edge_overlap"]] <- leading_edge_path()
  cat("[DRY-RUN] WGCNA candidate-protein shortlist; no outputs written.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Would write to results/{tables,source_data,reports,logs}/",
      "10_biological_integration/", SUBSTEP, "/<dataset>/ and /global/.\n", sep = "")
  quit(save = "no", status = 0L)
}

read_required <- function(path, label) {
  if (!file.exists(path)) {
    stop("missing_required_input: ", label, ": ", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE, guess_max = Inf)
}

relpath <- function(path) relative_to(normalizePath(path, winslash = "/", mustWork = FALSE))

# --------------------------------------------------- canonical input loaders

# Frozen Stage 01 membership. One row per ProteinGroupID; module identity,
# kME and the frozen hub/core flags are consumed verbatim.
load_membership <- function(dataset) {
  path <- membership_path(dataset)
  raw <- read_required(path, "frozen WGCNA module membership")
  keep <- c(
    "dataset", "ModuleID", "ModuleLegacyID", "ModuleColor", "ModuleLabel_Final",
    "final_label", "primary_label", "ProteinGroupID", "GeneSymbol",
    "GeneSymbols", "UniProt", "RepresentativeUniProt", "representative_gene_symbol",
    "mapping_status", "mapping_strategy", "gene_level_claim_allowed",
    "protein_level_claim_allowed", "protein_group_ambiguity_class",
    "same_gene_group", "n_members_canonical", "n_gene_symbols",
    "n_unmapped_members", "kME", "abs_kME", "GeneSignificanceP",
    "GeneSignificanceFDR", "is_core_kME_0.6", "is_top_hub_25"
  )
  out <- raw[, intersect(keep, names(raw)), drop = FALSE]
  out$membership_source_file <- relpath(path)
  out
}

# Canonical Stage 07 module-level three-contrast handoff. Selected with the
# frozen tier contract: primary global SUS - RES uses primary_wgcna_global;
# the other two named contrasts use their own tier rows.
load_module_effects <- function(dataset) {
  path <- handoff_path(dataset)
  raw <- read_required(path, "Stage 07 WGCNA inferential handoff")
  out <- raw %>%
    filter(
      .data$entity_level == "module",
      .data$contrast %in% wcp_contrast_levels(),
      .data$spatial_scope == "spatial_adjusted_global"
    ) %>%
    transmute(
      dataset = .data$dataset,
      ModuleID = .data$module_id,
      contrast = .data$contrast,
      module_supermodule_id = .data$module_supermodule_id,
      module_supermodule_label = .data$module_supermodule_label,
      module_display_label = .data$display_label,
      module_analysis_tier = .data$analysis_tier,
      module_estimate = suppressWarnings(as.numeric(.data$estimate)),
      module_CI_low = suppressWarnings(as.numeric(.data$CI_low)),
      module_CI_high = suppressWarnings(as.numeric(.data$CI_high)),
      module_p_value = suppressWarnings(as.numeric(.data$p_value)),
      module_tier_specific_fdr = suppressWarnings(as.numeric(.data$tier_specific_fdr)),
      module_tier_specific_family_id = .data$tier_specific_family_id,
      module_tier_specific_family_size = .data$tier_specific_family_size,
      module_support_class = .data$support_class,
      module_statistical_support_status = .data$statistical_support_status,
      module_claim_gate = .data$claim_gate,
      module_model_valid = .data$model_valid,
      module_model_stability_status = .data$model_stability_status,
      module_direction = .data$direction,
      module_safe_interpretation = .data$safe_interpretation,
      module_effect_source_file = relpath(path),
      module_effect_source_key = .data$source_key
    )
  dup <- duplicated(out[c("dataset", "ModuleID", "contrast")])
  if (any(dup)) {
    stop("Stage 07 handoff yielded duplicate dataset x ModuleID x contrast rows; ",
         "the tier selection contract is not being honoured.", call. = FALSE)
  }
  out
}

# Manifest-selected protein-level DA, exactly as the manuscript Figure 3
# renderer resolves it. Values are copied; only the orientation multiplier
# recorded in the canonical filename convention is applied so that the sign
# matches the canonical contrast label.
load_protein_da <- function(dataset) {
  manifest_path <- canonical_clusterprofiler_manifest_path(dataset)
  if (!file.exists(manifest_path)) {
    stop("missing_required_input: clusterProfiler manifest: ", manifest_path,
         call. = FALSE)
  }
  manifest <- readr::read_csv(
    manifest_path, show_col_types = FALSE, progress = FALSE, guess_max = Inf,
    col_types = readr::cols(.default = readr::col_character())
  ) %>%
    filter(.data$result_type == "GSEA_GO",
           .data$route_category == "phenotype_within_unit") %>%
    mutate(contrast = wcp_contrast_from_comparison(.data$comparison)) %>%
    filter(.data$contrast %in% wcp_contrast_levels()) %>%
    distinct(.data$comparison, .keep_all = TRUE)
  if (!nrow(manifest)) {
    stop("missing_required_input: no phenotype-within-unit GSEA_GO manifest rows for ",
         dataset, call. = FALSE)
  }

  bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, , drop = FALSE]
    input <- sus_res_resolve_manifest_input(row$input_gene_file[[1]], dataset,
                                            repo_root())
    z <- read_required(input, "manifest-selected protein DA input")
    required <- c("ProteinGroupID", "log2fc", "pval", "padj")
    missing <- setdiff(required, names(z))
    if (length(missing)) {
      stop("Manifest-selected DA input is missing: ",
           paste(missing, collapse = ", "), " (", input, ")", call. = FALSE)
    }
    orientation <- wcp_comparison_orientation(row$comparison[[1]])
    multiplier <- orientation$formal_effect_multiplier[[1]]
    tibble(
      dataset = dataset,
      spatial_unit = row$route_unit[[1]],
      contrast = row$contrast[[1]],
      da_comparison = row$comparison[[1]],
      ProteinGroupID = as.character(z$ProteinGroupID),
      da_gene_symbol = if ("gene_symbol" %in% names(z)) as.character(z$gene_symbol) else NA_character_,
      log2FC = suppressWarnings(as.numeric(z$log2fc)) * multiplier,
      raw_p = suppressWarnings(as.numeric(z$pval)),
      BH_FDR = suppressWarnings(as.numeric(z$padj)),
      serialized_effect_definition = orientation$serialized_effect_definition[[1]],
      effect_sign_was_flipped = orientation$sign_was_flipped[[1]],
      da_source_file = relpath(input),
      da_manifest_input_hash = row$input_hash[[1]],
      da_manifest_file = relpath(manifest_path)
    )
  }))
}

# Existing GSEA/WGCNA program-specific leading-edge overlap, already expressed
# in ProteinGroupID space. Descriptive context only - see guardrails above.
load_leading_edge <- function() {
  path <- leading_edge_path()
  status <- read_csv_optional(
    path, dataset = "global", evidence_domain = "gsea_wgcna_concordance",
    input_type = "program_specific_leading_edge_module_overlap.csv",
    required = FALSE
  )
  if (is.null(status$data) || !nrow(status$data)) {
    return(list(
      data = tibble(dataset = character(), ModuleID = character(),
                    contrast = character(), ProteinGroupID = character(),
                    gsea_leading_edge_programs = character(),
                    gsea_leading_edge_source_file = character()),
      status = status$status
    ))
  }
  le <- status$data %>%
    filter(is.finite(suppressWarnings(as.numeric(.data$overlap_FDR))),
           suppressWarnings(as.numeric(.data$overlap_FDR)) <= wcp_fdr_threshold(),
           .data$contrast %in% wcp_contrast_levels()) %>%
    select(dataset = "dataset", ModuleID = "entity_id", contrast = "contrast",
           biological_program = "biological_program",
           overlap_proteins = "overlap_proteins") %>%
    tidyr::separate_rows("overlap_proteins", sep = ";") %>%
    mutate(ProteinGroupID = trimws(.data$overlap_proteins)) %>%
    filter(nzchar(.data$ProteinGroupID)) %>%
    group_by(.data$dataset, .data$ModuleID, .data$contrast, .data$ProteinGroupID) %>%
    summarise(
      gsea_leading_edge_programs = paste(sort(unique(.data$biological_program)),
                                         collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(gsea_leading_edge_source_file = relpath(path))
  list(data = le, status = status$status)
}

# ------------------------------------------------------------ per dataset

build_dataset <- function(dataset, leading_edge) {
  members <- load_membership(dataset)
  module_effects <- load_module_effects(dataset)
  da <- load_protein_da(dataset)

  members <- wcp_rank_module_members(members)
  members$clean_mapping <- wcp_clean_mapping(members)
  members$module_label <- members$ModuleLabel_Final %||% members$final_label

  # ---- long evidence table: protein x spatial unit x contrast
  da <- wcp_flag_large_effect(da)
  da$fdr_supported <- is.finite(da$BH_FDR) & da$BH_FDR <= wcp_fdr_threshold()

  long <- members %>%
    inner_join(da, by = c("dataset", "ProteinGroupID")) %>%
    left_join(module_effects, by = c("dataset", "ModuleID", "contrast")) %>%
    left_join(leading_edge, by = c("dataset", "ModuleID", "contrast", "ProteinGroupID")) %>%
    mutate(
      gsea_leading_edge = !is.na(.data$gsea_leading_edge_programs),
      direction_matches_module = wcp_direction_matches_module(
        .data$log2FC, .data$module_estimate
      ),
      contract_version = wcp_contract_version(),
      Source = SCRIPT_ID
    )

  # ---- spatial consistency, per contrast, computed on de-duplicated rows
  consistency <- wcp_spatial_consistency(
    long %>% distinct(.data$dataset, .data$ProteinGroupID, .data$contrast,
                      .data$spatial_unit, .keep_all = TRUE) %>%
      select("dataset", "ProteinGroupID", "contrast", "spatial_unit",
             "log2FC", "BH_FDR")
  )
  # Protein-level effect magnitude, thresholded within dataset x contrast.
  consistency <- wcp_flag_typical_large_effect(consistency)

  primary <- wcp_primary_contrast()

  # ---- per-protein SUS - RES rollup
  sus_res <- long %>%
    filter(.data$contrast == primary) %>%
    group_by(.data$dataset, .data$ModuleID, .data$ProteinGroupID) %>%
    summarise(
      sus_res_fdr05_any_context = any(.data$fdr_supported %in% TRUE),
      sus_res_n_contexts_fdr05 = sum(.data$fdr_supported %in% TRUE),
      sus_res_large_effect_any_context =
        any(.data$large_effect_within_context %in% TRUE),
      sus_res_n_contexts_large_effect =
        sum(.data$large_effect_within_context %in% TRUE),
      sus_res_min_BH_FDR = suppressWarnings(min(.data$BH_FDR, na.rm = TRUE)),
      sus_res_n_contexts_direction_matches_module =
        sum(.data$direction_matches_module %in% TRUE),
      sus_res_n_contexts_direction_comparable =
        sum(!is.na(.data$direction_matches_module)),
      sus_res_gsea_leading_edge = any(.data$gsea_leading_edge %in% TRUE),
      sus_res_gsea_leading_edge_programs = {
        hit <- sort(unique(stats::na.omit(.data$gsea_leading_edge_programs)))
        if (!length(hit)) NA_character_ else paste(hit, collapse = ";")
      },
      .groups = "drop"
    ) %>%
    mutate(sus_res_min_BH_FDR = ifelse(is.finite(.data$sus_res_min_BH_FDR),
                                       .data$sus_res_min_BH_FDR, NA_real_))

  # Strongest single spatial context for SUS - RES, for convenient inspection.
  best_context <- long %>%
    filter(.data$contrast == primary, is.finite(.data$log2FC)) %>%
    arrange(.data$dataset, .data$ProteinGroupID,
            desc(abs(.data$log2FC)), .data$spatial_unit) %>%
    group_by(.data$dataset, .data$ProteinGroupID) %>%
    slice_head(n = 1L) %>%
    ungroup() %>%
    transmute(
      dataset = .data$dataset, ProteinGroupID = .data$ProteinGroupID,
      sus_res_strongest_spatial_unit = .data$spatial_unit,
      sus_res_strongest_log2FC = .data$log2FC,
      sus_res_strongest_raw_p = .data$raw_p,
      sus_res_strongest_BH_FDR = .data$BH_FDR
    )

  sus_res_consistency <- consistency %>%
    filter(.data$contrast == primary) %>%
    select(-"contrast") %>%
    rename_with(~ paste0("sus_res_", .x),
                .cols = -c("dataset", "ProteinGroupID"))

  module_effect_primary <- module_effects %>%
    filter(.data$contrast == primary,
           .data$module_analysis_tier == "primary_wgcna_global") %>%
    select(-"contrast")

  gsea_any <- long %>%
    group_by(.data$dataset, .data$ProteinGroupID) %>%
    summarise(gsea_leading_edge_any_contrast = any(.data$gsea_leading_edge %in% TRUE),
              .groups = "drop")

  summary_tbl <- members %>%
    left_join(module_effect_primary, by = c("dataset", "ModuleID")) %>%
    left_join(sus_res, by = c("dataset", "ModuleID", "ProteinGroupID")) %>%
    left_join(sus_res_consistency, by = c("dataset", "ProteinGroupID")) %>%
    left_join(best_context, by = c("dataset", "ProteinGroupID")) %>%
    left_join(gsea_any, by = c("dataset", "ProteinGroupID")) %>%
    mutate(
      has_protein_level_da = !is.na(.data$sus_res_n_spatial_contexts_tested) &
        .data$sus_res_n_spatial_contexts_tested > 0L,
      sus_res_fdr05_any_context = .data$sus_res_fdr05_any_context %in% TRUE,
      sus_res_large_effect_any_context = .data$sus_res_large_effect_any_context %in% TRUE,
      sus_res_large_effect_typical = .data$sus_res_large_effect_typical %in% TRUE,
      sus_res_spatially_consistent = .data$sus_res_spatially_consistent %in% TRUE,
      sus_res_gsea_leading_edge = .data$sus_res_gsea_leading_edge %in% TRUE,
      gsea_leading_edge_any_contrast = .data$gsea_leading_edge_any_contrast %in% TRUE,
      sus_res_direction_matches_module = wcp_direction_matches_module(
        .data$sus_res_median_log2FC, .data$module_estimate
      ),
      contract_version = wcp_contract_version(),
      Source = SCRIPT_ID
    )

  summary_tbl <- wcp_assign_candidate_flags(summary_tbl)
  # Reason strings read from the canonical flag names.
  summary_tbl$gsea_leading_edge <- summary_tbl$sus_res_gsea_leading_edge
  summary_tbl$sus_res_direction_matches_module_flag <-
    summary_tbl$sus_res_direction_matches_module %in% TRUE
  summary_tbl$phenotype_network_class <- wcp_phenotype_network_class(summary_tbl)
  summary_tbl$candidate_tier <- wcp_candidate_tier(summary_tbl)
  summary_tbl$candidate_tier_all <- wcp_candidate_tier_all(summary_tbl)
  summary_tbl$candidate_reason <- wcp_candidate_reason(summary_tbl)
  summary_tbl <- wcp_order_candidates(summary_tbl)

  # Invariant: a claimable protein with SUS - RES FDR support must land in A1,
  # A2 or D. Only the frozen mapping contract may exclude one, and that case is
  # reported rather than silently dropped.
  unclassified <- wcp_unclassified_fdr_support(summary_tbl)
  if (nrow(unclassified)) {
    stop(
      "FDR-supported SUS - RES proteins left unclassified in ", dataset, ": ",
      nrow(unclassified), " (first: ", unclassified$ProteinGroupID[[1]],
      "). Every claimable FDR-supported protein must be Tier A1, A2 or D.",
      call. = FALSE
    )
  }

  list(long = long, summary = summary_tbl, consistency = consistency)
}

# ----------------------------------------------------------------- ordering

summary_column_order <- function(df) {
  preferred <- c(
    # identity / annotation
    "dataset", "ModuleID", "module_label", "ModuleLegacyID", "ModuleColor",
    "module_supermodule_id", "module_supermodule_label", "module_display_label",
    "ProteinGroupID", "GeneSymbol", "representative_gene_symbol",
    "RepresentativeUniProt", "UniProt", "mapping_status",
    "gene_level_claim_allowed", "protein_level_claim_allowed",
    "protein_group_ambiguity_class", "same_gene_group", "n_gene_symbols",
    "n_members_canonical", "n_unmapped_members", "clean_mapping",
    # candidate classification
    "candidate_tier", "candidate_tier_all", "phenotype_network_class",
    "candidate_reason", "is_candidate",
    "is_tier_A1", "is_tier_A2", "is_tier_B", "is_tier_C", "is_tier_D", "is_tier_A",
    # WGCNA network information
    "kME", "abs_kME", "abs_kME_rank_in_module",
    "abs_kME_rank_fraction_in_module", "n_module_members",
    "is_core_kME_0.6", "is_top_hub_25", "is_top5_module_representative",
    "is_top10_module_hub", "high_kME", "very_high_kME", "top5_hub", "top10_hub",
    # omnibus condition evidence
    "GeneSignificanceP", "GeneSignificanceFDR",
    # protein-level DA, SUS - RES rollup
    "has_protein_level_da", "sus_res_fdr05_any_context",
    "sus_res_n_contexts_fdr05", "sus_res_min_BH_FDR",
    "sus_res_large_effect_typical", "sus_res_typical_large_effect_threshold",
    "sus_res_large_effect_any_context", "sus_res_n_contexts_large_effect",
    "sus_res_strongest_spatial_unit", "sus_res_strongest_log2FC",
    "sus_res_strongest_raw_p", "sus_res_strongest_BH_FDR",
    # spatial consistency
    "sus_res_n_spatial_contexts_tested", "sus_res_n_spatial_contexts_positive",
    "sus_res_n_spatial_contexts_negative", "sus_res_majority_direction",
    "sus_res_n_matching_majority_direction",
    "sus_res_fraction_matching_majority_direction",
    "sus_res_n_spatial_contexts_fdr05", "sus_res_max_abs_log2FC",
    "sus_res_median_log2FC", "sus_res_median_abs_log2FC",
    "sus_res_spatially_consistent",
    # module-level inference and concordance
    "module_analysis_tier", "module_estimate", "module_CI_low", "module_CI_high",
    "module_p_value", "module_tier_specific_fdr",
    "module_tier_specific_family_id", "module_tier_specific_family_size",
    "module_support_class", "module_statistical_support_status",
    "module_claim_gate", "module_model_valid", "module_model_stability_status",
    "module_direction", "module_safe_interpretation",
    "sus_res_direction_matches_module",
    "sus_res_n_contexts_direction_matches_module",
    "sus_res_n_contexts_direction_comparable",
    # GSEA-derived context
    "gsea_leading_edge", "sus_res_gsea_leading_edge",
    "sus_res_gsea_leading_edge_programs", "gsea_leading_edge_any_contrast",
    # provenance
    "membership_source_file", "module_effect_source_file",
    "module_effect_source_key", "contract_version", "Source"
  )
  df[, c(intersect(preferred, names(df)), setdiff(names(df), preferred)),
     drop = FALSE]
}

# The long table is one row per protein x spatial unit x contrast, so roughly
# 280k rows.  It therefore carries ONLY the evidence that actually varies per
# row plus the identity needed to join back.  Protein-level annotation
# (GeneSymbol aside), module prose and constant file paths live in
# wgcna_candidate_proteins_all.csv and in the run manifest; repeating them here
# inflated the file to ~550 MB with no added information.
#
# Provenance is preserved in a stronger, compact form: `da_comparison` names the
# exact manifest comparison the numbers came from and `da_manifest_input_hash`
# pins that input file's content.  The resolved paths are in the run manifest.
long_column_order <- function(df) {
  keep <- c(
    # identity, sufficient to join back to wgcna_candidate_proteins_all.csv
    "dataset", "ModuleID", "ProteinGroupID", "GeneSymbol",
    # per-context protein-level differential abundance
    "spatial_unit", "contrast", "log2FC", "raw_p", "BH_FDR", "fdr_supported",
    "large_effect_within_context", "large_effect_abs_log2FC_threshold",
    "n_tested_in_effect_family", "serialized_effect_definition",
    "effect_sign_was_flipped",
    # module-level inference for the SAME contrast
    "module_analysis_tier", "module_estimate", "module_p_value",
    "module_tier_specific_fdr", "module_support_class", "module_claim_gate",
    "module_model_valid", "module_direction", "direction_matches_module",
    # GSEA-derived context
    "gsea_leading_edge", "gsea_leading_edge_programs",
    # compact provenance
    "da_comparison", "da_manifest_input_hash"
  )
  df[, intersect(keep, names(df)), drop = FALSE]
}

# ------------------------------------------------------------------- report

write_readme <- function(paths, scope, summary_tbl) {
  dict <- wcp_flag_dictionary()
  lines <- c(
    paste0("# WGCNA candidate-protein shortlist (", scope, ")"),
    "",
    paste0("Generated by `", SCRIPT_ID, "` (contract ", wcp_contract_version(), ")."),
    "",
    "## Scope",
    "",
    "Additive, descriptive downstream layer. No WGCNA, differential-abundance or",
    "module-level model is recomputed here; every statistic is copied unchanged from",
    "its frozen upstream artifact and carries a provenance column.",
    "",
    "There is deliberately **no weighted composite importance score**. Each criterion",
    "is a separate boolean so it can be inspected and rejected independently.",
    "",
    "## Candidate tiers",
    "",
    "Tiers are **not mutually exclusive**; `candidate_tier_all` lists every tier a",
    "protein satisfies and `candidate_tier` is only a convenience label.",
    "",
    "The three phenotype-linked classes (A1, A2, D) share the SAME protein-level",
    "evidence requirement and differ only in WGCNA topology. The precedence",
    "A1 > A2 > D > B > C is a reading order, **not** a ranking of evidence: a Tier D",
    "protein routinely carries a smaller FDR and a larger effect than a Tier A1 one.",
    "",
    "- **Tier A1 - FDR-supported phenotype-linked module hub**: unambiguous mapping,",
    "  SUS - RES BH FDR <= 0.05 in at least one spatial context, |kME| >= 0.60, and in",
    "  the module's top 10 by |kME|. Module-direction concordance is reported as an",
    "  additional flag rather than required, because no module-level SUS - RES effect",
    "  is FDR-supported in any dataset.",
    "- **Tier A2 - FDR-supported phenotype-linked module member**: the same evidence",
    "  and |kME| >= 0.60, but outside the module's top 10 hubs. A bona fide core",
    "  member by the frozen |kME| >= 0.60 definition. No top-25 requirement is",
    "  imposed; the frozen `is_top_hub_25` flag is reported separately.",
    "- **Tier B - network + phenotype candidate**: |kME| >= 0.60 AND a large",
    "  typical SUS - RES effect (`sus_res_large_effect_typical`: median |log2FC|",
    "  across spatial contexts in the top decile for this dataset x contrast).",
    "  Protein-level FDR support is not required. Spatial directional agreement is",
    "  NOT part of this definition - see the caveats below - and the per-context",
    "  effect flag is not used either, because 'top decile in any context' gives a",
    "  protein one chance per context and flags 35-65% of all proteins.",
    "- **Tier C - canonical module representative**: top 5 within its module by",
    "  |kME|, regardless of differential abundance.",
    "- **Tier D - FDR-supported phenotype-linked peripheral module member**:",
    "  SUS - RES BH FDR <= 0.05 but |kME| < 0.60, so strongly stress-responsive",
    "  proteins are not discarded for being peripheral to a module. Unlike A1/A2 this",
    "  class does not additionally require clean mapping - a pre-existing asymmetry",
    "  left unchanged; on the current frozen inputs every FDR-supported protein is",
    "  cleanly mapped, so it currently has no effect.",
    "",
    "Every claimable protein with SUS - RES BH FDR <= 0.05 is therefore Tier A1, A2 or",
    "D. The run fails loudly if any such protein is left unclassified.",
    "`is_tier_A` is retained only as a compatibility aggregate (A1 OR A2); A1/A2 carry",
    "the actual distinction and are what the workbook sheets use.",
    "",
    "## Evidence classes",
    "",
    "| Evidence class | Meaning | Independent of protein-level DA? |",
    "| --- | --- | --- |"
  )
  notes <- wcp_evidence_class_notes()
  lines <- c(lines, paste0(
    "| ", notes[["Evidence class"]], " | ", notes[["Meaning"]], " | ",
    notes[["Independent of protein-level DA?"]], " |"
  ))
  lines <- c(
    lines, "",
    "## Interpretation caveats",
    "",
    "- `GeneSignificanceP`/`GeneSignificanceFDR` is an **omnibus** CON/RES/SUS ANOVA.",
    "  It is unsigned and is not a SUS-vs-RES phenotype correlation.",
    "- GSEA leading-edge membership is **not** independent evidence from the DA",
    "  columns: that GSEA ranks on the moderated t statistic of the same limma fit.",
    "- Hippocampal spatial contexts are repeated measurements from the **same**",
    "  animals. Agreement across them is spatial consistency, not replication.",
    "  They are positively correlated, so their signs are **not** independent",
    "  Bernoulli trials: unanimous direction arises far more often than any",
    "  independence-based calculation would imply. No such sign probability is",
    "  computed, exported, or used as inferential evidence anywhere in this layer,",
    "  and spatial agreement does not by itself create a candidate tier - it is",
    "  descriptive prioritisation information only.",
    "- `direction_matches_module` is descriptive concordance between a protein and a",
    "  summary of the same proteins in the same animals, not independent support.",
    "",
    "## Field dictionary",
    "",
    "| Field | Evidence class | Meaning | Provenance |",
    "| --- | --- | --- | --- |",
    paste0("| `", dict$Field, "` | ", dict$`Evidence class`, " | ", dict$Meaning,
           " | ", dict$Provenance, " |"),
    "",
    "## Counts",
    "",
    paste0("- WGCNA proteins examined: ", nrow(summary_tbl)),
    paste0("- Cleanly mapped: ", sum(summary_tbl$clean_mapping %in% TRUE)),
    paste0("- SUS - RES FDR-supported proteins: ",
           sum(summary_tbl$sus_res_fdr05_any_context %in% TRUE),
           " (unclassified: ", nrow(wcp_unclassified_fdr_support(summary_tbl)), ")"),
    paste0("- Tier A1: ", sum(summary_tbl$is_tier_A1 %in% TRUE)),
    paste0("- Tier A2: ", sum(summary_tbl$is_tier_A2 %in% TRUE)),
    paste0("- Tier B: ", sum(summary_tbl$is_tier_B %in% TRUE)),
    paste0("- Tier C: ", sum(summary_tbl$is_tier_C %in% TRUE)),
    paste0("- Tier D: ", sum(summary_tbl$is_tier_D %in% TRUE)),
    paste0("- Any candidate flag: ", sum(summary_tbl$is_candidate %in% TRUE)),
    "",
    "## Files",
    "",
    "- `wgcna_candidate_proteins_all.csv` - one row per dataset x ModuleID x",
    "  ProteinGroupID, with all annotation, flags and tiers.",
    "- `wgcna_candidate_proteins_long.csv` - one row per protein x spatial unit x",
    "  contrast, carrying only the evidence that varies per row. Join back to",
    "  `_all.csv` on `dataset` + `ProteinGroupID` for annotation. Written per",
    "  dataset only; the cross-dataset scope would just concatenate these.",
    "- `wgcna_candidate_proteins_shortlist.csv` - rows satisfying >= 1 tier.",
    "- `wgcna_top10_per_module.csv` - top 10 per module by |kME|.",
    "- `wgcna_candidate_proteins_shortlist.xlsx` - the same content for review.",
    "",
    "## Inspecting one module",
    "",
    "```",
    paste0("Rscript ", SCRIPT_ID, " --dataset neuron_neuropil --module WGCNA_m02"),
    "```",
    "",
    "or filter the exported table directly:",
    "",
    "```r",
    'x <- readr::read_csv("wgcna_candidate_proteins_all.csv")',
    'subset(x, ModuleID == "WGCNA_m02" & !is.na(candidate_tier))',
    "```",
    ""
  )
  path <- file.path(paths$reports, "README.md")
  dir_create(dirname(path))
  writeLines(lines, path)
  path
}

# ------------------------------------------------------------------ workbook

write_workbook <- function(path, scope, summary_tbl, top10, by_dataset) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("NOTE openxlsx is not installed; skipping ", basename(path))
    return(NA_character_)
  }
  wb <- openxlsx::createWorkbook(creator = "analysis/08_integration/build_candidate_protein_shortlist.R")

  title_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, fontColour = "#1F2933", textDecoration = "bold", halign = "left", valign = "center")
  note_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#5B6770", textDecoration = "italic", wrapText = TRUE, valign = "top")
  header_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#1F2933", fgFill = "#E9EDF0", textDecoration = "bold", wrapText = TRUE, halign = "center", valign = "center", border = "Bottom", borderColour = "#5B6770", borderStyle = "thin")
  body_base <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, valign = "top")
  body_wrap <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, wrapText = TRUE, valign = "top")
  fdr_style <- openxlsx::createStyle(numFmt = "0.00E+00")
  num_style <- openxlsx::createStyle(numFmt = "0.000")
  fraction_style <- openxlsx::createStyle(numFmt = "0.0%")
  tier_a_fill <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#E6F0E8", textDecoration = "bold")
  tier_d_fill <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#E8F2FA")

  sci_cols <- c("GeneSignificanceP", "GeneSignificanceFDR", "sus_res_min_BH_FDR",
                "sus_res_strongest_raw_p", "sus_res_strongest_BH_FDR",
                "module_p_value", "module_tier_specific_fdr", "BH_FDR", "raw_p")
  num_cols <- c("kME", "abs_kME", "sus_res_strongest_log2FC", "sus_res_max_abs_log2FC",
                "sus_res_median_log2FC", "sus_res_median_abs_log2FC",
                "module_estimate", "module_CI_low", "module_CI_high", "log2FC")
  frac_cols <- c("abs_kME_rank_fraction_in_module",
                 "sus_res_fraction_matching_majority_direction")

  add_table_sheet <- function(sheet, title, note, data) {
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE, tabColour = "#5B6770")
    ncol_data <- max(1L, ncol(data))
    openxlsx::mergeCells(wb, sheet, cols = seq_len(ncol_data), rows = 1)
    openxlsx::writeData(wb, sheet, title, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = seq_len(ncol_data), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 1, heights = 23)
    openxlsx::mergeCells(wb, sheet, cols = seq_len(ncol_data), rows = 2)
    openxlsx::writeData(
      wb, sheet,
      if (nrow(data)) note else paste0(note, "  [NONE IN THIS SCOPE - 0 proteins.]"),
      startRow = 2, startCol = 1
    )
    openxlsx::addStyle(wb, sheet, note_style, rows = 2, cols = seq_len(ncol_data), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 2, heights = 42)
    # openxlsx::writeDataTable on a zero-row frame emits a table whose range
    # spans one data row, so an empty sheet reads back as a single all-NA
    # record. Write just the header instead, plus an explicit "none" marker, so
    # an empty class cannot be mistaken for one candidate.
    if (nrow(data)) {
      openxlsx::writeDataTable(
        wb, sheet, data, startRow = 4, startCol = 1,
        tableName = gsub("[^A-Za-z0-9]", "", sheet),
        tableStyle = "TableStyleLight9", withFilter = TRUE
      )
    } else {
      # Header only: no data row at all, so the sheet reads back as 0 records.
      # The "none" marker lives in the row-2 note, above the header.
      openxlsx::writeData(wb, sheet, data[0, , drop = FALSE], startRow = 4,
                          startCol = 1, colNames = TRUE)
    }
    openxlsx::addStyle(wb, sheet, header_style, rows = 4, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 4, heights = 46)
    if (nrow(data)) {
      rows <- 5:(nrow(data) + 4L)
      openxlsx::addStyle(wb, sheet, body_base, rows = rows, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
      apply_fmt <- function(cols, style) {
        hit <- match(intersect(cols, names(data)), names(data))
        if (length(hit)) openxlsx::addStyle(wb, sheet, style, rows = rows, cols = hit, gridExpand = TRUE, stack = TRUE)
      }
      apply_fmt(sci_cols, fdr_style)
      apply_fmt(num_cols, num_style)
      apply_fmt(frac_cols, fraction_style)
      if ("candidate_reason" %in% names(data)) {
        openxlsx::addStyle(wb, sheet, body_wrap, rows = rows,
                           cols = match("candidate_reason", names(data)),
                           gridExpand = TRUE, stack = TRUE)
      }
      # Shade the phenotype-linked rows. The label column is `Tier` on the
      # renamed review sheets and `candidate_tier` on the full audit sheets.
      tier_col <- intersect(c("candidate_tier", "Tier"), names(data))
      if (length(tier_col)) {
        tier <- as.character(data[[tier_col[[1]]]])
        a_rows <- rows[tier %in% c("A1", "A2")]
        d_rows <- rows[tier %in% "D"]
        if (length(a_rows)) openxlsx::addStyle(wb, sheet, tier_a_fill, rows = a_rows, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
        if (length(d_rows)) openxlsx::addStyle(wb, sheet, tier_d_fill, rows = d_rows, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
      }
    }
    widths <- pmin(46, pmax(11, nchar(names(data)) + 2))
    if ("candidate_reason" %in% names(data)) widths[match("candidate_reason", names(data))] <- 58
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = widths)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 5, firstActiveCol = 2)
  }

  # ---- README sheet
  openxlsx::addWorksheet(wb, "README", gridLines = FALSE, tabColour = "#23384D")
  openxlsx::writeData(wb, "README", paste0("WGCNA candidate-protein shortlist (", scope, ")"), startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, "README", title_style, rows = 1, cols = 1:4, gridExpand = TRUE, stack = TRUE)
  intro <- c(
    "Additive descriptive layer. No WGCNA, differential-abundance or module-level model is recomputed; every statistic is copied unchanged from a frozen upstream artifact and carries a provenance column.",
    "There is deliberately NO weighted composite importance score. Every criterion is a separate boolean so it can be inspected and rejected independently.",
    "Candidate tiers are NOT mutually exclusive: candidate_tier_all lists every tier a protein satisfies; candidate_tier is a convenience label only.",
    "Row order is a deterministic SORT (tier, SUS-RES FDR support, |kME|, |SUS-RES log2FC|, spatial consistency, ProteinGroupID), never a score."
  )
  openxlsx::writeData(wb, "README", data.frame(Notes = intro), startRow = 3, startCol = 1)
  openxlsx::addStyle(wb, "README", body_wrap, rows = 3:(3 + length(intro)), cols = 1, gridExpand = TRUE, stack = TRUE)

  tier_defs <- data.frame(
    Tier = c("A1", "A2", "D", "B", "C"),
    Name = unname(wcp_tier_display_names()[c("A1", "A2", "D", "B", "C")]),
    Definition = c(
      "clean_mapping AND SUS-RES BH FDR <= 0.05 in >=1 spatial context AND |kME| >= 0.60 AND in the module's top 10 by |kME|.",
      "clean_mapping AND SUS-RES BH FDR <= 0.05 in >=1 spatial context AND |kME| >= 0.60 AND NOT in the module's top 10.",
      "SUS-RES BH FDR <= 0.05 in >=1 spatial context AND |kME| < 0.60.",
      "|kME| >= 0.60 AND top-decile MEDIAN |SUS-RES log2FC| within dataset x contrast (sus_res_large_effect_typical). Protein-level FDR support NOT required. Spatial directional agreement is NOT part of this definition.",
      "Top 5 within module by |kME|, regardless of differential abundance. Answers: which proteins best represent this module?"
    ),
    Note = c(
      "Extreme network centrality PLUS protein-level support. Module-direction concordance is an additional flag, not a requirement: no module-level SUS-RES effect is FDR-supported in any dataset.",
      "Bona fide core module member by the frozen |kME| >= 0.60 definition, outside the top 10. No top-25 requirement is imposed; is_top_hub_25 is reported separately.",
      "Protein-level inferential evidence with weak module centrality. Unlike A1/A2 this class does not additionally require clean_mapping - a pre-existing asymmetry left unchanged; on current frozen inputs every FDR-supported protein is cleanly mapped, so it has no effect.",
      "Descriptive prioritisation only; no inferential support is claimed.",
      "Network-structural claim only; carries no phenotype evidence."
    ),
    stringsAsFactors = FALSE
  )
  tier_start <- 3 + length(intro) + 2L
  openxlsx::writeData(wb, "README", "Candidate tiers", startRow = tier_start, startCol = 1)
  openxlsx::addStyle(wb, "README", openxlsx::createStyle(textDecoration = "bold", fontSize = 11, fontColour = "#23384D"), rows = tier_start, cols = 1)
  openxlsx::writeDataTable(wb, "README", tier_defs, startRow = tier_start + 1L, startCol = 1, tableName = "TierDefinitions", tableStyle = "TableStyleLight9")

  notes <- wcp_evidence_class_notes()
  class_start <- tier_start + nrow(tier_defs) + 4L
  openxlsx::writeData(wb, "README", "Evidence classes: what each family of columns can and cannot support", startRow = class_start, startCol = 1)
  openxlsx::addStyle(wb, "README", openxlsx::createStyle(textDecoration = "bold", fontSize = 11, fontColour = "#23384D"), rows = class_start, cols = 1)
  openxlsx::writeDataTable(wb, "README", notes, startRow = class_start + 1L, startCol = 1, tableName = "EvidenceClasses", tableStyle = "TableStyleLight9")

  dict <- wcp_flag_dictionary()
  dict_start <- class_start + nrow(notes) + 4L
  openxlsx::writeData(wb, "README", "Field dictionary", startRow = dict_start, startCol = 1)
  openxlsx::addStyle(wb, "README", openxlsx::createStyle(textDecoration = "bold", fontSize = 11, fontColour = "#23384D"), rows = dict_start, cols = 1)
  openxlsx::writeDataTable(wb, "README", dict, startRow = dict_start + 1L, startCol = 1, tableName = "FieldDictionary", tableStyle = "TableStyleLight9")

  caveats <- data.frame(Caveat = c(
    "Tier A1 / A2 / D all require the SAME protein-level evidence (SUS-RES BH FDR <= 0.05 in >=1 spatial context). They differ ONLY in WGCNA topology - top-10 hub, module member, peripheral member. A1 is NOT more statistically significant than A2 or D; a Tier D protein routinely has a smaller FDR and a larger effect.",
    "The precedence A1 > A2 > D > B > C used by candidate_tier and by the review sort is a reading order, not a numerical or statistical score.",
    "GeneSignificanceP/FDR is an OMNIBUS CON/RES/SUS ANOVA: unsigned, and not a SUS-vs-RES phenotype correlation.",
    "GSEA leading-edge membership is NOT independent of the DA columns; that GSEA ranks on the moderated t statistic of the same limma fit.",
    "Hippocampal spatial contexts are repeated measurements from the SAME animals: agreement across them is spatial consistency, not independent replication.",
    "Those contexts are positively correlated, so their signs are NOT independent Bernoulli trials. Unanimous direction arises far more often than an independence calculation implies. No sign probability is computed or used anywhere in this layer, and spatial agreement alone never creates a candidate tier.",
    "Tier B is defined by module membership plus a large TYPICAL SUS-RES effect (top-decile median |log2FC| for the dataset). Spatial directional agreement is descriptive prioritisation information only.",
    "direction_matches_module compares a protein with a summary of the same proteins in the same animals. It is descriptive concordance, not independent support.",
    "No module-level SUS-RES effect is FDR-supported in any dataset; module estimates are context, not confirmation.",
    "This table selects CANDIDATES for manual inspection. It does not establish that any protein is a driver, is causal, is validated, or has been replicated."
  ), stringsAsFactors = FALSE)
  cav_start <- dict_start + nrow(dict) + 4L
  openxlsx::writeData(wb, "README", "Interpretation caveats", startRow = cav_start, startCol = 1)
  openxlsx::addStyle(wb, "README", openxlsx::createStyle(textDecoration = "bold", fontSize = 11, fontColour = "#23384D"), rows = cav_start, cols = 1)
  openxlsx::writeDataTable(wb, "README", caveats, startRow = cav_start + 1L, startCol = 1, tableName = "Caveats", tableStyle = "TableStyleLight9")
  openxlsx::setColWidths(wb, "README", cols = 1:4, widths = c(30, 42, 92, 66))
  openxlsx::addStyle(wb, "README", body_wrap, rows = 3:(cav_start + nrow(caveats) + 2L), cols = 1:4, gridExpand = TRUE, stack = TRUE)
  openxlsx::freezePane(wb, "README", firstActiveRow = 2)

  shortlist <- summary_tbl[summary_tbl$is_candidate %in% TRUE, , drop = FALSE]

  # ---- compact human-facing review sheets, placed before the audit tables
  add_table_sheet(
    "Protein_review", "Protein review (compact)",
    paste(
      "Compact view of every candidate for manual biological review.",
      "Reading order: Tier A1, then Tier A2, then Tier D, then Tier B, then Tier C only;",
      "within each block by SUS-RES FDR, large-effect flag, |kME|, median |log2FC|, ProteinGroupID.",
      "This is a deterministic sort for reading, NOT a statistical ranking.",
      "Source paths and internal contract fields are in All_candidates and the run manifest."
    ),
    wcp_protein_review_table(shortlist)
  )
  add_table_sheet(
    "SUS_RES_FDR_hits", "Proteins with SUS-RES BH FDR <= 0.05",
    paste(
      "Every protein reaching SUS-RES BH FDR <= 0.05 in at least one spatial context,",
      "regardless of tier or module centrality. This is the only per-protein inferential",
      "evidence in the workbook. Same columns and sort as Protein_review."
    ),
    wcp_fdr_hits_table(summary_tbl)
  )
  add_table_sheet(
    "Module_review", "Module summary (descriptive)",
    paste(
      "One row per module. DESCRIPTIVE ONLY: candidate fraction is not a measure of",
      "biological strength, because Tier C contributes exactly 5 proteins to every module",
      "by construction, so small modules necessarily show a higher fraction."
    ),
    wcp_module_review_table(summary_tbl)
  )

  add_table_sheet("All_candidates", "All candidates (any tier flag)",
                  "Every WGCNA protein satisfying at least one candidate flag. Sorted by tier, then SUS-RES FDR support, |kME|, |SUS-RES log2FC|, spatial consistency, ProteinGroupID. The order is a sort, not a score.",
                  shortlist)
  add_table_sheet("Tier_A1", "Tier A1 - FDR-supported phenotype-linked module hubs",
                  "Clean mapping, SUS-RES BH FDR <= 0.05 in at least one spatial context, |kME| >= 0.60, AND in the module's top 10 by |kME|. A1, A2 and D share the same protein-level evidence requirement and differ only in WGCNA topology; A1 is not 'more significant' than A2 or D.",
                  summary_tbl[summary_tbl$is_tier_A1 %in% TRUE, , drop = FALSE])
  add_table_sheet("Tier_A2", "Tier A2 - FDR-supported phenotype-linked module members",
                  "Same protein-level evidence and |kME| >= 0.60 as Tier A1, but NOT in the module's top 10 hubs. These are bona fide core module members with inferential SUS-RES support. No top-25 requirement is imposed; the frozen is_top_hub_25 flag is reported separately.",
                  summary_tbl[summary_tbl$is_tier_A2 %in% TRUE, , drop = FALSE])
  add_table_sheet("Tier_B", "Tier B - network + phenotype candidates",
                  "|kME| >= 0.60 AND a large typical SUS-RES effect (median |log2FC| across spatial contexts in the top decile for this dataset x contrast). No protein-level FDR support is claimed. Spatial directional agreement is reported alongside but is not part of the definition.",
                  summary_tbl[summary_tbl$is_tier_B %in% TRUE, , drop = FALSE])
  add_table_sheet("Tier_C", "Tier C - canonical module representatives",
                  "Top 5 members of each module by |kME|, regardless of differential abundance. A network-structural statement only; it carries no phenotype evidence.",
                  summary_tbl[summary_tbl$is_tier_C %in% TRUE, , drop = FALSE])
  add_table_sheet("Tier_D", "Tier D - FDR-supported phenotype-linked peripheral module members",
                  "SUS-RES BH FDR <= 0.05 with |kME| < 0.60. Present so that stress-responsive proteins are not discarded merely for being peripheral to a WGCNA module. A Tier D protein routinely has a SMALLER FDR and a LARGER effect than a Tier A1 protein: the A1/A2/D split describes network position, not evidence strength.",
                  summary_tbl[summary_tbl$is_tier_D %in% TRUE, , drop = FALSE])
  add_table_sheet("Top10_per_module", "Top 10 module representatives by |kME|",
                  "Exactly the top 10 members of each module by |kME|, ties broken by ProteinGroupID ascending. Selected by module centrality only, never by differential-abundance significance.",
                  top10)

  # Per-dataset sheets only help in the cross-dataset workbook. In a
  # single-dataset workbook DS_<that dataset> would just duplicate
  # All_candidates row for row, so it is omitted.
  if (length(by_dataset) > 1L) {
    for (dataset in names(by_dataset)) {
      sheet <- substr(paste0("DS_", dataset), 1, 31)
      d <- by_dataset[[dataset]]
      add_table_sheet(sheet, paste0("All candidates - ", dataset),
                      paste0("Candidate proteins for ", dataset,
                             ". Same columns and same sort order as All_candidates."),
                      d[d$is_candidate %in% TRUE, , drop = FALSE])
    }
  }

  # saveWorkbook alone emits a drawing relationship and a content-type override
  # for every sheet without writing the drawing part, which leaves a package
  # that strict OOXML readers reject. See R/utilities/xlsx_package_utils.R.
  xlsx_save_valid_workbook(wb, path)
  path
}

# --------------------------------------------------------------------- main

leading_edge <- load_leading_edge()
input_status <- list(leading_edge$status)

built <- list()
for (dataset in datasets) {
  message("[", SCRIPT_ID, "] building ", dataset)
  built[[dataset]] <- build_dataset(dataset, leading_edge$data)
}

emit <- function(scope, summary_tbl, long_tbl, by_dataset) {
  paths <- integration_paths(SUBSTEP, scope)
  summary_tbl <- summary_column_order(summary_tbl)
  long_tbl <- long_column_order(long_tbl)

  shortlist <- summary_tbl[summary_tbl$is_candidate %in% TRUE, , drop = FALSE]
  top10 <- summary_tbl[
    !is.na(summary_tbl$abs_kME_rank_in_module) &
      summary_tbl$abs_kME_rank_in_module <= wcp_module_hub_n(), , drop = FALSE]
  top10 <- top10[order(top10$dataset, top10$ModuleID,
                       top10$abs_kME_rank_in_module,
                       top10$ProteinGroupID, method = "radix"), , drop = FALSE]

  written <- list(
    all = write_integration_table(summary_tbl, paths, "wgcna_candidate_proteins_all.csv"),
    shortlist = write_integration_table(shortlist, paths, "wgcna_candidate_proteins_shortlist.csv"),
    top10 = write_integration_table(top10, paths, "wgcna_top10_per_module.csv")
  )
  # The long evidence table is a machine-readable archive, not figure-backing
  # data, so it is written once to tables/ and not mirrored into source_data.
  # In the cross-dataset scope it is skipped entirely: it would be a byte-for-byte
  # concatenation of the per-dataset files and roughly doubles the tree on disk.
  written$long <- if (identical(scope, "global")) {
    list(table = NA_character_, source = NULL)
  } else {
    list(
      table = write_csv_safe(
        long_tbl, file.path(paths$tables, "wgcna_candidate_proteins_long.csv")
      ),
      source = NULL
    )
  }
  xlsx <- write_workbook(
    file.path(paths$tables, "wgcna_candidate_proteins_shortlist.xlsx"),
    scope, summary_tbl, top10, by_dataset
  )
  readme <- write_readme(paths, scope, summary_tbl)

  # Resolved input paths are recorded here rather than repeated on every long
  # row; the long table keeps the compact da_comparison + input-hash pair.
  scope_datasets <- unique(as.character(summary_tbl$dataset))
  inputs <- c(
    unique(stats::na.omit(as.character(summary_tbl$membership_source_file))),
    unique(stats::na.omit(as.character(summary_tbl$module_effect_source_file))),
    vapply(scope_datasets,
           function(d) relpath(canonical_clusterprofiler_manifest_path(d)),
           character(1), USE.NAMES = FALSE),
    relpath(leading_edge_path())
  )
  inputs <- unique(inputs[nzchar(inputs)])
  write_integration_manifest(
    paths,
    inputs = as.list(stats::setNames(inputs, make.unique(basename(inputs)))),
    outputs = list(
      tables = unlist(lapply(written, `[[`, "table"), use.names = FALSE),
      source_data = unlist(lapply(written, `[[`, "source"), use.names = FALSE),
      workbook = xlsx, report = readme
    ),
    parameters = list(
      contract_version = wcp_contract_version(),
      scope = scope,
      datasets = datasets,
      contrasts = wcp_contrast_levels(),
      primary_contrast = wcp_primary_contrast(),
      high_kME_threshold = wcp_high_kme_threshold(),
      very_high_kME_threshold = wcp_very_high_kme_threshold(),
      fdr_threshold = wcp_fdr_threshold(),
      large_effect_quantile = wcp_large_effect_quantile(),
      large_effect_family = "dataset x spatial_unit x contrast",
      min_spatial_contexts = wcp_min_spatial_contexts(),
      module_representative_n = wcp_module_representative_n(),
      module_hub_n = wcp_module_hub_n()
    ),
    notes = paste(
      "Additive descriptive candidate-protein layer. No WGCNA, differential",
      "abundance or module-level model is recomputed; all statistics are copied",
      "unchanged from frozen upstream artifacts. No weighted composite score is",
      "produced. GSEA leading-edge membership is derived from the same ranked DA",
      "statistics and is not independent evidence. Spatial contexts are repeated",
      "measurements from the same animals and are not independent replication."
    )
  )
  list(paths = paths, summary = summary_tbl, shortlist = shortlist, top10 = top10)
}

per_dataset_summary <- lapply(built, `[[`, "summary")
emitted <- list()
for (dataset in datasets) {
  emitted[[dataset]] <- emit(
    dataset, built[[dataset]]$summary, built[[dataset]]$long,
    stats::setNames(list(built[[dataset]]$summary), dataset)
  )
}

# The `global` scope is the cross-dataset table. It is only rewritten when the
# run actually covers every dataset; a single --dataset run must not silently
# replace the combined table with one dataset's rows.
write_global <- setequal(datasets, valid_datasets())
if (write_global) {
  combined_summary <- wcp_order_candidates(bind_rows(per_dataset_summary))
  combined_long <- bind_rows(lapply(built, `[[`, "long"))
  global <- emit("global", combined_summary, combined_long, per_dataset_summary)
} else {
  global <- emitted[[datasets[[1]]]]
  message("[", SCRIPT_ID, "] partial dataset selection (",
          paste(datasets, collapse = ", "),
          "); leaving the cross-dataset 'global' outputs untouched.")
}

# ------------------------------------------------------------ console report

fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)
cat("\n===== WGCNA candidate-protein shortlist =====\n")
cat("Contract    : ", wcp_contract_version(), "\n", sep = "")
cat("Datasets    : ", paste(datasets, collapse = ", "), "\n", sep = "")
cat("Primary     : ", wcp_primary_contrast(),
    "  (positive log2FC = higher in SUS)\n", sep = "")

for (dataset in datasets) {
  s <- emitted[[dataset]]$summary
  cat("\n--- ", dataset, " ---\n", sep = "")
  cat("  WGCNA proteins examined : ", fmt_int(nrow(s)),
      " in ", fmt_int(dplyr::n_distinct(s$ModuleID)), " modules\n", sep = "")
  cat("  Cleanly mapped          : ", fmt_int(sum(s$clean_mapping %in% TRUE)),
      " (", fmt_int(sum(!(s$clean_mapping %in% TRUE))), " gene-ambiguous, excluded from Tier A1/A2)\n", sep = "")
  cat("  SUS-RES FDR<=0.05       : ", fmt_int(sum(s$sus_res_fdr05_any_context %in% TRUE)),
      " protein(s) in >=1 spatial context",
      "  [A1 ", fmt_int(sum(s$is_tier_A1 %in% TRUE)),
      " + A2 ", fmt_int(sum(s$is_tier_A2 %in% TRUE)),
      " + D ", fmt_int(sum(s$is_tier_D %in% TRUE)),
      "; unclassified ", fmt_int(nrow(wcp_unclassified_fdr_support(s))), "]\n", sep = "")
  cat("  Tier A1/A2/B/C/D        : ",
      fmt_int(sum(s$is_tier_A1 %in% TRUE)), " / ",
      fmt_int(sum(s$is_tier_A2 %in% TRUE)), " / ",
      fmt_int(sum(s$is_tier_B %in% TRUE)), " / ",
      fmt_int(sum(s$is_tier_C %in% TRUE)), " / ",
      fmt_int(sum(s$is_tier_D %in% TRUE)), "\n", sep = "")
  cat("  Any candidate flag      : ", fmt_int(sum(s$is_candidate %in% TRUE)), "\n", sep = "")

  focus <- s
  if (nzchar(module_focus)) {
    focus <- s[as.character(s$ModuleID) == module_focus, , drop = FALSE]
    cat("  [--module ", module_focus, "] ", fmt_int(nrow(focus)),
        " member(s); ", fmt_int(sum(focus$is_candidate %in% TRUE)),
        " flagged\n", sep = "")
    if (!nrow(focus)) next
  }

  label <- function(df) {
    sym <- as.character(df$GeneSymbol)
    sym[is.na(sym) | !nzchar(sym)] <- as.character(df$ProteinGroupID)[is.na(sym) | !nzchar(sym)]
    sym
  }
  top_per_module <- focus[focus$is_candidate %in% TRUE, , drop = FALSE]
  if (nrow(top_per_module)) {
    cat("  Top candidate per module (tier, |kME|):\n")
    split_mod <- split(seq_len(nrow(top_per_module)), top_per_module$ModuleID)
    for (mid in sort(names(split_mod))) {
      row <- top_per_module[split_mod[[mid]][1], , drop = FALSE]
      cat(sprintf("    %-12s %-14s tier %-4s |kME|=%.2f  rank %s/%s\n",
                  mid, label(row), row$candidate_tier_all,
                  row$abs_kME, row$abs_kME_rank_in_module, row$n_module_members))
    }
  }
}

top_overall <- global$shortlist
if (nzchar(module_focus)) {
  top_overall <- top_overall[as.character(top_overall$ModuleID) == module_focus, , drop = FALSE]
}
top_overall <- utils::head(top_overall, 10L)
if (nrow(top_overall)) {
  cat("\n--- First 10 rows in display order (deterministic sort, NOT a statistical ranking) ---\n")
  sym <- as.character(top_overall$GeneSymbol)
  sym[is.na(sym) | !nzchar(sym)] <- as.character(top_overall$ProteinGroupID)[is.na(sym) | !nzchar(sym)]
  for (i in seq_len(nrow(top_overall))) {
    cat(sprintf("  %2d. %-14s %-16s %-10s tier %-4s |kME|=%.2f  SUS-RES FDR=%s\n",
                i, sym[[i]], top_overall$dataset[[i]], top_overall$ModuleID[[i]],
                top_overall$candidate_tier_all[[i]], top_overall$abs_kME[[i]],
                ifelse(is.na(top_overall$sus_res_min_BH_FDR[[i]]), "NA",
                       formatC(top_overall$sus_res_min_BH_FDR[[i]], format = "g", digits = 3))))
  }
}

cat("\nOutputs written under:\n")
for (dataset in datasets) cat("  ", relpath(emitted[[dataset]]$paths$tables), "\n", sep = "")
if (write_global) {
  cat("  ", relpath(global$paths$tables), "\n", sep = "")
} else {
  cat("  (cross-dataset 'global' outputs not rewritten: partial dataset selection)\n")
}
cat("\nInspect one module, e.g.:\n")
cat("  Rscript ", SCRIPT_ID, " --dataset neuron_neuropil --module WGCNA_m02\n", sep = "")
cat("  R: subset(readr::read_csv(\"wgcna_candidate_proteins_all.csv\"), ModuleID == \"WGCNA_m02\")\n")
cat("\nCaveats: GSEA leading-edge membership is NOT independent of the DA columns;\n")
cat("spatial contexts are repeated measurements from the same animals, not replication;\n")
cat("GeneSignificance is an omnibus CON/RES/SUS ANOVA, not a SUS-vs-RES correlation.\n")
