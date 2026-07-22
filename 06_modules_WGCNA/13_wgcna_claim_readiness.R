#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/13_wgcna_claim_readiness.R
# Stage: optional manuscript-readiness handoff
# Scope: microglia
# Consumes: Stage 05/06/07 plus optional Stage 09/11 and additive Stage 12.
# Produces: one non-circular 22-row claim-readiness table.
# Does not alter Stages 01-12, WGCNA state, labels, or model results.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_labeling_utils.R"))
source(repo_path("R", "wgcna_reviewed_label_registry.R"))
source(repo_path("R", "schema_validation.R"))

required_pkgs <- c("dplyr", "readr", "tibble", "yaml")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

run <- wgcna_cli()
if (!identical(run$dataset, "microglia")) stop("13_wgcna_claim_readiness.R is restricted to microglia.", call. = FALSE)
DATASET <- "microglia"
OUT <- create_module_dirs("06_modules_WGCNA", file.path("claim_readiness", DATASET))
FILES <- resolve_wgcna_files(DATASET)
TABLE_DIR <- path_results("tables", "06_modules_WGCNA")
AUDIT_DIR <- path_results("reviewer_audit", "microglia_wgcna_nature_readiness")
inputs <- list(
  stage05_module = file.path(TABLE_DIR, "group_effects", DATASET, "module_group_effects.csv"),
  stage05_supermodule = file.path(TABLE_DIR, "group_effects", DATASET, "supermodule_group_effects.csv"),
  stage06_module = file.path(TABLE_DIR, "module_annotation", DATASET, "WGCNA_module_biological_annotation.csv"),
  stage06_supermodule = file.path(TABLE_DIR, "module_annotation", DATASET, "WGCNA_supermodule_biological_annotation.csv"),
  stage07_lookup = file.path(TABLE_DIR, "interpretable_summary", DATASET, "WGCNA_final_label_lookup.csv"),
  stage09_neuropil = file.path(TABLE_DIR, "microglia_neuropil_independence", DATASET, "microglia_module_neuropil_independence_classification.csv"),
  stage11_robustness = file.path(TABLE_DIR, "module_robustness_sensitivity", DATASET, "WGCNA_claim_gate_audit.csv"),
  stage12_modules = file.path(AUDIT_DIR, "module_robustness_consensus.csv"),
  stage12_blocks = file.path(AUDIT_DIR, "higher_order_block_readiness_summary.csv"),
  stage01_supermodule_annotation = FILES$supermodule_annotation
)
if (run$dry_run) {
  required_names <- c("stage05_module", "stage05_supermodule", "stage06_module", "stage06_supermodule", "stage07_lookup", "stage12_modules", "stage12_blocks", "stage01_supermodule_annotation")
  for (nm in required_names) dry_run_line(nm, inputs[[nm]], if (file.exists(inputs[[nm]])) "PASS" else "FAIL")
  for (nm in c("stage09_neuropil", "stage11_robustness")) dry_run_line(paste0(nm, " (optional)"), inputs[[nm]], if (file.exists(inputs[[nm]])) "PASS" else "INFO")
  dry_run_line("Table output", file.path(OUT$tables, "WGCNA_entity_claim_readiness.csv"))
  dry_run_line("Source-data output", file.path(OUT$source_data, "WGCNA_entity_claim_readiness_source.csv"))
  quit(status = 0)
}
required <- unlist(inputs[c("stage05_module", "stage05_supermodule", "stage06_module", "stage06_supermodule", "stage07_lookup", "stage12_modules", "stage12_blocks", "stage01_supermodule_annotation")], use.names = FALSE)
if (any(!file.exists(required))) stop("Required claim-readiness inputs missing: ", paste(required[!file.exists(required)], collapse = ", "), call. = FALSE)
read_csv <- function(path) readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
repo_relative_path <- function(path) {
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (startsWith(resolved, paste0(root, "/"))) substring(resolved, nchar(root) + 2L) else resolved
}
lookup <- read_csv(inputs$stage07_lookup)
member_map <- wgcna_normalize_current_member_map(read_csv(inputs$stage01_supermodule_annotation), DATASET)
wgcna_validate_canonical_lookup(lookup, DATASET, member_map)
if (nrow(lookup) != 22L) stop("Stage 07 lookup must contain exactly 22 current microglia identities.", call. = FALSE)

modules <- read_csv(inputs$stage12_modules)
blocks <- read_csv(inputs$stage12_blocks)
if (nrow(modules) != 13L || anyDuplicated(modules$ModuleID)) stop("Stage 12 module readiness must contain exactly 13 unique ModuleIDs.", call. = FALSE)
if (!setequal(blocks$SupermoduleID, c("SM01", "SM03", "SM09"))) stop("Stage 12 block readiness must contain SM01, SM03 and SM09 only.", call. = FALSE)
module_ready <- modules |>
  dplyr::transmute(level = "module", entity_id = as.character(.data$ModuleID), primary_architecture_status, spatial_dependence_class,
                   animal_stability_status, strict_nonspatial_sensitivity_status,
                   conventional_preservation_status = .data$preservation_interpretation,
                   conventional_preservation_claim_gate_eligible = .data$conventional_preservation_claim_gate_eligible,
                   allowed_claim_scope = .data$claim_scope)
super_members <- member_map |>
  dplyr::group_by(.data$SupermoduleID) |>
  dplyr::summarise(n_member_modules = dplyr::n(), member_ModuleID = dplyr::first(.data$ModuleID), .groups = "drop")
singleton_members <- super_members |> dplyr::filter(.data$n_member_modules == 1L)
if (nrow(singleton_members) != 6L || anyDuplicated(singleton_members$member_ModuleID) ||
    any(!singleton_members$member_ModuleID %in% module_ready$entity_id)) {
  stop("Every singleton compatibility identity must map to exactly one distinct current module.", call. = FALSE)
}
super_ready <- super_members |>
  dplyr::left_join(blocks, by = c("SupermoduleID"), relationship = "many-to-one") |>
  dplyr::left_join(module_ready |> dplyr::select(member_ModuleID = "entity_id", member_primary = "primary_architecture_status", member_spatial = "spatial_dependence_class", member_animal = "animal_stability_status", member_strict = "strict_nonspatial_sensitivity_status", member_preservation = "conventional_preservation_status", member_preservation_gate = "conventional_preservation_claim_gate_eligible", member_scope = "allowed_claim_scope"), by = "member_ModuleID", relationship = "many-to-one") |>
  dplyr::transmute(level = "supermodule", entity_id = as.character(.data$SupermoduleID),
    primary_architecture_status = dplyr::if_else(.data$n_member_modules == 1L, .data$member_primary, as.character(.data$higher_order_block_classification)),
    spatial_dependence_class = dplyr::if_else(.data$n_member_modules == 1L, .data$member_spatial, dplyr::if_else(.data$spatial_adjusted_support == "fail", "region_organized_block_or_limited_support", "spatially_adjusted_block_support")),
    animal_stability_status = dplyr::if_else(.data$n_member_modules == 1L, .data$member_animal, as.character(.data$leave_one_animal_out_support)),
    strict_nonspatial_sensitivity_status = dplyr::if_else(.data$n_member_modules == 1L, .data$member_strict, as.character(.data$strict_nonspatial_support)),
    conventional_preservation_status = dplyr::if_else(.data$n_member_modules == 1L, .data$member_preservation, "not_applicable_higher_order_block"),
    conventional_preservation_claim_gate_eligible = dplyr::if_else(.data$n_member_modules == 1L, .data$member_preservation_gate, FALSE),
    allowed_claim_scope = dplyr::if_else(.data$n_member_modules == 1L, .data$member_scope, as.character(.data$claim_scope)))
readiness <- dplyr::bind_rows(module_ready, super_ready)
claim_identity <- dplyr::bind_rows(
  module_ready |> dplyr::transmute(
    level = "module", entity_id,
    canonical_claim_entity_id = .data$entity_id,
    claim_entity_role = "canonical_module",
    compatibility_target_level = NA_character_,
    compatibility_target_entity_id = NA_character_
  ),
  super_members |> dplyr::transmute(
    level = "supermodule", entity_id = as.character(.data$SupermoduleID),
    canonical_claim_entity_id = dplyr::if_else(.data$n_member_modules == 1L, .data$member_ModuleID, as.character(.data$SupermoduleID)),
    claim_entity_role = dplyr::if_else(.data$n_member_modules == 1L, "compatibility_alias", "higher_order_block"),
    compatibility_target_level = dplyr::if_else(.data$n_member_modules == 1L, "module", NA_character_),
    compatibility_target_entity_id = dplyr::if_else(.data$n_member_modules == 1L, .data$member_ModuleID, NA_character_)
  )
)

context_from <- function(path, level, ids) {
  x <- read_csv(path); id_col <- intersect(c("ModuleID", "SupermoduleID", "module_id", "supermodule_id"), names(x))[[1]]
  context_col <- intersect(c("dominant_microenvironment_class", "microenvironment_class", "dominant_context_class"), names(x))
  confidence_col <- intersect(c("annotation_confidence", "label_confidence", "context_confidence"), names(x))
  tibble::tibble(level = level, entity_id = as.character(x[[id_col]]),
    context_assignment_class = if (length(context_col)) as.character(x[[context_col[[1]]]]) else "not_available",
    context_assignment_confidence = if (length(confidence_col)) as.character(x[[confidence_col[[1]]]]) else "not_available") |>
    dplyr::filter(.data$entity_id %in% ids) |>
    dplyr::distinct(.data$level, .data$entity_id, .keep_all = TRUE)
}
context <- dplyr::bind_rows(context_from(inputs$stage06_module, "module", module_ready$entity_id), context_from(inputs$stage06_supermodule, "supermodule", super_ready$entity_id))

effect_from <- function(path, level, id_col, expected_ids) {
  selected <- read_csv(path) |>
    dplyr::filter(
      .data$dataset == DATASET,
      .data$contrast == "SUS - RES",
      .data$effect_scope == "spatial_adjusted_global",
      .data$spatial_unit == "global_spatial_adjusted"
    )
  selected_ids <- as.character(selected[[id_col]])
  id_counts <- table(selected_ids)
  if (nrow(selected) != length(expected_ids) || !setequal(selected_ids, expected_ids) ||
      any(id_counts != 1L)) {
    stop(
      "Stage 05 selected endpoint must contain exactly one SUS - RES / spatial_adjusted_global / global_spatial_adjusted row per current ",
      level, " identity; selected rows=", nrow(selected), ".",
      call. = FALSE
    )
  }
  source_file <- repo_relative_path(path)
  selected |>
    dplyr::transmute(
      level = level,
      entity_id = as.character(.data[[id_col]]),
      selected_contrast = as.character(.data$contrast),
      selected_effect_scope = as.character(.data$effect_scope),
      selected_spatial_unit = as.character(.data$spatial_unit),
      group_effect_estimate = suppressWarnings(as.numeric(.data$estimate)),
      group_effect_SE = suppressWarnings(as.numeric(.data$SE)),
      group_effect_CI_low = .data$group_effect_estimate - 1.96 * .data$group_effect_SE,
      group_effect_CI_high = .data$group_effect_estimate + 1.96 * .data$group_effect_SE,
      group_effect_p_value = suppressWarnings(as.numeric(.data$p_value)),
      group_effect_FDR_within_dataset_level = suppressWarnings(as.numeric(.data$FDR_within_dataset_level)),
      group_effect_FDR_global = suppressWarnings(as.numeric(.data$FDR_global)),
      group_effect_FDR = .data$group_effect_FDR_global,
      group_effect_adjustment_scope = "BH across combined current module and supermodule Stage 05 effect rows",
      group_effect_evidence_status = as.character(.data$evidence_status),
      group_effect_model_formula = as.character(.data$formula_used),
      group_effect_model_type = as.character(.data$model_type),
      group_effect_n_animals = suppressWarnings(as.integer(.data$n_animals_total)),
      group_effect_n_samples = suppressWarnings(as.integer(.data$n_samples_total)),
      group_effect_biological_replicate_unit = as.character(.data$biological_replicate_unit),
      group_effect_source_file = source_file,
      group_effect_source_key = paste(.data$dataset, level, .data[[id_col]], .data$contrast, .data$effect_scope, .data$spatial_unit, sep = "|"),
      model_stability = dplyr::if_else(.data$primary_model_stable %in% TRUE & !(.data$singular_model %in% TRUE), "stable", "singular_or_unstable_diagnostic_only")
    )
}
effects <- dplyr::bind_rows(
  effect_from(inputs$stage05_module, "module", "module_id", module_ready$entity_id),
  effect_from(inputs$stage05_supermodule, "supermodule", "supermodule_id", super_ready$entity_id)
)
neuropil <- if (file.exists(inputs$stage09_neuropil)) {
  x <- read_csv(inputs$stage09_neuropil)
  id_candidates <- intersect(c("ModuleID", "module_id"), names(x)); status_candidates <- intersect(c("neuropil_independence_status", "claim_gate_status", "classification"), names(x))
  if (!length(id_candidates) || !length(status_candidates)) {
    tibble::tibble(level = character(), entity_id = character(), neuropil_independence_status = character())
  } else {
    tibble::tibble(level = "module", entity_id = as.character(x[[id_candidates[[1]]]]), neuropil_independence_status = as.character(x[[status_candidates[[1]]]])) |>
      dplyr::distinct(.data$level, .data$entity_id, .keep_all = TRUE)
  }
} else tibble::tibble(level = character(), entity_id = character(), neuropil_independence_status = character())

out <- lookup |>
  dplyr::transmute(dataset, level, entity_id, canonical_biological_label, canonical_short_label, structural_status,
    biological_process_confidence = .data$biological_label_confidence) |>
  dplyr::left_join(claim_identity, by = c("level", "entity_id"), relationship = "one-to-one") |>
  dplyr::left_join(readiness, by = c("level", "entity_id"), relationship = "one-to-one") |>
  dplyr::left_join(context, by = c("level", "entity_id"), relationship = "many-to-one") |>
  dplyr::left_join(neuropil, by = c("level", "entity_id"), relationship = "many-to-one") |>
  dplyr::left_join(effects, by = c("level", "entity_id"), relationship = "many-to-one") |>
  dplyr::mutate(
    neuropil_independence_status = dplyr::if_else(.data$level == "supermodule", "not_available_supermodule_compatibility_identity", dplyr::coalesce(.data$neuropil_independence_status, "not_available")),
    group_effect_status = dplyr::case_when(.data$model_stability != "stable" ~ "diagnostic_only_model", is.finite(.data$group_effect_FDR_global) & .data$group_effect_FDR_global <= .05 ~ "FDR_supported", is.finite(.data$group_effect_FDR_global) ~ "not_FDR_supported", TRUE ~ "not_available"),
    separate_manuscript_claim_allowed = dplyr::case_when(
      .data$claim_entity_role == "compatibility_alias" ~ FALSE,
      .data$claim_entity_role == "higher_order_block" ~ .data$primary_architecture_status %in% c("robust_spatially_independent_block", "reproducible_spatially_organized_block"),
      .data$claim_entity_role == "canonical_module" ~ .data$primary_architecture_status %in% c("robust_spatially_independent", "reproducible_spatially_organized") & !is.na(.data$context_assignment_class) & .data$context_assignment_class != "not_available",
      TRUE ~ FALSE
    ),
    prohibited_claim_scope = "No purified-microglia, cell-intrinsic, causal, independent-replication, or unsupported stress-effect claim.",
    manuscript_placement = dplyr::case_when(
      .data$claim_entity_role == "compatibility_alias" ~ "compatibility_only",
      .data$separate_manuscript_claim_allowed ~ "main-text biological co-variation candidate",
      TRUE ~ "Supplementary or diagnostic only"
    ),
    allowed_wording = dplyr::case_when(
      .data$claim_entity_role == "compatibility_alias" ~ paste0("Compatibility alias for ", .data$compatibility_target_entity_id, "; not a separate biological claim."),
      grepl("reproducible_spatially_organized", .data$primary_architecture_status) ~ "A region-organized microglia-associated ROI co-variation program.",
      grepl("robust", .data$primary_architecture_status) ~ "A spatially robust microglia-associated ROI co-variation program.",
      TRUE ~ "A descriptive fixed-membership microglia-associated ROI co-variation pattern."
    ),
    readiness_contract_version = "microglia_wgcna_claim_readiness_v2"
  ) |>
  dplyr::select(dataset, level, entity_id, canonical_claim_entity_id, claim_entity_role, separate_manuscript_claim_allowed, compatibility_target_level, compatibility_target_entity_id, canonical_biological_label, canonical_short_label, structural_status, primary_architecture_status, spatial_dependence_class, animal_stability_status, strict_nonspatial_sensitivity_status, conventional_preservation_status, conventional_preservation_claim_gate_eligible, biological_process_confidence, context_assignment_class, context_assignment_confidence, neuropil_independence_status, selected_contrast, selected_effect_scope, selected_spatial_unit, group_effect_estimate, group_effect_SE, group_effect_CI_low, group_effect_CI_high, group_effect_p_value, group_effect_FDR_within_dataset_level, group_effect_FDR_global, group_effect_FDR, group_effect_adjustment_scope, group_effect_evidence_status, group_effect_status, group_effect_model_formula, group_effect_model_type, group_effect_n_animals, group_effect_n_samples, group_effect_biological_replicate_unit, group_effect_source_file, group_effect_source_key, model_stability, allowed_claim_scope, prohibited_claim_scope, manuscript_placement, allowed_wording, readiness_contract_version) |>
  dplyr::arrange(factor(.data$level, levels = c("module", "supermodule")), .data$entity_id)
if (nrow(out) != 22L || anyDuplicated(out[c("dataset", "level", "entity_id")])) stop("Claim-readiness join did not yield exactly one row per 22 stable identities.", call. = FALSE)
if (any(!out$entity_id %in% c(sprintf("WGCNA_m%02d", 1:13), sprintf("SM%02d", 1:9)))) stop("Claim-readiness table contains stale IDs.", call. = FALSE)
if (any(out$conventional_preservation_claim_gate_eligible %in% TRUE)) stop("Conventional unblocked WGCNA preservation cannot be claim-gating.", call. = FALSE)
aliases <- out |> dplyr::filter(.data$claim_entity_role == "compatibility_alias")
if (nrow(aliases) != 6L || any(aliases$separate_manuscript_claim_allowed) ||
    any(aliases$compatibility_target_level != "module") ||
    any(aliases$canonical_claim_entity_id != aliases$compatibility_target_entity_id) ||
    any(!aliases$compatibility_target_entity_id %in% out$entity_id[out$level == "module"])) {
  stop("Compatibility-alias validation failed.", call. = FALSE)
}
claimable_blocks <- out$entity_id[out$claim_entity_role == "higher_order_block" & out$separate_manuscript_claim_allowed]
if (!identical(claimable_blocks, "SM09")) stop("SM09 must be the only separately claimable higher-order block.", call. = FALSE)
if (anyDuplicated(out$canonical_claim_entity_id[out$separate_manuscript_claim_allowed])) stop("Canonical claim entities are duplicated among independently claimable rows.", call. = FALSE)
if (!isTRUE(all.equal(out$group_effect_FDR, out$group_effect_FDR_global, check.attributes = FALSE))) stop("group_effect_FDR must remain an exact compatibility alias for FDR_global.", call. = FALSE)
validate_table_schema(out, "wgcna_entity_claim_readiness", strict = TRUE)
readr::write_csv(out, file.path(OUT$tables, "WGCNA_entity_claim_readiness.csv"), na = "")
readr::write_csv(out, file.path(OUT$source_data, "WGCNA_entity_claim_readiness_source.csv"), na = "")
write_run_manifest(file.path(OUT$logs, "run_manifest.yml"), inputs = inputs, outputs = list(claim_readiness = file.path(OUT$tables, "WGCNA_entity_claim_readiness.csv"), claim_readiness_source = file.path(OUT$source_data, "WGCNA_entity_claim_readiness_source.csv")), parameters = list(dataset = DATASET, readiness_contract_version = "microglia_wgcna_claim_readiness_v2", selected_endpoint = "SUS - RES | spatial_adjusted_global | global_spatial_adjusted", FDR_claim_status_source = "FDR_global"), notes = "Non-circular manuscript handoff. Singleton supermodule IDs are compatibility aliases, not separate claims. Stages 05-07 and 12 do not read this output.")
message("Microglia WGCNA claim-readiness handoff complete.")
