#!/usr/bin/env Rscript
#
# Protein-level spatial / cell-affinity atlas, plus the SUS-RES overlay.
#
# SCOPE is the UNION of every WGCNA candidate-shortlist protein and every
# protein with SUS-RES FDR <= 0.05 in at least one canonical spatial unit,
# deduplicated by dataset x ProteinGroupID. The full candidate atlas is NOT
# narrowed to the FDR-supported subset.
#
# ORDER OF OPERATIONS MATTERS
#   Baseline spatial identity is established from CON animals ONLY. The SUS-RES
#   overlay is joined afterwards from canonical DA outputs and never contributes
#   to identity. No differential statistic is recomputed here.
#
# USAGE
#   Rscript 11_spatial_systems/09_protein_spatial_cell_atlas.R
#   Rscript 11_spatial_systems/09_protein_spatial_cell_atlas.R --dry-run

source("R/paths.R")
source("R/dataset_config.R")
source("R/integration_utils.R")
source("R/qc_exploration_utils.R")
source("R/protigy_input_utils.R")
source("R/spatial_systems_data_utils.R")
source("R/spatial_atlas_utils.R")

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

SCRIPT_ID <- "11_spatial_systems/09_protein_spatial_cell_atlas.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function() {
  d <- path_results("tables", "11_spatial_systems", "atlas"); dir_create(d); d
}
CAND <- path_results("tables", "10_biological_integration",
                     "wgcna_candidate_protein_shortlist", "global",
                     "wgcna_candidate_proteins_all.csv")
EMP <- path_results("tables", "03_qc_exploration", "05_empirical_roi_marker_discovery",
                    "empirical_roi_marker_sets.csv")
REF <- repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")
DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Protein spatial/cell atlas and SUS-RES overlay.\n")
  dry_run_inputs(SCRIPT_ID, list(candidate_shortlist = CAND,
                                 empirical_marker_sets = EMP,
                                 reference_marker_panels = REF))
  cat("[DRY-RUN] DA statistics are JOINED from canonical outputs, never recomputed.\n")
  quit(save = "no", status = 0L)
}

read_req <- function(p, lbl) {
  if (!file.exists(p)) stop("missing_required_input: ", lbl, ": ", p, call. = FALSE)
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
}

# ==================================================== PART 13: atlas scope

cand <- read_req(CAND, "WGCNA candidate shortlist")
cand$is_sus_res_fdr_supported <- cand$sus_res_fdr05_any_context %in% TRUE
cand$is_wgcna_candidate <- !is.na(cand$candidate_tier) & nzchar(cand$candidate_tier)

scope <- cand[cand$is_sus_res_fdr_supported | cand$is_wgcna_candidate, , drop = FALSE]
scope <- scope[!duplicated(paste(scope$dataset, scope$ProteinGroupID)), , drop = FALSE]
n_fdr <- sum(scope$is_sus_res_fdr_supported)
message(sprintf("Atlas scope: %d dataset x ProteinGroupID rows (%d FDR-supported)",
                nrow(scope), n_fdr))

# ============================== PART 14: CON baseline protein spatial identity

message("Building CON baseline protein spatial profiles")
baseline_rows <- list(); metric_rows <- list()
for (ds in DATASETS) {
  want <- unique(scope$ProteinGroupID[scope$dataset == ds])
  if (!length(want)) next
  inputs <- resolve_dataset_inputs(ds, purpose = "wgcna", script = SCRIPT_ID,
                                   stage = "networks")
  md <- path_processed("01_preprocessing", "06_merged_metadata_module_score", ds,
                       "sample_metadata_merged_clean_for_module_scores.xlsx")
  canonical <- qc_load_canonical_expression(inputs$expression_file, md,
                                            dataset = ds, strict = TRUE)
  lv <- sps_build_spatial_levels(ds, canonical = canonical)

  # CON animals only, using the VALIDATED bilateral level-2 matrix
  con_cols <- lv$level2$meta$StressGroup == "CON"
  b <- lv$level2$mat[intersect(want, rownames(lv$level2$mat)), con_cols, drop = FALSE]
  meta <- lv$level2$meta[con_cols, , drop = FALSE]
  units <- sort(unique(meta$SpatialUnit))
  # mean across CON animals per spatial unit
  prof <- vapply(units, function(u)
    rowMeans(b[, meta$SpatialUnit == u, drop = FALSE], na.rm = TRUE),
    numeric(nrow(b)))
  if (is.null(dim(prof))) prof <- matrix(prof, nrow = nrow(b))
  dimnames(prof) <- list(rownames(b), tolower(units))

  baseline_rows[[ds]] <- cbind(
    data.frame(dataset = ds, ProteinGroupID = rownames(prof),
               n_con_animals = length(unique(meta$AnimalID)),
               stringsAsFactors = FALSE),
    as.data.frame(prof, check.names = FALSE))
  metric_rows[[ds]] <- dplyr::bind_rows(lapply(rownames(prof), function(r)
    cbind(data.frame(dataset = ds, ProteinGroupID = r, stringsAsFactors = FALSE),
          sat_profile_metrics(prof[r, ], colnames(prof)))))
  message(sprintf("  %-16s %d proteins x %d units (CON n=%d)", ds, nrow(prof),
                  ncol(prof), length(unique(meta$AnimalID))))
}
baseline <- dplyr::bind_rows(baseline_rows)
metrics <- dplyr::bind_rows(metric_rows)

# ============== PART 14b: contrast-level bilateral support (already validated)

bsi_path <- path_results("tables", "11_spatial_systems", "bilateral",
                         "bilateral_spatial_identity_protein_level.csv")
contrast_support <- NULL
if (file.exists(bsi_path)) {
  bsi <- read_req(bsi_path, "bilateral spatial identity")
  contrast_support <- bsi %>%
    dplyr::group_by(.data$dataset, .data$ProteinGroupID) %>%
    dplyr::summarise(
      contrast_level_bilateral_support_n_contrasts = dplyr::n(),
      contrast_level_bilateral_sign_agreement = mean(.data$sign_agreement, na.rm = TRUE),
      contrast_level_bilateral_median_abs_LR_difference =
        stats::median(.data$abs_L_minus_R, na.rm = TRUE),
      .groups = "drop")
}

# ==================================== PART 15 + 16: compartment and reference

message("Attaching empirical compartment and reference context")
emp <- read_req(EMP, "empirical marker sets")
eff_cols <- grep("^logFC_", names(emp), value = TRUE)
fdr_cols <- grep("^FDR_", names(emp), value = TRUE)
emp_ctx_rows <- list()
for (ds in DATASETS) {
  id_col <- intersect(c(paste0("ProteinGroupID_", ds), "ProteinGroupID"), names(emp))[1]
  z <- emp[, c(id_col, "marker_set", eff_cols, fdr_cols), drop = FALSE]
  names(z)[1] <- "ProteinGroupID"
  z <- z[!is.na(z$ProteinGroupID) & nzchar(z$ProteinGroupID), , drop = FALSE]
  # continuous effects, not marker yes/no only; a protein may be partially or
  # mixedly assigned and is NOT forced into one exclusive compartment
  agg <- z %>%
    dplyr::group_by(.data$ProteinGroupID) %>%
    dplyr::summarise(
      empirical_marker_sets = paste(sort(unique(.data$marker_set)), collapse = ";"),
      n_empirical_marker_sets = dplyr::n_distinct(.data$marker_set),
      dplyr::across(dplyr::all_of(c(eff_cols, fdr_cols)),
                    ~ suppressWarnings(mean(.x, na.rm = TRUE))),
      .groups = "drop")
  agg$dataset <- ds
  emp_ctx_rows[[ds]] <- agg
}
emp_ctx <- dplyr::bind_rows(emp_ctx_rows)

ref <- read_req(REF, "reference marker panels")
ref_lut <- ref %>%
  dplyr::mutate(sym = toupper(trimws(.data$gene_symbol))) %>%
  dplyr::group_by(.data$sym) %>%
  dplyr::summarise(reference_marker_panels = paste(sort(unique(.data$marker_set)),
                                                   collapse = ";"),
                   reference_cell_classes = paste(sort(unique(stats::na.omit(.data$cell_class))),
                                                  collapse = ";"),
                   .groups = "drop")

# ============================== PART 17 + 18: SUS-RES overlay (JOINED only)

message("Joining canonical SUS-RES statistics")
overlay <- cand %>%
  dplyr::select(dplyr::any_of(c(
    "dataset", "ProteinGroupID",
    "sus_res_fdr05_any_context", "sus_res_n_spatial_contexts_fdr05",
    "sus_res_n_spatial_contexts_tested", "sus_res_strongest_spatial_unit",
    "sus_res_strongest_log2FC", "sus_res_strongest_raw_p",
    "sus_res_strongest_BH_FDR", "sus_res_min_BH_FDR",
    "sus_res_max_abs_log2FC", "sus_res_median_log2FC",
    "sus_res_majority_direction", "sus_res_spatially_consistent")))
overlay$sus_res_strongest_spatial_unit_canonical <- sat_canonical_spatial_unit(
  overlay$sus_res_strongest_spatial_unit, overlay$dataset)
overlay$da_provenance <- paste0(
  "joined from ", relative_to(CAND), "; no differential statistic recomputed")

# ------------------------------------------------------------- assemble

atlas <- scope %>%
  dplyr::select(dplyr::any_of(c(
    "dataset", "ProteinGroupID", "GeneSymbol", "ModuleID", "module_display_label",
    "module_supermodule_id", "abs_kME", "gene_level_claim_allowed",
    "candidate_tier", "candidate_tier_all", "phenotype_network_class",
    "is_tier_A1", "is_tier_A2", "is_tier_B", "is_tier_C", "is_tier_D",
    "is_sus_res_fdr_supported", "is_wgcna_candidate"))) %>%
  dplyr::left_join(metrics, by = c("dataset", "ProteinGroupID")) %>%
  dplyr::left_join(emp_ctx, by = c("dataset", "ProteinGroupID")) %>%
  dplyr::left_join(overlay, by = c("dataset", "ProteinGroupID"))
atlas$sym <- toupper(trimws(atlas$GeneSymbol))
atlas <- dplyr::left_join(atlas, ref_lut, by = "sym")
atlas$sym <- NULL
if (!is.null(contrast_support)) {
  atlas <- dplyr::left_join(atlas, contrast_support, by = c("dataset", "ProteinGroupID"))
}

# PART 17: where does the phenotype effect sit relative to baseline identity?
#
# "High affinity" is defined PROSPECTIVELY and within-protein: the upper third
# of that protein's own spatial ranking. It is never tuned per protein.
HIGH_AFFINITY_FRACTION <- 1 / 3
rank_of_unit <- function(rank_order, unit) {
  if (is.na(rank_order) || is.na(unit)) return(NA_integer_)
  u <- strsplit(rank_order, ">", fixed = TRUE)[[1]]
  m <- match(unit, u)
  if (is.na(m)) NA_integer_ else m
}
atlas$baseline_rank_of_strongest_effect_unit <- vapply(
  seq_len(nrow(atlas)),
  function(i) rank_of_unit(atlas$spatial_rank_order[i],
                           atlas$sus_res_strongest_spatial_unit_canonical[i]),
  integer(1))
atlas$baseline_high_affinity_cutoff_rank <- ceiling(
  atlas$n_spatial_units * HIGH_AFFINITY_FRACTION)

atlas$effect_identity_relationship <- dplyr::case_when(
  !(atlas$is_sus_res_fdr_supported %in% TRUE) ~ "not_fdr_supported",
  is.na(atlas$baseline_rank_of_strongest_effect_unit) ~ "insufficient_context",
  atlas$sus_res_n_spatial_contexts_fdr05 >= 3 ~ "spatially_broad",
  atlas$baseline_rank_of_strongest_effect_unit == 1L ~ "effect_at_baseline_peak",
  atlas$baseline_rank_of_strongest_effect_unit <= 2L ~ "effect_in_top2_baseline_units",
  atlas$baseline_rank_of_strongest_effect_unit <=
    atlas$baseline_high_affinity_cutoff_rank ~ "effect_in_high_affinity_unit",
  TRUE ~ "effect_outside_baseline_affinity")
atlas$high_affinity_rule <- sprintf(
  "upper %.0f%% of the protein's own baseline spatial ranking (prespecified, not tuned)",
  100 * HIGH_AFFINITY_FRACTION)

# PART 18: baseline location and phenotype direction are reported SEPARATELY.
atlas$directional_statement <- ifelse(
  atlas$is_sus_res_fdr_supported %in% TRUE & !is.na(atlas$peak_unit),
  sprintf("baseline peak %s (tau=%.2f); SUS-RES %s in %s (log2FC=%.2f, FDR=%.3g)",
          atlas$peak_unit, atlas$spatial_tau,
          ifelse(atlas$sus_res_strongest_log2FC < 0, "lower in SUS", "higher in SUS"),
          atlas$sus_res_strongest_spatial_unit_canonical,
          atlas$sus_res_strongest_log2FC, atlas$sus_res_strongest_BH_FDR),
  NA_character_)
atlas$interpretation_caveat <- paste0(
  "Baseline abundance and phenotype effect size are separate quantities. ",
  "A high baseline value does not imply a large effect, and a protein should ",
  "not be described as belonging to a unit unless its spatial specificity ",
  "supports it.")
atlas$contract_version <- sat_contract_version()

root <- OUT()
write_csv_safe(atlas, file.path(root, "protein_spatial_cell_affinity.csv"))
write_csv_safe(baseline, file.path(root, "protein_baseline_spatial_profile.csv"))
write_csv_safe(atlas[atlas$is_sus_res_fdr_supported %in% TRUE, , drop = FALSE],
               file.path(root, "protein_sus_res_fdr_supported_atlas.csv"))

cat("\n===== Protein spatial / cell atlas =====\n")
cat(sprintf("scope: %d dataset x ProteinGroupID (%d FDR-supported, %d candidates)\n",
            nrow(atlas), sum(atlas$is_sus_res_fdr_supported %in% TRUE),
            sum(atlas$is_wgcna_candidate %in% TRUE)))
print(table(atlas$dataset))
cat("\ncandidate tiers:\n"); print(table(atlas$candidate_tier, useNA = "ifany"))
cat("\neffect-at-identity classification (FDR-supported only):\n")
print(table(atlas$effect_identity_relationship[atlas$is_sus_res_fdr_supported %in% TRUE]))
cat("\nOutputs:", relative_to(root), "\n")
cat("Baseline identity from CON only; SUS-RES joined as a downstream overlay.\n")
