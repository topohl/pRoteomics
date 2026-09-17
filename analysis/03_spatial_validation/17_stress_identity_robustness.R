#!/usr/bin/env Rscript
#
# Spatial specificity, module distribution, and the stress-vs-baseline-identity
# result recomputed on explicit robustness subsets.
#
# WHAT THIS DOES NOT DO
#   It does not recompute any differential statistic, does not recompute any
#   FDR, and does not remove any protein from the canonical atlas. The
#   classification rule for effect-vs-baseline-identity is NOT reimplemented
#   here: the canonical per-protein labels produced by
#   analysis/03_spatial_validation/09_protein_spatial_cell_atlas.R are read and tabulated
#   over different row subsets, so the rule cannot drift between the headline
#   number and the robustness-qualified number.
#
# THE QUESTION FOR PART 16
#   The headline was 35 of 37 effects outside baseline affinity, with 15 of 37
#   strongest effects at baseline rank 10. Twenty-eight of those 37 are
#   CA2-SLM, and CA2-SLM carries the lowest baseline abundance, the highest
#   pre-imputation missingness, directional SUS missingness, and two QC-failed
#   SUS acquisitions. So: does the "effects preferentially occur outside
#   baseline affinity" pattern survive once imputation- and QC-sensitive
#   CA2-SLM findings are excluded?
#
# USAGE
#   Rscript analysis/03_spatial_validation/17_stress_identity_robustness.R
#   Rscript analysis/03_spatial_validation/17_stress_identity_robustness.R --dry-run

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/spatial/spatial_atlas_utils.R")
source("R/spatial/ca2_slm_robustness_utils.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/03_spatial_validation/17_stress_identity_robustness.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

NEUROPIL_UNITS <- c("CA1_slm", "CA1_so", "CA1_sr", "CA2_slm", "CA2_so", "CA2_sr",
                    "CA3_so", "CA3_sr", "DG_mo", "DG_po")
DA_DIR <- path_processed("02_id_mapping", "mapped", "neuron_neuropil", "forward", "per_file")
RAW <- repo_path("data", "raw", "pg_matrix", "quicksearch.pg_matrix.tsv")
META <- repo_path("data", "metadata", "TPE9_sample_metadata_males.xlsx")
ATLAS <- path_results("tables", "11_spatial_systems", "atlas",
                      "protein_spatial_cell_affinity.csv")
ROB <- path_results("tables", "11_spatial_systems", "ca2_slm_robustness",
                    "CA2_SLM_DAP_robustness.csv")
OUT <- function(...) {
  d <- path_results("tables", "11_spatial_systems", "ca2_slm_robustness")
  dir_create(d); file.path(d, ...)
}

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Spatial specificity, module distribution and stress identity.\n")
  dry_run_inputs(SCRIPT_ID, list(protein_atlas = ATLAS, ca2_slm_robustness = ROB,
                                 neuropil_da_dir = DA_DIR, raw_pg_matrix = RAW,
                                 sample_metadata = META))
  cat("[DRY-RUN] Canonical labels are TABULATED over subsets, never recomputed.\n")
  quit(save = "no", status = 0L)
}

for (p in c(ATLAS, ROB, RAW, META)) {
  if (!file.exists(p)) stop("missing_required_input: ", p, call. = FALSE)
}

atlas <- as.data.frame(readr::read_csv(ATLAS, show_col_types = FALSE,
                                       progress = FALSE, guess_max = Inf))
rob <- as.data.frame(readr::read_csv(ROB, show_col_types = FALSE,
                                     progress = FALSE, guess_max = Inf))
# This script must never write back over its own input: doing so made the
# second run join already-joined columns and produce .x/.y suffixes. The
# augmented table goes to a separate file and script 16 stays canonical.
rob <- rob[, setdiff(names(rob), c("CA2_SLM_log2FC", "next_largest_abs_log2FC",
  "next_largest_unit", "CA2_SLM_rank_by_abs_effect", "CA2_SLM_uniquely_strongest",
  "n_other_units_same_direction", "n_other_units_comparable_magnitude",
  "n_other_units_tested", "spatial_specificity_class", "spatial_specificity_rule")),
  drop = FALSE]
hits <- atlas[atlas$is_sus_res_fdr_supported %in% TRUE, , drop = FALSE]
message("FDR-supported hits in the canonical atlas: ", nrow(hits))

# ================ PART 14: spatial specificity of the robust subset

message("Comparing each protein's effect across all neuropil spatial units")
per_unit <- list()
for (u in NEUROPIL_UNITS) {
  tok <- gsub("_", "", u)
  f <- file.path(DA_DIR, sprintf("%ssus_%sres.csv", tok, tok))
  if (!file.exists(f)) { warning("missing contrast file: ", f); next }
  d <- as.data.frame(readr::read_csv(f, show_col_types = FALSE, progress = FALSE,
                                     guess_max = Inf))
  per_unit[[u]] <- data.frame(original_identifier = d$original_identifier,
                              spatial_unit = u, log2fc = d$log2fc,
                              padj = d$padj, stringsAsFactors = FALSE)
}
unit_long <- dplyr::bind_rows(per_unit)

spec_rows <- list()
for (pid in rob$original_identifier) {
  z <- unit_long[unit_long$original_identifier == pid, , drop = FALSE]
  if (!nrow(z)) next
  ca2 <- z$log2fc[z$spatial_unit == "CA2_slm"]
  others <- z[z$spatial_unit != "CA2_slm", , drop = FALSE]
  others <- others[order(-abs(others$log2fc)), , drop = FALSE]
  z2 <- z[order(-abs(z$log2fc)), , drop = FALSE]
  rank_ca2 <- match("CA2_slm", z2$spatial_unit)
  same_dir <- sum(sign(others$log2fc) == sign(ca2), na.rm = TRUE)
  # prespecified: "comparable elsewhere" means another unit reaches at least
  # half the CA2-SLM magnitude in the same direction
  comparable <- sum(abs(others$log2fc) >= 0.5 * abs(ca2) &
                      sign(others$log2fc) == sign(ca2), na.rm = TRUE)
  cls <- if (rank_ca2 == 1L && comparable == 0L) "CA2_SLM_selective" else
    if (rank_ca2 == 1L && comparable <= 2L) "CA2_SLM_strongest_but_broader" else
      if (rank_ca2 == 1L) "broad_same_direction" else "not_CA2_SLM_specific"
  spec_rows[[length(spec_rows) + 1L]] <- data.frame(
    original_identifier = pid,
    CA2_SLM_log2FC = ca2,
    next_largest_abs_log2FC = others$log2fc[1],
    next_largest_unit = others$spatial_unit[1],
    CA2_SLM_rank_by_abs_effect = as.integer(rank_ca2),
    CA2_SLM_uniquely_strongest = rank_ca2 == 1L,
    n_other_units_same_direction = as.integer(same_dir),
    n_other_units_comparable_magnitude = as.integer(comparable),
    n_other_units_tested = nrow(others),
    spatial_specificity_class = cls,
    spatial_specificity_rule = paste0(
      "comparable means another unit reaches at least half the CA2-SLM ",
      "magnitude in the same direction; prespecified, descriptive only"),
    stringsAsFactors = FALSE)
}
spec <- dplyr::bind_rows(spec_rows)
rob2 <- dplyr::left_join(rob, spec, by = "original_identifier")

# =================== observation status for ALL 37 hits at their own unit

message("Deriving observation status for every FDR-supported hit")
meta <- as.data.frame(readxl::read_excel(META))
meta$AnimalID <- as.character(meta$AnimalID)
meta <- meta[!(meta$exclude %in% TRUE), , drop = FALSE]
meta$StressGroup <- csr_expgroup_to_stress(meta$ExpGroup)

# The atlas is keyed on ProteinGroupID and carries no raw-matrix identifier.
# The canonical DAP membership table carries BOTH, plus the spatial unit each
# hit was called in, so it is the bridge to the raw matrix.
MEMBER <- path_results("source_data", "04_differential_expression_enrichment",
                       "sus_res_spatial_dap_atlas", "global",
                       "sus_res_dap_membership.csv")
if (!file.exists(MEMBER)) {
  stop("missing_required_input: DAP membership bridge: ", MEMBER, call. = FALSE)
}
member <- as.data.frame(readr::read_csv(MEMBER, show_col_types = FALSE,
                                        progress = FALSE, guess_max = Inf))
hits$original_identifier <- member$original_identifier[
  match(paste(hits$dataset, hits$ProteinGroupID),
        paste(member$dataset, member$ProteinGroupID))]
hits$dap_spatial_unit <- member$spatial_unit[
  match(paste(hits$dataset, hits$ProteinGroupID),
        paste(member$dataset, member$ProteinGroupID))]
if (any(is.na(hits$original_identifier))) {
  stop("unbridged FDR-supported hit(s): ",
       sum(is.na(hits$original_identifier)), call. = FALSE)
}

rawm <- utils::read.delim(RAW, check.names = FALSE, stringsAsFactors = FALSE)
obs_rows <- list()
for (ds in unique(hits$dataset)) {
  md <- meta[meta$celltype_layer == ds, , drop = FALSE]
  cols <- intersect(md$sample_id, names(rawm))
  m <- as.matrix(rawm[, cols, drop = FALSE]); storage.mode(m) <- "double"
  keepr <- rowMeans(is.na(m)) <= 0.7
  m <- m[keepr, , drop = FALSE]
  rownames(m) <- rawm$Protein.Names[keepr]
  md <- md[match(cols, md$sample_id), , drop = FALSE]
  md$region_layer <- paste0(md$region, "_", md$layer)

  h <- hits[hits$dataset == ds, , drop = FALSE]
  for (i in seq_len(nrow(h))) {
    key <- h$original_identifier[i]
    if (!key %in% rownames(m)) {
      obs_rows[[length(obs_rows) + 1L]] <- data.frame(
        dataset = ds, ProteinGroupID = h$ProteinGroupID[i],
        n_observed_SUS = NA_integer_, n_observed_RES = NA_integer_,
        n_samples_SUS = NA_integer_, n_samples_RES = NA_integer_,
        hit_fully_observed = NA, stringsAsFactors = FALSE)
      next
    }
    # neuropil units are Region_Layer; soma and microglia units are bare
    # Region, because their layer token is the compartment, not an axis
    u <- h$dap_spatial_unit[i]
    sel <- if (grepl("_", u, fixed = TRUE)) md$region_layer == u else md$region == u
    if (!any(sel)) {
      stop("no sample matched spatial unit '", u, "' in ", ds, call. = FALSE)
    }
    sus <- md$sample_id[sel & md$StressGroup == "SUS"]
    res <- md$sample_id[sel & md$StressGroup == "RES"]
    ov <- !is.na(m[key, , drop = TRUE])
    obs_rows[[length(obs_rows) + 1L]] <- data.frame(
      dataset = ds, ProteinGroupID = h$ProteinGroupID[i],
      n_observed_SUS = as.integer(sum(ov[sus])),
      n_observed_RES = as.integer(sum(ov[res])),
      n_samples_SUS = length(sus), n_samples_RES = length(res),
      hit_fully_observed = sum(ov[sus]) == length(sus) &&
        sum(ov[res]) == length(res), stringsAsFactors = FALSE)
  }
}
hit_obs <- dplyr::bind_rows(obs_rows)
hits <- dplyr::left_join(hits, hit_obs, by = c("dataset", "ProteinGroupID"))

# ================== PART 16: stress identity over explicit subsets

qualified <- rob2$ProteinGroupID[rob2$CA2_SLM_robustness_class %in%
  c("robust_to_missingness_and_QC", "supported_but_QC_sensitive")]
robust_only <- rob2$ProteinGroupID[rob2$CA2_SLM_robustness_class ==
                                     "robust_to_missingness_and_QC"]
ca2_ids <- rob2$ProteinGroupID

summarise_subset <- function(label, sel, note) {
  z <- hits[sel, , drop = FALSE]
  rel <- z$effect_identity_relationship
  rk <- z$baseline_rank_of_strongest_effect_unit
  data.frame(
    subset = label, n_hits = nrow(z),
    n_CA2_SLM = sum(z$sus_res_strongest_spatial_unit_canonical == "ca2_slm", na.rm = TRUE),
    effect_at_baseline_peak = sum(rel == "effect_at_baseline_peak", na.rm = TRUE),
    effect_in_top2_baseline_units = sum(rel == "effect_in_top2_baseline_units", na.rm = TRUE),
    effect_in_high_affinity_unit = sum(rel == "effect_in_high_affinity_unit", na.rm = TRUE),
    effect_outside_baseline_affinity = sum(rel == "effect_outside_baseline_affinity", na.rm = TRUE),
    fraction_outside_baseline_affinity = if (nrow(z))
      mean(rel == "effect_outside_baseline_affinity", na.rm = TRUE) else NA_real_,
    n_at_baseline_rank_10 = sum(rk == 10L, na.rm = TRUE),
    median_baseline_rank = if (nrow(z)) stats::median(rk, na.rm = TRUE) else NA_real_,
    classification_source = "canonical labels from 09_protein_spatial_cell_atlas.R, tabulated not recomputed",
    subset_definition = note, stringsAsFactors = FALSE)
}

is_ca2 <- hits$ProteinGroupID %in% ca2_ids
identity_cmp <- dplyr::bind_rows(
  summarise_subset("all_canonical_FDR_supported_hits", rep(TRUE, nrow(hits)),
                   "every FDR-supported protein in the canonical atlas; the published headline"),
  summarise_subset("fully_observed_hits", hits$hit_fully_observed %in% TRUE,
                   "no imputed value in any SUS or RES sample of the protein's own strongest unit"),
  summarise_subset("CA2_SLM_robustness_qualified",
                   (!is_ca2) | hits$ProteinGroupID %in% qualified,
                   "CA2-SLM hits restricted to robust_to_missingness_and_QC or supported_but_QC_sensitive; non-CA2-SLM hits retained"),
  summarise_subset("CA2_SLM_robust_only", (!is_ca2) | hits$ProteinGroupID %in% robust_only,
                   "CA2-SLM hits restricted to robust_to_missingness_and_QC only; non-CA2-SLM hits retained"),
  summarise_subset("non_CA2_SLM_hits_only", !is_ca2,
                   "the 9 FDR-supported hits whose strongest unit is not CA2-SLM"),
  summarise_subset("claimable_across_datasets",
                   (!is_ca2 & hits$hit_fully_observed %in% TRUE) |
                     hits$ProteinGroupID %in% robust_only,
                   "robust CA2-SLM hits plus fully observed non-CA2-SLM hits")
)

# ====================== PART 15: module distribution before and after

module_cmp <- dplyr::bind_rows(
  rob2 |> dplyr::count(ModuleID, name = "n_all_28_CA2_SLM_DAPs") |>
    dplyr::mutate(subset = "A_all_canonical_28"),
  rob2[rob2$ProteinGroupID %in% qualified, ] |>
    dplyr::count(ModuleID, name = "n_all_28_CA2_SLM_DAPs") |>
    dplyr::mutate(subset = "B_robustness_qualified"),
  rob2[rob2$ProteinGroupID %in% robust_only, ] |>
    dplyr::count(ModuleID, name = "n_all_28_CA2_SLM_DAPs") |>
    dplyr::mutate(subset = "C_robust_only")
) |>
  dplyr::rename(n_proteins = "n_all_28_CA2_SLM_DAPs") |>
  tidyr::pivot_wider(names_from = "subset", values_from = "n_proteins",
                     values_fill = 0L) |>
  as.data.frame()
module_cmp$enrichment_test <- paste0(
  "NOT TESTED: no prespecified universe exists for a post hoc module ",
  "enrichment test on a phenotype-selected subset of 28 proteins")

# ================= PART 17: claimability annotation for the atlas

claim <- rob2 |>
  dplyr::transmute(
    dataset = "neuron_neuropil",
    ProteinGroupID = .data$ProteinGroupID,
    CA2_SLM_robustness_class = .data$CA2_SLM_robustness_class,
    imputation_dependence = .data$imputation_dependence_class,
    LOO_sign_stability = .data$loo_sign_stable,
    fully_observed = .data$fully_observed,
    spatial_specificity_class = .data$spatial_specificity_class,
    QC_claimability = dplyr::case_when(
      .data$CA2_SLM_robustness_class == "robust_to_missingness_and_QC" ~ "claimable",
      .data$CA2_SLM_robustness_class == "supported_but_QC_sensitive" ~ "claimable_with_caveat",
      .data$CA2_SLM_robustness_class == "insufficient_observed_data" ~ "not_evaluable",
      TRUE ~ "not_claimable"),
    claimable_for_biological_interpretation =
      .data$CA2_SLM_robustness_class == "robust_to_missingness_and_QC",
    robustness_audit_status = "audited_ca2_slm_v1",
    robustness_reason = .data$classification_reason)

# ----------------------------------------------------------------- write

write_csv_safe(rob2, OUT("CA2_SLM_DAP_robustness_annotated.csv"))
write_csv_safe(spec, OUT("CA2_SLM_spatial_specificity.csv"))
write_csv_safe(unit_long[unit_long$original_identifier %in% rob$original_identifier, ],
               OUT("CA2_SLM_effect_across_neuropil_units.csv"))
write_csv_safe(identity_cmp, OUT("stress_identity_robustness_comparison.csv"))
write_csv_safe(module_cmp, OUT("CA2_SLM_module_distribution_comparison.csv"))
write_csv_safe(hit_obs, OUT("fdr_supported_hit_observation_status.csv"))
write_csv_safe(claim, path_results("tables", "11_spatial_systems", "atlas",
                                   "protein_claimability_annotation.csv"))

cat("\n===== Stress identity and specificity =====\n")
cat("\n--- spatial specificity of the 28 CA2-SLM DAPs ---\n")
print(table(spec$spatial_specificity_class))

cat("\n--- module distribution ---\n")
print(module_cmp[order(-module_cmp$A_all_canonical_28),
                 c("ModuleID", "A_all_canonical_28", "B_robustness_qualified",
                   "C_robust_only")], row.names = FALSE)

cat("\n--- stress vs baseline identity, by subset ---\n")
for (i in seq_len(nrow(identity_cmp))) {
  cat(sprintf("  %-34s n=%-3d outside=%-3d (%.0f%%)  rank10=%-3d  medRank=%.1f\n",
              substr(identity_cmp$subset[i], 1, 34), identity_cmp$n_hits[i],
              identity_cmp$effect_outside_baseline_affinity[i],
              100 * identity_cmp$fraction_outside_baseline_affinity[i],
              identity_cmp$n_at_baseline_rank_10[i],
              identity_cmp$median_baseline_rank[i]))
}
cat("\n--- claimability ---\n")
print(table(claim$QC_claimability))
cat("\nOutputs:", relative_to(dirname(OUT("x"))), "\n")
