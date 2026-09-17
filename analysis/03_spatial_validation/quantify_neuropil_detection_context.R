#!/usr/bin/env Rscript
#
#
# 28 of the 31 neuropil SUS-RES DAPs fall in CA2_slm, and 7 of 10 neuropil units
# have none at all. Three explanations must be separated before CA2_slm is
# called a biological hotspot:
#
#   A genuine spatial concentration of phenotype effects
#   B simply higher baseline abundance / detectability in CA2_slm
#   C broader context-specific statistical power (smaller residual variance,
#     so the same true effect clears FDR there and nowhere else)
#
# WHAT IS USED AS A DENOMINATOR
#   All proteins TESTED in each unit, taken from the canonical DAP-count table -
#   never the candidate subset, which is itself phenotype-selected and would
#   build the conclusion into the denominator.
#
# Every statistic here is READ from canonical DA outputs. Nothing is refitted.
#
# USAGE
#   Rscript analysis/03_spatial_validation/quantify_neuropil_detection_context.R
# Script: analysis/03_spatial_validation/quantify_neuropil_detection_context.R
# Stage: networks
# Scope: global
# Consumes: required results/tables/04_differential_expression_enrichment/sus_res_spatial_dap_atlas/global/sus_res_dap_counts.csv; optional results/tables/10_biological_integration/wgcna_candidate_protein_shortlist/neuron_neuropil/wgcna_candidate_proteins_long.csv
# Produces: results/tables/11_spatial_systems/atlas/neuropil_spatial_detection_context.csv
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Neuropil spatial detection context: is the CA2-SLM hit burden biological?

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/data_contracts/spatial_systems_data_utils.R")
source("R/spatial/spatial_atlas_utils.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/03_spatial_validation/quantify_neuropil_detection_context.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function() {
  d <- path_results("tables", "11_spatial_systems", "atlas"); dir_create(d); d
}
DAP <- path_results("tables", "04_differential_expression_enrichment",
                    "sus_res_spatial_dap_atlas", "global", "sus_res_dap_counts.csv")
LONG <- function(ds) path_results("tables", "10_biological_integration",
                                  "wgcna_candidate_protein_shortlist", ds,
                                  "wgcna_candidate_proteins_long.csv")
DS <- "neuron_neuropil"

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Neuropil spatial detection context (CA2-SLM audit).\n")
  dry_run_inputs(SCRIPT_ID, list(dap_counts = DAP, candidate_long = LONG(DS)))
  cat("[DRY-RUN] Denominator is ALL tested proteins, never the candidate subset.\n")
  quit(save = "no", status = 0L)
}

read_req <- function(p, lbl) {
  if (!file.exists(p)) stop("missing_required_input: ", lbl, ": ", p, call. = FALSE)
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
}

# ---- canonical per-unit DAP counts -----------------------------------------
# JOIN ON route_unit. The sibling `spatial_unit` column is collapsed to bare
# Region for soma and microglia, so joining on it silently loses 8 of 18 units.
dap <- read_req(DAP, "canonical DAP counts")
dap <- dap[dap$dataset == DS, , drop = FALSE]

# ---- canonical per-unit DA statistics --------------------------------------
long <- read_req(LONG(DS), "candidate per-unit DA")
sus_res <- long[long$contrast == "SUS - RES", , drop = FALSE]
if (!nrow(sus_res)) stop("No 'SUS - RES' rows found; contrast label changed?", call. = FALSE)

da_stats <- sus_res %>%
  dplyr::group_by(spatial_unit = .data$spatial_unit) %>%
  dplyr::summarise(
    n_candidate_proteins_with_da = dplyr::n_distinct(.data$ProteinGroupID),
    median_abs_log2FC = stats::median(abs(.data$log2FC), na.rm = TRUE),
    p90_abs_log2FC = unname(stats::quantile(abs(.data$log2FC), 0.9, na.rm = TRUE)),
    median_raw_p = stats::median(.data$raw_p, na.rm = TRUE),
    fraction_raw_p_lt_005 = mean(.data$raw_p < 0.05, na.rm = TRUE),
    min_BH_FDR = suppressWarnings(min(.data$BH_FDR, na.rm = TRUE)),
    # A precision proxy derived from PUBLISHED statistics only: for a two-sided
    # test, |t| = qnorm(1 - p/2) approximately, so SE ~ |log2FC| / |t|. Used
    # only to COMPARE units, never as a reported standard error.
    median_precision_proxy_SE = stats::median(
      abs(.data$log2FC) / pmax(stats::qnorm(1 - pmin(pmax(.data$raw_p, 1e-12), 1) / 2), 1e-6),
      na.rm = TRUE),
    .groups = "drop")

# ---- CON baseline abundance per unit, over ALL measured proteins ------------
message("Measuring CON baseline abundance per spatial unit (all measured proteins)")
inputs <- resolve_dataset_inputs(DS, purpose = "wgcna", script = SCRIPT_ID,
                                 stage = "networks")
md <- path_processed("01_preprocessing", "06_merged_metadata_module_score", DS,
                     "sample_metadata_merged_clean_for_module_scores.xlsx")
canonical <- qc_load_canonical_expression(inputs$expression_file, md,
                                          dataset = DS, strict = TRUE)
lv <- sps_build_spatial_levels(DS, canonical = canonical)
con <- lv$level2$meta$StressGroup == "CON"
b <- lv$level2$mat[, con, drop = FALSE]
meta <- lv$level2$meta[con, , drop = FALSE]

base_rows <- lapply(sort(unique(meta$SpatialUnit)), function(u) {
  cols <- meta$SpatialUnit == u
  v <- b[, cols, drop = FALSE]
  per_protein <- rowMeans(v, na.rm = TRUE)
  # between-animal spread at this unit, averaged over proteins: the quantity a
  # per-unit power difference would show up in
  spread <- apply(v, 1, stats::sd, na.rm = TRUE)
  data.frame(
    # normalise here too: the level-2 unit is "CA1_slm" while the canonical DA
    # tables use "ca1_slm", and an unnormalised join returns all-NA silently
    spatial_unit_canonical = sat_canonical_spatial_unit(u, DS),
    n_measured_proteins = nrow(v),
    median_baseline_abundance = stats::median(per_protein, na.rm = TRUE),
    mean_baseline_abundance = mean(per_protein, na.rm = TRUE),
    q10_baseline_abundance = unname(stats::quantile(per_protein, 0.1, na.rm = TRUE)),
    q90_baseline_abundance = unname(stats::quantile(per_protein, 0.9, na.rm = TRUE)),
    median_between_animal_sd = stats::median(spread, na.rm = TRUE),
    fraction_missing = mean(!is.finite(v)),
    stringsAsFactors = FALSE)
})
baseline <- dplyr::bind_rows(base_rows)

# ---- assemble ---------------------------------------------------------------
dap$spatial_unit_canonical <- sat_canonical_spatial_unit(dap$route_unit, DS)
da_stats$spatial_unit_canonical <- sat_canonical_spatial_unit(da_stats$spatial_unit, DS)

ctx <- dap %>%
  dplyr::select("route_unit", "spatial_unit_canonical",
                n_proteins_tested = "n_tested_ProteinGroupID",
                n_DAP_FDR05 = "n_DAP_FDR05",
                n_higher_in_SUS = "n_higher_in_SUS",
                n_higher_in_RES = "n_higher_in_RES",
                fraction_DAP_of_tested = "fraction_DAP_of_tested") %>%
  dplyr::left_join(da_stats, by = "spatial_unit_canonical") %>%
  dplyr::left_join(baseline, by = "spatial_unit_canonical") %>%
  dplyr::arrange(dplyr::desc(.data$n_DAP_FDR05))
ctx$dataset <- DS
ctx$denominator_definition <- "all ProteinGroupIDs tested in the unit (canonical DAP counts), not the candidate subset"
ctx$precision_proxy_note <- "SE proxy = |log2FC| / qnorm(1 - p/2), derived from published statistics for cross-unit comparison only"

root <- OUT()
write_csv_safe(ctx, file.path(root, "neuropil_spatial_detection_context.csv"))

cat("\n===== Neuropil spatial detection context =====\n")
cat(sprintf("%-9s %8s %6s %9s %9s %9s %9s %9s\n", "unit", "tested", "DAPs",
            "medAbsFC", "medRawP", "fracP<.05", "medBaseAb", "medBtwSD"))
for (i in seq_len(nrow(ctx))) {
  cat(sprintf("%-9s %8d %6d %9.4f %9.4f %9.4f %9.4f %9.4f\n",
              ctx$spatial_unit_canonical[i], ctx$n_proteins_tested[i],
              ctx$n_DAP_FDR05[i], ctx$median_abs_log2FC[i], ctx$median_raw_p[i],
              ctx$fraction_raw_p_lt_005[i], ctx$median_baseline_abundance[i],
              ctx$median_between_animal_sd[i]))
}
cat("\nInterpretation inputs:\n")
cat(sprintf("  proteins tested identical across units: %s\n",
            length(unique(ctx$n_proteins_tested)) == 1L))
ca2 <- ctx[ctx$spatial_unit_canonical == "ca2_slm", , drop = FALSE]
oth <- ctx[ctx$spatial_unit_canonical != "ca2_slm", , drop = FALSE]
if (nrow(ca2)) {
  cat(sprintf("  CA2_slm baseline abundance rank: %d of %d (higher = more abundant)\n",
              rank(-ctx$median_baseline_abundance)[ctx$spatial_unit_canonical == "ca2_slm"],
              nrow(ctx)))
  cat(sprintf("  CA2_slm between-animal SD rank : %d of %d (1 = least variable)\n",
              rank(ctx$median_between_animal_sd)[ctx$spatial_unit_canonical == "ca2_slm"],
              nrow(ctx)))
  cat(sprintf("  CA2_slm median |log2FC| vs other units: %.4f vs %.4f\n",
              ca2$median_abs_log2FC, stats::median(oth$median_abs_log2FC, na.rm = TRUE)))
  cat(sprintf("  CA2_slm fraction raw p<0.05 vs other units: %.4f vs %.4f\n",
              ca2$fraction_raw_p_lt_005,
              stats::median(oth$fraction_raw_p_lt_005, na.rm = TRUE)))
}
cat("\nOutput:", relative_to(file.path(root, "neuropil_spatial_detection_context.csv")), "\n")
cat("CA2-SLM is NOT called a biological hotspot by this script; the evidence is tabulated.\n")
