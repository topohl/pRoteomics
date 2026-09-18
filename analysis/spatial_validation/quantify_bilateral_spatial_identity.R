#!/usr/bin/env Rscript
#
#
# Runs the SAME anatomical contrast definitions - from the single registry in
# R/spatial/control_spatial_identity_utils.R - three times per dataset: on left-only
# samples, on right-only samples, and on the bilateral data. The canonical
# Stage-09 model is not rewritten and not replaced; this sits beside it.
#
# WHAT IS BEING ASKED
#   Not "is the effect significant on both sides?". Requiring independent
#   significance in each hemisphere halves the samples behind each estimate and
#   then penalises the result for it. The question is whether the anatomical
#   effect REPRODUCES in magnitude and direction on the opposite side of the
#   same brain. The both-sides-significant fraction is emitted only as a
#   sensitivity descriptor.
#
# THIS IS NOT INDEPENDENT REPLICATION
#   Both sides come from the same CON animals. n is the number of animals.
#
# PHENOTYPE-BLIND: CON animals only, exactly as the canonical workflow.
#
# USAGE
#   Rscript analysis/spatial_validation/quantify_bilateral_spatial_identity.R
#   Rscript analysis/spatial_validation/quantify_bilateral_spatial_identity.R --dry-run
# Script: analysis/spatial_validation/quantify_bilateral_spatial_identity.R
# Stage: networks
# Scope: per_dataset
# Consumes: required data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx; optional results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_modules_long.csv
# Produces: results/spatial_validation/quantify_bilateral_spatial_identity/global/tables/bilateral_spatial_identity_protein_level.csv; results/spatial_validation/quantify_bilateral_spatial_identity/global/tables/bilateral_spatial_identity_summary.csv
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Bilateral validation of anatomical spatial identity.
#  
#  

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/spatial/control_spatial_identity_utils.R")
source("R/data_contracts/spatial_systems_data_utils.R")
source("R/spatial/spatial_systems_bilateral_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))

# Phase 6G.3: destinations resolve through the normalized output contract,
# addressed by this analysis's own identity rather than by the historical
# 11_spatial_systems stage directory. Outputs already written there stay
# exactly where they are and are read, never rewritten.
ANALYSIS_ID <- "quantify_bilateral_spatial_identity"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

SCRIPT_ID <- "analysis/spatial_validation/quantify_bilateral_spatial_identity.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function() {
  d <- CANONICAL_PATHS$tables; dir_create(d); d
}
meta_path <- function(ds) {
  preprocessing_module_score_metadata(ds)
}
DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  inputs <- stats::setNames(lapply(DATASETS, meta_path),
                            paste0("canonical_metadata__", DATASETS))
  cat("[DRY-RUN] Bilateral spatial identity validation (CON only).\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Contrast definitions come from the single registry in ",
      "R/spatial/control_spatial_identity_utils.R.\n", sep = "")
  quit(save = "no", status = 0L)
}

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("missing_required_input: limma is required.", call. = FALSE)
}

# Prepare CON metadata. The canonical helper derives hemisphere from the sample
# name, which does not resolve for microglia; the ReplicateGroup-derived value
# was verified element-wise identical on both neuronal datasets, so it is used
# as the fallback rather than adding a second hemisphere authority.
prepare_con <- function(canonical, dataset) {
  meta <- canonical$meta
  x <- meta[meta$StressGroup == "CON" &
              (is.na(meta$Exclude) | !meta$Exclude), , drop = FALSE]
  if (!nrow(x)) stop("No CON samples for ", dataset, ".", call. = FALSE)
  hemi <- tryCatch(as.character(control_spatial_hemisphere(x$Sample)),
                   error = function(e) NULL)
  if (is.null(hemi)) hemi <- as.character(sps_normalize_hemisphere(x[[sps_hemisphere_source_field()]]))
  x$hemisphere <- factor(hemi, levels = c("L", "R"))
  x$anatomical_unit <- if (identical(dataset, "neuron_neuropil")) {
    paste(toupper(x$Region), toupper(x$Layer), sep = "_")
  } else toupper(x$Region)
  x
}

# Fit the canonical model on a given sample subset and return per-protein
# estimates for every registry contrast.
fit_side <- function(mat, meta, dataset, side_label) {
  d <- control_spatial_design(meta)
  design <- d$design
  corfit <- limma::duplicateCorrelation(mat, design, block = meta$AnimalID)
  fit <- limma::lmFit(mat, design, block = meta$AnimalID,
                      correlation = corfit$consensus)
  unit_levels <- sub("^anatomical_unit_", "",
                     colnames(design)[grepl("^anatomical_unit_", colnames(design))])
  specs <- control_spatial_contrast_registry(dataset, unit_levels)
  out <- list()
  for (nm in names(specs)) {
    v <- control_spatial_contrast_vector(specs[[nm]], colnames(design))
    cf <- limma::eBayes(
      limma::contrasts.fit(fit, matrix(v, ncol = 1,
                                       dimnames = list(names(v), nm))),
      robust = TRUE, trend = TRUE)
    tt <- limma::topTable(cf, number = Inf, sort.by = "none")
    out[[nm]] <- data.frame(
      dataset = dataset, contrast = nm, side = side_label,
      manuscript_locked = specs[[nm]]$manuscript_locked,
      ProteinGroupID = rownames(tt),
      estimate = tt$logFC,
      SE = sqrt(cf$s2.post) * cf$stdev.unscaled[, 1],
      P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
      n_samples = ncol(mat), n_animals = length(unique(meta$AnimalID)),
      hemisphere_in_model = d$hemisphere_included,
      stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

all_fits <- list()
for (ds in DATASETS) {
  inputs <- resolve_dataset_inputs(ds, purpose = "wgcna", script = SCRIPT_ID,
                                   stage = "networks")
  canonical <- qc_load_canonical_expression(inputs$expression_file, meta_path(ds),
                                            dataset = ds, strict = TRUE)
  meta <- prepare_con(canonical, ds)
  keep <- colnames(canonical$mat) %in% meta$Sample
  mat <- canonical$mat[, keep, drop = FALSE]
  meta <- meta[match(colnames(mat), meta$Sample), , drop = FALSE]

  for (side in c("L", "R", "bilateral")) {
    sel <- if (identical(side, "bilateral")) rep(TRUE, nrow(meta)) else
      meta$hemisphere == side
    m_side <- meta[sel, , drop = FALSE]
    x_side <- mat[, sel, drop = FALSE]
    # SIDE PURITY: a one-sided fit must contain only that hemisphere.
    if (!identical(side, "bilateral") &&
        !all(as.character(m_side$hemisphere) == side)) {
      stop("Side-specific fit for ", side, " contains other-side samples.",
           call. = FALSE)
    }
    message(sprintf("%-16s %-9s samples=%-4d animals=%d units=%d",
                    ds, side, ncol(x_side), length(unique(m_side$AnimalID)),
                    length(unique(m_side$anatomical_unit))))
    all_fits[[length(all_fits) + 1L]] <- fit_side(x_side, m_side, ds, side)
  }
}
fits <- dplyr::bind_rows(all_fits)

# --------------------------------------------------------- protein level

wide <- fits %>%
  dplyr::select("dataset", "contrast", "manuscript_locked", "ProteinGroupID",
                "side", "estimate", "SE", "adj.P.Val") %>%
  tidyr::pivot_wider(names_from = "side",
                     values_from = c("estimate", "SE", "adj.P.Val")) %>%
  dplyr::rename(estimate_L = "estimate_L", estimate_R = "estimate_R",
                estimate_bilateral = "estimate_bilateral",
                SE_L = "SE_L", SE_R = "SE_R", SE_bilateral = "SE_bilateral") %>%
  dplyr::mutate(
    sign_L = sign(.data$estimate_L),
    sign_R = sign(.data$estimate_R),
    sign_agreement = .data$sign_L == .data$sign_R,
    abs_L_minus_R = abs(.data$estimate_L - .data$estimate_R))

# gene symbol where a gene-level claim is allowed
sym_lut <- dplyr::bind_rows(lapply(DATASETS, function(ds) {
  p <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds, "modules",
                    "WGCNA_modules_long.csv")
  if (!file.exists(p)) return(NULL)
  z <- as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                     guess_max = Inf))
  cols <- intersect(c("ProteinGroupID", "GeneSymbol", "gene_level_claim_allowed"),
                    names(z))
  unique(cbind(data.frame(dataset = ds, stringsAsFactors = FALSE), z[, cols, drop = FALSE]))
}))
if (!is.null(sym_lut) && nrow(sym_lut)) {
  wide <- dplyr::left_join(wide, sym_lut, by = c("dataset", "ProteinGroupID"))
  wide$GeneSymbol[!(wide$gene_level_claim_allowed %in% TRUE)] <- NA_character_
}

# ------------------------------------------------------------- summary

summary_tbl <- wide %>%
  dplyr::group_by(.data$dataset, .data$contrast, .data$manuscript_locked) %>%
  dplyr::group_modify(~ {
    ag <- sps_paired_agreement(.x$estimate_L, .x$estimate_R)
    ag$n_evaluable_proteins <- ag$n_pairs
    ag$both_sides_FDR05_fraction <- mean(
      .x$adj.P.Val_L < 0.05 & .x$adj.P.Val_R < 0.05, na.rm = TRUE)
    ag$bilateral_FDR05_fraction <- mean(.x$adj.P.Val_bilateral < 0.05, na.rm = TRUE)
    ag
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    sensitivity_note = paste0(
      "both_sides_FDR05_fraction is a SENSITIVITY descriptor only. Each ",
      "one-sided fit uses half the samples, so requiring independent ",
      "significance on both sides measures power, not reproducibility."))

root <- OUT()
write_csv_safe(wide, file.path(root, "bilateral_spatial_identity_protein_level.csv"))
write_csv_safe(summary_tbl, file.path(root, "bilateral_spatial_identity_summary.csv"))

cat("\n===== Bilateral spatial identity =====\n")
for (ds in DATASETS) {
  s <- summary_tbl[summary_tbl$dataset == ds, , drop = FALSE]
  if (!nrow(s)) next
  cat(sprintf("\n--- %s ---\n", ds))
  cat(sprintf("  %-42s %6s %6s %6s %7s %8s\n",
              "contrast", "n", "r", "rho", "slope", "signAgr"))
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-42s %6d %6.3f %6.3f %7.3f %8.3f%s\n",
                substr(s$contrast[i], 1, 42), s$n_evaluable_proteins[i],
                s$pearson_r[i], s$spearman_rho[i], s$regression_slope[i],
                s$sign_agreement_fraction[i],
                if (!s$manuscript_locked[i]) "  [context]" else ""))
  }
}
cat("\nOutputs:", relative_to(root), "\n")
cat("Contrast definitions came from one registry; L, R and bilateral fits share them.\n")
