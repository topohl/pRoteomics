#!/usr/bin/env Rscript
#
#
# Inhibitory presynaptic share of every neuropil lamina (Figure 2l,
# 2026-10-28; 2026-10-29: carry-over rule shared with j and k, animals named in
# the outputs), control animals; stressed animals as replication.
#
# INHIBITORY SHARE: in every neuropil acquisition, the mean log2 abundance of
# the GABAergic presynaptic proteins GAD1, GAD2 and VGAT (Slc32a1) minus the
# log2 abundance of VGLUT1 (Slc17a7). A within-acquisition ratio, so it does
# not depend on normalisation; label-free, so it compares laminae, it is not a
# synapse count. All four proteins must be observed in the acquisition.
# Animal mean over retained acquisitions (hemispheres). Pre-declared pairs
# (strata within each region; the same stratum across regions): each pair is
# tested on the per-animal differences of the CON animals (paired t, df = 2),
# BH over the pairs. The pooled-error model (lamina + animal, one residual for
# every pair) is kept as a sensitivity only: the between-animal spread of the
# differences varies more than twentyfold between pairs (independent review
# 2026-10-28). The stressed animals are tested the same way as a replication,
# never pooled: every stressed animal with both laminae of a pair (six; five
# where DG po is missing).
# Sensitivities and context: GAD2 + VGAT only (GAD1 also marks interneuron
# somata and runs opposite in the hilus); VGLUT2 (Slc17a6, detected mainly in
# SLM and DG mo, mainly reuniens and supramammillary terminals; the VGLUT1
# denominator misses it, so the share there is an upper bound);
# GPHN and NLGN2 vs PSD-95 (Dlg4) and GluA2 (Gria2) as the postsynaptic
# counterpart (descriptive).
#
# DG po acquisitions with somatic carry-over are excluded by the rule of
# map_input_territories (R/spatial/neuropil_lamina_utils.R).
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/compare_lamina_synapses.R
#   Rscript analysis/spatial_validation/compare_lamina_synapses.R --dry-run
# Script: analysis/spatial_validation/compare_lamina_synapses.R
# Stage: networks
# Scope: global
# Consumes: required results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds
# Produces: results/spatial_validation/compare_lamina_synapses/global/tables/inhibitory_share_by_animal.csv; results/spatial_validation/compare_lamina_synapses/global/tables/inhibitory_share_contrasts.csv; +1 more
# Dataset behavior: global; reads the neuropil acquisitions of the joint bundle.
# Notes: Laminar inhibitory presynaptic share (GAD1/GAD2/VGAT vs VGLUT1), pre-declared lamina pairs paired over CON animals (Figure 2l); stressed animals as replication.
#
#

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/qc/compartment_abundance_utils.R")
source(repo_path("R", "preprocessing_paths.R"))
source(repo_path("R", "spatial_systems_paths.R"))
source("R/spatial/neuropil_lamina_utils.R")

ANALYSIS_ID <- "compare_lamina_synapses"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

SCRIPT_ID <- "analysis/spatial_validation/compare_lamina_synapses.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

JOINT_BUNDLE <- preprocessing_find("joint_compartment_qc_matrices.rds", owner = "build_joint_protigy_input",
                                   legacy_stage = "01_preprocessing", legacy_substep = "joint_compartment_qc",
                                   scope = "global", child = "models", legacy_family = "processed")
OUT <- list(
  share_acquisition = file.path(CANONICAL_PATHS$tables, "inhibitory_share_by_acquisition.csv"),
  share_animal      = file.path(CANONICAL_PATHS$tables, "inhibitory_share_by_animal.csv"),
  share_contrasts   = file.path(CANONICAL_PATHS$tables, "inhibitory_share_contrasts.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")
CONTRAST_GROUP <- "CON"
LAMINAE <- c("CA1_so", "CA1_sr", "CA1_slm", "CA2_so", "CA2_sr", "CA2_slm", "CA3_so", "CA3_sr", "DG_mo", "DG_po")
SOMA_UNITS <- c("CA1_sp", "CA2_sp", "CA3_sp", "DG_sg")
INHIBITORY <- c("Gad1", "Gad2", "Slc32a1")
EXCITATORY <- "Slc17a7"
POST_INHIBITORY <- c("Gphn", "Nlgn2")
POST_EXCITATORY <- c("Dlg4", "Gria2")
SHARE_PAIRS <- list(c("CA1_slm", "CA1_sr"), c("CA1_slm", "CA1_so"), c("CA1_sr", "CA1_so"),
                    c("CA2_slm", "CA2_sr"), c("CA2_slm", "CA2_so"), c("CA2_sr", "CA2_so"),
                    c("CA3_sr", "CA3_so"), c("DG_mo", "DG_po"),
                    c("CA2_sr", "CA1_sr"), c("CA3_sr", "CA1_sr"), c("CA2_sr", "CA3_sr"),
                    c("CA2_so", "CA1_so"), c("CA3_so", "CA1_so"), c("CA2_so", "CA3_so"),
                    c("CA2_slm", "CA1_slm"))

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Laminar inhibitory share (CON animals).\n")
  dry_run_inputs(SCRIPT_ID, list(joint_bundle = JOINT_BUNDLE))
  quit(save = "no", status = 0L)
}
if (!file.exists(JOINT_BUNDLE)) stop("missing_required_input: ", JOINT_BUNDLE, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))

rec <- ca_validate_bundle_reconstruction(readRDS(JOINT_BUNDLE))
meta <- rec$metadata
cols <- rec$columns
samples <- as.character(meta[[cols$sample]])
info <- data.frame(Sample = samples, AnimalID = as.character(meta[[cols$animal]]),
                   group = ca_normalize_group(meta[[cols$group]]),
                   hemisphere = ca_normalize_hemisphere(meta[[cols$hemisphere]]),
                   unit = paste(meta[[cols$region]], meta[[cols$layer]], sep = "_"), stringsAsFactors = FALSE)
mat <- rec$normalized_nonimputed[, samples, drop = FALSE]
obs <- rec$observed[rownames(mat), samples, drop = FALSE] & is.finite(mat)
ft <- rec$feature_table[match(rownames(mat), rec$feature_table$ProteinGroupID), , drop = FALSE]
symbol <- as.character(ft$official_gene_symbol)
use <- ft$joint_qc_eligible %in% TRUE & ft$gene_level_claim_allowed %in% TRUE & !is.na(symbol) & nzchar(symbol) &
  !duplicated(symbol)
mat <- mat[use, , drop = FALSE]
obs <- obs[use, , drop = FALSE]
rownames(mat) <- rownames(obs) <- symbol[use]
# soma carry-over into DG po: the rule of map_input_territories (R/spatial/neuropil_lamina_utils.R)
carry <- dg_po_carryover(info, obs, LAMINAE, SOMA_UNITS)
info <- info[info$unit %in% LAMINAE & !info$Sample %in% carry$Sample[carry$excluded], c("Sample", "AnimalID", "group",
                                                                                        "hemisphere", "unit"), drop = FALSE]
names(info)[names(info) == "unit"] <- "lamina"
mat <- mat[, info$Sample, drop = FALSE]
obs <- obs[, info$Sample, drop = FALSE]
need <- c(INHIBITORY, EXCITATORY, POST_INHIBITORY, POST_EXCITATORY)
if (!all(need %in% rownames(mat))) stop("Synapse proteins absent: ", paste(setdiff(need, rownames(mat)), collapse = ", "), call. = FALSE)

# ---- inhibitory share --------------------------------------------------------------
ratio <- function(num, den) colMeans(mat[num, , drop = FALSE]) - colMeans(mat[den, , drop = FALSE])
complete <- function(g) colSums(!obs[g, , drop = FALSE]) == 0
acq <- info
acq$inhibitory_presynaptic_share_log2 <- ifelse(complete(c(INHIBITORY, EXCITATORY)), ratio(INHIBITORY, EXCITATORY), NA_real_)
acq$inhibitory_postsynaptic_share_log2 <- ifelse(complete(c(POST_INHIBITORY, POST_EXCITATORY)),
                                                 ratio(POST_INHIBITORY, POST_EXCITATORY), NA_real_)
acq$gad2_vgat_share_log2 <- ifelse(complete(c("Gad2", "Slc32a1", EXCITATORY)), ratio(c("Gad2", "Slc32a1"), EXCITATORY), NA_real_)
acq$vglut2_minus_vglut1_log2 <- if ("Slc17a6" %in% rownames(mat))
  ifelse(complete(c("Slc17a6", EXCITATORY)), ratio("Slc17a6", EXCITATORY), NA_real_) else NA_real_
acq$definition <- paste0("presynaptic: mean log2 ", paste(INHIBITORY, collapse = "/"), " minus log2 ", EXCITATORY,
                         "; postsynaptic: mean log2 ", paste(POST_INHIBITORY, collapse = "/"), " minus mean log2 ",
                         paste(POST_EXCITATORY, collapse = "/"), "; observed values only")
key <- paste(acq$AnimalID, acq$lamina)
animal <- do.call(rbind, lapply(split(acq, key), function(d) {
  data.frame(AnimalID = d$AnimalID[1], group = d$group[1], lamina = d$lamina[1],
             n_acquisitions = sum(is.finite(d$inhibitory_presynaptic_share_log2)),
             inhibitory_presynaptic_share_log2 = mean(d$inhibitory_presynaptic_share_log2, na.rm = TRUE),
             inhibitory_postsynaptic_share_log2 = mean(d$inhibitory_postsynaptic_share_log2, na.rm = TRUE),
             gad2_vgat_share_log2 = mean(d$gad2_vgat_share_log2, na.rm = TRUE),
             vglut2_minus_vglut1_log2 = if (any(is.finite(d$vglut2_minus_vglut1_log2)))
               mean(d$vglut2_minus_vglut1_log2, na.rm = TRUE) else NA_real_,
             vglut2_acquisitions_observed = sum(is.finite(d$vglut2_minus_vglut1_log2)),
             stringsAsFactors = FALSE)
}))
rownames(animal) <- NULL
animal$presynaptic_share_centred <- animal$inhibitory_presynaptic_share_log2 -
  stats::ave(animal$inhibitory_presynaptic_share_log2, animal$AnimalID, FUN = function(x) mean(x, na.rm = TRUE))
animal$postsynaptic_share_centred <- animal$inhibitory_postsynaptic_share_log2 -
  stats::ave(animal$inhibitory_postsynaptic_share_log2, animal$AnimalID, FUN = function(x) mean(x, na.rm = TRUE))
con <- animal[animal$group == CONTRAST_GROUP, , drop = FALSE]
if (!all(table(factor(con$lamina, levels = LAMINAE)) == length(unique(con$AnimalID))) ||
    !all(is.finite(con$inhibitory_presynaptic_share_log2)))
  stop("Every CON animal needs a finite share in every lamina.", call. = FALSE)
# pooled-error model, kept as a sensitivity
fit <- stats::lm(inhibitory_presynaptic_share_log2 ~ lamina + AnimalID, data = con)
dfe <- fit$df.residual
se_pooled <- sqrt(2 * sum(stats::residuals(fit)^2) / dfe / length(unique(con$AnimalID)))
paired <- function(d, a, b, value = "inhibitory_presynaptic_share_log2") {
  x <- d[d$lamina == a, , drop = FALSE]
  y <- d[d$lamina == b, , drop = FALSE]
  y <- y[match(x$AnimalID, y$AnimalID), , drop = FALSE]
  v <- x[[value]] - y[[value]]
  v[is.finite(v)]
}
share_contrasts <- list()
for (set_name in c(CONTRAST_GROUP, "RES+SUS")) {
  gs <- if (set_name == CONTRAST_GROUP) CONTRAST_GROUP else c("RES", "SUS")
  d <- animal[animal$group %in% gs, , drop = FALSE]
  for (p in SHARE_PAIRS) {
    v <- paired(d, p[1], p[2])
    tt <- stats::t.test(v)
    g2 <- paired(d, p[1], p[2], "gad2_vgat_share_log2")
    row <- data.frame(animal_set = set_name, lamina_a = p[1], lamina_b = p[2],
                      comparison_type = if (sub("_.*", "", p[1]) == sub("_.*", "", p[2])) "strata within a region" else
                        "same stratum across regions",
                      estimate = mean(v), std_error = stats::sd(v) / sqrt(length(v)), ci_low = tt$conf.int[1],
                      ci_high = tt$conf.int[2], df = unname(tt$parameter), t_value = unname(tt$statistic),
                      p_value = tt$p.value, n_animals = length(v), animals_positive = sum(v > 0),
                      estimate_gad2_vgat_only = mean(g2),
                      animals = paste(sort(intersect(d$AnimalID[d$lamina == p[1] & is.finite(d$inhibitory_presynaptic_share_log2)],
                                                     d$AnimalID[d$lamina == p[2] & is.finite(d$inhibitory_presynaptic_share_log2)])),
                                      collapse = "; "), stringsAsFactors = FALSE)
    if (set_name == CONTRAST_GROUP) {
      row$pooled_std_error <- se_pooled
      row$pooled_df <- dfe
      row$pooled_p_value <- 2 * stats::pt(-abs(row$estimate / se_pooled), dfe)
    } else {
      row$pooled_std_error <- NA_real_
      row$pooled_df <- NA_real_
      row$pooled_p_value <- NA_real_
    }
    share_contrasts[[length(share_contrasts) + 1L]] <- row
  }
}
share_contrasts <- do.call(rbind, share_contrasts)
share_contrasts$p_adjust <- NA_real_
share_contrasts$pooled_p_adjust <- NA_real_
for (set_name in unique(share_contrasts$animal_set)) {
  i <- share_contrasts$animal_set == set_name
  share_contrasts$p_adjust[i] <- stats::p.adjust(share_contrasts$p_value[i], method = "BH")
  if (set_name == CONTRAST_GROUP)
    share_contrasts$pooled_p_adjust[i] <- stats::p.adjust(share_contrasts$pooled_p_value[i], method = "BH")
}
share_contrasts$group <- share_contrasts$animal_set
share_contrasts$role <- ifelse(share_contrasts$animal_set == CONTRAST_GROUP, "primary (control animals)",
                               "replication (stressed animals, never pooled with the controls)")
share_contrasts$method <- paste0("per-animal difference of the lamina shares (a - b), one-sample t test over animals; BH over ",
                                 "the ", length(SHARE_PAIRS), " pre-declared pairs of each animal set; pooled_*: the ",
                                 "lamina + animal model with one pooled residual (sensitivity, control animals)")

utils::write.csv(acq, OUT$share_acquisition, row.names = FALSE)
utils::write.csv(animal[order(animal$group, animal$AnimalID, match(animal$lamina, LAMINAE)), ], OUT$share_animal, row.names = FALSE)
utils::write.csv(share_contrasts, OUT$share_contrasts, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = list(joint_bundle = JOINT_BUNDLE), outputs = OUT,
                   parameters = list(scale = "joint shared-core normalisation, observed values only (no imputation)",
                                     inhibitory = INHIBITORY, excitatory = EXCITATORY,
                                     postsynaptic_inhibitory = POST_INHIBITORY, postsynaptic_excitatory = POST_EXCITATORY,
                                     share_pairs = vapply(SHARE_PAIRS, paste, "", collapse = " - "),
                                     comparison_group = CONTRAST_GROUP,
                                     comparison_animals = sort(unique(con$AnimalID)),
                                     replication_animals = sort(unique(animal$AnimalID[animal$group != CONTRAST_GROUP])),
                                     carryover_threshold = carry$carryover_threshold[1],
                                     excluded_acquisitions = carry$Sample[carry$excluded], inferential_p_values = TRUE,
                                     provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "The animal is the unit; retained hemisphere acquisitions are averaged within animal; the per-pair animals are in the contrasts table.")
message("Laminar inhibitory share complete")
