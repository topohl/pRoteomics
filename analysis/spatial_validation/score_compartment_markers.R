#!/usr/bin/env Rscript
#
#
# Reference marker classes scored in every acquisition of the three
# measurement compartments, on the joint cross-compartment scale.
#
# MARKER CLASSES: the external fidelity registry as resolved by the
# authoritative 04e workflow (render_compartment_abundance_figures):
# Microglia/PVM markers (Allen cell types), Soma markers (GO, MGI) and Neuropil
# markers (SynGO and GO, MGI); a protein counts when it is a primary-score
# eligible canonical marker of exactly one class (no conflicting classes) that is
# detected in its intended compartment in >= 2 of 3 CON animals
# (v2_04 intended_primary_eligible).
#
# SCALE: raw positive log2 abundance minus the acquisition's joint shared-core
# median offset (the 04e / joint-QC normalisation), expressed relative to the
# joint target median, so 0 is the typical shared-core protein. Observed values
# only: a marker that is not detected in an acquisition is left out of that
# acquisition's mean, never imputed; the number of contributing markers is kept.
#
# UNITS: every acquisition of all nine animals (group is carried), and the
# mean of each animal's acquisitions per compartment.
#
# COMPARISONS (2026-10-28): control animals only. The animal is the unit and
# every compartment is measured in every animal, so each marker class is a
# repeated-measures design on the animal means: a linear model with compartment
# and animal (block) as factors, each pair of compartments compared with the
# pooled error (df = 4 for 3 animals x 3 compartments), BH over all comparisons
# in the table. The within-protein-centred marker comparison stays in 04e.
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/score_compartment_markers.R
#   Rscript analysis/spatial_validation/score_compartment_markers.R --dry-run
# Script: analysis/spatial_validation/score_compartment_markers.R
# Stage: networks
# Scope: global
# Consumes: required results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds; results/qc/render_compartment_abundance_figures/global/tables/source_data/v2_04_all_marker_canonical_mapping_provenance.csv
# Produces: results/spatial_validation/score_compartment_markers/global/tables/marker_class_scores_by_acquisition.csv; results/spatial_validation/score_compartment_markers/global/tables/marker_class_scores_by_animal.csv; +4 more
# Dataset behavior: global; reads all three compartments from the joint bundle.
# Notes: Marker-class abundance per acquisition and animal on the joint scale, and CON-only compartment comparisons (Figure 2d).
#
#

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/qc/compartment_abundance_utils.R")
source(repo_path("R", "preprocessing_paths.R"))
source(repo_path("R", "qc_result_paths.R"))
source(repo_path("R", "spatial_systems_paths.R"))

ANALYSIS_ID <- "score_compartment_markers"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

SCRIPT_ID <- "analysis/spatial_validation/score_compartment_markers.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

JOINT_BUNDLE <- preprocessing_find("joint_compartment_qc_matrices.rds", owner = "build_joint_protigy_input",
                                   legacy_stage = "01_preprocessing", legacy_substep = "joint_compartment_qc",
                                   scope = "global", child = "models", legacy_family = "processed")
MARKERS <- qc_find("v2_04_all_marker_canonical_mapping_provenance.csv", owner = "render_compartment_abundance_figures",
                   legacy_substep = "04e_control_compartment_abundance_publication_figures", kind = "source_data")
OUT <- list(
  acquisition = file.path(CANONICAL_PATHS$tables, "marker_class_scores_by_acquisition.csv"),
  animal      = file.path(CANONICAL_PATHS$tables, "marker_class_scores_by_animal.csv"),
  summary     = file.path(CANONICAL_PATHS$tables, "marker_class_summary.csv"),
  proteins    = file.path(CANONICAL_PATHS$tables, "marker_class_proteins.csv"),
  contrasts   = file.path(CANONICAL_PATHS$tables, "marker_class_compartment_contrasts.csv")
)
CONTRAST_GROUP <- "CON"
COMPARTMENTS <- c("neuron_neuropil", "neuron_soma", "microglia")
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Marker classes scored per acquisition on the joint scale.\n")
  dry_run_inputs(SCRIPT_ID, list(joint_bundle = JOINT_BUNDLE, marker_classes = MARKERS))
  quit(save = "no", status = 0L)
}
for (f in c(JOINT_BUNDLE, MARKERS)) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))

rec <- ca_validate_bundle_reconstruction(readRDS(JOINT_BUNDLE))
meta <- rec$metadata
cols <- rec$columns
mat <- rec$normalized_nonimputed - rec$normalization_target_median

mk <- utils::read.csv(MARKERS, stringsAsFactors = FALSE, check.names = FALSE)
mk <- mk[mk$intended_primary_eligible %in% TRUE & mk$primary_score_eligible %in% TRUE &
           mk$joint_qc_eligible %in% TRUE & !is.na(mk$ProteinGroupID) & nzchar(mk$ProteinGroupID), , drop = FALSE]
if (any(mk$n_fidelity_marker_classes > 1L, na.rm = TRUE)) stop("A marker protein belongs to more than one class.", call. = FALSE)
proteins <- unique(mk[, c("fidelity_marker_class", "ProteinGroupID", "gene_symbol", "source_name",
                          "fidelity_subpanel", "intended_dataset")])
proteins <- proteins[!duplicated(paste(proteins$fidelity_marker_class, proteins$ProteinGroupID)), , drop = FALSE]
if (!all(proteins$ProteinGroupID %in% rownames(mat))) stop("Marker proteins absent from the joint matrix.", call. = FALSE)
source_label <- c(allen_celltypes = "Allen cell types", GO_MGI = "GO (MGI)", SynGO = "SynGO")

samples <- as.character(meta[[cols$sample]])
info <- data.frame(Sample = samples, dataset = as.character(meta[[cols$dataset]]),
                   AnimalID = as.character(meta[[cols$animal]]),
                   group = ca_normalize_group(meta[[cols$group]]),
                   region = as.character(meta[[cols$region]]), layer = as.character(meta[[cols$layer]]),
                   hemisphere = ca_normalize_hemisphere(meta[[cols$hemisphere]]), stringsAsFactors = FALSE)
if (length(unique(info$AnimalID)) != 9L || !setequal(info$dataset, c("neuron_soma", "neuron_neuropil", "microglia"))) {
  stop("Expected nine animals and three compartments in the joint bundle.", call. = FALSE)
}

rows <- list()
for (cl in sort(unique(proteins$fidelity_marker_class))) {
  ids <- proteins$ProteinGroupID[proteins$fidelity_marker_class == cl]
  x <- mat[ids, samples, drop = FALSE]
  n_obs <- colSums(is.finite(x))
  rows[[cl]] <- cbind(info, marker_class = cl,
                      marker_source = paste(source_label[sort(unique(proteins$source_name[proteins$fidelity_marker_class == cl]))], collapse = " + "),
                      score_log2_vs_typical_protein = ifelse(n_obs > 0, colSums(x, na.rm = TRUE) / pmax(n_obs, 1), NA_real_),
                      n_markers_observed = as.integer(n_obs), n_markers_in_class = length(ids),
                      stringsAsFactors = FALSE)
}
acq <- do.call(rbind, rows)
rownames(acq) <- NULL
acq$score_definition <- "mean over observed class markers of joint-normalised log2 abundance minus the joint target median (0 = typical shared-core protein); not imputed"

key <- paste(acq$dataset, acq$marker_class, acq$AnimalID)
animal <- do.call(rbind, lapply(split(acq, key), function(d) {
  data.frame(dataset = d$dataset[1], marker_class = d$marker_class[1], AnimalID = d$AnimalID[1], group = d$group[1],
             n_acquisitions = sum(is.finite(d$score_log2_vs_typical_protein)),
             animal_mean_score = mean(d$score_log2_vs_typical_protein, na.rm = TRUE),
             median_markers_observed = stats::median(d$n_markers_observed), stringsAsFactors = FALSE)
}))
rownames(animal) <- NULL
summary <- do.call(rbind, lapply(split(animal, paste(animal$dataset, animal$marker_class)), function(d) {
  a <- acq[acq$dataset == d$dataset[1] & acq$marker_class == d$marker_class[1], , drop = FALSE]
  data.frame(dataset = d$dataset[1], marker_class = d$marker_class[1], n_acquisitions = sum(d$n_acquisitions),
             n_animals = nrow(d), mean_of_animal_means = mean(d$animal_mean_score),
             min_animal_mean = min(d$animal_mean_score), max_animal_mean = max(d$animal_mean_score),
             median_fraction_markers_observed = stats::median(a$n_markers_observed / a$n_markers_in_class),
             stringsAsFactors = FALSE)
}))
rownames(summary) <- NULL

# compartment comparisons on the control animals' means (see COMPARISONS above)
con <- animal[animal$group == CONTRAST_GROUP, , drop = FALSE]
contrasts <- do.call(rbind, lapply(split(con, con$marker_class), function(d) {
  d$compartment <- factor(d$dataset, levels = COMPARTMENTS)
  if (any(table(d$compartment, d$AnimalID) != 1L)) {
    stop("Every control animal needs one mean per compartment for ", d$marker_class[1], call. = FALSE)
  }
  fit <- stats::lm(animal_mean_score ~ compartment + AnimalID, data = d)
  dfe <- fit$df.residual
  se <- sqrt(2 * sum(stats::residuals(fit)^2) / dfe / length(unique(d$AnimalID)))
  m <- tapply(d$animal_mean_score, d$compartment, mean)
  pr <- utils::combn(COMPARTMENTS, 2)
  do.call(rbind, lapply(seq_len(ncol(pr)), function(i) {
    est <- m[[pr[1, i]]] - m[[pr[2, i]]]
    data.frame(marker_class = d$marker_class[1], compartment_a = pr[1, i], compartment_b = pr[2, i],
               mean_a = m[[pr[1, i]]], mean_b = m[[pr[2, i]]], estimate = est, std_error = se,
               ci_low = est - stats::qt(0.975, dfe) * se, ci_high = est + stats::qt(0.975, dfe) * se,
               df = dfe, t_value = est / se, p_value = 2 * stats::pt(-abs(est / se), dfe),
               n_animals = length(unique(d$AnimalID)), group = CONTRAST_GROUP, stringsAsFactors = FALSE)
  }))
}))
rownames(contrasts) <- NULL
contrasts$p_adjust <- stats::p.adjust(contrasts$p_value, method = "BH")
contrasts$method <- paste0("animal means of ", CONTRAST_GROUP, " animals; linear model compartment + animal (block); ",
                           "pairwise difference (a - b) with pooled error; BH over the ", nrow(contrasts), " comparisons")
detected <- do.call(rbind, lapply(c("neuron_soma", "neuron_neuropil", "microglia"), function(ds) {
  s <- info$Sample[info$dataset == ds]
  data.frame(ProteinGroupID = proteins$ProteinGroupID, dataset = ds,
             fraction_acquisitions_observed = rowMeans(is.finite(mat[proteins$ProteinGroupID, s, drop = FALSE])))
}))
detected <- stats::reshape(detected, idvar = "ProteinGroupID", timevar = "dataset", direction = "wide")
proteins_out <- merge(proteins, detected, by = "ProteinGroupID")

utils::write.csv(acq, OUT$acquisition, row.names = FALSE)
utils::write.csv(animal, OUT$animal, row.names = FALSE)
utils::write.csv(summary, OUT$summary, row.names = FALSE)
utils::write.csv(proteins_out[order(proteins_out$fidelity_marker_class, proteins_out$gene_symbol), ], OUT$proteins, row.names = FALSE)
utils::write.csv(contrasts, OUT$contrasts, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = list(joint_bundle = JOINT_BUNDLE, marker_classes = MARKERS), outputs = OUT,
                   parameters = list(scale = "joint shared-core normalisation, relative to the joint target median",
                                     imputation = "none (observed values only)",
                                     marker_rule = "04e v2_04: intended_primary_eligible & primary_score_eligible & joint_qc_eligible; one class per protein",
                                     animals = sort(unique(info$AnimalID)),
                                     comparison_animals = sort(unique(con$AnimalID)), comparison_group = CONTRAST_GROUP,
                                     comparison_model = unique(contrasts$method), inferential_p_values = TRUE),
                   notes = "Acquisitions are repeated measurements within animals; the animal is the unit of every comparison.")
message("Compartment marker scores complete")
