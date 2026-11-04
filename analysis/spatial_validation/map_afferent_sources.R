#!/usr/bin/env Rscript
#
#
# Do lamina proteomes carry the proteins of the neurons that project there?
# (Figure 2k, 2026-10-28; classes revised the same day after an independent
# review: the afferent class was dominated by a 22-cell locus coeruleus type,
# the control class held proteins hippocampal cells make, classes overlapped.
# Later on 2026-10-28: carry-over rule shared with j and l, pinned reference hash,
# animals named in the outputs)
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/map_afferent_sources.R
#   Rscript analysis/spatial_validation/map_afferent_sources.R --dry-run
# Script: analysis/spatial_validation/map_afferent_sources.R
# Stage: networks
# Scope: global
# Consumes: required results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds; data/external/allen_abc_wmb_10x/abc_wmb_source_type_mean_log2cpm.csv
# Produces: results/spatial_validation/map_afferent_sources/global/tables/afferent_source_scores_by_animal.csv; results/spatial_validation/map_afferent_sources/global/tables/afferent_source_contrasts.csv; +1 more
# Dataset behavior: global; reads the neuropil acquisitions of the joint bundle.
# Notes: Laminar location of proteins by the cell type whose mRNA they are enriched in (entorhinal, other afferent, hippocampal, glial, control), CON animals (Figure 2k); stressed animals as replication.
#
# ORIGIN CLASSES from the whole-mouse-brain 10x reference (ABC atlas, mean
# log2(CPM+1) per source type; data/external/allen_abc_wmb_10x). mRNA of other
# animals is the proxy for where a protein is made, so every class is "mRNA
# enriched in", not "made only by":
#   entorhinal     - entorhinal layer II or III neurons >= 3 and >= 1.5 above
#                    every hippocampal cell type (other afferents may share it);
#   other afferent - nucleus reuniens, supramammillary, septal, raphe or locus
#                    coeruleus neurons, same rule, not entorhinal (descriptive);
#   CA3, granule, CA1, CA2, interneurons (Sst, Pvalb, Lamp5, Vip, Sncg) -
#                    that type >= 3, >= 0.75 above every other hippocampal type
#                    and above every afferent and control type, up to 30 each;
#   astrocytes, vascular (endothelial, pericytes, VLMC), oligodendrocytes -
#                    the same (>= 3) with a margin >= 1, up to 30 each;
#   control        - striatal medium spiny, cerebellar or olfactory-bulb
#                    granule neurons >= 3 while every hippocampal and afferent
#                    type stays < 2 (tested only with >= 5 such proteins).
# The classes are mutually exclusive. Only proteins observed in every retained
# acquisition of all ten laminae of every analysed animal enter (control and
# stressed animals, so both sets are tested on the same classes; observed
# values, nothing imputed; joint-QC eligible groups with an allowed gene-level
# claim, official symbols). The reference table is pinned by its SHA-256
# (README of data/external/allen_abc_wmb_10x).
#
# ACQUISITIONS: DG po acquisitions with somatic carry-over are excluded by the
# rule of map_input_territories (R/spatial/neuropil_lamina_utils.R); an animal
# without a retained acquisition of every lamina is left out (manifest:
# animals_left_out).
#
# SCORE: per animal and lamina, the mean over a class's proteins of the
# protein-centred laminar profile (as in map_input_territories), against 1000
# abundance-decile-matched random sets of proteins outside every class; the
# excess is the score minus the random-set mean.
#
# CONTRASTS (within a region so run order and plate cancel): entorhinal CA1
# SLM - SR, CA2 SLM - SR, DG mo - po (perforant path); CA3 CA1 SR - SLM, CA2
# SR - SLM (Schaffer collaterals); granule CA3 SR - SO (mossy fibres);
# astrocytes and vascular CA1 SLM - SR (context: the SLM borders the fissure);
# with a control class, its three SLM / DG mo contrasts and the entorhinal
# minus control difference for each. Per animal: the excess difference between
# the two laminae; across the control animals a one-sample t test, BH over the
# contrasts of each animal set (8 without a control class). Stressed animals:
# replication only, tested the same way, never pooled.
#
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

ANALYSIS_ID <- "map_afferent_sources"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

SCRIPT_ID <- "analysis/spatial_validation/map_afferent_sources.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

JOINT_BUNDLE <- preprocessing_find("joint_compartment_qc_matrices.rds", owner = "build_joint_protigy_input",
                                   legacy_stage = "01_preprocessing", legacy_substep = "joint_compartment_qc",
                                   scope = "global", child = "models", legacy_family = "processed")
ABC <- path_external("allen_abc_wmb_10x", "abc_wmb_source_type_mean_log2cpm.csv")
ABC_SHA256 <- "778e6a2dc82f55cbbbdd74060e56b579600e874b4be5add9227d4fc3e4e7dd88"   # data/external/allen_abc_wmb_10x/README.md
OUT <- list(
  sets      = file.path(CANONICAL_PATHS$tables, "afferent_source_sets.csv"),
  scores    = file.path(CANONICAL_PATHS$tables, "afferent_source_scores_by_animal.csv"),
  contrasts = file.path(CANONICAL_PATHS$tables, "afferent_source_contrasts.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")
CONTRAST_GROUP <- "CON"
N_DRAWS <- 1000L
SEED <- 20261028L
LAMINAE <- c("CA1_so", "CA1_sr", "CA1_slm", "CA2_so", "CA2_sr", "CA2_slm", "CA3_so", "CA3_sr", "DG_mo", "DG_po")
SOMA_UNITS <- c("CA1_sp", "CA2_sp", "CA3_sp", "DG_sg")
SOURCES <- list(
  "EC L2" = c("EC layer II - all MEC+LEC L2 excitatory", "EC layer II - MEC stellate (ocean) cells",
              "EC layer II - MEC pyramidal island cells", "EC layer II - LEC fan cells"),
  "EC L3" = c("EC layer III - MEC", "EC layer III - LEC"),
  "Reuniens" = "Nucleus reuniens glutamatergic", "Supramammillary" = "Supramammillary glutamatergic",
  "MS/DB cholinergic" = "MS/DB cholinergic", "MS/DB GABA" = "MS/DB GABAergic", "MS/DB glutamatergic" = "MS/DB glutamatergic",
  "Raphe" = c("Dorsal raphe serotonergic", "Median raphe serotonergic"), "LC" = "Locus coeruleus noradrenergic",
  "CA1" = "CA1 pyramidal cells", "CA2" = "CA2 pyramidal cells", "CA3" = "CA3 pyramidal cells",
  "DG granule" = "Dentate granule cells", "Mossy cells" = "Hilar mossy cells",
  "Sst" = "Hippocampal Sst (O-LM) interneurons", "Pvalb" = "Hippocampal Pvalb interneurons",
  "Lamp5" = "Hippocampal Lamp5/neurogliaform interneurons", "Vip" = "Hippocampal Vip interneurons",
  "Sncg" = "Hippocampal Sncg (CCK) interneurons", "Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
  "OPC" = "OPC", "Microglia" = "Microglia", "Vascular" = c("Endothelial cells", "Pericytes", "VLMC"),
  "Striatal MSN" = c("Negative control - striatal D1 MSN", "Negative control - striatal D2 MSN"),
  "Cerebellar granule" = "Negative control - cerebellar granule cells",
  "Olfactory bulb granule" = "Negative control - olfactory bulb granule cells")
RESIDENT <- c("CA1", "CA2", "CA3", "DG granule", "Mossy cells", "Sst", "Pvalb", "Lamp5", "Vip", "Sncg",
              "Astrocytes", "Oligodendrocytes", "OPC", "Microglia", "Vascular")
AFFERENT <- c("EC L2", "EC L3", "Reuniens", "Supramammillary", "MS/DB cholinergic", "MS/DB GABA",
              "MS/DB glutamatergic", "Raphe", "LC")
NON_PROJECTING <- c("Striatal MSN", "Cerebellar granule", "Olfactory bulb granule")
ENTORHINAL <- c("EC L2", "EC L3")
OTHER_AFFERENT <- setdiff(AFFERENT, ENTORHINAL)
CLASSES <- c(entorhinal = "Entorhinal neurons", other_afferent = "Other afferent neurons", ca3 = "CA3 pyramidal cells",
             granule = "Dentate granule cells", ca1 = "CA1 pyramidal cells", ca2 = "CA2 pyramidal cells",
             interneuron = "Interneurons", astrocyte = "Astrocytes", vascular = "Vascular cells",
             oligodendrocyte = "Oligodendrocytes", control = "Neurons not projecting to the hippocampus")
CONTRASTS <- list(
  c("entorhinal", "CA1_slm", "CA1_sr", "perforant path"), c("entorhinal", "CA2_slm", "CA2_sr", "perforant path"),
  c("entorhinal", "DG_mo", "DG_po", "perforant path"), c("ca3", "CA1_sr", "CA1_slm", "Schaffer collaterals"),
  c("ca3", "CA2_sr", "CA2_slm", "Schaffer collaterals"), c("granule", "CA3_sr", "CA3_so", "mossy fibres"),
  c("astrocyte", "CA1_slm", "CA1_sr", "context (fissure)"), c("vascular", "CA1_slm", "CA1_sr", "context (fissure)"))
CONTROL_PAIRS <- list(c("CA1_slm", "CA1_sr"), c("CA2_slm", "CA2_sr"), c("DG_mo", "DG_po"))

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Laminar location of proteins by the cell type whose mRNA they are enriched in (CON animals).\n")
  dry_run_inputs(SCRIPT_ID, list(joint_bundle = JOINT_BUNDLE, abc_reference = ABC))
  quit(save = "no", status = 0L)
}
for (f in c(JOINT_BUNDLE, ABC)) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
if (!identical(unname(tools::sha256sum(ABC)), ABC_SHA256))
  stop("The ABC reference table is not the pinned version (data/external/allen_abc_wmb_10x/README.md): ", ABC, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))
set.seed(SEED)

# ---- data -----------------------------------------------------------------------------
rec <- ca_validate_bundle_reconstruction(readRDS(JOINT_BUNDLE))
meta <- rec$metadata
cols <- rec$columns
samples <- as.character(meta[[cols$sample]])
info <- data.frame(Sample = samples, AnimalID = as.character(meta[[cols$animal]]),
                   group = ca_normalize_group(meta[[cols$group]]),
                   hemisphere = ca_normalize_hemisphere(meta[[cols$hemisphere]]),
                   unit = paste(meta[[cols$region]], meta[[cols$layer]], sep = "_"), stringsAsFactors = FALSE)
info$lamina <- info$unit
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
info <- info[info$lamina %in% LAMINAE & !info$Sample %in% carry$Sample[carry$excluded], , drop = FALSE]
mat <- mat[, info$Sample, drop = FALSE]
obs <- obs[, info$Sample, drop = FALSE]

profiles <- list()
abundance <- list()
for (a in sort(unique(info$AnimalID))) {
  ia <- info[info$AnimalID == a, , drop = FALSE]
  if (!all(LAMINAE %in% ia$lamina)) next
  P <- vapply(LAMINAE, function(k) {
    s <- ia$Sample[ia$lamina == k]
    v <- rowMeans(mat[, s, drop = FALSE])
    v[rowSums(!obs[, s, drop = FALSE]) > 0] <- NA
    v
  }, numeric(nrow(mat)))
  P <- P[stats::complete.cases(P), , drop = FALSE]
  abundance[[a]] <- rowMeans(P)
  profiles[[a]] <- P - rowMeans(P)
}
group_of <- stats::setNames(info$group[match(names(profiles), info$AnimalID)], names(profiles))
left_out <- setdiff(sort(unique(info$AnimalID)), names(profiles))

# ---- origin classes ---------------------------------------------------------------------
ref <- utils::read.csv(ABC, check.names = FALSE, stringsAsFactors = FALSE)
ref <- ref[!duplicated(ref$gene_symbol), , drop = FALSE]
if (!all(unlist(SOURCES) %in% names(ref))) stop("Source types absent from the ABC table: ",
                                                paste(setdiff(unlist(SOURCES), names(ref)), collapse = ", "), call. = FALSE)
universe <- Reduce(intersect, lapply(profiles, rownames))
universe <- intersect(universe, ref$gene_symbol)
X <- as.matrix(ref[match(universe, ref$gene_symbol), unlist(SOURCES), drop = FALSE])
rownames(X) <- universe
GX <- vapply(SOURCES, function(cc) apply(X[, cc, drop = FALSE], 1, max), numeric(length(universe)))
rownames(GX) <- universe
rowmax <- function(g) apply(GX[, g, drop = FALSE], 1, max)
resident_max <- rowmax(RESIDENT)
afferent_max <- rowmax(AFFERENT)
non_projecting_max <- rowmax(NON_PROJECTING)
entorhinal_max <- rowmax(ENTORHINAL)
other_afferent_max <- rowmax(OTHER_AFFERENT)
every_type <- names(SOURCES)
made_by <- function(g, min_margin, n = 30L) {
  margin <- GX[, g] - rowmax(setdiff(RESIDENT, g))
  above_outside <- GX[, g] > rowmax(setdiff(every_type, RESIDENT))
  ok <- GX[, g] >= 3 & margin >= min_margin & above_outside
  utils::head(names(sort(margin[ok], decreasing = TRUE)), n)
}
is_entorhinal <- entorhinal_max >= 3 & entorhinal_max - resident_max >= 1.5
sets <- list(
  entorhinal = universe[is_entorhinal],
  other_afferent = universe[!is_entorhinal & other_afferent_max >= 3 & other_afferent_max - resident_max >= 1.5],
  ca3 = made_by("CA3", 0.75), granule = made_by("DG granule", 0.75), ca1 = made_by("CA1", 0.75),
  ca2 = made_by("CA2", 0.75),
  interneuron = unique(unlist(lapply(c("Sst", "Pvalb", "Lamp5", "Vip", "Sncg"), made_by, min_margin = 0.75))),
  astrocyte = made_by("Astrocytes", 1), vascular = made_by("Vascular", 1), oligodendrocyte = made_by("Oligodendrocytes", 1),
  control = universe[non_projecting_max >= 3 & resident_max < 2 & afferent_max < 2])
if (anyDuplicated(unlist(sets))) stop("Origin classes overlap: ", paste(unique(unlist(sets)[duplicated(unlist(sets))]),
                                                                       collapse = ", "), call. = FALSE)
control_tested <- length(sets$control) >= 5L
main <- setdiff(names(sets), "control")
if (any(lengths(sets[main]) < 5L)) stop("Origin class with fewer than five proteins: ",
                                         paste(main[lengths(sets[main]) < 5L], collapse = ", "), call. = FALSE)
top_source <- colnames(GX)[apply(GX, 1, which.max)]
names(top_source) <- rownames(GX)
set_rows <- do.call(rbind, lapply(names(sets), function(k) if (length(sets[[k]])) data.frame(
  origin_class = k, origin_label = CLASSES[[k]], gene_symbol = sets[[k]], top_source_type = unname(top_source[sets[[k]]]),
  top_source_log2cpm = unname(apply(GX[sets[[k]], , drop = FALSE], 1, max)),
  max_hippocampal_log2cpm = unname(resident_max[sets[[k]]]), max_entorhinal_log2cpm = unname(entorhinal_max[sets[[k]]]),
  max_other_afferent_log2cpm = unname(other_afferent_max[sets[[k]]]),
  max_non_projecting_log2cpm = unname(non_projecting_max[sets[[k]]]),
  tested = k != "control" || control_tested, stringsAsFactors = FALSE)))
if (!control_tested) sets$control <- NULL
in_any <- unique(unlist(sets))

# ---- scores against matched random sets ----------------------------------------------------
scores <- list()
for (a in names(profiles)) {
  A <- profiles[[a]][universe, , drop = FALSE]
  ab <- abundance[[a]][universe]
  bin <- as.character(cut(ab, stats::quantile(ab, 0:10 / 10), include.lowest = TRUE))
  names(bin) <- universe
  pool <- setdiff(universe, in_any)
  for (k in names(sets)) {
    ids <- sets[[k]]
    nul <- matrix(NA_real_, N_DRAWS, length(LAMINAE))
    for (b in seq_len(N_DRAWS)) {
      r <- character(0)
      for (bb in unique(bin[ids])) {
        n <- sum(bin[ids] == bb)
        pl <- pool[bin[pool] == bb]
        r <- c(r, sample(pl, n, replace = length(pl) < n))
      }
      nul[b, ] <- colMeans(A[r, , drop = FALSE])
    }
    sc <- colMeans(A[ids, , drop = FALSE])
    scores[[length(scores) + 1L]] <- data.frame(
      AnimalID = a, group = group_of[[a]], origin_class = k, origin_label = CLASSES[[k]], lamina = LAMINAE,
      n_proteins = length(ids), score = unname(sc), null_mean = colMeans(nul), null_sd = apply(nul, 2, stats::sd),
      stringsAsFactors = FALSE)
  }
}
scores <- do.call(rbind, scores)
rownames(scores) <- NULL
scores$z <- (scores$score - scores$null_mean) / scores$null_sd

# ---- contrasts ----------------------------------------------------------------------------
excess_diff <- function(s, k, a, b) {
  x <- s[s$origin_class == k, , drop = FALSE]
  xa <- x[x$lamina == a, , drop = FALSE]
  xb <- x[x$lamina == b, , drop = FALSE]
  xb <- xb[match(xa$AnimalID, xb$AnimalID), , drop = FALSE]
  stats::setNames((xa$score - xb$score) - (xa$null_mean - xb$null_mean), xa$AnimalID)
}
specs <- lapply(CONTRASTS, function(cc) list(class = cc[1], minus = NA_character_, a = cc[2], b = cc[3], expectation = cc[4]))
if (control_tested) specs <- c(specs,
  lapply(CONTROL_PAIRS, function(pp) list(class = "control", minus = NA_character_, a = pp[1], b = pp[2],
                                          expectation = "control (not expected)")),
  lapply(CONTROL_PAIRS, function(pp) list(class = "entorhinal", minus = "control", a = pp[1], b = pp[2],
                                          expectation = "entorhinal minus control")))
contrasts <- list()
for (set_name in c(CONTRAST_GROUP, "RES+SUS")) {
  gs <- if (set_name == CONTRAST_GROUP) CONTRAST_GROUP else c("RES", "SUS")
  s <- scores[scores$group %in% gs, , drop = FALSE]
  for (sp in specs) {
    v <- excess_diff(s, sp$class, sp$a, sp$b)
    if (!is.na(sp$minus)) v <- v - excess_diff(s, sp$minus, sp$a, sp$b)[names(v)]
    tt <- stats::t.test(v)
    contrasts[[length(contrasts) + 1L]] <- data.frame(
      animal_set = set_name, origin_class = sp$class, minus_class = sp$minus, origin_label = CLASSES[[sp$class]],
      lamina_a = sp$a, lamina_b = sp$b, expectation = sp$expectation, estimate = mean(v),
      std_error = stats::sd(v) / sqrt(length(v)), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
      df = unname(tt$parameter), t_value = unname(tt$statistic), p_value = tt$p.value, n_animals = length(v),
      animals_positive = sum(v > 0), n_proteins = length(sets[[sp$class]]),
      animals = paste(sort(names(v)), collapse = "; "), stringsAsFactors = FALSE)
  }
}
contrasts <- do.call(rbind, contrasts)
contrasts$p_adjust <- NA_real_
for (set_name in unique(contrasts$animal_set)) {
  i <- contrasts$animal_set == set_name
  contrasts$p_adjust[i] <- stats::p.adjust(contrasts$p_value[i], method = "BH")
}
contrasts$role <- ifelse(contrasts$animal_set == CONTRAST_GROUP, "primary (control animals)",
                         "replication (stressed animals, never pooled with the controls)")
contrasts$method <- paste0("per animal: excess difference between the laminae (class score minus the abundance-matched ",
                           "random-set mean, ", N_DRAWS, " draws); one-sample t test over animals; BH over the ",
                           length(specs), " contrasts of each animal set")
# Context in the DG (added 2026-11-04): astrocyte and vascular proteins DG mo - po, and the entorhinal excess minus each,
# the comparison the Results make for the perforant-path target outside the fissure. Computed after the BH step, so the
# pre-declared contrasts and their adjusted P are unchanged; post hoc, not corrected (p_adjust NA).
CONTEXT_SPECS <- list(list(class = "astrocyte", minus = NA_character_, a = "DG_mo", b = "DG_po", expectation = "context (DG, post hoc)"),
                      list(class = "vascular", minus = NA_character_, a = "DG_mo", b = "DG_po", expectation = "context (DG, post hoc)"),
                      list(class = "entorhinal", minus = "astrocyte", a = "DG_mo", b = "DG_po", expectation = "entorhinal minus context (DG, post hoc)"),
                      list(class = "entorhinal", minus = "vascular", a = "DG_mo", b = "DG_po", expectation = "entorhinal minus context (DG, post hoc)"))
context <- list()
for (set_name in c(CONTRAST_GROUP, "RES+SUS")) {
  gs <- if (set_name == CONTRAST_GROUP) CONTRAST_GROUP else c("RES", "SUS")
  s <- scores[scores$group %in% gs, , drop = FALSE]
  for (sp in CONTEXT_SPECS) {
    v <- excess_diff(s, sp$class, sp$a, sp$b)
    if (!is.na(sp$minus)) v <- v - excess_diff(s, sp$minus, sp$a, sp$b)[names(v)]
    tt <- stats::t.test(v)
    context[[length(context) + 1L]] <- data.frame(
      animal_set = set_name, origin_class = sp$class, minus_class = sp$minus, origin_label = CLASSES[[sp$class]],
      lamina_a = sp$a, lamina_b = sp$b, expectation = sp$expectation, estimate = mean(v),
      std_error = stats::sd(v) / sqrt(length(v)), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
      df = unname(tt$parameter), t_value = unname(tt$statistic), p_value = tt$p.value, n_animals = length(v),
      animals_positive = sum(v > 0), n_proteins = length(sets[[sp$class]]),
      animals = paste(sort(names(v)), collapse = "; "), p_adjust = NA_real_,
      role = if (set_name == CONTRAST_GROUP) "context (control animals; post hoc, not corrected)" else "context (stressed animals; post hoc, not corrected)",
      method = paste0("per animal: excess difference between the laminae (class score minus the abundance-matched random-set mean, ",
                      N_DRAWS, " draws); one-sample t test over animals; post hoc context, not corrected"), stringsAsFactors = FALSE)
  }
}
contrasts <- rbind(contrasts, do.call(rbind, context))

utils::write.csv(set_rows, OUT$sets, row.names = FALSE)
utils::write.csv(scores, OUT$scores, row.names = FALSE)
utils::write.csv(contrasts, OUT$contrasts, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = list(joint_bundle = JOINT_BUNDLE, abc_reference = ABC), outputs = OUT,
                   parameters = list(scale = "joint shared-core normalisation, observed values only (no imputation)",
                                     laminae = LAMINAE, resident_types = RESIDENT, afferent_types = AFFERENT,
                                     non_projecting_types = NON_PROJECTING, class_sizes = as.list(table(set_rows$origin_class)),
                                     control_class_tested = control_tested,
                                     universe_proteins = length(universe), null = paste(N_DRAWS, "abundance-decile-matched sets"),
                                     seed = SEED, comparison_group = CONTRAST_GROUP,
                                     control_animals = sort(names(profiles)[group_of == CONTRAST_GROUP]),
                                     replication_animals = sort(names(profiles)[group_of != CONTRAST_GROUP]),
                                     animals_left_out = left_out, carryover_threshold = carry$carryover_threshold[1],
                                     excluded_acquisitions = carry$Sample[carry$excluded], abc_reference_sha256 = ABC_SHA256,
                                     inferential_p_values = TRUE, provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "mRNA of a source cell type (other animals, ABC atlas) is a proxy for where a protein is made; the animal is the unit.")
message("Afferent sources complete: ", length(profiles), " animals, ", length(universe), " proteins")
