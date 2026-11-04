#!/usr/bin/env Rscript
#
#
# The nine proteomic animals and their acquisitions, for Extended Data (added 2026-10-30, figure review: show the
# proteomic design openly, quality control across every unit, and every animal's values of the exemplar proteins).
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/describe_proteomic_animals.R
#   Rscript analysis/spatial_validation/describe_proteomic_animals.R --dry-run
# Script: analysis/spatial_validation/describe_proteomic_animals.R
# Stage: networks
# Scope: global
# Consumes: required results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds; data/processed/01_preprocessing/protigy_input_animal_level/<dataset>/<dataset>_animal_level.gct; data/external/behaviour_outcomes/behaviour_outcomes_ov4.csv; +1 more
# Produces: results/spatial_validation/describe_proteomic_animals/global/tables/acquisition_missingness.csv; results/spatial_validation/describe_proteomic_animals/global/tables/compartment_pca_terms.csv; +4 more
# Dataset behavior: global; reads the joint acquisition bundle and the three animal-level matrices.
# Notes: Descriptive: no test. Missingness before imputation per acquisition, within-compartment PCA with the share of
#        each PC explained by anatomy, animal, group, plate and hemisphere, the design of the nine animals with their
#        behaviour, and the animal-level values of the 21 exemplar proteins of Figure 3e-g. (2026-11-04) CA2 SLM
#        missingness by group, with exact permutation tests of the group labels over the animals.
#
# MISSINGNESS: the raw quantities of the joint feature table (5,741 protein groups x 323 acquisitions). Per compartment
# the analysed proteins are those missing in at most 70% of its acquisitions (the imputation filter; gate: 5,054
# neuropil, 5,538 soma, 5,229 microglia-enriched ROI). Per acquisition: the fraction of those missing before imputation.
#
# PCA PER COMPARTMENT: the 4,242 protein groups quantified in every compartment (the matrix of Figure 2c), restricted
# to one compartment's acquisitions, each protein centred; prcomp. For PC1-PC5 the share of the PC's variance explained
# (eta squared of a one-way fit) by the spatial unit, the animal, the group, the MS plate and the hemisphere.
#
# DESIGN: per animal the group, the behaviour batch, the plate of its acquisitions, the median fraction missing, and the
# behaviour outcomes of data/external/behaviour_outcomes (ISS = CombZ, its rank among the 59 males with an ISS, 1 = lowest;
# early movement and its rank among the 53 RFID-tracked males, 1 = most active).
#
# EXEMPLAR PROTEINS: the seven leading-edge proteins of each Figure 3e-g exemplar (declared below from the Figure 3
# export), their animal-level values (hemispheres averaged, imputed; the matrices of the Figure 3 model) in the
# exemplar's unit, also centred on the mean of the three controls.
#
#

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source(repo_path("R", "preprocessing_paths.R"))
source(repo_path("R", "spatial_systems_paths.R"))
source("R/spatial/neuropil_lamina_utils.R")
source("R/spatial/stress_effect_utils.R")

ANALYSIS_ID <- "describe_proteomic_animals"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)
SCRIPT_ID <- "analysis/spatial_validation/describe_proteomic_animals.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

JOINT_BUNDLE <- preprocessing_find("joint_compartment_qc_matrices.rds", owner = "build_joint_protigy_input",
                                   legacy_stage = "01_preprocessing", legacy_substep = "joint_compartment_qc",
                                   scope = "global", child = "models", legacy_family = "processed")
OUTCOMES <- path_external("behaviour_outcomes", "behaviour_outcomes_ov4.csv")
OUTCOME_INFO <- path_external("behaviour_outcomes", "outcome_info.csv")
OUT <- list(
  missingness = file.path(CANONICAL_PATHS$tables, "acquisition_missingness.csv"),
  pca_terms   = file.path(CANONICAL_PATHS$tables, "compartment_pca_terms.csv"),
  pca_scores  = file.path(CANONICAL_PATHS$tables, "compartment_pca_scores.csv"),
  design      = file.path(CANONICAL_PATHS$tables, "proteomic_animal_design.csv"),
  exemplars   = file.path(CANONICAL_PATHS$tables, "exemplar_proteins_by_animal.csv"),
  ca2_groups  = file.path(CANONICAL_PATHS$tables, "ca2_slm_missingness_by_group.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")
DATASETS <- c(neuron_neuropil = "Neuropil", neuron_soma = "Neuronal soma", microglia = "Microglia-enriched ROI")
FILTER_N <- c(neuron_neuropil = 5054L, neuron_soma = 5538L, microglia = 5229L)
MISSING_FILTER <- 0.70
N_PC <- 5L
PCA_TERMS <- c(unit = "region_layer", animal = "animal_id", group = "experimental_group", plate = "plate", hemisphere = "ReplicateGroup")
# The Figure 3e-g exemplars and their seven top-ranked leading-edge proteins (Figure 3 export protein_zoom_values.csv,
# pRoteomics 83ace9b); unit = the animal-level matrix unit (region + layer).
EXEMPLARS <- data.frame(
  exemplar = rep(1:3, each = 7L),
  dataset = rep(c("neuron_neuropil", "neuron_soma", "microglia"), each = 7L),
  unit = rep(c("CA3sr", "CA2sp", "CA1microglia"), each = 7L),
  gene = c("Synpo", "Cnr1", "Ly6h", "Plppr4", "Eif4ebp2", "App", "Dbi",
           "Csdc2", "Lsm8", "Ddx23", "Cirbp", "Lsm3", "Dcps", "Rbm8a",
           "Cox6b1", "Iscu", "Ndufv3", "Ndufs8", "Ndufv2", "Cox5b", "Uqcrh"),
  ProteinGroupID = c("PG:neuron_neuropil:990B3F0D5E6E260C", "PG:neuron_neuropil:BA55EB377B089760", "PG:neuron_neuropil:249706BC590DAD1B",
                     "PG:neuron_neuropil:493191D0CD9A0E44", "PG:neuron_neuropil:C14B5979D8B7C678", "PG:neuron_neuropil:9ED9E1D8CDD4F38F",
                     "PG:neuron_neuropil:35B164D36402204D",
                     "PG:neuron_soma:F6775AD7570E667B", "PG:neuron_soma:864D151CF63FEFEA", "PG:neuron_soma:2BA3E3CFECCD74B1",
                     "PG:neuron_soma:3EB10AF9FA367139", "PG:neuron_soma:8DC1FF0CD676934A", "PG:neuron_soma:58ADD4989D3BFE69",
                     "PG:neuron_soma:EDA7A0D8962ACEBA",
                     "PG:microglia:54CE0440EE8DFE4A", "PG:microglia:539DC4664703DF4D", "PG:microglia:BA774D5D3E2C4456",
                     "PG:microglia:38848678964070E2", "PG:microglia:B525F1EFE43C35B7", "PG:microglia:14EA55C14EBD3B63",
                     "PG:microglia:AC490AEACC00E281"),
  stringsAsFactors = FALSE)

gct_paths <- vapply(names(DATASETS), animal_level_gct_path, "")
if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] The nine proteomic animals: missingness, compartment PCA terms, design with behaviour, exemplar proteins.\n")
  dry_run_inputs(SCRIPT_ID, c(as.list(gct_paths), list(joint_bundle = JOINT_BUNDLE, outcomes = OUTCOMES, outcome_info = OUTCOME_INFO)))
  quit(save = "no", status = 0L)
}
for (f in c(gct_paths, JOINT_BUNDLE, OUTCOMES, OUTCOME_INFO)) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))
info <- utils::read.csv(OUTCOME_INFO, stringsAsFactors = FALSE)
if (!identical(info$sha256[info$item == "output"], unname(tools::sha256sum(OUTCOMES))))
  stop("behaviour_outcomes_ov4.csv differs from outcome_info.csv; rebuild it with its script.", call. = FALSE)

jb <- readRDS(JOINT_BUNDLE); md <- jb$metadata; ft <- jb$feature_table
S <- md$.joint_sample_id
if (length(S) != 323L || !all(S %in% names(ft)) || !identical(colnames(jb$primary$matrix), S))
  stop("The joint bundle's acquisitions are not aligned.", call. = FALSE)
md$group <- unname(STRESS_GROUPS[as.character(md$experimental_group)])
md$animal <- sub("^A0*", "", md$animal_id)
if (anyNA(md$group)) stop("Unknown experimental group code.", call. = FALSE)

# ---- missingness before imputation --------------------------------------------------------------------------------
V <- as.matrix(ft[, S]); miss_rows <- list()
for (ds in names(DATASETS)) {
  s <- md$dataset == ds
  keep <- rowMeans(is.na(V[, s, drop = FALSE])) <= MISSING_FILTER
  if (sum(keep) != FILTER_N[[ds]]) stop(ds, ": the imputation filter keeps ", sum(keep), " proteins, not ", FILTER_N[[ds]], call. = FALSE)
  nm <- colSums(is.na(V[keep, s, drop = FALSE]))
  miss_rows[[ds]] <- data.frame(dataset = ds, compartment = DATASETS[[ds]], acquisition = md$Sample[s], unit = md$region_layer[s],
                                region = md$region[s], layer = md$layer[s], animal_id = md$animal[s], group = md$group[s],
                                hemisphere = md$ReplicateGroup[s], plate = md$plate[s], run_order = md$run_order[s],
                                n_proteins_analysed = sum(keep), n_missing = unname(nm), fraction_missing = unname(nm) / sum(keep),
                                n_quantified_all_features = unname(colSums(!is.na(V[, s, drop = FALSE]))), stringsAsFactors = FALSE)
}
missingness <- do.call(rbind, miss_rows)

# ---- CA2 SLM missingness by group (added 2026-11-04) ----------------------------------------------------------------
# Per animal the mean fraction missing over its CA2 SLM acquisitions (the unit of the robustness screen); per group the
# mean and range of the animal means. SUS - RES by exact permutation of the group labels over those six animals
# (two-sided on |difference|, 20 assignments, smallest attainable P 0.1) and the three groups by the between-group sum
# of squares over all 1,680 assignments of the nine animals (upper tail). A check of the statement that missingness was
# unequal between the groups.
ca2 <- missingness[missingness$unit == "CA2_slm", , drop = FALSE]
am <- stats::aggregate(fraction_missing ~ animal_id + group, data = ca2, FUN = mean)
if (nrow(am) != 9L || any(table(am$group) != 3L)) stop("CA2 SLM missingness: expected three animals per group.", call. = FALSE)
ss_between <- function(v, g) sum(tapply(v, g, function(x) length(x) * (mean(x) - mean(v))^2))
v <- am$fraction_missing; idx <- seq_along(v); obs_ss <- ss_between(v, am$group); null_ss <- numeric()
for (a in utils::combn(idx, 3L, simplify = FALSE)) for (b in utils::combn(setdiff(idx, a), 3L, simplify = FALSE)) {
  g <- character(length(v)); g[a] <- "A"; g[b] <- "B"; g[setdiff(idx, c(a, b))] <- "C"; null_ss <- c(null_ss, ss_between(v, g)) }
sr <- am[am$group %in% c("SUS", "RES"), , drop = FALSE]
d_obs <- mean(sr$fraction_missing[sr$group == "SUS"]) - mean(sr$fraction_missing[sr$group == "RES"])
null_d <- vapply(utils::combn(nrow(sr), 3L, simplify = FALSE), function(k) mean(sr$fraction_missing[k]) - mean(sr$fraction_missing[-k]), 0)
ca2_groups <- rbind(
  do.call(rbind, lapply(c("CON", "RES", "SUS"), function(g) { x <- am$fraction_missing[am$group == g]
    data.frame(row_type = "group", comparison = g, n_animals = length(x), n_acquisitions = sum(ca2$group == g), mean_fraction_missing = mean(x),
               min_animal_mean = min(x), max_animal_mean = max(x), estimate = NA_real_, statistic = NA_real_, p_exact = NA_real_,
               n_assignments = NA_integer_, test = NA_character_, stringsAsFactors = FALSE) })),
  data.frame(row_type = "test", comparison = c("SUS - RES", "CON / RES / SUS"), n_animals = c(nrow(sr), nrow(am)),
             n_acquisitions = c(sum(ca2$group %in% c("SUS", "RES")), nrow(ca2)), mean_fraction_missing = NA_real_, min_animal_mean = NA_real_,
             max_animal_mean = NA_real_, estimate = c(d_obs, NA_real_), statistic = c(d_obs, obs_ss),
             p_exact = c(mean(abs(null_d) >= abs(d_obs) - 1e-12), mean(null_ss >= obs_ss - 1e-12)), n_assignments = c(length(null_d), length(null_ss)),
             test = c("exact permutation of the group labels over the six animals, two-sided on |difference| (smallest attainable P 0.1)",
                      "exact permutation of the group labels over the nine animals, between-group sum of squares, upper tail"),
             stringsAsFactors = FALSE))
ca2_groups$unit <- "CA2_slm"
ca2_groups$measure <- paste0("animal mean of the fraction of the ", FILTER_N[["neuron_neuropil"]], " analysed neuropil protein groups missing before imputation")
if (length(null_ss) != 1680L || length(null_d) != 20L) stop("CA2 SLM missingness: wrong number of assignments.", call. = FALSE)

# ---- PCA within each compartment -----------------------------------------------------------------------------------
eta2 <- function(y, g) { g <- factor(g); if (nlevels(g) < 2L) return(NA_real_); ss <- sum((y - mean(y))^2)
  sum(tapply(y, g, function(v) length(v) * (mean(v) - mean(y))^2)) / ss }
pca_terms <- list(); pca_scores <- list()
for (ds in names(DATASETS)) {
  s <- md$dataset == ds
  X <- t(jb$primary$matrix[, s, drop = FALSE]); X <- scale(X, center = TRUE, scale = FALSE)
  pc <- stats::prcomp(X, center = FALSE, scale. = FALSE)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  for (k in seq_len(N_PC)) for (tm in names(PCA_TERMS))
    pca_terms[[length(pca_terms) + 1L]] <- data.frame(dataset = ds, compartment = DATASETS[[ds]], PC = paste0("PC", k),
      variance_explained = ve[[k]], term = tm, eta_squared = eta2(pc$x[, k], md[[PCA_TERMS[[tm]]]][s]), stringsAsFactors = FALSE)
  pca_scores[[ds]] <- data.frame(dataset = ds, compartment = DATASETS[[ds]], acquisition = md$Sample[s], unit = md$region_layer[s],
                                 animal_id = md$animal[s], group = md$group[s], plate = md$plate[s], hemisphere = md$ReplicateGroup[s],
                                 PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC1_variance = ve[[1]], PC2_variance = ve[[2]], stringsAsFactors = FALSE)
}
pca_terms <- do.call(rbind, pca_terms); pca_scores <- do.call(rbind, pca_scores)

# ---- design of the nine animals, with their behaviour --------------------------------------------------------------
beh <- utils::read.csv(OUTCOMES, stringsAsFactors = FALSE, colClasses = c(animal_id = "character", AnimalNum = "character"))
males <- beh[beh$Sex == "Male", , drop = FALSE]
males$iss_rank <- rank(males$CombZ, ties.method = "min")
tm <- males[!is.na(males$early_movement), , drop = FALSE]; tm$movement_rank <- rank(-tm$early_movement, ties.method = "min")
animals <- sort(unique(md$animal))
design <- do.call(rbind, lapply(animals, function(a) {
  s <- md$animal == a; b <- males[males$animal_id == a, , drop = FALSE]
  if (nrow(b) != 1L) stop("Proteomic animal ", a, " has no unique behaviour row among the males.", call. = FALSE)
  if (!identical(b$Group, unique(md$group[s]))) stop("Animal ", a, ": proteomic group ", unique(md$group[s]), " differs from behaviour ", b$Group, call. = FALSE)
  mm <- missingness[missingness$animal_id == a, , drop = FALSE]
  data.frame(animal_id = a, behaviour_id = b$AnimalNum, group = b$Group, batch = b$Batch,
             plate = paste(sort(unique(md$plate[s])), collapse = "+"), n_acquisitions = sum(s),
             median_fraction_missing = stats::median(mm$fraction_missing),
             ISS = b$CombZ, ISS_rank_among_males = b$iss_rank, n_males_with_ISS = nrow(males),
             early_movement = b$early_movement, early_movement_rank_among_tracked_males = tm$movement_rank[match(a, tm$animal_id)],
             n_tracked_males = nrow(tm), cort_delta = b$cort_delta, suc_pref = b$suc_pref, stringsAsFactors = FALSE)
}))
design <- design[order(match(design$group, c("CON", "RES", "SUS")), design$ISS), , drop = FALSE]

# ---- the exemplar proteins of Figure 3e-g, animal by animal ---------------------------------------------------------
ex_rows <- list()
for (ds in unique(EXEMPLARS$dataset)) {
  g <- read_animal_level_gct(ds)
  for (u in unique(EXEMPLARS$unit[EXEMPLARS$dataset == ds])) {
    ident <- canonical_protein_identity(ds, u, nrow(g$m))
    cols <- which(g$unit == u)
    if (length(cols) != 9L) stop(ds, " ", u, ": expected nine animals, found ", length(cols), call. = FALSE)
    e <- EXEMPLARS[EXEMPLARS$dataset == ds & EXEMPLARS$unit == u, , drop = FALSE]
    for (i in seq_len(nrow(e))) {
      r <- which(ident$ProteinGroupID == e$ProteinGroupID[i])
      if (length(r) != 1L || !identical(ident$official_gene_symbol[r], e$gene[i]))
        stop("Exemplar protein ", e$gene[i], " (", e$ProteinGroupID[i], ") not found once in ", ds, " ", u, call. = FALSE)
      v <- g$m[r, cols]; grp <- g$group[cols]
      ex_rows[[length(ex_rows) + 1L]] <- data.frame(exemplar = e$exemplar[i], dataset = ds, unit = u, gene = e$gene[i],
        ProteinGroupID = e$ProteinGroupID[i], animal_id = sub("^A0*", "", g$animal[cols]), group = grp, log2_abundance = unname(v),
        log2_vs_con_mean = unname(v - mean(v[grp == "CON"])), stringsAsFactors = FALSE)
    }
  }
}
exemplars <- do.call(rbind, ex_rows)
if (nrow(exemplars) != 21L * 9L) stop("Exemplar table must hold 21 proteins x 9 animals.", call. = FALSE)

# ---- write -----------------------------------------------------------------------------------------------------------
utils::write.csv(missingness, OUT$missingness, row.names = FALSE)
utils::write.csv(pca_terms, OUT$pca_terms, row.names = FALSE)
utils::write.csv(pca_scores, OUT$pca_scores, row.names = FALSE)
utils::write.csv(design, OUT$design, row.names = FALSE)
utils::write.csv(exemplars, OUT$exemplars, row.names = FALSE)
utils::write.csv(ca2_groups, OUT$ca2_groups, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = c(as.list(gct_paths), list(joint_bundle = JOINT_BUNDLE, outcomes = OUTCOMES, outcome_info = OUTCOME_INFO),
                                        as.list(stats::setNames(vapply(unique(EXEMPLARS$unit), function(u)
                                          mapped_unit_table_path(EXEMPLARS$dataset[match(u, EXEMPLARS$unit)], u), ""), unique(EXEMPLARS$unit)))),
                   outputs = OUT,
                   parameters = list(missing_filter = MISSING_FILTER, filter_n = as.list(FILTER_N), n_pc = N_PC, pca_terms = as.list(PCA_TERMS),
                                     exemplars = "Figure 3e-g, seven top-ranked leading-edge proteins each (Figure 3 export, pRoteomics 83ace9b)",
                                     provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = paste("Descriptive; the one test is the CA2 SLM missingness by group (exact permutation, added 2026-11-04).",
                                 "Behaviour outcomes from the frozen behaviour bundles (data/external/behaviour_outcomes)."))
message("Proteomic animals described: ", nrow(missingness), " acquisitions, ", nrow(design), " animals, ", nrow(exemplars) / 9L, " exemplar proteins")
