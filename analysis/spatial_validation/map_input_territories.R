#!/usr/bin/env Rscript
#
#
# Input territories of the hippocampal neuropil: are dendritic laminae that
# receive the same afferent pathway alike across subfields, and is that
# likeness carried by presynaptic proteins? (Figure 2j, 2026-10-28)
#
# TERRITORIES (pre-declared): perforant path / temporoammonic input = SLM of
# CA1 and CA2 and the DG molecular layer. Same-input pairs CA1 SLM-CA2 SLM,
# CA1 SLM-DG mo and CA2 SLM-DG mo are compared with the same-region pairs of
# those laminae, whose inputs differ (CA1 SLM-SR, CA1 SLM-SO, CA2 SLM-SR,
# CA2 SLM-SO, DG mo-po). Mossy fibres are the positive control: DG po
# (hilus)-CA3 SR against DG po-mo and CA3 SR-SO. The hilus also holds CA3c and
# the CA3 SR dissection includes stratum lucidum.
#
# STATISTIC: per animal, a lamina's profile is the mean over its retained
# acquisitions (hemispheres) of the joint-normalised observed log2 abundance.
# A protein enters when it is observed in every retained acquisition of all ten
# laminae of that animal, and it is centred over the ten laminae. For a protein
# class, the similarity of two laminae is the Spearman correlation of their
# centred profiles (Fisher z). Input organisation d = mean z of the same-input
# pairs minus mean z of the same-region pairs.
#
# PROTEIN CLASSES: SynGO cellular-component annotations (cached export),
# presynaptic = GO:0098793 and its SynGO descendants, postsynaptic = GO:0098794
# and its descendants, each without the genes of the other; human genes mapped
# to mouse by MGI homology. Genes enriched in glial or vascular cells (Allen
# cortex-hippocampus 10x trimmed means: highest glial/vascular minus highest
# neuronal log2 value > 1) and mitochondrial genes (GO:0005739 and its is_a /
# part_of descendants, MGI GAF, curated or experimental evidence only - not
# IEA or high-throughput codes, which tag many synaptic proteins) are removed
# from every class. Background = detected proteins outside SynGO, with the same
# removals. Protein groups: joint-QC eligible with an allowed gene-level claim,
# official gene symbols.
#
# NULL: per animal and class, random background sets of the class's size,
# matched on the protein's mean abundance decile; one-sided P = (k + 1) / (B + 1).
# Presynaptic vs postsynaptic: permutations of the class labels within the
# union of both classes. The excess of a class is d minus its own null median.
#
# INFERENCE: control animals only (n = 3). Per animal: the null and
# permutation P values. Across animals: one-sample t test (df = n - 1) of the
# excess of each class, and of the presynaptic minus the postsynaptic excess;
# BH over the control comparisons. The stressed animals (RES, SUS) are a replication only, tested
# the same way and never pooled with the controls.
#
# QC: DG po borders the granule cell layer. Soma-only proteins (observed in
# >= 90% of the acquisitions of every soma unit and < 50% of every other unit
# except DG po itself) index carry-over; a DG po acquisition with more soma-only detections than
# any acquisition of the other nine laminae is excluded. An animal without a
# retained acquisition of every lamina is left out.
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/map_input_territories.R
#   Rscript analysis/spatial_validation/map_input_territories.R --dry-run
# Script: analysis/spatial_validation/map_input_territories.R
# Stage: networks
# Scope: global
# Consumes: required results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds; data/external/reference_markers/syngo/raw/annotations.xlsx; data/external/reference_markers/syngo/raw/ontologies.xlsx; +4 more
# Produces: results/spatial_validation/map_input_territories/global/tables/input_territory_scores_by_animal.csv; results/spatial_validation/map_input_territories/global/tables/input_territory_tests.csv; +3 more
# Dataset behavior: global; reads the neuropil acquisitions of the joint bundle.
# Notes: Input-territory organisation of presynaptic vs postsynaptic proteins across the ten neuropil laminae, CON animals (Figure 2j); stressed animals as replication.
#
#

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/qc/compartment_abundance_utils.R")
source(repo_path("R", "preprocessing_paths.R"))
source(repo_path("R", "spatial_systems_paths.R"))

ANALYSIS_ID <- "map_input_territories"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

SCRIPT_ID <- "analysis/spatial_validation/map_input_territories.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

JOINT_BUNDLE <- preprocessing_find("joint_compartment_qc_matrices.rds", owner = "build_joint_protigy_input",
                                   legacy_stage = "01_preprocessing", legacy_substep = "joint_compartment_qc",
                                   scope = "global", child = "models", legacy_family = "processed")
REF <- list(
  syngo_annotations = path_external("reference_markers", "syngo", "raw", "annotations.xlsx"),
  syngo_ontologies  = path_external("reference_markers", "syngo", "raw", "ontologies.xlsx"),
  homology          = path_external("reference_markers", "go_mgi", "processed", "mgi_human_mouse_homology_map.csv"),
  allen_10x         = path_external("reference_markers", "allen_mouse_ctx_hip_10x", "allen_mouse_ctx_hip_10x_trimmed_means.csv"),
  mgi_gaf           = path_external("reference_markers", "go_mgi", "raw", "mgi.gaf.gz"),
  go_obo            = path_external("reference_markers", "go_mgi", "raw", "go-basic.obo")
)
OUT <- list(
  scores     = file.path(CANONICAL_PATHS$tables, "input_territory_scores_by_animal.csv"),
  tests      = file.path(CANONICAL_PATHS$tables, "input_territory_tests.csv"),
  pairs      = file.path(CANONICAL_PATHS$tables, "input_territory_pair_similarity.csv"),
  classes    = file.path(CANONICAL_PATHS$tables, "input_territory_protein_classes.csv"),
  carryover  = file.path(CANONICAL_PATHS$tables, "input_territory_carryover_qc.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")
CONTRAST_GROUP <- "CON"
N_DRAWS <- 1000L
SEED <- 20261028L
LAMINAE <- c("CA1_so", "CA1_sr", "CA1_slm", "CA2_so", "CA2_sr", "CA2_slm", "CA3_so", "CA3_sr", "DG_mo", "DG_po")
SOMA_UNITS <- c("CA1_sp", "CA2_sp", "CA3_sp", "DG_sg")
TERRITORIES <- list(
  perforant_path = list(
    label = "perforant path / temporoammonic (CA1 SLM, CA2 SLM, DG mo)",
    same_input = list(c("CA1_slm", "CA2_slm"), c("CA1_slm", "DG_mo"), c("CA2_slm", "DG_mo")),
    same_region = list(c("CA1_slm", "CA1_sr"), c("CA1_slm", "CA1_so"), c("CA2_slm", "CA2_sr"),
                       c("CA2_slm", "CA2_so"), c("DG_mo", "DG_po"))),
  mossy_fibre = list(
    label = "mossy fibre (DG po, CA3 SR) - positive control",
    same_input = list(c("DG_po", "CA3_sr")),
    same_region = list(c("DG_po", "DG_mo"), c("CA3_sr", "CA3_so")))
)

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Input-territory organisation of presynaptic and postsynaptic proteins (CON animals).\n")
  dry_run_inputs(SCRIPT_ID, c(list(joint_bundle = JOINT_BUNDLE), REF))
  quit(save = "no", status = 0L)
}
for (f in c(JOINT_BUNDLE, unlist(REF))) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))
set.seed(SEED)

# ---- data: joint bundle, observed values only --------------------------------
rec <- ca_validate_bundle_reconstruction(readRDS(JOINT_BUNDLE))
meta <- rec$metadata
cols <- rec$columns
samples <- as.character(meta[[cols$sample]])
info <- data.frame(Sample = samples, dataset = as.character(meta[[cols$dataset]]),
                   AnimalID = as.character(meta[[cols$animal]]), group = ca_normalize_group(meta[[cols$group]]),
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
symbol <- symbol[use]
pg_of <- stats::setNames(rownames(mat), symbol)
rownames(mat) <- rownames(obs) <- symbol
if (!all(LAMINAE %in% info$unit) || !all(SOMA_UNITS %in% info$unit)) stop("Expected neuropil and soma units are missing.", call. = FALSE)

# ---- QC: soma carry-over into DG po -----------------------------------------
units_all <- sort(unique(info$unit))
det <- vapply(units_all, function(k) rowMeans(obs[, info$unit == k, drop = FALSE]), numeric(nrow(obs)))
soma_only <- rownames(det)[apply(det[, SOMA_UNITS, drop = FALSE] >= 0.9, 1, all) &
                           apply(det[, setdiff(units_all, c(SOMA_UNITS, "DG_po")), drop = FALSE] < 0.5, 1, all)]
np_info <- info[info$unit %in% LAMINAE, , drop = FALSE]
np_info$soma_only_detected <- colSums(obs[soma_only, np_info$Sample, drop = FALSE])
threshold <- max(np_info$soma_only_detected[np_info$unit != "DG_po"])
np_info$carryover_threshold <- threshold
np_info$excluded <- np_info$unit == "DG_po" & np_info$soma_only_detected > threshold
carry <- np_info[, c("Sample", "AnimalID", "group", "hemisphere", "unit", "soma_only_detected",
                     "carryover_threshold", "excluded")]
carry$n_soma_only_proteins <- length(soma_only)
kept <- np_info[!np_info$excluded, , drop = FALSE]

# ---- protein classes ------------------------------------------------------------
ann <- as.data.frame(readxl::read_excel(REF$syngo_annotations))
ann <- ann[ann$go_domain == "CC", , drop = FALSE]
ont <- as.data.frame(readxl::read_excel(REF$syngo_ontologies))
syngo_desc <- function(root) {
  s <- root
  repeat {
    k <- setdiff(ont$id[ont$parent_id %in% s], s)
    if (!length(k)) break
    s <- c(s, k)
  }
  s
}
hom <- utils::read.csv(REF$homology, stringsAsFactors = FALSE)
to_mouse <- function(h) unique(hom$mouse_symbol[hom$hgnc_symbol %in% h])
pre_h <- unique(ann$hgnc_symbol[ann$go_id %in% syngo_desc("GO:0098793")])
post_h <- unique(ann$hgnc_symbol[ann$go_id %in% syngo_desc("GO:0098794")])
syngo_any <- to_mouse(unique(ann$hgnc_symbol))
pre_all <- to_mouse(setdiff(pre_h, post_h))
post_all <- to_mouse(setdiff(post_h, pre_h))

tm <- utils::read.csv(REF$allen_10x, check.names = FALSE, stringsAsFactors = FALSE)
X <- log2(as.matrix(tm[, -1, drop = FALSE]) + 1)
rownames(X) <- tm[[1]]
cell_type <- sub("^[0-9]+_", "", colnames(X))
glial_types <- c("Astro", "Oligo", "Micro-PVM", "Endo", "VLMC", "SMC-Peri")
if (!all(glial_types %in% cell_type)) stop("Allen 10x reference lacks a glial/vascular type.", call. = FALSE)
glial <- rownames(X)[apply(X[, cell_type %in% glial_types, drop = FALSE], 1, max) -
                     apply(X[, !cell_type %in% glial_types, drop = FALSE], 1, max) > 1]

obo <- readLines(REF$go_obo, warn = FALSE)
st <- grep("^\\[Term\\]", obo)
en <- c(st[-1] - 1L, length(obo))
go_par <- do.call(rbind, lapply(seq_along(st), function(i) {
  l <- obo[st[i]:en[i]]
  id <- sub("^id: ", "", l[grepl("^id: ", l)][1])
  p <- c(sub("^is_a: (GO:[0-9]+).*", "\\1", l[grepl("^is_a: ", l)]),
         sub("^relationship: part_of (GO:[0-9]+).*", "\\1", l[grepl("^relationship: part_of ", l)]))
  if (length(p)) data.frame(id = id, parent = p, stringsAsFactors = FALSE) else NULL
}))
mito_terms <- "GO:0005739"
repeat {
  k <- setdiff(go_par$id[go_par$parent %in% mito_terms], mito_terms)
  if (!length(k)) break
  mito_terms <- c(mito_terms, k)
}
gaf <- readr::read_tsv(REF$mgi_gaf, comment = "!", col_names = FALSE, col_types = readr::cols(.default = "c"),
                       progress = FALSE)
mito <- unique(gaf$X3[gaf$X9 == "C" & !grepl("NOT", gaf$X4) & gaf$X5 %in% mito_terms &
                       !gaf$X7 %in% c("IEA", "HDA", "HTP", "HEP", "HGI", "HMP")])

removed <- union(glial, mito)
classes <- data.frame(gene_symbol = symbol, ProteinGroupID = unname(pg_of[symbol]),
                      syngo_presynaptic_only = symbol %in% pre_all, syngo_postsynaptic_only = symbol %in% post_all,
                      syngo_any = symbol %in% syngo_any, glial_or_vascular = symbol %in% glial,
                      mitochondrial = symbol %in% mito, stringsAsFactors = FALSE)
classes$protein_class <- ifelse(classes$gene_symbol %in% removed, "removed",
                         ifelse(classes$syngo_presynaptic_only, "presynaptic",
                         ifelse(classes$syngo_postsynaptic_only, "postsynaptic",
                         ifelse(classes$syngo_any, "synaptic_other", "background"))))
set_of <- function(cl) classes$gene_symbol[classes$protein_class == cl]

# ---- per-animal laminar profiles ------------------------------------------------
animals <- sort(unique(kept$AnimalID))
profiles <- list()
for (a in animals) {
  ka <- kept[kept$AnimalID == a, , drop = FALSE]
  if (!all(LAMINAE %in% ka$unit)) next
  P <- vapply(LAMINAE, function(k) {
    s <- ka$Sample[ka$unit == k]
    v <- rowMeans(mat[, s, drop = FALSE])
    v[rowSums(!obs[, s, drop = FALSE]) > 0] <- NA
    v
  }, numeric(nrow(mat)))
  profiles[[a]] <- P[stats::complete.cases(P), , drop = FALSE]
}
left_out <- setdiff(animals, names(profiles))
group_of <- stats::setNames(info$group[match(names(profiles), info$AnimalID)], names(profiles))

zmat <- function(A) {
  A <- A - rowMeans(A)
  atanh(stats::cor(A, method = "spearman"))
}
pair_z <- function(C, pairs) vapply(pairs, function(p) C[p[1], p[2]], numeric(1))
dstat <- function(A, terr) {
  C <- zmat(A)
  c(same_input = mean(pair_z(C, terr$same_input)), same_region = mean(pair_z(C, terr$same_region)))
}

scores <- list()
pairs_out <- list()
for (a in names(profiles)) {
  A <- profiles[[a]]
  ip <- intersect(set_of("presynaptic"), rownames(A))
  iq <- intersect(set_of("postsynaptic"), rownames(A))
  ib <- intersect(set_of("background"), rownames(A))
  bin <- as.character(cut(rowMeans(A), stats::quantile(rowMeans(A), 0:10 / 10), include.lowest = TRUE))
  names(bin) <- rownames(A)
  draw <- function(ids) {
    out <- character(0)
    for (b in unique(bin[ids])) {
      n <- sum(bin[ids] == b)
      pool <- ib[bin[ib] == b]
      out <- c(out, sample(pool, n, replace = length(pool) < n))
    }
    out
  }
  for (tn in names(TERRITORIES)) {
    terr <- TERRITORIES[[tn]]
    dp <- dstat(A[ip, , drop = FALSE], terr)
    dq <- dstat(A[iq, , drop = FALSE], terr)
    null_p <- matrix(NA_real_, N_DRAWS, 2)
    null_q <- matrix(NA_real_, N_DRAWS, 2)
    perm <- numeric(N_DRAWS)
    both <- c(ip, iq)
    for (b in seq_len(N_DRAWS)) {
      null_p[b, ] <- dstat(A[draw(ip), , drop = FALSE], terr)
      null_q[b, ] <- dstat(A[draw(iq), , drop = FALSE], terr)
      sh <- sample(both)
      perm[b] <- diff(rev(dstat(A[sh[seq_along(ip)], , drop = FALSE], terr))) -
        diff(rev(dstat(A[sh[-seq_along(ip)], , drop = FALSE], terr)))
    }
    d_p <- dp[[1]] - dp[[2]]
    d_q <- dq[[1]] - dq[[2]]
    nd_p <- null_p[, 1] - null_p[, 2]
    nd_q <- null_q[, 1] - null_q[, 2]
    p_up <- function(x, nd) (sum(nd >= x) + 1) / (length(nd) + 1)
    p_lo <- function(x, nd) (sum(nd <= x) + 1) / (length(nd) + 1)
    base <- data.frame(AnimalID = a, group = group_of[[a]], territory = tn, territory_label = terr$label,
                       stringsAsFactors = FALSE)
    scores[[length(scores) + 1L]] <- rbind(
      cbind(base, protein_class = "presynaptic", n_proteins = length(ip), same_input_z = dp[[1]],
            same_region_z = dp[[2]], input_organisation_d = d_p, null_median_d = stats::median(nd_p),
            null_q025_d = stats::quantile(nd_p, 0.025, names = FALSE), null_q975_d = stats::quantile(nd_p, 0.975, names = FALSE),
            p_greater_than_random = p_up(d_p, nd_p), p_less_than_random = p_lo(d_p, nd_p),
            p_presynaptic_minus_postsynaptic = p_up(d_p - d_q, perm)),
      cbind(base, protein_class = "postsynaptic", n_proteins = length(iq), same_input_z = dq[[1]],
            same_region_z = dq[[2]], input_organisation_d = d_q, null_median_d = stats::median(nd_q),
            null_q025_d = stats::quantile(nd_q, 0.025, names = FALSE), null_q975_d = stats::quantile(nd_q, 0.975, names = FALSE),
            p_greater_than_random = p_up(d_q, nd_q), p_less_than_random = p_lo(d_q, nd_q),
            p_presynaptic_minus_postsynaptic = p_up(d_p - d_q, perm)),
      cbind(base, protein_class = "matched_random", n_proteins = length(ip), same_input_z = mean(null_p[, 1]),
            same_region_z = mean(null_p[, 2]), input_organisation_d = stats::median(nd_p), null_median_d = stats::median(nd_p),
            null_q025_d = stats::quantile(nd_p, 0.025, names = FALSE), null_q975_d = stats::quantile(nd_p, 0.975, names = FALSE),
            p_greater_than_random = NA_real_, p_less_than_random = NA_real_, p_presynaptic_minus_postsynaptic = NA_real_))
    for (cl in c("presynaptic", "postsynaptic")) {
      ids <- if (cl == "presynaptic") ip else iq
      C <- zmat(A[ids, , drop = FALSE])
      pl <- c(lapply(terr$same_input, function(p) c(p, "same_input")), lapply(terr$same_region, function(p) c(p, "same_region")))
      pairs_out[[length(pairs_out) + 1L]] <- data.frame(
        AnimalID = a, group = group_of[[a]], territory = tn, protein_class = cl,
        lamina_a = vapply(pl, `[`, "", 1), lamina_b = vapply(pl, `[`, "", 2), pair_class = vapply(pl, `[`, "", 3),
        fisher_z = vapply(pl, function(p) C[p[1], p[2]], numeric(1)), n_proteins = length(ids), stringsAsFactors = FALSE)
    }
  }
}
scores <- do.call(rbind, scores)
rownames(scores) <- NULL
pairs_out <- do.call(rbind, pairs_out)
pairs_out$spearman_r <- tanh(pairs_out$fisher_z)

# ---- across-animal tests: controls; stressed animals as replication -------------
one_sample <- function(x) {
  tt <- stats::t.test(x)
  data.frame(estimate = mean(x), std_error = stats::sd(x) / sqrt(length(x)), ci_low = tt$conf.int[1],
             ci_high = tt$conf.int[2], df = unname(tt$parameter), t_value = unname(tt$statistic),
             p_value = tt$p.value, n_animals = length(x), animals_positive = sum(x > 0))
}
tests <- list()
for (set_name in c(CONTRAST_GROUP, "RES+SUS")) {
  gs <- if (set_name == CONTRAST_GROUP) CONTRAST_GROUP else c("RES", "SUS")
  for (tn in names(TERRITORIES)) {
    s <- scores[scores$group %in% gs & scores$territory == tn, , drop = FALSE]
    pre <- s[s$protein_class == "presynaptic", , drop = FALSE]
    post <- s[s$protein_class == "postsynaptic", , drop = FALSE]
    post <- post[match(pre$AnimalID, post$AnimalID), , drop = FALSE]
    rows <- list(
      cbind(comparison = "presynaptic vs matched random", one_sample(pre$input_organisation_d - pre$null_median_d)),
      cbind(comparison = "postsynaptic vs matched random", one_sample(post$input_organisation_d - post$null_median_d)),
      cbind(comparison = "presynaptic vs postsynaptic",
            one_sample((pre$input_organisation_d - pre$null_median_d) - (post$input_organisation_d - post$null_median_d))))
    tests[[length(tests) + 1L]] <- cbind(animal_set = set_name, territory = tn, territory_label = TERRITORIES[[tn]]$label,
                                         do.call(rbind, rows), stringsAsFactors = FALSE)
  }
}
tests <- do.call(rbind, tests)
rownames(tests) <- NULL
tests$p_adjust <- NA_real_
for (set_name in unique(tests$animal_set)) {
  i <- tests$animal_set == set_name
  tests$p_adjust[i] <- stats::p.adjust(tests$p_value[i], method = "BH")
}
tests$role <- ifelse(tests$animal_set == CONTRAST_GROUP, "primary (control animals)",
                     "replication (stressed animals, never pooled with the controls)")
tests$method <- paste0("per animal: d = mean Fisher z (Spearman, protein-centred laminar profiles) of same-input pairs minus ",
                       "same-region pairs; excess = d minus the class's abundance-matched random-set median (", N_DRAWS,
                       " draws); across animals: one-sample t test of the excess, or of the presynaptic minus the ",
                       "postsynaptic excess; BH over the ", sum(tests$animal_set == CONTRAST_GROUP),
                       " comparisons of each animal set")
scores$excess_over_random <- scores$input_organisation_d - scores$null_median_d
scores$null_q025_excess <- scores$null_q025_d - scores$null_median_d
scores$null_q975_excess <- scores$null_q975_d - scores$null_median_d

classes$used_in_animals <- vapply(classes$gene_symbol, function(g)
  sum(vapply(profiles, function(A) g %in% rownames(A), logical(1))), integer(1))

utils::write.csv(scores, OUT$scores, row.names = FALSE)
utils::write.csv(tests, OUT$tests, row.names = FALSE)
utils::write.csv(pairs_out, OUT$pairs, row.names = FALSE)
utils::write.csv(classes[order(classes$protein_class, classes$gene_symbol), ], OUT$classes, row.names = FALSE)
utils::write.csv(carry[order(carry$unit, carry$AnimalID, carry$hemisphere), ], OUT$carryover, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = c(list(joint_bundle = JOINT_BUNDLE), REF), outputs = OUT,
                   parameters = list(scale = "joint shared-core normalisation, observed values only (no imputation)",
                                     laminae = LAMINAE, territories = lapply(TERRITORIES, function(t) t[c("same_input", "same_region")]),
                                     similarity = "Spearman correlation of protein-centred laminar profiles, Fisher z",
                                     protein_classes = "SynGO CC presynaptic (GO:0098793) / postsynaptic (GO:0098794) subtrees, mutually exclusive; glial/vascular-enriched (Allen 10x) and mitochondrial (GO:0005739, curated or experimental evidence) genes removed; joint-QC eligible protein groups with an allowed gene-level claim, official symbols",
                                     null = paste(N_DRAWS, "abundance-decile-matched background sets per animal and class"),
                                     seed = SEED, comparison_group = CONTRAST_GROUP,
                                     control_animals = sort(names(profiles)[group_of == CONTRAST_GROUP]),
                                     replication_animals = sort(names(profiles)[group_of != CONTRAST_GROUP]),
                                     animals_left_out = left_out, carryover_threshold = threshold,
                                     excluded_acquisitions = carry$Sample[carry$excluded], inferential_p_values = TRUE),
                   notes = "The animal is the unit of every across-animal comparison; hemispheres are averaged within animal.")
message("Input territories complete: ", length(profiles), " animals (", sum(group_of == CONTRAST_GROUP), " CON)")
