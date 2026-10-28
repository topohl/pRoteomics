#!/usr/bin/env Rscript
#
#
# Input territories of the hippocampal neuropil: are laminae that share an
# input territory alike across subfields, and is that likeness carried by
# presynaptic proteins? (Figure 2j, 2026-10-28. 2026-10-29: the Allen trimmed
# means are read on their own log2 scale - they were log-transformed twice
# before, which loosened the glial filter -; carry-over rule shared with k and
# l; post hoc sensitivity analyses)
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/map_input_territories.R
#   Rscript analysis/spatial_validation/map_input_territories.R --dry-run
# Script: analysis/spatial_validation/map_input_territories.R
# Stage: networks
# Scope: global
# Consumes: required results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds; data/external/reference_markers/syngo/raw/annotations.xlsx; data/external/reference_markers/syngo/raw/ontologies.xlsx; +4 more
# Produces: results/spatial_validation/map_input_territories/global/tables/input_territory_scores_by_animal.csv; results/spatial_validation/map_input_territories/global/tables/input_territory_tests.csv; +5 more
# Dataset behavior: global; reads the neuropil acquisitions of the joint bundle.
# Notes: Input-territory organisation of presynaptic vs postsynaptic proteins across the ten neuropil laminae, CON animals (Figure 2j); stressed animals as replication; post hoc sensitivity analyses.
#
# TERRITORIES (pre-declared). Perforant-path territory: the SLM of CA1 and CA2
# and the DG molecular layer (DG mo = the whole molecular layer; its inner,
# middle and outer thirds are not separated). Entorhinal layer II cells project
# to DG mo and CA2 SLM, layer III (and layer II island) cells to CA1 SLM, so the
# laminae share an input pathway, not the same neurons; thalamic (reuniens),
# supramammillary and neurogliaform terminals also end there. Same-input pairs
# CA1 SLM-CA2 SLM, CA1 SLM-DG mo and CA2 SLM-DG mo are compared with the
# same-region pairs of those laminae, whose inputs differ (CA1 SLM-SR, CA1
# SLM-SO, CA2 SLM-SR, CA2 SLM-SO, DG mo-po). Mossy fibres are the positive
# control (one source, the granule cells): DG po (hilus)-CA3 SR against DG
# po-mo and CA3 SR-SO. The hilus also holds CA3c and mossy cells; the CA3 SR
# dissection includes stratum lucidum.
#
# STATISTIC (per animal and protein class): a lamina's profile is the mean over
# its retained acquisitions (hemispheres) of the joint-normalised observed log2
# abundance. A protein enters when it is observed in every retained acquisition
# of all ten laminae of that animal (nothing imputed) and is centred over the
# ten laminae, so a lamina-wide shift (e.g. synapse density) cancels and the
# similarity reflects relative composition. The similarity of two laminae is
# the Spearman correlation of their centred profiles over the class's proteins
# (Fisher z). Input organisation d = mean z of the same-input pairs minus mean z
# of the same-region pairs; d > 0: laminae that share an input are more alike
# than laminae that share a region.
#
# PROTEIN CLASSES: SynGO cellular-component annotations (cached export),
# presynaptic = GO:0098793 and its SynGO descendants, postsynaptic = GO:0098794
# and its descendants, each without the genes of the other (genes annotated to
# both, e.g. Grm2, enter neither class); human genes mapped to mouse by MGI
# homology. Genes enriched in glial or vascular cells (Allen cortex-hippocampus
# 10x trimmed means, which are log2(CPM + 1): highest glial/vascular minus
# highest neuronal value > 1) and mitochondrial genes (GO:0005739 and its is_a /
# part_of descendants, MGI GAF, curated or experimental evidence only - not
# IEA or high-throughput codes, which tag many synaptic proteins) are removed
# from every class. Background = the proteins of the animal's complete-case
# profile without a SynGO cellular-component annotation, with the same
# removals. Protein groups: joint-QC eligible with an allowed gene-level claim,
# official gene symbols.
#
# NULL AND EXCESS: per animal and class, 1,000 random background sets of the
# class's size, matched on the protein's mean abundance decile. Random proteins
# already give d away from zero (the laminae are alike in general), so the
# excess - d minus the median d of the random sets - is plotted and tested.
# Per-animal P: one-sided (k + 1) / (B + 1) against the random sets, and
# presynaptic vs postsynaptic by permuting the class labels within the union
# of both classes; both treat proteins as independent and are descriptive.
#
# INFERENCE: control animals only (n = 3). Across animals: one-sample t test
# (df = n - 1) of the excess of each class, and of the presynaptic minus the
# postsynaptic excess (a paired test); BH over the 6 control comparisons. The
# stressed animals (RES, SUS) are a replication only, tested the same way (BH
# separately) and never pooled with the controls. An animal enters only with a
# retained acquisition of every lamina: 5 of the 6 stressed animals (A0003 has
# no retained DG po; see QC and SENSITIVITY).
#
# QC: DG po borders the granule cell layer. Soma-only proteins (observed in
# >= 90% of the acquisitions of every soma unit and < 50% of every other unit
# except DG po itself) index carry-over; a DG po acquisition with more soma-only
# detections than any acquisition of the other nine laminae is excluded
# (R/spatial/neuropil_lamina_utils.R; the same rule in k and l). An animal
# without a retained acquisition of every lamina is left out.
#
# SENSITIVITY (post hoc, 2026-10-29; input_territory_sensitivity.csv; separate
# random stream, so the primary results are unchanged): (1) lamina-label
# permutation - 20,000 relabellings of the ten laminae, the same for every
# animal of a set: is the excess specific to the declared pairs? (statistic: d
# of the class Fisher z minus the mean Fisher z of 300 matched random sets); (2)
# cross-host same-input pairs only (CA1/CA2 SLM vs DG mo); (3) without the DG
# mo-po reference pair that both territories share; (4) a null matched on
# abundance and laminar SD quintiles; (5) without GABAergic-terminal proteins
# (Allen 10x: highest GABAergic cluster > 1 above every other neuronal
# cluster; GABAergic = Gad2 >= 6); (6) an animal left out for lack of DG po
# re-enters with its least-contaminated DG po acquisition (replication set).
# input_territory_protein_contrasts.csv lists the control-animal territory
# contrast of every synaptic protein (descriptive: which proteins carry it).
#

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/qc/compartment_abundance_utils.R")
source(repo_path("R", "preprocessing_paths.R"))
source(repo_path("R", "spatial_systems_paths.R"))
source("R/spatial/neuropil_lamina_utils.R")

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
  carryover  = file.path(CANONICAL_PATHS$tables, "input_territory_carryover_qc.csv"),
  sensitivity = file.path(CANONICAL_PATHS$tables, "input_territory_sensitivity.csv"),
  contrasts  = file.path(CANONICAL_PATHS$tables, "input_territory_protein_contrasts.csv")
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

# ---- QC: soma carry-over into DG po (R/spatial/neuropil_lamina_utils.R) ----------
np_info <- dg_po_carryover(info, obs, LAMINAE, SOMA_UNITS)
threshold <- np_info$carryover_threshold[1]
carry <- np_info[, c("Sample", "AnimalID", "group", "hemisphere", "unit", "soma_only_detected",
                     "carryover_threshold", "excluded", "n_soma_only_proteins")]
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
X <- as.matrix(tm[, -1, drop = FALSE])   # the trimmed means are log2(CPM + 1) already
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
                                         do.call(rbind, rows), animals = paste(sort(pre$AnimalID), collapse = "; "),
                                         stringsAsFactors = FALSE)
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

# ---- sensitivity analyses (post hoc, 2026-10-29) ----------------------------------
# A separate random stream, so the primary draws above are reproduced exactly. Every
# analysis keeps the classes, the per-animal statistic and the excess definition of
# the primary analysis and changes one thing (SENSITIVITY in the header).
set.seed(SEED + 1L)
SENS_DRAWS <- 300L
N_LAMINA_PERM <- 20000L
SENS_PAIRS <- list(
  list(territory = "perforant_path", pairs = "as declared", terr = TERRITORIES$perforant_path),
  list(territory = "perforant_path", pairs = "cross-host same-input pairs only (CA1 SLM-DG mo, CA2 SLM-DG mo)",
       terr = list(same_input = list(c("CA1_slm", "DG_mo"), c("CA2_slm", "DG_mo")),
                   same_region = TERRITORIES$perforant_path$same_region)),
  list(territory = "perforant_path", pairs = "without the DG mo-po reference pair",
       terr = list(same_input = TERRITORIES$perforant_path$same_input,
                   same_region = list(c("CA1_slm", "CA1_sr"), c("CA1_slm", "CA1_so"), c("CA2_slm", "CA2_sr"),
                                      c("CA2_slm", "CA2_so")))),
  list(territory = "mossy_fibre", pairs = "as declared", terr = TERRITORIES$mossy_fibre),
  list(territory = "mossy_fibre", pairs = "without the DG po-mo reference pair",
       terr = list(same_input = TERRITORIES$mossy_fibre$same_input, same_region = list(c("CA3_sr", "CA3_so"))))
)
sens_d <- function(Z) vapply(SENS_PAIRS, function(s) mean(pair_z(Z, s$terr$same_input)) - mean(pair_z(Z, s$terr$same_region)),
                             numeric(1))
declared <- which(vapply(SENS_PAIRS, function(s) s$pairs == "as declared", logical(1)))
permuted <- which(vapply(SENS_PAIRS, function(s) s$pairs == "as declared" || grepl("^cross-host", s$pairs), logical(1)))

# GABAergic-terminal genes: Allen 10x clusters with Gad2 >= 6 log2(CPM + 1) are GABAergic; a gene is
# GABAergic-terminal when its highest GABAergic value exceeds every other neuronal cluster by > 1.
neuronal <- !cell_type %in% glial_types
is_gaba <- neuronal & X["Gad2", ] >= 6
gaba_terminal <- rownames(X)[apply(X[, is_gaba, drop = FALSE], 1, max) -
                             apply(X[, neuronal & !is_gaba, drop = FALSE], 1, max) > 1]

# An animal left out for lack of a retained DG po acquisition re-enters with its least-contaminated one.
sens_profiles <- profiles
for (a in left_out) {
  na <- np_info[np_info$AnimalID == a, , drop = FALSE]
  po <- na[na$unit == "DG_po", , drop = FALSE]
  if (!nrow(po) || !all(LAMINAE %in% na$unit)) next
  ka <- rbind(na[!na$excluded, , drop = FALSE], po[which.min(po$soma_only_detected), , drop = FALSE])
  P <- vapply(LAMINAE, function(k) {
    s <- ka$Sample[ka$unit == k]
    v <- rowMeans(mat[, s, drop = FALSE])
    v[rowSums(!obs[, s, drop = FALSE]) > 0] <- NA
    v
  }, numeric(nrow(mat)))
  sens_profiles[[a]] <- P[stats::complete.cases(P), , drop = FALSE]
}
sens_group <- stats::setNames(info$group[match(names(sens_profiles), info$AnimalID)], names(sens_profiles))

sens_rows <- list()
zr_mean <- list()
add_row <- function(a, analysis, k, cl, d, nd) sens_rows[[length(sens_rows) + 1L]] <<- data.frame(
  AnimalID = a, group = sens_group[[a]], territory = SENS_PAIRS[[k]]$territory, pairs = SENS_PAIRS[[k]]$pairs,
  analysis = analysis, protein_class = cl, d = d, null_median_d = stats::median(nd), excess = d - stats::median(nd),
  stringsAsFactors = FALSE)
REENTER <- "animal left out for lack of DG po re-enters (its least-contaminated DG po acquisition)"
for (a in names(sens_profiles)) {
  A <- sens_profiles[[a]]
  mu <- rowMeans(A)
  lsd <- apply(A - mu, 1, stats::sd)
  ib <- intersect(set_of("background"), rownames(A))
  bin_ab <- stats::setNames(as.character(cut(mu, stats::quantile(mu, 0:10 / 10), include.lowest = TRUE)), rownames(A))
  bin_as <- stats::setNames(paste(cut(mu, stats::quantile(mu, 0:5 / 5), include.lowest = TRUE, labels = FALSE),
                                  cut(lsd, stats::quantile(lsd, 0:5 / 5), include.lowest = TRUE, labels = FALSE)),
                            rownames(A))
  sdraw <- function(ids, bin) {
    out <- character(0)
    for (b in unique(bin[ids])) {
      n <- sum(bin[ids] == b)
      pool <- ib[bin[ib] == b]
      if (!length(pool)) pool <- ib
      out <- c(out, pool[sample.int(length(pool), n, replace = length(pool) < n)])
    }
    out
  }
  for (cl in c("presynaptic", "postsynaptic")) {
    ids <- intersect(set_of(cl), rownames(A))
    Zo <- zmat(A[ids, , drop = FALSE])
    do <- sens_d(Zo)
    # declared null (abundance deciles): pair variants, and the mean random z for the permutation and the pair table
    D <- matrix(NA_real_, SENS_DRAWS, length(SENS_PAIRS))
    Zs <- 0
    for (b in seq_len(SENS_DRAWS)) {
      Z <- zmat(A[sdraw(ids, bin_ab), , drop = FALSE])
      D[b, ] <- sens_d(Z)
      Zs <- Zs + Z
    }
    zr_mean[[a]][[cl]] <- list(Zo = Zo, Zr = Zs / SENS_DRAWS)
    if (a %in% left_out) {
      for (k in declared) add_row(a, REENTER, k, cl, do[k], D[, k])
      next
    }
    for (k in setdiff(seq_along(SENS_PAIRS), declared)) add_row(a, "lamina pairs changed", k, cl, do[k], D[, k])
    # null matched on abundance and laminar SD quintiles
    D2 <- t(replicate(SENS_DRAWS, sens_d(zmat(A[sdraw(ids, bin_as), , drop = FALSE]))))
    for (k in declared) add_row(a, "null matched on abundance and laminar SD", k, cl, do[k], D2[, k])
    # without GABAergic-terminal proteins
    ids_ng <- setdiff(ids, gaba_terminal)
    do_ng <- sens_d(zmat(A[ids_ng, , drop = FALSE]))
    D3 <- t(replicate(SENS_DRAWS, sens_d(zmat(A[sdraw(ids_ng, bin_ab), , drop = FALSE]))))
    for (k in declared) add_row(a, "without GABAergic-terminal proteins", k, cl, do_ng[k], D3[, k])
  }
}
sens_rows <- do.call(rbind, sens_rows)

sens_tests <- list()
for (set_name in c(CONTRAST_GROUP, "RES+SUS")) {
  gs <- if (set_name == CONTRAST_GROUP) CONTRAST_GROUP else c("RES", "SUS")
  s <- sens_rows[sens_rows$group %in% gs, , drop = FALSE]
  keys <- unique(s[, c("analysis", "territory", "pairs")])
  for (r in seq_len(nrow(keys))) {
    ss <- s[s$analysis == keys$analysis[r] & s$territory == keys$territory[r] & s$pairs == keys$pairs[r], , drop = FALSE]
    if (keys$analysis[r] == REENTER) {
      # the set with the re-entered animal: the other animals' primary excess plus its own
      base <- scores[scores$group %in% gs & scores$territory == keys$territory[r] &
                     scores$protein_class %in% c("presynaptic", "postsynaptic"), , drop = FALSE]
      ss <- rbind(ss, data.frame(AnimalID = base$AnimalID, group = base$group, territory = base$territory,
                                 pairs = keys$pairs[r], analysis = REENTER, protein_class = base$protein_class,
                                 d = base$input_organisation_d, null_median_d = base$null_median_d,
                                 excess = base$input_organisation_d - base$null_median_d, stringsAsFactors = FALSE))
    }
    pre <- ss[ss$protein_class == "presynaptic", , drop = FALSE]
    post <- ss[ss$protein_class == "postsynaptic", , drop = FALSE]
    post <- post[match(pre$AnimalID, post$AnimalID), , drop = FALSE]
    per <- function(x) paste(paste0(pre$AnimalID, ":", sprintf("%.3f", x)), collapse = "; ")
    lab <- function(statistic, x) cbind(animal_set = set_name, territory = keys$territory[r], pairs = keys$pairs[r],
                                        analysis = keys$analysis[r], statistic = statistic, one_sample(x),
                                        per_animal = per(x), stringsAsFactors = FALSE)
    sens_tests[[length(sens_tests) + 1L]] <- rbind(
      lab("presynaptic excess", pre$excess),
      lab("postsynaptic excess", post$excess),
      lab("presynaptic minus postsynaptic excess", pre$excess - post$excess))
  }
}
sens_tests <- do.call(rbind, sens_tests)
sens_tests$p_lamina_permutation <- NA_real_

# lamina-label permutation: relabel the ten laminae (the same relabelling for every animal of a set) and
# recompute d of (class z minus mean matched random z); P in the direction of the observed value
perm_idx <- rbind(seq_along(LAMINAE), t(replicate(N_LAMINA_PERM, sample.int(length(LAMINAE)))))
side <- function(M, pl) {
  ij <- cbind(match(vapply(pl, `[`, "", 1), LAMINAE), match(vapply(pl, `[`, "", 2), LAMINAE))
  rowMeans(vapply(seq_len(nrow(ij)), function(r) M[cbind(perm_idx[, ij[r, 1]], perm_idx[, ij[r, 2]])],
                  numeric(nrow(perm_idx))))
}
p_dir <- function(v) if (v[1] >= 0) mean(v >= v[1]) else mean(v <= v[1])
perm_rows <- list()
for (set_name in c(CONTRAST_GROUP, "RES+SUS")) {
  gs <- if (set_name == CONTRAST_GROUP) CONTRAST_GROUP else c("RES", "SUS")
  an <- names(profiles)[group_of[names(profiles)] %in% gs]
  for (k in permuted) {
    terr <- SENS_PAIRS[[k]]$terr
    for (stat in c("presynaptic excess", "postsynaptic excess", "presynaptic minus postsynaptic excess")) {
      V <- vapply(an, function(a) {
        dp <- zr_mean[[a]]$presynaptic$Zo - zr_mean[[a]]$presynaptic$Zr
        dq <- zr_mean[[a]]$postsynaptic$Zo - zr_mean[[a]]$postsynaptic$Zr
        M <- switch(stat, "presynaptic excess" = dp, "postsynaptic excess" = dq, dp - dq)
        diag(M) <- 0
        side(M, terr$same_input) - side(M, terr$same_region)
      }, numeric(nrow(perm_idx)))
      g <- rowMeans(V)
      perm_rows[[length(perm_rows) + 1L]] <- data.frame(
        animal_set = set_name, territory = SENS_PAIRS[[k]]$territory, pairs = SENS_PAIRS[[k]]$pairs,
        analysis = "lamina-label permutation", statistic = stat, estimate = g[1], std_error = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_, df = NA_real_, t_value = NA_real_, p_value = NA_real_,
        n_animals = length(an), animals_positive = sum(V[1, ] > 0),
        per_animal = paste(paste0(an, ":", sprintf("%.4f", apply(V, 2, p_dir))), collapse = "; "),
        p_lamina_permutation = p_dir(g), stringsAsFactors = FALSE)
    }
  }
}
sens_tests <- rbind(sens_tests, do.call(rbind, perm_rows))
rownames(sens_tests) <- NULL
sens_tests$role <- "sensitivity (post hoc, 2026-10-29); the primary statistics are in input_territory_tests.csv"
sens_tests$ci_note <- ifelse(is.na(sens_tests$ci_low), "no interval: a permutation P for the observed statistic", "")
sens_tests$method <- ifelse(
  sens_tests$analysis == "lamina-label permutation",
  paste0("d of (class Fisher z minus the mean Fisher z of ", SENS_DRAWS, " abundance-matched random sets), mean over the ",
         "animals of the set; P = share of ", N_LAMINA_PERM, " relabellings of the ten laminae plus the observed one (the ",
         "same relabelling for every animal) at least as extreme in the direction of the observed value; per_animal: ",
         "each animal's own permutation P"),
  paste0("excess = d minus the median d of ", SENS_DRAWS, " random sets; one-sample t test over the animals of the set ",
         "(df = n - 1); per_animal: the animals' values"))

# matched random similarity of every territory pair (mean Fisher z of the declared null), for the pair table
pairs_out$matched_random_fisher_z <- vapply(seq_len(nrow(pairs_out)), function(i) {
  zr_mean[[pairs_out$AnimalID[i]]][[pairs_out$protein_class[i]]]$Zr[pairs_out$lamina_a[i], pairs_out$lamina_b[i]]
}, numeric(1))

# the proteins that carry the result: control-animal means of the centred laminar profiles (descriptive)
con_an <- names(profiles)[group_of == CONTRAST_GROUP]
common <- Reduce(intersect, lapply(profiles[con_an], rownames))
Cm <- Reduce(`+`, lapply(profiles[con_an], function(A) (A - rowMeans(A))[common, , drop = FALSE])) / length(con_an)
mean_of <- function(u) rowMeans(Cm[, u, drop = FALSE])
contrasts_out <- data.frame(
  gene_symbol = common, ProteinGroupID = unname(pg_of[common]),
  protein_class = classes$protein_class[match(common, classes$gene_symbol)],
  gabaergic_terminal = common %in% gaba_terminal, n_control_animals = length(con_an),
  perforant_territory_contrast = mean_of(c("CA1_slm", "CA2_slm", "DG_mo")) -
    mean_of(c("CA1_sr", "CA1_so", "CA2_sr", "CA2_so", "DG_po")),
  mossy_territory_contrast = mean_of(c("DG_po", "CA3_sr")) - mean_of(c("DG_mo", "CA3_so")),
  stringsAsFactors = FALSE)
contrasts_out <- cbind(contrasts_out, stats::setNames(as.data.frame(Cm[, LAMINAE, drop = FALSE]), paste0("centred_", LAMINAE)))
contrasts_out <- contrasts_out[contrasts_out$protein_class %in% c("presynaptic", "postsynaptic", "synaptic_other"), , drop = FALSE]
contrasts_out <- contrasts_out[order(contrasts_out$protein_class, -contrasts_out$perforant_territory_contrast), , drop = FALSE]
rownames(contrasts_out) <- NULL

utils::write.csv(scores, OUT$scores, row.names = FALSE)
utils::write.csv(tests, OUT$tests, row.names = FALSE)
utils::write.csv(pairs_out, OUT$pairs, row.names = FALSE)
utils::write.csv(classes[order(classes$protein_class, classes$gene_symbol), ], OUT$classes, row.names = FALSE)
utils::write.csv(carry[order(carry$unit, carry$AnimalID, carry$hemisphere), ], OUT$carryover, row.names = FALSE)
utils::write.csv(sens_tests, OUT$sensitivity, row.names = FALSE)
utils::write.csv(contrasts_out, OUT$contrasts, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = c(list(joint_bundle = JOINT_BUNDLE), REF), outputs = OUT,
                   parameters = list(scale = "joint shared-core normalisation, observed values only (no imputation)",
                                     laminae = LAMINAE, territories = lapply(TERRITORIES, function(t) t[c("same_input", "same_region")]),
                                     similarity = "Spearman correlation of protein-centred laminar profiles, Fisher z",
                                     protein_classes = "SynGO CC presynaptic (GO:0098793) / postsynaptic (GO:0098794) subtrees, mutually exclusive; glial/vascular-enriched (Allen 10x trimmed means, log2(CPM + 1): highest glial/vascular minus highest neuronal > 1) and mitochondrial (GO:0005739, curated or experimental evidence) genes removed; joint-QC eligible protein groups with an allowed gene-level claim, official symbols",
                                     null = paste(N_DRAWS, "background sets per animal and class (complete-case proteins without a SynGO CC annotation, same removals), matched on abundance decile"),
                                     seed = SEED, comparison_group = CONTRAST_GROUP,
                                     control_animals = sort(names(profiles)[group_of == CONTRAST_GROUP]),
                                     replication_animals = sort(names(profiles)[group_of != CONTRAST_GROUP]),
                                     animals_left_out = left_out, carryover_threshold = threshold,
                                     excluded_acquisitions = carry$Sample[carry$excluded], inferential_p_values = TRUE,
                                     sensitivity = list(seed = SEED + 1L, random_sets = SENS_DRAWS, lamina_relabellings = N_LAMINA_PERM,
                                                        analyses = unique(sens_tests$analysis),
                                                        gabaergic_terminal_rule = "Allen 10x: highest GABAergic cluster (Gad2 >= 6) > 1 above every other neuronal cluster",
                                                        reentered_animals = intersect(left_out, names(sens_profiles))),
                                     provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "The animal is the unit of every across-animal comparison; retained hemisphere acquisitions are averaged within animal.")
message("Input territories complete: ", length(profiles), " animals (", sum(group_of == CONTRAST_GROUP), " CON)")
