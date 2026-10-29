#!/usr/bin/env Rscript
#
#
# Stress exposure: stressed animals (RES + SUS) against controls, per unit, at
# protein, programme and declared-panel level (added 2026-10-29).
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/test_stress_exposure.R
#   Rscript analysis/spatial_validation/test_stress_exposure.R --dry-run
# Script: analysis/spatial_validation/test_stress_exposure.R
# Stage: networks
# Scope: global
# Consumes: required data/processed/01_preprocessing/protigy_input_animal_level/<dataset>/<dataset>_animal_level.gct; results/preprocessing/map_protein_identifiers/<dataset>/tables/mapped/forward/per_file/<unit>sus_<unit>res.csv; +2 more
# Produces: results/spatial_validation/test_stress_exposure/global/tables/stress_exposure_protein_tests.csv; results/spatial_validation/test_stress_exposure/global/tables/stress_exposure_unit_summary.csv; +6 more
# Dataset behavior: global; reads the animal-level matrices of neuron_neuropil, neuron_soma and microglia and the microglia acquisitions of the joint bundle.
# Notes: Stressed (RES + SUS) vs CON per unit: protein-level contrast in the Figure 3 model, GO BP GSEA as in Figure 3 with a theme summary, and a declared animal-level microglia panel with exact label permutation.
#
# WHY: Figure 3 compares SUS with RES and each with CON. The pooled exposure
# contrast - every stressed animal against the controls - is not reported, yet
# it is the axis the behaviour uses (a wider spread with no mean shift) and,
# with 6 against 3 animals, the one with the most label permutations (84).
#
# PROTEIN LEVEL: the Figure 3 model (fit_animal_level_limma.R): animal-level
# matrices (hemispheres averaged, imputed), one fit per dataset with
# ~0 + unit x group cell means, eBayes(robust = TRUE, trend = FALSE), the
# residual variance pooled over the units of a dataset. Contrast per unit:
# (RES + SUS) / 2 - CON; BH within the unit. Protein-level inference (pooled
# variance), labelled as such.
#
# PROGRAMMES: GO BP GSEA of each unit's contrast exactly as Figure 3
# (run_clusterprofiler_enrichment): protein groups collapsed to genes by the
# median of the moderated t, gseGO BP with gene sets of 10-800, nPermSimple
# 100000, the seed derived from the comparison name with base 20260824, BH
# within the comparison. Theme summary per manuscript theme
# (config/manuscript_go_theme_registry.tsv), as Figure 3b: median NES of the
# theme's terms and the number with FDR < 0.05.
#
# DECLARED PANEL (animal level): within-acquisition log2 ratios on the observed
# joint bundle (normalisation-free; every protein of a ratio observed in the
# acquisition), animal mean over acquisitions, in each microglia-enriched ROI
# and as CA2 minus the mean of the other three ROIs:
#   iba1_vs_homeostatic       Aif1 - mean(P2ry12, Csf1r)
#   complement_vs_homeostatic mean(C1qa, C1qb, C1qc) - mean(P2ry12, Csf1r)
#   homeostatic_vs_neuronal   mean(P2ry12, Csf1r) - mean(Map2, Snap25)  (microglial content)
#   choroid_vs_neuronal       Ttr - mean(Map2, Snap25)                    (composition control)
# The Aif1 ratio was seen in exploration of these animals on 2026-10-28, so its
# test here is not an independent confirmation. Stressed vs CON: exact
# permutation of the group labels over animals (84 labellings, smallest P
# 1/84), difference of animal means with a pooled two-sample 95% CI; SUS vs
# RES the same (20 labellings). BH over the panel (4 ratios x 5 locations)
# within each comparison. Checks: left-only and right-only acquisitions, and
# without A765 (the control animal with the fewest detected proteins).
#
# POST HOC LEAD (2026-10-29): Npy, lower in the three susceptible animals across
# the units where it is observed in every acquisition; selected after a probe of
# these animals, so described, not confirmed (see the block before the writes).
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
source("R/spatial/stress_effect_utils.R")
source("R/enrichment/protein_group_enrichment_utils.R")
source("R/enrichment/clusterprofiler_reproducibility.R")
source("R/enrichment/manuscript_go_theme_utils.R")
suppressPackageStartupMessages(library(limma))

ANALYSIS_ID <- "test_stress_exposure"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

SCRIPT_ID <- "analysis/spatial_validation/test_stress_exposure.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

JOINT_BUNDLE <- preprocessing_find("joint_compartment_qc_matrices.rds", owner = "build_joint_protigy_input",
                                   legacy_stage = "01_preprocessing", legacy_substep = "joint_compartment_qc",
                                   scope = "global", child = "models", legacy_family = "processed")
THEMES <- repo_path("config", "manuscript_go_theme_registry.tsv")
OUT <- list(
  proteins = file.path(CANONICAL_PATHS$tables, "stress_exposure_protein_tests.csv"),
  units    = file.path(CANONICAL_PATHS$tables, "stress_exposure_unit_summary.csv"),
  go       = file.path(CANONICAL_PATHS$tables, "stress_exposure_go_bp.csv"),
  themes   = file.path(CANONICAL_PATHS$tables, "stress_exposure_theme_summary.csv"),
  panel    = file.path(CANONICAL_PATHS$tables, "stress_exposure_panel_by_animal.csv"),
  tests    = file.path(CANONICAL_PATHS$tables, "stress_exposure_panel_tests.csv"),
  lead_animals = file.path(CANONICAL_PATHS$tables, "stress_exposure_npy_by_animal.csv"),
  lead_tests   = file.path(CANONICAL_PATHS$tables, "stress_exposure_npy_tests.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")
GSEA <- list(seed_base = 20260824L, n_perm = 100000L, min_size = 10L, max_size = 800L)
RATIOS <- list(
  iba1_vs_homeostatic       = list(num = "Aif1", den = c("P2ry12", "Csf1r")),
  complement_vs_homeostatic = list(num = c("C1qa", "C1qb", "C1qc"), den = c("P2ry12", "Csf1r")),
  homeostatic_vs_neuronal   = list(num = c("P2ry12", "Csf1r"), den = c("Map2", "Snap25")),
  choroid_vs_neuronal       = list(num = "Ttr", den = c("Map2", "Snap25"))
)
ROIS <- c("CA1microglia", "CA2microglia", "CA3microglia", "DGmicroglia")
LOW_DEPTH_CON <- "A765"

gct_paths <- vapply(STRESS_DATASETS, animal_level_gct_path, "")
if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Stressed (RES + SUS) vs CON per unit: proteins, GO BP programmes, declared microglia panel.\n")
  dry_run_inputs(SCRIPT_ID, c(as.list(gct_paths), list(joint_bundle = JOINT_BUNDLE, theme_registry = THEMES)))
  quit(save = "no", status = 0L)
}
for (f in c(gct_paths, JOINT_BUNDLE, THEMES)) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests", "work")], dir_create))
suppressPackageStartupMessages(library(org.Mm.eg.db))
# 4 workers for fgsea (this machine's limit); the worker count does not change the results
# (R/enrichment/clusterprofiler_reproducibility.R)
BiocParallel::register(BiocParallel::SnowParam(workers = 4L), default = TRUE)

# ---- protein level: the Figure 3 model with the exposure contrast -------------------------------
protein_rows <- list(); unit_rows <- list(); gsea_inputs <- list(); mapped_inputs <- character()
for (ds in STRESS_DATASETS) {
  g <- read_animal_level_gct(ds)
  cell <- factor(g$cell, levels = unique(g$cell))
  design <- model.matrix(~ 0 + cell); colnames(design) <- levels(cell)
  fit <- lmFit(g$m, design)
  for (u in unique(g$unit)) {
    cn <- paste0(u, "_", c("CON", "RES", "SUS"))
    if (!all(cn %in% colnames(design))) stop("Unit without three groups: ", ds, " ", u, call. = FALSE)
    w <- setNames(numeric(ncol(design)), colnames(design))
    cm <- cbind(stressed_vs_con = replace(w, cn, c(-1, 0.5, 0.5)),
                res_vs_con = replace(w, cn, c(-1, 1, 0)), sus_vs_con = replace(w, cn, c(-1, 0, 1)))
    eb <- eBayes(contrasts.fit(fit, cm), robust = TRUE, trend = FALSE)
    tt <- topTable(eb, coef = "stressed_vs_con", number = Inf, sort.by = "none", adjust.method = "BH")
    ident <- canonical_protein_identity(ds, u, nrow(g$m))
    mapped_inputs <- c(mapped_inputs, mapped_unit_table_path(ds, u))
    protein_rows[[length(protein_rows) + 1L]] <- data.frame(
      dataset = ds, unit = u, ProteinGroupID = ident$ProteinGroupID, official_gene_symbol = ident$official_gene_symbol,
      gene_level_claim_allowed = ident$gene_level_claim_allowed, ave_expr = tt$AveExpr,
      log2fc_stressed_minus_con = tt$logFC, t = tt$t, p_value = tt$P.Value, p_adjust = tt$adj.P.Val,
      log2fc_res_minus_con = eb$coefficients[, "res_vs_con"], log2fc_sus_minus_con = eb$coefficients[, "sus_vs_con"],
      stringsAsFactors = FALSE)
    unit_rows[[length(unit_rows) + 1L]] <- data.frame(
      dataset = ds, unit = u, n_proteins = nrow(tt), n_fdr_005 = sum(tt$adj.P.Val < 0.05), n_fdr_010 = sum(tt$adj.P.Val < 0.1),
      n_up_fdr_005 = sum(tt$adj.P.Val < 0.05 & tt$logFC > 0), n_down_fdr_005 = sum(tt$adj.P.Val < 0.05 & tt$logFC < 0),
      share_p_below_005 = mean(tt$P.Value < 0.05), prop_true_null = limma::propTrueNull(tt$P.Value),
      median_t = stats::median(tt$t), n_con = sum(g$unit == u & g$group == "CON"),
      n_stressed = sum(g$unit == u & g$group != "CON"), stringsAsFactors = FALSE)
    # canonical input of the Figure 3 GSEA, with this contrast's statistics
    df <- ident
    df$log2fc <- tt$logFC; df$t <- tt$t; df$pval <- tt$P.Value; df$padj <- tt$adj.P.Val
    gsea_inputs[[paste(ds, u, sep = "|")]] <- list(dataset = ds, unit = u, df = df)
  }
}
protein_tests <- do.call(rbind, protein_rows)
unit_summary <- do.call(rbind, unit_rows)

# ---- programmes: GO BP GSEA as Figure 3, cached per comparison -----------------------------------
# A crash of a long GSEA series loses at most the comparison in progress; a cached result is used only for the
# identical ranking, seed and settings, so it equals a recomputation.
GSEA_CACHE <- file.path(CANONICAL_PATHS$work, "gsea_cache")
dir_create(GSEA_CACHE)
go_rows <- list()
for (k in names(gsea_inputs)) {
  x <- gsea_inputs[[k]]
  comparison <- paste0(x$unit, "stressed_", x$unit, "con")
  inputs <- build_enrichment_gene_inputs(x$df, strict = TRUE)
  repro <- clusterprofiler_gsea_reproducibility_table(comparison, "BP", GSEA$seed_base, GSEA$n_perm)
  repro <- repro[repro$analysis_type == "gseGO_BP", , drop = FALSE]
  key <- list(ranked = inputs$ranked, seed = repro$gsea_seed, settings = GSEA, comparison = comparison,
              packages = c(clusterProfiler = as.character(utils::packageVersion("clusterProfiler")),
                           fgsea = as.character(utils::packageVersion("fgsea")),
                           GO.db = as.character(utils::packageVersion("GO.db"))))
  cache_file <- file.path(GSEA_CACHE, paste0(digest::digest(key, algo = "sha256"), ".rds"))
  hit <- if (file.exists(cache_file)) tryCatch(readRDS(cache_file), error = function(e) NULL) else NULL
  if (!is.null(hit) && identical(hit$key, key)) {
    go <- hit$go
  } else {
    message("GSEA ", x$dataset, " ", comparison)
    go <- as.data.frame(run_seeded_clusterprofiler_gsea(
      clusterProfiler::gseGO, gsea_seed = repro$gsea_seed, n_perm_simple = repro$n_perm_simple,
      geneList = inputs$ranked, ont = "BP", keyType = "SYMBOL", minGSSize = GSEA$min_size, maxGSSize = GSEA$max_size,
      pvalueCutoff = 1, verbose = FALSE, OrgDb = org.Mm.eg.db, pAdjustMethod = "BH"))
    tmp <- paste0(cache_file, ".tmp")
    saveRDS(list(key = key, go = go), tmp)
    file.rename(tmp, cache_file)
  }
  if (nrow(go)) {
    go_rows[[k]] <- data.frame(dataset = x$dataset, unit = x$unit, comparison = comparison, ID = go$ID,
                               Description = go$Description, setSize = go$setSize, NES = go$NES, p_value = go$pvalue,
                               p_adjust = go$p.adjust, core_enrichment = go$core_enrichment, gsea_seed = repro$gsea_seed,
                               n_perm_simple = repro$n_perm_simple, stringsAsFactors = FALSE)
  }
}
go_bp <- do.call(rbind, go_rows)
rownames(go_bp) <- NULL

registry <- read_manuscript_go_theme_registry(THEMES)
# a term counts in every theme it belongs to (terms in two themes are rare and count twice)
mapping <- map_go_terms_to_manuscript_themes(unique(go_bp[, c("ID", "Description")]), registry)$assignments
mapping <- unique(mapping[, c("GO_ID", "theme_id"), drop = FALSE])
labels <- unique(registry[, c("theme_id", "display_label", "theme_role", "display_order")])
labels <- labels[!duplicated(labels$theme_id), , drop = FALSE]
th <- merge(go_bp, mapping, by.x = "ID", by.y = "GO_ID")
theme_summary <- do.call(rbind, lapply(split(th, paste(th$theme_id, th$dataset, th$unit)), function(d) data.frame(
  theme_id = d$theme_id[1], dataset = d$dataset[1], unit = d$unit[1], n_terms_tested = length(unique(d$ID)),
  n_terms_fdr_005 = length(unique(d$ID[d$p_adjust < 0.05])),
  n_positive_fdr_005 = length(unique(d$ID[d$p_adjust < 0.05 & d$NES > 0])),
  n_negative_fdr_005 = length(unique(d$ID[d$p_adjust < 0.05 & d$NES < 0])),
  median_NES = stats::median(d$NES), stringsAsFactors = FALSE)))
theme_summary <- merge(labels, theme_summary, by = "theme_id")
theme_summary <- theme_summary[order(theme_summary$display_order, theme_summary$dataset, theme_summary$unit), , drop = FALSE]
rownames(theme_summary) <- NULL

# ---- declared panel: within-acquisition ratios in the microglia-enriched ROIs ------------------------
rec <- ca_validate_bundle_reconstruction(readRDS(JOINT_BUNDLE))
meta <- rec$metadata
cols <- rec$columns
samples <- as.character(meta[[cols$sample]])
info <- data.frame(Sample = samples, AnimalID = as.character(meta[[cols$animal]]),
                   group = ca_normalize_group(meta[[cols$group]]), hemisphere = ca_normalize_hemisphere(meta[[cols$hemisphere]]),
                   unit = paste0(meta[[cols$region]], meta[[cols$layer]]), stringsAsFactors = FALSE)
mat <- rec$normalized_nonimputed[, samples, drop = FALSE]
obs <- rec$observed[rownames(mat), samples, drop = FALSE] & is.finite(mat)
ft <- rec$feature_table[match(rownames(mat), rec$feature_table$ProteinGroupID), , drop = FALSE]
symbol <- as.character(ft$official_gene_symbol)
use <- ft$joint_qc_eligible %in% TRUE & ft$gene_level_claim_allowed %in% TRUE & !is.na(symbol) & nzchar(symbol) & !duplicated(symbol)
mat <- mat[use, , drop = FALSE]; obs <- obs[use, , drop = FALSE]
rownames(mat) <- rownames(obs) <- symbol[use]
need <- unique(unlist(RATIOS))
if (!all(need %in% rownames(mat))) stop("Panel proteins absent: ", paste(setdiff(need, rownames(mat)), collapse = ", "), call. = FALSE)
acq <- info[info$unit %in% ROIS, , drop = FALSE]
for (r in names(RATIOS)) {
  p <- RATIOS[[r]]
  ok <- colSums(!obs[c(p$num, p$den), acq$Sample, drop = FALSE]) == 0
  acq[[r]] <- ifelse(ok, colMeans(mat[p$num, acq$Sample, drop = FALSE]) - colMeans(mat[p$den, acq$Sample, drop = FALSE]), NA_real_)
}
animal_means <- function(a, hemi = NULL) {
  if (!is.null(hemi)) a <- a[a$hemisphere == hemi, , drop = FALSE]
  do.call(rbind, lapply(split(a, paste(a$AnimalID, a$unit)), function(d) {
    out <- data.frame(AnimalID = d$AnimalID[1], group = d$group[1], location = d$unit[1], stringsAsFactors = FALSE)
    for (r in names(RATIOS)) { v <- d[[r]][is.finite(d[[r]])]; out[[r]] <- if (length(v)) mean(v) else NA_real_ }
    out$n_acquisitions <- nrow(d)
    out
  }))
}
with_ca2_contrast <- function(am) {
  ca2 <- do.call(rbind, lapply(unique(am$AnimalID), function(a) {
    x <- am[am$AnimalID == a, , drop = FALSE]
    out <- data.frame(AnimalID = a, group = x$group[1], location = "CA2_minus_other_ROIs", stringsAsFactors = FALSE)
    for (r in names(RATIOS)) {
      v <- stats::setNames(x[[r]], x$location)
      out[[r]] <- if (all(is.finite(v[ROIS]))) v[["CA2microglia"]] - mean(v[setdiff(ROIS, "CA2microglia")]) else NA_real_
    }
    out$n_acquisitions <- sum(x$n_acquisitions)
    out
  }))
  rbind(am, ca2)
}
panel <- list()
for (check in c("all", "left", "right")) {
  am <- with_ca2_contrast(animal_means(acq, if (check == "all") NULL else c(left = "Left", right = "Right")[[check]]))
  am$acquisitions <- check
  panel[[check]] <- am
}
panel <- do.call(rbind, panel)
rownames(panel) <- NULL

ci_two_sample <- function(v, g) {
  ok <- is.finite(v); v <- v[ok]; g <- g[ok]
  if (sum(g) < 2L || sum(!g) < 2L) return(c(NA_real_, NA_real_))
  stats::t.test(v[g], v[!g], var.equal = TRUE)$conf.int[1:2]
}
test_rows <- list()
for (check in c("all", "left", "right", "without_A765")) {
  am <- panel[panel$acquisitions == if (check == "without_A765") "all" else check, , drop = FALSE]
  if (check == "without_A765") am <- am[am$AnimalID != LOW_DEPTH_CON, , drop = FALSE]
  for (cmp in c("stressed_vs_con", "sus_vs_res")) {
    if (cmp == "sus_vs_res" && check == "without_A765") next
    a <- if (cmp == "stressed_vs_con") am else am[am$group %in% c("RES", "SUS"), , drop = FALSE]
    for (loc in c(ROIS, "CA2_minus_other_ROIs")) for (r in names(RATIOS)) {
      d <- a[a$location == loc, , drop = FALSE]
      g <- if (cmp == "stressed_vs_con") d$group != "CON" else d$group == "SUS"
      et <- exact_label_test(d[[r]], g)
      ci <- ci_two_sample(d[[r]], g)
      ok <- is.finite(d[[r]])
      sep <- any(ok & g) && any(ok & !g) &&
        (max(d[[r]][ok & !g]) < min(d[[r]][ok & g]) || min(d[[r]][ok & !g]) > max(d[[r]][ok & g]))
      test_rows[[length(test_rows) + 1L]] <- data.frame(
        comparison = cmp, acquisitions = check, location = loc, ratio = r, estimate = et$estimate, ci_low = ci[1],
        ci_high = ci[2], n_group = et$n1, n_reference = et$n0, p_permutation = et$p, n_labellings = et$n_labellings,
        p_floor = if (et$n_labellings) 1 / et$n_labellings else NA_real_, complete_separation = sep,
        animals = paste(sort(d$AnimalID[ok]), collapse = "; "), stringsAsFactors = FALSE)
    }
  }
}
panel_tests <- do.call(rbind, test_rows)
panel_tests$p_adjust <- NA_real_
for (cmp in unique(panel_tests$comparison)) {
  i <- panel_tests$comparison == cmp & panel_tests$acquisitions == "all"
  panel_tests$p_adjust[i] <- stats::p.adjust(panel_tests$p_permutation[i], method = "BH")
}
panel_tests$role <- ifelse(panel_tests$acquisitions == "all", "declared panel (BH over 4 ratios x 5 locations per comparison)",
                           "check (left-only / right-only acquisitions; without the lowest-depth control)")
panel_tests$method <- paste0("animal mean of within-acquisition log2 ratios (every protein observed); difference of group means ",
                             "(stressed - CON or SUS - RES); exact permutation of the group labels over animals, two-sided on ",
                             "|difference|; pooled two-sample 95% CI")

# ---- post hoc lead: Npy lower in the susceptible animals (added 2026-10-29) ----------------------------
# Npy had the most SUS-RES separations in a probe of glucocorticoid and neuropeptide proteins on 2026-10-29,
# so this is a description of a lead, not a test of a hypothesis. Score per animal: the mean, over the units
# where the protein is observed in every acquisition, of its animal-level value (Figure 3 matrices; one protein
# group per symbol) centred on the unit's mean over the nine animals; the same on observed values (animal mean
# of the joint observed values). Tests: SUS vs the other six animals (84 labellings, smallest P 1/84), SUS vs
# RES, SUS vs CON and RES vs CON (20 labellings each, two-sided smallest P 0.10). Background: the same official
# score for every protein observed in every acquisition of the same units gives the lead's rank.
LEAD <- "Npy"
frames <- lapply(STRESS_DATASETS, stress_gene_frame); names(frames) <- STRESS_DATASETS
unit_ds <- unlist(lapply(STRESS_DATASETS, function(ds) stats::setNames(rep(ds, length(frames[[ds]]$units)), frames[[ds]]$units)))
full_obs <- function(sym, u) { s <- info$Sample[info$unit == u]; sym %in% rownames(obs) && length(s) > 0 && all(obs[sym, s]) }
lead_units <- names(unit_ds)[vapply(names(unit_ds), function(u) full_obs(LEAD, u), TRUE)]
official_value <- function(sym, u) {
  d <- frames[[unit_ds[[u]]]]
  j <- which(d$unit == u)
  stats::setNames(d$m[d$rows[match(sym, d$genes)], j], d$animal[j])
}
animals <- sort(unique(info$AnimalID))
grp <- stats::setNames(info$group[match(animals, info$AnimalID)], animals)
lead_rows <- do.call(rbind, lapply(lead_units, function(u) {
  v <- official_value(LEAD, u)[animals]
  s <- info[info$unit == u, , drop = FALSE]
  o <- tapply(mat[LEAD, s$Sample], s$AnimalID, mean)[animals]
  data.frame(protein = LEAD, unit = u, AnimalID = animals, group = unname(grp[animals]), official_value = unname(v),
             official_centred = unname(v - mean(v)), observed_value = unname(o), observed_centred = unname(o - mean(o)),
             stringsAsFactors = FALSE)
}))
lead_score <- function(col) tapply(lead_rows[[col]], lead_rows$AnimalID, mean)[animals]
lead_scores <- data.frame(protein = LEAD, unit = "score over units", AnimalID = animals, group = unname(grp[animals]),
                          official_value = NA_real_, official_centred = unname(lead_score("official_centred")),
                          observed_value = NA_real_, observed_centred = unname(lead_score("observed_centred")), stringsAsFactors = FALSE)
# proteome background on the official score over the same units
bg_genes <- Reduce(intersect, lapply(lead_units, function(u) {
  d <- frames[[unit_ds[[u]]]]
  g <- intersect(d$genes, rownames(obs))
  s <- info$Sample[info$unit == u]
  g[rowSums(!obs[g, s, drop = FALSE]) == 0]
}))
bg_score <- sapply(lead_units, function(u) {
  d <- frames[[unit_ds[[u]]]]; j <- which(d$unit == u)
  x <- d$m[d$rows[match(bg_genes, d$genes)], j, drop = FALSE]
  colnames(x) <- d$animal[j]
  x <- x[, animals, drop = FALSE]
  x - rowMeans(x)
}, simplify = "array")
bg_score <- apply(bg_score, c(1, 2), mean)
rownames(bg_score) <- bg_genes
is_sus <- grp[animals] == "SUS"
bg_diff <- rowMeans(bg_score[, is_sus, drop = FALSE]) - rowMeans(bg_score[, !is_sus, drop = FALSE])
lead_tests <- do.call(rbind, lapply(c("official_centred", "observed_centred"), function(col) {
  sc <- lead_scores[[col]]
  do.call(rbind, lapply(list(c("SUS", "CON+RES"), c("SUS", "RES"), c("SUS", "CON"), c("RES", "CON")), function(cmp) {
    keep <- grp %in% c(cmp[1], strsplit(cmp[2], "+", fixed = TRUE)[[1]])
    g <- grp[keep] == cmp[1]
    et <- exact_label_test(sc[keep], g)
    ci <- ci_two_sample(sc[keep], g)
    data.frame(protein = LEAD, values = col, comparison = paste(cmp[1], "vs", cmp[2]), units = paste(lead_units, collapse = "; "),
               estimate = et$estimate, ci_low = ci[1], ci_high = ci[2], p_permutation = et$p, n_labellings = et$n_labellings,
               complete_separation = max(sc[keep][g]) < min(sc[keep][!g]) || min(sc[keep][g]) > max(sc[keep][!g]),
               proteome_percentile = if (col == "official_centred" && cmp[2] == "CON+RES") 100 * mean(bg_diff <= bg_diff[[LEAD]]) else NA_real_,
               n_background_proteins = if (col == "official_centred" && cmp[2] == "CON+RES") length(bg_diff) else NA_integer_,
               role = "post hoc lead (selected after a probe of these animals); description, not confirmation",
               stringsAsFactors = FALSE)
  }))
}))
lead_by_animal <- rbind(lead_rows, lead_scores)
utils::write.csv(lead_by_animal, OUT$lead_animals, row.names = FALSE)
utils::write.csv(lead_tests, OUT$lead_tests, row.names = FALSE)

utils::write.csv(protein_tests, OUT$proteins, row.names = FALSE)
utils::write.csv(unit_summary, OUT$units, row.names = FALSE)
utils::write.csv(go_bp, OUT$go, row.names = FALSE)
utils::write.csv(theme_summary, OUT$themes, row.names = FALSE)
utils::write.csv(panel, OUT$panel, row.names = FALSE)
utils::write.csv(panel_tests, OUT$tests, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = c(as.list(gct_paths), list(joint_bundle = JOINT_BUNDLE, theme_registry = THEMES),
                                        as.list(stats::setNames(unique(mapped_inputs), basename(unique(mapped_inputs))))),
                   outputs = OUT,
                   parameters = list(protein_model = "Figure 3 model: ~0 + unit x group cell means per dataset, eBayes(robust = TRUE, trend = FALSE), contrast (RES + SUS)/2 - CON per unit, BH within unit",
                                     gsea = GSEA, ratios = RATIOS, rois = ROIS, low_depth_control = LOW_DEPTH_CON,
                                     panel_note = "the iba1_vs_homeostatic ratio was seen in exploration of these animals on 2026-10-28",
                                     provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "Protein and programme results are protein-level (pooled variance, gene ranks); the panel is animal-level.")
message("Stress exposure complete: ", nrow(unit_summary), " units, ", sum(unit_summary$n_fdr_005), " proteins at FDR < 0.05")
