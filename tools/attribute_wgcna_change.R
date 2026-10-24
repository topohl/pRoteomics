#!/usr/bin/env Rscript
# Why did the WGCNA modules change after the 2026-10-19 re-normalisation?
#
# Old network: the frozen state in the pre-renormalisation archive. New network: this checkout's
# rebuilt Stage 01 state. The new expression matrix is the old one plus one constant per sample
# (the centring offset), so protein inclusion, missingness and the imputed values cannot differ
# except through that offset; this is checked first. The remaining change is decomposed by
# rebuilding the network with the build_wgcna_modules.R settings (signed bicor, maxPOutliers 0.05,
# TOM, average linkage, cutreeDynamic deepSplit 2 / minClusterSize 30, mergeCloseModules 0.25):
#   A  old data, old samples, old power            reproduces the old modules (check)
#   B  new data, old samples, old power            effect of the offset alone
#   C  new data, new samples, old power            + samples retained by the sample-tree rule
#   D  new data, new samples, new power            + rule-selected soft power = new network (check)
#   E  old data, old samples, new power            soft power alone, for comparison
# Agreement is the adjusted Rand index over shared proteins (unassigned counts as a module).
# Also written: per-sample offsets by group, the change in protein-protein bicor on the shared
# universe, and how strongly each old module eigengene tracked the offset.
#
# Writes results/wgcna/attribute_wgcna_change/<dataset>/tables/ and .../global/tables/.
# USAGE  Rscript tools/attribute_wgcna_change.R [--dataset all|<dataset>] [--threads 8]

source(if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R"))
source(repo_path("R", "wgcna_paths.R"))
suppressPackageStartupMessages({ library(WGCNA); library(dplyr); library(readr) })

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1L]] }
DATASETS <- { d <- arg("--dataset", "all"); if (identical(d, "all")) c("neuron_neuropil", "neuron_soma", "microglia") else d }
OLD_ROOT <- arg("--old-root", "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19")
enableWGCNAThreads(nThreads = as.integer(arg("--threads", "8")))

build <- function(x, power) {
  adj <- adjacency(x, power = power, type = "signed", corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.05))
  tom <- TOMsimilarity(adj, verbose = 0); rm(adj); invisible(gc())
  diss <- 1 - tom; rm(tom); invisible(gc())
  tree <- hclust(as.dist(diss), method = "average")
  mods <- cutreeDynamic(dendro = tree, distM = diss, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30, verbose = 0)
  rm(diss); invisible(gc())
  cols <- labels2colors(mods)
  merged <- mergeCloseModules(x, cols, cutHeight = 0.25, verbose = 0)$colors
  setNames(as.character(merged), colnames(x))
}
ari <- function(a, b) {
  ids <- intersect(names(a), names(b)); a <- a[ids]; b <- b[ids]
  tab <- table(a, b); n <- sum(tab); c2 <- function(v) sum(v * (v - 1) / 2)
  sij <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / (n * (n - 1) / 2)
  (sij - e) / ((sa + sb) / 2 - e)
}
n_mod <- function(x) length(setdiff(unique(x), c("grey", "unassigned")))
group_of <- function(info, ids) {
  info <- as.data.frame(info)
  key <- if (all(ids %in% rownames(info))) rownames(info) else {
    hit <- names(info)[vapply(info, function(v) all(ids %in% as.character(v)), logical(1))]
    if (length(hit)) as.character(info[[hit[1]]]) else NULL
  }
  if (is.null(key) || !"ExpGroup" %in% names(info)) return(rep(NA_character_, length(ids)))
  toupper(as.character(info$ExpGroup[match(ids, key)]))
}

overview <- list()
for (ds in DATASETS) {
  cat("== ", ds, "\n", sep = "")
  out <- wgcna_dirs("attribute_wgcna_change", ds, create = TRUE)$tables
  so <- readRDS(file.path(OLD_ROOT, "data", "processed", "06_modules_WGCNA", "01_WGCNA", ds, "wgcna_final_model_state.rds"))
  sn <- readRDS(wgcna_modules_artifact("wgcna_final_model_state.rds", ds, child = "models"))
  xo <- as.matrix(so$expression.data); storage.mode(xo) <- "double"
  xn <- as.matrix(sn$expression.data); storage.mode(xn) <- "double"
  col_o <- setNames(as.character(so$mergedColors), colnames(xo))
  col_n <- setNames(as.character(sn$mergedColors), colnames(xn))
  p_o <- as.numeric(so$softPower); p_n <- as.numeric(sn$softPower)
  P <- intersect(colnames(xo), colnames(xn)); S <- intersect(rownames(xo), rownames(xn))

  # 1. The matrices differ only by one constant per sample.
  d <- xn[S, P, drop = FALSE] - xo[S, P, drop = FALSE]
  offset <- apply(d, 1, median)
  residual <- max(abs(sweep(d, 1, offset)))
  info_n <- as.data.frame(sn$sample_info); info_o <- as.data.frame(so$sample_info)
  offsets <- tibble(dataset = ds, sample = S, group = group_of(info_n, S), offset = unname(offset))
  write_csv(offsets, file.path(out, "sample_offsets.csv"))
  inclusion <- tibble(dataset = ds,
    proteins_old = ncol(xo), proteins_new = ncol(xn), proteins_shared = length(P),
    samples_old = nrow(xo), samples_new = nrow(xn), samples_shared = length(S),
    samples_only_new = paste(setdiff(rownames(xn), rownames(xo)), collapse = "; "),
    samples_only_old = paste(setdiff(rownames(xo), rownames(xn)), collapse = "; "),
    max_residual_after_removing_offset = residual,
    old_power = p_o, new_power = p_n)
  write_csv(inclusion, file.path(out, "inclusion_and_offset_check.csv"))

  # 2. Protein-protein correlation change on the shared universe.
  ro <- bicor(xo[S, P], maxPOutliers = 0.05); rn <- bicor(xn[S, P], maxPOutliers = 0.05)
  ut <- upper.tri(ro); dr <- rn[ut] - ro[ut]
  corr_change <- tibble(dataset = ds, n_pairs = length(dr),
    median_abs_change = median(abs(dr)), p95_abs_change = unname(quantile(abs(dr), 0.95)),
    mean_change = mean(dr), r_old_vs_new = as.numeric(stats::cor(ro[ut], rn[ut])),
    mean_r_old = mean(ro[ut]), mean_r_new = mean(rn[ut]))
  rm(ro, rn, dr, ut); invisible(gc())
  write_csv(corr_change, file.path(out, "correlation_change_summary.csv"))

  # 3. Old module eigengenes versus the offset.
  me <- moduleEigengenes(xo[S, ], colors = col_o)$eigengenes
  me_assoc <- tibble(dataset = ds, old_module_eigengene = colnames(me),
    n_proteins = as.integer(table(col_o)[sub("^ME", "", colnames(me))]),
    r_with_offset = vapply(me, function(v) as.numeric(stats::cor(v, offset)), numeric(1))) |> arrange(desc(abs(r_with_offset)))
  if (!is.null(so$WGCNA_modules_long)) {
    key <- distinct(transmute(as.data.frame(so$WGCNA_modules_long), old_module_eigengene = paste0("ME", WGCNAInternalColor), ModuleID))
    me_assoc <- left_join(me_assoc, key, by = "old_module_eigengene")
  }
  write_csv(me_assoc, file.path(out, "old_module_offset_association.csv"))

  # 4. Stepwise rebuilds.
  run <- function(step, x, power, what) {
    t0 <- Sys.time(); cl <- build(x, power)
    tibble(dataset = ds, step = step, description = what, power = power, n_samples = nrow(x), n_proteins = ncol(x),
           n_modules = n_mod(cl), ari_vs_old_network = ari(cl, col_o), ari_vs_new_network = ari(cl, col_n),
           minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  }
  steps <- bind_rows(
    run("A", xo, p_o, "old data, old samples, old power (reproduces old network)"),
    run("B", xn[rownames(xo)[rownames(xo) %in% rownames(xn)], colnames(xo)[colnames(xo) %in% colnames(xn)]], p_o, "new data, old samples, old power: offset only"),
    run("C", xn[, colnames(xn) %in% colnames(xo)], p_o, "new data, new samples, old power: + sample inclusion"),
    run("D", xn, p_n, "new data, new samples, new power (reproduces new network)"),
    run("E", xo, p_n, "old data, old samples, new power: soft power only")
  )
  steps <- bind_rows(tibble(dataset = ds, step = "old", description = "frozen old network", power = p_o,
                            n_samples = nrow(xo), n_proteins = ncol(xo), n_modules = n_mod(col_o),
                            ari_vs_old_network = 1, ari_vs_new_network = ari(col_o, col_n), minutes = NA_real_), steps)
  write_csv(steps, file.path(out, "attribution_steps.csv"))
  print(as.data.frame(steps[, c("step", "power", "n_samples", "n_modules", "ari_vs_old_network", "ari_vs_new_network")]), row.names = FALSE)
  overview[[ds]] <- left_join(inclusion, corr_change, by = "dataset")
}
g <- wgcna_dirs("attribute_wgcna_change", "global", create = TRUE)$tables
write_csv(bind_rows(overview), file.path(g, "attribution_overview.csv"))
