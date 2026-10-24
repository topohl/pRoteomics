# Leave-one-animal-out rebuilds of an archived (original-centring) WGCNA network from its own stored input
# (normalisation review 2026-10-22). Used where the review's offset-shifted rebuilds do not reproduce the archived
# network (neuropil: ARI 0.54, because the archived build retained other samples). Same build as
# validate_sample_centring_followup.R part C (signed bicor, maxPOutliers 0.05, TOM, average linkage, cutreeDynamic
# deepSplit 2, minClusterSize 30, merge 0.25) with the archived soft power. Writes
# validation/wgcna_builds/<dataset>_archived_original_{full,loao_<animal>}.rds; the full rebuild must reproduce the
# archived network (ARI >= 0.99), otherwise the script stops.
# USAGE  Rscript tools/normalisation_review/rebuild_archived_network_loao.R neuron_neuropil [--threads 8]

suppressPackageStartupMessages({ library(WGCNA); library(dplyr) })
source(file.path("R", "paths.R"))
args <- commandArgs(TRUE); ds <- args[1]
th <- if ("--threads" %in% args) as.integer(args[match("--threads", args) + 1L]) else 8L
enableWGCNAThreads(nThreads = th)
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
cache <- path_results("reviewer_audit", "normalisation_choice_2026-10-21", "validation", "wgcna_builds")
so <- readRDS(file.path(ARCHIVE, "data", "processed", "06_modules_WGCNA", "01_WGCNA", ds, "wgcna_final_model_state.rds"))
x <- as.matrix(so$expression.data); storage.mode(x) <- "double"
pw <- so$softPower; col_old <- setNames(as.character(so$mergedColors), colnames(x))
# the animal is part of the acquisition name (..._Tobias_<animal>_<hemisphere>_<unit>_...)
animal <- sub("_.*", "", sub("^.*Tobias_", "", rownames(x)))
stopifnot(!anyNA(animal), all(grepl("^A[0-9]+$", animal)))
cat(ds, ": ", nrow(x), " samples, ", ncol(x), " proteins, power ", pw, ", animals ", paste(sort(unique(animal)), collapse = ","), "\n", sep = "")
build <- function(x, power) {
  adj <- adjacency(x, power = power, type = "signed", corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.05))
  tom <- TOMsimilarity(adj, verbose = 0); rm(adj); invisible(gc())
  diss <- 1 - tom; rm(tom); invisible(gc())
  tree <- hclust(as.dist(diss), method = "average")
  mods <- cutreeDynamic(dendro = tree, distM = diss, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30, verbose = 0)
  rm(diss); invisible(gc())
  merged <- mergeCloseModules(x, labels2colors(mods), cutHeight = 0.25, verbose = 0)$colors
  setNames(as.character(merged), colnames(x))
}
ari <- function(a, b) {
  ids <- intersect(names(a), names(b)); a <- a[ids]; b <- b[ids]
  tab <- table(a, b); n <- sum(tab); c2 <- function(v) sum(v * (v - 1) / 2)
  sij <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / (n * (n - 1) / 2)
  (sij - e) / ((sa + sb) / 2 - e)
}
key <- function(k) file.path(cache, sprintf("%s_archived_original_%s.rds", ds, k))
full <- if (file.exists(key("full"))) readRDS(key("full")) else { f <- build(x, pw); saveRDS(f, key("full")); f }
a <- ari(full, col_old); cat("full rebuild vs archived network: ARI", round(a, 4), "\n")
if (a < 0.99) stop("the rebuild does not reproduce the archived network (ARI ", round(a, 3), "); stability would not refer to it.")
for (an in sort(unique(animal))) {
  if (file.exists(key(paste0("loao_", an)))) next
  cl <- build(x[animal != an, , drop = FALSE], pw); saveRDS(cl, key(paste0("loao_", an)))
  cat("  without", an, ": ARI vs full", round(ari(cl, full), 3), "\n")
}
cat("done\n")
