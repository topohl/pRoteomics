#!/usr/bin/env Rscript
#
#
# Do published hippocampal stress datasets show the same SUS / RES / CON
# differences as our spatial proteome? (Extended Data, added 2026-10-29)
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/compare_external_stress_studies.R
#   Rscript analysis/spatial_validation/compare_external_stress_studies.R --dry-run
# Script: analysis/spatial_validation/compare_external_stress_studies.R
# Stage: networks
# Scope: global
# Consumes: required data/processed/01_preprocessing/protigy_input_animal_level/<dataset>/<dataset>_animal_level.gct; data/external/stress_signatures/signature_info.csv; data/external/stress_signatures/signatures/<id>.csv; config/manuscript_go_theme_registry.tsv
# Produces: results/spatial_validation/compare_external_stress_studies/global/tables/external_concordance_signatures.csv; results/spatial_validation/compare_external_stress_studies/global/tables/external_concordance_omnibus.csv; +2 more
# Dataset behavior: global; reads the animal-level matrices of neuron_neuropil, neuron_soma and microglia.
# Notes: Rank concordance of our per-unit SUS-CON, RES-CON, SUS-RES and stressed-CON t with published signatures (gene, GO-BP, Figure 3 theme terms), exact null over all 1,680 assignments of the nine animals to three groups.
#
# SIGNATURES: data/external/stress_signatures (README there): one statistic per
# mouse gene and contrast for each published dataset, rebuilt from the deposited
# data, with the units each should be compared with (match_units). The set
# 'primary' enters the tests and the omnibus; 'exploratory' is reported only.
#
# OUR SIDE: the Figure 3 model (fit_animal_level_limma.R inputs and design:
# animal-level values, ~0 + unit x group, eBayes(robust = TRUE)), gene-level t
# (one protein group per symbol) for SUS - CON, RES - CON, SUS - RES and
# stressed (RES + SUS)/2 - CON in every unit.
#
# STATISTIC: for a signature and contrast, the Spearman correlation of our t with
# the published statistic over shared genes, in each matched unit, averaged over
# the matched units. Three levels: genes; GO biological-process terms (mean
# statistic of the term's genes, 10-800 shared genes); the Figure 3 theme terms
# only (descriptive). Omnibus per contrast and level: the mean over the primary
# signatures that carry the contrast.
#
# NULL: our t recomputed for every assignment of the animals to the three groups with
# the observed group sizes (1,680 for nine animals in groups of three, 560 for the
# 3/3/2 design); one-sided p = share of assignments with a correlation at least as
# large. Contrast-specific p: only the assignments that keep the third group fixed
# (20 with 3/3/3; 10 or 20 with 3/3/2), and all splits for stressed - CON. BH over
# signature x contrast within each level (primary set). The p columns keep their
# names (p_all_1680) for compatibility; n_relabellings gives the actual count.
#
# AGREEMENT: Spearman correlation between every pair of published signatures,
# per contrast (shared genes) - how much the published studies agree with each other.
#
#
#

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))
source("R/spatial/neuropil_lamina_utils.R")
source("R/spatial/stress_effect_utils.R")
source("R/enrichment/manuscript_go_theme_utils.R")
suppressPackageStartupMessages({ library(limma); library(Matrix); library(parallel) })

ANALYSIS_ID <- "compare_external_stress_studies"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)
SCRIPT_ID <- "analysis/spatial_validation/compare_external_stress_studies.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

SIG_DIR <- path_external("stress_signatures")
SIG_INFO <- file.path(SIG_DIR, "signature_info.csv")
THEME_REGISTRY <- repo_path("config", "manuscript_go_theme_registry.tsv")
OUT <- list(
  signatures = file.path(CANONICAL_PATHS$tables, "external_concordance_signatures.csv"),
  omnibus    = file.path(CANONICAL_PATHS$tables, "external_concordance_omnibus.csv"),
  units      = file.path(CANONICAL_PATHS$tables, "external_concordance_units.csv"),
  agreement  = file.path(CANONICAL_PATHS$tables, "external_signature_agreement.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")
N_WORKERS <- 4L          # host limit for parallel R workers
GO_MIN <- 10L; GO_MAX <- 800L
THEMES <- c("rna_processing_splicing_rnp", "ribosome_translation", "chromatin_organization", "mitochondrial_respiration_oxphos",
            "synaptic_signaling_vesicle", "neuron_projection_development", "autophagy_lysosome_endosome")
ALL_UNITS <- c("CA1so", "CA1sr", "CA1slm", "CA2so", "CA2sr", "CA2slm", "CA3so", "CA3sr", "DGmo", "DGpo",
               "CA1sp", "CA2sp", "CA3sp", "DGsg", "CA1microglia", "CA2microglia", "CA3microglia", "DGmicroglia")

gct_paths <- vapply(STRESS_DATASETS, animal_level_gct_path, "")
if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Concordance with published stress signatures (all relabellings with the observed group sizes).\n")
  dry_run_inputs(SCRIPT_ID, c(as.list(gct_paths), list(signature_info = SIG_INFO, theme_registry = THEME_REGISTRY)))
  quit(save = "no", status = 0L)
}
for (f in c(gct_paths, SIG_INFO, THEME_REGISTRY)) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))

# ---- signatures (checksums pinned in signature_info.csv) ----
info <- utils::read.csv(SIG_INFO, stringsAsFactors = FALSE, check.names = FALSE)
sig_files <- file.path(SIG_DIR, "signatures", paste0(info$id, ".csv"))
for (i in seq_len(nrow(info))) {
  if (!file.exists(sig_files[i])) stop("missing_required_input: ", sig_files[i], call. = FALSE)
  h <- digest::digest(file = sig_files[i], algo = "sha256")
  if (!identical(h, info$sha256[i])) stop("Signature checksum differs from signature_info.csv: ", info$id[i], call. = FALSE)
}
sigs <- lapply(seq_len(nrow(info)), function(i) {
  x <- utils::read.csv(sig_files[i], stringsAsFactors = FALSE, check.names = FALSE)
  st <- as.matrix(x[, intersect(STRESS_CONTRASTS, names(x)), drop = FALSE]); rownames(st) <- x$gene
  mu <- strsplit(info$match_units[i], " ", fixed = TRUE)[[1]]
  list(id = info$id[i], set = info$set[i], stat = st, match_units = if (identical(mu, "all")) ALL_UNITS else mu)
})
names(sigs) <- info$id

# ---- GO-BP term membership and the Figure 3 theme terms ----
go_bp_sets <- function(genes) {
  a <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = genes, keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL")))
  a <- unique(a[a$ONTOLOGYALL %in% "BP" & !is.na(a$GOALL), c("SYMBOL", "GOALL")])
  split(a$SYMBOL, a$GOALL)
}

run_dataset <- function(ds) {
  d <- stress_gene_frame(ds)
  L <- stress_all_labellings(d$animals, d$true_group)
  true_idx <- which(apply(L, 1, function(l) identical(unname(l), unname(d$true_group))))
  stopifnot(length(true_idx) == 1L)
  sets_all <- go_bp_sets(d$genes)
  tm <- map_go_terms_to_manuscript_themes(data.frame(ID = names(sets_all), Description = AnnotationDbi::Term(names(sets_all))),
                                          THEME_REGISTRY)$assignments
  theme_terms <- unique(tm$GO_ID[tm$theme_id %in% THEMES])
  E <- list()
  for (s in sigs) {
    mu <- intersect(s$match_units, d$units)
    ctr <- colnames(s$stat)
    if (!length(mu) || !length(ctr)) next
    st <- s$stat[apply(is.finite(s$stat), 1, all), , drop = FALSE]
    g <- intersect(d$genes, rownames(st))
    if (length(g) < 500L) stop("Fewer than 500 shared genes: ", s$id, " / ", ds, call. = FALSE)
    sets <- lapply(sets_all, function(x) intersect(x, g))
    sets <- sets[lengths(sets) >= GO_MIN & lengths(sets) <= GO_MAX]
    A <- Matrix::sparseMatrix(i = rep(seq_along(sets), lengths(sets)), j = match(unlist(sets), g), x = rep(1 / lengths(sets), lengths(sets)),
                              dims = c(length(sets), length(g)), dimnames = list(names(sets), g))
    cols <- as.vector(outer(mu, ctr, paste, sep = "|"))
    E[[s$id]] <- list(id = s$id, g = g, cols = cols, ctr = sub(".*[|]", "", cols), r = apply(st[g, , drop = FALSE], 2, rank),
                      A = A, es_r = apply(as.matrix(A %*% st[g, , drop = FALSE]), 2, rank), th = which(rownames(A) %in% theme_terms))
  }
  score <- function(t) unlist(unname(lapply(E, function(e) {
    tg <- t[e$g, e$cols, drop = FALSE]
    rg <- apply(tg, 2, rank); os <- as.matrix(e$A %*% tg)
    ro <- apply(os, 2, rank); rth <- apply(os[e$th, , drop = FALSE], 2, rank)
    v <- c(vapply(seq_along(e$cols), function(j) stats::cor(rg[, j], e$r[, e$ctr[j]]), 0),
           vapply(seq_along(e$cols), function(j) stats::cor(ro[, j], e$es_r[, e$ctr[j]]), 0),
           vapply(seq_along(e$cols), function(j) stats::cor(rth[, j], e$es_r[e$th, e$ctr[j]]), 0))
    names(v) <- paste(rep(c("gene", "go_bp", "fig3_themes"), each = length(e$cols)), e$id, e$cols, sep = "|")
    v })))
  cl <- makeCluster(N_WORKERS); on.exit(stopCluster(cl))
  clusterCall(cl, setwd, getwd())
  clusterEvalQ(cl, { source("R/paths.R"); source("R/spatial/stress_effect_utils.R"); suppressPackageStartupMessages(library(Matrix)); NULL })
  clusterExport(cl, c("d", "E", "score"), envir = environment())
  vals <- parLapply(cl, seq_len(nrow(L)), function(i) score(stress_contrast_t(d, L[i, ])))
  list(values = do.call(rbind, vals), true_idx = true_idx, labellings = L)
}

t0 <- Sys.time()
runs <- lapply(STRESS_DATASETS, run_dataset); names(runs) <- STRESS_DATASETS
stopifnot(all(vapply(runs, function(r) identical(r$labellings, runs[[1]]$labellings), TRUE)),
          length(unique(vapply(runs, `[[`, 0L, "true_idx"))) == 1L)
lab <- runs[[1]]$labellings; ti <- runs[[1]]$true_idx
V <- do.call(cbind, lapply(runs, `[[`, "values"))
comp <- rep(names(runs), vapply(runs, function(r) ncol(r$values), 0L))
parts <- do.call(rbind, strsplit(colnames(V), "|", fixed = TRUE)); colnames(parts) <- c("level", "signature", "unit", "contrast")

same <- function(g) apply(lab, 1, function(l) all((l == g) == (lab[ti, ] == g)))
keep <- list(SUS_vs_RES = same("CON"), SUS_vs_CON = same("RES"), RES_vs_CON = same("SUS"), stressed_vs_CON = rep(TRUE, nrow(lab)))
summ <- function(S, ctr) { o <- S[ti]
  c(rho = o, null_median = stats::median(S), null_q95 = unname(stats::quantile(S, 0.95)), p_all = mean(S >= o - 1e-12),
    p_specific = mean(S[keep[[ctr]]] >= o - 1e-12), n_specific = sum(keep[[ctr]])) }

key <- paste(parts[, "level"], parts[, "signature"], parts[, "contrast"], sep = "|")
series <- list(); rows <- list()
for (k in unique(key)) {
  sel <- key == k; lv <- parts[sel, "level"][1]; sg <- parts[sel, "signature"][1]; ct <- parts[sel, "contrast"][1]
  S <- rowMeans(V[, sel, drop = FALSE]); series[[k]] <- S; s <- summ(S, ct)
  byc <- vapply(STRESS_DATASETS, function(d) { w <- sel & comp == d; if (any(w)) mean(V[ti, w]) else NA_real_ }, 0)
  rows[[k]] <- data.frame(level = lv, signature = sg, set = sigs[[sg]]$set, contrast = ct, n_units = sum(sel),
                          units_positive = sum(V[ti, sel] > 0), mean_rho = s[["rho"]], null_median = s[["null_median"]],
                          null_q95 = s[["null_q95"]], p_all_1680 = s[["p_all"]], n_relabellings = nrow(lab),
                          p_contrast_specific = s[["p_specific"]],
                          n_contrast_specific = s[["n_specific"]], rho_neuropil = byc[["neuron_neuropil"]],
                          rho_soma = byc[["neuron_soma"]], rho_microglia = byc[["microglia"]], stringsAsFactors = FALSE)
}
per_sig <- do.call(rbind, rows); rownames(per_sig) <- NULL
per_sig$q_bh <- NA_real_
for (lv in unique(per_sig$level)) { w <- per_sig$level == lv & per_sig$set == "primary"; per_sig$q_bh[w] <- stats::p.adjust(per_sig$p_all_1680[w], "BH") }

primary <- info$id[info$set == "primary"]
om <- list()
for (lv in unique(parts[, "level"])) for (ct in STRESS_CONTRASTS) {
  ks <- paste(lv, primary, ct, sep = "|"); ks <- ks[ks %in% names(series)]
  if (!length(ks)) next
  S <- rowMeans(do.call(cbind, series[ks])); s <- summ(S, ct)
  om[[length(om) + 1L]] <- data.frame(level = lv, contrast = ct, n_signatures = length(ks), signatures = paste(sub("^[^|]*[|]([^|]*)[|].*$", "\\1", ks), collapse = "; "),
                                       mean_rho = s[["rho"]], null_median = s[["null_median"]], null_q95 = s[["null_q95"]],
                                       p_all_1680 = s[["p_all"]], n_relabellings = nrow(lab), p_contrast_specific = s[["p_specific"]],
                                       stringsAsFactors = FALSE)
}
omnibus <- do.call(rbind, om)

units_tab <- data.frame(parts, compartment = comp, rho = V[ti, ], p_all_1680 = colMeans(sweep(V, 2, V[ti, ], ">=")),
                        stringsAsFactors = FALSE, row.names = NULL)
units_tab$set <- vapply(units_tab$signature, function(s) sigs[[s]]$set, "")

agree <- list()
for (ct in STRESS_CONTRASTS) {
  h <- names(sigs)[vapply(sigs, function(s) ct %in% colnames(s$stat), TRUE)]
  for (a in seq_along(h)) for (b in seq_along(h)) if (a < b) {
    g <- intersect(rownames(sigs[[h[a]]]$stat), rownames(sigs[[h[b]]]$stat))
    agree[[length(agree) + 1L]] <- data.frame(contrast = ct, signature_a = h[a], signature_b = h[b], set_a = sigs[[h[a]]]$set,
      set_b = sigs[[h[b]]]$set, n_genes = length(g), rho = stats::cor(sigs[[h[a]]]$stat[g, ct], sigs[[h[b]]]$stat[g, ct], method = "spearman"),
      stringsAsFactors = FALSE)
  }
}
agreement <- do.call(rbind, agree)

utils::write.csv(per_sig, OUT$signatures, row.names = FALSE)
utils::write.csv(omnibus, OUT$omnibus, row.names = FALSE)
utils::write.csv(units_tab, OUT$units, row.names = FALSE)
utils::write.csv(agreement, OUT$agreement, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = c(as.list(gct_paths), list(signature_info = SIG_INFO, theme_registry = THEME_REGISTRY),
                                        as.list(stats::setNames(sig_files, basename(sig_files)))),
                   outputs = OUT,
                   parameters = list(our_model = "Figure 3 per-unit model on the animal-level values; gene level = one protein group per symbol",
                                     statistic = "Spearman of our t with the published statistic over shared genes; mean over matched units",
                                     levels = c("gene", sprintf("GO-BP term means (%d-%d shared genes)", GO_MIN, GO_MAX), "Figure 3 theme terms (descriptive)"),
                                     null = sprintf("all %d assignments of the %d animals to the groups with the observed sizes (%s); contrast-specific subsets keep the third group fixed",
                                                    nrow(lab), ncol(lab), paste(names(table(lab[ti, ])), table(lab[ti, ]), sep = " ", collapse = ", ")),
                                     primary_signatures = primary, workers = N_WORKERS,
                                     runtime_minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
                                     provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "One-sided (agreement) exact relabelling p values; BH over signature x contrast within level, primary set.")
message("External concordance complete: ", nrow(per_sig), " signature x contrast x level rows")
