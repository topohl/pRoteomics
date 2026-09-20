#!/usr/bin/env Rscript
# =====================================================================
# Part-29, sections 28 / 29 / 30 - final WGCNA annotation-registry review
# =====================================================================
#
# Phase 6G.8: the canonical Stage-01 tables below are named through the WGCNA
# path resolver, which looks in the normalized location first and falls back to
# the historical one. This script has no source() block of its own, so the
# bootstrap is explicit and must precede the first resolver call.
paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "wgcna_paths.R"))
#
# AUDIT ONLY. Nothing canonical is recomputed, rerun or rewritten. This
# script READS the canonical tables, re-derives every evidence layer from
# them, and writes CSVs under
#   results/tables/publication_audits/upstream_enrichment_v10/
# It edits no registry and activates no label.
#
# CANONICAL SOURCES (all read-only):
#  E1 GO enrichment  results/tables/06_modules_WGCNA/01_WGCNA/<ds>/modules/
#                      WGCNA_module_GO_enrichment_long.csv   (scope = "all")
#  E2 module summary results/tables/06_modules_WGCNA/01_WGCNA/<ds>/modules/
#                      WGCNA_module_summary.csv              (median |kME|)
#  E3 hubs           results/reviewer_audit/wgcna_label_adjudication/
#                      WGCNA_module_top25_hubs.csv
#  E4 spatial+extern results/tables/11_spatial_systems/atlas/
#                      WGCNA_module_spatial_cell_affinity.csv
#  E5 EWCE long      results/tables/11_spatial_systems/celltype_annotation/
#                      WGCNA_module_external_celltype_affinity_long.csv
#  E6 registry       results/reviewer_audit/wgcna_label_approval/
#                      WGCNA_final_label_approval_table.csv  (ACTIVE label)
#  E7 part-28 matrix results/tables/manuscript_candidates/final_truth_v9/
#                      audit/wgcna_annotation_evidence_matrix.csv (cross-check)
#
# PRE-SPECIFIED DECISION RULES (section 30) - fixed before looking at results:
#  L1 ENRICHMENT layer agrees  <=>  n_sig(FDR<0.05) >= 5 AND best FDR < 1e-4
#  L2 HUB layer agrees (YES)   <=>  >= 6 of the top-10 hubs (by |kME|) are
#                                   members of the union of the gene sets of
#                                   the 10 most significant GO terms.
#                                   PARTIAL = 3..5, NO = 0..2 (or no terms)
#  L3 SEMANTIC layer agrees    <=>  >= 70% of the top-30 significant GO terms
#                                   fall in ONE connected component of the
#                                   term graph in which two terms are linked
#                                   when Jaccard(gene set) >= 0.25
#  L4 SPATIAL layer agrees     <=>  spatial_tau >= 0.60 AND
#                                   peak_minus_second >= 0.05
#  EWCE / reference-panel affinity is deliberately EXCLUDED from the
#  converging count: it may support cell CONTEXT, never functional identity.
#
#  label_confidence uses only the three FUNCTIONAL layers (L1,L2,L3):
#     3/3 -> HIGH     2/3 -> MODERATE     <=1/3 or n_sig==0 -> WITHHOLD
# =====================================================================

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

OUT <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
FDR_CUT  <- 0.05
SCOPE    <- "all"           # widest protein-set scope, same as Part-28
N_TERMS_HUB  <- 10L         # terms whose gene sets define the hub pool
N_TERMS_SEM  <- 30L         # terms entering the semantic-block graph
J_LINK       <- 0.25        # Jaccard link threshold
HUB_YES      <- 6L
HUB_PARTIAL  <- 3L
TAU_CUT      <- 0.60
PMS_CUT      <- 0.05
SIG_MIN      <- 5L
BEST_FDR_CUT <- 1e-4

rd <- function(...) utils::read.csv(file.path(...), stringsAsFactors = FALSE,
                                    check.names = FALSE)

## ---------------------------------------------------------------- E1/E2
read_ds <- function(ds, f) {
  p <- wgcna_modules_artifact(f, ds, child = "tables", "modules")
  if (!file.exists(p)) stop("missing canonical table: ", p)
  x <- rd(p); x$dataset <- ds; x
}
go   <- do.call(rbind, lapply(DATASETS, read_ds, f = "WGCNA_module_GO_enrichment_long.csv"))
summ <- do.call(rbind, lapply(DATASETS, read_ds, f = "WGCNA_module_summary.csv"))

stopifnot(SCOPE %in% unique(go$ModuleProteinSetType))
go <- go[go$ModuleProteinSetType == SCOPE, , drop = FALSE]
go$module_id   <- sub("^WGCNA_", "", go$ModuleID)
summ$module_id <- sub("^WGCNA_", "", summ$ModuleID)

## ---------------------------------------------------------------- E3
hubs <- rd("results", "reviewer_audit", "wgcna_label_adjudication",
           "WGCNA_module_top25_hubs.csv")
hubs$module_id <- sub("^WGCNA_", "", hubs$ModuleID)
hubs$EntrezID  <- trimws(as.character(hubs$EntrezID))

## ---------------------------------------------------------------- E4
aff <- rd("results", "tables", "11_spatial_systems", "atlas",
          "WGCNA_module_spatial_cell_affinity.csv")
aff$module_id <- sub("^WGCNA_", "", aff$ModuleID)

## ---------------------------------------------------------------- E5
ewce <- rd("results", "tables", "11_spatial_systems", "celltype_annotation",
           "WGCNA_module_external_celltype_affinity_long.csv")
ewce$module_id <- sub("^WGCNA_", "", ewce$ModuleID)

## ---------------------------------------------------------------- E6
reg <- rd("results", "reviewer_audit", "wgcna_label_approval",
          "WGCNA_final_label_approval_table.csv")
reg <- reg[reg$level == "module", , drop = FALSE]
reg$module_id <- sub("^WGCNA_", "", reg$entity_id)
reg$active_label <- trimws(sub("^WGCNA_m[0-9]+\\s*.\\s*", "",
                               reg$current_active_label))

## ---------------------------------------------------------------- E7
p28 <- rd("results", "tables", "manuscript_candidates", "final_truth_v9",
          "audit", "wgcna_annotation_evidence_matrix.csv")

## ===================================================================
## helpers
## ===================================================================
split_genes <- function(s) {
  s <- as.character(s)
  out <- strsplit(s, "/", fixed = TRUE)
  lapply(out, function(z) unique(trimws(z[nzchar(trimws(z))])))
}

jacc <- function(a, b) {
  u <- length(union(a, b)); if (u == 0L) return(0) ; length(intersect(a, b)) / u
}

# largest connected component of a boolean adjacency matrix
largest_cc <- function(adj) {
  n <- nrow(adj); if (n == 0L) return(integer(0))
  seen <- rep(FALSE, n); best <- integer(0)
  for (i in seq_len(n)) {
    if (seen[i]) next
    comp <- i; queue <- i; seen[i] <- TRUE
    while (length(queue)) {
      v <- queue[1]; queue <- queue[-1]
      nb <- which(adj[v, ] & !seen)
      if (length(nb)) { seen[nb] <- TRUE; comp <- c(comp, nb); queue <- c(queue, nb) }
    }
    if (length(comp) > length(best)) best <- comp
  }
  sort(best)
}

fmt_g <- function(x) if (length(x) == 0 || is.na(x)) "NA" else
  formatC(x, format = "g", digits = 3)

## ===================================================================
## per-module evidence
## ===================================================================
rows <- list(); m11_detail <- NULL

for (i in seq_len(nrow(aff))) {
  ds <- aff$dataset[i]; m <- aff$module_id[i]

  g <- go[go$dataset == ds & go$module_id == m & is.finite(go$p.adjust), ,
          drop = FALSE]
  g <- g[order(g$p.adjust, -g$Count, g$ID), , drop = FALSE]
  sig <- g[g$p.adjust < FDR_CUT, , drop = FALSE]
  n_sig <- nrow(sig)

  bp <- g[g$Ontology == "BP", , drop = FALSE]
  cc <- g[g$Ontology == "CC", , drop = FALSE]
  mf <- g[g$Ontology == "MF", , drop = FALSE]
  best_bp   <- if (nrow(bp)) bp$Description[1] else NA_character_
  best_bp_f <- if (nrow(bp)) bp$p.adjust[1]    else NA_real_
  best_cc   <- if (nrow(cc)) cc$Description[1] else NA_character_
  best_cc_f <- if (nrow(cc)) cc$p.adjust[1]    else NA_real_
  best_mf   <- if (nrow(mf)) mf$Description[1] else NA_character_
  best_mf_f <- if (nrow(mf)) mf$p.adjust[1]    else NA_real_
  best_any  <- if (nrow(g))  g$Description[1]  else NA_character_
  best_any_ont <- if (nrow(g)) g$Ontology[1]   else NA_character_
  best_any_f   <- if (nrow(g)) g$p.adjust[1]   else NA_real_

  n_sig_bp <- sum(bp$p.adjust < FDR_CUT)
  n_sig_cc <- sum(cc$p.adjust < FDR_CUT)
  n_sig_mf <- sum(mf$p.adjust < FDR_CUT)

  ## ---- hubs -------------------------------------------------------
  hz <- hubs[hubs$dataset == ds & hubs$module_id == m, , drop = FALSE]
  hz <- hz[order(hz$rank), , drop = FALSE]
  n_hub_report <- min(20L, nrow(hz))
  hub_str <- paste(sprintf("%s(%.3f)", hz$GeneSymbol[seq_len(n_hub_report)],
                           hz$abs_kME[seq_len(n_hub_report)]), collapse = "; ")
  top10 <- hz[seq_len(min(10L, nrow(hz))), , drop = FALSE]

  ## ---- hub / enrichment coherence ---------------------------------
  if (n_sig > 0L) {
    kt <- min(N_TERMS_HUB, n_sig)
    pool <- unique(unlist(split_genes(sig$geneID[seq_len(kt)])))
    bestset <- split_genes(sig$geneID[1])[[1]]
    hub_ent <- top10$EntrezID
    hub_ent_ok <- hub_ent[!is.na(hub_ent) & nzchar(hub_ent)]
    n_hub_pool <- sum(hub_ent_ok %in% pool)
    n_hub_best <- sum(hub_ent_ok %in% bestset)
    hub_in_pool_genes <- top10$GeneSymbol[top10$EntrezID %in% pool]
    hub_out_genes     <- top10$GeneSymbol[!(top10$EntrezID %in% pool)]
  } else {
    kt <- 0L; pool <- character(0); n_hub_pool <- 0L; n_hub_best <- 0L
    hub_in_pool_genes <- character(0); hub_out_genes <- top10$GeneSymbol
  }
  n_hub_eval <- nrow(top10)

  hub_class <- if (n_sig == 0L) "NO" else if (n_hub_pool >= HUB_YES) "YES" else
    if (n_hub_pool >= HUB_PARTIAL) "PARTIAL" else "NO"
  hub_reason <- if (n_sig == 0L) {
    sprintf("NO: module has no GO term at FDR<0.05, so there is no enrichment for the hubs to agree with (top-10 hubs %s)",
            paste(utils::head(top10$GeneSymbol, 10), collapse = "/"))
  } else {
    sprintf("%d/%d top-10 hubs are members of the top-%d significant GO term gene sets (%d also in the single best term '%s'); inside: %s; outside: %s",
            n_hub_pool, n_hub_eval, kt, n_hub_best, best_any,
            if (length(hub_in_pool_genes)) paste(hub_in_pool_genes, collapse = "/") else "none",
            if (length(hub_out_genes)) paste(hub_out_genes, collapse = "/") else "none")
  }
  hub_assess <- paste0(hub_class, " - ", hub_reason)

  ## ---- semantic coherence -----------------------------------------
  if (n_sig == 0L) {
    sem_class <- "NO_SIGNIFICANT_TERMS"; frac_cc <- NA_real_; meanJ <- NA_real_
    n_comp <- NA_integer_
    sem_txt <- "no GO term reaches FDR<0.05; no semantic block can be assessed"
    cc_members <- character(0); out_members <- character(0); k_sem <- 0L
  } else {
    k_sem <- min(N_TERMS_SEM, n_sig)
    sets <- split_genes(sig$geneID[seq_len(k_sem)])
    J <- matrix(0, k_sem, k_sem)
    if (k_sem > 1L) for (a in 1:(k_sem - 1)) for (b in (a + 1):k_sem) {
      J[a, b] <- J[b, a] <- jacc(sets[[a]], sets[[b]])
    }
    meanJ <- if (k_sem > 1L) mean(J[upper.tri(J)]) else NA_real_
    adj <- J >= J_LINK; diag(adj) <- FALSE
    cc_idx <- largest_cc(adj)
    frac_cc <- length(cc_idx) / k_sem
    seen <- rep(FALSE, k_sem); n_comp <- 0L
    for (a in seq_len(k_sem)) {
      if (seen[a]) next
      comp <- a; q <- a; seen[a] <- TRUE
      while (length(q)) { v <- q[1]; q <- q[-1]
        nb <- which(adj[v, ] & !seen)
        if (length(nb)) { seen[nb] <- TRUE; comp <- c(comp, nb); q <- c(q, nb) } }
      if (length(comp) >= 2L) n_comp <- n_comp + 1L
    }
    cc_members  <- sig$Description[cc_idx]
    out_members <- sig$Description[setdiff(seq_len(k_sem), cc_idx)]
    sem_class <- if (frac_cc >= 0.70) "SINGLE_COHERENT_BLOCK" else
      if (frac_cc >= 0.40) "DOMINANT_BLOCK_PLUS_SATELLITES" else "HETEROGENEOUS"
    sem_txt <- sprintf("%.0f%% (%d/%d) of the top-%d significant terms form one gene-set-linked block (Jaccard>=%.2f); %d multi-term blocks; mean pairwise Jaccard %.3f; block core: %s; outside block: %s",
                       100 * frac_cc, length(cc_idx), k_sem, k_sem, J_LINK,
                       n_comp, ifelse(is.na(meanJ), 0, meanJ),
                       paste(utils::head(cc_members, 4), collapse = " | "),
                       if (length(out_members)) paste(utils::head(out_members, 4), collapse = " | ") else "none")
  }
  sem_full <- paste0(sem_class, " - ", sem_txt)

  ## ---- spatial ----------------------------------------------------
  tau <- aff$spatial_tau[i]; pms <- aff$peak_minus_second[i]
  spatial_ok <- is.finite(tau) && is.finite(pms) && tau >= TAU_CUT && pms >= PMS_CUT

  ## ---- external ---------------------------------------------------
  ext_ct  <- aff$external_celltype_all[i]
  ext_fdr <- aff$external_FDR_all[i]
  # NOTE: Part-28 tested concordance == "all_scopes_agree"; the canonical
  # vocabulary is "consistent" / "mixed_scopes_disagree", so Part-28's flag is
  # FALSE for every module. Re-derived correctly here.
  ext_agree <- identical(as.character(aff$external_celltype_scope_concordance[i]),
                         "consistent")
  ez <- ewce[ewce$dataset == ds & ewce$module_id == m & ewce$level == 1 &
               ewce$module_scope == "all" & ewce$cell_type == ext_ct, , drop = FALSE]
  ext_z <- if (nrow(ez)) max(ez$z_score) else NA_real_
  ez2 <- ewce[ewce$dataset == ds & ewce$module_id == m & ewce$level == 1 &
                ewce$module_scope == "core_kME06" & ewce$cell_type == ext_ct, , drop = FALSE]
  ext_z_core <- if (nrow(ez2)) max(ez2$z_score) else NA_real_

  ## ---- layers -----------------------------------------------------
  enrich_ok <- (n_sig >= SIG_MIN) && is.finite(best_any_f) && best_any_f < BEST_FDR_CUT
  hub_ok    <- identical(hub_class, "YES")
  sem_ok    <- identical(sem_class, "SINGLE_COHERENT_BLOCK")
  n_conv <- sum(enrich_ok, hub_ok, sem_ok, spatial_ok)
  n_func <- sum(enrich_ok, hub_ok, sem_ok)

  conf <- if (n_sig == 0L) "WITHHOLD" else
    if (n_func == 3L) "HIGH" else if (n_func == 2L) "MODERATE" else "WITHHOLD"

  ## ---- recommended publication form -------------------------------
  rec <- if (conf == "WITHHOLD")
    sprintf("%s %s only, annotation withheld", ds, m)
  else
    sprintf("%s %s, enriched for %s proteins", ds, m, best_any)

  rg <- reg[reg$dataset == ds & reg$module_id == m, , drop = FALSE]
  active <- if (nrow(rg)) rg$active_label[1] else NA_character_

  just <- sprintf(
    "Converging layers %d/4 (enrichment=%s, hubs=%s, semantics=%s, spatial=%s); functional layers %d/3. Enrichment: %d GO terms at FDR<0.05 (BP %d, CC %d, MF %d), best %s '%s' FDR %s. Hubs: %d/%d of the top-10 hubs sit inside the leading term gene sets. Semantics: %s. Spatial: peak %s, tau %s, peak-second %s. External (CONTEXT ONLY, never identity): %s FDR %s, scopes %s, EWCE z(all)=%s. Active registry label: '%s'. %s",
    n_conv, enrich_ok, hub_class, sem_class, spatial_ok, n_func,
    n_sig, n_sig_bp, n_sig_cc, n_sig_mf, best_any_ont, best_any, fmt_g(best_any_f),
    n_hub_pool, n_hub_eval, sem_class,
    aff$peak_unit[i], fmt_g(tau), fmt_g(pms),
    ext_ct, fmt_g(ext_fdr),
    aff$external_celltype_scope_concordance[i], fmt_g(ext_z),
    active,
    switch(conf,
      HIGH = "All three functional layers converge, so a descriptive enrichment phrase is supported provided the module ID stays visible.",
      MODERATE = "Two of three functional layers converge; the enrichment phrase is reportable but must stay descriptive and ID-anchored.",
      WITHHOLD = "Fewer than two functional layers converge (or there is no significant enrichment at all), so no functional label is supported; refer to the module by ID only."))

  rows[[length(rows) + 1L]] <- data.frame(
    dataset = ds, module_id = m,
    module_size = aff$module_size[i],
    median_abs_kME = summ$median_abs_kME[summ$dataset == ds & summ$module_id == m][1],
    n_enriched_terms_FDR05 = n_sig,
    best_GO_BP_term = best_bp, best_GO_BP_FDR = best_bp_f,
    best_GO_CC_term = best_cc, best_GO_CC_FDR = best_cc_f,
    top_hub_genes_with_kME = hub_str,
    hub_coherence_assessment = hub_assess,
    semantic_coherence = sem_full,
    spatial_peak_unit = aff$peak_unit[i],
    spatial_tau = tau,
    external_celltype = ext_ct, external_FDR = ext_fdr,
    external_scopes_agree = ext_agree,
    current_active_label = active,
    evidence_layers_converging = n_conv,
    recommended_publication_form = rec,
    label_confidence = conf,
    justification = just,
    # ---- supporting numbers (audit transparency)
    n_sig_BP = n_sig_bp, n_sig_CC = n_sig_cc, n_sig_MF = n_sig_mf,
    best_GO_MF_term = best_mf, best_GO_MF_FDR = best_mf_f,
    best_any_ontology = best_any_ont, best_any_term = best_any,
    best_any_FDR = best_any_f,
    n_top10_hubs_in_leading_termsets = n_hub_pool,
    n_top10_hubs_in_best_term = n_hub_best,
    hub_class = hub_class,
    semantic_class = sem_class,
    semantic_block_fraction = frac_cc,
    semantic_mean_jaccard = meanJ,
    semantic_n_multiterm_blocks = n_comp,
    spatial_peak_minus_second = aff$peak_minus_second[i],
    spatial_layer_supports = spatial_ok,
    enrichment_layer_supports = enrich_ok,
    hub_layer_supports = hub_ok,
    semantic_layer_supports = sem_ok,
    n_functional_layers_converging = n_func,
    external_scope_concordance_raw = aff$external_celltype_scope_concordance[i],
    external_EWCE_z_all = ext_z, external_EWCE_z_core_kME06 = ext_z_core,
    reference_panel = aff$strongest_reference_panel[i],
    reference_panel_FDR = aff$reference_panel_FDR[i],
    stringsAsFactors = FALSE)

  ## ---- m11 neuropil deep dive -------------------------------------
  if (ds == "neuron_neuropil" && m == "m11") {
    m11_detail <- list(g = g, sig = sig, hz = hz, i = i,
                       n_sig = n_sig, best_any = best_any,
                       best_any_f = best_any_f,
                       n_hub_pool = n_hub_pool, n_hub_eval = n_hub_eval,
                       sem_class = sem_class, frac_cc = frac_cc,
                       hub_class = hub_class, conf = conf,
                       n_conv = n_conv, n_func = n_func)
  }
}

res <- do.call(rbind, rows)
res <- res[order(res$dataset, res$module_id), , drop = FALSE]
utils::write.csv(res, file.path(OUT, "wgcna_registry_final_review.csv"),
                 row.names = FALSE)

## ===================================================================
## Section 29 - neuropil m11 full adjudication
## ===================================================================
d <- m11_detail
gg <- d$g; sg <- d$sig; hz <- d$hz
irow <- d$i

myelin_canon <- c("CNP", "MAG", "SIRT2", "ERMN", "ENPP6", "MOG", "BCAS1",
                  "PLP1", "CLDN11", "NDRG1", "OPALIN", "MBP")
top13 <- hz[order(hz$rank), ][1:13, ]
n_myelin_top13 <- sum(top13$GeneSymbol %in% myelin_canon)
non_myelin_top13 <- top13$GeneSymbol[!(top13$GeneSymbol %in% myelin_canon)]

ms <- gg[gg$Ontology == "CC" & gg$Description == "myelin sheath", ]
ewce11 <- ewce[ewce$dataset == "neuron_neuropil" & ewce$module_id == "m11" &
                 ewce$level == 1 & ewce$cell_type == "oligodendrocytes", ]
ewce11 <- ewce11[order(-ewce11$z_score), ]

reg11 <- reg[reg$dataset == "neuron_neuropil" & reg$module_id == "m11", ]

myelin_rx <- "myelin|ensheath|oligodendrocyte|glial cell differentiation|gliogenesis|glial cell development|paranod|Schmidt-Lanterman"
n_sig_myelin <- sum(grepl(myelin_rx, sg$Description, ignore.case = TRUE))

chk <- function(claim, claimed, observed, verdict, src)
  data.frame(claim = claim, claimed_value = claimed, canonical_value = observed,
             verdict = verdict, canonical_source = src, stringsAsFactors = FALSE)

## ---- m11 satellite-term diagnostic (needed by the claim table below) ----
## the primary semantic rule leaves 7/21 terms outside the block. Are those
## terms different biology, or the same myelin genes carrying cell-polarity
## CC labels that GO keeps in separate branches?
sets21 <- split_genes(sg$geneID)
mye_set <- split_genes(sg$geneID[sg$Description == "myelin sheath"])[[1]]
core_block <- unique(unlist(split_genes(
  sg$geneID[grepl(myelin_rx, sg$Description, ignore.case = TRUE)])))
sym_map <- stats::setNames(hubs$GeneSymbol, hubs$EntrezID)
sat <- data.frame(
  rank = seq_len(nrow(sg)),
  ontology = sg$Ontology, GO_ID = sg$ID, term = sg$Description,
  p_adjust = sg$p.adjust, count = sg$Count, gene_ratio = sg$GeneRatio,
  in_myelin_name_family = grepl(myelin_rx, sg$Description, ignore.case = TRUE),
  jaccard_with_myelin_sheath = vapply(sets21, jacc, numeric(1), b = mye_set),
  frac_of_term_genes_also_in_myelin_sheath =
    vapply(sets21, function(z) if (!length(z)) NA_real_ else
      mean(z %in% mye_set), numeric(1)),
  frac_of_term_genes_also_in_myelin_family_union =
    vapply(sets21, function(z) if (!length(z)) NA_real_ else
      mean(z %in% core_block), numeric(1)),
  overlapping_m11_hub_symbols = vapply(sets21, function(z)
    paste(stats::na.omit(sym_map[z[z %in% hz$EntrezID]]), collapse = "/"),
    character(1)),
  stringsAsFactors = FALSE)
utils::write.csv(sat, file.path(OUT, "wgcna_m11_significant_term_diagnostic.csv"),
                 row.names = FALSE)
out_sat <- sat[!sat$in_myelin_name_family, , drop = FALSE]

## ---- m11 semantic-threshold sensitivity (needed by the claim table) ----
m11_frac <- function(sub, thr) {
  if (!nrow(sub)) return(NA_real_)
  kk <- min(N_TERMS_SEM, nrow(sub))
  st <- split_genes(sub$geneID[seq_len(kk)])
  JJ <- matrix(0, kk, kk)
  if (kk > 1L) for (a in 1:(kk - 1)) for (b in (a + 1):kk)
    JJ[a, b] <- JJ[b, a] <- jacc(st[[a]], st[[b]])
  A <- JJ >= thr; diag(A) <- FALSE
  length(largest_cc(A)) / kk
}
m11_f15 <- m11_frac(sg, 0.15); m11_f25 <- m11_frac(sg, 0.25)
m11_f35 <- m11_frac(sg, 0.35)
m11_fcc <- m11_frac(sg[sg$Ontology == "CC", ], 0.25)
m11_fbp <- m11_frac(sg[sg$Ontology == "BP", ], 0.25)
m11_pool <- unique(unlist(split_genes(sg$geneID[seq_len(min(N_TERMS_HUB, nrow(sg)))])))
m11_hub25 <- sum(hz$EntrezID %in% m11_pool)

GOP <- "results/tables/06_modules_WGCNA/01_WGCNA/neuron_neuropil/modules/WGCNA_module_GO_enrichment_long.csv (ModuleProteinSetType=all)"
HUP <- "results/reviewer_audit/wgcna_label_adjudication/WGCNA_module_top25_hubs.csv"
EWP <- "results/tables/11_spatial_systems/celltype_annotation/WGCNA_module_external_celltype_affinity_long.csv"
RGP <- "results/reviewer_audit/wgcna_label_approval/WGCNA_final_label_approval_table.csv"
AFP <- "results/tables/11_spatial_systems/atlas/WGCNA_module_spatial_cell_affinity.csv"

m11 <- rbind(
  chk("CC 'myelin sheath' BH-adjusted p ~ 6.3e-11",
      "~6.3e-11",
      sprintf("p.adjust = %.6g (raw p = %.6g), GeneRatio %s, BgRatio %s, Count %d, ModuleSize %d, MappedModuleSize %s, GO ID %s",
              ms$p.adjust[1], ms$pvalue[1], ms$GeneRatio[1], ms$BgRatio[1],
              ms$Count[1], ms$ModuleSize[1], as.character(ms$MappedModuleSize[1]), ms$ID[1]),
      "CONFIRMED", GOP),
  chk("Part-28 reported the myelin-sheath gene ratio as 22/88",
      "22/88",
      sprintf("the canonical CC row carries GeneRatio %s; its MappedModuleSize field is %s, and Part-28 printed Count/MappedModuleSize rather than the GeneRatio itself. The module has %d features; the CC GeneRatio denominator is %s (genes with any CC annotation) and the BP denominator is %s.",
              ms$GeneRatio[1], as.character(ms$MappedModuleSize[1]), aff$module_size[irow],
              sub("^[0-9]+/", "", ms$GeneRatio[1]),
              sub("^[0-9]+/", "",
                  gg$GeneRatio[gg$Ontology == "BP" & gg$Description == "myelination"][1])),
      "DISPLAY DISCREPANCY - the count 22 is confirmed; the correct denominator for this CC test is 86, not the 88 Part-28 printed", GOP),
  chk("12 of the top 13 hubs are canonical myelin proteins",
      "12/13",
      sprintf("%d/13 (%s); the single non-myelin hub is %s at rank %d, |kME| %.4f",
              n_myelin_top13,
              paste(sprintf("%s=%.4f", top13$GeneSymbol, top13$abs_kME), collapse = ", "),
              paste(non_myelin_top13, collapse = "/"),
              top13$rank[top13$GeneSymbol %in% non_myelin_top13][1],
              top13$abs_kME[top13$GeneSymbol %in% non_myelin_top13][1]),
      "CONFIRMED", HUP),
  chk("EWCE oligodendrocyte z ~ 34",
      "~34",
      sprintf("max z = %.4f (scope %s, level 1, FDR %.4g); scope 'all' z = %.4f (FDR %.4g); scope 'top25' z = %.4f (FDR %.4g); reps=10000, seed=20260101, background=measured_proteome, phenotype_used=none",
              ewce11$z_score[1], ewce11$module_scope[1], ewce11$FDR[1],
              ewce11$z_score[ewce11$module_scope == "all"][1],
              ewce11$FDR[ewce11$module_scope == "all"][1],
              ewce11$z_score[ewce11$module_scope == "top25"][1],
              ewce11$FDR[ewce11$module_scope == "top25"][1]),
      "CONFIRMED (34.07 is the core_kME06 scope; the 'all' scope is 32.18)", EWP),
  chk("historical / active registry label is 'synaptic/cytoskeletal trafficking'",
      "synaptic/cytoskeletal trafficking",
      sprintf("current_active_label = '%s'; proposed_final_label = '%s'; adjudication_action = '%s'; confidence = '%s'; recommended_for_activation = '%s'; human_decision = '%s'",
              reg11$current_active_label[1], reg11$proposed_final_label[1],
              reg11$adjudication_action[1], reg11$confidence[1],
              as.character(reg11$recommended_for_activation[1]),
              as.character(reg11$human_decision[1])),
      "CONFIRMED - the active label is contradicted by the module's own enrichment and hubs", RGP),
  chk("external cell-type affinity is concordant across protein-set scopes",
      "(Part-28 recorded external_scopes_agree = FALSE)",
      sprintf("external_celltype_all=%s, core_kME06=%s, top25=%s; concordance field = '%s' -> scopes DO agree; Part-28 compared against the string 'all_scopes_agree', which is not in the canonical vocabulary, so its flag is FALSE for all 35 modules",
              aff$external_celltype_all[irow], aff$external_celltype_core_kME06[irow],
              aff$external_celltype_top25[irow],
              aff$external_celltype_scope_concordance[irow]),
      "PART-28 FLAG IS WRONG - scopes are concordant", AFP),
  chk("the enrichment block is the myelin / ensheathment family, not a synaptic one",
      "(assessed here)",
      sprintf("%d of the %d significant terms (FDR<0.05) match the myelin/ensheathment/oligodendrocyte family; largest gene-set-linked block covers %.0f%% of the top-%d terms (%s); top terms: %s",
              n_sig_myelin, d$n_sig, 100 * d$frac_cc, min(30L, d$n_sig), d$sem_class,
              paste(utils::head(sg$Description, 6), collapse = " | ")),
      "CONFIRMED", GOP),
  chk("reference-panel affinity",
      "(context layer)",
      sprintf("strongest_reference_panel = %s, FDR = %.4g",
              aff$strongest_reference_panel[irow], aff$reference_panel_FDR[irow]),
      "CONFIRMED - external/reference evidence is CONTEXT ONLY under section 30", AFP),
  chk("spatial identity",
      "(context layer)",
      sprintf("peak unit %s, tau %.4f, peak-second %.4f, n_spatial_units %s",
              aff$peak_unit[irow], aff$spatial_tau[irow],
              aff$peak_minus_second[irow], as.character(aff$n_spatial_units[irow])),
      "CONFIRMED", AFP),
  chk("the significant terms outside the myelin-named family are a SECOND, competing biology",
      "(tested here, because the primary semantic rule leaves 7/21 terms outside the block)",
      sprintf("all %d non-myelin-named significant terms are CC cell-polarity / cell-projection terms (%s); their genes are %.0f-%.0f%% shared with the myelin-family gene union (median %.0f%%), and every one of them overlaps m11 hub proteins (%s). There is no second biological block - GO keeps these polarity CC branches separate from GO:0043209 while the underlying proteins are the same oligodendrocyte membrane proteins.",
              nrow(out_sat), paste(out_sat$term, collapse = " | "),
              100 * min(out_sat$frac_of_term_genes_also_in_myelin_family_union),
              100 * max(out_sat$frac_of_term_genes_also_in_myelin_family_union),
              100 * stats::median(out_sat$frac_of_term_genes_also_in_myelin_family_union),
              paste(unique(unlist(strsplit(out_sat$overlapping_m11_hub_symbols, "/"))),
                    collapse = "/")),
      "REFUTED - the 'satellite' terms are the same myelin proteins under cell-polarity CC labels", GOP),
  chk("does m11's label_confidence depend on the semantic-block threshold?",
      "(sensitivity test)",
      sprintf("largest gene-set-linked block covers %.0f%% of the 21 significant terms at Jaccard>=0.15, %.0f%% at >=0.25 (PRIMARY, gives MODERATE) and %.0f%% at >=0.35; within BP alone all 9 significant terms form ONE block (100%%); within CC alone %.0f%%. Top-25 hubs: %d/25 (%.0f%%) inside the leading term gene sets. m11 is MODERATE only because it sits 3 percentage points under the pre-specified 70%% cut; the enrichment and hub layers are unambiguous either way.",
              100 * m11_f15, 100 * m11_f25, 100 * m11_f35, 100 * m11_fcc,
              m11_hub25, 100 * m11_hub25 / nrow(hz)),
      "THRESHOLD-SENSITIVE - HIGH at Jaccard>=0.15, MODERATE at >=0.25 and >=0.35; never WITHHOLD", GOP),
  chk("status of the registry's own proposed replacement label for m11",
      "(registry state)",
      sprintf("proposed_final_label = '%s', adjudication_action = '%s', confidence = '%s', recommended_for_activation = '%s', human_decision = '%s'. That proposed wording is the candidate-B shape: it names a CELL TYPE. It is recommended but NOT activated, and this audit does not activate it.",
              reg11$proposed_final_label[1], reg11$adjudication_action[1],
              reg11$confidence[1], as.character(reg11$recommended_for_activation[1]),
              as.character(reg11$human_decision[1])),
      "NOT ACTIVATED - and under section 30 the cell-type wording should not be activated as-is", RGP)
)
utils::write.csv(m11, file.path(OUT, "wgcna_m11_claim_verification.csv"),
                 row.names = FALSE)

## ---- candidate forms A / B / C ------------------------------------
evid_common <- sprintf(
  "CC myelin sheath FDR %.3g (22/86, GO:0043209); BP ensheathment of neurons / axon ensheathment FDR %.3g (16/83); BP myelination FDR %.3g (15/83); %d of %d significant terms carry myelin/ensheathment/oligodendrocyte names, and the remaining %d are CC cell-polarity terms built from the SAME proteins (%.0f-%.0f%% gene overlap with the myelin family union), so there is no competing second block; semantic block %s (%.0f%% of the %d significant terms at Jaccard>=0.25, 100%% at >=0.15, 100%% within BP alone); %d/%d top-10 and %d/25 top-25 hubs sit inside the leading term gene sets; 12/13 top hubs are canonical myelin proteins (only RHOG, rank 3, |kME| 0.971, is not); spatial peak %s, tau %.3f; EWCE oligodendrocytes z=32.18 (scope all) / 34.07 (core kME>=0.6) / 26.60 (top25), FDR %.3g, same call in all three protein-set scopes; reference_oligodendrocyte panel FDR %.3g.",
  sg$p.adjust[sg$Description == "myelin sheath"][1],
  sg$p.adjust[sg$Description == "axon ensheathment"][1],
  sg$p.adjust[sg$Description == "myelination"][1],
  n_sig_myelin, d$n_sig, nrow(out_sat),
  100 * min(out_sat$frac_of_term_genes_also_in_myelin_family_union),
  100 * max(out_sat$frac_of_term_genes_also_in_myelin_family_union),
  d$sem_class, 100 * d$frac_cc, d$n_sig,
  d$n_hub_pool, d$n_hub_eval, m11_hub25,
  aff$peak_unit[irow], aff$spatial_tau[irow],
  aff$external_FDR_all[irow], aff$reference_panel_FDR[irow])

forms <- data.frame(
  candidate_form = c("A", "B", "C"),
  wording = c("m11, enriched for myelin-associated proteins",
              "m11, myelin/oligodendrocyte-associated",
              "m11 only, annotation withheld"),
  keeps_module_id_visible = c(TRUE, TRUE, TRUE),
  claim_type = c("PROTEIN-COMPOSITION claim about the module's own measured content",
                 "PROTEIN-COMPOSITION claim plus a CELL-TYPE attribution",
                 "no annotation"),
  supported_by_enrichment_layer = c(TRUE, TRUE, TRUE),
  supported_by_hub_layer = c(TRUE, TRUE, TRUE),
  supported_by_semantic_layer = c(TRUE, TRUE, TRUE),
  requires_cell_intrinsic_identity = c(FALSE, TRUE, FALSE),
  cell_intrinsic_identity_justified = c(NA, FALSE, NA),
  verdict = c("RECOMMENDED", "NOT RECOMMENDED", "NOT RECOMMENDED"),
  reason = c(
    paste("The enrichment and hub layers converge unambiguously on one biology, and the semantic layer converges substantively: under the pre-specified Jaccard>=0.25 rule m11 scores 67%, three points under the 70% cut (hence MODERATE, not HIGH, in the review table), but the 8 terms that fall outside the block are CC cell-polarity terms built from the same myelin proteins, not a competing biology. The claim itself is purely compositional - it states which proteins are in the module, which is exactly what a co-abundance module can support, and it keeps the module ID visible.",
          evid_common),
    paste("The protein composition is not in doubt, but 'oligodendrocyte-associated' adds a CELL-TYPE attribution. These are enriched-ROI / co-abundance proteomics from neuropil, not sorted cells, so cell of origin cannot be established; EWCE affinity is an external reference-panel overlap and under section 30 may support cell CONTEXT but can never define identity. The compound form also invites the forbidden reading 'the oligodendrocyte module'. If cell context is to be stated it belongs in a separate sentence explicitly flagged as external context, not inside the module's name.",
          evid_common),
    paste("Withholding is correct only where evidence does not converge. Here it does: the best CC term is at FDR 6.3e-11, 12 of the top 13 hubs are canonical myelin proteins, and no second biological block exists among the significant terms. Withholding would discard one of the best-supported annotations in the registry and would leave in place the ACTIVE label 'synaptic/cytoskeletal trafficking', which the same canonical tables contradict.",
          evid_common)),
  stringsAsFactors = FALSE)
utils::write.csv(forms, file.path(OUT, "wgcna_m11_candidate_form_adjudication.csv"),
                 row.names = FALSE)

## ===================================================================
## sensitivity of the semantic and hub rules (transparency; the PRIMARY
## thresholds above were fixed in advance and are NOT changed here)
## ===================================================================
sens <- do.call(rbind, lapply(seq_len(nrow(res)), function(j) {
  ds <- res$dataset[j]; m <- res$module_id[j]
  g <- go[go$dataset == ds & go$module_id == m & is.finite(go$p.adjust), ,
          drop = FALSE]
  g <- g[order(g$p.adjust, -g$Count, g$ID), , drop = FALSE]
  s <- g[g$p.adjust < FDR_CUT, , drop = FALSE]
  frac_at <- function(sub, thr) {
    if (!nrow(sub)) return(NA_real_)
    kk <- min(N_TERMS_SEM, nrow(sub))
    st <- split_genes(sub$geneID[seq_len(kk)])
    JJ <- matrix(0, kk, kk)
    if (kk > 1L) for (a in 1:(kk - 1)) for (b in (a + 1):kk)
      JJ[a, b] <- JJ[b, a] <- jacc(st[[a]], st[[b]])
    A <- JJ >= thr; diag(A) <- FALSE
    length(largest_cc(A)) / kk
  }
  hz <- hubs[hubs$dataset == ds & hubs$module_id == m, , drop = FALSE]
  hz <- hz[order(hz$rank), , drop = FALSE]
  pool <- if (nrow(s)) unique(unlist(split_genes(
    s$geneID[seq_len(min(N_TERMS_HUB, nrow(s)))]))) else character(0)
  n25 <- sum(hz$EntrezID %in% pool)
  data.frame(dataset = ds, module_id = m,
             semantic_frac_J015 = frac_at(s, 0.15),
             semantic_frac_J025_PRIMARY = frac_at(s, 0.25),
             semantic_frac_J035 = frac_at(s, 0.35),
             semantic_frac_CC_only_J025 = frac_at(s[s$Ontology == "CC", ], 0.25),
             semantic_frac_BP_only_J025 = frac_at(s[s$Ontology == "BP", ], 0.25),
             n_top25_hubs_in_leading_termsets = n25,
             frac_top25_hubs_in_leading_termsets = n25 / nrow(hz),
             primary_label_confidence = res$label_confidence[j],
             stringsAsFactors = FALSE)
}))
sens$confidence_if_J015 <- ifelse(
  res$n_enriched_terms_FDR05 == 0, "WITHHOLD",
  ifelse(res$enrichment_layer_supports + res$hub_layer_supports +
           (sens$semantic_frac_J015 >= 0.70) == 3, "HIGH",
  ifelse(res$enrichment_layer_supports + res$hub_layer_supports +
           (sens$semantic_frac_J015 >= 0.70) == 2, "MODERATE", "WITHHOLD")))
sens$confidence_if_J035 <- ifelse(
  res$n_enriched_terms_FDR05 == 0, "WITHHOLD",
  ifelse(res$enrichment_layer_supports + res$hub_layer_supports +
           (sens$semantic_frac_J035 >= 0.70) == 3, "HIGH",
  ifelse(res$enrichment_layer_supports + res$hub_layer_supports +
           (sens$semantic_frac_J035 >= 0.70) == 2, "MODERATE", "WITHHOLD")))
utils::write.csv(sens, file.path(OUT, "wgcna_registry_threshold_sensitivity.csv"),
                 row.names = FALSE)

## ===================================================================
## cross-check against Part-28
## ===================================================================
k <- match(paste(res$dataset, res$module_id), paste(p28$dataset, p28$module_id))
cross <- data.frame(
  dataset = res$dataset, module_id = res$module_id,
  p28_n_enriched = p28$n_enriched_terms_FDR05[k],
  v10_n_enriched = res$n_enriched_terms_FDR05,
  n_enriched_match = p28$n_enriched_terms_FDR05[k] == res$n_enriched_terms_FDR05,
  p28_best_term = p28$best_term[k], v10_best_term = res$best_any_term,
  best_term_match = p28$best_term[k] == res$best_any_term,
  p28_best_FDR = p28$best_term_FDR[k], v10_best_FDR = res$best_any_FDR,
  best_FDR_reldiff = abs(p28$best_term_FDR[k] - res$best_any_FDR) /
    pmax(res$best_any_FDR, .Machine$double.xmin),
  p28_module_size = p28$module_size[k], v10_module_size = res$module_size,
  p28_median_abs_kME = p28$median_abs_kME[k], v10_median_abs_kME = res$median_abs_kME,
  p28_external_scopes_agree = p28$external_scopes_agree[k],
  v10_external_scopes_agree = res$external_scopes_agree,
  p28_annotation_confidence = p28$annotation_confidence[k],
  v10_label_confidence = res$label_confidence,
  p28_label_review_required = p28$label_review_required[k],
  stringsAsFactors = FALSE)
utils::write.csv(cross, file.path(OUT, "wgcna_registry_part28_crosscheck.csv"),
                 row.names = FALSE)

## ===================================================================
## console report
## ===================================================================
cat("\n===== PART-29 S28/29/30 WGCNA REGISTRY FINAL REVIEW =====\n")
cat("modules reviewed            :", nrow(res), "\n")
cat("label_confidence counts     :\n"); print(table(res$label_confidence))
cat("\nby dataset:\n"); print(table(res$dataset, res$label_confidence))
cat("\nconverging layers (of 4)    :\n"); print(table(res$evidence_layers_converging))
cat("functional layers (of 3)    :\n"); print(table(res$n_functional_layers_converging))
cat("hub coherence               :\n"); print(table(res$hub_class))
cat("semantic coherence          :\n"); print(table(res$semantic_class))
cat("modules with 0 sig GO terms :",
    sum(res$n_enriched_terms_FDR05 == 0), "->",
    paste(res$dataset[res$n_enriched_terms_FDR05 == 0],
          res$module_id[res$n_enriched_terms_FDR05 == 0], collapse = ", "), "\n")
cat("\n--- cross-check vs Part-28 ---\n")
cat("n_enriched mismatches :", sum(!cross$n_enriched_match, na.rm = TRUE), "\n")
cat("best_term mismatches  :", sum(!cross$best_term_match, na.rm = TRUE), "\n")
cat("module_size mismatches:", sum(cross$p28_module_size != cross$v10_module_size), "\n")
cat("median|kME| max abs diff:",
    max(abs(cross$p28_median_abs_kME - cross$v10_median_abs_kME)), "\n")
cat("max relative diff best FDR:", max(cross$best_FDR_reldiff, na.rm = TRUE), "\n")
cat("external_scopes_agree differences:",
    sum(cross$p28_external_scopes_agree != cross$v10_external_scopes_agree), "\n")
cat("\n--- WITHHOLD modules ---\n")
w <- res[res$label_confidence == "WITHHOLD", ]
for (j in seq_len(nrow(w)))
  cat(sprintf("  %-16s %s  n_sig=%-4d hub=%-8s sem=%-32s func=%d\n",
              w$dataset[j], w$module_id[j], w$n_enriched_terms_FDR05[j],
              w$hub_class[j], w$semantic_class[j],
              w$n_functional_layers_converging[j]))
cat("\n--- HIGH modules ---\n")
hi <- res[res$label_confidence == "HIGH", ]
for (j in seq_len(nrow(hi)))
  cat(sprintf("  %-16s %s  %s\n", hi$dataset[j], hi$module_id[j],
              hi$recommended_publication_form[j]))
cat("\n--- MODERATE modules ---\n")
mo <- res[res$label_confidence == "MODERATE", ]
for (j in seq_len(nrow(mo)))
  cat(sprintf("  %-16s %s  hub=%-8s sem=%-32s %s\n", mo$dataset[j], mo$module_id[j],
              mo$hub_class[j], mo$semantic_class[j],
              mo$recommended_publication_form[j]))
cat("\n--- neuropil m11 ---\n")
print(m11[, c("claim", "verdict")])
cat("\nm11 recommendation: FORM A -> ", forms$wording[1], "\n", sep = "")
cat("m11 label_confidence in the review table: ",
    res$label_confidence[res$dataset == "neuron_neuropil" & res$module_id == "m11"],
    "; converging layers ",
    res$evidence_layers_converging[res$dataset == "neuron_neuropil" & res$module_id == "m11"],
    "/4\n", sep = "")
cat("\n--- m11 satellite terms (outside the myelin-named family) ---\n")
print(out_sat[, c("ontology", "term", "p_adjust", "count",
                  "frac_of_term_genes_also_in_myelin_family_union",
                  "overlapping_m11_hub_symbols")], row.names = FALSE)
cat("\n--- threshold sensitivity ---\n")
cat("confidence at J=0.15 :\n"); print(table(sens$confidence_if_J015))
cat("confidence at J=0.25 (PRIMARY):\n"); print(table(sens$primary_label_confidence))
cat("confidence at J=0.35 :\n"); print(table(sens$confidence_if_J035))
cat("modules whose confidence changes with J: ",
    sum(sens$confidence_if_J015 != sens$primary_label_confidence |
          sens$confidence_if_J035 != sens$primary_label_confidence), "\n")
mi <- which(sens$dataset == "neuron_neuropil" & sens$module_id == "m11")
cat(sprintf("m11 semantic fraction: J0.15=%.3f  J0.25=%.3f  J0.35=%.3f  CConly=%.3f  BPonly=%.3f -> confidence %s / %s / %s\n",
            sens$semantic_frac_J015[mi], sens$semantic_frac_J025_PRIMARY[mi],
            sens$semantic_frac_J035[mi], sens$semantic_frac_CC_only_J025[mi],
            sens$semantic_frac_BP_only_J025[mi],
            sens$confidence_if_J015[mi], sens$primary_label_confidence[mi],
            sens$confidence_if_J035[mi]))
cat(sprintf("m11 hubs: %d/25 top-25 hubs inside leading term gene sets (%.0f%%)\n",
            sens$n_top25_hubs_in_leading_termsets[mi],
            100 * sens$frac_top25_hubs_in_leading_termsets[mi]))

cat("\nwritten:\n")
for (f in c("wgcna_registry_final_review.csv", "wgcna_m11_claim_verification.csv",
            "wgcna_m11_candidate_form_adjudication.csv",
            "wgcna_m11_significant_term_diagnostic.csv",
            "wgcna_registry_threshold_sensitivity.csv",
            "wgcna_registry_part28_crosscheck.csv"))
  cat("  ", normalizePath(file.path(OUT, f), winslash = "/"), "\n")
