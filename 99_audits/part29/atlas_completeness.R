#!/usr/bin/env Rscript

# Part-29 sections 7-10: is the six-row atlas complete enough to be called an
# atlas, and is any omitted supported program a coherent recurrent biology?
#
# AUDIT ONLY. Nothing canonical is rerun or modified. Every number is a count or
# a semantic-similarity statistic over GO IDs the canonical workflow already
# tested. All clustering is PHENOTYPE-BLIND: NES, FDR, dataset and direction are
# never inputs to a cluster, only descriptions of one afterwards.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({
  library(GOSemSim); library(GO.db); library(simplifyEnrichment)
})

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
REP <- file.path("results", "reports", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(AUD, recursive = TRUE, showWarnings = FALSE)
dir.create(REP, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
SEM_CUTOFF <- 0.70          # the repository's own pinned Wang cutoff
MIN_RECURRENCE_UNITS <- 3L  # a cluster is "recurrent" only across >=3 units
MIN_CLUSTER_TERMS <- 3L     # and only if it holds >=3 distinct GO terms

TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)

PRIMARY <- sort(unique(TH$theme_id[TH$theme_claim_eligible %in% TRUE &
                                     nzchar(TH$theme_id)]))
stopifnot(length(PRIMARY) == 6L)

# ------------------------------------------------------- S7 coverage classes
sup <- TH[is.finite(TH$GSEA_FDR) & TH$GSEA_FDR < ALPHA, , drop = FALSE]
sup$class <- with(sup, ifelse(
  theme_claim_eligible %in% TRUE & theme_id %in% PRIMARY,
  ifelse(assignment_status == "multi_theme", "MULTI_THEME", "PRIMARY_THEME"),
  ifelse(theme_role == "qc_review" | assignment_status == "qc_review",
         "QC_REVIEW_THEME",
         ifelse(nzchar(theme_id), "SUPPORTING_THEME", "UNCLASSIFIED"))))
sup$direction <- ifelse(sup$NES > 0, "positive", "negative")
sup$in_primary <- sup$class %in% c("PRIMARY_THEME", "MULTI_THEME")

cov_block <- function(d, label) {
  if (!nrow(d)) return(NULL)
  ug <- unique(d$GO_ID)
  ugp <- unique(d$GO_ID[d$in_primary])
  data.frame(
    scope = label,
    n_supported_occurrences = nrow(d),
    n_occurrences_primary = sum(d$in_primary),
    occurrence_coverage = sum(d$in_primary) / nrow(d),
    n_unique_supported_GO = length(ug),
    n_unique_GO_primary = length(ugp),
    unique_term_coverage = length(ugp) / length(ug),
    n_PRIMARY_THEME = sum(d$class == "PRIMARY_THEME"),
    n_MULTI_THEME = sum(d$class == "MULTI_THEME"),
    n_SUPPORTING_THEME = sum(d$class == "SUPPORTING_THEME"),
    n_QC_REVIEW_THEME = sum(d$class == "QC_REVIEW_THEME"),
    n_UNCLASSIFIED = sum(d$class == "UNCLASSIFIED"),
    stringsAsFactors = FALSE)
}

cov <- rbind(
  cov_block(sup, "ALL_CONTRASTS"),
  do.call(rbind, lapply(sort(unique(sup$contrast)), function(ct)
    cov_block(sup[sup$contrast == ct, ], paste0("contrast=", ct)))),
  do.call(rbind, lapply(sort(unique(sup$dataset)), function(ds)
    cov_block(sup[sup$dataset == ds & sup$contrast == "SUS - RES", ],
              paste0("SUS-RES dataset=", ds)))),
  do.call(rbind, lapply(sort(unique(sup$spatial_unit[sup$contrast == "SUS - RES"])),
    function(u) cov_block(sup[sup$spatial_unit == u & sup$contrast == "SUS - RES", ],
                          paste0("SUS-RES unit=", u)))),
  do.call(rbind, lapply(c("positive", "negative"), function(dr)
    cov_block(sup[sup$direction == dr & sup$contrast == "SUS - RES", ],
              paste0("SUS-RES direction=", dr)))))
utils::write.csv(cov, file.path(AUD, "atlas_supported_term_coverage.csv"),
                 row.names = FALSE)

# ------------------------------ S8 semantic structure of the OMITTED supported
#
# The omitted set is every FDR-supported GO term that no primary theme carries.
# It is clustered on ontology alone. NES, FDR, dataset and direction are
# attached afterwards purely to describe a cluster, never to form one.
omit_ids <- sort(setdiff(unique(sup$GO_ID),
                         unique(sup$GO_ID[sup$in_primary])))
ont <- AnnotationDbi::select(GO.db::GO.db, keys = omit_ids,
                             columns = c("ONTOLOGY", "TERM"), keytype = "GOID")
omit_ids <- ont$GOID[!is.na(ont$ONTOLOGY) & ont$ONTOLOGY == "BP"]
term_of <- stats::setNames(ont$TERM, ont$GOID)

sd_bp <- GOSemSim::godata(annoDb = "org.Mm.eg.db", ont = "BP",
                          computeIC = FALSE)
S <- as.matrix(GOSemSim::mgoSim(omit_ids, omit_ids, semData = sd_bp,
                                measure = "Wang", combine = NULL))
S[!is.finite(S)] <- 0
diag(S) <- 1
S <- S[omit_ids, omit_ids, drop = FALSE]

hc <- stats::hclust(stats::as.dist(1 - S), method = "average")
cl <- stats::cutree(hc, h = 1 - SEM_CUTOFF)
# binary_cut is recorded as a cross-check but is NOT used to form the clusters:
# on 934 broadly-connected BP terms it returns a single block, which is a fact
# about the ontology graph rather than a usable grouping. The clusters below
# come from average linkage at Wang >= 0.30 - roughly four times the median
# off-diagonal similarity of this very set (0.079), so a cluster is far above
# the relatedness any two BP terms share by graph structure alone, and the
# level at which the set resolves into 67 groups with only 15 singletons.
set.seed(20260914L)
BINARY_CUT_BLOCKS <- length(unique(simplifyEnrichment::binary_cut(S)))
CLUSTER_CUTOFF <- 0.30
blk <- stats::cutree(hc, h = 1 - CLUSTER_CUTOFF)

medoid <- function(g) {
  if (length(g) == 1L) return(g)
  a <- vapply(g, function(x) mean(S[x, setdiff(g, x)]), numeric(1))
  g[order(-a, g)][1]
}

om <- sup[!sup$in_primary & sup$GO_ID %in% omit_ids, , drop = FALSE]
clusters <- do.call(rbind, lapply(sort(unique(blk)), function(k) {
  g <- names(blk)[blk == k]
  z <- om[om$GO_ID %in% g, , drop = FALSE]
  med <- medoid(g)
  W <- S[g, g, drop = FALSE]
  data.frame(
    semantic_cluster = sprintf("OMIT_B%02d", k),
    n_unique_GO_terms = length(g),
    n_FDR_supported_occurrences = nrow(z),
    datasets_represented = paste(sort(unique(z$dataset)), collapse = "; "),
    n_datasets = length(unique(z$dataset)),
    spatial_units_represented = paste(sort(unique(z$spatial_unit)),
                                      collapse = "; "),
    n_spatial_units = length(unique(z$spatial_unit)),
    contrasts_represented = paste(sort(unique(z$contrast)), collapse = "; "),
    sign_directions = paste(sort(unique(ifelse(z$NES > 0, "+", "-"))),
                            collapse = "/"),
    semantic_medoid = med, semantic_medoid_term = unname(term_of[med]),
    median_within_similarity = if (length(g) > 1L)
      stats::median(W[upper.tri(W)]) else NA_real_,
    representative_GO_terms = paste(utils::head(
      unname(term_of[g[order(g)]]), 8), collapse = "; "),
    n_duplicate_clusters = length(unique(cl[g])),
    stringsAsFactors = FALSE)
}))

# S10 candidate-row rule: ALL of coherence, distinctness, recurrence, not a
# parent/child duplicate of an existing row, not a QC artefact
prim_ids <- sort(unique(TH$GO_ID[TH$theme_claim_eligible %in% TRUE]))
Sp <- as.matrix(GOSemSim::mgoSim(omit_ids, prim_ids, semData = sd_bp,
                                 measure = "Wang", combine = NULL))
Sp[!is.finite(Sp)] <- 0
clusters$max_similarity_to_a_primary_term <- vapply(
  sort(unique(blk)), function(k) max(Sp[names(blk)[blk == k], , drop = FALSE]),
  numeric(1))
# a single close term does not make a whole cluster redundant; what matters is
# what FRACTION of the cluster an existing row already covers
clusters$fraction_terms_near_a_primary_term <- vapply(
  sort(unique(blk)), function(k) {
    g <- names(blk)[blk == k]
    mean(apply(Sp[g, , drop = FALSE], 1, max) >= 0.50)
  }, numeric(1))

# SPECIFICITY. A row must name a biological PROGRAM, not a generic regulatory
# or metabolic parent. Depth in the BP DAG is an objective, phenotype-blind
# proxy, and the bar is set by the atlas itself: a candidate medoid must be at
# least as deep as the median term already displayed, so a new row would be no
# vaguer than the six it would join.
anc <- as.list(GO.db::GOBPANCESTOR)
go_depth <- function(id) {
  a <- anc[[id]]
  if (is.null(a)) return(NA_integer_)
  sum(a != "all")
}
PRIMARY_DEPTH <- stats::median(vapply(prim_ids, go_depth, integer(1)),
                               na.rm = TRUE)
clusters$medoid_GO_depth <- vapply(clusters$semantic_medoid, go_depth,
                                   integer(1))
clusters$primary_theme_median_depth <- PRIMARY_DEPTH
clusters$specific_enough <- !is.na(clusters$medoid_GO_depth) &
  clusters$medoid_GO_depth >= PRIMARY_DEPTH
qc_ids <- unique(TH$GO_ID[TH$theme_role == "qc_review"])
clusters$overlaps_QC_theme <- vapply(sort(unique(blk)), function(k)
  any(names(blk)[blk == k] %in% qc_ids), logical(1))

clusters$coherent <- is.na(clusters$median_within_similarity) |
  clusters$median_within_similarity >= 0.30
clusters$recurrent <- clusters$n_spatial_units >= MIN_RECURRENCE_UNITS
clusters$substantial <- clusters$n_unique_GO_terms >= MIN_CLUSTER_TERMS
clusters$distinct_from_primary <-
  clusters$fraction_terms_near_a_primary_term < 0.50
clusters$verdict <- ifelse(
  clusters$substantial & clusters$coherent & clusters$recurrent &
    clusters$distinct_from_primary & clusters$specific_enough &
    !clusters$overlaps_QC_theme,
  "ADD", "DO_NOT_ADD")
clusters$verdict_reason <- with(clusters, paste0(
  ifelse(substantial, "", "too few terms; "),
  ifelse(coherent, "", "not semantically coherent; "),
  ifelse(recurrent, "", "not recurrent across >=3 spatial units; "),
  ifelse(distinct_from_primary, "",
         "mostly already covered by an existing row; "),
  ifelse(specific_enough, "",
         "medoid is a generic GO parent, shallower than the displayed terms; "),
  ifelse(overlaps_QC_theme, "overlaps a QC-review theme; ", "")))
clusters$verdict_reason[clusters$verdict == "ADD"] <-
  "meets every candidate-row criterion"
clusters <- clusters[order(-clusters$n_FDR_supported_occurrences), ,
                     drop = FALSE]
utils::write.csv(clusters, file.path(AUD, "atlas_omitted_semantic_clusters.csv"),
                 row.names = FALSE)
utils::write.csv(
  clusters[, c("semantic_cluster", "n_unique_GO_terms",
               "n_FDR_supported_occurrences", "n_spatial_units",
               "semantic_medoid_term", "median_within_similarity",
               "fraction_terms_near_a_primary_term", "medoid_GO_depth",
               "primary_theme_median_depth", "verdict", "verdict_reason")],
  file.path(AUD, "atlas_candidate_new_rows.csv"), row.names = FALSE)

a <- cov[cov$scope == "ALL_CONTRASTS", ]
s <- cov[cov$scope == "contrast=SUS - RES", ]
cat("\n===== PART-29 ATLAS COMPLETENESS =====\n")
cat("FDR-supported occurrences (all contrasts):", a$n_supported_occurrences,
    "| in a primary theme:", a$n_occurrences_primary,
    sprintf("(%.1f%%)\n", 100 * a$occurrence_coverage))
cat("unique supported GO IDs:", a$n_unique_supported_GO, "| in a primary theme:",
    a$n_unique_GO_primary, sprintf("(%.1f%%)\n", 100 * a$unique_term_coverage))
cat("SUS-RES occurrence coverage:", sprintf("%.1f%%", 100 * s$occurrence_coverage),
    "| unique-term coverage:", sprintf("%.1f%%\n", 100 * s$unique_term_coverage))
cat("classes:", paste(sprintf("%s=%d", c("PRIMARY", "MULTI", "SUPPORTING",
    "QC", "UNCLASSIFIED"), c(a$n_PRIMARY_THEME, a$n_MULTI_THEME,
    a$n_SUPPORTING_THEME, a$n_QC_REVIEW_THEME, a$n_UNCLASSIFIED)),
    collapse = "  "), "\n")
cat("\nomitted supported BP terms:", length(omit_ids), "in",
    nrow(clusters), "semantic clusters\n")
cat("candidate rows ADD:", sum(clusters$verdict == "ADD"), "\n")
print(utils::head(clusters[, c("semantic_cluster", "n_unique_GO_terms",
  "n_FDR_supported_occurrences", "n_spatial_units", "semantic_medoid_term",
  "verdict")], 12), row.names = FALSE)
cat("\nwritten to:", AUD, "\n")
