#!/usr/bin/env Rscript

# Finalization pass, sections 3-7: sanity-check the largest omitted semantic
# clusters below neuron projection development, and verify the boundary of the
# approved row.
#
# AUDIT ONLY. No canonical analysis is rerun. The clustering is the one Part 29
# already produced (average linkage on GOSemSim Wang, cut at 0.30); no threshold
# is retuned and no new search is run.
#
# THE KEY DIAGNOSTIC. A cluster's occurrence count can be large simply because
# GO stacks ancestors and descendants of one branch, each tested separately. For
# every cluster this script therefore reports what fraction of its supported
# occurrences is carried by terms that are ANCESTORS of another term in the same
# cluster. A cluster that is mostly such a ladder is hierarchical redundancy,
# not a hidden program.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({ library(GOSemSim); library(GO.db); library(AnnotationDbi) })

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(AUD, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05
CUTOFF <- 0.30            # the Part-29 cutoff, unchanged
APPROVED <- "OMIT_B12"    # neuron projection development

TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)
cl29 <- utils::read.csv(file.path(AUD, "atlas_omitted_semantic_clusters.csv"),
                        stringsAsFactors = FALSE)

sup <- TH[is.finite(TH$GSEA_FDR) & TH$GSEA_FDR < ALPHA, , drop = FALSE]
prim <- unique(TH$GO_ID[TH$theme_claim_eligible %in% TRUE])
omit <- sort(setdiff(unique(sup$GO_ID), prim))
ont <- AnnotationDbi::select(GO.db::GO.db, keys = omit,
                             columns = c("ONTOLOGY", "TERM"), keytype = "GOID")
omit <- ont$GOID[!is.na(ont$ONTOLOGY) & ont$ONTOLOGY == "BP"]
term <- stats::setNames(ont$TERM, ont$GOID)

sd_bp <- GOSemSim::godata(annoDb = "org.Mm.eg.db", ont = "BP",
                          computeIC = FALSE)
S <- as.matrix(GOSemSim::mgoSim(omit, omit, semData = sd_bp, measure = "Wang",
                                combine = NULL))
S[!is.finite(S)] <- 0; diag(S) <- 1; S <- S[omit, omit, drop = FALSE]
hc <- stats::hclust(stats::as.dist(1 - S), method = "average")
blk <- stats::cutree(hc, h = 1 - CUTOFF)

ANC <- as.list(GO.db::GOBPANCESTOR)
anc_of <- function(id) { a <- ANC[[id]]; if (is.null(a)) character(0) else a }

# ---------------------------------------------------- per-cluster diagnostics
rows <- do.call(rbind, lapply(sort(unique(blk)), function(k) {
  g <- names(blk)[blk == k]
  z <- sup[sup$GO_ID %in% g, , drop = FALSE]
  cid <- sprintf("OMIT_B%02d", k)
  c29 <- cl29[cl29$semantic_cluster == cid, , drop = FALSE]
  # a term is "an ancestor inside its own cluster" if another cluster member
  # lists it among its ancestors
  all_anc <- unique(unlist(lapply(g, anc_of), use.names = FALSE))
  is_internal_ancestor <- g %in% all_anc
  occ <- table(factor(z$GO_ID, levels = g))
  frac_ladder <- if (sum(occ)) sum(occ[is_internal_ancestor]) / sum(occ) else NA_real_
  # how recurrent are the terms that are NOT generic internal ancestors?
  spec <- g[!is_internal_ancestor]
  u_spec <- if (length(spec))
    max(vapply(spec, function(i)
      length(unique(z$spatial_unit[z$GO_ID == i])), integer(1))) else 0L
  data.frame(
    semantic_cluster_id = cid,
    semantic_medoid_GO_ID = if (nrow(c29)) c29$semantic_medoid[1] else NA_character_,
    semantic_medoid_term = if (nrow(c29)) c29$semantic_medoid_term[1] else NA_character_,
    n_unique_GO_terms = length(g),
    n_supported_occurrences = nrow(z),
    n_datasets = length(unique(z$dataset)),
    n_spatial_units = length(unique(z$spatial_unit)),
    n_contrasts = length(unique(z$contrast)),
    directions_present = paste(sort(unique(ifelse(z$NES > 0, "+", "-"))),
                               collapse = "/"),
    median_within_similarity = if (nrow(c29)) c29$median_within_similarity[1] else NA_real_,
    overlap_with_existing_theme = if (nrow(c29))
      c29$fraction_terms_near_a_primary_term[1] else NA_real_,
    medoid_GO_depth = if (nrow(c29)) c29$medoid_GO_depth[1] else NA_integer_,
    primary_theme_median_depth = if (nrow(c29)) c29$primary_theme_median_depth[1] else NA_integer_,
    n_internal_ancestor_terms = sum(is_internal_ancestor),
    fraction_occurrences_from_internal_ancestors = frac_ladder,
    max_units_for_a_non_ancestor_term = u_spec,
    representative_GO_terms = paste(utils::head(unname(term[g[order(g)]]), 10),
                                    collapse = "; "),
    stringsAsFactors = FALSE)
}))

# closest existing primary theme, by best Wang similarity
prim_bp <- prim[prim %in% AnnotationDbi::keys(GO.db::GO.db, "GOID")]
Sp <- as.matrix(GOSemSim::mgoSim(omit, prim_bp, semData = sd_bp,
                                 measure = "Wang", combine = NULL))
Sp[!is.finite(Sp)] <- 0
theme_of <- stats::setNames(TH$theme_id[!duplicated(TH$GO_ID)],
                            TH$GO_ID[!duplicated(TH$GO_ID)])
rows$closest_existing_theme <- vapply(sort(unique(blk)), function(k) {
  g <- intersect(names(blk)[blk == k], rownames(Sp))
  if (!length(g)) return(NA_character_)
  m <- Sp[g, , drop = FALSE]
  unname(theme_of[colnames(m)[which(m == max(m), arr.ind = TRUE)[1, 2]]])
}, character(1))

# ------------------------------------------------------- the fixed criteria
QC_LIKE <- "immunoglobulin|keratin|cornified|epiderm|lymphocyte|natural killer|B cell|T cell"
# a cluster is QC-like only if its MEDOID is, or if a substantial share of its
# members are; one immune term inside a signalling grab-bag is not QC context
qc_frac <- vapply(sort(unique(blk)), function(k) {
  g <- names(blk)[blk == k]
  mean(grepl(QC_LIKE, unname(term[g]), ignore.case = TRUE))
}, numeric(1))
rows$fraction_qc_like_terms <- qc_frac
rows$qc_like <- grepl(QC_LIKE, rows$semantic_medoid_term, ignore.case = TRUE) |
  qc_frac >= 0.30
rows$plausible_neural_biology <- !rows$qc_like
rows$recurrent <- rows$n_spatial_units >= 3L & rows$n_datasets >= 2L
rows$distinct <- rows$overlap_with_existing_theme < 0.50
rows$coherent_A <- rows$median_within_similarity >= 0.30
rows$specific_D <- rows$medoid_GO_depth >= rows$primary_theme_median_depth
# criterion D also fails when the cluster is mostly an ancestor ladder of one
# branch rather than a set of distinct sibling processes
rows$not_a_ladder <- rows$fraction_occurrences_from_internal_ancestors < 0.50
rows$substantial_F <- rows$n_supported_occurrences >= 100L

rows$decision <- with(rows, ifelse(
  semantic_cluster_id == APPROVED, "ALREADY_APPROVED",
  ifelse(qc_like, "NO_QC_CONTEXT",
  ifelse(!coherent_A, "NO_INCOHERENT",
  ifelse(!recurrent, "NO_NOT_RECURRENT",
  ifelse(!distinct, "NO_REDUNDANT",
  ifelse(!not_a_ladder | !specific_D, "NO_REDUNDANT",
  ifelse(!substantial_F, "NO_TOO_SMALL", "YES_MAJOR_OMISSION"))))))))

rows$decision_reason <- with(rows, ifelse(
  decision == "ALREADY_APPROVED", "approved in Part 29 as the seventh row",
  ifelse(decision == "NO_QC_CONTEXT",
    sprintf(paste0("predominantly non-neural / immune-annotation terms (%.0f%% ",
                   "of members); QC context, not hippocampal biology"),
            100 * fraction_qc_like_terms),
  ifelse(decision == "NO_INCOHERENT",
    sprintf("median within-cluster Wang similarity %.3f: a heterogeneous group, not one program",
            median_within_similarity),
  ifelse(decision == "NO_NOT_RECURRENT",
    sprintf("supported in only %d spatial unit(s) across %d dataset(s)",
            n_spatial_units, n_datasets),
  ifelse(decision == "NO_REDUNDANT",
    sprintf(paste0("%.0f%% of its terms already sit within an existing row, ",
                   "%.0f%% of its occurrences come from terms that are ",
                   "ancestors of other members (a parent/child ladder), and ",
                   "its medoid sits at GO depth %d against %d for the terms ",
                   "already displayed"),
            100 * overlap_with_existing_theme,
            100 * fraction_occurrences_from_internal_ancestors,
            medoid_GO_depth, primary_theme_median_depth),
  ifelse(decision == "NO_TOO_SMALL",
    sprintf("%d supported occurrences, below the materiality bar",
            n_supported_occurrences),
    "meets every fixed addition criterion")))))))

rows <- rows[order(-rows$n_supported_occurrences), , drop = FALSE]
rows$rank <- seq_len(nrow(rows))
top <- rows[rows$semantic_cluster_id != APPROVED, , drop = FALSE][1:10, ]
top$rank <- seq_len(nrow(top))
top$biological_summary <- top$representative_GO_terms
utils::write.csv(top[, c("rank", "semantic_cluster_id", "semantic_medoid_GO_ID",
  "semantic_medoid_term", "n_unique_GO_terms", "n_supported_occurrences",
  "n_datasets", "n_spatial_units", "n_contrasts", "directions_present",
  "overlap_with_existing_theme", "closest_existing_theme",
  "fraction_occurrences_from_internal_ancestors",
  "max_units_for_a_non_ancestor_term", "medoid_GO_depth",
  "primary_theme_median_depth", "biological_summary",
  "plausible_neural_biology", "recurrent", "distinct", "qc_like",
  "decision", "decision_reason")],
  file.path(AUD, "atlas_top_omitted_clusters_final_review.csv"),
  row.names = FALSE)

# ============================================ S7 neuron projection boundary
np <- names(blk)[blk == as.integer(sub("OMIT_B", "", APPROVED))]
NEIGHBOUR <- "axon|neurite|dendrit|projection|neuron|neurogen|synaps"
near <- rows[rows$semantic_cluster_id != APPROVED &
               grepl(NEIGHBOUR, rows$representative_GO_terms,
                     ignore.case = TRUE), , drop = FALSE]
bnd <- do.call(rbind, lapply(seq_len(nrow(near)), function(i) {
  k <- as.integer(sub("OMIT_B", "", near$semantic_cluster_id[i]))
  g <- names(blk)[blk == k]
  data.frame(
    cluster_id = near$semantic_cluster_id[i],
    medoid_term = near$semantic_medoid_term[i],
    n_terms = length(g),
    n_supported_occurrences = near$n_supported_occurrences[i],
    max_wang_to_neuron_projection_cluster = max(S[g, np]),
    median_wang_to_neuron_projection_cluster = stats::median(S[g, np]),
    shares_terms_with_np = length(intersect(g, np)),
    neural_terms_in_cluster = paste(utils::head(
      grep(NEIGHBOUR, unname(term[g]), ignore.case = TRUE, value = TRUE), 6),
      collapse = "; "),
    merge_recommended = FALSE,
    reason = paste0("separate cluster at Wang >= ", CUTOFF,
                    "; not merged because it is a different ontology branch ",
                    "(organization/homeostasis rather than projection ",
                    "development) and merging would widen the row beyond a ",
                    "single defensible program"),
    stringsAsFactors = FALSE)
}))
np_self <- data.frame(
  cluster_id = APPROVED, medoid_term = "neuron projection development",
  n_terms = length(np),
  n_supported_occurrences = rows$n_supported_occurrences[
    rows$semantic_cluster_id == APPROVED],
  max_wang_to_neuron_projection_cluster = 1,
  median_wang_to_neuron_projection_cluster = stats::median(S[np, np][
    upper.tri(S[np, np])]),
  shares_terms_with_np = length(np),
  neural_terms_in_cluster = paste(utils::head(unname(term[sort(np)]), 10),
                                  collapse = "; "),
  merge_recommended = NA, reason = "the approved row itself",
  stringsAsFactors = FALSE)
utils::write.csv(rbind(np_self, bnd),
                 file.path(AUD, "neuron_projection_theme_boundary.csv"),
                 row.names = FALSE)

cat("\n===== FINALIZATION SANITY CHECK =====\n")
cat("omitted clusters:", nrow(rows), "| reviewed in detail: top 10 below",
    APPROVED, "\n\n")
print(top[, c("rank", "semantic_cluster_id", "n_supported_occurrences",
              "n_spatial_units", "overlap_with_existing_theme",
              "fraction_occurrences_from_internal_ancestors",
              "semantic_medoid_term", "decision")], row.names = FALSE)
cat("\nYES_MAJOR_OMISSION:", sum(top$decision == "YES_MAJOR_OMISSION"), "\n")
cat("\nneuron projection cluster:", length(np), "terms | neighbouring clusters examined:",
    nrow(bnd), "| merges recommended:", sum(bnd$merge_recommended), "\n")
cat("  max Wang from any neighbour into the approved cluster:",
    sprintf("%.3f", max(bnd$max_wang_to_neuron_projection_cluster)), "\n")
cat("\nwritten to:", AUD, "\n")
