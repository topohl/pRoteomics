#!/usr/bin/env Rscript

# Re-audit of what the atlas still omits, against the SEVEN-row atlas.
#
# This is not a repeat of the earlier six-row check. Neuron projection
# development is now a displayed row, so the omitted set is genuinely different
# and the question is asked again from scratch.
#
# The criteria are the ones already fixed for this project and are not retuned:
# semantic coherence, biological distinctness, recurrence, not a parent/child
# ladder, at least as specific as the displayed terms, not QC context, and
# material size. Clustering is phenotype-blind.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({ library(GOSemSim); library(GO.db); library(AnnotationDbi) })

OUT <- file.path("results", "tables", "publication_audits", "program_evidence")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05
CUTOFF <- 0.30
MIN_TERMS <- 3L
MIN_UNITS <- 3L
MIN_OCC <- 100L
LADDER_MAX <- 0.50
QC_RE <- "immunoglobulin|keratin|cornified|epiderm|lymphocyte|natural killer"

TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)
sup <- TH[is.finite(TH$GSEA_FDR) & TH$GSEA_FDR < ALPHA, , drop = FALSE]
displayed <- unique(TH$GO_ID[TH$theme_claim_eligible %in% TRUE])
omit <- sort(setdiff(unique(sup$GO_ID), displayed))
ont <- AnnotationDbi::select(GO.db::GO.db, keys = omit,
                             columns = c("ONTOLOGY", "TERM"), keytype = "GOID")
omit <- ont$GOID[!is.na(ont$ONTOLOGY) & ont$ONTOLOGY == "BP"]
term <- stats::setNames(ont$TERM, ont$GOID)

sd_bp <- GOSemSim::godata(annoDb = "org.Mm.eg.db", ont = "BP",
                          computeIC = FALSE)
S <- as.matrix(GOSemSim::mgoSim(omit, omit, semData = sd_bp, measure = "Wang",
                                combine = NULL))
S[!is.finite(S)] <- 0; diag(S) <- 1; S <- S[omit, omit, drop = FALSE]
blk <- stats::cutree(stats::hclust(stats::as.dist(1 - S), method = "average"),
                     h = 1 - CUTOFF)

Sd <- as.matrix(GOSemSim::mgoSim(omit, intersect(displayed, rownames(S)),
                                 semData = sd_bp, measure = "Wang",
                                 combine = NULL))
if (!is.matrix(Sd)) Sd <- matrix(Sd, nrow = length(omit))
Sd[!is.finite(Sd)] <- 0
disp_all <- AnnotationDbi::select(GO.db::GO.db,
  keys = intersect(displayed, AnnotationDbi::keys(GO.db::GO.db, "GOID")),
  columns = "ONTOLOGY", keytype = "GOID")
disp_bp <- disp_all$GOID[disp_all$ONTOLOGY == "BP"]
Sd2 <- as.matrix(GOSemSim::mgoSim(omit, disp_bp, semData = sd_bp,
                                  measure = "Wang", combine = NULL))
Sd2[!is.finite(Sd2)] <- 0

ANC <- as.list(GO.db::GOBPANCESTOR)
anc_of <- function(id) { a <- ANC[[id]]; if (is.null(a)) character(0) else a }
depth_of <- function(id) { a <- anc_of(id); if (!length(a)) NA_integer_ else sum(a != "all") }
DISP_DEPTH <- stats::median(vapply(disp_bp, depth_of, integer(1)), na.rm = TRUE)

rows <- do.call(rbind, lapply(sort(unique(blk)), function(k) {
  g <- names(blk)[blk == k]
  z <- sup[sup$GO_ID %in% g, , drop = FALSE]
  W <- S[g, g, drop = FALSE]
  internal_anc <- g %in% unique(unlist(lapply(g, anc_of), use.names = FALSE))
  occ <- table(factor(z$GO_ID, levels = g))
  med <- if (length(g) == 1L) g else {
    a <- vapply(g, function(x) mean(W[x, setdiff(g, x)]), numeric(1))
    g[order(-a, g)][1]
  }
  data.frame(
    cluster_id = sprintf("OMIT%02d", k),
    medoid_GO_ID = med, medoid_term = unname(term[med]),
    n_terms = length(g), n_supported_occurrences = nrow(z),
    n_datasets = length(unique(z$dataset)),
    n_spatial_units = length(unique(z$spatial_unit)),
    median_within_similarity = if (length(g) > 1L)
      stats::median(W[upper.tri(W)]) else NA_real_,
    frac_terms_near_displayed = mean(apply(Sd2[g, , drop = FALSE], 1, max) >= 0.50),
    frac_occ_from_internal_ancestors = if (sum(occ)) sum(occ[internal_anc]) / sum(occ) else NA_real_,
    medoid_depth = depth_of(med), displayed_median_depth = DISP_DEPTH,
    frac_qc_terms = mean(grepl(QC_RE, unname(term[g]), ignore.case = TRUE)),
    example_terms = paste(utils::head(unname(term[g[order(g)]]), 8), collapse = "; "),
    stringsAsFactors = FALSE)
}))

rows$coherent <- is.na(rows$median_within_similarity) |
  rows$median_within_similarity >= CUTOFF
rows$recurrent <- rows$n_spatial_units >= MIN_UNITS & rows$n_datasets >= 2L
rows$substantial <- rows$n_terms >= MIN_TERMS &
  rows$n_supported_occurrences >= MIN_OCC
rows$distinct <- rows$frac_terms_near_displayed < 0.50
rows$not_a_ladder <- rows$frac_occ_from_internal_ancestors < LADDER_MAX
rows$specific_enough <- !is.na(rows$medoid_depth) &
  rows$medoid_depth >= rows$displayed_median_depth
rows$qc_like <- grepl(QC_RE, rows$medoid_term, ignore.case = TRUE) |
  rows$frac_qc_terms >= 0.30

rows$rejection_class <- with(rows, ifelse(
  qc_like, "QC_OR_NON_NEURAL",
  ifelse(!coherent, "TOO_DIFFUSE",
  ifelse(!recurrent, "TOO_RARE",
  ifelse(!distinct, "ALREADY_COVERED",
  ifelse(!not_a_ladder | !specific_enough, "ANCESTOR_REDUNDANCY",
  ifelse(!substantial, "TOO_SMALL", "GENUINE_MISSING_CANDIDATE")))))))
rows$is_true_omission <- rows$rejection_class == "GENUINE_MISSING_CANDIDATE"
rows$rejection_reason <- with(rows, sprintf(
  paste0("%.0f%% of terms already within a displayed row; %.0f%% of occurrences ",
         "from internal ancestors; medoid depth %s vs %s displayed; %d terms, ",
         "%d occurrences, %d units, %d datasets; within-similarity %.3f"),
  100 * frac_terms_near_displayed, 100 * frac_occ_from_internal_ancestors,
  medoid_depth, displayed_median_depth, n_terms, n_supported_occurrences,
  n_spatial_units, n_datasets, median_within_similarity))

rows <- rows[order(-rows$n_supported_occurrences), , drop = FALSE]
rows$rank <- seq_len(nrow(rows))
utils::write.csv(rows, file.path(OUT, "omitted_program_recheck.csv"),
                 row.names = FALSE)

cat("\n===== OMITTED-PROGRAM RE-AUDIT (seven-row atlas) =====\n")
cat("supported terms outside the seven rows:", length(omit),
    "| clusters:", nrow(rows), "\n")
cat("GENUINE_MISSING_CANDIDATE:", sum(rows$is_true_omission), "\n\n")
print(utils::head(rows[, c("rank", "cluster_id", "n_terms",
  "n_supported_occurrences", "n_spatial_units", "frac_terms_near_displayed",
  "frac_occ_from_internal_ancestors", "medoid_term", "rejection_class")], 10),
  row.names = FALSE)
cat("\nwritten to:", OUT, "\n")
