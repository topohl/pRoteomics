#!/usr/bin/env Rscript

# Part-29 sections 11-15: can a phenotype-blind registry recover a clean
# mitochondrial theme, and how much do the six themes overlap each other?
#
# AUDIT ONLY. Nothing canonical is rerun or modified. The candidate registry is
# written to the isolated audit layer and is NOT activated.
#
# THE V3 RULE IS ONE ONTOLOGY OPERATION, NOT A BLACKLIST. A term leaves the
# mitochondrial theme if and only if it lies in the glycolysis sub-DAG, i.e. it
# IS GO:0006096 (glycolytic process) or has GO:0006096 among its BP ancestors.
# No term is removed for being inconvenient, non-significant or wrongly signed,
# and no NES or FDR is consulted while the membership is decided.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({ library(GO.db); library(AnnotationDbi) })

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
REP <- file.path("results", "reports", "publication_audits",
                 "upstream_enrichment_v10")
CFG <- file.path("config")
dir.create(AUD, recursive = TRUE, showWarnings = FALSE)
dir.create(REP, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
GLYCOLYSIS_ROOT <- "GO:0006096"
MITO <- "mitochondrial_respiration_oxphos"

TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)
elig <- TH[TH$theme_claim_eligible %in% TRUE, , drop = FALSE]
map <- unique(elig[, c("GO_ID", "GO_description", "theme_id")])
themes <- sort(unique(map$theme_id))
members <- lapply(themes, function(t) sort(unique(map$GO_ID[map$theme_id == t])))
names(members) <- themes

# =============================================== S11 the two candidate registries
anc <- as.list(GO.db::GOBPANCESTOR)
in_glycolysis <- function(id)
  identical(id, GLYCOLYSIS_ROOT) ||
  (!is.null(anc[[id]]) && GLYCOLYSIS_ROOT %in% anc[[id]])

V2 <- members[[MITO]]
is_gly <- vapply(V2, in_glycolysis, logical(1))
V3 <- V2[!is_gly]

term_of <- stats::setNames(map$GO_description[!duplicated(map$GO_ID)],
                           map$GO_ID[!duplicated(map$GO_ID)])

cand <- data.frame(
  theme_id = MITO, GO_ID = V2, GO_description = unname(term_of[V2]),
  in_V2_broad = TRUE, in_V3_mito_specific = !is_gly,
  removal_rule = ifelse(is_gly,
    paste0("lies in the ", GLYCOLYSIS_ROOT, " glycolytic-process sub-DAG"), ""),
  stringsAsFactors = FALSE)
cand <- cand[order(cand$GO_description), , drop = FALSE]
utils::write.csv(cand, file.path(AUD,
  "mitochondrial_registry_v3_candidate_membership.csv"), row.names = FALSE)

# ============================================ S12 what the swap does to the atlas
cells <- function(ids, label) {
  z <- elig[elig$theme_id == MITO & elig$GO_ID %in% ids, , drop = FALSE]
  z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID)), ,
         drop = FALSE]
  k <- paste(z$dataset, z$spatial_unit, z$contrast, sep = "\r")
  do.call(rbind, lapply(split(seq_len(nrow(z)), k), function(ix) {
    w <- z[ix, , drop = FALSE]
    sup <- is.finite(w$GSEA_FDR) & w$GSEA_FDR < ALPHA
    data.frame(dataset = w$dataset[1], spatial_unit = w$spatial_unit[1],
               contrast = w$contrast[1], registry = label, n_terms = nrow(w),
               median_NES = stats::median(w$NES, na.rm = TRUE),
               n_supported = sum(sup), dot = sum(sup) > 0,
               stringsAsFactors = FALSE)
  }))
}
c2 <- cells(V2, "V2_BROAD"); c3 <- cells(V3, "V3_MITO_SPECIFIC")
key <- function(d) paste(d$dataset, d$spatial_unit, d$contrast, sep = "\r")
m <- merge(c2, c3, by = c("dataset", "spatial_unit", "contrast"),
           suffixes = c("_V2", "_V3"))
m$difference <- m$median_NES_V2 - m$median_NES_V3
m$sign_changed <- sign(m$median_NES_V2) != sign(m$median_NES_V3)
m$support_dot_changed <- m$dot_V2 != m$dot_V3
m$n_terms_lost <- m$n_terms_V2 - m$n_terms_V3
m$supported_occurrences_lost <- m$n_supported_V2 - m$n_supported_V3
m <- m[order(m$contrast, m$dataset, m$spatial_unit), , drop = FALSE]
utils::write.csv(m, file.path(AUD,
  "mitochondrial_registry_v2_v3_comparison.csv"), row.names = FALSE)

# ================================================== S14 theme overlap matrix
#
# Gene-level overlap uses the pinned mouse BP annotation, so two themes can be
# non-overlapping in GO IDs and still share most of their genes.
g2go <- AnnotationDbi::mapIds(
  org.Mm.eg.db::org.Mm.eg.db, keys = unique(map$GO_ID), keytype = "GOALL",
  column = "SYMBOL", multiVals = "list")
genes_of <- function(ids) sort(unique(unlist(g2go[ids], use.names = FALSE)))

supported_pairs <- unique(paste(elig$GO_ID, elig$dataset, elig$spatial_unit,
                                elig$contrast)[is.finite(elig$GSEA_FDR) &
                                                 elig$GSEA_FDR < ALPHA])
ov <- do.call(rbind, lapply(seq_along(themes), function(i)
  do.call(rbind, lapply(seq_along(themes), function(j) {
    if (j <= i) return(NULL)
    a <- members[[themes[i]]]; b <- members[[themes[j]]]
    ga <- genes_of(a); gb <- genes_of(b)
    shared <- intersect(a, b)
    sup_shared <- sum(vapply(strsplit(supported_pairs, " "), function(p)
      p[1] %in% shared, logical(1)))
    data.frame(
      theme_a = themes[i], theme_b = themes[j],
      n_GO_a = length(a), n_GO_b = length(b),
      n_shared_GO = length(shared),
      jaccard_GO = length(shared) / length(union(a, b)),
      n_genes_a = length(ga), n_genes_b = length(gb),
      n_shared_genes = length(intersect(ga, gb)),
      jaccard_genes = length(intersect(ga, gb)) / length(union(ga, gb)),
      n_shared_supported_occurrences = sup_shared,
      shared_GO_terms = paste(unname(term_of[shared]), collapse = "; "),
      stringsAsFactors = FALSE)
  }))))
ov <- ov[order(-ov$jaccard_genes), , drop = FALSE]
utils::write.csv(ov, file.path(AUD, "atlas_theme_overlap_matrix.csv"),
                 row.names = FALSE)

# ============================================== S15 overlap sensitivity
#
# Descriptive only: recompute each theme cell using ONLY the GO terms unique to
# that theme, to see whether shared terms are what make two rows look alike.
shared_any <- unique(unlist(lapply(seq_along(themes), function(i)
  unlist(lapply(setdiff(seq_along(themes), i), function(j)
    intersect(members[[themes[i]]], members[[themes[j]]]))))))

theme_cells <- function(th, ids, label) {
  z <- elig[elig$theme_id == th & elig$GO_ID %in% ids, , drop = FALSE]
  if (!nrow(z)) return(NULL)
  z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID)), ,
         drop = FALSE]
  k <- paste(z$dataset, z$spatial_unit, z$contrast, sep = "\r")
  do.call(rbind, lapply(split(seq_len(nrow(z)), k), function(ix) {
    w <- z[ix, , drop = FALSE]
    sup <- is.finite(w$GSEA_FDR) & w$GSEA_FDR < ALPHA
    data.frame(theme = th, dataset = w$dataset[1],
               spatial_unit = w$spatial_unit[1], contrast = w$contrast[1],
               variant = label, n_terms = nrow(w),
               median_NES = stats::median(w$NES, na.rm = TRUE),
               dot = sum(sup) > 0, stringsAsFactors = FALSE)
  }))
}
sens <- do.call(rbind, lapply(themes, function(th) {
  full <- theme_cells(th, members[[th]], "all_terms")
  uniq <- theme_cells(th, setdiff(members[[th]], shared_any), "unique_terms_only")
  if (is.null(full) || is.null(uniq)) return(NULL)
  z <- merge(full, uniq, by = c("theme", "dataset", "spatial_unit", "contrast"),
             suffixes = c("_all", "_unique"))
  z$difference <- z$median_NES_all - z$median_NES_unique
  z$sign_changed <- sign(z$median_NES_all) != sign(z$median_NES_unique)
  z$support_dot_changed <- z$dot_all != z$dot_unique
  z
}))
sens <- sens[order(sens$theme, sens$contrast, sens$dataset, sens$spatial_unit), ,
             drop = FALSE]
utils::write.csv(sens, file.path(AUD, "atlas_theme_overlap_sensitivity.csv"),
                 row.names = FALSE)

cat("\n===== PART-29 REGISTRY AND OVERLAP =====\n")
cat("mitochondrial V2:", length(V2), "terms | V3 (glycolysis sub-DAG removed):",
    length(V3), "| removed:", sum(is_gly), "\n")
cat("  removed terms:", paste(unname(term_of[V2[is_gly]]), collapse = "; "), "\n")
cat("  cells:", nrow(m), "| sign changes:", sum(m$sign_changed),
    "| support-dot changes:", sum(m$support_dot_changed),
    "| median |diff|:", sprintf("%.4f", stats::median(abs(m$difference))),
    "| max |diff|:", sprintf("%.4f", max(abs(m$difference))), "\n")
ca1 <- m[m$dataset == "microglia" & m$spatial_unit == "CA1", ]
cat("  CA1 microglia:", paste(sprintf("%s %.3f->%.3f", ca1$contrast,
    ca1$median_NES_V2, ca1$median_NES_V3), collapse = " | "), "\n")
cat("  supported occurrences lost overall:", sum(m$supported_occurrences_lost),
    "\n")
cat("\ntheme overlap, top 4 by gene Jaccard:\n")
print(utils::head(ov[, c("theme_a", "theme_b", "n_shared_GO", "jaccard_GO",
                         "n_shared_genes", "jaccard_genes")], 4),
      row.names = FALSE)
cat("\noverlap sensitivity: cells", nrow(sens), "| sign changes",
    sum(sens$sign_changed), "| dot changes", sum(sens$support_dot_changed),
    "| median |diff|", sprintf("%.4f", stats::median(abs(sens$difference))), "\n")
cat("\nwritten to:", AUD, "\n")
