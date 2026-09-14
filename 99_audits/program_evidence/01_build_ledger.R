#!/usr/bin/env Rscript

# Program-evidence audit of the curated spatial GO-program atlas.
#
# USER-REQUESTED VERIFICATION, not a new discovery pass. Nothing canonical is
# rerun or modified: differential abundance, GSEA, CAMERA and WGCNA are all read
# as-is. This layer only reconstructs program membership from the registry,
# joins it to the canonical enrichment results, and summarises the leading-edge
# proteins that already drive those results.
#
# WHAT A READER MUST BE ABLE TO TRACE BY HAND:
#   program -> registry anchor -> GO term -> supported GO term -> leading-edge
#   protein
# Every step below is one join on an identifier, so any row of the ledger can be
# checked against the canonical table with a text search.
#
# PHENOTYPE-BLIND BY CONSTRUCTION. Membership comes from the ontology registry.
# Semantic clusters come from GO similarity alone. NES, FDR and direction are
# attached afterwards to DESCRIBE a term, never to decide whether it belongs.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({ library(GOSemSim); library(GO.db); library(AnnotationDbi) })

OUT <- file.path("results", "tables", "publication_audits", "program_evidence")
REP <- file.path("results", "reports", "publication_audits", "program_evidence")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(REP, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05
EPS_FLOOR <- 1e-10        # clusterProfiler::gseGO default, never overridden
SEM_CUTOFF <- 0.30        # the cutoff this project already uses for program-level
                          # semantic grouping; ~4x the median similarity of an
                          # arbitrary BP pair, and not retuned here
OFF_THEME_SIM <- 0.30     # a term is off-theme if it is no closer than that to
                          # ANY of its own program's registry anchors

TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)
REG <- utils::read.delim(file.path("config", "manuscript_go_theme_registry.tsv"),
                         stringsAsFactors = FALSE, quote = "")
stopifnot(unique(REG$registry_version) == "manuscript_go_themes_v3")

# the seven displayed rows, in registry display order
LABEL <- c(rna_processing_splicing_rnp = "RNA processing",
           ribosome_translation = "Translation / ribosome",
           chromatin_organization = "Chromatin",
           mitochondrial_respiration_oxphos = "Mitochondrial respiration",
           synaptic_signaling_vesicle = "Synaptic signalling / vesicle",
           neuron_projection_development = "Neuron projection development",
           autophagy_lysosome_endosome = "Autophagy / endolysosomal")
PROGS <- names(LABEL)

prim <- TH[TH$theme_claim_eligible %in% TRUE & TH$theme_id %in% PROGS, ,
           drop = FALSE]
stopifnot(nrow(prim) > 0)
prim$supported <- is.finite(prim$GSEA_FDR) & prim$GSEA_FDR < ALPHA
prim$at_eps_floor <- is.finite(prim$raw_p) & prim$raw_p <= EPS_FLOOR

members <- lapply(PROGS, function(p) sort(unique(prim$GO_ID[prim$theme_id == p])))
names(members) <- PROGS
all_ids <- sort(unique(unlist(members)))
term_of <- stats::setNames(prim$GO_description[!duplicated(prim$GO_ID)],
                           prim$GO_ID[!duplicated(prim$GO_ID)])

# ------------------------------------------------- semantic structure, per program
sd_bp <- GOSemSim::godata(annoDb = "org.Mm.eg.db", ont = "BP",
                          computeIC = FALSE)
S <- as.matrix(GOSemSim::mgoSim(all_ids, all_ids, semData = sd_bp,
                                measure = "Wang", combine = NULL))
S[!is.finite(S)] <- 0; diag(S) <- 1; S <- S[all_ids, all_ids, drop = FALSE]

ANC <- as.list(GO.db::GOBPANCESTOR)
anc_of <- function(id) { a <- ANC[[id]]; if (is.null(a)) character(0) else a }

clust <- list(); medoid <- list()
for (p in PROGS) {
  g <- members[[p]]
  W <- S[g, g, drop = FALSE]
  cl <- if (length(g) == 1L) stats::setNames(1L, g) else
    stats::cutree(stats::hclust(stats::as.dist(1 - W), method = "average"),
                  h = 1 - SEM_CUTOFF)
  clust[[p]] <- cl
  medoid[[p]] <- stats::setNames(vapply(sort(unique(cl)), function(k) {
    gg <- names(cl)[cl == k]
    if (length(gg) == 1L) return(gg)
    a <- vapply(gg, function(x) mean(W[x, setdiff(gg, x)]), numeric(1))
    gg[order(-a, gg)][1]
  }, character(1)), sort(unique(cl)))
}

# similarity of every term to its own program's registry anchors
anchor_sim <- function(p, id) {
  a <- intersect(REG$anchor_go_id[REG$theme_id == p &
    REG$match_scope != "exclude_anchor_and_descendants"], rownames(S))
  if (!length(a)) return(NA_real_)
  max(S[id, a])
}

# ---------------------------------------------------------- the master ledger
prim$membership_class <- ifelse(prim$match_type == "exact_anchor",
                                "DEFINITION_ANCHOR", "CLOSURE_DESCENDANT")
prim$program_label <- unname(LABEL[prim$theme_id])

# per (program, term) descriptors
key <- paste(prim$theme_id, prim$GO_ID, sep = "|")
uk <- !duplicated(key)
pt <- prim[uk, c("theme_id", "GO_ID"), drop = FALSE]
pt$semantic_cluster <- vapply(seq_len(nrow(pt)), function(i)
  sprintf("%s_C%02d", pt$theme_id[i], clust[[pt$theme_id[i]]][[pt$GO_ID[i]]]),
  character(1))
pt$semantic_cluster_medoid <- vapply(seq_len(nrow(pt)), function(i) {
  k <- clust[[pt$theme_id[i]]][[pt$GO_ID[i]]]
  unname(term_of[medoid[[pt$theme_id[i]]][[as.character(k)]]])
}, character(1))
pt$similarity_to_own_anchor <- vapply(seq_len(nrow(pt)), function(i)
  anchor_sim(pt$theme_id[i], pt$GO_ID[i]), numeric(1))
pt$off_theme_flag <- is.finite(pt$similarity_to_own_anchor) &
  pt$similarity_to_own_anchor < OFF_THEME_SIM
# is this term an ancestor of another term in the SAME program?
pt$redundant_ancestor_flag <- vapply(seq_len(nrow(pt)), function(i) {
  g <- setdiff(members[[pt$theme_id[i]]], pt$GO_ID[i])
  pt$GO_ID[i] %in% unique(unlist(lapply(g, anc_of), use.names = FALSE))
}, logical(1))
# does it also belong to another displayed program?
pt$shared_with_programs <- vapply(seq_len(nrow(pt)), function(i) {
  o <- setdiff(PROGS[vapply(PROGS, function(p) pt$GO_ID[i] %in% members[[p]],
                            logical(1))], pt$theme_id[i])
  if (length(o)) paste(o, collapse = "; ") else ""
}, character(1))

# per-term recurrence over supported contexts only
sup <- prim[prim$supported, , drop = FALSE]
rec <- stats::aggregate(
  list(n_supported_contexts = sup$source_comparison),
  by = list(theme_id = sup$theme_id, GO_ID = sup$GO_ID),
  FUN = function(z) length(unique(z)))
rec2 <- stats::aggregate(
  list(n_supported_units = sup$spatial_unit),
  by = list(theme_id = sup$theme_id, GO_ID = sup$GO_ID),
  FUN = function(z) length(unique(z)))
rec3 <- stats::aggregate(
  list(n_supported_datasets = sup$dataset),
  by = list(theme_id = sup$theme_id, GO_ID = sup$GO_ID),
  FUN = function(z) length(unique(z)))
pt <- merge(pt, rec, by = c("theme_id", "GO_ID"), all.x = TRUE)
pt <- merge(pt, rec2, by = c("theme_id", "GO_ID"), all.x = TRUE)
pt <- merge(pt, rec3, by = c("theme_id", "GO_ID"), all.x = TRUE)
for (k in c("n_supported_contexts", "n_supported_units", "n_supported_datasets"))
  pt[[k]][is.na(pt[[k]])] <- 0L

led <- merge(prim, pt, by = c("theme_id", "GO_ID"), all.x = TRUE)
led$why_in_program <- sprintf("%s (%s) via %s, path length %s",
  led$anchor_GO_ID, led$anchor_label, led$match_type,
  led$shortest_path_length)
led <- led[, c("theme_id", "program_label", "GO_ID", "GO_description",
  "why_in_program", "membership_class", "anchor_GO_ID", "anchor_label",
  "shortest_path_length", "source_comparison", "dataset", "spatial_unit",
  "contrast", "NES", "raw_p", "GSEA_FDR", "supported", "at_eps_floor",
  "leading_edge_genes", "leading_edge_proteins", "n_supported_contexts",
  "n_supported_units", "n_supported_datasets", "semantic_cluster",
  "semantic_cluster_medoid", "similarity_to_own_anchor", "off_theme_flag",
  "redundant_ancestor_flag", "shared_with_programs")]
led <- led[order(match(led$theme_id, PROGS), led$GO_ID, led$source_comparison), ,
           drop = FALSE]
utils::write.csv(led, file.path(OUT, "program_evidence_ledger.csv"),
                 row.names = FALSE)

# =========================================================== protein evidence
#
# Counted ONLY over FDR-supported (term x context) rows. A gene that appears in
# one supported term in one context is not evidence that the program is about
# that gene; a gene recurring across many supported terms AND many contexts is.
split_genes <- function(x) {
  z <- trimws(unlist(strsplit(as.character(x), "[;,/]")))
  z[nzchar(z)]
}
prot <- do.call(rbind, lapply(PROGS, function(p) {
  z <- sup[sup$theme_id == p, , drop = FALSE]
  if (!nrow(z)) return(NULL)
  rows <- do.call(rbind, lapply(seq_len(nrow(z)), function(i) {
    g <- split_genes(z$leading_edge_genes[i])
    if (!length(g)) return(NULL)
    data.frame(theme_id = p, gene = g, GO_ID = z$GO_ID[i],
               comparison = z$source_comparison[i], dataset = z$dataset[i],
               spatial_unit = z$spatial_unit[i], stringsAsFactors = FALSE)
  }))
  if (is.null(rows)) return(NULL)
  agg <- do.call(rbind, lapply(split(seq_len(nrow(rows)), rows$gene), function(ix) {
    r <- rows[ix, , drop = FALSE]
    data.frame(theme_id = p, program_label = unname(LABEL[p]),
               gene = r$gene[1],
               n_supported_terms = length(unique(r$GO_ID)),
               n_supported_contexts = length(unique(r$comparison)),
               n_datasets = length(unique(r$dataset)),
               n_spatial_units = length(unique(r$spatial_unit)),
               example_terms = paste(utils::head(
                 unname(term_of[sort(unique(r$GO_ID))]), 3), collapse = "; "),
               stringsAsFactors = FALSE)
  }))
  agg$n_program_supported_terms <- length(unique(z$GO_ID))
  agg$frac_supported_terms <- agg$n_supported_terms / agg$n_program_supported_terms
  # a gene counts as CORE only if it recurs across several supported terms AND
  # several contexts - the distinction the brief asks for between a recurrent
  # protein and one carried in by a single broad annotation
  agg$evidence_class <- ifelse(
    agg$n_supported_terms >= 3L & agg$n_supported_contexts >= 3L, "RECURRENT_CORE",
    ifelse(agg$n_supported_terms == 1L & agg$n_supported_contexts == 1L,
           "SINGLE_APPEARANCE", "INTERMEDIATE"))
  agg[order(-agg$n_supported_terms, -agg$n_supported_contexts, agg$gene), ,
      drop = FALSE]
}))
utils::write.csv(prot, file.path(OUT, "program_leading_edge_protein_evidence.csv"),
                 row.names = FALSE)

saveRDS(list(S = S, members = members, clust = clust, medoid = medoid,
             LABEL = LABEL, PROGS = PROGS, term_of = term_of, pt = pt,
             prim = prim, sup = sup, prot = prot, REG = REG),
        file.path(OUT, "program_evidence_state.rds"))

cat("\n===== PROGRAM EVIDENCE LEDGER =====\n")
cat("ledger rows:", nrow(led), "| programs:", length(PROGS),
    "| GO terms:", length(all_ids), "\n")
cat("supported (term x context) rows:", sum(led$supported), "\n")
cat("leading-edge gene records:", nrow(prot), "\n\n")
s <- do.call(rbind, lapply(PROGS, function(p) {
  q <- pt[pt$theme_id == p, ]
  z <- prot[prot$theme_id == p, ]
  data.frame(program = unname(LABEL[p]), terms = nrow(q),
             supported_terms = sum(q$n_supported_contexts > 0),
             clusters = length(unique(q$semantic_cluster)),
             off_theme = sum(q$off_theme_flag),
             redundant_ancestors = sum(q$redundant_ancestor_flag),
             shared = sum(nzchar(q$shared_with_programs)),
             core_proteins = sum(z$evidence_class == "RECURRENT_CORE"),
             stringsAsFactors = FALSE)
}))
print(s, row.names = FALSE)
cat("\nwritten to:", OUT, "\n")
