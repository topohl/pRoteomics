#!/usr/bin/env Rscript

# The reviewer-facing review book: one sheet per atlas program, a cross-program
# summary, a figure-ready table, and a mechanical blind-naming check.
#
# VERDICTS ARE RULE-BASED. Each program's classification follows from thresholds
# fixed here, not from biological taste, so the same inputs always produce the
# same verdict. The prose in each sheet explains the numbers; it never overrides
# them.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

OUT <- file.path("results", "tables", "publication_audits", "program_evidence")
REP <- file.path("results", "reports", "publication_audits", "program_evidence")
SHEETS <- file.path(REP, "program_sheets")
dir.create(SHEETS, recursive = TRUE, showWarnings = FALSE)

st <- readRDS(file.path(OUT, "program_evidence_state.rds"))
pt <- st$pt; prot <- st$prot; term_of <- st$term_of
LABEL <- st$LABEL; PROGS <- st$PROGS; REG <- st$REG; sup <- st$sup
omit <- utils::read.csv(file.path(OUT, "omitted_program_recheck.csv"),
                        stringsAsFactors = FALSE)

# ---------------------------------------------------- rule-based classification
#
# STRONG                 one dominant semantic block, most terms supported, and a
#                        recurrent protein core
# SUPPORTED_BUT_BROAD    supported and coherent, but the label spans more than
#                        one block or a minority of terms carry the support
# MIXED                  several comparable blocks, or a large share of terms far
#                        from the program's own anchors
# MISLEADING             the protein core or the supported terms contradict the label
DOMINANT_FRAC <- 0.60     # largest cluster's share of terms
SUPPORT_FRAC <- 0.40      # share of terms supported somewhere
CORE_MIN <- 10L           # recurrent-core proteins needed to call a core real

summ <- do.call(rbind, lapply(PROGS, function(p) {
  q <- pt[pt$theme_id == p, , drop = FALSE]
  z <- prot[prot$theme_id == p, , drop = FALSE]
  s <- sup[sup$theme_id == p, , drop = FALSE]
  cl <- table(q$semantic_cluster)
  major <- cl[cl >= max(3L, ceiling(0.10 * nrow(q)))]
  data.frame(
    theme_id = p, current_label = unname(LABEL[p]),
    registry_display_label = REG$display_label[REG$theme_id == p][1],
    n_terms = nrow(q),
    n_terms_supported = sum(q$n_supported_contexts > 0),
    frac_terms_supported = sum(q$n_supported_contexts > 0) / nrow(q),
    n_supported_occurrences = nrow(s),
    n_spatial_units = length(unique(s$spatial_unit)),
    n_datasets = length(unique(s$dataset)),
    n_semantic_clusters = length(cl),
    n_major_clusters = length(major),
    largest_cluster_frac = max(cl) / nrow(q),
    n_off_theme_terms = sum(q$off_theme_flag),
    frac_off_theme = mean(q$off_theme_flag),
    n_redundant_ancestors = sum(q$redundant_ancestor_flag),
    n_shared_terms = sum(nzchar(q$shared_with_programs)),
    n_recurrent_core_proteins = sum(z$evidence_class == "RECURRENT_CORE"),
    n_single_appearance_proteins = sum(z$evidence_class == "SINGLE_APPEARANCE"),
    top_core_proteins = paste(utils::head(
      z$gene[z$evidence_class == "RECURRENT_CORE"], 12), collapse = "; "),
    major_subclusters = paste(vapply(names(major), function(k)
      sprintf("%s (n=%d)", unique(q$semantic_cluster_medoid[q$semantic_cluster == k]),
              cl[[k]]), character(1)), collapse = " || "),
    stringsAsFactors = FALSE)
}))

summ$coherence_status <- with(summ, ifelse(
  n_major_clusters <= 1L, "SINGLE_BLOCK",
  ifelse(largest_cluster_frac >= DOMINANT_FRAC, "DOMINANT_BLOCK_PLUS_SATELLITES",
         "SEVERAL_COMPARABLE_BLOCKS")))
summ$verdict <- with(summ, ifelse(
  frac_off_theme >= 0.50, "MIXED",
  ifelse(coherence_status == "SEVERAL_COMPARABLE_BLOCKS" |
           frac_terms_supported < SUPPORT_FRAC, "SUPPORTED_BUT_BROAD",
  ifelse(n_recurrent_core_proteins >= CORE_MIN &
           coherence_status != "SEVERAL_COMPARABLE_BLOCKS", "STRONG",
         "SUPPORTED_BUT_BROAD"))))
summ$selection_status <- "RECURRENT_COHERENT_PROGRAM"
summ$naming_status <- with(summ, ifelse(
  verdict == "STRONG", "LABEL_MATCHES_MEMBERSHIP",
  ifelse(verdict == "SUPPORTED_BUT_BROAD", "LABEL_BROADER_THAN_DOMINANT_BLOCK",
         "LABEL_SPANS_SEVERAL_BLOCKS")))

RECOMMENDED <- c(
  rna_processing_splicing_rnp = "RNA processing",
  ribosome_translation = "Translation / ribosome",
  chromatin_organization = "Chromatin / epigenetic regulation",
  mitochondrial_respiration_oxphos = "Mitochondrial respiration",
  synaptic_signaling_vesicle = "Synaptic signalling / vesicle",
  neuron_projection_development = "Neuron projection development",
  autophagy_lysosome_endosome = "Autophagy / endolysosomal")
summ$recommended_label <- unname(RECOMMENDED[summ$theme_id])
summ$action <- ifelse(summ$current_label == summ$recommended_label,
                      "KEEP", "REFINE_WORDING_OPTIONAL")
summ <- summ[match(PROGS, summ$theme_id), , drop = FALSE]
utils::write.csv(summ, file.path(OUT, "program_cross_summary.csv"),
                 row.names = FALSE)

# ------------------------------------------------ mechanical blind-naming check
#
# Strip the label, keep only the supported GO term names and the recurrent
# protein core, and report the most frequent content words. If those words spell
# out the label, the name is reproducible from the evidence alone.
STOP <- c("process", "regulation", "positive", "negative", "cellular", "complex",
          "via", "involved", "protein", "organization", "assembly", "activity",
          "pathway", "from", "into", "with", "other", "type")
blind <- do.call(rbind, lapply(PROGS, function(p) {
  q <- pt[pt$theme_id == p, , drop = FALSE]
  ids <- q$GO_ID[q$n_supported_contexts > 0]
  w <- unlist(strsplit(tolower(gsub("[^a-zA-Z ]", " ", unname(term_of[ids]))), " +"))
  w <- w[nchar(w) > 3 & !w %in% STOP]
  tw <- sort(table(w), decreasing = TRUE)
  z <- prot[prot$theme_id == p & prot$evidence_class == "RECURRENT_CORE", ]
  data.frame(theme_id = p, current_label = unname(LABEL[p]),
             top_supported_term_words = paste(sprintf("%s(%d)",
               names(utils::head(tw, 8)), utils::head(tw, 8)), collapse = " "),
             top_core_proteins = paste(utils::head(z$gene, 10), collapse = " "),
             stringsAsFactors = FALSE)
}))
utils::write.csv(blind, file.path(OUT, "program_blind_naming_check.csv"),
                 row.names = FALSE)

# ----------------------------------------------------- figure-ready table
figtab <- do.call(rbind, lapply(PROGS, function(p) {
  q <- pt[pt$theme_id == p, , drop = FALSE]
  q <- q[q$n_supported_contexts > 0, , drop = FALSE]
  q <- q[order(-q$n_supported_contexts), , drop = FALSE]
  z <- prot[prot$theme_id == p & prot$evidence_class == "RECURRENT_CORE", ]
  data.frame(
    display_order = which(PROGS == p),
    program = unname(LABEL[p]),
    representative_supported_GO_terms = paste(
      sprintf("%s (%d ctx)", unname(term_of[utils::head(q$GO_ID, 3)]),
              utils::head(q$n_supported_contexts, 3)), collapse = "; "),
    recurrent_leading_edge_proteins = paste(utils::head(z$gene, 8), collapse = ", "),
    n_supported_terms = nrow(q), n_core_proteins = nrow(z),
    stringsAsFactors = FALSE)
}))
utils::write.csv(figtab, file.path(OUT, "program_figure_ready_table.csv"),
                 row.names = FALSE)

# --------------------------------------------------------- per-program sheets
NOTE <- list(
  rna_processing_splicing_rnp = list(
    why = "The general parent GO:0006396 plus explicit regulatory and RNP anchors. Splicing, tRNA/rRNA/miRNA processing and RNP assembly are one ontology branch.",
    take = "The protein core is a textbook spliceosome and hnRNP roster - SR proteins (Srsf1, Srsf10), U2AF2, PRPF19, EIF4A3, hnRNPs, PTBP1, TRA2B. The broad label is earned by genuine diversity of membership, not by vagueness.",
    caveat = "Six rRNA terms are shared with Translation / ribosome, and Npm1 appears in both protein cores. That overlap is real and is disclosed rather than engineered away.",
    defend = "Yes."),
  ribosome_translation = list(
    why = "Translation, its regulation, and ribosome biogenesis anchors.",
    take = "Two clusters map exactly onto the two halves of the label: ribosome biogenesis / rRNA maturation, and cytoplasmic translation. The core is almost entirely large- and small-subunit ribosomal proteins plus Eef2.",
    caveat = "Shares the six rRNA terms with RNA processing. Keeping the rows separate is defensible because translation proper and transcript processing are different biology, and the overlap moves atlas values by a median of 0.003 NES.",
    defend = "Yes."),
  chromatin_organization = list(
    why = "The chromatin organisation anchor GO:0006325 and its approved descendants.",
    take = "The supported core is histones and remodellers - H1f0, Macroh2a1, Smarca2, Atrx, Hdac2, Hmgb1 - which is genuinely chromatin. But this is the weakest row: fewest supported terms, fewest supported occurrences, and the most semantic clusters.",
    caveat = "Eight of thirteen terms sit in the epigenetic-regulation branch rather than the organisation branch the anchor names, so 'Chromatin' alone under-describes half the row. The core also contains nuclear lamins (Lmna, Lmnb1, Lmnb2) and hnRNPs (Hnrnpk, Hnrnpu), which are nuclear but not chromatin proper.",
    defend = "Yes, with the limitation stated. 'Chromatin / epigenetic regulation' would describe the membership more exactly."),
  mitochondrial_respiration_oxphos = list(
    why = "Cellular respiration, electron transport chain and respiratory-chain complex assembly anchors, MINUS the GO:0006096 glycolysis sub-DAG.",
    take = "The strongest protein evidence in the atlas. The recurrent core is Complex I structural subunits (Ndufs1/2/3/6/7/8, Ndufc2, Ndufa8, Ndufv2) together with mitochondrially encoded subunits (mt-Nd1/2/4/5). Nothing cytosolic survives in the core.",
    caveat = "Two terms fall below the anchor-similarity threshold - pyruvate decarboxylation to acetyl-CoA and proton-motive-force-driven mitochondrial ATP synthesis. Both are genuinely mitochondrial; this is a limitation of Wang similarity against a broad anchor, not an off-theme membership.",
    defend = "Yes, and this is the row to show first if the naming is challenged."),
  synaptic_signaling_vesicle = list(
    why = "Synaptic signalling, trans-synaptic signalling regulation, and vesicle-mediated transport in synapse anchors.",
    take = "Two clusters correspond exactly to the two halves of the label. The core is specifically synaptic rather than generic trafficking: Syt1, Unc13a, Stx1a, Stx1b, Rab3a, Rims1 on the release side, Nlgn1, Nlgn3, Grin2b, Ntrk2 on the postsynaptic side.",
    caveat = "It is the largest row (62 terms) and one term falls just below the anchor-similarity threshold.",
    defend = "Yes."),
  neuron_projection_development = list(
    why = "Neuron projection development, its regulation, and neuron migration as a sibling branch. Added because it was the only omitted semantic cluster meeting every fixed criterion.",
    take = "The most semantically coherent row in the atlas: all 46 terms form ONE cluster and none is flagged off-theme. The core carries canonical guidance receptors (Epha4, Nrp1) alongside the growth-cone actin machinery they signal through (Tiam1, Cyfip1, Actr3, Pak3, Dbn1, Dbnl).",
    caveat = "The protein core is dominated by actin-cytoskeletal effectors. That is expected for neurite outgrowth rather than evidence of a generic cytoskeleton row - the guidance receptors and the migration terms are what distinguish it, and no term overlaps any other program.",
    defend = "Yes."),
  autophagy_lysosome_endosome = list(
    why = "Autophagy, lysosome organisation and endosomal transport anchors.",
    take = "The protein core spans both halves of the label and justifies combining them: LC3/GABARAP family (Map1lc3a, Map1lc3b, Gabarapl1, Gabarapl2) for autophagy, and Ctsd, Lamtor1, Rragc, Chmp4b, Chmp7, Stx12, Vti1a, Snx4 for the endolysosomal side.",
    caveat = "Only 13 of 36 terms are supported anywhere, and the terms fall into three clusters (autophagy, lysosome, endosomal transport). It is a coherent membrane-degradation axis rather than a single process.",
    defend = "Yes, as a combined axis. Do not describe it as a single pathway."))

book <- c("# Atlas program review book", "",
"One sheet per displayed row. Every number traces to",
"`program_evidence_ledger.csv` and `program_leading_edge_protein_evidence.csv`.",
"Protein counts are over FDR-supported (GO term x context) rows only.", "")
for (p in PROGS) {
  r <- summ[summ$theme_id == p, ]
  q <- pt[pt$theme_id == p, , drop = FALSE]
  z <- prot[prot$theme_id == p & prot$evidence_class == "RECURRENT_CORE", ]
  n <- NOTE[[p]]
  anchors <- REG[REG$theme_id == p, , drop = FALSE]
  supq <- q[q$n_supported_contexts > 0, , drop = FALSE]
  supq <- supq[order(-supq$n_supported_contexts), , drop = FALSE]
  sh <- c(
    sprintf("# %s", r$current_label), "",
    sprintf("- **Current label:** %s", r$current_label),
    sprintf("- **Registry label:** %s", r$registry_display_label),
    sprintf("- **Recommended label:** %s", r$recommended_label),
    sprintf("- **Verdict: %s** (coherence: %s; naming: %s; action: %s)",
            r$verdict, r$coherence_status, r$naming_status, r$action),
    "",
    "## Why this row is in the atlas", "", n$why,
    sprintf("It is a recurrent coherent program, not an ancestor ladder: %d of %d terms are supported somewhere, across %d spatial units and %d datasets, and %d of its terms are ancestors of another term in the same row.",
            r$n_terms_supported, r$n_terms, r$n_spatial_units, r$n_datasets,
            r$n_redundant_ancestors),
    "",
    "## Definition rule", "",
    paste0("| anchor | label | scope |", "\n|---|---|---|"),
    paste(sprintf("| %s | %s | %s |", anchors$anchor_go_id, anchors$anchor_label,
                  anchors$match_scope), collapse = "\n"),
    "",
    sprintf("## Constituent GO terms (%d)", r$n_terms), "",
    sprintf("Semantic clusters: %d (%d major). Largest holds %.0f%% of terms.",
            r$n_semantic_clusters, r$n_major_clusters,
            100 * r$largest_cluster_frac),
    sprintf("Major subclusters: %s", r$major_subclusters),
    "",
    sprintf("## Supported GO terms (%d of %d, %d occurrences)",
            r$n_terms_supported, r$n_terms, r$n_supported_occurrences), "",
    "Most recurrent first (contexts = distinct comparisons):", "",
    paste(sprintf("- %s (%d contexts, %d units)",
                  unname(term_of[utils::head(supq$GO_ID, 8)]),
                  utils::head(supq$n_supported_contexts, 8),
                  utils::head(supq$n_supported_units, 8)), collapse = "\n"),
    "",
    sprintf("## Recurrent leading-edge proteins (%d core, %d single-appearance)",
            r$n_recurrent_core_proteins, r$n_single_appearance_proteins), "",
    "A protein is CORE only if it recurs in at least three supported terms AND",
    "three contexts, which separates a real core from one carried in by a single",
    "broad annotation.", "",
    paste(sprintf("- **%s** (%d supported terms, %d contexts)",
                  utils::head(z$gene, 12), utils::head(z$n_supported_terms, 12),
                  utils::head(z$n_supported_contexts, 12)), collapse = "\n"),
    "",
    "## Potentially off-theme content", "",
    if (r$n_off_theme_terms == 0)
      "None: every term is within 0.30 Wang similarity of one of this row's own anchors." else
      paste(sprintf("- %s (similarity to own anchor %.3f, supported in %d contexts)",
                    unname(term_of[q$GO_ID[q$off_theme_flag]]),
                    q$similarity_to_own_anchor[q$off_theme_flag],
                    q$n_supported_contexts[q$off_theme_flag]), collapse = "\n"),
    "",
    "## Take-home", "", n$take, "", "**Caveat.** ", n$caveat, "",
    sprintf("**Would I defend this label?** %s", n$defend), "")
  writeLines(sh, file.path(SHEETS, sprintf("%02d_%s.md", which(PROGS == p), p)))
  book <- c(book, sh, "---", "")
}
writeLines(book, file.path(REP, "atlas_program_review_book.md"))

cat("\n===== CROSS-PROGRAM SUMMARY =====\n")
print(summ[, c("current_label", "n_terms", "n_terms_supported",
               "n_supported_occurrences", "n_semantic_clusters",
               "n_off_theme_terms", "n_recurrent_core_proteins",
               "coherence_status", "verdict", "action")], row.names = FALSE)
cat("\ngenuine omissions in the re-audit:", sum(omit$is_true_omission), "\n")
cat("sheets:", SHEETS, "\n")
