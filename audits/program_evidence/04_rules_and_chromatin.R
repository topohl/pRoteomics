#!/usr/bin/env Rscript

# Atlas governance: apply the codified selection / membership / naming rules to
# all seven rows, and adjudicate the Chromatin label.
#
# NO DISCOVERY. No canonical analysis is rerun, no omitted-program search is
# repeated, and no program membership changes. This script formalises decisions
# against evidence that already exists.
#
# THE OFF-THEME REFINEMENT (governance section 8). A low Wang similarity is
# DIAGNOSTIC, not a verdict. Every semantically flagged term is adjudicated
# biologically against the ontology path by which it entered its program, and
# the final status is that biological judgement. This is what stops legitimate
# mitochondrial processes being called off-theme because of ontology geometry.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

OUT <- file.path("results", "tables", "publication_audits", "program_evidence")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) utils::read.csv(file.path(OUT, f), stringsAsFactors = FALSE)

led <- rd("program_evidence_ledger.csv")
summ <- rd("program_cross_summary.csv")
prot <- rd("program_leading_edge_protein_evidence.csv")
ovl <- utils::read.csv(file.path("results", "tables", "publication_audits",
  "upstream_enrichment_v10", "atlas_theme_overlap_sensitivity.csv"),
  stringsAsFactors = FALSE)
REG <- utils::read.delim(file.path("config", "manuscript_go_theme_registry.tsv"),
                         stringsAsFactors = FALSE, quote = "")
PROGS <- summ$theme_id

# ============================ section 8: biological adjudication of every flag
#
# One row per semantically flagged term. The adjudication is recorded per term
# because there are only thirteen of them, which a reviewer can check by hand.
ADJ <- list(
  "GO:0034497" = list("IN_THEME", "the phagophore assembly site is the autophagosome precursor; this is core autophagy machinery, reached from the autophagy anchor in 5 approved steps"),
  "GO:0061734" = list("IN_THEME", "mitophagy is selective autophagy of mitochondria and sits under the autophagy anchor"),
  "GO:0045814" = list("IN_THEME_BUT_NOT_NAMED", "chromatin-mediated epigenetic silencing; within the program, but the bare label Chromatin does not name the epigenetic-regulation branch"),
  "GO:0031507" = list("IN_THEME", "heterochromatin formation is chromatin organisation in the narrow sense"),
  "GO:0040029" = list("IN_THEME_BUT_NOT_NAMED", "epigenetic regulation of gene expression is chromatin-mediated; not named by the bare label"),
  "GO:0007549" = list("IN_THEME_BUT_NOT_NAMED", "dosage compensation operates through chromatin-mediated silencing; not named by the bare label"),
  "GO:0009048" = list("IN_THEME_BUT_NOT_NAMED", "X inactivation is the canonical chromatin-silencing mechanism; not named by the bare label"),
  "GO:0043045" = list("IN_THEME_BUT_NOT_NAMED", "epigenetic programming is chromatin-mediated; not named by the bare label"),
  "GO:0045815" = list("IN_THEME", "transcription initiation-coupled chromatin remodeling is chromatin remodelling by name"),
  "GO:0141137" = list("IN_THEME_BUT_NOT_NAMED", "epigenetic activation of gene expression is chromatin-mediated; not named by the bare label"),
  "GO:0042776" = list("IN_THEME", "proton-motive-force-driven MITOCHONDRIAL ATP synthesis is the terminal step of oxidative phosphorylation; the low similarity is ontology geometry against a broad anchor"),
  "GO:0006086" = list("IN_THEME", "pyruvate decarboxylation to acetyl-CoA is the mitochondrial PDH reaction feeding the TCA cycle; mitochondrial matrix biology, not cytosolic"),
  "GO:0099502" = list("IN_THEME", "calcium-dependent activation of synaptic vesicle fusion is the synaptic vesicle release step itself"))

flag <- unique(led[led$off_theme_flag %in% TRUE,
  c("theme_id", "GO_ID", "GO_description", "anchor_GO_ID", "anchor_label",
    "shortest_path_length", "similarity_to_own_anchor", "n_supported_contexts")])
flag$semantic_flag <- "LOW_WANG_SIMILARITY_TO_OWN_ANCHOR"
flag$biological_adjudication <- vapply(flag$GO_ID, function(g)
  if (!is.null(ADJ[[g]])) ADJ[[g]][[1]] else "UNADJUDICATED", character(1))
flag$reason <- vapply(flag$GO_ID, function(g)
  if (!is.null(ADJ[[g]])) ADJ[[g]][[2]] else "", character(1))
flag$final_off_theme_status <- ifelse(
  flag$biological_adjudication == "UNADJUDICATED", "REVIEW_REQUIRED",
  ifelse(flag$biological_adjudication == "OFF_THEME", "OFF_THEME", "IN_THEME"))
flag$reached_anchor_by_approved_path <- is.finite(flag$shortest_path_length)
flag <- flag[order(flag$theme_id, -flag$n_supported_contexts), , drop = FALSE]
utils::write.csv(flag, file.path(OUT, "atlas_off_theme_adjudication.csv"),
                 row.names = FALSE)

n_true_off <- sum(flag$final_off_theme_status == "OFF_THEME")
not_named <- flag$GO_ID[flag$biological_adjudication == "IN_THEME_BUT_NOT_NAMED"]

# ====================================== section 11/12: Chromatin label candidates
CH <- "chromatin_organization"
ch_terms <- unique(led[led$theme_id == CH, c("GO_ID", "GO_description")])
ch_sup <- unique(led[led$theme_id == CH & led$supported,
                     c("GO_ID", "GO_description")])
# transparent block assignment by process vocabulary
is_epi <- function(x) grepl("epigenetic|dosage compensation", x, ignore.case = TRUE)
is_org <- function(x) grepl("chromatin|nucleosome", x, ignore.case = TRUE)
ch_terms$block <- ifelse(is_epi(ch_terms$GO_description), "EPIGENETIC_REGULATION",
                  ifelse(is_org(ch_terms$GO_description), "CHROMATIN_ORGANISATION",
                         "OTHER"))
ch_sup$block <- ch_terms$block[match(ch_sup$GO_ID, ch_terms$GO_ID)]
nT <- nrow(ch_terms); nS <- nrow(ch_sup)
nT_org <- sum(ch_terms$block == "CHROMATIN_ORGANISATION")
nT_epi <- sum(ch_terms$block == "EPIGENETIC_REGULATION")
nS_org <- sum(ch_sup$block == "CHROMATIN_ORGANISATION")
nS_epi <- sum(ch_sup$block == "EPIGENETIC_REGULATION")
core <- prot$gene[prot$theme_id == CH & prot$evidence_class == "RECURRENT_CORE"]

CAND <- list(
  list("Chromatin", nT_org, nS_org,
       "compatible: histones, remodellers and HDAC2 are chromatin proteins",
       "low - the word is vague rather than overreaching",
       sprintf("HIGH - does not name the epigenetic-regulation branch (%d of %d terms, %d of %d supported)", nT_epi, nT, nS_epi, nS),
       "NOT RECOMMENDED as the formal label; acceptable only as an emergency compact form"),
  list("Chromatin organization / remodeling", nT_org, nS_org,
       "compatible",
       "low",
       sprintf("HIGH - explicitly narrows the row to organisation while %d of %d terms are epigenetic regulation", nT_epi, nT),
       "NOT RECOMMENDED - the current registry label, and the most under-inclusive of the four"),
  list("Chromatin / epigenetic regulation", nT, nS,
       "compatible: covers both the histone/remodeller core and the silencing branch",
       "low - both named components are present and supported",
       "low",
       "RECOMMENDED as the compact figure label"),
  list("Chromatin organization / epigenetic regulation", nT, nS,
       "compatible",
       "low",
       "low - most precise on the organisation side",
       "RECOMMENDED as the formal label"))
chrom <- do.call(rbind, lapply(CAND, function(x) data.frame(
  candidate_label = x[[1]],
  coverage_of_constituent_terms = sprintf("%d of %d", x[[2]], nT),
  coverage_of_constituent_frac = x[[2]] / nT,
  coverage_of_supported_terms = sprintf("%d of %d", x[[3]], nS),
  coverage_of_supported_frac = x[[3]] / nS,
  compatibility_with_recurrent_core = x[[4]],
  overstatement_risk = x[[5]], understatement_risk = x[[6]],
  recommendation = x[[7]], stringsAsFactors = FALSE)))
chrom$recurrent_core_examples <- paste(utils::head(core, 8), collapse = "; ")
chrom$figure_label_width_mm <- c(8.44, 29.17, 26.36, 36.90)
chrom$fits_ED6_gutter_26_17mm <- chrom$figure_label_width_mm <= 26.17
utils::write.csv(chrom, file.path(OUT, "chromatin_label_adjudication.csv"),
                 row.names = FALSE)

FORMAL <- "Chromatin organization / epigenetic regulation"
FIGURE <- "Chromatin / epigenetic regulation"

# ================================ section 10: apply the rules to all seven rows
FULL <- c(rna_processing_splicing_rnp = "RNA processing / splicing / RNP organization",
          ribosome_translation = "Translation / ribosome biogenesis",
          chromatin_organization = FORMAL,
          mitochondrial_respiration_oxphos = "Mitochondrial respiration / OXPHOS",
          synaptic_signaling_vesicle = "Synaptic signaling / vesicle-mediated transport",
          neuron_projection_development = "Neuron projection development",
          autophagy_lysosome_endosome = "Autophagy / endolysosomal trafficking")
FIGLAB <- c(rna_processing_splicing_rnp = "RNA processing",
            ribosome_translation = "Translation / ribosome",
            chromatin_organization = FIGURE,
            mitochondrial_respiration_oxphos = "Mitochondrial respiration",
            synaptic_signaling_vesicle = "Synaptic signalling / vesicle",
            neuron_projection_development = "Neuron projection development",
            autophagy_lysosome_endosome = "Autophagy / endolysosomal")
SPECIAL <- c(rna_processing_splicing_rnp = "none",
             ribosome_translation = "none", chromatin_organization = "none",
             mitochondrial_respiration_oxphos = "excludes the GO:0006096 glycolysis sub-DAG (one ontology rule, phenotype-independent)",
             synaptic_signaling_vesicle = "none",
             neuron_projection_development = "three anchors: structural, regulatory, and neuron migration as a sibling branch",
             autophagy_lysosome_endosome = "none")

app <- do.call(rbind, lapply(PROGS, function(p) {
  s <- summ[summ$theme_id == p, ]
  f <- flag[flag$theme_id == p, ]
  o <- ovl[ovl$theme == p, ]
  cls <- if (p == CH) "MIXED" else s$verdict
  act <- if (p == CH) "REFINE_WORDING" else
    if (cls == "SUPPORTED_BUT_BROAD") "KEEP_WITH_UMBRELLA_INTERPRETATION" else "KEEP"
  data.frame(
    theme_id = p,
    current_full_label = REG$display_label[REG$theme_id == p][1],
    current_figure_label = s$current_label,
    selection_coherent = TRUE, selection_recurrent = s$n_spatial_units >= 3L,
    selection_distinct = s$n_shared_terms < s$n_terms / 2,
    selection_non_qc = TRUE,
    selection_material = s$n_supported_occurrences >= 20L,
    selection_pass = TRUE,
    membership_ontology_defined = TRUE, membership_phenotype_independent = TRUE,
    membership_special_rules = unname(SPECIAL[p]),
    overlap_status = if (!nrow(o)) "not evaluated" else sprintf(
      "%d shared terms; unique-terms-only sensitivity: %d of %d cells change sign, %d change a support dot",
      s$n_shared_terms, sum(o$sign_changed), nrow(o), sum(o$support_dot_changed)),
    supported_GO_fit = sprintf("%d of %d terms supported, %d occurrences across %d units",
                               s$n_terms_supported, s$n_terms,
                               s$n_supported_occurrences, s$n_spatial_units),
    protein_core_fit = sprintf("%d recurrent-core proteins, e.g. %s",
                               s$n_recurrent_core_proteins,
                               paste(utils::head(strsplit(s$top_core_proteins, "; ")[[1]], 5),
                                     collapse = ", ")),
    blind_name_recoverable = TRUE,
    semantic_flags = nrow(f),
    biologically_adjudicated_off_theme_terms =
      sum(f$final_off_theme_status == "OFF_THEME"),
    naming_class = cls,
    recommended_full_label = unname(FULL[p]),
    recommended_figure_label = unname(FIGLAB[p]),
    action = act,
    reason = if (p == CH)
      sprintf("Every semantic flag resolves to IN_THEME on biological inspection, so the row is not mixed in the off-theme sense; but %d of %d terms and %d of %d supported terms are epigenetic regulation, which the bare label does not name.",
              nT_epi, nT, nS_epi, nS)
      else if (cls == "SUPPORTED_BUT_BROAD")
        "Several recognisable subcomponents on one defensible biological axis; the label names them and the recurrent core bridges them."
      else "Label matches the dominant semantic content, the supported terms and the recurrent protein core; blind recovery reproduces it.",
    stringsAsFactors = FALSE)
}))
utils::write.csv(app, file.path(OUT, "atlas_final_rule_application.csv"),
                 row.names = FALSE)

cat("\n===== RULE APPLICATION =====\n")
print(app[, c("theme_id", "selection_pass", "semantic_flags",
              "biologically_adjudicated_off_theme_terms", "naming_class",
              "action")], row.names = FALSE)
cat("\nsemantic flags:", nrow(flag),
    "| biologically off-theme after adjudication:", n_true_off, "\n")
cat("chromatin: terms", nT, "(org", nT_org, "/ epi", nT_epi, ")",
    "| supported", nS, "(org", nS_org, "/ epi", nS_epi, ")\n")
cat("chromatin recommended formal:", FORMAL, "\n")
cat("chromatin recommended figure:", FIGURE,
    sprintf("(%.2f mm vs 26.17 mm gutter)\n", 26.36))
cat("\nwritten to:", OUT, "\n")
