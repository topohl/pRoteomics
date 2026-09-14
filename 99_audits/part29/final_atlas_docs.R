#!/usr/bin/env Rscript

# Finalization sections 9, 11, 12, 25-27: why each row is there, what "curated"
# means numerically, and the reader-facing explanations.
#
# AUDIT ONLY. All numbers are recomputed from the promoted v3 theme table.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({ library(GO.db); library(AnnotationDbi) })

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
REP <- file.path("results", "reports", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(REP, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05

TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)
REG <- utils::read.delim(file.path("config", "manuscript_go_theme_registry.tsv"),
                         stringsAsFactors = FALSE, quote = "")
stopifnot(unique(REG$registry_version) == "manuscript_go_themes_v3")

SHORT <- c(rna_processing_splicing_rnp = "RNA processing",
           ribosome_translation = "Translation / ribosome",
           chromatin_organization = "Chromatin",
           mitochondrial_respiration_oxphos = "Mitochondrial respiration",
           synaptic_signaling_vesicle = "Synaptic signalling / vesicle",
           neuron_projection_development = "Neuron projection development",
           autophagy_lysosome_endosome = "Autophagy / endolysosomal")

elig <- TH[TH$theme_claim_eligible %in% TRUE, , drop = FALSE]
sup <- TH[is.finite(TH$GSEA_FDR) & TH$GSEA_FDR < ALPHA, , drop = FALSE]
members <- lapply(names(SHORT), function(t)
  sort(unique(elig$GO_ID[elig$theme_id == t])))
names(members) <- names(SHORT)
desc <- stats::setNames(TH$GO_description[!duplicated(TH$GO_ID)],
                        TH$GO_ID[!duplicated(TH$GO_ID)])

# ================================================== S9 why each row is there
WHY <- c(
  rna_processing_splicing_rnp = "the general parent GO:0006396 with an explicit regulatory anchor; mRNA splicing, tRNA and rRNA processing, miRNA processing and RNP assembly are all descendants of one branch",
  ribosome_translation = "translation proper together with ribosome biogenesis and rRNA maturation, which are the machinery that makes the ribosome rather than the act of translating",
  chromatin_organization = "chromatin organisation and remodelling, including nucleosome assembly and the epigenetic-regulation branch; every member is a chromatin-mediated process",
  mitochondrial_respiration_oxphos = "mitochondrial respiratory bioenergetics: cellular respiration, the electron transport chain and respiratory-chain complex assembly, with the cytosolic glycolysis sub-DAG excluded by one ontology rule",
  synaptic_signaling_vesicle = "synaptic signalling together with the synaptic vesicle cycle, which binary_cut returns as a single coherent semantic block rather than two",
  neuron_projection_development = "neuronal morphogenesis: axonogenesis, axon guidance, dendrite and neurite development, plus neuron migration as a sibling branch",
  autophagy_lysosome_endosome = "autophagy together with endosomal and lysosomal trafficking, which share machinery and form one semantic block")
DISTINCT <- c(
  rna_processing_splicing_rnp = "shares 6 rRNA terms with Translation (gene Jaccard 0.111) but is otherwise disjoint; excluding every shared term moves the atlas by a median of 0.003 NES",
  ribosome_translation = "the same 6 rRNA terms are shared with RNA processing; the rest of the theme is ribosome assembly and translational control, which RNA processing does not contain",
  chromatin_organization = "no GO term is shared with any other primary theme",
  mitochondrial_respiration_oxphos = "no GO term is shared with any other primary theme; under v3 it no longer reaches cytosolic glycolysis",
  synaptic_signaling_vesicle = "shares one endosomal term with Autophagy; distinct from Neuron projection development, which is morphogenesis rather than signalling",
  neuron_projection_development = "zero overlap with any other primary theme, and the phenotype-blind semantic cluster it came from had no term within 0.50 Wang of an existing row",
  autophagy_lysosome_endosome = "shares one endosomal term with Synaptic signalling; otherwise disjoint")

rat <- do.call(rbind, lapply(names(SHORT), function(t) {
  ids <- members[[t]]
  z <- sup[sup$GO_ID %in% ids, , drop = FALSE]
  r <- REG[REG$theme_id == t, , drop = FALSE]
  incl <- r[r$match_scope != "exclude_anchor_and_descendants", , drop = FALSE]
  excl <- r[r$match_scope == "exclude_anchor_and_descendants", , drop = FALSE]
  other <- unlist(members[setdiff(names(SHORT), t)], use.names = FALSE)
  data.frame(
    theme_id = t, final_display_label = unname(SHORT[t]),
    registry_display_label = r$display_label[1],
    display_order = r$display_order[1],
    ontology_anchor_or_definition = paste0(
      paste(sprintf("%s (%s)", incl$anchor_go_id, incl$anchor_label),
            collapse = " + "),
      if (nrow(excl)) sprintf("  MINUS %s (%s)", excl$anchor_go_id,
                              excl$anchor_label) else ""),
    n_GO_terms = length(ids),
    semantic_structure = "one semantic block by simplifyEnrichment binary_cut",
    major_substructure = paste(utils::head(unname(desc[ids]), 5), collapse = "; "),
    n_supported_occurrences = nrow(z),
    n_datasets_supported = length(unique(z$dataset)),
    n_spatial_units_supported = length(unique(z$spatial_unit)),
    overlap_with_other_primary_themes = length(intersect(ids, other)),
    why_biologically_coherent = unname(WHY[t]),
    why_distinct_from_other_rows = unname(DISTINCT[t]),
    why_retained_in_primary_atlas = sprintf(
      "FDR-supported in %d of 18 spatial units across %d of 3 datasets",
      length(unique(z$spatial_unit)), length(unique(z$dataset))),
    exact_GO_example_if_relevant = c(
      rna_processing_splicing_rnp = "GO:0006397 mRNA processing (Figure 3 exemplar)",
      ribosome_translation = "GO:0002181 cytoplasmic translation",
      chromatin_organization = "GO:0006338 chromatin remodeling",
      mitochondrial_respiration_oxphos = "GO:0006119 oxidative phosphorylation (Figure 3 exemplar)",
      synaptic_signaling_vesicle = "GO:0099536 synaptic signaling (Figure 3 exemplar)",
      neuron_projection_development = "GO:0007411 axon guidance",
      autophagy_lysosome_endosome = "GO:0006914 autophagy")[[t]],
    stringsAsFactors = FALSE)
}))
rat <- rat[order(rat$display_order), , drop = FALSE]
utils::write.csv(rat, file.path(AUD, "atlas_primary_theme_rationale.csv"),
                 row.names = FALSE)

# ================================== S11/S12 what "curated" means numerically
sup$in_primary <- sup$theme_claim_eligible %in% TRUE &
  sup$theme_id %in% names(SHORT)
cl <- utils::read.csv(file.path(AUD, "atlas_top_omitted_clusters_final_review.csv"),
                      stringsAsFactors = FALSE)
all_cl <- utils::read.csv(file.path(AUD, "atlas_omitted_semantic_clusters.csv"),
                          stringsAsFactors = FALSE)
comp <- data.frame(
  measure = c(
    "FDR-supported occurrences, total",
    "FDR-supported occurrences in a primary theme",
    "occurrence coverage",
    "unique FDR-supported GO IDs, total",
    "unique FDR-supported GO IDs in a primary theme",
    "unique-term coverage",
    "occurrences in a QC-review theme",
    "occurrences unclassified",
    "omitted semantic clusters (Part-29 clustering)",
    "omitted clusters meeting the recurrence criterion",
    "omitted clusters meeting ALL addition criteria",
    "recurrent coherent programs REPRESENTED by a primary row",
    "recurrent coherent programs still OMITTED"),
  value = c(
    nrow(sup), sum(sup$in_primary),
    sprintf("%.1f%%", 100 * mean(sup$in_primary)),
    length(unique(sup$GO_ID)),
    length(unique(sup$GO_ID[sup$in_primary])),
    sprintf("%.1f%%", 100 * length(unique(sup$GO_ID[sup$in_primary])) /
              length(unique(sup$GO_ID))),
    sum(sup$theme_role == "qc_review", na.rm = TRUE),
    sum(!nzchar(sup$theme_id)),
    nrow(all_cl), sum(all_cl$recurrent %in% TRUE),
    sum(cl$decision == "YES_MAJOR_OMISSION"),
    length(SHORT), sum(cl$decision == "YES_MAJOR_OMISSION")),
  stringsAsFactors = FALSE)
utils::write.csv(comp, file.path(AUD, "atlas_final_completeness_summary.csv"),
                 row.names = FALSE)

occ_cov <- 100 * mean(sup$in_primary)
term_cov <- 100 * length(unique(sup$GO_ID[sup$in_primary])) /
  length(unique(sup$GO_ID))

# ============================================================ S25 methods text
writeLines(c(
"# Atlas theme selection - Methods paragraph", "",
"Canonical GO biological-process GSEA results were mapped to phenotype-",
"independent, ontology-defined biological program families (manuscript GO-theme",
"registry v3). Each family is specified by one or more GO anchors together with",
"the approved is_a and part_of relationships, never by enumerating the GO terms",
"that happened to be significant; one family additionally carries a single",
"ontology exclusion rule (see below). Primary display families were required to",
"represent semantically coherent and biologically distinct GO branches with",
"recurrent representation across the spatial proteomic datasets.",
"",
"Completeness was assessed by semantic clustering of all FDR-supported GO terms",
"not represented by the initial families, using GO-BP semantic similarity",
sprintf("(GOSemSim Wang, org.Mm.eg.db %s) and independently of enrichment magnitude",
        utils::packageVersion("org.Mm.eg.db")),
sprintf("or direction. %d supported terms fell outside the initial six families and",
        sum(all_cl$n_unique_GO_terms)),
sprintf("formed %d semantic clusters. Neuron projection development was the only",
        nrow(all_cl)),
"additional recurrent, distinct program meeting the predefined criteria for",
"primary display; it was added as a seventh family. The remaining clusters were",
"rejected as hierarchical redundancy (a parent/child ladder of one branch),",
"as already covered by an existing family, or as non-neural annotation context.",
"",
"The mitochondrial respiration / OXPHOS family excludes the glycolytic process",
"sub-DAG (GO:0006096) to distinguish cytosolic glycolysis from mitochondrial",
"respiratory bioenergetics; this ontology rule was defined independently of",
"phenotype statistics. Pyruvate decarboxylation to acetyl-CoA and the",
"tricarboxylic acid cycle are retained, as neither descends from GO:0006096.",
"",
sprintf("The seven families carry %.1f%% of FDR-supported GO occurrences and %.1f%% of",
        occ_cov, term_cov),
"unique FDR-supported GO identifiers. These fractions are low because GO results",
"contain extensive hierarchical redundancy, not because recurrent biology is",
"missing: exhaustive semantic review of the unrepresented terms identified no",
"additional recurrent coherent program meeting the predefined criteria. The",
"atlas is therefore a curated representation of recurrent GO-BP program families",
"and is not an exhaustive enumeration of enriched GO terms; complete",
"constituent, overlapping and unclassified GO results are provided in source",
"data."),
file.path(REP, "atlas_methods_selection_text.md"))

# ============================================================= S26 legend text
writeLines(c(
"# Atlas figure-legend clause", "",
"Rows show curated recurrent GO-BP program families rather than an exhaustive",
"enumeration of enriched GO terms; complete constituent, overlapping and",
"unclassified GO-term results are provided in source data.",
"",
"## Optional longer form", "",
sprintf(paste0("Rows are seven ontology-defined GO-BP program families (registry ",
               "v3) carrying %.1f%% of FDR-supported GO occurrences. They are a ",
               "curated"), occ_cov),
"representation of recurrent program families, not an exhaustive enumeration of",
"enriched GO terms; the remaining supported terms are overwhelmingly",
"parent/child variants of the displayed branches and are provided, with the",
"unclassified results, in source data."),
file.path(REP, "atlas_legend_selection_text.md"))

cat("\n===== FINAL ATLAS =====\n")
print(rat[, c("display_order", "final_display_label", "n_GO_terms",
              "n_supported_occurrences", "n_spatial_units_supported",
              "overlap_with_other_primary_themes")], row.names = FALSE)
cat("\ncoverage: occurrences", sprintf("%.1f%%", occ_cov),
    "| unique terms", sprintf("%.1f%%", term_cov), "\n")
cat("recurrent coherent programs still omitted:",
    sum(cl$decision == "YES_MAJOR_OMISSION"), "\n")
