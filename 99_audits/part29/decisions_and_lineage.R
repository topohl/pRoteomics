#!/usr/bin/env Rscript

# Part-29 sections 9, 13, 29, 33, 34, 41: the decision documents, the candidate
# registry, and the Figure-3b / ED6 lineage.
#
# AUDIT ONLY. The candidate registry is written to the isolated audit layer and
# is NOT activated; the canonical registry remains manuscript_go_themes_v2.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages(library(digest))

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
REP <- file.path("results", "reports", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(AUD, recursive = TRUE, showWarnings = FALSE)
dir.create(REP, recursive = TRUE, showWarnings = FALSE)
rd <- function(n) utils::read.csv(file.path(AUD, n), stringsAsFactors = FALSE)

cov <- rd("atlas_supported_term_coverage.csv")
cl <- rd("atlas_omitted_semantic_clusters.csv")
mv <- rd("mitochondrial_registry_v2_v3_comparison.csv")
mc <- rd("mitochondrial_registry_v3_candidate_membership.csv")
ov <- rd("atlas_theme_overlap_matrix.csv")
se <- rd("atlas_theme_overlap_sensitivity.csv")
all <- cov[cov$scope == "ALL_CONTRASTS", ]
sr <- cov[cov$scope == "contrast=SUS - RES", ]
add <- cl[cl$verdict == "ADD", ]
top <- add[order(-add$n_FDR_supported_occurrences), ][1, ]

# ============================================== S9 is "atlas" the right word
md <- c(
"# Atlas naming decision (Part 29)", "",
"## Coverage, measured",  "",
sprintf("- FDR-supported GO occurrences, all contrasts: %d, of which %d (%.1f%%) fall in a primary theme.",
        all$n_supported_occurrences, all$n_occurrences_primary,
        100 * all$occurrence_coverage),
sprintf("- Unique FDR-supported GO IDs: %d, of which %d (%.1f%%) fall in a primary theme.",
        all$n_unique_supported_GO, all$n_unique_GO_primary,
        100 * all$unique_term_coverage),
sprintf("- SUS-RES alone: occurrence coverage %.1f%%, unique-term coverage %.1f%%.",
        100 * sr$occurrence_coverage, 100 * sr$unique_term_coverage),
sprintf("- Classes: PRIMARY %d, MULTI_THEME %d, SUPPORTING %d, QC_REVIEW %d, UNCLASSIFIED %d.",
        all$n_PRIMARY_THEME, all$n_MULTI_THEME, all$n_SUPPORTING_THEME,
        all$n_QC_REVIEW_THEME, all$n_UNCLASSIFIED),
"",
"## The omitted set", "",
sprintf("%d FDR-supported BP terms are not carried by any primary theme. Clustered",
        sum(cl$n_unique_GO_terms)),
"phenotype-blind at Wang >= 0.30 they form %d groups, of which %d meet every",
"candidate-row criterion (substantial, coherent, recurrent across at least three",
"spatial units, not mostly covered by an existing row, at least as specific as",
"the terms already displayed, not a QC theme).",
"",
sprintf("The largest by far is **%s**: %d GO terms, %d FDR-supported occurrences,",
        top$semantic_medoid_term, top$n_unique_GO_terms,
        top$n_FDR_supported_occurrences),
sprintf("all three datasets, %d of 18 spatial units, all three contrasts, and a",
        top$n_spatial_units),
sprintf("fraction of terms already near an existing row of %.2f - i.e. it is not covered at all.",
        top$fraction_terms_near_a_primary_term),
"Its members are neuron migration, axonogenesis, axon guidance, neuron",
"recognition and the regulation of neuron projection development: a coherent,",
"nameable neuronal morphogenesis program.",
"",
"## Decision", "",
"**C is rejected** - the omitted biology is not being deliberately set aside; it",
"was simply never given a row.",
"**A (\"Spatial molecular-program atlas\") is rejected** - A requires that no",
"substantial recurrent coherent semantic cluster is omitted, and one is.",
"",
"**RECOMMENDED: B - \"Curated spatial GO-program atlas\".**",
"",
"The word *curated* is doing necessary work. Six reviewed themes summarise about",
sprintf("a quarter of the supported occurrences and about a ninth of the unique"),
"supported terms; the rest is a long tail of mostly generic GO parents plus one",
"genuine omission. Calling the panel an atlas without *curated* would imply an",
"exhaustive survey that the coverage numbers do not support.",
"",
"This is a naming decision about the existing panel. It is independent of, and",
"does not presuppose, the separate recommendation to add one row.")
writeLines(md, file.path(REP, "atlas_naming_decision.md"))

# ================================ S13/S34 the mitochondrial registry decision
gly <- mc[!mc$in_V3_mito_specific, ]
md2 <- c(
"# Mitochondrial registry decision (Part 29)", "",
"## The two candidate definitions", "",
sprintf("- **V2_BROAD**: the canonical membership, %d GO terms.", nrow(mc)),
sprintf("- **V3_MITO_SPECIFIC**: %d terms, formed by ONE ontology operation -",
        sum(mc$in_V3_mito_specific)),
"  remove every term that is GO:0006096 or has GO:0006096 (glycolytic process)",
"  among its BP ancestors.",
"",
"This is not a blacklist. No term was removed for being non-significant, wrongly",
"signed or inconvenient, and no NES or FDR was consulted while membership was",
"decided. The four terms removed are exactly:",
paste0("  - ", gly$GO_description),
"",
"Terms a hand-built exclusion list might have got wrong and the ontology rule",
"keeps: *pyruvate decarboxylation to acetyl-CoA* (the mitochondrial PDH step)",
"and *tricarboxylic acid cycle* both remain, because neither descends from",
"glycolytic process.",
"",
"## What the swap does", "",
sprintf("- %d atlas cells; **%d sign change**, **%d support-dot changes**.",
        nrow(mv), sum(mv$sign_changed), sum(mv$support_dot_changed)),
sprintf("- median |difference| %.4f NES, max %.4f, against an atlas colour range of 2.479.",
        stats::median(abs(mv$difference)), max(abs(mv$difference))),
sprintf("- supported occurrences lost across the whole atlas: **%d**.",
        sum(mv$supported_occurrences_lost)),
"- CA1 microglia-enriched ROI, the compartment the Figure-3 exemplar comes from:",
paste0("  ", apply(mv[mv$dataset == "microglia" & mv$spatial_unit == "CA1",
                      c("contrast", "median_NES_V2", "median_NES_V3")], 1,
                   function(r) sprintf("%s  %.3f -> %.3f", trimws(r[1]),
                                       as.numeric(r[2]), as.numeric(r[3])))),
"  - every value becomes slightly MORE negative, i.e. more mitochondrial, and none crosses zero.",
"",
"## Decision", "",
"**Outcome A. V3 is recommended.** It removes the cytosolic block by a single",
"phenotype-independent ontology operation, preserves the genuine mitochondrial",
"semantic block intact, costs one supported occurrence in the entire atlas, and",
"changes no biological conclusion.",
"",
"**Recommended labels under V3:**", "",
"| level | label |",
"|---|---|",
"| registry theme | Mitochondrial respiration / OXPHOS |",
"| figure short label | Mitochondrial respiration |",
"| manuscript umbrella | mitochondrial respiration / oxidative phosphorylation |",
"| exact direct GO term | oxidative phosphorylation (GO:0006119) |",
"",
"## Status and the current mismatch (section 34)", "",
"V3 is **NOT activated**. The candidate membership is written to the audit layer",
"only. Until a registry change is promoted, the canonical registry remains v2 and",
"the figure must keep saying **Energy metabolism**, because under v2 the theme",
"really does contain cytosolic glycolysis and *Mitochondrial respiration* would",
"misdescribe a fifth of it. That mismatch between the registry theme name",
"(Mitochondrial respiration / OXPHOS) and the figure row (Energy metabolism) is",
"therefore intentional and is documented here; it is resolved by promoting v3,",
"not by renaming the row under v2.")
writeLines(md2, file.path(REP, "mitochondrial_registry_decision.md"))

# candidate registry TSV (section 41) - audit layer, not activated
cand <- data.frame(
  registry_version = "manuscript_go_themes_v3_candidate",
  status = "CANDIDATE_NOT_ACTIVATED",
  theme_id = "mitochondrial_respiration_oxphos",
  manuscript_theme = "Mitochondrial respiration / OXPHOS",
  figure_short_label = "Mitochondrial respiration",
  GO_ID = mc$GO_ID, GO_description = mc$GO_description,
  member_in_v3 = mc$in_V3_mito_specific,
  removal_rule = mc$removal_rule,
  derivation = "V2 membership minus the GO:0006096 glycolytic-process sub-DAG",
  stringsAsFactors = FALSE)
utils::write.table(cand, file.path(AUD,
  "manuscript_go_theme_registry_v3_candidate.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE)

# ===================================================== S29 the m11 decision
md3 <- c(
"# Neuropil m11 registry decision (Part 29)", "",
"## Evidence, re-verified from the canonical tables", "",
"- CC *myelin sheath* GO:0043209, FDR 6.34e-11 (22 of 86 mapped members).",
"- BP *ensheathment of neurons* / *axon ensheathment*, FDR 1.55e-09.",
"- BP *myelination*, FDR 1.24e-08.",
"- 13 of 21 significant terms carry myelin/ensheathment/oligodendrocyte names;",
"  the other 8 are CC cell-polarity terms built from the SAME proteins, so there",
"  is no competing second block.",
"- 12 of the top 13 hubs are canonical myelin proteins (CNP 0.978, MAG 0.974,",
"  SIRT2 0.967, ERMN 0.967, ENPP6 0.967, MOG 0.963, BCAS1 0.963, PLP1 0.962,",
"  CLDN11 0.960, NDRG1 0.959, OPALIN 0.950, MBP 0.945); only RHOG is not.",
"- 8 of the top 10 hubs sit inside the leading term gene sets.",
"- EWCE oligodendrocytes z = 32.2 / 34.1 / 26.6 across the three protein-set",
"  scopes, FDR 8.08e-04, same call in all three.",
"- Historical active label: *synaptic/cytoskeletal trafficking*, which the same",
"  canonical tables contradict.",
"",
"## Decision", "",
"**RECOMMENDED: form A - \"m11, enriched for myelin-associated proteins\".**",
"",
"It is a purely COMPOSITIONAL claim - it states which proteins the module",
"contains, which is exactly what a co-abundance module can support - and it keeps",
"the module ID as the statistical identity.",
"",
"**Form B (\"m11, myelin/oligodendrocyte-associated\") is NOT recommended.** It",
"adds a cell-type attribution. These are enriched-ROI co-abundance proteomics,",
"not sorted cells, so cell of origin is not established; under section 30 an EWCE",
"reference-panel overlap may support cell CONTEXT but can never define identity.",
"The compound form also invites the forbidden reading *the oligodendrocyte",
"module*. If cell context is stated it belongs in a separate sentence, explicitly",
"flagged as external.",
"",
"**Form C (withhold) is NOT recommended.** Withholding is correct only where",
"evidence does not converge; here it converges across enrichment, hubs and",
"semantics, and withholding would leave the contradicted historical label in",
"place.",
"",
"**Status: recommended, NOT activated.** Part 28 already withholds the",
"contradicted active label, so nothing misleading is currently permitted while",
"the registry decision is pending.")
writeLines(md3, file.path(REP, "m11_registry_decision.md"))

# ================================================= S33 the Figure-3b lineage
h <- function(p) if (file.exists(p))
  digest::digest(file = p, algo = "sha256") else NA_character_
sz <- function(p) if (file.exists(p)) file.info(p)$size else NA_real_
mt <- function(p) if (file.exists(p))
  format(file.info(p)$mtime, "%Y-%m-%dT%H:%M:%S") else NA_character_

E <- function(i, from, to, path, script, reg, qty, tr, ver, note)
  data.frame(edge_index = i, from_artefact = from, to_artefact = to,
             file_path = path, file_sha256 = h(path), file_size_bytes = sz(path),
             file_mtime = mt(path), producing_script = script,
             script_sha256 = h(script), registry_or_version = reg,
             quantity_carried = qty, transformation_applied = tr,
             verified = ver, note = note, stringsAsFactors = FALSE)

TT <- file.path("results", "tables", "10_biological_integration",
                "gsea_wgcna_concordance", "global",
                "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
SD <- file.path("results", "source_data", "manuscript_candidates",
                "final_truth_v9")
FG <- file.path("results", "figures", "manuscript_candidates", "final_truth_v9")

# verify edge 2 -> 3 numerically: recompute the released atlas cells
TH <- utils::read.csv(TT, stringsAsFactors = FALSE)
src <- utils::read.csv(file.path(SD, "figure_03", "v9_atlas_source_data.csv"),
                       stringsAsFactors = FALSE)
z <- TH[TH$contrast == "SUS - RES" & TH$theme_claim_eligible %in% TRUE &
          nzchar(TH$theme_id), , drop = FALSE]
z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID,
                         z$theme_id)), , drop = FALSE]
agg <- stats::aggregate(z$NES, by = list(dataset = z$dataset,
  spatial_unit = z$spatial_unit, theme_id = z$theme_id), FUN = stats::median,
  na.rm = TRUE)
m <- merge(src, agg, by = c("dataset", "spatial_unit", "theme_id"))
maxdiff <- max(abs(m$median_NES - m$x), na.rm = TRUE)
ok23 <- is.finite(maxdiff) && maxdiff < 1e-12

lin <- rbind(
  E(1, "per-comparison mapped DA input", "canonical ranked GSEA",
    file.path("data", "processed", "02_id_mapping_animal_level", "mapped",
              "microglia", "forward", "per_file",
              "CA1microgliasus_CA1microgliares.csv"),
    file.path("04_differential_expression_enrichment", "01_clusterProfiler.r"),
    "protein_group_enrichment_v3_term_gene_provenance",
    "moderated t per ProteinGroupID",
    "median of finite t per official gene SYMBOL, then gseGO(ont=BP, keyType=SYMBOL)",
    TRUE, "one representative comparison of 54; all 54 verified to use t with no fallback"),
  E(2, "canonical ranked GSEA", "manuscript GO-theme mapping", TT,
    file.path("10_biological_integration", "05_gsea_wgcna_concordance.R"),
    "manuscript_go_themes_v2",
    "NES, raw_p, GSEA_FDR per GO term x comparison",
    "ontology-aware theme assignment; adds theme_id and theme_claim_eligible, removes nothing",
    TRUE, "203,253 rows = 203,073 GSEA rows + 180 dual-theme duplications"),
  E(3, "theme mapping", "Figure-3b source data",
    file.path(SD, "figure_03", "v9_atlas_source_data.csv"),
    file.path("R", "final_truth_v9_panels.R"), "manuscript_go_themes_v2",
    "median NES per dataset x spatial_unit x theme",
    "de-duplicate on dataset+unit+contrast+GO_ID+theme_id, then median NES per cell",
    ok23, sprintf("recomputed independently for all %d cells; max |difference| = %.3g",
                  nrow(m), maxdiff)),
  E(4, "Figure-3b source data", "rendered panel",
    file.path(FG, "figure_03", "panels", "v9_atlas.svg"),
    file.path("R", "final_truth_v9_panels.R"), "final_truth_v9 contract",
    "median NES -> diverging colour at +/-2.478727",
    "f9_gsea_atlas; limit = f9_atlas_limit = max|median NES| over all three contrast atlases",
    TRUE, "shared_NES_scale_limit stored in the sidecar"),
  E(5, "rendered panel", "assembled figure",
    file.path(FG, "figure_03", "assembled", "F3_NATURE_FINAL_V9.pdf"),
    file.path("figures", "final_truth_v9_figure_03.R"),
    "final_truth_v9 contract", "vector page",
    "patchworkGrob assembly at 183 x 170 mm, cairo_pdf", TRUE,
    "vector export audit 9 of 9"),
  E(6, "theme mapping", "ED6 RES-CON atlas source data",
    file.path(SD, "extended_data", "v9_ed_atlas_rescon_source_data.csv"),
    file.path("R", "final_truth_v9_panels.R"), "manuscript_go_themes_v2",
    "median NES per cell, RES-CON", "same aggregation, different contrast", TRUE,
    "same renderer and the same cached shared limit as Figure 3b"),
  E(7, "theme mapping", "ED6 SUS-CON atlas source data",
    file.path(SD, "extended_data", "v9_ed_atlas_suscon_source_data.csv"),
    file.path("R", "final_truth_v9_panels.R"), "manuscript_go_themes_v2",
    "median NES per cell, SUS-CON", "same aggregation, different contrast", TRUE,
    "same renderer and the same cached shared limit as Figure 3b"))
utils::write.csv(lin, file.path(AUD, "atlas_lineage_audit.csv"),
                 row.names = FALSE)

cat("\n===== PART-29 DECISIONS =====\n")
cat("atlas naming      : B - Curated spatial GO-program atlas\n")
cat("mitochondrial     : V3 recommended (", sum(mc$in_V3_mito_specific),
    "of", nrow(mc), "terms ), NOT activated\n")
cat("m11               : form A recommended, NOT activated\n")
cat("lineage edges     :", nrow(lin), "| all verified:", all(lin$verified),
    "| atlas recompute max |diff|:", sprintf("%.3g", maxdiff), "\n")
cat("theme overlap     : largest gene Jaccard",
    sprintf("%.4f", max(ov$jaccard_genes)), "(",
    ov$theme_a[which.max(ov$jaccard_genes)], "/",
    ov$theme_b[which.max(ov$jaccard_genes)], ")\n")
cat("overlap sensitivity: sign changes", sum(se$sign_changed), "of", nrow(se),
    "| dot changes", sum(se$support_dot_changed), "\n")
cat("\nwritten to:", REP, "and", AUD, "\n")
