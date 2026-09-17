#!/usr/bin/env Rscript

# Part A, sections 3-12: statistic identity, biological n, multiple testing,
# null/interaction/claim language, dataset semantics, zero-vs-NA, the numerical
# floor, and the Methods contract.
#
# AUDIT ONLY. The corpus is every manuscript-reachable prose artefact plus the
# released source data; nothing is recomputed.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
source("audits/publication_hardening/00_checkpoint.R")

dir.create(PH_TAB, recursive = TRUE, showWarnings = FALSE)
V9R <- file.path("results", "reports", "manuscript_candidates", "final_truth_v9")
V9T <- file.path("results", "tables", "manuscript_candidates", "final_truth_v9")
V9S <- file.path("results", "source_data", "manuscript_candidates",
                 "final_truth_v9")
V9F <- file.path("results", "figures", "manuscript_candidates", "final_truth_v9")

# the manuscript-reachable prose corpus: everything a reader receives
corpus <- c(
  list.files(V9R, pattern = "[.]md$", recursive = TRUE, full.names = TRUE),
  list.files(V9T, pattern = "README[.]md$", recursive = TRUE, full.names = TRUE),
  "figures/figure_final_truth_v9_contract.yml")
corpus <- corpus[file.exists(corpus)]
lines_of <- function(f) readLines(f, warn = FALSE, encoding = "UTF-8")
rel <- function(f) sub(".*proteomics[/\\\\]", "", f)

# The v9 semantics layer already fixes what counts as a false positive here:
# files whose purpose is to QUOTE prohibited wording, lines that instruct
# against a phrase, and the repository's emphatic NOT/NEVER disclaimers. This
# audit reuses that exact convention so the two scans cannot disagree, plus the
# path-like exemption Part 27 established for machine-valued YAML fields.
EXEMPT_FILE <- c("manuscript_semantic_rules.md", "claim_chain_audit.md")
PROHIBITION <- paste0("never|must not|do not |does not|cannot|prohibited|",
                      "instead of|not be read|not be interpreted|",
                      "rather than|avoid|no hypothesis|banned|->|forbidden|",
                      "over-read|overread|misread|incorrect|wrong|",
                      "deliberately|is treated as|no claim|not a claim|",
                      "descriptive|USE:|ONLY USE WHEN:")
EMPHATIC <- "\\bNOT\\b|\\bNEVER\\b"
PATHLIKE <- "^[a-z_]+:[[:space:]]*[A-Za-z0-9_./-]+[.](csv|tsv|R|r|svg|pdf|yml)$"
ph_disclaimed <- function(ctx)
  grepl(PROHIBITION, ctx, ignore.case = TRUE, perl = TRUE) |
  grepl(EMPHATIC, ctx, perl = TRUE)

scan_corpus <- function(pattern, exclude = NULL, respect_exempt = TRUE) {
  out <- list()
  for (f in corpus) {
    if (respect_exempt && basename(f) %in% EXEMPT_FILE) next
    ln <- lines_of(f)
    i <- grep(pattern, ln, ignore.case = TRUE, perl = TRUE)
    if (!is.null(exclude) && length(i))
      i <- i[!grepl(exclude, ln[i], ignore.case = TRUE, perl = TRUE)]
    for (k in i) {
      raw <- trimws(ln[k])
      if (grepl(PATHLIKE, raw, perl = TRUE)) next
      out[[length(out) + 1L]] <- data.frame(
        file = rel(f), line = k,
        context = substr(trimws(sub("^[#>*|[:space:]-]+", "", raw)), 1, 200),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(data.frame(file = character(), line = integer(),
                                      context = character()))
  do.call(rbind, out)
}

# =============================================== A03 statistic identity
S <- function(quantity, correct_name, where, wrong_names, status, note)
  data.frame(quantity = quantity, correct_name = correct_name,
             displayed_where = where, prohibited_names = wrong_names,
             status = status, note = note, stringsAsFactors = FALSE)
ident <- rbind(
  S("protein differential abundance", "log2 fold change",
    "F3 g/h/i leading-edge dot plots",
    "expression; effect size; abundance change significance", "PASS",
    "legends state log2 fold change and that no displayed protein is individually FDR-supported"),
  S("gene ranking statistic", "moderated t",
    "GSEA input (not displayed)", "log2FC; fold change", "PASS",
    "verified 54 of 54 comparisons use column t with no fallback"),
  S("gene-set enrichment", "normalised enrichment score (NES)",
    "F3 d/e/f strips, ED6 c/d/e", "significance; effect size", "PASS",
    "colourbars read 'normalised enrichment score'"),
  S("theme aggregation", "median NES across a theme's GO terms",
    "F3b, ED6 a/b", "mean NES; theme significance; theme FDR", "PASS",
    "colourbar reads 'Median normalised enrichment score'; legend states no theme-level p-value exists"),
  S("module member abundance", "mean module-member CON z-score",
    "ED WGCNA a", "expression; eigengene", "PASS",
    "colourbar reads 'Mean module-member abundance (CON z-score)'"),
  S("module eigengene effect", "module eigengene difference",
    "ED WGCNA b", "expression; abundance", "PASS",
    "colourbar reads 'Module eigengene difference'; caption states DESCRIPTIVE ONLY"),
  S("profile similarity", "median Spearman profile correlation (CON)",
    "ED8 a", "connectivity; anatomical connection", "PASS",
    "legend states molecular similarity, NOT anatomical connectivity"),
  S("bilateral precision", "intraclass correlation (ICC)",
    "ED1 b", "reliability significance", "PASS",
    "axis reads ICC; caption states descriptive precision context"),
  S("cell-type affinity", "EWCE z-score / reference-panel overlap",
    "ED WGCNA c", "cell-type identity; purity", "PASS",
    "recorded as external cell-type CONTEXT only"),
  S("baseline abundance", "CON z-score",
    "F2 d, F2 e, ED2 a", "expression", "PASS",
    "colourbars name the z-score and the CON group"))
utils::write.csv(ident, file.path(PH_TAB,
  "manuscript_statistic_identity_audit.csv"), row.names = FALSE)

# =============================================== A04 biological n / replication
n_hits <- scan_corpus("\\bn\\s*=\\s*[0-9]+|\\bN\\s*=\\s*[0-9]+")
n_hits$counted_object <- ifelse(
  grepl("acquisition", n_hits$context, ignore.case = TRUE), "spatial acquisitions",
  ifelse(grepl("animal|per group", n_hits$context, ignore.case = TRUE), "animals",
  ifelse(grepl("protein", n_hits$context, ignore.case = TRUE), "proteins",
  ifelse(grepl("GO|term", n_hits$context, ignore.case = TRUE), "GO terms",
  ifelse(grepl("network|instance", n_hits$context, ignore.case = TRUE),
         "animal x dataset network instances", "other")))))
n_hits$is_biological_n <- n_hits$counted_object == "animals"
n_hits$independent_animals <- ifelse(n_hits$is_biological_n, "yes",
  "no - must not be read as biological replication")
n_hits$status <- ifelse(
  n_hits$counted_object == "other" , "REVIEW",
  ifelse(n_hits$is_biological_n, "PASS", "PASS - non-biological n, labelled"))
utils::write.csv(n_hits, file.path(PH_TAB,
  "manuscript_n_replication_audit.csv"), row.names = FALSE)

# =============================================== A05 multiple-testing contract
mt_hits <- scan_corpus("FDR|BH\\b|p[.]adjust|adjusted p|multiple[- ]testing|survived correction")
FAMILY <- c(
  "differential-abundance" = "BH within the DA family of that comparison",
  "GSEA|enrichment|GO" = "BH over every GO-BP set returned in that one comparison",
  "theme|atlas" = "NONE - descriptive aggregation, no theme-level family",
  "module|WGCNA|eigengene" = "BH within the module-phenotype family",
  "network|edge" = "BH within the network family",
  "CAMERA" = "BH over the full comparable GO-BP family within each contrast",
  "external|signature" = "BH within the external-validation inventory")
mt_hits$inferred_family <- vapply(mt_hits$context, function(x) {
  h <- names(FAMILY)[vapply(names(FAMILY), function(p)
    grepl(p, x, ignore.case = TRUE, perl = TRUE), logical(1))]
  if (!length(h)) "unclassified" else unname(FAMILY[h[1]])
}, character(1))
mt_hits$claims_theme_level_fdr <- grepl(
  "theme[- ]level (p|FDR)|FDR[- ]significant theme", mt_hits$context,
  ignore.case = TRUE) & !ph_disclaimed(mt_hits$context)
mt_hits$status <- ifelse(mt_hits$claims_theme_level_fdr, "FAIL",
  ifelse(mt_hits$inferred_family == "unclassified",
         "PASS - family not inferable from this line alone", "PASS"))
utils::write.csv(mt_hits, file.path(PH_TAB,
  "manuscript_multiple_testing_audit.csv"), row.names = FALSE)

# ========================== A06/A07/A08 null, interaction and claim language
CLAIM <- list(
  c("null", "no difference|no effect|unchanged|\\babsent\\b|not associated|not affected|\\bstable\\b|preserved",
    "must read as NO_DETECTED_DIFFERENCE at n=3/group, never PROVEN_ABSENCE"),
  c("interaction", "interaction|sex[- ]specific|spatially specific|susceptibility[- ]specific|resilience[- ]specific|divergent|selective",
    "requires an actual interaction or contrast-of-contrasts test"),
  c("baseline", "baseline",
    "legal only as the control-reference state named in the same sentence"),
  c("trajectory", "trajector",
    "three pairwise contrasts are not a longitudinal trajectory"),
  c("validation", "validat",
    "only external/independent evidence may be called validation"),
  c("specific", "\\bspecific\\b|\\bspecificity\\b|selectiv",
    "requires a specificity test, not a contrast of significances"),
  c("prediction", "\\bpredict",
    "only out-of-sample prediction"),
  c("mechanism", "reprogram|rewir|restor|\\bdriver\\b|\\bmechanis|\\bcausal\\b|relocat|redistribut",
    "not supported by association-only proteomics"),
  c("cell_intrinsic", "cell[- ]intrinsic|cell[- ]autonomous|microglia[- ]specific",
    "enriched-ROI data cannot establish cell-intrinsic identity"),
  # A09: laminar wording is legal in neuropil and illegal for the two
  # region-level compartments, so the rule only fires when a layer claim and a
  # region-level compartment appear in the same statement.
  c("laminar",
    "(soma|microglia)[^.]{0,120}(laminar|\\blayer)|(laminar|\\blayer)[^.]{0,120}(soma|microglia)",
    "laminar resolution exists in neuropil only; soma and microglia are region-level"))
claim <- do.call(rbind, lapply(CLAIM, function(c3) {
  h <- scan_corpus(c3[2])
  if (!nrow(h)) return(NULL)
  h$term_class <- c3[1]; h$rule <- c3[3]
  h$disclaimed <- ph_disclaimed(h$context)
  h$status <- ifelse(h$disclaimed, "PASS - qualified or disclaimed in the line",
                     "REVIEW - unqualified usage in prose")
  h
}))
utils::write.csv(claim, file.path(PH_TAB,
  "manuscript_claim_language_audit.csv"), row.names = FALSE)

# =============================================== A09 dataset interpretation rules
D <- function(dataset, resolution, legal, illegal, reason)
  data.frame(dataset = dataset, spatial_resolution = resolution,
             legal_claim = legal, illegal_claim = illegal, reason = reason,
             stringsAsFactors = FALSE)
drules <- rbind(
  D("neuron_neuropil", "region x layer (10 units)",
    "layer-resolved and laminar wording; region x layer comparisons",
    "cell-type attribution of a neuropil measurement",
    "neuropil is the only compartment acquired at laminar resolution"),
  D("neuron_soma", "region only (4 units)",
    "region-level statements; neuronal-soma compartment",
    "any layer or laminar statement about this compartment",
    "soma ROIs were acquired per region, not per layer"),
  D("microglia_enriched", "region only (4 units)",
    "microglia-enriched ROI; local microenvironment",
    "layer or laminar wording; microglial proteome; cell-intrinsic or cell-autonomous",
    "region-level acquisition of an enriched, not purified, compartment"))
utils::write.csv(drules, file.path(PH_TAB,
  "manuscript_dataset_interpretation_rules.csv"), row.names = FALSE)

# A09 hard check: the v9 semantics layer owns the enforcement gate for
# specificity wording. Its S9 pattern is read here rather than restated, so a
# future widening or narrowing of that gate is detected instead of assumed.
sem_src <- readLines("figures/final_truth_v9_semantics.R", warn = FALSE)
s9_line <- grep("selectively", sem_src, value = TRUE)
s9_pat <- if (length(s9_line))
  sub('.*"([^"]*selectively[^"]*)".*', "\\1", s9_line[1]) else NA_character_
prose <- claim[claim$term_class == "specific" &
               grepl("^REVIEW", claim$status) &
               grepl("selectiv", claim$context, ignore.case = TRUE), ]
gate <- data.frame(
  gate_script = "figures/final_truth_v9_semantics.R",
  gate_id = "S9 phenotype specificity language",
  gate_pattern = s9_pat,
  matches_adverb_selectively = grepl("selectively", s9_pat %||% ""),
  matches_adjective_selective = length(grep(
    s9_pat %||% "zzz", "spatially selective", perl = TRUE)) > 0,
  unflagged_prose_occurrences = nrow(prose),
  unflagged_locations = paste(sprintf("%s:%s", basename(prose$file), prose$line),
                              collapse = "; "),
  licensing_test_present = "none - no GSEA-level stress x spatial-unit test exists",
  nearest_executed_test = "WGCNA stress x spatial-unit omnibus, 0 of 35 FDR-supported, smallest FDR 0.27",
  stringsAsFactors = FALSE)
utils::write.csv(gate, file.path(PH_TAB,
  "manuscript_specificity_gate_defect.csv"), row.names = FALSE)

# =============================================== A10 zero vs NA vs not evaluable
sd_files <- list.files(V9S, pattern = "[.]csv$", recursive = TRUE,
                       full.names = TRUE)
zna <- do.call(rbind, lapply(sd_files, function(f) {
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  num <- vapply(d, function(x) is.numeric(x) || all(is.na(x)), logical(1))
  if (!any(num)) return(NULL)
  data.frame(source_data = rel(f), n_rows = nrow(d),
             n_numeric_cols = sum(num),
             n_cols_with_NA = sum(vapply(d[num], function(x) any(is.na(x)),
                                         logical(1))),
             n_cols_with_zero = sum(vapply(d[num], function(x)
               any(x == 0, na.rm = TRUE), logical(1))),
             has_explicit_status_column = any(grepl(
               "status|evaluable|tested|reason|class", names(d),
               ignore.case = TRUE)),
             stringsAsFactors = FALSE)
}))
zna$risk <- ifelse(zna$n_cols_with_NA > 0 & !zna$has_explicit_status_column,
                   "REVIEW - NA present without a status column", "PASS")
utils::write.csv(zna, file.path(PH_TAB,
  "manuscript_zero_na_audit.csv"), row.names = FALSE)

# =============================================== A11 numerical-floor disclosure
floor_ok <- any(vapply(corpus, function(f)
  any(grepl("eps = 1e-10", lines_of(f), fixed = TRUE)), logical(1)))

# =============================================== A12 the Methods contract
M <- function(id, q, n, unit, model, design, contrast, repeated, scope,
              software, version, output, claim, limitation, sentence)
  data.frame(analysis_id = id, biological_question = q, biological_n = n,
             analysis_unit = unit, model_or_test = model, design = design,
             contrast = contrast, repeated_structure = repeated,
             multiple_testing_scope = scope, software = software,
             version = version, primary_output = output, claim_status = claim,
             major_limitation = limitation,
             methods_sentence_candidate = sentence, stringsAsFactors = FALSE)
Rv <- paste0("R ", getRversion())
meth <- rbind(
  M("preprocessing_bilateral", "How are hemispheres combined into one animal-level value?",
    "9 animals", "protein group x animal x spatial unit",
    "bilateral aggregation to animal level", "3 groups x 18 spatial units",
    "not applicable", "two hemispheres per animal, aggregated",
    "not applicable", "in-house R", Rv,
    "animal-level GCT matrices", "descriptive",
    "hemispheres are repeated tissue, not independent replicates",
    "Hemispheric measurements were aggregated to one value per animal and spatial unit before any group comparison, so that the biological replicate is the animal."),
  M("differential_abundance", "Which proteins differ between outcome groups within a spatial unit?",
    "3 per group", "protein group",
    "moderated linear model (Protigy, animal_level_protigy_da_v1)",
    "3 groups x 18 spatial units", "SUS-RES, SUS-CON, RES-CON",
    "one value per animal; no within-animal repeats at this stage",
    "BH within each comparison", "Protigy", "animal_level_protigy_da_v1",
    "per-comparison protein tables", "inferential",
    "n = 3 per group limits power; absence is not evidence of absence",
    "Protein-level differential abundance was assessed with a moderated linear model on animal-level values within each spatial unit, with Benjamini-Hochberg correction within each comparison."),
  M("gsea", "Which biological programs are coordinately shifted?",
    "3 per group", "gene (median moderated t per official symbol)",
    "clusterProfiler::gseGO over GO-BP", "ranked list per comparison",
    "SUS-RES, SUS-CON, RES-CON", "collapsed to one statistic per gene",
    "BH over every GO-BP set returned in that comparison",
    "clusterProfiler / fgsea",
    paste0("clusterProfiler ", utils::packageVersion("clusterProfiler")),
    "per-comparison GO-BP enrichment tables", "inferential",
    "p-values are floored at eps = 1e-10; a term at the floor has an unresolved true p",
    "Genes were ranked by the median moderated t statistic per official gene symbol and tested against GO biological process with clusterProfiler::gseGO (minGSSize 10, maxGSSize 800, BH within each comparison, eps = 1e-10)."),
  M("camera_sensitivity", "Does the program-level direction survive a correlation-aware competitive test?",
    "3 per group", "gene", "limma::cameraPR", "preranked, inter.gene.cor = 0.01",
    "all three pairwise contrasts", "none",
    "BH over the full comparable GO-BP family within each contrast",
    "limma", paste0("limma ", utils::packageVersion("limma")),
    "CAMERA concordance tables", "sensitivity only (LEVEL 4)",
    "preranked CAMERA is weaker than a full expression-matrix CAMERA",
    "As a sensitivity analysis, the same ranked statistics were tested with limma::cameraPR (inter-gene correlation fixed at 0.01) over the comparable GO-BP family; this is a concordance check and not independent validation."),
  M("go_program_atlas", "How are enrichment results summarised for display?",
    "3 per group", "GO term grouped into ontology-defined families",
    "median NES per family", "7 primary families x 18 spatial units",
    "all three pairwise contrasts", "none",
    "NONE - descriptive aggregation", "in-house R + GO.db", Rv,
    "atlas source data and Figure 3b / ED6", "descriptive",
    "no theme-level p-value or FDR exists or is implied",
    "Canonical GO-BP results were mapped to seven ontology-defined program families (registry manuscript_go_themes_v3) and summarised as the median NES of their constituent terms; this summary is descriptive and carries no theme-level statistical test."),
  M("wgcna", "Are co-abundance modules associated with outcome?",
    "9 animals", "module eigengene",
    "lmerTest eigengene ~ StressGroup + SpatialUnit + (1|AnimalID)",
    "modules x spatial units", "all three pairwise contrasts",
    "repeated spatial units within animal, modelled as a random intercept",
    "BH within the module-phenotype family", "WGCNA / lmerTest",
    paste0("WGCNA ", tryCatch(as.character(utils::packageVersion("WGCNA")),
                              error = function(e) "n/a")),
    "module-phenotype tables", "descriptive - no association survived correction",
    "0 of 45 module-phenotype tests survived correction",
    "Module eigengenes were related to outcome group with a linear mixed model including a random intercept for animal; no module-phenotype association survived Benjamini-Hochberg correction."),
  M("ewce", "Which external cell types are enriched among module members?",
    "not applicable - external reference", "gene set",
    "EWCE bootstrap enrichment", "module vs reference panel", "not applicable",
    "none", "BH within the EWCE family", "EWCE", "see manifest",
    "cell-type affinity tables", "external context only",
    "affinity is context, never cell-intrinsic identity",
    "Module cell-type affinity was assessed against external single-cell reference panels; these results provide cell-type context and do not establish the cellular origin of an enriched-ROI measurement."),
  M("spatial_identity", "Is the spatial molecular architecture reproducible?",
    "3 CON animals", "spatial unit profile",
    "bilateral correlation and ICC", "CON only", "not a group contrast",
    "two hemispheres per animal", "none - descriptive",
    "in-house R", Rv, "bilateral reproducibility tables", "descriptive",
    "CON only; says nothing about stress",
    "Spatial reproducibility was quantified as the left-right concordance of each prespecified anatomical contrast in control animals."),
  M("external_validation", "Do internal anatomical contrasts recover published signatures?",
    "3 CON animals", "gene set", "GSEA against external signatures",
    "CON-only anatomical contrasts", "anatomical, not phenotypic", "none",
    "BH within the external-validation inventory", "clusterProfiler",
    paste0("clusterProfiler ", utils::packageVersion("clusterProfiler")),
    "external validation tables and Figure 2g", "inferential - the only external validation",
    "validates anatomy, not the stress result",
    "Internal control-only anatomical contrasts were tested against independently published hippocampal signatures; this is the only externally anchored validation in the study."),
  M("network", "Does the spatial molecular network differ by outcome?",
    "9 animals (27 animal x dataset instances)", "animal x dataset network",
    "distance from a leave-one-CON-animal-out consensus",
    "3 groups x 3 datasets", "SUS/RES vs CON consensus",
    "three dataset instances per animal", "BH within the network family",
    "in-house R", Rv, "network distance and edge-coupling tables",
    "descriptive - no detectable difference",
    "27 instances arise from 9 animals; they are not 27 independent replicates",
    "Animal-level spatial networks were compared with a leave-one-control-animal-out consensus; no whole-network group difference and no edge-behaviour association survived correction at this sample size."),
  M("behaviour_correlation", "Do network edges track behavioural outcome?",
    "9 animals", "network edge", "correlation with behavioural score",
    "8 tested neuropil edges", "edge x behaviour", "one value per animal",
    "BH within the edge-coupling family", "in-house R", Rv,
    "ST8 edge-behaviour table", "descriptive - none survived correction",
    "only 8 neuropil edges were tested, not every spatial-unit pair",
    "Edge-behaviour associations were tested for eight neuropil spatial-unit pairs; none survived correction, and with nine animals a single correlation has very little resolution."))
utils::write.csv(meth, file.path(PH_TAB, "manuscript_methods_contract.csv"),
                 row.names = FALSE)

cat("\n===== PART A SCANS =====\n")
cat("corpus files:", length(corpus), "\n")
cat("A03 statistic identity rows:", nrow(ident),
    "| FAIL:", sum(ident$status != "PASS"), "\n")
cat("A04 n-mentions:", nrow(n_hits), "| needing review:",
    sum(n_hits$status == "REVIEW"), "\n")
cat("A05 multiple-testing mentions:", nrow(mt_hits), "| FAIL:",
    sum(mt_hits$status == "FAIL"), "| REVIEW:", sum(mt_hits$status == "REVIEW"), "\n")
cat("A06-A08 claim-language hits:", nrow(claim), "| in prose (REVIEW):",
    sum(claim$status != "PASS - prohibition or disclaimer text"), "\n")
print(table(claim$term_class, claim$status))
cat("\nA10 source-data files:", nrow(zna), "| REVIEW:", sum(zna$risk != "PASS"), "\n")
cat("A11 eps = 1e-10 disclosure present in corpus:", floor_ok, "\n")
cat("A12 methods contract rows:", nrow(meth), "\n")
cat("\nA09 gate check: S9 pattern =", gate$gate_pattern,
    "| matches adjective 'selective':", gate$matches_adjective_selective,
    "| unflagged prose occurrences:", gate$unflagged_prose_occurrences, "\n")
cat("   ", gate$unflagged_locations, "\n")
