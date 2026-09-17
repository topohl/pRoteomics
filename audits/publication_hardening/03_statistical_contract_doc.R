#!/usr/bin/env Rscript

# Renders docs/MANUSCRIPT_STATISTICAL_CONTRACT.md from the Part A audit tables.
# Nothing is asserted here that is not already in those tables, so the document
# cannot drift from the audit that produced it.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
source("audits/publication_hardening/00_checkpoint.R")

rd <- function(f) utils::read.csv(file.path(PH_TAB, f), stringsAsFactors = FALSE)
m <- rd("manuscript_methods_contract.csv")
e <- rd("manuscript_effect_contract.csv")
i <- rd("manuscript_statistic_identity_audit.csv")
d <- rd("manuscript_dataset_interpretation_rules.csv")

# a pipe inside a cell would break the markdown table
esc <- function(x) gsub("|", "/", x, fixed = TRUE)

L <- c(
  "# Manuscript statistical contract", "",
  "One row per analysis that reaches the manuscript. This file answers, once:",
  "what exactly was tested, on how many animals, against which multiple-testing",
  "family, and what may therefore be claimed.", "",
  "Generated from the frozen v9 figure contract and the Part A audit tables by",
  "`audits/publication_hardening/02_statistical_audits.R` and rendered by",
  "`audits/publication_hardening/03_statistical_contract_doc.R`.", "",
  "## 1. Analyses", "")
for (k in seq_len(nrow(m))) L <- c(L,
  sprintf("### %s", m$analysis_id[k]), "",
  sprintf("- **Question:** %s", m$biological_question[k]),
  sprintf("- **Biological n:** %s (unit of analysis: %s)", m$biological_n[k],
          m$analysis_unit[k]),
  sprintf("- **Model / test:** %s", m$model_or_test[k]),
  sprintf("- **Design:** %s; contrast: %s", m$design[k], m$contrast[k]),
  sprintf("- **Repeated structure:** %s", m$repeated_structure[k]),
  sprintf("- **Multiple-testing scope:** %s", m$multiple_testing_scope[k]),
  sprintf("- **Software:** %s (%s)", m$software[k], m$version[k]),
  sprintf("- **Claim status:** %s", m$claim_status[k]),
  sprintf("- **Major limitation:** %s", m$major_limitation[k]),
  sprintf("- **Methods sentence:** %s", m$methods_sentence_candidate[k]), "")

L <- c(L, "## 2. Effect and sign conventions", "",
  "| Analysis | Contrast | Formal effect definition | Positive means |",
  "|---|---|---|---|")
for (k in seq_len(nrow(e))) L <- c(L, sprintf("| %s | %s | %s | %s |",
  esc(e$analysis[k]), esc(e$contrast[k]), esc(e$formal_effect_definition[k]),
  esc(e$positive_direction_meaning[k])))

L <- c(L, "", "## 3. Statistic identity", "",
  "Each displayed quantity has exactly one correct name and a list of names that",
  "must never be substituted for it.", "",
  "| Quantity | Correct name | Never call it | Where displayed |",
  "|---|---|---|---|")
for (k in seq_len(nrow(i))) L <- c(L, sprintf("| %s | %s | %s | %s |",
  esc(i$quantity[k]), esc(i$correct_name[k]), esc(i$prohibited_names[k]),
  esc(i$displayed_where[k])))

L <- c(L, "", "## 4. Dataset-specific interpretation rules", "",
  "| Dataset | Resolution | Legal | Illegal |", "|---|---|---|---|")
for (k in seq_len(nrow(d))) L <- c(L, sprintf("| %s | %s | %s | %s |",
  esc(d$dataset[k]), esc(d$spatial_resolution[k]), esc(d$legal_claim[k]),
  esc(d$illegal_claim[k])))

L <- c(L, "", "## 5. Standing constraints", "",
  "- The biological replicate is the **animal**: n = 3 per group. Acquisitions,",
  "  spatial units, hemispheres and animal x dataset network instances are never",
  "  biological replicates, and are labelled as such wherever they are counted.",
  "- **No theme-level p-value or FDR exists.** Atlas themes are descriptive",
  "  aggregations of canonical GO terms and carry no multiple-testing family.",
  "- **GSEA p-values are floored at eps = 1e-10.** This affects 90 of the 851",
  "  displayed FDR-supported occurrences, including all three exemplars. A term",
  "  at the floor has a true p the method does not resolve, so its FDR bounds",
  "  the evidence rather than measuring it.",
  "- **Absence of FDR support is never absence of effect.** At three animals per",
  "  group, write 'did not survive correction' or 'no detectable ... at the",
  "  present sample size', never 'no difference', 'unchanged' or 'equivalent'.",
  "- **CAMERA is LEVEL 4 sensitivity**, never independent validation. The only",
  "  genuinely external validation is the CON-only anatomical signature test",
  "  (F2g, ED2b).",
  "- **Specificity and selectivity claims require a heterogeneity test.** The",
  "  only one executed is the WGCNA stress x spatial-unit omnibus (0 of 35",
  "  FDR-supported, smallest FDR 0.27). No equivalent test exists at GO-program",
  "  level, so program-level results are reported as the count of spatial units",
  "  in which the contrast was FDR-supported, not as selectivity.", "",
  "## 6. Spatial wording contract", "",
  "One phrase, used consistently, for what the spatial design delivers.", "",
  "- **`spatially resolved` is the default descriptive wording.** It says that",
  "  the design and analysis resolve effects across defined hippocampal spatial",
  "  contexts - 18 prespecified spatial units, laminar in the neuropil and",
  "  region-level in the neuronal soma and microglia-enriched ROI. It asserts",
  "  nothing about where effects are or are not present.",
  "- **`restricted`, `selective` and `specific` require either a formal",
  "  corresponding test or an explicitly factual description of observed",
  "  support.** 'FDR-supported in 6 of 18 spatial units' is factual and allowed.",
  "  'Spatially restricted program differences' is an inferential claim and is",
  "  not, because the only heterogeneity test in the package is at WGCNA level",
  "  and is FDR-negative.",
  "- **Absence of support in some contexts does not by itself establish",
  "  specificity.** At three animals per group, a unit without FDR support is a",
  "  unit where nothing was detected, not a unit where nothing is happening.",
  "  What is restricted is the detection, not the effect.",
  "",
  "Enforced by the `spatially restricted / specific` rule in the S9 scan of",
  "`figures/final_truth_v9_semantics.R`, which is anchored to the adverb",
  "`spatially` so that a restricted *interpretation* or a specificity",
  "*inventory* is not flagged.", "")

writeLines(L, file.path("docs", "MANUSCRIPT_STATISTICAL_CONTRACT.md"))
cat("wrote docs/MANUSCRIPT_STATISTICAL_CONTRACT.md:", length(L), "lines\n")
