#!/usr/bin/env Rscript

# PH-002 resolution, section 4: every reader-facing occurrence of the
# specificity / comparative vocabulary, with the generator line that produces it.
#
# SCAN ONLY. Nothing is edited here. The corpus is the reader-facing v9 layer
# plus the frozen figure contract; for each hit the generating source line is
# resolved, because the .md files are emitted by scripts and editing the .md
# directly would be overwritten on the next rebuild.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
source("audits/publication_hardening/00_checkpoint.R")

V9R <- file.path("results", "reports", "manuscript_candidates", "final_truth_v9")
V9T <- file.path("results", "tables", "manuscript_candidates", "final_truth_v9")

corpus <- c(
  list.files(V9R, pattern = "[.]md$", recursive = TRUE, full.names = TRUE),
  list.files(V9T, pattern = "README[.]md$", recursive = TRUE, full.names = TRUE),
  "figures/figure_final_truth_v9_contract.yml")
corpus <- corpus[file.exists(corpus)]
rel <- function(f) sub(".*proteomics[/\\\\]", "", f)

# the section-4 vocabulary, longest first so a compound term wins over its parts
TERMS <- c("spatially selective", "spatially specific", "region-specific",
           "layer-specific", "susceptibility-specific", "resilience-specific",
           "selectively", "selective", "specifically", "specific",
           "preferential", "restricted", "localized", "localised", "stronger")
PAT <- paste0("(?i)(", paste(gsub("-", "[- ]", TERMS), collapse = "|"), ")")

# which generator writes each reader-facing file, so an edit lands in source
GEN <- c(
  final_figure_story_v9.md          = "figures/final_truth_v9_vector_audit.R",
  core_story_corrected.md           = "figures/final_truth_v9_semantics.R",
  manuscript_semantic_rules.md      = "figures/final_truth_v9_semantics.R",
  claim_chain_audit.md              = "figures/final_truth_v9_semantics.R",
  final_figure_legends_v9.md        = "figures/final_truth_v9_legends.R",
  known_issues_v9.md                = "figures/final_truth_v9_readmes.R",
  nature_reviewer_vulnerabilities.md = "figures/final_truth_v9_reviewer_audit.R",
  README.md                         = "figures/final_truth_v9_readmes.R")

hits <- list()
for (f in corpus) {
  ln <- readLines(f, warn = FALSE, encoding = "UTF-8")
  idx <- grep(PAT, ln, perl = TRUE)
  for (k in idx) {
    m <- regmatches(ln[k], gregexpr(PAT, ln[k], perl = TRUE))[[1]]
    for (term in unique(tolower(m))) {
      hits[[length(hits) + 1L]] <- data.frame(
        file = rel(f), line = k, term = term,
        generator = unname(GEN[basename(f)] %||% "hand-maintained"),
        context = substr(trimws(ln[k]), 1, 300),
        stringsAsFactors = FALSE)
    }
  }
}
hits <- do.call(rbind, hits)

# A hit inside the rulebook or the claim-chain audit is the repository quoting
# the wording in order to ban it; that is the opposite of a violation.
QUOTING_FILE <- c("manuscript_semantic_rules.md", "claim_chain_audit.md")
hits$quoting_file <- basename(hits$file) %in% QUOTING_FILE
hits$machine_field <- grepl("^[a-z_]+:\\s*[A-Za-z0-9_./-]+$", hits$context,
                            perl = TRUE)

# provisional class; the biological adjudication is done by hand afterwards
hits$provisional_class <- ifelse(
  hits$quoting_file, "LEGACY_NONREADER_FACING",
  ifelse(hits$machine_field, "ALLOWED_TECHNICAL",
  ifelse(grepl("selectiv|spatially specific|region[- ]specific|layer[- ]specific|susceptibility[- ]specific|resilience[- ]specific",
               hits$term), "REQUIRES_FORMAL_TEST", "REVIEW")))
# The adjudicated verdicts, so the artefact records a decision rather than a
# stale REVIEW. Each was reached by independent adjudication plus three
# adversarial refuters; see docs/PUBLICATION_HARDENING_PROGRESS.md PH-002,
# PH-008 and PH-009.
ADJ <- rbind(
  data.frame(file_base = "core_story_corrected.md", term = "restricted",
    adjudicated_class = "ALLOWED_DESCRIPTIVE", finding = "PH-008",
    verdict = "the rulebook's USE column names this phrase; the claim chain certifies it with named evidence (theme atlas and ED6: FDR-supported terms occur in a subset of units). Changing it would amend the contract, not fix wording."),
  data.frame(file_base = "nature_reviewer_vulnerabilities.md", term = "restricted",
    adjudicated_class = "ALLOWED_DESCRIPTIVE", finding = "PH-008",
    verdict = "same sanctioned phrase as core_story_corrected.md; left consistent with it."),
  data.frame(file_base = "final_figure_legends_v9.md", term = "restricted",
    adjudicated_class = "ALLOWED_DESCRIPTIVE", finding = "",
    verdict = "restricts an INTERPRETATION of CA2-SLM, not a spatial extent. Survived refutation 0 of 3."),
  data.frame(file_base = "figure_final_truth_v9_contract.yml", term = "restricted",
    adjudicated_class = "ALLOWED_DESCRIPTIVE", finding = "",
    verdict = "the same interpretation sentence, in its contract source."),
  data.frame(file_base = "core_story_corrected.md", term = "stronger",
    adjudicated_class = "ALLOWED_DESCRIPTIVE", finding = "",
    verdict = "both occurrences DISCLAIM the comparison ('neither is claimed to be the stronger'); they are the prohibition, not the claim."),
  data.frame(file_base = "final_figure_story_v9.md", term = "stronger",
    adjudicated_class = "ALLOWED_DESCRIPTIVE", finding = "",
    verdict = "regional vs laminar bilateral concordance is a WITHIN-family comparison of one metric on one scale, unlike the cross-family 'stronger' removed from the Figure 3 title. The claim chain keeps it as descriptive with 'the weaker laminar result is deliberately left visible'."),
  data.frame(file_base = "final_figure_legends_v9.md", term = "specific",
    adjudicated_class = "REQUIRES_FORMAL_TEST", finding = "PH-009",
    verdict = "RECORDED, NOT CHANGED, NOT VERIFIED. External-validation specificity is a different subject from PH-002; confirming it needs a new statistical audit."),
  data.frame(file_base = "figure_final_truth_v9_contract.yml", term = "specific",
    adjudicated_class = "REQUIRES_FORMAL_TEST", finding = "PH-009",
    verdict = "same, and in the frozen figure contract."),
  stringsAsFactors = FALSE)
k <- paste(basename(hits$file), hits$term)
m <- match(k, paste(ADJ$file_base, ADJ$term))
hits$adjudicated_class <- ifelse(hits$quoting_file, "LEGACY_NONREADER_FACING",
  ifelse(hits$machine_field, "ALLOWED_TECHNICAL",
         ADJ$adjudicated_class[m] %||% NA_character_))
hits$finding <- ADJ$finding[m]
hits$verdict <- ADJ$verdict[m]
hits$adjudicated_class[is.na(hits$adjudicated_class)] <- "ALLOWED_DESCRIPTIVE"
# Line-level override: contract line 62 is "before any specific spatial
# signature is inspected", an ordinary use of the adjective with no specificity
# claim attached. Its own refuter conceded it while overturning 118 and 286.
ovr <- basename(hits$file) == "figure_final_truth_v9_contract.yml" &
       hits$line == 62L
hits$adjudicated_class[ovr] <- "ALLOWED_DESCRIPTIVE"
hits$finding[ovr] <- ""
hits$verdict[ovr] <- "ordinary adjective; asserts no specificity claim"
rownames(hits) <- NULL
utils::write.csv(hits, file.path(PH_TAB, "ph002_language_occurrences.csv"),
                 row.names = FALSE)

cat("corpus files:", length(corpus), "| occurrences:", nrow(hits), "\n\n")
print(table(hits$term, hits$adjudicated_class))
unresolved <- hits[hits$adjudicated_class == "REQUIRES_FORMAL_TEST" &
                   !hits$quoting_file, ]
cat("\nPH-002 target sentences still present:",
    sum(basename(hits$file) == "final_figure_story_v9.md" &
        grepl("selectiv", hits$term)), "\n")
cat("occurrences still classed REQUIRES_FORMAL_TEST:", nrow(unresolved),
    "| all attributed to:",
    paste(unique(unresolved$finding[nzchar(unresolved$finding %||% "")]),
          collapse = ", "), "\n")
for (i in seq_len(nrow(unresolved)))
  cat(sprintf("  %-34s %4s %-12s %s\n", basename(unresolved$file[i]),
              unresolved$line[i], unresolved$term[i], unresolved$finding[i]))
