#!/usr/bin/env Rscript

# Finalization sections 13-16: promote the validated candidate registry to
# canonical, and record the provenance of the change.
#
# This changes a REGISTRY, not a result. No differential abundance, no GSEA and
# no CAMERA is rerun; only the mapping of already-computed GO terms to display
# themes changes.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages(library(digest))

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
REP <- file.path("results", "reports", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(REP, recursive = TRUE, showWarnings = FALSE)

REG <- file.path("config", "manuscript_go_theme_registry.tsv")
CAND <- file.path("config", "manuscript_go_theme_registry_v3_candidate.tsv")
val <- utils::read.csv(file.path(AUD,
  "manuscript_go_theme_registry_v3_validation.csv"), stringsAsFactors = FALSE)

get <- function(k) val$value[val$check == k][1]
stopifnot(identical(get("registry parses under the v3 reader"), "PASS"))
stopifnot(identical(get("exactly one exclusion rule, on the glycolysis sub-DAG"), "1"))
stopifnot(identical(get("PDH retained"), "TRUE"),
          identical(get("TCA retained"), "TRUE"))
stopifnot(identical(get("mitochondrial terms added"), "none"))
stopifnot(identical(get("neuron projection terms overlapping another primary theme"),
                    "none"))
stopifnot(identical(get("primary themes in v3"), "7"))

old_hash <- digest::digest(file = REG, algo = "sha256")
new_hash <- digest::digest(file = CAND, algo = "sha256")
old <- utils::read.delim(REG, stringsAsFactors = FALSE, quote = "")
new <- utils::read.delim(CAND, stringsAsFactors = FALSE, quote = "")

file.copy(CAND, REG, overwrite = TRUE)
stopifnot(identical(digest::digest(file = REG, algo = "sha256"), new_hash))

md <- c(
"# Manuscript GO-theme registry v3 - promotion provenance", "",
sprintf("- previous version: %s", unique(old$registry_version)),
sprintf("- new version: %s", unique(new$registry_version)),
sprintf("- previous file sha256: `%s`", old_hash),
sprintf("- new file sha256: `%s`", new_hash),
sprintf("- rows: %d -> %d", nrow(old), nrow(new)),
sprintf("- primary themes: %d -> %d",
        length(unique(old$theme_id[old$theme_role == "primary"])),
        length(unique(new$theme_id[new$theme_role == "primary"]))),
"",
"## Change 1 - the mitochondrial theme excludes the glycolysis sub-DAG", "",
"The mitochondrial respiration / OXPHOS theme excludes the glycolytic process",
"sub-DAG to distinguish cytosolic glycolysis from mitochondrial respiratory",
"bioenergetics; this ontology rule was defined independently of phenotype",
"statistics.",
"",
"Cytosolic glycolysis had entered the theme only because GO:0045333 cellular",
"respiration is an approved ancestor of it under is_a / part_of. The fix is ONE",
"rule - a new `exclude_anchor_and_descendants` row on GO:0006096 - and not a",
"list of individual terms. Four terms leave the theme:",
"",
"- glycolytic process", "- canonical glycolysis",
"- glycolytic process through fructose-6-phosphate",
"- glycolytic process through glucose-6-phosphate",
"",
"Pyruvate decarboxylation to acetyl-CoA (GO:0006086) and the tricarboxylic acid",
"cycle (GO:0006099) are RETAINED, because neither descends from GO:0006096.",
"A hand-built exclusion list could easily have removed them by mistake.",
"",
"Glycolysis was NOT removed for being non-significant. The rule is ontological",
"and was fixed before any enrichment value was consulted; the measured effect on",
"the atlas was one sign change in 54 cells, no change to any support dot, and",
"one supported occurrence.",
"",
"## Change 2 - a seventh primary theme", "",
"Neuron projection development is added, defined by ontology anchors rather than",
"by enumerating the GO IDs that happened to be supported:",
"",
"| anchor | label | scope |",
"|---|---|---|",
"| GO:0031175 | neuron projection development | anchor_and_descendants |",
"| GO:0010975 | regulation of neuron projection development | anchor_and_descendants |",
"| GO:0001764 | neuron migration | anchor_and_descendants |",
"",
"The regulatory anchor is explicit because regulation edges are not traversed",
"from a structural anchor, and neuron migration is a sibling branch rather than",
"a descendant. This follows the convention the RNA-processing theme already uses.",
"",
sprintf("The definition admits %s tested GO terms, of which %s overlap another",
        get("neuron projection theme size (tested GO terms)"),
        get("neuron projection terms overlapping another primary theme")),
sprintf("primary theme, carrying %s FDR-supported occurrences. Its full ontology",
        get("neuron projection supported occurrences")),
sprintf("closure is %s, so most of the branch was simply never tested in these",
        sub(".* of ", "", get("unintended admission: theme size against its ontology closure"))),
"proteomes - the anchors do not drag in a large unrelated branch.",
"",
"## Change 3 - display order", "",
"Neuron projection development takes position 6 and autophagy / endolysosomal",
"trafficking moves from 6 to 7, so the new row sits between synaptic signalling",
"and endolysosomal trafficking. The order is fixed and phenotype-independent: it",
"follows a molecular-to-cellular progression (RNA, translation, chromatin,",
"bioenergetics, synaptic signalling, neuronal morphogenesis, degradation) and is",
"never sorted by NES, FDR, direction or number of supported cells.",
"",
"## What did not change", "",
"- No differential abundance, GSEA or CAMERA result was recomputed.",
"- The other five primary themes keep their membership exactly.",
"- The QC-review and supporting themes are untouched and remain in source data.",
"- The three Figure-3 exemplar terms are unchanged.")
writeLines(md, file.path(REP, "manuscript_go_theme_registry_v3_provenance.md"))

cat("promoted", basename(CAND), "->", basename(REG), "\n")
cat("  old sha256:", old_hash, "\n  new sha256:", new_hash, "\n")
cat("  primary themes:", length(unique(new$theme_id[new$theme_role == "primary"])), "\n")
