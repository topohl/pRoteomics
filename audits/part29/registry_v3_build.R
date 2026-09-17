#!/usr/bin/env Rscript

# Finalization sections 13-17: build and validate the candidate manuscript
# GO-theme registry v3.
#
# Two changes, both phenotype-independent and both expressed as ontology rules:
#
#  1. The mitochondrial theme gains ONE exclusion row: the GO:0006096
#     glycolytic-process sub-DAG. The glycolytic terms entered only because
#     GO:0045333 cellular respiration reaches them through approved edges. This
#     is a single ontology rule, not a hand list of individual terms, and it was
#     chosen without reference to any NES or FDR.
#
#  2. A seventh primary theme, neuron projection development, is defined from
#     ontology anchors rather than by enumerating the supported GO IDs Part 29
#     happened to find.
#
# This script only WRITES AND VALIDATES the candidate. Promotion is a separate
# step and happens only if every validation check passes.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({ library(GO.db); library(AnnotationDbi); library(digest) })
source(file.path("R", "paths.R"))
source(repo_path("R", "manuscript_go_theme_utils.R"))

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(AUD, recursive = TRUE, showWarnings = FALSE)

REG <- file.path("config", "manuscript_go_theme_registry.tsv")
CAND <- file.path("config", "manuscript_go_theme_registry_v3_candidate.tsv")
GLYCOLYSIS <- "GO:0006096"
NEW_THEME <- "neuron_projection_development"

old <- utils::read.delim(REG, stringsAsFactors = FALSE, check.names = FALSE,
                         quote = "")
stopifnot(all(old$registry_version == "manuscript_go_themes_v2"))

R <- function(theme_id, display_label, anchor_go_id, anchor_label, theme_role,
              display_order, match_scope, rationale)
  data.frame(theme_id = theme_id, display_label = display_label,
             anchor_go_id = anchor_go_id, anchor_label = anchor_label,
             theme_role = theme_role, display_order = display_order,
             match_scope = match_scope, rationale = rationale,
             registry_version = "manuscript_go_themes_v3",
             stringsAsFactors = FALSE)

new <- old
new$registry_version <- "manuscript_go_themes_v3"

# ---- 1. the single glycolysis exclusion ----------------------------------
new <- rbind(new, R(
  "mitochondrial_respiration_oxphos", "Mitochondrial respiration / OXPHOS",
  GLYCOLYSIS, "glycolytic process", "primary", 4L,
  "exclude_anchor_and_descendants",
  paste0("Cytosolic glycolysis reaches this theme only because cellular ",
         "respiration (GO:0045333) is an approved ancestor of it. The theme ",
         "is mitochondrial respiratory bioenergetics, so the glycolytic ",
         "sub-DAG is excluded by one ontology rule. Defined independently of ",
         "any enrichment statistic; pyruvate decarboxylation and the TCA ",
         "cycle are retained because neither descends from GO:0006096.")))

# ---- 2. the seventh theme -------------------------------------------------
# Anchors follow the registry's own convention: a structural anchor plus an
# explicit regulatory anchor, because regulation edges are not traversed from a
# process anchor. Neuron migration is a sibling branch, not a descendant of
# neuron projection development, so it needs its own anchor.
NP_ORDER <- 6L
new <- rbind(new, R(
  NEW_THEME, "Neuron projection development", "GO:0031175",
  "neuron projection development", "primary", NP_ORDER,
  "anchor_and_descendants",
  paste0("Structural anchor covering axonogenesis, axon guidance, dendrite ",
         "and neurite development; descendants are followed only through ",
         "is_a and part_of.")))
new <- rbind(new, R(
  NEW_THEME, "Neuron projection development", "GO:0010975",
  "regulation of neuron projection development", "primary", NP_ORDER,
  "anchor_and_descendants",
  "Explicit regulatory anchor; regulation edges are not traversed from the structural anchor."))
new <- rbind(new, R(
  NEW_THEME, "Neuron projection development", "GO:0001764", "neuron migration",
  "primary", NP_ORDER, "anchor_and_descendants",
  paste0("Neuron migration is a sibling branch of projection development, not ",
         "a descendant of it, and was part of the same phenotype-blind ",
         "semantic cluster; it needs its own anchor.")))

# Autophagy moves from 6 to 7 so the new row sits between synaptic signalling
# and endolysosomal trafficking. Order is fixed and phenotype-independent.
new$display_order[new$theme_id == "autophagy_lysosome_endosome"] <- 7L

new <- new[order(new$display_order, new$theme_id, new$anchor_go_id), ,
           drop = FALSE]
utils::write.table(new, CAND, sep = "\t", row.names = FALSE, quote = FALSE)
cat("candidate registry written:", CAND, "|", nrow(new), "rows\n")

# ========================================================== validation
TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)
tested <- unique(TH[, c("GO_ID", "GO_description")])
sup <- TH[is.finite(TH$GSEA_FDR) & TH$GSEA_FDR < 0.05, , drop = FALSE]

reg_v3 <- read_manuscript_go_theme_registry(CAND)
desc_of_go <- stats::setNames(tested$GO_description, tested$GO_ID)
names(tested) <- c("ID", "Description")
asg <- map_go_terms_to_manuscript_themes(tested, reg_v3)$assignments
cat("assignments:", nrow(asg), "| themes:",
    length(unique(asg$theme_id[asg$theme_role == "primary"])), "\n")

old_members <- unique(TH[TH$theme_claim_eligible %in% TRUE,
                         c("theme_id", "GO_ID")])
new_members <- unique(asg[asg$theme_role == "primary", c("theme_id", "GO_ID")])

# --- mitochondrial: the exclusion did exactly one thing
mito_old <- sort(old_members$GO_ID[old_members$theme_id == "mitochondrial_respiration_oxphos"])
mito_new <- sort(new_members$GO_ID[new_members$theme_id == "mitochondrial_respiration_oxphos"])
removed <- setdiff(mito_old, mito_new); added_m <- setdiff(mito_new, mito_old)

# --- neuron projection: what the anchors admit beyond the Part-29 cluster
np_new <- sort(new_members$GO_ID[new_members$theme_id == NEW_THEME])
cl_terms <- NULL
bnd <- file.path(AUD, "neuron_projection_theme_boundary.csv")
np_cluster_n <- 42L
np_supported <- sup$GO_ID[sup$GO_ID %in% np_new]
overlap_other <- new_members$GO_ID[new_members$theme_id != NEW_THEME]
np_overlap <- intersect(np_new, overlap_other)

desc_of <- function(id) {
  d <- go_bp_allowed_descendants(id)
  if (is.null(d) || !nrow(d)) character(0) else d$descendant_GO_ID
}
np_all_ontology <- unique(unlist(lapply(
  c("GO:0031175", "GO:0010975", "GO:0001764"), desc_of)))

val <- data.frame(
  check = c(
    "registry parses under the v3 reader",
    "exactly one exclusion rule, on the glycolysis sub-DAG",
    "mitochondrial terms removed",
    "mitochondrial terms added",
    "PDH retained", "TCA retained",
    "no other theme's membership changed",
    "neuron projection theme size (tested GO terms)",
    "neuron projection terms overlapping another primary theme",
    "neuron projection supported occurrences",
    "unintended admission: theme size against its ontology closure",
    "primary themes in v3"),
  value = c(
    "PASS",
    as.character(sum(new$match_scope == "exclude_anchor_and_descendants")),
    paste(sort(unname(desc_of_go[removed])), collapse = "; "),
    if (length(added_m)) paste(added_m, collapse = "; ") else "none",
    as.character("GO:0006086" %in% mito_new),
    as.character("GO:0006099" %in% mito_new),
    "verified below",
    as.character(length(np_new)),
    if (length(np_overlap)) paste(np_overlap, collapse = "; ") else "none",
    as.character(sum(sup$GO_ID %in% np_new)),
    sprintf("%d tested terms of %d in the full ontology closure",
            length(np_new), length(np_all_ontology)),
    as.character(length(unique(new$theme_id[new$theme_role == "primary"])))),
  stringsAsFactors = FALSE)

others <- setdiff(unique(old_members$theme_id), "mitochondrial_respiration_oxphos")
unchanged <- vapply(others, function(t)
  identical(sort(old_members$GO_ID[old_members$theme_id == t]),
            sort(new_members$GO_ID[new_members$theme_id == t])), logical(1))
val$value[val$check == "no other theme's membership changed"] <-
  sprintf("%d of %d unchanged: %s", sum(unchanged), length(others),
          paste(names(unchanged)[!unchanged], collapse = "; "))

CL <- utils::read.csv(file.path(AUD, "atlas_omitted_semantic_clusters.csv"),
                      stringsAsFactors = FALSE)
np_cluster_terms <- unique(sup$GO_ID[!sup$GO_ID %in%
  unique(TH$GO_ID[TH$theme_claim_eligible %in% TRUE])])
# the 42 terms of OMIT_B12, recovered from the Part-29 boundary file
BND <- utils::read.csv(file.path(AUD, "neuron_projection_theme_boundary.csv"),
                       stringsAsFactors = FALSE)
np_terms_listed <- trimws(unlist(strsplit(
  BND$neural_terms_in_cluster[BND$cluster_id == "OMIT_B12"], ";")))
captured <- sum(np_terms_listed %in% unname(desc_of_go[np_new]))
val <- rbind(val, data.frame(
  check = c("cluster terms listed in the Part-29 boundary file",
            "of those captured by the ontology anchors",
            "terms admitted by the anchors but never tested"),
  value = c(as.character(length(np_terms_listed)),
            as.character(captured),
            as.character(length(setdiff(np_all_ontology, tested$ID)))),
  stringsAsFactors = FALSE))

utils::write.csv(val, file.path(AUD,
  "manuscript_go_theme_registry_v3_validation.csv"), row.names = FALSE)

cat("\n===== REGISTRY v3 VALIDATION =====\n")
print(val, row.names = FALSE)
cat("\nold registry sha256:", digest::digest(file = REG, algo = "sha256"), "\n")
cat("candidate sha256   :", digest::digest(file = CAND, algo = "sha256"), "\n")
saveRDS(list(np_new = np_new, removed = removed, mito_new = mito_new),
        file.path(AUD, "registry_v3_sets.rds"))
